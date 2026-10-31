#include "../galerkin/regular_qualification.hpp"
#include "device_gmres.hpp"
#include "support_fmm_solve.hpp"
#include "workspace_size.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cublas_v2.h>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
namespace hornsculpt::wideband {
namespace {
using cuda_detail::checkedSizeProduct;
using C = std::complex<double>;
using Z = cuDoubleComplex;
using Clock = std::chrono::steady_clock;
void check(cudaError_t code, const char *context) {
    if (code != cudaSuccess)
        throw std::runtime_error(std::string(context) + ": " + cudaGetErrorString(code));
}
void check(cublasStatus_t code, const char *context) {
    if (code != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error(std::string(context) + ": cuBLAS " + std::to_string(code));
}
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("BEM calculation cancelled.");
}
Z device(C value) {
    return make_cuDoubleComplex(value.real(), value.imag());
}
class Buffer {
  public:
    explicit Buffer(std::size_t count) : count_(count) {
        check(cudaMalloc(reinterpret_cast<void **>(&data_), checkedSizeProduct(count, sizeof(Z))),
              "GMRES allocation");
    }
    ~Buffer() {
        cudaFree(data_);
    }
    Buffer(const Buffer &) = delete;
    Buffer &operator=(const Buffer &) = delete;
    Z *data() {
        return data_;
    }
    std::size_t bytes() const {
        return count_ * sizeof(Z);
    }

  private:
    Z *data_{};
    std::size_t count_{};
};
__global__ void difference(const Z *right, const Z *left, Z *output, int count) {
    const int i = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (i < count)
        output[i] = cuCsub(right[i], left[i]);
}
double elapsed(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
std::string scientific(double value) {
    std::ostringstream output;
    output << std::scientific << std::setprecision(6) << value;
    return output.str();
}
BemNumerics innerKrylovOptions(BemNumerics value) {
    value.relativeResidualTolerance *= supportFmmKrylovSettingsScale;
    return value;
}
} // namespace

bool preferSupportDirectAudit(const SupportTree &tree,
                              std::span<const BemFastMultipoleLevelEvidence> levels,
                              unsigned regularPoints, unsigned singleToDoubleThroughputRatio) {
    const auto panels = tree.fluxLeaf.size();
    if (!panels || panels > maximumDirectAuditPanels || !regularPoints || levels.empty())
        return false;
    const long double direct =
        static_cast<long double>(panels) * panels *
        (9.L + static_cast<long double>(regularPoints) * regularPoints *
                   (13.L + 1.L / std::max(1u, singleToDoubleThroughputRatio)));
    long double multipole = 0;
    for (const auto &level : levels) {
        if (level.level >= tree.levels.size())
            throw std::invalid_argument("Audit cost estimate has an invalid level index.");
        const auto &inventory = tree.levels[level.level];
        const long double pairs =
            static_cast<long double>(inventory.pressurePressureFarInteractions) +
            static_cast<long double>(inventory.pressureFluxFarInteractions);
        const long double p = level.momentBandlimit, n = p + 1;
        if (level.representation == BemFastMultipoleRepresentation::regularSpherical) {
            const long double rotation = n * (4 * p * p + 8 * p + 3) / 3;
            const long double axial = n * n + p * n * (2 * p + 1) / 3;
            multipole += pairs * (2 * rotation + axial);
            // A polar projection of each distinct axial distance class is
            // O(p^4); cached directions do not multiply this preparation work.
            multipole +=
                static_cast<long double>(inventory.squaredOffsetClasses.size()) * n * n * n * n;
        } else {
            const long double q =
                2.L * (level.quadratureOrder + 1.L) * (level.quadratureOrder + 1.L);
            multipole += pairs * q;
        }
    }
    // One original-operator audit replaces two independently seeded plans.
    return direct <= 2 * multipole;
}

class SupportFmmLinearSolver::Implementation {
  public:
    Implementation(SupportFmmOperator &op, std::size_t size, const BemNumerics &settings,
                   std::stop_token token, SupportDirectAuditPolicy policy,
                   SupportFmmProgress callback, unsigned quadratureDegree)
        : sourceOperation(op), operation(&op), krylovOperation(&op), n(static_cast<int>(size)),
          options(settings), stream(op.stream()), stop(token), progress(std::move(callback)),
          krylov(
              [this, token](const Z *in, Z *out) {
                  graphAction(*krylovOperation,
                              [&] { krylovOperation->applyPressureDevice(in, out, token); });
                  captureOperatorBytes(*krylovOperation);
                  ++pressureApplications;
                  report(BemSolvePhase::krylovSolve);
              },
              [&op, token](Z *value) { op.preconditionDevice(value, token); }, size,
              innerKrylovOptions(settings), op.stream(), token,
              [this] {
                  ++completedIterations;
                  report(BemSolvePhase::krylovSolve);
              }),
          rhs(size), product(size), residual(size), correction(size), comparisonProduct(size),
          comparisonRhs(size), temporary(size), directPolicy(std::move(policy)),
          panelDegree(quadratureDegree), sourcePanelDegree(quadratureDegree) {
        if (panelDegree != 4 && panelDegree != 6 && panelDegree != 8 && panelDegree != 12)
            throw std::invalid_argument("FMM endpoint comparison needs degree 4, 6, 8 or 12.");
        const auto expectedRule = galerkin::triangleRule(panelDegree);
        const auto actualRule = op.panelRule();
        bool identical = expectedRule.size() == actualRule.size();
        for (std::size_t i = 0; identical && i < expectedRule.size(); ++i) {
            identical &= expectedRule[i].weight == actualRule[i].weight;
            for (unsigned j = 0; j < 3; ++j)
                identical &= expectedRule[i].barycentric[j] == actualRule[i].barycentric[j];
        }
        if (!identical)
            throw std::invalid_argument(
                "FMM endpoint degree does not describe its actual panel rule.");
        if (static_cast<bool>(directPolicy.makeDirect) !=
            static_cast<bool>(directPolicy.preferDirect))
            throw std::invalid_argument("Direct audit policy needs both a decision and a factory.");
        if (op.activeFaceIds().size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::invalid_argument("FMM active face count exceeds the cuBLAS vector limit.");
        check(cublasCreate(&blas), "FMM audit cuBLAS creation");
        try {
            check(cublasSetStream(blas, stream), "FMM audit cuBLAS stream");
        } catch (...) {
            cublasDestroy(blas);
            throw;
        }
        captureOperatorBytes(op);
    }
    ~Implementation() {
        cublasDestroy(blas);
    }
    void report(BemSolvePhase phase) {
        currentPhase = phase;
        if (progress)
            progress(phase, completedIterations, pressureApplications, fluxApplications);
    }
    void observeGraphResources(const SupportFmmWorkStatistics &stats) {
        graphEvidence.maximumCachedExecutables =
            std::max(graphEvidence.maximumCachedExecutables,
                     static_cast<std::uint64_t>(stats.m2lGraphMaximumCachedExecutables));
        graphEvidence.maximumCachedNodes =
            std::max(graphEvidence.maximumCachedNodes,
                     static_cast<std::uint64_t>(stats.m2lGraphMaximumCachedNodes));
        graphEvidence.peakOwnedHostMetadataBytes =
            std::max(graphEvidence.peakOwnedHostMetadataBytes,
                     static_cast<std::uint64_t>(stats.m2lGraphPeakOwnedHostPayloadBytes));
        if (stats.m2lGraphMemoryQueries) {
            const auto freeBytes =
                static_cast<std::uint64_t>(stats.m2lGraphMinimumObservedFreeBytes);
            if (!graphEvidence.memoryObservationAvailable ||
                freeBytes < graphEvidence.minimumObservedFreeDeviceBytes) {
                graphEvidence.minimumObservedFreeDeviceBytes = freeBytes;
                graphEvidence.totalDeviceBytes = stats.m2lGraphTotalDeviceBytes;
                graphEvidence.requiredFreeReserveBytes = stats.m2lGraphRequiredFreeReserveBytes;
            }
            graphEvidence.memoryObservationAvailable = true;
        }
    }
    template <class Operator, class Action> void graphAction(Operator &op, Action action) {
        if constexpr (std::is_same_v<Operator, SupportFmmOperator>) {
            // Values, never object addresses: serial audits may reuse an address.
            const auto before = op.workStatistics();
            action();
            const auto &after = op.workStatistics();
            if (after.m2lGraphCaptureAttempts < before.m2lGraphCaptureAttempts)
                throw std::logic_error("FMM graph statistics regressed across an action.");
            graphEvidence.captureAttempts +=
                after.m2lGraphCaptureAttempts - before.m2lGraphCaptureAttempts;
            if (after.m2lGraphInstantiations < before.m2lGraphInstantiations)
                throw std::logic_error("FMM graph statistics regressed across an action.");
            graphEvidence.instantiations +=
                after.m2lGraphInstantiations - before.m2lGraphInstantiations;
            if (after.m2lGraphLaunches < before.m2lGraphLaunches)
                throw std::logic_error("FMM graph statistics regressed across an action.");
            graphEvidence.launches += after.m2lGraphLaunches - before.m2lGraphLaunches;
            if (after.m2lGraphReuses < before.m2lGraphReuses)
                throw std::logic_error("FMM graph statistics regressed across an action.");
            graphEvidence.reuses += after.m2lGraphReuses - before.m2lGraphReuses;
            if (after.m2lGraphExecutedBatches < before.m2lGraphExecutedBatches)
                throw std::logic_error("FMM graph statistics regressed across an action.");
            graphEvidence.executedBatches +=
                after.m2lGraphExecutedBatches - before.m2lGraphExecutedBatches;
            if (after.m2lGraphDirectChunks < before.m2lGraphDirectChunks)
                throw std::logic_error("FMM graph statistics regressed across an action.");
            graphEvidence.canonicalChunks +=
                after.m2lGraphDirectChunks - before.m2lGraphDirectChunks;
            if (after.m2lGraphDirectBatches < before.m2lGraphDirectBatches)
                throw std::logic_error("FMM graph statistics regressed across an action.");
            graphEvidence.canonicalBatches +=
                after.m2lGraphDirectBatches - before.m2lGraphDirectBatches;
            if (after.m2lGraphInvalidatedExecutables < before.m2lGraphInvalidatedExecutables)
                throw std::logic_error("FMM graph statistics regressed across an action.");
            graphEvidence.invalidatedExecutables +=
                after.m2lGraphInvalidatedExecutables - before.m2lGraphInvalidatedExecutables;
            if (after.m2lGraphBudgetBypasses < before.m2lGraphBudgetBypasses)
                throw std::logic_error("FMM graph statistics regressed across an action.");
            graphEvidence.budgetBypasses +=
                after.m2lGraphBudgetBypasses - before.m2lGraphBudgetBypasses;
            if (after.m2lGraphAdmissionClosures < before.m2lGraphAdmissionClosures)
                throw std::logic_error("FMM graph statistics regressed across an action.");
            graphEvidence.admissionClosures +=
                after.m2lGraphAdmissionClosures - before.m2lGraphAdmissionClosures;
            graphEvidence.captureMilliseconds +=
                after.m2lGraphCaptureMilliseconds - before.m2lGraphCaptureMilliseconds;
            graphEvidence.instantiationMilliseconds +=
                after.m2lGraphInstantiationMilliseconds - before.m2lGraphInstantiationMilliseconds;
            observeGraphResources(after);
        } else {
            action();
        }
    }
    void releaseOperator(SupportFmmOperator &op) {
        graphAction(op, [&] { op.releaseDeviceState(); });
    }
    void captureOperatorBytes(const SupportFmmOperator &op) {
        stationaryOperatorBytes = std::max(stationaryOperatorBytes, op.stationaryDeviceBytes());
        peakOperatorBytes = std::max(peakOperatorBytes, op.peakDeviceBytes());
    }
    void captureOperatorBytes(const SupportDirectAudit &op) {
        // Direct buffers coexist with the retained family, but never with an
        // active FMM hierarchy. Include any actual remaining own buffers rather
        // than relying on that serial policy for the ownership measurement.
        const auto live = sourceOperation.sharedDeviceBytes() + sourceOperation.ownDeviceBytes() +
                          (refinedOperation ? refinedOperation->ownDeviceBytes() : 0) +
                          op.ownDeviceBytes();
        stationaryOperatorBytes = std::max(stationaryOperatorBytes, live);
        peakOperatorBytes = std::max(peakOperatorBytes, live);
    }
    template <class Operator> void applyPressure(Operator &op, const Z *in, Z *out) {
        report(BemSolvePhase::residualAudit);
        graphAction(op, [&] { op.applyPressureDevice(in, out, stop); });
        captureOperatorBytes(op);
        ++pressureApplications;
        report(BemSolvePhase::residualAudit);
    }
    template <class Operator> void applyFlux(Operator &op, const Z *in, Z *out) {
        graphAction(op, [&] { op.applyFaceFluxDevice(in, out, stop); });
        captureOperatorBytes(op);
        ++fluxApplications;
        report(currentPhase);
    }
    void synchronize() {
        check(cudaStreamSynchronize(stream), "GMRES synchronization");
    }
    void copy(const Z *source, Z *destination) {
        check(cudaMemcpyAsync(destination, source, static_cast<std::size_t>(n) * sizeof(Z),
                              cudaMemcpyDeviceToDevice, stream),
              "GMRES vector copy");
    }
    void clear(Z *value) {
        check(cudaMemsetAsync(value, 0, static_cast<std::size_t>(n) * sizeof(Z), stream),
              "GMRES vector clear");
    }
    void axpy(C coefficient, const Z *input, Z *output) {
        const Z alpha = device(coefficient);
        check(cublasZaxpy(blas, n, &alpha, input, 1, output, 1), "GMRES vector update");
    }
    double norm(const Z *value) {
        double result{};
        check(cublasDznrm2(blas, n, value, 1, &result), "GMRES vector norm");
        synchronize();
        if (!std::isfinite(result))
            throw std::runtime_error("GMRES encountered a non-finite vector.");
        return result;
    }
    void subtract(const Z *right, const Z *left, Z *output) {
        difference<<<(n + 255) / 256, 256, 0, stream>>>(right, left, output, n);
        check(cudaGetLastError(), "GMRES residual subtraction");
    }

    double ratio(double numerator, double denominator) const {
        return denominator > 0  ? numerator / denominator
               : numerator == 0 ? 0
                                : std::numeric_limits<double>::infinity();
    }
    double evaluate(SupportFmmOperator &audit, const Z *pressure, SupportFmmSolution &result) {
        const auto start = Clock::now();
        applyPressure(audit, pressure, product.data());
        subtract(rhs.data(), product.data(), residual.data());
        const double rightNorm = norm(rhs.data());
        result.boundaryResidual = ratio(norm(residual.data()), rightNorm);
        result.residualCheckMilliseconds += elapsed(start);
        captureOperatorBytes(audit);
        return rightNorm;
    }
    void improve(SupportFmmOperator &audit, Z *pressure, SupportFmmSolution &result) {
        krylovOperation = &audit;
        unsigned corrections = 0;
        while (result.boundaryResidual > options.relativeResidualTolerance) {
            cancelled(stop);
            if (corrections == options.maximumRefinementSteps)
                throw std::runtime_error("FMM candidate seed " + scientific(audit.accuracySeed()) +
                                         " residual " + scientific(result.boundaryResidual) +
                                         " exceeds tolerance " +
                                         scientific(options.relativeResidualTolerance) +
                                         "; bounded defect corrections did not converge.");
            report(BemSolvePhase::krylovSolve);
            result.iterations +=
                krylov.solve(residual.data(), correction.data(), result.krylovResidual);
            axpy(1., correction.data(), pressure);
            ++corrections;
            ++result.refinementSteps;
            evaluate(audit, pressure, result);
        }
    }
    bool preferDirect(std::span<const BemFastMultipoleLevelEvidence> levels) const {
        return directPolicy.preferDirect && directPolicy.preferDirect(levels, panelDegree);
    }
    SupportFmmSolution directSolve(const Z *flux, Z *pressure, SupportFmmSolution result) {
        captureOperatorBytes(*operation);
        releaseOperator(*operation);
        auto start = Clock::now();
        report(BemSolvePhase::residualAudit);
        auto audit = directPolicy.makeDirect(panelDegree);
        if (!audit)
            throw std::runtime_error("Direct audit factory returned no operator.");
        applyFlux(*audit, flux, rhs.data());
        result.residualCheckMilliseconds += elapsed(start);
        auto evaluateDirect = [&] {
            const auto checkStart = Clock::now();
            applyPressure(*audit, pressure, product.data());
            subtract(rhs.data(), product.data(), residual.data());
            result.boundaryResidual = ratio(norm(residual.data()), norm(rhs.data()));
            result.residualCheckMilliseconds += elapsed(checkStart);
            captureOperatorBytes(*audit);
        };
        evaluateDirect();
        krylovOperation = operation;
        for (unsigned corrections = 0; result.boundaryResidual > options.relativeResidualTolerance;
             ++corrections) {
            cancelled(stop);
            if (corrections == options.maximumRefinementSteps)
                throw std::runtime_error(
                    "Direct Galerkin audit residual " + scientific(result.boundaryResidual) +
                    " exceeds tolerance " + scientific(options.relativeResidualTolerance) +
                    "; bounded FMM defect corrections did not converge.");
            // Keep the direct buffers and a rebuilt Krylov hierarchy serial,
            // just as for two successive seed operators. Only vectors survive.
            audit.reset();
            report(BemSolvePhase::krylovSolve);
            result.iterations +=
                krylov.solve(residual.data(), correction.data(), result.krylovResidual);
            axpy(1., correction.data(), pressure);
            captureOperatorBytes(*operation);
            releaseOperator(*operation);
            const auto preparationStart = Clock::now();
            report(BemSolvePhase::residualAudit);
            audit = directPolicy.makeDirect(panelDegree);
            if (!audit)
                throw std::runtime_error("Direct audit factory returned no operator.");
            result.residualCheckMilliseconds += elapsed(preparationStart);
            ++result.refinementSteps;
            evaluateDirect();
        }
        result.residualEvaluation = BemResidualEvaluation::directGalerkin;
        result.operatorDiscrepancy = 0;
        result.comparisonAccuracySeed = result.verificationAccuracySeed = 0;
        result.comparisonLevels.clear();
        result.verificationLevels.clear();
        releaseOperator(*operation);
        return result;
    }
    SupportFmmSolution solveSeed(const Z *flux, Z *pressure, bool initial) {
        krylovOperation = operation;
        report(BemSolvePhase::krylovSolve);
        cancelled(stop);
        if (!pressure || (!flux && !operation->activeFaceIds().empty()) || flux == pressure)
            throw std::invalid_argument(
                "Support FMM solve needs valid distinct face-flux and pressure arrays.");
        SupportFmmSolution result;
        result.krylovTarget = options.relativeResidualTolerance * supportFmmKrylovSettingsScale *
                              cuda_detail::deviceGmresResidualFraction;
        captureOperatorBytes(*operation);
        result.accuracySeed = operation->accuracySeed();
        // Zero Neumann data and zero pressure solve the original discrete
        // equation exactly by linearity, including zero-flow coupled ports.
        double fluxNorm = 0;
        if (!operation->activeFaceIds().empty()) {
            check(cublasDznrm2(blas, static_cast<int>(operation->activeFaceIds().size()), flux, 1,
                               &fluxNorm),
                  "FMM face-flux norm");
            synchronize();
            if (!std::isfinite(fluxNorm))
                throw std::runtime_error("FMM face flux is non-finite.");
        }
        if (fluxNorm == 0) {
            homogeneous = true;
            clear(pressure);
            synchronize();
            result.residualEvaluation = BemResidualEvaluation::directGalerkin;
            releaseOperator(*operation);
            return result;
        }
        homogeneous = false;
        if (!initial || operation->levelEvidence().empty())
            applyFlux(*operation, flux, rhs.data());
        if (!initial) {
            krylovOperation = operation;
            report(BemSolvePhase::krylovSolve);
            result.iterations += krylov.solve(rhs.data(), correction.data(), result.krylovResidual);
            copy(correction.data(), pressure);
        }
        if (operation->levelEvidence().empty()) {
            result.residualEvaluation = BemResidualEvaluation::directGalerkin;
            evaluate(*operation, pressure, result);
            improve(*operation, pressure, result);
            return result;
        }
        if (preferDirect(operation->levelEvidence()))
            return directSolve(flux, pressure, std::move(result));
        double seed = operation->accuracySeed() / 10;
        captureOperatorBytes(*operation);
        releaseOperator(*operation);
        auto start = Clock::now();
        report(BemSolvePhase::residualAudit);
        auto audit = std::make_unique<SupportFmmOperator>(*operation, seed, stop);
        captureOperatorBytes(*audit);
        if (preferDirect(audit->levelEvidence())) {
            releaseOperator(*audit);
            audit.reset();
            return directSolve(flux, pressure, std::move(result));
        }
        applyFlux(*audit, flux, rhs.data());
        result.residualCheckMilliseconds += elapsed(start);
        evaluate(*audit, pressure, result);
        for (;;) {
            improve(*audit, pressure, result);
            result.comparisonAccuracySeed = seed;
            result.comparisonLevels.assign(audit->levelEvidence().begin(),
                                           audit->levelEvidence().end());
            copy(product.data(), comparisonProduct.data());
            copy(rhs.data(), comparisonRhs.data());
            synchronize();
            releaseOperator(*audit);
            audit.reset();
            for (;;) {
                seed /= 10;
                start = Clock::now();
                report(BemSolvePhase::residualAudit);
                audit = std::make_unique<SupportFmmOperator>(*operation, seed, stop);
                captureOperatorBytes(*audit);
                result.residualCheckMilliseconds += elapsed(start);
                if (preferDirect(audit->levelEvidence())) {
                    releaseOperator(*audit);
                    audit.reset();
                    return directSolve(flux, pressure, std::move(result));
                }
                const auto next = audit->levelEvidence();
                if (next.size() != result.comparisonLevels.size() ||
                    !std::equal(next.begin(), next.end(), result.comparisonLevels.begin()))
                    break;
                if (result.accuracyRefinements == maximumSupportFmmAccuracyRefinements)
                    throw std::runtime_error("FMM seed refinement exhausted without a distinct "
                                             "actual LF/HF level plan.");
                ++result.accuracyRefinements;
                releaseOperator(*audit);
                audit.reset();
            }
            start = Clock::now();
            applyFlux(*audit, flux, rhs.data());
            subtract(rhs.data(), comparisonRhs.data(), temporary.data());
            const double fluxDiscrepancy = norm(temporary.data());
            result.residualCheckMilliseconds += elapsed(start);
            const double rightNorm = evaluate(*audit, pressure, result);
            start = Clock::now();
            subtract(product.data(), comparisonProduct.data(), temporary.data());
            const double pressureDiscrepancy = norm(temporary.data());
            result.operatorDiscrepancy = ratio(pressureDiscrepancy + fluxDiscrepancy, rightNorm);
            result.verificationAccuracySeed = seed;
            result.verificationLevels.assign(audit->levelEvidence().begin(),
                                             audit->levelEvidence().end());
            result.residualCheckMilliseconds += elapsed(start);
            if (result.boundaryResidual <= options.relativeResidualTolerance &&
                result.operatorDiscrepancy <= .25 * options.relativeResidualTolerance) {
                releaseOperator(*audit);
                return result;
            }
            if (result.accuracyRefinements == maximumSupportFmmAccuracyRefinements)
                throw std::runtime_error(
                    "FMM accuracy-seed refinement exhausted: verification residual " +
                    scientific(result.boundaryResidual) + ", separate A/B discrepancy " +
                    scientific(result.operatorDiscrepancy) + "; requested tolerance " +
                    scientific(options.relativeResidualTolerance) + " and comparison budget " +
                    scientific(.25 * options.relativeResidualTolerance) + " were not both met.");
            ++result.accuracyRefinements;
        }
    }
    bool compareQuadrature(const Z *flux, const Z *pressure, SupportFmmSolution &result) {
        const auto nextDegree = galerkin::refinedEndpointDegree(panelDegree);
        if (!nextDegree)
            throw std::runtime_error("FMM endpoint quadrature has no independent refinement rule.");
        result.farQuadratureComparisonDegree = panelDegree;
        result.farQuadratureVerificationDegree = nextDegree;
        copy(product.data(), comparisonProduct.data());
        copy(rhs.data(), comparisonRhs.data());
        synchronize();
        captureOperatorBytes(*operation);
        releaseOperator(*operation);
        report(BemSolvePhase::residualAudit);
        const auto preparation = Clock::now();
        const auto apply = [&](auto &audit) {
            result.farQuadraturePreparationMilliseconds += elapsed(preparation);
            auto start = Clock::now();
            applyFlux(audit, flux, rhs.data());
            const auto rightNorm = norm(rhs.data());
            result.farQuadratureFluxMilliseconds += elapsed(start);
            start = Clock::now();
            applyPressure(audit, pressure, product.data());
            synchronize();
            result.farQuadraturePressureMilliseconds += elapsed(start);
            subtract(rhs.data(), comparisonRhs.data(), temporary.data());
            const double bDifference = norm(temporary.data());
            subtract(product.data(), comparisonProduct.data(), temporary.data());
            const double aDifference = norm(temporary.data());
            result.farQuadratureDiscrepancy = ratio(aDifference + bDifference, rightNorm);
            subtract(rhs.data(), product.data(), residual.data());
            result.boundaryResidual = ratio(norm(residual.data()), rightNorm);
            captureOperatorBytes(audit);
        };
        if (result.residualEvaluation == BemResidualEvaluation::directGalerkin) {
            auto audit = directPolicy.makeDirect(nextDegree);
            if (!audit)
                throw std::runtime_error("Direct quadrature factory returned no operator.");
            apply(*audit);
        } else {
            // Identical seed and support partition: this measures the panel
            // integration change, independently of the earlier seed change.
            SupportFmmOperator audit(*operation, result.verificationAccuracySeed,
                                     galerkin::triangleRule(nextDegree),
                                     SupportFmmOperatorRole::nativeAudit, stop);
            apply(audit);
            releaseOperator(audit);
        }
        return result.farQuadratureDiscrepancy <= .025 * options.relativeResidualTolerance &&
               result.boundaryResidual <= options.relativeResidualTolerance;
    }
    SupportFmmSolution solve(const Z *flux, Z *pressure, bool initial) {
        unsigned iterations = 0, corrections = 0, promotions = 0;
        double residualMilliseconds = 0, preparationMilliseconds = 0, pressureMilliseconds = 0,
               fluxMilliseconds = 0;
        for (;;) {
            auto result = solveSeed(flux, pressure, initial);
            iterations += result.iterations;
            corrections += result.refinementSteps;
            residualMilliseconds += result.residualCheckMilliseconds;
            const auto quadratureStart = Clock::now();
            const bool needsComparison = !homogeneous && !operation->levelEvidence().empty();
            const bool accepted = !needsComparison || compareQuadrature(flux, pressure, result);
            if (needsComparison)
                residualMilliseconds += elapsed(quadratureStart);
            preparationMilliseconds += result.farQuadraturePreparationMilliseconds;
            pressureMilliseconds += result.farQuadraturePressureMilliseconds;
            fluxMilliseconds += result.farQuadratureFluxMilliseconds;
            if (accepted) {
                result.iterations = iterations;
                result.refinementSteps = corrections;
                result.residualCheckMilliseconds = residualMilliseconds;
                result.farQuadratureRefinements = promotions;
                result.farQuadraturePreparationMilliseconds = preparationMilliseconds;
                result.farQuadraturePressureMilliseconds = pressureMilliseconds;
                result.farQuadratureFluxMilliseconds = fluxMilliseconds;
                krylovOperation = operation;
                return result;
            }
            if (panelDegree == 12)
                throw std::runtime_error(
                    "BEM endpoint quadrature 12/14 did not converge: separate A/B action "
                    "difference " +
                    scientific(result.farQuadratureDiscrepancy) + " exceeds budget " +
                    scientific(.025 * options.relativeResidualTolerance) +
                    " or higher-rule residual " + scientific(result.boundaryResidual) +
                    " exceeds tolerance " + scientific(options.relativeResidualTolerance) +
                    "; refine the mesh. Seed refinement cannot certify panel quadrature.");
            const auto nextDegree = galerkin::refinedEndpointDegree(panelDegree);
            ++promotions;
            releaseOperator(*operation);
            // Restore a complete reusable state before any operation that can
            // throw, including cancellation during the replacement constructor.
            operation = &sourceOperation;
            krylovOperation = &sourceOperation;
            refinedOperation.reset();
            panelDegree = sourcePanelDegree;
            // The existing pressure is reused through the genuine defect-
            // correction solve. Both seed evidences are then recomputed at the
            // NEW rule; an old seed indicator is never attached to changed p.
            const auto rebuilding = Clock::now();
            refinedOperation = std::make_unique<SupportFmmOperator>(
                sourceOperation, sourceOperation.accuracySeed(), galerkin::triangleRule(nextDegree),
                sourceOperation.role(), stop);
            preparationMilliseconds += elapsed(rebuilding);
            captureOperatorBytes(*refinedOperation);
            operation = refinedOperation.get();
            panelDegree = nextDegree;
            initial = true;
        }
    }
    std::size_t bytes() const {
        return krylov.deviceBytes() + rhs.bytes() + product.bytes() + residual.bytes() +
               correction.bytes() + comparisonProduct.bytes() + comparisonRhs.bytes() +
               temporary.bytes();
    }
    SupportFmmOperator &sourceOperation;
    SupportFmmOperator *operation;
    std::unique_ptr<SupportFmmOperator> refinedOperation;
    SupportFmmOperator *krylovOperation;
    int n;
    BemNumerics options;
    cudaStream_t stream;
    std::stop_token stop;
    SupportFmmProgress progress;
    BemSolvePhase currentPhase{BemSolvePhase::preparing};
    unsigned completedIterations{};
    std::size_t pressureApplications{}, fluxApplications{};
    cuda_detail::DeviceGmresSolver krylov;
    Buffer rhs, product, residual, correction, comparisonProduct, comparisonRhs, temporary;
    SupportDirectAuditPolicy directPolicy;
    unsigned panelDegree;
    const unsigned sourcePanelDegree;
    bool homogeneous{};
    cublasHandle_t blas{};
    std::size_t peakOperatorBytes{}, stationaryOperatorBytes{};
    BemFastMultipoleGraphEvidence graphEvidence;
};
SupportFmmLinearSolver::SupportFmmLinearSolver(SupportFmmOperator &operation, std::size_t n,
                                               const BemNumerics &options, std::stop_token stop,
                                               SupportDirectAuditPolicy directPolicy,
                                               SupportFmmProgress progress,
                                               unsigned panelQuadratureDegree)
    : implementation_(std::make_unique<Implementation>(operation, n, options, stop,
                                                       std::move(directPolicy), std::move(progress),
                                                       panelQuadratureDegree)) {}
SupportFmmLinearSolver::~SupportFmmLinearSolver() = default;
SupportFmmSolution SupportFmmLinearSolver::solve(const Z *flux, Z *pressure, bool initial) {
    return implementation_->solve(flux, pressure, initial);
}
std::size_t SupportFmmLinearSolver::deviceBytes() const noexcept {
    return implementation_->bytes();
}
std::size_t SupportFmmLinearSolver::peakOperatorDeviceBytes() const noexcept {
    return implementation_->peakOperatorBytes;
}
std::size_t SupportFmmLinearSolver::stationaryOperatorDeviceBytes() const noexcept {
    return implementation_->stationaryOperatorBytes;
}
BemFastMultipoleGraphEvidence SupportFmmLinearSolver::graphExecutionEvidence() const noexcept {
    return implementation_->graphEvidence;
}
} // namespace hornsculpt::wideband
