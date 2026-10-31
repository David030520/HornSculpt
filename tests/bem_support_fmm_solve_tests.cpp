#include "bem_hf_test_oracle.hpp"
#include "bem_support_operator_oracle.hpp"
#include <cuda/support_fmm_solve.hpp>
#include <iostream>

namespace {
using namespace hornsculpt;
using namespace hornsculpt::wideband;
using C = std::complex<double>;
using Z = cuDoubleComplex;
void check(cudaError_t status) {
    if (status != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(status));
}
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
class Device {
  public:
    explicit Device(std::size_t count) : count_(count) {
        check(cudaMalloc(reinterpret_cast<void **>(&data), count * sizeof(Z)));
    }
    ~Device() {
        cudaFree(data);
    }
    Device(const Device &) = delete;
    Device &operator=(const Device &) = delete;
    void upload(std::span<const C> values, cudaStream_t stream) {
        require(values.size() == count_, "Solve fixture device size mismatch.");
        std::vector<Z> packed;
        for (const auto v : values)
            packed.push_back(make_cuDoubleComplex(v.real(), v.imag()));
        check(cudaMemcpyAsync(data, packed.data(), count_ * sizeof(Z), cudaMemcpyHostToDevice,
                              stream));
        check(cudaStreamSynchronize(stream));
    }
    std::vector<C> download(cudaStream_t stream) {
        std::vector<Z> packed(count_);
        check(cudaMemcpyAsync(packed.data(), data, count_ * sizeof(Z), cudaMemcpyDeviceToHost,
                              stream));
        check(cudaStreamSynchronize(stream));
        std::vector<C> result;
        for (const auto v : packed)
            result.emplace_back(v.x, v.y);
        return result;
    }
    Z *data{};

  private:
    std::size_t count_{};
};
std::vector<C> directSolve(std::vector<C> a, std::vector<C> b) {
    const auto n = b.size();
    for (std::size_t j = 0; j < n; ++j) {
        std::size_t pivot = j;
        for (auto i = j + 1; i < n; ++i)
            if (std::abs(a[i * n + j]) > std::abs(a[pivot * n + j]))
                pivot = i;
        require(std::abs(a[pivot * n + j]) > 0, "Direct Galerkin fixture is singular.");
        for (std::size_t column = 0; column < n; ++column)
            std::swap(a[j * n + column], a[pivot * n + column]);
        std::swap(b[j], b[pivot]);
        for (auto row = j + 1; row < n; ++row) {
            const C factor = a[row * n + j] / a[j * n + j];
            for (auto column = j + 1; column < n; ++column)
                a[row * n + column] -= factor * a[j * n + column];
            b[row] -= factor * b[j];
        }
    }
    std::vector<C> result(n);
    for (auto remaining = n; remaining; --remaining) {
        const auto row = remaining - 1;
        C value = b[row];
        for (auto column = row + 1; column < n; ++column)
            value -= a[row * n + column] * result[column];
        result[row] = value / a[row * n + row];
    }
    return result;
}
void checkSolve(bool hasFar, bool directAudit = false) {
    constexpr double k = 20;
    auto fixture = test::makeWidebandFarBlockFixture(k);
    for (auto &panel : fixture.space.panels)
        panel.boundary = BoundaryCondition::driven;
    SupportTreeSettings treeSettings;
    treeSettings.leafCapacity = hasFar ? 2 : 64;
    treeSettings.nearDistanceRatio = .5;
    auto tree = std::make_shared<SupportTree>(makeSupportTree(fixture.space, treeSettings));
    const std::array<std::uint32_t, 2> faces{5, 1};
    auto nearPlan = std::make_shared<SupportNearPlan>(
        makeSupportNearPlan(fixture.space, *tree, faces, 16 * 1024 * 1024));
    SupportNearOperatorSettings nearSettings;
    nearSettings.quadrature.nearDistanceRatio = .5;
    auto near =
        std::make_shared<SupportNearOperator>(fixture.space, *tree, *nearPlan, k, nearSettings);
    auto pc = std::make_shared<SupportBlockJacobi>(near->pressureMatrix(), tree->pressureOrder, 32,
                                                   std::stop_token{}, near->stream());
    BemNumerics options;
    options.relativeResidualTolerance = 5e-5;
    options.fastMultipole.restartLength = 4;
    options.fastMultipole.maximumIterations = 80;
    SupportFmmOperatorSettings settings;
    settings.role = SupportFmmOperatorRole::distancePrimary;
    settings.accuracySeed = .1 * options.relativeResidualTolerance;
    settings.maximumClassCacheBytes = 64 * 1024 * 1024;
    settings.panelRule = galerkin::triangleRule(4);
    SupportFmmOperator operation(fixture.space, tree, nearPlan, near, pc, k, settings);
    require(operation.levelEvidence().empty() != hasFar,
            "Solve fixture has incorrect actual near/far plan.");
    const auto reference =
        test::support_operator_oracle::directReference(fixture.space, k, nearSettings.quadrature);
    SupportDirectAuditPolicy directPolicy;
    if (directAudit) {
        directPolicy.preferDirect = [&](auto levels, unsigned degree) {
            return preferSupportDirectAudit(
                *tree, levels, static_cast<unsigned>(galerkin::triangleRule(degree).size()), 64);
        };
        require(directPolicy.preferDirect(operation.levelEvidence(), 4),
                "Small actual far fixture should choose the bounded direct audit.");
        directPolicy.makeDirect = [&](unsigned degree) {
            return std::make_unique<SupportDirectAudit>(fixture.space, *tree, *nearPlan, near, k,
                                                        degree);
        };
    }
    unsigned reportedIterations{}, totalIterations{};
    std::size_t reportedPressure{}, reportedFlux{};
    bool sawKrylov{}, sawAudit{};
    SupportFmmLinearSolver solver(
        operation, 8, options, {}, std::move(directPolicy),
        [&](BemSolvePhase phase, unsigned iterations, std::size_t pressureApplications,
            std::size_t fluxApplications) {
            require(iterations >= reportedIterations && pressureApplications >= reportedPressure &&
                        fluxApplications >= reportedFlux,
                    "Actual solve progress counters regressed across RHS/audit phases.");
            reportedIterations = iterations;
            reportedPressure = pressureApplications;
            reportedFlux = fluxApplications;
            sawKrylov |= phase == BemSolvePhase::krylovSolve;
            sawAudit |= phase == BemSolvePhase::residualAudit;
        });
    Device flux(2), pressure(8);
    BemFastMultipoleGraphEvidence previousGraph;
    for (unsigned rhs = 0; rhs < 2; ++rhs) {
        const std::vector<C> q =
            rhs ? std::vector<C>{{-.4, .2}, {.6, .3}} : std::vector<C>{{.3, -.7}, {0, 0}};
        flux.upload(q, operation.stream());
        const auto evidence = solver.solve(flux.data, pressure.data);
        const auto graph = solver.graphExecutionEvidence();
        require(graph.launches >= previousGraph.launches && graph.reuses >= previousGraph.reuses &&
                    graph.canonicalBatches >= previousGraph.canonicalBatches &&
                    graph.invalidatedExecutables >= previousGraph.invalidatedExecutables &&
                    graph.captureMilliseconds >= previousGraph.captureMilliseconds &&
                    graph.maximumCachedExecutables >= previousGraph.maximumCachedExecutables,
                "Frequency graph evidence regressed across serial RHS/audit lifetimes.");
        require(graph.reuses <= graph.launches && graph.launches <= graph.executedBatches &&
                    graph.instantiations <= graph.captureAttempts &&
                    graph.invalidatedExecutables <= graph.instantiations,
                "Frequency graph evidence duplicated serial audit counters.");
        if (!hasFar)
            require(graph == BemFastMultipoleGraphEvidence{},
                    "NoFar solve invented graph execution or memory observations.");
        previousGraph = graph;
        require(solver.stationaryOperatorDeviceBytes() >= operation.sharedDeviceBytes() &&
                    solver.peakOperatorDeviceBytes() >= solver.stationaryOperatorDeviceBytes() &&
                    solver.peakOperatorDeviceBytes() >= operation.peakDeviceBytes(),
                "Serial audits lost retained family/shared payload or a simultaneous peak.");
        totalIterations += evidence.iterations;
        require(reportedIterations == totalIterations && reportedPressure >= totalIterations &&
                    reportedFlux >= rhs + 1 && sawKrylov && sawAudit,
                "Solve progress omitted actual Krylov/operator work or audit phases.");
        const auto actual = pressure.download(operation.stream());
        std::vector<C> right(8);
        for (unsigned row = 0; row < 8; ++row)
            for (unsigned column = 0; column < 2; ++column)
                right[row] += reference.b[row * 2 + column] * q[column];
        const auto expected = directSolve(reference.a, right);
        double norm{}, difference{}, residualNorm{}, rightNorm{};
        for (unsigned row = 0; row < 8; ++row) {
            C residual = right[row];
            for (unsigned column = 0; column < 8; ++column)
                residual -= reference.a[row * 8 + column] * actual[column];
            residualNorm += std::norm(residual);
            rightNorm += std::norm(right[row]);
            difference += std::norm(actual[row] - expected[row]);
            norm += std::norm(expected[row]);
        }
        const double originalResidual = std::sqrt(residualNorm / rightNorm),
                     parity = std::sqrt(difference / norm);
        std::cout << "Full CUDA support solve far=" << hasFar << ", rhs=" << rhs
                  << ", iterations=" << evidence.iterations
                  << ", verification residual=" << evidence.boundaryResidual
                  << ", direct residual=" << originalResidual << ", pressure parity=" << parity
                  << ", discrepancy=" << evidence.operatorDiscrepancy << '\n';
        require(originalResidual < options.relativeResidualTolerance && parity < 2e-4,
                "Complete support solve fails independent direct Galerkin residual/parity.");
        require(evidence.boundaryResidual <= options.relativeResidualTolerance,
                "Complete support solve accepted a failed residual gate.");
        if (hasFar)
            require(evidence.farQuadratureComparisonDegree >= 4 &&
                        evidence.farQuadratureVerificationDegree >
                            evidence.farQuadratureComparisonDegree &&
                        evidence.farQuadratureDiscrepancy <=
                            .025 * options.relativeResidualTolerance,
                    "Complete support solve omitted the separate fixed-seed quadrature gate.");
        else
            require(evidence.farQuadratureComparisonDegree == 0 &&
                        evidence.farQuadratureVerificationDegree == 0,
                    "A no-far operator must not invent endpoint quadrature evidence.");
        if (hasFar && !directAudit) {
            require(
                evidence.residualEvaluation == BemResidualEvaluation::higherOrderFastMultipole &&
                    evidence.comparisonAccuracySeed > evidence.verificationAccuracySeed &&
                    !evidence.comparisonLevels.empty() && !evidence.verificationLevels.empty() &&
                    evidence.comparisonLevels != evidence.verificationLevels &&
                    evidence.operatorDiscrepancy <= .25 * options.relativeResidualTolerance,
                "Full support solve omitted distinct seed/level convergence evidence.");
            require(operation.ownDeviceBytes() == 0,
                    "Base hierarchy remained resident during serial seed audits.");
        } else {
            require(evidence.residualEvaluation == BemResidualEvaluation::directGalerkin &&
                        evidence.comparisonLevels.empty() && evidence.verificationLevels.empty() &&
                        evidence.comparisonAccuracySeed == 0 &&
                        evidence.verificationAccuracySeed == 0,
                    "Exact no-far solve fabricated seed refinement evidence.");
        }
        require(operation.activeFaceIds().size() == 2 && operation.activeFaceIds()[0] == 5 &&
                    operation.activeFaceIds()[1] == 1,
                "A zero-flow port changed the reusable active face basis.");
        require(operation.role() == SupportFmmOperatorRole::distancePrimary,
                "Serial audits changed the primary operator's immutable numerical role.");
    }
    flux.upload(std::vector<C>(2), operation.stream());
    const auto beforeZero = solver.graphExecutionEvidence();
    const auto zero = solver.solve(flux.data, pressure.data, true);
    const auto afterZero = solver.graphExecutionEvidence();
    require(afterZero.launches == beforeZero.launches && afterZero.reuses == beforeZero.reuses &&
                afterZero.canonicalBatches == beforeZero.canonicalBatches,
            "Zero RHS was counted as a graph or canonical GPU application.");
    require(zero.iterations == 0 && zero.boundaryResidual == 0 && zero.operatorDiscrepancy == 0,
            "Zero-flow support solve did not use exact linearity.");
    for (const auto p : pressure.download(operation.stream()))
        require(p == C{}, "Zero-flow solve retained prior pressure.");
}
} // namespace
int main() {
    int count{};
    if (cudaGetDeviceCount(&count) != cudaSuccess || !count)
        return 77;
    try {
        checkSolve(true);
        checkSolve(false);
        checkSolve(true, true);
        SupportTree large;
        large.fluxLeaf.resize(maximumDirectAuditPanels + 1);
        const std::array<BemFastMultipoleLevelEvidence, 1> levels{
            {{0, BemFastMultipoleRepresentation::regularSpherical, 128, 256, 256}}};
        require(!preferSupportDirectAudit(large, levels, 6, 64),
                "Direct audit cost policy ignored the hard large-mesh ceiling.");
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
