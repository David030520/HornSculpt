#pragma once

#include "support_direct_audit.hpp"
#include "support_fmm_operator.hpp"
#include <functional>

namespace hornsculpt::wideband {
inline constexpr unsigned maximumSupportFmmAccuracyRefinements = 4;
// Only the inner left-preconditioned solve is tightened. Outer residual,
// expansion seeds and quadrature qualification keep the requested numerics.
inline constexpr double supportFmmKrylovSettingsScale = .05;
struct SupportDirectAuditPolicy {
    // Optional bounded small-mesh alternative. The decision is re-evaluated
    // for every proposed seed plan; it never changes the Krylov operator.
    std::function<bool(std::span<const BemFastMultipoleLevelEvidence>, unsigned regularDegree)>
        preferDirect;
    std::function<std::unique_ptr<SupportDirectAudit>(unsigned regularDegree)> makeDirect;
};
// Heuristic operation-count comparison, not a timing guarantee. Direct work
// includes nine weak P1 coefficients, thirteen FP64 sample accumulators and
// the regular Green samples. RCR work
// includes two rotations and one axial m-block multiply per far interaction;
// diagonal work uses the actual full angular grid. Device FP32/FP64 throughput
// weights the regular-sample term. The hard mesh ceiling always takes priority.
bool preferSupportDirectAudit(const SupportTree &tree,
                              std::span<const BemFastMultipoleLevelEvidence> levels,
                              unsigned regularPoints, unsigned singleToDoubleThroughputRatio = 1);
struct SupportFmmSolution {
    BemResidualEvaluation residualEvaluation{BemResidualEvaluation::higherOrderFastMultipole};
    double boundaryResidual{}, krylovResidual{}, operatorDiscrepancy{};
    double krylovTarget{};
    double accuracySeed{}, comparisonAccuracySeed{}, verificationAccuracySeed{};
    double residualCheckMilliseconds{};
    unsigned farQuadratureComparisonDegree{}, farQuadratureVerificationDegree{},
        farQuadratureRefinements{};
    double farQuadratureDiscrepancy{}, farQuadraturePreparationMilliseconds{},
        farQuadraturePressureMilliseconds{}, farQuadratureFluxMilliseconds{};
    unsigned iterations{}, refinementSteps{}, accuracyRefinements{};
    std::vector<BemFastMultipoleLevelEvidence> comparisonLevels, verificationLevels;
};
// Bounded empirical convergence check, not a certified original-operator bound:
// candidate/verification use successively ten times smaller accuracy seeds.
// They share the exact near operator and geometric support partition. Only one
// audit state exists at a time; actual LF/HF plans must differ before acceptance.
using SupportFmmProgress = std::function<void(BemSolvePhase, unsigned, std::size_t, std::size_t)>;
class SupportFmmLinearSolver {
  public:
    SupportFmmLinearSolver(SupportFmmOperator &operation, std::size_t pressureCount,
                           const BemNumerics &options, std::stop_token stop = {},
                           SupportDirectAuditPolicy directPolicy = {},
                           SupportFmmProgress progress = {}, unsigned panelQuadratureDegree = 4);
    ~SupportFmmLinearSolver();
    SupportFmmLinearSolver(const SupportFmmLinearSolver &) = delete;
    SupportFmmLinearSolver &operator=(const SupportFmmLinearSolver &) = delete;
    SupportFmmSolution solve(const cuDoubleComplex *faceFlux, cuDoubleComplex *pressure,
                             bool useInitialGuess = false);
    std::size_t deviceBytes() const noexcept;
    // Maximum simultaneously owned operator payload, including shared near,
    // preconditioner and retained moment-family allocations exactly once.
    // Solver vectors are separate in deviceBytes(); these are not driver-wide
    // cudaMemGetInfo measurements or sums of independent historical maxima.
    std::size_t peakOperatorDeviceBytes() const noexcept;
    std::size_t stationaryOperatorDeviceBytes() const noexcept;
    // Cumulative across all solves on this frequency-local solver instance.
    BemFastMultipoleGraphEvidence graphExecutionEvidence() const noexcept;

  private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};
} // namespace hornsculpt::wideband
