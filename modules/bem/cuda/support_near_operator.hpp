#pragma once

#include "../fmm/support_near_plan.hpp"
#include "../fmm/wideband_transfer.hpp"
#include "../galerkin/pair_quadrature.hpp"
#include <array>
#include <cuComplex.h>
#include <cuda_runtime_api.h>
#include <functional>
#include <memory>

namespace hornsculpt::wideband {
struct SupportNearOperatorSettings {
    galerkin::PairQuadratureSettings quadrature;
    std::size_t maximumDeviceBytes{std::numeric_limits<std::size_t>::max()};
    // Per-batch job/patch/mask capacity. Two pinned host slots share an explicit
    // total bound of twice this value; the device batch exists only once.
    // The adaptive quadrature's own bounded working state is separate.
    std::size_t maximumConstructionBatchBytes{16 * 1024 * 1024};
    std::size_t maximumIteratorWorkspaceBytes{64 * 1024 * 1024};
    // Temporary O(T) centroid-sphere qualification bounds, separate from jobs.
    std::size_t maximumClassificationWorkspaceBytes{64 * 1024 * 1024};
    unsigned maximumBatchJobs{32768};
    // Bounds actual Green samples, including singular sectors and adaptive
    // patch rules. Byte capacity alone does not bound a WDDM kernel's work.
    std::size_t maximumBatchGreenSamples{8 * 1024 * 1024};
    bool profileConstruction{};
    // Internal CUDA estimator capacity, not a second CPU near backend.
    unsigned maximumActiveAdaptivePairs{512};
    std::size_t maximumAdaptiveWorkspaceBytes{256 * 1024 * 1024};
};
struct SupportNearConstructionStats {
    std::size_t maskedPanelJobs{}, quadratureJobs{}, batches{}, patches{};
    std::size_t maximumBatchJobs{}, maximumBatchPatches{};
    std::size_t submittedBatches{}, maximumResidentBatches{};
    // A query observed unfinished assembly when host preparation began. This
    // counts starts, not elapsed overlap or GPU utilization.
    std::size_t preparationsStartedWhileAssemblyPending{};
    std::size_t pinnedConstructionHostBytes{}, maximumPinnedConstructionHostBytes{},
        deviceConstructionBatchBytes{};
    std::size_t iteratorWorkspaceBytes{};
    std::size_t maximumBatchGreenSamples{}, liveDeviceBytes{}, peakDeviceBytes{};
    double maximumKernelMilliseconds{};
    double kernelMilliseconds{}, batchUploadAndSyncWallMilliseconds{};
    // Host time inside submission/retirement only; intervening preparation is
    // excluded. These intervals overlap GPU event time and are not additive.
    double batchSubmissionMilliseconds{}, batchRetirementMilliseconds{};
    double iteratorMilliseconds{}, classificationMilliseconds{}, adaptiveDecisionMilliseconds{};
    double estimatorKernelMilliseconds{}, estimatorWallMilliseconds{};
    std::size_t adaptiveQuadratureCalls{};
    std::size_t classificationWorkspaceBytes{};
    std::array<std::size_t, 3> regularRuleJobs{}; // Actual unique jobs at degrees 4/6/8.
    std::size_t iteratorBatchBytes{}, adaptiveWorkspacePeakBytes{},
        maximumSinglePairWorkspaceBytes{};
    std::size_t adaptiveRounds{}, adaptivePatchEvaluations{}, maximumEstimatorBatchPatches{};
};
struct SupportDeviceCsr {
    const std::uint32_t *rowOffsets{}, *columns{};
    const cuDoubleComplex *values{};
    std::uint32_t rows{}, columnCount{}, nonzeros{};
};

// Exact near part of the same quadrature-discretized weak Burton-Miller form.
// Full panel supports are integrated once per owned mask; adaptive patch chunks
// carry disjoint integration domains. No wall/inactive flux column is stored.
class SupportNearOperator {
  public:
    SupportNearOperator(const galerkin::Space &space, const SupportTree &tree,
                        const SupportNearPlan &plan, double wavenumber,
                        SupportNearOperatorSettings settings = {}, std::stop_token stop = {},
                        cudaStream_t stream = nullptr,
                        std::function<void(const SupportNearConstructionStats &)> progress = {});
    ~SupportNearOperator();
    SupportNearOperator(const SupportNearOperator &) = delete;
    SupportNearOperator &operator=(const SupportNearOperator &) = delete;
    std::vector<Complex> applyPressure(std::span<const Complex> values, std::stop_token stop = {});
    // Face inputs are Neumann flux q=-i*rho*omega*v, in activeFaceIds order.
    // The solver converts public RMS normal velocities exactly once.
    std::vector<Complex> applyFaceFlux(std::span<const Complex> values, std::stop_token stop = {});
    void applyPressureDevice(const cuDoubleComplex *values, cuDoubleComplex *output,
                             std::stop_token stop = {});
    void applyFaceFluxDevice(const cuDoubleComplex *values, cuDoubleComplex *output,
                             std::stop_token stop = {});
    SupportDeviceCsr pressureMatrix() const noexcept;
    SupportDeviceCsr faceFluxMatrix() const noexcept;
    cudaStream_t stream() const noexcept;
    std::size_t deviceBytes() const noexcept;
    std::size_t peakDeviceBytes() const noexcept;
    SupportNearConstructionStats constructionStats() const noexcept;

  private:
    class Implementation;
    std::unique_ptr<Implementation> impl_;
};
} // namespace hornsculpt::wideband
