#pragma once

#include "../fmm/support_near_plan.hpp"
#include "../galerkin/pair_quadrature.hpp"
#include "../galerkin/space.hpp"
#include <array>
#include <complex>
#include <cuda_runtime_api.h>
#include <functional>
#include <memory>
#include <span>

namespace hornsculpt::wideband {
struct AdaptivePatchTask {
    std::uint32_t testPanel{}, sourcePanel{};
    galerkin::PairQuadraturePatch patch;
};
using AdaptivePatchMoments = galerkin::PairPatchEstimate;
// Batched FP64 evaluation only. CPU subdivision/heap decisions remain separate;
// A call evaluates independent patches in stream-ordered chunks bounded by
// both embedded rules' actual Green-sample count, with cancellation per chunk.
// degree8 uses degree6/8 (12/16 points), degree14 uses degree12/14 (49/64).
// Projected patches additionally compare inner Gauss orders 14/20 at the fine
// outer rule, preserving separate inner and outer error/refinement decisions.
// Space is borrowed for pair validation and must outlive this evaluator.
class AdaptivePatchEstimator {
  public:
    static constexpr std::size_t maximumChunkGreenSamples = 8 * 1024 * 1024;
    AdaptivePatchEstimator(const galerkin::Space &space, double wavenumber,
                           std::size_t maximumBatchTasks = 4096,
                           std::size_t maximumDeviceBytes = std::numeric_limits<std::size_t>::max(),
                           cudaStream_t stream = nullptr);
    ~AdaptivePatchEstimator();
    AdaptivePatchEstimator(const AdaptivePatchEstimator &) = delete;
    AdaptivePatchEstimator &operator=(const AdaptivePatchEstimator &) = delete;
    std::vector<AdaptivePatchMoments> evaluate(std::span<const AdaptivePatchTask> tasks,
                                               std::stop_token stop = {});
    // Borrowed result storage remains valid until the next evaluation or destruction.
    std::span<const AdaptivePatchMoments> evaluateView(std::span<const AdaptivePatchTask> tasks,
                                                       std::stop_token stop = {});
    std::size_t deviceBytes() const noexcept;
    // Owned pinned task/evaluation/status mirrors and the converted result
    // vector; borrowed geometry and immutable rule tables are excluded.
    std::size_t hostBytes() const noexcept;
    std::size_t lastChunkCount() const noexcept;
    std::size_t lastMaximumChunkGreenSamples() const noexcept;
    double lastKernelMilliseconds() const noexcept;

  private:
    struct Implementation;
    std::unique_ptr<Implementation> impl_;
};

struct AdaptivePairBatchStats {
    std::size_t rounds{}, evaluatedPatches{}, completedPairs{}, maximumRoundPatches{};
    std::size_t liveWorkspaceBytes{}, peakWorkspaceBytes{}, maximumSinglePairWorkspaceBytes{};
    double decisionMilliseconds{}, estimatorKernelMilliseconds{}, estimatorWallMilliseconds{};
};

// Many independently advancing pair heaps share a hard host-payload budget.
// Every round evaluates all pending patches in bounded CUDA chunks. Completed
// rules are consumed before their slots are reused; no global patch/job pool.
class AdaptivePairBatch {
  public:
    using Consumer =
        std::function<void(const SupportNearJob &, const galerkin::PairQuadratureRule &)>;
    AdaptivePairBatch(const galerkin::Space &space, double wavenumber,
                      galerkin::PairQuadratureSettings settings, std::size_t maximumHostBytes,
                      std::size_t maximumDeviceBytes = std::numeric_limits<std::size_t>::max(),
                      std::size_t maximumActivePairs = 512, cudaStream_t stream = nullptr,
                      std::stop_token stop = {});
    ~AdaptivePairBatch();
    bool submit(const SupportNearJob &job);
    void advance(const Consumer &consume);
    void drain(const Consumer &consume);
    std::size_t deviceBytes() const noexcept;
    AdaptivePairBatchStats statistics() const noexcept;

  private:
    struct Implementation;
    std::unique_ptr<Implementation> impl_;
};
} // namespace hornsculpt::wideband
