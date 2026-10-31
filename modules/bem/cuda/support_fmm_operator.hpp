#pragma once

#include "../fmm/translation_plan_cache.hpp"
#include "support_near_operator.hpp"
#include "support_preconditioner.hpp"
#include "wideband_operator.hpp"
#include <hornsculpt/bem.hpp>

namespace hornsculpt::wideband {
struct SupportFmmWorkStatistics {
    // Class-kind order: rotation, axial, angular, diagonal. Counts include
    // requests during complete applications, not hypothetical retained plans.
    std::array<std::size_t, 4> classCacheHits{}, classCacheMisses{};
    std::size_t applications{}, azimuthUploads{};
    std::size_t rcrBatchesQueued{}, rcrCompletionChecks{}, maximumRcrChunkWork{};
    unsigned maximumRcrChunkBatches{};
    std::size_t leafMomentCacheBytes{}, leafCachedPressureDofs{}, leafCachedFaceSlots{};
    double leafMomentConstructionMilliseconds{};
    // Graph recording is not execution. Direct and graph batches partition
    // actual LF-M2L work; opaque driver graph storage is never owned payload.
    std::size_t m2lGraphCaptureAttempts{}, m2lGraphInstantiations{}, m2lGraphLaunches{},
        m2lGraphReuses{}, m2lGraphExecutedBatches{}, m2lGraphDirectChunks{},
        m2lGraphDirectBatches{}, m2lGraphInvalidatedExecutables{}, m2lGraphBudgetBypasses{},
        m2lGraphAdmissionClosures{};
    std::size_t m2lGraphCachedEntries{}, m2lGraphCachedExecutables{}, m2lGraphCachedNodes{},
        m2lGraphMaximumCachedNodes{}, m2lGraphMaximumCachedExecutables{},
        m2lGraphOwnedHostPayloadBytes{}, m2lGraphPeakOwnedHostPayloadBytes{};
    bool m2lGraphAdmissionClosed{};
    double m2lGraphCaptureMilliseconds{}, m2lGraphInstantiationMilliseconds{};
    // cudaMemGetInfo observations include other allocations and driver state;
    // they are not an attribution of exact bytes to graph executables.
    std::size_t m2lGraphMemoryQueries{}, m2lGraphLastFreeBeforeBytes{},
        m2lGraphLastFreeAfterBytes{}, m2lGraphMinimumObservedFreeBytes{},
        m2lGraphTotalDeviceBytes{}, m2lGraphRequiredFreeReserveBytes{};
    double cpuPlanPreparationMilliseconds{}, classPreparationMilliseconds{},
        azimuthUploadMilliseconds{};
    // Synchronous wall time of the last application; class construction is
    // included in the phase that actually needed that class.
    double nearMilliseconds{}, emissionMilliseconds{}, upwardMilliseconds{},
        translationMilliseconds{}, downwardMilliseconds{}, testingMilliseconds{};
};
// Internal numerical roles; audits never inherit the primary's endpoint cutoffs.
enum class SupportFmmOperatorRole { nativeAudit, distancePrimary };
struct SupportFmmOperatorSettings {
    SupportFmmOperatorRole role{SupportFmmOperatorRole::nativeAudit};
    double accuracySeed{5e-6};
    unsigned classBatchNodes{32};
    // Shared across axial distance classes, directional rotations and HF
    // diagonal classes. Evicted classes are reconstructed, never retained in
    // an unbounded host mirror. Preparation/reconstruction time is real work.
    std::size_t maximumClassCacheBytes{256 * 1024 * 1024};
    // One shared immutable rule bank across serial seed hierarchies, including
    // genuinely integrated base/candidate/verifier degrees. Actual free VRAM,
    // largest planned state and reserved class space also bound it.
    std::size_t maximumLeafMomentCacheBytes{std::size_t{2} * 1024 * 1024 * 1024};
    // Admission only on a second exact chunk contact. Reaching any graph cap
    // closes admission for this hierarchy; existing executable hits survive.
    std::size_t maximumM2lGraphEntries{16384}, maximumM2lGraphExecutables{8192},
        maximumM2lGraphNodes{524288}, maximumM2lGraphHostBytes{8 * 1024 * 1024};
    std::size_t maximumHostPlanBytes{512 * 1024 * 1024};
    unsigned maximumPlanWorkers{4};
    std::size_t maximumDeviceBytes{std::numeric_limits<std::size_t>::max()};
    std::vector<galerkin::TriangleQuadraturePoint> panelRule;
};

// Complete matrix-free weak P1/face-DP0 Burton-Miller operator. Exact near,
// support partition and spatial preconditioner are immutable shared objects;
// each accuracy seed owns its complete LF/HF level plan and device states.
// The space/tree/plan snapshot must outlive all operators that reference it.
class SupportFmmOperator {
  public:
    SupportFmmOperator(const galerkin::Space &space, std::shared_ptr<const SupportTree> tree,
                       std::shared_ptr<const SupportNearPlan> nearPlan,
                       std::shared_ptr<SupportNearOperator> near,
                       std::shared_ptr<SupportBlockJacobi> preconditioner, double wavenumber,
                       SupportFmmOperatorSettings settings, std::stop_token stop = {});
    // A seed-only successor is always an uncut native audit.
    SupportFmmOperator(SupportFmmOperator &source, double accuracySeed, std::stop_token stop = {});
    // A quadrature action comparison keeps the seed/partition fixed and only
    // replaces the canonical leaf integration rule. As with seed refinement,
    // callers release the previous hierarchy before creating its successor.
    SupportFmmOperator(SupportFmmOperator &source, double accuracySeed,
                       std::vector<galerkin::TriangleQuadraturePoint> panelRule,
                       SupportFmmOperatorRole role, std::stop_token stop = {});
    ~SupportFmmOperator();
    SupportFmmOperator(const SupportFmmOperator &) = delete;
    SupportFmmOperator &operator=(const SupportFmmOperator &) = delete;
    // Arrays use global P1 numbering and nearPlan.activeFaceIds order. They
    // must belong to this stream's execution order and may not alias.
    void applyPressureDevice(const cuDoubleComplex *pressure, cuDoubleComplex *output,
                             std::stop_token stop = {});
    void applyFaceFluxDevice(const cuDoubleComplex *faceFlux, cuDoubleComplex *output,
                             std::stop_token stop = {});
    void preconditionDevice(cuDoubleComplex *values, std::stop_token stop = {});
    // Serial seed comparison keeps only one hierarchy resident. Immutable
    // shared state and evidence survive; the next application rebuilds this
    // seed's hierarchy before use. Completed peak evidence is retained.
    void releaseDeviceState();
    cudaStream_t stream() const noexcept;
    std::span<const std::uint32_t> activeFaceIds() const noexcept;
    const FastMultipoleStatistics &statistics() const noexcept;
    std::span<const BemFastMultipoleLevelEvidence> levelEvidence() const noexcept;
    double accuracySeed() const noexcept;
    SupportFmmOperatorRole role() const noexcept;
    // Actual immutable M2L endpoint policy. Requires a resident hierarchy;
    // level states and M2M/L2L retain their native expansion degrees. Only the
    // primary's separated LF classes may shorten both endpoint prefixes.
    unsigned m2lSourceOrder(unsigned level, double centerDistance) const;
    unsigned m2lTargetOrder(unsigned level, double centerDistance) const;
    bool hasDistanceTargetGraphs() const noexcept;
    std::span<const galerkin::TriangleQuadraturePoint> panelRule() const noexcept;
    // Explicit buffers/class payload only. Opaque driver graph allocations
    // are observed via separate cudaMemGetInfo statistics, never invented bytes.
    std::size_t deviceBytes() const noexcept;
    // Live allocations unique to this seed. sharedDeviceBytes() comprises
    // exact near, spatial preconditioner and the bounded LF moment family,
    // counted once when serial audit objects coexist. A rule/degree miss can
    // replace the idle family, so callers sample shared bytes after preparation.
    std::size_t ownDeviceBytes() const noexcept;
    std::size_t sharedDeviceBytes() const noexcept;
    std::size_t ownPeakDeviceBytes() const noexcept;
    // Maximum completed preparation/application residency; excludes a
    // transient replacement plan coexisting with the cache it will evict.
    std::size_t ownStationaryDeviceBytes() const noexcept;
    const SupportFmmWorkStatistics &workStatistics() const noexcept;
    TranslationPlanCacheStatistics hostPlanStatistics() const;
    // Maxima of simultaneously resident own+shared allocations. Unlike adding
    // separate historical maxima, these remain valid across family/rule changes.
    std::size_t stationaryDeviceBytes() const noexcept;
    std::size_t peakDeviceBytes() const noexcept;

  private:
    void ensureState(std::stop_token stop);
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
} // namespace hornsculpt::wideband
