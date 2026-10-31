#pragma once

#include "scaled_translation.hpp"
#include <memory>
#include <stop_token>

namespace hornsculpt::wideband {
struct AxialPlanRequest {
    ScaledAxialKind kind{};
    unsigned sourceOrder{}, targetOrder{}, polarCount{};
    double wavenumber{}, sourceRadius{}, targetRadius{}, distance{}, projectionRadius{};
};
struct TranslationPlanCacheStatistics {
    std::size_t hits{}, plansBuilt{};
    // Sum of elapsed independent builder calls (parallel sums can exceed wall
    // time); no timing from GPU upload or application is included.
    double builderMilliseconds{};
};

// A solve and its serial audits share this bounded repository of immutable
// CPU plans. There is no global cache. Independent distance classes can be
// prepared by a fixed worker pool; only completed plans enter the LRU.
class TranslationPlanCache {
  public:
    explicit TranslationPlanCache(std::size_t maximumBytes, unsigned maximumWorkers = 4);
    ~TranslationPlanCache();
    TranslationPlanCache(const TranslationPlanCache &) = delete;
    TranslationPlanCache &operator=(const TranslationPlanCache &) = delete;
    std::shared_ptr<const ScaledAxialPlan> axial(const AxialPlanRequest &request,
                                                 std::stop_token stop = {});
    void prepare(std::span<const AxialPlanRequest> requests, std::stop_token stop = {});
    std::size_t liveBytes() const noexcept;
    // Conservative high-water bound including reserved in-flight builder
    // scratch. liveBytes() counts actual completed-plan vector capacities.
    std::size_t peakBytes() const noexcept;
    TranslationPlanCacheStatistics statistics() const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
} // namespace hornsculpt::wideband
