#include "translation_plan_cache.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <condition_variable>
#include <map>
#include <mutex>
#include <stdexcept>
#include <thread>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace hornsculpt::wideband {
namespace {
using Key = std::array<std::uint64_t, 10>;
Key key(const AxialPlanRequest &r) {
    return {static_cast<std::uint64_t>(r.kind),
            r.sourceOrder,
            r.targetOrder,
            r.polarCount,
            std::bit_cast<std::uint64_t>(r.wavenumber),
            std::bit_cast<std::uint64_t>(r.sourceRadius),
            std::bit_cast<std::uint64_t>(r.targetRadius),
            std::bit_cast<std::uint64_t>(r.distance),
            std::bit_cast<std::uint64_t>(r.projectionRadius),
            0};
}
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Support translation plan preparation cancelled.");
}
} // namespace
struct TranslationPlanCache::Impl {
    struct Entry {
        std::shared_ptr<const ScaledAxialPlan> plan;
        std::size_t bytes{}, touch{};
    };
    mutable std::mutex mutex;
    std::condition_variable changed;
    std::map<Key, Entry> entries;
    std::size_t limit{}, bytes{}, reserved{}, peak{}, clock{};
    unsigned workers{};
    TranslationPlanCacheStatistics stats;
    Impl(std::size_t maximumBytes, unsigned maximumWorkers)
        : limit(maximumBytes), workers(maximumWorkers) {
        if (!limit || !workers || workers > 32)
            throw std::invalid_argument("Support translation CPU cache limits are invalid.");
        workers = std::min(workers, std::max(1u, std::thread::hardware_concurrency()));
#ifdef _WIN32
        MEMORYSTATUSEX memory{};
        memory.dwLength = sizeof(memory);
        if (!GlobalMemoryStatusEx(&memory))
            throw std::runtime_error("Support translation host memory query failed.");
        limit = std::min(limit, static_cast<std::size_t>(memory.ullAvailPhys / 4));
        if (!limit)
            throw std::runtime_error("No host memory remains for support translation plans.");
#endif
    }
    std::shared_ptr<const ScaledAxialPlan> axial(const AxialPlanRequest &r, std::stop_token stop) {
        cancelled(stop);
        const auto requestKey = key(r);
        {
            std::lock_guard lock(mutex);
            const auto found = entries.find(requestKey);
            if (found != entries.end()) {
                ++stats.hits;
                found->second.touch = ++clock;
                return found->second.plan;
            }
        }
        if (r.sourceOrder > maximumScaledSphericalOrder ||
            r.targetOrder > maximumScaledSphericalOrder || r.polarCount > 1024)
            throw std::invalid_argument("Axial CPU request exceeds its numerical domain.");
        const auto maximumOrder = std::max(r.sourceOrder, r.targetOrder);
        std::size_t entriesCount = 0;
        for (unsigned m = 0; m <= std::min(r.sourceOrder, r.targetOrder); ++m)
            entriesCount +=
                static_cast<std::size_t>(r.sourceOrder - m + 1) * (r.targetOrder - m + 1);
        // The builder allocates one GL sphere and bounded source/target wave
        // tables. Count their conservative capacities before any worker starts.
        const std::size_t reservation =
            sizeof(ScaledAxialPlan) + entriesCount * sizeof(Complex) +
            2 * (std::min(r.sourceOrder, r.targetOrder) + 2) * sizeof(std::size_t) +
            2 * static_cast<std::size_t>(r.polarCount) * r.polarCount * sizeof(SphereSample) +
            8 * static_cast<std::size_t>(maximumOrder + 2) * (maximumOrder + 2) *
                sizeof(ScalarWave) +
            64 * 1024;
        if (reservation > limit)
            throw std::runtime_error("One axial CPU builder exceeds the shared host budget.");
        {
            std::unique_lock lock(mutex);
            while (bytes + reserved > limit - reservation) {
                auto oldest = entries.end();
                for (auto it = entries.begin(); it != entries.end(); ++it)
                    if (it->second.plan.use_count() == 1 &&
                        (oldest == entries.end() || it->second.touch < oldest->second.touch))
                        oldest = it;
                if (oldest != entries.end()) {
                    bytes -= oldest->second.bytes;
                    entries.erase(oldest);
                    continue;
                }
                if (!reserved)
                    throw std::runtime_error("Pinned CPU plans leave insufficient builder memory.");
                cancelled(stop);
                changed.wait_for(lock, std::chrono::milliseconds(25));
            }
            reserved += reservation;
            peak = std::max(peak, bytes + reserved);
        }
        std::shared_ptr<const ScaledAxialPlan> plan;
        const auto buildStart = std::chrono::steady_clock::now();
        try {
            plan = std::make_shared<const ScaledAxialPlan>(makeScaledAxialPlan(
                r.kind, r.sourceOrder, r.targetOrder, r.wavenumber, r.sourceRadius, r.targetRadius,
                r.distance, r.polarCount, r.projectionRadius, stop));
            cancelled(stop);
        } catch (...) {
            std::lock_guard lock(mutex);
            reserved -= reservation;
            changed.notify_all();
            throw;
        }
        const auto payload = sizeof(ScaledAxialPlan) +
                             plan->blockOffsets.capacity() * sizeof(std::size_t) +
                             plan->coefficients.capacity() * sizeof(Complex);
        std::lock_guard lock(mutex);
        ++stats.plansBuilt;
        stats.builderMilliseconds +=
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - buildStart)
                .count();
        reserved -= reservation;
        changed.notify_all();
        if (payload > reservation)
            throw std::runtime_error("Axial CPU builder exceeded its reserved payload capacity.");
        // Concurrent duplicate requests return the first completed object.
        if (auto found = entries.find(requestKey); found != entries.end()) {
            ++stats.hits;
            found->second.touch = ++clock;
            return found->second.plan;
        }
        entries.emplace(requestKey, Entry{plan, payload, ++clock});
        bytes += payload;
        return plan;
    }
};
TranslationPlanCache::TranslationPlanCache(std::size_t maximumBytes, unsigned maximumWorkers)
    : impl_(std::make_unique<Impl>(maximumBytes, maximumWorkers)) {}
TranslationPlanCache::~TranslationPlanCache() = default;
std::shared_ptr<const ScaledAxialPlan> TranslationPlanCache::axial(const AxialPlanRequest &request,
                                                                   std::stop_token stop) {
    return impl_->axial(request, stop);
}
void TranslationPlanCache::prepare(std::span<const AxialPlanRequest> requests,
                                   std::stop_token stop) {
    cancelled(stop);
    // Unique requests prevent duplicate in-flight plans within this pool.
    std::map<Key, AxialPlanRequest> unique;
    for (const auto &request : requests)
        unique.emplace(key(request), request);
    std::vector<AxialPlanRequest> work;
    for (const auto &[unused, request] : unique) {
        (void)unused;
        work.push_back(request);
    }
    std::atomic_size_t next{};
    std::atomic_bool failed{};
    std::mutex failureMutex;
    std::exception_ptr failure;
    const auto count = std::min<std::size_t>(impl_->workers, work.size());
    std::vector<std::jthread> workers;
    for (std::size_t worker = 0; worker < count; ++worker)
        workers.emplace_back([&] {
            while (!failed.load(std::memory_order_relaxed)) {
                const auto index = next.fetch_add(1);
                if (index >= work.size())
                    break;
                try {
                    (void)impl_->axial(work[index], stop);
                } catch (...) {
                    std::lock_guard lock(failureMutex);
                    if (!failure)
                        failure = std::current_exception();
                    failed.store(true, std::memory_order_relaxed);
                }
            }
        });
    workers.clear(); // jthread joins: no background work escapes this call.
    if (failure)
        std::rethrow_exception(failure);
    cancelled(stop);
}
std::size_t TranslationPlanCache::liveBytes() const noexcept {
    std::lock_guard lock(impl_->mutex);
    return impl_->bytes;
}
std::size_t TranslationPlanCache::peakBytes() const noexcept {
    std::lock_guard lock(impl_->mutex);
    return impl_->peak;
}
TranslationPlanCacheStatistics TranslationPlanCache::statistics() const {
    std::lock_guard lock(impl_->mutex);
    return impl_->stats;
}
} // namespace hornsculpt::wideband
