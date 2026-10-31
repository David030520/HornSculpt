#include "interaction_plan.hpp"
#include "regular_qualification.hpp"
#include "space.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <compare>
#include <exception>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>

namespace hornsculpt::galerkin {
namespace {
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("BEM interaction planning cancelled.");
}

std::uint32_t index32(std::size_t value) {
    if (value > std::numeric_limits<std::uint32_t>::max())
        throw std::runtime_error("BEM interaction plan exceeds 32-bit indexing limits.");
    return static_cast<std::uint32_t>(value);
}

struct Cell {
    std::int64_t x{}, y{}, z{};
    auto operator<=>(const Cell &) const = default;
};
using PanelPair = std::pair<std::uint32_t, std::uint32_t>;

template <class T>
void reserveTracked(std::vector<T> &values, std::size_t count,
                    AdaptiveQuadratureWorkspace &budget) {
    if (count <= values.capacity())
        return;
    const auto capacity = std::max(count, values.capacity() + values.capacity() / 2 + 1);
    if (capacity > std::numeric_limits<std::size_t>::max() / sizeof(T))
        throw std::runtime_error("Dense quadrature plan capacity overflow.");
    const auto oldBytes = values.capacity() * sizeof(T), newBytes = capacity * sizeof(T);
    budget.acquire(newBytes);
    try {
        values.reserve(capacity);
    } catch (...) {
        budget.release(newBytes);
        throw;
    }
    budget.release(oldBytes);
    if (values.capacity() > capacity)
        budget.acquire((values.capacity() - capacity) * sizeof(T));
}

std::vector<PanelPair> candidates(const Space &space, double ratio, std::stop_token stop,
                                  AdaptiveQuadratureWorkspace &budget, const std::vector<bool> &active) {
    double maximumEdge = 0;
    for (const auto &panel : space.panels) {
        if (!(panel.maximumEdgeLength > 0) || !std::isfinite(panel.maximumEdgeLength))
            throw std::invalid_argument("Interaction plan has an invalid panel size.");
        maximumEdge = std::max(maximumEdge, panel.maximumEdgeLength);
    }
    // A panel point is at most its maximum edge length from vertex 0. Thus every
    // true near pair has vertex-0 separation <= (ratio+2)*maximumEdge, including
    // differently sized panels. Only 27 neighboring cells can contain it.
    const auto origin = space.panels.front().vertices[0];
    double extent = 0;
    for (const auto &panel : space.panels) {
        cancelled(stop);
        const auto point = panel.vertices[0];
        extent = std::max({extent, std::abs(point.x - origin.x), std::abs(point.y - origin.y),
                           std::abs(point.z - origin.z)});
    }
    // This margin also bounds the normalized grid coordinates, so floating
    // subtraction/division cannot separate a real near pair by two cells.
    const double width =
        (ratio + 2) * maximumEdge + 128 * std::numeric_limits<double>::epsilon() * extent;
    if (!(width > 0) || !std::isfinite(width))
        throw std::invalid_argument("Interaction-plan grid scale is not finite.");
    auto coordinate = [width](double value, double base) {
        const double cell = std::floor((value - base) / width);
        // Keep a margin for neighbor offsets, avoiding signed integer overflow.
        constexpr double limit = static_cast<double>(std::numeric_limits<std::int64_t>::max()) / 2;
        if (!std::isfinite(cell) || cell < -limit || cell > limit)
            throw std::invalid_argument(
                "Interaction-plan spatial index exceeds its coordinate range.");
        return static_cast<std::int64_t>(cell);
    };
    using IndexedCell = std::pair<Cell, std::uint32_t>;
    std::vector<IndexedCell> grid;
    reserveTracked(grid, space.panels.size(), budget);
    std::vector<Cell> locations;
    reserveTracked(locations, space.panels.size(), budget);
    for (std::size_t index = 0; index < space.panels.size(); ++index) {
        cancelled(stop);
        const auto point = space.panels[index].vertices[0];
        const Cell cell{coordinate(point.x, origin.x), coordinate(point.y, origin.y),
                        coordinate(point.z, origin.z)};
        locations.push_back(cell);
        grid.emplace_back(cell, index32(index));
    }
    std::sort(grid.begin(), grid.end());
    std::vector<PanelPair> result;
    for (std::size_t index = 0; index < space.panels.size(); ++index) {
        cancelled(stop);
        const auto cell = locations[index];
        for (int x = -1; x <= 1; ++x)
            for (int y = -1; y <= 1; ++y)
                for (int z = -1; z <= 1; ++z) {
                    const Cell neighbor{cell.x + x, cell.y + y, cell.z + z};
                    auto found = std::lower_bound(grid.begin(), grid.end(), neighbor,
                                                  [](const IndexedCell &entry, const Cell &key) {
                                                      return entry.first < key;
                                                  });
                    for (; found != grid.end() && found->first == neighbor; ++found) {
                        const auto other = found->second;
                        if (other <= index)
                            continue;
                        if (!active[index] && !active[other])
                            continue;
                        const auto &a = space.panels[index], &b = space.panels[other];
                        const double bound =
                            ratio * std::max(a.maximumEdgeLength, b.maximumEdgeLength) +
                            a.maximumEdgeLength + b.maximumEdgeLength;
                        if (std::abs(a.vertices[0].x - b.vertices[0].x) > bound ||
                            std::abs(a.vertices[0].y - b.vertices[0].y) > bound ||
                            std::abs(a.vertices[0].z - b.vertices[0].z) > bound)
                            continue;
                        if (classifyPairIds(space.panels[index].nodes.data(),
                                            space.panels[other].nodes.data())
                                .topology == PairTopology::disjoint) {
                            reserveTracked(result, result.size() + 1, budget);
                            result.emplace_back(index32(index), other);
                        }
                    }
                }
    }
    std::sort(result.begin(), result.end());
    const auto temporaryBytes =
        grid.capacity() * sizeof(IndexedCell) + locations.capacity() * sizeof(Cell);
    std::vector<IndexedCell>().swap(grid);
    std::vector<Cell>().swap(locations);
    budget.release(temporaryBytes);
    return result;
}

struct PlannedPair {
    PanelPair pair;
    PairQuadratureRule rule;
};
} // namespace

InteractionPlan makeInteractionPlan(const Space &space, const PairQuadratureSettings &settings,
                                    double maximumWavenumber, std::stop_token stop,
                                    unsigned maximumThreads, std::size_t maximumHostBytes,
                                    std::span<const std::uint32_t> testPanels) {
    cancelled(stop);
    if (space.panels.empty() || space.vertices.empty() || !std::isfinite(maximumWavenumber) ||
        maximumWavenumber < 0 || !std::isfinite(settings.nearDistanceRatio) ||
        settings.nearDistanceRatio < 0 || !std::isfinite(settings.relativeTolerance) ||
        settings.relativeTolerance <= 0 || settings.maximumSubdivisionDepth > 16 ||
        (settings.regularOrder != 2 && settings.regularOrder != 4) || settings.singularOrder < 2 ||
        settings.singularOrder > 20 || space.incidenceOffsets.size() != space.vertices.size() + 1 ||
        maximumHostBytes < 4)
        throw std::invalid_argument("Invalid BEM interaction-plan input or quadrature settings.");
    index32(space.panels.size());
    std::vector<bool> active(space.panels.size(), testPanels.empty());
    for (auto panel : testPanels) {
        if (panel >= active.size() || active[panel])
            throw std::invalid_argument("BEM target-panel list has invalid or duplicate indices.");
        active[panel] = true;
    }
    InteractionPlan result;
    // Independent pools bound simultaneously live classification, worker and
    // output payloads, including old+new capacities during vector growth.
    AdaptiveQuadratureWorkspace fixedBudget(maximumHostBytes / 4),
        outputBudget(maximumHostBytes / 2);
    std::vector<RegularPanelBounds> bounds;
    reserveTracked(bounds, space.panels.size(), fixedBudget);
    // Use the qualifier's multiplication/comparison order, including at the
    // rounded boundary; dividing the threshold can disagree by one ulp.
    bool allPairs = !((1e-6 + 4 * 3.5e-8) <= .2 * settings.relativeTolerance);
    for (const auto &panel : space.panels) {
        bounds.push_back(makeRegularPanelBounds(panel));
        allPairs |=
            !bounds.back().qualifiedShape || maximumWavenumber * panel.maximumEdgeLength > 2;
    }
    reserveTracked(result.nearOffsets, space.panels.size() + 1, fixedBudget);
    result.nearOffsets.resize(space.panels.size() + 1, 0);
    for (std::size_t test = 0; test < space.panels.size(); ++test) {
        cancelled(stop);
        if (!active[test])
            continue;
        std::vector<std::uint32_t> adjacent;
        for (const auto node : space.panels[test].nodes)
            for (auto offset = space.incidenceOffsets[node];
                 offset < space.incidenceOffsets[node + 1]; ++offset) {
                reserveTracked(adjacent, adjacent.size() + 1, fixedBudget);
                adjacent.push_back(space.incidences[offset].panel);
            }
        std::sort(adjacent.begin(), adjacent.end());
        adjacent.erase(std::unique(adjacent.begin(), adjacent.end()), adjacent.end());
        for (const auto trial : adjacent) {
            reserveTracked(result.touching, result.touching.size() + 1, fixedBudget);
            result.touching.push_back(
                {index32(test), trial, classifyPair(space.panels[test], space.panels[trial])});
        }
        const auto adjacentBytes = adjacent.capacity() * sizeof(std::uint32_t);
        std::vector<std::uint32_t>().swap(adjacent);
        fixedBudget.release(adjacentBytes);
    }
    // With qualified shapes, k*h<=2 and a sufficient arithmetic budget, every
    // unqualified pair has centroid-ball gap<h. Its vertex-zero distance is
    // <(1+8/3)*h<4*h, so ratio2 is a conservative spatial candidate search.
    const auto pairs =
        allPairs ? std::vector<PanelPair>{}
                 : candidates(space, std::max(2., settings.nearDistanceRatio), stop, fixedBudget, active);
    result.candidatePairs = pairs.size();
    unsigned count =
        maximumThreads ? maximumThreads : std::max(1U, std::thread::hardware_concurrency());
    count = static_cast<unsigned>(
        std::min<std::size_t>({count, 8, allPairs ? space.panels.size() : pairs.size()}));
    if (!count)
        return result;
    std::vector<std::vector<PlannedPair>> perWorker(count);
    std::atomic<std::size_t> streamedCandidates{};
    std::atomic<std::size_t> next{};
    std::atomic<bool> failed{};
    std::exception_ptr failure;
    std::mutex failureMutex;
    std::stop_source pairStop;
    std::stop_callback propagateStop(stop, [&] { pairStop.request_stop(); });
    std::vector<std::jthread> workers;
    workers.reserve(count);
    for (unsigned worker = 0; worker < count; ++worker)
        workers.emplace_back([&, worker](std::stop_token localStop) {
            try {
                AdaptiveQuadratureWorkspace workspace(maximumHostBytes / 4 / count);
                const auto process = [&](PanelPair pair) {
                    cancelled(pairStop.get_token());
                    if (!active[pair.first] && !active[pair.second])
                        return;
                    const auto &test = space.panels[pair.first], &trial = space.panels[pair.second];
                    if (classifyPair(test, trial).topology != PairTopology::disjoint)
                        return;
                    if (allPairs)
                        ++streamedCandidates;
                    const auto degree = qualifiedRegularDegree(
                        bounds[pair.first], bounds[pair.second], maximumWavenumber,
                        settings.relativeTolerance, settings.regularOrder);
                    const double h = std::max(test.maximumEdgeLength, trial.maximumEdgeLength),
                                 limit = bounds[pair.first].radius + bounds[pair.second].radius +
                                         settings.nearDistanceRatio * h,
                                 dx = test.centroid.x - trial.centroid.x,
                                 dy = test.centroid.y - trial.centroid.y,
                                 dz = test.centroid.z - trial.centroid.z;
                    const bool geometricNear = degree && settings.nearDistanceRatio > 0 &&
                                               dx * dx + dy * dy + dz * dz < limit * limit &&
                                               isNearPair(test, trial, settings.nearDistanceRatio);
                    if (degree && !geometricNear)
                        return;
                    try {
                        PairQuadratureRule retained;
                        {
                            AdaptivePairQuadratureState state(test, trial, settings,
                                                              maximumWavenumber, workspace,
                                                              pairStop.get_token());
                            while (!state.done()) {
                                const auto tasks = state.requests();
                                std::array<PairPatchEstimate, 4> estimates;
                                for (std::size_t i = 0; i < tasks.size(); ++i) {
                                    cancelled(pairStop.get_token());
                                    estimates[i] = estimatePairPatchCpu(test, trial, tasks[i],
                                                                        maximumWavenumber, pairStop.get_token());
                                }
                                state.commit(std::span(estimates).first(tasks.size()));
                            }
                            const auto &rule = state.result();
                            // The copy is intentional: both old and retained
                            // capacities are charged BEFORE allocating it.
                            reserveTracked(retained.patches, rule.patches.size(), workspace);
                            retained.patches.assign(rule.patches.begin(), rule.patches.end());
                            retained.topology = rule.topology;
                            retained.estimatedRelativeError = rule.estimatedRelativeError;
                            retained.maximumDepth = rule.maximumDepth;
                            retained.peakWorkspaceBytes = rule.peakWorkspaceBytes;
                        }
                        reserveTracked(perWorker[worker], perWorker[worker].size() + 1, workspace);
                        perWorker[worker].push_back({pair, std::move(retained)});
                    } catch (const std::exception &error) {
                        if (pairStop.stop_requested())
                            throw;
                        std::ostringstream context;
                        context << "Near pair panels " << pair.first << '/' << pair.second
                                << " (faces " << test.faceId << '/' << trial.faceId
                                << ", k=" << maximumWavenumber
                                << ", gap=" << panelDistance(test, trial)
                                << " m, depth=" << settings.maximumSubdivisionDepth
                                << "): " << error.what();
                        throw std::runtime_error(context.str());
                    }
                };
                for (;;) {
                    cancelled(stop);
                    if (failed.load(std::memory_order_relaxed) || localStop.stop_requested())
                        return;
                    const auto index = next.fetch_add(1, std::memory_order_relaxed);
                    if (index >= (allPairs ? space.panels.size() : pairs.size()))
                        return;
                    if (allPairs) {
                        // Only a row counter is shared. No all-pair job list is
                        // allocated even for needles, high kh or strict budgets.
                        for (std::size_t source = index + 1; source < space.panels.size(); ++source)
                            process({index32(index), index32(source)});
                    } else
                        process(pairs[index]);
                }
            } catch (...) {
                std::lock_guard lock(failureMutex);
                if (!failure)
                    failure = std::current_exception();
                failed.store(true, std::memory_order_relaxed);
                pairStop.request_stop();
            }
        });
    for (auto &worker : workers)
        worker.join();
    if (failure)
        std::rethrow_exception(failure);
    if (allPairs)
        result.candidatePairs = streamedCandidates.load();
    cancelled(stop);
    std::vector<PlannedPair> planned;
    for (auto &worker : perWorker)
        for (auto &pair : worker) {
            reserveTracked(planned, planned.size() + 1, outputBudget);
            planned.push_back(std::move(pair));
        }
    std::sort(planned.begin(), planned.end(),
              [](const auto &a, const auto &b) { return a.pair < b.pair; });
    for (const auto &entry : planned) {
        cancelled(stop);
        const auto first = index32(result.patches.size()),
                   size = index32(entry.rule.patches.size());
        index32(result.patches.size() + entry.rule.patches.size());
        reserveTracked(result.patches, result.patches.size() + entry.rule.patches.size(),
                       outputBudget);
        result.patches.insert(result.patches.end(), entry.rule.patches.begin(),
                              entry.rule.patches.end());
        reserveTracked(result.near, result.near.size() + 2, outputBudget);
        if (active[entry.pair.first])
            result.near.push_back({entry.pair.first, entry.pair.second, first, size, 0});
        if (active[entry.pair.second])
            result.near.push_back({entry.pair.second, entry.pair.first, first, size, 1});
        result.maximumDepth = std::max(result.maximumDepth, entry.rule.maximumDepth);
        result.maximumEstimatedRelativeError =
            std::max(result.maximumEstimatedRelativeError, entry.rule.estimatedRelativeError);
    }
    index32(result.near.size());
    reserveTracked(result.nearTrials, result.near.size(), outputBudget);
    std::sort(result.near.begin(), result.near.end(), [](const auto &a, const auto &b) {
        return std::pair{a.test, a.trial} < std::pair{b.test, b.trial};
    });
    for (const auto &pair : result.near) {
        ++result.nearOffsets[pair.test + 1];
        result.nearTrials.push_back(pair.trial);
    }
    for (std::size_t row = 1; row < result.nearOffsets.size(); ++row)
        result.nearOffsets[row] += result.nearOffsets[row - 1];
    return result;
}
} // namespace hornsculpt::galerkin
