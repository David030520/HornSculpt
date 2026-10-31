#include "interaction_plan.hpp"
#include "space.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
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
    bool operator==(const Cell &) const = default;
};
struct CellHash {
    std::size_t operator()(const Cell &cell) const {
        auto hash = std::hash<std::int64_t>{}(cell.x);
        for (const auto coordinate : {cell.y, cell.z})
            hash ^= std::hash<std::int64_t>{}(coordinate) + 0x9e3779b97f4a7c15ULL + (hash << 6) +
                    (hash >> 2);
        return hash;
    }
};
using PanelPair = std::pair<std::uint32_t, std::uint32_t>;

std::vector<PanelPair> candidates(const Space &space, double ratio, std::stop_token stop,
                                  const std::vector<bool> &active) {
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
    std::unordered_map<Cell, std::vector<std::uint32_t>, CellHash> grid;
    grid.reserve(space.panels.size());
    std::vector<Cell> locations;
    locations.reserve(space.panels.size());
    for (std::size_t index = 0; index < space.panels.size(); ++index) {
        cancelled(stop);
        const auto point = space.panels[index].vertices[0];
        const Cell cell{coordinate(point.x, origin.x), coordinate(point.y, origin.y),
                        coordinate(point.z, origin.z)};
        locations.push_back(cell);
        grid[cell].push_back(index32(index));
    }
    std::vector<PanelPair> result;
    for (std::size_t index = 0; index < space.panels.size(); ++index) {
        cancelled(stop);
        const auto cell = locations[index];
        for (int x = -1; x <= 1; ++x)
            for (int y = -1; y <= 1; ++y)
                for (int z = -1; z <= 1; ++z) {
                    const auto found = grid.find({cell.x + x, cell.y + y, cell.z + z});
                    if (found == grid.end())
                        continue;
                    for (const auto other : found->second) {
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
                                .topology == PairTopology::disjoint)
                            result.emplace_back(index32(index), other);
                    }
                }
    }
    std::sort(result.begin(), result.end());
    return result;
}

struct PlannedPair {
    PanelPair pair;
    PairQuadratureRule rule;
};
} // namespace

InteractionPlan makeInteractionPlan(const Space &space, const PairQuadratureSettings &settings,
                                    double maximumWavenumber, std::stop_token stop,
                                    unsigned maximumThreads,
                                    std::span<const std::uint32_t> testPanels) {
    cancelled(stop);
    if (space.panels.empty() || space.vertices.empty() || !std::isfinite(maximumWavenumber) ||
        maximumWavenumber < 0 || !std::isfinite(settings.nearDistanceRatio) ||
        settings.nearDistanceRatio < 0 || !std::isfinite(settings.relativeTolerance) ||
        settings.relativeTolerance <= 0 || settings.maximumSubdivisionDepth > 16 ||
        (settings.regularOrder != 2 && settings.regularOrder != 4) || settings.singularOrder < 2 ||
        settings.singularOrder > 20 || space.incidenceOffsets.size() != space.vertices.size() + 1)
        throw std::invalid_argument("Invalid BEM interaction-plan input or quadrature settings.");
    index32(space.panels.size());
    std::vector<bool> active(space.panels.size(), testPanels.empty());
    for (auto panel : testPanels) {
        if (panel >= active.size() || active[panel])
            throw std::invalid_argument("BEM target-panel list has invalid or duplicate indices.");
        active[panel] = true;
    }
    InteractionPlan result;
    result.nearOffsets.resize(space.panels.size() + 1, 0);
    for (std::size_t test = 0; test < space.panels.size(); ++test) {
        cancelled(stop);
        if (!active[test])
            continue;
        std::vector<std::uint32_t> adjacent;
        for (const auto node : space.panels[test].nodes)
            for (auto offset = space.incidenceOffsets[node];
                 offset < space.incidenceOffsets[node + 1]; ++offset)
                adjacent.push_back(space.incidences[offset].panel);
        std::sort(adjacent.begin(), adjacent.end());
        adjacent.erase(std::unique(adjacent.begin(), adjacent.end()), adjacent.end());
        for (const auto trial : adjacent)
            result.touching.push_back(
                {index32(test), trial, classifyPair(space.panels[test], space.panels[trial])});
    }
    if (settings.nearDistanceRatio == 0)
        return result;
    const auto pairs = candidates(space, settings.nearDistanceRatio, stop, active);
    result.candidatePairs = pairs.size();
    unsigned count =
        maximumThreads ? maximumThreads : std::max(1U, std::thread::hardware_concurrency());
    count = static_cast<unsigned>(std::min<std::size_t>({count, 8, pairs.size()}));
    if (!count)
        return result;
    std::vector<std::vector<PlannedPair>> perWorker(count);
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
                for (;;) {
                    cancelled(stop);
                    if (failed.load(std::memory_order_relaxed) || localStop.stop_requested())
                        return;
                    const auto index = next.fetch_add(1, std::memory_order_relaxed);
                    if (index >= pairs.size())
                        return;
                    const auto pair = pairs[index];
                    const auto &test = space.panels[pair.first], &trial = space.panels[pair.second];
                    if (!isNearPair(test, trial, settings.nearDistanceRatio))
                        continue;
                    PairQuadratureRule rule;
                    try {
                        rule = makePairQuadrature(test, trial, settings, maximumWavenumber,
                                                  pairStop.get_token());
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
                    if (rule.patches.empty())
                        throw std::logic_error("Near interaction has no adaptive patch plan.");
                    perWorker[worker].push_back({pair, std::move(rule)});
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
    cancelled(stop);
    std::vector<PlannedPair> planned;
    for (auto &worker : perWorker)
        for (auto &pair : worker)
            planned.push_back(std::move(pair));
    std::sort(planned.begin(), planned.end(),
              [](const auto &a, const auto &b) { return a.pair < b.pair; });
    for (const auto &entry : planned) {
        cancelled(stop);
        const auto first = index32(result.patches.size()),
                   size = index32(entry.rule.patches.size());
        index32(result.patches.size() + entry.rule.patches.size());
        result.patches.insert(result.patches.end(), entry.rule.patches.begin(),
                              entry.rule.patches.end());
        if (active[entry.pair.first])
            result.near.push_back({entry.pair.first, entry.pair.second, first, size, 0});
        if (active[entry.pair.second])
            result.near.push_back({entry.pair.second, entry.pair.first, first, size, 1});
        result.maximumDepth = std::max(result.maximumDepth, entry.rule.maximumDepth);
        result.maximumEstimatedRelativeError =
            std::max(result.maximumEstimatedRelativeError, entry.rule.estimatedRelativeError);
    }
    index32(result.near.size());
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
