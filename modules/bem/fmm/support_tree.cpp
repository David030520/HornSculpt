#include "support_tree.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>

namespace hornsculpt::wideband {
namespace {
SupportBounds emptyBounds() {
    const auto inf = std::numeric_limits<double>::infinity();
    return {{inf, inf, inf}, {-inf, -inf, -inf}, 0};
}
void include(SupportBounds &a, const SupportBounds &b) {
    a.lower = {std::min(a.lower.x, b.lower.x), std::min(a.lower.y, b.lower.y),
               std::min(a.lower.z, b.lower.z)};
    a.upper = {std::max(a.upper.x, b.upper.x), std::max(a.upper.y, b.upper.y),
               std::max(a.upper.z, b.upper.z)};
    a.maximumPanelEdge = std::max(a.maximumPanelEdge, b.maximumPanelEdge);
}
double gap(const SupportBounds &a, const SupportBounds &b) {
    const double x = std::max({0., a.lower.x - b.upper.x, b.lower.x - a.upper.x});
    const double y = std::max({0., a.lower.y - b.upper.y, b.lower.y - a.upper.y});
    const double z = std::max({0., a.lower.z - b.upper.z, b.lower.z - a.upper.z});
    return std::hypot(x, y, z);
}
template <class T> std::size_t bytes(const std::vector<T> &v) {
    return v.capacity() * sizeof(T);
}

class Builder {
  public:
    Builder(const galerkin::Space &s, SupportTreeSettings config, std::stop_token token)
        : space(s), settings(config), stop(token) {}
    SupportTree build() {
        check();
        tree.separationRatio = settings.separationRatio;
        tree.nearDistanceRatio = settings.nearDistanceRatio;
        if (space.vertices.empty() || space.panels.empty() ||
            space.incidenceOffsets.size() != space.vertices.size() + 1 ||
            space.vertices.size() >= noSupportNode || space.panels.size() >= noSupportNode)
            throw std::invalid_argument("Wideband support geometry is empty or invalid.");
        tree.fluxSupports.reserve(space.panels.size());
        auto global = emptyBounds();
        for (const auto &panel : space.panels) {
            auto b = emptyBounds();
            for (const auto p : panel.vertices) {
                if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
                    throw std::invalid_argument("Wideband support coordinates are not finite.");
                include(b, {p, p, panel.maximumEdgeLength});
            }
            tree.fluxSupports.push_back({panel.centroid, b});
            include(global, b);
        }
        tree.pressureSupports.reserve(space.vertices.size());
        for (std::size_t i = 0; i < space.vertices.size(); ++i) {
            check();
            auto b = emptyBounds();
            const auto begin = space.incidenceOffsets[i], end = space.incidenceOffsets[i + 1];
            if (begin >= end || end > space.incidences.size())
                throw std::invalid_argument("Wideband pressure support has invalid incidences.");
            for (auto j = begin; j < end; ++j) {
                const auto panel = space.incidences[j].panel;
                if (panel >= tree.fluxSupports.size())
                    throw std::invalid_argument("Wideband support panel index is invalid.");
                include(b, tree.fluxSupports[panel].bounds);
            }
            tree.pressureSupports.push_back({space.vertices[i], b});
        }
        tree.pressureOrder.resize(space.vertices.size());
        tree.fluxOrder.resize(space.panels.size());
        tree.pressureLeaf.resize(space.vertices.size());
        tree.fluxLeaf.resize(space.panels.size());
        std::iota(tree.pressureOrder.begin(), tree.pressureOrder.end(), 0u);
        std::iota(tree.fluxOrder.begin(), tree.fluxOrder.end(), 0u);
        const Vec3 center{std::midpoint(global.lower.x, global.upper.x),
                          std::midpoint(global.lower.y, global.upper.y),
                          std::midpoint(global.lower.z, global.upper.z)};
        const double half =
            .5 * std::max({global.upper.x - global.lower.x, global.upper.y - global.lower.y,
                           global.upper.z - global.lower.z});
        if (!(half > 0) || !std::isfinite(half))
            throw std::invalid_argument("Wideband root bounds are invalid.");
        node(center, std::nextafter(half, std::numeric_limits<double>::infinity()), {}, 0,
             noSupportNode, 0, static_cast<std::uint32_t>(space.vertices.size()), 0,
             static_cast<std::uint32_t>(space.panels.size()));
        nestExpansionRadii();
        partition(0, 0, SupportKind::pressure, tree.pressurePressure);
        partition(0, 0, SupportKind::flux, tree.pressureFlux);
        activate(tree.pressurePressure, pressureOutgoing);
        activate(tree.pressureFlux, fluxOutgoing);
        countOffsets();
        for (const auto &n : tree.nodes) {
            auto &level = tree.levels[n.level];
            if (n.activeChannels & pressureOutgoing)
                ++level.pressureOutgoingNodes;
            if (n.activeChannels & fluxOutgoing)
                ++level.fluxOutgoingNodes;
            if (n.activeChannels & pressureIncoming)
                ++level.pressureIncomingNodes;
            if (n.activeChannels & (pressureOutgoing | pressureIncoming))
                level.maximumActiveSupportRadius =
                    std::max(level.maximumActiveSupportRadius, n.pressureRadius);
            if (n.activeChannels & fluxOutgoing)
                level.maximumActiveSupportRadius =
                    std::max(level.maximumActiveSupportRadius, n.fluxRadius);
        }
        return std::move(tree);
    }

  private:
    void nestExpansionRadii() {
        for (std::size_t i = tree.nodes.size(); i-- > 1;) {
            check();
            const auto &child = tree.nodes[i];
            auto &parent = tree.nodes[child.parent];
            const auto distance =
                std::hypot(child.center.x - parent.center.x, child.center.y - parent.center.y,
                           child.center.z - parent.center.z);
            if (child.pressureCount)
                parent.pressureRadius = std::max(
                    parent.pressureRadius, std::nextafter(distance + child.pressureRadius,
                                                          std::numeric_limits<double>::infinity()));
            if (child.fluxCount)
                parent.fluxRadius = std::max(
                    parent.fluxRadius, std::nextafter(distance + child.fluxRadius,
                                                      std::numeric_limits<double>::infinity()));
        }
        for (auto &level : tree.levels)
            level.maximumSupportRadius = 0;
        for (const auto &node : tree.nodes)
            tree.levels[node.level].maximumSupportRadius =
                std::max({tree.levels[node.level].maximumSupportRadius, node.pressureRadius,
                          node.fluxRadius});
    }
    void countOffsets() {
        using Offset = std::array<std::int64_t, 3>;
        std::vector<std::set<Offset>> pressure(tree.levels.size()), flux(tree.levels.size());
        const auto add = [&](const SupportPartition &partition, auto &sets) {
            for (const auto pair : partition.farInteractions) {
                check();
                const auto &t = tree.nodes[pair.target], &s = tree.nodes[pair.source];
                Offset offset{};
                for (unsigned a = 0; a < 3; ++a)
                    offset[a] = std::int64_t(t.coordinates[a]) - s.coordinates[a];
                sets[t.level].insert(offset);
            }
        };
        add(tree.pressurePressure, pressure);
        add(tree.pressureFlux, flux);
        for (std::size_t l = 0; l < tree.levels.size(); ++l) {
            auto &level = tree.levels[l];
            level.pressurePressureOffsets = pressure[l].size();
            level.pressureFluxOffsets = flux[l].size();
            const auto distances = [](const auto &offsets) {
                std::set<std::uint64_t> result;
                for (const auto &v : offsets)
                    result.insert(
                        static_cast<std::uint64_t>(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
                return result;
            };
            const auto polarClasses = [](const auto &offsets) {
                std::set<std::array<std::uint64_t, 3>> result;
                for (const auto &v : offsets) {
                    const auto z2 = static_cast<std::uint64_t>(v[2] * v[2]);
                    const auto norm2 = static_cast<std::uint64_t>(v[0] * v[0] + v[1] * v[1]) + z2;
                    const auto divisor = std::gcd(z2, norm2);
                    result.insert({v[2] < 0   ? 0u
                                   : v[2] > 0 ? 2u
                                              : 1u,
                                   z2 / divisor, norm2 / divisor});
                }
                return result;
            };
            auto pressurePolar = polarClasses(pressure[l]);
            const auto fluxPolar = polarClasses(flux[l]);
            level.pressurePressurePolarClasses = pressurePolar.size();
            level.pressureFluxPolarClasses = fluxPolar.size();
            pressurePolar.insert(fluxPolar.begin(), fluxPolar.end());
            level.unionPolarClasses = pressurePolar.size();
            auto pressureDistances = distances(pressure[l]);
            const auto fluxDistances = distances(flux[l]);
            level.pressurePressureDistanceClasses = pressureDistances.size();
            level.pressureFluxDistanceClasses = fluxDistances.size();
            pressureDistances.insert(fluxDistances.begin(), fluxDistances.end());
            level.squaredOffsetClasses.assign(pressureDistances.begin(), pressureDistances.end());
            pressure[l].insert(flux[l].begin(), flux[l].end());
            level.unionOffsets = pressure[l].size();
        }
    }
    void activateSubtree(std::uint32_t index, unsigned channel) {
        auto &n = tree.nodes[index];
        const bool nonempty = channel == fluxOutgoing ? n.fluxCount != 0 : n.pressureCount != 0;
        if (!nonempty || (n.activeChannels & channel))
            return;
        n.activeChannels |= channel;
        for (const auto child : n.children)
            if (child != noSupportNode)
                activateSubtree(child, channel);
    }
    void activate(SupportPartition &partition, unsigned outgoing) {
        std::map<std::pair<std::uint32_t, std::uint32_t>, SupportLevelPair> pairs;
        for (const auto pair : partition.farInteractions) {
            check();
            activateSubtree(pair.target, pressureIncoming);
            activateSubtree(pair.source, outgoing);
            auto &level = tree.levels[tree.nodes[pair.target].level];
            if (outgoing == pressureOutgoing)
                ++level.pressurePressureFarInteractions;
            else
                ++level.pressureFluxFarInteractions;
            const auto &t = tree.nodes[pair.target], &s = tree.nodes[pair.source];
            auto &statistics = pairs[{t.level, s.level}];
            statistics.targetLevel = t.level;
            statistics.sourceLevel = s.level;
            ++statistics.interactions;
            statistics.representedPairs +=
                std::uint64_t(t.pressureCount) *
                (outgoing == pressureOutgoing ? s.pressureCount : s.fluxCount);
        }
        for (const auto &entry : pairs)
            partition.farLevelPairs.push_back(entry.second);
    }
    void check() const {
        if (stop.stop_requested())
            throw std::runtime_error("Wideband support construction cancelled.");
    }
    std::array<std::uint32_t, 9> split(std::vector<std::uint32_t> &order,
                                       const std::vector<Support> &supports, std::uint32_t first,
                                       std::uint32_t count, Vec3 center) {
        const auto octant = [&](std::uint32_t id) {
            const auto p = supports[id].anchor;
            return unsigned(p.x >= center.x) | (unsigned(p.y >= center.y) << 1) |
                   (unsigned(p.z >= center.z) << 2);
        };
        std::array<std::uint32_t, 9> offsets{};
        for (std::size_t i = first; i < std::size_t(first) + count; ++i)
            ++offsets[octant(order[i]) + 1];
        offsets[0] = first;
        for (unsigned i = 1; i < 9; ++i)
            offsets[i] += offsets[i - 1];
        auto cursor = offsets;
        std::vector<std::uint32_t> sorted(count);
        for (std::size_t i = first; i < std::size_t(first) + count; ++i)
            sorted[cursor[octant(order[i])]++ - first] = order[i];
        std::copy(sorted.begin(), sorted.end(), order.begin() + first);
        return offsets;
    }
    double completeSupportRadius(Vec3 center, std::uint32_t first, std::uint32_t count,
                                 SupportKind kind) const {
        double result = 0;
        const auto includePanel = [&](std::uint32_t panel) {
            for (const auto p : space.panels[panel].vertices)
                result =
                    std::max(result, std::hypot(p.x - center.x, p.y - center.y, p.z - center.z));
        };
        for (std::size_t i = first; i < std::size_t(first) + count; ++i) {
            check();
            if (kind == SupportKind::flux)
                includePanel(tree.fluxOrder[i]);
            else {
                const auto vertex = tree.pressureOrder[i];
                for (auto j = space.incidenceOffsets[vertex];
                     j < space.incidenceOffsets[std::size_t(vertex) + 1]; ++j)
                    includePanel(space.incidences[j].panel);
            }
        }
        // A triangle is contained in the ball enclosing its three vertices.
        // Artificial AABB corners are not part of the P1/DP0 support and can
        // substantially inflate the electrical radius of tilted/curved leaves.
        return count ? std::nextafter(result, std::numeric_limits<double>::infinity()) : 0;
    }
    std::uint32_t node(Vec3 center, double half, std::array<std::uint32_t, 3> coords,
                       std::uint32_t level, std::uint32_t parent, std::uint32_t fp,
                       std::uint32_t np, std::uint32_t fq, std::uint32_t nq) {
        check();
        if (tree.nodes.size() >= noSupportNode)
            throw std::overflow_error("Wideband node index overflow.");
        SupportNode value;
        value.center = center;
        value.halfWidth = half;
        value.coordinates = coords;
        value.level = level;
        value.parent = parent;
        value.firstPressure = fp;
        value.pressureCount = np;
        value.firstFlux = fq;
        value.fluxCount = nq;
        value.pressureBounds = emptyBounds();
        value.fluxBounds = emptyBounds();
        for (std::size_t i = fp; i < std::size_t(fp) + np; ++i)
            include(value.pressureBounds, tree.pressureSupports[tree.pressureOrder[i]].bounds);
        for (std::size_t i = fq; i < std::size_t(fq) + nq; ++i)
            include(value.fluxBounds, tree.fluxSupports[tree.fluxOrder[i]].bounds);
        const bool leaf = (np <= settings.leafCapacity && nq <= settings.leafCapacity) ||
                          level == settings.maximumDepth;
        if (leaf) {
            value.pressureRadius = completeSupportRadius(center, fp, np, SupportKind::pressure);
            value.fluxRadius = completeSupportRadius(center, fq, nq, SupportKind::flux);
        }
        const auto id = static_cast<std::uint32_t>(tree.nodes.size());
        tree.nodes.push_back(value);
        if (tree.levels.size() <= level)
            tree.levels.resize(std::size_t(level) + 1);
        tree.levels[level].nodes.push_back(id);
        tree.levels[level].maximumSupportRadius = std::max(
            {tree.levels[level].maximumSupportRadius, value.pressureRadius, value.fluxRadius});
        if (leaf) {
            for (std::size_t i = fp; i < std::size_t(fp) + np; ++i)
                tree.pressureLeaf[tree.pressureOrder[i]] = id;
            for (std::size_t i = fq; i < std::size_t(fq) + nq; ++i)
                tree.fluxLeaf[tree.fluxOrder[i]] = id;
            return id;
        }
        const auto po = split(tree.pressureOrder, tree.pressureSupports, fp, np, center);
        const auto qo = split(tree.fluxOrder, tree.fluxSupports, fq, nq, center);
        for (unsigned oct = 0; oct < 8; ++oct) {
            if (po[oct] == po[oct + 1] && qo[oct] == qo[oct + 1])
                continue;
            const double childHalf = half * .5;
            const Vec3 childCenter{center.x + ((oct & 1) ? childHalf : -childHalf),
                                   center.y + ((oct & 2) ? childHalf : -childHalf),
                                   center.z + ((oct & 4) ? childHalf : -childHalf)};
            auto childCoords = coords;
            for (unsigned a = 0; a < 3; ++a)
                childCoords[a] = 2 * coords[a] + ((oct >> a) & 1);
            const auto child = node(childCenter, childHalf, childCoords, level + 1, id, po[oct],
                                    po[oct + 1] - po[oct], qo[oct], qo[oct + 1] - qo[oct]);
            tree.nodes[id].children[oct] = child;
        }
        return id;
    }
    void partition(std::uint32_t ti, std::uint32_t si, SupportKind kind, SupportPartition &out) {
        check();
        const auto &t = tree.nodes[ti], &s = tree.nodes[si];
        const auto ns = kind == SupportKind::pressure ? s.pressureCount : s.fluxCount;
        if (!t.pressureCount || !ns)
            return;
        const auto &sb = kind == SupportKind::pressure ? s.pressureBounds : s.fluxBounds;
        const auto sr = kind == SupportKind::pressure ? s.pressureRadius : s.fluxRadius;
        const double distance =
            std::hypot(t.center.x - s.center.x, t.center.y - s.center.y, t.center.z - s.center.z);
        const double supportGap = gap(t.pressureBounds, sb);
        const bool far =
            t.level == s.level && supportGap > 0 &&
            supportGap > settings.nearDistanceRatio *
                             std::max(t.pressureBounds.maximumPanelEdge, sb.maximumPanelEdge) &&
            distance > settings.separationRatio * (t.pressureRadius + sr);
        if (far || (t.leaf() && s.leaf())) {
            if (interactions == settings.maximumInteractions)
                throw std::runtime_error("Wideband interaction memory budget exceeded.");
            ++interactions;
            const auto count = std::uint64_t(t.pressureCount) * ns;
            if (far) {
                out.farInteractions.push_back({ti, si});
                out.representedFarPairs += count;
            } else {
                out.nearInteractions.push_back({ti, si});
                out.representedNearPairs += count;
            }
            return;
        }
        if (!t.leaf() && (s.leaf() || t.level < s.level))
            for (auto child : t.children) {
                if (child != noSupportNode)
                    partition(child, si, kind, out);
            }
        // A coarse leaf cannot become a far endpoint against finer nodes.
        // Unequal-level leaf products remain exact near work, with no implicit
        // angular resampling or duplicated virtual support ownership.
        else if (!s.leaf() && (t.leaf() || s.level < t.level))
            for (auto child : s.children) {
                if (child != noSupportNode)
                    partition(ti, child, kind, out);
            }
        else
            for (auto tc : t.children)
                if (tc != noSupportNode)
                    for (auto sc : s.children)
                        if (sc != noSupportNode)
                            partition(tc, sc, kind, out);
    }
    const galerkin::Space &space;
    SupportTreeSettings settings;
    std::stop_token stop;
    SupportTree tree;
    std::size_t interactions{};
};
} // namespace

bool SupportNode::leaf() const noexcept {
    return std::all_of(children.begin(), children.end(),
                       [](auto id) { return id == noSupportNode; });
}
std::size_t SupportTree::allocatedBytes() const noexcept {
    std::size_t result =
        bytes(pressureSupports) + bytes(fluxSupports) + bytes(pressureOrder) + bytes(fluxOrder) +
        bytes(pressureLeaf) + bytes(fluxLeaf) + bytes(nodes) + bytes(levels) +
        bytes(pressurePressure.nearInteractions) + bytes(pressurePressure.farInteractions) +
        bytes(pressurePressure.farLevelPairs) + bytes(pressureFlux.nearInteractions) +
        bytes(pressureFlux.farInteractions) + bytes(pressureFlux.farLevelPairs);
    for (const auto &level : levels)
        result += bytes(level.nodes) + bytes(level.squaredOffsetClasses);
    return result;
}
SupportTree makeSupportTree(const galerkin::Space &space, SupportTreeSettings settings,
                            std::stop_token stop) {
    if (!settings.leafCapacity || settings.maximumDepth > 30 ||
        !std::isfinite(settings.separationRatio) || settings.separationRatio <= 1 ||
        !std::isfinite(settings.nearDistanceRatio) || settings.nearDistanceRatio < 0)
        throw std::invalid_argument("Wideband support tree settings are invalid.");
    return Builder(space, settings, stop).build();
}
} // namespace hornsculpt::wideband
