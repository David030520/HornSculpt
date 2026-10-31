#include "fixture_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <fmm/support_tree.hpp>
#include <galerkin/pair_quadrature.hpp>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <stop_token>

using namespace hornsculpt;
using namespace hornsculpt::wideband;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
template <class F> void rejects(F f) {
    bool failed = false;
    try {
        f();
    } catch (const std::exception &) {
        failed = true;
    }
    require(failed, "Invalid or cancelled support construction was accepted.");
}
SurfaceMesh separated() {
    auto mesh = testing::sphere(2);
    const auto other = testing::sphere(1);
    const auto first = static_cast<std::uint32_t>(mesh.vertices.size());
    for (auto p : other.vertices)
        mesh.vertices.push_back({p.x + 2, p.y + .03, p.z - .02});
    for (auto t : other.triangles) {
        for (auto &v : t.vertices)
            v += first;
        mesh.triangles.push_back(t);
    }
    return mesh;
}
bool contains(const SupportBounds &b, Vec3 p) {
    constexpr double eps = 1e-13;
    return p.x >= b.lower.x - eps && p.x <= b.upper.x + eps && p.y >= b.lower.y - eps &&
           p.y <= b.upper.y + eps && p.z >= b.lower.z - eps && p.z <= b.upper.z + eps;
}
void hierarchy(const galerkin::Space &space, const SupportTree &tree) {
    std::vector<unsigned> pm(space.vertices.size()), qm(space.panels.size()),
        levels(tree.nodes.size());
    bool crosses = false;
    for (std::size_t l = 0; l < tree.levels.size(); ++l)
        for (auto ni : tree.levels[l].nodes) {
            require(ni < tree.nodes.size() && tree.nodes[ni].level == l,
                    "Invalid support level inventory.");
            ++levels[ni];
            require(tree.levels[l].maximumSupportRadius >=
                        std::max(tree.nodes[ni].pressureRadius, tree.nodes[ni].fluxRadius),
                    "Level angular radius fails to bound all complete supports.");
        }
    for (std::size_t ni = 0; ni < tree.nodes.size(); ++ni) {
        const auto &n = tree.nodes[ni];
        require(levels[ni] == 1, "Node lost or duplicated in support levels.");
        if (ni) {
            require(n.parent < ni && n.level == tree.nodes[n.parent].level + 1,
                    "Invalid dyadic parent.");
            const auto &p = tree.nodes[n.parent];
            const auto centerDistance = std::hypot(n.center.x - p.center.x, n.center.y - p.center.y,
                                                   n.center.z - p.center.z);
            require((!n.pressureCount || p.pressureRadius >= centerDistance + n.pressureRadius) &&
                        (!n.fluxCount || p.fluxRadius >= centerDistance + n.fluxRadius),
                    "Parent expansion sphere does not contain a complete child expansion sphere.");
            require(n.halfWidth == p.halfWidth * .5, "Child is not a dyadic half cell.");
            for (unsigned axis = 0; axis < 3; ++axis)
                require(n.coordinates[axis] / 2 == p.coordinates[axis],
                        "Integer dyadic coordinate is inconsistent.");
        }
        if (!n.leaf()) {
            auto pressureEnd = n.firstPressure, fluxEnd = n.firstFlux;
            for (auto child : n.children)
                if (child != noSupportNode) {
                    const auto &c = tree.nodes[child];
                    require(c.parent == ni && c.firstPressure == pressureEnd &&
                                c.firstFlux == fluxEnd,
                            "Child spans do not partition their parent in octant order.");
                    pressureEnd += c.pressureCount;
                    fluxEnd += c.fluxCount;
                }
            require(pressureEnd == n.firstPressure + n.pressureCount &&
                        fluxEnd == n.firstFlux + n.fluxCount,
                    "Child spans lose parent supports.");
        }
        for (unsigned kind = 0; kind < 2; ++kind) {
            const auto &order = kind ? tree.fluxOrder : tree.pressureOrder;
            const auto &supports = kind ? tree.fluxSupports : tree.pressureSupports;
            const auto &bounds = kind ? n.fluxBounds : n.pressureBounds;
            const auto rad = kind ? n.fluxRadius : n.pressureRadius;
            const auto first = kind ? n.firstFlux : n.firstPressure,
                       count = kind ? n.fluxCount : n.pressureCount;
            require(std::size_t(first) + count <= order.size(), "Invalid node support span.");
            for (std::size_t j = first; j < std::size_t(first) + count; ++j) {
                const auto id = order[j];
                require(id < supports.size(), "Unknown DOF in support span.");
                const auto &s = supports[id];
                require(contains(bounds, s.bounds.lower) && contains(bounds, s.bounds.upper),
                        "Incomplete support bounds.");
                // The physical support is the union of complete triangles;
                // corners of its AABB generally do not belong to that union.
                const auto checkPanel = [&](std::uint32_t panel) {
                    for (const auto v : space.panels[panel].vertices)
                        require(std::hypot(v.x - n.center.x, v.y - n.center.y, v.z - n.center.z) <=
                                    rad + 1e-12,
                                "Support radius loses an actual incident panel vertex.");
                };
                if (kind)
                    checkPanel(id);
                else
                    for (auto incidence = space.incidenceOffsets[id];
                         incidence < space.incidenceOffsets[std::size_t(id) + 1]; ++incidence)
                        checkPanel(space.incidences[incidence].panel);
                if (n.leaf()) {
                    auto &members = kind ? qm : pm;
                    ++members[id];
                    require((kind ? tree.fluxLeaf[id] : tree.pressureLeaf[id]) == ni,
                            "Support leaf mapping inconsistent.");
                    crosses = crosses ||
                              std::abs(s.bounds.lower.x - n.center.x) > n.halfWidth + 1e-12 ||
                              std::abs(s.bounds.upper.x - n.center.x) > n.halfWidth + 1e-12;
                }
            }
        }
    }
    require(crosses || tree.nodes.size() == 1,
            "Fixture did not exercise support crossing dyadic cell boundaries.");
    for (auto n : pm)
        require(n == 1, "Pressure support lost or duplicated.");
    for (auto n : qm)
        require(n == 1, "Flux support lost or duplicated.");
    for (std::size_t v = 0; v < space.vertices.size(); ++v)
        for (auto i = space.incidenceOffsets[v]; i < space.incidenceOffsets[v + 1]; ++i)
            for (auto p : space.panels[space.incidences[i].panel].vertices)
                require(contains(tree.pressureSupports[v].bounds, p),
                        "P1 support omits part of its vertex star.");
}
void activation(const SupportTree &tree) {
    const auto descends = [&](std::uint32_t child, std::uint32_t ancestor) {
        while (child != noSupportNode) {
            if (child == ancestor)
                return true;
            child = tree.nodes[child].parent;
        }
        return false;
    };
    std::vector<SupportLevel> expected(tree.levels.size());
    using Offset = std::array<std::int64_t, 3>;
    std::vector<std::set<Offset>> pressureOffsets(tree.levels.size()),
        fluxOffsets(tree.levels.size());
    const auto recordOffsets = [&](const auto &partition, auto &offsets) {
        for (const auto pair : partition.farInteractions) {
            const auto &t = tree.nodes[pair.target], &s = tree.nodes[pair.source];
            offsets[t.level].insert({std::int64_t(t.coordinates[0]) - s.coordinates[0],
                                     std::int64_t(t.coordinates[1]) - s.coordinates[1],
                                     std::int64_t(t.coordinates[2]) - s.coordinates[2]});
        }
    };
    recordOffsets(tree.pressurePressure, pressureOffsets);
    recordOffsets(tree.pressureFlux, fluxOffsets);
    for (std::size_t i = 0; i < tree.nodes.size(); ++i) {
        const auto &node = tree.nodes[i];
        unsigned flags = 0;
        for (const auto pair : tree.pressurePressure.farInteractions) {
            if (node.pressureCount && descends(static_cast<std::uint32_t>(i), pair.target))
                flags |= pressureIncoming;
            if (node.pressureCount && descends(static_cast<std::uint32_t>(i), pair.source))
                flags |= pressureOutgoing;
        }
        for (const auto pair : tree.pressureFlux.farInteractions) {
            if (node.pressureCount && descends(static_cast<std::uint32_t>(i), pair.target))
                flags |= pressureIncoming;
            if (node.fluxCount && descends(static_cast<std::uint32_t>(i), pair.source))
                flags |= fluxOutgoing;
        }
        require(flags == node.activeChannels,
                "Expansion activation omits required descendants or includes unused ancestors.");
        auto &e = expected[node.level];
        if (flags & pressureOutgoing)
            ++e.pressureOutgoingNodes;
        if (flags & fluxOutgoing)
            ++e.fluxOutgoingNodes;
        if (flags & pressureIncoming)
            ++e.pressureIncomingNodes;
        if (flags & (pressureOutgoing | pressureIncoming))
            e.maximumActiveSupportRadius =
                std::max(e.maximumActiveSupportRadius, node.pressureRadius);
        if (flags & fluxOutgoing)
            e.maximumActiveSupportRadius = std::max(e.maximumActiveSupportRadius, node.fluxRadius);
    }
    for (auto pair : tree.pressurePressure.farInteractions)
        ++expected[tree.nodes[pair.target].level].pressurePressureFarInteractions;
    for (auto pair : tree.pressureFlux.farInteractions)
        ++expected[tree.nodes[pair.target].level].pressureFluxFarInteractions;
    for (std::size_t l = 0; l < tree.levels.size(); ++l) {
        const auto &e = expected[l], &a = tree.levels[l];
        require(e.pressureOutgoingNodes == a.pressureOutgoingNodes &&
                    e.fluxOutgoingNodes == a.fluxOutgoingNodes &&
                    e.pressureIncomingNodes == a.pressureIncomingNodes &&
                    e.maximumActiveSupportRadius == a.maximumActiveSupportRadius &&
                    e.pressurePressureFarInteractions == a.pressurePressureFarInteractions &&
                    e.pressureFluxFarInteractions == a.pressureFluxFarInteractions,
                "Active-level memory statistics disagree with independent ancestry oracle.");
        require(a.pressurePressureOffsets == pressureOffsets[l].size() &&
                    a.pressureFluxOffsets == fluxOffsets[l].size(),
                "Per-operator translation-class statistics disagree with actual dyadic offsets.");
        const auto distances = [](const auto &offsets) {
            std::set<std::uint64_t> values;
            for (const auto &v : offsets)
                values.insert(static_cast<std::uint64_t>(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
            return values;
        };
        auto pressureDistances = distances(pressureOffsets[l]);
        const auto fluxDistances = distances(fluxOffsets[l]);
        require(a.pressurePressureDistanceClasses == pressureDistances.size() &&
                    a.pressureFluxDistanceClasses == fluxDistances.size(),
                "Axial distance classes disagree with actual integer offsets.");
        pressureDistances.insert(fluxDistances.begin(), fluxDistances.end());
        require(a.squaredOffsetClasses ==
                    std::vector<std::uint64_t>(pressureDistances.begin(), pressureDistances.end()),
                "Axial distance-class list is incomplete or uses an inexact floating-point key.");
        pressureOffsets[l].insert(fluxOffsets[l].begin(), fluxOffsets[l].end());
        require(a.unionOffsets == pressureOffsets[l].size(),
                "Shared translation-class cache statistics lose or duplicate offsets.");
        const auto polarClasses = [](const auto &offsets) {
            std::set<std::array<std::int64_t, 3>> keys;
            for (const auto &offset : offsets) {
                const auto z2 = offset[2] * offset[2];
                const auto r2 = offset[0] * offset[0] + offset[1] * offset[1] + z2;
                const auto g = std::gcd(z2, r2);
                keys.insert({(offset[2] > 0) - (offset[2] < 0), z2 / g, r2 / g});
            }
            return keys;
        };
        require(a.pressureFluxPolarClasses == polarClasses(fluxOffsets[l]).size() &&
                    a.unionPolarClasses == polarClasses(pressureOffsets[l]).size() &&
                    a.pressurePressurePolarClasses <= a.unionPolarClasses,
                "Polar rotation classes disagree with exact rational cos(beta) keys.");
    }
    require(tree.nodes.front().activeChannels == 0 &&
                tree.levels.front().maximumActiveSupportRadius == 0,
            "Unused root must not allocate high-frequency spectra.");
}
void partition(const galerkin::Space &space, const SupportTree &tree, bool flux, double nearRatio,
               double separation) {
    const auto &part = flux ? tree.pressureFlux : tree.pressurePressure;
    const auto &order = flux ? tree.fluxOrder : tree.pressureOrder;
    const auto columns = order.size();
    std::vector<unsigned char> seen(space.vertices.size() * columns);
    std::uint64_t nearCount = 0, farCount = 0;
    auto visit = [&](SupportInteraction pair, bool far) {
        require(pair.target < tree.nodes.size() && pair.source < tree.nodes.size(),
                "Unknown interaction node.");
        const auto &a = tree.nodes[pair.target], &b = tree.nodes[pair.source];
        require(far ? a.level == b.level : (a.leaf() && b.leaf()),
                "Incorrect interaction level/leaf contract.");
        if (far) {
            const auto br = flux ? b.fluxRadius : b.pressureRadius;
            require(std::hypot(a.center.x - b.center.x, a.center.y - b.center.y,
                               a.center.z - b.center.z) > separation * (a.pressureRadius + br),
                    "Far pair violates complete-support separation.");
        }
        const auto first = flux ? b.firstFlux : b.firstPressure,
                   count = flux ? b.fluxCount : b.pressureCount;
        for (std::size_t i = a.firstPressure; i < std::size_t(a.firstPressure) + a.pressureCount;
             ++i)
            for (std::size_t j = first; j < std::size_t(first) + count; ++j) {
                const auto ti = tree.pressureOrder[i], si = order[j];
                require(++seen[std::size_t(ti) * columns + si] == 1, "Support pair is duplicated.");
                if (far) {
                    ++farCount;
                    for (auto x = space.incidenceOffsets[ti]; x < space.incidenceOffsets[ti + 1];
                         ++x) {
                        const auto tp = space.incidences[x].panel;
                        const auto checkPanel = [&](std::uint32_t sp) {
                            require(galerkin::classifyPair(space.panels[tp], space.panels[sp])
                                                .topology == galerkin::PairTopology::disjoint &&
                                        !galerkin::isNearPair(space.panels[tp], space.panels[sp],
                                                              nearRatio),
                                    "Far support pair contains singular, touching or near panels.");
                        };
                        if (flux)
                            checkPanel(si);
                        else
                            for (auto y = space.incidenceOffsets[si];
                                 y < space.incidenceOffsets[si + 1]; ++y)
                                checkPanel(space.incidences[y].panel);
                    }
                } else
                    ++nearCount;
            }
    };
    for (auto p : part.nearInteractions)
        visit(p, false);
    for (auto p : part.farInteractions)
        visit(p, true);
    std::uint64_t statisticPairs = 0;
    std::size_t statisticInteractions = 0;
    for (const auto &entry : part.farLevelPairs) {
        require(entry.targetLevel == entry.sourceLevel,
                "Unequal-level far statistics violate the angular-grid contract.");
        std::size_t count = 0;
        std::uint64_t pairs = 0;
        for (const auto interaction : part.farInteractions) {
            const auto &t = tree.nodes[interaction.target], &s = tree.nodes[interaction.source];
            if (t.level == entry.targetLevel && s.level == entry.sourceLevel) {
                ++count;
                pairs += std::uint64_t(t.pressureCount) * (flux ? s.fluxCount : s.pressureCount);
            }
        }
        require(count == entry.interactions && pairs == entry.representedPairs,
                "Far-level statistics disagree with pair enumeration.");
        statisticPairs += pairs;
        statisticInteractions += count;
    }
    require(statisticPairs == farCount && statisticInteractions == part.farInteractions.size(),
            "Far-level inventory loses interactions.");
    for (auto count : seen)
        require(count == 1, "Support partition omits a DOF pair.");
    require(nearCount == part.representedNearPairs && farCount == part.representedFarPairs &&
                nearCount + farCount == seen.size(),
            "Support pair statistics disagree with exact enumeration.");
}
std::array<std::uint64_t, 2> boxBoundedNearPairs(const SupportTree &tree) {
    // Independent conservative comparator: enclose the full support boxes,
    // including their artificial corners, then nest those balls bottom-up.
    // It uses the same tree and physical gap contract, so the comparison
    // isolates the effect of the physical triangle radius.
    std::vector<std::array<double, 2>> radii(tree.nodes.size());
    for (std::size_t i = tree.nodes.size(); i-- > 0;) {
        const auto &node = tree.nodes[i];
        const auto ball = [&](const SupportBounds &b) {
            return std::hypot(
                std::max(std::abs(b.lower.x - node.center.x), std::abs(b.upper.x - node.center.x)),
                std::max(std::abs(b.lower.y - node.center.y), std::abs(b.upper.y - node.center.y)),
                std::max(std::abs(b.lower.z - node.center.z), std::abs(b.upper.z - node.center.z)));
        };
        if (node.pressureCount)
            radii[i][0] = std::max(radii[i][0], ball(node.pressureBounds));
        if (node.fluxCount)
            radii[i][1] = std::max(radii[i][1], ball(node.fluxBounds));
        if (node.parent != noSupportNode) {
            const auto &parent = tree.nodes[node.parent];
            const auto shift =
                std::hypot(node.center.x - parent.center.x, node.center.y - parent.center.y,
                           node.center.z - parent.center.z);
            for (unsigned kind = 0; kind < 2; ++kind)
                if (kind ? node.fluxCount : node.pressureCount)
                    radii[node.parent][kind] =
                        std::max(radii[node.parent][kind],
                                 std::nextafter(shift + radii[i][kind],
                                                std::numeric_limits<double>::infinity()));
        }
        const auto rounding = 1 + 8 * std::numeric_limits<double>::epsilon();
        require(node.pressureRadius <= radii[i][0] * rounding &&
                    node.fluxRadius <= radii[i][1] * rounding,
                "Physical support balls exceed the conservative box-ball bound.");
    }
    const auto nearCount = [&](auto &&self, std::uint32_t ti, std::uint32_t si,
                               unsigned kind) -> std::uint64_t {
        const auto &t = tree.nodes[ti], &s = tree.nodes[si];
        const auto ns = kind ? s.fluxCount : s.pressureCount;
        if (!t.pressureCount || !ns)
            return 0;
        const auto &b = kind ? s.fluxBounds : s.pressureBounds;
        const auto &a = t.pressureBounds;
        const auto gap = std::hypot(std::max({0., a.lower.x - b.upper.x, b.lower.x - a.upper.x}),
                                    std::max({0., a.lower.y - b.upper.y, b.lower.y - a.upper.y}),
                                    std::max({0., a.lower.z - b.upper.z, b.lower.z - a.upper.z}));
        const auto distance =
            std::hypot(t.center.x - s.center.x, t.center.y - s.center.y, t.center.z - s.center.z);
        if (t.level == s.level && gap > 0 &&
            gap > tree.nearDistanceRatio * std::max(a.maximumPanelEdge, b.maximumPanelEdge) &&
            distance > tree.separationRatio * (radii[ti][0] + radii[si][kind]))
            return 0;
        if (t.leaf() && s.leaf())
            return std::uint64_t(t.pressureCount) * ns;
        std::uint64_t count = 0;
        if (!t.leaf() && (s.leaf() || t.level < s.level)) {
            for (const auto child : t.children)
                if (child != noSupportNode)
                    count += self(self, child, si, kind);
        } else if (!s.leaf() && (t.leaf() || s.level < t.level)) {
            for (const auto child : s.children)
                if (child != noSupportNode)
                    count += self(self, ti, child, kind);
        } else
            for (const auto tc : t.children)
                if (tc != noSupportNode)
                    for (const auto sc : s.children)
                        if (sc != noSupportNode)
                            count += self(self, tc, sc, kind);
        return count;
    };
    return {nearCount(nearCount, 0, 0, 0), nearCount(nearCount, 0, 0, 1)};
}
void run() {
    auto mesh = separated();
    auto space = galerkin::buildSpace(mesh);
    SupportTreeSettings settings;
    settings.leafCapacity = 4;
    auto tree = makeSupportTree(space, settings);
    hierarchy(space, tree);
    activation(tree);
    partition(space, tree, false, settings.nearDistanceRatio, settings.separationRatio);
    partition(space, tree, true, settings.nearDistanceRatio, settings.separationRatio);
    require(!tree.pressurePressure.farInteractions.empty() &&
                !tree.pressureFlux.farInteractions.empty(),
            "Separated bodies failed to exercise far partitions.");
    require(tree.allocatedBytes() > tree.nodes.size() * sizeof(SupportNode),
            "Support RAM accounting omits owned arrays.");
    auto gradedMesh = testing::sphere(3);
    const auto coarseMesh = testing::sphere(0);
    const auto base = static_cast<std::uint32_t>(gradedMesh.vertices.size());
    for (auto p : coarseMesh.vertices)
        gradedMesh.vertices.push_back({p.x + 2, p.y, p.z});
    for (auto t : coarseMesh.triangles) {
        for (auto &v : t.vertices)
            v += base;
        gradedMesh.triangles.push_back(t);
    }
    const auto gradedSpace = galerkin::buildSpace(gradedMesh);
    auto gradedSettings = settings;
    gradedSettings.leafCapacity = 8;
    const auto graded = makeSupportTree(gradedSpace, gradedSettings);
    hierarchy(gradedSpace, graded);
    activation(graded);
    partition(gradedSpace, graded, false, gradedSettings.nearDistanceRatio,
              gradedSettings.separationRatio);
    partition(gradedSpace, graded, true, gradedSettings.nearDistanceRatio,
              gradedSettings.separationRatio);
    const auto boxNear = boxBoundedNearPairs(graded);
    std::cout << "Graded physical/box-ball direct pairs: P1/P1 "
              << graded.pressurePressure.representedNearPairs << "/" << boxNear[0] << "; P1/DP0 "
              << graded.pressureFlux.representedNearPairs << "/" << boxNear[1] << ".\n";
    require(graded.pressurePressure.representedNearPairs < boxNear[0] &&
                graded.pressureFlux.representedNearPairs < boxNear[1],
            "Actual triangle balls did not reduce direct work on the graded curved fixture.");
    require(std::any_of(graded.pressureFlux.nearInteractions.begin(),
                        graded.pressureFlux.nearInteractions.end(),
                        [&](auto pair) {
                            return graded.nodes[pair.target].level !=
                                   graded.nodes[pair.source].level;
                        }),
            "Graded fixture must exercise exact near products between unequal levels.");
    settings.maximumDepth = 0;
    auto shallow = makeSupportTree(space, settings);
    hierarchy(space, shallow);
    activation(shallow);
    partition(space, shallow, false, settings.nearDistanceRatio, settings.separationRatio);
    partition(space, shallow, true, settings.nearDistanceRatio, settings.separationRatio);
    require(shallow.nodes.size() == 1 && shallow.pressureFlux.farInteractions.empty(),
            "Depth bound did not retain exact near work.");
    double exactRadius = 0;
    const auto center = shallow.nodes.front().center;
    for (const auto &panel : space.panels)
        for (const auto p : panel.vertices)
            exactRadius =
                std::max(exactRadius, std::hypot(p.x - center.x, p.y - center.y, p.z - center.z));
    require(shallow.nodes.front().pressureRadius <= exactRadius * (1 + 1e-14) &&
                shallow.nodes.front().fluxRadius <= exactRadius * (1 + 1e-14),
            "Leaf ball encloses artificial support-box corners instead of physical triangles.");
    const auto &rootBounds = shallow.nodes.front().pressureBounds;
    const auto boxRadius = std::hypot(rootBounds.upper.x - center.x, rootBounds.upper.y - center.y,
                                      rootBounds.upper.z - center.z);
    require(exactRadius < boxRadius,
            "Shallow fixture does not distinguish physical vertices from artificial box corners.");
    for (auto &p : mesh.vertices) {
        p.x += 1e6;
        p.y -= 2e6;
        p.z += 3e6;
    }
    const auto movedSpace = galerkin::buildSpace(mesh);
    settings.maximumDepth = 20;
    const auto moved = makeSupportTree(movedSpace, settings);
    hierarchy(movedSpace, moved);
    activation(moved);
    partition(movedSpace, moved, false, settings.nearDistanceRatio, settings.separationRatio);
    partition(movedSpace, moved, true, settings.nearDistanceRatio, settings.separationRatio);
    settings.maximumDepth = 20;
    settings.maximumInteractions = 1;
    rejects([&] { makeSupportTree(space, settings); });
    settings.maximumInteractions = std::numeric_limits<std::size_t>::max();
    std::stop_source stop;
    stop.request_stop();
    rejects([&] { makeSupportTree(space, settings, stop.get_token()); });
    settings.leafCapacity = 0;
    rejects([&] { makeSupportTree(space, settings); });
    settings.leafCapacity = 4;
    settings.separationRatio = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { makeSupportTree(space, settings); });
    std::cout << "Support partition: " << tree.nodes.size() << " nodes, " << tree.allocatedBytes()
              << " allocated bytes; P1/P1 " << tree.pressurePressure.representedFarPairs
              << " far DOF pairs; P1/DP0 " << tree.pressureFlux.representedFarPairs << ".\n";
}
} // namespace
int main() {
    try {
        run();
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
