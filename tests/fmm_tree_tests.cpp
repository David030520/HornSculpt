#include "fixture_geometry.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fmm/tree.hpp>
#include <galerkin/pair_quadrature.hpp>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <stop_token>
#include <string>
#include <vector>

using namespace hornsculpt;

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

void checkInterpolation() {
    for (unsigned order = 3; order <= 8; ++order) {
        const auto nodes = fmm::chebyshevNodes(order);
        require(nodes.size() == order, "Chebyshev node count is incorrect.");
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            require(std::isfinite(nodes[i]) && std::abs(nodes[i]) < 1,
                    "First-kind Chebyshev node lies outside its interval.");
            const auto values = fmm::interpolationWeights(nodes, nodes[i]);
            for (std::size_t j = 0; j < values.size(); ++j)
                require(std::abs(values[j] - (i == j ? 1.0 : 0.0)) < 1e-13,
                        "Chebyshev interpolation fails at its own nodes.");
        }
        for (const double x : std::array{-1.0, -.8321, -.1, 0.0, .347, .999, 1.0}) {
            const auto values = fmm::interpolationWeights(nodes, x);
            for (unsigned degree = 0; degree < order; ++degree) {
                double interpolated = 0;
                for (std::size_t j = 0; j < values.size(); ++j)
                    interpolated += values[j] * std::pow(nodes[j], static_cast<int>(degree));
                require(std::abs(interpolated - std::pow(x, static_cast<int>(degree))) < 2e-12,
                        "Chebyshev interpolation fails polynomial reproduction.");
            }
        }
    }
}

void checkHierarchy(const galerkin::Space &space, const fmm::Tree &tree, unsigned capacity) {
    require(!tree.nodes.empty() && tree.nodes.front().parent == fmm::noNode &&
                tree.panelOrder.size() == space.panels.size() &&
                tree.panelLeaf.size() == space.panels.size(),
            "FMM tree has inconsistent root or panel mapping.");
    std::vector<unsigned> membership(space.panels.size());
    std::vector<unsigned> levels(tree.nodes.size());
    for (std::size_t depth = 0; depth < tree.levels.size(); ++depth)
        for (const auto node : tree.levels[depth]) {
            require(node < tree.nodes.size() && tree.nodes[node].depth == depth,
                    "FMM level contains a node from a different depth.");
            ++levels[node];
        }
    for (std::size_t index = 0; index < tree.nodes.size(); ++index) {
        const auto &node = tree.nodes[index];
        require(levels[index] == 1 && node.panelCount > 0 &&
                    static_cast<std::size_t>(node.firstPanel) + node.panelCount <=
                        tree.panelOrder.size() &&
                    std::isfinite(node.halfWidth) && node.halfWidth > 0,
                "FMM node has invalid level, bounds or panel interval.");
        for (std::size_t p = node.firstPanel; p < node.firstPanel + node.panelCount; ++p) {
            const auto panel = tree.panelOrder[p];
            require(panel < space.panels.size(), "FMM tree refers to an unknown panel.");
            for (const auto vertex : space.panels[panel].vertices)
                require(std::abs(vertex.x - node.center.x) <= node.halfWidth + 1e-13 &&
                            std::abs(vertex.y - node.center.y) <= node.halfWidth + 1e-13 &&
                            std::abs(vertex.z - node.center.z) <= node.halfWidth + 1e-13,
                        "FMM interpolation cube fails to enclose a complete panel.");
        }
        if (node.leaf()) {
            require(node.panelCount <= capacity, "FMM leaf exceeds its panel capacity.");
            for (std::size_t p = node.firstPanel; p < node.firstPanel + node.panelCount; ++p) {
                const auto panel = tree.panelOrder[p];
                ++membership[panel];
                require(tree.panelLeaf[panel] == index, "FMM panel has an inconsistent owner.");
            }
        } else {
            require(node.children[0] < tree.nodes.size() && node.children[1] < tree.nodes.size(),
                    "FMM parent has an invalid child index.");
            const auto &left = tree.nodes[node.children[0]], &right = tree.nodes[node.children[1]];
            require(left.parent == index && right.parent == index && left.depth == node.depth + 1 &&
                        right.depth == node.depth + 1 && left.firstPanel == node.firstPanel &&
                        right.firstPanel == left.firstPanel + left.panelCount &&
                        left.panelCount + right.panelCount == node.panelCount,
                    "FMM children do not partition their parent's complete panel interval.");
        }
    }
    require(
        std::all_of(membership.begin(), membership.end(), [](auto count) { return count == 1; }),
        "FMM leaves lose or duplicate panels.");
    require(tree.leafNodes.size() == static_cast<std::size_t>(std::count_if(
                                         tree.nodes.begin(), tree.nodes.end(),
                                         [](const auto &node) { return node.leaf(); })),
            "FMM leaf inventory is inconsistent.");
}

void checkPartition(const galerkin::Space &space, const fmm::Tree &tree, double nearRatio, double k,
                    unsigned order, double separation) {
    const auto count = space.panels.size();
    std::vector<unsigned char> membership(count * count);
    std::size_t represented = 0;
    const auto add = [&](std::uint32_t test, std::uint32_t trial, bool far) {
        require(test < count && trial < count, "FMM interaction uses an unknown panel.");
        auto &visited = membership[static_cast<std::size_t>(test) * count + trial];
        require(visited == 0, "FMM interaction partition duplicates a panel pair.");
        visited = far ? 2 : 1;
        if (far) {
            require(galerkin::classifyPair(space.panels[test], space.panels[trial]).topology ==
                            galerkin::PairTopology::disjoint &&
                        !galerkin::isNearPair(space.panels[test], space.panels[trial], nearRatio),
                    "A singular or near panel pair was incorrectly assigned to M2L.");
            ++represented;
        }
    };
    for (const auto pair : tree.directInteractions)
        add(pair.test, pair.trial, false);
    for (const auto pair : tree.farInteractions) {
        require(pair.test < tree.nodes.size() && pair.trial < tree.nodes.size(),
                "M2L interaction uses an unknown node.");
        const auto &test = tree.nodes[pair.test], &trial = tree.nodes[pair.trial];
        const double radiusTest = std::sqrt(3.0) * test.halfWidth,
                     radiusTrial = std::sqrt(3.0) * trial.halfWidth;
        const double distance =
            std::hypot(test.center.x - trial.center.x, test.center.y - trial.center.y,
                       test.center.z - trial.center.z);
        require(distance > separation * (radiusTest + radiusTrial) &&
                    k * std::max(radiusTest, radiusTrial) <= static_cast<double>(order) / 2,
                "M2L interaction violates spatial or electrical-size admissibility.");
        for (std::size_t t = test.firstPanel; t < test.firstPanel + test.panelCount; ++t)
            for (std::size_t s = trial.firstPanel; s < trial.firstPanel + trial.panelCount; ++s)
                add(tree.panelOrder[t], tree.panelOrder[s], true);
    }
    require(
        std::all_of(membership.begin(), membership.end(), [](auto count) { return count != 0; }),
        "FMM interaction partition omits a panel pair.");
    require(represented == tree.representedFarPanelPairs &&
                represented + tree.directInteractions.size() == count * count,
            "FMM represented-pair statistics do not describe the actual partition.");
    for (std::size_t t = 0; t < count; ++t)
        require(membership[t * count + t] == 1, "A self interaction was assigned to M2L.");
}

void checkTree() {
    const auto space = galerkin::buildSpace(testing::sphere(3));
    constexpr unsigned capacity = 4, order = 5;
    constexpr double separation = 2.0, nearRatio = 1.0;
    const auto low = fmm::makeTree(space, capacity, order, separation, nearRatio, 10);
    checkHierarchy(space, low, capacity);
    checkPartition(space, low, nearRatio, 10, order, separation);
    require(!low.farInteractions.empty(), "Tree test did not exercise M2L interactions.");
    const auto high = fmm::makeTree(space, capacity, order, separation, nearRatio, 500);
    checkHierarchy(space, high, capacity);
    checkPartition(space, high, nearRatio, 500, order, separation);
    for (const auto leaf : high.leafNodes)
        require(high.nodes[leaf].panelCount == 1 ||
                    500 * std::sqrt(3.0) * high.nodes[leaf].halfWidth <= .5 * order,
                "An electrically unresolved cluster stopped splitting before a single panel.");
    require(high.representedFarPanelPairs < low.representedFarPanelPairs,
            "Helmholtz tree ignored the wavelength when admitting far interactions.");
    const auto fineSpace = galerkin::buildSpace(testing::sphere(4));
    const auto fine = fmm::makeTree(fineSpace, capacity, order, separation, nearRatio, 10);
    checkHierarchy(fineSpace, fine, capacity);
    const double lowFraction = static_cast<double>(low.directInteractions.size()) /
                               static_cast<double>(space.panels.size() * space.panels.size());
    const double fineFraction =
        static_cast<double>(fine.directInteractions.size()) /
        static_cast<double>(fineSpace.panels.size() * fineSpace.panels.size());
    require(fineFraction < lowFraction,
            "Refining the mesh did not reduce the direct pair fraction at fixed electrical size.");
    require(std::any_of(fine.farInteractions.begin(), fine.farInteractions.end(),
                        [&](const auto pair) {
                            return !fine.nodes[pair.test].leaf() || !fine.nodes[pair.trial].leaf();
                        }),
            "Fine FMM tree has no interaction above the leaf level.");
    std::cout << "Direct pair fraction 512T/2048T: " << lowFraction << '/' << fineFraction
              << "; far blocks " << low.farInteractions.size() << '/' << fine.farInteractions.size()
              << '\n';
    std::stop_source cancelled;
    cancelled.request_stop();
    bool rejected = false;
    try {
        (void)fmm::makeTree(space, capacity, order, separation, nearRatio, 10,
                            cancelled.get_token());
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    require(rejected, "FMM tree construction ignored cancellation.");
    rejected = false;
    try {
        (void)fmm::makeTree(space, capacity, order, separation, nearRatio, 10, {},
                            {.maximumDirectPairs = 7});
    } catch (const std::runtime_error &error) {
        rejected = std::string(error.what()).find("memory budget") != std::string::npos;
    }
    require(rejected, "FMM tree construction exceeded its explicit direct-pair memory budget.");
    std::size_t auditCount = 0;
    for (const auto pair : low.farInteractions)
        auditCount += low.nodes[pair.test].panelCount;
    require(low.auditJobCount == auditCount,
            "Counted audit jobs disagree with the exact far partition.");
    fmm::TreeMemoryBudget budget;
    budget.nodeBytes = 3;
    budget.directPairBytes = 5;
    budget.farPairBytes = 7;
    budget.auditJobBytes = 11;
    budget.maximumBytes = 3 * low.nodes.size() + 5 * low.directInteractions.size() +
                          7 * low.farInteractions.size() + 11 * auditCount;
    const auto bounded =
        fmm::makeTree(space, capacity, order, separation, nearRatio, 10, {}, budget);
    const auto equal = [](const auto &a, const auto &b) {
        return a.test == b.test && a.trial == b.trial;
    };
    require(bounded.constructionMinimumBytes == budget.maximumBytes &&
                std::equal(low.directInteractions.begin(), low.directInteractions.end(),
                           bounded.directInteractions.begin(), bounded.directInteractions.end(),
                           equal) &&
                std::equal(low.farInteractions.begin(), low.farInteractions.end(),
                           bounded.farInteractions.begin(), bounded.farInteractions.end(), equal),
            "Count/reserve/materialize changed the ordered interaction partition.");
    for (unsigned component = 0; component < 4; ++component) {
        fmm::TreeMemoryBudget limit;
        if (component == 0) {
            limit.nodeBytes = 1;
            limit.maximumBytes = low.nodes.size() - 1;
        }
        if (component == 1) {
            limit.directPairBytes = 1;
            limit.maximumBytes = low.directInteractions.size() - 1;
        }
        if (component == 2) {
            limit.farPairBytes = 1;
            limit.maximumBytes = low.farInteractions.size() - 1;
        }
        if (component == 3) {
            limit.auditJobBytes = 1;
            limit.maximumBytes = auditCount - 1;
        }
        rejected = false;
        try {
            (void)fmm::makeTree(space, capacity, order, separation, nearRatio, 10, {}, limit);
        } catch (const std::runtime_error &error) {
            rejected = std::string(error.what()).find("before interaction materialization") !=
                       std::string::npos;
        }
        require(rejected, "A mandatory FMM storage component escaped the counting budget.");
    }
    budget.fixedBytes = std::numeric_limits<std::size_t>::max();
    rejected = false;
    try {
        (void)fmm::makeTree(space, capacity, order, separation, nearRatio, 10, {}, budget);
    } catch (const std::overflow_error &) {
        rejected = true;
    }
    require(rejected, "FMM mandatory storage accounting overflowed silently.");
}
} // namespace

int main() {
    try {
        checkInterpolation();
        checkTree();
        std::cout << "FMM interpolation, hierarchy, near protection, pair partition and scaling "
                     "passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
