#include "tree.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace hornsculpt::fmm {
namespace {
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Fast Multipole construction cancelled.");
}
double component(Vec3 value, unsigned axis) {
    return axis == 0 ? value.x : axis == 1 ? value.y : value.z;
}
double cubeGap(const Node &a, const Node &b) {
    const double h = a.halfWidth + b.halfWidth;
    const double x = std::max(0.0, std::abs(a.center.x - b.center.x) - h);
    const double y = std::max(0.0, std::abs(a.center.y - b.center.y) - h);
    const double z = std::max(0.0, std::abs(a.center.z - b.center.z) - h);
    return std::sqrt(x * x + y * y + z * z);
}
class Builder {
  public:
    Builder(const galerkin::Space &geometry, unsigned capacity, unsigned order, double separation,
            double nearRatio, double k, std::stop_token token, TreeMemoryBudget memoryBudget)
        : space(geometry), leafCapacity(capacity), polynomialOrder(order),
          separationRatio(separation), nearDistanceRatio(nearRatio), wavenumber(k), stop(token),
          budget(memoryBudget) {}

    Tree build() {
        tree.panelOrder.resize(space.panels.size());
        tree.panelLeaf.resize(space.panels.size());
        std::iota(tree.panelOrder.begin(), tree.panelOrder.end(), 0u);
        makeNode(0, static_cast<std::uint32_t>(space.panels.size()), noNode, 0);
        checkBudget();
        // First traversal counts exact storage, including the expanded audit
        // jobs. A rejected plan never materializes an enormous interaction list.
        visit(0, 0);
        tree.directInteractions.reserve(directCount);
        tree.farInteractions.reserve(farCount);
        materialize = true;
        visit(0, 0);
        if (tree.directInteractions.size() != directCount ||
            tree.farInteractions.size() != farCount)
            throw std::logic_error("FMM counted and materialized interaction partitions disagree.");
        const auto less = [](const Interaction &a, const Interaction &b) {
            return a.test == b.test ? a.trial < b.trial : a.test < b.test;
        };
        std::sort(tree.directInteractions.begin(), tree.directInteractions.end(), less);
        std::sort(tree.farInteractions.begin(), tree.farInteractions.end(), less);
        return std::move(tree);
    }

  private:
    void checkBudget() {
        std::size_t bytes = budget.fixedBytes;
        const auto add = [&bytes](std::size_t count, std::size_t width) {
            if (width && count > (std::numeric_limits<std::size_t>::max() - bytes) / width)
                throw std::overflow_error("Fast Multipole construction byte count overflow.");
            bytes += count * width;
        };
        add(tree.nodes.size(), budget.nodeBytes);
        add(directCount, budget.directPairBytes);
        add(farCount, budget.farPairBytes);
        add(tree.auditJobCount, budget.auditJobBytes);
        tree.constructionMinimumBytes = bytes;
        if (bytes > budget.maximumBytes || directCount > budget.maximumDirectPairs) {
            std::ostringstream message;
            message << "Fast Multipole counting exceeds its memory budget before interaction "
                       "materialization: nodes="
                    << tree.nodes.size() << ", direct=" << directCount << ", far=" << farCount
                    << ", audit=" << tree.auditJobCount << ", minimumBytes=" << bytes
                    << ", budgetBytes=" << budget.maximumBytes << ". Mesh was not changed.";
            throw std::runtime_error(message.str());
        }
    }
    std::uint32_t makeNode(std::uint32_t first, std::uint32_t count, std::uint32_t parent,
                           std::uint32_t depth) {
        cancelled(stop);
        if (tree.nodes.size() >= noNode)
            throw std::overflow_error("Fast Multipole node indices overflow.");
        const auto id = static_cast<std::uint32_t>(tree.nodes.size());
        Node node;
        node.parent = parent;
        node.firstPanel = first;
        node.panelCount = count;
        node.depth = depth;
        const double inf = std::numeric_limits<double>::infinity();
        Vec3 lower{inf, inf, inf}, upper{-inf, -inf, -inf};
        for (std::uint32_t i = first; i < first + count; ++i) {
            const auto &panel = space.panels[tree.panelOrder[i]];
            node.maximumPanelEdge = std::max(node.maximumPanelEdge, panel.maximumEdgeLength);
            for (const auto &v : panel.vertices) {
                lower.x = std::min(lower.x, v.x);
                lower.y = std::min(lower.y, v.y);
                lower.z = std::min(lower.z, v.z);
                upper.x = std::max(upper.x, v.x);
                upper.y = std::max(upper.y, v.y);
                upper.z = std::max(upper.z, v.z);
            }
        }
        node.center = {(lower.x + upper.x) * .5, (lower.y + upper.y) * .5,
                       (lower.z + upper.z) * .5};
        const double extents[3]{upper.x - lower.x, upper.y - lower.y, upper.z - lower.z};
        const auto axis = static_cast<unsigned>(std::max_element(extents, extents + 3) - extents);
        // Round outward, so all leaf quadrature points stay within [-1,1].
        node.halfWidth = std::nextafter(.5 * extents[axis], inf);
        if (!(node.halfWidth > 0) || !std::isfinite(node.halfWidth))
            throw std::invalid_argument("Fast Multipole panel bounds are invalid.");
        tree.nodes.push_back(node);
        if (tree.levels.size() <= depth)
            tree.levels.resize(static_cast<std::size_t>(depth) + 1);
        tree.levels[depth].push_back(id);
        // Capacity is an upper bound. Electrically large clusters must keep
        // subdividing, otherwise even far-away panels would become permanently
        // direct merely because a leaf was reached before the Helmholtz
        // interpolation resolution criterion. A single panel is indivisible;
        // its unresolved interactions correctly remain direct.
        const bool electricallyResolved =
            wavenumber * std::sqrt(3.0) * node.halfWidth <= .5 * polynomialOrder;
        if (count == 1 || (count <= leafCapacity && electricallyResolved)) {
            tree.leafNodes.push_back(id);
            for (std::uint32_t i = first; i < first + count; ++i)
                tree.panelLeaf[tree.panelOrder[i]] = id;
        } else {
            const auto middle = first + count / 2;
            std::nth_element(tree.panelOrder.begin() + first, tree.panelOrder.begin() + middle,
                             tree.panelOrder.begin() + first + count,
                             [&](std::uint32_t a, std::uint32_t b) {
                                 const double x = component(space.panels[a].centroid, axis),
                                              y = component(space.panels[b].centroid, axis);
                                 return x == y ? a < b : x < y;
                             });
            const auto left = makeNode(first, count / 2, id, depth + 1);
            const auto right = makeNode(middle, count - count / 2, id, depth + 1);
            tree.nodes[id].children = {left, right};
        }
        return id;
    }
    bool admissible(const Node &a, const Node &b) const {
        const double dx = a.center.x - b.center.x, dy = a.center.y - b.center.y,
                     dz = a.center.z - b.center.z;
        const double radiusA = std::sqrt(3.0) * a.halfWidth, radiusB = std::sqrt(3.0) * b.halfWidth;
        // Low-frequency interpolation FMM: refine electrically large clusters.
        // This is an explicit resolution criterion, not an a-posteriori error bound.
        if (wavenumber * std::max(radiusA, radiusB) > .5 * polynomialOrder)
            return false;
        const double minimumDistance = separationRatio * (radiusA + radiusB);
        return dx * dx + dy * dy + dz * dz > minimumDistance * minimumDistance &&
               cubeGap(a, b) > nearDistanceRatio * std::max(a.maximumPanelEdge, b.maximumPanelEdge);
    }
    void visit(std::uint32_t test, std::uint32_t trial) {
        cancelled(stop);
        const auto &a = tree.nodes[test];
        const auto &b = tree.nodes[trial];
        if (test != trial && admissible(a, b)) {
            if (materialize) {
                tree.farInteractions.push_back({test, trial});
                return;
            }
            if (farCount == std::numeric_limits<std::size_t>::max() ||
                a.panelCount > std::numeric_limits<std::size_t>::max() - tree.auditJobCount)
                throw std::overflow_error("Fast Multipole audit count overflow.");
            ++farCount;
            tree.auditJobCount += a.panelCount;
            const auto represented = static_cast<std::size_t>(a.panelCount) * b.panelCount;
            if (represented >
                std::numeric_limits<std::size_t>::max() - tree.representedFarPanelPairs)
                throw std::overflow_error("Fast Multipole interaction count overflow.");
            tree.representedFarPanelPairs += represented;
            checkBudget();
        } else if (a.leaf() && b.leaf()) {
            const auto additional = static_cast<std::size_t>(a.panelCount) * b.panelCount;
            if (!materialize) {
                if (additional > tree.directInteractions.max_size() - directCount)
                    throw std::overflow_error("Fast Multipole direct interaction count overflow.");
                directCount += additional;
                checkBudget();
                return;
            }
            for (std::uint32_t i = a.firstPanel; i < a.firstPanel + a.panelCount; ++i)
                for (std::uint32_t j = b.firstPanel; j < b.firstPanel + b.panelCount; ++j)
                    tree.directInteractions.push_back({tree.panelOrder[i], tree.panelOrder[j]});
        } else if (!a.leaf() && (b.leaf() || a.halfWidth >= b.halfWidth)) {
            visit(a.children[0], trial);
            visit(a.children[1], trial);
        } else {
            visit(test, b.children[0]);
            visit(test, b.children[1]);
        }
    }
    const galerkin::Space &space;
    unsigned leafCapacity, polynomialOrder;
    double separationRatio, nearDistanceRatio, wavenumber;
    std::stop_token stop;
    TreeMemoryBudget budget;
    std::size_t directCount{}, farCount{};
    bool materialize{};
    Tree tree;
};
} // namespace

Tree makeTree(const galerkin::Space &space, unsigned leafPanelCapacity, unsigned interpolationOrder,
              double separationRatio, double nearDistanceRatio, double wavenumber,
              std::stop_token stop, TreeMemoryBudget budget) {
    if (space.panels.empty() || space.panels.size() > UINT32_MAX || leafPanelCapacity < 4 ||
        leafPanelCapacity > 256 || interpolationOrder < 3 || interpolationOrder > 8 ||
        !std::isfinite(separationRatio) || separationRatio < 1.5 || separationRatio > 10 ||
        !std::isfinite(nearDistanceRatio) || nearDistanceRatio < 0 || !std::isfinite(wavenumber) ||
        wavenumber < 0)
        throw std::invalid_argument("Invalid Fast Multipole hierarchy settings.");
    return Builder(space, leafPanelCapacity, interpolationOrder, separationRatio, nearDistanceRatio,
                   wavenumber, stop, budget)
        .build();
}

std::vector<double> chebyshevNodes(unsigned order) {
    if (order < 2 || order > 32)
        throw std::invalid_argument("Invalid Chebyshev interpolation order.");
    std::vector<double> nodes(order);
    for (unsigned i = 0; i < order; ++i)
        nodes[i] = std::cos(std::numbers::pi * (2 * i + 1) / (2 * order));
    return nodes;
}
std::vector<double> interpolationWeights(const std::vector<double> &nodes, double x) {
    if (nodes.size() < 2 || !std::isfinite(x))
        throw std::invalid_argument("Invalid Chebyshev interpolation argument.");
    std::vector<double> weights(nodes.size());
    double total = 0;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (std::abs(x - nodes[i]) <= 8 * std::numeric_limits<double>::epsilon()) {
            weights[i] = 1;
            for (std::size_t j = 0; j < i; ++j)
                weights[j] = 0;
            return weights;
        }
        const double sine = std::sin(std::numbers::pi * (2 * i + 1) / (2 * nodes.size()));
        weights[i] = (i % 2 ? -sine : sine) / (x - nodes[i]);
        total += weights[i];
    }
    if (!std::isfinite(total) || total == 0)
        throw std::runtime_error("Chebyshev interpolation is not finite.");
    for (auto &weight : weights)
        weight /= total;
    return weights;
}
} // namespace hornsculpt::fmm
