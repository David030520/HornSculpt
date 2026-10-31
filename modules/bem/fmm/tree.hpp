#pragma once

#include "../galerkin/space.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stop_token>
#include <vector>

namespace hornsculpt::fmm {

inline constexpr std::uint32_t noNode = UINT32_MAX;

struct Node {
    Vec3 center{};
    double halfWidth{}, maximumPanelEdge{};
    std::uint32_t parent{noNode};
    std::array<std::uint32_t, 2> children{noNode, noNode};
    std::uint32_t firstPanel{}, panelCount{}, depth{};
    bool leaf() const noexcept {
        return children[0] == noNode;
    }
};

struct Interaction {
    std::uint32_t test{}, trial{};
};

struct Tree {
    std::vector<Node> nodes;
    std::vector<std::uint32_t> panelOrder;
    std::vector<std::uint32_t> panelLeaf;
    std::vector<std::uint32_t> leafNodes;
    std::vector<std::vector<std::uint32_t>> levels;
    // The Cartesian products of these cluster pairs and direct panel pairs
    // partition every ordered panel pair exactly once, including the diagonal.
    std::vector<Interaction> farInteractions;
    std::vector<Interaction> directInteractions;
    std::size_t representedFarPanelPairs{};
    std::size_t auditJobCount{}, constructionMinimumBytes{};
};

// Caller-supplied mandatory storage costs. Counting never allocates pair lists.
struct TreeMemoryBudget {
    std::size_t maximumDirectPairs{std::numeric_limits<std::size_t>::max()};
    std::size_t maximumBytes{std::numeric_limits<std::size_t>::max()};
    std::size_t fixedBytes{}, nodeBytes{}, directPairBytes{}, farPairBytes{}, auditJobBytes{};
};

// A balanced spatial hierarchy encloses complete panels in interpolation cubes.
// Disjoint cubes alone are insufficient: both the separation and the electrical
// box size must pass before an interaction is admitted to M2L. All panel pairs
// inside the configured near-distance contract remain direct.
// Leaf capacity is a maximum: electrically unresolved leaves split further,
// down to individual complete panels, without changing the mesh or quadrature.
// The optional budget counts mandatory node, direct, far and expanded audit
// storage before materializing interaction lists; it never changes the mesh.
Tree makeTree(const galerkin::Space &space, unsigned leafPanelCapacity, unsigned interpolationOrder,
              double separationRatio, double nearDistanceRatio, double wavenumber,
              std::stop_token stop = {}, TreeMemoryBudget budget = {});

// First-kind Chebyshev interpolation; values are L_i(x), x in [-1,1].
std::vector<double> chebyshevNodes(unsigned order);
std::vector<double> interpolationWeights(const std::vector<double> &nodes, double x);

} // namespace hornsculpt::fmm
