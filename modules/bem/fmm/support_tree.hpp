#pragma once

#include "../galerkin/space.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stop_token>
#include <vector>

namespace hornsculpt::wideband {

inline constexpr std::uint32_t noSupportNode = UINT32_MAX;
enum class SupportKind { pressure, flux };
inline constexpr unsigned pressureOutgoing = 1, fluxOutgoing = 2, pressureIncoming = 4;

struct SupportBounds {
    Vec3 lower{}, upper{};
    double maximumPanelEdge{};
};

struct Support {
    Vec3 anchor{};
    SupportBounds bounds{};
};

struct SupportNode {
    Vec3 center{};
    double halfWidth{};
    std::array<std::uint32_t, 3> coordinates{};
    std::uint32_t level{}, parent{noSupportNode};
    std::array<std::uint32_t, 8> children{noSupportNode, noSupportNode, noSupportNode,
                                          noSupportNode, noSupportNode, noSupportNode,
                                          noSupportNode, noSupportNode};
    std::uint32_t firstPressure{}, pressureCount{}, firstFlux{}, fluxCount{};
    SupportBounds pressureBounds{}, fluxBounds{};
    // Radius about the dyadic cell center encloses every complete support
    // triangle and each nonempty child expansion ball (distance to child center
    // + child radius). Leaf balls use actual triangle vertices, not artificial
    // AABB corners. Panels may cross cell boundaries; this is not a cell radius.
    double pressureRadius{}, fluxRadius{};
    unsigned activeChannels{};
    bool leaf() const noexcept;
};

struct SupportInteraction {
    std::uint32_t target{}, source{};
};

struct SupportLevelPair {
    std::uint32_t targetLevel{}, sourceLevel{};
    std::size_t interactions{};
    std::uint64_t representedPairs{};
};

struct SupportPartition {
    // Each pair denotes a Cartesian product of DOFs in the two node spans.
    // Near entries always refer to leaves; far entries have equal levels.
    std::vector<SupportInteraction> nearInteractions, farInteractions;
    std::vector<SupportLevelPair> farLevelPairs;
    std::uint64_t representedNearPairs{}, representedFarPairs{};
};

struct SupportLevel {
    std::vector<std::uint32_t> nodes;
    // Complete geometrical level inventory, including inactive ancestors.
    double maximumSupportRadius{};
    // Only subtrees rooted at actual far endpoints require expansion storage.
    // Ancestors above the coarsest far level remain inactive. The angular-order
    // policy uses this active radius, with one electrical order per level.
    double maximumActiveSupportRadius{};
    std::size_t pressureOutgoingNodes{}, fluxOutgoingNodes{}, pressureIncomingNodes{};
    std::size_t pressurePressureFarInteractions{}, pressureFluxFarInteractions{};
    std::size_t pressurePressureOffsets{}, pressureFluxOffsets{}, unionOffsets{};
    std::size_t pressurePressureDistanceClasses{}, pressureFluxDistanceClasses{};
    // Wigner small-d depends only on polar angle. Count exact reduced integer
    // keys (sign(z), z^2/gcd(z^2,|offset|^2), |offset|^2/gcd(...)); z=0 is one class.
    std::size_t pressurePressurePolarClasses{}, pressureFluxPolarClasses{}, unionPolarClasses{};
    // Exact squared integer dyadic offsets. Axial translation builders depend
    // on this distance class; directional rotations remain separate.
    std::vector<std::uint64_t> squaredOffsetClasses;
};

struct SupportTree {
    double separationRatio{}, nearDistanceRatio{};
    std::vector<Support> pressureSupports, fluxSupports;
    std::vector<std::uint32_t> pressureOrder, fluxOrder, pressureLeaf, fluxLeaf;
    std::vector<SupportNode> nodes;
    std::vector<SupportLevel> levels;
    SupportPartition pressurePressure, pressureFlux;
    std::size_t allocatedBytes() const noexcept;
};

struct SupportTreeSettings {
    unsigned leafCapacity{32}, maximumDepth{20};
    double separationRatio{1.5}, nearDistanceRatio{2};
    std::size_t maximumInteractions{std::numeric_limits<std::size_t>::max()};
};

// P1 support is the entire incident panel star; DP0 support is one whole panel.
// Anchors select dyadic cells, while complete support bounds determine the MAC.
// There is no low-frequency kR cutoff: angular resolution is a separate policy.
SupportTree makeSupportTree(const galerkin::Space &space, SupportTreeSettings settings = {},
                            std::stop_token stop = {});

} // namespace hornsculpt::wideband
