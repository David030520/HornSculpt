#pragma once

#include "pair_quadrature.hpp"
#include <cstddef>
#include <cstdint>
#include <span>
#include <stop_token>
#include <vector>

namespace hornsculpt::galerkin {
struct Space;

struct TouchingInteraction {
    std::uint32_t test{}, trial{};
    PairPermutation permutation;
};

struct NearInteraction {
    std::uint32_t test{}, trial{}, firstPatch{}, patchCount{};
    // Reverse jobs reuse the same patches, exchanging test/trial barycentrics.
    std::uint32_t swapBarycentric{};
};

struct InteractionPlan {
    // Ordered (test,trial), including self and both directions of touching pairs.
    std::vector<TouchingInteraction> touching;
    // Ordered (test,trial), disjoint pairs satisfying isNearPair only.
    std::vector<NearInteraction> near;
    std::vector<PairQuadraturePatch> patches;
    // Exact near membership: nearTrials[nearOffsets[t]..nearOffsets[t+1]) is
    // strictly increasing. All other disjoint pairs use regular quadrature.
    std::vector<std::uint32_t> nearOffsets, nearTrials;
    std::size_t candidatePairs{}; // Unordered spatial candidates, excluding touching.
    unsigned maximumDepth{};
    double maximumEstimatedRelativeError{};
};

// Space must come from buildSpace. Candidate generation is conservative for
// variable panel sizes. Each unordered near pair is integrated once; the result
// is deterministic across worker counts. Cancellation is checked between pairs
// and during indexing; quadrature failures propagate without a coarser fallback.
// maximumThreads=0 uses at most eight hardware threads; explicit values are also
// capped at eight. nearDistanceRatio=0 builds touching interactions only.
// A nonempty testPanels list restricts target rows without changing source
// integration. This is used for one sector of the explicit rotational space.
InteractionPlan makeInteractionPlan(const Space &space, const PairQuadratureSettings &settings,
                                    double maximumWavenumber, std::stop_token stop = {},
                                    unsigned maximumThreads = 0,
                                    std::span<const std::uint32_t> testPanels = {});
} // namespace hornsculpt::galerkin
