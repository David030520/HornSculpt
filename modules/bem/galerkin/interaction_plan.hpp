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
    // Ordered disjoint pairs: explicit geometric near OR not qualified regular.
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
// capped at eight. nearDistanceRatio=0 disables geometric correction only;
// numerical regular qualification still applies. The payload budget includes
// candidates, retained patches and bounded concurrent adaptive workspaces.
InteractionPlan makeInteractionPlan(const Space &space, const PairQuadratureSettings &settings,
                                    double maximumWavenumber, std::stop_token stop = {},
                                    unsigned maximumThreads = 0,
                                    std::size_t maximumHostBytes = std::size_t{1024} * 1024 * 1024,
                                    std::span<const std::uint32_t> testPanels = {});
} // namespace hornsculpt::galerkin
