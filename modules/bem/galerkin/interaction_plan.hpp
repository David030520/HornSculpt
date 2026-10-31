#pragma once

#include "pair_quadrature.hpp"
#include <cstddef>
#include <cstdint>
#include <span>
#include <stop_token>
#include <vector>

namespace hornsculpt::galerkin {
struct Space;

#if defined(__CUDACC__)
#define HS_INTERACTION_HD __host__ __device__
#else
#define HS_INTERACTION_HD
#endif
// Packed unordered distinct pair, independent of the mesh's panel count.
HS_INTERACTION_HD inline std::size_t unorderedPanelPairIndex(std::uint32_t a, std::uint32_t b) {
    const std::size_t low = a < b ? a : b, high = a < b ? b : a;
    return high * (high - 1) / 2 + low;
}
#undef HS_INTERACTION_HD

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
    // For streamed all-pair planning, one byte per unordered pair replaces
    // identical full-panel patches and their two directed near jobs. Nonzero
    // degree 8/14 was accepted by the SAME adaptive moment error estimator;
    // zero means no override. Empty for spatially indexed plans.
    std::vector<std::uint8_t> wholePanelDegrees;
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
// Completed rules use a shared output pool; worker scratch holds only the
// currently integrated pair. After workers finish, flattening reuses the
// classification and scratch allowances under the same total payload limit,
// reserves exact final sizes, and releases consumed rules and metadata promptly.
InteractionPlan makeInteractionPlan(const Space &space, const PairQuadratureSettings &settings,
                                    double maximumWavenumber, std::stop_token stop = {},
                                    unsigned maximumThreads = 0,
                                    std::size_t maximumHostBytes = std::size_t{1024} * 1024 * 1024,
                                    std::span<const std::uint32_t> testPanels = {});
} // namespace hornsculpt::galerkin
