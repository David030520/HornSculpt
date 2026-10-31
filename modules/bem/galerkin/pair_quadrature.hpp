#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <stop_token>
#include <vector>

#if defined(__CUDACC__)
#define HS_QUAD_HD __host__ __device__
#else
#define HS_QUAD_HD
#endif

namespace hornsculpt::galerkin {
struct Panel;

enum class PairTopology : std::uint8_t { disjoint, shared_vertex, shared_edge, coincident };

struct PairPermutation {
    PairTopology topology{PairTopology::disjoint};
    std::uint8_t test[3]{0, 1, 2};
    std::uint8_t trial[3]{0, 1, 2};
};

struct PairQuadraturePoint {
    double test[3]{};
    double trial[3]{};
    // Reference triangle measure: sum(weight)=1/4; physical factor=4*A_test*A_trial.
    double weight{};
};

struct TriangleQuadraturePoint {
    double barycentric[3]{};
    double weight{};
};
struct PairQuadraturePatch {
    // Each consecutive triple is one child vertex in parent barycentric coordinates.
    double test[9]{};
    double trial[9]{};
    double areaFraction{};
    // Spatial depth and polynomial accuracy are independent. Critical leaves
    // can use degree 14 without changing their geometric subdivision limit.
    unsigned degree{8};
    // A projected patch subdivides only the outer test triangle. The inner
    // triangle is integrated by an oriented projection fan and sinh maps.
    unsigned projection{};
    // Composite Gauss intervals in each transformed inner coordinate. Their
    // error is estimated separately from outer spatial/polynomial refinement.
    unsigned projectionIntervals{1};
};

HS_QUAD_HD inline PairPermutation classifyPairIds(const std::uint32_t *test,
                                                  const std::uint32_t *trial) {
    PairPermutation result;
    std::uint32_t shared[3]{};
    unsigned count = 0;
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            if (test[i] == trial[j] && count < 3)
                shared[count++] = test[i];
    for (unsigned i = 0; i < count; ++i)
        for (unsigned j = i + 1; j < count; ++j)
            if (shared[j] < shared[i]) {
                const auto value = shared[i];
                shared[i] = shared[j];
                shared[j] = value;
            }
    result.topology = count == 3   ? PairTopology::coincident
                      : count == 2 ? PairTopology::shared_edge
                      : count == 1 ? PairTopology::shared_vertex
                                   : PairTopology::disjoint;
    unsigned tx = 0, ty = 0;
    for (unsigned s = 0; s < count; ++s)
        for (unsigned i = 0; i < 3; ++i) {
            if (test[i] == shared[s])
                result.test[tx++] = static_cast<std::uint8_t>(i);
            if (trial[i] == shared[s])
                result.trial[ty++] = static_cast<std::uint8_t>(i);
        }
    for (unsigned i = 0; i < 3; ++i) {
        bool testShared = false, trialShared = false;
        for (unsigned s = 0; s < count; ++s) {
            testShared = testShared || test[i] == shared[s];
            trialShared = trialShared || trial[i] == shared[s];
        }
        if (!testShared)
            result.test[tx++] = static_cast<std::uint8_t>(i);
        if (!trialShared)
            result.trial[ty++] = static_cast<std::uint8_t>(i);
    }
    return result;
}

HS_QUAD_HD inline PairQuadraturePoint mapPairPatch(const TriangleQuadraturePoint &x,
                                                   const TriangleQuadraturePoint &y,
                                                   const PairQuadraturePatch &patch) {
    PairQuadraturePoint point{};
    point.weight = x.weight * y.weight * patch.areaFraction;
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            point.test[i] += x.barycentric[j] * patch.test[3 * j + i];
            point.trial[i] += y.barycentric[j] * patch.trial[3 * j + i];
        }
    return point;
}

HS_QUAD_HD inline unsigned pairSectorCount(PairTopology topology) {
    return topology == PairTopology::coincident      ? 6
           : topology == PairTopology::shared_edge   ? 5
           : topology == PairTopology::shared_vertex ? 2
                                                     : 1;
}

// Sauter-Schwab domains on {(s,t):0<=t<=s<=1}. Weight is the absolute
// four-dimensional Jacobian only; multiply by the four Gauss weights.
HS_QUAD_HD inline PairQuadraturePoint canonicalPairTransform(PairTopology topology, unsigned sector,
                                                             double r, double u, double v,
                                                             double w) {
    double xs{}, xt{}, ys{}, yt{}, jacobian{};
    if (topology == PairTopology::disjoint) {
        xs = r;
        xt = r * u;
        ys = v;
        yt = v * w;
        jacobian = r * v;
    } else if (topology == PairTopology::shared_vertex) {
        xs = r;
        xt = r * u;
        ys = r * v;
        yt = r * v * w;
        jacobian = r * r * r * v;
        if (sector == 1) {
            const double s = xs, t = xt;
            xs = ys;
            xt = yt;
            ys = s;
            yt = t;
        }
    } else if (topology == PairTopology::shared_edge) {
        jacobian = r * r * r * u * u;
        if (sector == 0) {
            xs = r;
            xt = r * u * w;
            ys = r * (1 - u * v);
            yt = r * u * (1 - v);
        } else if (sector == 1) {
            xs = r;
            xt = r * u;
            ys = r * (1 - u * v * w);
            yt = r * u * v * (1 - w);
            jacobian *= v;
        } else if (sector == 2) {
            xs = r * (1 - u * v);
            xt = r * u * (1 - v);
            ys = r;
            yt = r * u * v * w;
            jacobian *= v;
        } else if (sector == 3) {
            xs = r * (1 - u * v * w);
            xt = r * u * v * (1 - w);
            ys = r;
            yt = r * u;
            jacobian *= v;
        } else {
            xs = r * (1 - u * v * w);
            xt = r * u * (1 - v * w);
            ys = r;
            yt = r * u * v;
            jacobian *= v;
        }
    } else {
        jacobian = r * r * r * u * u * v;
        const unsigned domain = sector % 3;
        if (domain == 0) {
            xs = r;
            xt = r * (1 - u + u * v);
            ys = r * (1 - u * v * w);
            yt = r * (1 - u);
        } else if (domain == 1) {
            xs = r;
            xt = r * u * (1 - v + v * w);
            ys = r * (1 - u * v);
            yt = r * u * (1 - v);
        } else {
            xs = r * (1 - u * v * w);
            xt = r * u * (1 - v * w);
            ys = r;
            yt = r * u * (1 - v);
        }
        if (sector >= 3) {
            const double s = xs, t = xt;
            xs = ys;
            xt = yt;
            ys = s;
            yt = t;
        }
    }
    return {{1 - xs, xs - xt, xt}, {1 - ys, ys - yt, yt}, jacobian};
}

HS_QUAD_HD inline PairQuadraturePoint permutePairPoint(const PairQuadraturePoint &point,
                                                       const PairPermutation &permutation) {
    PairQuadraturePoint result{};
    result.weight = point.weight;
    for (unsigned i = 0; i < 3; ++i) {
        result.test[permutation.test[i]] = point.test[i];
        result.trial[permutation.trial[i]] = point.trial[i];
    }
    return result;
}

struct GaussPoint {
    double position{}, weight{};
};
std::vector<GaussPoint> gaussLegendre01(unsigned order);
// Positive-weight Dunavant rules of degree 2/4/6/8 and tensor Gauss rules of
// degree 12/14 on the reference triangle; weights sum to 1/2.
std::vector<TriangleQuadraturePoint> triangleRule(unsigned degree);
PairPermutation classifyPair(const Panel &test, const Panel &trial);
bool isNearPair(const Panel &test, const Panel &trial, double distanceRatio = 1.0);
double panelDistance(const Panel &test, const Panel &trial);
std::vector<PairQuadraturePoint> referencePairRule(PairTopology topology, unsigned order);

struct PairQuadratureSettings {
    unsigned regularOrder{4};      // Polynomial-exact triangle rule: 2 or 4.
    unsigned singularOrder{4};     // One-dimensional Gauss order, 2..20.
    double nearDistanceRatio{1.0}; // Zero explicitly disables near correction for diagnostics.
    double relativeTolerance{1e-5};
    unsigned maximumSubdivisionDepth{8};
    // Owned vector payloads, including transient old+new capacity during growth.
    // Shared immutable rule tables and stack/allocator bookkeeping are excluded.
    std::size_t maximumWorkspaceBytes{std::numeric_limits<std::size_t>::max()};
};

struct PairQuadratureRule {
    PairTopology topology{PairTopology::disjoint};
    std::vector<PairQuadraturePoint> points;
    // Near-disjoint rules use compact leaf patches instead of explicit point pairs.
    std::vector<PairQuadraturePatch> patches;
    double estimatedRelativeError{-1.0}; // Negative means no embedded estimate (far/singular rule).
    unsigned maximumDepth{};
    std::size_t peakWorkspaceBytes{};
};

struct PairPatchEstimate {
    std::array<std::complex<double>, 36> integral{};
    double absoluteError{}, innerError{}, outerError{};
};

// A shared, single-threaded bound for concurrently advancing adaptive pairs.
// Payload reservations include transient vector growth; allocator bookkeeping
// and immutable quadrature tables are excluded. No pair owns an unbounded pool.
class AdaptiveQuadratureWorkspace {
  public:
    explicit AdaptiveQuadratureWorkspace(std::size_t maximumBytes);
    void acquire(std::size_t bytes);
    void release(std::size_t bytes) noexcept;
    std::size_t liveBytes() const noexcept;
    std::size_t peakBytes() const noexcept;

  private:
    std::size_t maximum_, live_{}, peak_{};
};

// One identical largest-error state machine is driven by either CPU estimates
// (dense integration/reference) or batched CUDA estimates (support near build).
// A commit supplies all currently requested patches in their original order.
// The stopping norm is the sum of leaf errors relative to the complete 36-
// moment integral, with the original .2 safety factor and depth/patch limits.
class AdaptivePairQuadratureState {
  public:
    AdaptivePairQuadratureState(const Panel &test, const Panel &trial,
                                const PairQuadratureSettings &settings, double wavenumber,
                                AdaptiveQuadratureWorkspace &workspace, std::stop_token stop = {});
    ~AdaptivePairQuadratureState();
    AdaptivePairQuadratureState(const AdaptivePairQuadratureState &) = delete;
    AdaptivePairQuadratureState &operator=(const AdaptivePairQuadratureState &) = delete;
    std::span<const PairQuadraturePatch> requests() const noexcept;
    void commit(std::span<const PairPatchEstimate> estimates);
    bool done() const noexcept;
    const PairQuadratureRule &result() const;
    PairQuadratureRule takeResult();

  private:
    struct Implementation;
    std::unique_ptr<Implementation> impl_;
};

PairPatchEstimate estimatePairPatchCpu(const Panel &test, const Panel &trial,
                                       const PairQuadraturePatch &patch, double wavenumber,
                                       std::stop_token stop = {});

// Returns points in each panel's ORIGINAL local vertex order. Near-disjoint
// pairs use embedded degree-6/8 comparisons for G and its gradient. The sum of
// patch errors is controlled relative to the complete pair integral; refinement
// targets the largest remaining contribution. At the spatial depth limit,
// degree-12/14 estimates can promote individual leaves. Unresolved limits fail.
// forceAdaptive bypasses the geometric-far shortcut for pairs that fail a
// separate frequency-dependent regular-quadrature qualification.
PairQuadratureRule makePairQuadrature(const Panel &test, const Panel &trial,
                                      const PairQuadratureSettings &settings = {},
                                      double wavenumber = 0.0, std::stop_token stop = {},
                                      bool forceAdaptive = false);
} // namespace hornsculpt::galerkin

#undef HS_QUAD_HD
