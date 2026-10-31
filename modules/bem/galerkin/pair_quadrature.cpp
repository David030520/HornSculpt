#include "pair_quadrature.hpp"
#include "near_projection.hpp"
#include "space.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace hornsculpt::galerkin {
namespace {
Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3 multiply(Vec3 a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
double squareLength(Vec3 a) {
    return dot(a, a);
}
Vec3 localPosition(const Panel &panel, const double *barycentric) {
    return add(multiply(subtract(panel.vertices[1], panel.vertices[0]), barycentric[1]),
               multiply(subtract(panel.vertices[2], panel.vertices[0]), barycentric[2]));
}
Vec3 displacement(const Panel &test, const double *testBarycentric, const Panel &trial,
                  const double *trialBarycentric) {
    // Form differences before interpolation. Absolute barycentric sums amplify
    // partition-of-unity roundoff at remote CAD placements and can erase gaps.
    return add(
        subtract(test.vertices[0], trial.vertices[0]),
        subtract(localPosition(test, testBarycentric), localPosition(trial, trialBarycentric)));
}

double pointTriangleDistanceSquared(Vec3 p, const std::array<Vec3, 3> &triangle) {
    const auto a = triangle[0], b = triangle[1], c = triangle[2];
    const auto ab = subtract(b, a), ac = subtract(c, a), ap = subtract(p, a);
    const double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0)
        return squareLength(ap);
    const auto bp = subtract(p, b);
    const double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3)
        return squareLength(bp);
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0)
        return squareLength(subtract(p, add(a, multiply(ab, d1 / (d1 - d3)))));
    const auto cp = subtract(p, c);
    const double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6)
        return squareLength(cp);
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0)
        return squareLength(subtract(p, add(a, multiply(ac, d2 / (d2 - d6)))));
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        const double t = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return squareLength(subtract(p, add(b, multiply(subtract(c, b), t))));
    }
    const double inverse = 1 / (va + vb + vc);
    return squareLength(
        subtract(p, add(a, add(multiply(ab, vb * inverse), multiply(ac, vc * inverse)))));
}

double segmentDistanceSquared(Vec3 p, Vec3 q, Vec3 a, Vec3 b) {
    const auto u = subtract(q, p), v = subtract(b, a), w = subtract(p, a);
    const double aa = dot(u, u), bb = dot(u, v), cc = dot(v, v), dd = dot(u, w), ee = dot(v, w);
    const double determinant = aa * cc - bb * bb;
    double s =
        determinant > aa * cc * 1e-14 ? std::clamp((bb * ee - cc * dd) / determinant, 0.0, 1.0) : 0;
    double t = (bb * s + ee) / cc;
    if (t < 0) {
        t = 0;
        s = std::clamp(-dd / aa, 0.0, 1.0);
    } else if (t > 1) {
        t = 1;
        s = std::clamp((bb - dd) / aa, 0.0, 1.0);
    }
    return squareLength(subtract(add(w, multiply(u, s)), multiply(v, t)));
}

double triangleDistance(const std::array<Vec3, 3> &a, const std::array<Vec3, 3> &b) {
    double distanceSquared = std::numeric_limits<double>::infinity();
    for (unsigned i = 0; i < 3; ++i) {
        distanceSquared = std::min({distanceSquared, pointTriangleDistanceSquared(a[i], b),
                                    pointTriangleDistanceSquared(b[i], a)});
        for (unsigned j = 0; j < 3; ++j)
            distanceSquared =
                std::min(distanceSquared,
                         segmentDistanceSquared(a[i], a[(i + 1) % 3], b[j], b[(j + 1) % 3]));
    }
    return std::sqrt(distanceSquared);
}

using Barycentric = std::array<double, 3>;
struct Subtriangle {
    std::array<Barycentric, 3> vertices{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    double fraction{1};
};
Barycentric midpoint(Barycentric a, Barycentric b) {
    return {(a[0] + b[0]) * .5, (a[1] + b[1]) * .5, (a[2] + b[2]) * .5};
}
std::array<Subtriangle, 4> subdivide(const Subtriangle &triangle) {
    const auto a = triangle.vertices[0], b = triangle.vertices[1], c = triangle.vertices[2];
    const auto ab = midpoint(a, b), bc = midpoint(b, c), ca = midpoint(c, a);
    const double fraction = triangle.fraction * .25;
    return {{{{a, ab, ca}, fraction},
             {{ab, b, bc}, fraction},
             {{ca, bc, c}, fraction},
             {{ab, bc, ca}, fraction}}};
}
std::array<Vec3, 3> physical(const Panel &panel, const Subtriangle &triangle, Vec3 origin) {
    const auto offset = subtract(panel.vertices[0], origin);
    return {add(offset, localPosition(panel, triangle.vertices[0].data())),
            add(offset, localPosition(panel, triangle.vertices[1].data())),
            add(offset, localPosition(panel, triangle.vertices[2].data()))};
}
double diameter(const std::array<Vec3, 3> &p) {
    return std::sqrt(
        std::max({squareLength(subtract(p[0], p[1])), squareLength(subtract(p[1], p[2])),
                  squareLength(subtract(p[2], p[0]))}));
}
PairQuadraturePoint mappedPoint(const TriangleQuadraturePoint &x, const TriangleQuadraturePoint &y,
                                const Subtriangle &test, const Subtriangle &trial) {
    PairQuadraturePoint result{};
    result.weight = x.weight * y.weight * test.fraction * trial.fraction;
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            result.test[i] += x.barycentric[j] * test.vertices[j][i];
            result.trial[i] += y.barycentric[j] * trial.vertices[j][i];
        }
    return result;
}

// Nine P1/P1 moments of G and its three dimensionless gradient components
// control S, D, K' and the weakly regularized H integrand with one rule.
using Moments = std::array<std::complex<double>, 36>;
inline std::array<std::complex<double>, 4> pointKernels(Vec3 delta, double k, double scale) {
    const double radius = std::sqrt(squareLength(delta));
    if (!(radius > 0))
        throw std::invalid_argument(
            "Disjoint quadrature panels intersect or contain duplicate geometric nodes.");
    const auto green = std::exp(std::complex<double>(0, -k * radius)) / radius;
    const auto derivative =
        -green * std::complex<double>(1, k * radius) * (scale / (radius * radius));
    return {green, derivative * delta.x, derivative * delta.y, derivative * delta.z};
}
void accumulateMoment(Moments &result, const Panel &test, const Panel &trial,
                      const PairQuadraturePoint &point, double k, double scale) {
    const auto kernels = pointKernels(displacement(test, point.test, trial, point.trial), k, scale);
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            const double weight = point.weight * point.test[i] * point.trial[j];
            for (unsigned component = 0; component < 4; ++component)
                result[component * 9 + i * 3 + j] += kernels[component] * weight;
        }
}
Moments moments(const Panel &test, const Panel &trial, const Subtriangle &tx, const Subtriangle &ty,
                const std::vector<TriangleQuadraturePoint> &rule, double k, double scale) {
    Moments result{};
    for (const auto &x : rule)
        for (const auto &y : rule) {
            accumulateMoment(result, test, trial, mappedPoint(x, y, tx, ty), k, scale);
        }
    return result;
}

Moments projectionMoments(const Panel &test, const Panel &trial, const Subtriangle &tx,
                          const std::vector<TriangleQuadraturePoint> &outer, unsigned order,
                          unsigned intervals, double k, double scale, std::stop_token stop) {
    PairQuadraturePatch patch{};
    patch.projection = 1;
    patch.areaFraction = tx.fraction;
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            patch.test[3 * i + j] = tx.vertices[i][j];
    static const auto gauss14 = gaussLegendre01(14), gauss20 = gaussLegendre01(20);
    const auto &base = order == 14 ? gauss14 : gauss20;
    std::vector<GaussPoint> gauss;
    gauss.reserve(base.size() * intervals);
    for (unsigned interval = 0; interval < intervals; ++interval)
        for (const auto &g : base)
            gauss.push_back({(interval + g.position) / intervals, g.weight / intervals});
    Moments result{};
    const auto offset = subtract(test.vertices[0], trial.vertices[0]);
    for (const auto &x : outer) {
        const auto projected = makeProjectionOuter(x, patch, test, trial);
        const auto localTest = localPosition(test, projected.test);
        std::array<std::complex<double>, 12> inner{};
        for (unsigned edge = 0; edge < 3; ++edge) {
            const auto geometry = makeProjectionEdge(projected, edge, trial);
            for (const auto &u : gauss) {
                if (stop.stop_requested())
                    throw std::runtime_error("Near-pair quadrature cancelled.");
                const auto ray = makeProjectionRay(projected, geometry, u);
                for (const auto &v : gauss) {
                    const auto point = mapProjectionRay(projected, ray, v);
                    const auto delta =
                        add(offset, subtract(localTest, localPosition(trial, point.trial)));
                    const auto kernels = pointKernels(delta, k, scale);
                    for (unsigned j = 0; j < 3; ++j) {
                        const double weight = point.weight * point.trial[j];
                        for (unsigned c = 0; c < 4; ++c)
                            inner[3 * c + j] += kernels[c] * weight;
                    }
                }
            }
        }
        // The test basis is constant throughout the complete inner integral.
        // Attaching it here replaces 36 per-point complex accumulations by 12,
        // while retaining every original integration point and Jacobian.
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j)
                for (unsigned c = 0; c < 4; ++c)
                    result[9 * c + 3 * i + j] += projected.test[i] * inner[3 * c + j];
    }
    return result;
}

// The inner near field changes most rapidly where the projected trial edges
// cross the outer triangle. Split the integration domain on these affine lines,
// keeping BOTH half-polygons. This does not clip or modify the physical mesh.
std::vector<Subtriangle> projectionPartition(const Panel &test, const Panel &trial) {
    std::array<Barycentric, 3> projected{};
    for (unsigned i = 0; i < 3; ++i) {
        const auto relative = subtract(test.vertices[i], trial.vertices[0]);
        projected[i][1] = dot(relative, trial.gradients[1]);
        projected[i][2] = dot(relative, trial.gradients[2]);
        projected[i][0] = 1 - projected[i][1] - projected[i][2];
    }
    using Polygon = std::vector<Barycentric>;
    std::vector<Polygon> polygons{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    for (unsigned edge = 0; edge < 3; ++edge) {
        auto value = [&](const Barycentric &point) {
            return point[0] * projected[0][edge] + point[1] * projected[1][edge] +
                   point[2] * projected[2][edge];
        };
        std::vector<Polygon> children;
        for (const auto &polygon : polygons) {
            bool positive = false, negative = false;
            for (const auto &p : polygon) {
                positive |= value(p) > 0;
                negative |= value(p) < 0;
            }
            if (!(positive && negative)) {
                children.push_back(polygon);
                continue;
            }
            Polygon plus, minus;
            for (unsigned i = 0; i < polygon.size(); ++i) {
                const auto a = polygon[i], b = polygon[(i + 1) % polygon.size()];
                const double va = value(a), vb = value(b);
                if (va >= 0)
                    plus.push_back(a);
                if (va <= 0)
                    minus.push_back(a);
                if ((va < 0 && vb > 0) || (va > 0 && vb < 0)) {
                    const double t = va / (va - vb);
                    Barycentric crossing{};
                    for (unsigned j = 0; j < 3; ++j)
                        crossing[j] = a[j] + t * (b[j] - a[j]);
                    plus.push_back(crossing);
                    minus.push_back(crossing);
                }
            }
            if (plus.size() >= 3)
                children.push_back(std::move(plus));
            if (minus.size() >= 3)
                children.push_back(std::move(minus));
        }
        polygons = std::move(children);
    }
    std::vector<Subtriangle> result;
    double area{};
    for (const auto &polygon : polygons)
        for (unsigned i = 1; i + 1 < polygon.size(); ++i) {
            const auto a = polygon[0], b = polygon[i], c = polygon[i + 1];
            const double fraction =
                std::abs((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]));
            if (fraction > 0) {
                result.push_back({{a, b, c}, fraction});
                area += fraction;
            }
        }
    if (std::abs(area - 1) > 1024 * std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Projected near-pair partition lost reference area.");
    return result;
}

bool useProjection(const Panel &test, const Panel &trial, double scale) {
    // Below 1/40 panel diameter a Cartesian grid resolves a narrow 4D ridge.
    // A projected fan is nonsingular for every OPEN-interior test point when
    // the triangle is on one side of the trial plane. Zero heights at vertices
    // or an edge are admissible; no parallel-normal condition is needed.
    // This selects integration coordinates, not an accuracy or geometry change.
    if (panelDistance(test, trial) >= .025 * scale)
        return false;
    const double roundoff = 128 * std::numeric_limits<double>::epsilon() * scale;
    double minimum = std::numeric_limits<double>::infinity(), maximum = -minimum;
    for (unsigned i = 0; i < 3; ++i) {
        const double height = dot(subtract(test.vertices[i], trial.vertices[0]), trial.normal);
        minimum = std::min(minimum, height);
        maximum = std::max(maximum, height);
    }
    // This tolerance classifies plane arithmetic only. Actual heights and
    // integration points are never clipped, shifted or replaced by a gap.
    return (minimum >= -roundoff && maximum > roundoff) ||
           (maximum <= roundoff && minimum < -roundoff);
}
double maxNorm(const Moments &values) {
    double result = 0;
    for (auto value : values) {
        const double magnitude = std::abs(value);
        if (!std::isfinite(magnitude))
            throw std::runtime_error("Near-pair quadrature produced a non-finite moment.");
        result = std::max(result, magnitude);
    }
    return result;
}

struct AdaptiveRule {
    const Panel &test;
    const Panel &trial;
    PairQuadratureSettings settings;
    double k;
    bool projection{};
    std::stop_token stop;
    AdaptiveQuadratureWorkspace &workspace;
    PairQuadratureRule result;
    Moments integral{};
    double absoluteError{};
    struct Leaf {
        Subtriangle test, trial;
        Moments integral{};
        double error{}, testDiameter{}, trialDiameter{}, innerError{}, outerError{};
        unsigned depth{}, degree{8}, innerIntervals{1};
    };
    struct LargestError {
        bool operator()(const Leaf &a, const Leaf &b) const {
            return a.error < b.error;
        }
    };
    std::vector<Leaf> pending, atLimit;
    std::array<Leaf, 4> requestedLeaves;
    std::array<PairQuadraturePatch, 4> requestedPatches;
    // Three projected edge lines split a triangle into at most seven polygons;
    // their fan triangulation has at most sixteen triangles. Requests remain
    // bounded to four so CPU and CUDA batch consumers share the same contract.
    std::array<Subtriangle, 16> initialParts;
    unsigned initialCount{}, nextInitial{}, requestedCount{};
    bool finished{};
    std::size_t charged{sizeof(AdaptiveRule)};
    PairQuadraturePatch patch(const Leaf &leaf) const {
        PairQuadraturePatch mappedPatch{};
        mappedPatch.degree = leaf.degree;
        mappedPatch.projection = projection ? 1 : 0;
        mappedPatch.projectionIntervals = leaf.innerIntervals;
        mappedPatch.areaFraction = leaf.test.fraction * leaf.trial.fraction;
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j) {
                mappedPatch.test[3 * i + j] = leaf.test.vertices[i][j];
                mappedPatch.trial[3 * i + j] = leaf.trial.vertices[i][j];
            }
        return mappedPatch;
    }
    AdaptiveRule(const Panel &a, const Panel &b, PairQuadratureSettings q, double frequency,
                 AdaptiveQuadratureWorkspace &budget, std::stop_token token)
        : test(a), trial(b), settings(q), k(frequency), stop(token), workspace(budget) {
        account(sizeof(AdaptiveRule));
        projection =
            useProjection(test, trial, std::max(test.maximumEdgeLength, trial.maximumEdgeLength));
        if (projection) {
            const auto parts = projectionPartition(test, trial);
            if (parts.size() > initialParts.size())
                throw std::logic_error("Projected line partition exceeds its triangle bound.");
            for (const auto &part : parts)
                initialParts[initialCount++] = part;
        } else
            initialParts[initialCount++] = {};
        requestInitial();
    }
    void requestInitial() {
        requestedCount = 0;
        while (nextInitial < initialCount && requestedCount < requestedLeaves.size()) {
            requestedLeaves[requestedCount] = leaf(initialParts[nextInitial++], {}, 0, 8);
            requestedPatches[requestedCount] = patch(requestedLeaves[requestedCount]);
            ++requestedCount;
        }
    }
    ~AdaptiveRule() {
        workspace.release(charged);
    }
    std::size_t ownedBytes() const {
        return sizeof(AdaptiveRule) + (pending.capacity() + atLimit.capacity()) * sizeof(Leaf) +
               result.patches.capacity() * sizeof(PairQuadraturePatch);
    }
    void account(std::size_t value) {
        if (value > settings.maximumWorkspaceBytes)
            throw std::runtime_error("Pair quadrature workspace budget exceeded.");
        result.peakWorkspaceBytes = std::max(result.peakWorkspaceBytes, value);
    }
    template <class T> void grow(std::vector<T> &values, std::size_t count) {
        if (count <= values.capacity())
            return;
        const auto capacity = std::max(count, values.capacity() + values.capacity() / 2 + 1);
        const auto current = ownedBytes(), oldBytes = values.capacity() * sizeof(T);
        if (capacity >
            (settings.maximumWorkspaceBytes - std::min(current, settings.maximumWorkspaceBytes)) /
                sizeof(T))
            throw std::runtime_error(
                "Pair quadrature workspace budget exceeded before allocation.");
        const auto added = capacity * sizeof(T);
        account(current + added);
        workspace.acquire(added);
        charged += added;
        try {
            values.reserve(capacity);
        } catch (...) {
            workspace.release(added);
            charged -= added;
            throw;
        }
        workspace.release(oldBytes);
        charged -= oldBytes;
        // Standard vector implementations may reserve more than requested.
        // Charge that excess immediately and reject it rather than hide it.
        if (values.capacity() > capacity) {
            const auto extra = (values.capacity() - capacity) * sizeof(T);
            workspace.acquire(extra);
            charged += extra;
        }
        account(ownedBytes());
    }
    void cancelled() const {
        if (stop.stop_requested())
            throw std::runtime_error("Near-pair quadrature cancelled.");
    }
    Leaf leaf(const Subtriangle &a, const Subtriangle &b, unsigned depth, unsigned degree,
              unsigned innerIntervals = 1) const {
        return {a,
                b,
                {},
                0,
                diameter(physical(test, a, test.vertices[0])),
                diameter(physical(trial, b, test.vertices[0])),
                0,
                0,
                depth,
                degree,
                innerIntervals};
    }
    void include(Leaf value) {
        absoluteError += value.error;
        if (!std::isfinite(absoluteError))
            throw std::runtime_error("Non-finite accumulated quadrature error.");
        for (unsigned i = 0; i < integral.size(); ++i)
            integral[i] += value.integral[i];
        const bool outerAtLimit =
            value.depth == settings.maximumSubdivisionDepth && value.degree == 14;
        const unsigned nextInner = value.innerIntervals * 2;
        const bool innerAtLimit =
            nextInner > (1u << settings.maximumSubdivisionDepth) ||
            std::uint64_t(192) * 400 * nextInner * nextInner > std::numeric_limits<unsigned>::max();
        if ((!projection && outerAtLimit) ||
            (projection && (value.innerError > value.outerError ? innerAtLimit : outerAtLimit))) {
            grow(atLimit, atLimit.size() + 1);
            atLimit.push_back(std::move(value));
        } else {
            grow(pending, pending.size() + 1);
            pending.push_back(std::move(value));
            std::push_heap(pending.begin(), pending.end(), LargestError{});
        }
    }
    void append(const Leaf &value) {
        result.maximumDepth = std::max(result.maximumDepth, value.depth);
        grow(result.patches, result.patches.size() + 1);
        result.patches.push_back(patch(value));
    }
    void advance() {
        cancelled();
        requestedCount = 0;
        if (nextInitial < initialCount) {
            requestInitial();
            return;
        }
        if (absoluteError <= settings.relativeTolerance * .2 * maxNorm(integral)) {
            for (const auto &value : atLimit)
                append(value);
            while (!pending.empty()) {
                append(pending.front());
                std::pop_heap(pending.begin(), pending.end(), LargestError{});
                pending.pop_back();
            }
            const auto reference = maxNorm(integral);
            if (!(reference > 0))
                throw std::runtime_error(
                    "Near-pair quadrature has no finite nonzero reference scale.");
            result.estimatedRelativeError = absoluteError / reference;
            if (!std::isfinite(result.estimatedRelativeError) ||
                result.estimatedRelativeError > settings.relativeTolerance)
                throw std::runtime_error(
                    "Near-pair quadrature global error estimate exceeds the requested tolerance.");
            finished = true;
            return;
        }
        if (pending.empty())
            throw std::runtime_error("Near-pair quadrature did not attain the requested error "
                                     "tolerance before the subdivision limit.");
        const auto value = pending.front();
        std::pop_heap(pending.begin(), pending.end(), LargestError{});
        pending.pop_back();
        absoluteError = std::max(0.0, absoluteError - value.error);
        for (unsigned i = 0; i < integral.size(); ++i)
            integral[i] -= value.integral[i];
        if (projection && value.innerError > value.outerError) {
            requestedLeaves[0] =
                leaf(value.test, value.trial, value.depth, value.degree, value.innerIntervals * 2);
            requestedCount = 1;
        } else if (value.depth == settings.maximumSubdivisionDepth) {
            requestedLeaves[0] =
                leaf(value.test, value.trial, value.depth, 14, value.innerIntervals);
            requestedCount = 1;
        } else {
            if (pending.size() + atLimit.size() + 4 > 50000)
                throw std::runtime_error(
                    "Near-pair quadrature exceeds the fifty-thousand-patch safety limit.");
            if (projection || value.testDiameter >= value.trialDiameter)
                for (const auto &child : subdivide(value.test))
                    requestedLeaves[requestedCount++] =
                        leaf(child, value.trial, value.depth + 1, 8, value.innerIntervals);
            else
                for (const auto &child : subdivide(value.trial))
                    requestedLeaves[requestedCount++] = leaf(value.test, child, value.depth + 1, 8);
        }
        for (unsigned i = 0; i < requestedCount; ++i)
            requestedPatches[i] = patch(requestedLeaves[i]);
    }
    void commit(std::span<const PairPatchEstimate> estimates) {
        cancelled();
        if (finished || estimates.size() != requestedCount)
            throw std::invalid_argument("Adaptive patch commit count/state mismatch.");
        for (const auto &estimate : estimates) {
            if (!(estimate.absoluteError >= 0) || !std::isfinite(estimate.absoluteError) ||
                !(estimate.innerError >= 0) || !std::isfinite(estimate.innerError) ||
                !(estimate.outerError >= 0) || !std::isfinite(estimate.outerError))
                throw std::runtime_error(
                    "Near-pair quadrature produced a non-finite error estimate.");
            (void)maxNorm(estimate.integral);
        }
        for (unsigned i = 0; i < requestedCount; ++i) {
            auto value = requestedLeaves[i];
            value.integral = estimates[i].integral;
            value.error = estimates[i].absoluteError;
            value.innerError = estimates[i].innerError;
            value.outerError = estimates[i].outerError;
            include(std::move(value));
        }
        advance();
    }
};
} // namespace

AdaptiveQuadratureWorkspace::AdaptiveQuadratureWorkspace(std::size_t maximumBytes)
    : maximum_(maximumBytes) {}
void AdaptiveQuadratureWorkspace::acquire(std::size_t bytes) {
    if (bytes > maximum_ - std::min(live_, maximum_))
        throw std::runtime_error(
            "Adaptive pair batch workspace budget exceeded before allocation.");
    live_ += bytes;
    peak_ = std::max(peak_, live_);
}
void AdaptiveQuadratureWorkspace::release(std::size_t bytes) noexcept {
    live_ -= bytes;
}
std::size_t AdaptiveQuadratureWorkspace::liveBytes() const noexcept {
    return live_;
}
std::size_t AdaptiveQuadratureWorkspace::peakBytes() const noexcept {
    return peak_;
}
struct AdaptivePairQuadratureState::Implementation {
    AdaptiveRule rule;
    Implementation(const Panel &test, const Panel &trial, const PairQuadratureSettings &settings,
                   double k, AdaptiveQuadratureWorkspace &workspace, std::stop_token stop)
        : rule(test, trial, settings, k, workspace, stop) {}
};
AdaptivePairQuadratureState::AdaptivePairQuadratureState(const Panel &test, const Panel &trial,
                                                         const PairQuadratureSettings &settings,
                                                         double k,
                                                         AdaptiveQuadratureWorkspace &workspace,
                                                         std::stop_token stop) {
    if (!std::isfinite(k) || k < 0 || !std::isfinite(settings.relativeTolerance) ||
        settings.relativeTolerance <= 0 || settings.maximumSubdivisionDepth > 16 ||
        classifyPair(test, trial).topology != PairTopology::disjoint)
        throw std::invalid_argument("Invalid adaptive pair settings or non-disjoint topology.");
    if (sizeof(Implementation) > settings.maximumWorkspaceBytes)
        throw std::runtime_error("Pair quadrature workspace budget exceeded before allocation.");
    workspace.acquire(sizeof(Implementation));
    try {
        impl_ = std::make_unique<Implementation>(test, trial, settings, k, workspace, stop);
    } catch (...) {
        workspace.release(sizeof(Implementation));
        throw;
    }
}
AdaptivePairQuadratureState::~AdaptivePairQuadratureState() = default;
std::span<const PairQuadraturePatch> AdaptivePairQuadratureState::requests() const noexcept {
    return std::span(impl_->rule.requestedPatches).first(impl_->rule.requestedCount);
}
void AdaptivePairQuadratureState::commit(std::span<const PairPatchEstimate> estimates) {
    impl_->rule.commit(estimates);
}
bool AdaptivePairQuadratureState::done() const noexcept {
    return impl_->rule.finished;
}
const PairQuadratureRule &AdaptivePairQuadratureState::result() const {
    if (!done())
        throw std::logic_error("Adaptive quadrature result is not complete.");
    return impl_->rule.result;
}
PairQuadratureRule AdaptivePairQuadratureState::takeResult() {
    if (!done())
        throw std::logic_error("Adaptive quadrature result is not complete.");
    auto &rule = impl_->rule;
    rule.workspace.release(rule.result.patches.capacity() * sizeof(PairQuadraturePatch));
    rule.charged -= rule.result.patches.capacity() * sizeof(PairQuadraturePatch);
    return std::move(rule.result);
}
PairPatchEstimate estimatePairPatchCpu(const Panel &test, const Panel &trial,
                                       const PairQuadraturePatch &patch, double k,
                                       std::stop_token stop) {
    if ((patch.degree != 8 && patch.degree != 14) || !std::isfinite(k) || k < 0)
        throw std::invalid_argument("Invalid CPU adaptive patch estimate.");
    static const auto low6 = triangleRule(6), high8 = triangleRule(8), low12 = triangleRule(12),
                      high14 = triangleRule(14);
    Subtriangle a, b;
    a.fraction = patch.areaFraction;
    b.fraction = 1;
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            a.vertices[i][j] = patch.test[3 * i + j];
            b.vertices[i][j] = patch.trial[3 * i + j];
        }
    const auto scale = std::max(test.maximumEdgeLength, trial.maximumEdgeLength);
    const auto &outerLow = patch.degree == 8 ? low6 : low12;
    const auto &outerHigh = patch.degree == 8 ? high8 : high14;
    const auto low = patch.projection ? projectionMoments(test, trial, a, outerLow, 20,
                                                          patch.projectionIntervals, k, scale, stop)
                                      : moments(test, trial, a, b, outerLow, k, scale);
    const auto high = patch.projection
                          ? projectionMoments(test, trial, a, outerHigh, 20,
                                              patch.projectionIntervals, k, scale, stop)
                          : moments(test, trial, a, b, outerHigh, k, scale);
    Moments difference{};
    for (unsigned i = 0; i < 36; ++i)
        difference[i] = high[i] - low[i];
    (void)maxNorm(high);
    const double outerError = maxNorm(difference);
    double innerError{};
    if (patch.projection) {
        const auto innerLow = projectionMoments(test, trial, a, outerHigh, 14,
                                                patch.projectionIntervals, k, scale, stop);
        for (unsigned i = 0; i < 36; ++i)
            difference[i] = high[i] - innerLow[i];
        innerError = maxNorm(difference);
    }
    return {high, outerError + innerError, innerError, outerError};
}

std::vector<GaussPoint> gaussLegendre01(unsigned order) {
    if (order < 2 || order > 20)
        throw std::invalid_argument("Gauss order must be between 2 and 20.");
    std::vector<GaussPoint> result(order);
    for (unsigned i = 0; i < (order + 1) / 2; ++i) {
        double x = std::cos(std::numbers::pi * (static_cast<double>(i) + .75) /
                            (static_cast<double>(order) + .5)),
               derivative = 0;
        for (unsigned iteration = 0; iteration < 32; ++iteration) {
            double p0 = 1, p1 = x;
            for (unsigned degree = 2; degree <= order; ++degree) {
                const double p = ((2 * degree - 1) * x * p1 - (degree - 1) * p0) / degree;
                p0 = p1;
                p1 = p;
            }
            derivative = order * (x * p1 - p0) / (x * x - 1);
            const double change = p1 / derivative;
            x -= change;
            if (std::abs(change) < 4 * std::numeric_limits<double>::epsilon())
                break;
            if (iteration == 31)
                throw std::runtime_error("Gauss-Legendre roots did not converge.");
        }
        const double weight = 1 / ((1 - x * x) * derivative * derivative);
        result[i] = {(1 - x) * .5, weight};
        result[order - 1 - i] = {(1 + x) * .5, weight};
    }
    return result;
}

std::vector<TriangleQuadraturePoint> triangleRule(unsigned degree) {
    if (degree == 12 || degree == 14) {
        // x=u, y=(1-u)v: the affine-triangle Jacobian raises the
        // polynomial degree in u by one, hence order=(degree+2)/2.
        const auto gauss = gaussLegendre01((degree + 2) / 2);
        std::vector<TriangleQuadraturePoint> result;
        result.reserve(gauss.size() * gauss.size());
        for (const auto &u : gauss)
            for (const auto &v : gauss)
                result.push_back({{(1 - u.position) * (1 - v.position), u.position,
                                   (1 - u.position) * v.position},
                                  u.weight * v.weight * (1 - u.position)});
        return result;
    }
    if (degree == 2)
        return {{{2.0 / 3, 1.0 / 6, 1.0 / 6}, 1.0 / 6},
                {{1.0 / 6, 2.0 / 3, 1.0 / 6}, 1.0 / 6},
                {{1.0 / 6, 1.0 / 6, 2.0 / 3}, 1.0 / 6}};
    std::vector<TriangleQuadraturePoint> result;
    // Fixed rule payloads are budgeted before construction by the adaptive
    // state; avoid implicit vector growth and its transient extra allocation.
    if (degree == 4 || degree == 6 || degree == 8)
        result.reserve(degree == 4 ? 6 : degree == 6 ? 12 : 16);
    auto three = [&](double repeated, double other, double weight) {
        for (unsigned i = 0; i < 3; ++i) {
            TriangleQuadraturePoint point{{repeated, repeated, repeated}, weight};
            point.barycentric[i] = other;
            result.push_back(point);
        }
    };
    auto six = [&](double a, double b, double c, double weight) {
        const double coordinates[3]{a, b, c};
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j)
                if (i != j)
                    result.push_back(
                        {{coordinates[i], coordinates[j], coordinates[3 - i - j]}, weight});
    };
    if (degree == 4) {
        three(.445948490915965, .108103018168070, .1116907948390055);
        three(.091576213509771, .816847572980459, .054975871827661);
    } else if (degree == 6) {
        three(.249286745170910, .501426509658179, .0583931378631895);
        three(.063089014491502, .873821971016996, .0254224531851035);
        six(.310352451033785, .636502499121399, .053145049844816, .041425537809187);
    } else if (degree == 8) {
        result.push_back({{1.0 / 3, 1.0 / 3, 1.0 / 3}, .0721578038388935});
        three(.459292588292723, .081414823414554, .0475458171336425);
        three(.170569307751760, .658861384496480, .051608685267359);
        three(.050547228317031, .898905543365938, .016229248811599);
        six(.263112829634638, .728492392955404, .008394777409958, .0136151570872175);
    } else
        throw std::invalid_argument("Triangle quadrature degree must be 2, 4, 6, 8, 12 or 14.");
    return result;
}

PairPermutation classifyPair(const Panel &test, const Panel &trial) {
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = i + 1; j < 3; ++j)
            if (test.nodes[i] == test.nodes[j] || trial.nodes[i] == trial.nodes[j])
                throw std::invalid_argument("Panel contains repeated node IDs.");
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            if (test.nodes[i] == trial.nodes[j]) {
                const auto delta = subtract(test.vertices[i], trial.vertices[j]);
                if (squareLength(delta) != 0)
                    throw std::invalid_argument("A shared node ID has inconsistent coordinates.");
            }
    return classifyPairIds(test.nodes.data(), trial.nodes.data());
}

bool isNearPair(const Panel &test, const Panel &trial, double distanceRatio) {
    if (!std::isfinite(distanceRatio) || distanceRatio < 0)
        throw std::invalid_argument("Near distance ratio must be finite and nonnegative.");
    if (distanceRatio == 0)
        return false;
    if (classifyPairIds(test.nodes.data(), trial.nodes.data()).topology != PairTopology::disjoint)
        return false;
    const double size = std::max(test.maximumEdgeLength, trial.maximumEdgeLength);
    const double bound = distanceRatio * size + test.maximumEdgeLength + trial.maximumEdgeLength;
    if (squareLength(subtract(test.centroid, trial.centroid)) >= bound * bound)
        return false;
    return panelDistance(test, trial) < distanceRatio * size;
}

double panelDistance(const Panel &test, const Panel &trial) {
    std::array<Vec3, 3> localTest{}, localTrial{};
    for (unsigned i = 0; i < 3; ++i) {
        localTest[i] = subtract(test.vertices[i], test.vertices[0]);
        localTrial[i] = subtract(trial.vertices[i], test.vertices[0]);
    }
    return triangleDistance(localTest, localTrial);
}

std::vector<PairQuadraturePoint> referencePairRule(PairTopology topology, unsigned order) {
    const auto gauss = gaussLegendre01(order);
    std::vector<PairQuadraturePoint> result;
    result.reserve(pairSectorCount(topology) * order * order * order * order);
    for (unsigned sector = 0; sector < pairSectorCount(topology); ++sector)
        for (const auto &r : gauss)
            for (const auto &u : gauss)
                for (const auto &v : gauss)
                    for (const auto &w : gauss) {
                        auto point = canonicalPairTransform(topology, sector, r.position,
                                                            u.position, v.position, w.position);
                        point.weight *= r.weight * u.weight * v.weight * w.weight;
                        result.push_back(point);
                    }
    return result;
}

PairQuadratureRule makePairQuadrature(const Panel &test, const Panel &trial,
                                      const PairQuadratureSettings &settings, double wavenumber,
                                      std::stop_token stop, bool forceAdaptive) {
    if (stop.stop_requested())
        throw std::runtime_error("Near-pair quadrature cancelled.");
    if ((settings.regularOrder != 2 && settings.regularOrder != 4) || settings.singularOrder < 2 ||
        settings.singularOrder > 20 || !std::isfinite(wavenumber) || wavenumber < 0 ||
        !std::isfinite(settings.nearDistanceRatio) || settings.nearDistanceRatio < 0 ||
        !std::isfinite(settings.relativeTolerance) || settings.relativeTolerance <= 0 ||
        settings.maximumSubdivisionDepth > 16)
        throw std::invalid_argument("Invalid triangle-pair quadrature settings.");
    const auto permutation = classifyPair(test, trial);
    PairQuadratureRule result;
    result.topology = permutation.topology;
    if (permutation.topology != PairTopology::disjoint) {
        const auto count = std::size_t(pairSectorCount(permutation.topology)) *
                           settings.singularOrder * settings.singularOrder *
                           settings.singularOrder * settings.singularOrder;
        const auto expected =
            count * sizeof(PairQuadraturePoint) + settings.singularOrder * sizeof(GaussPoint);
        if (expected > settings.maximumWorkspaceBytes)
            throw std::runtime_error(
                "Pair quadrature workspace budget exceeded before allocation.");
        result.points = referencePairRule(permutation.topology, settings.singularOrder);
        result.peakWorkspaceBytes = result.points.capacity() * sizeof(PairQuadraturePoint) +
                                    settings.singularOrder * sizeof(GaussPoint);
        if (result.peakWorkspaceBytes > settings.maximumWorkspaceBytes)
            throw std::runtime_error("Pair quadrature actual workspace exceeds budget.");
        for (auto &point : result.points)
            point = permutePairPoint(point, permutation);
        return result;
    }
    if (!forceAdaptive && !isNearPair(test, trial, settings.nearDistanceRatio)) {
        const auto count = settings.regularOrder == 2 ? 3u : 6u;
        const auto expected = count * sizeof(TriangleQuadraturePoint) +
                              std::size_t(count) * count * sizeof(PairQuadraturePoint);
        if (expected > settings.maximumWorkspaceBytes)
            throw std::runtime_error(
                "Pair quadrature workspace budget exceeded before allocation.");
        const auto regular = triangleRule(settings.regularOrder);
        result.points.reserve(std::size_t(count) * count);
        result.peakWorkspaceBytes = regular.capacity() * sizeof(TriangleQuadraturePoint) +
                                    result.points.capacity() * sizeof(PairQuadraturePoint);
        if (result.peakWorkspaceBytes > settings.maximumWorkspaceBytes)
            throw std::runtime_error("Pair quadrature actual workspace exceeds budget.");
        for (const auto &x : regular)
            for (const auto &y : regular)
                result.points.push_back(mappedPoint(x, y, {}, {}));
        return result;
    }
    AdaptiveQuadratureWorkspace workspace(settings.maximumWorkspaceBytes);
    AdaptivePairQuadratureState state(test, trial, settings, wavenumber, workspace, stop);
    while (!state.done()) {
        std::array<PairPatchEstimate, 4> estimates;
        const auto requests = state.requests();
        for (std::size_t i = 0; i < requests.size(); ++i) {
            if (stop.stop_requested())
                throw std::runtime_error("Near-pair quadrature cancelled.");
            estimates[i] = estimatePairPatchCpu(test, trial, requests[i], wavenumber, stop);
        }
        state.commit(std::span(estimates).first(requests.size()));
    }
    return state.takeResult();
}
} // namespace hornsculpt::galerkin
