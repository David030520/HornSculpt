#include "surface_charts.hpp"
#include "surface_sizing_field.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/intersections.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

namespace hornsculpt::detail {
namespace {
constexpr auto absent = std::numeric_limits<std::uint32_t>::max();
using EdgeKey = std::array<std::uint32_t, 2>;
using PredicateKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
PredicateKernel::Point_2 predicatePoint(ChartPoint p) {
    return {p.u, p.v};
}
bool intersect(ChartPoint a, ChartPoint b, ChartPoint c, ChartPoint d) {
    return CGAL::do_intersect(PredicateKernel::Segment_2(predicatePoint(a), predicatePoint(b)),
                              PredicateKernel::Segment_2(predicatePoint(c), predicatePoint(d)));
}
Vec3 sub(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3 mul(Vec3 a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double length(Vec3 a) {
    return std::hypot(a.x, a.y, a.z);
}
EdgeKey key(std::uint32_t a, std::uint32_t b) {
    return {std::min(a, b), std::max(a, b)};
}
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Meshing cancelled.");
}

bool finitePosition(Vec3 p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

SurfaceProjection checkedProjection(SurfaceProjection result) {
    const double normalLength = length(result.normal);
    if (!finitePosition(result.position) || !finitePosition(result.normal) ||
        !std::isfinite(normalLength) || std::abs(normalLength - 1) > 1e-6)
        throw std::runtime_error(
            "Exact surface projection returned invalid position or unit normal.");
    result.normal = mul(result.normal, 1 / normalLength);
    return result;
}

SurfaceProjection exactSupportingSurface(const SurfaceGeometry &geometry, std::uint32_t faceId,
                                         Vec3 p) {
    return checkedProjection(geometry.projectSupportingSurface(faceId, p));
}

SurfaceChartMetric tangentMetric(const SurfaceChartSample &reference, Vec3 normal) {
    if (!(dot(cross(reference.du, reference.dv), normal) > 0))
        throw std::runtime_error("Exact surface orientation disagrees with its source chart.");
    const auto du = sub(reference.du, mul(normal, dot(reference.du, normal)));
    const auto dv = sub(reference.dv, mul(normal, dot(reference.dv, normal)));
    if (!(length(cross(du, dv)) > 0))
        throw std::runtime_error("Exact surface projection collapsed its chart tangent metric.");
    return {du, dv, reference.targetEdgeLength};
}

Vec3 exactSharedPosition(const SurfaceGeometry &geometry, std::span<const std::uint32_t> faces,
                         Vec3 p) {
    if (faces.empty())
        throw std::runtime_error("An exact atlas node has no physical face.");
    if (faces.size() == 1)
        return exactSupportingSurface(geometry, faces.front(), p).position;
    if (faces.size() == 2) {
        const auto result = geometry.projectBoundary(faces, p);
        if (!finitePosition(result))
            throw std::runtime_error("Exact boundary projection returned a non-finite position.");
        return result;
    }
    const auto result = geometry.projectJunction(faces, p);
    if (!finitePosition(result))
        throw std::runtime_error("Exact junction projection returned a non-finite position.");
    return result;
}
double orient(ChartPoint a, ChartPoint b, ChartPoint c) {
    return (b.u - a.u) * (c.v - a.v) - (b.v - a.v) * (c.u - a.u);
}
std::array<double, 3> barycentric(const ChartReferenceTriangle &t, ChartPoint p) {
    const double d = orient(t.uv[0], t.uv[1], t.uv[2]);
    if (!(d > 0) || !std::isfinite(d))
        throw std::runtime_error("Surface chart contains a singular or inverted UV triangle.");
    const double b = orient(t.uv[0], p, t.uv[2]) / d;
    const double c = orient(t.uv[0], t.uv[1], p) / d;
    return {1 - b - c, b, c};
}
bool within(const std::array<double, 3> &b) {
    constexpr double tolerance = 2e-10;
    return b[0] >= -tolerance && b[1] >= -tolerance && b[2] >= -tolerance;
}
std::optional<std::array<double, 3>> roundoffBoundaryWeights(const ChartReferenceTriangle &triangle,
                                                             ChartPoint p) {
    double scale = std::max(std::abs(p.u), std::abs(p.v));
    for (const auto uv : triangle.uv)
        scale = std::max({scale, std::abs(uv.u), std::abs(uv.v)});
    const double coordinateRoundoff = 16 * std::numeric_limits<double>::epsilon() * scale;
    double distance = std::numeric_limits<double>::infinity();
    std::array<double, 3> weights{};
    for (unsigned i = 0; i < 3; ++i) {
        const auto a = triangle.uv[i], b = triangle.uv[(i + 1) % 3];
        // A boundary midpoint can round outside an extremely thin reference
        // lens. Bound coordinate subtraction and determinant arithmetic in UV
        // area units; a fixed barycentric bound is ill-conditioned here.
        // This is unrelated to the requested geometric or quality tolerances.
        const double du = b.u - a.u, dv = b.v - a.v;
        const double roundoff = coordinateRoundoff * (std::abs(du) + std::abs(dv));
        if (orient(a, b, p) < -roundoff)
            return std::nullopt;
        const double t =
            std::clamp(((p.u - a.u) * du + (p.v - a.v) * dv) / (du * du + dv * dv), 0., 1.);
        const double error = std::hypot(p.u - (a.u + t * du), p.v - (a.v + t * dv));
        if (error < distance) {
            distance = error;
            weights = {};
            weights[i] = 1 - t;
            weights[(i + 1) % 3] = t;
        }
    }
    return distance <= coordinateRoundoff ? std::optional{weights} : std::nullopt;
}
struct Edge {
    EdgeKey vertices;
    std::array<std::uint32_t, 2> faces{absent, absent};
    std::array<std::uint32_t, 2> slots{};
    bool glued{};
};
struct DisjointSet {
    std::vector<std::uint32_t> parent;
    std::vector<std::vector<std::uint32_t>> members;
    explicit DisjointSet(std::size_t n) : parent(n), members(n) {
        std::iota(parent.begin(), parent.end(), 0u);
        for (std::uint32_t i = 0; i < n; ++i)
            members[i].push_back(i);
    }
    std::uint32_t find(std::uint32_t v) {
        auto r = v;
        while (parent[r] != r)
            r = parent[r];
        while (parent[v] != v) {
            const auto next = parent[v];
            parent[v] = r;
            v = next;
        }
        return r;
    }
    std::vector<std::uint32_t> join(std::uint32_t a, std::uint32_t b) {
        a = find(a);
        b = find(b);
        if (a == b)
            return {};
        if (members[a].size() < members[b].size())
            std::swap(a, b);
        parent[b] = a;
        auto changed = std::move(members[b]);
        members[a].insert(members[a].end(), changed.begin(), changed.end());
        return changed;
    }
};
struct RawHalfedge {
    std::uint32_t from, to, edge;
    std::uint32_t sourceTriangle{}, sourceSlot{};
};
struct RawChart {
    SurfaceChart chart;
    std::vector<std::uint32_t> sourceVertices;
    std::vector<std::array<std::uint32_t, 3>> triangles;
    std::vector<ChartPoint> uv;
    std::vector<std::vector<RawHalfedge>> loops;
    bool projected{};
};

struct EdgeSample {
    double t;
    std::uint32_t atlasVertex;
};
struct CurvePoint {
    Vec3 position;
    ChartPoint uv;
    std::uint32_t atlasVertex{absent}, sourceVertex{absent};
    std::optional<CurveUseInterval> curve;
    double curveParameter{};
};
using CurvePiece = std::vector<CurvePoint>;

ChartPoint nativePeriodLift(const SurfaceChart &chart, std::uint32_t useId, double parameter,
                            ChartPoint sourceUv) {
    const auto use = chart.parameterization->curveUse(useId, parameter);
    const auto domain = chart.parameterization->domain(chart.faceId);
    const auto component = [](double delta, const std::optional<double> &period) {
        if (!period)
            return 0.;
        if (!std::isfinite(*period) || *period <= 0)
            throw std::runtime_error("Native periodic surface has an invalid period.");
        return std::nearbyint(delta / *period) * *period;
    };
    const ChartPoint lift{
        component(sourceUv.u - use.uv.u, domain.uPeriod),
        component(sourceUv.v - chart.parameterOrientation * use.uv.v, domain.vPeriod)};
    const auto original =
        chart.parameterization
            ->surface(chart.faceId, {sourceUv.u, chart.parameterOrientation * sourceUv.v})
            .position;
    const auto bound = chart.parameterization
                           ->surface(chart.faceId, {use.uv.u + lift.u,
                                                    use.uv.v + chart.parameterOrientation * lift.v})
                           .position;
    if (!finitePosition(original) || !finitePosition(bound) ||
        !std::isfinite(use.modelToleranceMetres) || use.modelToleranceMetres < 0)
        throw std::runtime_error("Native periodic-use binding contains invalid data.");
    const double roundoff =
        256 * std::numeric_limits<double>::epsilon() * std::max(length(original), length(bound));
    if (use.vertexToleranceMetres &&
        (!std::isfinite(*use.vertexToleranceMetres) || *use.vertexToleranceMetres < 0))
        throw std::runtime_error("Native curve use returned invalid source vertex tolerance.");
    const double residual = length(sub(original, bound));
    if (residual > chart.parameterAccuracy + roundoff) {
        std::ostringstream message;
        message << "Flaeche " << chart.faceId << ": UV-/Randrest " << residual
                << " m ueberschreitet geforderte Genauigkeit " << chart.parameterAccuracy
                << " m (Quelltoleranz " << use.modelToleranceMetres << " m";
        if (use.vertexToleranceMetres)
            message << ", Vertex " << *use.vertexToleranceMetres << " m";
        message << ").";
        throw std::runtime_error(message.str());
    }
    return lift;
}
ChartPoint nativeCurveUv(const SurfaceChart &chart, const CurveUseInterval &curve, double parameter,
                         ChartPoint sourceUv) {
    const auto lift = nativePeriodLift(chart, curve.useId, curve.t0, sourceUv);
    const auto point = chart.parameterization->curveUse(curve.useId, parameter);
    if (!std::isfinite(point.uv.u) || !std::isfinite(point.uv.v) ||
        !std::isfinite(point.uvDerivative.u) || !std::isfinite(point.uvDerivative.v) ||
        !std::isfinite(point.modelToleranceMetres) || point.modelToleranceMetres < 0)
        throw std::runtime_error("Native curve-use evaluation returned invalid parameters.");
    // Only an integer true period is added. Source UV/pcurve representation
    // residuals are checked above, never transported as a continuous offset.
    return {point.uv.u + lift.u, chart.parameterOrientation * point.uv.v + lift.v};
}

std::vector<CurvePiece> boundaryPieces(const RawChart &raw, const std::vector<RawHalfedge> &loop,
                                       const SurfaceMesh &source, const std::vector<Edge> &edges,
                                       const std::vector<std::vector<EdgeSample>> &samples,
                                       const std::vector<Vec3> &vertices) {
    CurvePiece path;
    for (const auto &h : loop) {
        const auto &edge = edges[h.edge];
        const bool forward = raw.sourceVertices[h.from] == edge.vertices[0];
        std::optional<CurveUseInterval> nativeCurve;
        if (raw.chart.parameterization)
            nativeCurve =
                raw.chart.parameterization->triangles()[h.sourceTriangle].edges[h.sourceSlot];
        path.push_back({source.vertices[raw.sourceVertices[h.from]], raw.uv[h.from], absent,
                        raw.sourceVertices[h.from], nativeCurve,
                        nativeCurve ? nativeCurve->t0 : 0});
        auto append = [&](const EdgeSample &node) {
            const double t = forward ? node.t : 1 - node.t;
            if (t >= 1)
                return;
            const auto a = raw.uv[h.from], b = raw.uv[h.to];
            if (t == 0)
                path.back().atlasVertex = node.atlasVertex;
            else {
                const double parameter =
                    nativeCurve ? nativeCurve->t0 + t * (nativeCurve->t1 - nativeCurve->t0) : 0;
                const auto uv = nativeCurve
                                    ? nativeCurveUv(raw.chart, *nativeCurve, parameter, a)
                                    : ChartPoint{a.u + t * (b.u - a.u), a.v + t * (b.v - a.v)};
                path.push_back({vertices[node.atlasVertex], uv, node.atlasVertex, absent,
                                nativeCurve, parameter});
            }
        };
        if (forward)
            for (const auto &node : samples[h.edge])
                append(node);
        else
            for (auto it = samples[h.edge].rbegin(); it != samples[h.edge].rend(); ++it)
                append(*it);
    }
    std::vector<std::size_t> controls;
    for (std::size_t i = 0; i < path.size(); ++i)
        if (path[i].atlasVertex != absent)
            controls.push_back(i);
    if (controls.size() < 3)
        throw std::runtime_error("Surface atlas boundary has fewer than three resampled nodes.");
    std::vector<CurvePiece> pieces;
    for (std::size_t i = 0; i < controls.size(); ++i) {
        CurvePiece piece;
        for (auto j = controls[i];; j = (j + 1) % path.size()) {
            piece.push_back(path[j]);
            if (j == controls[(i + 1) % controls.size()])
                break;
        }
        pieces.push_back(std::move(piece));
    }
    return pieces;
}

enum class LensKind { none, extension, retainKnots };
LensKind referenceLens(const RawChart &raw, const CurvePiece &piece, double tolerance,
                       std::vector<ChartReferenceTriangle> *triangles = nullptr) {
    if (raw.chart.planar || !raw.projected || piece.size() < 3)
        return LensKind::none;
    bool left = false, right = false;
    for (std::size_t i = 0; i < piece.size(); ++i) {
        const auto orientation = CGAL::orientation(
            predicatePoint(piece[i].uv), predicatePoint(piece[(i + 1) % piece.size()].uv),
            predicatePoint(piece[(i + 2) % piece.size()].uv));
        left = left || orientation == CGAL::LEFT_TURN;
        right = right || orientation == CGAL::RIGHT_TURN;
    }
    if (left && right)
        return LensKind::retainKnots;
    if (!right)
        return LensKind::none;
    const auto closingA = predicatePoint(piece.front().uv),
               closingB = predicatePoint(piece.back().uv);
    const PredicateKernel::Segment_2 closing(closingA, closingB);
    for (const auto &loop : raw.loops)
        for (const auto &edge : loop) {
            const auto a = predicatePoint(raw.uv[edge.from]), b = predicatePoint(raw.uv[edge.to]);
            const auto first = CGAL::orientation(closingA, closingB, a),
                       second = CGAL::orientation(closingA, closingB, b);
            const auto third = CGAL::orientation(a, b, closingA),
                       fourth = CGAL::orientation(a, b, closingB);
            if (first != CGAL::COLLINEAR && second != CGAL::COLLINEAR && first != second &&
                third != CGAL::COLLINEAR && fourth != CGAL::COLLINEAR && third != fourth)
                return LensKind::retainKnots;
            if (a != closingA && a != closingB && closing.has_on(a)) {
                const bool original = std::any_of(piece.begin(), piece.end(), [&](const auto &p) {
                    return predicatePoint(p.uv) == a;
                });
                if (!original)
                    return LensKind::retainKnots;
            }
            bool interior = true;
            for (std::size_t i = 0; i < piece.size(); ++i)
                if (CGAL::orientation(predicatePoint(piece[i].uv),
                                      predicatePoint(piece[(i + 1) % piece.size()].uv),
                                      a) != CGAL::RIGHT_TURN) {
                    interior = false;
                    break;
                }
            if (interior)
                return LensKind::retainKnots;
        }
    const auto chord = sub(piece.back().position, piece.front().position);
    const double squared = dot(chord, chord);
    if (!(squared > 0))
        return LensKind::retainKnots;
    for (const auto &point : piece) {
        const auto q = sub(point.position, piece.front().position);
        const double t = std::clamp(dot(q, chord) / squared, 0.0, 1.0);
        if (length(sub(q, mul(chord, t))) > tolerance / 2 * (1 + 1e-12))
            return LensKind::retainKnots;
    }
    std::vector<ChartReferenceTriangle> added;
    for (std::size_t i = 1; i + 1 < piece.size(); ++i) {
        const std::array<std::size_t, 3> indices{0, i + 1, i};
        ChartReferenceTriangle t;
        for (std::size_t j = 0; j < 3; ++j) {
            t.uv[j] = piece[indices[j]].uv;
            t.positions[j] = piece[indices[j]].position;
        }
        const auto predicate = CGAL::orientation(predicatePoint(t.uv[0]), predicatePoint(t.uv[1]),
                                                 predicatePoint(t.uv[2]));
        if (predicate == CGAL::COLLINEAR)
            continue;
        const double determinant = orient(t.uv[0], t.uv[1], t.uv[2]);
        if (!(determinant > 0))
            return LensKind::retainKnots;
        const ChartPoint center{(t.uv[0].u + t.uv[1].u + t.uv[2].u) / 3,
                                (t.uv[0].v + t.uv[1].v + t.uv[2].v) / 3};
        bool oldDomain = false;
        for (const auto &loop : raw.loops)
            for (const auto &edge : loop) {
                const auto a = raw.uv[edge.from], b = raw.uv[edge.to];
                if ((a.v > center.v) != (b.v > center.v) &&
                    center.u < (b.u - a.u) * (center.v - a.v) / (b.v - a.v) + a.u)
                    oldDomain = !oldDomain;
            }
        if (oldDomain)
            return LensKind::retainKnots;
        const auto a = sub(t.positions[1], t.positions[0]), b = sub(t.positions[2], t.positions[0]);
        const double u1 = t.uv[1].u - t.uv[0].u, v1 = t.uv[1].v - t.uv[0].v;
        const double u2 = t.uv[2].u - t.uv[0].u, v2 = t.uv[2].v - t.uv[0].v;
        const auto du = mul(sub(mul(a, v2), mul(b, v1)), 1 / determinant);
        const auto dv = mul(sub(mul(b, u1), mul(a, u2)), 1 / determinant);
        if (raw.chart.geometry) {
            const auto centroid =
                mul(add(t.positions[0], add(t.positions[1], t.positions[2])), 1. / 3);
            SurfaceProjection exact;
            try {
                exact = exactSupportingSurface(*raw.chart.geometry, raw.chart.faceId, centroid);
            } catch (const SurfaceProjectionDomainError &) {
                // An optional fan may leave a finite supporting domain. The
                // unchanged source-curve knots avoid that extension entirely.
                return LensKind::retainKnots;
            }
            const auto normal = cross(du, dv);
            // A curve fan may span a plane normal to the actual surface and
            // collapse under closest-point lifting. Keep its original curve
            // knots rather than inventing a tangent metric for that mapping.
            if (!(dot(normal, exact.normal) >= .5 * length(normal)))
                return LensKind::retainKnots;
        }
        const double g00 = dot(du, du), g01 = dot(du, dv), g11 = dot(dv, dv);
        const double jacobian = length(cross(du, dv));
        const double maximum = (g00 + g11 + std::hypot(g00 - g11, 2 * g01)) / 2;
        if (!(jacobian > 0) || !std::isfinite(maximum) || maximum / jacobian > 2 * (1 + 1e-10))
            return LensKind::retainKnots;
        added.push_back(t);
    }
    // Every added point lies in the convex hull of a source-curve arc whose
    // distance to its chord is <=tol/2. Continuity of the arc then bounds the
    // extension's distance to that original curve by tol (two chord distances).
    if (triangles)
        triangles->insert(triangles->end(), added.begin(), added.end());
    return LensKind::extension;
}

bool simpleProjectedBoundary(const RawChart &raw) {
    struct Segment {
        ChartPoint a, b;
        std::size_t loop, index, count;
    };
    std::vector<Segment> segments;
    std::vector<std::vector<PredicateKernel::Point_2>> polygons;
    for (std::size_t li = 0; li < raw.loops.size(); ++li) {
        const auto &loop = raw.loops[li];
        std::vector<PredicateKernel::Point_2> polygon;
        for (std::size_t i = 0; i < loop.size(); ++i) {
            const auto a = raw.uv[loop[i].from], b = raw.uv[loop[i].to];
            polygon.push_back(predicatePoint(a));
            segments.push_back({a, b, li, i, loop.size()});
        }
        polygons.push_back(std::move(polygon));
    }
    for (std::size_t i = 0; i < segments.size(); ++i)
        for (std::size_t j = i + 1; j < segments.size(); ++j) {
            const auto &s = segments[i], &t = segments[j];
            if (s.loop == t.loop &&
                (t.index == s.index + 1 || (s.index == 0 && t.index + 1 == s.count)))
                continue;
            const auto a = s.a, b = s.b, c = t.a, d = t.b;
            if (std::max(a.u, b.u) < std::min(c.u, d.u) ||
                std::max(c.u, d.u) < std::min(a.u, b.u) ||
                std::max(a.v, b.v) < std::min(c.v, d.v) || std::max(c.v, d.v) < std::min(a.v, b.v))
                continue;
            if (intersect(a, b, c, d))
                return false;
        }
    std::size_t outer = 0;
    for (const auto &polygon : polygons) {
        const auto orientation =
            CGAL::orientation_2(polygon.begin(), polygon.end(), PredicateKernel{});
        if (orientation == CGAL::COLLINEAR)
            return false;
        if (orientation == CGAL::COUNTERCLOCKWISE)
            ++outer;
    }
    return outer == 1;
}

bool validComponentProjection(std::span<const std::uint32_t> component, Vec3 direction,
                              const SurfaceMesh &source, const std::vector<Edge> &edges,
                              const std::vector<std::array<std::uint32_t, 3>> &faceEdges,
                              std::stop_token stop) {
    const auto &first = source.triangles[component.front()];
    const auto origin = source.vertices[first.vertices[0]],
               initial = sub(source.vertices[first.vertices[1]], origin);
    auto u = sub(initial, mul(direction, dot(initial, direction)));
    u = mul(u, 1 / length(u));
    const auto v = cross(direction, u);
    std::set<std::uint32_t> members(component.begin(), component.end());
    RawChart projected;
    std::map<std::uint32_t, std::uint32_t> local;
    auto vertex = [&](std::uint32_t vi) {
        auto [it, added] = local.emplace(vi, static_cast<std::uint32_t>(local.size()));
        if (added) {
            const auto p = sub(source.vertices[vi], origin);
            projected.uv.push_back({dot(p, u), dot(p, v)});
        }
        return it->second;
    };
    std::vector<RawHalfedge> boundary;
    for (auto f : component) {
        cancelled(stop);
        const auto &t = source.triangles[f];
        for (std::size_t i = 0; i < 3; ++i) {
            const auto &edge = edges[faceEdges[f][i]];
            const auto other = edge.faces[0] == f ? edge.faces[1] : edge.faces[0];
            if (other == absent || !members.contains(other))
                boundary.push_back(
                    {vertex(t.vertices[i]), vertex(t.vertices[(i + 1) % 3]), faceEdges[f][i]});
        }
    }
    std::map<std::uint32_t, std::size_t> outgoing;
    for (std::size_t i = 0; i < boundary.size(); ++i)
        if (!outgoing.emplace(boundary[i].from, i).second)
            return false;
    std::vector<bool> visited(boundary.size());
    for (std::size_t i = 0; i < boundary.size(); ++i)
        if (!visited[i]) {
            std::vector<RawHalfedge> loop;
            auto next = i;
            do {
                if (visited[next])
                    return false;
                visited[next] = true;
                loop.push_back(boundary[next]);
                const auto found = outgoing.find(boundary[next].to);
                if (found == outgoing.end())
                    return false;
                next = found->second;
            } while (next != i);
            if (loop.size() < 3)
                return false;
            projected.loops.push_back(std::move(loop));
        }
    return !projected.loops.empty() && simpleProjectedBoundary(projected);
}

std::optional<Vec3> fitProjectionCone(std::span<const std::uint32_t> component,
                                      const std::vector<Vec3> &normals, Vec3 convexMean,
                                      double cosine, std::stop_token stop) {
    // The closest point of the normals' convex hull gives the smallest enclosing
    // spherical cap when that cap is a hemisphere. Frank-Wolfe uses only convex
    // combinations; acceptance still checks every actual face normal explicitly.
    auto center = convexMean;
    for (unsigned iteration = 0; iteration < 512; ++iteration) {
        cancelled(stop);
        const double magnitude = length(center);
        if (!(magnitude > 1e-12))
            return std::nullopt;
        double minimum = std::numeric_limits<double>::infinity();
        Vec3 worst{};
        for (auto f : component) {
            const double value = dot(center, normals[f]);
            if (value < minimum) {
                minimum = value;
                worst = normals[f];
            }
        }
        if (minimum / magnitude >= cosine)
            return mul(center, 1 / magnitude);
        const double gap = dot(center, center) - minimum;
        const auto delta = sub(worst, center);
        const double squared = dot(delta, delta);
        if (!(squared > 0) || gap <= 1e-13 * magnitude * magnitude)
            return std::nullopt;
        center = add(center, mul(delta, std::clamp(gap / squared, 0.0, 1.0)));
    }
    return std::nullopt;
}

// A polygonal cylinder/cone is intrinsically developable even though a single
// 3D planar projection fails. Accept its unfolding only after checking every
// shared edge, closure assignment and the complete disk boundary for overlap.
bool unfoldDevelopable(RawChart &raw, const SurfaceMesh &source, std::stop_token stop) {
    if (raw.loops.size() != 1)
        return false;
    const auto vertexCount = raw.sourceVertices.size();
    std::vector<ChartPoint> uv(vertexCount);
    std::vector<bool> assigned(vertexCount), visited(raw.triangles.size());
    std::map<EdgeKey, std::vector<std::uint32_t>> incident;
    double extent = 0;
    const auto origin = source.vertices[raw.sourceVertices.front()];
    for (auto vi : raw.sourceVertices)
        extent = std::max(extent, length(sub(source.vertices[vi], origin)));
    const double rounding = 512 * std::numeric_limits<double>::epsilon() * extent;
    for (std::uint32_t f = 0; f < raw.triangles.size(); ++f)
        for (std::size_t i = 0; i < 3; ++i)
            incident[key(raw.triangles[f][i], raw.triangles[f][(i + 1) % 3])].push_back(f);
    const auto first = raw.triangles.front();
    uv[first[0]] = {0, 0};
    uv[first[1]] = {length(sub(source.vertices[raw.sourceVertices[first[1]]],
                               source.vertices[raw.sourceVertices[first[0]]])),
                    0};
    assigned[first[0]] = assigned[first[1]] = true;
    std::vector<std::uint32_t> queue{0};
    for (std::size_t next = 0; next < queue.size(); ++next) {
        if (next % 128 == 0)
            cancelled(stop);
        const auto f = queue[next];
        if (visited[f])
            continue;
        const auto &t = raw.triangles[f];
        std::size_t edge = 0;
        while (edge < 3 && !(assigned[t[edge]] && assigned[t[(edge + 1) % 3]]))
            ++edge;
        if (edge == 3)
            return false;
        const auto a = t[edge], b = t[(edge + 1) % 3], c = t[(edge + 2) % 3];
        const auto ab =
            sub(source.vertices[raw.sourceVertices[b]], source.vertices[raw.sourceVertices[a]]);
        const auto ac =
            sub(source.vertices[raw.sourceVertices[c]], source.vertices[raw.sourceVertices[a]]);
        const double physicalLength = length(ab);
        const double du = uv[b].u - uv[a].u, dv = uv[b].v - uv[a].v;
        const double planarLength = std::hypot(du, dv);
        if (std::abs(planarLength - physicalLength) > 1e-9 * physicalLength + rounding)
            return false;
        const double x = dot(ab, ac) / physicalLength;
        const double y = length(cross(ab, ac)) / physicalLength;
        const ChartPoint candidate{uv[a].u + (du * x - dv * y) / planarLength,
                                   uv[a].v + (dv * x + du * y) / planarLength};
        if (assigned[c] && std::hypot(candidate.u - uv[c].u, candidate.v - uv[c].v) >
                               1e-9 * std::max(physicalLength, length(ac)) + rounding)
            return false;
        uv[c] = candidate;
        assigned[c] = true;
        visited[f] = true;
        for (std::size_t i = 0; i < 3; ++i)
            for (auto neighbor : incident.at(key(t[i], t[(i + 1) % 3])))
                if (!visited[neighbor])
                    queue.push_back(neighbor);
    }
    if (std::find(visited.begin(), visited.end(), false) != visited.end())
        return false;
    for (const auto &t : raw.triangles) {
        if (!(orient(uv[t[0]], uv[t[1]], uv[t[2]]) > 0))
            return false;
        for (std::size_t i = 0; i < 3; ++i) {
            const auto a = t[i], b = t[(i + 1) % 3];
            const double expected = length(sub(source.vertices[raw.sourceVertices[a]],
                                               source.vertices[raw.sourceVertices[b]]));
            if (std::abs(std::hypot(uv[a].u - uv[b].u, uv[a].v - uv[b].v) - expected) >
                1e-9 * expected + rounding)
                return false;
        }
    }
    const auto &loop = raw.loops.front();
    for (std::size_t i = 0; i < loop.size(); ++i) {
        if (i % 128 == 0)
            cancelled(stop);
        const auto a = uv[loop[i].from], b = uv[loop[i].to];
        for (std::size_t j = i + 2; j < loop.size(); ++j) {
            if (i == 0 && j + 1 == loop.size())
                continue;
            const auto c = uv[loop[j].from], d = uv[loop[j].to];
            if (std::max(a.u, b.u) < std::min(c.u, d.u) ||
                std::max(c.u, d.u) < std::min(a.u, b.u) ||
                std::max(a.v, b.v) < std::min(c.v, d.v) || std::max(c.v, d.v) < std::min(a.v, b.v))
                continue;
            if (intersect(a, b, c, d))
                return false;
        }
    }
    raw.uv = std::move(uv);
    return true;
}

// The uniform graph Laplacian is symmetric positive definite after elimination
// of the convex Dirichlet border (Tutte embedding). MVC weights are not symmetric.
void embedDisk(RawChart &raw, const SurfaceMesh &source, std::stop_token stop) {
    if (raw.loops.size() != 1 || raw.loops.front().size() < 3)
        throw std::runtime_error("A non-planar chart was not cut into a topological disk.");
    const auto n = raw.sourceVertices.size();
    std::vector<std::set<std::uint32_t>> graph(n);
    for (const auto &t : raw.triangles)
        for (std::size_t i = 0; i < 3; ++i) {
            graph[t[i]].insert(t[(i + 1) % 3]);
            graph[t[(i + 1) % 3]].insert(t[i]);
        }
    raw.uv.resize(n);
    std::vector<bool> boundary(n);
    double perimeter = 0;
    for (const auto &h : raw.loops.front())
        perimeter += length(sub(source.vertices[raw.sourceVertices[h.to]],
                                source.vertices[raw.sourceVertices[h.from]]));
    if (!(perimeter > 0))
        throw std::runtime_error("A surface chart has a zero-length border.");
    double arc = 0;
    for (const auto &h : raw.loops.front()) {
        const double angle = 2 * std::numbers::pi * arc / perimeter;
        raw.uv[h.from] = {perimeter * std::cos(angle) / (2 * std::numbers::pi),
                          perimeter * std::sin(angle) / (2 * std::numbers::pi)};
        boundary[h.from] = true;
        arc += length(sub(source.vertices[raw.sourceVertices[h.to]],
                          source.vertices[raw.sourceVertices[h.from]]));
    }
    std::vector<std::uint32_t> index(n, absent), interior;
    for (std::uint32_t i = 0; i < n; ++i)
        if (!boundary[i]) {
            index[i] = static_cast<std::uint32_t>(interior.size());
            interior.push_back(i);
        }
    const auto m = interior.size();
    auto solve = [&](bool vertical) {
        std::vector<double> x(m), r(m), p(m), z(m), ap(m), rhs(m);
        for (std::size_t i = 0; i < m; ++i)
            for (auto j : graph[interior[i]])
                if (boundary[j])
                    rhs[i] += vertical ? raw.uv[j].v : raw.uv[j].u;
        auto product = [&](const std::vector<double> &a, std::vector<double> &b) {
            for (std::size_t i = 0; i < m; ++i) {
                b[i] = static_cast<double>(graph[interior[i]].size()) * a[i];
                for (auto j : graph[interior[i]])
                    if (!boundary[j])
                        b[i] -= a[index[j]];
            }
        };
        auto inner = [](const auto &a, const auto &b) {
            return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
        };
        r = rhs;
        for (std::size_t i = 0; i < m; ++i)
            z[i] = r[i] / static_cast<double>(graph[interior[i]].size());
        p = z;
        double rz = inner(r, z);
        const double tolerance = 1e-24 * std::max(inner(rhs, rhs), perimeter * perimeter * 1e-12);
        bool done = inner(r, r) <= tolerance;
        for (std::size_t iteration = 0; !done && iteration < 4 * m + 100; ++iteration) {
            if (iteration % 32 == 0)
                cancelled(stop);
            product(p, ap);
            const double denominator = inner(p, ap);
            if (!(denominator > 0) || !std::isfinite(denominator))
                throw std::runtime_error("The surface chart's Dirichlet system is singular.");
            const double alpha = rz / denominator;
            for (std::size_t i = 0; i < m; ++i) {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }
            done = inner(r, r) <= tolerance;
            if (done)
                break;
            for (std::size_t i = 0; i < m; ++i)
                z[i] = r[i] / static_cast<double>(graph[interior[i]].size());
            const double next = inner(r, z), beta = next / rz;
            for (std::size_t i = 0; i < m; ++i)
                p[i] = z[i] + beta * p[i];
            rz = next;
        }
        if (!done)
            throw std::runtime_error("The surface chart parameterization did not converge.");
        for (std::size_t i = 0; i < m; ++i) {
            if (vertical)
                raw.uv[interior[i]].v = x[i];
            else
                raw.uv[interior[i]].u = x[i];
        }
    };
    solve(false);
    solve(true);
    for (const auto &t : raw.triangles) {
        const auto a = raw.uv[t[0]], b = raw.uv[t[1]], c = raw.uv[t[2]];
        const double determinant = orient(a, b, c);
        const double scale =
            std::max({std::hypot(b.u - a.u, b.v - a.v), std::hypot(c.u - a.u, c.v - a.v),
                      std::hypot(c.u - b.u, c.v - b.v)});
        if (!(determinant > 64 * std::numeric_limits<double>::epsilon() * scale * scale))
            throw std::runtime_error(
                "The surface chart parameterization is inverted or numerically singular.");
    }
}
} // namespace

struct ChartLookup {
    double lowerU{}, lowerV{}, cellU{1}, cellV{1};
    std::uint32_t countU{1}, countV{1};
    std::vector<std::vector<std::uint32_t>> cells;
    std::size_t cell(ChartPoint p) const {
        const auto u = static_cast<std::uint32_t>(
            std::clamp(std::floor((p.u - lowerU) / cellU), 0.0, static_cast<double>(countU - 1)));
        const auto v = static_cast<std::uint32_t>(
            std::clamp(std::floor((p.v - lowerV) / cellV), 0.0, static_cast<double>(countV - 1)));
        return static_cast<std::size_t>(v) * countU + u;
    }
};

struct AtlasCurveRegistry {
    struct Side {
        std::size_t chart{}, loop{};
        std::vector<ChartPoint> uv;
        std::uint32_t curveUse{};
        ChartPoint parameterLift;
    };
    struct ParameterPath {
        // Physical curves use their canonical C(t). Artificial cuts use the
        // native UV polyline on this face, evaluated directly through S/D1.
        std::uint32_t curveId{}, faceId{};
        std::vector<double> parameters;
        std::vector<double> integrationBreaks;
        std::vector<ParameterPoint> uv;
    };
    struct Segment {
        std::uint32_t a{}, b{};
        std::vector<Vec3> positions;
        std::vector<Side> sides;
        std::optional<ParameterPath> parameterPath;
    };
    std::map<EdgeKey, Segment> segments;
    std::set<EdgeKey> retained;
    std::shared_ptr<const SurfaceGeometry> geometry;
    std::shared_ptr<const SurfaceParameterization> parameterization;
    double geometricTolerance{};
};

namespace {
AtlasCurveRegistry::ParameterPath parameterPathForPiece(const SurfaceChart &chart,
                                                        const CurvePiece &piece) {
    AtlasCurveRegistry::ParameterPath result;
    result.faceId = chart.faceId;
    result.curveId = piece.front().curve ? piece.front().curve->curveId : 0;
    for (std::size_t i = 0; i < piece.size(); ++i) {
        const auto &point = piece[i];
        result.uv.push_back({point.uv.u, chart.parameterOrientation * point.uv.v});
        if (result.curveId) {
            if (i + 1 == piece.size() && point.sourceVertex != absent && i > 0 &&
                piece[i - 1].curve)
                result.parameters.push_back(piece[i - 1].curve->t1);
            else if (point.curve && point.curve->curveId == result.curveId)
                result.parameters.push_back(point.curveParameter);
            else
                throw std::runtime_error(
                    "A native shared segment crosses an unretained curve junction.");
        } else {
            const double next =
                i ? result.parameters.back() + length(sub(point.position, piece[i - 1].position))
                  : 0;
            result.parameters.push_back(next);
        }
    }
    if (result.parameters.size() < 2 || result.parameters.front() == result.parameters.back())
        throw std::runtime_error("A native shared segment has a collapsed parameter interval.");
    const bool increasing = result.parameters.back() > result.parameters.front();
    for (std::size_t i = 1; i < result.parameters.size(); ++i)
        if (increasing ? result.parameters[i] <= result.parameters[i - 1]
                       : result.parameters[i] >= result.parameters[i - 1])
            throw std::runtime_error("Native shared-curve parameters are not strictly monotone.");
    if (result.curveId) {
        const auto partition = chart.parameterization->curvePartition(result.curveId);
        if (partition.size() < 2)
            throw std::runtime_error("Native shared curve has no complete span partition.");
        for (std::size_t i = 0; i < partition.size(); ++i)
            if (!std::isfinite(partition[i]) || (i && partition[i] <= partition[i - 1]))
                throw std::runtime_error("Native shared-curve spans are not finite and ordered.");
        const auto [low, high] = std::minmax(result.parameters.front(), result.parameters.back());
        if (low < partition.front() || high > partition.back())
            throw std::runtime_error("Native shared-curve spans do not cover its true trim.");
        result.integrationBreaks.assign(partition.begin(), partition.end());
    }
    return result;
}
AtlasCurveRegistry::Side parameterSideForPiece(const SurfaceChart &chart, const CurvePiece &piece,
                                               std::size_t chartIndex, std::size_t loopIndex,
                                               std::vector<ChartPoint> uv, double tolerance) {
    AtlasCurveRegistry::Side side{chartIndex, loopIndex, std::move(uv)};
    if (!piece.front().curve)
        return side;
    side.curveUse = piece.front().curve->useId;
    side.parameterLift =
        nativePeriodLift(chart, side.curveUse, piece.front().curveParameter, piece.front().uv);
    const auto path = parameterPathForPiece(chart, piece);
    for (std::size_t i = 0; i < piece.size(); ++i) {
        const auto use = chart.parameterization->curveUse(side.curveUse, path.parameters[i]);
        if (use.vertexToleranceMetres &&
            (!std::isfinite(*use.vertexToleranceMetres) || *use.vertexToleranceMetres < 0))
            throw std::runtime_error("Native curve use has invalid source vertex tolerance.");
        const auto surface = chart.parameterization->surface(chart.faceId, use.uv);
        const auto curve = chart.parameterization->curve(path.curveId, path.parameters[i]);
        if (!finitePosition(surface.position) || !finitePosition(curve.position) ||
            !std::isfinite(use.modelToleranceMetres) || use.modelToleranceMetres < 0)
            throw std::runtime_error("Native curve representation contains invalid source data.");
        const double residual = length(sub(surface.position, curve.position));
        if (residual > tolerance) {
            std::ostringstream message;
            message << std::setprecision(17) << "Native surface/curve representation on face "
                    << chart.faceId << ", curve " << path.curveId << ", use " << side.curveUse
                    << ", t=" << path.parameters[i] << ": residual " << residual
                    << " m exceeds requested geometric accuracy " << tolerance << " m.";
            throw std::runtime_error(message.str());
        }
    }
    return side;
}
} // namespace

namespace {
ParametricCurveSample parameterPathValue(const AtlasCurveRegistry &registry,
                                         const AtlasCurveRegistry::ParameterPath &path,
                                         std::size_t interval, double fraction) {
    const double a = path.parameters[interval], b = path.parameters[interval + 1];
    if (path.curveId)
        return registry.parameterization->curve(path.curveId, a + fraction * (b - a));
    const auto u = path.uv[interval], v = path.uv[interval + 1];
    const auto value = registry.parameterization->surface(
        path.faceId, {u.u + fraction * (v.u - u.u), u.v + fraction * (v.v - u.v)});
    return {value.position,
            add(mul(value.du, (v.u - u.u) / (b - a)), mul(value.dv, (v.v - u.v) / (b - a)))};
}
double parameterArc(const AtlasCurveRegistry &registry,
                    const AtlasCurveRegistry::ParameterPath &path, std::size_t interval, double low,
                    double high, double tolerance, std::stop_token stop) {
    const auto integrate = [&](double a, double b) {
        constexpr double nodes[]{0, .5384693101056831, .9061798459386640};
        constexpr double weights[]{.5688888888888889, .4786286704993665, .2369268850561891};
        const double midpoint = (a + b) / 2, half = (b - a) / 2;
        const auto speed = [&](double x) {
            const auto value = parameterPathValue(registry, path, interval, x);
            if (!finitePosition(value.position) || !finitePosition(value.tangent))
                throw std::runtime_error(
                    "Native shared curve has non-finite position or derivative.");
            return length(value.tangent) *
                   std::abs(path.parameters[interval + 1] - path.parameters[interval]);
        };
        double sum = weights[0] * speed(midpoint);
        for (unsigned i = 1; i < 3; ++i)
            sum += weights[i] *
                   (speed(midpoint - half * nodes[i]) + speed(midpoint + half * nodes[i]));
        return half * sum;
    };
    struct Panel {
        double low, high, value, error;
    };
    const auto panel = [&](double a, double b) {
        cancelled(stop);
        const double midpoint = (a + b) / 2;
        const double whole = integrate(a, b);
        const double split = integrate(a, midpoint) + integrate(midpoint, b);
        if (!std::isfinite(split) || !std::isfinite(whole))
            throw std::runtime_error("Native shared-curve arclength exceeds the numerical range.");
        return Panel{a, b, split, std::abs(split - whole)};
    };
    const auto less = [](const Panel &a, const Panel &b) {
        return a.error != b.error ? a.error < b.error : a.low > b.low;
    };
    // A whole-versus-halves error estimate can miss a narrow polynomial span
    // in both estimates. Seed every native span before adapting, preserving
    // the same total numerical error budget and the actual directed trim.
    std::vector<double> breaks{low, high};
    const double start = path.parameters[interval];
    const double delta = path.parameters[interval + 1] - start;
    for (const double parameter : path.integrationBreaks) {
        const double fraction = (parameter - start) / delta;
        if (fraction > low && fraction < high)
            breaks.push_back(fraction);
    }
    std::sort(breaks.begin(), breaks.end());
    if (breaks.size() > 4097)
        throw std::runtime_error("Native shared-curve spans exceed the integration panel budget.");
    std::vector<Panel> pending;
    for (std::size_t i = 1; i < breaks.size(); ++i)
        pending.push_back(panel(breaks[i - 1], breaks[i]));
    std::make_heap(pending.begin(), pending.end(), less);
    for (;;) {
        cancelled(stop);
        // The segment keeps one total error budget. Splitting a panel does
        // not force an unrelated smooth span into a needlessly tiny subbudget.
        double error = 0, value = 0;
        for (const auto &part : pending) {
            error += part.error;
            value += part.value;
        }
        if (error <= tolerance)
            return value;
        if (pending.size() == 4096)
            throw std::runtime_error(
                "Native shared-curve arclength exhausted its panel budget on face " +
                std::to_string(path.faceId) + ".");
        std::pop_heap(pending.begin(), pending.end(), less);
        const auto worst = pending.back();
        pending.pop_back();
        const double middle = (worst.low + worst.high) / 2;
        if (!(middle > worst.low && middle < worst.high))
            throw std::runtime_error(
                "Native shared-curve integration reached its representable limit.");
        pending.push_back(panel(worst.low, middle));
        std::push_heap(pending.begin(), pending.end(), less);
        pending.push_back(panel(middle, worst.high));
        std::push_heap(pending.begin(), pending.end(), less);
    }
}
struct ParameterMidpoint {
    std::size_t interval{};
    double fraction{}, parameter{};
    ParameterPoint uv;
    Vec3 position;
};
ParameterMidpoint parameterMidpoint(const AtlasCurveRegistry &registry,
                                    const AtlasCurveRegistry::Segment &segment,
                                    std::stop_token stop) {
    const auto &path = *segment.parameterPath;
    const double scale = length(sub(segment.positions.back(), segment.positions.front())) +
                         registry.geometricTolerance;
    const double tolerance = 1e-8 * scale; // numerical arc integration, in metres
    std::vector<double> arc(path.parameters.size());
    for (std::size_t i = 1; i < arc.size(); ++i)
        arc[i] =
            arc[i - 1] + parameterArc(registry, path, i - 1, 0, 1, tolerance / arc.size(), stop);
    if (!(arc.back() > 0))
        throw std::runtime_error("Native shared curve has zero arclength.");
    const double half = arc.back() / 2;
    const auto upper = std::upper_bound(arc.begin(), arc.end(), half);
    const auto interval =
        std::min(static_cast<std::size_t>(upper - arc.begin() - 1), arc.size() - 2);
    double low = 0, high = 1, fraction = .5;
    bool converged = false;
    for (unsigned iteration = 0; iteration < 48; ++iteration) {
        cancelled(stop);
        fraction = (low + high) / 2;
        const double current =
            arc[interval] + parameterArc(registry, path, interval, 0, fraction, tolerance, stop);
        if (std::abs(current - half) <= tolerance) {
            converged = true;
            break;
        }
        if (current < half)
            low = fraction;
        else
            high = fraction;
    }
    if (!converged)
        throw std::runtime_error("Native shared-curve arclength midpoint did not converge.");
    const auto a = path.uv[interval], b = path.uv[interval + 1];
    return {interval,
            fraction,
            path.parameters[interval] +
                fraction * (path.parameters[interval + 1] - path.parameters[interval]),
            {a.u + fraction * (b.u - a.u), a.v + fraction * (b.v - a.v)},
            parameterPathValue(registry, path, interval, fraction).position};
}
} // namespace

bool SurfaceChart::contains(ChartPoint p) const {
    if (!std::isfinite(p.u) || !std::isfinite(p.v))
        return false;
    bool inside = false;
    for (const auto &loop : loops) {
        const auto n = loop.vertices.size();
        for (std::size_t i = 0; i < n; ++i) {
            const auto a = loop.vertices[i].uv, b = loop.vertices[(i + 1) % n].uv;
            const double scale = std::hypot(b.u - a.u, b.v - a.v);
            if (std::abs(orient(a, b, p)) <= 1e-12 * scale * scale &&
                (p.u - a.u) * (p.u - b.u) + (p.v - a.v) * (p.v - b.v) <= 1e-12 * scale * scale)
                return true;
            if ((a.v > p.v) != (b.v > p.v) && p.u < (b.u - a.u) * (p.v - a.v) / (b.v - a.v) + a.u)
                inside = !inside;
        }
    }
    return inside;
}

bool SurfaceChart::containsParameterDomain(ChartPoint p) const {
    if (!std::isfinite(p.u) || !std::isfinite(p.v))
        return false;
    if (!parameterization)
        return true;
    const auto domain = parameterization->domain(faceId);
    const auto admitted = [](double value, const std::optional<std::array<double, 2>> &bounds,
                             const std::optional<double> &period) {
        if (period)
            return true;
        return !bounds || (value >= (*bounds)[0] && value <= (*bounds)[1]);
    };
    // No geometric containment tolerance extends a nonperiodic source domain.
    // This predicate rejects optional points; evaluations of retained/used
    // geometry keep the provider's own strict validation and diagnostics.
    return admitted(p.u, domain.uBounds, domain.uPeriod) &&
           admitted(parameterOrientation * p.v, domain.vBounds, domain.vPeriod);
}

SurfaceChartSample SurfaceChart::referenceSample(ChartPoint p) const {
    if (!std::isfinite(p.u) || !std::isfinite(p.v))
        throw std::invalid_argument("Surface chart coordinates must be finite.");
    const ChartReferenceTriangle *triangle = nullptr;
    std::array<double, 3> weights{};
    auto inspect = [&](std::uint32_t i) {
        const auto b = barycentric(reference[i], p);
        if (within(b)) {
            triangle = &reference[i];
            weights = b;
            return true;
        }
        if (const auto boundary = roundoffBoundaryWeights(reference[i], p)) {
            triangle = &reference[i];
            weights = *boundary;
            return true;
        }
        return false;
    };
    if (lookup) {
        for (auto i : lookup->cells[lookup->cell(p)])
            if (inspect(i))
                break;
    } else {
        for (std::uint32_t i = 0; i < reference.size(); ++i)
            if (inspect(i))
                break;
    }
    SurfaceChartSample result;
    if (planar && !parameterization) {
        result.position = add(planeOrigin, add(mul(planeU, p.u), mul(planeV, p.v)));
        result.du = planeU;
        result.dv = planeV;
        result.targetEdgeLength = triangle ? weights[0] * triangle->sizes[0] +
                                                 weights[1] * triangle->sizes[1] +
                                                 weights[2] * triangle->sizes[2]
                                           : reference.front().sizes[0];
    } else {
        if (!triangle) {
            if (!parameterization)
                throw std::runtime_error("A surface chart evaluation left its reference domain.");
            const auto value = parameterization->surface(faceId, {p.u, parameterOrientation * p.v});
            if (!finitePosition(value.position))
                throw std::runtime_error(
                    "Native scalar-size query returned a non-finite position.");
            const std::array<std::uint32_t, 1> patch{faceId};
            return {value.position, value.du, mul(value.dv, parameterOrientation),
                    sizing ? sizing->at(value.position, patch) : reference.front().sizes[0]};
        }
        const auto &t = *triangle;
        const auto a = sub(t.positions[1], t.positions[0]), b = sub(t.positions[2], t.positions[0]);
        const double u1 = t.uv[1].u - t.uv[0].u, v1 = t.uv[1].v - t.uv[0].v;
        const double u2 = t.uv[2].u - t.uv[0].u, v2 = t.uv[2].v - t.uv[0].v;
        const double determinant = u1 * v2 - v1 * u2;
        result.position = add(t.positions[0], add(mul(a, weights[1]), mul(b, weights[2])));
        result.du = mul(sub(mul(a, v2), mul(b, v1)), 1 / determinant);
        result.dv = mul(sub(mul(b, u1), mul(a, u2)), 1 / determinant);
        result.targetEdgeLength =
            weights[0] * t.sizes[0] + weights[1] * t.sizes[1] + weights[2] * t.sizes[2];
    }
    if (sizing) {
        const std::array<std::uint32_t, 1> patch{faceId};
        result.targetEdgeLength = sizing->at(result.position, patch);
    }
    if (!(result.targetEdgeLength > 0) || !std::isfinite(result.targetEdgeLength))
        throw std::runtime_error("Surface chart has an invalid scalar size.");
    return result;
}

double SurfaceChart::targetSize(ChartPoint p) const {
    return referenceSample(p).targetEdgeLength;
}

SurfaceChartMetric SurfaceChart::boundaryMetric(ChartPoint p) const {
    if (parameterization) {
        const auto value = sample(p);
        return {value.du, value.dv, value.targetEdgeLength};
    }
    const auto reference = referenceSample(p);
    if (!geometry)
        return {reference.du, reference.dv, reference.targetEdgeLength};
    // Only the local physical tangent plane is used here. The trimmed
    // projection must never be mistaken for a full-rank interior position map.
    const auto projection = checkedProjection(geometry->projectSurface(faceId, reference.position));
    return tangentMetric(reference, projection.normal);
}

SurfaceChartSample SurfaceChart::sample(ChartPoint p) const {
    if (parameterization) {
        const auto value = parameterization->surface(faceId, {p.u, parameterOrientation * p.v});
        const auto dv = mul(value.dv, static_cast<double>(parameterOrientation));
        if (!finitePosition(value.position) || !finitePosition(value.du) || !finitePosition(dv) ||
            !(length(cross(value.du, dv)) > 0))
            throw std::runtime_error(
                "Native surface evaluation is non-finite or singular on face " +
                std::to_string(faceId) + ".");
        const std::array<std::uint32_t, 1> patch{faceId};
        const double h =
            sizing ? sizing->at(value.position, patch) : referenceSample(p).targetEdgeLength;
        return {value.position, value.du, dv, h};
    }
    auto result = referenceSample(p);
    if (geometry) {
        // A chordal chart can extend across a curved trim. Clamping its
        // interior to that trim collapses a two-dimensional strip to a curve.
        const auto projection = exactSupportingSurface(*geometry, faceId, result.position);
        const auto metric = tangentMetric(result, projection.normal);
        result.position = projection.position;
        // The PL chart supplies coordinates, while the exact normal supplies a
        // local tangent metric. This is not an exact CAD-parameter derivative.
        result.du = metric.du;
        result.dv = metric.dv;
    }
    return result;
}

namespace {
class SurfaceAtlasBuilder {
  public:
    SurfaceAtlasBuilder(const SurfaceMesh &sourceMesh, const MeshingOptions &settings,
                        std::stop_token cancellation)
        : source(sourceMesh), options(settings), stop(cancellation) {}
    SurfaceAtlas build() {
        readSourceTopology();
        readParameters();
        selectProjectionCharts();
        parameterizeCharts();
        buildScalarField();
        initializeFixedAtlas();
        sampleSharedCurves();
        publishCharts();
        projectAtlasVertices();
        return std::move(atlas);
    }

  private:
    static constexpr double projectionCosine = .5;
    const SurfaceMesh &source;
    const MeshingOptions &options;
    std::stop_token stop;
    std::uint32_t count{};
    std::vector<bool> fixed;
    std::vector<Edge> edges;
    std::vector<std::array<std::uint32_t, 3>> faceEdges;
    std::vector<Vec3> normals;
    std::map<std::uint32_t, BoundaryCondition> conditions;
    std::vector<std::vector<std::uint32_t>> adjacency, components;
    std::vector<std::uint32_t> componentOf, physicalComponentOf;
    std::size_t physicalComponentCount{};
    std::vector<Vec3> chartNormals;
    std::vector<bool> retained;
    std::vector<RawChart> rawCharts;
    std::shared_ptr<SurfaceSizingField> sizing;
    SurfaceAtlas atlas;
    std::vector<std::uint32_t> originalAtlas;
    std::vector<std::vector<EdgeSample>> samples;
    std::vector<std::vector<std::uint32_t>> boundaryGraph, chartComponents;
    std::vector<bool> edgeVisited, requiredCurveKnot;
    std::vector<std::array<ChartPoint, 3>> parameterCorners;

    void readParameters() {
        if (!options.parameterization)
            return;
        const auto &snapshot = options.parameterization->sourceMesh();
        if (snapshot.vertices != source.vertices ||
            snapshot.fixedVertices != source.fixedVertices ||
            snapshot.triangles.size() != source.triangles.size() ||
            options.parameterization->triangles().size() != source.triangles.size())
            throw std::invalid_argument(
                "Surface parameterization belongs to a different source snapshot.");
        parameterCorners.resize(source.triangles.size());
        for (std::size_t f = 0; f < source.triangles.size(); ++f) {
            cancelled(stop);
            if (source.triangles[f].vertices != snapshot.triangles[f].vertices ||
                source.triangles[f].faceId != snapshot.triangles[f].faceId)
                throw std::invalid_argument(
                    "Surface parameterization source indices or face IDs differ.");
            const auto domain = options.parameterization->domain(source.triangles[f].faceId);
            if (domain.orientation != 1 && domain.orientation != -1)
                throw std::invalid_argument("Native surface orientation must be +1 or -1.");
            const auto validBound = [](const std::optional<std::array<double, 2>> &bound) {
                return !bound || (std::isfinite((*bound)[0]) && std::isfinite((*bound)[1]) &&
                                  (*bound)[0] < (*bound)[1]);
            };
            const auto validPeriod = [](const std::optional<double> &period) {
                return !period || (std::isfinite(*period) && *period > 0);
            };
            if (!validBound(domain.uBounds) || !validBound(domain.vBounds) ||
                !validPeriod(domain.uPeriod) || !validPeriod(domain.vPeriod))
                throw std::invalid_argument("Invalid native surface parameter domain.");
            const auto &triangle = options.parameterization->triangles()[f];
            for (std::size_t i = 0; i < 3; ++i) {
                const auto uv = triangle.uv[i];
                if (!std::isfinite(uv.u) || !std::isfinite(uv.v))
                    throw std::invalid_argument("Native source corner parameters must be finite.");
                parameterCorners[f][i] = {uv.u, domain.orientation * uv.v};
                if (triangle.edges[i]) {
                    const auto &edge = *triangle.edges[i];
                    if (!edge.curveId || !edge.useId || !std::isfinite(edge.t0) ||
                        !std::isfinite(edge.t1) || edge.t0 == edge.t1)
                        throw std::invalid_argument("Invalid native source curve-use interval.");
                }
            }
            const auto &uv = parameterCorners[f];
            if (!(orient(uv[0], uv[1], uv[2]) > 0))
                throw std::runtime_error(
                    "Native source parameters are inverted or singular on face " +
                    std::to_string(source.triangles[f].faceId) + ".");
        }
        // Physical source components were captured before this computational
        // cut. Different periodic UV uses must not be glued through 3D welding.
        for (std::uint32_t f = 0; f < count; ++f)
            std::erase_if(adjacency[f], [&](std::uint32_t other) {
                for (auto ei : faceEdges[f]) {
                    const auto &edge = edges[ei];
                    if ((edge.faces[0] == other || edge.faces[1] == other) &&
                        !parameterEdgeContinuous(edge))
                        return true;
                }
                return false;
            });
    }

    bool parameterEdgeContinuous(const Edge &edge) const {
        const auto a = edge.faces[0], b = edge.faces[1], i = edge.slots[0], j = edge.slots[1];
        if (b == absent)
            return false;
        const auto same = [](ChartPoint p, ChartPoint q) {
            const double scale =
                std::max({1., std::abs(p.u), std::abs(p.v), std::abs(q.u), std::abs(q.v)});
            return std::hypot(p.u - q.u, p.v - q.v) <=
                   128 * std::numeric_limits<double>::epsilon() * scale;
        };
        return same(parameterCorners[a][i], parameterCorners[b][(j + 1) % 3]) &&
               same(parameterCorners[a][(i + 1) % 3], parameterCorners[b][j]);
    }

    std::uint32_t original(std::uint32_t vi) {
        if (originalAtlas[vi] == absent) {
            if (atlas.vertices.size() >= absent)
                throw std::length_error("Surface atlas vertex count exceeds its index range.");
            originalAtlas[vi] = static_cast<std::uint32_t>(atlas.vertices.size());
            atlas.vertices.push_back(source.vertices[vi]);
        }
        return originalAtlas[vi];
    }

    void readSourceTopology() {
        cancelled(stop);
        if (source.triangles.empty() || source.vertices.empty() ||
            source.triangles.size() > absent / 3)
            throw std::invalid_argument(
                "A surface atlas requires a nonempty representable triangle mesh.");
        if (!(options.targetEdgeLength > 0) || !std::isfinite(options.targetEdgeLength) ||
            !(options.adaptive.minimumEdgeLength > 0) ||
            options.adaptive.minimumEdgeLength > options.targetEdgeLength ||
            !(options.adaptive.geometricTolerance > 0) ||
            !std::isfinite(options.adaptive.geometricTolerance) ||
            !(options.adaptive.growthRate > 0) || !std::isfinite(options.adaptive.growthRate))
            throw std::invalid_argument("Invalid surface atlas size settings.");
        count = static_cast<std::uint32_t>(source.triangles.size());
        fixed.assign(source.vertices.size(), false);
        for (auto v : source.fixedVertices) {
            if (v >= fixed.size() || fixed[v])
                throw std::invalid_argument("Invalid fixed surface vertex.");
            fixed[v] = true;
        }
        std::map<EdgeKey, std::uint32_t> edgeIndex;

        faceEdges.resize(count);
        normals.resize(count);

        for (std::uint32_t f = 0; f < count; ++f) {
            const auto &t = source.triangles[f];
            const auto [entry, added] = conditions.emplace(t.faceId, t.boundary);
            if (!added && entry->second != t.boundary)
                throw std::invalid_argument(
                    "A physical face has inconsistent boundary conditions.");
            for (auto v : t.vertices) {
                if (v >= source.vertices.size())
                    throw std::invalid_argument("Surface atlas has an invalid source index.");
                const auto p = source.vertices[v];
                if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
                    throw std::invalid_argument("Surface atlas source coordinates must be finite.");
            }
            normals[f] = cross(sub(source.vertices[t.vertices[1]], source.vertices[t.vertices[0]]),
                               sub(source.vertices[t.vertices[2]], source.vertices[t.vertices[0]]));
            const double magnitude = length(normals[f]);
            if (!(magnitude > 0))
                throw std::invalid_argument("Surface atlas source has a degenerate triangle.");
            normals[f] = mul(normals[f], 1 / magnitude);
            for (std::uint32_t i = 0; i < 3; ++i) {
                const auto a = t.vertices[i], b = t.vertices[(i + 1) % 3];
                const auto k = key(a, b);
                auto [found, created] =
                    edgeIndex.emplace(k, static_cast<std::uint32_t>(edges.size()));
                if (created)
                    edges.push_back({k});
                auto &e = edges[found->second];
                const auto side = e.faces[0] == absent ? 0u : 1u;
                if (e.faces[side] != absent)
                    throw std::invalid_argument("Surface atlas source has a nonmanifold edge.");
                if (side == 1 && source.triangles[e.faces[0]].vertices[e.slots[0]] == a)
                    throw std::invalid_argument(
                        "Surface atlas source has inconsistent edge orientation.");
                e.faces[side] = f;
                e.slots[side] = i;
                faceEdges[f][i] = found->second;
            }
        }
        for (const auto &[id, h] : options.faceTargetEdgeLengths)
            if (!conditions.contains(id) || !std::isfinite(h) ||
                h < options.adaptive.minimumEdgeLength || h > options.targetEdgeLength)
                throw std::invalid_argument("Invalid local surface atlas target.");
        adjacency.resize(count);
        for (const auto &e : edges)
            if (e.faces[1] != absent &&
                source.triangles[e.faces[0]].faceId == source.triangles[e.faces[1]].faceId) {
                adjacency[e.faces[0]].push_back(e.faces[1]);
                adjacency[e.faces[1]].push_back(e.faces[0]);
            }

        componentOf.assign(count, absent);
        for (std::uint32_t f = 0; f < count; ++f)
            if (componentOf[f] == absent) {
                const auto ci = static_cast<std::uint32_t>(components.size());
                std::vector<std::uint32_t> component{f};
                componentOf[f] = ci;
                for (std::size_t j = 0; j < component.size(); ++j)
                    for (auto next : adjacency[component[j]])
                        if (componentOf[next] == absent) {
                            componentOf[next] = ci;
                            component.push_back(next);
                        }
                components.push_back(std::move(component));
            }
        physicalComponentOf = componentOf;
        physicalComponentCount = components.size();
    }

    void selectProjectionCharts() {
        splitSmoothComponents();
        if (options.parameterization) {
            // Native UVs already provide an injective chart, including separate
            // periodic seam uses. Preview normals must not introduce projection
            // cuts or pinch points into that exact parameter domain.
            chartNormals.clear();
            for (std::uint32_t ci = 0; ci < components.size(); ++ci) {
                chartNormals.push_back(normals[components[ci].front()]);
                for (const auto f : components[ci])
                    componentOf[f] = ci;
            }
            return;
        }
        selectNormalRegions();
        partitionProjectedRegions();
        mergeCompatibleRegions();
        conditionArtificialCorners();
    }

    void splitSmoothComponents() {
        const double coneCosine = std::cos(std::numbers::pi / 4);
        // A strong source crease is geometric, even when both sides carry the same
        // physical face tag. It must not be bridged by a new triangle.
        for (std::uint32_t f = 0; f < count; ++f)
            std::erase_if(adjacency[f], [&](std::uint32_t other) {
                return !options.parameterization && dot(normals[f], normals[other]) < coneCosine;
            });
        std::vector<std::vector<std::uint32_t>> smoothComponents;
        std::vector<bool> smoothVisited(count);
        for (std::uint32_t start = 0; start < count; ++start)
            if (!smoothVisited[start]) {
                std::vector<std::uint32_t> component{start};
                smoothVisited[start] = true;
                for (std::size_t i = 0; i < component.size(); ++i)
                    for (auto next : adjacency[component[i]])
                        if (!smoothVisited[next]) {
                            smoothVisited[next] = true;
                            component.push_back(next);
                        }
                smoothComponents.push_back(std::move(component));
            }
        components = std::move(smoothComponents);
    }

    void selectNormalRegions() {
        const double coneCosine = std::cos(std::numbers::pi / 4);

        // Whole patches may use up to 60 degrees (physical metric condition <=2).
        // This avoids introducing artificial seams in nearly planar CAD patches.

        std::vector<Vec3> directions;
        for (double a : {-1.0, 1.0})
            for (double b : {-1.0, 1.0})
                for (double c : {-1.0, 1.0})
                    directions.push_back({a, b, c});
        for (auto &direction : directions)
            direction = mul(direction, 1 / length(direction));
        for (const auto &component : components) {
            cancelled(stop);
            auto normal = normals[component.front()];
            Vec3 areaNormal{};
            double summedArea = 0;
            for (auto f : component) {
                const auto &t = source.triangles[f];
                const auto n =
                    cross(sub(source.vertices[t.vertices[1]], source.vertices[t.vertices[0]]),
                          sub(source.vertices[t.vertices[2]], source.vertices[t.vertices[0]]));
                areaNormal = add(areaNormal, n);
                summedArea += length(n);
            }
            if (length(areaNormal) > 1e-10 * summedArea)
                normal = mul(areaNormal, 1 / length(areaNormal));
            bool projected = true;
            for (auto f : component)
                projected = projected && dot(normal, normals[f]) >= projectionCosine;
            std::map<std::uint32_t, double> angleSums;
            std::set<std::uint32_t> border;
            for (auto f : component) {
                const auto &t = source.triangles[f];
                for (std::size_t i = 0; i < 3; ++i) {
                    const auto a = sub(source.vertices[t.vertices[(i + 1) % 3]],
                                       source.vertices[t.vertices[i]]);
                    const auto b = sub(source.vertices[t.vertices[(i + 2) % 3]],
                                       source.vertices[t.vertices[i]]);
                    angleSums[t.vertices[i]] += std::atan2(length(cross(a, b)), dot(a, b));
                    const auto &e = edges[faceEdges[f][i]];
                    if (e.faces[1] == absent ||
                        physicalComponentOf[e.faces[0]] != physicalComponentOf[e.faces[1]])
                        border.insert(e.vertices.begin(), e.vertices.end());
                }
            }
            bool developable = true;
            for (const auto &[vi, angle] : angleSums)
                if (!border.contains(vi) && std::abs(angle - 2 * std::numbers::pi) > 1e-8) {
                    developable = false;
                    break;
                }
            if (projected || developable) {
                chartComponents.push_back(component);
                chartNormals.push_back(normal);
                continue;
            }
            // Normal octants are traversal-independent and have a 54.736-degree
            // projection bound. Their meridian cuts avoid oblique raster staircases
            // on axis-aligned horns while keeping physical face identities intact.
            std::map<std::uint32_t, std::uint32_t> labels;
            for (auto f : component) {
                // The octant dot-product maximum is exactly its component signs.
                // A score epsilon creates artificial islands near flat meridians.
                const std::uint32_t label = (normals[f].x > 0 ? 4u : 0u) |
                                            (normals[f].y > 0 ? 2u : 0u) |
                                            (normals[f].z > 0 ? 1u : 0u);
                labels.emplace(f, label);
            }
            std::set<std::uint32_t> visited;
            for (auto start : component)
                if (visited.insert(start).second) {
                    std::vector<std::uint32_t> region{start};
                    const auto label = labels.at(start);
                    for (std::size_t i = 0; i < region.size(); ++i)
                        for (auto next : adjacency[region[i]])
                            if (labels.at(next) == label && visited.insert(next).second)
                                region.push_back(next);
                    if (region.size() > 8)
                        continue;
                    std::map<std::uint32_t, double> alternatives;
                    for (auto f : region)
                        for (auto ei : faceEdges[f]) {
                            const auto &e = edges[ei];
                            if (e.faces[1] == absent ||
                                physicalComponentOf[e.faces[0]] != physicalComponentOf[e.faces[1]])
                                continue;
                            const auto other = e.faces[0] == f ? e.faces[1] : e.faces[0];
                            // UV seams and geometric creases split a physical face
                            // into separate computational components. An adjacent
                            // triangle across that cut has no label in this region.
                            const auto neighbor = labels.find(other);
                            if (neighbor == labels.end())
                                continue;
                            const auto candidate = neighbor->second;
                            if (candidate != label)
                                alternatives[candidate] +=
                                    length(sub(source.vertices[e.vertices[0]],
                                               source.vertices[e.vertices[1]]));
                        }
                    std::uint32_t replacement = label;
                    double shared = 0;
                    for (const auto &[candidate, boundaryLength] : alternatives) {
                        bool valid = true;
                        for (auto f : region)
                            valid = valid && dot(normals[f], directions[candidate]) >= coneCosine;
                        if (valid && boundaryLength > shared) {
                            replacement = candidate;
                            shared = boundaryLength;
                        }
                    }
                    for (auto f : region)
                        labels[f] = replacement;
                }
            visited.clear();
            for (auto start : component)
                if (visited.insert(start).second) {
                    std::vector<std::uint32_t> region{start};
                    const auto label = labels.at(start);
                    for (std::size_t i = 0; i < region.size(); ++i)
                        for (auto next : adjacency[region[i]])
                            if (labels.at(next) == label && visited.insert(next).second)
                                region.push_back(next);
                    chartComponents.push_back(std::move(region));
                    chartNormals.push_back(directions[label]);
                }
        }
    }

    // Validate global injectivity before merging and corner conditioning. A
    // later median split would introduce new artificial corners after their
    // conditioning and can isolate a thin but smooth source facet indefinitely.
    void partitionProjectedRegions() {
        if (options.parameterization)
            return;
        for (std::size_t ci = 0; ci < chartComponents.size(); ++ci) {
            cancelled(stop);
            const auto region = chartComponents[ci];
            const auto normal = chartNormals[ci];
            if (region.size() < 2 ||
                !std::all_of(region.begin(), region.end(),
                             [&](auto f) { return dot(normal, normals[f]) >= projectionCosine; }) ||
                validComponentProjection(region, normal, source, edges, faceEdges, stop))
                continue;
            const auto &t = source.triangles[region.front()];
            const auto origin = source.vertices[t.vertices[0]];
            const auto initial = sub(source.vertices[t.vertices[1]], origin);
            auto u = sub(initial, mul(normal, dot(initial, normal)));
            u = mul(u, 1 / length(u));
            const auto v = cross(normal, u);
            double lowU = std::numeric_limits<double>::infinity(), lowV = lowU, highU = -lowU,
                   highV = -lowU;
            for (auto f : region)
                for (auto vi : source.triangles[f].vertices) {
                    const auto d = sub(source.vertices[vi], origin);
                    const auto pu = dot(d, u), pv = dot(d, v);
                    lowU = std::min(lowU, pu);
                    highU = std::max(highU, pu);
                    lowV = std::min(lowV, pv);
                    highV = std::max(highV, pv);
                }
            const auto axis = highU - lowU >= highV - lowV ? u : v;
            std::vector<std::pair<double, std::uint32_t>> order;
            for (auto f : region) {
                double center = 0;
                for (auto vi : source.triangles[f].vertices)
                    center += dot(sub(source.vertices[vi], origin), axis) / 3;
                order.emplace_back(center, f);
            }
            std::sort(order.begin(), order.end());
            std::map<std::uint32_t, bool> side;
            for (std::size_t j = 0; j < order.size(); ++j)
                side.emplace(order[j].second, j >= order.size() / 2);
            std::set<std::uint32_t> visited;
            for (const auto &[position, start] : order) {
                (void)position;
                if (!visited.insert(start).second)
                    continue;
                std::vector<std::uint32_t> piece{start};
                for (std::size_t j = 0; j < piece.size(); ++j)
                    for (auto next : adjacency[piece[j]]) {
                        const auto found = side.find(next);
                        if (found != side.end() && found->second == side.at(start) &&
                            visited.insert(next).second)
                            piece.push_back(next);
                    }
                chartComponents.push_back(std::move(piece));
                chartNormals.push_back(normal);
            }
            chartComponents[ci].clear();
        }
        std::vector<std::vector<std::uint32_t>> regions;
        std::vector<Vec3> directions;
        for (std::size_t i = 0; i < chartComponents.size(); ++i)
            if (!chartComponents[i].empty()) {
                regions.push_back(std::move(chartComponents[i]));
                directions.push_back(chartNormals[i]);
            }
        chartComponents = std::move(regions);
        chartNormals = std::move(directions);
    }

    void mergeCompatibleRegions() {
        const double coneCosine = std::cos(std::numbers::pi / 4);
        std::vector<std::uint32_t> regionOf(count);
        for (std::uint32_t i = 0; i < chartComponents.size(); ++i)
            for (auto f : chartComponents[i])
                regionOf[f] = i;
        for (;;) {
            cancelled(stop);
            std::map<EdgeKey, double> sharedLengths;
            for (const auto &edge : edges) {
                if (edge.faces[1] == absent ||
                    physicalComponentOf[edge.faces[0]] != physicalComponentOf[edge.faces[1]] ||
                    dot(normals[edge.faces[0]], normals[edge.faces[1]]) < coneCosine)
                    continue;
                const auto a = regionOf[edge.faces[0]], b = regionOf[edge.faces[1]];
                if (a != b)
                    sharedLengths[key(a, b)] += length(
                        sub(source.vertices[edge.vertices[0]], source.vertices[edge.vertices[1]]));
            }
            std::vector<std::pair<EdgeKey, double>> candidates(sharedLengths.begin(),
                                                               sharedLengths.end());
            std::sort(candidates.begin(), candidates.end(), [](const auto &a, const auto &b) {
                if (a.second != b.second)
                    return a.second > b.second;
                return a.first < b.first;
            });
            bool merged = false;
            for (const auto &[pair, sharedLength] : candidates) {
                (void)sharedLength;
                auto combined = chartComponents[pair[0]];
                combined.insert(combined.end(), chartComponents[pair[1]].begin(),
                                chartComponents[pair[1]].end());
                Vec3 areaNormal{};
                double area = 0;
                for (auto f : combined) {
                    const auto &t = source.triangles[f];
                    const auto normal =
                        cross(sub(source.vertices[t.vertices[1]], source.vertices[t.vertices[0]]),
                              sub(source.vertices[t.vertices[2]], source.vertices[t.vertices[0]]));
                    areaNormal = add(areaNormal, normal);
                    area += length(normal);
                }
                if (length(areaNormal) <= 1e-10 * area)
                    continue;
                auto direction = mul(areaNormal, 1 / length(areaNormal));
                bool valid = true;
                for (auto f : combined)
                    valid = valid && dot(direction, normals[f]) >= projectionCosine;
                if (!valid) {
                    const auto fitted = fitProjectionCone(
                        combined, normals, mul(areaNormal, 1 / area), projectionCosine, stop);
                    if (fitted) {
                        direction = *fitted;
                        valid = true;
                    }
                }
                if (!valid ||
                    !validComponentProjection(combined, direction, source, edges, faceEdges, stop))
                    continue;
                chartComponents[pair[0]] = std::move(combined);
                chartComponents[pair[1]].clear();
                chartNormals[pair[0]] = direction;
                for (auto f : chartComponents[pair[0]])
                    regionOf[f] = pair[0];
                merged = true;
                break;
            }
            if (!merged)
                break;
        }
        components.clear();
        std::vector<Vec3> mergedNormals;
        for (std::size_t i = 0; i < chartComponents.size(); ++i)
            if (!chartComponents[i].empty()) {
                components.push_back(std::move(chartComponents[i]));
                mergedNormals.push_back(chartNormals[i]);
            }
        chartNormals = std::move(mergedNormals);
        for (std::uint32_t ci = 0; ci < components.size(); ++ci)
            for (auto f : components[ci])
                componentOf[f] = ci;
    }

    std::uint32_t adjacentFace(std::uint32_t face, std::size_t side) const {
        const auto &edge = edges[faceEdges[face][side]];
        return edge.faces[0] == face ? edge.faces[1] : edge.faces[0];
    }

    bool smoothNeighbors(std::uint32_t face, std::uint32_t other) const {
        // The trimmed adjacency excludes physical boundaries and genuine source
        // creases, including creases whose two sides have the same face ID.
        return other != absent && std::find(adjacency[face].begin(), adjacency[face].end(),
                                            other) != adjacency[face].end();
    }

    bool connectedRegion(std::span<const std::uint32_t> region) const {
        if (region.empty())
            return false;
        std::set<std::uint32_t> remaining(region.begin(), region.end());
        std::vector<std::uint32_t> queue{region.front()};
        remaining.erase(region.front());
        for (std::size_t i = 0; i < queue.size(); ++i)
            for (auto next : adjacency[queue[i]])
                if (remaining.erase(next))
                    queue.push_back(next);
        return remaining.empty();
    }

    std::optional<Vec3> conditionedProjection(std::span<const std::uint32_t> region,
                                              Vec3 direction) const {
        const auto valid = [&](Vec3 candidate) {
            return std::all_of(region.begin(), region.end(),
                               [&](auto face) {
                                   return dot(candidate, normals[face]) >= projectionCosine;
                               }) &&
                   validComponentProjection(region, candidate, source, edges, faceEdges, stop);
        };
        if (valid(direction))
            return direction;
        Vec3 areaNormal{};
        double area = 0;
        for (auto face : region) {
            const auto &t = source.triangles[face];
            const auto normal =
                cross(sub(source.vertices[t.vertices[1]], source.vertices[t.vertices[0]]),
                      sub(source.vertices[t.vertices[2]], source.vertices[t.vertices[0]]));
            areaNormal = add(areaNormal, normal);
            area += length(normal);
        }
        const auto fitted =
            fitProjectionCone(region, normals, mul(areaNormal, 1 / area), projectionCosine, stop);
        if (fitted && valid(*fitted))
            return fitted;
        return std::nullopt;
    }

    struct CornerScore {
        std::size_t rejected{};
        double deficit{};
    };

    std::optional<CornerScore> artificialCornerScore(std::span<const std::uint32_t> region,
                                                     Vec3 direction,
                                                     std::span<const std::uint32_t> labels) const {
        struct Corner {
            std::uint32_t before{absent}, after{absent};
            bool incomingArtificial{}, outgoingArtificial{};
        };
        std::map<std::uint32_t, Corner> corners;
        for (auto face : region) {
            const auto &t = source.triangles[face];
            for (std::size_t i = 0; i < 3; ++i) {
                const auto other = adjacentFace(face, i);
                if (other != absent && labels[other] == labels[face])
                    continue;
                auto &from = corners[t.vertices[i]];
                auto &to = corners[t.vertices[(i + 1) % 3]];
                if (from.after != absent || to.before != absent)
                    return std::nullopt;
                from.after = t.vertices[(i + 1) % 3];
                to.before = t.vertices[i];
                from.outgoingArtificial = to.incomingArtificial = smoothNeighbors(face, other);
            }
        }
        const auto &first = source.triangles[region.front()];
        const auto edge =
            sub(source.vertices[first.vertices[1]], source.vertices[first.vertices[0]]);
        auto u = sub(edge, mul(direction, dot(edge, direction)));
        u = mul(u, 1 / length(u));
        const auto v = cross(direction, u);
        CornerScore score;
        for (const auto &[vertex, corner] : corners) {
            if (corner.before == absent || corner.after == absent)
                return std::nullopt;
            if (!corner.incomingArtificial && !corner.outgoingArtificial)
                continue;
            const auto a = sub(source.vertices[corner.before], source.vertices[vertex]);
            const auto b = sub(source.vertices[corner.after], source.vertices[vertex]);
            if (CGAL::orientation(
                    PredicateKernel::Point_2(dot(a, u), dot(a, v)), PredicateKernel::Point_2(0, 0),
                    PredicateKernel::Point_2(dot(b, u), dot(b, v))) != CGAL::LEFT_TURN)
                continue;
            const double angle =
                std::atan2(length(cross(a, b)), dot(a, b)) * 180 / std::numbers::pi;
            if (angle < options.quality.minimumAngleDegrees) {
                ++score.rejected;
                score.deficit += options.quality.minimumAngleDegrees - angle;
            }
        }
        return score;
    }

    bool transferArtificialFan(std::uint32_t vertex, std::span<const std::uint32_t> star) {
        std::set<std::uint32_t> donors;
        for (auto face : star)
            donors.insert(componentOf[face]);
        if (donors.size() < 2)
            return false;
        for (auto donor : donors) {
            std::vector<std::uint32_t> moved;
            double angle = 0;
            std::set<std::uint32_t> recipients;
            unsigned boundarySides = 0;
            for (auto face : star) {
                if (componentOf[face] != donor)
                    continue;
                moved.push_back(face);
                const auto &t = source.triangles[face];
                std::size_t c = 0;
                while (t.vertices[c] != vertex)
                    ++c;
                const auto a =
                    sub(source.vertices[t.vertices[(c + 2) % 3]], source.vertices[vertex]);
                const auto b =
                    sub(source.vertices[t.vertices[(c + 1) % 3]], source.vertices[vertex]);
                angle += std::atan2(length(cross(a, b)), dot(a, b)) * 180 / std::numbers::pi;
                for (const auto side : {c, (c + 2) % 3}) {
                    const auto other = adjacentFace(face, side);
                    if (other == absent || componentOf[other] != donor) {
                        ++boundarySides;
                        if (smoothNeighbors(face, other))
                            recipients.insert(componentOf[other]);
                    }
                }
            }
            if (moved.empty() || moved.size() == components[donor].size() || boundarySides != 2 ||
                recipients.size() != 1 || angle >= options.quality.minimumAngleDegrees)
                continue;
            const auto recipient = *recipients.begin();
            auto reduced = components[donor], enlarged = components[recipient];
            std::erase_if(reduced, [&](auto face) {
                return std::find(moved.begin(), moved.end(), face) != moved.end();
            });
            enlarged.insert(enlarged.end(), moved.begin(), moved.end());
            if (!connectedRegion(reduced) || !connectedRegion(enlarged))
                continue;
            bool crease = false;
            for (auto face : moved)
                for (std::size_t side = 0; side < 3; ++side) {
                    const auto other = adjacentFace(face, side);
                    crease = crease || (other != absent && componentOf[other] == recipient &&
                                        !smoothNeighbors(face, other));
                }
            if (crease)
                continue;
            const auto oldProjection = conditionedProjection(reduced, chartNormals[donor]);
            const auto newProjection = conditionedProjection(enlarged, chartNormals[recipient]);
            if (!oldProjection || !newProjection)
                continue;
            auto labels = componentOf;
            for (auto face : moved)
                labels[face] = recipient;
            const auto oldFirst =
                artificialCornerScore(components[donor], chartNormals[donor], componentOf);
            const auto oldSecond =
                artificialCornerScore(components[recipient], chartNormals[recipient], componentOf);
            const auto newFirst = artificialCornerScore(reduced, *oldProjection, labels);
            const auto newSecond = artificialCornerScore(enlarged, *newProjection, labels);
            if (!oldFirst || !oldSecond || !newFirst || !newSecond)
                continue;
            const auto oldCount = oldFirst->rejected + oldSecond->rejected;
            const auto newCount = newFirst->rejected + newSecond->rejected;
            const auto oldDeficit = oldFirst->deficit + oldSecond->deficit;
            const auto newDeficit = newFirst->deficit + newSecond->deficit;
            const auto roundoff =
                64 * std::numeric_limits<double>::epsilon() * std::max(1., oldDeficit);
            if (newCount > oldCount ||
                (newCount == oldCount && !(newDeficit < oldDeficit - roundoff)))
                continue;
            components[donor] = std::move(reduced);
            components[recipient] = std::move(enlarged);
            for (auto face : moved)
                componentOf[face] = recipient;
            chartNormals[donor] = *oldProjection;
            chartNormals[recipient] = *newProjection;
            return true;
        }
        return false;
    }

    void conditionArtificialCorners() {
        std::vector<std::vector<std::uint32_t>> stars(source.vertices.size());
        for (std::uint32_t f = 0; f < count; ++f)
            for (auto vi : source.triangles[f].vertices)
                stars[vi].push_back(f);
        // An artificial cut can bound an acute source-facet fan. The angle
        // survives any boundary subdivision, so move the complete fan across
        // that cut while preserving physical boundaries and source creases.
        for (std::uint32_t iteration = 0; iteration < count; ++iteration) {
            cancelled(stop);
            bool changed = false;
            for (std::uint32_t vi = 0; vi < stars.size() && !changed; ++vi)
                changed = transferArtificialFan(vi, stars[vi]);
            if (!changed)
                return;
        }
        throw std::runtime_error("Surface chart corner conditioning exceeded its work budget.");
    }

    std::optional<std::set<std::uint32_t>> annularCut(
        std::span<const std::uint32_t> component) const {
        std::set<std::uint32_t> vertices, componentEdges;
        std::map<std::uint32_t, std::vector<std::uint32_t>> boundary;
        const auto ci = componentOf[component.front()];
        for (const auto face : component)
            for (const auto edgeId : faceEdges[face]) {
                if (!componentEdges.insert(edgeId).second)
                    continue;
                const auto &edge = edges[edgeId];
                vertices.insert(edge.vertices.begin(), edge.vertices.end());
                if (edge.faces[1] == absent || componentOf[edge.faces[0]] != ci ||
                    componentOf[edge.faces[1]] != ci) {
                    boundary[edge.vertices[0]].push_back(edge.vertices[1]);
                    boundary[edge.vertices[1]].push_back(edge.vertices[0]);
                }
            }
        // An oriented connected annulus has Euler characteristic zero and two
        // simple boundary cycles. Other surface topologies use the disk cutter.
        if (vertices.size() + component.size() != componentEdges.size() || boundary.empty())
            return std::nullopt;
        std::map<std::uint32_t, unsigned> boundaryLoop;
        unsigned loopCount = 0;
        for (const auto &[start, neighbors] : boundary) {
            if (neighbors.size() != 2)
                return std::nullopt;
            if (boundaryLoop.contains(start))
                continue;
            std::vector<std::uint32_t> pending{start};
            boundaryLoop.emplace(start, loopCount);
            for (std::size_t i = 0; i < pending.size(); ++i)
                for (const auto next : boundary.at(pending[i]))
                    if (boundaryLoop.emplace(next, loopCount).second)
                        pending.push_back(next);
            ++loopCount;
        }
        if (loopCount != 2)
            return std::nullopt;
        std::map<std::uint32_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>> graph;
        for (const auto edgeId : componentEdges) {
            const auto &edge = edges[edgeId];
            graph[edge.vertices[0]].emplace_back(edge.vertices[1], edgeId);
            graph[edge.vertices[1]].emplace_back(edge.vertices[0], edgeId);
        }
        using Candidate = std::pair<double, std::uint32_t>;
        std::priority_queue<Candidate, std::vector<Candidate>, std::greater<Candidate>> pending;
        std::map<std::uint32_t, double> distance;
        std::map<std::uint32_t, std::pair<std::uint32_t, std::uint32_t>> predecessor;
        for (const auto vertex : vertices)
            distance[vertex] = std::numeric_limits<double>::infinity();
        for (const auto &[vertex, loop] : boundaryLoop)
            if (loop == 0) {
                distance[vertex] = 0;
                pending.emplace(0, vertex);
            }
        while (!pending.empty()) {
            cancelled(stop);
            const auto [travel, vertex] = pending.top();
            pending.pop();
            if (travel != distance.at(vertex))
                continue;
            if (const auto found = boundaryLoop.find(vertex);
                found != boundaryLoop.end() && found->second == 1) {
                std::set<std::uint32_t> cut;
                auto current = vertex;
                while (predecessor.contains(current)) {
                    const auto [previous, edgeId] = predecessor.at(current);
                    cut.insert(edgeId);
                    current = previous;
                }
                return cut;
            }
            for (const auto &[next, edgeId] : graph.at(vertex)) {
                const double candidate =
                    travel + length(sub(source.vertices[next], source.vertices[vertex]));
                if (candidate < distance.at(next)) {
                    distance[next] = candidate;
                    predecessor[next] = {vertex, edgeId};
                    pending.emplace(candidate, next);
                }
            }
        }
        throw std::runtime_error("Annular surface boundaries have no connecting source path.");
    }

    void parameterizeCharts() {
        DisjointSet corners(3 * count);
        retained.assign(count, false);

        for (std::size_t ci = 0; ci < components.size(); ++ci) {
            const auto component = components[ci];
            cancelled(stop);
            RawChart raw;
            const auto &first = source.triangles[component.front()];
            raw.chart.faceId = first.faceId;
            raw.chart.boundary = first.boundary;
            raw.chart.geometry = options.geometry;
            raw.chart.parameterization = options.parameterization;
            raw.chart.parameterAccuracy = options.adaptive.geometricTolerance;
            if (options.parameterization)
                raw.chart.parameterOrientation =
                    options.parameterization->domain(first.faceId).orientation;
            const auto origin = source.vertices[first.vertices[0]];
            auto direction = chartNormals[ci];
            double extent = 0, deflection = 0;
            for (auto f : component)
                for (auto vi : source.triangles[f].vertices) {
                    const auto delta = sub(source.vertices[vi], origin);
                    extent = std::max(extent, length(delta));
                    deflection =
                        std::max(deflection, std::abs(dot(delta, normals[component.front()])));
                }
            bool planar = deflection <= 512 * std::numeric_limits<double>::epsilon() * extent;
            if (planar)
                for (auto f : component)
                    planar = planar && dot(normals[f], normals[component.front()]) > 0;
            if (planar)
                direction = normals[component.front()];
            const auto initialEdge = sub(source.vertices[first.vertices[1]], origin);
            const auto projectedEdge =
                sub(initialEdge, mul(direction, dot(initialEdge, direction)));
            const auto u = mul(projectedEdge, 1 / length(projectedEdge));
            const auto v = cross(direction, u);
            raw.chart.planar = planar;
            raw.chart.planeOrigin = origin;
            raw.chart.planeU = u;
            raw.chart.planeV = v;
            bool projectable = true;
            for (auto f : component)
                projectable = projectable && dot(direction, normals[f]) >= projectionCosine;
            raw.projected = !options.parameterization && (planar || projectable);
            auto glue = [&](Edge &e) {
                const auto a = e.faces[0], b = e.faces[1], i = e.slots[0], j = e.slots[1];
                if (options.parameterization && !parameterEdgeContinuous(e))
                    return;
                corners.join(3 * a + i, 3 * b + (j + 1) % 3);
                corners.join(3 * a + (i + 1) % 3, 3 * b + j);
                e.glued = true;
            };
            bool fullyFixed = true;
            for (auto f : component)
                for (auto vi : source.triangles[f].vertices)
                    fullyFixed = fullyFixed && fixed[vi];
            if (fullyFixed) {
                for (auto f : component) {
                    retained[f] = true;
                    for (auto ei : faceEdges[f]) {
                        auto &e = edges[ei];
                        if (!e.glued && e.faces[1] != absent &&
                            componentOf[e.faces[0]] == componentOf[e.faces[1]])
                            glue(e);
                    }
                }
                continue;
            }
            if (options.parameterization || planar || projectable) {
                for (auto f : component)
                    for (auto ei : faceEdges[f]) {
                        auto &e = edges[ei];
                        if (!e.glued && e.faces[1] != absent &&
                            componentOf[e.faces[0]] == componentOf[e.faces[1]])
                            glue(e);
                    }
            } else if (const auto cut = annularCut(component)) {
                // Open an annulus along its shortest boundary-to-boundary path.
                // A face-tree cut can instead wind obliquely around a thin ring
                // and create acute artificial corners that subdivision cannot
                // improve. On an extruded strip this path is a generator, so
                // the unfolded ends meet the unchanged physical rims normally.
                for (const auto face : component)
                    for (const auto edgeId : faceEdges[face]) {
                        auto &edge = edges[edgeId];
                        if (!edge.glued && !cut->contains(edgeId) && edge.faces[1] != absent &&
                            componentOf[edge.faces[0]] == componentOf[edge.faces[1]])
                            glue(edge);
                    }
            } else {
                std::set<std::uint32_t> visited{component.front()};
                std::vector<std::uint32_t> queue{component.front()};
                for (std::size_t j = 0; j < queue.size(); ++j)
                    for (auto ei : faceEdges[queue[j]]) {
                        auto &e = edges[ei];
                        if (e.faces[1] == absent ||
                            componentOf[e.faces[0]] != componentOf[e.faces[1]])
                            continue;
                        const auto next = e.faces[0] == queue[j] ? e.faces[1] : e.faces[0];
                        if (visited.insert(next).second) {
                            glue(e);
                            queue.push_back(next);
                        }
                    }
                std::size_t borderCount = 3 * component.size() - 2 * (component.size() - 1);
                std::deque<std::uint32_t> candidates;
                for (auto f : component)
                    for (auto ei : faceEdges[f])
                        if (!edges[ei].glued)
                            candidates.push_back(ei);
                while (!candidates.empty()) {
                    const auto ei = candidates.front();
                    candidates.pop_front();
                    auto &e = edges[ei];
                    if (e.glued || e.faces[1] == absent ||
                        componentOf[e.faces[0]] != componentOf[e.faces[1]] || borderCount <= 4)
                        continue;
                    const auto a = e.faces[0], b = e.faces[1], i = e.slots[0], j = e.slots[1];
                    const bool firstSame =
                        corners.find(3 * a + i) == corners.find(3 * b + (j + 1) % 3);
                    const bool secondSame =
                        corners.find(3 * a + (i + 1) % 3) == corners.find(3 * b + j);
                    // Gluing boundary edges with exactly one identified endpoint
                    // removes a leaf of the cut graph while retaining a disk.
                    if (firstSame == secondSame)
                        continue;
                    auto changed = corners.join(3 * a + i, 3 * b + (j + 1) % 3);
                    auto changed2 = corners.join(3 * a + (i + 1) % 3, 3 * b + j);
                    changed.insert(changed.end(), changed2.begin(), changed2.end());
                    e.glued = true;
                    borderCount -= 2;
                    for (auto corner : changed)
                        for (auto incident : faceEdges[corner / 3])
                            candidates.push_back(incident);
                }
            }
            std::map<std::uint32_t, std::uint32_t> local;
            auto localIndex = [&](std::uint32_t corner) {
                const auto root = corners.find(corner);
                auto [it, added] = local.emplace(root, static_cast<std::uint32_t>(local.size()));
                if (added) {
                    raw.sourceVertices.push_back(source.triangles[corner / 3].vertices[corner % 3]);
                    if (options.parameterization)
                        raw.uv.push_back(parameterCorners[corner / 3][corner % 3]);
                }
                return it->second;
            };
            std::vector<RawHalfedge> border;
            for (auto f : component) {
                const std::array<std::uint32_t, 3> localTriangle{
                    localIndex(3 * f), localIndex(3 * f + 1), localIndex(3 * f + 2)};
                raw.triangles.push_back(localTriangle);
                for (std::uint32_t i = 0; i < 3; ++i)
                    if (!edges[faceEdges[f][i]].glued)
                        border.push_back(
                            {localTriangle[i], localTriangle[(i + 1) % 3], faceEdges[f][i], f, i});
            }
            std::map<std::uint32_t, std::size_t> outgoing;
            for (std::size_t i = 0; i < border.size(); ++i)
                if (!outgoing.emplace(border[i].from, i).second)
                    throw std::runtime_error("Surface chart cut has a branching boundary.");
            std::vector<bool> traced(border.size());
            for (std::size_t i = 0; i < border.size(); ++i)
                if (!traced[i]) {
                    std::vector<RawHalfedge> loop;
                    auto next = i;
                    do {
                        if (traced[next])
                            throw std::runtime_error(
                                "Surface chart cut has a broken boundary cycle.");
                        traced[next] = true;
                        loop.push_back(border[next]);
                        const auto found = outgoing.find(border[next].to);
                        if (found == outgoing.end())
                            throw std::runtime_error("Surface chart cut has an open boundary.");
                        next = found->second;
                    } while (next != i);
                    raw.loops.push_back(std::move(loop));
                }
            if (options.parameterization) {
                if (!simpleProjectedBoundary(raw))
                    throw std::runtime_error(
                        "Native surface chart has overlapping parameter boundaries on face " +
                        std::to_string(first.faceId) + ".");
            } else if (planar || projectable) {
                for (auto vi : raw.sourceVertices) {
                    const auto delta = sub(source.vertices[vi], origin);
                    raw.uv.push_back({dot(delta, u), dot(delta, v)});
                }
                if (!simpleProjectedBoundary(raw))
                    throw std::runtime_error(
                        "Conditioned projected surface chart has overlapping boundaries.");
            } else if (!unfoldDevelopable(raw, source, stop))
                embedDisk(raw, source, stop);
            rawCharts.push_back(std::move(raw));
        }
    }

    void buildScalarField() {
        std::vector<double> vertexSizes(source.vertices.size(), options.targetEdgeLength);
        std::map<EdgeKey, double> featureSizes;
        std::vector<double> areas(physicalComponentCount), perimeters(physicalComponentCount);
        for (std::uint32_t f = 0; f < count; ++f) {
            const auto &t = source.triangles[f];
            areas[physicalComponentOf[f]] +=
                .5 *
                length(cross(sub(source.vertices[t.vertices[1]], source.vertices[t.vertices[0]]),
                             sub(source.vertices[t.vertices[2]], source.vertices[t.vertices[0]])));
        }
        for (const auto &e : edges)
            if (e.faces[1] == absent ||
                physicalComponentOf[e.faces[0]] != physicalComponentOf[e.faces[1]]) {
                const double d =
                    length(sub(source.vertices[e.vertices[1]], source.vertices[e.vertices[0]]));
                for (auto f : e.faces)
                    if (f != absent)
                        perimeters[physicalComponentOf[f]] += d;
            }
        auto feature = [&](EdgeKey edge, double size) {
            size = std::clamp(size, options.adaptive.minimumEdgeLength, options.targetEdgeLength);
            auto [it, created] = featureSizes.emplace(edge, size);
            if (!created)
                it->second = std::min(it->second, size);
            for (auto vi : edge)
                vertexSizes[vi] = std::min(vertexSizes[vi], size);
        };
        auto chordSize = [&](double radius) {
            const double e = std::min(radius, options.adaptive.geometricTolerance);
            return 2 * std::sqrt(e * (2 * radius - e));
        };
        using SeamKey = std::array<std::uint32_t, 2>;
        std::map<SeamKey, std::map<std::uint32_t, std::vector<std::uint32_t>>> physicalSeams;
        for (const auto &e : edges) {
            const auto a = source.vertices[e.vertices[0]], b = source.vertices[e.vertices[1]];
            const bool physical = e.faces[1] == absent || physicalComponentOf[e.faces[0]] !=
                                                              physicalComponentOf[e.faces[1]];
            // Mandatory physical edge endpoints cannot be coarsened away.
            // Grade their actual spacing into the surface, independently of
            // optional shape limits. Free preview edges and artificial chart
            // cuts do not become sizing features.
            if (fixed[e.vertices[0]] && fixed[e.vertices[1]] &&
                (physical || !smoothNeighbors(e.faces[0], e.faces[1])))
                feature(e.vertices, length(sub(b, a)));
            if (physical) {
                double h = options.targetEdgeLength;
                for (auto f : e.faces)
                    if (f != absent && !retained[f] && perimeters[physicalComponentOf[f]] > 0)
                        h = std::min(h, 4 * areas[physicalComponentOf[f]] /
                                            perimeters[physicalComponentOf[f]]);
                feature(e.vertices, h);
                const auto left = source.triangles[e.faces[0]].faceId;
                const auto right =
                    e.faces[1] == absent ? absent : source.triangles[e.faces[1]].faceId;
                auto &seam = physicalSeams[key(left, right)];
                seam[e.vertices[0]].push_back(e.vertices[1]);
                seam[e.vertices[1]].push_back(e.vertices[0]);
            } else if (options.adaptive.enabled) {
                const auto n1 = normals[e.faces[0]], n2 = normals[e.faces[1]];
                const double sine = length(cross(n1, n2));
                if (sine <= 1e-7)
                    continue;
                const double angle = std::atan2(sine, std::clamp(dot(n1, n2), -1.0, 1.0));
                const auto c =
                    source.vertices[source.triangles[e.faces[0]].vertices[(e.slots[0] + 2) % 3]];
                const auto d =
                    source.vertices[source.triangles[e.faces[1]].vertices[(e.slots[1] + 2) % 3]];
                const double dual =
                    (length(cross(sub(c, a), sub(b, a))) + length(cross(sub(d, a), sub(b, a)))) /
                    (2 * length(sub(b, a)));
                const double radius = dual / (2 * std::sin(angle / 2));
                if (radius > 0 && std::isfinite(radius))
                    feature(e.vertices, chordSize(radius));
            }
        }
        if (options.adaptive.enabled)
            for (const auto &[seamId, seam] : physicalSeams) {
                (void)seamId;
                std::set<std::uint32_t> visited;
                for (const auto &[start, neighbors] : seam) {
                    (void)neighbors;
                    if (visited.contains(start))
                        continue;
                    std::vector<std::uint32_t> connected{start};
                    visited.insert(start);
                    bool closed = true;
                    for (std::size_t i = 0; i < connected.size(); ++i) {
                        const auto &next = seam.at(connected[i]);
                        closed = closed && next.size() == 2;
                        for (auto vi : next)
                            if (visited.insert(vi).second)
                                connected.push_back(vi);
                    }
                    for (auto vi : connected) {
                        const auto &next = seam.at(vi);
                        if (next.size() != 2)
                            continue;
                        const auto a = sub(source.vertices[next[0]], source.vertices[vi]),
                                   b = sub(source.vertices[next[1]], source.vertices[vi]);
                        const double denominator = 2 * length(cross(a, b));
                        if (denominator <= 1e-7 * length(a) * length(b))
                            continue;
                        const double radius =
                            length(a) * length(b) * length(sub(a, b)) / denominator;
                        const double h =
                            closed ? std::min(chordSize(radius),
                                              2 * radius * std::sin(std::numbers::pi / 12))
                                   : chordSize(radius);
                        feature(key(vi, next[0]), h);
                        feature(key(vi, next[1]), h);
                    }
                }
            }
        std::vector<SurfaceEdgeSeed> seeds;
        for (const auto &[edge, h] : featureSizes)
            seeds.push_back({edge, h});
        sizing = std::make_shared<SurfaceSizingField>(
            source, std::move(vertexSizes), options.adaptive.growthRate, options.targetEdgeLength,
            seeds, stop, options.faceTargetEdgeLengths);
    }

    void initializeFixedAtlas() {
        originalAtlas.assign(source.vertices.size(), absent);
        for (std::uint32_t vi = 0; vi < fixed.size(); ++vi)
            if (fixed[vi])
                atlas.fixedVertices.push_back(original(vi));
        for (std::uint32_t f = 0; f < count; ++f)
            if (retained[f]) {
                auto t = source.triangles[f];
                for (auto &vi : t.vertices)
                    vi = original(vi);
                atlas.fixedTriangles.push_back(t);
            }
    }

    void initializeSharedCurves() {

        samples.resize(edges.size());
        boundaryGraph.resize(source.vertices.size());
        for (std::uint32_t ei = 0; ei < edges.size(); ++ei)
            if (!edges[ei].glued)
                for (auto vi : edges[ei].vertices)
                    boundaryGraph[vi].push_back(ei);
        edgeVisited.assign(edges.size(), false);
        requiredCurveKnot.assign(source.vertices.size(), false);
        for (const auto &raw : rawCharts)
            if (!options.parameterization && !raw.chart.planar)
                for (const auto &loop : raw.loops) {
                    std::vector<int> turns(loop.size());
                    for (std::size_t i = 0; i < loop.size(); ++i) {
                        const auto a = raw.uv[loop[(i + loop.size() - 1) % loop.size()].from];
                        const auto b = raw.uv[loop[i].from], c = raw.uv[loop[i].to];
                        const double scale =
                            std::hypot(b.u - a.u, b.v - a.v) * std::hypot(c.u - b.u, c.v - b.v);
                        const double value = orient(a, b, c);
                        if (std::abs(value) > 256 * std::numeric_limits<double>::epsilon() * scale)
                            turns[i] = value > 0 ? 1 : -1;
                        if (!raw.projected && turns[i] < 0)
                            requiredCurveKnot[raw.sourceVertices[loop[i].from]] = true;
                    }
                    if (raw.projected) {
                        int previous = 0;
                        for (auto it = turns.rbegin(); it != turns.rend(); ++it)
                            if (*it != 0) {
                                previous = *it;
                                break;
                            }
                        for (std::size_t i = 0; i < turns.size(); ++i)
                            if (turns[i] != 0) {
                                if (previous != 0 && turns[i] != previous)
                                    requiredCurveKnot[raw.sourceVertices[loop[i].from]] = true;
                                previous = turns[i];
                            }
                    }
                }
    }

    bool mandatory(std::uint32_t vi) const {
        if (fixed[vi] || requiredCurveKnot[vi] || boundaryGraph[vi].size() != 2)
            return true;
        if (options.parameterization) {
            std::optional<std::uint32_t> firstCurve;
            for (auto ei : boundaryGraph[vi]) {
                const auto &edge = edges[ei];
                const auto &curve =
                    options.parameterization->triangles()[edge.faces[0]].edges[edge.slots[0]];
                const auto id = curve ? curve->curveId : 0;
                if (firstCurve && *firstCurve != id)
                    return true;
                firstCurve = id;
            }
        }
        return false;
    }

    Vec3 parameterEdgePosition(std::uint32_t edgeId, double fraction) const {
        const auto &edge = edges[edgeId];
        const auto f = edge.faces[0], slot = edge.slots[0];
        if (source.triangles[f].vertices[slot] != edge.vertices[0])
            fraction = 1 - fraction;
        const auto &triangle = options.parameterization->triangles()[f];
        if (triangle.edges[slot]) {
            const auto &curve = *triangle.edges[slot];
            return options.parameterization
                ->curve(curve.curveId, curve.t0 + fraction * (curve.t1 - curve.t0))
                .position;
        }
        const auto a = triangle.uv[slot], b = triangle.uv[(slot + 1) % 3];
        return options.parameterization
            ->surface(source.triangles[f].faceId,
                      {a.u + fraction * (b.u - a.u), a.v + fraction * (b.v - a.v)})
            .position;
    }

    std::uint32_t addVertexSample(std::uint32_t vi) {
        const auto id = original(vi);
        for (auto ei : boundaryGraph[vi])
            samples[ei].push_back({edges[ei].vertices[0] == vi ? 0.0 : 1.0, id});
        return id;
    }

    void sampleChain(std::uint32_t start, std::uint32_t firstEdge) {
        std::vector<std::uint32_t> chainVertices{start}, chainEdges;
        auto current = start, ei = firstEdge;
        for (;;) {
            if (edgeVisited[ei])
                throw std::runtime_error("Surface atlas boundary chain overlaps itself.");
            edgeVisited[ei] = true;
            chainEdges.push_back(ei);
            const auto &e = edges[ei];
            const auto next = e.vertices[0] == current ? e.vertices[1] : e.vertices[0];
            chainVertices.push_back(next);
            if (next == start || mandatory(next))
                break;
            const auto &choices = boundaryGraph[next];
            ei = choices[0] == ei ? choices[1] : choices[0];
            current = next;
        }
        const bool closed = chainVertices.front() == chainVertices.back();
        std::vector<double> arc(chainVertices.size());
        for (std::size_t i = 1; i < arc.size(); ++i)
            arc[i] = arc[i - 1] + length(sub(source.vertices[chainVertices[i]],
                                             source.vertices[chainVertices[i - 1]]));
        const double total = arc.back();
        if (!(total > 0))
            throw std::runtime_error("Surface atlas has a zero-length shared curve.");
        auto locate = [&](double s) {
            auto it = std::upper_bound(arc.begin(), arc.end(), s);
            return std::min(static_cast<std::size_t>(it - arc.begin() - 1), chainEdges.size() - 1);
        };
        auto position = [&](double s) {
            const auto i = locate(s);
            const double t = (s - arc[i]) / (arc[i + 1] - arc[i]);
            return add(
                source.vertices[chainVertices[i]],
                mul(sub(source.vertices[chainVertices[i + 1]], source.vertices[chainVertices[i]]),
                    t));
        };
        auto sizeAt = [&](double s) {
            const auto &e = edges[chainEdges[locate(s)]];
            std::array<std::uint32_t, 2> ids{source.triangles[e.faces[0]].faceId, 0};
            std::size_t n = 1;
            if (e.faces[1] != absent && source.triangles[e.faces[1]].faceId != ids[0])
                ids[n++] = source.triangles[e.faces[1]].faceId;
            return sizing->at(position(s), std::span<const std::uint32_t>(ids.data(), n));
        };
        std::vector<double> knots{0, total};
        // A skipped source corner is permitted only within the geometric chord
        // budget. Artificial cuts use the same budget, never curvature seeds.
        std::vector<std::pair<double, double>> intervals;
        intervals.emplace_back(0, total);
        while (!intervals.empty()) {
            cancelled(stop);
            const auto [left, right] = intervals.back();
            intervals.pop_back();
            const auto a = position(left), b = position(right), ab = sub(b, a);
            const double denominator = dot(ab, ab);
            double largest = 0, worst = 0;
            for (auto it = std::upper_bound(arc.begin(), arc.end(), left);
                 it != arc.end() && *it < right; ++it) {
                const auto q = sub(
                    source.vertices[chainVertices[static_cast<std::size_t>(it - arc.begin())]], a);
                const double t =
                    denominator > 0 ? std::clamp(dot(q, ab) / denominator, 0.0, 1.0) : 0;
                const double deviation = length(sub(q, mul(ab, t)));
                if (deviation > largest) {
                    largest = deviation;
                    worst = *it;
                }
            }
            if (largest > options.adaptive.geometricTolerance / 2) {
                knots.push_back(worst);
                intervals.emplace_back(left, worst);
                intervals.emplace_back(worst, right);
            }
        }
        std::sort(knots.begin(), knots.end());
        knots.erase(std::unique(knots.begin(), knots.end()), knots.end());
        if (closed && knots.size() < 4) {
            knots = {0, total / 3, 2 * total / 3, total};
        }
        const auto geometricKnots = knots;
        knots.clear();
        knots.push_back(0);
        for (std::size_t interval = 1; interval < geometricKnots.size(); ++interval) {
            const double end = geometricKnots[interval];
            const auto beginIndex = knots.size() - 1;
            for (double s = knots.back(); s < end;) {
                cancelled(stop);
                const double h = .8 * std::min(sizeAt(s), sizeAt(std::min(end, s + sizeAt(s))));
                const double next = std::min(end, s + h);
                if (!(next > s))
                    throw std::runtime_error("Surface atlas boundary resampling made no progress.");
                s = next;
                knots.push_back(s);
            }
            if (knots.size() > beginIndex + 2) {
                const auto n = knots.size();
                const double last = knots[n - 1] - knots[n - 2],
                             previous = knots[n - 2] - knots[n - 3];
                if (last < .5 * previous)
                    knots.erase(knots.end() - 2);
            }
        }
        for (double s : knots) {
            auto near = std::lower_bound(arc.begin(), arc.end(), s);
            if (near != arc.end() &&
                std::abs(*near - s) <= 64 * std::numeric_limits<double>::epsilon() * total) {
                addVertexSample(chainVertices[static_cast<std::size_t>(near - arc.begin())]);
                continue;
            }
            if (near != arc.begin() && std::abs(*std::prev(near) - s) <=
                                           64 * std::numeric_limits<double>::epsilon() * total) {
                addVertexSample(
                    chainVertices[static_cast<std::size_t>(std::prev(near) - arc.begin())]);
                continue;
            }
            const auto i = locate(s);
            const auto edgeId = chainEdges[i];
            double t = (s - arc[i]) / (arc[i + 1] - arc[i]);
            if (edges[edgeId].vertices[0] != chainVertices[i])
                t = 1 - t;
            if (atlas.vertices.size() >= absent)
                throw std::length_error("Surface atlas vertex count exceeds its index range.");
            const auto id = static_cast<std::uint32_t>(atlas.vertices.size());
            const auto p =
                options.parameterization ? parameterEdgePosition(edgeId, t) : position(s);
            if (!finitePosition(p))
                throw std::runtime_error("Native boundary evaluation returned a non-finite point.");
            atlas.vertices.push_back(p);
            samples[edgeId].push_back({t, id});
        }
    }

    void sortSamples() {
        for (auto &list : samples) {
            std::sort(list.begin(), list.end(),
                      [](const auto &a, const auto &b) { return a.t < b.t; });
            list.erase(std::unique(list.begin(), list.end(),
                                   [](const auto &a, const auto &b) {
                                       return a.atlasVertex == b.atlasVertex;
                                   }),
                       list.end());
        }
    }

    void sampleBoundaries() {
        std::fill(edgeVisited.begin(), edgeVisited.end(), false);
        for (auto &list : samples)
            list.clear();
        for (std::uint32_t ei = 0; ei < edges.size(); ++ei) {
            const auto &e = edges[ei];
            if (!e.glued &&
                (retained[e.faces[0]] || (e.faces[1] != absent && retained[e.faces[1]]))) {
                edgeVisited[ei] = true;
                addVertexSample(e.vertices[0]);
                addVertexSample(e.vertices[1]);
            }
        }
        // Required source knots divide chains before distributing samples, so
        // a later geometric guard cannot create an arbitrarily short segment.
        for (std::uint32_t vi = 0; vi < boundaryGraph.size(); ++vi)
            if (mandatory(vi))
                for (auto ei : boundaryGraph[vi])
                    if (!edgeVisited[ei])
                        sampleChain(vi, ei);
        for (std::uint32_t ei = 0; ei < edges.size(); ++ei)
            if (!edges[ei].glued && !edgeVisited[ei])
                sampleChain(edges[ei].vertices[0], ei);
        sortSamples();
    }

    void sampleSharedCurves() {
        initializeSharedCurves();

        for (;;) {
            sampleBoundaries();
            bool changed = false;
            for (const auto &raw : rawCharts)
                if (raw.projected && !raw.chart.planar)
                    for (const auto &loop : raw.loops)
                        for (const auto &piece :
                             boundaryPieces(raw, loop, source, edges, samples, atlas.vertices))
                            if (referenceLens(raw, piece, options.adaptive.geometricTolerance) ==
                                LensKind::retainKnots)
                                for (const auto &point : piece)
                                    if (point.sourceVertex != absent &&
                                        !requiredCurveKnot[point.sourceVertex]) {
                                        requiredCurveKnot[point.sourceVertex] = true;
                                        changed = true;
                                    }
            if (!changed)
                break;
        }
    }

    void publishCharts() {
        auto curves = std::make_shared<AtlasCurveRegistry>();
        curves->geometry = options.geometry;
        curves->parameterization = options.parameterization;
        curves->geometricTolerance = options.adaptive.geometricTolerance;
        for (const auto &triangle : atlas.fixedTriangles)
            for (std::size_t i = 0; i < 3; ++i)
                curves->retained.insert(key(triangle.vertices[i], triangle.vertices[(i + 1) % 3]));
        for (auto &raw : rawCharts) {
            cancelled(stop);
            auto &chart = raw.chart;
            chart.sizing = sizing;
            chart.geometry = options.geometry;
            std::set<std::uint32_t> boundaryOriginal;
            std::vector<ChartReferenceTriangle> extensions;
            for (const auto &rawLoop : raw.loops) {
                SurfaceChartLoop loop;
                for (const auto &halfedge : rawLoop) {
                    boundaryOriginal.insert(raw.sourceVertices[halfedge.from]);
                    boundaryOriginal.insert(raw.sourceVertices[halfedge.to]);
                }
                for (const auto &piece :
                     boundaryPieces(raw, rawLoop, source, edges, samples, atlas.vertices)) {
                    loop.vertices.push_back({piece.front().uv, piece.front().atlasVertex});
                    if (referenceLens(raw, piece, options.adaptive.geometricTolerance,
                                      &extensions) == LensKind::retainKnots)
                        throw std::runtime_error(
                            "A shared curve cannot be covered by a bounded reference lens.");
                    const auto a = piece.front().atlasVertex, b = piece.back().atlasVertex;
                    if (a == b)
                        throw std::runtime_error(
                            "A shared curve segment closes onto the same atlas node.");
                    std::vector<Vec3> positions;
                    std::vector<ChartPoint> uv;
                    for (const auto &point : piece) {
                        positions.push_back(point.position);
                        uv.push_back(point.uv);
                    }
                    std::optional<AtlasCurveRegistry::ParameterPath> nativePath;
                    if (options.parameterization)
                        nativePath = parameterPathForPiece(chart, piece);
                    const auto segmentKey = key(a, b);
                    auto found = curves->segments.find(segmentKey);
                    if (found == curves->segments.end())
                        found = curves->segments
                                    .emplace(segmentKey,
                                             AtlasCurveRegistry::Segment{
                                                 a, b, positions, {}, nativePath})
                                    .first;
                    else {
                        if (found->second.a != a) {
                            std::reverse(positions.begin(), positions.end());
                            std::reverse(uv.begin(), uv.end());
                            if (nativePath) {
                                std::reverse(nativePath->parameters.begin(),
                                             nativePath->parameters.end());
                                std::reverse(nativePath->uv.begin(), nativePath->uv.end());
                            }
                        }
                        if (nativePath) {
                            const auto &existing = found->second.parameterPath;
                            if (!existing || existing->curveId != nativePath->curveId ||
                                existing->parameters.size() != nativePath->parameters.size())
                                throw std::runtime_error("Native shared-curve sides disagree on "
                                                         "their canonical binding.");
                            if (nativePath->curveId)
                                for (std::size_t j = 0; j < nativePath->parameters.size(); ++j) {
                                    const double t = nativePath->parameters[j],
                                                 other = existing->parameters[j];
                                    if (std::abs(t - other) >
                                        256 * std::numeric_limits<double>::epsilon() *
                                            std::max({1., std::abs(t), std::abs(other)}))
                                        throw std::runtime_error("Native shared-curve sides have "
                                                                 "different spatial parameters.");
                                }
                        }
                        if (positions.size() != found->second.positions.size())
                            throw std::runtime_error(
                                "Shared curve sides disagree on their original source polyline.");
                        for (std::size_t j = 0; j < positions.size(); ++j)
                            if (length(sub(positions[j], found->second.positions[j])) >
                                1e-12 * options.targetEdgeLength)
                                throw std::runtime_error("Shared curve sides disagree on their "
                                                         "original source positions.");
                    }
                    if (options.parameterization)
                        found->second.sides.push_back(parameterSideForPiece(
                            chart, piece, atlas.charts.size(), chart.loops.size(), std::move(uv),
                            options.adaptive.geometricTolerance));
                    else
                        found->second.sides.push_back(
                            {atlas.charts.size(), chart.loops.size(), std::move(uv)});
                }
                chart.loops.push_back(std::move(loop));
            }
            for (std::size_t i = 0; i < raw.sourceVertices.size(); ++i)
                if (fixed[raw.sourceVertices[i]] &&
                    !boundaryOriginal.contains(raw.sourceVertices[i]))
                    chart.fixedPoints.push_back({raw.uv[i], original(raw.sourceVertices[i])});
            const std::array<std::uint32_t, 1> patch{chart.faceId};
            for (const auto &t : raw.triangles) {
                ChartReferenceTriangle reference;
                for (std::size_t i = 0; i < 3; ++i) {
                    reference.uv[i] = raw.uv[t[i]];
                    reference.positions[i] = source.vertices[raw.sourceVertices[t[i]]];
                    reference.sizes[i] = sizing->at(reference.positions[i], patch);
                }
                chart.reference.push_back(reference);
            }
            for (auto &extension : extensions) {
                for (std::size_t i = 0; i < 3; ++i)
                    extension.sizes[i] = sizing->at(extension.positions[i], patch);
                chart.reference.push_back(std::move(extension));
            }
            auto lookup = std::make_shared<ChartLookup>();
            double upperU = -std::numeric_limits<double>::infinity(), upperV = upperU;
            lookup->lowerU = lookup->lowerV = std::numeric_limits<double>::infinity();
            for (const auto &p : raw.uv) {
                lookup->lowerU = std::min(lookup->lowerU, p.u);
                lookup->lowerV = std::min(lookup->lowerV, p.v);
                upperU = std::max(upperU, p.u);
                upperV = std::max(upperV, p.v);
            }
            const double width = upperU - lookup->lowerU, height = upperV - lookup->lowerV;
            if (!(width > 0 && height > 0))
                throw std::runtime_error("Surface chart has a collapsed bounding box.");
            const double ideal = std::sqrt(static_cast<double>(chart.reference.size()) / 4);
            lookup->countU = static_cast<std::uint32_t>(
                std::clamp(ideal * std::sqrt(width / height), 1.0, 256.0));
            lookup->countV = static_cast<std::uint32_t>(
                std::clamp(ideal * std::sqrt(height / width), 1.0, 256.0));
            lookup->cellU = width / lookup->countU;
            lookup->cellV = height / lookup->countV;
            lookup->cells.resize(static_cast<std::size_t>(lookup->countU) * lookup->countV);
            for (std::uint32_t i = 0; i < chart.reference.size(); ++i) {
                const auto &uv = chart.reference[i].uv;
                const ChartPoint lo{std::min({uv[0].u, uv[1].u, uv[2].u}),
                                    std::min({uv[0].v, uv[1].v, uv[2].v})};
                const ChartPoint hi{std::max({uv[0].u, uv[1].u, uv[2].u}),
                                    std::max({uv[0].v, uv[1].v, uv[2].v})};
                const auto low = lookup->cell(lo), high = lookup->cell(hi);
                for (std::size_t y = low / lookup->countU; y <= high / lookup->countU; ++y)
                    for (std::size_t x = low % lookup->countU; x <= high % lookup->countU; ++x)
                        lookup->cells[y * lookup->countU + x].push_back(i);
            }
            chart.lookup = std::move(lookup);
            atlas.charts.push_back(std::move(chart));
        }
        atlas.curves = std::move(curves);
    }

    void projectAtlasVertices() {
        if (options.parameterization || !options.geometry)
            return;
        std::vector<std::set<std::uint32_t>> incident(atlas.vertices.size());
        for (const auto &chart : atlas.charts) {
            for (const auto &loop : chart.loops)
                for (const auto &vertex : loop.vertices)
                    incident[vertex.atlasVertex].insert(chart.faceId);
            for (const auto &vertex : chart.fixedPoints)
                incident[vertex.atlasVertex].insert(chart.faceId);
        }
        for (const auto &triangle : atlas.fixedTriangles)
            for (auto vertex : triangle.vertices)
                incident[vertex].insert(triangle.faceId);
        const std::set<std::uint32_t> pinned(atlas.fixedVertices.begin(),
                                             atlas.fixedVertices.end());
        for (std::uint32_t i = 0; i < atlas.vertices.size(); ++i) {
            cancelled(stop);
            // Boundary resampling can retire earlier samples during lens repair.
            if (incident[i].empty())
                continue;
            const std::vector<std::uint32_t> faces(incident[i].begin(), incident[i].end());
            const auto projected = exactSharedPosition(*options.geometry, faces, atlas.vertices[i]);
            if (pinned.contains(i)) {
                if (length(sub(projected, atlas.vertices[i])) >
                    1e-10 * std::max(1.0, length(atlas.vertices[i])))
                    throw std::runtime_error("A fixed source vertex is not on its exact geometry.");
            } else
                atlas.vertices[i] = projected;
        }
    }
};
} // namespace

SurfaceAtlas makeSurfaceAtlas(const SurfaceMesh &source, const MeshingOptions &options,
                              std::stop_token stop) {
    return SurfaceAtlasBuilder(source, options, stop).build();
}

namespace {
struct CornerLineEvidence {
    const BoundaryCorner &corner;
    double roundoff;

    bool same(Vec3 a, Vec3 b) const {
        return length(sub(a, b)) <= roundoff;
    }
    bool onInterval(Vec3 p, Vec3 end) const {
        const auto direction = sub(end, corner.position), delta = sub(p, corner.position);
        const double extent = length(direction);
        const auto unit = mul(direction, 1 / extent);
        const double along = dot(delta, unit);
        return along >= -roundoff && along <= extent + roundoff &&
               length(sub(delta, mul(unit, along))) <= roundoff;
    }
};
CornerLineEvidence lineEvidence(const BoundaryCorner &corner) {
    if (!finitePosition(corner.position))
        throw std::runtime_error("Geometrieanbieter meldet einen nicht endlichen Eckpunkt.");
    const auto &constraint = *corner.constraint;
    double scale = length(corner.position);
    if (!std::isfinite(scale))
        throw std::runtime_error("Geometrieanbieter-Eckpunkt ueberschreitet den Zahlenbereich.");
    for (const auto end : constraint.boundaryEnds) {
        const double extent = length(sub(end, corner.position));
        if (!finitePosition(end) || !std::isfinite(length(end)) || !std::isfinite(extent) ||
            !(extent > 0))
            throw std::runtime_error("Geometrieanbieter meldet ein ungueltiges Randintervall.");
        scale = std::max({scale, length(end), length(sub(end, corner.position))});
    }
    const auto normal = constraint.loopNormal;
    if (!finitePosition(normal) ||
        std::abs(length(normal) - 1) > 128 * std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Geometrieanbieter meldet eine ungueltige Ecknormalenrichtung.");
    // Only arithmetic error in point/line comparisons is covered. Neither CAD
    // modelling tolerance nor the requested mesh deviation enters this proof.
    CornerLineEvidence evidence{corner, 128 * std::numeric_limits<double>::epsilon() * scale};
    for (const auto end : constraint.boundaryEnds)
        if (std::abs(dot(sub(end, corner.position), normal)) > evidence.roundoff)
            throw std::runtime_error("Geometrieanbieter meldet keine koplanaren Randintervalle.");
    return evidence;
}
bool retainedLine(const SurfaceAtlas &atlas, std::size_t chartIndex, std::size_t loopIndex,
                  std::size_t edgeIndex, Vec3 end, const CornerLineEvidence &evidence) {
    if (!atlas.curves)
        return false;
    const auto &loop = atlas.charts[chartIndex].loops[loopIndex].vertices;
    const auto &a = loop[edgeIndex], &b = loop[(edgeIndex + 1) % loop.size()];
    const auto found = atlas.curves->segments.find(key(a.atlasVertex, b.atlasVertex));
    if (found == atlas.curves->segments.end())
        return false;
    const auto &segment = found->second;
    const auto sameUv = [](ChartPoint x, ChartPoint y) { return x.u == y.u && x.v == y.v; };
    bool bound = false;
    for (const auto &side : segment.sides) {
        if (side.chart != chartIndex || side.loop != loopIndex ||
            side.uv.size() != segment.positions.size() || side.uv.empty())
            continue;
        bound |= (segment.a == a.atlasVertex && segment.b == b.atlasVertex &&
                  sameUv(side.uv.front(), a.uv) && sameUv(side.uv.back(), b.uv)) ||
                 (segment.b == a.atlasVertex && segment.a == b.atlasVertex &&
                  sameUv(side.uv.back(), a.uv) && sameUv(side.uv.front(), b.uv));
    }
    return bound && segment.positions.size() >= 2 &&
           evidence.onInterval(atlas.vertices.at(a.atlasVertex), end) &&
           evidence.onInterval(atlas.vertices.at(b.atlasVertex), end) &&
           std::all_of(segment.positions.begin(), segment.positions.end(),
                       [&](Vec3 p) { return evidence.onInterval(p, end); });
}
std::optional<double> retainedCornerAngle(const SurfaceAtlas &atlas,
                                          const CornerLineEvidence &evidence) {
    const auto &corner = evidence.corner;
    const auto &constraint = *corner.constraint;
    std::optional<double> result;
    for (std::size_t ci = 0; ci < atlas.charts.size(); ++ci) {
        const auto &chart = atlas.charts[ci];
        if (chart.faceId != corner.faceId || !chart.planar)
            continue;
        for (std::size_t li = 0; li < chart.loops.size(); ++li) {
            const auto &loop = chart.loops[li].vertices;
            for (std::size_t i = 0; i < loop.size(); ++i) {
                const auto p = atlas.vertices.at(loop[i].atlasVertex);
                if (!evidence.same(p, corner.position))
                    continue;
                const auto previous = (i + loop.size() - 1) % loop.size();
                const bool forward =
                    retainedLine(atlas, ci, li, previous, constraint.boundaryEnds[0], evidence) &&
                    retainedLine(atlas, ci, li, i, constraint.boundaryEnds[1], evidence);
                const bool reverse =
                    retainedLine(atlas, ci, li, previous, constraint.boundaryEnds[1], evidence) &&
                    retainedLine(atlas, ci, li, i, constraint.boundaryEnds[0], evidence);
                if (!forward && !reverse)
                    continue;
                // A computational cut through this junction does not establish
                // the complete physical corner in any one current chart.
                if (result)
                    return std::nullopt;
                const auto a = sub(p, atlas.vertices.at(loop[previous].atlasVertex));
                const auto b = sub(atlas.vertices.at(loop[(i + 1) % loop.size()].atlasVertex), p);
                if (!(length(a) > evidence.roundoff && length(b) > evidence.roundoff))
                    continue;
                const auto product = cross(a, b);
                const double turn =
                    std::copysign(std::atan2(length(product), dot(a, b)),
                                  (forward ? 1 : -1) * dot(constraint.loopNormal, product));
                result = std::numbers::pi - turn;
            }
        }
    }
    return result;
}
void checkCornerQuality(std::uint32_t face, double angle, const MeshQualityLimits &quality) {
    const double angleLimit = quality.minimumAngleDegrees * (std::numbers::pi / 180.);
    // A planar fan has a triangle angle <= its retained convex corner. For
    // theta<=pi/3, aspect=longest/shortestAltitude >= 1/sin(theta). This is a
    // necessary condition only; it does not certify a feasible triangulation.
    const double aspectLimit = std::asin(1 / quality.maximumAspectRatio);
    const double required = std::max(angleLimit, aspectLimit);
    const double roundoff =
        8 * std::numeric_limits<double>::epsilon() * (std::abs(angle) + required);
    if (!(angle < required && required - angle > roundoff))
        return;
    std::ostringstream message;
    message << std::setprecision(7) << "Flaeche " << face << ": erhaltener ebener Eckwinkel "
            << angle * (180. / std::numbers::pi) << " Grad ist mit ";
    if (angleLimit >= aspectLimit)
        message << "Mindestwinkel " << quality.minimumAngleDegrees << " Grad";
    else
        message << "Seitenverhaeltnis " << quality.maximumAspectRatio;
    message << " unvereinbar. Qualitaetsgrenze anpassen oder Geometrie aendern.";
    throw std::invalid_argument(message.str());
}
} // namespace
void validateAtlasBoundaryCorners(const SurfaceAtlas &atlas, const MeshingOptions &options,
                                  std::stop_token stop) {
    if (!options.quality.enabled || !options.geometry)
        return;
    std::set<std::uint32_t> active;
    for (const auto &chart : atlas.charts)
        active.insert(chart.faceId);
    for (const auto &triangle : atlas.fixedTriangles)
        active.insert(triangle.faceId);
    for (const auto &corner : options.geometry->boundaryCorners()) {
        cancelled(stop);
        if (!active.contains(corner.faceId))
            continue;
        if (!finitePosition(corner.position))
            throw std::runtime_error("Geometrieanbieter meldet einen nicht endlichen Eckpunkt.");
        if (!corner.constraint)
            continue;
        const auto evidence = lineEvidence(corner);
        if (const auto angle = retainedCornerAngle(atlas, evidence))
            checkCornerQuality(corner.faceId, *angle, options.quality);
    }
}

std::vector<std::size_t> refineSurfaceAtlas(SurfaceAtlas &atlas,
                                            std::span<const std::array<std::uint32_t, 2>> segments,
                                            std::stop_token stop) {
    cancelled(stop);
    if (segments.empty())
        return {};
    if (!atlas.curves)
        throw std::invalid_argument("The atlas has no original shared-curve registry.");
    std::set<EdgeKey> requested;
    for (const auto &segment : segments) {
        const auto k = key(segment[0], segment[1]);
        if (atlas.curves->retained.contains(k))
            throw std::runtime_error(
                "Shared curve " + std::to_string(k[0]) + "-" + std::to_string(k[1]) +
                " belongs to a retained fixed triangle; conforming refinement is infeasible.");
        if (k[0] == k[1] || !atlas.curves->segments.contains(k))
            throw std::invalid_argument(
                "Requested shared curve segment is unknown or already retired.");
        requested.insert(k);
    }
    // Immutable lookup and sizing data are shared. Chart vectors and the curve
    // registry are copied and committed together after all requests.
    auto next = atlas;
    auto registry = std::make_shared<AtlasCurveRegistry>(*atlas.curves);
    std::set<std::size_t> changedCharts;
    auto same = [](ChartPoint a, ChartPoint b) { return a.u == b.u && a.v == b.v; };
    for (const auto &segmentKey : requested) {
        cancelled(stop);
        const auto parent = registry->segments.at(segmentKey);
        std::vector<double> arc(parent.positions.size());
        for (std::size_t i = 1; i < arc.size(); ++i)
            arc[i] = arc[i - 1] + length(sub(parent.positions[i], parent.positions[i - 1]));
        if (!(arc.back() > 0))
            throw std::runtime_error("The original shared curve has zero length.");
        const double middle = arc.back() / 2;
        const auto upper = std::upper_bound(arc.begin(), arc.end(), middle);
        auto interval = std::min(static_cast<std::size_t>(upper - arc.begin() - 1), arc.size() - 2);
        double fraction = (middle - arc[interval]) / (arc[interval + 1] - arc[interval]);
        auto position =
            add(parent.positions[interval],
                mul(sub(parent.positions[interval + 1], parent.positions[interval]), fraction));
        std::optional<ParameterMidpoint> native;
        if (parent.parameterPath) {
            native = parameterMidpoint(*registry, parent, stop);
            interval = native->interval;
            fraction = native->fraction;
            position = native->position;
        }
        if (next.vertices.size() >= absent)
            throw std::length_error("Surface atlas vertex count exceeds its index range.");
        const auto id = static_cast<std::uint32_t>(next.vertices.size());
        auto outputPosition = position;
        if (!registry->parameterization && registry->geometry) {
            std::set<std::uint32_t> uniqueFaces;
            for (const auto &side : parent.sides)
                uniqueFaces.insert(next.charts.at(side.chart).faceId);
            const std::vector<std::uint32_t> faces(uniqueFaces.begin(), uniqueFaces.end());
            outputPosition = exactSharedPosition(*registry->geometry, faces, position);
        }
        // A projected midpoint must advance along the physical curve. Keeping
        // a second UV node at an endpoint creates degenerate cells and can
        // trigger an endless front cascade; the copied atlas stays uncommitted.
        if (outputPosition == next.vertices.at(parent.a) ||
            outputPosition == next.vertices.at(parent.b))
            throw std::runtime_error("Shared curve " + std::to_string(parent.a) + "-" +
                                     std::to_string(parent.b) +
                                     " midpoint projection collapsed onto a physical endpoint; "
                                     "no atlas change was committed.");
        next.vertices.push_back(outputPosition);
        AtlasCurveRegistry::Segment left{parent.a, id, {}, {}}, right{id, parent.b, {}, {}};
        std::vector<std::size_t> leftIndices, rightIndices;
        for (std::size_t i = 0; i <= interval; ++i)
            leftIndices.push_back(i);
        for (std::size_t i = interval + 1; i < parent.positions.size(); ++i)
            rightIndices.push_back(i);
        const bool leftHasMiddle = length(sub(parent.positions[interval], position)) == 0;
        const bool rightHasMiddle = length(sub(parent.positions[interval + 1], position)) == 0;
        for (auto i : leftIndices)
            left.positions.push_back(parent.positions[i]);
        if (!leftHasMiddle)
            left.positions.push_back(position);
        right.positions.push_back(position);
        for (auto i : rightIndices)
            if (!(rightHasMiddle && i == interval + 1))
                right.positions.push_back(parent.positions[i]);
        if (left.positions.size() < 2 || right.positions.size() < 2)
            throw std::runtime_error("Shared-curve subdivision reached its representable limit.");
        if (native) {
            auto a = *parent.parameterPath, b = a;
            a.parameters.clear();
            a.uv.clear();
            b.parameters.clear();
            b.uv.clear();
            for (auto i : leftIndices) {
                a.parameters.push_back(parent.parameterPath->parameters[i]);
                a.uv.push_back(parent.parameterPath->uv[i]);
            }
            if (!leftHasMiddle) {
                a.parameters.push_back(native->parameter);
                a.uv.push_back(native->uv);
            }
            b.parameters.push_back(native->parameter);
            b.uv.push_back(native->uv);
            for (auto i : rightIndices)
                if (!(rightHasMiddle && i == interval + 1)) {
                    b.parameters.push_back(parent.parameterPath->parameters[i]);
                    b.uv.push_back(parent.parameterPath->uv[i]);
                }
            left.parameterPath = std::move(a);
            right.parameterPath = std::move(b);
        }
        for (const auto &side : parent.sides) {
            cancelled(stop);
            if (side.uv.size() != parent.positions.size() || side.chart >= next.charts.size() ||
                side.loop >= next.charts[side.chart].loops.size())
                throw std::runtime_error(
                    "The shared curve registry has inconsistent chart references.");
            const auto a = side.uv[interval], b = side.uv[interval + 1];
            ChartPoint uv{a.u + fraction * (b.u - a.u), a.v + fraction * (b.v - a.v)};
            if (native && side.curveUse) {
                const auto use =
                    registry->parameterization->curveUse(side.curveUse, native->parameter);
                const auto &chart = next.charts[side.chart];
                uv = {use.uv.u + side.parameterLift.u,
                      chart.parameterOrientation * use.uv.v + side.parameterLift.v};
                const auto point =
                    registry->parameterization->surface(chart.faceId, use.uv).position;
                if (!finitePosition(point) ||
                    length(sub(point, outputPosition)) > registry->geometricTolerance)
                    throw std::runtime_error(
                        "Native split surface/curve representation exceeds geometric accuracy.");
            }
            auto &loop = next.charts[side.chart].loops[side.loop].vertices;
            std::optional<std::size_t> insertion;
            for (std::size_t i = 0; i < loop.size(); ++i) {
                const auto &first = loop[i], &second = loop[(i + 1) % loop.size()];
                const bool forward =
                    first.atlasVertex == parent.a && second.atlasVertex == parent.b &&
                    same(first.uv, side.uv.front()) && same(second.uv, side.uv.back());
                const bool reverse =
                    first.atlasVertex == parent.b && second.atlasVertex == parent.a &&
                    same(first.uv, side.uv.back()) && same(second.uv, side.uv.front());
                if (forward || reverse) {
                    insertion = i + 1;
                    break;
                }
            }
            if (!insertion)
                throw std::runtime_error("A shared curve side is no longer a chart boundary.");
            // The last-to-first edge is represented by appending before wrap.
            loop.insert(loop.begin() + static_cast<std::ptrdiff_t>(*insertion), {uv, id});
            changedCharts.insert(side.chart);
            AtlasCurveRegistry::Side leftSide{
                side.chart, side.loop, {}, side.curveUse, side.parameterLift},
                rightSide{side.chart, side.loop, {}, side.curveUse, side.parameterLift};
            for (auto i : leftIndices)
                leftSide.uv.push_back(side.uv[i]);
            if (!leftHasMiddle)
                leftSide.uv.push_back(uv);
            rightSide.uv.push_back(uv);
            for (auto i : rightIndices)
                if (!(rightHasMiddle && i == interval + 1))
                    rightSide.uv.push_back(side.uv[i]);
            left.sides.push_back(std::move(leftSide));
            right.sides.push_back(std::move(rightSide));
        }
        registry->segments.erase(segmentKey);
        if (!registry->segments.emplace(key(left.a, left.b), std::move(left)).second ||
            !registry->segments.emplace(key(right.a, right.b), std::move(right)).second)
            throw std::runtime_error("Shared-curve subdivision produced a duplicate segment.");
    }
    cancelled(stop);
    next.curves = std::move(registry);
    atlas = std::move(next);
    return {changedCharts.begin(), changedCharts.end()};
}
} // namespace hornsculpt::detail
