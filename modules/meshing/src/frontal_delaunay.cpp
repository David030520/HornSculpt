#include "frontal_delaunay.hpp"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <algorithm>
#include <cmath>
#include <deque>
#include <iomanip>
#include <list>
#include <map>
#include <numbers>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace hornsculpt::detail {
namespace {
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
struct Quality {
    double angle{180}, aspect{}, priority{}, shortest{std::numeric_limits<double>::max()};
    double sizeViolationEnergy{}, geometryViolationEnergy{};
    bool accepted{true};
};
struct VertexData {
    FrontalVertex value;
    std::uint64_t insertionId{};
    std::uint64_t positionRevision{};
    std::uint32_t output{std::numeric_limits<std::uint32_t>::max()};
};
struct FaceData {
    int nesting{-1};
    std::uint64_t serial{};
    bool blocked{};
    // A face handle may be recycled by the CDT. Identity includes its current
    // corners and their positions, including temporary line-search positions.
    std::array<std::pair<std::uint64_t, std::uint64_t>, 3> qualityIdentity{};
    std::optional<Quality> cachedQuality;
    bool inside() const {
        return nesting >= 0 && nesting % 2 == 1;
    }
};
using VertexBase = CGAL::Triangulation_vertex_base_with_info_2<VertexData, Kernel>;
using ConstraintFaceBase = CGAL::Constrained_triangulation_face_base_2<Kernel>;
using FaceBase = CGAL::Triangulation_face_base_with_info_2<FaceData, Kernel, ConstraintFaceBase>;
using DataStructure = CGAL::Triangulation_data_structure_2<VertexBase, FaceBase>;
using ConstraintTriangulation =
    CGAL::Constrained_triangulation_2<Kernel, DataStructure, CGAL::Exact_predicates_tag>;
using Triangulation =
    CGAL::Constrained_Delaunay_triangulation_2<Kernel, DataStructure, CGAL::Exact_predicates_tag>;
using Vertex = Triangulation::Vertex_handle;
using Face = Triangulation::Face_handle;
using Edge = Triangulation::Edge;
// A non-convergent chart must fail before it can consume unbounded memory.
// This is a work budget, not a relaxed geometric or triangle-quality criterion.
constexpr std::size_t maximumInteriorInsertionsPerChart = 1000000;
constexpr std::size_t maximumCandidateEvaluationsPerChart = 5000000;
constexpr std::size_t maximumQualityVertexTrialsPerChart = 100000;
Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double norm(Vec3 a) {
    return std::hypot(a.x, a.y, a.z);
}
bool finite(Vec3 a) {
    return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z);
}
ChartPoint midpoint(ChartPoint a, ChartPoint b) {
    return {std::midpoint(a.u, b.u), std::midpoint(a.v, b.v)};
}
void checkStop(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Frontal-Delaunay meshing cancelled.");
}
struct Metric {
    double a{}, b{}, c{};
    static Metric from(const SurfaceChartMetric &sample) {
        const double uu = dot(sample.du, sample.du), uv = dot(sample.du, sample.dv);
        const double normalSquared = dot(cross(sample.du, sample.dv), cross(sample.du, sample.dv));
        if (!finite(sample.du) || !finite(sample.dv) || !std::isfinite(uu) || uu <= 0 ||
            !std::isfinite(normalSquared) || normalSquared <= 0)
            throw std::invalid_argument("Singular surface chart metric.");
        // Gram-Schmidt via cross products avoids cancellation in vv-uv*uv/uu.
        Metric result{std::sqrt(uu), uv / std::sqrt(uu), std::sqrt(normalSquared / uu)};
        if (!std::isfinite(result.b) || !std::isfinite(result.c) || result.c <= 0)
            throw std::invalid_argument("Non-finite surface chart metric.");
        return result;
    }
    ChartPoint forward(ChartPoint uv) const {
        return {a * uv.u + b * uv.v, c * uv.v};
    }
    ChartPoint inverse(ChartPoint xy) const {
        const double v = xy.v / c;
        return {(xy.u - b * v) / a, v};
    }
};
struct BoundarySegment {
    ChartPoint middle, a, b;
    Metric metric;
    std::array<std::uint32_t, 2> atlas;
    ChartPoint start, end;
    bool encroached(ChartPoint uv) const {
        const auto p = metric.forward({uv.u - middle.u, uv.v - middle.v});
        const double au = p.u - a.u, av = p.v - a.v, bu = p.u - b.u, bv = p.v - b.v;
        const double product = au * bu + av * bv;
        const double l2 = (a.u - b.u) * (a.u - b.u) + (a.v - b.v) * (a.v - b.v);
        const double roundoff = 64 * std::numeric_limits<double>::epsilon() *
                                (std::abs(au * bu) + std::abs(av * bv) + l2);
        return product < -roundoff;
    }
};
struct ParameterFrame {
    ChartPoint origin;
    Metric metric;
    Kernel::Point_2 point(ChartPoint uv) const {
        const auto p = metric.forward({uv.u - origin.u, uv.v - origin.v});
        return {p.u, p.v};
    }
    ChartPoint uv(const Kernel::Point_2 &p) const {
        const auto delta = metric.inverse({p.x(), p.y()});
        return {origin.u + delta.u, origin.v + delta.v};
    }
};
Metric representativeMetric(const SurfaceChart &chart, std::stop_token stop) {
    struct FacetMetric {
        double area;
        Vec3 du, dv;
    };
    std::vector<FacetMetric> facets;
    facets.reserve(chart.reference.size());
    double area = 0, uu = 0, uv = 0;
    for (const auto &triangle : chart.reference) {
        checkStop(stop);
        const auto a = subtract(triangle.positions[1], triangle.positions[0]);
        const auto b = subtract(triangle.positions[2], triangle.positions[0]);
        const double weight = norm(cross(a, b));
        if (!std::isfinite(weight) || weight <= 0)
            throw std::invalid_argument("Invalid reference facet in the chart metric.");
        Vec3 du, dv;
        if (chart.parameterization) {
            // Native D1 defines the metric. Reconstructing it from a thin
            // preview triangle amplifies coordinate roundoff near poles/seams
            // and can make a regular chart's global CDT frame almost singular.
            ChartPoint center{};
            for (const auto uv : triangle.uv) {
                center.u += uv.u / 3;
                center.v += uv.v / 3;
            }
            const auto exact = chart.parameterization->surface(
                chart.faceId, {center.u, chart.parameterOrientation * center.v});
            du = exact.du;
            dv = {chart.parameterOrientation * exact.dv.x, chart.parameterOrientation * exact.dv.y,
                  chart.parameterOrientation * exact.dv.z};
        } else {
            const double u1 = triangle.uv[1].u - triangle.uv[0].u;
            const double v1 = triangle.uv[1].v - triangle.uv[0].v;
            const double u2 = triangle.uv[2].u - triangle.uv[0].u;
            const double v2 = triangle.uv[2].v - triangle.uv[0].v;
            const double determinant = u1 * v2 - v1 * u2;
            if (!std::isfinite(determinant) || determinant == 0)
                throw std::invalid_argument("Invalid reference facet in the chart metric.");
            du = {(a.x * v2 - b.x * v1) / determinant, (a.y * v2 - b.y * v1) / determinant,
                  (a.z * v2 - b.z * v1) / determinant};
            dv = {(b.x * u1 - a.x * u2) / determinant, (b.y * u1 - a.y * u2) / determinant,
                  (b.z * u1 - a.z * u2) / determinant};
        }
        facets.push_back({weight, du, dv});
        area += weight;
        uu += weight * dot(du, du);
        uv += weight * dot(du, dv);
    }
    // Average the Gram tensor by physical area rather than selecting one
    // arbitrary source facet. This remains covariant under affine UV changes
    // and avoids amplifying metric variation across a curved chart.
    const double projection = uv / uu;
    double residual = 0;
    for (const auto &facet : facets) {
        checkStop(stop);
        const Vec3 orthogonal{facet.dv.x - projection * facet.du.x,
                              facet.dv.y - projection * facet.du.y,
                              facet.dv.z - projection * facet.du.z};
        residual += facet.area * dot(orthogonal, orthogonal);
    }
    // Compute the Schur complement as a sum of squares, avoiding cancellation
    // in vv - uv*uv/uu for strongly scaled or skewed chart coordinates.
    const double scale = std::sqrt(uu / area);
    const Metric result{scale, uv / area / scale, std::sqrt(residual / area)};
    if (!std::isfinite(result.a) || result.a <= 0 || !std::isfinite(result.b) ||
        !std::isfinite(result.c) || result.c <= 0)
        throw std::invalid_argument("Singular representative chart metric.");
    return result;
}
SurfaceChartSample sample(const SurfaceChart &chart, ChartPoint uv) {
    auto s = chart.sample(uv);
    if (!finite(s.position) || !std::isfinite(s.targetEdgeLength) || s.targetEdgeLength <= 0)
        throw std::invalid_argument("Invalid scalar size or position on a surface chart.");
    return s;
}
void floodDomain(Triangulation &tri, Face initial, int nesting, std::list<Edge> &borders) {
    if (initial->info().nesting != -1)
        return;
    std::deque<Face> todo{initial};
    while (!todo.empty()) {
        const auto f = todo.front();
        todo.pop_front();
        if (f->info().nesting != -1)
            continue;
        f->info().nesting = nesting;
        for (int e = 0; e < 3; ++e) {
            if (f->neighbor(e)->info().nesting != -1)
                continue;
            if (tri.is_constrained({f, e}))
                borders.push_back({f, e});
            else
                todo.push_back(f->neighbor(e));
        }
    }
}
void markDomain(Triangulation &tri) {
    for (auto f = tri.all_faces_begin(); f != tri.all_faces_end(); ++f)
        f->info().nesting = -1;
    std::list<Edge> borders;
    floodDomain(tri, tri.infinite_face(), 0, borders);
    while (!borders.empty()) {
        const auto edge = borders.front();
        borders.pop_front();
        floodDomain(tri, edge.first->neighbor(edge.second), edge.first->info().nesting + 1,
                    borders);
    }
}
double normalDeviation(const SurfaceChart &chart, std::span<const Vertex> vertices);
double sampledProjectionDeviation(const SurfaceChart &chart, std::span<const Vertex> vertices,
                                  double upperBoundLimit = 0);
Quality quality(const SurfaceChart &chart, const std::array<Vertex, 3> &vertices,
                const MeshingOptions &options) {
    Quality q;
    std::array<Vec3, 3> p;
    for (int i = 0; i < 3; ++i)
        p[i] = vertices[i]->info().value.position;
    const auto n = cross(subtract(p[1], p[0]), subtract(p[2], p[0]));
    const double area2 = norm(n);
    if (!std::isfinite(area2))
        throw std::runtime_error("Non-finite physical chart triangle area.");
    double longest = 0;
    for (int i = 0; i < 3; ++i) {
        const auto edge = subtract(p[(i + 1) % 3], p[i]), other = subtract(p[(i + 2) % 3], p[i]);
        const double length = norm(edge);
        q.shortest = std::min(q.shortest, length);
        longest = std::max(longest, length);
        q.angle = std::min(
            q.angle, area2 > 0 ? std::atan2(area2, dot(edge, other)) * 180 / std::numbers::pi : 0.);
        const auto target = chart.targetSize(
            midpoint(vertices[i]->info().value.uv, vertices[(i + 1) % 3]->info().value.uv));
        const double ratio = length / target;
        q.priority = std::max(q.priority, ratio);
        q.accepted = q.accepted && ratio <= 4. / 3.;
        const double excess = std::max(0., ratio / (4. / 3.) - 1);
        q.sizeViolationEnergy += excess * excess;
        if (options.adaptive.enabled) {
            const std::array<Vertex, 2> segment{vertices[i], vertices[(i + 1) % 3]};
            const double deviation =
                normalDeviation(chart, segment) / options.adaptive.geometricTolerance;
            q.priority = std::max(q.priority, deviation);
            q.accepted = q.accepted && deviation <= 1;
            q.geometryViolationEnergy += std::pow(std::max(0., deviation - 1), 2);
        }
    }
    if (options.adaptive.enabled) {
        const double deviation =
            normalDeviation(chart, vertices) / options.adaptive.geometricTolerance;
        q.priority = std::max(q.priority, deviation);
        q.accepted = q.accepted && deviation <= 1;
        q.geometryViolationEnergy += std::pow(std::max(0., deviation - 1), 2);
    }
    // A coarse initial CDT can connect two UV copies of the same 3D seam
    // vertex. Its chordal area is then zero despite a valid local PL metric.
    // Such cells still require front refinement; only the final mesh is gated.
    q.aspect = area2 > 0 ? longest * longest / area2 : std::numeric_limits<double>::infinity();
    // Disabling shape limits must not allow a degenerate physical boundary.
    q.accepted = q.accepted && area2 > 0;
    if (!options.quality.enabled) {
        if (area2 <= 0)
            q.priority = std::max(q.priority, 1e12);
    } else {
        q.priority = std::max(
            {q.priority, q.angle > 0 ? options.quality.minimumAngleDegrees / q.angle : 1e12,
             std::isfinite(q.aspect) ? q.aspect / options.quality.maximumAspectRatio : 1e12});
        q.accepted = q.accepted && q.angle >= options.quality.minimumAngleDegrees &&
                     q.aspect <= options.quality.maximumAspectRatio;
    }
    // Cheap shape/size/normal tests reject early. A cell that would otherwise
    // pass must also satisfy the same bounded-face distance used by final QA,
    // so a real geometry defect participates in front and repair decisions.
    if (q.accepted && options.adaptive.enabled && chart.geometry) {
        const double deviation =
            sampledProjectionDeviation(chart, vertices, options.adaptive.geometricTolerance) /
            options.adaptive.geometricTolerance;
        q.priority = std::max(q.priority, deviation);
        q.accepted = deviation <= 1;
        q.geometryViolationEnergy += std::pow(std::max(0., deviation - 1), 2);
    }
    return q;
}
Quality quality(const SurfaceChart &chart, Face face, const MeshingOptions &options) {
    std::array<std::pair<std::uint64_t, std::uint64_t>, 3> identity;
    for (int i = 0; i < 3; ++i) {
        const auto &vertex = face->vertex(i)->info();
        identity[i] = {vertex.insertionId, vertex.positionRevision};
    }
    auto &data = face->info();
    if (!data.cachedQuality || data.qualityIdentity != identity) {
        data.cachedQuality =
            quality(chart, {{face->vertex(0), face->vertex(1), face->vertex(2)}}, options);
        data.qualityIdentity = identity;
    }
    return *data.cachedQuality;
}
struct FrontItem {
    std::array<Vertex, 3> vertices;
    std::uint64_t serial{};
    double priority{};
    bool operator<(const FrontItem &other) const {
        if (priority != other.priority)
            return priority < other.priority;
        return serial > other.serial;
    }
};
bool activeEdge(ConstraintTriangulation &tri, const SurfaceChart &chart, Face f, int edge,
                const MeshingOptions &options) {
    const auto neighbor = f->neighbor(edge);
    return tri.is_constrained({f, edge}) || !neighbor->info().inside() ||
           quality(chart, neighbor, options).accepted;
}
ChartPoint nativeTangentCandidate(const SurfaceChart &chart, Vec3 target, ChartPoint uv) {
    // The tangent construction is physical; this bounded damped inversion only
    // proposes native coordinates. Domain, constraints and physical triangle
    // quality are checked independently before and after any insertion.
    for (unsigned iteration = 0; iteration < 8; ++iteration) {
        if (!chart.containsParameterDomain(uv))
            break;
        SurfaceChartSample value;
        try {
            value = sample(chart, uv);
        } catch (const std::out_of_range &) {
            // An optional off-centre can be outside a source modelling halo
            // despite lying in its UV rectangle. Keep it for the caller's
            // ordinary domain/boundary test; never enlarge that halo here.
            return uv;
        }
        const auto residual = subtract(target, value.position);
        const auto metric = Metric::from({value.du, value.dv, value.targetEdgeLength});
        const Vec3 first{value.du.x / metric.a, value.du.y / metric.a, value.du.z / metric.a};
        const Vec3 second{(value.dv.x - metric.b * first.x) / metric.c,
                          (value.dv.y - metric.b * first.y) / metric.c,
                          (value.dv.z - metric.b * first.z) / metric.c};
        const auto step = metric.inverse({dot(residual, first), dot(residual, second)});
        bool moved = false;
        for (unsigned backtrack = 0; backtrack < 8; ++backtrack) {
            const double fraction = std::ldexp(1., -static_cast<int>(backtrack));
            const ChartPoint candidate{uv.u + fraction * step.u, uv.v + fraction * step.v};
            if (!chart.containsParameterDomain(candidate))
                continue;
            try {
                if (!(norm(subtract(sample(chart, candidate).position, target)) < norm(residual)))
                    continue;
            } catch (const std::out_of_range &) {
                continue;
            }
            uv = candidate;
            moved = true;
            break;
        }
        if (!moved)
            break;
    }
    return uv;
}
// Rebay (JCP 106, 1993, doi:10.1006/jcph.1993.1097): front points lie on
// the dual edge normal to an active edge. The off-centre is limited by the
// circumcentre so the insertion retains Delaunay-refinement separation.
std::vector<ChartPoint> candidates(const SurfaceChart &chart, Face face, int edge, bool constrained,
                                   const MeshingOptions &options) {
    const auto a = face->vertex((edge + 1) % 3)->info().value.uv;
    const auto b = face->vertex((edge + 2) % 3)->info().value.uv;
    const auto c = face->vertex(edge)->info().value.uv;
    const auto middle = midpoint(a, b);
    const auto s = [&] {
        if (constrained)
            return chart.boundaryMetric(middle);
        const auto interior = sample(chart, middle);
        return SurfaceChartMetric{interior.du, interior.dv, interior.targetEdgeLength};
    }();
    const auto metric = Metric::from(s);
    auto ea = metric.forward({a.u - middle.u, a.v - middle.v});
    auto eb = metric.forward({b.u - middle.u, b.v - middle.v});
    auto ec = metric.forward({c.u - middle.u, c.v - middle.v});
    Vec3 physicalMiddle{}, first{}, second{};
    if (chart.parameterization) {
        const auto pa = face->vertex((edge + 1) % 3)->info().value.position;
        const auto pb = face->vertex((edge + 2) % 3)->info().value.position;
        physicalMiddle = {std::midpoint(pa.x, pb.x), std::midpoint(pa.y, pb.y),
                          std::midpoint(pa.z, pb.z)};
        first = {s.du.x / metric.a, s.du.y / metric.a, s.du.z / metric.a};
        second = {(s.dv.x - metric.b * first.x) / metric.c,
                  (s.dv.y - metric.b * first.y) / metric.c,
                  (s.dv.z - metric.b * first.z) / metric.c};
        const auto tangent = [&](Vec3 p) {
            const auto delta = subtract(p, physicalMiddle);
            return ChartPoint{dot(delta, first), dot(delta, second)};
        };
        ea = tangent(pa);
        eb = tangent(pb);
        ec = tangent(face->vertex(edge)->info().value.position);
    }
    const double length = std::hypot(eb.u - ea.u, eb.v - ea.v);
    // A regular curved edge can project to zero in this midpoint tangent
    // plane. Subdivide its native interval instead of inventing an off-centre
    // or misreporting a zero physical edge. Protected intervals stay shared.
    if (chart.parameterization && length == 0)
        return {middle};
    if (!std::isfinite(length) || length <= 0)
        throw std::runtime_error("Zero-length frontal edge.");
    const double tx = (eb.u - ea.u) / length, ty = (eb.v - ea.v) / length;
    double nx = -ty, ny = tx;
    if (ec.u * nx + ec.v * ny < 0) {
        nx = -nx;
        ny = -ny;
    }
    const double x = ec.u * tx + ec.v * ty, y = ec.u * nx + ec.v * ny;
    if (!std::isfinite(y) || y <= 0)
        return {};
    const double centreHeight = (x * x + y * y - length * length / 4) / (2 * y);
    if (!std::isfinite(centreHeight))
        return {};
    // A small mandatory boundary segment grows gradually, instead of producing
    // one long skinny triangle directly to the global target size.
    const double side = std::max(.75 * length, std::min(s.targetEdgeLength, 1.25 * length));
    const double desired = std::sqrt(std::max(0., side * side - length * length / 4));
    // For an obtuse/right cell the circumcentre can be across/on the front
    // edge. It remains a valid Delaunay candidate inside the overall domain;
    // rejecting its sign would strand cells next to an accepted neighbor.
    const double height = std::min(desired, centreHeight);
    const auto toUV = [&](double h) {
        const auto delta = metric.inverse({nx * h, ny * h});
        const ChartPoint initial{middle.u + delta.u, middle.v + delta.v};
        if (!chart.parameterization || !chart.containsParameterDomain(initial))
            return initial;
        const Vec3 target{physicalMiddle.x + h * (nx * first.x + ny * second.x),
                          physicalMiddle.y + h * (nx * first.y + ny * second.y),
                          physicalMiddle.z + h * (nx * first.z + ny * second.z)};
        return nativeTangentCandidate(chart, target, initial);
    };
    std::vector<ChartPoint> result;
    if (options.adaptive.enabled) {
        // A physically well-shaped cell can still underresolve curvature.
        // Its circumcentre/off-centre need not sample the geometric defect.
        // Split the failing chord or cell in its native domain instead.
        const std::array<Vertex, 3> vertices{face->vertex(0), face->vertex(1), face->vertex(2)};
        for (int i = 0; i < 3; ++i) {
            const std::array<Vertex, 2> segment{vertices[i], vertices[(i + 1) % 3]};
            if (normalDeviation(chart, segment) > options.adaptive.geometricTolerance ||
                (chart.geometry &&
                 sampledProjectionDeviation(chart, segment, options.adaptive.geometricTolerance) >
                     options.adaptive.geometricTolerance))
                result.push_back(
                    midpoint(segment[0]->info().value.uv, segment[1]->info().value.uv));
        }
        if (normalDeviation(chart, vertices) > options.adaptive.geometricTolerance ||
            (chart.geometry &&
             sampledProjectionDeviation(chart, vertices, options.adaptive.geometricTolerance) >
                 options.adaptive.geometricTolerance)) {
            ChartPoint center{};
            for (const auto vertex : vertices) {
                center.u += vertex->info().value.uv.u / 3;
                center.v += vertex->info().value.uv.v / 3;
            }
            result.push_back(center);
        }
    }
    result.push_back(toUV(height));
    if (centreHeight > height * (1 + 1e-12))
        result.push_back(toUV(centreHeight));
    return result;
}
void flipConstrainedTopology(ConstraintTriangulation &tri, Face face, int edge) {
    const auto neighbor = face->neighbor(edge);
    using Key = std::array<std::uint64_t, 2>;
    const auto key = [](Face f, int e) {
        Key result{f->vertex((e + 1) % 3)->info().insertionId,
                   f->vertex((e + 2) % 3)->info().insertionId};
        if (result[1] < result[0])
            std::swap(result[0], result[1]);
        return result;
    };
    std::map<Key, bool> constraints;
    for (const auto f : {face, neighbor})
        for (int e = 0; e < 3; ++e)
            constraints[key(f, e)] = tri.is_constrained({f, e});
    // The base triangulation's flip updates topology. Constraints belong to
    // geometric edges rather than face slots and must follow those edges.
    tri.flip(face, edge);
    for (const auto f : {face, neighbor})
        for (int e = 0; e < 3; ++e) {
            const auto found = constraints.find(key(f, e));
            const bool constrained = found != constraints.end() && found->second;
            f->set_constraint(e, constrained);
            f->neighbor(e)->set_constraint(tri.mirror_index(f, e), constrained);
        }
}
double violationEnergy(const Quality &quality, const MeshingOptions &options) {
    if (quality.accepted)
        return 0;
    if (!options.quality.enabled)
        return std::isfinite(quality.aspect)
                   ? quality.sizeViolationEnergy + quality.geometryViolationEnergy
                   : std::numeric_limits<double>::infinity();
    const double angle = std::max(0., options.quality.minimumAngleDegrees / quality.angle - 1);
    const double aspect = std::max(0., quality.aspect / options.quality.maximumAspectRatio - 1);
    return angle * angle + aspect * aspect + quality.sizeViolationEnergy +
           quality.geometryViolationEnergy;
}
double normalDeviation(const SurfaceChart &chart, std::span<const Vertex> vertices) {
    ChartPoint uv{};
    const auto origin = vertices.front()->info().value.position;
    Vec3 center{};
    for (const auto vertex : vertices) {
        const double weight = 1. / static_cast<double>(vertices.size());
        uv.u += vertex->info().value.uv.u * weight;
        uv.v += vertex->info().value.uv.v * weight;
        const auto delta = subtract(vertex->info().value.position, origin);
        center.x += delta.x * weight;
        center.y += delta.y * weight;
        center.z += delta.z * weight;
    }
    const auto reference = sample(chart, uv);
    const auto normal = vertices.size() == 2
                            ? cross(reference.du, reference.dv)
                            : cross(subtract(vertices[1]->info().value.position, origin),
                                    subtract(vertices[2]->info().value.position, origin));
    const double length = norm(normal);
    if (!(length > 0))
        return std::numeric_limits<double>::infinity();
    // Surface-normal error measures geometric bending, independently of the
    // tangential change of speed induced by native parameter coordinates.
    return std::abs(dot(subtract(subtract(reference.position, origin), center), normal)) / length;
}
std::optional<double> surfaceDistanceWitness(const SurfaceChart &chart, Vec3 position,
                                             ChartPoint uv, double limit) {
    auto reference = sample(chart, uv);
    double distance = norm(subtract(reference.position, position));
    // A witness need not be stationary or unique: any certified physical face
    // point proves d_min <= its distance. Bounded Gauss-Newton steps remove
    // tangential UV-speed errors without replacing the authoritative MIN query.
    constexpr unsigned maximumSteps = 6, maximumBacktracks = 6;
    for (unsigned iteration = 0; iteration <= maximumSteps; ++iteration) {
        if (distance <= limit && chart.parameterization->isInsideTrim(
                                     chart.faceId, {uv.u, chart.parameterOrientation * uv.v}))
            return distance;
        if (distance <= limit || iteration == maximumSteps)
            break;
        const auto residual = subtract(position, reference.position);
        const double uu = dot(reference.du, reference.du);
        if (!(uu > 0) || !std::isfinite(uu))
            break;
        const double projection = dot(reference.du, reference.dv) / uu;
        const Vec3 transverse{reference.dv.x - projection * reference.du.x,
                              reference.dv.y - projection * reference.du.y,
                              reference.dv.z - projection * reference.du.z};
        const double vv = dot(transverse, transverse);
        if (!(vv > 0) || !std::isfinite(vv))
            break;
        const double dv = dot(residual, transverse) / vv;
        const double du = dot(residual, reference.du) / uu - projection * dv;
        bool improved = false;
        for (unsigned backtrack = 0; backtrack < maximumBacktracks; ++backtrack) {
            const double fraction = std::ldexp(1., -static_cast<int>(backtrack));
            const ChartPoint candidate{uv.u + fraction * du, uv.v + fraction * dv};
            if (!chart.containsParameterDomain(candidate))
                continue;
            SurfaceChartSample next;
            try {
                next = sample(chart, candidate);
            } catch (const std::out_of_range &) {
                // A source modelling halo can be narrower in physical space
                // than its UV bounds. Optional witness steps cannot enlarge it.
                continue;
            }
            const double nextDistance = norm(subtract(next.position, position));
            if (!(nextDistance < distance))
                continue;
            uv = candidate;
            reference = next;
            distance = nextDistance;
            improved = true;
            break;
        }
        if (!improved)
            break;
    }
    return {};
}
double sampledProjectionDeviation(const SurfaceChart &chart, std::span<const Vertex> vertices,
                                  double upperBoundLimit) {
    if (!chart.geometry)
        return 0;
    double maximum = 0;
    const auto measure = [&](Vec3 position, ChartPoint uv) {
        if (upperBoundLimit > 0 && chart.parameterization) {
            if (const auto upper = surfaceDistanceWitness(chart, position, uv, upperBoundLimit))
                return *upper;
        }
        // This is distance to the physical bounded face, not an interior map.
        // A valid boundary chord can have its unrestricted stationary minimum
        // beyond a finite natural spline domain. The nearest trim still gives
        // the correct geometric approximation error without inventing a map.
        const auto distance = chart.geometry->distanceToSurface(chart.faceId, position);
        if (!std::isfinite(distance) || distance < 0)
            throw std::runtime_error("Non-finite physical surface deviation projection.");
        return distance;
    };
    if (vertices.size() == 2) {
        const auto a = vertices[0]->info().value.position, b = vertices[1]->info().value.position;
        return measure({std::midpoint(a.x, b.x), std::midpoint(a.y, b.y), std::midpoint(a.z, b.z)},
                       midpoint(vertices[0]->info().value.uv, vertices[1]->info().value.uv));
    }
    const auto origin = vertices[0]->info().value.position;
    Vec3 center{};
    ChartPoint centerUv{};
    for (int i = 0; i < 3; ++i) {
        const auto a = vertices[i]->info().value.position;
        const auto b = vertices[(i + 1) % 3]->info().value.position;
        maximum = std::max(maximum, measure({std::midpoint(a.x, b.x), std::midpoint(a.y, b.y),
                                             std::midpoint(a.z, b.z)},
                                            midpoint(vertices[i]->info().value.uv,
                                                     vertices[(i + 1) % 3]->info().value.uv)));
        const auto delta = subtract(a, origin);
        center.x += delta.x / 3;
        center.y += delta.y / 3;
        center.z += delta.z / 3;
        centerUv.u += vertices[i]->info().value.uv.u / 3;
        centerUv.v += vertices[i]->info().value.uv.v / 3;
    }
    maximum =
        std::max(maximum, measure({origin.x + center.x, origin.y + center.y, origin.z + center.z},
                                  centerUv));
    return maximum;
}
std::size_t optimizePhysicalQuality(ConstraintTriangulation &tri, const SurfaceChart &chart,
                                    const MeshingOptions &options, std::stop_token stop) {
    std::size_t flips = 0;
    const auto aligned = [&](const std::array<Vertex, 3> &vertices, Vec3 previousNormal) {
        const auto p = vertices[0]->info().value.position;
        const auto n = cross(subtract(vertices[1]->info().value.position, p),
                             subtract(vertices[2]->info().value.position, p));
        ChartPoint uv{};
        for (auto v : vertices) {
            uv.u += v->info().value.uv.u / 3;
            uv.v += v->info().value.uv.v / 3;
        }
        const auto local = sample(chart, uv);
        return dot(n, previousNormal) > 0 && dot(n, cross(local.du, local.dv)) > 0;
    };
    // Final quality swaps strictly reduce dimensionless physical violations.
    // An improving intermediate triangle can still fail a hard criterion;
    // the final mesh must pass every criterion. Swaps need not preserve
    // Euclidean Delaunayhood; no CDT insertion follows this bounded stage.
    for (unsigned pass = 0; pass < 5; ++pass) {
        checkStop(stop);
        using Key = std::array<std::uint64_t, 2>;
        std::map<Key, std::array<Vertex, 2>> candidates;
        for (auto f = tri.finite_faces_begin(); f != tri.finite_faces_end(); ++f) {
            checkStop(stop);
            if (f->info().inside() && !quality(chart, f, options).accepted)
                for (int e = 0; e < 3; ++e)
                    if (!tri.is_constrained({f, e}) && f->neighbor(e)->info().inside()) {
                        auto a = f->vertex((e + 1) % 3), b = f->vertex((e + 2) % 3);
                        Key key{a->info().insertionId, b->info().insertionId};
                        if (key[1] < key[0])
                            std::swap(key[0], key[1]);
                        candidates.emplace(key, std::array<Vertex, 2>{a, b});
                    }
        }
        bool changed = false;
        for (const auto &entry : candidates) {
            checkStop(stop);
            Face f;
            int e;
            if (!tri.is_edge(entry.second[0], entry.second[1], f, e) ||
                tri.is_constrained({f, e}) || !f->info().inside() ||
                !f->neighbor(e)->info().inside())
                continue;
            const auto other = f->neighbor(e);
            const int opposite = tri.mirror_index(f, e);
            const auto a = f->vertex(e), b = f->vertex((e + 1) % 3), c = f->vertex((e + 2) % 3),
                       d = other->vertex(opposite);
            const std::array<Vertex, 3> first{a, b, d}, second{a, d, c};
            // Convexity is checked in the actual CDT coordinates. A flip must
            // retain the exact polygon and all shared boundary constraints.
            if (CGAL::orientation(a->point(), b->point(), d->point()) != CGAL::LEFT_TURN ||
                CGAL::orientation(a->point(), d->point(), c->point()) != CGAL::LEFT_TURN)
                continue;
            const auto oldFirst = quality(chart, f, options);
            const auto oldSecond = quality(chart, other, options);
            if (oldFirst.accepted && oldSecond.accepted)
                continue;
            const double oldEnergy =
                violationEnergy(oldFirst, options) + violationEnergy(oldSecond, options);
            const double newEnergy = violationEnergy(quality(chart, first, options), options) +
                                     violationEnergy(quality(chart, second, options), options);
            // A relative roundoff guard prevents numerically reversible swaps.
            const double roundoff = 64 * std::numeric_limits<double>::epsilon();
            if (!std::isfinite(newEnergy) || !(newEnergy < oldEnergy) ||
                (std::isfinite(oldEnergy) && oldEnergy - newEnergy <= roundoff * oldEnergy))
                continue;
            const auto normal = [&](Face face) {
                return cross(subtract(face->vertex(1)->info().value.position,
                                      face->vertex(0)->info().value.position),
                             subtract(face->vertex(2)->info().value.position,
                                      face->vertex(0)->info().value.position));
            };
            const auto n0 = normal(f), n1 = normal(other);
            const Vec3 direction{n0.x + n1.x, n0.y + n1.y, n0.z + n1.z};
            if (!aligned(first, direction) || !aligned(second, direction))
                continue;
            flipConstrainedTopology(tri, f, e);
            f->info().nesting = other->info().nesting = 1;
            ++flips;
            changed = true;
        }
        if (!changed)
            break;
    }
    return flips;
}

struct InteriorRing {
    std::vector<Face> faces;
    std::vector<Vec3> normals;
    std::map<std::uint64_t, Vertex> neighbors;
    double energy{};
};
std::optional<InteriorRing> interiorRing(ConstraintTriangulation &tri, Vertex vertex,
                                         const SurfaceChart &chart, const MeshingOptions &options) {
    InteriorRing ring;
    auto current = tri.incident_faces(vertex), end = current;
    do {
        const Face face = current;
        if (tri.is_infinite(face) || !face->info().inside())
            return {};
        for (int i = 0; i < 3; ++i) {
            if (face->vertex(i) != vertex)
                ring.neighbors.emplace(face->vertex(i)->info().insertionId, face->vertex(i));
            if (face->vertex(i) != vertex && tri.is_constrained({face, i}))
                return {};
        }
        ring.faces.push_back(face);
        ring.energy += violationEnergy(quality(chart, face, options), options);
        const auto p = face->vertex(0)->info().value.position;
        ring.normals.push_back(cross(subtract(face->vertex(1)->info().value.position, p),
                                     subtract(face->vertex(2)->info().value.position, p)));
        ++current;
    } while (current != end);
    return ring;
}
bool preservesUvOrientation(const InteriorRing &ring, Vertex vertex, const Kernel::Point_2 &point) {
    for (const auto face : ring.faces) {
        std::array<Kernel::Point_2, 3> points;
        for (int i = 0; i < 3; ++i)
            points[i] = face->vertex(i) == vertex ? point : face->vertex(i)->point();
        if (CGAL::orientation(points[0], points[1], points[2]) != CGAL::LEFT_TURN)
            return false;
    }
    return true;
}
bool preservesPhysicalRing(const InteriorRing &ring, const SurfaceChart &chart) {
    for (std::size_t j = 0; j < ring.faces.size(); ++j) {
        const auto face = ring.faces[j];
        const auto origin = face->vertex(0)->info().value.position;
        const auto n = cross(subtract(face->vertex(1)->info().value.position, origin),
                             subtract(face->vertex(2)->info().value.position, origin));
        ChartPoint center{};
        for (int i = 0; i < 3; ++i) {
            const auto &data = face->vertex(i)->info().value;
            center.u += data.uv.u / 3;
            center.v += data.uv.v / 3;
        }
        const auto reference = sample(chart, center);
        if (dot(n, ring.normals[j]) <= 0 || dot(n, cross(reference.du, reference.dv)) <= 0)
            return false;
    }
    return true;
}
class TrialVertexPosition {
  public:
    TrialVertexPosition(Vertex vertex, ChartPoint uv, Vec3 position)
        : vertex_(vertex), previous_(vertex->info().value) {
        vertex_->info().value.uv = uv;
        vertex_->info().value.position = position;
        ++vertex_->info().positionRevision;
    }
    ~TrialVertexPosition() {
        if (!committed_) {
            vertex_->info().value = previous_;
            ++vertex_->info().positionRevision;
        }
    }
    void commit(const Kernel::Point_2 &point) {
        vertex_->set_point(point);
        committed_ = true;
    }

  private:
    Vertex vertex_;
    FrontalVertex previous_;
    bool committed_{};
};
std::optional<double> ringEnergyAt(const InteriorRing &ring, Vertex vertex, ChartPoint uv,
                                   const SurfaceChart &chart, const MeshingOptions &options,
                                   const ParameterFrame &frame) {
    if (!chart.containsParameterDomain(uv) || !chart.contains(uv) ||
        !preservesUvOrientation(ring, vertex, frame.point(uv)))
        return {};
    TrialVertexPosition trial(vertex, uv, sample(chart, uv).position);
    double energy = 0;
    for (const auto face : ring.faces)
        energy += violationEnergy(quality(chart, face, options), options);
    return std::isfinite(energy) ? std::optional{energy} : std::nullopt;
}
std::optional<ChartPoint> qualityDescentGoal(const InteriorRing &ring, Vertex vertex,
                                             const SurfaceChart &chart,
                                             const MeshingOptions &options,
                                             const ParameterFrame &frame, std::stop_token stop) {
    const auto original = vertex->info().value;
    const auto local = sample(chart, original.uv);
    const auto metric = Metric::from({local.du, local.dv, local.targetEdgeLength});
    double shortest = std::numeric_limits<double>::infinity();
    for (const auto &entry : ring.neighbors)
        shortest = std::min(shortest,
                            norm(subtract(entry.second->info().value.position, original.position)));
    if (!(shortest > 0) || !std::isfinite(shortest))
        return {};
    // Central differences balance O(h^2) truncation and O(epsilon/h) rounding.
    // Directions and the trust radius use physical metres, independent of UV
    // units. This is a derivative step, never a triangle-quality tolerance.
    const double difference = std::cbrt(std::numeric_limits<double>::epsilon()) * shortest;
    std::array<double, 2> gradient{};
    for (unsigned axis = 0; axis < 2; ++axis) {
        std::array<double, 2> energies{};
        for (unsigned side = 0; side < 2; ++side) {
            checkStop(stop);
            const double step = side ? difference : -difference;
            const auto delta = metric.inverse({axis == 0 ? step : 0, axis == 1 ? step : 0});
            const auto energy =
                ringEnergyAt(ring, vertex, {original.uv.u + delta.u, original.uv.v + delta.v},
                             chart, options, frame);
            if (!energy)
                return {};
            energies[side] = *energy;
        }
        gradient[axis] = (energies[1] - energies[0]) / (2 * difference);
    }
    const double magnitude = std::hypot(gradient[0], gradient[1]);
    if (!std::isfinite(magnitude) || magnitude == 0)
        return {};
    // A trial move is bounded by the shortest existing incident edge; the
    // dyadic search below accepts only actual surface/energy improvements.
    const auto step =
        metric.inverse({-shortest * gradient[0] / magnitude, -shortest * gradient[1] / magnitude});
    return ChartPoint{original.uv.u + step.u, original.uv.v + step.v};
}
bool relaxInteriorVertex(ConstraintTriangulation &tri, Vertex vertex, const SurfaceChart &chart,
                         const MeshingOptions &options, const ParameterFrame &frame,
                         std::stop_token stop, bool fineSearch) {
    const auto ring = interiorRing(tri, vertex, chart, options);
    if (!ring || !(ring->energy > 0))
        return false;
    const auto original = vertex->info().value;
    const auto goal = qualityDescentGoal(*ring, vertex, chart, options, frame, stop);
    if (!goal)
        return false;
    double bestEnergy = ring->energy;
    std::optional<FrontalVertex> best;
    // 2^-17 remains above the cbrt(epsilon) finite-difference length. Finer
    // trials would resolve less than the direction estimate can distinguish.
    for (unsigned step = 0; step < (fineSearch ? 18u : 8u); ++step) {
        checkStop(stop);
        const double fraction = std::ldexp(1., -static_cast<int>(step));
        const ChartPoint uv{original.uv.u + fraction * (goal->u - original.uv.u),
                            original.uv.v + fraction * (goal->v - original.uv.v)};
        if (!chart.containsParameterDomain(uv) || !chart.contains(uv))
            continue;
        const auto point = frame.point(uv);
        if (!preservesUvOrientation(*ring, vertex, point))
            continue;
        const auto value = sample(chart, uv);
        TrialVertexPosition trial(vertex, uv, value.position);
        double energy = 0;
        for (const auto face : ring->faces)
            energy += violationEnergy(quality(chart, face, options), options);
        // Numerical progress guard only; acceptance still uses the exact
        // configured angle, aspect, size and geometric thresholds.
        const double roundoff = std::sqrt(std::numeric_limits<double>::epsilon());
        if (!std::isfinite(energy) || !(energy < bestEnergy) ||
            (std::isfinite(bestEnergy) && bestEnergy - energy <= roundoff * bestEnergy) ||
            !preservesPhysicalRing(*ring, chart))
            continue;
        best = vertex->info().value;
        bestEnergy = energy;
    }
    if (!best)
        return false;
    // The first improving full step can overshoot the useful local minimum
    // and cause oscillating sweeps. Choose the best checked step instead.
    TrialVertexPosition accepted(vertex, best->uv, best->position);
    accepted.commit(frame.point(best->uv));
    return true;
}
std::size_t relaxPhysicalQuality(ConstraintTriangulation &tri, const SurfaceChart &chart,
                                 const MeshingOptions &options, const ParameterFrame &frame,
                                 std::size_t &vertexTrials, std::stop_token stop,
                                 bool fineSearch = false) {
    std::map<std::uint64_t, Vertex> vertices;
    for (auto face = tri.finite_faces_begin(); face != tri.finite_faces_end(); ++face) {
        checkStop(stop);
        if (face->info().inside() && !quality(chart, face, options).accepted)
            for (int i = 0; i < 3; ++i)
                if (face->vertex(i)->info().value.atlasVertex ==
                    std::numeric_limits<std::uint32_t>::max())
                    vertices.emplace(face->vertex(i)->info().insertionId, face->vertex(i));
    }
    std::size_t moved = 0;
    for (const auto &entry : vertices) {
        checkStop(stop);
        if (++vertexTrials > maximumQualityVertexTrialsPerChart)
            throw std::runtime_error("Frontal-Delaunay face " + std::to_string(chart.faceId) +
                                     " exhausted its quality-repair work budget; the current "
                                     "triangles are not accepted.");
        moved += relaxInteriorVertex(tri, entry.second, chart, options, frame, stop, fineSearch);
    }
    return moved;
}
std::vector<AtlasSegment> outsideBoundaryRequests(Face face, ChartPoint candidate,
                                                  const std::vector<BoundarySegment> &segments) {
    ChartPoint center{};
    for (int i = 0; i < 3; ++i) {
        center.u += face->vertex(i)->info().value.uv.u / 3;
        center.v += face->vertex(i)->info().value.uv.v / 3;
    }
    // An outside circumcentre can encroach a different boundary than the bad
    // cell's own edges. Only the first crossed barrier can block this route;
    // unrelated segments outside the domain must not trigger extra refinement.
    const Kernel::Segment_2 ray({center.u, center.v}, {candidate.u, candidate.v});
    const auto intersects = Kernel{}.do_intersect_2_object();
    double nearest = std::numeric_limits<double>::infinity();
    std::vector<const BoundarySegment *> crossed;
    const double dx = candidate.u - center.u, dy = candidate.v - center.v;
    constexpr double roundoff = 128 * std::numeric_limits<double>::epsilon();
    for (const auto &segment : segments) {
        const auto a = segment.start, b = segment.end;
        if (!intersects(ray, Kernel::Segment_2({a.u, a.v}, {b.u, b.v})))
            continue;
        const double denominator = dx * (b.v - a.v) - dy * (b.u - a.u);
        if (denominator == 0)
            continue;
        const double t =
            ((a.u - center.u) * (b.v - a.v) - (a.v - center.v) * (b.u - a.u)) / denominator;
        if (!std::isfinite(t) || t < 0 || t > 1)
            continue;
        if (t < nearest - roundoff) {
            nearest = t;
            crossed.clear();
        }
        if (std::abs(t - nearest) <= roundoff)
            crossed.push_back(&segment);
    }
    std::vector<AtlasSegment> result;
    for (const auto segment : crossed)
        if (segment->encroached(candidate))
            result.push_back(segment->atlas);
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}
std::vector<AtlasSegment>
currentBoundaryRequests(ConstraintTriangulation &tri, const SurfaceChart &chart,
                        const ParameterFrame &frame, const std::vector<BoundarySegment> &segments,
                        const MeshingOptions &options, std::stop_token stop) {
    std::set<AtlasSegment> requests;
    for (auto f = tri.finite_faces_begin(); f != tri.finite_faces_end(); ++f) {
        checkStop(stop);
        if (!f->info().inside())
            continue;
        const auto q = quality(chart, f, options);
        if (q.accepted)
            continue;
        // An obtuse boundary ear can have every circumcentre outside the
        // domain. Its existing opposite vertex already encroaches a shared
        // segment; requiring a viable new point first would miss that standard
        // Delaunay-refinement case entirely.
        std::vector<AtlasSegment> existingEncroachment;
        for (int e = 0; e < 3; ++e) {
            if (!tri.is_constrained({f, e}))
                continue;
            AtlasSegment key{f->vertex((e + 1) % 3)->info().value.atlasVertex,
                             f->vertex((e + 2) % 3)->info().value.atlasVertex};
            if (key[1] < key[0])
                std::swap(key[0], key[1]);
            for (const auto &segment : segments)
                if (segment.atlas == key && segment.encroached(f->vertex(e)->info().value.uv))
                    existingEncroachment.push_back(key);
        }
        if (!existingEncroachment.empty()) {
            requests.insert(existingEncroachment.begin(), existingEncroachment.end());
            continue;
        }
        std::vector<AtlasSegment> best;
        for (int e = 0; e < 3; ++e) {
            if (!activeEdge(tri, chart, f, e, options))
                continue;
            for (const auto uv : candidates(chart, f, e, tri.is_constrained({f, e}), options)) {
                if (!std::isfinite(uv.u) || !std::isfinite(uv.v))
                    continue;
                if (!chart.containsParameterDomain(uv) || !chart.contains(uv)) {
                    auto outside = outsideBoundaryRequests(f, uv, segments);
                    if (!outside.empty() && (best.empty() || outside.size() < best.size() ||
                                             (outside.size() == best.size() && outside < best)))
                        best = std::move(outside);
                    continue;
                }
                Triangulation::Locate_type type;
                int index;
                auto located = tri.locate(frame.point(uv), type, index, f);
                if (type == Triangulation::VERTEX)
                    continue;
                if (!located->info().inside()) {
                    if (type != Triangulation::EDGE || !located->neighbor(index)->info().inside())
                        continue;
                }
                const auto s = sample(chart, uv);
                const double separation = .2 * std::min(s.targetEdgeLength, q.shortest);
                bool close = false;
                // Physical quality swaps need not preserve UV Delaunayhood, so
                // a final diagnostic must not use a Delaunay conflict traversal.
                for (auto v = tri.finite_vertices_begin(); v != tri.finite_vertices_end(); ++v)
                    if (norm(subtract(s.position, v->info().value.position)) < separation) {
                        close = true;
                        break;
                    }
                if (close)
                    continue;
                std::vector<AtlasSegment> blocked;
                for (const auto &segment : segments)
                    if (segment.encroached(uv))
                        blocked.push_back(segment.atlas);
                if (blocked.empty())
                    continue;
                std::sort(blocked.begin(), blocked.end());
                blocked.erase(std::unique(blocked.begin(), blocked.end()), blocked.end());
                // One viable candidate per currently rejected cell is enough;
                // splitting the entire historical candidate set wastes elements.
                if (best.empty() || blocked.size() < best.size() ||
                    (blocked.size() == best.size() && blocked < best))
                    best = std::move(blocked);
            }
        }
        requests.insert(best.begin(), best.end());
    }
    return {requests.begin(), requests.end()};
}
} // namespace

BoundaryRefinementRequired::BoundaryRefinementRequired(std::string message,
                                                       std::vector<AtlasSegment> segments,
                                                       std::optional<double> violationEnergy)
    : std::runtime_error(std::move(message)), segments_(std::move(segments)),
      violationEnergy_(violationEnergy) {
    if (violationEnergy_ && (std::isnan(*violationEnergy_) || *violationEnergy_ <= 0))
        throw std::invalid_argument("Boundary refinement needs a positive quality violation.");
    for (auto &segment : segments_) {
        if (segment[0] == segment[1])
            throw std::invalid_argument("A boundary refinement segment has equal endpoints.");
        if (segment[1] < segment[0])
            std::swap(segment[0], segment[1]);
    }
    std::sort(segments_.begin(), segments_.end());
    segments_.erase(std::unique(segments_.begin(), segments_.end()), segments_.end());
    if (segments_.empty())
        throw std::invalid_argument("Boundary refinement requires a current blocking segment.");
}

namespace {
class FrontalTriangulator {
  public:
    FrontalTriangulator(const SurfaceChart &inputChart, std::span<const Vec3> inputVertices,
                        const MeshingOptions &inputOptions, std::stop_token inputStop)
        : chart(inputChart), atlasVertices(inputVertices), options(inputOptions), stop(inputStop) {
        checkStop(stop);
        if (chart.loops.empty() || chart.reference.empty() ||
            !std::isfinite(options.adaptive.geometricTolerance) ||
            options.adaptive.geometricTolerance <= 0 ||
            !std::isfinite(options.quality.minimumAngleDegrees) ||
            options.quality.minimumAngleDegrees <= 0 || options.quality.minimumAngleDegrees >= 60 ||
            !std::isfinite(options.quality.maximumAspectRatio) ||
            options.quality.maximumAspectRatio < 2 / std::sqrt(3.))
            throw std::invalid_argument("Missing chart boundary or invalid frontal criteria.");
        const auto &seed = chart.reference.front().uv;
        const ChartPoint origin{(seed[0].u + seed[1].u + seed[2].u) / 3,
                                (seed[0].v + seed[1].v + seed[2].v) / 3};
        frame = {origin, representativeMetric(chart, stop)};
        // Reject an already inconsistent exact-source boundary before trying
        // to repair it by adding interior nodes. Such refinement cannot move
        // a mandatory atlas point onto a different provider's physical face.
        if (chart.geometry && options.adaptive.enabled && !chart.loops.front().vertices.empty()) {
            const auto id = chart.loops.front().vertices.front().atlasVertex;
            if (id < atlasVertices.size()) {
                const double deviation =
                    chart.geometry->distanceToSurface(chart.faceId, atlasVertices[id]);
                if (!std::isfinite(deviation) || deviation < 0 ||
                    deviation > options.adaptive.geometricTolerance) {
                    std::ostringstream message;
                    message << "Source chart " << chart.faceId
                            << " mandatory boundary has sampled physical deviation " << deviation
                            << " m; interior refinement cannot repair that source mismatch.";
                    throw std::runtime_error(message.str());
                }
            }
        }
    }

    FrontalChartMesh run() {
        initializeBoundary();
        initializeFront();
        refineFront();
        return optimizeAndExport();
    }

  private:
    static std::array<std::uint64_t, 3> faceKey(Face f) {
        std::array<std::uint64_t, 3> key{f->vertex(0)->info().insertionId,
                                         f->vertex(1)->info().insertionId,
                                         f->vertex(2)->info().insertionId};
        std::sort(key.begin(), key.end());
        return key;
    }

    Vertex addMandatory(const SurfaceChartVertex &v) {
        if (v.atlasVertex >= atlasVertices.size() || !std::isfinite(v.uv.u) ||
            !std::isfinite(v.uv.v) || !finite(atlasVertices[v.atlasVertex]))
            throw std::invalid_argument("Invalid shared chart vertex.");
        const auto count = triangulation.number_of_vertices();
        const auto node = triangulation.insert(frame.point(v.uv));
        if (triangulation.number_of_vertices() == count) {
            if (node->info().value.atlasVertex != v.atlasVertex)
                throw std::invalid_argument("Different atlas vertices share one chart coordinate.");
        } else {
            node->info().value = {v.uv, atlasVertices[v.atlasVertex], v.atlasVertex};
            node->info().insertionId = ++vertexSerial;
        }
        return node;
    }

    void initializeBoundary() {
        for (const auto &loop : chart.loops) {
            if (loop.vertices.size() < 3)
                throw std::invalid_argument("A chart loop has fewer than three vertices.");
            std::vector<Vertex> vertices;
            for (const auto &v : loop.vertices)
                vertices.push_back(addMandatory(v));
            for (size_t i = 0; i < vertices.size(); ++i) {
                if (vertices[i] == vertices[(i + 1) % vertices.size()])
                    throw std::invalid_argument("Repeated chart-loop vertex.");
                const auto &a = vertices[i]->info().value;
                const auto &b = vertices[(i + 1) % vertices.size()]->info().value;
                const double target = chart.targetSize(midpoint(a.uv, b.uv));
                if (norm(subtract(a.position, b.position)) > (4. / 3.) * target * (1 + 1e-12)) {
                    std::ostringstream message;
                    message << "Shared chart boundary " << a.atlasVertex << "-" << b.atlasVertex
                            << " exceeds its scalar size target on face " << chart.faceId << ".";
                    throw BoundaryRefinementRequired(message.str(),
                                                     {{a.atlasVertex, b.atlasVertex}});
                }
                const auto middle = midpoint(a.uv, b.uv);
                const auto metric = Metric::from(chart.boundaryMetric(middle));
                auto ids = std::array<std::uint32_t, 2>{a.atlasVertex, b.atlasVertex};
                if (ids[1] < ids[0])
                    std::swap(ids[0], ids[1]);
                boundarySegments.push_back({middle,
                                            metric.forward({a.uv.u - middle.u, a.uv.v - middle.v}),
                                            metric.forward({b.uv.u - middle.u, b.uv.v - middle.v}),
                                            metric, ids, a.uv, b.uv});
                const auto before = triangulation.number_of_vertices();
                triangulation.insert_constraint(vertices[i], vertices[(i + 1) % vertices.size()]);
                if (before != triangulation.number_of_vertices())
                    throw std::invalid_argument("Intersecting chart constraints.");
            }
            checkStop(stop);
        }
        for (const auto &v : chart.fixedPoints)
            addMandatory(v);
        markDomain(triangulation);
    }

    void enqueue(Face f) {
        if (!f->info().inside() || f->info().blocked)
            return;
        const auto q = quality(chart, f, options);
        if (q.accepted)
            return;
        for (int e = 0; e < 3; ++e)
            if (activeEdge(triangulation, chart, f, e, options)) {
                front.push(
                    {{f->vertex(0), f->vertex(1), f->vertex(2)}, f->info().serial, q.priority});
                return;
            }
    }

    void initializeFront() {
        std::map<std::array<std::uint64_t, 3>, Face> initialFaces;
        for (auto f = triangulation.all_faces_begin(); f != triangulation.all_faces_end(); ++f)
            initialFaces.emplace(faceKey(f), f);
        for (const auto &entry : initialFaces)
            entry.second->info().serial = ++serial;
        for (const auto &entry : initialFaces) {
            checkStop(stop);
            if (!triangulation.is_infinite(entry.second))
                enqueue(entry.second);
        }
    }

    std::vector<int> activeEdges(Face face) {
        std::vector<int> edges;
        for (int e = 0; e < 3; ++e)
            if (activeEdge(triangulation, chart, face, e, options))
                edges.push_back(e);
        std::sort(edges.begin(), edges.end(), [&](int a, int b) {
            const auto length = [&](int e) {
                return norm(subtract(face->vertex((e + 1) % 3)->info().value.position,
                                     face->vertex((e + 2) % 3)->info().value.position));
            };
            return length(a) < length(b);
        });
        return edges;
    }

    void updateFront(Vertex vertex) {
        // CGAL handles contain allocator addresses; they must not order
        // front updates or identical jobs would create different meshes.
        std::map<std::array<std::uint64_t, 3>, Face> changed;
        auto incident = triangulation.incident_faces(vertex), end = incident;
        do {
            if (!triangulation.is_infinite(incident)) {
                incident->info().nesting = 1;
                changed.emplace(faceKey(incident), incident);
            }
            ++incident;
        } while (incident != end);
        // An unchanged neighbor can become part of the front when the
        // adjacent newly inserted triangles satisfy their size/shape.
        const auto direct = changed;
        for (const auto &entry : direct)
            for (int e = 0; e < 3; ++e)
                if (!triangulation.is_infinite(entry.second->neighbor(e))) {
                    const auto f = entry.second->neighbor(e);
                    changed.emplace(faceKey(f), f);
                }
        for (const auto &entry : changed) {
            entry.second->info().serial = ++serial;
            entry.second->info().blocked = false;
        }
        for (const auto &entry : changed)
            enqueue(entry.second);
    }

    bool insertCandidate(Face face, ChartPoint uv) {
        ++result.candidateEvaluations;
        if (result.candidateEvaluations > maximumCandidateEvaluationsPerChart ||
            result.insertedInteriorVertices >= maximumInteriorInsertionsPerChart) {
            std::ostringstream message;
            message << "Frontal-Delaunay chart " << chart.faceId
                    << " exhausted its work budget (1000000 interior insertions / "
                       "5000000 candidate evaluations). Check chart distortion and "
                       "boundary constraints; no mesh was accepted.";
            throw std::runtime_error(message.str());
        }
        if (!chart.containsParameterDomain(uv) || !chart.contains(uv))
            return false;
        const auto point = frame.point(uv);
        Triangulation::Locate_type location;
        int localIndex;
        const auto located = triangulation.locate(point, location, localIndex, face);
        if (location == Triangulation::VERTEX || !located->info().inside() ||
            (location == Triangulation::EDGE &&
             triangulation.is_constrained({located, localIndex})))
            return false;
        // A constrained CDT permits encroached segments, but inserting
        // inside their diametral disks creates skinny triangles which
        // cannot be repaired while those shared segments stay fixed.
        // Use each segment's physical tangent metric, not raw UV units.
        bool encroached = false;
        for (const auto &segment : boundarySegments)
            if (segment.encroached(uv)) {
                encroached = true;
                break;
            }
        if (encroached)
            return false;
        const auto s = sample(chart, uv);
        std::vector<Face> conflicts;
        triangulation.get_conflicts(point, std::back_inserter(conflicts), located);
        // A point constructed in a local surface metric need not enter the
        // selected cell's global CDT cavity. Inserting it elsewhere leaves
        // that bad cell untouched and can cause an endless refinement cascade.
        if (std::find(conflicts.begin(), conflicts.end(), face) == conflicts.end())
            return false;
        const double separation =
            .2 * std::min(s.targetEdgeLength, quality(chart, face, options).shortest);
        bool close = false;
        for (const auto other : conflicts)
            for (int i = 0; i < 3; ++i)
                if (!triangulation.is_infinite(other->vertex(i)) &&
                    norm(subtract(s.position, other->vertex(i)->info().value.position)) <
                        separation)
                    close = true;
        if (close)
            return false;
        const auto before = triangulation.number_of_vertices();
        const auto vertex = triangulation.insert(point, located);
        if (triangulation.number_of_vertices() == before)
            return false;
        vertex->info().value = {uv, s.position, std::numeric_limits<std::uint32_t>::max()};
        vertex->info().insertionId = ++vertexSerial;
        ++result.insertedInteriorVertices;
        updateFront(vertex);
        return true;
    }

    void refineFront() {
        while (!front.empty()) {
            checkStop(stop);
            const auto item = front.top();
            front.pop();
            Face face;
            // Handles to removed faces are never dereferenced. Vertices survive all
            // insertions, so a queued triple can be checked against the current CDT.
            if (!triangulation.is_face(item.vertices[0], item.vertices[1], item.vertices[2],
                                       face) ||
                face->info().serial != item.serial || !face->info().inside() ||
                face->info().blocked || quality(chart, face, options).accepted)
                continue;
            bool inserted = false;
            for (const int edge : activeEdges(face)) {
                for (const auto uv : candidates(
                         chart, face, edge, triangulation.is_constrained({face, edge}), options)) {
                    if (insertCandidate(face, uv)) {
                        inserted = true;
                        break;
                    }
                }
                if (inserted)
                    break;
            }
            if (!inserted)
                face->info().blocked = true;
        }
    }

    void validateFinal(ConstraintTriangulation &finished) {
        double worstAngle = 180, worstAspect = 0, maximumDeviation = 0;
        std::size_t sizeFailures = 0, geometryFailures = 0;
        double maximumSizeEnergy = 0, maximumGeometryEnergy = 0;
        double totalViolationEnergy = 0;
        size_t rejected = 0;
        Face worstFace;
        double worstPriority = -1;
        for (auto f = finished.finite_faces_begin(); f != finished.finite_faces_end(); ++f) {
            checkStop(stop);
            if (f->info().inside()) {
                auto q = quality(chart, f, options);
                // Accepted quality already certifies every physical-distance
                // sample for this exact corner identity and position revision.
                // Rejected cells still need distances for diagnostic priority.
                if (!q.accepted && options.adaptive.enabled && chart.geometry) {
                    const std::array<Vertex, 3> vertices{f->vertex(0), f->vertex(1), f->vertex(2)};
                    // A certified point on the trimmed face proves the same
                    // distance gate: d_min <= |chord - S(uv)| <= tolerance.
                    // Only inconclusive witnesses require a global MIN search.
                    const double deviation = sampledProjectionDeviation(
                        chart, vertices, options.adaptive.geometricTolerance);
                    maximumDeviation = std::max(maximumDeviation, deviation);
                    q.accepted = q.accepted && deviation <= options.adaptive.geometricTolerance;
                    q.priority =
                        std::max(q.priority, deviation / options.adaptive.geometricTolerance);
                }
                if (!q.accepted) {
                    totalViolationEnergy += violationEnergy(q, options);
                    sizeFailures += q.sizeViolationEnergy > 0;
                    geometryFailures += q.geometryViolationEnergy > 0;
                    maximumSizeEnergy = std::max(maximumSizeEnergy, q.sizeViolationEnergy);
                    maximumGeometryEnergy =
                        std::max(maximumGeometryEnergy, q.geometryViolationEnergy);
                    ++rejected;
                    worstAngle = std::min(worstAngle, q.angle);
                    worstAspect = std::max(worstAspect, q.aspect);
                    if (q.priority > worstPriority ||
                        (q.priority == worstPriority && faceKey(f) < faceKey(worstFace))) {
                        worstPriority = q.priority;
                        worstFace = f;
                    }
                }
            }
        }
        if (rejected) {
            std::ostringstream message;
            message << std::setprecision(17) << "Frontal-Delaunay chart " << chart.faceId
                    << " cannot meet shared-boundary size/quality constraints: " << rejected
                    << " triangles, min angle " << worstAngle << ", max aspect " << worstAspect
                    << ", sampled physical deviation " << maximumDeviation << " m"
                    << ", size failures " << sizeFailures << " (energy " << maximumSizeEnergy
                    << "), geometric failures " << geometryFailures << " (energy "
                    << maximumGeometryEnergy << ")"
                    << ". Refine the common atlas boundary or improve its parametrization.";
            message << " Worst rejected cell:";
            for (int i = 0; i < 3; ++i) {
                const auto &v = worstFace->vertex(i)->info().value;
                message << " [atlas=" << v.atlasVertex << " uv=" << v.uv.u << "," << v.uv.v
                        << " xyz=" << v.position.x << "," << v.position.y << "," << v.position.z
                        << "]";
            }
            auto requests =
                currentBoundaryRequests(finished, chart, frame, boundarySegments, options, stop);
            if (!requests.empty()) {
                message << " Current candidates are blocked by shared-segment encroachment:";
                for (std::size_t i = 0; i < std::min<std::size_t>(requests.size(), 8); ++i)
                    message << " " << requests[i][0] << "-" << requests[i][1];
                message
                    << ". Coordinated atlas-boundary refinement is required; no segment was split "
                       "locally.";
                throw BoundaryRefinementRequired(message.str(), std::move(requests),
                                                 totalViolationEnergy);
            }
            throw std::runtime_error(message.str());
        }
    }

    FrontalChartMesh exportMesh(ConstraintTriangulation &finished) {
        std::map<std::array<std::uint64_t, 3>, Face> outputFaces;
        for (auto f = finished.finite_faces_begin(); f != finished.finite_faces_end(); ++f)
            if (f->info().inside())
                outputFaces.emplace(faceKey(f), f);
        for (const auto &entry : outputFaces) {
            const auto f = entry.second;
            std::array<std::uint32_t, 3> indices;
            for (int i = 0; i < 3; ++i) {
                auto &v = f->vertex(i)->info();
                if (v.output == std::numeric_limits<std::uint32_t>::max()) {
                    if (result.vertices.size() >= std::numeric_limits<std::uint32_t>::max())
                        throw std::runtime_error("Frontal mesh exceeds the vertex index range.");
                    v.output = static_cast<std::uint32_t>(result.vertices.size());
                    result.vertices.push_back(v.value);
                }
                indices[i] = v.output;
            }
            result.triangles.push_back(indices);
        }
        if (result.triangles.empty())
            throw std::runtime_error("The chart contains no interior triangles.");
        checkStop(stop);
        return std::move(result);
    }

    bool
    preservesInsertionOrientation(std::span<const Face> affected, ChartPoint uv,
                                  const Kernel::Point_2 &point, Vec3 position,
                                  const std::optional<std::array<Vertex, 2>> &splitEdge) const {
        for (const auto original : affected) {
            const auto p0 = original->vertex(0)->info().value.position;
            const auto previousNormal =
                cross(subtract(original->vertex(1)->info().value.position, p0),
                      subtract(original->vertex(2)->info().value.position, p0));
            for (int e = 0; e < 3; ++e) {
                const auto a = original->vertex((e + 1) % 3), b = original->vertex((e + 2) % 3);
                if (splitEdge && ((a == (*splitEdge)[0] && b == (*splitEdge)[1]) ||
                                  (b == (*splitEdge)[0] && a == (*splitEdge)[1])))
                    continue;
                if (CGAL::orientation(a->point(), b->point(), point) != CGAL::LEFT_TURN)
                    return false;
                const auto normal =
                    cross(subtract(b->info().value.position, a->info().value.position),
                          subtract(position, a->info().value.position));
                const auto reference =
                    sample(chart, {(a->info().value.uv.u + b->info().value.uv.u + uv.u) / 3,
                                   (a->info().value.uv.v + b->info().value.uv.v + uv.v) / 3});
                if (dot(normal, previousNormal) <= 0 ||
                    dot(normal, cross(reference.du, reference.dv)) <= 0)
                    return false;
            }
        }
        return true;
    }

    bool insertPhysicalCandidate(ConstraintTriangulation &finished, Face face, ChartPoint uv,
                                 std::optional<int> bisectedEdge = {}) {
        checkStop(stop);
        if (++result.candidateEvaluations > maximumCandidateEvaluationsPerChart ||
            result.insertedInteriorVertices >= maximumInteriorInsertionsPerChart)
            throw std::runtime_error("Physical chart refinement exhausted its work budget.");
        auto point = frame.point(uv);
        ConstraintTriangulation::Locate_type location;
        int index;
        Face located;
        std::optional<std::array<Vertex, 2>> splitEdge;
        std::vector<Face> affected;
        if (bisectedEdge) {
            if (finished.is_constrained({face, *bisectedEdge}) ||
                !face->neighbor(*bisectedEdge)->info().inside())
                return false;
            splitEdge = {
                {face->vertex((*bisectedEdge + 1) % 3), face->vertex((*bisectedEdge + 2) % 3)}};
            affected = {face, face->neighbor(*bisectedEdge)};
            point = CGAL::midpoint((*splitEdge)[0]->point(), (*splitEdge)[1]->point());
            uv = frame.uv(point);
        }
        located = finished.locate(point, location, index, face);
        if (!chart.containsParameterDomain(uv) || !chart.contains(uv))
            return false;
        if (location == ConstraintTriangulation::VERTEX || !located->info().inside() ||
            (location == ConstraintTriangulation::EDGE &&
             (finished.is_constrained({located, index}) ||
              !located->neighbor(index)->info().inside())))
            return false;
        // A constrained insertion has no Delaunay conflict cavity. It must
        // actually replace this rejected cell, rather than repeatedly adding
        // unrelated points around an unchanged bad triangle.
        const bool bisectedCavity =
            bisectedEdge && (located == affected[0] || located == affected[1]);
        if (!bisectedCavity && located != face &&
            !(location == ConstraintTriangulation::EDGE && located->neighbor(index) == face))
            return false;
        const Kernel::Point_2 center{(face->vertex(0)->point().x() + face->vertex(1)->point().x() +
                                      face->vertex(2)->point().x()) /
                                         3,
                                     (face->vertex(0)->point().y() + face->vertex(1)->point().y() +
                                      face->vertex(2)->point().y()) /
                                         3};
        const Kernel::Segment_2 route(center, point);
        const auto intersects = Kernel{}.do_intersect_2_object();
        bool blocked = false;
        for (const auto &segment : boundarySegments)
            blocked = blocked || segment.encroached(uv) ||
                      intersects(route, Kernel::Segment_2(frame.point(segment.start),
                                                          frame.point(segment.end)));
        if (blocked)
            return false;
        const auto value = sample(chart, uv);
        const auto q = quality(chart, face, options);
        const double separation = .2 * std::min(value.targetEdgeLength, q.shortest);
        for (auto v = finished.finite_vertices_begin(); v != finished.finite_vertices_end(); ++v)
            if (norm(subtract(value.position, v->info().value.position)) < separation) {
                blocked = true;
                break;
            }
        if (blocked)
            return false;
        if (!bisectedEdge)
            affected = {located};
        if (!bisectedEdge && location == ConstraintTriangulation::EDGE) {
            affected.push_back(located->neighbor(index));
            splitEdge = {{located->vertex((index + 1) % 3), located->vertex((index + 2) % 3)}};
        }
        if (!preservesInsertionOrientation(affected, uv, point, value.position, splitEdge))
            return false;
        std::vector<std::array<Vertex, 3>> replaced;
        for (const auto old : affected)
            replaced.push_back({old->vertex(0), old->vertex(1), old->vertex(2)});
        const auto vertex = finished.insert(point, location, located, index);
        vertex->info().value = {uv, value.position, std::numeric_limits<std::uint32_t>::max()};
        vertex->info().insertionId = ++vertexSerial;
        if (bisectedEdge) {
            Face oldEdge;
            int oldIndex;
            if (finished.is_edge((*splitEdge)[0], (*splitEdge)[1], oldEdge, oldIndex)) {
                // An inexact midpoint can locate one ulp beside the old edge.
                // Insert with its true location, then remove that old diagonal;
                // the four final fan triangles were checked before mutation.
                flipConstrainedTopology(finished, oldEdge, oldIndex);
                ++result.qualityFlips;
            }
        }
        for (const auto &old : replaced) {
            Face retained;
            if (finished.is_face(old[0], old[1], old[2], retained))
                throw std::runtime_error(
                    "Physical refinement did not replace its selected cavity.");
        }
        auto incident = finished.incident_faces(vertex), end = incident;
        do {
            incident->info().nesting = 1;
            ++incident;
        } while (incident != end);
        ++result.insertedInteriorVertices;
        return true;
    }

    bool refinePhysicalConstraints(ConstraintTriangulation &finished, bool sizeOnly = false) {
        std::map<std::array<std::uint64_t, 3>, Face> rejected;
        for (auto f = finished.finite_faces_begin(); f != finished.finite_faces_end(); ++f) {
            checkStop(stop);
            if (f->info().inside()) {
                const auto q = quality(chart, f, options);
                if (!q.accepted && (!sizeOnly || q.sizeViolationEnergy > 0))
                    rejected.emplace(faceKey(f), f);
            }
        }
        for (const auto &entry : rejected) {
            const auto face = entry.second;
            std::optional<int> longestEdge;
            double longestLength = 0;
            for (int edge = 0; edge < 3; ++edge)
                if (!finished.is_constrained({face, edge}) &&
                    face->neighbor(edge)->info().inside()) {
                    const auto a = face->vertex((edge + 1) % 3)->info().value.position;
                    const auto b = face->vertex((edge + 2) % 3)->info().value.position;
                    const double length = norm(subtract(a, b));
                    if (length > longestLength) {
                        longestLength = length;
                        longestEdge = edge;
                    }
                }
            if (longestEdge && insertPhysicalCandidate(finished, face, {}, longestEdge))
                return true;
            for (int edge = 0; edge < 3; ++edge)
                for (const auto uv :
                     candidates(chart, face, edge, finished.is_constrained({face, edge}), options))
                    if (insertPhysicalCandidate(finished, face, uv))
                        return true;
        }
        return false;
    }

    FrontalChartMesh optimizeAndExport() {
        // Transfer the base topology without copying vertices or constraints.
        // Quality optimization is a constrained triangulation operation; keeping
        // it in a CDT object would violate that object's Euclidean invariant.
        ConstraintTriangulation finished(
            std::move(static_cast<ConstraintTriangulation &>(triangulation)));
        result.qualityFlips = optimizePhysicalQuality(finished, chart, options, stop);
        // A varying native metric can leave the affine CDT at a physical
        // edge-flip minimum. Relax only free interior nodes on the same surface;
        // every accepted step strictly reduces the complete one-ring energy.
        // Further insertion uses the constrained topology, never CDT methods;
        // all final quality and embedding gates remain mandatory.
        // An encroaching free interior node is repairable. Requesting a shared
        // split before trying that repair can reproduce the same encroachment
        // at ever smaller scales after each fresh chart triangulation.
        std::size_t vertexTrials = 0;
        double bestEnergy = std::numeric_limits<double>::infinity();
        unsigned stagnantSweeps = 0;
        for (;;) {
            double energy = 0;
            for (auto f = finished.finite_faces_begin(); f != finished.finite_faces_end(); ++f) {
                checkStop(stop);
                if (f->info().inside())
                    energy += violationEnergy(quality(chart, f, options), options);
            }
            if (energy == 0)
                break;
            // This is an optimization work limit, not a mesh tolerance.
            // Insignificant oscillation or insertion without net progress
            // must return to the unchanged final validity checks.
            if (energy < .99 * bestEnergy) {
                bestEnergy = energy;
                stagnantSweeps = 0;
            } else if (++stagnantSweeps >= 32) {
                // Tiny improving displacements can stall above a size
                // constraint. Resolve that missing resolution by replacing
                // an offending cell/cavity, rather than waiting for an
                // exactly zero gradient or accepting an oversized edge.
                if (refinePhysicalConstraints(finished, true)) {
                    do {
                        result.qualityFlips +=
                            optimizePhysicalQuality(finished, chart, options, stop);
                    } while (refinePhysicalConstraints(finished, true));
                    bestEnergy = std::numeric_limits<double>::infinity();
                    stagnantSweeps = 0;
                    continue;
                }
                break;
            }
            const auto moved =
                relaxPhysicalQuality(finished, chart, options, frame, vertexTrials, stop);
            const auto flips = optimizePhysicalQuality(finished, chart, options, stop);
            result.qualityFlips += flips;
            // Stop at an accepted mesh or a genuine fixed point, not after
            // an arbitrary sweep count while substantial descent continues.
            // A stalled chart is then checked for newly exposed blockers.
            if (!moved && !flips) {
                if (relaxPhysicalQuality(finished, chart, options, frame, vertexTrials, stop, true))
                    continue;
                if (refinePhysicalConstraints(finished))
                    result.qualityFlips += optimizePhysicalQuality(finished, chart, options, stop);
                else
                    break;
            }
        }
        if (!finished.is_valid())
            throw std::runtime_error(
                "Physical quality optimization invalidated constrained chart topology.");
        validateFinal(finished);
        return exportMesh(finished);
    }

    const SurfaceChart &chart;
    std::span<const Vec3> atlasVertices;
    const MeshingOptions &options;
    std::stop_token stop;
    ParameterFrame frame{};
    Triangulation triangulation;
    std::vector<BoundarySegment> boundarySegments;
    std::uint64_t vertexSerial{}, serial{};
    std::priority_queue<FrontItem> front;
    FrontalChartMesh result;
};
} // namespace

FrontalChartMesh meshFrontalChart(const SurfaceChart &chart, std::span<const Vec3> atlasVertices,
                                  const MeshingOptions &options, std::stop_token stop) {
    return FrontalTriangulator(chart, atlasVertices, options, stop).run();
}
} // namespace hornsculpt::detail
