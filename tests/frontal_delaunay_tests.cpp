#include "../modules/meshing/src/frontal_delaunay.hpp"
#include "fixture_geometry.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
#include <stdexcept>
#include <stop_token>
#include <thread>
using namespace hornsculpt;
using namespace hornsculpt::detail;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class Function> void rejects(Function f) {
    try {
        f();
    } catch (const std::exception &) {
        return;
    }
    throw std::runtime_error("Invalid frontal chart was accepted.");
}
Vec3 sub(Vec3 a, Vec3 b) {
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
struct Fixture {
    SurfaceChart chart;
    std::vector<Vec3> vertices;
    double expectedArea{};
};
void loopPoint(Fixture &f, SurfaceChartLoop &loop, double x, double y) {
    loop.vertices.push_back({{x, y}, static_cast<std::uint32_t>(f.vertices.size())});
    f.vertices.push_back({x, y, 0});
}
void reference(Fixture &f, std::array<ChartPoint, 3> p, std::array<double, 3> sizes) {
    ChartReferenceTriangle t;
    t.uv = p;
    t.sizes = sizes;
    for (int i = 0; i < 3; ++i)
        t.positions[i] = {p[i].u, p[i].v, 0};
    f.chart.reference.push_back(t);
}
Fixture plate(bool graded = false) {
    Fixture f;
    f.chart.faceId = 11;
    SurfaceChartLoop loop;
    std::vector<double> xs;
    if (graded) {
        const double first = std::log(5.) / .3, total = first + (.1 - .008 / .3) / .01;
        const unsigned n = static_cast<unsigned>(std::ceil(total / .8));
        for (unsigned i = 0; i <= n; ++i) {
            const double s = total * i / n;
            xs.push_back(s <= first ? .002 * std::expm1(.3 * s) / .3
                                    : .008 / .3 + (s - first) * .01);
        }
        xs.back() = .1;
    } else
        for (unsigned i = 0; i <= 10; ++i)
            xs.push_back(.01 * i);
    const unsigned left = graded ? 40 : 8, right = 8;
    for (auto x : xs)
        loopPoint(f, loop, x, 0);
    for (unsigned i = 1; i <= right; ++i)
        loopPoint(f, loop, .1, .08 * i / right);
    for (size_t i = xs.size() - 1; i > 0; --i)
        loopPoint(f, loop, xs[i - 1], .08);
    for (unsigned i = left - 1; i > 0; --i)
        loopPoint(f, loop, 0, .08 * i / left);
    f.chart.loops.push_back(std::move(loop));
    if (graded) {
        for (const auto bounds :
             std::array<std::array<double, 2>, 2>{{{0, .008 / .3}, {.008 / .3, .1}}}) {
            auto h = [](double x) { return std::min(.01, .002 + .3 * x); };
            const auto a = bounds[0], b = bounds[1];
            reference(f, {{{a, 0}, {b, 0}, {b, .08}}}, {{h(a), h(b), h(b)}});
            reference(f, {{{a, 0}, {b, .08}, {a, .08}}}, {{h(a), h(b), h(a)}});
        }
    } else {
        reference(f, {{{0, 0}, {.1, 0}, {.1, .08}}}, {{.01, .01, .01}});
        reference(f, {{{0, 0}, {.1, .08}, {0, .08}}}, {{.01, .01, .01}});
    }
    f.expectedArea = .008;
    return f;
}
class VariableSpeedPlane final : public SurfaceParameterization {
  public:
    explicit VariableSpeedPlane(double strength, bool trimEvidence = false)
        : strength_(strength), trimEvidence_(trimEvidence) {}
    double parameter(double x) const {
        return .1 * std::log1p(x / .1 * std::expm1(strength_)) / strength_;
    }
    const SurfaceMesh &sourceMesh() const noexcept override {
        return source_;
    }
    std::span<const ParameterTriangle> triangles() const noexcept override {
        return {};
    }
    ParameterDomain domain(std::uint32_t) const override {
        return {{std::array{0., .1}}, {std::array{0., .08}}, {}, {}, 1};
    }
    ParametricSurfaceSample surface(std::uint32_t, ParameterPoint p) const override {
        const double x = .1 * std::expm1(strength_ * p.u / .1) / std::expm1(strength_);
        const double dx = strength_ * std::exp(strength_ * p.u / .1) / std::expm1(strength_);
        return {{x, p.v, 0}, {dx, 0, 0}, {0, 1, 0}};
    }
    ParametricCurveSample curve(std::uint32_t, double) const override {
        throw std::logic_error("Chart fixture has no source curves.");
    }
    bool isInsideTrim(std::uint32_t, ParameterPoint p) const override {
        return trimEvidence_ && p.u > 0 && p.u < .1 && p.v > 0 && p.v < .08;
    }
    CurveUseSample curveUse(std::uint32_t, double) const override {
        throw std::logic_error("Chart fixture has no source curve uses.");
    }
    std::span<const double> curvePartition(std::uint32_t) const override {
        throw std::logic_error("Chart fixture has no source curves.");
    }

  private:
    double strength_;
    bool trimEvidence_{};
    SurfaceMesh source_;
};
class PlaneDeviationGeometry final : public SurfaceGeometry {
  public:
    explicit PlaneDeviationGeometry(double offset = 0) : offset_(offset) {}
    mutable std::size_t evaluations{};
    std::span<const BoundaryCorner> boundaryCorners() const noexcept override {
        return {};
    }
    SurfaceProjection projectSurface(std::uint32_t face, Vec3 point) const override {
        return projectSupportingSurface(face, point);
    }
    SurfaceProjection projectSupportingSurface(std::uint32_t, Vec3 point) const override {
        ++evaluations;
        return {{point.x, point.y, offset_}, {0, 0, 1}};
    }
    Vec3 projectJunction(std::span<const std::uint32_t>, Vec3) const override {
        throw std::logic_error("Chart fixture has no source junctions.");
    }
    Vec3 projectBoundary(std::span<const std::uint32_t>, Vec3) const override {
        throw std::logic_error("Chart fixture has no source curves.");
    }

  private:
    double offset_;
};
class OptionalCandidatePatch final : public SurfaceParameterization {
  public:
    explicit OptionalCandidatePatch(bool turning) : turning_(turning) {}
    const SurfaceMesh &sourceMesh() const noexcept override {
        return source_;
    }
    std::span<const ParameterTriangle> triangles() const noexcept override {
        return {};
    }
    ParameterDomain domain(std::uint32_t) const override {
        return turning_
                   ? ParameterDomain{{std::array{0., 1.}}, {std::array{0., .2}}, {}, {}, 1}
                   : ParameterDomain{{std::array{-.02, .02}}, {std::array{-.02, .02}}, {}, {}, 1};
    }
    ParametricSurfaceSample surface(std::uint32_t, ParameterPoint p) const override {
        if (turning_) {
            const double t = p.u - .5;
            return {{4 * t * t * t + .5, t * t * t - .25 * t, p.v},
                    {12 * t * t, 3 * t * t - .25, 0},
                    {0, 0, 1}};
        }
        // The natural rectangle is larger than the immutable physical halo.
        if (p.u < 0 || p.v < 0 || p.v > .25 * p.u + 1e-15 || p.v > .01 - p.u + 1e-15)
            throw std::out_of_range("Optional off-centre is outside the source halo.");
        return {{p.u, p.v, 0}, {1, 0, 0}, {0, 1, 0}};
    }
    ParametricCurveSample curve(std::uint32_t, double) const override {
        throw std::logic_error("Chart fixture has no source curves.");
    }
    CurveUseSample curveUse(std::uint32_t, double) const override {
        throw std::logic_error("Chart fixture has no source curve uses.");
    }
    std::span<const double> curvePartition(std::uint32_t) const override {
        throw std::logic_error("Chart fixture has no source curves.");
    }

  private:
    bool turning_;
    SurfaceMesh source_;
};
void optionalNativeCandidates() {
    for (const bool turning : {false, true}) {
        Fixture fixture;
        fixture.chart.faceId = 77;
        const auto provider = std::make_shared<OptionalCandidatePatch>(turning);
        fixture.chart.parameterization = provider;
        const std::array<ChartPoint, 3> uv =
            turning ? std::array<ChartPoint, 3>{{{0, 0}, {1, 0}, {.5, .2}}}
                    : std::array<ChartPoint, 3>{{{0, 0}, {.01, 0}, {.008, .002}}};
        SurfaceChartLoop loop;
        for (const auto point : uv) {
            const auto id = static_cast<std::uint32_t>(fixture.vertices.size());
            fixture.vertices.push_back(provider->surface(77, {point.u, point.v}).position);
            loop.vertices.push_back({point, id});
        }
        fixture.chart.loops.push_back(loop);
        const double target = turning ? .8 : .02;
        reference(fixture, uv, {{target, target, target}});
        for (unsigned i = 0; i < 3; ++i)
            fixture.chart.reference.front().positions[i] = fixture.vertices[i];
        MeshingOptions options{target};
        options.adaptive.enabled = false;
        bool requested = false;
        try {
            (void)meshFrontalChart(fixture.chart, fixture.vertices, options);
        } catch (const BoundaryRefinementRequired &) {
            requested = true;
        }
        require(requested, "An optional native probe hid a required physical boundary split.");
    }
}
class FiniteConePatch final : public SurfaceGeometry, public SurfaceParameterization {
  public:
    static constexpr double radius = .02, height = .004, angle = .3;
    mutable std::size_t trimmedQueries{}, supportingQueries{};
    const SurfaceMesh &sourceMesh() const noexcept override {
        return source_;
    }
    std::span<const ParameterTriangle> triangles() const noexcept override {
        return {};
    }
    ParameterDomain domain(std::uint32_t) const override {
        return {{std::array{0., angle}}, {std::array{0., height}}, {}, {}, 1};
    }
    ParametricSurfaceSample surface(std::uint32_t, ParameterPoint p) const override {
        const double r = radius + p.v, c = std::cos(p.u), s = std::sin(p.u);
        return {{r * c, r * s, p.v}, {-r * s, r * c, 0}, {c, s, 1}};
    }
    ParametricCurveSample curve(std::uint32_t, double) const override {
        throw std::logic_error("Cone chart fixture has no source curves.");
    }
    CurveUseSample curveUse(std::uint32_t, double) const override {
        throw std::logic_error("Cone chart fixture has no source curve uses.");
    }
    std::span<const double> curvePartition(std::uint32_t) const override {
        throw std::logic_error("Cone chart fixture has no source curves.");
    }
    std::span<const BoundaryCorner> boundaryCorners() const noexcept override {
        return {};
    }
    SurfaceProjection projectSurface(std::uint32_t, Vec3 p) const override {
        ++trimmedQueries;
        const double u = std::clamp(std::atan2(p.y, p.x), 0., angle);
        const double v =
            std::clamp((p.x * std::cos(u) + p.y * std::sin(u) - radius + p.z) / 2, 0., height);
        return projection(u, v);
    }
    SurfaceProjection projectSupportingSurface(std::uint32_t, Vec3 p) const override {
        ++supportingQueries;
        const double u = std::atan2(p.y, p.x), v = (std::hypot(p.x, p.y) - radius + p.z) / 2;
        if (u < 0 || u > angle || v < 0 || v > height)
            throw SurfaceProjectionDomainError("No stationary root in the finite cone domain.");
        return projection(u, v);
    }
    Vec3 projectJunction(std::span<const std::uint32_t>, Vec3) const override {
        throw std::logic_error("Cone chart fixture has no source junctions.");
    }
    Vec3 projectBoundary(std::span<const std::uint32_t>, Vec3) const override {
        throw std::logic_error("Cone chart fixture has no source curves.");
    }

  private:
    SurfaceProjection projection(double u, double v) const {
        return {surface(1, {u, v}).position,
                {std::cos(u) / std::sqrt(2.), std::sin(u) / std::sqrt(2.), -1 / std::sqrt(2.)}};
    }
    SurfaceMesh source_;
};
void finiteTrimChordDeviation() {
    const auto provider = std::make_shared<FiniteConePatch>();
    Fixture fixture;
    fixture.chart.faceId = 1;
    fixture.chart.geometry = provider;
    fixture.chart.parameterization = provider;
    const std::array<ChartPoint, 4> corners{
        {{0, 0}, {provider->angle, 0}, {provider->angle, provider->height}, {0, provider->height}}};
    SurfaceChartLoop loop;
    for (const auto uv : corners) {
        loop.vertices.push_back({uv, static_cast<std::uint32_t>(fixture.vertices.size())});
        fixture.vertices.push_back(provider->surface(1, {uv.u, uv.v}).position);
    }
    fixture.chart.loops.push_back(loop);
    for (const auto triangle : {std::array{0, 1, 2}, std::array{0, 2, 3}}) {
        ChartReferenceTriangle record;
        for (std::size_t i = 0; i < triangle.size(); ++i) {
            record.uv[i] = corners[triangle[i]];
            record.positions[i] = fixture.vertices[triangle[i]];
            record.sizes[i] = .01;
        }
        fixture.chart.reference.push_back(record);
    }
    const auto a = fixture.vertices[0], b = fixture.vertices[1];
    const Vec3 chord{std::midpoint(a.x, b.x), std::midpoint(a.y, b.y), 0};
    rejects([&] { (void)provider->projectSupportingSurface(1, chord); });
    const double measured = norm(sub(provider->projectSurface(1, chord).position, chord));
    const double expected = provider->radius * (1 - std::cos(provider->angle / 2));
    require(std::abs(measured - expected) < 1e-12,
            "The finite cone chord did not measure distance to its physical trim.");
    const auto supportingBefore = provider->supportingQueries;
    const auto trimmedBefore = provider->trimmedQueries;
    MeshingOptions options{.01};
    options.quality.enabled = false;
    options.adaptive.geometricTolerance = .0003;
    const auto mesh = meshFrontalChart(fixture.chart, fixture.vertices, options);
    require(!mesh.triangles.empty() && provider->supportingQueries == supportingBefore &&
                provider->trimmedQueries >= trimmedBefore + 4 * mesh.triangles.size(),
            "Physical chord validation required an unavailable stationary interior map.");
    // The normal component is sag/sqrt(2), below this budget; the physical
    // distance to the finite rim is above it. That distinct error must request
    // coordinated refinement, rather than first appearing as a final failure.
    options.adaptive.geometricTolerance = .00021;
    bool requested = false;
    try {
        (void)meshFrontalChart(fixture.chart, fixture.vertices, options);
    } catch (const BoundaryRefinementRequired &request) {
        requested = !request.segments().empty();
    }
    require(requested, "Physical trim distance did not trigger adaptive boundary refinement.");
}
Fixture annulus() {
    Fixture f;
    f.chart.faceId = 23;
    constexpr unsigned n = 64;
    constexpr double outer = .05, inner = .015;
    SurfaceChartLoop outside, inside;
    for (unsigned i = 0; i < n; ++i) {
        const double a = 2 * std::numbers::pi * i / n;
        loopPoint(f, outside, outer * std::cos(a), outer * std::sin(a));
    }
    for (unsigned i = 0; i < n; ++i) {
        const double a = -2 * std::numbers::pi * i / n;
        loopPoint(f, inside, inner * std::cos(a), inner * std::sin(a));
    }
    f.chart.loops = {outside, inside};
    for (unsigned i = 0; i < n; ++i) {
        const double a = 2 * std::numbers::pi * i / n, b = 2 * std::numbers::pi * (i + 1) / n;
        const ChartPoint p{outer * std::cos(a), outer * std::sin(a)},
            q{outer * std::cos(b), outer * std::sin(b)},
            r{inner * std::cos(a), inner * std::sin(a)},
            s{inner * std::cos(b), inner * std::sin(b)};
        reference(f, {{p, q, r}}, {{.008, .008, .008}});
        reference(f, {{q, s, r}}, {{.008, .008, .008}});
    }
    f.expectedArea = .5 * n * std::sin(2 * std::numbers::pi / n) * (outer * outer - inner * inner);
    return f;
}
Fixture obtuseStrip() {
    Fixture f;
    f.chart.faceId = 34;
    f.chart.planar = true;
    f.chart.planeU = {1, 0, 0};
    f.chart.planeV = {0, 1, 0};
    SurfaceChartLoop loop;
    for (double x : {.45, .46240911187315553, .483, .49})
        loopPoint(f, loop, x, -.35);
    for (double x : {.49, .47699355798661924, .46874894906677228, .46, .45})
        loopPoint(f, loop, x, -.342);
    f.chart.loops.push_back(std::move(loop));
    reference(f, {{{.45, -.35}, {.49, -.35}, {.49, -.342}}}, {{.04, .04, .04}});
    reference(f, {{{.45, -.35}, {.49, -.342}, {.45, -.342}}}, {{.04, .04, .04}});
    f.expectedArea = .04 * .008;
    return f;
}
Fixture mirroredStripNeighbor(Fixture &first) {
    auto second = first;
    second.chart.faceId = 35;
    auto &loop = second.chart.loops.front().vertices;
    for (auto &v : loop) {
        if (v.atlasVertex < 4)
            continue;
        v.uv.v = -.7 - v.uv.v;
        v.atlasVertex = static_cast<std::uint32_t>(first.vertices.size());
        first.vertices.push_back({v.uv.u, v.uv.v, 0});
    }
    std::reverse(loop.begin(), loop.end());
    for (auto &t : second.chart.reference) {
        for (int i = 0; i < 3; ++i) {
            t.uv[i].v = -.7 - t.uv[i].v;
            t.positions[i].y = t.uv[i].v;
        }
        std::swap(t.uv[1], t.uv[2]);
        std::swap(t.positions[1], t.positions[2]);
    }
    second.vertices = first.vertices;
    return second;
}
void splitStripBoundary(Fixture &first, Fixture &second, AtlasSegment segment) {
    const auto start = first.vertices[segment[0]], end = first.vertices[segment[1]];
    const Vec3 middle{std::midpoint(start.x, end.x), std::midpoint(start.y, end.y), 0};
    const auto id = static_cast<std::uint32_t>(first.vertices.size());
    first.vertices.push_back(middle);
    second.vertices = first.vertices;
    for (auto *fixture : {&first, &second}) {
        bool found = false;
        for (auto &loop : fixture->chart.loops)
            for (std::size_t i = 0; i < loop.vertices.size(); ++i) {
                const auto a = loop.vertices[i].atlasVertex;
                const auto b = loop.vertices[(i + 1) % loop.vertices.size()].atlasVertex;
                if (AtlasSegment{std::min(a, b), std::max(a, b)} != segment)
                    continue;
                loop.vertices.insert(loop.vertices.begin() + i + 1, {{middle.x, middle.y}, id});
                found = true;
                break;
            }
        require(found, "Shared refinement did not reach both neighboring strip boundaries.");
    }
}

struct Metrics {
    double area{}, angle{180}, aspect{}, leftEdge{}, rightEdge{};
    size_t leftCount{}, rightCount{};
};
Metrics inspect(const Fixture &f, const FrontalChartMesh &m) {
    Metrics q;
    std::map<std::pair<unsigned, unsigned>, unsigned> edges;
    std::set<unsigned> atlas;
    for (const auto &v : m.vertices) {
        require(std::isfinite(v.position.x) && std::isfinite(v.position.y) &&
                    std::isfinite(v.position.z),
                "Non-finite frontal vertex.");
        if (v.atlasVertex != std::numeric_limits<unsigned>::max()) {
            atlas.insert(v.atlasVertex);
            const auto p = f.vertices[v.atlasVertex];
            require(p.x == v.position.x && p.y == v.position.y && p.z == v.position.z,
                    "A shared atlas vertex moved.");
        }
    }
    for (const auto &t : m.triangles) {
        std::array<Vec3, 3> p;
        for (int i = 0; i < 3; ++i) {
            require(t[i] < m.vertices.size(), "Invalid triangle index.");
            p[i] = m.vertices[t[i]].position;
        }
        const auto n = cross(sub(p[1], p[0]), sub(p[2], p[0]));
        const double area2 = norm(n);
        require(n.z > 0, "Frontal triangle orientation changed.");
        q.area += area2 / 2;
        double longest = 0;
        for (int i = 0; i < 3; ++i) {
            const auto a = sub(p[(i + 1) % 3], p[i]), b = sub(p[(i + 2) % 3], p[i]);
            const double length = norm(a);
            longest = std::max(longest, length);
            q.angle = std::min(q.angle, std::atan2(area2, dot(a, b)) * 180 / std::numbers::pi);
            ++edges[std::minmax(t[i], t[(i + 1) % 3])];
            const auto ua = m.vertices[t[i]].uv, ub = m.vertices[t[(i + 1) % 3]].uv;
            require(
                length <=
                    4. / 3 *
                        f.chart.sample({(ua.u + ub.u) / 2, (ua.v + ub.v) / 2}).targetEdgeLength *
                        (1 + 1e-12),
                "An output edge exceeds its scalar target.");
            const double x = (p[i].x + p[(i + 1) % 3].x) / 2;
            if (x < .015) {
                q.leftEdge += length;
                ++q.leftCount;
            }
            if (x > .06) {
                q.rightEdge += length;
                ++q.rightCount;
            }
        }
        q.aspect = std::max(q.aspect, longest * longest / area2);
    }
    size_t boundary = 0;
    for (auto [e, count] : edges) {
        require(count == 1 || count == 2, "Non-manifold frontal chart edge.");
        boundary += count == 1;
    }
    size_t expectedBoundary = 0;
    for (auto &loop : f.chart.loops) {
        expectedBoundary += loop.vertices.size();
        for (auto &v : loop.vertices)
            require(atlas.contains(v.atlasVertex), "A shared boundary point was dropped.");
    }
    for (auto &v : f.chart.fixedPoints)
        require(atlas.contains(v.atlasVertex), "A fixed chart interior point was dropped.");
    require(boundary == expectedBoundary, "The frontal mesh has a crack or filled a hole.");
    require(std::abs(q.area - f.expectedArea) < 1e-11 * f.expectedArea,
            "The physical chart area changed.");
    require(q.angle >= 20 && q.aspect <= 4,
            "The frontal chart violates hard 20-degree/aspect-four gates.");
    return q;
}
void variableSpeedPlane() {
    for (const auto [strength, graded] : std::array<std::pair<double, bool>, 4>{
             {{1., false}, {2., false}, {3., false}, {3., true}}}) {
        auto fixture = plate(graded);
        auto provider = std::make_shared<VariableSpeedPlane>(strength);
        fixture.chart.parameterization = provider;
        auto geometry = std::make_shared<PlaneDeviationGeometry>();
        fixture.chart.geometry = geometry;
        for (auto &loop : fixture.chart.loops)
            for (auto &vertex : loop.vertices)
                vertex.uv.u = provider->parameter(vertex.uv.u);
        for (auto &triangle : fixture.chart.reference)
            for (auto &point : triangle.uv)
                point.u = provider->parameter(point.u);
        MeshingOptions options{.01};
        // This plane has exactly zero geometric approximation error. Its UV
        // midpoint can differ from its chord midpoint by millimetres; that
        // tangential speed difference must not lock physical quality repair.
        options.adaptive.geometricTolerance = 1e-8;
        FrontalChartMesh mesh;
        for (;;) {
            try {
                mesh = meshFrontalChart(fixture.chart, fixture.vertices, options);
                break;
            } catch (const BoundaryRefinementRequired &request) {
                require(fixture.vertices.size() < 10000,
                        "Variable-speed plane exhausted shared refinement.");
                for (const auto segment : request.segments()) {
                    auto &loop = fixture.chart.loops.front().vertices;
                    for (std::size_t i = 0; i < loop.size(); ++i) {
                        const auto a = loop[i].atlasVertex,
                                   b = loop[(i + 1) % loop.size()].atlasVertex;
                        if (AtlasSegment{std::min(a, b), std::max(a, b)} != segment)
                            continue;
                        const auto p = fixture.vertices[a], q = fixture.vertices[b];
                        const Vec3 middle{std::midpoint(p.x, q.x), std::midpoint(p.y, q.y), 0};
                        const auto id = static_cast<std::uint32_t>(fixture.vertices.size());
                        fixture.vertices.push_back(middle);
                        loop.insert(loop.begin() + static_cast<std::ptrdiff_t>(i + 1),
                                    {{provider->parameter(middle.x), middle.y}, id});
                        break;
                    }
                }
            }
        }
        const auto q = inspect(fixture, mesh);
        for (const auto &vertex : mesh.vertices)
            require(vertex.position.z == 0, "Planar quality repair left its source surface.");
        require(geometry->evaluations >= 4 * mesh.triangles.size(),
                "Final validation did not project physical chord midpoints and cell centroids.");
        auto certifiedGeometry = std::make_shared<PlaneDeviationGeometry>();
        fixture.chart.geometry = certifiedGeometry;
        fixture.chart.parameterization = std::make_shared<VariableSpeedPlane>(strength, true);
        const auto certified = meshFrontalChart(fixture.chart, fixture.vertices, options);
        require(certified.triangles == mesh.triangles &&
                    certified.vertices.size() == mesh.vertices.size(),
                "Certified native distance witnesses changed variable-speed plane topology.");
        for (std::size_t i = 0; i < mesh.vertices.size(); ++i)
            require(certified.vertices[i].position == mesh.vertices[i].position,
                    "Certified native distance witnesses changed physical mesh positions.");
        require(certifiedGeometry->evaluations < geometry->evaluations / 2,
                "Variable-speed native witnesses did not avoid redundant global projections.");
        fixture.chart.parameterization = provider;
        // Exact projection remains authoritative even when the native position
        // map reports zero normal deviation. Inconsistent providers must fail.
        fixture.chart.geometry = std::make_shared<PlaneDeviationGeometry>(.001);
        bool rejectedDeviation = false;
        try {
            (void)meshFrontalChart(fixture.chart, fixture.vertices, options);
        } catch (const std::runtime_error &error) {
            rejectedDeviation =
                std::string(error.what()).find("sampled physical deviation") != std::string::npos;
        }
        require(rejectedDeviation, "Exact physical deviation did not reject an inaccurate mesh.");
        std::cout << "variable-speed plane " << strength << (graded ? " graded" : "") << ": "
                  << mesh.triangles.size() << " triangles, " << q.angle << " deg / " << q.aspect
                  << " aspect\n";
    }
}
void sameMesh(const FrontalChartMesh &a, const FrontalChartMesh &b) {
    require(a.vertices.size() == b.vertices.size() && a.triangles == b.triangles,
            "Identical frontal inputs produced different topology or output ordering.");
    for (size_t i = 0; i < a.vertices.size(); ++i) {
        const auto &p = a.vertices[i];
        const auto &q = b.vertices[i];
        require(p.uv.u == q.uv.u && p.uv.v == q.uv.v && p.position.x == q.position.x &&
                    p.position.y == q.position.y && p.position.z == q.position.z &&
                    p.atlasVertex == q.atlasVertex,
                "Identical frontal inputs produced different vertex positions.");
    }
}
SurfaceMesh revolvedFixture(unsigned n, bool tapered) {
    SurfaceMesh source;
    constexpr unsigned levels = 4;
    const double height = tapered ? .006 : .06;
    for (unsigned z = 0; z <= levels; ++z)
        for (unsigned i = 0; i < n; ++i) {
            const double phi = 2 * std::numbers::pi * i / n;
            const double radius = tapered ? .014 - .01 * z / levels : .02;
            source.vertices.push_back(
                {radius * std::cos(phi), radius * std::sin(phi), height * z / levels});
        }
    for (unsigned z = 0; z < levels; ++z)
        for (unsigned i = 0; i < n; ++i) {
            const auto a = z * n + i, b = z * n + (i + 1) % n, c = (z + 1) * n + i,
                       d = (z + 1) * n + (i + 1) % n;
            source.triangles.push_back({{a, b, c}, 1, BoundaryCondition::wall});
            source.triangles.push_back({{b, d, c}, 1, BoundaryCondition::wall});
        }
    const auto bottom = static_cast<unsigned>(source.vertices.size());
    source.vertices.push_back({0, 0, 0});
    const auto top = static_cast<unsigned>(source.vertices.size());
    source.vertices.push_back({0, 0, height});
    for (unsigned i = 0; i < n; ++i) {
        source.triangles.push_back({{bottom, (i + 1) % n, i}, 2, BoundaryCondition::driven});
        source.triangles.push_back(
            {{top, levels * n + i, levels * n + (i + 1) % n}, 3, BoundaryCondition::wall});
    }
    return source;
}
void revolvedAtlas(bool tapered = false) {
    const auto source = revolvedFixture(tapered ? 48 : 32, tapered);
    MeshingOptions options;
    options.targetEdgeLength = tapered ? .004 : .008;
    options.adaptive.enabled = false;
    std::stop_source timeout;
    std::jthread deadline([&](std::stop_token finished) {
        const auto end = std::chrono::steady_clock::now() + std::chrono::seconds(15);
        while (!finished.stop_requested() && std::chrono::steady_clock::now() < end)
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        if (!finished.stop_requested())
            timeout.request_stop();
    });
    const auto begin = std::chrono::steady_clock::now();
    const auto atlas = makeSurfaceAtlas(source, options, timeout.get_token());
    SurfaceMesh output;
    output.vertices = atlas.vertices;
    output.triangles = atlas.fixedTriangles;
    for (const auto &chart : atlas.charts) {
        std::cout << (tapered ? "frustum face " : "cylinder face ") << chart.faceId << " boundary "
                  << chart.loops.front().vertices.size() << std::endl;
        const auto mesh = meshFrontalChart(chart, atlas.vertices, options, timeout.get_token());
        std::vector<unsigned> mapping;
        for (const auto &v : mesh.vertices) {
            if (v.atlasVertex != std::numeric_limits<unsigned>::max())
                mapping.push_back(v.atlasVertex);
            else {
                mapping.push_back(static_cast<unsigned>(output.vertices.size()));
                output.vertices.push_back(v.position);
            }
        }
        for (const auto &t : mesh.triangles)
            output.triangles.push_back(
                {{mapping[t[0]], mapping[t[1]], mapping[t[2]]}, chart.faceId, chart.boundary});
    }
    const auto quality = inspectMesh(output);
    require(quality.watertight() && quality.eulerCharacteristic == 2 &&
                quality.minimumAngleDegrees >= 20 && quality.maximumAspectRatio <= 4,
            "The genuine revolved atlas failed final closed-surface quality.");
    std::cout << (tapered ? "frustum atlas: " : "cylinder atlas: ") << output.triangles.size()
              << " triangles, " << quality.minimumAngleDegrees << " deg / "
              << quality.maximumAspectRatio << " aspect, "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count()
              << " s\n";
}
void fixedPhysicalBoundaryGrading() {
    constexpr unsigned count = 128;
    auto source = revolvedFixture(count, false);
    for (unsigned i = 0; i < count; ++i)
        source.fixedVertices.push_back(i);
    require(inspectMesh(source).watertight(), "The fixed-rim fixture is not closed.");
    MeshingOptions options{.008};
    options.adaptive.enabled = false;
    options.quality.enabled = false;
    const auto atlas = makeSurfaceAtlas(source, options);
    const auto found = std::find_if(atlas.charts.begin(), atlas.charts.end(),
                                    [](const auto &chart) { return chart.faceId == 2; });
    require(found != atlas.charts.end(), "The fixed-rim cap is missing.");
    const auto &chart = *found;
    const double rimLength = norm(sub(source.vertices[1], source.vertices[0]));
    std::set<std::array<unsigned, 2>> boundary;
    ChartPoint center{};
    for (const auto &loop : chart.loops)
        for (std::size_t i = 0; i < loop.vertices.size(); ++i) {
            const auto a = loop.vertices[i], b = loop.vertices[(i + 1) % loop.vertices.size()];
            boundary.insert(
                {std::min(a.atlasVertex, b.atlasVertex), std::max(a.atlasVertex, b.atlasVertex)});
            require(chart.sample(a.uv).targetEdgeLength <= rimLength * (1 + 1e-10),
                    "The scalar field ignores mandatory physical boundary spacing.");
            center.u += a.uv.u / count;
            center.v += a.uv.v / count;
        }
    require(chart.sample(center).targetEdgeLength > 3 * rimLength,
            "Mandatory boundary spacing became uniform interior refinement.");
    const auto cap = meshFrontalChart(chart, atlas.vertices, options);
    double maximumRimFanEdge = 0, maximumEdge = 0;
    for (const auto &triangle : cap.triangles) {
        bool touchesRim = false;
        double longest = 0;
        for (std::size_t i = 0; i < 3; ++i) {
            const auto &a = cap.vertices[triangle[i]];
            const auto &b = cap.vertices[triangle[(i + 1) % 3]];
            touchesRim = touchesRim || boundary.contains({std::min(a.atlasVertex, b.atlasVertex),
                                                          std::max(a.atlasVertex, b.atlasVertex)});
            longest = std::max(longest, norm(sub(a.position, b.position)));
        }
        maximumEdge = std::max(maximumEdge, longest);
        if (touchesRim)
            maximumRimFanEdge = std::max(maximumRimFanEdge, longest);
    }
    require(maximumRimFanEdge < 4 * rimLength && maximumEdge > 2.5 * rimLength,
            "Fixed physical spacing did not grade from the rim into the cap interior.");
    require(cap.triangles.size() > count && cap.triangles.size() < 12 * count,
            "Fixed-rim sizing produced a star fan or unnecessary uniform fine mesh.");
    require(atlas.fixedVertices.size() == count,
            "Boundary grading changed mandatory physical point metadata.");
    std::cout << "fixed physical rim, shape limits off: " << cap.triangles.size()
              << " cap triangles, rim-adjacent edge/rim " << maximumRimFanEdge / rimLength
              << ", interior edge/rim " << maximumEdge / rimLength << '\n';
}
void sphericalBoundaryValence() {
    // This closed curved fixture previously left a 14.7-degree triangle beside
    // a fixed shared seam. Encroachment prevention and bounded physical swaps
    // must meet the same quality gates without changing a single seam position.
    auto source = testing::sphere(3);
    for (std::size_t i = 0; i < source.triangles.size(); ++i) {
        source.triangles[i].faceId = static_cast<std::uint32_t>(i / 64 + 1);
        source.triangles[i].boundary = i < 64 ? BoundaryCondition::driven : BoundaryCondition::wall;
    }
    MeshingOptions options{.025};
    options.adaptive.geometricTolerance = .005;
    const auto atlas = makeSurfaceAtlas(source, options);
    SurfaceMesh output;
    output.vertices = atlas.vertices;
    output.triangles = atlas.fixedTriangles;
    std::size_t swaps = 0;
    for (const auto &chart : atlas.charts) {
        const auto mesh = meshFrontalChart(chart, atlas.vertices, options);
        swaps += mesh.qualityFlips;
        std::vector<std::uint32_t> mapping;
        for (const auto &v : mesh.vertices) {
            if (v.atlasVertex != std::numeric_limits<std::uint32_t>::max()) {
                mapping.push_back(v.atlasVertex);
                const auto p = atlas.vertices[v.atlasVertex];
                require(p.x == v.position.x && p.y == v.position.y && p.z == v.position.z,
                        "A spherical boundary point moved during physical quality optimization.");
            } else {
                mapping.push_back(static_cast<std::uint32_t>(output.vertices.size()));
                output.vertices.push_back(v.position);
            }
        }
        for (const auto &t : mesh.triangles)
            output.triangles.push_back(
                {{mapping[t[0]], mapping[t[1]], mapping[t[2]]}, chart.faceId, chart.boundary});
    }
    const auto q = inspectMesh(output);
    require(q.watertight() && q.eulerCharacteristic == 2 && q.signedVolume > 0 &&
                q.minimumAngleDegrees >= 20 && q.maximumAspectRatio <= 4,
            "The spherical shared-boundary valence regression violates closed-surface quality.");
    std::set<std::uint32_t> tags;
    for (const auto &t : output.triangles) {
        tags.insert(t.faceId);
        require((t.faceId == 1) == (t.boundary == BoundaryCondition::driven),
                "Spherical boundary assignments changed.");
    }
    require(tags.size() == 8, "Spherical physical patches were lost.");
    std::cout << "spherical boundary: " << output.triangles.size() << " triangles, "
              << q.minimumAngleDegrees << " deg / " << q.maximumAspectRatio << " aspect, " << swaps
              << " physical quality swaps\n";
}
void outsideBoundaryRefinement() {
    auto first = obtuseStrip();
    auto second = mirroredStripNeighbor(first);
    MeshingOptions options;
    options.targetEdgeLength = .04;
    // The obtuse cell has a legal minimum angle but violates aspect 4. Its
    // circumcentre lies outside the domain and encroaches a different boundary.
    for (const auto *fixture : {&first, &second}) {
        bool requested = false;
        try {
            (void)meshFrontalChart(fixture->chart, fixture->vertices, options);
        } catch (const BoundaryRefinementRequired &error) {
            requested = true;
            require(error.segments().size() == 1 && error.segments()[0] == AtlasSegment{1, 2},
                    "Outside candidate requested an unrelated or unnecessary boundary split.");
        }
        require(requested, "Outside circumcentre hid its crossed boundary encroachment.");
    }
    splitStripBoundary(first, second, {1, 2});
    const auto firstMesh = meshFrontalChart(first.chart, first.vertices, options);
    const auto secondMesh = meshFrontalChart(second.chart, second.vertices, options);
    inspect(first, firstMesh);
    inspect(second, secondMesh);
    sameMesh(firstMesh, meshFrontalChart(first.chart, first.vertices, options));
    const auto middle = static_cast<std::uint32_t>(first.vertices.size() - 1);
    std::map<AtlasSegment, std::pair<unsigned, int>> uses;
    for (const auto *mesh : {&firstMesh, &secondMesh})
        for (const auto &triangle : mesh->triangles)
            for (int i = 0; i < 3; ++i) {
                const auto a = mesh->vertices[triangle[i]].atlasVertex;
                const auto b = mesh->vertices[triangle[(i + 1) % 3]].atlasVertex;
                const AtlasSegment key{std::min(a, b), std::max(a, b)};
                if (key != AtlasSegment{1, middle} && key != AtlasSegment{2, middle})
                    continue;
                ++uses[key].first;
                uses[key].second += a < b ? 1 : -1;
            }
    require(uses.size() == 2, "Shared child edges were lost after boundary refinement.");
    for (const auto &[key, use] : uses) {
        (void)key;
        require(use.first == 2 && use.second == 0,
                "Refined neighboring boundaries are not conforming with opposite orientation.");
    }
    options.quality.enabled = false;
    const auto unrestricted = obtuseStrip();
    (void)meshFrontalChart(unrestricted.chart, unrestricted.vertices, options);
}
} // namespace
int main() {
    try {
        MeshingOptions options;
        options.targetEdgeLength = .01;
        auto uniform = plate();
        uniform.chart.fixedPoints.push_back(
            {{.05, .04}, static_cast<unsigned>(uniform.vertices.size())});
        uniform.vertices.push_back({.05, .04, 0});
        const auto first = meshFrontalChart(uniform.chart, uniform.vertices, options);
        const auto q = inspect(uniform, first);
        for (int repeat = 0; repeat < 3; ++repeat)
            sameMesh(first, meshFrontalChart(uniform.chart, uniform.vertices, options));
        require(first.insertedInteriorVertices > 0,
                "Frontal test did not perform real interior insertion.");
        std::cout << "plate: " << first.triangles.size() << " triangles, " << q.angle << " deg / "
                  << q.aspect << " aspect\n";
        auto stretched = uniform;
        for (auto &loop : stretched.chart.loops)
            for (auto &v : loop.vertices) {
                v.uv.u *= 1000;
                v.uv.v *= .001;
            }
        for (auto &v : stretched.chart.fixedPoints) {
            v.uv.u *= 1000;
            v.uv.v *= .001;
        }
        for (auto &t : stretched.chart.reference)
            for (auto &v : t.uv) {
                v.u *= 1000;
                v.v *= .001;
            }
        const auto transformed = meshFrontalChart(stretched.chart, stretched.vertices, options);
        inspect(stretched, transformed);
        require(transformed.triangles.size() < 2 * first.triangles.size(),
                "Affine UV units caused excessive refinement.");
        const auto ring = annulus();
        const auto ringMesh = meshFrontalChart(ring.chart, ring.vertices, options);
        const auto ringQuality = inspect(ring, ringMesh);
        std::cout << "annulus: " << ringMesh.triangles.size() << " triangles, " << ringQuality.angle
                  << " deg / " << ringQuality.aspect << " aspect\n";
        const auto graded = plate(true);
        const auto gradedMesh = meshFrontalChart(graded.chart, graded.vertices, options);
        const auto gradedQuality = inspect(graded, gradedMesh);
        require(gradedQuality.leftCount && gradedQuality.rightCount &&
                    gradedQuality.leftEdge / gradedQuality.leftCount <
                        .6 * gradedQuality.rightEdge / gradedQuality.rightCount,
                "The scalar size gradient did not create local refinement.");
        std::cout << "graded: " << gradedMesh.triangles.size() << " triangles, left/right edge "
                  << gradedQuality.leftEdge / gradedQuality.leftCount << " / "
                  << gradedQuality.rightEdge / gradedQuality.rightCount << " m\n";
        std::stop_source stop;
        stop.request_stop();
        rejects(
            [&] { meshFrontalChart(uniform.chart, uniform.vertices, options, stop.get_token()); });
        auto invalid = uniform;
        invalid.chart.loops[0].vertices[0].atlasVertex = 999999;
        rejects([&] { meshFrontalChart(invalid.chart, invalid.vertices, options); });
        auto coarse = uniform;
        for (auto &t : coarse.chart.reference)
            t.sizes = {.001, .001, .001};
        bool requested = false;
        try {
            (void)meshFrontalChart(coarse.chart, coarse.vertices, options);
        } catch (const BoundaryRefinementRequired &error) {
            requested = true;
            require(error.segments().size() == 1 && error.segments()[0][0] == 0 &&
                        error.segments()[0][1] == 1,
                    "A coarse shared boundary did not request its exact canonical segment.");
        }
        require(requested,
                "An oversized shared boundary was not distinguished from a numerical failure.");
        // Every circumcentre of this obtuse boundary ear lies outside its
        // domain. The existing opposite vertex still encroaches segment 0--1.
        Fixture ear;
        ear.chart.faceId = 31;
        SurfaceChartLoop earLoop;
        loopPoint(ear, earLoop, 0, 0);
        loopPoint(ear, earLoop, .01, 0);
        loopPoint(ear, earLoop, .008, .002);
        ear.chart.loops.push_back(earLoop);
        reference(ear, {{{0, 0}, {.01, 0}, {.008, .002}}}, {{.02, .02, .02}});
        requested = false;
        try {
            (void)meshFrontalChart(ear.chart, ear.vertices, options);
        } catch (const BoundaryRefinementRequired &error) {
            requested = true;
            require(error.segments().size() == 1 && error.segments()[0] == AtlasSegment{0, 1},
                    "An existing boundary-ear vertex did not identify its encroached segment.");
        }
        require(requested, "An outside circumcentre hid existing boundary encroachment.");
        options.quality.enabled = false;
        require(meshFrontalChart(ear.chart, ear.vertices, options).triangles.size() == 1,
                "Disabled shape limits still refine a valid small boundary ear.");
        options.quality = {true, 5, 10};
        require(meshFrontalChart(ear.chart, ear.vertices, options).triangles.size() == 1,
                "The front ignores custom angle/aspect limits.");
        options.quality = {};
        variableSpeedPlane();
        optionalNativeCandidates();
        finiteTrimChordDeviation();
        revolvedAtlas();
        revolvedAtlas(true);
        fixedPhysicalBoundaryGrading();
        sphericalBoundaryValence();
        outsideBoundaryRefinement();
        std::cout
            << "Frontal-Delaunay domains, holes, metric, sizing and cancellation tests passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
