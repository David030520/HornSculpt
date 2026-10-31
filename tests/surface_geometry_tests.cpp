#include "../modules/meshing/src/surface_charts.hpp"
#include <hornsculpt/meshing.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>

using namespace hornsculpt;
using namespace hornsculpt::detail;
namespace {
constexpr double radius = .01, height = .02;
void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
template <typename F> void rejects(F operation) {
    bool rejected = false;
    try {
        operation();
    } catch (const std::exception &) {
        rejected = true;
    }
    require(rejected, "Invalid exact geometry was silently accepted.");
}
class CylinderGeometry final : public SurfaceGeometry {
  public:
    std::span<const BoundaryCorner> boundaryCorners() const noexcept override {
        return {};
    }
    mutable std::atomic<std::size_t> surfaceCalls{};
    SurfaceProjection projectSurface(std::uint32_t face, Vec3 p) const override {
        ++surfaceCalls;
        if (face == 1) {
            const double r = std::hypot(p.x, p.y);
            if (!(r > 0))
                throw std::runtime_error("Cylinder axis has no unique side projection.");
            return {{radius * p.x / r, radius * p.y / r, std::clamp(p.z, 0., height)},
                    {p.x / r, p.y / r, 0}};
        }
        if (face == 2 || face == 3) {
            const double r = std::hypot(p.x, p.y), scale = r > radius ? radius / r : 1.;
            return {{scale * p.x, scale * p.y, face == 2 ? 0. : height},
                    {0, 0, face == 2 ? -1. : 1.}};
        }
        throw std::invalid_argument("Unknown analytic face.");
    }
    SurfaceProjection projectSupportingSurface(std::uint32_t face, Vec3 p) const override {
        ++surfaceCalls;
        if (face == 1) {
            const double r = std::hypot(p.x, p.y);
            if (!(r > 0))
                throw std::runtime_error("Cylinder axis has no unique supporting projection.");
            return {{radius * p.x / r, radius * p.y / r, p.z}, {p.x / r, p.y / r, 0}};
        }
        if (face == 2 || face == 3)
            return {{p.x, p.y, face == 2 ? 0. : height}, {0, 0, face == 2 ? -1. : 1.}};
        throw std::invalid_argument("Unknown analytic supporting face.");
    }
    Vec3 projectJunction(std::span<const std::uint32_t>, Vec3) const override {
        throw std::runtime_error("This analytic source has no physical triple junction.");
    }
    Vec3 projectBoundary(std::span<const std::uint32_t> faces, Vec3 p) const override {
        require(faces.size() == 2 && faces[0] == 1 && (faces[1] == 2 || faces[1] == 3),
                "Artificial seams were sent to the common-curve evaluator.");
        auto result = projectSurface(1, p).position;
        result.z = faces[1] == 2 ? 0 : height;
        return result;
    }
};
SurfaceMesh cylinder() {
    constexpr unsigned n = 24;
    SurfaceMesh mesh;
    for (unsigned layer = 0; layer < 2; ++layer)
        for (unsigned i = 0; i < n; ++i) {
            const double a = 2 * std::numbers::pi * i / n;
            mesh.vertices.push_back({radius * std::cos(a), radius * std::sin(a), layer * height});
        }
    mesh.vertices.push_back({0, 0, 0});
    mesh.vertices.push_back({0, 0, height});
    for (unsigned i = 0; i < n; ++i) {
        const unsigned j = (i + 1) % n;
        mesh.triangles.push_back({{i, j, j + n}, 1, BoundaryCondition::wall});
        mesh.triangles.push_back({{i, j + n, i + n}, 1, BoundaryCondition::wall});
        mesh.triangles.push_back({{2 * n, j, i}, 2, BoundaryCondition::driven});
        mesh.triangles.push_back({{2 * n + 1, i + n, j + n}, 3, BoundaryCondition::wall});
    }
    mesh.fixedVertices = {0};
    return mesh;
}
MeshingOptions options() {
    MeshingOptions result;
    result.targetEdgeLength = .006;
    result.adaptive.enabled = false;
    result.geometry = std::make_shared<CylinderGeometry>();
    return result;
}
void exactAtlas() {
    const auto source = cylinder();
    auto atlas = makeSurfaceAtlas(source, options());
    bool lifted = false, artificial = false;
    std::array<std::uint32_t, 2> physical{}, artificialSegment{};
    for (const auto &chart : atlas.charts) {
        for (const auto &t : chart.reference) {
            const ChartPoint uv{(t.uv[0].u + t.uv[1].u + t.uv[2].u) / 3,
                                (t.uv[0].v + t.uv[1].v + t.uv[2].v) / 3};
            const auto sample = chart.sample(uv);
            if (chart.faceId == 1) {
                require(std::abs(std::hypot(sample.position.x, sample.position.y) - radius) < 1e-14,
                        "Interior samples remained on the polygonal cylinder.");
                require(std::abs(sample.du.x * sample.position.x +
                                 sample.du.y * sample.position.y) < 1e-14 &&
                            std::abs(sample.dv.x * sample.position.x +
                                     sample.dv.y * sample.position.y) < 1e-14,
                        "Chart differential is not tangent to the exact cylinder.");
                lifted = true;
            }
        }
        for (const auto &loop : chart.loops)
            for (std::size_t i = 0; i < loop.vertices.size(); ++i) {
                const auto a = loop.vertices[i].atlasVertex;
                const auto b = loop.vertices[(i + 1) % loop.vertices.size()].atlasVertex;
                const auto p = atlas.vertices[a], q = atlas.vertices[b];
                if (chart.faceId == 2)
                    physical = {a, b};
                if (chart.faceId == 1 && std::abs(p.z - q.z) > .001) {
                    artificial = true;
                    artificialSegment = {a, b};
                }
            }
    }
    require(lifted && artificial && physical[0] != physical[1],
            "Exact atlas fixture was incomplete.");
    const auto newId = atlas.vertices.size();
    refineSurfaceAtlas(atlas, std::span{&physical, 1});
    const auto p = atlas.vertices[newId];
    require(std::abs(std::hypot(p.x, p.y) - radius) < 1e-14 && std::abs(p.z) < 1e-14,
            "Shared refinement inserted a polygonal chord midpoint.");
    unsigned incidences = 0;
    for (const auto &chart : atlas.charts)
        for (const auto &loop : chart.loops)
            for (const auto &vertex : loop.vertices)
                incidences += vertex.atlasVertex == newId;
    require(incidences == 2, "Exact curve refinement is not conforming.");
    require(atlas.vertices[atlas.fixedVertices.front()] == source.vertices.front(),
            "Exact projection moved a fixed source point.");
    const auto seamId = atlas.vertices.size();
    refineSurfaceAtlas(atlas, std::span{&artificialSegment, 1});
    const auto seam = atlas.vertices[seamId];
    require(std::abs(std::hypot(seam.x, seam.y) - radius) < 1e-14,
            "Artificial cut refinement left the exact physical surface.");
}
void exactMesh() {
    const auto source = cylinder();
    const auto result = remeshIsotropic(source, options());
    require(result.quality.watertight() && result.quality.minimumAngleDegrees >= 20 &&
                result.quality.maximumAspectRatio <= 4,
            "Exact geometric lifting failed the normal output quality gate.");
    std::set<std::uint32_t> tags;
    for (const auto &t : result.mesh.triangles) {
        tags.insert(t.faceId);
        for (auto index : t.vertices) {
            const auto p = result.mesh.vertices[index];
            if (t.faceId == 1)
                require(std::abs(std::hypot(p.x, p.y) - radius) < 1e-14,
                        "Output side vertex is not on the exact cylinder.");
            else
                require(std::abs(p.z - (t.faceId == 2 ? 0 : height)) < 1e-14,
                        "Output cap vertex left its exact plane.");
        }
        require((t.faceId == 2) == (t.boundary == BoundaryCondition::driven),
                "Exact meshing changed boundary assignments.");
    }
    require(tags == std::set<std::uint32_t>{1, 2, 3}, "Exact meshing lost physical faces.");
    std::cout << "exact cylinder: " << result.mesh.triangles.size() << " triangles, "
              << result.quality.minimumAngleDegrees << " deg / "
              << result.quality.maximumAspectRatio << " aspect\n";
}
void sizingDoesNotProject() {
    const auto geometry = std::make_shared<CylinderGeometry>();
    auto settings = options();
    settings.geometry = geometry;
    const auto atlas = makeSurfaceAtlas(cylinder(), settings);
    geometry->surfaceCalls = 0;
    for (const auto &chart : atlas.charts) {
        auto referenceOnly = chart;
        referenceOnly.geometry.reset();
        for (const auto &t : chart.reference) {
            const ChartPoint uv{(t.uv[0].u + t.uv[1].u + t.uv[2].u) / 3,
                                (t.uv[0].v + t.uv[1].v + t.uv[2].v) / 3};
            require(chart.targetSize(uv) == referenceOnly.sample(uv).targetEdgeLength,
                    "Size-only evaluation changed the source scalar field.");
        }
    }
    require(geometry->surfaceCalls == 0,
            "Size-only evaluation invoked the exact geometry provider.");
    const auto &chart = atlas.charts.front();
    const auto &t = chart.reference.front();
    (void)chart.sample(
        {(t.uv[0].u + t.uv[1].u + t.uv[2].u) / 3, (t.uv[0].v + t.uv[1].v + t.uv[2].v) / 3});
    require(geometry->surfaceCalls == 1, "Physical sampling bypassed the exact geometry provider.");
}
class SphereGeometry final : public SurfaceGeometry {
  public:
    std::span<const BoundaryCorner> boundaryCorners() const noexcept override {
        return {};
    }
    SurfaceProjection projectSurface(std::uint32_t, Vec3 p) const override {
        const double r = std::hypot(p.x, p.y, p.z);
        if (!(r > 0))
            throw std::runtime_error("Sphere centre has no unique projection.");
        return {{.03 * p.x / r, .03 * p.y / r, .03 * p.z / r}, {p.x / r, p.y / r, p.z / r}};
    }
    SurfaceProjection projectSupportingSurface(std::uint32_t face, Vec3 p) const override {
        return projectSurface(face, p);
    }
    Vec3 projectJunction(std::span<const std::uint32_t>, Vec3) const override {
        throw std::runtime_error("This analytic source has no physical triple junction.");
    }
    Vec3 projectBoundary(std::span<const std::uint32_t>, Vec3) const override {
        throw std::runtime_error("The spherical patch has no physical face junction.");
    }
};
SurfaceMesh sphericalPatch() {
    SurfaceMesh source;
    constexpr unsigned n = 8;
    for (unsigned row = 0; row <= n; ++row)
        for (unsigned column = 0; column <= n; ++column) {
            const double theta = std::numbers::pi / 3 * (1 + double(row) / n);
            const double phi = std::numbers::pi / 2 * double(column) / n;
            source.vertices.push_back({.03 * std::sin(theta) * std::cos(phi),
                                       .03 * std::sin(theta) * std::sin(phi),
                                       .03 * std::cos(theta)});
        }
    for (unsigned row = 0; row < n; ++row)
        for (unsigned column = 0; column < n; ++column) {
            const auto a = row * (n + 1) + column, b = a + n + 1;
            source.triangles.push_back({{a, b, b + 1}, 1, BoundaryCondition::wall});
            source.triangles.push_back({{a, b + 1, a + 1}, 1, BoundaryCondition::wall});
        }
    return source;
}
void curvedBoundaryLens() {
    const auto source = sphericalPatch();
    auto settings = options();
    settings.geometry = std::make_shared<SphereGeometry>();
    const auto atlas = makeSurfaceAtlas(source, settings);
    for (const auto &chart : atlas.charts)
        for (const auto &t : chart.reference) {
            const ChartPoint uv{(t.uv[0].u + t.uv[1].u + t.uv[2].u) / 3,
                                (t.uv[0].v + t.uv[1].v + t.uv[2].v) / 3};
            const auto lifted = chart.sample(uv);
            require(std::abs(std::hypot(lifted.position.x, lifted.position.y, lifted.position.z) -
                             .03) < 1e-14,
                    "A curved boundary reference left its exact sphere.");
            const Vec3 normal{lifted.du.y * lifted.dv.z - lifted.du.z * lifted.dv.y,
                              lifted.du.z * lifted.dv.x - lifted.du.x * lifted.dv.z,
                              lifted.du.x * lifted.dv.y - lifted.du.y * lifted.dv.x};
            require(normal.x * lifted.position.x + normal.y * lifted.position.y +
                            normal.z * lifted.position.z >
                        0,
                    "A curve-only lens collapsed or reversed the exact surface lifting.");
        }
}
// A provider may expose a regular map only in its declared local source
// neighbourhood. Optional fan extensions must not turn that limit into a
// complete-atlas failure or quietly substitute a different projection.
class LocalSphereGeometry final : public SurfaceGeometry {
  public:
    std::span<const BoundaryCorner> boundaryCorners() const noexcept override {
        return {};
    }
    explicit LocalSphereGeometry(const SurfaceMesh &source) : source_(source) {}
    mutable std::atomic<unsigned> domainRejections{};
    SurfaceProjection projectSurface(std::uint32_t face, Vec3 p) const override {
        return exact_.projectSurface(face, p);
    }
    SurfaceProjection projectSupportingSurface(std::uint32_t face, Vec3 p) const override {
        if (!inSourceNeighbourhood(p)) {
            ++domainRejections;
            throw SurfaceProjectionDomainError("Outside the local supporting-sphere domain.");
        }
        return exact_.projectSupportingSurface(face, p);
    }
    Vec3 projectBoundary(std::span<const std::uint32_t> ids, Vec3 p) const override {
        return exact_.projectBoundary(ids, p);
    }
    Vec3 projectJunction(std::span<const std::uint32_t> ids, Vec3 p) const override {
        return exact_.projectJunction(ids, p);
    }

  private:
    const SurfaceMesh source_;
    SphereGeometry exact_;
    bool inSourceNeighbourhood(Vec3 p) const {
        const auto sub = [](Vec3 a, Vec3 b) { return Vec3{a.x - b.x, a.y - b.y, a.z - b.z}; };
        const auto dot = [](Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; };
        for (const auto &t : source_.triangles) {
            const auto a = source_.vertices[t.vertices[0]];
            const auto e = sub(source_.vertices[t.vertices[1]], a);
            const auto f = sub(source_.vertices[t.vertices[2]], a), q = sub(p, a);
            const Vec3 n{e.y * f.z - e.z * f.y, e.z * f.x - e.x * f.z, e.x * f.y - e.y * f.x};
            if (std::abs(dot(q, n)) > 1e-12 * std::sqrt(dot(n, n)))
                continue;
            const double aa = dot(e, e), bb = dot(e, f), cc = dot(f, f), d = aa * cc - bb * bb;
            const double u = (cc * dot(q, e) - bb * dot(q, f)) / d;
            const double v = (aa * dot(q, f) - bb * dot(q, e)) / d;
            if (u >= -1e-10 && v >= -1e-10 && u + v <= 1 + 1e-10)
                return true;
        }
        return false;
    }
};
SurfaceMesh sphericalAnnulus() {
    SurfaceMesh source;
    constexpr unsigned count = 32;
    for (unsigned row = 0; row < 3; ++row)
        for (unsigned i = 0; i < count; ++i) {
            const double theta = std::numbers::pi / 6 + row * std::numbers::pi / 12;
            const double phi = 2 * std::numbers::pi * i / count;
            source.vertices.push_back({.03 * std::sin(theta) * std::cos(phi),
                                       .03 * std::sin(theta) * std::sin(phi),
                                       .03 * std::cos(theta)});
        }
    for (unsigned row = 0; row < 2; ++row)
        for (unsigned i = 0; i < count; ++i) {
            const auto a = row * count + i, b = row * count + (i + 1) % count;
            source.triangles.push_back({{a, a + count, b + count}, 1, BoundaryCondition::wall});
            source.triangles.push_back({{a, b + count, b}, 1, BoundaryCondition::wall});
        }
    return source;
}
void optionalSupportingDomain() {
    const auto source = sphericalAnnulus();
    auto geometry = std::make_shared<LocalSphereGeometry>(source);
    auto settings = options();
    settings.geometry = geometry;
    settings.targetEdgeLength = .012;
    settings.adaptive.geometricTolerance = .001;
    auto atlas = makeSurfaceAtlas(source, settings);
    require(geometry->domainRejections > 0,
            "Optional supporting-domain regression did not exercise a rejected lens.");
    std::set<std::uint32_t> active;
    for (const auto &chart : atlas.charts)
        for (const auto &loop : chart.loops)
            for (const auto &v : loop.vertices)
                active.insert(v.atlasVertex);
    for (unsigned i = 0; i < 32; ++i) {
        const auto p = source.vertices[i];
        require(std::any_of(active.begin(), active.end(),
                            [&](auto id) {
                                const auto q = atlas.vertices[id];
                                return std::hypot(p.x - q.x, p.y - q.y, p.z - q.z) < 1e-14;
                            }),
                "An unavailable optional lens discarded an active original hole knot.");
    }
    auto chart = atlas.charts.front();
    chart.reference = {chart.reference.front()};
    chart.lookup.reset();
    chart.sizing.reset();
    chart.planar = false;
    for (auto &p : chart.reference.front().positions)
        p.x += .1;
    const auto &uv = chart.reference.front().uv;
    bool rejected = false;
    try {
        (void)chart.sample({(uv[0].u + uv[1].u + uv[2].u) / 3, (uv[0].v + uv[1].v + uv[2].v) / 3});
    } catch (const SurfaceProjectionDomainError &) {
        rejected = true;
    }
    require(rejected, "A used chart interior silently replaced an unavailable supporting map.");
}
class AnnulusGeometry final : public SurfaceGeometry {
  public:
    std::span<const BoundaryCorner> boundaryCorners() const noexcept override {
        return {};
    }
    SurfaceProjection projectSurface(std::uint32_t, Vec3 p) const override {
        const double r = std::hypot(p.x, p.y);
        if (!(r > 0))
            throw std::runtime_error("Hole centre has no unique trimmed projection.");
        const double scale = r < 1 ? 1 / r : 1.;
        return {{p.x * scale, p.y * scale, 0}, {0, 0, 1}};
    }
    SurfaceProjection projectSupportingSurface(std::uint32_t, Vec3 p) const override {
        return {{p.x, p.y, 0}, {0, 0, 1}};
    }
    Vec3 projectJunction(std::span<const std::uint32_t>, Vec3) const override {
        throw std::runtime_error("No triple junction in this chart fixture.");
    }
    Vec3 projectBoundary(std::span<const std::uint32_t>, Vec3) const override {
        throw std::runtime_error("No shared boundary query in this chart fixture.");
    }
};
void curvedTrimInteriorRetainsRank() {
    SurfaceChart chart;
    chart.faceId = 1;
    chart.geometry = std::make_shared<AnnulusGeometry>();
    ChartReferenceTriangle reference;
    reference.uv = {{{1, 0}, {1, 1}, {0, 1}}};
    reference.positions = {{{1, 0, 0}, {1, 1, 0}, {0, 1, 0}}};
    reference.sizes = {{.2, .2, .2}};
    chart.reference.push_back(reference);
    const auto a = chart.sample({.55, .55});
    const auto b = chart.sample({.56, .55});
    const auto c = chart.sample({.55, .56});
    require(std::abs(a.position.x - .55) < 1e-15 && std::abs(a.position.y - .55) < 1e-15 &&
                std::abs(b.position.x - .56) < 1e-15 && std::abs(b.position.y - .55) < 1e-15 &&
                std::abs(c.position.x - .55) < 1e-15 && std::abs(c.position.y - .56) < 1e-15 &&
                a.position.z == 0 && b.position.z == 0 && c.position.z == 0,
            "Chart interior was clamped onto its curved physical trim.");
    const double area = (b.position.x - a.position.x) * (c.position.y - a.position.y) -
                        (b.position.y - a.position.y) * (c.position.x - a.position.x);
    require(std::abs(area - .0001) < 1e-16 && a.du == Vec3{1, 0, 0} && a.dv == Vec3{0, 1, 0},
            "Supporting-surface sampling lost its two-dimensional interior map.");
    const auto trimmed = chart.geometry->projectSurface(1, {.55, .55, 0});
    require(std::abs(std::hypot(trimmed.position.x, trimmed.position.y) - 1) < 1e-14,
            "Supporting query changed the independent trimmed-face contract.");
}
class InvalidGeometry final : public SurfaceGeometry {
  public:
    std::span<const BoundaryCorner> boundaryCorners() const noexcept override {
        return {};
    }
    SurfaceProjection projectSurface(std::uint32_t, Vec3 p) const override {
        return {p, {std::numeric_limits<double>::quiet_NaN(), 0, 0}};
    }
    SurfaceProjection projectSupportingSurface(std::uint32_t face, Vec3 p) const override {
        return projectSurface(face, p);
    }
    Vec3 projectJunction(std::span<const std::uint32_t>, Vec3) const override {
        throw std::runtime_error("This analytic source has no physical triple junction.");
    }
    Vec3 projectBoundary(std::span<const std::uint32_t>, Vec3) const override {
        throw std::runtime_error("No valid common curve.");
    }
};
} // namespace
int main() {
    try {
        exactAtlas();
        exactMesh();
        sizingDoesNotProject();
        curvedBoundaryLens();
        curvedTrimInteriorRetainsRank();
        optionalSupportingDomain();
        auto invalid = options();
        invalid.geometry = std::make_shared<InvalidGeometry>();
        rejects([&] { makeSurfaceAtlas(cylinder(), invalid); });
        auto atlas = makeSurfaceAtlas(cylinder(), options());
        auto &chart = atlas.charts.front();
        chart.geometry = invalid.geometry;
        const auto &t = chart.reference.front();
        const ChartPoint sample{(t.uv[0].u + t.uv[1].u + t.uv[2].u) / 3,
                                (t.uv[0].v + t.uv[1].v + t.uv[2].v) / 3};
        rejects([&] { chart.sample(sample); });
        auto displaced = cylinder();
        displaced.vertices.front().x *= .99;
        rejects([&] { makeSurfaceAtlas(displaced, options()); });
        std::cout << "Exact surface geometry tests passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
