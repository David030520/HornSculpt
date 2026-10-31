#include "../modules/meshing/src/surface_sizing_field.hpp"
#include <hornsculpt/draft_geometry.hpp>
#include <hornsculpt/meshing.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class F> void rejects(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const std::exception &) {
        rejected = true;
    }
    require(rejected, "Invalid mesh was accepted.");
}
void surfaceFieldContinuity() {
    SurfaceMesh mesh;
    mesh.vertices = {{0, 0, 0},     {1, 0, 0},     {0, 1, 0},     {0, 0, 1},
                     {0, 0, 1.001}, {1, 0, 1.001}, {0, 1, 1.001}, {0, 0, 2.001}};
    constexpr std::array<std::array<std::uint32_t, 3>, 4> tetra{
        {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}}};
    for (const auto &t : tetra) {
        mesh.triangles.push_back({t, 1, BoundaryCondition::wall});
        mesh.triangles.push_back({{t[0] + 4, t[1] + 4, t[2] + 4}, 2, BoundaryCondition::wall});
    }
    std::vector<double> sizes{.001, .001, .001, .001, .1, .1, .1, .1};
    detail::SurfaceSizingField field(mesh, sizes, .3, .1, {}, {});
    const std::array<std::uint32_t, 1> first{1}, second{2};
    require(std::abs(field.at({.01, .01, 1.001}, second) - .1) < 1e-12,
            "Sizing bled across the air gap between disconnected closed bodies.");
    require(std::abs(field.at({0, 0, 1}, first) - .001) < 1e-12,
            "The small feature's own surface target was lost.");
    SurfaceMesh joined;
    joined.vertices.assign(mesh.vertices.begin(), mesh.vertices.begin() + 4);
    for (std::size_t i = 0; i < tetra.size(); ++i)
        joined.triangles.push_back({tetra[i], i == 0 ? 1u : 2u, BoundaryCondition::wall});
    detail::SurfaceSizingField shared(joined, {.001, .05, .05, .05}, .1, .1, {}, {});
    const Vec3 seam{.25, 0, 0};
    require(std::abs(shared.at(seam, first) - shared.at(seam, second)) < 1e-12,
            "The surface size field is discontinuous across a physical patch seam.");
    SurfaceMesh fan;
    fan.vertices = {{0, 0, 0}, {.08, 0, 0}, {0, .08, 0}};
    fan.triangles = {{{0, 1, 2}, 1, BoundaryCondition::wall}};
    detail::SurfaceSizingField broad(fan, {.002, .002, .002}, .3, .008, {}, {});
    require(broad.at({.025, .025, 0}, first) > .007,
            "Small boundary seeds prevented coarse sizing inside a broad CAD fan.");
    const std::array<detail::SurfaceEdgeSeed, 1> feature{{{{0, 1}, .002}}};
    detail::SurfaceSizingField seeded(fan, {.002, .002, .002}, .3, .008, feature, {});
    require(seeded.at({.04, 0, 0}, first) <= .002 + 1e-12,
            "Field subdivision forgot a physical feature along a long original edge.");
    const std::map<std::uint32_t, double> patchTarget{{1, .002}};
    detail::SurfaceSizingField selected(fan, {.008, .008, .008}, .3, .008, {}, {}, patchTarget);
    require(selected.at({.025, .025, 0}, first) <= .002 + 1e-12,
            "The local face target was lost at new interior field samples.");
    SurfaceMesh skinny;
    skinny.vertices = {{0, 0, 0}, {.01, 0, 0}, {.01, 1e-11, 0}};
    skinny.triangles = {{{0, 1, 2}, 1, BoundaryCondition::wall}};
    detail::SurfaceSizingField thin(skinny, {.01, .01, .01}, .1, .01, {}, {});
    require(std::abs(thin.at({.005, 2.5e-12, 0}, first) - .01) < 1e-12,
            "Barycentric sizing lost a valid thin source triangle to cancellation.");
}
void optionalQualityLimits() {
    // A valid closed boundary with deliberately skinny, fixed physical facets.
    SurfaceMesh source;
    source.vertices = {{0, 0, 0}, {.01, 0, 0}, {.01, .0001, 0}, {0, 0, .01}};
    source.fixedVertices = {0, 1, 2, 3};
    source.triangles = {{{0, 2, 1}, 1, BoundaryCondition::driven},
                        {{0, 1, 3}, 2, BoundaryCondition::wall},
                        {{0, 3, 2}, 3, BoundaryCondition::wall},
                        {{1, 2, 3}, 4, BoundaryCondition::wall}};
    require(inspectMesh(source).watertight(), "The quality-policy fixture must be closed.");
    MeshingOptions options{.02};
    rejects([&] { remeshIsotropic(source, options); });
    options.quality.enabled = false;
    const auto unrestricted = remeshIsotropic(source, options);
    require(unrestricted.quality.watertight() && unrestricted.quality.minimumAngleDegrees < 1 &&
                unrestricted.quality.maximumAspectRatio > 50 &&
                unrestricted.mesh.triangles.size() == 4,
            "Disabling shape limits still rejects or refines fixed skinny elements.");
    options.quality = {true, .01, 1000};
    require(remeshIsotropic(source, options).quality.watertight(),
            "Custom shape limits are not used by the final boundary check.");
    options.quality.minimumAngleDegrees = 1;
    rejects([&] { remeshIsotropic(source, options); });
    options.quality = {true, .01, 50};
    rejects([&] { remeshIsotropic(source, options); });
    options.quality.enabled = false;
    auto open = source;
    open.triangles.pop_back();
    rejects([&] { remeshIsotropic(open, options); });
    auto degenerate = source;
    degenerate.vertices[2] = degenerate.vertices[1];
    rejects([&] { remeshIsotropic(degenerate, options); });
}
SurfaceMesh sphere() {
    SurfaceMesh mesh;
    mesh.vertices = {{.1, 0, 0}, {0, .1, 0}, {-.1, 0, 0}, {0, -.1, 0}, {0, 0, .1}, {0, 0, -.1}};
    const std::array<std::array<std::uint32_t, 3>, 8> faces = {
        {{4, 0, 1}, {4, 1, 2}, {4, 2, 3}, {4, 3, 0}, {5, 1, 0}, {5, 2, 1}, {5, 3, 2}, {5, 0, 3}}};
    for (std::uint32_t i = 0; i < faces.size(); ++i)
        mesh.triangles.push_back(
            {faces[i], i + 1, i == 0 ? BoundaryCondition::driven : BoundaryCondition::wall});
    for (int refinement = 0; refinement < 3; ++refinement) {
        SurfaceMesh next;
        next.vertices = mesh.vertices;
        std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> midpoints;
        auto midpoint = [&](std::uint32_t a, std::uint32_t b) {
            const auto edge = std::minmax(a, b);
            const auto key = std::pair{edge.first, edge.second};
            if (auto it = midpoints.find(key); it != midpoints.end())
                return it->second;
            Vec3 p{mesh.vertices[a].x + mesh.vertices[b].x, mesh.vertices[a].y + mesh.vertices[b].y,
                   mesh.vertices[a].z + mesh.vertices[b].z};
            const double scale = .1 / std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
            p = {p.x * scale, p.y * scale, p.z * scale};
            const auto index = static_cast<std::uint32_t>(next.vertices.size());
            next.vertices.push_back(p);
            midpoints.emplace(key, index);
            return index;
        };
        for (const auto &t : mesh.triangles) {
            const auto a = t.vertices[0], b = t.vertices[1], c = t.vertices[2], ab = midpoint(a, b),
                       bc = midpoint(b, c), ca = midpoint(c, a);
            for (const auto triangle : std::array<std::array<std::uint32_t, 3>, 4>{
                     {{a, ab, ca}, {ab, b, bc}, {ca, bc, c}, {ab, bc, ca}}})
                next.triangles.push_back({triangle, t.faceId, t.boundary});
        }
        mesh = std::move(next);
    }
    return mesh;
}
SurfaceMesh flatPlate() {
    SurfaceMesh mesh;
    mesh.vertices = {{-.05, -.05, 0},    {.05, -.05, 0},    {.05, .05, 0},    {-.05, .05, 0},
                     {-.05, -.05, .005}, {.05, -.05, .005}, {.05, .05, .005}, {-.05, .05, .005}};
    const std::array<std::array<std::uint32_t, 4>, 6> faces{
        {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}}};
    for (std::uint32_t i = 0; i < faces.size(); ++i) {
        const auto q = faces[i];
        const auto boundary = i == 0 ? BoundaryCondition::driven : BoundaryCondition::wall;
        mesh.triangles.push_back({{q[0], q[1], q[2]}, i + 1, boundary});
        mesh.triangles.push_back({{q[0], q[2], q[3]}, i + 1, boundary});
    }
    return mesh;
}
void sourceDensityIndependence(const MeshingResult &coarse) {
    auto source = flatPlate();
    // Subdivision changes only the provider's tessellation, not the physical
    // piecewise-planar surface, its face IDs, or the requested scalar size.
    for (unsigned pass = 0; pass < 4; ++pass) {
        SurfaceMesh refined;
        refined.vertices = source.vertices;
        std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> midpoints;
        const auto midpoint = [&](std::uint32_t a, std::uint32_t b) {
            const auto edge = std::minmax(a, b);
            const auto key = std::pair{edge.first, edge.second};
            if (const auto found = midpoints.find(key); found != midpoints.end())
                return found->second;
            const auto p = source.vertices[a], q = source.vertices[b];
            const auto index = static_cast<std::uint32_t>(refined.vertices.size());
            refined.vertices.push_back(
                {std::midpoint(p.x, q.x), std::midpoint(p.y, q.y), std::midpoint(p.z, q.z)});
            midpoints.emplace(key, index);
            return index;
        };
        for (const auto &triangle : source.triangles) {
            const auto a = triangle.vertices[0], b = triangle.vertices[1], c = triangle.vertices[2],
                       ab = midpoint(a, b), bc = midpoint(b, c), ca = midpoint(c, a);
            for (const auto vertices : std::array<std::array<std::uint32_t, 3>, 4>{
                     {{a, ab, ca}, {ab, b, bc}, {ca, bc, c}, {ab, bc, ca}}})
                refined.triangles.push_back({vertices, triangle.faceId, triangle.boundary});
        }
        source = std::move(refined);
    }
    const auto result = remeshIsotropic(source);
    require(result.quality.watertight() && result.quality.minimumAngleDegrees >= 20 &&
                result.quality.maximumAspectRatio <= 4,
            "Source-density independence bypassed the mesh quality gates.");
    require(result.mesh.triangles.size() < coarse.mesh.triangles.size() * 1.25 &&
                result.mesh.triangles.size() > coarse.mesh.triangles.size() * .75,
            "Provider tessellation density unnecessarily changed the BEM mesh density.");
    require(std::abs(result.quality.signedVolume / coarse.quality.signedVolume - 1) < 1e-10,
            "Resampling an equivalent planar source changed its volume.");
    std::cout << "Equivalent source tessellation: " << source.triangles.size()
              << " source triangles, " << coarse.mesh.triangles.size() << " -> "
              << result.mesh.triangles.size() << " BEM triangles\n";
}
SurfaceMesh longBore() {
    SurfaceMesh mesh;
    constexpr std::uint32_t count = 64;
    for (int ring = 0; ring < 4; ++ring)
        for (std::uint32_t i = 0; i < count; ++i) {
            const double angle = 2 * std::numbers::pi * i / count;
            const double c = std::cos(angle), s = std::sin(angle);
            const double radius = ring < 2 ? .00225 : .025 / std::max(std::abs(c), std::abs(s));
            mesh.vertices.push_back({radius * c, radius * s, ring % 2 ? .06 : 0});
        }
    auto quad = [&](std::array<std::uint32_t, 4> v, std::uint32_t face) {
        const auto condition = face == 6 ? BoundaryCondition::driven : BoundaryCondition::wall;
        mesh.triangles.push_back({{v[0], v[1], v[2]}, face, condition});
        mesh.triangles.push_back({{v[0], v[2], v[3]}, face, condition});
    };
    for (std::uint32_t i = 0; i < count; ++i) {
        const auto j = (i + 1) % count;
        const auto side = 2 + ((i + count / 8) / (count / 4)) % 4;
        quad({i, i + count, j + count, j}, 1);
        quad({i + 2 * count, j + 2 * count, j + 3 * count, i + 3 * count}, side);
        quad({i, j, j + 2 * count, i + 2 * count}, 6);
        quad({i + count, i + 3 * count, j + 3 * count, j + count}, 7);
    }
    return mesh;
}
void check(const MeshingResult &result, std::size_t expectedFaces) {
    const auto &q = result.quality;
    require(q.watertight() && q.connectedComponents == 1,
            "Generated boundary is not one closed manifold.");
    require(q.minimumAngleDegrees >= 20 && q.maximumAspectRatio <= 4,
            "Elongated elements escaped the quality gate.");
    std::set<std::uint32_t> faces;
    bool driven = false;
    for (const auto &t : result.mesh.triangles) {
        faces.insert(t.faceId);
        driven = driven || t.boundary == BoundaryCondition::driven;
    }
    require(faces.size() == expectedFaces && driven,
            "Physical face IDs or driven condition were lost.");
    std::cout << result.mesh.triangles.size() << " triangles; minimum angle "
              << q.minimumAngleDegrees << " deg; maximum aspect " << q.maximumAspectRatio
              << "; edges " << q.minimumEdgeLength << " / " << q.medianEdgeLength << " / "
              << q.maximumEdgeLength << " m\n";
}
void localFaceRefinement(const MeshingResult &coarse) {
    const auto source = flatPlate();
    MeshingOptions options;
    options.adaptive.enabled = false;
    options.faceTargetEdgeLengths = {{2, .002}};
    const auto refined = remeshIsotropic(source, options);
    check(refined, 6);
    const auto statistics = [](const SurfaceMesh &mesh, std::uint32_t face) {
        std::array<double, 3> result{}; // face count, sum of edge lengths, maximum edge
        for (const auto &t : mesh.triangles) {
            if (t.faceId != face)
                continue;
            ++result[0];
            for (std::size_t i = 0; i < 3; ++i) {
                const auto a = mesh.vertices[t.vertices[i]],
                           b = mesh.vertices[t.vertices[(i + 1) % 3]];
                const double edge = std::hypot(a.x - b.x, a.y - b.y, a.z - b.z);
                result[1] += edge;
                result[2] = std::max(result[2], edge);
            }
        }
        result[1] /= 3 * result[0];
        return result;
    };
    const auto before = statistics(coarse.mesh, 2), after = statistics(refined.mesh, 2);
    const auto opposite = statistics(refined.mesh, 1);
    require(after[0] > before[0] * 4 && after[1] < before[1] * .5 && after[2] < .0032,
            "Manual local sizing did not refine the requested face with automatic adaptation "
            "disabled.");
    require(opposite[1] > after[1] * 2,
            "A local face target unnecessarily became a uniform global refinement.");
    constexpr std::array<Vec3, 6> directions{
        {{0, 0, -1}, {0, 0, 1}, {0, -1, 0}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0}}};
    for (const auto &t : refined.mesh.triangles) {
        require(t.faceId >= 1 && t.faceId <= 6, "A local target changed physical surface IDs.");
        require((t.faceId == 1) == (t.boundary == BoundaryCondition::driven),
                "A local target changed a physical boundary condition.");
        const auto a = refined.mesh.vertices[t.vertices[0]],
                   b = refined.mesh.vertices[t.vertices[1]],
                   c = refined.mesh.vertices[t.vertices[2]];
        const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
        const Vec3 n{u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
        const auto expected = directions[t.faceId - 1];
        require(n.x * expected.x + n.y * expected.y + n.z * expected.z > 0,
                "A locally refined surface normal was reversed.");
    }
    std::cout << "Local face target: " << before[0] << " -> " << after[0]
              << " triangles, mean edge " << before[1] << " -> " << after[1]
              << " m, opposite mean edge " << opposite[1] << " m\n";
    for (const auto &[id, size] : std::map<std::uint32_t, double>{{0, .002}, {99, .002}}) {
        options.faceTargetEdgeLengths = {{id, size}};
        rejects([&] { remeshIsotropic(source, options); });
    }
    for (const double size : {0., -.001, .00001, .009, std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::infinity()}) {
        options.faceTargetEdgeLengths = {{2, size}};
        rejects([&] { remeshIsotropic(source, options); });
    }
}
} // namespace
int main() {
    std::cout << std::unitbuf;
    try {
        surfaceFieldContinuity();
        const auto plate = remeshIsotropic(flatPlate());
        check(plate, 6);
        sourceDensityIndependence(plate);
        localFaceRefinement(plate);
        double longestSideEdge = 0;
        for (const auto &triangle : plate.mesh.triangles) {
            if (triangle.faceId < 3)
                continue;
            for (std::size_t i = 0; i < 3; ++i) {
                const auto a = plate.mesh.vertices[triangle.vertices[i]];
                const auto b = plate.mesh.vertices[triangle.vertices[(i + 1) % 3]];
                if (std::abs(a.z - b.z) < 1e-12)
                    longestSideEdge = std::max(longestSideEdge, std::hypot(a.x - b.x, a.y - b.y));
            }
        }
        require(longestSideEdge > .006,
                "Thin flat CAD side forced uniformly small straight-edge elements.");
        auto fixedPlate = flatPlate();
        const auto fixedPoint = Vec3{-.017, .013, 0};
        fixedPlate.vertices.push_back(fixedPoint);
        const auto original = fixedPlate.triangles.front();
        fixedPlate.triangles.front().vertices = {original.vertices[0], original.vertices[1], 8};
        fixedPlate.triangles.push_back(
            {{original.vertices[1], original.vertices[2], 8}, original.faceId, original.boundary});
        fixedPlate.triangles.push_back(
            {{original.vertices[2], original.vertices[0], 8}, original.faceId, original.boundary});
        fixedPlate.fixedVertices = {8};
        const auto fixedResult = remeshIsotropic(fixedPlate);
        check(fixedResult, 6);
        require(fixedResult.mesh.fixedVertices.size() == 1,
                "Fixed source vertex metadata was lost.");
        const auto retained = fixedResult.mesh.vertices[fixedResult.mesh.fixedVertices.front()];
        require(retained.x == fixedPoint.x && retained.y == fixedPoint.y &&
                    retained.z == fixedPoint.z,
                "A fixed source position changed during remeshing.");
        fixedPlate.fixedVertices = {8, 8};
        rejects([&] { remeshIsotropic(fixedPlate); });
        fixedPlate.fixedVertices = {static_cast<std::uint32_t>(fixedPlate.vertices.size())};
        rejects([&] { remeshIsotropic(fixedPlate); });
        optionalQualityLimits();
        const auto ball = sphere();
        require(inspectMesh(ball).watertight(), "Sphere fixture is open.");
        auto open = ball;
        open.triangles.pop_back();
        rejects([&] { remeshIsotropic(open); });
        auto reversed = ball;
        std::swap(reversed.triangles[0].vertices[1], reversed.triangles[0].vertices[2]);
        rejects([&] { remeshIsotropic(reversed); });
        auto duplicate = ball;
        duplicate.triangles.push_back(duplicate.triangles[0]);
        rejects([&] { inspectMesh(duplicate); });
        auto conflicting = ball;
        conflicting.triangles[0].boundary = BoundaryCondition::wall;
        rejects([&] { remeshIsotropic(conflicting); });
        std::stop_source cancelled;
        cancelled.request_stop();
        rejects([&] { remeshIsotropic(ball, {}, cancelled.get_token()); });
        MeshingOptions coarseOptions{.025};
        coarseOptions.adaptive.geometricTolerance = .005;
        const auto coarse = remeshIsotropic(ball, coarseOptions);
        check(coarse, 8);
        auto fineOptions = coarseOptions;
        fineOptions.targetEdgeLength = .017;
        const auto fine = remeshIsotropic(ball, fineOptions);
        check(fine, 8);
        require(fine.mesh.triangles.size() > coarse.mesh.triangles.size() * 1.5,
                "Refining the requested size did not refine the mesh.");
        const auto bore = longBore();
        const auto boreInputQuality = inspectMesh(bore);
        require(boreInputQuality.watertight() && boreInputQuality.eulerCharacteristic == 0,
                "Long-bore fixture must have one preserved through-hole.");
        MeshingOptions boreOptions{.012};
        boreOptions.adaptive.geometricTolerance = .00005;
        const auto adaptedBore = remeshIsotropic(bore, boreOptions);
        check(adaptedBore, 7);
        require(adaptedBore.quality.eulerCharacteristic == 0,
                "Remeshing erased the bore topology.");
        double maximumBoreEdge = 0, maximumFlatEdge = 0;
        for (const auto &t : adaptedBore.mesh.triangles)
            for (int j = 0; j < 3; ++j) {
                const auto a = adaptedBore.mesh.vertices[t.vertices[j]],
                           b = adaptedBore.mesh.vertices[t.vertices[(j + 1) % 3]];
                const double length =
                    std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) +
                              (a.z - b.z) * (a.z - b.z));
                if (t.faceId == 1)
                    maximumBoreEdge = std::max(maximumBoreEdge, length);
                else
                    maximumFlatEdge = std::max(maximumFlatEdge, length);
            }
        require(maximumBoreEdge < .002 && maximumFlatEdge > maximumBoreEdge * 2,
                "Small bore curvature was not refined locally with scalar size grading.");
        MeshingOptions invalid{.008};
        invalid.adaptive.growthRate = 0;
        rejects([&] { remeshIsotropic(bore, invalid); });
        const auto horn = DraftGeometry{}.previewMesh();
        const auto hornResult = remeshIsotropic(horn, {.008});
        check(hornResult, 5);
        const auto repeatedHorn = remeshIsotropic(horn, {.008});
        require(
            hornResult.mesh.fixedVertices == repeatedHorn.mesh.fixedVertices &&
                hornResult.mesh.vertices.size() == repeatedHorn.mesh.vertices.size() &&
                hornResult.mesh.triangles.size() == repeatedHorn.mesh.triangles.size() &&
                std::equal(hornResult.mesh.vertices.begin(), hornResult.mesh.vertices.end(),
                           repeatedHorn.mesh.vertices.begin(),
                           [](Vec3 a, Vec3 b) { return a.x == b.x && a.y == b.y && a.z == b.z; }) &&
                std::equal(hornResult.mesh.triangles.begin(), hornResult.mesh.triangles.end(),
                           repeatedHorn.mesh.triangles.begin(),
                           [](const Triangle &a, const Triangle &b) {
                               return a.vertices == b.vertices && a.faceId == b.faceId &&
                                      a.boundary == b.boundary;
                           }),
            "Identical horn inputs produced a different mesh or physical boundary tags.");
        for (const auto &t : hornResult.mesh.triangles)
            if (t.boundary == BoundaryCondition::driven) {
                require(t.faceId == 5, "Driven face identity changed.");
                for (auto index : t.vertices)
                    require(std::abs(hornResult.mesh.vertices[index].z) < 1e-10,
                            "Driven face left the throat plane.");
                const auto a = hornResult.mesh.vertices[t.vertices[0]],
                           b = hornResult.mesh.vertices[t.vertices[1]],
                           c = hornResult.mesh.vertices[t.vertices[2]];
                require((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) > 0,
                        "Driven face normal reversed.");
            }
        const double volumeRelativeError =
            std::abs(hornResult.quality.signedVolume / inspectMesh(horn).signedVolume - 1);
        require(volumeRelativeError < .025,
                "Remeshing significantly changed the physical shell volume.");
        DraftParameters boxy;
        boxy.crossSectionExponent = 4;
        check(remeshIsotropic(DraftGeometry{boxy}.previewMesh(), {.008}), 5);
        DraftParameters formula;
        formula.guideEnabled = true;
        formula.guide.type = CurveType::superformula;
        formula.guide.superformula.n1 = 4;
        formula.guide.superformula.n2 = 4;
        formula.guide.superformula.n3 = 4;
        check(remeshIsotropic(DraftGeometry{formula}.previewMesh(), {.008}), 5);
        formula.guide.rotation = .3;
        check(remeshIsotropic(DraftGeometry{formula}.previewMesh(), {.008}), 5);
        std::cout << "Meshing tests passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
