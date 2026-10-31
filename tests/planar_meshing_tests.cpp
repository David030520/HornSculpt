#include <hornsculpt/meshing.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <stop_token>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
SurfaceMesh box() {
    constexpr unsigned n = 16, stride = n + 1, offset = stride * stride;
    constexpr double step = .0025;
    SurfaceMesh mesh;
    for (unsigned layer = 0; layer < 2; ++layer)
        for (unsigned y = 0; y <= n; ++y)
            for (unsigned x = 0; x <= n; ++x)
                mesh.vertices.push_back({x * step, y * step, layer == 0 ? 0.0 : -.005});
    for (unsigned y = 0; y < n; ++y)
        for (unsigned x = 0; x < n; ++x) {
            const auto a = y * stride + x, b = a + 1, c = a + stride, d = c + 1;
            mesh.triangles.push_back({{a, b, d}, 5, BoundaryCondition::driven});
            mesh.triangles.push_back({{a, d, c}, 5, BoundaryCondition::driven});
            mesh.triangles.push_back({{a + offset, d + offset, b + offset}, 4});
            mesh.triangles.push_back({{a + offset, c + offset, d + offset}, 4});
        }
    std::vector<std::uint32_t> rim;
    for (unsigned i = 0; i < n; ++i)
        rim.push_back(i);
    for (unsigned i = 0; i < n; ++i)
        rim.push_back(i * stride + n);
    for (unsigned i = n; i > 0; --i)
        rim.push_back(n * stride + i);
    for (unsigned i = n; i > 0; --i)
        rim.push_back(i * stride);
    constexpr std::uint32_t sides[]{1, 2, 3, 6};
    for (unsigned i = 0; i < rim.size(); ++i) {
        const auto a = rim[i], b = rim[(i + 1) % rim.size()], face = sides[i / n];
        mesh.triangles.push_back({{a, a + offset, b}, face});
        mesh.triangles.push_back({{b, a + offset, b + offset}, face});
    }
    mesh.fixedVertices = {0, 8 * stride + 8};
    return mesh;
}
std::map<std::uint32_t, double> areas(const SurfaceMesh &mesh) {
    std::map<std::uint32_t, double> result;
    for (const auto &t : mesh.triangles) {
        const auto a = mesh.vertices[t.vertices[0]], b = mesh.vertices[t.vertices[1]],
                   c = mesh.vertices[t.vertices[2]];
        const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
        result[t.faceId] +=
            .5 * std::hypot(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
        require((t.faceId == 5) == (t.boundary == BoundaryCondition::driven),
                "A physical boundary condition changed.");
    }
    return result;
}
void planarReduction() {
    const auto source = box();
    const auto before = inspectMesh(source);
    require(before.watertight(), "The dense-box test fixture is not closed.");
    MeshingOptions options;
    options.targetEdgeLength = .01;
    const auto result = remeshIsotropic(source, options);
    const auto &mesh = result.mesh;
    require(mesh.triangles.size() < .8 * source.triangles.size(),
            "The planar interior was not meaningfully coarsened.");
    const auto &after = result.quality;
    require(after.watertight() && after.eulerCharacteristic == before.eulerCharacteristic &&
                after.connectedComponents == before.connectedComponents,
            "Coarsening changed the mesh topology or embedding.");
    require(after.minimumAngleDegrees >= options.quality.minimumAngleDegrees &&
                after.maximumAspectRatio <= options.quality.maximumAspectRatio &&
                after.maximumEdgeLength <= 4.0 / 3.0 * options.targetEdgeLength,
            "Coarsening violated an isotropic quality gate.");
    const auto originalAreas = areas(source), finalAreas = areas(mesh);
    require(originalAreas.size() == finalAreas.size(), "A physical patch disappeared.");
    for (const auto &[face, area] : originalAreas)
        require(std::abs(finalAreas.at(face) - area) <= 1e-12 * area,
                "The area of a physical patch changed.");
    // Sampling on straight physical seams may change; their complete areas and
    // the closed topology above must survive, as must explicitly fixed points.
    require(mesh.fixedVertices.size() == source.fixedVertices.size(), "A fixed index disappeared.");
    for (unsigned i = 0; i < source.fixedVertices.size(); ++i) {
        const auto p = mesh.vertices.at(mesh.fixedVertices[i]);
        const auto original = source.vertices[source.fixedVertices[i]];
        require(p.x == original.x && p.y == original.y && p.z == original.z,
                "A remapped fixed point changed.");
    }
    std::cout << "Dense planar box: " << source.triangles.size() << " -> " << mesh.triangles.size()
              << " triangles, minAngle=" << after.minimumAngleDegrees
              << ", maxAspect=" << after.maximumAspectRatio << '\n';
    options.faceTargetEdgeLengths = {{5, .004}};
    const auto local = remeshIsotropic(source, options);
    const auto &locallySized = local.mesh;
    const auto countFace = [](const SurfaceMesh &m, std::uint32_t face) {
        return std::count_if(m.triangles.begin(), m.triangles.end(),
                             [&](const Triangle &t) { return t.faceId == face; });
    };
    require(countFace(locallySized, 5) > countFace(mesh, 5) * 1.5,
            "Planar meshing discarded a manually refined physical surface.");
    for (const auto &t : locallySized.triangles) {
        if (t.faceId != 5)
            continue;
        for (unsigned i = 0; i < 3; ++i) {
            const auto a = locallySized.vertices[t.vertices[i]],
                       b = locallySized.vertices[t.vertices[(i + 1) % 3]];
            require(std::hypot(a.x - b.x, a.y - b.y, a.z - b.z) <= 4.0 / 3.0 * .004,
                    "Planar meshing exceeded a local face target.");
        }
    }
    const auto &localQuality = local.quality;
    require(localQuality.watertight() &&
                localQuality.eulerCharacteristic == before.eulerCharacteristic &&
                localQuality.minimumAngleDegrees >= options.quality.minimumAngleDegrees &&
                localQuality.maximumAspectRatio <= options.quality.maximumAspectRatio,
            "Locally bounded meshing damaged mesh topology or quality.");
    options.faceTargetEdgeLengths.clear();
    options.targetEdgeLength = .001;
    const auto limited = remeshIsotropic(source, options);
    require(limited.mesh.triangles.size() > source.triangles.size() &&
                limited.quality.maximumEdgeLength <= 4.0 / 3.0 * options.targetEdgeLength,
            "The maximum output edge-length bound was ignored.");
    std::stop_source cancelled;
    cancelled.request_stop();
    bool rejected = false;
    try {
        remeshIsotropic(source, options, cancelled.get_token());
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    require(rejected, "Cancellation returned a mesh.");
}
void nonplanarPatchIsNotFlattened() {
    SurfaceMesh mesh;
    mesh.vertices = {{.1, 0, 0}, {0, .1, 0}, {-.1, 0, 0}, {0, -.1, 0}, {0, 0, .1}, {0, 0, -.1}};
    constexpr std::array<std::array<std::uint32_t, 3>, 8> faces{
        {{4, 0, 1}, {4, 1, 2}, {4, 2, 3}, {4, 3, 0}, {5, 1, 0}, {5, 2, 1}, {5, 3, 2}, {5, 0, 3}}};
    for (const auto vertices : faces)
        mesh.triangles.push_back({vertices, 1});
    const auto original = mesh;
    MeshingOptions options;
    options.targetEdgeLength = .03;
    options.adaptive.enabled = false;
    const auto result = remeshIsotropic(mesh, options);
    const double volumeError = result.quality.signedVolume / inspectMesh(original).signedVolume - 1;
    std::cout << "Bent same-ID surface: " << result.mesh.triangles.size()
              << " triangles, relative volume error=" << volumeError << '\n';
    require(result.quality.watertight() && result.quality.minimumAngleDegrees >= 20 &&
                result.quality.maximumAspectRatio <= 4 && std::abs(volumeError) < 1e-10,
            "A nonplanar physical patch was flattened or lost its polyhedral shape.");
}
} // namespace
int main() {
    try {
        planarReduction();
        nonplanarPatchIsNotFlattened();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
