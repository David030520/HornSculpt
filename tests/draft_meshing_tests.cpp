#include <hornsculpt/draft_geometry.hpp>
#include <hornsculpt/meshing.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <vector>

using namespace hornsculpt;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
Vec3 sub(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double norm(Vec3 p) {
    return std::hypot(p.x, p.y, p.z);
}
using Edge = std::pair<std::uint32_t, std::uint32_t>;
struct BoundaryMeasures {
    double area{};
    double perimeter{};
    std::size_t rimEdges{};
    std::size_t drivenTriangles{};
};
BoundaryMeasures measure(const SurfaceMesh &mesh) {
    BoundaryMeasures result;
    std::map<Edge, unsigned> drivenEdges;
    for (const auto &t : mesh.triangles) {
        if (t.boundary != BoundaryCondition::driven)
            continue;
        ++result.drivenTriangles;
        const auto a = mesh.vertices[t.vertices[0]];
        const auto b = mesh.vertices[t.vertices[1]];
        const auto c = mesh.vertices[t.vertices[2]];
        result.area += .5 * norm(cross(sub(b, a), sub(c, a)));
        for (unsigned i = 0; i < 3; ++i)
            ++drivenEdges[std::minmax(t.vertices[i], t.vertices[(i + 1) % 3])];
    }
    std::set<Edge> rim;
    for (const auto &[edge, count] : drivenEdges) {
        if (count != 1)
            continue;
        const auto a = mesh.vertices[edge.first], b = mesh.vertices[edge.second];
        result.perimeter += norm(sub(b, a));
        rim.insert(edge);
    }
    result.rimEdges = rim.size();
    std::size_t joinedEdges = 0;
    for (const auto &t : mesh.triangles) {
        require(t.faceId >= 1 && t.faceId <= 5, "Unexpected computational patch in Draft.");
        if (t.boundary != BoundaryCondition::wall)
            continue;
        for (unsigned i = 0; i < 3; ++i)
            if (rim.contains(std::minmax(t.vertices[i], t.vertices[(i + 1) % 3]))) {
                require(t.faceId == static_cast<std::uint32_t>(DraftFace::innerWall),
                        "Driven perimeter is not joined directly to the physical inner wall.");
                ++joinedEdges;
            }
    }
    require(joinedEdges == rim.size(), "A driven boundary edge has no matching wall.");
    return result;
}
void check(const DraftParameters &parameters, double targetEdgeLength, bool qualityLimits = true) {
    const DraftGeometry geometry(parameters);
    const auto source = geometry.previewMesh();
    const auto before = measure(source);
    MeshingOptions options;
    options.targetEdgeLength = targetEdgeLength;
    options.quality.enabled = qualityLimits;
    const auto result = remeshIsotropic(source, options);
    const auto after = measure(result.mesh);
    require(std::abs(after.area - before.area) <= 1e-10 * before.area,
            "Remeshing changed the driven polygon area.");
    require(std::abs(after.perimeter - before.perimeter) <= 1e-10 * before.perimeter,
            "Remeshing changed the driven polygon perimeter.");
    require(after.drivenTriangles < 6 * after.rimEdges,
            "The flat driven disk retains excessive interior refinement.");
    require(result.mesh.fixedVertices.size() == source.fixedVertices.size(),
            "Remeshing lost fixed physical-boundary vertices.");
    for (const auto sourceIndex : source.fixedVertices) {
        const auto original = source.vertices[sourceIndex];
        require(original.z == 0,
                "Draft pins an artificial wall layer instead of its physical boundary.");
        const bool found =
            std::any_of(result.mesh.fixedVertices.begin(), result.mesh.fixedVertices.end(),
                        [&](std::uint32_t index) {
                            require(index < result.mesh.vertices.size(),
                                    "Invalid remapped fixed-vertex index.");
                            const auto p = result.mesh.vertices[index];
                            return p.x == original.x && p.y == original.y && p.z == original.z;
                        });
        require(found, "A fixed boundary coordinate changed or disappeared.");
    }
    require(result.quality.watertight() && result.quality.connectedComponents == 1,
            "Boundary preservation damaged mesh closure.");
    require(result.quality.selfIntersectionPairs == 0,
            "Remeshing introduced intersecting surface triangles.");
    if (qualityLimits)
        require(result.quality.minimumAngleDegrees >= options.quality.minimumAngleDegrees &&
                    result.quality.maximumAspectRatio <= options.quality.maximumAspectRatio,
                "Boundary preservation bypassed the selected quality limits.");
    std::cout << "Entry=" << parameters.entry.width * 1000.0 << "x"
              << parameters.entry.height * 1000.0 << "mm n=" << parameters.entry.exponent
              << " alpha0=" << parameters.profile.throatAngle * 180.0 / std::numbers::pi
              << "deg h=" << targetEdgeLength * 1000.0 << "mm qualityLimits=" << qualityLimits
              << ": areaRelativeError=" << after.area / before.area - 1.0
              << " perimeterRelativeError=" << after.perimeter / before.perimeter - 1.0
              << " rimEdges=" << after.rimEdges << " triangles=" << result.mesh.triangles.size()
              << " drivenTriangles=" << after.drivenTriangles
              << " minAngle=" << result.quality.minimumAngleDegrees
              << " maxAspect=" << result.quality.maximumAspectRatio << '\n';
}
} // namespace
int main() {
    try {
        DraftParameters parameters;
        check(parameters, .008);
        check(parameters, .008, false);
        check(parameters, .004);
        parameters.profile.throatAngle = .2;
        parameters.profile.k = .7;
        check(parameters, .008);
        parameters.entry = {.032, .008, 8.0};
        parameters.profile.throatAngle = .12;
        parameters.verticalProfile.throatAngle = .21;
        check(parameters, .008);
        parameters = {};
        parameters.crossSectionExponent = 12;
        check(parameters, .016);
        parameters = {};
        parameters.crossSectionExponent = 8;
        for (double exponent : {4., 8., 12.}) {
            parameters.entry = {.1, .01, exponent};
            check(parameters, .016);
        }
        parameters.entry = {.04, .01, 12.};
        check(parameters, .016);
        parameters = {};
        parameters.mouthMorph.enabled = true;
        check(parameters, .008);
        parameters = {};
        parameters.asymmetricEnabled = true;
        parameters.leftAngle = .5;
        parameters.bottomAngle = .3;
        parameters.leftProfile.throatAngle = .1;
        parameters.bottomProfile.terminationScale = .18;
        check(parameters, .008);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
