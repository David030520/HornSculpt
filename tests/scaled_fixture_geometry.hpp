#pragma once

#include "fixture_geometry.hpp"
#include <hornsculpt/meshing.hpp>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <vector>

namespace hornsculpt::testing {

inline SurfaceMesh subdividePolyhedralFixture(const SurfaceMesh &source) {
    SurfaceMesh result;
    result.vertices = source.vertices;
    result.fixedVertices = source.fixedVertices;
    result.triangles.reserve(4 * source.triangles.size());
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> midpoints;
    auto midpoint = [&](std::uint32_t a, std::uint32_t b) {
        const auto key = std::minmax(a, b);
        if (const auto found = midpoints.find(key); found != midpoints.end())
            return found->second;
        const auto x = source.vertices[a], y = source.vertices[b];
        const auto index = static_cast<std::uint32_t>(result.vertices.size());
        result.vertices.push_back({(x.x + y.x) * .5, (x.y + y.y) * .5, (x.z + y.z) * .5});
        midpoints.emplace(key, index);
        return index;
    };
    for (const auto &triangle : source.triangles) {
        const auto a = triangle.vertices[0], b = triangle.vertices[1], c = triangle.vertices[2],
                   ab = midpoint(a, b), bc = midpoint(b, c), ac = midpoint(a, c);
        for (const auto vertices : std::array<std::array<std::uint32_t, 3>, 4>{
                 {{a, ab, ac}, {ab, b, bc}, {ac, bc, c}, {ab, bc, ac}}})
            result.triangles.push_back({vertices, triangle.faceId, triangle.boundary});
    }
    return result;
}

// The projector is identity for a fixed polyhedral surface. A sphere fixture
// supplies radial projection explicitly instead of changing physical CAD geometry.
template <class ProjectMidpoint>
inline SurfaceMesh refineFixture(SurfaceMesh mesh, std::size_t triangleCount,
                                 ProjectMidpoint projectMidpoint,
                                 const MeshQualityLimits &limits = {}) {
    if (triangleCount < mesh.triangles.size() || (triangleCount - mesh.triangles.size()) % 2 != 0)
        throw std::invalid_argument("Refinement needs an even nonnegative panel increment.");
    auto acceptable = [&](Vec3 a, Vec3 b, Vec3 c) {
        const Vec3 ab{b.x - a.x, b.y - a.y, b.z - a.z}, ac{c.x - a.x, c.y - a.y, c.z - a.z};
        const double area2 = std::hypot(ab.y * ac.z - ab.z * ac.y, ab.z * ac.x - ab.x * ac.z,
                                        ab.x * ac.y - ab.y * ac.x);
        const std::array lengths{std::hypot(ab.x, ab.y, ab.z), std::hypot(ac.x, ac.y, ac.z),
                                 std::hypot(c.x - b.x, c.y - b.y, c.z - b.z)};
        if (!(area2 > 0))
            return false;
        if (!limits.enabled)
            return true;
        const double longest = *std::max_element(lengths.begin(), lengths.end());
        if (longest * longest / area2 > limits.maximumAspectRatio)
            return false;
        for (unsigned i = 0; i < 3; ++i) {
            const double opposite = lengths[i], first = lengths[(i + 1) % 3],
                         second = lengths[(i + 2) % 3];
            const double angle = std::acos(std::clamp(
                (first * first + second * second - opposite * opposite) / (2 * first * second),
                -1.0, 1.0));
            if (angle * 180 / std::numbers::pi < limits.minimumAngleDegrees)
                return false;
        }
        return true;
    };
    using Edge = std::pair<std::uint32_t, std::uint32_t>;
    while (mesh.triangles.size() < triangleCount) {
        std::map<Edge, std::vector<std::size_t>> uses;
        for (std::size_t t = 0; t < mesh.triangles.size(); ++t)
            for (unsigned e = 0; e < 3; ++e) {
                const auto &v = mesh.triangles[t].vertices;
                uses[std::minmax(v[e], v[(e + 1) % 3])].push_back(t);
            }
        std::vector<std::pair<double, Edge>> candidates;
        for (const auto &[edge, adjacent] : uses) {
            if (adjacent.size() != 2)
                throw std::logic_error("Fixture edge must have exactly two neighbors.");
            const auto a = mesh.vertices[edge.first], b = mesh.vertices[edge.second];
            candidates.push_back({std::hypot(a.x - b.x, a.y - b.y, a.z - b.z), edge});
        }
        std::sort(candidates.begin(), candidates.end(), std::greater<>());
        std::vector<std::optional<Edge>> selected(mesh.triangles.size());
        std::map<Edge, std::uint32_t> midpoints;
        const auto wanted = (triangleCount - mesh.triangles.size()) / 2;
        for (const auto &[length, edge] : candidates) {
            (void)length;
            const auto &adjacent = uses.at(edge);
            if (selected[adjacent[0]] || selected[adjacent[1]])
                continue;
            const auto a = mesh.vertices[edge.first], b = mesh.vertices[edge.second];
            const auto midpoint =
                projectMidpoint(Vec3{(a.x + b.x) * .5, (a.y + b.y) * .5, (a.z + b.z) * .5});
            bool valid = true;
            for (const auto t : adjacent) {
                const auto &v = mesh.triangles[t].vertices;
                for (const auto vertex : v)
                    if (vertex != edge.first && vertex != edge.second)
                        valid &= acceptable(a, midpoint, mesh.vertices[vertex]) &&
                                 acceptable(midpoint, b, mesh.vertices[vertex]);
            }
            if (!valid)
                continue;
            selected[adjacent[0]] = selected[adjacent[1]] = edge;
            midpoints.emplace(edge, static_cast<std::uint32_t>(mesh.vertices.size()));
            mesh.vertices.push_back(midpoint);
            if (midpoints.size() == wanted)
                break;
        }
        if (midpoints.empty())
            throw std::runtime_error(
                "Cannot refine the fixture while retaining its quality limits.");
        std::vector<Triangle> next;
        next.reserve(mesh.triangles.size() + 2 * midpoints.size());
        for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
            const auto &triangle = mesh.triangles[t];
            if (!selected[t]) {
                next.push_back(triangle);
                continue;
            }
            const auto &v = triangle.vertices;
            for (unsigned e = 0; e < 3; ++e)
                if (Edge(std::minmax(v[e], v[(e + 1) % 3])) == *selected[t]) {
                    const auto m = midpoints.at(*selected[t]);
                    next.push_back({{v[e], m, v[(e + 2) % 3]}, triangle.faceId, triangle.boundary});
                    next.push_back(
                        {{m, v[(e + 1) % 3], v[(e + 2) % 3]}, triangle.faceId, triangle.boundary});
                    break;
                }
        }
        mesh.triangles = std::move(next);
    }
    return mesh;
}

// Geodesic octahedron lattice, followed by conforming longest-edge bisection.
// The exact even panel count avoids misleading "100k" evidence at smaller meshes.
inline SurfaceMesh scaledSphere(std::size_t triangleCount) {
    if (triangleCount < 8 || triangleCount % 2 != 0)
        throw std::invalid_argument("Sphere triangle count must be even and at least eight.");
    const auto divisions = static_cast<int>(std::sqrt(static_cast<double>(triangleCount) / 8));
    const auto seed = sphere(0);
    SurfaceMesh mesh;
    std::map<std::array<int, 3>, std::uint32_t> nodes;
    auto normalized = [](Vec3 p) {
        const auto scale = .1 / std::hypot(p.x, p.y, p.z);
        return Vec3{p.x * scale, p.y * scale, p.z * scale};
    };
    for (const auto &face : seed.triangles) {
        const auto a = seed.vertices[face.vertices[0]], b = seed.vertices[face.vertices[1]],
                   c = seed.vertices[face.vertices[2]];
        auto node = [&](int i, int j) {
            const auto k = divisions - i - j;
            const std::array<int, 3> key{
                static_cast<int>(std::lround(10 * (k * a.x + i * b.x + j * c.x))),
                static_cast<int>(std::lround(10 * (k * a.y + i * b.y + j * c.y))),
                static_cast<int>(std::lround(10 * (k * a.z + i * b.z + j * c.z)))};
            const auto [found, inserted] =
                nodes.emplace(key, static_cast<std::uint32_t>(mesh.vertices.size()));
            if (inserted)
                mesh.vertices.push_back(
                    normalized({double(key[0]), double(key[1]), double(key[2])}));
            return found->second;
        };
        for (int i = 0; i < divisions; ++i)
            for (int j = 0; j < divisions - i; ++j) {
                mesh.triangles.push_back(
                    {{node(i, j), node(i + 1, j), node(i, j + 1)}, 1, BoundaryCondition::driven});
                if (i + j + 1 < divisions)
                    mesh.triangles.push_back({{node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)},
                                              1,
                                              BoundaryCondition::driven});
            }
    }
    return refineFixture(std::move(mesh), triangleCount, normalized);
}
inline SurfaceMesh refinePolyhedralFixture(SurfaceMesh source, std::size_t triangleCount) {
    // Uniform subdivision refines the initially least favorable elements
    // without changing their angles. The remaining increment is conforming.
    while (source.triangles.size() <= triangleCount / 4)
        source = subdividePolyhedralFixture(source);
    return refineFixture(std::move(source), triangleCount, [](Vec3 p) { return p; });
}
} // namespace hornsculpt::testing
