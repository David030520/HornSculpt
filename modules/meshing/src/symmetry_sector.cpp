#include <algorithm>
#include <cmath>
#include <hornsculpt/meshing.hpp>
#include <hornsculpt/symmetry_sector.hpp>
#include <map>
#include <set>
#include <stdexcept>

namespace hornsculpt {
namespace {
Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
bool finite(Vec3 p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
struct Edge {
    unsigned count{};
};
void validateOpenVertexLinks(const SurfaceMesh &mesh) {
    std::vector<std::vector<std::array<std::uint32_t, 2>>> links(mesh.vertices.size());
    for (const auto &triangle : mesh.triangles)
        for (unsigned i = 0; i < 3; ++i)
            links[triangle.vertices[i]].push_back(
                {triangle.vertices[(i + 1) % 3], triangle.vertices[(i + 2) % 3]});
    for (const auto &link : links) {
        if (link.empty())
            continue;
        std::map<std::uint32_t, std::vector<std::uint32_t>> graph;
        for (const auto &edge : link) {
            graph[edge[0]].push_back(edge[1]);
            graph[edge[1]].push_back(edge[0]);
        }
        unsigned ends = 0;
        for (const auto &[vertex, neighbors] : graph) {
            (void)vertex;
            if (neighbors.size() == 1)
                ++ends;
            else if (neighbors.size() != 2)
                throw std::invalid_argument("Symmetry sector has a non-manifold vertex link.");
        }
        if (ends != 0 && ends != 2)
            throw std::invalid_argument("Symmetry sector has a non-manifold boundary vertex.");
        std::set<std::uint32_t> visited;
        std::vector<std::uint32_t> pending{graph.begin()->first};
        while (!pending.empty()) {
            const auto vertex = pending.back();
            pending.pop_back();
            if (!visited.insert(vertex).second)
                continue;
            for (auto neighbor : graph.at(vertex))
                pending.push_back(neighbor);
        }
        if (visited.size() != graph.size())
            throw std::invalid_argument("Symmetry sector has a disconnected vertex link.");
    }
}

} // namespace
SymmetrySectorValidation validateSymmetrySector(const SurfaceMesh &mesh,
                                                const OrthogonalSymmetryPlanes &planes,
                                                const ExteriorBaffle *baffle) {
    // Dimensionless unit/orthogonality budget; does not alter caller geometry.
    constexpr double directionTolerance = 1e-10;
    if (!finite(planes.origin) || !finite(planes.firstNormal) || !finite(planes.secondNormal) ||
        !std::isfinite(planes.planeDistanceTolerance) || planes.planeDistanceTolerance <= 0 ||
        std::abs(dot(planes.firstNormal, planes.firstNormal) - 1) > directionTolerance ||
        std::abs(dot(planes.secondNormal, planes.secondNormal) - 1) > directionTolerance ||
        std::abs(dot(planes.firstNormal, planes.secondNormal)) > directionTolerance)
        throw std::invalid_argument("Symmetry planes require finite orthogonal unit normals and a "
                                    "positive metre tolerance.");
    if (baffle &&
        (!finite(baffle->origin) || !finite(baffle->normal) ||
         !std::isfinite(baffle->planeDistanceTolerance) || baffle->planeDistanceTolerance <= 0 ||
         std::abs(dot(baffle->normal, baffle->normal) - 1) > directionTolerance ||
         std::abs(dot(baffle->normal, planes.firstNormal)) > directionTolerance ||
         std::abs(dot(baffle->normal, planes.secondNormal)) > directionTolerance))
        throw std::invalid_argument(
            "Exterior baffle must have a finite unit normal orthogonal to both symmetry normals.");
    const auto quality = inspectMesh(mesh);
    if (quality.degenerateTriangles || quality.nonManifoldEdges || quality.inconsistentEdges ||
        quality.selfIntersectionPairs)
        throw std::invalid_argument(
            "Symmetry sector is not an oriented, nondegenerate, intersection-free manifold.");
    validateOpenVertexLinks(mesh);
    std::vector<std::array<double, 2>> distance;
    distance.reserve(mesh.vertices.size());
    for (const auto p : mesh.vertices) {
        const auto x = subtract(p, planes.origin);
        const double d1 = dot(x, planes.firstNormal), d2 = dot(x, planes.secondNormal);
        if (!std::isfinite(d1) || !std::isfinite(d2) || d1 < -planes.planeDistanceTolerance ||
            d2 < -planes.planeDistanceTolerance)
            throw std::invalid_argument("Symmetry sector extends outside a retained half-space.");
        if (baffle &&
            dot(subtract(p, baffle->origin), baffle->normal) < -baffle->planeDistanceTolerance)
            throw std::invalid_argument("Exterior boundary lies behind the baffle plane.");
        distance.push_back({d1, d2});
    }
    std::map<std::pair<std::uint32_t, std::uint32_t>, Edge> edges;
    for (const auto &triangle : mesh.triangles) {
        for (unsigned plane = 0; plane < 2; ++plane) {
            bool coplanar = true;
            for (auto v : triangle.vertices)
                coplanar &= std::abs(distance[v][plane]) <= planes.planeDistanceTolerance;
            if (coplanar)
                throw std::invalid_argument("A physical triangle lies on a symmetry plane.");
        }
        for (unsigned i = 0; i < 3; ++i)
            ++edges[std::minmax(triangle.vertices[i], triangle.vertices[(i + 1) % 3])].count;
    }
    SymmetrySectorValidation result;
    for (const auto &[edge, use] : edges) {
        if (use.count != 1)
            continue;
        const bool first = std::abs(distance[edge.first][0]) <= planes.planeDistanceTolerance &&
                           std::abs(distance[edge.second][0]) <= planes.planeDistanceTolerance;
        const bool second = std::abs(distance[edge.first][1]) <= planes.planeDistanceTolerance &&
                            std::abs(distance[edge.second][1]) <= planes.planeDistanceTolerance;
        const bool onBaffle =
            baffle &&
            std::abs(dot(subtract(mesh.vertices[edge.first], baffle->origin), baffle->normal)) <=
                baffle->planeDistanceTolerance &&
            std::abs(dot(subtract(mesh.vertices[edge.second], baffle->origin), baffle->normal)) <=
                baffle->planeDistanceTolerance;
        if (!first && !second && !onBaffle)
            throw std::invalid_argument("An open edge is not on a symmetry or baffle plane.");
        if (first)
            ++result.firstPlaneEdges;
        if (second)
            ++result.secondPlaneEdges;
        if (first && second)
            ++result.intersectionAxisEdges;
    }
    return result;
}
void validateBaffleBoundary(const SurfaceMesh &mesh, const ExteriorBaffle &baffle) {
    constexpr double directionTolerance = 1e-10;
    if (!finite(baffle.origin) || !finite(baffle.normal) ||
        !std::isfinite(baffle.planeDistanceTolerance) || baffle.planeDistanceTolerance <= 0 ||
        std::abs(dot(baffle.normal, baffle.normal) - 1) > directionTolerance)
        throw std::invalid_argument(
            "Baffle requires a finite unit normal and positive metre tolerance.");
    const auto quality = inspectMesh(mesh);
    if (quality.degenerateTriangles || quality.nonManifoldEdges || quality.inconsistentEdges ||
        quality.selfIntersectionPairs)
        throw std::invalid_argument(
            "Baffle boundary is not an oriented intersection-free manifold.");
    validateOpenVertexLinks(mesh);
    std::vector<double> distances;
    for (const auto point : mesh.vertices) {
        const double d = dot(subtract(point, baffle.origin), baffle.normal);
        if (!std::isfinite(d) || d < -baffle.planeDistanceTolerance)
            throw std::invalid_argument("Exterior boundary lies behind the baffle plane.");
        distances.push_back(d);
    }
    std::map<std::pair<std::uint32_t, std::uint32_t>, Edge> edges;
    for (const auto &triangle : mesh.triangles)
        for (unsigned i = 0; i < 3; ++i)
            ++edges[std::minmax(triangle.vertices[i], triangle.vertices[(i + 1) % 3])].count;
    for (const auto &[edge, use] : edges)
        if (use.count == 1 && (std::abs(distances[edge.first]) > baffle.planeDistanceTolerance ||
                               std::abs(distances[edge.second]) > baffle.planeDistanceTolerance))
            throw std::invalid_argument("An exterior open edge is not on the baffle plane.");
}
} // namespace hornsculpt
