#include "mesh_validation.hpp"

#ifdef _MSC_VER
#pragma warning(push)
// Boost's terminating output-iterator path is instantiated by CGAL here.
#pragma warning(disable : 4702)
#endif
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif
#include <algorithm>
#include <array>
#include <iterator>
#include <map>
#include <set>
#include <vector>

namespace hornsculpt::detail {
EmbeddingValidation
validateEmbedding(const SurfaceMesh &mesh,
                  std::vector<std::pair<std::size_t, std::size_t>> *intersections) {
    EmbeddingValidation result;
    std::vector<std::vector<std::array<std::uint32_t, 2>>> links(mesh.vertices.size());
    for (const auto &triangle : mesh.triangles) {
        const auto &ids = triangle.vertices;
        for (int j = 0; j < 3; ++j)
            links[ids[j]].push_back({ids[(j + 1) % 3], ids[(j + 2) % 3]});
    }
    // The link of every vertex of a closed 2-manifold is one cycle.
    for (const auto &link : links) {
        if (link.empty())
            continue;
        std::map<std::uint32_t, std::vector<std::uint32_t>> graph;
        for (const auto edge : link) {
            graph[edge[0]].push_back(edge[1]);
            graph[edge[1]].push_back(edge[0]);
        }
        const bool invalid = std::any_of(graph.begin(), graph.end(),
                                         [](const auto &item) { return item.second.size() != 2; });
        std::set<std::uint32_t> visited;
        std::vector<std::uint32_t> pending{graph.begin()->first};
        while (!pending.empty()) {
            const auto vertex = pending.back();
            pending.pop_back();
            if (!visited.insert(vertex).second)
                continue;
            for (const auto neighbor : graph.at(vertex))
                pending.push_back(neighbor);
        }
        if (invalid || visited.size() != graph.size())
            ++result.nonManifoldVertices;
    }
    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    std::vector<Kernel::Point_3> points;
    points.reserve(mesh.vertices.size());
    for (const auto p : mesh.vertices)
        points.emplace_back(p.x, p.y, p.z);
    std::vector<std::array<std::size_t, 3>> triangles;
    triangles.reserve(mesh.triangles.size());
    for (const auto &t : mesh.triangles)
        triangles.push_back({t.vertices[0], t.vertices[1], t.vertices[2]});
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    // Exact predicates distinguish real overlap from permitted shared edges and
    // vertices, including nearly coplanar tangent triangles at flat curve extrema.
    CGAL::Polygon_mesh_processing::triangle_soup_self_intersections(points, triangles,
                                                                    std::back_inserter(pairs));
    result.selfIntersectionPairs = pairs.size();
    if (intersections)
        *intersections = std::move(pairs);
    return result;
}
} // namespace hornsculpt::detail
