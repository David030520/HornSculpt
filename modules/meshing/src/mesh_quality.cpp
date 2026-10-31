#include "mesh_validation.hpp"
#include <hornsculpt/meshing.hpp>
#include <hornsculpt/triangle_metrics.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
#include <stdexcept>

namespace hornsculpt {
namespace {
Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
double length(Vec3 a) {
    return std::sqrt(dot(a, a));
}
struct EdgeUse {
    std::size_t count{};
    int orientation{};
};
} // namespace

bool MeshQuality::watertight() const noexcept {
    return boundaryEdges == 0 && nonManifoldEdges == 0 && inconsistentEdges == 0 &&
           degenerateTriangles == 0 && connectedComponents > 0 && selfIntersectionPairs == 0 &&
           nonManifoldVertices == 0;
}

MeshQuality inspectMesh(const SurfaceMesh &mesh) {
    if (mesh.vertices.empty() || mesh.triangles.empty())
        throw std::invalid_argument("The surface mesh is empty.");
    Vec3 lower = mesh.vertices.front(), upper = lower;
    for (const auto p : mesh.vertices) {
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
            throw std::invalid_argument("The mesh contains non-finite coordinates.");
        lower = {std::min(lower.x, p.x), std::min(lower.y, p.y), std::min(lower.z, p.z)};
        upper = {std::max(upper.x, p.x), std::max(upper.y, p.y), std::max(upper.z, p.z)};
    }
    // Relative to the model's diagonal, rather than an implicit metre-sized object.
    const double areaTolerance = dot(subtract(upper, lower), subtract(upper, lower)) * 1e-14;
    std::map<std::pair<std::uint32_t, std::uint32_t>, EdgeUse> edges;
    std::vector<std::size_t> parent(mesh.vertices.size());
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&parent](std::size_t i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    };
    std::set<std::array<std::uint32_t, 3>> uniqueTriangles;
    std::set<std::uint32_t> used;
    MeshQuality q;
    q.minimumAngleDegrees = 180.0;
    q.minimumEdgeLength = std::numeric_limits<double>::infinity();
    for (const auto &t : mesh.triangles) {
        for (auto v : t.vertices) {
            if (v >= mesh.vertices.size())
                throw std::invalid_argument("Triangle index outside vertex array.");
            used.insert(v);
        }
        auto canonical = t.vertices;
        std::sort(canonical.begin(), canonical.end());
        if (!uniqueTriangles.insert(canonical).second)
            throw std::invalid_argument("The mesh contains duplicate triangles.");
        const Vec3 a = mesh.vertices[t.vertices[0]], b = mesh.vertices[t.vertices[1]],
                   c = mesh.vertices[t.vertices[2]];
        const auto metrics = triangleMetrics(a, b, c);
        if (metrics.twiceArea <= areaTolerance) {
            ++q.degenerateTriangles;
            q.maximumAspectRatio = std::numeric_limits<double>::infinity();
            q.minimumAngleDegrees = 0;
        } else {
            q.maximumAspectRatio = std::max(q.maximumAspectRatio, metrics.maximumAspectRatio);
            q.minimumAngleDegrees = std::min(q.minimumAngleDegrees, metrics.minimumAngleDegrees);
        }
        // Translation by a common origin reduces cancellation for remote CAD placements.
        q.signedVolume +=
            dot(subtract(a, lower), cross(subtract(b, lower), subtract(c, lower))) / 6;
        for (int j = 0; j < 3; ++j) {
            const auto i = t.vertices[j], k = t.vertices[(j + 1) % 3];
            const auto edge = std::minmax(i, k);
            auto &use = edges[{edge.first, edge.second}];
            ++use.count;
            use.orientation += i < k ? 1 : -1;
            const auto pi = find(i), pk = find(k);
            if (pi != pk)
                parent[pi] = pk;
        }
    }
    std::vector<double> edgeLengths;
    edgeLengths.reserve(edges.size());
    for (const auto &[edge, use] : edges) {
        if (use.count == 1)
            ++q.boundaryEdges;
        if (use.count > 2)
            ++q.nonManifoldEdges;
        if (use.count == 2 && use.orientation != 0)
            ++q.inconsistentEdges;
        const double l = length(subtract(mesh.vertices[edge.first], mesh.vertices[edge.second]));
        edgeLengths.push_back(l);
        q.minimumEdgeLength = std::min(q.minimumEdgeLength, l);
        q.maximumEdgeLength = std::max(q.maximumEdgeLength, l);
    }
    std::sort(edgeLengths.begin(), edgeLengths.end());
    q.medianEdgeLength = edgeLengths[edgeLengths.size() / 2];
    std::set<std::size_t> components;
    for (auto v : used)
        components.insert(find(v));
    q.connectedComponents = components.size();
    q.eulerCharacteristic = static_cast<std::int64_t>(used.size()) -
                            static_cast<std::int64_t>(edges.size()) +
                            static_cast<std::int64_t>(mesh.triangles.size());
    const auto embedding = detail::validateEmbedding(mesh);
    q.selfIntersectionPairs = embedding.selfIntersectionPairs;
    q.nonManifoldVertices = embedding.nonManifoldVertices;
    return q;
}
} // namespace hornsculpt
