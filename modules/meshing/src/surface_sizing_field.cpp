#include "surface_sizing_field.hpp"

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4701 4702)
#endif
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/iterator.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif
#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <utility>

namespace hornsculpt::detail {
namespace {
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using Tree = CGAL::AABB_tree<CGAL::AABB_traits_3<Kernel, Primitive>>;
void checkStop(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Meshing cancelled.");
}
std::array<double, 3> barycentric(const Kernel::Point_3 &point, const Kernel::Point_3 &a,
                                  const Kernel::Point_3 &b, const Kernel::Point_3 &c) {
    const auto u = b - a, w = c - a, q = point - a;
    // Dominant-plane determinants avoid squaring the conditioning of thin faces.
    const auto normal = CGAL::cross_product(u, w);
    const std::array<double, 3> n{normal.x(), normal.y(), normal.z()};
    const auto axis = static_cast<std::size_t>(
        std::max_element(n.begin(), n.end(),
                         [](double x, double y) { return std::abs(x) < std::abs(y); }) -
        n.begin());
    const auto i = (axis + 1) % 3, j = (axis + 2) % 3;
    const std::array<double, 3> uv{u.x(), u.y(), u.z()}, wv{w.x(), w.y(), w.z()},
        qv{q.x(), q.y(), q.z()};
    const double denominator = uv[i] * wv[j] - uv[j] * wv[i];
    if (denominator == 0 || !std::isfinite(denominator))
        throw std::runtime_error("Sizing projection encountered a degenerate source triangle.");
    const double s = (qv[i] * wv[j] - qv[j] * wv[i]) / denominator;
    const double t = (uv[i] * qv[j] - uv[j] * qv[i]) / denominator;
    return {1 - s - t, s, t};
}
} // namespace

struct SurfaceSizingField::Data {
    Mesh mesh;
    std::vector<std::array<std::uint32_t, 3>> triangles;
    std::vector<double> sizes;
    std::map<std::uint32_t, std::unique_ptr<Tree>> patches;
    std::stop_token stop;

    auto closest(Vec3 position, std::span<const std::uint32_t> ids) const {
        const Kernel::Point_3 query(position.x, position.y, position.z);
        std::pair<Kernel::Point_3, Mesh::Face_index> result;
        double best = std::numeric_limits<double>::infinity();
        for (auto id : ids) {
            const auto found = patches.find(id);
            if (found == patches.end())
                throw std::runtime_error("Sizing requested an unknown physical patch.");
            const auto candidate = found->second->closest_point_and_primitive(query);
            const double squared = CGAL::to_double(CGAL::squared_distance(query, candidate.first));
            if (squared < best) {
                best = squared;
                result = candidate;
            }
        }
        if (!std::isfinite(best))
            throw std::runtime_error("Sizing requires a finite point on a physical patch.");
        return result;
    }
};

SurfaceSizingField::SurfaceSizingField(const SurfaceMesh &source, std::vector<double> vertexSizes,
                                       double growthRate, double targetEdgeLength,
                                       std::span<const SurfaceEdgeSeed> edgeSeeds,
                                       std::stop_token stop,
                                       const std::map<std::uint32_t, double> &faceTargetEdgeLengths)
    : data_(std::make_unique<Data>()) {
    auto &data = *data_;
    data.sizes = std::move(vertexSizes);
    data.stop = stop;
    if (data.sizes.size() != source.vertices.size())
        throw std::invalid_argument("Surface sizing seed count does not match the mesh.");
    std::vector<Mesh::Vertex_index> vertices;
    for (const auto p : source.vertices)
        vertices.push_back(data.mesh.add_vertex({p.x, p.y, p.z}));
    auto physical =
        data.mesh.add_property_map<Mesh::Face_index, std::uint32_t>("f:sizing-physical", 0).first;
    for (const auto &triangle : source.triangles) {
        checkStop(stop);
        const auto &v = triangle.vertices;
        const auto face = data.mesh.add_face(vertices[v[0]], vertices[v[1]], vertices[v[2]]);
        if (face == Mesh::null_face())
            throw std::invalid_argument("Sizing requires an oriented source surface.");
        if (face.idx() != data.triangles.size())
            throw std::runtime_error("Source triangle indexing is not contiguous.");
        data.triangles.push_back(v);
        physical[face] = triangle.faceId;
    }
    auto edgeSizes = data.mesh
                         .add_property_map<Mesh::Edge_index, double>(
                             "e:feature-size", std::numeric_limits<double>::infinity())
                         .first;
    for (const auto &seed : edgeSeeds) {
        checkStop(stop);
        if (seed.vertices[0] >= vertices.size() || seed.vertices[1] >= vertices.size() ||
            !std::isfinite(seed.size) || seed.size <= 0)
            throw std::invalid_argument("Invalid source edge sizing seed.");
        const auto found =
            CGAL::edge(vertices[seed.vertices[0]], vertices[seed.vertices[1]], data.mesh);
        if (!found.second)
            throw std::invalid_argument("Sizing seed does not identify a source edge.");
        edgeSizes[found.first] = std::min(edgeSizes[found.first], seed.size);
    }
    auto triangleVertices = [&](Mesh::Face_index f) {
        const auto h = data.mesh.halfedge(f);
        return std::array<Mesh::Vertex_index, 3>{data.mesh.source(h), data.mesh.target(h),
                                                 data.mesh.target(data.mesh.next(h))};
    };
    auto splitEdge = [&](Mesh::Halfedge_index h, const Kernel::Point_3 &position) {
        if (data.mesh.face(h) == Mesh::null_face())
            h = data.mesh.opposite(h);
        const auto opposite = data.mesh.opposite(h);
        const auto left = data.mesh.face(h), right = data.mesh.face(opposite);
        const auto leftId = physical[left];
        const auto rightId = right == Mesh::null_face() ? 0 : physical[right];
        const auto leftCorner = data.mesh.target(data.mesh.next(h));
        const auto rightCorner = right == Mesh::null_face()
                                     ? Mesh::null_vertex()
                                     : data.mesh.target(data.mesh.next(opposite));
        const double inheritedSize = edgeSizes[data.mesh.edge(h)];
        const auto added = CGAL::Euler::split_edge(h, data.mesh);
        edgeSizes[data.mesh.edge(h)] = inheritedSize;
        edgeSizes[data.mesh.edge(added)] = inheritedSize;
        const auto inserted = data.mesh.target(added);
        data.mesh.point(inserted) = position;
        auto triangulate = [&](Mesh::Face_index face, Mesh::Vertex_index corner, std::uint32_t id) {
            if (face == Mesh::null_face())
                return;
            auto toNew = Mesh::null_halfedge(), toCorner = Mesh::null_halfedge();
            for (const auto q : CGAL::halfedges_around_face(data.mesh.halfedge(face), data.mesh)) {
                if (data.mesh.target(q) == inserted)
                    toNew = q;
                if (data.mesh.target(q) == corner)
                    toCorner = q;
            }
            const auto diagonal = CGAL::Euler::split_face(toNew, toCorner, data.mesh);
            physical[data.mesh.face(diagonal)] = id;
            physical[data.mesh.face(data.mesh.opposite(diagonal))] = id;
        };
        triangulate(left, leftCorner, leftId);
        triangulate(right, rightCorner, rightId);
        return std::pair{inserted, inheritedSize};
    };
    // Original CAD fans may have only small boundary seeds on a broad plane.
    // Add unconstrained field samples so their P1 interpolant can recover the
    // requested coarse size away from those boundaries. These subdivisions
    // preserve the source piecewise-planar surface and are not BEM elements.
    if (!std::isfinite(targetEdgeLength) || targetEdgeLength <= 0)
        throw std::invalid_argument("Sizing requires a finite positive target edge length.");
    std::deque<Mesh::Edge_index> pending;
    for (const auto edge : data.mesh.edges())
        pending.push_back(edge);
    const double maximumSquared = targetEdgeLength * targetEdgeLength;
    while (!pending.empty()) {
        checkStop(stop);
        const auto edge = pending.front();
        pending.pop_front();
        const auto h = data.mesh.halfedge(edge);
        const auto a = data.mesh.point(data.mesh.source(h)),
                   b = data.mesh.point(data.mesh.target(h));
        if (CGAL::to_double(CGAL::squared_distance(a, b)) <= maximumSquared)
            continue;
        const auto [inserted, inheritedSize] = splitEdge(h, a + .5 * (b - a));
        if (inserted.idx() != data.sizes.size())
            throw std::runtime_error("Sizing sample vertex indexing is not contiguous.");
        data.sizes.push_back(std::min(targetEdgeLength, inheritedSize));
        for (const auto q : CGAL::halfedges_around_target(data.mesh.halfedge(inserted), data.mesh))
            pending.push_back(data.mesh.edge(q));
    }
    std::map<std::uint32_t, std::vector<Mesh::Face_index>> faces;
    data.triangles.resize(data.mesh.num_faces());
    std::vector<std::vector<std::pair<std::uint32_t, double>>> neighbors(data.mesh.num_vertices());
    for (const auto f : data.mesh.faces()) {
        checkStop(stop);
        const auto v = triangleVertices(f);
        data.triangles[f.idx()] = {v[0].idx(), v[1].idx(), v[2].idx()};
        faces[physical[f]].push_back(f);
        // Seed the complete requested patch, including new field samples.
        // The shared surface graph then grades into neighboring patches.
        if (const auto local = faceTargetEdgeLengths.find(physical[f]);
            local != faceTargetEdgeLengths.end())
            for (const auto vertex : v)
                data.sizes[vertex.idx()] = std::min(data.sizes[vertex.idx()], local->second);
    }
    for (const auto edge : data.mesh.edges()) {
        const auto h = data.mesh.halfedge(edge);
        const auto a = data.mesh.source(h), b = data.mesh.target(h);
        const double length = std::sqrt(
            CGAL::to_double(CGAL::squared_distance(data.mesh.point(a), data.mesh.point(b))));
        neighbors[a.idx()].push_back({b.idx(), length});
        neighbors[b.idx()].push_back({a.idx(), length});
    }
    for (const auto &[id, list] : faces) {
        auto tree = std::make_unique<Tree>(list.begin(), list.end(), data.mesh);
        tree->accelerate_distance_queries();
        data.patches.emplace(id, std::move(tree));
    }
    using Entry = std::pair<double, std::uint32_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
    for (std::size_t i = 0; i < data.sizes.size(); ++i) {
        if (!std::isfinite(data.sizes[i]) || data.sizes[i] <= 0)
            throw std::invalid_argument("Surface sizing requires finite positive seeds.");
        queue.push({data.sizes[i], static_cast<std::uint32_t>(i)});
    }
    while (!queue.empty()) {
        checkStop(stop);
        const auto [size, vertex] = queue.top();
        queue.pop();
        if (size != data.sizes[vertex])
            continue;
        for (const auto [other, length] : neighbors[vertex]) {
            const double candidate = size + growthRate * length;
            if (candidate < data.sizes[other]) {
                data.sizes[other] = candidate;
                queue.push({candidate, other});
            }
        }
    }
}

SurfaceSizingField::~SurfaceSizingField() = default;

double SurfaceSizingField::at(Vec3 position, std::span<const std::uint32_t> ids) const {
    checkStop(data_->stop);
    const auto [projected, face] = data_->closest(position, ids);
    const auto &v = data_->triangles[face.idx()];
    const auto a = data_->mesh.point(Mesh::Vertex_index(v[0]));
    const auto b = data_->mesh.point(Mesh::Vertex_index(v[1]));
    const auto c = data_->mesh.point(Mesh::Vertex_index(v[2]));
    auto weights = barycentric(projected, a, b, c);
    for (auto &weight : weights)
        weight = std::clamp(weight, 0.0, 1.0);
    const double sum = weights[0] + weights[1] + weights[2];
    double result = 0;
    for (std::size_t j = 0; j < 3; ++j)
        result += weights[j] / sum * data_->sizes[v[j]];
    if (!std::isfinite(result) || result <= 0)
        throw std::runtime_error("Surface sizing produced a non-finite or non-positive target.");
    return result;
}
} // namespace hornsculpt::detail
