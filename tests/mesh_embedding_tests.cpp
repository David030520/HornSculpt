#include <hornsculpt/meshing.hpp>
#include <iostream>
#include <stdexcept>

using namespace hornsculpt;

namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
SurfaceMesh tetrahedron(Vec3 shift = {}) {
    SurfaceMesh mesh;
    mesh.vertices = {{shift.x, shift.y, shift.z},
                     {shift.x + 1, shift.y, shift.z},
                     {shift.x, shift.y + 1, shift.z},
                     {shift.x, shift.y, shift.z + 1}};
    mesh.triangles = {{{0, 2, 1}, 1}, {{0, 1, 3}, 1}, {{0, 3, 2}, 1}, {{1, 2, 3}, 1}};
    return mesh;
}
void append(SurfaceMesh &mesh, const SurfaceMesh &other) {
    const auto offset = static_cast<std::uint32_t>(mesh.vertices.size());
    mesh.vertices.insert(mesh.vertices.end(), other.vertices.begin(), other.vertices.end());
    for (auto triangle : other.triangles) {
        for (auto &v : triangle.vertices)
            v += offset;
        mesh.triangles.push_back(triangle);
    }
}
} // namespace

int main() {
    try {
        auto mesh = tetrahedron();
        auto q = inspectMesh(mesh);
        require(q.watertight() && q.selfIntersectionPairs == 0, "Valid tetrahedron rejected");
        append(mesh, tetrahedron({2, 0, 0}));
        q = inspectMesh(mesh);
        require(q.watertight() && q.connectedComponents == 2, "Separate closed bodies rejected");
        mesh = tetrahedron();
        append(mesh, tetrahedron({0.2, 0.2, 0.2}));
        q = inspectMesh(mesh);
        require(q.boundaryEdges == 0 && q.selfIntersectionPairs > 0 && !q.watertight(),
                "Crossing closed bodies accepted");
        mesh = tetrahedron();
        mesh.vertices.insert(mesh.vertices.end(), {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}});
        mesh.triangles.insert(mesh.triangles.end(),
                              {{{0, 4, 5}, 2}, {{0, 6, 4}, 2}, {{0, 5, 6}, 2}, {{4, 6, 5}, 2}});
        q = inspectMesh(mesh);
        require(q.boundaryEdges == 0 && q.nonManifoldEdges == 0 && q.nonManifoldVertices == 1 &&
                    !q.watertight(),
                "Bow-tie vertex accepted");
        mesh.vertices = {{0, 0, 0}, {2, 0, 0}, {0, 2, 0}, {1, 0, 0}, {2, 1, 0}, {0, 1, 0}};
        mesh.triangles = {{{0, 1, 2}, 1}, {{3, 4, 5}, 2}};
        require(inspectMesh(mesh).selfIntersectionPairs == 1, "Coplanar overlap missed");
        mesh.vertices = {{0, 0, 0}, {1, 0, 0}, {1, 1, 1e-12}, {1, 1 + 1e-7, 1e-12}, {0, 1, 0}};
        mesh.triangles = {{{0, 1, 2}, 1}, {{0, 3, 4}, 1}};
        require(inspectMesh(mesh).selfIntersectionPairs == 0,
                "Nearly coplanar shared-vertex contact falsely reported as a crossing");
        mesh = tetrahedron({1e5, -1e5, 1e5});
        require(inspectMesh(mesh).watertight(), "Translated valid body rejected");
        for (auto &p : mesh.vertices) {
            p.x *= 1e-7;
            p.y *= 1e-7;
            p.z *= 1e-7;
        }
        require(inspectMesh(mesh).watertight(), "Scaled valid body rejected");
        std::cout << "Embedding validation: closed, crossing, coplanar, bow-tie, scale and "
                     "translation passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
