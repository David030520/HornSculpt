#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <galerkin/exterior_domain.hpp>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace hornsculpt;
using hornsculpt::galerkin::validateExteriorComponents;

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

void appendTetrahedron(SurfaceMesh &mesh, Vec3 origin, double size, bool inverted = false) {
    const auto offset = static_cast<std::uint32_t>(mesh.vertices.size());
    mesh.vertices.push_back(origin);
    mesh.vertices.push_back({origin.x + size, origin.y, origin.z});
    mesh.vertices.push_back({origin.x, origin.y + size, origin.z});
    mesh.vertices.push_back({origin.x, origin.y, origin.z + size});
    for (auto face : std::array<std::array<std::uint32_t, 3>, 4>{
             {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}}}) {
        for (auto &node : face)
            node += offset;
        if (inverted)
            std::swap(face[1], face[2]);
        mesh.triangles.push_back({face, 1, BoundaryCondition::wall});
    }
}

void translate(SurfaceMesh &mesh, Vec3 value) {
    for (auto &vertex : mesh.vertices) {
        vertex.x += value.x;
        vertex.y += value.y;
        vertex.z += value.z;
    }
}

double globalSignedVolume(const SurfaceMesh &mesh) {
    double sixVolumes = 0;
    for (const auto &face : mesh.triangles) {
        const auto a = mesh.vertices[face.vertices[0]], b = mesh.vertices[face.vertices[1]],
                   c = mesh.vertices[face.vertices[2]];
        sixVolumes += a.x * (b.y * c.z - b.z * c.y) + a.y * (b.z * c.x - b.x * c.z) +
                      a.z * (b.x * c.y - b.y * c.x);
    }
    return sixVolumes / 6;
}

void rejected(const SurfaceMesh &mesh, const std::string &reason) {
    try {
        validateExteriorComponents(mesh);
    } catch (const std::invalid_argument &error) {
        require(std::string(error.what()).find(reason) != std::string::npos,
                "Exterior-domain validation failed for an unrelated reason");
        return;
    }
    throw std::runtime_error("Invalid exterior-domain surface was accepted");
}

void checkSeparateBodies() {
    SurfaceMesh mesh;
    appendTetrahedron(mesh, {0, 0, 0}, 4);
    validateExteriorComponents(mesh);
    appendTetrahedron(mesh, {10, 0, 0}, .5);
    validateExteriorComponents(mesh);
    // Its bounding box lies inside the first body's bounding box, but this
    // component is outside the tetrahedron. A bounding-box test is insufficient.
    appendTetrahedron(mesh, {3, 3, 3}, .5);
    validateExteriorComponents(mesh);
    std::reverse(mesh.triangles.begin(), mesh.triangles.end());
    validateExteriorComponents(mesh);
    translate(mesh, {1e12, -1e12, 1e12});
    validateExteriorComponents(mesh);
}

void checkInvertedComponents() {
    SurfaceMesh onlyInverted;
    appendTetrahedron(onlyInverted, {0, 0, 0}, 1, true);
    rejected(onlyInverted, "outward");

    // Total signed volume is positive in both cases: (4^3-.5^3)/6.
    SurfaceMesh separated;
    appendTetrahedron(separated, {0, 0, 0}, 4);
    appendTetrahedron(separated, {10, 0, 0}, .5, true);
    require(globalSignedVolume(separated) > 0,
            "Inverted-component regression fixture must have positive total volume");
    rejected(separated, "outward");
    translate(separated, {1e12, -1e12, 1e12});
    rejected(separated, "outward");

    SurfaceMesh cavity;
    appendTetrahedron(cavity, {0, 0, 0}, 4);
    appendTetrahedron(cavity, {.5, .5, .5}, .5, true);
    require(globalSignedVolume(cavity) > 0,
            "Cavity regression fixture must have positive total volume");
    rejected(cavity, "outward");
}

void checkNestedOutwardComponents() {
    SurfaceMesh nested;
    appendTetrahedron(nested, {0, 0, 0}, 4);
    appendTetrahedron(nested, {.5, .5, .5}, .5);
    rejected(nested, "Nested");
    std::reverse(nested.triangles.begin(), nested.triangles.end());
    rejected(nested, "Nested");
    translate(nested, {1e12, -1e12, 1e12});
    rejected(nested, "Nested");

    // A third body outside both shells does not make the nested pair admissible.
    appendTetrahedron(nested, {1e12 + 10, -1e12, 1e12}, 1);
    rejected(nested, "Nested");
}

void checkDistantPlacements() {
    SurfaceMesh mesh;
    appendTetrahedron(mesh, {1e12, -1e12, 1e12}, 1);
    appendTetrahedron(mesh, {-1e12, 1e12, -1e12}, .5);
    validateExteriorComponents(mesh);
    auto inverted = mesh;
    for (std::size_t panel = 4; panel < inverted.triangles.size(); ++panel)
        std::swap(inverted.triangles[panel].vertices[1], inverted.triangles[panel].vertices[2]);
    rejected(inverted, "outward");
}
} // namespace

int main() {
    try {
        checkSeparateBodies();
        checkInvertedComponents();
        checkNestedOutwardComponents();
        checkDistantPlacements();
        std::cout << "Exterior component orientation, nesting and CAD-placement tests passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
