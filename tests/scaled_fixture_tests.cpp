#include "scaled_fixture_geometry.hpp"
#include <cmath>
#include <hornsculpt/meshing.hpp>
#include <iostream>
#include <stdexcept>

using namespace hornsculpt;
int main() {
    try {
        auto polyhedron = testing::sphere(0);
        for (std::size_t i = 0; i < polyhedron.triangles.size(); ++i) {
            polyhedron.triangles[i].faceId = static_cast<std::uint32_t>(i + 1);
            polyhedron.triangles[i].boundary =
                i % 2 ? BoundaryCondition::wall : BoundaryCondition::driven;
        }
        const auto originalQuality = inspectMesh(polyhedron);
        const auto refined = testing::refinePolyhedralFixture(polyhedron, 500);
        const auto refinedQuality = inspectMesh(refined);
        if (!refinedQuality.watertight() || refinedQuality.selfIntersectionPairs != 0 ||
            std::abs(refinedQuality.signedVolume / originalQuality.signedVolume - 1) > 1e-13)
            throw std::runtime_error(
                "Polyhedral refinement changed the enclosed volume or topology.");
        std::map<std::uint32_t, double> area;
        for (const auto &triangle : refined.triangles) {
            const auto &parent = polyhedron.triangles.at(triangle.faceId - 1);
            if (triangle.boundary != parent.boundary)
                throw std::runtime_error("Polyhedral refinement changed a boundary condition.");
            const auto a = polyhedron.vertices[parent.vertices[0]],
                       b = polyhedron.vertices[parent.vertices[1]],
                       c = polyhedron.vertices[parent.vertices[2]];
            const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
            const Vec3 normal{u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
            for (const auto vertex : triangle.vertices) {
                const auto p = refined.vertices[vertex];
                if (std::abs(normal.x * (p.x - a.x) + normal.y * (p.y - a.y) +
                             normal.z * (p.z - a.z)) > 1e-16)
                    throw std::runtime_error(
                        "A refined corner left its original polyhedral facet.");
            }
            const auto x = refined.vertices[triangle.vertices[0]],
                       y = refined.vertices[triangle.vertices[1]],
                       z = refined.vertices[triangle.vertices[2]];
            const Vec3 xy{y.x - x.x, y.y - x.y, y.z - x.z}, xz{z.x - x.x, z.y - x.y, z.z - x.z};
            area[triangle.faceId] +=
                .5 * std::hypot(xy.y * xz.z - xy.z * xz.y, xy.z * xz.x - xy.x * xz.z,
                                xy.x * xz.y - xy.y * xz.x);
        }
        for (const auto &[face, value] : area) {
            (void)face;
            if (std::abs(value / (.005 * std::sqrt(3.0)) - 1) > 1e-13)
                throw std::runtime_error("Polyhedral refinement changed a physical face area.");
        }
        for (const auto count : {50000u, 100000u}) {
            const auto mesh = testing::scaledSphere(count);
            const auto quality = inspectMesh(mesh);
            if (mesh.triangles.size() != count || mesh.vertices.size() != count / 2 + 2 ||
                !quality.watertight() || quality.selfIntersectionPairs != 0 ||
                quality.connectedComponents != 1 || quality.eulerCharacteristic != 2 ||
                quality.signedVolume <= 0 || quality.minimumAngleDegrees < 20 ||
                quality.maximumAspectRatio > 4)
                throw std::runtime_error("Scaled sphere violates the benchmark mesh contract.");
            for (const auto p : mesh.vertices)
                if (std::abs(std::hypot(p.x, p.y, p.z) - .1) > 1e-14)
                    throw std::runtime_error("Scaled sphere vertex is not on the analytic sphere.");
            std::cout << count << "T " << mesh.vertices.size()
                      << "V min-angle=" << quality.minimumAngleDegrees
                      << " aspect=" << quality.maximumAspectRatio
                      << " intersections=" << quality.selfIntersectionPairs << '\n';
        }
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
