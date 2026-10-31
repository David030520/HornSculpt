#include "space.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hornsculpt::galerkin {
namespace {
Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
Vec3 divided(Vec3 value, double divisor) {
    return {value.x / divisor, value.y / divisor, value.z / divisor};
}
double length(Vec3 value) {
    return std::hypot(value.x, value.y, value.z);
}
bool finite(Vec3 value) {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}
} // namespace

Space buildSpace(const SurfaceMesh &mesh) {
    if (mesh.vertices.empty() || mesh.triangles.empty())
        throw std::invalid_argument("Galerkin space requires a nonempty surface mesh.");
    if (mesh.vertices.size() > std::numeric_limits<std::uint32_t>::max() ||
        mesh.triangles.size() > std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument("Galerkin mesh exceeds the 32-bit panel index limit.");
    for (const auto vertex : mesh.vertices)
        if (!finite(vertex))
            throw std::invalid_argument("Galerkin mesh has non-finite coordinates.");

    Space result;
    result.vertices = mesh.vertices;
    result.panels.reserve(mesh.triangles.size());
    result.incidenceOffsets.resize(mesh.vertices.size() + 1, 0);
    result.massDiagonal.resize(mesh.vertices.size(), 0.0);
    result.lumpedMass.resize(mesh.vertices.size(), 0.0);
    double areaCompensation = 0.0;
    for (const auto &triangle : mesh.triangles) {
        Panel panel;
        panel.nodes = triangle.vertices;
        panel.faceId = triangle.faceId;
        panel.boundary = triangle.boundary;
        for (std::size_t local = 0; local < 3; ++local) {
            if (panel.nodes[local] >= result.vertices.size())
                throw std::invalid_argument("Galerkin triangle has an invalid vertex index.");
            panel.vertices[local] = result.vertices[panel.nodes[local]];
        }
        const auto a = panel.vertices[0], b = panel.vertices[1], c = panel.vertices[2];
        const auto ab = subtract(b, a), ac = subtract(c, a), bc = subtract(c, b);
        const auto areaVector = cross(ab, ac);
        const double twiceArea = length(areaVector);
        panel.maximumEdgeLength = std::max({length(ab), length(ac), length(bc)});
        // Relative to panel size so unit changes do not alter admissibility.
        const double degeneracyThreshold = 128.0 * std::numeric_limits<double>::epsilon() *
                                           panel.maximumEdgeLength * panel.maximumEdgeLength;
        if (!std::isfinite(twiceArea) || !std::isfinite(panel.maximumEdgeLength) ||
            twiceArea <= degeneracyThreshold)
            throw std::invalid_argument(
                "Galerkin triangle is degenerate or numerically unresolved.");
        panel.area = 0.5 * twiceArea;
        panel.normal = divided(areaVector, twiceArea);
        // This form avoids adding three possibly large translated coordinates.
        panel.centroid = {a.x + (ab.x + ac.x) / 3.0, a.y + (ab.y + ac.y) / 3.0,
                          a.z + (ab.z + ac.z) / 3.0};
        const std::array oppositeEdges{bc, subtract(a, c), ab};
        for (std::size_t local = 0; local < 3; ++local) {
            panel.gradients[local] = divided(cross(panel.normal, oppositeEdges[local]), twiceArea);
            panel.curls[local] = cross(panel.normal, panel.gradients[local]);
            if (!finite(panel.gradients[local]) || !finite(panel.curls[local]))
                throw std::invalid_argument("Galerkin basis gradient is not finite.");
            panel.massP1DP0[local] = panel.area / 3.0;
            for (std::size_t trial = 0; trial < 3; ++trial)
                panel.massP1P1[3 * local + trial] =
                    panel.area * (local == trial ? 2.0 : 1.0) / 12.0;
            ++result.incidenceOffsets[static_cast<std::size_t>(panel.nodes[local]) + 1];
            result.massDiagonal[panel.nodes[local]] += panel.area / 6.0;
            result.lumpedMass[panel.nodes[local]] += panel.area / 3.0;
        }
        const double correctedArea = panel.area - areaCompensation;
        const double nextArea = result.surfaceArea + correctedArea;
        areaCompensation = (nextArea - result.surfaceArea) - correctedArea;
        result.surfaceArea = nextArea;
        result.panels.push_back(std::move(panel));
    }
    for (std::size_t vertex = 0; vertex < mesh.vertices.size(); ++vertex) {
        if (result.incidenceOffsets[vertex + 1] == 0)
            throw std::invalid_argument("Galerkin mesh contains an orphan pressure vertex.");
        result.incidenceOffsets[vertex + 1] += result.incidenceOffsets[vertex];
    }
    result.incidences.resize(result.incidenceOffsets.back());
    auto cursor = result.incidenceOffsets;
    for (std::size_t triangle = 0; triangle < result.panels.size(); ++triangle) {
        const auto &panel = result.panels[triangle];
        for (std::size_t local = 0; local < 3; ++local) {
            result.incidences[cursor[panel.nodes[local]]++] = {static_cast<std::uint32_t>(triangle),
                                                               static_cast<std::uint8_t>(local)};
        }
    }
    return result;
}

} // namespace hornsculpt::galerkin
