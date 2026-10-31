#include "axisymmetric_reduction.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>

namespace hornsculpt::bem_detail {
namespace {
Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3 scaled(Vec3 a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
double length(Vec3 a) {
    return std::hypot(a.x, a.y, a.z);
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
bool finite(Vec3 a) {
    return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z);
}
std::array<std::uint32_t, 3> cyclic(std::array<std::uint32_t, 3> nodes) {
    if (nodes[1] < nodes[0] && nodes[1] < nodes[2])
        return {nodes[1], nodes[2], nodes[0]};
    if (nodes[2] < nodes[0] && nodes[2] < nodes[1])
        return {nodes[2], nodes[0], nodes[1]};
    return nodes;
}
std::size_t dofCount(const std::vector<std::uint32_t> &map, std::size_t expected) {
    if (map.empty() || map.size() != expected)
        throw std::invalid_argument("Axisymmetric DOF map does not match the mesh.");
    const auto maximum = *std::max_element(map.begin(), map.end());
    if (maximum >= expected)
        throw std::invalid_argument("Axisymmetric DOF index exceeds the mesh size.");
    std::vector<bool> present(static_cast<std::size_t>(maximum) + 1);
    for (auto index : map)
        present[index] = true;
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("Axisymmetric DOF indices must be contiguous from zero.");
    return present.size();
}
} // namespace

AxisymmetricLayout validateAxisymmetricReduction(const SurfaceMesh &mesh, Vec3 worldOrigin,
                                                 const AxisymmetricBem &reduction) {
    if (!finite(worldOrigin) || !finite(reduction.origin) || !finite(reduction.axis) ||
        mesh.vertices.size() > std::numeric_limits<std::uint32_t>::max() ||
        mesh.triangles.size() > std::numeric_limits<std::uint32_t>::max() ||
        std::abs(length(reduction.axis) - 1) > 1e-12 || reduction.sectors < 8 ||
        reduction.sectors % 4 != 0 || reduction.sectors > mesh.vertices.size() ||
        reduction.testPanels.empty() || mesh.triangles.size() % reduction.sectors != 0 ||
        reduction.testPanels.size() != mesh.triangles.size() / reduction.sectors)
        throw std::invalid_argument("Invalid axisymmetric frame, sector count or target sector.");
    const AxisymmetricLayout layout{dofCount(reduction.pressureDofs, mesh.vertices.size()),
                                    dofCount(reduction.fluxDofs, mesh.triangles.size())};
    const auto origin = subtract(reduction.origin, worldOrigin), axis = reduction.axis;
    double extent{};
    for (auto p : mesh.vertices) {
        if (!finite(p))
            throw std::invalid_argument("Axisymmetric mesh has a non-finite vertex.");
        extent = std::max(extent, length(subtract(p, origin)));
    }
    // Coordinate roundoff only: an explicitly generated periodic mesh must
    // close to FP64 accuracy. World magnitudes account for the centering loss.
    const double tolerance = 1024 * std::numeric_limits<double>::epsilon() *
                             std::max({extent, length(worldOrigin), length(reduction.origin)});
    if (!(tolerance > 0) || !std::isfinite(tolerance))
        throw std::invalid_argument("Axisymmetric mesh has an invalid coordinate scale.");
    std::vector<std::vector<std::uint32_t>> orbits(layout.pressureCount);
    for (std::uint32_t vertex = 0; vertex < mesh.vertices.size(); ++vertex)
        orbits[reduction.pressureDofs[vertex]].push_back(vertex);
    Vec3 horizontal{};
    for (auto p : mesh.vertices) {
        const auto d = subtract(p, origin), radial = subtract(d, scaled(axis, dot(d, axis)));
        if (length(radial) > tolerance) {
            horizontal = scaled(radial, 1 / length(radial));
            break;
        }
    }
    if (length(horizontal) == 0)
        throw std::invalid_argument("Axisymmetric mesh contains no off-axis vertices.");
    const auto vertical = cross(axis, horizontal);
    constexpr auto missing = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::vector<std::uint32_t>> rotation(
        layout.pressureCount, std::vector<std::uint32_t>(reduction.sectors, missing));
    std::vector<unsigned> vertexSector(mesh.vertices.size());
    std::vector<bool> pole(layout.pressureCount);
    for (std::size_t dof = 0; dof < orbits.size(); ++dof) {
        const auto first = subtract(mesh.vertices[orbits[dof].front()], origin);
        const double z = dot(first, axis), radius = length(subtract(first, scaled(axis, z)));
        pole[dof] = radius <= tolerance;
        if (orbits[dof].size() != (pole[dof] ? 1 : reduction.sectors))
            throw std::invalid_argument(
                "Axisymmetric pressure orbit has incomplete angular coverage.");
        for (auto vertex : orbits[dof]) {
            const auto d = subtract(mesh.vertices[vertex], origin);
            double angle = std::atan2(dot(d, vertical), dot(d, horizontal));
            if (angle < 0)
                angle += 2 * std::numbers::pi;
            const auto sector = static_cast<unsigned>(std::llround(angle * reduction.sectors /
                                                                   (2 * std::numbers::pi))) %
                                reduction.sectors;
            const double theta = 2 * std::numbers::pi * sector / reduction.sectors;
            const auto expected =
                add(scaled(axis, z), add(scaled(horizontal, radius * std::cos(theta)),
                                         scaled(vertical, radius * std::sin(theta))));
            if (length(subtract(d, expected)) > tolerance ||
                (!pole[dof] && rotation[dof][sector] != missing))
                throw std::invalid_argument(
                    "Axisymmetric vertex orbit is not an exact periodic rotation.");
            vertexSector[vertex] = sector;
            if (pole[dof])
                std::fill(rotation[dof].begin(), rotation[dof].end(), vertex);
            else
                rotation[dof][sector] = vertex;
        }
    }
    std::map<std::array<std::uint32_t, 3>, std::uint32_t> triangles;
    for (std::uint32_t t = 0; t < mesh.triangles.size(); ++t) {
        for (auto vertex : mesh.triangles[t].vertices)
            if (vertex >= mesh.vertices.size())
                throw std::invalid_argument("Axisymmetric triangle has an invalid vertex.");
        if (!triangles.emplace(cyclic(mesh.triangles[t].vertices), t).second)
            throw std::invalid_argument("Axisymmetric mesh repeats an oriented triangle.");
    }
    std::vector<bool> covered(mesh.triangles.size());
    std::vector<std::set<std::uint32_t>> stripNodes(layout.fluxCount);
    std::vector<unsigned> stripTargets(layout.fluxCount);
    unsigned targetSector = reduction.sectors;
    for (auto panel : reduction.testPanels) {
        if (panel >= mesh.triangles.size())
            throw std::invalid_argument("Axisymmetric target panel index is invalid.");
        const auto &test = mesh.triangles[panel];
        const auto strip = reduction.fluxDofs[panel];
        ++stripTargets[strip];
        std::set<unsigned> angles;
        for (auto vertex : test.vertices) {
            const auto dof = reduction.pressureDofs[vertex];
            stripNodes[strip].insert(dof);
            if (!pole[dof])
                angles.insert(vertexSector[vertex]);
        }
        if (angles.size() != 2)
            throw std::invalid_argument("Axisymmetric target triangle does not span one sector.");
        auto first = *angles.begin(), last = *angles.rbegin();
        if (first == 0 && last == reduction.sectors - 1)
            first = last;
        else if (last != first + 1)
            throw std::invalid_argument("Axisymmetric target triangle skips angular sectors.");
        if (targetSector == reduction.sectors)
            targetSector = first;
        if (targetSector != first)
            throw std::invalid_argument(
                "Axisymmetric targets must belong to the same angular sector.");
        for (unsigned sector = 0; sector < reduction.sectors; ++sector) {
            std::array<std::uint32_t, 3> nodes{};
            for (unsigned i = 0; i < 3; ++i) {
                const auto vertex = test.vertices[i];
                nodes[i] = rotation[reduction.pressureDofs[vertex]]
                                   [(vertexSector[vertex] + sector) % reduction.sectors];
            }
            const auto found = triangles.find(cyclic(nodes));
            if (found == triangles.end() || covered[found->second] ||
                reduction.fluxDofs[found->second] != strip ||
                mesh.triangles[found->second].faceId != test.faceId ||
                mesh.triangles[found->second].boundary != test.boundary)
                throw std::invalid_argument(
                    "Axisymmetric rotated triangles, strip tags or normals disagree.");
            covered[found->second] = true;
        }
    }
    for (std::size_t strip = 0; strip < layout.fluxCount; ++strip)
        if (stripNodes[strip].size() != 2 || stripTargets[strip] < 1 || stripTargets[strip] > 2)
            throw std::invalid_argument(
                "Axisymmetric flux DOF must span one complete meridian strip.");
    if (std::find(covered.begin(), covered.end(), false) != covered.end())
        throw std::invalid_argument(
            "Axisymmetric target sector does not cover the full rotation mesh.");
    return layout;
}
void validateAxisymmetricTrace(const AxisymmetricBem &reduction,
                               std::span<const std::complex<double>> pressure,
                               std::span<const std::complex<double>> normalDerivative) {
    const auto check = [](const std::vector<std::uint32_t> &dofs,
                          std::span<const std::complex<double>> values) {
        if (values.empty() || values.size() != dofs.size())
            throw std::invalid_argument("Axisymmetric trace dimensions disagree with its DOFs.");
        std::vector<std::complex<double>> representatives(values.size());
        std::vector<bool> seen(values.size());
        for (std::size_t i = 0; i < values.size(); ++i) {
            const auto id = dofs[i];
            if (id >= values.size() || !std::isfinite(values[i].real()) ||
                !std::isfinite(values[i].imag()))
                throw std::invalid_argument("Invalid axisymmetric trace value or DOF.");
            if (seen[id] && representatives[id] != values[i])
                throw std::invalid_argument(
                    "Expanded axisymmetric trace aliases contain different values.");
            representatives[id] = values[i];
            seen[id] = true;
        }
    };
    check(reduction.pressureDofs, pressure);
    check(reduction.fluxDofs, normalDerivative);
}
} // namespace hornsculpt::bem_detail
