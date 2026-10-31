#include "axisymmetric_space.hpp"
#include "../axisymmetric_reduction.hpp"

#include <cmath>
#include <map>
#include <numbers>
#include <stdexcept>

namespace hornsculpt::galerkin {
namespace {
Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 scaled(Vec3 a, double factor) {
    return {a.x * factor, a.y * factor, a.z * factor};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double norm(Vec3 v) {
    return std::hypot(v.x, v.y, v.z);
}
SurfaceMesh meshOf(const Space &space) {
    SurfaceMesh mesh;
    mesh.vertices = space.vertices;
    mesh.triangles.reserve(space.panels.size());
    for (const auto &panel : space.panels)
        mesh.triangles.push_back({panel.nodes, panel.faceId, panel.boundary});
    return mesh;
}
AxisymmetricSpace extract(const SurfaceMesh &mesh, Vec3 worldOrigin,
                          const AxisymmetricBem &reduction) {
    if (reduction.surfaceContract != AxisymmetricSurfaceContract::continuousRevolution)
        throw std::invalid_argument(
            "Continuous rotational space requires its explicit surface contract.");
    const auto layout = bem_detail::validateAxisymmetricReduction(mesh, worldOrigin, reduction);
    AxisymmetricSpace out;
    auto &m = out.meridian;
    m.origin = reduction.origin;
    m.axis = scaled(reduction.axis, 1 / norm(reduction.axis));
    m.points.resize(layout.pressureCount);
    out.elements.resize(layout.fluxCount);
    out.lumpedMass.resize(layout.pressureCount);
    const auto localOrigin = subtract(reduction.origin, worldOrigin);
    std::vector<std::size_t> representatives(layout.pressureCount, mesh.vertices.size()),
        orbitSizes(layout.pressureCount);
    for (std::size_t vertex = 0; vertex < mesh.vertices.size(); ++vertex) {
        const auto dof = reduction.pressureDofs[vertex];
        if (!orbitSizes[dof]++)
            representatives[dof] = vertex;
    }
    bool radialFound = false;
    for (std::size_t dof = 0; dof < layout.pressureCount; ++dof) {
        const auto d = subtract(mesh.vertices[representatives[dof]], localOrigin);
        const double z = dot(d, m.axis);
        const auto radial = subtract(d, scaled(m.axis, z));
        const double radius = orbitSizes[dof] == 1 ? 0 : norm(radial);
        m.points[dof] = {radius, z};
        if (radius > 0 && !radialFound) {
            m.radial = scaled(radial, 1 / radius);
            radialFound = true;
        }
    }
    if (!radialFound)
        throw std::invalid_argument("Rotational space has no off-axis meridian point.");
    std::vector<bool> found(layout.fluxCount);
    for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
        const auto &triangle = mesh.triangles[t];
        const auto element = reduction.fluxDofs[t];
        std::uint32_t first{}, second{};
        bool oriented = false;
        for (unsigned edge = 0; edge < 3; ++edge) {
            const auto a = triangle.vertices[edge], b = triangle.vertices[(edge + 1) % 3],
                       c = triangle.vertices[(edge + 2) % 3];
            const auto same = reduction.pressureDofs[a], other = reduction.pressureDofs[c];
            if (same != reduction.pressureDofs[b] || same == other)
                continue;
            const auto relativeA = subtract(mesh.vertices[a], localOrigin),
                       relativeB = subtract(mesh.vertices[b], localOrigin);
            const auto radialA = subtract(relativeA, scaled(m.axis, dot(relativeA, m.axis))),
                       radialB = subtract(relativeB, scaled(m.axis, dot(relativeB, m.axis)));
            const auto unitA = scaled(radialA, 1 / norm(radialA)),
                       unitB = scaled(radialB, 1 / norm(radialB));
            const double turn = dot(m.axis, cross(unitA, unitB));
            if (!std::isfinite(turn) || turn == 0)
                throw std::invalid_argument("Rotational triangle has no directed angular edge.");
            first = turn > 0 ? same : other;
            second = turn > 0 ? other : same;
            oriented = true;
            break;
        }
        if (!oriented)
            throw std::invalid_argument("Rotational triangle does not span one meridian segment.");
        if (found[element]) {
            const auto &e = out.elements[element];
            if (e.first != first || e.second != second || e.faceId != triangle.faceId ||
                e.boundary != triangle.boundary)
                throw std::invalid_argument(
                    "Rotational strip orientation or physical tags disagree.");
        } else {
            found[element] = true;
            auto &e = out.elements[element];
            e.first = first;
            e.second = second;
            e.faceId = triangle.faceId;
            e.boundary = triangle.boundary;
        }
    }
    std::map<std::uint32_t, BoundaryCondition> faceConditions;
    std::vector<unsigned> incoming(layout.pressureCount), outgoing(layout.pressureCount);
    std::vector<std::vector<std::uint32_t>> neighbours(layout.pressureCount);
    long double volume = 0;
    for (const auto &e : out.elements) {
        if (!e.faceId ||
            (e.boundary != BoundaryCondition::wall && e.boundary != BoundaryCondition::driven))
            throw std::invalid_argument("Rotational element has an invalid physical face.");
        const auto [entry, inserted] = faceConditions.emplace(e.faceId, e.boundary);
        if (!inserted && entry->second != e.boundary)
            throw std::invalid_argument(
                "Rotational face group has inconsistent boundary conditions.");
        ++outgoing[e.first];
        ++incoming[e.second];
        neighbours[e.first].push_back(e.second);
        neighbours[e.second].push_back(e.first);
    }
    for (std::size_t dof = 0; dof < layout.pressureCount; ++dof)
        if (m.points[dof].radius == 0 ? incoming[dof] + outgoing[dof] != 1
                                      : incoming[dof] != 1 || outgoing[dof] != 1)
            throw std::invalid_argument(
                "Rotational meridian is open, branching or inconsistently oriented.");
    std::vector<bool> visited(layout.pressureCount);
    std::vector<std::uint32_t> pending{0};
    visited[0] = true;
    for (std::size_t index = 0; index < pending.size(); ++index)
        for (auto next : neighbours[pending[index]])
            if (!visited[next]) {
                visited[next] = true;
                pending.push_back(next);
            }
    if (pending.size() != layout.pressureCount)
        throw std::invalid_argument("Rotational meridian has disconnected components.");
    for (auto &e : out.elements) {
        const auto a = m.points[e.first], b = m.points[e.second];
        const double dr = b.radius - a.radius, dz = b.axial - a.axial,
                     radiusSum = a.radius + b.radius;
        e.length = std::hypot(dr, dz);
        e.area = std::numbers::pi * radiusSum * e.length;
        if (!(e.length > 0) || !(e.area > 0) || !std::isfinite(e.length) || !std::isfinite(e.area))
            throw std::invalid_argument(
                "Rotational meridian has a degenerate or non-finite element.");
        e.radialNormal = dz / e.length;
        e.axialNormal = -dr / e.length;
        const double u = a.radius / radiusSum, v = b.radius / radiusSum;
        // Integrate N_i N_j and N_i against 2*pi*L*r(s) ds exactly.
        e.massP1P1 = {e.area * ((1 + 2 * u) / 6), e.area / 6, e.area / 6,
                      e.area * ((1 + 2 * v) / 6)};
        e.massP1DP0 = {e.area * ((1 + u) / 3), e.area * ((1 + v) / 3)};
        out.lumpedMass[e.first] += e.massP1DP0[0];
        out.lumpedMass[e.second] += e.massP1DP0[1];
        m.segments.push_back({e.first, e.second, e.faceId, e.boundary});
        const long double r0 = a.radius, r1 = b.radius;
        volume += static_cast<long double>(dz) * (r0 * r0 + r0 * r1 + r1 * r1);
    }
    if (!(volume > 0) || !std::isfinite(volume))
        throw std::invalid_argument("Rotational air volume is not outward oriented.");
    for (double mass : out.lumpedMass)
        if (!(mass > 0) || !std::isfinite(mass))
            throw std::invalid_argument("Rotational pressure mass is not finite and positive.");
    return out;
}
} // namespace

AxisymmetricSpace buildAxisymmetricSpace(const Space &space, Vec3 worldOrigin,
                                         const AxisymmetricBem &reduction) {
    return extract(meshOf(space), worldOrigin, reduction);
}
AxisymmetricSpace buildAxisymmetricSpace(const SurfaceMesh &mesh, Vec3 worldOrigin,
                                         const AxisymmetricBem &reduction) {
    // buildSpace establishes the same finite/nondegenerate P1 mesh contract as
    // the overload consumed by the numerical assembly.
    (void)buildSpace(mesh);
    return extract(mesh, worldOrigin, reduction);
}
} // namespace hornsculpt::galerkin
