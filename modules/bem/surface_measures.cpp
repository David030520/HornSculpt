#include "axisymmetric_reduction.hpp"
#include "galerkin/axisymmetric_space.hpp"
#include "galerkin/space.hpp"
#include <hornsculpt/surface_measures.hpp>

#include <cmath>
#include <stdexcept>

namespace hornsculpt {
SurfaceMeasures surfaceMeasures(const SurfaceMesh &mesh,
                                const std::optional<AxisymmetricBem> &axisymmetric,
                                Vec3 worldOrigin) {
    if (!std::isfinite(worldOrigin.x) || !std::isfinite(worldOrigin.y) ||
        !std::isfinite(worldOrigin.z))
        throw std::invalid_argument("Surface measures require a finite coordinate origin.");
    const auto cartesian = galerkin::buildSpace(mesh);
    SurfaceMeasures out;
    out.areas.resize(mesh.triangles.size());
    out.pressureIntegrals.resize(mesh.triangles.size());
    if (!axisymmetric ||
        axisymmetric->surfaceContract == AxisymmetricSurfaceContract::triangulatedOrbit) {
        if (axisymmetric)
            bem_detail::validateAxisymmetricReduction(mesh, worldOrigin, *axisymmetric);
        for (std::size_t t = 0; t < cartesian.panels.size(); ++t) {
            out.areas[t] = cartesian.panels[t].area;
            out.pressureIntegrals[t] = cartesian.panels[t].massP1DP0;
        }
        return out;
    }
    const auto space = galerkin::buildAxisymmetricSpace(cartesian, worldOrigin, *axisymmetric);
    for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
        const auto &triangle = mesh.triangles[t];
        const auto &e = space.elements[axisymmetric->fluxDofs[t]];
        const auto a = space.meridian.points[e.first], b = space.meridian.points[e.second];
        const double sectorArea = e.area / axisymmetric->sectors;
        auto &weights = out.pressureIntegrals[t];
        if (a.radius == 0 || b.radius == 0) {
            out.areas[t] = sectorArea;
            weights.fill(sectorArea / 3);
            continue;
        }
        unsigned atFirst = 0;
        for (auto vertex : triangle.vertices)
            atFirst += axisymmetric->pressureDofs[vertex] == e.first ? 1 : 0;
        if (atFirst != 1 && atFirst != 2)
            throw std::invalid_argument(
                "Continuous surface cell has invalid endpoint multiplicities.");
        const double ratio = (atFirst == 2 ? a.radius : b.radius) / (a.radius + b.radius);
        // Integrate on the two affine halves of [0,1] x [0,dphi],
        // with physical Jacobian L*r(s); the diagonal is s+phi/dphi=1.
        out.areas[t] = sectorArea * ((1 + ratio) / 3);
        for (unsigned local = 0; local < 3; ++local) {
            const bool repeated =
                (axisymmetric->pressureDofs[triangle.vertices[local]] == e.first) == (atFirst == 2);
            weights[local] = repeated ? sectorArea * ((1 + 2 * ratio) / 12) : sectorArea / 6;
        }
    }
    return out;
}
} // namespace hornsculpt
