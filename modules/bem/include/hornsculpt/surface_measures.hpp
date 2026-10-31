#pragma once

#include <array>
#include <hornsculpt/bem.hpp>

namespace hornsculpt {
struct SurfaceMeasures {
    std::vector<double> areas;
    // Exact integrals of the three local nodal basis functions, in the
    // original mesh triangle's node order. Weights sum to areas[t], in m^2.
    std::vector<std::array<double, 3>> pressureIntegrals;
};

// Cartesian and archived triangulated-orbit traces use planar triangle area
// and A/3. Continuous m=0 uses the exact revolved generator, with each displayed
// ring triangle interpreted as its (s,phi) half-cell; pole fans use their full
// angular sector. This preserves the radial P1 pressure/volume-flow duality.
// mesh coordinates + worldOrigin are global; omission means world coordinates.
SurfaceMeasures surfaceMeasures(const SurfaceMesh &mesh,
                                const std::optional<AxisymmetricBem> &axisymmetric = {},
                                Vec3 worldOrigin = {});
} // namespace hornsculpt
