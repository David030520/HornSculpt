#pragma once

#include "space.hpp"
#include <array>
#include <hornsculpt/axisymmetric.hpp>
#include <hornsculpt/bem.hpp>

namespace hornsculpt::galerkin {
struct AxisymmetricElement {
    std::uint32_t first{}, second{}, faceId{};
    BoundaryCondition boundary{BoundaryCondition::wall};
    double length{};
    double radialNormal{}, axialNormal{};
    double area{};
    // Full 2*pi integration, local endpoint order first/second, in m^2.
    std::array<double, 4> massP1P1{};
    std::array<double, 2> massP1DP0{};
};

struct AxisymmetricSpace {
    // World-space frame; points are (radius, axial) relative to that frame.
    // Point/element indices equal the public pressure/flux DOF indices.
    AxisymmetricMeridian meridian;
    std::vector<AxisymmetricElement> elements;
    std::vector<double> lumpedMass;
};

// Extracts an oriented piecewise linear generator from the already represented
// periodic surface; it never fits or obtains another CAD source. The caller's
// space must come from buildSpace. local coordinates + worldOrigin are global.
// Requires the continuous-revolution contract; old polygonal traces are read
// through surfaceMeasures rather than passed to a second numerical solver.
AxisymmetricSpace buildAxisymmetricSpace(const Space &space, Vec3 worldOrigin,
                                         const AxisymmetricBem &reduction);
AxisymmetricSpace buildAxisymmetricSpace(const SurfaceMesh &mesh, Vec3 worldOrigin,
                                         const AxisymmetricBem &reduction);
} // namespace hornsculpt::galerkin
