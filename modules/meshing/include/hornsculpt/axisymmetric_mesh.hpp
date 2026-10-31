#pragma once

#include <hornsculpt/axisymmetric.hpp>
#include <stop_token>

namespace hornsculpt {

struct AxisymmetricMesh {
    SurfaceMesh mesh;
    std::vector<std::uint32_t> pressureDofs;
    std::vector<std::uint32_t> fluxDofs;
    std::vector<std::uint32_t> testPanels;
};

// Integration tessellation for the m=0 Galerkin space, not a 3D isotropic
// remeshing result. Unknowns live only on the generator. All sector copies
// retain identical physical face IDs and boundary conditions.
AxisymmetricMesh buildAxisymmetricMesh(const AxisymmetricMeridian &meridian, unsigned sectors = 64,
                                       std::stop_token stop = {});

} // namespace hornsculpt
