#pragma once
#include "galerkin/space.hpp"
#include <hornsculpt/bem.hpp>
#include <span>

namespace hornsculpt::bem_detail {
// All coordinates are in the same local frame. cross[e,a] and self[e,j]
// are double surface integrals of G divided by the exterior test-cell area.
struct BaffleExteriorOperator {
    std::vector<std::complex<double>> cross;
    std::vector<std::complex<double>> self;
};
void validateBaffleExteriorMesh(const SurfaceMesh &cavity, const SurfaceMesh &exterior,
                                Vec3 planeOrigin, Vec3 normal);
BaffleExteriorOperator baffleExteriorOperator(const galerkin::Space &cavity,
                                              std::span<const unsigned> apertures,
                                              const galerkin::Space &exterior,
                                              const BemNumerics &numerics, double k,
                                              std::stop_token stop = {});
SurfaceMesh appendBaffleExterior(const SurfaceMesh &cavity, const SurfaceMesh &exterior);
} // namespace hornsculpt::bem_detail
