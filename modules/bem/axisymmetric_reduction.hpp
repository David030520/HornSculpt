#pragma once

#include <hornsculpt/bem.hpp>
#include <span>

namespace hornsculpt::bem_detail {
struct AxisymmetricLayout {
    std::size_t pressureCount{}, fluxCount{};
};

// Checks the public reduction against the actual full mesh, including angular
// closure, orientation and complete coverage by rotated target-sector panels.
AxisymmetricLayout validateAxisymmetricReduction(const SurfaceMesh &localMesh, Vec3 worldOrigin,
                                                 const AxisymmetricBem &reduction);
// Expanded traces are aliases of shared meridional DOFs, not independently
// sampled azimuthal measurements. Their aliases must retain the same value.
void validateAxisymmetricTrace(const AxisymmetricBem &reduction,
                               std::span<const std::complex<double>> pressure,
                               std::span<const std::complex<double>> normalDerivative);
} // namespace hornsculpt::bem_detail
