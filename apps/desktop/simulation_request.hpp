#pragma once

#include "project_file.hpp"
#include <span>

namespace hornsculpt {
// Result-storage estimate, not a complete process peak. Boundary
// arrays are shared; per-snapshot observer and diagnostic arrays are independent.
struct ResultMemoryDimensions {
    std::size_t triangles{}, vertices{}, exteriorTriangles{}, observers{}, faces{};
    std::size_t frequencies{}, excitations{}, snapshots{};
    std::size_t domains{};
};
long double estimateResultBytes(const ResultMemoryDimensions &dimensions);

// Bound virtual ports and unresolved previous assignments must never silently
// become material walls when coupling is disabled.
bool hasBoundInterfacePorts(const lem::Network &network);
void validateAxisymmetricSamplingSettings(const AxisymmetricSamplingSettings &settings);
OrthogonalSymmetryPlanes inferXYSymmetry(const SurfaceMesh &mesh, double toleranceMetres);
AcousticFrame simulationFrame(const SimulationSettings &settings, Vec3 sourceOrigin = {},
                              Vec3 sourceAxis = {0, 0, 1});

// Structural checks only; the solver additionally checks the cavity and aperture mesh.
void validateInfiniteBaffleSettings(const InfiniteBaffle &baffle);

std::vector<std::uint8_t> exteriorObserverValidity(const BemRequest &request);
inline bool observerValid(const BemRequest &request, std::size_t index) {
    return request.observerValidity.empty() ||
           (index < request.observerValidity.size() && request.observerValidity[index] != 0);
}
// The desktop and parameter studies share the same 1 m, 5 degree H/V sampling.
BemRequest makeSimulationRequest(const SimulationSettings &settings,
                                 std::span<const FaceVelocity> velocities, Vec3 origin = {},
                                 Vec3 axis = {0, 0, 1});
} // namespace hornsculpt
