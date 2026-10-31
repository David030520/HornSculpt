#pragma once
#include <cstddef>
#include <hornsculpt/exterior_baffle.hpp>
#include <hornsculpt/mesh.hpp>
#include <hornsculpt/symmetry_planes.hpp>

namespace hornsculpt {
struct SymmetrySectorValidation {
    std::size_t firstPlaneEdges{};
    std::size_t secondPlaneEdges{};
    std::size_t intersectionAxisEdges{};
};
// Checks one oriented domain-boundary view, not a combined interface complex.
// Each open edge must have both endpoints on at least one plane within the
// given tolerance; all vertices must lie in the retained half-spaces. Rejects
// degenerate, intersecting or non-manifold meshes and coplanar triangles.
// Original coordinates are never changed. A valid result says nothing about
// excitation symmetry or solver support.
SymmetrySectorValidation validateSymmetrySector(const SurfaceMesh &mesh,
                                                const OrthogonalSymmetryPlanes &planes,
                                                const ExteriorBaffle *baffle = nullptr);
// Open edges must lie on the baffle; geometry remains in its exterior half-space.
void validateBaffleBoundary(const SurfaceMesh &mesh, const ExteriorBaffle &baffle);
} // namespace hornsculpt
