#pragma once

#include <hornsculpt/mesh.hpp>

namespace hornsculpt::galerkin {

// Requires a closed, orientable, non-self-intersecting surface mesh whose
// components do not touch. This supplements, rather than replaces, inspectMesh.
// Every component must point outward and belong to the same unbounded air
// domain. Separate bodies are supported; nested shells and sealed air cavities
// are rejected without changing winding or removing surfaces.
// Throws std::invalid_argument when these exterior-domain conditions fail.
void validateExteriorComponents(const SurfaceMesh &mesh);

} // namespace hornsculpt::galerkin
