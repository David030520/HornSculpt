#pragma once
#include "project_file.hpp"
#include <stop_token>
namespace hornsculpt {
// Shared meshing of the local Draft followed by rigid copies. CAD retains its
// exact provider and independent scalar mesh size; the final union is checked.
SurfaceMesh assemblyPreview(const SurfaceMesh &local, const Project &project);
AcousticFrame assemblyAcousticFrame(const DraftAssembly &assembly, double hornLength);
std::shared_ptr<const CadSnapshot> appendAssemblyCad(const CadSnapshot &current,
                                                     const CadSnapshot &added);
MeshingResult meshDraftAssembly(const Project &project, const SurfaceMesh &local,
                                std::stop_token stop = {});
} // namespace hornsculpt
