#pragma once
#include <hornsculpt/gmsh_mesh.hpp>

namespace hornsculpt {
struct MeshImportSelection {
    GmshSurfaceImport imported;
    double lengthScaleToMetres{};
};
} // namespace hornsculpt
