#pragma once

#include <hornsculpt/mesh.hpp>
#include <set>
#include <stdexcept>
#include <utility>

namespace hornsculpt {

// Apply a source-relative physical-face orientation mask. Coordinates, tags,
// velocities and fixed vertices remain unchanged. This is not a validity repair:
// callers must validate the resulting closed exterior before meshing/BEM use.
inline SurfaceMesh applyFaceOrientation(const SurfaceMesh &source,
                                        const std::set<std::uint32_t> &reversedFaceIds) {
    auto missing = reversedFaceIds;
    for (const auto &triangle : source.triangles)
        missing.erase(triangle.faceId);
    if (!missing.empty())
        throw std::invalid_argument("Orientierungsmaske enthaelt unbekannte Flaechen.");
    auto result = source;
    for (auto &triangle : result.triangles)
        if (reversedFaceIds.contains(triangle.faceId))
            std::swap(triangle.vertices[1], triangle.vertices[2]);
    return result;
}

} // namespace hornsculpt
