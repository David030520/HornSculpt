#pragma once

#include <cstddef>
#include <hornsculpt/mesh.hpp>
#include <utility>
#include <vector>

namespace hornsculpt::detail {
struct EmbeddingValidation {
    std::size_t selfIntersectionPairs{};
    std::size_t nonManifoldVertices{};
};

// Requires finite vertices and valid triangle indices, checked by inspectMesh.
EmbeddingValidation
validateEmbedding(const SurfaceMesh &mesh,
                  std::vector<std::pair<std::size_t, std::size_t>> *intersections = nullptr);
} // namespace hornsculpt::detail
