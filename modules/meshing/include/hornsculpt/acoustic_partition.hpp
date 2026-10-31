#pragma once

#include <hornsculpt/acoustic_domains.hpp>
#include <hornsculpt/meshing.hpp>
#include <span>
#include <string>
#include <vector>

namespace hornsculpt {
// Physical surfaces bounding one finite air volume. Original winding points
// from material into air. Planar opening rims explicitly request transparent
// interfaces; they do not authorize repair of an invalid physical source.
struct AcousticRegionDefinition {
    std::uint32_t id{};
    std::string name;
    std::vector<std::uint32_t> faceIds;
    // Existing transparent mesh patches. Their winding points into the first
    // region listing them; a second listing names the other interior side.
    // A single listing connects to the exterior. Never also a physical face.
    std::vector<std::uint32_t> interfaceFaceIds;
    // A changed provider identity must be reviewed even if numeric IDs collide.
    bool needsAssignmentReview{};
    bool operator==(const AcousticRegionDefinition &) const = default;
};
struct AcousticPartition {
    SurfaceMesh mesh;
    BemSubdomainModel model;
};
// Validates ownership, shared interface topology, domain orientation and
// embedded closed boundaries. Independent of optional triangle-shape limits.
void validateAcousticPartition(const SurfaceMesh &mesh, const BemSubdomainModel &model);
// Leaves physical coordinates and triangles unchanged and meshes each unique
// planar opening once with the ordinary scalar Frontal-Delaunay contract.
// Each resulting oriented domain boundary is checked independently.
AcousticPartition partitionAcousticSurface(const SurfaceMesh &physical,
                                           std::span<const AcousticRegionDefinition> regions,
                                           const MeshingOptions &options,
                                           std::stop_token stop = {});
} // namespace hornsculpt
