#pragma once

#include <hornsculpt/acoustic_domains.hpp>
#include <hornsculpt/exterior_baffle.hpp>
#include <hornsculpt/meshing.hpp>
#include <hornsculpt/symmetry_planes.hpp>
#include <optional>
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
    // Explicit physical surfaces shared by this domain and another domain. A single
    // listing uses the exterior; two listings connect two interior domains.
    std::vector<std::uint32_t> twoSidedFaceIds;
    bool operator==(const AcousticRegionDefinition &) const = default;
};
// Presentation groups are disjoint subsets of the owning solver face lists.
struct AcousticElementGroup {
    std::uint32_t id{};
    std::string name;
    std::vector<std::uint32_t> faceIds;
    bool operator==(const AcousticElementGroup &) const = default;
};
struct AcousticDomainDefinition {
    std::uint32_t id{};
    std::string name;
    BemDomainKind kind{BemDomainKind::interior};
    std::vector<std::uint32_t> faceIds;
    std::vector<std::uint32_t> twoSidedFaceIds;
    bool needsAssignmentReview{};
    std::vector<AcousticElementGroup> groups;
    bool operator==(const AcousticDomainDefinition &) const = default;
};
struct AcousticInterfaceDefinition {
    std::uint32_t id{};
    std::string name;
    std::uint32_t firstDomainId{};
    std::uint32_t secondDomainId{};
    std::vector<std::uint32_t> faceIds;
    bool needsAssignmentReview{};
    std::vector<AcousticElementGroup> groups;
    bool operator==(const AcousticInterfaceDefinition &) const = default;
};
struct AcousticModelDefinition {
    std::vector<AcousticDomainDefinition> domains;
    std::vector<AcousticInterfaceDefinition> interfaces;
    std::optional<ExteriorBaffle> exteriorBaffle;
    bool operator==(const AcousticModelDefinition &) const = default;
};
struct AcousticPartition {
    SurfaceMesh mesh;
    BemSubdomainModel model;
};
// Validates ownership, shared interface topology, domain orientation and
// embedded closed boundaries. Independent of optional triangle-shape limits.
void validateAcousticPartition(const SurfaceMesh &mesh, const BemSubdomainModel &model,
                               const ExteriorBaffle *baffle = nullptr);
void validateSymmetryAcousticPartition(const SurfaceMesh &mesh, const BemSubdomainModel &model,
                                       const OrthogonalSymmetryPlanes &planes,
                                       const ExteriorBaffle *baffle = nullptr);
// Leaves physical coordinates and triangles unchanged and meshes each unique
// planar opening once with the ordinary scalar Frontal-Delaunay contract.
// Each resulting oriented domain boundary is checked independently.
AcousticPartition partitionAcousticSurface(const SurfaceMesh &physical,
                                           std::span<const AcousticRegionDefinition> regions,
                                           const MeshingOptions &options, std::stop_token stop = {},
                                           const OrthogonalSymmetryPlanes *symmetry = nullptr,
                                           const ExteriorBaffle *baffle = nullptr);
// Explicit domain and interface identities are retained in the result.
AcousticPartition partitionAcousticSurface(const SurfaceMesh &physical,
                                           const AcousticModelDefinition &definition,
                                           const MeshingOptions &options, std::stop_token stop = {},
                                           const OrthogonalSymmetryPlanes *symmetry = nullptr);
} // namespace hornsculpt
