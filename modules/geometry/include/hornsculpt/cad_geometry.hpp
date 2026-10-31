#pragma once

#include <hornsculpt/mesh.hpp>
#include <hornsculpt/surface_geometry.hpp>
#include <hornsculpt/surface_parameters.hpp>
#include <memory>
#include <span>
#include <stop_token>
#include <string>
#include <vector>

namespace hornsculpt {
struct CadBrepPart {
    std::string data;
    // Export-local face traversal order within this exact BRep artifact.
    std::vector<std::uint32_t> faceIds;
};
struct CadFaceMeasures {
    std::uint32_t faceId{};
    double area{};
    Vec3 center;
    Vec3 normal;
    std::string surfaceType;
    std::uint32_t edgeCount{};
    std::uint32_t wireCount{};
};

// Immutable exact CAD source. Public coordinates and measures use SI units.
class CadGeometry final : public SurfaceGeometry {
  public:
    static std::shared_ptr<const CadGeometry> fromBrepParts(std::span<const CadBrepPart> parts,
                                                            double metresPerUnit);
    static std::shared_ptr<const CadGeometry> fromStep(const std::string &path);
    ~CadGeometry() override;
    std::span<const CadFaceMeasures> faces() const noexcept;
    // Canonical BRep artifacts in metres, for self-contained project persistence.
    std::span<const CadBrepPart> brepParts() const noexcept;
    std::span<const BoundaryCorner> boundaryCorners() const noexcept override;
    MeshingSource meshingSource(double geometricTolerance, std::stop_token stop = {}) const;
    SurfaceMesh previewMesh(double geometricTolerance, std::stop_token stop = {}) const;
    SurfaceProjection projectSurface(std::uint32_t faceId, Vec3 approximate) const override;
    double distanceToSurface(std::uint32_t faceId, Vec3 position) const override;
    SurfaceProjection projectSupportingSurface(std::uint32_t faceId,
                                               Vec3 approximate) const override;
    Vec3 projectJunction(std::span<const std::uint32_t> faceIds, Vec3 approximate) const override;
    Vec3 projectBoundary(std::span<const std::uint32_t> faceIds, Vec3 approximate) const override;

  private:
    struct Data;
    explicit CadGeometry(std::unique_ptr<Data> data);
    std::unique_ptr<Data> data_;
};
} // namespace hornsculpt
