#pragma once

#include <cstddef>
#include <hornsculpt/mesh.hpp>
#include <hornsculpt/surface_geometry.hpp>
#include <hornsculpt/surface_parameters.hpp>
#include <map>
#include <memory>
#include <stop_token>

namespace hornsculpt {

struct MeshQuality {
    std::size_t boundaryEdges{};
    std::size_t nonManifoldEdges{};
    std::size_t inconsistentEdges{};
    std::size_t degenerateTriangles{};
    std::size_t connectedComponents{};
    std::size_t selfIntersectionPairs{};
    std::size_t nonManifoldVertices{};
    std::int64_t eulerCharacteristic{}; // Used vertices minus edges plus triangles.
    double minimumAngleDegrees{};
    double maximumAspectRatio{}; // longest edge / shortest altitude; equilateral = 2/sqrt(3).
    double minimumEdgeLength{};
    double medianEdgeLength{};
    double maximumEdgeLength{};
    double signedVolume{};
    bool watertight() const noexcept;
};

struct AdaptiveSizing {
    bool enabled{true};
    double minimumEdgeLength{0.0001};   // metres; lower bound for the local target.
    double geometricTolerance{0.00025}; // metres; curvature-based chord deviation target.
    double growthRate{0.3};             // maximum target increase per metre along source edges.
};

struct MeshQualityLimits {
    bool enabled{true}; // Optional angle/aspect constraints, independent of boundary validity.
    double minimumAngleDegrees{20.0};
    double maximumAspectRatio{4.0};
};

struct MeshingOptions {
    double targetEdgeLength{0.008}; // metres; isotropic scalar size field.
    MeshQualityLimits quality;
    AdaptiveSizing adaptive;
    // Optional refinement targets in metres for existing physical face IDs.
    // Values lie within [adaptive.minimumEdgeLength, targetEdgeLength].
    std::map<std::uint32_t, double> faceTargetEdgeLengths;
    // Optional exact geometry; null uses the source piecewise-linear surface.
    std::shared_ptr<const SurfaceGeometry> geometry;
    // Exact source-corner and trim bindings captured with this input mesh.
    std::shared_ptr<const SurfaceParameterization> parameterization;
};

struct MeshingResult {
    SurfaceMesh mesh;
    MeshQuality quality;
};

// Invalid connectivity/coordinates are rejected. No holes are filled or faces invented.
MeshQuality inspectMesh(const SurfaceMesh &mesh);

// Rebuilds an oriented closed boundary with a shared atlas and Frontal-Delaunay fronts.
// Face IDs, conditions and fixed coordinates survive. Cancellation returns no partial mesh.
MeshingResult remeshIsotropic(const SurfaceMesh &boundary, const MeshingOptions &options = {},
                              std::stop_token stop = {});

} // namespace hornsculpt
