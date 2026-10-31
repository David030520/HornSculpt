#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <hornsculpt/meshing.hpp>
#include <memory>
#include <span>
#include <stop_token>
#include <vector>

namespace hornsculpt::detail {
class SurfaceSizingField;
struct ChartLookup;
struct AtlasCurveRegistry;

struct ChartPoint {
    double u{}, v{};
};

struct SurfaceChartVertex {
    ChartPoint uv;
    // Shared physical boundaries and both sides of virtual cuts use the same ID.
    std::uint32_t atlasVertex{};
};

struct SurfaceChartLoop {
    // No repeated final vertex. Domain lies to the left of each directed edge.
    std::vector<SurfaceChartVertex> vertices;
};

struct ChartReferenceTriangle {
    std::array<ChartPoint, 3> uv;
    std::array<Vec3, 3> positions;
    std::array<double, 3> sizes;
};

struct SurfaceChartSample {
    Vec3 position, du, dv;
    double targetEdgeLength{};
};
struct SurfaceChartMetric {
    Vec3 du, dv;
    double targetEdgeLength{};
};

struct SurfaceChart {
    std::uint32_t faceId{};
    BoundaryCondition boundary{BoundaryCondition::wall};
    std::vector<SurfaceChartLoop> loops;
    std::vector<SurfaceChartVertex> fixedPoints;
    // Geometry lookup only; these vertices are not mandatory output nodes.
    std::vector<ChartReferenceTriangle> reference;
    // Builder-owned immutable acceleration; manually constructed test charts may
    // omit it. The planar extension also evaluates coarser polygonal hole rims.
    std::shared_ptr<const ChartLookup> lookup;
    std::shared_ptr<const SurfaceSizingField> sizing;
    std::shared_ptr<const SurfaceGeometry> geometry;
    std::shared_ptr<const SurfaceParameterization> parameterization;
    // Native V is reflected only in chart coordinates for an outward UV frame.
    int parameterOrientation{1};
    double parameterAccuracy{}; // metres; requested representation accuracy
    bool planar{};
    Vec3 planeOrigin, planeU, planeV;

    bool contains(ChartPoint uv) const;
    // Optional front candidates must respect the provider's source domain even
    // when the polygon containment test accepts a nearby boundary point.
    bool containsParameterDomain(ChartPoint uv) const;
    SurfaceChartSample sample(ChartPoint uv) const;
    // Approximate tangent metric at a constrained chord. It supplies no
    // position map: a physical chord can leave a finite supporting domain.
    SurfaceChartMetric boundaryMetric(ChartPoint uv) const;
    // Scalar sizing belongs to the source field and does not project geometry.
    double targetSize(ChartPoint uv) const;

  private:
    SurfaceChartSample referenceSample(ChartPoint uv) const;
};

struct SurfaceAtlas {
    std::vector<Vec3> vertices;
    std::vector<SurfaceChart> charts;
    // Atlas indices. Fully fixed physical components bypass triangulation.
    std::vector<std::uint32_t> fixedVertices;
    std::vector<Triangle> fixedTriangles;
    std::shared_ptr<const AtlasCurveRegistry> curves;
};

// Geometry references and scalar fields remain immutable; shared boundary
// sampling may be refined between front passes. Artificial chart cuts create
// neither physical faces nor local sizing seeds.
SurfaceAtlas makeSurfaceAtlas(const SurfaceMesh &source, const MeshingOptions &options,
                              std::stop_token stop = {});

// Checks proven planar line corners only after their current loop and complete
// registry intervals have been retained. Curved/unknown evidence cannot block.
void validateAtlasBoundaryCorners(const SurfaceAtlas &atlas, const MeshingOptions &options,
                                  std::stop_token stop = {});

// Atomic, conforming splits at original curve-arclength midpoints. Retained
// fixed-triangle edges are explicitly infeasible and are never split on one side.
// Returns the sorted unique charts whose loops changed. Existing vertex IDs,
// positions, reference geometry and sizing remain unchanged.
std::vector<std::size_t> refineSurfaceAtlas(SurfaceAtlas &atlas,
                                            std::span<const std::array<std::uint32_t, 2>> segments,
                                            std::stop_token stop = {});

} // namespace hornsculpt::detail
