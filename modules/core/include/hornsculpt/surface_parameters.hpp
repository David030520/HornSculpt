#pragma once

#include <array>
#include <cstdint>
#include <hornsculpt/mesh.hpp>
#include <memory>
#include <optional>
#include <span>

namespace hornsculpt {

// Provider coordinates have their own units; positions and derivatives are
// metres and metres per provider parameter respectively.
struct ParameterPoint {
    double u{}, v{};
    bool operator==(const ParameterPoint &) const = default;
};

struct ParameterDomain {
    // Admissible bounds of this immutable source snapshot; null is unbounded.
    // Any source-defined modelling halo beyond nominal surface patch bounds
    // must be checked in world space without clamping positions or derivatives.
    std::optional<std::array<double, 2>> uBounds, vBounds;
    std::optional<double> uPeriod, vPeriod;
    // Outward source orientation relative to du x dv. Derivatives stay native.
    int orientation{1};
};

struct ParametricSurfaceSample {
    Vec3 position, du, dv;
};

struct ParametricCurveSample {
    Vec3 position, tangent;
};

struct CurveUseSample {
    ParameterPoint uv, uvDerivative;
    // Unmodified source modelling tolerance, not a mesh acceptance allowance.
    // Consumers measure |surface(uv) - curve(t)| against the requested accuracy.
    double modelToleranceMetres{};
    // Present only at an exact native endpoint; kept separate from edge/face
    // tolerance and never added to the user's geometric accuracy budget.
    std::optional<double> vertexToleranceMetres;
};

struct CurveUseInterval {
    std::uint32_t curveId{}, useId{};
    // Canonical spatial-curve parameters, directed from corner i to i+1.
    // Both uses of a periodic seam share curveId but retain distinct useId/UVs.
    double t0{}, t1{};
};

struct ParameterTriangle {
    std::array<ParameterPoint, 3> uv;
    std::array<std::optional<CurveUseInterval>, 3> edges;
};

// Immutable parameter records captured before source vertices are welded.
// Concurrent const calls must be safe. Native topology and corner coordinates
// belong to this one source snapshot; no nearest-point reconstruction or cache
// of a previous preview call is part of this contract.
class SurfaceParameterization {
  public:
    virtual ~SurfaceParameterization() = default;
    // Positions, indices, face IDs and fixed nodes identify the snapshot.
    // Boundary conditions may be assigned by the caller without changing it.
    virtual const SurfaceMesh &sourceMesh() const noexcept = 0;
    virtual std::span<const ParameterTriangle> triangles() const noexcept = 0;
    virtual ParameterDomain domain(std::uint32_t faceId) const = 0;
    virtual ParametricSurfaceSample surface(std::uint32_t faceId, ParameterPoint uv) const = 0;
    // Optional conservative evidence for a distance upper bound: strictly
    // inside the physical trim, beyond membership in the natural UV domain.
    // This potentially expensive query is separate from position/D1 evaluation.
    // False includes unknown; callers then perform the bounded-face projection.
    virtual bool isInsideTrim(std::uint32_t, ParameterPoint) const {
        return false;
    }
    virtual ParametricCurveSample curve(std::uint32_t curveId, double t) const = 0;
    virtual CurveUseSample curveUse(std::uint32_t useId, double t) const = 0;
    // Strictly increasing, finite partition of the complete native trim,
    // including both endpoints and every polynomial/continuity span boundary.
    // Quadrature must not infer these spans from geometric tessellation or D1
    // samples: narrow spans can be smooth but have strongly localized speed.
    virtual std::span<const double> curvePartition(std::uint32_t curveId) const = 0;
};

struct MeshingSource {
    SurfaceMesh mesh;
    std::shared_ptr<const SurfaceParameterization> parameters;
};

} // namespace hornsculpt
