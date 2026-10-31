#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <hornsculpt/mesh.hpp>
#include <optional>
#include <span>
#include <stdexcept>

namespace hornsculpt {

struct SurfaceProjection {
    Vec3 position; // metres
    Vec3 normal;   // oriented outward from the source solid; finite unit vector
};

struct PlanarLineConstraint {
    // Other endpoints of the complete incoming/outgoing analytic line trims.
    std::array<Vec3, 2> boundaryEnds;
    // Unit normal with the oriented face domain on the loop's left.
    Vec3 loopNormal;
};

struct BoundaryCorner {
    std::uint32_t faceId{};
    Vec3 position; // metres; physical trim junction
    // Evidence only for an exact plane bounded by two analytic line trims.
    // A quality precheck must additionally prove retention of this junction
    // and both ray intervals in the actual mesh boundary. Curved/singular
    // cases have no constraint; their tangent angle alone does not prove that
    // a finite triangulation is impossible.
    std::optional<PlanarLineConstraint> constraint;
};

// The requested supporting-surface neighbourhood has no unique regular map.
// Optional chart extensions may retain source knots instead; used interiors
// must propagate this error. Corrupt/non-finite provider data are different errors.
class SurfaceProjectionDomainError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

// Immutable geometric source. Concurrent const calls must be safe. Projection
// failures and ambiguous shared curves are errors, never tessellation fallbacks.
class SurfaceGeometry {
  public:
    virtual ~SurfaceGeometry() = default;
    // Immutable physical trim evidence, not a completeness or mesh-quality
    // guarantee. Providers without exact trim information return an empty span.
    virtual std::span<const BoundaryCorner> boundaryCorners() const noexcept = 0;
    virtual SurfaceProjection projectSurface(std::uint32_t faceId, Vec3 approximate) const = 0;
    // Unsigned distance in metres to the physical trimmed face. Unlike a
    // position/normal map, distance remains well-defined at multiple minima.
    // Providers with such trims override the unique-projection default.
    virtual double distanceToSurface(std::uint32_t faceId, Vec3 position) const {
        const auto projected = projectSurface(faceId, position).position;
        return std::hypot(position.x - projected.x, position.y - projected.y,
                          position.z - projected.z);
    }
    // Local regular projection onto the supporting surface, without clamping
    // to face trims. Chordal chart interiors can cross a curved trim slightly.
    // A provider may restrict its parameter search domain, but must reject
    // missing/ambiguous stationary solutions rather than return its boundary.
    // Coordinates use metres; physical shared curves use projectBoundary.
    virtual SurfaceProjection projectSupportingSurface(std::uint32_t faceId,
                                                       Vec3 approximate) const = 0;
    // Three or more distinct physical faces identify a common source junction.
    // The provider validates topological identity and its modelling tolerance;
    // this is a local source-vertex query, not a remote closest-point search.
    // Input/output positions use metres; ambiguous junctions are errors.
    virtual Vec3 projectJunction(std::span<const std::uint32_t> faceIds,
                                 Vec3 approximate) const = 0;
    // At least two distinct physical faces identify their common trimmed curve.
    // Artificial chart seams inside one face use projectSupportingSurface.
    virtual Vec3 projectBoundary(std::span<const std::uint32_t> faceIds,
                                 Vec3 approximate) const = 0;
};

} // namespace hornsculpt
