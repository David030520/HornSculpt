#pragma once
#include <hornsculpt/mesh.hpp>

namespace hornsculpt {
// Unit normals point into the retained half-spaces. Both planes pass through
// origin. Coordinates and the explicit distance tolerance are in metres.
struct OrthogonalSymmetryPlanes {
    Vec3 origin;
    Vec3 firstNormal;
    Vec3 secondNormal;
    double planeDistanceTolerance{};
    bool operator==(const OrthogonalSymmetryPlanes &) const = default;
};
} // namespace hornsculpt
