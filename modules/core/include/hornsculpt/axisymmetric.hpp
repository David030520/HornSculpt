#pragma once

#include <hornsculpt/mesh.hpp>

namespace hornsculpt {

struct MeridianPoint {
    double radius{}; // metres, nonnegative
    double axial{};  // metres, from the explicit origin along axis
    bool operator==(const MeridianPoint &) const = default;
};

struct MeridianSegment {
    std::uint32_t first{}, second{}, faceId{};
    BoundaryCondition boundary{BoundaryCondition::wall};
    bool operator==(const MeridianSegment &) const = default;
};

// Oriented generator of a closed surface of revolution. The planar outward
// normal of first->second is (dz,-dr). Axis endpoints are not joined by a
// degenerate surface. Alternatively a closed generator describes a toroidal
// boundary. No implicit surface or physical face is introduced.
struct AxisymmetricMeridian {
    Vec3 origin{};
    Vec3 axis{0, 0, 1};
    Vec3 radial{1, 0, 0};
    std::vector<MeridianPoint> points;
    std::vector<MeridianSegment> segments;
    bool operator==(const AxisymmetricMeridian &) const = default;
};

} // namespace hornsculpt
