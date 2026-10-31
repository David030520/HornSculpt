#pragma once
#include <hornsculpt/mesh.hpp>
namespace hornsculpt {
// The normal points into the exterior acoustic half-space; units are metres.
struct ExteriorBaffle {
    Vec3 origin;
    Vec3 normal;
    double planeDistanceTolerance{};
    bool operator==(const ExteriorBaffle &) const = default;
};
} // namespace hornsculpt
