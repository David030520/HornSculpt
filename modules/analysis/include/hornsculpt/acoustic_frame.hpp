#pragma once

#include <hornsculpt/mesh.hpp>

namespace hornsculpt {

// World-space origin in metres. Positive H rotates forward towards horizontal;
// positive V rotates forward towards forward x horizontal.
struct AcousticFrame {
    Vec3 origin{};
    Vec3 forward{0, 0, 1};
    Vec3 horizontal{1, 0, 0};
    bool operator==(const AcousticFrame &) const = default;
};

// Returns an orthonormal frame. Throws std::invalid_argument for non-finite
// values, zero directions or axes too nearly parallel to determine a stable roll.
AcousticFrame normalizedAcousticFrame(AcousticFrame frame);

// Geometry-derived frame: project world X into the transverse plane, using Y
// instead when forward is close to X. This preserves the original H/V convention.
AcousticFrame defaultAcousticFrame(Vec3 origin = {}, Vec3 forward = {0, 0, 1});

// Accepts the same inputs as normalizedAcousticFrame and returns the unit V axis.
Vec3 verticalDirection(const AcousticFrame &frame);

} // namespace hornsculpt
