#include <hornsculpt/acoustic_frame.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hornsculpt {
namespace {
bool finite(Vec3 value) {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

Vec3 normalized(Vec3 value) {
    if (!finite(value))
        throw std::invalid_argument("Akustische Richtung muss endlich sein.");
    const auto scale = std::max({std::abs(value.x), std::abs(value.y), std::abs(value.z)});
    if (scale == 0)
        throw std::invalid_argument("Akustische Richtung darf kein Nullvektor sein.");
    // Scale first so very large or small finite direction vectors remain valid.
    value = {value.x / scale, value.y / scale, value.z / scale};
    const auto length = std::hypot(value.x, value.y, value.z);
    return {value.x / length, value.y / length, value.z / length};
}

Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
} // namespace

AcousticFrame normalizedAcousticFrame(AcousticFrame frame) {
    if (!finite(frame.origin))
        throw std::invalid_argument("Akustischer Ursprung muss endlich sein.");
    frame.forward = normalized(frame.forward);
    const auto vertical = cross(frame.forward, normalized(frame.horizontal));
    // Cross-product magnitude is sin(angle). Below sqrt(epsilon), determining
    // roll amplifies direction rounding by more than half the double precision.
    if (std::hypot(vertical.x, vertical.y, vertical.z) <=
        std::sqrt(std::numeric_limits<double>::epsilon()))
        throw std::invalid_argument("Akustische Achsen dürfen nicht parallel sein.");
    frame.horizontal = normalized(cross(normalized(vertical), frame.forward));
    return frame;
}

AcousticFrame defaultAcousticFrame(Vec3 origin, Vec3 forward) {
    forward = normalized(forward);
    return normalizedAcousticFrame(
        {origin, forward, std::abs(forward.x) < .9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0}});
}

Vec3 verticalDirection(const AcousticFrame &frame) {
    const auto normalizedFrame = normalizedAcousticFrame(frame);
    return normalized(cross(normalizedFrame.forward, normalizedFrame.horizontal));
}
} // namespace hornsculpt
