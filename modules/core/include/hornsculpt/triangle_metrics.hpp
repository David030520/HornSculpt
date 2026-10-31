#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <hornsculpt/mesh.hpp>
#include <limits>
#include <numbers>

namespace hornsculpt {
struct TriangleMetrics {
    double twiceArea{};
    double minimumAngleDegrees{};
    double maximumAspectRatio{std::numeric_limits<double>::infinity()};
    double shortestEdge{};
    double medianEdge{};
    double longestEdge{};
};
// Shared geometric diagnostics; topology and scale-dependent degeneracy are checked separately.
inline TriangleMetrics triangleMetrics(Vec3 a, Vec3 b, Vec3 c) {
    const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
    const auto norm = [](Vec3 p) { return std::hypot(p.x, p.y, p.z); };
    const double area = norm({u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x});
    std::array<double, 3> edges{norm(u), norm(v), std::hypot(b.x - c.x, b.y - c.y, b.z - c.z)};
    std::sort(edges.begin(), edges.end());
    TriangleMetrics result{area,     0,        std::numeric_limits<double>::infinity(),
                           edges[0], edges[1], edges[2]};
    if (area <= 0 || edges[0] <= 0)
        return result;
    result.maximumAspectRatio = edges[2] * edges[2] / area;
    result.minimumAngleDegrees = 180;
    for (int i = 0; i < 3; ++i) {
        const auto opposite = edges[i], left = edges[(i + 1) % 3], right = edges[(i + 2) % 3];
        result.minimumAngleDegrees =
            std::min(result.minimumAngleDegrees,
                     std::acos(std::clamp((left * left + right * right - opposite * opposite) /
                                              (2 * left * right),
                                          -1.0, 1.0)) *
                         180 / std::numbers::pi);
    }
    return result;
}
} // namespace hornsculpt
