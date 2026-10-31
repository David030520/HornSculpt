#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>
#include <set>
#include <span>
#include <stdexcept>

namespace hornsculpt::bem_detail {
// An explicit planning raster retains the existing complete-sweep policy:
// plan at its maximum k, even when this call solves only a scheduled subset.
// It is not a claim of an error certificate uniform over a frequency interval.
inline double planningWavenumber(std::span<const double> solved, std::span<const double> planned,
                                 double soundSpeed) {
    if (solved.empty() || !std::isfinite(soundSpeed) || soundSpeed <= 0)
        throw std::invalid_argument("Interaction planning requires frequencies and sound speed.");
    std::set<double> raster;
    double maximum{};
    for (const auto frequency : planned) {
        if (!std::isfinite(frequency) || frequency <= 0 || !raster.insert(frequency).second)
            throw std::invalid_argument(
                "Interaction planning frequencies must be finite, positive and unique.");
        maximum = std::max(maximum, frequency);
    }
    for (const auto frequency : solved) {
        if (!std::isfinite(frequency) || frequency <= 0 ||
            (!planned.empty() && !raster.contains(frequency)))
            throw std::invalid_argument(
                "Solved frequency is invalid or absent from its interaction planning raster.");
        if (planned.empty())
            maximum = std::max(maximum, frequency);
    }
    const double result = 2 * std::numbers::pi * maximum / soundSpeed;
    if (!std::isfinite(result))
        throw std::invalid_argument("Interaction planning wavenumber overflows.");
    return result;
}
} // namespace hornsculpt::bem_detail
