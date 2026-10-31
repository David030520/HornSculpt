#include <hornsculpt/geometry_analysis.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace hornsculpt {
namespace {
// Periodic trapezoidal integration converges rapidly for smooth polar sections.
constexpr std::size_t angularSamples = 1024;
double sectionArea(const DraftGeometry &geometry, double z, std::stop_token stop) {
    double sum = 0.0;
    for (std::size_t i = 0; i < angularSamples; ++i) {
        if ((i % 64) == 0 && stop.stop_requested())
            throw std::runtime_error("Geometrieanalyse abgebrochen.");
        const double phi = 2.0 * std::numbers::pi * static_cast<double>(i) / angularSamples;
        const double r = geometry.radius(z, phi);
        sum += r * r;
    }
    return std::numbers::pi * sum / angularSamples;
}

double sectionAreaDerivative(const DraftGeometry &geometry, double z, AxialDerivativeSide side,
                             std::stop_token stop) {
    double sum = 0.0;
    for (std::size_t i = 0; i < angularSamples; ++i) {
        if ((i % 64) == 0 && stop.stop_requested())
            throw std::runtime_error("Geometrieanalyse abgebrochen.");
        const double phi = 2.0 * std::numbers::pi * static_cast<double>(i) / angularSamples;
        sum += geometry.radius(z, phi) * geometry.radiusDerivative(z, phi, side);
    }
    return 2.0 * std::numbers::pi * sum / angularSamples;
}
} // namespace

GeometryAnalysis analyzeGeometry(const DraftGeometry &geometry, std::size_t sampleCount, std::stop_token stop) {
    if (sampleCount < 2 || sampleCount > 100000)
        throw std::invalid_argument("Analyse benoetigt 2..100000 Stuetzstellen.");
    GeometryAnalysis result;
    const double dz = geometry.length() / static_cast<double>(sampleCount - 1);
    std::vector<double> positions;
    for (std::size_t i = 0; i < sampleCount; ++i) {
        positions.push_back(i + 1 == sampleCount ? geometry.length() : dz * static_cast<double>(i));
    }
    if (geometry.parameters().slot.enabled) {
        const auto breaks = geometry.axialBreakpoints();
        // Short physical inlet sections must remain resolved even when the
        // user's global analysis stations are much farther apart.
        for (std::size_t i = 1; i < breaks.size(); ++i)
            for (unsigned j = 0; j <= 16; ++j) {
                const double z = std::lerp(breaks[i - 1], breaks[i], j / 16.0);
                const auto nearby = std::find_if(positions.begin(), positions.end(), [&](double x) {
                    return std::abs(x - z) <=
                           16 * std::numeric_limits<double>::epsilon() * geometry.length();
                });
                if (nearby != positions.end())
                    *nearby = z;
                else
                    positions.push_back(z);
            }
        std::sort(positions.begin(), positions.end());
    }
    result.samples.reserve(positions.size() + 2);
    double volume = 0.0, previousArea = 0.0, previousZ = 0.0;
    for (const double z : positions) {
        const double area = sectionArea(geometry, z, stop);
        if (z != 0) {
            const double midpointArea = sectionArea(geometry, (z + previousZ) / 2.0, stop);
            volume += (z - previousZ) / 6.0 * (previousArea + 4.0 * midpointArea + area);
        }
        const auto sample = [&](AxialDerivativeSide side) {
            return GeometrySample{
                z,
                volume,
                area,
                sectionAreaDerivative(geometry, z, side, stop),
                std::atan(geometry.radiusDerivative(z, 0.0, side)),
                std::atan(geometry.radiusDerivative(z, std::numbers::pi / 2.0, side)),
                std::atan(geometry.radiusDerivative(z, std::numbers::pi, side)),
                std::atan(geometry.radiusDerivative(z, -std::numbers::pi / 2.0, side))};
        };
        const auto at = sample(AxialDerivativeSide::unspecified);
        const auto &slot = geometry.parameters().slot;
        const bool sharp = slot.enabled &&
                           z == slot.adapterLength + slot.straightLength &&
                           std::isnan(at.areaDerivative);
        if (sharp)
            result.samples.push_back(sample(AxialDerivativeSide::before));
        // NaN breaks derivative plots at a real crease. Volume and area stay
        // continuous; the two adjacent samples retain both one-sided limits.
        result.samples.push_back(at);
        if (sharp)
            result.samples.push_back(sample(AxialDerivativeSide::after));
        previousArea = area;
        previousZ = z;
    }
    return result;
}

} // namespace hornsculpt
