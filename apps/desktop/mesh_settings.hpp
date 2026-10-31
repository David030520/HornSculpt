#pragma once

#include <cmath>
#include <hornsculpt/meshing.hpp>
#include <stdexcept>

namespace hornsculpt {
struct MeshSettings {
    double targetEdgeLength{.008};
    AdaptiveSizing adaptive;
    MeshQualityLimits quality;
};

inline void validateMeshSettings(const MeshSettings &s) {
    const auto length = [](double value) {
        return std::isfinite(value) && value >= .000001 && value <= .1;
    };
    if (!length(s.targetEdgeLength) || s.targetEdgeLength < .00005 ||
        !length(s.adaptive.minimumEdgeLength) || !length(s.adaptive.geometricTolerance) ||
        s.adaptive.minimumEdgeLength > s.targetEdgeLength ||
        !std::isfinite(s.adaptive.growthRate) || s.adaptive.growthRate < .01 ||
        s.adaptive.growthRate > 1 || !std::isfinite(s.quality.minimumAngleDegrees) ||
        s.quality.minimumAngleDegrees < .001 || s.quality.minimumAngleDegrees > 59.999 ||
        !std::isfinite(s.quality.maximumAspectRatio) || s.quality.maximumAspectRatio < 1.155 ||
        s.quality.maximumAspectRatio > 1000000)
        throw std::invalid_argument("Ungültige gemeinsame Netzeinstellungen.");
}
} // namespace hornsculpt
