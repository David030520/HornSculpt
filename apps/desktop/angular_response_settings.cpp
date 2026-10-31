#include "angular_response_settings.hpp"

#include <cmath>
#include <set>
#include <stdexcept>

namespace hornsculpt {

void validateAngularResponseSettings(const AngularResponseSettings &settings) {
    if (settings.anglesDegrees.empty() || settings.anglesDegrees.size() > 12)
        throw std::invalid_argument("Für Winkelfrequenzgänge sind 1 bis 12 Winkel erforderlich.");
    const auto validAngle = [](double angle) {
        return std::isfinite(angle) && angle >= -180 && angle <= 180;
    };
    std::set<double> seen;
    for (const double angle : settings.anglesDegrees) {
        if (!validAngle(angle))
            throw std::invalid_argument(
                "Winkel muss endlich sein und zwischen −180° und 180° liegen.");
        if (!seen.insert(angle).second)
            throw std::invalid_argument("Winkel für Frequenzgänge müssen eindeutig sein.");
    }
    if (!validAngle(settings.referenceAngleDegrees))
        throw std::invalid_argument(
            "Referenzwinkel muss endlich sein und zwischen −180° und 180° liegen.");
}

} // namespace hornsculpt
