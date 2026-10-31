#pragma once

#include <vector>

namespace hornsculpt {

// Display-only selection; changing it does not change the acoustic request.
struct AngularResponseSettings {
    std::vector<double> anglesDegrees{0, 15, 30, 45, 60};
    bool normalizeToAngle{};
    double referenceAngleDegrees{};
    bool operator==(const AngularResponseSettings &) const = default;
};

// Angles are signed degrees in [-180, 180]. Throws invalid_argument unless
// the selection contains 1..12 distinct finite angles and a finite reference.
void validateAngularResponseSettings(const AngularResponseSettings &settings);

} // namespace hornsculpt
