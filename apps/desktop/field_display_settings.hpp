#pragma once

#include <cmath>
#include <initializer_list>
#include <stdexcept>

namespace hornsculpt {
struct FieldDisplaySettings {
    bool vertical{};
    double minimumLateral{-0.1}, maximumLateral{0.1}; // m, relative to result observation frame
    double minimumAxial{-0.05}, maximumAxial{0.2};    // m
    unsigned lateralCount{41}, axialCount{51};
    int quantity{}; // 0 SPL, 1 phase, 2 real RMS phasor rotated by phaseDegrees
    double minimumDb{40}, maximumDb{110}, phaseDegrees{};
    double exclusionDistance{0.0001}; // m
    double frequencyHz{}; // 0 selects the available middle frequency
    int excitationIndex{-1}; // -1 selects the first available excitation
    bool operator==(const FieldDisplaySettings &) const = default;
};
struct PhaseDisplaySettings {
    double angleDegrees{}, referenceAngleDegrees{}, delaySeconds{};
    bool relative{}, unwrap{};
    bool operator==(const PhaseDisplaySettings &) const = default;
};
inline void validateFieldDisplaySettings(const FieldDisplaySettings &s) {
    for (const double v : {s.minimumLateral, s.maximumLateral, s.minimumAxial, s.maximumAxial,
                           s.minimumDb, s.maximumDb, s.phaseDegrees, s.exclusionDistance, s.frequencyHz})
        if (!std::isfinite(v))
            throw std::invalid_argument("Ungültige Schallfeld-Einstellung.");
    if (s.maximumLateral <= s.minimumLateral || s.maximumAxial <= s.minimumAxial ||
        s.minimumLateral < -100 || s.maximumLateral > 100 || s.minimumAxial < -100 ||
        s.maximumAxial > 100 || s.lateralCount < 3 || s.lateralCount > 201 || s.axialCount < 3 ||
        s.axialCount > 201 || s.quantity < 0 || s.quantity > 2 || s.maximumDb <= s.minimumDb ||
        s.minimumDb < -300 || s.maximumDb > 300 || s.phaseDegrees < -180 || s.phaseDegrees > 180 ||
        s.exclusionDistance < 1e-6 || s.exclusionDistance > .01 || s.frequencyHz < 0 ||
        s.frequencyHz > 100000 || s.excitationIndex < -1 || s.excitationIndex > 10000)
        throw std::invalid_argument("Ungültige Schallfeld-Einstellung.");
}
inline void validatePhaseDisplaySettings(const PhaseDisplaySettings &s) {
    if (!std::isfinite(s.angleDegrees) || !std::isfinite(s.referenceAngleDegrees) ||
        !std::isfinite(s.delaySeconds) || std::abs(s.angleDegrees) > 180 ||
        std::abs(s.referenceAngleDegrees) > 180 || std::abs(s.delaySeconds) > 1)
        throw std::invalid_argument("Ungültige Phasen-Einstellung.");
}
} // namespace hornsculpt
