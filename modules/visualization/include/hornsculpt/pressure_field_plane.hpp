#pragma once

#include <QColor>
#include <QString>
#include <array>
#include <complex>
#include <cstdint>
#include <hornsculpt/acoustic_frame.hpp>
#include <vector>

namespace hornsculpt {
// Immutable finite-distance samples in the same world/SI coordinates as the
// CAD. Row/column positions must match the normalized frame and stated ranges,
// including masked nodes. A renderer must omit triangles touching an invalid sample. No acoustic
// reconstruction or solver execution belongs to this presentation contract.
struct PressureFieldPlane {
    unsigned columns{}, rows{};
    std::vector<Vec3> positions;
    std::vector<std::complex<double>> pressure;
    std::vector<bool> valid;
    AcousticFrame frame;
    bool vertical{};
    double minimumLateral{}, maximumLateral{}, minimumAxial{}, maximumAxial{};
    double frequency{};
    std::uint64_t revision{};
    bool current{};
    QString excitation;
};
struct PressureFieldAppearance {
    int quantity{}; // 0 SPL, 1 phase, 2 rotated real RMS phasor
    double minimumDb{40}, maximumDb{110}, phaseDegrees{};
    bool operator==(const PressureFieldAppearance &) const = default;
};
void validatePressureFieldPlane(const PressureFieldPlane &plane);
void validatePressureFieldAppearance(const PressureFieldAppearance &appearance);
std::vector<std::uint32_t> pressureFieldTriangleIndices(const PressureFieldPlane &plane);
std::array<double, 2> pressureFieldRange(PressureFieldAppearance appearance);
// Shared legend/reference colour convention; the GL shader evaluates the same
// map AFTER interpolating complex pressure, so phase wraps never blend RGB.
QColor pressureFieldColor(double value, PressureFieldAppearance appearance);
} // namespace hornsculpt
