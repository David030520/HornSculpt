#pragma once

#include <complex>
#include <optional>
#include <span>
#include <vector>

namespace hornsculpt {

// RMS pressure in Pa; a true zero has level -infinity, without an artificial floor.
double soundPressureLevelDecibels(std::complex<double> pressure, double referencePressure = 20e-6);

// Volume impedance in Pa s/m^3 normalized by rho*c/S, using SI inputs.
std::complex<double> normalizedVolumeImpedance(std::complex<double> volumeImpedance,
                                               double sourceArea, double density,
                                               double soundSpeed);

struct AcousticsAngularSample {
    double angleDegrees{};
    double levelDecibels{};
};

struct Cta2034Cuts {
    std::optional<double> onAxisDb;
    std::optional<double> listeningWindowDb;
    std::optional<double> earlyReflectionsDb;
    std::optional<double> earlyReflectionsDiDb;
    std::optional<double> soundPowerDb;
    std::optional<double> soundPowerDiDb;
    std::optional<double> predictedInRoomDb;
};

// CTA-2034 directional power averages from exact H/V samples. Positive vertical
// angles point up. Missing required directions leave the corresponding curve
// unavailable. The CTA sound-power curve is a 70-direction approximation, not
// an exact spherical integral; the in-room curve is an empirical estimate.
Cta2034Cuts cta2034Cuts(std::span<const AcousticsAngularSample> horizontal,
                        std::span<const AcousticsAngularSample> vertical);

// Level at a freely selected angle, linearly interpolated in dB between samples.
// No extrapolation or interpolation across missing/non-finite levels.
std::optional<double> angularLevelDecibels(std::span<const AcousticsAngularSample> samples,
                                           double angleDegrees);

enum class AngularResponseNormalization { absolute, referenceAngle };

struct AngularFrequencyRow {
    double frequencyHz{};
    std::span<const AcousticsAngularSample> samples;
};

struct AngularFrequencySample {
    double frequencyHz{};
    std::optional<double> levelDecibels;
};

// Extract one angular cut through actual frequency rows, retaining their gaps.
// Angles use the same linear dB interpolation as angularLevelDecibels. Relative
// levels subtract the reference at EACH frequency. Missing/non-finite target or
// reference levels (including a true pressure zero) produce no level, never zero
// or a numeric floor. Frequencies must be finite, positive and strictly rising;
// no frequency interpolation, angular extrapolation or periodic wrapping occurs.
std::vector<AngularFrequencySample> angularFrequencyResponse(
    std::span<const AngularFrequencyRow> rows, double angleDegrees,
    AngularResponseNormalization normalization = AngularResponseNormalization::absolute,
    double referenceAngleDegrees = 0.0);

struct Beamwidth {
    double leftDegrees{};
    double rightDegrees{};
    double widthDegrees{};
};

// First threshold crossings of the connected lobe containing 0 degrees, relative
// to its measured/interpolated on-axis level. Interpolation is linear in dB.
// Missing axis/crossings or non-finite levels on the path return no result;
// non-finite/unordered angles and non-positive attenuation are invalid arguments.
std::optional<Beamwidth> beamwidthDegrees(std::span<const AcousticsAngularSample> samples,
                                          double attenuationDb = 6.0);

} // namespace hornsculpt
