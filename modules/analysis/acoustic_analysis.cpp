#include <hornsculpt/acoustic_analysis.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hornsculpt {
namespace {
double fraction(double value, double a, double b) {
    const double scale = std::max({std::abs(value), std::abs(a), std::abs(b)});
    return (value / scale - a / scale) / (b / scale - a / scale);
}
} // namespace

Cta2034Cuts cta2034Cuts(std::span<const AcousticsAngularSample> horizontal,
                        std::span<const AcousticsAngularSample> vertical) {
    const auto find = [](std::span<const AcousticsAngularSample> samples, int angle)
        -> std::optional<double> {
        for (const auto &sample : samples)
            if (sample.angleDegrees == angle && !std::isnan(sample.levelDecibels) &&
                sample.levelDecibels != std::numeric_limits<double>::infinity())
                return sample.levelDecibels;
        return {};
    };
    struct Direction { bool vertical; int angle; double weight; };
    const auto average = [&](const std::vector<Direction> &directions)
        -> std::optional<double> {
        std::vector<std::pair<double, double>> levels;
        double totalWeight = 0;
        for (const auto &[isVertical, angle, weight] : directions) {
            const auto level = find(isVertical ? vertical : horizontal, angle);
            if (!level)
                return {};
            levels.emplace_back(*level, weight);
            totalWeight += weight;
        }
        const double maximum = std::max_element(levels.begin(), levels.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; })->first;
        if (!std::isfinite(maximum))
            return maximum;
        double sum = 0;
        for (const auto &[level, weight] : levels)
            sum += weight * std::pow(10., (level - maximum) / 10.);
        return maximum + 10 * std::log10(sum / totalWeight);
    };
    // ANSI/CEA-2034-A sections 5.2, 13 and appendices C and F. The five reflection
    // groups contribute equally; directions within each group are power averaged.
    const std::vector<Direction> listening{
        {false, 0, 1}, {false, -10, 1}, {false, 10, 1}, {false, -20, 1}, {false, 20, 1},
        {false, -30, 1}, {false, 30, 1}, {true, -10, 1}, {true, 10, 1}};
    std::vector<Direction> reflections{{true, -20, 1./3}, {true, -30, 1./3}, {true, -40, 1./3},
        {true, 40, 1./3}, {true, 50, 1./3}, {true, 60, 1./3},
        {false, 0, 1./7}, {false, 180, 1./3}};
    for (int angle = 10; angle <= 90; angle += 10) {
        const double weight = angle <= 30 ? 1./7 : angle <= 80 ? 1./10 : 1./3;
        reflections.push_back({false, -angle, weight});
        reflections.push_back({false, angle, weight});
    }
    constexpr double weights[] = {.000604486, .004730189, .008955027, .012387354,
        .014989611, .016868154, .018165962, .019006744, .019477787, .019629373,
        .019477787, .019006744, .018165962, .016868154, .014989611, .012387354,
        .008955027, .004730189, .000604486};
    std::vector<Direction> power{{false, 0, weights[0]}, {false, 180, weights[18]}};
    for (int angle = 10; angle < 180; angle += 10)
        for (bool isVertical : {false, true})
            for (int sign : {-1, 1})
                power.push_back({isVertical, sign * angle, weights[angle / 10]});
    Cta2034Cuts result{find(horizontal, 0), average(listening), average(reflections), {},
                       average(power), {}, {}};
    if (result.listeningWindowDb && result.earlyReflectionsDb &&
        std::isfinite(*result.listeningWindowDb) && std::isfinite(*result.earlyReflectionsDb))
        result.earlyReflectionsDiDb = *result.listeningWindowDb - *result.earlyReflectionsDb;
    if (result.listeningWindowDb && result.soundPowerDb &&
        std::isfinite(*result.listeningWindowDb) && std::isfinite(*result.soundPowerDb))
        result.soundPowerDiDb = *result.listeningWindowDb - *result.soundPowerDb;
    if (result.listeningWindowDb && result.earlyReflectionsDb && result.soundPowerDb) {
        const double maximum = std::max({*result.listeningWindowDb, *result.earlyReflectionsDb,
                                         *result.soundPowerDb});
        result.predictedInRoomDb = std::isfinite(maximum)
            ? maximum + 10 * std::log10(.12 * std::pow(10., (*result.listeningWindowDb - maximum) / 10.) +
                 .44 * std::pow(10., (*result.earlyReflectionsDb - maximum) / 10.) +
                 .44 * std::pow(10., (*result.soundPowerDb - maximum) / 10.))
            : maximum;
    }
    return result;
}

double soundPressureLevelDecibels(std::complex<double> pressure, double referencePressure) {
    if (!std::isfinite(pressure.real()) || !std::isfinite(pressure.imag()) ||
        !std::isfinite(referencePressure) || referencePressure <= 0.0)
        throw std::invalid_argument("SPL requires finite pressure and a positive reference.");
    const double scale = std::max(std::abs(pressure.real()), std::abs(pressure.imag()));
    if (scale == 0.0)
        return -std::numeric_limits<double>::infinity();
    const double scaledAmplitude = std::hypot(pressure.real() / scale, pressure.imag() / scale);
    return 20.0 * (std::log10(scale) + std::log10(scaledAmplitude) - std::log10(referencePressure));
}

std::complex<double> normalizedVolumeImpedance(std::complex<double> volumeImpedance,
                                               double sourceArea, double density,
                                               double soundSpeed) {
    if (!std::isfinite(volumeImpedance.real()) || !std::isfinite(volumeImpedance.imag()) ||
        !std::isfinite(sourceArea) || sourceArea <= 0.0 || !std::isfinite(density) ||
        density <= 0.0 || !std::isfinite(soundSpeed) || soundSpeed <= 0.0)
        throw std::invalid_argument("Impedance normalization requires finite positive SI data.");
    const auto result = volumeImpedance * ((sourceArea / density) / soundSpeed);
    if (!std::isfinite(result.real()) || !std::isfinite(result.imag()))
        throw std::invalid_argument("Normalized impedance exceeds the finite numeric range.");
    return result;
}

std::optional<double> angularLevelDecibels(std::span<const AcousticsAngularSample> samples,
                                           double angleDegrees) {
    if (!std::isfinite(angleDegrees))
        throw std::invalid_argument("Reference angle must be finite.");
    for (std::size_t i = 0; i < samples.size(); ++i)
        if (!std::isfinite(samples[i].angleDegrees) ||
            (i != 0 && samples[i].angleDegrees <= samples[i - 1].angleDegrees))
            throw std::invalid_argument("Angular samples must be finite and strictly increasing.");
    const auto right = std::lower_bound(
        samples.begin(), samples.end(), angleDegrees,
        [](const auto &sample, double angle) { return sample.angleDegrees < angle; });
    if (right == samples.end() || !std::isfinite(right->levelDecibels))
        return {};
    if (right->angleDegrees == angleDegrees)
        return right->levelDecibels;
    if (right == samples.begin() || !std::isfinite((right - 1)->levelDecibels))
        return {};
    const auto &left = *(right - 1);
    return std::lerp(left.levelDecibels, right->levelDecibels,
                     fraction(angleDegrees, left.angleDegrees, right->angleDegrees));
}

std::optional<Beamwidth> beamwidthDegrees(std::span<const AcousticsAngularSample> samples,
                                          double attenuationDb) {
    if (!std::isfinite(attenuationDb) || attenuationDb <= 0.0)
        throw std::invalid_argument("Beamwidth attenuation must be finite and positive.");
    const auto reference = angularLevelDecibels(samples, 0.0);
    if (!reference)
        return std::nullopt;
    const auto right =
        std::lower_bound(samples.begin(), samples.end(), 0.0, [](const auto &sample, double angle) {
            return sample.angleDegrees < angle;
        });
    if (right == samples.end())
        return std::nullopt;
    const auto index = static_cast<std::ptrdiff_t>(right - samples.begin());
    const bool exactAxis = right->angleDegrees == 0.0;
    const double axisLevel = *reference;
    const double threshold = axisLevel - attenuationDb;
    if (!std::isfinite(axisLevel) || !std::isfinite(threshold) || threshold >= axisLevel)
        return std::nullopt;
    const auto crossing = [&](std::ptrdiff_t first, std::ptrdiff_t step) -> std::optional<double> {
        AcousticsAngularSample previous{0.0, axisLevel};
        for (auto i = first; i >= 0 && i < static_cast<std::ptrdiff_t>(samples.size()); i += step) {
            const auto &current = samples[static_cast<std::size_t>(i)];
            // NaNs and infinite null depths cannot define a linear dB crossing.
            if (!std::isfinite(current.levelDecibels))
                return std::nullopt;
            if (current.levelDecibels <= threshold)
                return std::lerp(
                    previous.angleDegrees, current.angleDegrees,
                    fraction(threshold, previous.levelDecibels, current.levelDecibels));
            previous = current;
        }
        return std::nullopt;
    };
    const auto leftCrossing = crossing(index - 1, -1);
    const auto rightCrossing = crossing(exactAxis ? index + 1 : index, 1);
    if (!leftCrossing || !rightCrossing)
        return std::nullopt;
    const double width = *rightCrossing - *leftCrossing;
    if (!std::isfinite(width) || *leftCrossing >= 0.0 || *rightCrossing <= 0.0)
        return std::nullopt;
    return Beamwidth{*leftCrossing, *rightCrossing, width};
}

std::vector<AngularFrequencySample>
angularFrequencyResponse(std::span<const AngularFrequencyRow> rows, double angleDegrees,
                         AngularResponseNormalization normalization, double referenceAngleDegrees) {
    if (!std::isfinite(angleDegrees) || !std::isfinite(referenceAngleDegrees))
        throw std::invalid_argument("Response and reference angles must be finite.");
    if (normalization != AngularResponseNormalization::absolute &&
        normalization != AngularResponseNormalization::referenceAngle)
        throw std::invalid_argument("Unknown angular response normalization.");
    std::vector<AngularFrequencySample> response;
    response.reserve(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto &row = rows[i];
        if (!std::isfinite(row.frequencyHz) || row.frequencyHz <= 0 ||
            (i != 0 && row.frequencyHz <= rows[i - 1].frequencyHz))
            throw std::invalid_argument(
                "Response frequencies must be positive and strictly increasing.");
        auto level = angularLevelDecibels(row.samples, angleDegrees);
        if (level && normalization == AngularResponseNormalization::referenceAngle) {
            const auto reference = angularLevelDecibels(row.samples, referenceAngleDegrees);
            if (reference && std::isfinite(*level - *reference))
                level = *level - *reference;
            else
                level.reset();
        }
        response.push_back({row.frequencyHz, level});
    }
    return response;
}
} // namespace hornsculpt
