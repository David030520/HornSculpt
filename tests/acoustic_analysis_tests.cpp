#include <hornsculpt/acoustic_analysis.hpp>
#include <hornsculpt/acoustic_frame.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
void near(double actual, double expected) {
    require(std::isfinite(actual) && std::abs(actual - expected) < 1e-10,
            "The acoustic value does not match its analytical reference.");
}
template <class Function> void rejects(Function function) {
    bool rejected = false;
    try {
        function();
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "An invalid acoustic analysis input was accepted.");
}
void acousticFrames() {
    const auto nearVector = [](Vec3 actual, Vec3 expected) {
        near(actual.x, expected.x);
        near(actual.y, expected.y);
        near(actual.z, expected.z);
    };
    const auto standard = defaultAcousticFrame();
    nearVector(standard.origin, {});
    nearVector(standard.forward, {0, 0, 1});
    nearVector(standard.horizontal, {1, 0, 0});
    nearVector(verticalDirection(standard), {0, 1, 0});
    for (const double sign : {-1.0, 1.0}) {
        const auto alongX = defaultAcousticFrame({.3, -.2, .1}, {2 * sign, 0, 0});
        nearVector(alongX.origin, {.3, -.2, .1});
        nearVector(alongX.forward, {sign, 0, 0});
        nearVector(alongX.horizontal, {0, 1, 0});
        nearVector(verticalDirection(alongX), {0, 0, sign});
    }
    const auto oblique = normalizedAcousticFrame({{1, 2, 3}, {0, 3, 4}, {4, 3, 4}});
    nearVector(oblique.origin, {1, 2, 3});
    nearVector(oblique.forward, {0, .6, .8});
    nearVector(oblique.horizontal, {1, 0, 0});
    nearVector(verticalDirection(oblique), {0, .8, -.6});
    const auto repeated = normalizedAcousticFrame(oblique);
    nearVector(repeated.forward, oblique.forward);
    nearVector(repeated.horizontal, oblique.horizontal);
    for (const double scale :
         {std::numeric_limits<double>::denorm_min(), 1.0, std::numeric_limits<double>::max()}) {
        const auto rotated = normalizedAcousticFrame({{1, 2, 3}, {0, -scale, 0}, {0, 0, scale}});
        nearVector(rotated.forward, {0, -1, 0});
        nearVector(rotated.horizontal, {0, 0, 1});
        nearVector(verticalDirection(rotated), {-1, 0, 0});
    }
    for (const auto bad :
         {AcousticFrame{{}, {}, {1, 0, 0}}, AcousticFrame{{}, {0, 0, 1}, {}},
          AcousticFrame{{}, {0, 0, 1}, {0, 0, 2}}, AcousticFrame{{}, {0, 0, 1}, {1e-10, 0, -1}}}) {
        rejects([&] { (void)normalizedAcousticFrame(bad); });
        rejects([&] { (void)verticalDirection(bad); });
    }
    for (const double invalid :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        rejects([&] { (void)normalizedAcousticFrame({{invalid, 0, 0}}); });
        rejects([&] { (void)normalizedAcousticFrame({{}, {0, invalid, 1}}); });
        rejects([&] { (void)normalizedAcousticFrame({{}, {0, 0, 1}, {1, 0, invalid}}); });
        rejects([&] { (void)defaultAcousticFrame({}, {0, 0, invalid}); });
    }
    rejects([] { (void)defaultAcousticFrame({}, {}); });
}
void analyticThresholds() {
    std::vector<AcousticsAngularSample> triangle;
    for (int angle = -90; angle <= 90; angle += 15)
        triangle.push_back({static_cast<double>(angle), 80.0 - std::abs(angle) / 5.0});
    for (const double attenuation : {3.0, 6.0, 10.0}) {
        const auto result = beamwidthDegrees(triangle, attenuation);
        require(result.has_value(), "A closed analytical lobe was not measured.");
        near(result->leftDegrees, -5.0 * attenuation);
        near(result->rightDegrees, 5.0 * attenuation);
        near(result->widthDegrees, 10.0 * attenuation);
    }
    const auto defaultWidth = beamwidthDegrees(triangle);
    near(defaultWidth->widthDegrees, 60.0);
    for (auto &sample : triangle)
        sample.levelDecibels -= 123.0;
    near(beamwidthDegrees(triangle)->widthDegrees, defaultWidth->widthDegrees);
}
void physicalNormalization() {
    near(soundPressureLevelDecibels({20e-6, 0.0}), 0.0);
    near(soundPressureLevelDecibels({0.0, 20e-6}), 0.0);
    near(soundPressureLevelDecibels({1.0, 0.0}), 93.97940008672037);
    require(std::isinf(soundPressureLevelDecibels({0.0, 0.0})) &&
                soundPressureLevelDecibels({0.0, 0.0}) < 0.0,
            "A pressure null was replaced by an artificial SPL floor.");
    constexpr double density = 1.21, soundSpeed = 343.0, area = .0005;
    const double reference = density * soundSpeed / area;
    const auto impedance =
        normalizedVolumeImpedance({reference, -2.0 * reference}, area, density, soundSpeed);
    near(impedance.real(), 1.0);
    near(impedance.imag(), -2.0);
    rejects([] { (void)soundPressureLevelDecibels({1.0, 0.0}, 0.0); });
    rejects([] { (void)normalizedVolumeImpedance({1.0, 0.0}, 0.0, 1.21, 343.0); });
}
void arbitraryReferenceAngles() {
    const std::vector<AcousticsAngularSample> samples{{-180, 40}, {-30, 80}, {0, 90},
                                                      {10, 93},   {30, 85},  {180, 40}};
    near(*angularLevelDecibels(samples, 10), 93);
    near(*angularLevelDecibels(samples, 5), 91.5);
    near(*angularLevelDecibels(samples, -15), 85);
    near(*angularLevelDecibels(samples, -180), 40);
    near(*angularLevelDecibels(samples, 180), 40);
    require(!angularLevelDecibels(samples, 181), "Reference angle was extrapolated.");
    auto missing = samples;
    missing[3].levelDecibels = std::numeric_limits<double>::quiet_NaN();
    require(!angularLevelDecibels(missing, 5), "Reference interpolation crossed missing data.");
    near(*angularLevelDecibels(missing, 0), 90);
    rejects([&] { (void)angularLevelDecibels(samples, std::numeric_limits<double>::infinity()); });
}
void angularFrequencyResponses() {
    const std::vector<AcousticsAngularSample> low{{-10, 80}, {-5, 86}, {0, 90}, {5, 94}, {10, 84}};
    const std::vector<AcousticsAngularSample> middle{
        {-10, 60}, {-5, 70}, {0, 85}, {5, 81}, {10, 70}};
    const std::vector<AcousticsAngularSample> high{{-10, 75}, {-5, 80}, {0, 89}, {5, 80}, {10, 70}};
    const std::vector<AngularFrequencyRow> rows{{500, low}, {1000, middle}, {2000, high}};
    const auto absolute = angularFrequencyResponse(rows, 2.5);
    const auto relative =
        angularFrequencyResponse(rows, 2.5, AngularResponseNormalization::referenceAngle);
    const auto offAxis =
        angularFrequencyResponse(rows, 2.5, AngularResponseNormalization::referenceAngle, -2.5);
    const double expectedAbsolute[]{92, 83, 84.5}, expectedRelative[]{2, -2, -4.5},
        expectedOffAxis[]{4, 5.5, 0};
    require(absolute.size() == rows.size() && relative.size() == rows.size() &&
                offAxis.size() == rows.size(),
            "Angular response lost an actual frequency row.");
    for (std::size_t i = 0; i < rows.size(); ++i) {
        near(absolute[i].frequencyHz, rows[i].frequencyHz);
        require(absolute[i].levelDecibels && relative[i].levelDecibels && offAxis[i].levelDecibels,
                "A supported angular response sample was omitted.");
        near(*absolute[i].levelDecibels, expectedAbsolute[i]);
        near(*relative[i].levelDecibels, expectedRelative[i]);
        near(*offAxis[i].levelDecibels, expectedOffAxis[i]);
    }
    // Only angular interpolation is allowed; an absent frequency stays absent.
    const std::vector<AngularFrequencyRow> sparse{{500, low}, {2000, high}};
    const auto sparseResponse = angularFrequencyResponse(sparse, 0);
    require(sparseResponse.size() == 2 && sparseResponse[0].frequencyHz == 500 &&
                sparseResponse[1].frequencyHz == 2000,
            "Angular extraction invented frequency samples.");
    near(*sparseResponse[0].levelDecibels, 90);
    near(*sparseResponse[1].levelDecibels, 89);
    const auto exactReference =
        angularFrequencyResponse(rows, -2.5, AngularResponseNormalization::referenceAngle, -2.5);
    for (const auto &sample : exactReference)
        near(*sample.levelDecibels, 0);
    require(angularFrequencyResponse({}, 0).empty(), "Empty frequency data invented a response.");
}
void angularFrequencyGaps() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();
    const std::vector<AcousticsAngularSample> valid{{-5, 84}, {0, 90}, {5, 86}, {10, 80}};
    const std::vector<AcousticsAngularSample> missingReference{
        {-5, 84}, {0, nan}, {5, 86}, {10, 80}};
    const std::vector<AcousticsAngularSample> zeroReference{
        {-5, 84}, {0, -infinity}, {5, 86}, {10, 80}};
    const std::vector<AcousticsAngularSample> missingTarget{{-5, 84}, {0, 90}, {5, nan}, {10, 80}};
    const std::vector<AcousticsAngularSample> infiniteTarget{
        {-5, 84}, {0, 90}, {5, infinity}, {10, 80}};
    const std::vector<AngularFrequencyRow> rows{{100, valid},          {200, missingReference},
                                                {300, zeroReference},  {400, missingTarget},
                                                {500, infiniteTarget}, {600, {}},
                                                {700, valid}};
    const auto absolute = angularFrequencyResponse(rows, 5);
    const auto relative =
        angularFrequencyResponse(rows, 5, AngularResponseNormalization::referenceAngle);
    require(absolute.size() == 7 && relative.size() == 7,
            "Missing angular values removed frequency rows.");
    for (std::size_t i = 0; i < rows.size(); ++i) {
        near(relative[i].frequencyHz, rows[i].frequencyHz);
        if (i == 0 || i == 6) {
            require(relative[i].levelDecibels.has_value(), "Valid normalized data was lost.");
            near(*relative[i].levelDecibels, -4);
        } else
            require(!relative[i].levelDecibels,
                    "Missing target/reference fabricated a normalized response.");
    }
    near(*absolute[1].levelDecibels, 86);
    near(*absolute[2].levelDecibels, 86);
    for (std::size_t i : {3u, 4u, 5u})
        require(!absolute[i].levelDecibels, "Invalid absolute level was fabricated.");
    const auto crossedTarget = angularFrequencyResponse(rows, 7.5);
    const auto crossedReference =
        angularFrequencyResponse(rows, 5, AngularResponseNormalization::referenceAngle, 2.5);
    require(!crossedTarget[3].levelDecibels && !crossedReference[1].levelDecibels &&
                !crossedReference[2].levelDecibels,
            "Angular response interpolated across missing data.");
    for (const auto &sample : angularFrequencyResponse(rows, 15))
        require(!sample.levelDecibels, "An angular response was extrapolated.");
    for (const auto &sample :
         angularFrequencyResponse(rows, 5, AngularResponseNormalization::referenceAngle, 15))
        require(!sample.levelDecibels, "A missing reference angle fabricated normalized data.");
    const std::vector<AcousticsAngularSample> extreme{{0, -std::numeric_limits<double>::max()},
                                                      {5, std::numeric_limits<double>::max()}};
    const std::vector<AngularFrequencyRow> overflow{{1000, extreme}};
    require(!angularFrequencyResponse(overflow, 5, AngularResponseNormalization::referenceAngle)[0]
                 .levelDecibels,
            "Overflowing relative dB values escaped as infinity.");
    const std::vector<AngularFrequencyRow> zero{{1000, zeroReference}};
    require(!angularFrequencyResponse(zero, 0, AngularResponseNormalization::referenceAngle)[0]
                 .levelDecibels,
            "A null at the reference angle was converted into a fabricated zero-dB curve.");
}
void invalidAngularFrequencyInputs() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();
    const std::vector<AcousticsAngularSample> samples{{0, 90}, {5, 86}};
    for (double bad : {-1.0, 0.0, nan, infinity}) {
        const std::vector<AngularFrequencyRow> rows{{bad, samples}};
        rejects([&] { (void)angularFrequencyResponse(rows, 0); });
    }
    for (const std::vector<AngularFrequencyRow> rows :
         {std::vector<AngularFrequencyRow>{{100, samples}, {100, samples}},
          std::vector<AngularFrequencyRow>{{200, samples}, {100, samples}}})
        rejects([&] { (void)angularFrequencyResponse(rows, 0); });
    const std::vector<AngularFrequencyRow> rows{{100, samples}};
    for (double bad : {nan, infinity}) {
        rejects([&] { (void)angularFrequencyResponse(rows, bad); });
        rejects([&] {
            (void)angularFrequencyResponse(rows, 0, AngularResponseNormalization::referenceAngle,
                                           bad);
        });
    }
    rejects([&] {
        (void)angularFrequencyResponse(rows, 0, static_cast<AngularResponseNormalization>(100));
    });
    const std::vector<AcousticsAngularSample> unordered{{5, 86}, {0, 90}};
    const std::vector<AngularFrequencyRow> invalid{{100, unordered}};
    rejects([&] { (void)angularFrequencyResponse(invalid, 0); });
}
void asymmetryAndSideLobes() {
    const std::vector<AcousticsAngularSample> samples{{-90, 3}, {-60, -12}, {-20, -4}, {0, 0},
                                                      {15, -3}, {45, -12},  {70, 2},   {90, -12}};
    const auto result = beamwidthDegrees(samples);
    require(result.has_value(), "An asymmetric lobe was not measured.");
    near(result->leftDegrees, -30.0);
    near(result->rightDegrees, 25.0);
    near(result->widthDegrees, 55.0);
    const std::vector<AcousticsAngularSample> offAxisPeak{
        {-40, -12}, {-10, -2}, {0, 0}, {10, 2}, {40, -12}};
    const auto shifted = beamwidthDegrees(offAxisPeak);
    require(shifted.has_value(), "A peak shifted within the on-axis lobe was rejected.");
    near(shifted->leftDegrees, -22.0);
    near(shifted->rightDegrees, 10.0 + 30.0 * 8.0 / 14.0);
    const std::vector<AcousticsAngularSample> plateau{{-60, -12}, {-30, -6}, {-20, -6}, {0, 0},
                                                      {20, -6},   {30, -6},  {60, -12}};
    near(beamwidthDegrees(plateau)->widthDegrees, 40.0);
}
void axisInterpolationAndResampling() {
    const std::vector<AcousticsAngularSample> coarse{{-60, -12}, {-10, 0}, {10, 0}, {60, -12}};
    const auto original = beamwidthDegrees(coarse);
    require(original.has_value(), "A bracketed on-axis value was not interpolated.");
    near(original->widthDegrees, 70.0);
    std::vector<AcousticsAngularSample> refined;
    for (std::size_t i = 0; i + 1 < coarse.size(); ++i) {
        refined.push_back(coarse[i]);
        refined.push_back({.5 * (coarse[i].angleDegrees + coarse[i + 1].angleDegrees),
                           .5 * (coarse[i].levelDecibels + coarse[i + 1].levelDecibels)});
    }
    refined.push_back(coarse.back());
    near(beamwidthDegrees(refined)->widthDegrees, original->widthDegrees);
}
void missingAndInvalidData() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();
    require(!beamwidthDegrees({}), "Empty data fabricated a beamwidth.");
    for (const std::vector<AcousticsAngularSample> samples :
         {std::vector<AcousticsAngularSample>{{-20, -3}, {0, 0}, {20, -3}},
          std::vector<AcousticsAngularSample>{{10, 0}, {30, -12}},
          std::vector<AcousticsAngularSample>{{-20, -12}, {0, nan}, {20, -12}},
          std::vector<AcousticsAngularSample>{{-30, -12}, {-10, nan}, {0, 0}, {30, -12}},
          std::vector<AcousticsAngularSample>{{-30, -12}, {0, 0}, {10, -infinity}},
          std::vector<AcousticsAngularSample>{{-30, -12}, {-1, nan}, {1, 0}, {30, -12}}})
        require(!beamwidthDegrees(samples), "Missing data or an open lobe fabricated a crossing.");
    const std::vector<AcousticsAngularSample> outsideGap{
        {-90, nan}, {-30, -12}, {0, 0}, {30, -12}, {90, nan}};
    near(beamwidthDegrees(outsideGap)->widthDegrees, 30.0);
    const std::vector<AcousticsAngularSample> valid{{-30, -12}, {0, 0}, {30, -12}};
    for (const double attenuation : {0.0, -6.0, nan, infinity})
        rejects([&] { (void)beamwidthDegrees(valid, attenuation); });
    for (const std::vector<AcousticsAngularSample> invalid :
         {std::vector<AcousticsAngularSample>{{0, 0}, {-10, -12}},
          std::vector<AcousticsAngularSample>{{0, 0}, {0, -12}},
          std::vector<AcousticsAngularSample>{{-30, -12}, {nan, 0}},
          std::vector<AcousticsAngularSample>{{-infinity, -12}, {0, 0}}})
        rejects([&] { (void)beamwidthDegrees(invalid); });
}
void ctaCuts() {
    std::vector<AcousticsAngularSample> horizontal, vertical;
    for (int angle = -180; angle <= 180; angle += 5) {
        horizontal.push_back({double(angle), 80});
        vertical.push_back({double(angle), 80});
    }
    auto curves = cta2034Cuts(horizontal, vertical);
    near(*curves.onAxisDb, 80);
    near(*curves.listeningWindowDb, 80);
    near(*curves.earlyReflectionsDb, 80);
    near(*curves.earlyReflectionsDiDb, 0);
    near(*curves.soundPowerDb, 80);
    near(*curves.soundPowerDiDb, 0);
    near(*curves.predictedInRoomDb, 80);
    vertical.front().levelDecibels = 120;
    vertical.back().levelDecibels = 120;
    vertical[36].levelDecibels = 120;
    near(*cta2034Cuts(horizontal, vertical).soundPowerDb, 80);
    vertical.front().levelDecibels = 80;
    vertical.back().levelDecibels = 80;
    vertical[36].levelDecibels = 80;
    // One listening direction at +10 dB changes the nine-point power mean.
    horizontal[38].levelDecibels = 90;
    curves = cta2034Cuts(horizontal, vertical);
    near(*curves.listeningWindowDb, 80 + 10 * std::log10(2.));
    // One of seven front-wall samples contributes 1/35 of reflection power.
    near(*curves.earlyReflectionsDb, 80 + 10 * std::log10(1 + 9. / 35));
    // Independently calculated from the published Table 7 sum 0.996373088,
    // one 90 dB sample on a uniform 80 dB field.
    near(*curves.soundPowerDb, 80.181704804004);
    near(*curves.soundPowerDiDb, *curves.listeningWindowDb - *curves.soundPowerDb);
    near(*curves.predictedInRoomDb,
         80 + 10 * std::log10(.12 * std::pow(10., (*curves.listeningWindowDb - 80) / 10) +
                                .44 * std::pow(10., (*curves.earlyReflectionsDb - 80) / 10) +
                                .44 * std::pow(10., (*curves.soundPowerDb - 80) / 10)));
    horizontal.pop_back();
    curves = cta2034Cuts(horizontal, vertical);
    require(curves.listeningWindowDb.has_value() && !curves.earlyReflectionsDb &&
                !curves.earlyReflectionsDiDb && !curves.soundPowerDb && !curves.predictedInRoomDb,
            "Missing rear direction must not fabricate early reflections.");
}
void ctaPublishedWeighting() {
    // Table 7, one 80 dB direction against true zero at all other directions.
    // Golden values use the independent 70-point weight sum 0.996373088.
    constexpr double expected[] = {
        47.829642582, 56.764565051, 59.536449109, 60.945565601, 61.773683737,
        62.286455686, 62.608364126, 62.804857365, 62.911176235, 62.944844389,
        62.911176235, 62.804857365, 62.608364126, 62.286455686, 61.773683737,
        60.945565601, 59.536449109, 56.764565051, 47.829642582};
    const double zero = -std::numeric_limits<double>::infinity();
    for (int angle = 0; angle <= 180; angle += 10) {
        for (bool isVertical : {false, true}) {
            if (isVertical && (angle == 0 || angle == 180))
                continue;
            for (int sign : {-1, 1}) {
                if ((angle == 0 || angle == 180) && sign < 0)
                    continue;
                std::vector<AcousticsAngularSample> h, v;
                for (int a = -180; a <= 180; a += 5) {
                    h.push_back({double(a), zero});
                    v.push_back({double(a), zero});
                }
                auto &plane = isVertical ? v : h;
                plane[static_cast<std::size_t>((sign * angle + 180) / 5)].levelDecibels = 80;
                const auto result = cta2034Cuts(h, v);
                require(result.soundPowerDb && std::abs(*result.soundPowerDb - expected[angle / 10]) < 1e-8,
                        "CTA sound power differs from the independently derived Table 7 value.");
            }
        }
    }
}
} // namespace
int main() {
    try {
        acousticFrames();
        analyticThresholds();
        physicalNormalization();
        arbitraryReferenceAngles();
        angularFrequencyResponses();
        angularFrequencyGaps();
        invalidAngularFrequencyInputs();
        asymmetryAndSideLobes();
        axisInterpolationAndResampling();
        missingAndInvalidData();
        ctaCuts();
        ctaPublishedWeighting();
        std::cout << "Acoustic frames, angular frequency responses, beamwidth, interpolation and "
                     "gap tests passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
