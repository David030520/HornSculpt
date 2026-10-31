#include <hornsculpt/axisymmetric_geometry.hpp>
#include <hornsculpt/draft_geometry.hpp>
#include <hornsculpt/geometry_analysis.hpp>
#include <hornsculpt/meshing.hpp>

#ifdef HORNSCULPT_SLOT_CUDA_TEST
#include <hornsculpt/bem.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

using namespace hornsculpt;
namespace {
constexpr double pi = std::numbers::pi;
constexpr double dimensionTolerance = 1e-11; // metres, analytic evaluation only.

void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const char *message) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) > tolerance)
        throw std::runtime_error(std::string(message) + ": actual=" + std::to_string(actual) +
                                 ", expected=" + std::to_string(expected));
}

template <class Function> void rejects(Function operation, const char *message) {
    try {
        operation();
    } catch (const std::invalid_argument &) {
        return;
    }
    throw std::runtime_error(message);
}

DraftParameters fixture() {
    DraftParameters p;
    p.entry = {.012, .012, 2};
    p.length = .06;
    p.wallThickness = .002;
    p.horizontalAngle = .45;
    p.verticalAngle = .30;
    p.profile.throatAngle = .10;
    p.verticalProfile.throatAngle = .15;
    p.slot.enabled = true;
    p.slot.aperture = {.006, .018, 4};
    p.slot.adapterLength = .012;
    p.slot.straightLength = .003;
    p.slot.transitionLength = .004;
    p.crossSectionExponent = 4;
    p.transition.startFraction = .4;
    p.transition.endFraction = .8;
    return p;
}

DraftParameters doubleRadiusFixture() {
    auto p = fixture();
    p.slot.contour = SlotContour::roundedRectangle;
    p.slot.cornerRadius = .001;
    p.slot.exit = SlotExit::circularFillet;
    p.slot.exitRadius = .001;
    p.profile.k = p.verticalProfile.k = 0;
    p.profile.terminationScale = p.verticalProfile.terminationScale = 0;
    return p;
}

double slotExit(const DraftParameters &p) {
    return p.slot.adapterLength + p.slot.straightLength;
}

double apertureArea(const EntryParameters &aperture) {
    const double gamma = std::tgamma(1.0 + 1.0 / aperture.exponent);
    return aperture.width * aperture.height * gamma * gamma /
           std::tgamma(1.0 + 2.0 / aperture.exponent);
}

void elongatedSlotLargeFillet() {
    DraftParameters p;
    p.slot.enabled = true;
    p.slot.aperture = {.011234567891234, .036789123456789, 4.25};
    p.slot.adapterLength = .023456789123456;
    p.slot.straightLength = .010;
    p.slot.contour = SlotContour::roundedRectangle;
    p.slot.cornerRadius = .002;
    p.slot.exit = SlotExit::circularFillet;
    p.slot.exitRadius = .018;
    p.transition.startFraction = .3;
    p.profile.k = p.verticalProfile.k = 0;
    p.profile.terminationScale = p.verticalProfile.terminationScale = 0;
    p.horizontalAngle = .5;
    p.verticalAngle = .3;
    const DraftGeometry geometry(p);
    // Near a moving rounded corner, the tangent-circle equation can be much
    // steeper than the wall slope. Axial bracket width alone is insufficient.
    for (const double degrees : {0., 30., 64., 64.5, 65., 72., 90.}) {
        const double phi = degrees * pi / 180;
        const auto fillet = *geometry.slotExitFillet(phi);
        near(fillet.radius, p.slot.exitRadius, 0, "Large fillet changed its physical radius");
        near(geometry.radius(fillet.end, phi),
             fillet.centreRadius - std::sqrt(fillet.radius * fillet.radius -
                                            std::pow(fillet.end - fillet.start, 2)),
             dimensionTolerance, "Large fillet misses the downstream profile");
        constexpr double tangentStep = 1e-9;
        near(geometry.radiusDerivative(fillet.end - tangentStep, phi),
             geometry.radiusDerivative(fillet.end + tangentStep, phi), 2e-5,
             "Large fillet lost the downstream tangent");
    }
}

void trueRadii() {
    const double a = .006, b = .0175, corner = .002;
    const RoundedRectangle contour(a, b, corner);
    near(contour.area(), 4 * a * b - (4 - pi) * corner * corner, 1e-15,
         "Rounded-rectangle analytic area");
    near(contour.radius(0), a, dimensionTolerance, "Rounded-rectangle half-width");
    near(contour.radius(pi / 2), b, dimensionTolerance, "Rounded-rectangle half-height");
    const double cx = a - corner, cy = b - corner;
    for (unsigned i = 1; i < 16; ++i) {
        const double theta = pi * i / 32;
        const double x = cx + corner * std::cos(theta), y = cy + corner * std::sin(theta);
        const double phi = std::atan2(y, x);
        near(contour.radius(phi), std::hypot(x, y), dimensionTolerance,
             "Corner follows the requested physical circle");
        const double h = 1e-6;
        near(contour.angularDerivative(phi),
             (contour.radius(phi + h) - contour.radius(phi - h)) / (2 * h), 1e-8,
             "Rounded-rectangle polar derivative");
    }
    const RoundedRectangle square(a, b, 0);
    near(square.radius(std::atan2(b, a)), std::hypot(a, b), dimensionTolerance,
         "Zero corner radius is an exact rectangle");
    const RoundedRectangle circle(a, a, a);
    for (double phi : {0., .3, .7, 1.2})
        near(circle.radius(phi), a, dimensionTolerance, "Maximum square corner radius is a circle");
    rejects([&] { (void)RoundedRectangle(a, b, a + .0001); },
            "Oversized corner radius was silently clamped");

    DraftParameters p;
    p.entry = {.0254, .0254, 2};
    p.length = .18;
    p.slot.enabled = true;
    p.slot.contour = SlotContour::roundedRectangle;
    p.slot.aperture = {.012, .012, 4};
    p.slot.cornerRadius = .006;
    p.slot.adapterLength = .02;
    p.slot.straightLength = .012;
    p.slot.exit = SlotExit::circularFillet;
    p.slot.exitRadius = .018;
    p.profile.k = p.verticalProfile.k = 0;
    p.profile.terminationScale = p.verticalProfile.terminationScale = 0;
    p.horizontalAngle = p.verticalAngle = .45;
    p.transition.startFraction = .4;
    const DraftGeometry geometry(p);
    const double m = std::tan(.45), R = p.slot.exitRadius;
    const double expectedStart = slotExit(p) - R * std::tan(.45 / 2);
    const double expectedEnd = slotExit(p) + R * (1 - std::cos(.45)) / m;
    for (double phi : {0., .3, .7, pi / 2, 2.6, 4.1}) {
        const auto fillet = geometry.slotExitFillet(phi);
        require(fillet.has_value(), "Active physical exit radius has no arc");
        near(fillet->radius, R, 0, "Requested exit radius was modified");
        near(fillet->start, expectedStart, dimensionTolerance, "Cone fillet start");
        near(fillet->end, expectedEnd, dimensionTolerance, "Cone fillet end");
        for (unsigned i = 0; i <= 10; ++i) {
            const double z = std::lerp(fillet->start, fillet->end, i / 10.);
            const double radial = geometry.radius(z, phi);
            near(std::hypot(z - fillet->start, radial - fillet->centreRadius),
                 R, dimensionTolerance, "Exit follows a true constant-radius circle");
        }
        near(geometry.radiusDerivative(fillet->start, phi), 0, 1e-10,
             "Fillet is tangent to straight slot");
        near(geometry.radiusDerivative(fillet->end, phi), m, 1e-9,
             "Fillet is tangent to unchanged OS-SE profile");
    }
    near(geometry.slotStraightEnd(), expectedStart, dimensionTolerance,
         "Straight-slot face must end before the complete fillet");
    near(geometry.slotTransitionEnd(), expectedEnd, dimensionTolerance,
         "Composition end tracks the actual fillet");
    const auto breaks = geometry.axialBreakpoints();
    require(std::none_of(breaks.begin(), breaks.end(), [&](double z) {
                return z == slotExit(p);
            }), "Nominal sharp exit still splits the circular fillet");
    require(std::any_of(breaks.begin(), breaks.end(), [&](double z) {
                return std::abs(z - expectedStart) <= dimensionTolerance;
            }), "Fillet start is absent from axial sampling");
    auto invalid = p;
    invalid.slot.straightLength = .001;
    rejects([&] { (void)DraftGeometry(invalid); }, "Fillet consumed the adapter silently");
    invalid = p;
    invalid.transition.startFraction = slotExit(p) / p.length;
    rejects([&] { (void)DraftGeometry(invalid); }, "Form transition overlapped the actual fillet");
    invalid = p;
    invalid.profile.k = invalid.verticalProfile.k = 1;
    invalid.profile.throatAngle = invalid.verticalProfile.throatAngle = 0;
    invalid.slot.exitRadius = .001;
    rejects([&] { (void)DraftGeometry(invalid); },
            "Trivial tangent point was falsely accepted as a circular fillet");

    p.slot.aperture.height = .035;
    p.slot.cornerRadius = .002;
    p.verticalAngle = .3;
    p.slot.exitRadius = .002;
    const DraftGeometry rectangular(p);
    near(diffractionSlotArea(p.slot), contour.area(), 1e-15,
         "Slot-area helper ignores true corner radius");
    for (double phi : {0., .3, .7, 1.2, pi / 2}) {
        near(rectangular.radius(p.slot.adapterLength, phi), contour.radius(phi),
             dimensionTolerance, "Adapter misses rounded-rectangle aperture");
        const auto fillet = *rectangular.slotExitFillet(phi);
        const double z = (fillet.start + fillet.end) / 2;
        near(std::hypot(z - fillet.start,
                        rectangular.radius(z, phi) - fillet.centreRadius),
             p.slot.exitRadius, dimensionTolerance, "Directional fillet changed its radius");
    }
    invalid = p;
    invalid.slot.cornerRadius = 0;
    bool rejectedProfileJump = false;
    try {
        (void)DraftGeometry(invalid);
    } catch (const std::invalid_argument &error) {
        rejectedProfileJump =
            std::string_view(error.what()).find("Tangentialgleichung") != std::string_view::npos;
    }
    require(rejectedProfileJump,
            "A moving sharp rectangle corner was accepted as a tangent-circle root");

    auto osse = p;
    osse.slot.straightLength = .010;
    osse.horizontalAngle = .5;
    osse.profile.k = osse.verticalProfile.k = 1;
    osse.profile.throatAngle = .1;
    osse.verticalProfile.throatAngle = .15;
    osse.profile.terminationScale = osse.verticalProfile.terminationScale = .12;
    const DraftGeometry osseRounded(osse);
    constexpr double tangentStep = 1e-8;
    for (double phi : {0., .3, .7, 1.2, pi / 2, 2.6, 4.1}) {
        const auto fillet = *osseRounded.slotExitFillet(phi);
        const double z = (fillet.start + fillet.end) / 2;
        near(std::hypot(z - fillet.start,
                        osseRounded.radius(z, phi) - fillet.centreRadius),
             osse.slot.exitRadius, dimensionTolerance, "Nonconical OS-SE exit radius");
        near(osseRounded.radiusDerivative(fillet.end - tangentStep, phi),
             osseRounded.radiusDerivative(fillet.end + tangentStep, phi), 2e-5,
             "Nonconical OS-SE fillet lost its downstream tangent");
    }
    p.guideEnabled = true;
    p.guide.position = .8;
    const DraftGeometry guided(p);
    for (unsigned i = 0; i < 48; ++i) {
        const double phi = 2 * pi * i / 48;
        near(guided.radius(guided.guidePosition(), phi), guided.guideRadius(phi),
             dimensionTolerance, "True slot radii broke the guide hit");
    }
    p.guideEnabled = false;
    p.slot.exitRadius = 0;
    const DraftGeometry sharp(p);
    require(!sharp.slotExitFillet(0) &&
                std::isnan(sharp.radiusDerivative(slotExit(p), 0)),
            "Zero exit radius did not restore the physical tangent jump");
}

void disabledEquivalence() {
    const DraftGeometry original;
    auto p = original.parameters();
    p.slot.aperture = {.004, .042, 8};
    p.slot.adapterLength = .03;
    p.slot.straightLength = .009;
    p.slot.transitionLength = 0;
    const DraftGeometry disabled(p);
    require(!DraftParameters{}.slot.enabled, "The diffraction slot must be opt-in.");
    for (int axial = 0; axial <= 20; ++axial)
        for (int angle = 0; angle < 32; ++angle) {
            const double z = original.length() * axial / 20.0;
            const double phi = 2 * pi * angle / 32.0;
            require(original.radius(z, phi) == disabled.radius(z, phi),
                    "Inactive slot parameters changed the existing horn.");
            require(original.radiusDerivative(z, phi) == disabled.radiusDerivative(z, phi),
                    "Inactive slot parameters changed the existing tangents.");
        }
    require(disabled.axialBreakpoints() == std::vector<double>({0, p.length}),
            "Disabled slot introduced extra axial boundaries.");
}

void dimensionsAndSegments() {
    const auto p = fixture();
    const DraftGeometry geometry(p);
    const Superellipse entry(p.entry.width / 2, p.entry.height / 2, p.entry.exponent);
    const Superellipse aperture(p.slot.aperture.width / 2, p.slot.aperture.height / 2,
                                p.slot.aperture.exponent);
    const double exit = slotExit(p);
    near(geometry.length(), p.length, 0, "The slot changed total horn length");
    for (int angle = 0; angle < 64; ++angle) {
        const double phi = 2 * pi * angle / 64.0;
        near(geometry.radius(0, phi), entry.radius(phi), dimensionTolerance,
             "The source connection changed");
        for (int axial = 0; axial <= 10; ++axial) {
            const double z = p.slot.adapterLength + p.slot.straightLength * axial / 10.0;
            near(geometry.radius(z, phi), aperture.radius(phi), dimensionTolerance,
                 "The straight slot does not preserve its aperture");
        }
        near(geometry.radius(p.slot.adapterLength / 2, phi),
             (entry.radius(phi) + aperture.radius(phi)) / 2, dimensionTolerance,
             "The symmetric adapter midpoint is incorrect");
        const double lower = std::min(entry.radius(phi), aperture.radius(phi));
        const double upper = std::max(entry.radius(phi), aperture.radius(phi));
        for (int axial = 0; axial <= 20; ++axial) {
            const double radius = geometry.radius(p.slot.adapterLength * axial / 20.0, phi);
            require(radius >= lower - dimensionTolerance && radius <= upper + dimensionTolerance,
                    "The adapter overshoots its endpoint radii.");
        }
    }
    const std::vector<double> expected{0, p.slot.adapterLength, exit,
                                       exit + p.slot.transitionLength, p.length};
    const auto breaks = geometry.axialBreakpoints();
    require(breaks.size() == expected.size(), "A physical axial segment boundary is missing.");
    for (std::size_t i = 0; i < breaks.size(); ++i)
        near(breaks[i], expected[i], dimensionTolerance, "Unexpected axial segment boundary");

    auto remaining = p;
    remaining.slot.enabled = false;
    remaining.entry = p.slot.aperture;
    remaining.length = p.length - exit;
    const DraftGeometry flare(remaining);
    for (int angle = 0; angle < 16; ++angle) {
        const double phi = 2 * pi * angle / 16.0;
        for (int axial = 0; axial <= 10; ++axial) {
            const double u = axial == 10
                                 ? remaining.length
                                 : p.slot.transitionLength +
                                       (remaining.length - p.slot.transitionLength) * axial / 10.0;
            const double z = axial == 10 ? p.length : exit + u;
            near(geometry.radius(z, phi), flare.radius(u, phi), dimensionTolerance,
                 "OS-SE does not use the remaining physical length");
            near(geometry.radiusDerivative(z, phi), flare.radiusDerivative(u, phi), 1e-9,
                 "The downstream OS-SE tangent changed outside the local transition");
        }
    }
}

void tangentsAndAnalysis() {
    auto p = fixture();
    const DraftGeometry smooth(p);
    const double exit = slotExit(p);
    // These distances resolve the finite transition without relying on tessellation.
    constexpr double derivativeStep = 1e-7;
    for (const double phi : {0.0, pi / 4, pi / 2}) {
        near(smooth.radiusDerivative(0, phi), 0, 1e-10, "Adapter inlet tangent is not axial");
        near(smooth.radiusDerivative(p.slot.adapterLength, phi), 0, 1e-10,
             "Adapter and straight slot are not tangent");
        near(smooth.radiusDerivative(exit, phi), 0, 1e-10,
             "The smoothed exit does not meet the straight slot tangentially");
        for (const double z : {p.slot.adapterLength, exit, exit + p.slot.transitionLength}) {
            const double before = smooth.radiusDerivative(z - derivativeStep, phi);
            const double at = smooth.radiusDerivative(z, phi);
            const double after = smooth.radiusDerivative(z + derivativeStep, phi);
            near(before, after, 2e-4, "A smooth segment boundary has a tangent jump");
            // One-sided derivative differences estimate radial curvature in 1/m.
            near((at - before) / derivativeStep, (after - at) / derivativeStep, .1,
                 "A smooth segment boundary is not C2");
        }
        for (const double z : {.004, .013, .017, .028, .048})
            near(smooth.radiusDerivative(z, phi),
                 (smooth.radius(z + derivativeStep, phi) - smooth.radius(z - derivativeStep, phi)) /
                     (2 * derivativeStep),
                 2e-6, "The analytical wall tangent disagrees with the actual geometry");
    }

    p.slot.transitionLength = 0;
    const DraftGeometry sharp(p);
    for (const double phi : {0.0, pi / 2}) {
        require(std::isnan(sharp.radiusDerivative(exit, phi)),
                "An actual sharp exit must not report a single wall tangent.");
        near(sharp.radiusDerivative(exit, phi, AxialDerivativeSide::before), 0, 1e-10,
             "The sharp exit upstream tangent is not axial");
        const double throatAngle = phi == 0 ? p.profile.throatAngle : p.verticalProfile.throatAngle;
        near(sharp.radiusDerivative(exit, phi, AxialDerivativeSide::after), std::tan(throatAngle),
             1e-9, "The sharp exit downstream tangent is incorrect");
        near(sharp.radius(exit - derivativeStep, phi), sharp.radius(exit, phi), dimensionTolerance,
             "The sharp exit has a radial step");
    }
    require(sharp.axialBreakpoints().size() == 4,
            "Zero transition length left duplicated axial breakpoints.");

    const auto sharpAnalysis = analyzeGeometry(sharp, 121);
    std::vector<GeometrySample> edgeSamples;
    for (const auto &sample : sharpAnalysis.samples)
        if (std::abs(sample.z - exit) < 1e-12)
            edgeSamples.push_back(sample);
    require(edgeSamples.size() == 3,
            "Sharp-edge analysis must expose both tangents and their discontinuity.");
    near(edgeSamples[0].horizontalWallAngle, 0, 1e-10, "Missing upstream edge angle");
    require(std::isnan(edgeSamples[1].horizontalWallAngle) &&
                std::isnan(edgeSamples[1].verticalWallAngle) &&
                std::isnan(edgeSamples[1].areaDerivative),
            "The analysis connects incompatible derivatives through a sharp edge.");
    near(edgeSamples[2].horizontalWallAngle, p.profile.throatAngle, 1e-9,
         "Missing downstream horizontal edge angle");
    near(edgeSamples[2].verticalWallAngle, p.verticalProfile.throatAngle, 1e-9,
         "Missing downstream vertical edge angle");
    near(edgeSamples[0].volume, edgeSamples[1].volume, 0,
         "An edge marker introduced artificial acoustic volume");
    near(edgeSamples[1].volume, edgeSamples[2].volume, 0,
         "One-sided edge samples introduced artificial acoustic volume");

    auto tangent = p;
    tangent.profile.throatAngle = 0;
    tangent.verticalProfile.throatAngle = 0;
    near(DraftGeometry(tangent).radiusDerivative(exit, 0), 0, 1e-10,
         "A tangential zero-length transition was incorrectly marked as a sharp edge");
    tangent.slot.aperture.width = .001;
    tangent.slot.adapterLength = .02;
    tangent.horizontalAngle = 80 * pi / 180;
    near(DraftGeometry(tangent).radiusDerivative(slotExit(tangent), 0), 0, 0,
         "The exact tangent limit was replaced by a nearby nonzero-curvature sample");

    const auto analysis = analyzeGeometry(smooth, 121);
    const auto sampleAt = [&](double position) -> const GeometrySample & {
        const auto found =
            std::find_if(analysis.samples.begin(), analysis.samples.end(),
                         [&](const auto &sample) { return std::abs(sample.z - position) < 1e-12; });
        require(found != analysis.samples.end(), "Analysis omitted a requested axial sample.");
        return *found;
    };
    near(sampleAt(0).area, apertureArea(p.entry), 1e-10, "The physical source area changed");
    near(sampleAt(.013).area, apertureArea(p.slot.aperture), 1e-9,
         "Slot cross-sectional analysis is incorrect");
    near(sampleAt(.013).areaDerivative, 0, 1e-9, "The straight slot has nonzero area gradient");
    near(sampleAt(.014).volume - sampleAt(.013).volume, apertureArea(p.slot.aperture) * .001, 1e-12,
         "Volume accumulation does not use the actual straight-slot cross section");
}

void zeroLengthAndCovariance() {
    auto p = fixture();
    p.slot.straightLength = 0;
    const DraftGeometry zeroLength(p);
    require(zeroLength.axialBreakpoints().size() == 4,
            "A zero-length straight slot produced duplicate boundaries.");
    near(zeroLength.radius(p.slot.adapterLength, 0), p.slot.aperture.width / 2, dimensionTolerance,
         "The zero-length slot lost its aperture width");
    near(zeroLength.radius(p.slot.adapterLength, pi / 2), p.slot.aperture.height / 2,
         dimensionTolerance, "The zero-length slot lost its aperture height");
    p.slot.transitionLength = 0;
    require(DraftGeometry(p).axialBreakpoints().size() == 3,
            "Coincident adapter, slot and exit boundaries were not merged.");

    p = fixture();
    p.slot.aperture = p.entry;
    const DraftGeometry matching(p);
    for (int i = 0; i <= 30; ++i)
        for (const double phi : {0., pi / 4, pi / 2}) {
            const double z = slotExit(p) * i / 30.0;
            near(matching.radius(z, phi), p.entry.width / 2, dimensionTolerance,
                 "Identical inlet and slot contours created a spurious adapter deformation");
            near(matching.radiusDerivative(z, phi), 0, 1e-10,
                 "Identical inlet and slot contours created a spurious wall angle");
        }

    p = fixture();
    const DraftGeometry original(p);
    std::swap(p.entry.width, p.entry.height);
    std::swap(p.slot.aperture.width, p.slot.aperture.height);
    std::swap(p.profile, p.verticalProfile);
    std::swap(p.horizontalAngle, p.verticalAngle);
    const DraftGeometry swapped(p);
    for (int axial = 0; axial <= 20; ++axial)
        for (int angular = 0; angular < 32; ++angular) {
            const double z = p.length * axial / 20.;
            const double phi = 2 * pi * angular / 32.;
            near(swapped.radius(z, phi), original.radius(z, phi + pi / 2), dimensionTolerance,
                 "Exchanging X/Y changes the physical slot instead of rotating it");
            near(swapped.radiusDerivative(z, phi), original.radiusDerivative(z, phi + pi / 2), 1e-9,
                 "X/Y exchange does not rotate the wall tangents consistently");
        }
}

void compositionAndRejections() {
    auto touching = fixture();
    touching.length = .16;
    touching.slot.adapterLength = .02;
    touching.slot.straightLength = .003;
    touching.slot.transitionLength = .005;
    const DraftGeometry beforeMorph(touching);
    touching.mouthMorph.enabled = true;
    const double finish = slotExit(touching) + touching.slot.transitionLength;
    touching.mouthMorph.startFraction = finish / touching.length;
    const DraftGeometry touchingMorph(touching);
    near(touchingMorph.radius(finish, 0), beforeMorph.radius(finish, 0), dimensionTolerance,
         "A morph exactly at the inlet boundary was rejected or changed the inlet");
    auto p = fixture();
    p.guideEnabled = true;
    p.guide.position = .85;
    p.guide.width = .12;
    p.guide.aspectRatio = .8;
    p.guide.rotation = .2;
    p.guide.superellipseExponent = 4;
    const DraftGeometry guided(p);
    for (int i = 0; i < 64; ++i) {
        const double phi = 2 * pi * i / 64.0;
        near(guided.radius(guided.guidePosition(), phi), guided.guideRadius(phi), 1e-9,
             "The slot composition misses the requested guide contour");
    }

    auto invalid = [](auto change, const char *message) {
        auto parameters = fixture();
        change(parameters);
        rejects([&] { (void)DraftGeometry(parameters); }, message);
    };
    invalid([](auto &q) { q.slot.aperture.width = 0; }, "Zero slot width was accepted.");
    invalid([](auto &q) { q.slot.aperture.height = -.01; }, "Negative slot height was accepted.");
    invalid([](auto &q) { q.slot.aperture.exponent = .5; }, "A nonconvex slot was accepted.");
    invalid([](auto &q) { q.slot.adapterLength = -.01; }, "Negative adapter length was accepted.");
    invalid([](auto &q) { q.slot.straightLength = -.001; },
            "Negative straight length was accepted.");
    invalid([](auto &q) { q.slot.adapterLength = q.length; }, "A slot without flare was accepted.");
    invalid([](auto &q) { q.slot.transitionLength = q.length; },
            "A transition outside the remaining flare was accepted.");
    invalid([](auto &q) { q.slot.transitionLength = std::numeric_limits<double>::infinity(); },
            "An infinite slot transition was accepted.");
    invalid([](auto &q) { q.transition.startFraction = .1; },
            "A shape transition overlapping the slot was accepted.");
    invalid(
        [](auto &q) {
            q.transition.method = TransitionMethod::lut;
            q.transition.knots = {{0, 0}, {.2, .1}, {1, 1}};
        },
        "A LUT changing shape inside the adapter was accepted.");
    invalid(
        [](auto &q) {
            q.guideEnabled = true;
            q.guide.position = .2;
        },
        "A guide inside the adapter was accepted.");
    invalid(
        [](auto &q) {
            q.mouthMorph.enabled = true;
            q.mouthMorph.startFraction = .2;
        },
        "Mouth morphing overlapping the slot was accepted.");

    p = fixture();
    p.crossSectionExponent = 2;
    p.verticalAngle = p.horizontalAngle;
    p.verticalProfile = p.profile;
    rejects([&] { (void)makeAxisymmetricMeridian(DraftGeometry(p)); },
            "A noncircular slot was accepted as an axisymmetric horn.");
    p.slot.aperture = {.008, .008, 2};
    const auto meridian = makeAxisymmetricMeridian(DraftGeometry(p));
    std::set<std::uint32_t> meridianTags;
    for (const auto &segment : meridian.segments)
        meridianTags.insert(segment.faceId);
    require(meridianTags.contains(6) && meridianTags.contains(7) && meridianTags.contains(8),
            "A rotationally symmetric inlet lost its physical segment identities.");
}

Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double norm(Vec3 a) {
    return std::hypot(a.x, a.y, a.z);
}

struct SourceMeasures {
    double area{}, perimeter{};
};
SourceMeasures sourceMeasures(const SurfaceMesh &mesh, bool withSlot) {
    using Edge = std::pair<std::uint32_t, std::uint32_t>;
    SourceMeasures result;
    std::map<Edge, unsigned> drivenEdges;
    std::set<std::uint32_t> faces;
    for (const auto &triangle : mesh.triangles) {
        faces.insert(triangle.faceId);
        require((triangle.faceId == 5) == (triangle.boundary == BoundaryCondition::driven),
                "The slot introduced another source or changed a wall condition.");
        if (triangle.boundary != BoundaryCondition::driven)
            continue;
        const auto a = mesh.vertices.at(triangle.vertices[0]);
        const auto b = mesh.vertices.at(triangle.vertices[1]);
        const auto c = mesh.vertices.at(triangle.vertices[2]);
        require(a.z == 0 && b.z == 0 && c.z == 0, "A source face moved away from the inlet.");
        result.area += .5 * norm(cross(subtract(b, a), subtract(c, a)));
        for (unsigned i = 0; i < 3; ++i)
            ++drivenEdges[std::minmax(triangle.vertices[i], triangle.vertices[(i + 1) % 3])];
    }
    require(faces == (withSlot ? std::set<std::uint32_t>{1, 2, 3, 4, 5, 6, 7, 8}
                               : std::set<std::uint32_t>{1, 2, 3, 4, 5}),
            "Meshing lost or invented a physical face tag.");
    for (const auto &[edge, count] : drivenEdges)
        if (count == 1)
            result.perimeter +=
                norm(subtract(mesh.vertices[edge.first], mesh.vertices[edge.second]));
    require(result.area > 0 && result.perimeter > 0, "The source patch is absent.");
    return result;
}

MeshingResult checkedMesh(const DraftParameters &p, double edgeLength,
                          double geometricTolerance = .0002) {
    const DraftGeometry geometry(p);
    const auto preview = geometry.previewMesh(geometricTolerance);
    const auto before = sourceMeasures(preview, p.slot.enabled);
    if (p.slot.enabled)
        for (const auto &triangle : preview.triangles)
            if (triangle.faceId == 8)
                for (const auto index : triangle.vertices) {
                    const auto point = preview.vertices.at(index);
                    require(point.z <= geometry.slotStraightEnd() + dimensionTolerance,
                            "Upstream outer-wall face 8 extends beyond the slot assembly.");
                    if (point.z >= 0)
                        near(std::hypot(point.x, point.y),
                             geometry.radius(point.z, std::atan2(point.y, point.x)) +
                                 p.wallThickness,
                             dimensionTolerance, "Face 8 does not follow the assembly outer wall");
                }
    const auto sourceQuality = inspectMesh(preview);
    require(sourceQuality.watertight() && sourceQuality.selfIntersectionPairs == 0,
            "The physical slot preview is open or self-intersecting.");
    MeshingOptions options;
    options.targetEdgeLength = edgeLength;
    options.adaptive.geometricTolerance = geometricTolerance;
    const auto result = remeshIsotropic(preview, options);
    const auto after = sourceMeasures(result.mesh, p.slot.enabled);
    near(after.area, before.area, 1e-10 * before.area, "Remeshing changed the driver area");
    near(after.perimeter, before.perimeter, 1e-10 * before.perimeter,
         "Remeshing changed the physical connection perimeter");
    require(result.quality.watertight() && result.quality.connectedComponents == 1 &&
                result.quality.selfIntersectionPairs == 0,
            "The remeshed slot is not a closed intersection-free manifold.");
    require(result.quality.minimumAngleDegrees >= options.quality.minimumAngleDegrees &&
                result.quality.maximumAspectRatio <= options.quality.maximumAspectRatio,
            "The slot bypassed the enabled triangle-quality limits.");
    require(result.mesh.fixedVertices.size() == preview.fixedVertices.size(),
            "The mesh lost or added fixed wall vertices.");
    for (const auto index : preview.fixedVertices) {
        const auto original = preview.vertices.at(index);
        require(original.z == 0, "An internal slot wall layer was artificially pinned.");
        require(
            std::any_of(result.mesh.fixedVertices.begin(), result.mesh.fixedVertices.end(),
                        [&](auto target) { return result.mesh.vertices.at(target) == original; }),
            "The physical connection changed during remeshing.");
    }
    std::cout << "slot=" << p.slot.enabled << " transition=" << p.slot.transitionLength
              << " contour=" << static_cast<int>(p.slot.contour)
              << " exit=" << static_cast<int>(p.slot.exit)
              << " cornerRadius=" << p.slot.cornerRadius << " exitRadius=" << p.slot.exitRadius
              << " h=" << edgeLength << " geometricTolerance=" << geometricTolerance
              << " triangles=" << result.mesh.triangles.size()
              << " minAngle=" << result.quality.minimumAngleDegrees
              << " maxAspect=" << result.quality.maximumAspectRatio << " drivenArea=" << after.area
              << '\n';
    return result;
}

void meshChecks() {
    auto p = fixture();
    (void)checkedMesh(p, .004);
    (void)checkedMesh(p, .003);
    p.slot.transitionLength = 0;
    (void)checkedMesh(p, .004);
    (void)checkedMesh(doubleRadiusFixture(), .004);
}

void radiiDefaultMeshChecks() {
    // Reproduces the saved 180 mm OS-SE example with unmodified product
    // meshing defaults, including 20 degree / aspect 4 quality constraints.
    DraftParameters p;
    p.horizontalAngle = .5;
    p.verticalAngle = .3;
    p.profile.throatAngle = .1;
    p.verticalProfile.throatAngle = .15;
    p.crossSectionExponent = 4;
    p.slot.enabled = true;
    p.slot.contour = SlotContour::roundedRectangle;
    p.slot.cornerRadius = .002;
    p.slot.exit = SlotExit::circularFillet;
    p.slot.exitRadius = .002;
    p.slot.straightLength = .010;
    p.transition.startFraction = .3;
    p.transition.endFraction = .7;
    const MeshingOptions defaults;
    (void)checkedMesh(p, defaults.targetEdgeLength, defaults.adaptive.geometricTolerance);
}

void defaultMeshChecks() {
    const MeshingOptions defaults;
    std::string failures;
    for (const double width : {DraftParameters{}.slot.aperture.width, .014}) {
        DraftParameters p;
        p.slot.enabled = true;
        p.slot.aperture.width = width;
        // This is the explicit, earliest valid shape-transition start. All other
        // dimensions, quality limits and sizing controls retain product defaults.
        p.transition.startFraction =
            (p.slot.adapterLength + p.slot.straightLength + p.slot.transitionLength) / p.length;
        std::cout << "Default slot mesh: L=" << p.length << " aperture=" << width << 'x'
                  << p.slot.aperture.height << " wall=" << p.wallThickness
                  << " adapter=" << p.slot.adapterLength << " straight=" << p.slot.straightLength
                  << " transition=" << p.slot.transitionLength << " h=" << defaults.targetEdgeLength
                  << " geometricTolerance=" << defaults.adaptive.geometricTolerance
                  << " minAngle=" << defaults.quality.minimumAngleDegrees
                  << " maxAspect=" << defaults.quality.maximumAspectRatio << std::endl;
        try {
            (void)checkedMesh(p, defaults.targetEdgeLength, defaults.adaptive.geometricTolerance);
        } catch (const std::exception &error) {
            // Run both user-facing aperture choices, but retain every failure
            // and fail the command after reporting their individual diagnostics.
            const auto message =
                "Default slot width=" + std::to_string(width) + ": " + error.what();
            std::cerr << message << '\n';
            failures += message + '\n';
        }
    }
    if (!failures.empty())
        throw std::runtime_error(failures);
}

#ifdef HORNSCULPT_SLOT_CUDA_TEST
bool finite(std::complex<double> value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

struct AcousticDifference {
    double pressure{}, impedance{};
};

AcousticDifference difference(const BemResult &actual, const BemResult &reference,
                              const BemRequest &request) {
    require(actual.frequencies.size() == request.frequencies.size() &&
                reference.frequencies.size() == request.frequencies.size() &&
                !actual.device.empty() && !reference.device.empty() && actual.peakDeviceBytes > 0 &&
                reference.peakDeviceBytes > 0,
            "A convergence comparison lacks actual GPU frequency results.");
    AcousticDifference worst;
    for (std::size_t i = 0; i < request.frequencies.size(); ++i) {
        const auto &a = actual.frequencies[i].excitations.at(0);
        const auto &b = reference.frequencies[i].excitations.at(0);
        require(a.boundaryResidual <= request.numerics.relativeResidualTolerance &&
                    b.boundaryResidual <= request.numerics.relativeResidualTolerance,
                "A convergence solve failed the separate algebraic residual budget.");
        require(a.observerPressure.size() == request.observers.size() &&
                    b.observerPressure.size() == request.observers.size(),
                "A convergence comparison has incomplete field observations.");
        double errorSquared = 0, referenceSquared = 0;
        for (std::size_t j = 0; j < request.observers.size(); ++j) {
            require(finite(a.observerPressure[j]) && finite(b.observerPressure[j]),
                    "A convergence pressure is not finite.");
            errorSquared += std::norm(a.observerPressure[j] - b.observerPressure[j]);
            referenceSquared += std::norm(b.observerPressure[j]);
        }
        require(referenceSquared > 1e-20,
                "The acoustic comparison has no nonzero reference field.");
        const double pressure = std::sqrt(errorSquared / referenceSquared);
        const auto za = a.impedances.at(0).volumeImpedance;
        const auto zb = b.impedances.at(0).volumeImpedance;
        require(finite(za) && finite(zb) && std::abs(zb) > 1e-12,
                "The impedance comparison has no finite nonzero reference.");
        const double impedance = std::abs(za - zb) / std::abs(zb);
        std::cout << "convergence f=" << request.frequencies[i]
                  << " complexFieldRelativeL2=" << pressure
                  << " complexImpedanceRelative=" << impedance << '\n';
        worst.pressure = std::max(worst.pressure, pressure);
        worst.impedance = std::max(worst.impedance, impedance);
    }
    return worst;
}

void convergenceChecks(const DraftParameters &p, CudaBemSolver &mediumSolver,
                       const BemRequest &request, const BemResult &medium) {
    auto unit = request;
    unit.excitations.resize(1);
    const auto coarseMesh = checkedMesh(p, .006);
    const auto fineMesh = checkedMesh(p, .003);
    require(coarseMesh.mesh.triangles.size() < medium.triangleCount &&
                medium.triangleCount < fineMesh.mesh.triangles.size(),
            "The convergence study did not actually refine the acoustic mesh.");
    CudaBemSolver coarseSolver(coarseMesh.mesh);
    const auto coarse = coarseSolver.solve(unit);
    CudaBemSolver fineSolver(fineMesh.mesh);
    const auto fine = fineSolver.solve(unit);
    std::cout << "Mesh convergence 6/4/3 mm; differences are complex, including phase.\n";
    const auto coarseDifference = difference(coarse, fine, unit);
    const auto mediumDifference = difference(medium, fine, unit);
    // At 4 kHz even the coarse nominal edge is under lambda/14. These are
    // bounded regression budgets for this fixture, not a general BEM accuracy
    // claim. Require actual refinement improvement unless already below 1%.
    constexpr double meshDifferenceBudget = .08;
    constexpr double alreadyConvergedBudget = .01;
    require(mediumDifference.pressure < meshDifferenceBudget &&
                mediumDifference.impedance < meshDifferenceBudget,
            "The slot field or impedance did not settle within the mesh-refinement budget.");
    require(mediumDifference.pressure <=
                    std::max(alreadyConvergedBudget, coarseDifference.pressure) &&
                mediumDifference.impedance <=
                    std::max(alreadyConvergedBudget, coarseDifference.impedance),
            "Refining the slot mesh did not improve the finest-mesh comparison.");

    auto strictNear = unit;
    strictNear.numerics.nearRelativeTolerance = 1e-6;
    strictNear.numerics.nearMaximumDepth = 8;
    const auto strict = mediumSolver.solve(strictNear);
    std::cout << "Independent near-quadrature convergence on the unchanged 4 mm mesh: "
                 "tolerance/depth 1e-4/6 versus 1e-6/8.\n";
    const auto nearDifference = difference(medium, strict, unit);
    constexpr double nearQuadratureBudget = 1e-3;
    require(nearDifference.pressure < nearQuadratureBudget &&
                nearDifference.impedance < nearQuadratureBudget,
            "The slot field or impedance is not stable under stricter near quadrature.");
    std::cout << "meshRatioPressure="
              << mediumDifference.pressure / std::max(1e-15, coarseDifference.pressure)
              << " meshRatioImpedance="
              << mediumDifference.impedance / std::max(1e-15, coarseDifference.impedance)
              << " nearPressure=" << nearDifference.pressure
              << " nearImpedance=" << nearDifference.impedance << '\n';
}

void bemChecks() {
    for (const int variant : {0, 1, 2, 3}) {
        auto p = variant == 3 ? doubleRadiusFixture() : fixture();
        p.slot.enabled = variant != 0;
        if (variant == 1)
            p.slot.transitionLength = 0;
        const auto mesh = checkedMesh(p, .004);
        BemRequest request;
        request.frequencies = {1000, 4000};
        request.excitations = {{"unit", {{5, {1, 0}}}}, {"double", {{5, {2, 0}}}}};
        request.observers = {
            {0, 0, p.length + 1}, {.707, 0, p.length + .707}, {0, .707, p.length + .707}};
        request.farField = FarFieldObservation{{0, 0, p.length}, 1};
        CudaBemSolver solver(mesh.mesh);
        const auto result = solver.solve(request);
        require(!result.device.empty() && result.peakDeviceBytes > 0 &&
                    result.triangleCount == mesh.mesh.triangles.size(),
                "The slot acoustic test did not run a real GPU solve.");
        require(result.frequencies.size() == request.frequencies.size(),
                "An acoustic frequency result is absent.");
        for (const auto &frequency : result.frequencies) {
            require(frequency.excitations.size() == 2, "The linearity excitation is absent.");
            const auto &unit = frequency.excitations[0];
            const auto &twice = frequency.excitations[1];
            require(unit.boundaryResidual <= request.numerics.relativeResidualTolerance &&
                        twice.boundaryResidual <= request.numerics.relativeResidualTolerance,
                    "The slot BEM solve failed its original-operator residual tolerance.");
            require(unit.observerPressure.size() == request.observers.size() &&
                        twice.observerPressure.size() == request.observers.size(),
                    "The acoustic observation set is incomplete.");
            for (std::size_t i = 0; i < request.observers.size(); ++i) {
                require(finite(unit.observerPressure[i]) && finite(twice.observerPressure[i]),
                        "A slot acoustic pressure is not finite.");
                require(std::abs(twice.observerPressure[i] - 2.0 * unit.observerPressure[i]) <=
                            1e-4 * std::max(1e-12, std::abs(2.0 * unit.observerPressure[i])),
                        "The slot acoustic result violates source linearity.");
            }
            require(unit.impedances.size() == 1 && unit.impedances.front().faceId == 5 &&
                        finite(unit.impedances.front().volumeImpedance),
                    "The driver impedance was not calculated on the original source.");
            const auto &impedance = unit.impedances.front();
            near(impedance.area, sourceMeasures(mesh.mesh, p.slot.enabled).area, 1e-10,
                 "The acoustic source area differs from the geometric driver area");
            std::cout << "GPU=" << result.device << " variant=" << variant
                      << " f=" << frequency.frequency << " residual=" << unit.boundaryResidual
                      << " pAxis=" << unit.observerPressure[0]
                      << " pH45=" << unit.observerPressure[1]
                      << " pV45=" << unit.observerPressure[2]
                      << " Zin=" << impedance.volumeImpedance << '\n';
        }
        if (variant == 2)
            convergenceChecks(p, solver, request, result);
    }
}
#endif
} // namespace

int main(int argc, char **argv) {
    try {
        if (argc == 1) {
            elongatedSlotLargeFillet();
            trueRadii();
            disabledEquivalence();
            dimensionsAndSegments();
            tangentsAndAnalysis();
            zeroLengthAndCovariance();
            compositionAndRejections();
            std::cout
                << "Diffraction-slot geometry, derivatives, analysis and constraints passed.\n";
        } else if (argc == 2 && std::string_view(argv[1]) == "--elongated-fillet-only") {
            elongatedSlotLargeFillet();
        } else if (argc == 2 && std::string_view(argv[1]) == "--mesh") {
            meshChecks();
        } else if (argc == 2 && std::string_view(argv[1]) == "--default-mesh") {
            defaultMeshChecks();
        } else if (argc == 2 && std::string_view(argv[1]) == "--radii-default-mesh") {
            radiiDefaultMeshChecks();
        } else if (argc == 2 && std::string_view(argv[1]) == "--bem") {
#ifdef HORNSCULPT_SLOT_CUDA_TEST
            bemChecks();
#else
            throw std::runtime_error("This test was built without real CUDA BEM support.");
#endif
        } else {
            throw std::invalid_argument(
                "Usage: diffraction_slot_tests [--elongated-fillet-only|--mesh|--default-mesh|"
                "--radii-default-mesh|--bem]");
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
