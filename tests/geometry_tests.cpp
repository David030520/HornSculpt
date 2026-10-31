#include <hornsculpt/draft_geometry.hpp>
#include <hornsculpt/geometry_analysis.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

using namespace hornsculpt;

namespace {
void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const std::string &message) {
    if (std::abs(actual - expected) > tolerance || !std::isfinite(actual))
        throw std::runtime_error(message + ": expected " + std::to_string(expected) + ", got " +
                                 std::to_string(actual));
}

template <class Function> void rejects(Function function, const std::string &message) {
    bool rejected = false;
    try {
        function();
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, message);
}

void profileReferences() {
    OsseParameters p;
    p.terminationScale = 0.0;
    const OsseProfile os(p);
    const double angle = 0.61, z = 0.12;
    near(os.radius(0.0, angle), p.throatRadius, 1e-14, "OS throat");
    near(os.radius(z, angle), std::hypot(p.throatRadius, z * std::tan(angle)), 1e-14,
         "OS analytical radius");
    near(os.angleThrough(z, os.radius(z, angle)), angle, 1e-14, "OS inverse");
    p.k = 0.0;
    const OsseProfile cone(p);
    near(cone.radius(z, angle), p.throatRadius + z * std::tan(angle), 1e-14, "cone radius");
    near(cone.slope(0.0, angle), std::tan(angle), 1e-14, "cone throat slope");
    p.terminationScale = 0.3;
    p.terminationQ = 0.9;
    p.terminationExponent = 4.0;
    const OsseProfile terminated(p);
    const double term = p.terminationScale * p.length / p.terminationQ *
                        (1.0 - std::pow(1.0 - std::pow(p.terminationQ, 4.0), 0.25));
    near(terminated.radius(p.length, angle), p.throatRadius + p.length * std::tan(angle) + term,
         1e-14, "SE termination reference");
    const double step = 1e-7;
    near(terminated.slope(z, angle),
         (terminated.radius(z + step, angle) - terminated.radius(z - step, angle)) / (2.0 * step),
         1e-8, "SE analytical derivative");
    rejects([&] { (void)terminated.angleThrough(z, 1e-8); }, "infeasible guide rejected");
    p.terminationQ = 1.0;
    require(std::isinf(OsseProfile(p).slope(p.length, angle)), "q=1 vertical mouth tangent");
}

void curveReferences() {
    Superellipse ellipse(0.2, 0.1, 2.0);
    near(ellipse.area(), .02 * std::numbers::pi, 1e-14, "ellipse enclosed area");
    near(Superellipse(.2, .1, 1e9).area(), .08, 1e-9, "superellipse rectangular limit area");
    near(ellipse.radius(0.0), 0.2, 1e-14, "SE horizontal axis");
    near(ellipse.radius(std::numbers::pi / 2.0), 0.1, 1e-14, "SE vertical axis");
    near(Superellipse(1.0, 1.0, 4.0).radius(std::numbers::pi / 4.0), std::pow(2.0, 0.25), 1e-14,
         "SE diagonal reference");
    for (double exponent : {2.0, 4.0, 6.0, 2.75}) {
        const double diagonal = std::sqrt(2.0) * std::pow(2.0, -1.0 / exponent);
        for (double scale : {1e-100, 1.0, 1e100}) {
            const auto actual = Superellipse(scale, scale, exponent).radius(std::numbers::pi / 4.0);
            near(actual / scale, diagonal, 1e-14, "scaled superellipse diagonal radius");
        }
    }
    rejects([&] { (void)ellipse.radius(std::numeric_limits<double>::infinity()); },
            "nonfinite superellipse angle rejected");
    SuperformulaParameters p;
    const Superformula circle(p);
    for (int i = 0; i < 100; ++i)
        near(circle.radius(0.13 * i), 1.0, 1e-14, "SF circle specialization");
    p.m1 = 3.0;
    p.m2 = 4.0;
    p.n2 = 4.0;
    rejects([&] { (void)Superformula(p); }, "nonperiodic SF rejected");
    p = {};
    p.n1 = 0.0;
    rejects([&] { (void)Superformula(p); }, "invalid SF exponent rejected");

    // ATH 4.8.2, section 6.4.2, public six-parameter example. These reference
    // radii are before the independent physical width/aspect transformation.
    p = {.95, .82, 4, 4, .8, 7, 1.9};
    const Superformula ath(p);
    const std::array angles{0., std::numbers::pi / 6, std::numbers::pi / 4, std::numbers::pi / 2};
    const std::array radii{.6383833424176973, 1.1192185405034356, 1.1711514318004939,
                           .6241771981884541};
    for (std::size_t i = 0; i < angles.size(); ++i)
        near(ath.radius(angles[i]), radii[i], 1e-13, "ATH public superformula example");
    p.n1 = -.8;
    const Superformula negative(p);
    for (const auto phi : angles)
        near(negative.radius(phi) * ath.radius(phi), 1., 1e-13,
             "nonzero negative n1 has a finite reciprocal contour");
    p = {.95, .82, 8, 4, .8, 7, 1.9};
    const Superformula unequal(p);
    for (unsigned i = 0; i < 100; ++i)
        near(unequal.radius(i * .071), unequal.radius(i * .071 + 2 * std::numbers::pi), 1e-12,
             "independent periodic superformula frequencies");
    p = {.8, 2.7, 4.8, 1.9, 2.9, .3, 3.6};
    rejects([&] { (void)Superformula(p); }, "reported nonperiodic morph contour rejected");
    p.n2 = -1;
    rejects([&] { (void)Superformula(p); }, "nonpositive secondary exponent rejected");
}

void guideControlSensitivity() {
    DraftParameters p;
    p.guideEnabled = true;
    p.guide.type = CurveType::superformula;
    p.guide.width = .6;
    p.guide.position = .55;
    p.transition.startFraction = .1;
    p.transition.endFraction = .9;
    const DraftGeometry original(p);
    require(std::abs(original.guideRadius(.4) - original.guideRadius(0)) > .001,
            "New superformula contours start in a circular control degeneracy.");
    const auto checkHit = [](const DraftGeometry &g) {
        for (unsigned i = 0; i < 120; ++i) {
            const double phi = i * 2 * std::numbers::pi / 120;
            near(g.radius(g.guidePosition(), phi), g.guideRadius(phi), 1e-12,
                 "composed profile meets the entire guide, including off-axis directions");
        }
    };
    checkHit(original);
    double SuperformulaParameters::*members[] = {
        &SuperformulaParameters::a,  &SuperformulaParameters::b,  &SuperformulaParameters::m1,
        &SuperformulaParameters::m2, &SuperformulaParameters::n1, &SuperformulaParameters::n2,
        &SuperformulaParameters::n3};
    const double values[]{1.1, .95, 8, 2, 1.3, 6, 3};
    for (unsigned member = 0; member < std::size(members); ++member) {
        auto changed = p;
        changed.guide.superformula.*members[member] = values[member];
        const DraftGeometry variant(changed);
        checkHit(variant);
        double difference = 0;
        for (unsigned i = 0; i < 64; ++i)
            difference = std::max(difference, std::abs(variant.guideRadius(i * .093) -
                                                       original.guideRadius(i * .093)));
        require(difference > 1e-5, "A nondegenerate guide control has no geometric effect.");
    }
    auto rotated = p;
    rotated.guide.rotation = .31;
    const DraftGeometry rotation(rotated);
    for (unsigned i = 0; i < 20; ++i)
        near(rotation.guideRadius(.31 + i * .1), original.guideRadius(i * .1), 1e-12,
             "guide rotation is independent of physical scaling");
    rotated.guide.positionMode = GuidePositionMode::absolute;
    rotated.guide.position = .55 * rotated.length;
    rotated.transition.method = TransitionMethod::lut;
    rotated.transition.knots = {{0, 0}, {.3, .1}, {.55, .6}, {1, 1}};
    rotated.mouthMorph.enabled = true;
    rotated.mouthMorph.startFraction = .7;
    checkHit(DraftGeometry(rotated));
    rotated.guide.superformula.n1 = -.8;
    checkHit(DraftGeometry(rotated));

    // A narrow vertical guide can be valid even when its radius is below the
    // horizontal entry half-width. The inactive horizontal direction must not
    // impose its full throat radius on every polar ray.
    p = {};
    p.entry = {.032, .008, 8};
    p.profile.terminationScale = p.verticalProfile.terminationScale = 0;
    p.guideEnabled = true;
    p.guide.width = .06;
    p.guide.aspectRatio = .25;
    const DraftGeometry rectangular(p);
    checkHit(rectangular);
    const double z = .7 * rectangular.length(), step = 1e-7;
    near(rectangular.radiusDerivative(z, 1.1),
         (rectangular.radius(z + step, 1.1) - rectangular.radius(z - step, 1.1)) / (2 * step), 1e-8,
         "unequal-entry guide retains its analytic axial derivative");
    p.guide.width = .01;
    rejects([&] { (void)DraftGeometry(p); }, "genuinely undersized guide is still rejected");
}

void transitions() {
    ShapeTransitionParameters p;
    p.startFraction = 0.2;
    p.endFraction = 0.8;
    const ShapeTransition range(p);
    near(range.blend(0.1), 0.0, 0.0, "transition lower plateau");
    near(range.blend(0.5), 0.5, 1e-14, "transition midpoint");
    near(range.blend(0.9), 1.0, 0.0, "transition upper plateau");
    near(range.derivative(0.2), 0.0, 0.0, "transition smooth start");
    near(range.derivative(0.8), 0.0, 0.0, "transition smooth end");
    p.method = TransitionMethod::lut;
    p.knots = {{0.0, 0.0}, {0.25, 0.0}, {0.6, 0.7}, {0.8, 1.0}, {1.0, 1.0}};
    const ShapeTransition lut(p);
    near(lut.blend(0.6), 0.7, 1e-14, "LUT interpolation knot");
    double previous = 0.0;
    for (int i = 0; i <= 1000; ++i) {
        const double value = lut.blend(i / 1000.0);
        require(value >= previous - 1e-14 && value <= 1.0 + 1e-14, "LUT no overshoot");
        previous = value;
    }
    p.knots[2].fraction = p.knots[1].fraction;
    rejects([&] { (void)ShapeTransition(p); }, "duplicate LUT positions rejected");
}

void guidesAndComposition() {
    DraftParameters p;
    p.guideEnabled = true;
    p.guide.width = 0.2;
    p.guide.aspectRatio = 0.8;
    p.guide.position = 0.55;
    p.transition.startFraction = 0.1;
    p.transition.endFraction = 0.9;
    DraftGeometry geometry(p);
    near(geometry.guideRadius(0.0), 0.1, 1e-12, "physical guide width");
    near(geometry.guideRadius(std::numbers::pi / 2.0), 0.08, 1e-12, "physical guide height");
    for (int i = 0; i < 360; ++i) {
        const double phi = std::numbers::pi * i / 180.0;
        near(geometry.radius(geometry.guidePosition(), phi), geometry.guideRadius(phi), 1e-12,
             "guide exact inside form transition");
        near(geometry.radius(0.0, phi), .5 * p.entry.width, 1e-12, "circular throat");
    }
    auto changed = p;
    changed.horizontalAngle *= .7;
    const DraftGeometry changedBase(changed);
    require(std::abs(changedBase.radius(p.length * .25, 0) -
                     geometry.radius(p.length * .25, 0)) > 1e-6,
            "The base angle must remain effective during the guide transition");
    near(changedBase.radius(changedBase.guidePosition(), 0), geometry.guideRadius(0),
         1e-12, "base-angle editing retains the exact guide hit");
    changed = p;
    changed.guide.superellipseExponent = 6;
    const DraftGeometry changedCorners(changed);
    near(changedCorners.guideRadius(0), geometry.guideRadius(0), 1e-12,
         "guide exponent retains horizontal dimension");
    near(changedCorners.guideRadius(std::numbers::pi / 2), geometry.guideRadius(std::numbers::pi / 2),
         1e-12, "guide exponent retains vertical dimension");
    require(std::abs(changedCorners.guideRadius(std::numbers::pi / 4) -
                     geometry.guideRadius(std::numbers::pi / 4)) > 1e-4,
            "Guide exponent must change the corners");
    const double z = geometry.length() * 0.43, h = geometry.length() * 1e-6;
    near(geometry.radiusDerivative(z, 0.7),
         (geometry.radius(z + h, 0.7) - geometry.radius(z - h, 0.7)) / (2.0 * h), 1e-7,
         "composed derivative");
    p.guide.rotation = 0.37;
    const DraftGeometry rotated(p);
    near(rotated.guideRadius(0.37), geometry.guideRadius(0.0), 1e-12, "guide rotation");
    p.guide.positionMode = GuidePositionMode::absolute;
    p.guide.position = geometry.guidePosition();
    near(DraftGeometry(p).guidePosition(), geometry.guidePosition(), 1e-14,
         "absolute guide position");
    p.guide.type = CurveType::superformula;
    p.guide.superformula.n1 = 3.0;
    p.guide.superformula.n2 = p.guide.superformula.n3 = 4.0;
    p.transition.method = TransitionMethod::lut;
    p.transition.knots = {{0.0, 0.0}, {0.3, 0.1}, {0.55, 0.6}, {1.0, 1.0}};
    const DraftGeometry sf(p);
    for (int i = 0; i < 100; ++i) {
        const double phi = i * 0.071;
        near(sf.radius(sf.guidePosition(), phi), sf.guideRadius(phi), 1e-12, "SF guide with LUT");
    }
    p.guide.positionMode = GuidePositionMode::relative;
    p.guide.position = 0.02;
    p.transition.method = TransitionMethod::smoothRange;
    p.transition.startFraction = 0.2;
    rejects([&] { (void)DraftGeometry(p); }, "unreachable pre-transition guide rejected");
}

class IndependentCone final : public LongitudinalProfile {
  public:
    explicit IndependentCone(double throatRadius = .01) : throatRadius_(throatRadius) {}
    double radius(double z, double angle) const override {
        return throatRadius_ + z * std::tan(angle);
    }
    double angleThrough(double z, double r) const override {
        return std::atan((r - throatRadius_) / z);
    }
    double slope(double, double angle) const override {
        return std::tan(angle);
    }
    double length() const noexcept override {
        return 0.23;
    }

  private:
    double throatRadius_;
};

void analysisReferences() {
    DraftParameters p;
    p.length = .23;
    p.entry.width = p.entry.height = .02;
    p.horizontalAngle = p.verticalAngle = 0.4;
    const DraftGeometry geometry(p, std::make_shared<IndependentCone>());
    const auto analysis = analyzeGeometry(geometry, 31);
    const double tangent = std::tan(0.4), r0 = 0.01;
    near(geometry.length(), 0.23, 0.0, "replaceable profile length");
    for (const auto &sample : analysis.samples) {
        const double z = sample.z, r = r0 + z * tangent;
        near(sample.area, std::numbers::pi * r * r, 1e-12, "cone area");
        near(sample.areaDerivative, 2.0 * std::numbers::pi * r * tangent, 1e-11, "cone dA/dz");
        near(sample.volume,
             std::numbers::pi *
                 (r0 * r0 * z + r0 * tangent * z * z + tangent * tangent * z * z * z / 3.0),
             1e-12, "cone exact volume");
        near(sample.horizontalWallAngle, 0.4, 1e-14, "cone horizontal wall angle");
        near(sample.verticalWallAngle, 0.4, 1e-14, "cone vertical wall angle");
    }
    p.profile.terminationQ = 1.0;
    p.verticalProfile.terminationQ = 1.0;
    const auto verticalMouth = analyzeGeometry(DraftGeometry(p), 11).samples.back();
    require(std::isinf(verticalMouth.areaDerivative), "vertical endpoint area derivative explicit");
    near(verticalMouth.horizontalWallAngle, std::numbers::pi / 2.0, 1e-14,
         "vertical endpoint wall angle");
}

void directionalProfiles() {
    DraftParameters p;
    p.profile.throatAngle = .08;
    p.profile.k = 1.1;
    p.profile.terminationScale = .16;
    p.profile.terminationQ = .95;
    p.profile.terminationExponent = 4.2;
    p.verticalProfile.throatAngle = .19;
    p.verticalProfile.k = .7;
    p.verticalProfile.terminationScale = .24;
    p.verticalProfile.terminationQ = .9;
    p.verticalProfile.terminationExponent = 3.2;
    p.crossSectionExponent = 4.0;
    p.transition.startFraction = .1;
    p.transition.endFraction = .9;
    const DraftGeometry original(p);
    OsseParameters horizontalParameters;
    static_cast<OsseDirectionalParameters &>(horizontalParameters) = p.profile;
    horizontalParameters.throatRadius = .5 * p.entry.width;
    horizontalParameters.length = p.length;
    auto verticalParameters = horizontalParameters;
    verticalParameters.throatRadius = .5 * p.entry.height;
    static_cast<OsseDirectionalParameters &>(verticalParameters) = p.verticalProfile;
    const OsseProfile h(horizontalParameters), v(verticalParameters);
    for (double fraction : {0.0, .03, .35, .8, 1.0}) {
        const double z = original.length() * fraction;
        near(original.radius(z, 0.0), h.radius(z, p.horizontalAngle), 1e-13,
             "horizontal profile survives the complete shape transition");
        near(original.radius(z, std::numbers::pi / 2.0), v.radius(z, p.verticalAngle), 1e-13,
             "vertical profile survives the complete shape transition");
        near(original.radiusDerivative(z, 0.0), h.slope(z, p.horizontalAngle), 1e-12,
             "horizontal profile derivative");
        near(original.radiusDerivative(z, std::numbers::pi / 2.0), v.slope(z, p.verticalAngle),
             1e-12, "vertical profile derivative");
    }
    double OsseDirectionalParameters::*members[]{
        &OsseDirectionalParameters::throatAngle, &OsseDirectionalParameters::k,
        &OsseDirectionalParameters::terminationScale, &OsseDirectionalParameters::terminationQ,
        &OsseDirectionalParameters::terminationExponent};
    const double values[]{.27, 1.5, .35, .82, 6.0};
    for (unsigned i = 0; i < std::size(members); ++i) {
        auto changedH = p, changedV = p;
        // All five controls change their own axis without altering the other.
        const auto member = members[i];
        changedH.profile.*member = values[i];
        changedV.verticalProfile.*member = values[i];
        const DraftGeometry horizontalChange(changedH), verticalChange(changedV);
        const double z = .7 * original.length();
        near(horizontalChange.radius(z, std::numbers::pi / 2.0),
             original.radius(z, std::numbers::pi / 2.0), 1e-13,
             "horizontal parameter does not change vertical geometry");
        near(verticalChange.radius(z, 0.0), original.radius(z, 0.0), 1e-13,
             "vertical parameter does not change horizontal geometry");
        require(std::abs(horizontalChange.radius(z, 0.0) - original.radius(z, 0.0)) > 1e-7,
                "horizontal parameter is active");
        require(std::abs(verticalChange.radius(z, std::numbers::pi / 2.0) -
                         original.radius(z, std::numbers::pi / 2.0)) > 1e-7,
                "vertical parameter is active");
    }
    p.transition.method = TransitionMethod::lut;
    p.transition.knots = {{0.0, 0.0}, {.4, .0}, {.8, 1.0}, {1.0, 1.0}};
    const DraftGeometry lut(p);
    near(std::atan(lut.radiusDerivative(0.0, 0.0)), p.profile.throatAngle, 1e-14,
         "LUT preserves horizontal throat angle");
    near(std::atan(lut.radiusDerivative(0.0, std::numbers::pi / 2.0)),
         p.verticalProfile.throatAngle, 1e-14, "LUT preserves vertical throat angle");
    p.guideEnabled = true;
    p.guide.width = .24;
    p.guide.aspectRatio = .8;
    p.guide.position = .8;
    p.guide.rotation = .23;
    const DraftGeometry guide(p);
    for (unsigned i = 0; i < 100; ++i) {
        const double phi = i * .071;
        near(guide.radius(guide.guidePosition(), phi), guide.guideRadius(phi), 1e-12,
             "independent profiles retain rotated guide constraints");
        const double z = .53 * guide.length(), step = 1e-7 * guide.length();
        near(guide.radiusDerivative(z, phi),
             (guide.radius(z + step, phi) - guide.radius(z - step, phi)) / (2.0 * step), 2e-7,
             "independent guide profile derivative");
    }
    p.guideEnabled = false;
    p.profile.terminationQ = 1.0;
    const DraftGeometry verticalEndpoint(p);
    require(std::isinf(verticalEndpoint.radiusDerivative(verticalEndpoint.length(), 0.0)) &&
                std::isfinite(verticalEndpoint.radiusDerivative(verticalEndpoint.length(),
                                                                std::numbers::pi / 2.0)),
            "a vertical endpoint on one axis does not contaminate the other axis");
}

void closedPhysicalShell(const DraftGeometry &geometry = DraftGeometry{}) {
    const auto mesh = geometry.previewMesh();
    struct EdgeUse {
        unsigned count{};
        int orientation{};
    };
    std::map<std::pair<std::uint32_t, std::uint32_t>, EdgeUse> edges;
    std::map<std::uint32_t, unsigned> faces;
    double signedVolume = 0.0;
    for (const auto &triangle : mesh.triangles) {
        ++faces[triangle.faceId];
        require((triangle.boundary == BoundaryCondition::driven) == (triangle.faceId == 5),
                "driven tag");
        for (unsigned i = 0; i < 3; ++i) {
            const auto a = triangle.vertices[i], b = triangle.vertices[(i + 1) % 3];
            auto &use = edges[std::minmax(a, b)];
            ++use.count;
            use.orientation += a < b ? 1 : -1;
        }
        const auto &a = mesh.vertices[triangle.vertices[0]];
        const auto &b = mesh.vertices[triangle.vertices[1]];
        const auto &c = mesh.vertices[triangle.vertices[2]];
        const Vec3 ab{b.x - a.x, b.y - a.y, b.z - a.z};
        const Vec3 ac{c.x - a.x, c.y - a.y, c.z - a.z};
        const Vec3 normal{ab.y * ac.z - ab.z * ac.y, ab.z * ac.x - ab.x * ac.z,
                          ab.x * ac.y - ab.y * ac.x};
        require(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z > 1e-24,
                "nondegenerate shell triangle");
        signedVolume += (a.x * (b.y * c.z - b.z * c.y) + a.y * (b.z * c.x - b.x * c.z) +
                         a.z * (b.x * c.y - b.y * c.x)) /
                        6.0;
    }
    for (const auto &[key, use] : edges) {
        (void)key;
        require(use.count == 2 && use.orientation == 0, "closed consistently oriented shell edges");
    }
    require(faces.size() == 5, "five physical surface patches");
    require(signedVolume > 0.0, "outward shell orientation");
}

void asymmetricProfiles() {
    DraftParameters p;
    p.asymmetricEnabled = true;
    p.entry = {.032, .024, 3.0};
    p.leftAngle = .46;
    p.bottomAngle = .28;
    p.profile = {.04, 1.1, .1, .95, 4.0};
    p.verticalProfile = {.08, .9, .15, .94, 3.5};
    p.leftProfile = {.13, .7, .18, .91, 4.5};
    p.bottomProfile = {.18, 1.3, .13, .93, 3.2};
    p.crossSectionExponent = 4.0;
    const DraftGeometry geometry(p);
    const OsseDirectionalParameters profiles[]{p.profile, p.verticalProfile, p.leftProfile,
                                               p.bottomProfile};
    const double angles[]{p.horizontalAngle, p.verticalAngle, p.leftAngle, p.bottomAngle};
    for (unsigned i = 0; i < 4; ++i) {
        OsseParameters referenceParameters;
        static_cast<OsseDirectionalParameters &>(referenceParameters) = profiles[i];
        referenceParameters.length = p.length;
        referenceParameters.throatRadius = .5 * (i % 2 == 0 ? p.entry.width : p.entry.height);
        const OsseProfile reference(referenceParameters);
        const double phi = std::numbers::pi * i / 2;
        for (double fraction : {0., .2, .61, 1.}) {
            const double z = fraction * p.length;
            near(geometry.radius(z, phi), reference.radius(z, angles[i]), 1e-13,
                 "each cardinal profile preserves its independent OS-SE radius");
            near(geometry.radiusDerivative(z, phi), reference.slope(z, angles[i]), 1e-12,
                 "each cardinal profile preserves its independent OS-SE slope");
        }
    }
    const auto checkContinuity = [](const DraftGeometry &g) {
        const double z = .57 * g.length(), angularStep = 1e-6, axialStep = 1e-7;
        for (int i = -4; i <= 4; ++i) {
            const double phi = std::numbers::pi * i / 2;
            const double center = g.radius(z, phi);
            near(g.radius(z, phi + 2 * std::numbers::pi), center, 1e-13,
                 "quadrant composition is periodic, including negative azimuths");
            near((g.radius(z, phi + angularStep) - center) / angularStep,
                 (center - g.radius(z, phi - angularStep)) / angularStep, 3e-5,
                 "quadrant transitions retain their angular tangent");
        }
        for (unsigned i = 0; i < 80; ++i) {
            const double phi = i * .1 - 4;
            near(g.radiusDerivative(z, phi),
                 (g.radius(z + axialStep, phi) - g.radius(z - axialStep, phi)) / (2 * axialStep),
                 2e-7, "four-sided composition retains the analytic axial derivative");
        }
    };
    checkContinuity(geometry);
    closedPhysicalShell(geometry);
    auto changed = p;
    changed.leftAngle = .3;
    changed.leftProfile.terminationScale = .27;
    changed.bottomAngle = .65;
    changed.bottomProfile.k = .6;
    const DraftGeometry oppositeChanged(changed);
    for (unsigned i = 0; i <= 16; ++i) {
        const double phi = std::numbers::pi * i / 32;
        near(oppositeChanged.radius(.11, phi), geometry.radius(.11, phi), 1e-14,
             "left and bottom controls leave the right/top quadrant unchanged");
    }
    p.guideEnabled = true;
    p.guide.width = .3;
    p.guide.aspectRatio = .8;
    p.guide.position = .65;
    p.guide.rotation = .23;
    p.transition.method = TransitionMethod::lut;
    p.transition.knots = {{0., 0.}, {.2, .05}, {.65, .8}, {1., 1.}};
    p.mouthMorph.enabled = true;
    p.mouthMorph.startFraction = .75;
    p.mouthMorph.endFraction = 1.;
    p.mouthMorph.target.width = .36;
    p.mouthMorph.target.aspectRatio = .8;
    const DraftGeometry transformed(p);
    const PhysicalContour mouth(p.mouthMorph.target);
    checkContinuity(transformed);
    for (unsigned i = 0; i < 80; ++i) {
        const double phi = .1 * i - 4;
        near(transformed.radius(transformed.guidePosition(), phi), transformed.guideRadius(phi),
             1e-12, "all four sides hit the rotated guide through a LUT transition");
        near(transformed.radius(transformed.length(), phi), mouth.radius(phi), 1e-12,
             "mouth morph retains its independent target for all four sides");
    }
    auto invalid = p;
    invalid.leftAngle = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { (void)DraftGeometry(invalid); }, "active invalid left angle rejected");
    invalid = p;
    invalid.bottomProfile.terminationQ = 1.1;
    rejects([&] { (void)DraftGeometry(invalid); }, "active invalid bottom profile rejected");
    invalid.asymmetricEnabled = false;
    invalid.leftAngle = std::numeric_limits<double>::quiet_NaN();
    const DraftGeometry inactive(invalid);
    invalid.leftProfile = invalid.profile;
    invalid.bottomProfile = invalid.verticalProfile;
    invalid.leftAngle = invalid.horizontalAngle;
    invalid.bottomAngle = invalid.verticalAngle;
    invalid.asymmetricEnabled = true;
    const DraftGeometry mirrored(invalid);
    for (unsigned i = 0; i < 80; ++i)
        near(inactive.radius(.1, .1 * i), mirrored.radius(.1, .1 * i), 1e-14,
             "symmetric mode ignores dormant values and matches four equal opposite sides");
}

void asymmetricAnalysis() {
    DraftParameters p;
    p.asymmetricEnabled = true;
    p.entry.width = p.entry.height = .02;
    p.length = .23;
    p.horizontalAngle = .2;
    p.verticalAngle = .3;
    p.leftAngle = .4;
    p.bottomAngle = .5;
    const auto cone = std::make_shared<IndependentCone>();
    const DraftGeometry geometry(p, cone, cone, cone, cone);
    const double h = std::tan(p.horizontalAngle) + std::tan(p.leftAngle);
    const double v = std::tan(p.verticalAngle) + std::tan(p.bottomAngle);
    const double diameter = p.entry.width, quarterPi = std::numbers::pi / 4;
    for (const auto &sample : analyzeGeometry(geometry, 13).samples) {
        const double z = sample.z;
        // Four quarter ellipses give A = pi/4 * (right+left) * (top+bottom).
        near(sample.area, quarterPi * (diameter + h * z) * (diameter + v * z), 1e-10,
             "four-sided cross-sectional area");
        near(sample.areaDerivative, quarterPi * (diameter * (h + v) + 2 * h * v * z), 1e-9,
             "four-sided area derivative");
        near(sample.volume,
             quarterPi *
                 (diameter * diameter * z + diameter * (h + v) * z * z / 2 + h * v * z * z * z / 3),
             1e-11, "four-sided volume from final geometry");
        near(sample.horizontalWallAngle, p.horizontalAngle, 1e-14, "right wall angle");
        near(sample.verticalWallAngle, p.verticalAngle, 1e-14, "top wall angle");
        near(sample.leftWallAngle, p.leftAngle, 1e-14, "left wall angle");
        near(sample.bottomWallAngle, p.bottomAngle, 1e-14, "bottom wall angle");
    }
    rejects([&] { (void)DraftGeometry(p, cone, cone, nullptr, cone); },
            "missing active direction provider rejected");
}

double pointTriangleDistance(Vec3 p, Vec3 a, Vec3 b, Vec3 c) {
    const auto minus = [](Vec3 u, Vec3 v) { return Vec3{u.x - v.x, u.y - v.y, u.z - v.z}; };
    const auto product = [](Vec3 u, Vec3 v) { return u.x * v.x + u.y * v.y + u.z * v.z; };
    const auto cross = [](Vec3 u, Vec3 v) {
        return Vec3{u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
    };
    const auto norm = [&](Vec3 u) { return std::hypot(u.x, u.y, u.z); };
    const auto n = cross(minus(b, a), minus(c, a));
    if (product(cross(minus(b, a), minus(p, a)), n) >= 0 &&
        product(cross(minus(c, b), minus(p, b)), n) >= 0 &&
        product(cross(minus(a, c), minus(p, c)), n) >= 0)
        return std::abs(product(minus(p, a), n)) / norm(n);
    double distance = std::numeric_limits<double>::infinity();
    for (const auto edge : {std::pair{a, b}, std::pair{b, c}, std::pair{c, a}}) {
        const auto d = minus(edge.second, edge.first), q = minus(p, edge.first);
        const double t = std::clamp(product(d, q) / product(d, d), 0.0, 1.0);
        distance = std::min(distance, norm({q.x - t * d.x, q.y - t * d.y, q.z - t * d.z}));
    }
    return distance;
}
void generalTessellationReferences() {
    const auto check = [](const DraftGeometry &geometry) {
        constexpr double tolerance = .00025;
        const auto mesh = geometry.previewMesh(tolerance);
        for (auto index : mesh.fixedVertices) {
            require(index < mesh.vertices.size(), "valid physical fixed index");
            near(mesh.vertices[index].z, 0, 0, "only physical driven perimeter fixed");
        }
        for (const auto &t : mesh.triangles)
            require(t.faceId >= 1 && t.faceId <= 5, "no artificial wall patch");
        for (unsigned i = 0; i < 11; ++i)
            for (unsigned j = 0; j < 24; ++j) {
                const double z = geometry.length() * (i + .37) / 11.0;
                const double phi = 2 * std::numbers::pi * (j + .23) / 24.0;
                for (bool outer : {false, true}) {
                    const double r =
                        geometry.radius(z, phi) + (outer ? geometry.parameters().wallThickness : 0);
                    const Vec3 p{r * std::cos(phi), r * std::sin(phi), z};
                    double error = std::numeric_limits<double>::infinity();
                    for (const auto &t : mesh.triangles)
                        if (t.faceId == static_cast<std::uint32_t>(outer ? DraftFace::outerWall
                                                                         : DraftFace::innerWall))
                            error = std::min(error,
                                             pointTriangleDistance(p, mesh.vertices[t.vertices[0]],
                                                                   mesh.vertices[t.vertices[1]],
                                                                   mesh.vertices[t.vertices[2]]));
                    require(error <= tolerance, "full-surface geometric deviation away from axes");
                }
            }
        closedPhysicalShell(geometry);
    };
    const DraftGeometry defaultGeometry;
    check(defaultGeometry);
    // This is the actual new-project GUI action, including the published ATH
    // superformula defaults, not only an evaluable radial guide definition.
    DraftParameters defaultGuide;
    defaultGuide.guideEnabled = true;
    defaultGuide.guide.type = CurveType::superformula;
    check(DraftGeometry(defaultGuide));
    DraftParameters p;
    p.entry = {.032, .008, 8};
    p.profile.throatAngle = .12;
    p.verticalProfile.throatAngle = .21;
    p.mouthMorph.enabled = true;
    p.mouthMorph.target.width = .34;
    p.mouthMorph.target.aspectRatio = .65;
    p.mouthMorph.target.rotation = .17;
    p.mouthMorph.target.superellipseExponent = 4;
    check(DraftGeometry(p));
    p = {};
    p.asymmetricEnabled = true;
    p.leftAngle = p.horizontalAngle * .75;
    p.bottomAngle = p.verticalAngle * .65;
    check(DraftGeometry(p));
    p.transition.method = TransitionMethod::lut;
    p.transition.knots = {{0., 0.}, {.3, .1}, {.7, .8}, {1., 1.}};
    p.crossSectionExponent = 4.8;
    check(DraftGeometry(p));
    p = {};
    p.length = .23;
    p.entry.width = p.entry.height = .02;
    p.horizontalAngle = p.verticalAngle = .4;
    check(DraftGeometry(p, std::make_shared<IndependentCone>()));
    const auto loose = defaultGeometry.previewMesh(.001),
               strict = defaultGeometry.previewMesh(.00005);
    require(strict.triangles.size() > loose.triangles.size(),
            "global detail tolerance refines the geometry");
    rejects([&] { (void)defaultGeometry.previewMesh(0); }, "zero geometric tolerance rejected");
    rejects([&] { (void)defaultGeometry.previewMesh(std::numeric_limits<double>::quiet_NaN()); },
            "nonfinite geometric tolerance rejected");
    std::stop_source cancellation;
    cancellation.request_stop();
    bool cancelled = false;
    try {
        (void)defaultGeometry.previewMesh(.00025, cancellation.get_token());
    } catch (const std::runtime_error &) {
        cancelled = true;
    }
    require(cancelled, "pre-cancelled source tessellation stops before producing a mesh");
}

void localContourSampling() {
    DraftParameters p;
    const auto original = DraftGeometry(p).previewMesh(.00025);
    p.mouthMorph.enabled = true;
    p.mouthMorph.target.width = .6;
    p.mouthMorph.target.rotation = .3;
    p.mouthMorph.target.superellipseExponent = 6;
    const auto morphed = DraftGeometry(p).previewMesh(.00025);
    require(morphed.fixedVertices.size() <= original.fixedVertices.size() + 4,
            "downstream morph detail does not dictate entry contour density");
    const double radius = .5 * p.entry.width;
    const auto chordCount = static_cast<std::size_t>(
        std::ceil(std::numbers::pi / std::acos(1 - .00025 / (4 * radius))));
    require(original.fixedVertices.size() <= 2 * chordCount,
            "circular contour sampling scales with its own chord deviation");
    std::map<double, std::set<std::uint32_t>> rings;
    for (const auto &triangle : morphed.triangles)
        if (triangle.faceId == static_cast<std::uint32_t>(DraftFace::innerWall))
            for (auto index : triangle.vertices)
                rings[morphed.vertices[index].z].insert(index);
    require(rings.begin()->second.size() < rings.rbegin()->second.size(),
            "each axial contour has an independent angular resolution");
    closedPhysicalShell(DraftGeometry(p));
}

void entryReferences() {
    DraftParameters p;
    p.entry = {.032, .008, 8.0};
    p.profile.throatAngle = .12;
    p.verticalProfile.throatAngle = .21;
    p.crossSectionExponent = 3.0;
    p.transition.startFraction = .15;
    p.transition.endFraction = .8;
    const Superellipse entry(.016, .004, 8.0);
    const auto checkEntry = [&](const DraftGeometry &geometry) {
        for (unsigned i = 0; i < 360; ++i) {
            const double phi = i * std::numbers::pi / 180.0;
            near(geometry.radius(0.0, phi), entry.radius(phi), 1e-13,
                 "actual entry contour independent of target and transition");
        }
        near(std::atan(geometry.radiusDerivative(0.0, 0.0)), p.profile.throatAngle, 1e-14,
             "rectangular entry horizontal tangent");
        near(std::atan(geometry.radiusDerivative(0.0, std::numbers::pi / 2.0)),
             p.verticalProfile.throatAngle, 1e-14, "rectangular entry vertical tangent");
        const auto mesh = geometry.previewMesh();
        double drivenArea = 0.0;
        for (const auto &t : mesh.triangles) {
            if (t.boundary != BoundaryCondition::driven)
                continue;
            const auto a = mesh.vertices[t.vertices[0]], b = mesh.vertices[t.vertices[1]],
                       c = mesh.vertices[t.vertices[2]];
            drivenArea += .5 * std::abs((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
        }
        const double exactArea = p.entry.width * p.entry.height *
                                 std::pow(std::tgamma(1.0 + 1.0 / p.entry.exponent), 2) /
                                 std::tgamma(1.0 + 2.0 / p.entry.exponent);
        // The inscribed contour is governed by a length tolerance, not by a
        // hidden relative-area target. Convex perimeter is bounded by its box.
        const double areaBound = 2 * (p.entry.width + p.entry.height) * .00025;
        require(drivenArea <= exactArea && exactArea - drivenArea < areaBound,
                "driven polygon area follows the common geometric tolerance");
        closedPhysicalShell(geometry);
    };
    checkEntry(DraftGeometry(p));
    p.transition.method = TransitionMethod::lut;
    p.transition.knots = {{0.0, 0.0}, {.3, .1}, {.7, .8}, {1.0, 1.0}};
    p.guideEnabled = true;
    p.guide.position = .7;
    p.guide.width = .22;
    p.guide.aspectRatio = .7;
    p.guide.rotation = .18;
    const DraftGeometry guided(p);
    checkEntry(guided);
    for (unsigned i = 0; i < 100; ++i) {
        const double phi = i * .071;
        near(guided.radius(guided.guidePosition(), phi), guided.guideRadius(phi), 1e-12,
             "rectangular entry retains the rotated guide");
    }
    p.guideEnabled = false;
    p.entry.width = 0.0;
    rejects([&] { (void)DraftGeometry(p); }, "zero entry width rejected");
    p.entry.width = .032;
    p.entry.exponent = 1.0;
    rejects([&] { (void)DraftGeometry(p); }, "unsupported entry exponent rejected");
    p.entry.exponent = 8.0;
    p.length = .23;
    const DraftGeometry injected(p, std::make_shared<IndependentCone>(.016),
                                 std::make_shared<IndependentCone>(.004));
    near(injected.radius(0.0, .71), entry.radius(.71), 1e-13,
         "independent profile providers retain arbitrary entry axes");
    rejects([&] { (void)DraftGeometry(p, std::make_shared<IndependentCone>()); },
            "injected profile with inconsistent entry dimensions rejected");
}

void mouthMorphReferences() {
    DraftParameters p;
    p.entry = {.03, .012, 4.0};
    p.profile.throatAngle = .12;
    p.verticalProfile.throatAngle = .08;
    p.guideEnabled = true;
    p.guide.position = .45;
    p.guide.width = .15;
    p.guide.aspectRatio = .7;
    p.guide.rotation = .11;
    const DraftGeometry base(p);
    p.mouthMorph.enabled = true;
    p.mouthMorph.startFraction = .65;
    p.mouthMorph.endFraction = .95;
    p.mouthMorph.target.width = .26;
    p.mouthMorph.target.aspectRatio = .75;
    p.mouthMorph.target.rotation = .23;
    p.mouthMorph.target.superellipseExponent = 6.0;
    const DraftGeometry morphed(p);
    const PhysicalContour target(p.mouthMorph.target);
    for (unsigned i = 0; i < 100; ++i) {
        const double phi = i * .071;
        for (const double fraction : {0.0, p.guide.position, p.mouthMorph.startFraction}) {
            const double z = fraction * base.length();
            near(morphed.radius(z, phi), base.radius(z, phi), 1e-13,
                 "mouth morph preserves entry, guide and transition start");
            near(morphed.radiusDerivative(z, phi), base.radiusDerivative(z, phi), 1e-12,
                 "mouth morph preserves boundary tangents");
        }
        near(morphed.radius(morphed.length(), phi), target.radius(phi), 1e-13,
             "mouth morph reaches the exact physical target");
        near(morphed.radiusDerivative(morphed.length(), phi),
             base.radiusDerivative(base.length(), phi), 1e-12,
             "mouth morph preserves the original mouth tangent");
        for (const double fraction : {.69, .81, .93}) {
            const double z = fraction * morphed.length(), step = 1e-7 * morphed.length();
            near(morphed.radiusDerivative(z, phi),
                 (morphed.radius(z + step, phi) - morphed.radius(z - step, phi)) / (2.0 * step),
                 2e-7, "mouth morph analytical derivative");
        }
        const double settledZ = .98 * morphed.length();
        near(morphed.radius(settledZ, phi),
             base.radius(settledZ, phi) + target.radius(phi) - base.radius(base.length(), phi),
             1e-13, "full mouth correction retains the evolving longitudinal profile");
        near(morphed.radiusDerivative(settledZ, phi), base.radiusDerivative(settledZ, phi), 1e-12,
             "settled correction retains the longitudinal derivative");
    }
    near(morphed.radius(morphed.length(), p.mouthMorph.target.rotation), .13, 1e-13,
         "rotated physical mouth width");
    closedPhysicalShell(morphed);
    auto scaled = p;
    scaled.length *= 2.0;
    scaled.entry.width *= 2.0;
    scaled.entry.height *= 2.0;
    scaled.wallThickness *= 2.0;
    scaled.guide.width *= 2.0;
    scaled.mouthMorph.target.width *= 2.0;
    const DraftGeometry twice(scaled);
    for (const double fraction : {.1, .75, 1.0}) {
        near(twice.radius(fraction * twice.length(), .71),
             2.0 * morphed.radius(fraction * morphed.length(), .71), 1e-12,
             "mouth morph respects SI scaling");
        near(twice.radiusDerivative(fraction * twice.length(), .71),
             morphed.radiusDerivative(fraction * morphed.length(), .71), 1e-11,
             "mouth morph derivative respects SI scaling");
    }
    {
        auto invalid = p;
        invalid.mouthMorph.target.type = CurveType::superformula;
        invalid.mouthMorph.target.superformula = {.8, 2.7, 4.8, 1.9, 2.9, .3, 3.6};
        bool contextual = false;
        try {
            (void)DraftGeometry(invalid);
        } catch (const std::invalid_argument &error) {
            contextual = std::string(error.what()).starts_with("Mund-Morphing: ") &&
                         std::string(error.what()).find("m1/m2") != std::string::npos;
        }
        require(contextual, "Nonperiodic mouth controls must identify their input section.");
    }
    p.mouthMorph.target.type = CurveType::superformula;
    p.mouthMorph.target.superformula.n1 = 3.0;
    p.mouthMorph.target.superformula.n2 = p.mouthMorph.target.superformula.n3 = 4.0;
    const DraftGeometry sf(p);
    const PhysicalContour sfTarget(p.mouthMorph.target);
    for (unsigned i = 0; i < 100; ++i)
        near(sf.radius(sf.length(), i * .071), sfTarget.radius(i * .071), 1e-13,
             "superformula mouth target is separate from the guide");
    p.mouthMorph.startFraction = p.guide.position;
    const DraftGeometry startsAtGuide(p);
    near(startsAtGuide.radius(startsAtGuide.guidePosition(), .71), startsAtGuide.guideRadius(.71),
         1e-13, "mouth morph may start exactly at the guide");
    near(startsAtGuide.radiusDerivative(startsAtGuide.guidePosition(), .71),
         base.radiusDerivative(base.guidePosition(), .71), 1e-12,
         "guide tangent retained at shared morph boundary");
    p.mouthMorph.startFraction = .4;
    rejects([&] { (void)DraftGeometry(p); }, "mouth morph before guide rejected explicitly");
    p.guideEnabled = false;
    p.mouthMorph.startFraction = 0.0;
    p.mouthMorph.endFraction = .5;
    p.mouthMorph.target.width = .001;
    rejects([&] { (void)DraftGeometry(p); },
            "mouth morph with negative intermediate radius rejected");
    p = {};
    p.length = .23;
    p.entry.width = p.entry.height = .02;
    p.horizontalAngle = p.verticalAngle = .4;
    p.mouthMorph.enabled = true;
    const DraftGeometry injected(p, std::make_shared<IndependentCone>());
    near(injected.radiusDerivative(injected.length(), .71), std::tan(.4), 1e-13,
         "mouth transform remains independent of longitudinal profile type");
}
} // namespace

int main() {
    try {
        profileReferences();
        curveReferences();
        guideControlSensitivity();
        transitions();
        guidesAndComposition();
        analysisReferences();
        directionalProfiles();
        asymmetricProfiles();
        asymmetricAnalysis();
        closedPhysicalShell();
        generalTessellationReferences();
        localContourSampling();
        entryReferences();
        mouthMorphReferences();
        std::cout << "Geometry, guide, transition, analysis and closed-shell tests passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
