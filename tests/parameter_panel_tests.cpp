#include "parameter_panel.hpp"

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QPixmap>
#include <QPushButton>
#include <QScrollArea>
#include <QTabWidget>
#include <QTableWidget>
#include <cmath>
#include <hornsculpt/osse_profile.hpp>
#include <iostream>
#include <numbers>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
QDoubleSpinBox *number(ParameterPanel &panel, const char *name) {
    auto *result = panel.findChild<QDoubleSpinBox *>(QString::fromLatin1(name));
    require(result != nullptr, "The tested parameter input is missing.");
    return result;
}
DraftParameters exactHorn() {
    DraftParameters p;
    const double areaFactor = std::pow(std::tgamma(1.125), 2) / std::tgamma(1.25);
    p.entry = {.015306 / (.05 * areaFactor), .05, 8};
    p.length = .8;
    p.profile.terminationScale = p.verticalProfile.terminationScale = 0;
    p.crossSectionExponent = 8;
    OsseParameters horizontal;
    static_cast<OsseDirectionalParameters &>(horizontal) = p.profile;
    horizontal.throatRadius = p.entry.width / 2;
    horizontal.length = p.length;
    p.horizontalAngle = OsseProfile(horizontal).angleThrough(p.length, .75);
    p.verticalAngle = 0;
    p.wallThickness = .006123456789123;
    p.transition.startFraction = .012345678912345;
    p.transition.endFraction = .876543219876543;
    p.transition.knots = {{0, 0}, {.456789123456789, .234567891234567}, {1, 1}};
    p.guide.position = .678912345678912;
    p.guide.width = .301234567891234;
    p.guide.rotation = .123456789123456;
    p.guide.superformula.n1 = 4.123456789123456;
    p.mouthMorph.target.width = .456789123456789;
    p.mouthMorph.target.superformula.n2 = 3.234567891234567;
    return p;
}
void retainedValues(const DraftParameters &actual, const DraftParameters &expected) {
    require(actual.entry.width == expected.entry.width &&
                actual.entry.height == expected.entry.height &&
                actual.entry.exponent == expected.entry.exponent,
            "Displaying imported throat dimensions rounded the model.");
    require(actual.horizontalAngle == expected.horizontalAngle &&
                actual.verticalAngle == expected.verticalAngle,
            "Displaying the computed opening angle rounded the model.");
    require(actual.transition.startFraction == expected.transition.startFraction &&
                actual.transition.endFraction == expected.transition.endFraction &&
                actual.transition.knots.size() == expected.transition.knots.size(),
            "An unrelated edit changed exact shape-transition parameters.");
    for (std::size_t i = 0; i < expected.transition.knots.size(); ++i)
        require(actual.transition.knots[i].fraction == expected.transition.knots[i].fraction &&
                    actual.transition.knots[i].blend == expected.transition.knots[i].blend,
                "Displaying a LUT rounded an untouched knot.");
    require(actual.guide.width == expected.guide.width &&
                actual.guide.rotation == expected.guide.rotation &&
                actual.guide.superformula.n1 == expected.guide.superformula.n1 &&
                actual.mouthMorph.target.width == expected.mouthMorph.target.width &&
                actual.mouthMorph.target.superformula.n2 ==
                    expected.mouthMorph.target.superformula.n2,
            "An unrelated edit rounded guide or mouth parameters.");
}
void precisionAndIndependentEdits() {
    ParameterPanel panel;
    const auto original = exactHorn();
    int changes = 0;
    panel.onChanged = [&] { ++changes; };
    panel.setParameters(original);
    require(changes == 0, "Loading parameters emitted an edit.");
    auto result = panel.parameters();
    retainedValues(result, original);
    require(result.length == original.length && result.wallThickness == original.wallThickness &&
                result.guide.position == original.guide.position,
            "Untouched SI dimensions changed during display conversion.");
    require(number(panel, "entryWidth")->value() != original.entry.width * 1000 &&
                number(panel, "angleH")->value() !=
                    original.horizontalAngle * 180 / std::numbers::pi,
            "The fixture did not exercise display rounding.");
    const DraftGeometry before(original), after(result);
    require(before.radius(original.length, 0) == after.radius(result.length, 0),
            "An unchanged imported project would reconstruct a different mouth.");

    number(panel, "thickness")->setValue(7.25);
    result = panel.parameters();
    retainedValues(result, original);
    require(changes == 1 && result.wallThickness == 7.25 / 1000 && result.length == original.length,
            "Editing wall thickness did not affect exactly its own model value.");
    number(panel, "length")->setValue(850.12);
    result = panel.parameters();
    retainedValues(result, original);
    require(changes == 2 && result.length == 850.12 / 1000 && result.wallThickness == 7.25 / 1000,
            "Editing length changed another precise imported field.");
    panel.setParameters(original);
    result = panel.parameters();
    retainedValues(result, original);
    require(result.length == original.length && result.wallThickness == original.wallThickness,
            "Reloading exact parameters retained an earlier edited numeric value.");
}
void lutAndPositionMode() {
    ParameterPanel panel;
    const auto original = exactHorn();
    panel.setParameters(original);
    auto *table = panel.findChild<QTableWidget *>();
    require(table != nullptr, "The shape LUT is missing.");
    table->item(1, 1)->setText("37.25");
    auto result = panel.parameters();
    auto expected = original;
    expected.transition.knots[1].blend = .3725;
    retainedValues(result, expected);
    require(result.wallThickness == original.wallThickness,
            "Editing a LUT cell rounded an unrelated wall dimension.");

    auto *mode = panel.findChild<QComboBox *>("guidePositionMode");
    require(mode != nullptr, "The guide position mode is missing.");
    mode->setCurrentIndex(1);
    result = panel.parameters();
    require(result.guide.positionMode == GuidePositionMode::absolute &&
                result.guide.position == original.guide.position * original.length,
            "Relative-to-absolute conversion used rounded display values.");
    retainedValues(result, expected);
    mode->setCurrentIndex(0);
    result = panel.parameters();
    require(std::abs(result.guide.position - original.guide.position) < 2e-16,
            "Guide position conversion lost precision beyond floating-point arithmetic.");
    number(panel, "guidePosition")->setValue(72.123);
    number(panel, "length")->setValue(811.11);
    mode->setCurrentIndex(1);
    result = panel.parameters();
    require(result.guide.position == (72.123 / 100) * (811.11 / 1000),
            "Position-mode conversion ignored deliberately edited position or length.");
    retainedValues(result, expected);
}
void independentSides() {
    ParameterPanel panel;
    auto p = exactHorn();
    p.asymmetricEnabled = true;
    p.leftAngle = .345678912345678;
    p.bottomAngle = .234567891234567;
    p.leftProfile.throatAngle = .031234567891234;
    p.bottomProfile.throatAngle = .021234567891234;
    p.leftProfile.k = .612345678912345;
    p.bottomProfile.k = .712345678912345;
    p.leftProfile.terminationScale = .123456789123456;
    p.bottomProfile.terminationScale = .234567891234567;
    p.leftProfile.terminationQ = .912345678912345;
    p.bottomProfile.terminationQ = .923456789123456;
    p.leftProfile.terminationExponent = 3.123456789123456;
    p.bottomProfile.terminationExponent = 4.123456789123456;
    panel.setParameters(p);
    int changes = 0;
    panel.onChanged = [&] { ++changes; };
    auto *mode = panel.findChild<QCheckBox *>("asymmetricEnabled");
    auto *copy = panel.findChild<QPushButton *>("copyOppositeProfiles");
    auto *guide = panel.findChild<QCheckBox *>("guideEnabled");
    require(mode && copy && guide, "Asymmetry controls are missing.");
    auto result = panel.parameters();
    require(result.asymmetricEnabled && result.leftAngle == p.leftAngle &&
                result.bottomAngle == p.bottomAngle &&
                result.leftProfile.terminationQ == p.leftProfile.terminationQ &&
                result.bottomProfile.terminationExponent == p.bottomProfile.terminationExponent,
            "Displaying four sides rounded precise profile values.");
    number(panel, "angleL")->setValue(31.25);
    number(panel, "sB")->setValue(.314);
    result = panel.parameters();
    retainedValues(result, p);
    require(changes == 2 && result.leftAngle == 31.25 / (180 / std::numbers::pi) &&
                result.bottomAngle == p.bottomAngle &&
                result.bottomProfile.terminationScale == .314 &&
                result.leftProfile.terminationScale == p.leftProfile.terminationScale,
            "An opposite-side edit affected another direction.");
    mode->setChecked(false);
    require(!panel.parameters().asymmetricEnabled && !number(panel, "angleL")->isEnabled(),
            "Symmetric mode did not disable opposite-side controls.");
    mode->setChecked(true);
    require(panel.parameters().leftAngle == result.leftAngle &&
                panel.parameters().bottomProfile.terminationScale == .314,
            "Mode switching discarded the stored opposite profiles.");
    guide->setChecked(true);
    require(number(panel, "angleH")->isEnabled() && number(panel, "angleV")->isEnabled() &&
                number(panel, "angleL")->isEnabled() && number(panel, "angleB")->isEnabled() &&
                number(panel, "throatAngleB")->isEnabled() &&
                !number(panel, "crossExponent")->isEnabled() &&
                number(panel, "guideExponent")->isEnabled() &&
                number(panel, "transitionStart")->isEnabled(),
            "Guide controls must expose the active base profile and hide the replaced contour.");
    guide->setChecked(false);
    require(number(panel, "crossExponent")->isEnabled(),
            "Direct contour did not become editable after disabling the guide.");
    const int beforeCopy = changes;
    copy->click();
    result = panel.parameters();
    require(changes == beforeCopy + 1 && result.leftAngle == p.horizontalAngle &&
                result.bottomAngle == p.verticalAngle &&
                result.leftProfile.throatAngle == p.profile.throatAngle &&
                result.bottomProfile.terminationScale == p.verticalProfile.terminationScale,
            "Explicit side copy must preserve exact values and emit one edit.");
    retainedValues(result, p);
    panel.setParameters(p);
    require(panel.parameters().leftAngle == p.leftAngle &&
                panel.parameters().bottomProfile.terminationScale ==
                    p.bottomProfile.terminationScale,
            "Reload/undo retained edited values from the opposite sides.");
}

void sectionNavigation() {
    ParameterPanel panel;
    auto p = exactHorn();
    p.guideEnabled = true;
    panel.setParameters(p);
    panel.resize(340, 650);
    panel.show();
    QCoreApplication::processEvents();
    auto *tabs = panel.findChild<QTabWidget *>("draftSections");
    auto *guideLink = panel.findChild<QPushButton *>("guideAngleHint");
    require(tabs && tabs->count() == 5 && guideLink && guideLink->isVisible(),
            "Draft sections or guide-angle explanation are missing");
    int changes = 0;
    panel.onChanged = [&] { ++changes; };
    guideLink->click();
    require(tabs->currentIndex() == 2 && number(panel, "guideWidth")->isVisible(),
            "The angle explanation does not reveal the responsible guide inputs");
    auto *metrics = panel.findChild<QLabel *>("guideMetrics");
    auto *formLink = panel.findChild<QPushButton *>("guideFormLink");
    auto *shapeLink = panel.findChild<QPushButton *>("guideShapeHint");
    require(metrics && !metrics->text().isEmpty() && formLink && shapeLink,
            "Guide dimensions and shape navigation are missing");
    formLink->click();
    require(tabs->currentIndex() == 1 && shapeLink->isVisible(),
            "The guide does not link to its axial form schedule");
    shapeLink->click();
    require(tabs->currentIndex() == 2, "Direct contour does not link to the active guide");
    for (int i = 0; i < tabs->count(); ++i) {
        tabs->setCurrentIndex(i);
        require(qobject_cast<QScrollArea *>(tabs->currentWidget()) != nullptr,
                "A Draft section cannot scroll independently from its navigation");
    }
    require(changes == 0 && panel.parameters().guide.width == p.guide.width,
            "Section navigation must not edit or round model parameters");
    p.guideEnabled = false;
    p.mouthMorph.enabled = false;
    panel.setParameters(p);
    tabs->setCurrentIndex(2);
    require(!number(panel, "guideWidth")->isVisible(), "Inactive guide details take up space");
    tabs->setCurrentIndex(3);
    require(!number(panel, "mouthWidth")->isVisible(), "Inactive mouth details take up space");
    panel.findChild<QCheckBox *>("mouthMorphEnabled")->setChecked(true);
    require(number(panel, "mouthWidth")->isVisible(),
            "Enabling mouth morph does not reveal inputs");
}
void diffractionSlotControls() {
    ParameterPanel panel;
    DraftParameters p;
    p.slot.aperture = {.012345678912345, .034567891234567, 4.123456789123456};
    p.slot.adapterLength = .020123456789123;
    p.slot.straightLength = .003123456789123;
    p.slot.transitionLength = .005123456789123;
    panel.setParameters(p);
    auto *enabled = panel.findChild<QCheckBox *>("slotEnabled");
    auto *repair = panel.findChild<QPushButton *>("slotFormStart");
    auto *method = panel.findChild<QComboBox *>("transitionMethod");
    require(enabled && repair && method, "Diffraction slot controls are missing.");
    int changes = 0;
    panel.onChanged = [&] { ++changes; };
    enabled->setChecked(true);
    auto actual = panel.parameters();
    require(changes == 1 && actual.slot.enabled &&
                actual.transition.startFraction == p.transition.startFraction &&
                actual.slot.aperture.width == p.slot.aperture.width &&
                actual.slot.aperture.height == p.slot.aperture.height &&
                actual.slot.aperture.exponent == p.slot.aperture.exponent &&
                actual.slot.adapterLength == p.slot.adapterLength &&
                actual.slot.straightLength == p.slot.straightLength &&
                actual.slot.transitionLength == p.slot.transitionLength && repair->isEnabled(),
            "Slot activation rounded values or silently repaired the form transition.");
    repair->click();
    actual = panel.parameters();
    const double finish =
        (p.slot.adapterLength + p.slot.straightLength + p.slot.transitionLength) / p.length;
    require(changes == 2 && actual.transition.startFraction == finish &&
                actual.transition.endFraction == p.transition.endFraction &&
                actual.mouthMorph.startFraction == p.mouthMorph.startFraction &&
                actual.length == p.length && !repair->isEnabled(),
            "Explicit form repair must use exact SI values and change only the range start.");
    (void)DraftGeometry(actual);
    number(panel, "slotTransitionLength")->setValue(0);
    require(panel.parameters().slot.transitionLength == 0 &&
                panel.parameters().slot.adapterLength == p.slot.adapterLength,
            "A sharp slot exit cannot be entered independently.");
    enabled->setChecked(false);
    require(!number(panel, "slotWidth")->isEnabled(), "Disabled slot inputs remain active.");
    enabled->setChecked(true);
    require(panel.parameters().slot.transitionLength == 0 &&
                panel.parameters().slot.aperture.width == p.slot.aperture.width,
            "Toggling the slot discarded its retained settings.");
    method->setCurrentIndex(1);
    const auto lut = panel.parameters().transition.knots;
    require(!repair->isEnabled(), "The range repair must not replace a user's LUT.");
    repair->click();
    actual = panel.parameters();
    require(actual.transition.method == TransitionMethod::lut &&
                actual.transition.knots.size() == lut.size() &&
                actual.transition.knots[1].blend == lut[1].blend,
            "A slot action silently changed the selected LUT.");
    panel.setParameters(p);
    require(panel.parameters().slot.transitionLength == p.slot.transitionLength &&
                !panel.parameters().slot.enabled,
            "Reload/undo retained edited slot values.");
    DraftParameters ellipse;
    ellipse.slot.enabled = true;
    ellipse.slot.aperture.exponent = 2;
    panel.setParameters(ellipse);
    auto *metrics = panel.findChild<QLabel *>("slotMetrics");
    require(
        metrics && metrics->text().contains("329.9 mm²") && metrics->text().contains("157.000 mm"),
        "Slot summary must show the actual ellipse area and retain the blend within OS-SE length.");
}
void physicalSlotRadiusControls() {
    ParameterPanel panel;
    DraftParameters p;
    p.slot.enabled = true;
    p.slot.cornerRadius = .002345678912345;
    p.slot.exitRadius = .018123456789123;
    p.slot.transitionLength = .007123456789123;
    p.slot.aperture.exponent = 4.123456789123456;
    panel.setParameters(p);
    auto *contour = panel.findChild<QComboBox *>("slotContour");
    auto *exit = panel.findChild<QComboBox *>("slotExit");
    auto *repair = panel.findChild<QPushButton *>("slotFormStart");
    auto *metrics = panel.findChild<QLabel *>("slotMetrics");
    require(contour && exit && repair && metrics, "Physical slot radius controls are missing.");
    require(!number(panel, "slotCornerRadius")->isEnabled() &&
                !number(panel, "slotExitRadius")->isEnabled() &&
                number(panel, "slotExponent")->isEnabled() &&
                number(panel, "slotTransitionLength")->isEnabled() &&
                panel.parameters().slot.cornerRadius == p.slot.cornerRadius &&
                panel.parameters().slot.exitRadius == p.slot.exitRadius,
            "Dormant physical radii changed or replaced existing shape/length controls.");
    contour->setCurrentIndex(static_cast<int>(SlotContour::roundedRectangle));
    exit->setCurrentIndex(static_cast<int>(SlotExit::circularFillet));
    auto actual = panel.parameters();
    require(actual.slot.contour == SlotContour::roundedRectangle &&
                actual.slot.exit == SlotExit::circularFillet &&
                actual.slot.cornerRadius == p.slot.cornerRadius &&
                actual.slot.exitRadius == p.slot.exitRadius &&
                number(panel, "slotCornerRadius")->isEnabled() &&
                number(panel, "slotExitRadius")->isEnabled() &&
                !number(panel, "slotExponent")->isEnabled() &&
                !number(panel, "slotTransitionLength")->isEnabled(),
            "Selecting true radii changed retained SI values or enabled inactive parameters.");
    number(panel, "slotCornerRadius")->setValue(3.125);
    require(panel.parameters().slot.cornerRadius == 3.125 / 1000 &&
                panel.parameters().slot.exitRadius == p.slot.exitRadius,
            "Editing the corner radius changed the independent exit radius.");
    number(panel, "slotExitRadius")->setValue(18.25);
    require(panel.parameters().slot.exitRadius == 18.25 / 1000 &&
                panel.parameters().slot.cornerRadius == 3.125 / 1000,
            "A freely entered exit radius was clamped or changed the corner radius.");
    number(panel, "slotExitRadius")->setValue(0);
    require(panel.parameters().slot.exitRadius == 0 &&
                panel.parameters().slot.cornerRadius == 3.125 / 1000,
            "A sharp exit cannot be entered independently.");
    number(panel, "slotCornerRadius")->setValue(0);
    require(panel.parameters().slot.cornerRadius == 0 && panel.parameters().slot.exitRadius == 0,
            "A sharp rectangular cross-section cannot be entered.");
    number(panel, "slotCornerRadius")->setValue(18);
    require(panel.parameters().slot.cornerRadius == .018 &&
                metrics->text().contains("ungültig"),
            "An impossible corner radius was silently changed or crashed the input callback.");
    contour->setCurrentIndex(static_cast<int>(SlotContour::superellipse));
    exit->setCurrentIndex(static_cast<int>(SlotExit::axialBlend));
    actual = panel.parameters();
    require(actual.slot.cornerRadius == .018 && actual.slot.exitRadius == 0 &&
                actual.slot.aperture.exponent == p.slot.aperture.exponent &&
                actual.slot.transitionLength == p.slot.transitionLength,
            "Mode switching discarded dormant physical radii, exponent or smoothing length.");

    p.slot.contour = SlotContour::roundedRectangle;
    p.slot.cornerRadius = .002;
    p.slot.exit = SlotExit::circularFillet;
    p.slot.exitRadius = .002;
    p.slot.straightLength = .010;
    p.profile.k = p.verticalProfile.k = 0;
    p.profile.terminationScale = p.verticalProfile.terminationScale = 0;
    p.horizontalAngle = .5;
    p.verticalAngle = .3;
    p.transition.startFraction = .25;
    p.transition.endFraction = .65;
    const double actualJoin = DraftGeometry(p).slotTransitionEnd();
    p.transition.startFraction = 0;
    panel.setParameters(p);
    require(!repair->isEnabled(),
            "The radius repair guessed a join extent before geometry supplied it.");
    int changes = 0;
    panel.onChanged = [&] { ++changes; };
    panel.setSlotTransitionEnd(actualJoin);
    require(changes == 0 && repair->isEnabled(),
            "Computed join extent did not enable the explicit repair, or emitted an edit.");
    repair->click();
    actual = panel.parameters();
    require(changes == 1 && actual.transition.startFraction == actualJoin / p.length &&
                actual.slot.exitRadius == p.slot.exitRadius &&
                actual.slot.cornerRadius == p.slot.cornerRadius,
            "Form repair did not use the exact computed circular join.");
    (void)DraftGeometry(actual);
    panel.setSlotTransitionEnd(actualJoin);
    number(panel, "slotExitRadius")->setValue(2.5);
    require(!repair->isEnabled(), "An edit retained a stale circular join extent.");
    panel.setSlotTransitionEnd(actualJoin);
    panel.setParameters(p);
    require(!repair->isEnabled(), "Reload/undo retained a previous geometry's circular join.");
}
} // namespace
int main(int argc, char **argv) {
    QApplication application(argc, argv);
    try {
        if (argc == 2) {
            ParameterPanel preview;
            DraftParameters p;
            p.slot.enabled = true;
            p.transition.startFraction =
                (p.slot.adapterLength + p.slot.straightLength + p.slot.transitionLength) / p.length;
            preview.setParameters(p);
            preview.resize(440, 620);
            preview.findChild<QTabWidget *>("draftSections")->setCurrentIndex(4);
            preview.show();
            QApplication::processEvents();
            require(preview.grab().save(QString::fromLocal8Bit(argv[1])),
                    "The requested parameter-panel screenshot could not be saved.");
        }
        precisionAndIndependentEdits();
        lutAndPositionMode();
        independentSides();
        sectionNavigation();
        diffractionSlotControls();
        physicalSlotRadiusControls();
        std::cout << "Exact imported SI parameters, independent edits, LUT and guide conversion "
                     "passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
