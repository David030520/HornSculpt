#include "draft_wizard.hpp"

#include <QAbstractButton>
#include <QApplication>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QLocale>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
void near(double actual, double expected, const char *message) {
    require(std::abs(actual - expected) < 1e-12, message);
}
} // namespace

int main(int argc, char **argv) {
    QApplication app(argc, argv);
    try {
        using namespace hornsculpt;
        BasicDraftSpecification request{35.5, 240, 90, 60, 4, 4.5};
        const auto p = makeBasicDraft(request);
        near(p.entry.width, .0355, "Diameter must convert mm to m exactly once");
        near(p.entry.height, .0355, "Entry must remain circular");
        near(p.length, .24, "Length must use SI units");
        near(p.horizontalAngle, std::numbers::pi / 4, "Horizontal full angle must become half angle");
        near(p.verticalAngle, std::numbers::pi / 6, "Vertical full angle must become half angle");
        near(p.wallThickness, .0045, "Wall must use SI units");
        require(!p.guideEnabled && !p.mouthMorph.enabled && !p.slot.enabled && !p.asymmetricEnabled,
                "A basic draft must not contain hidden transformations");
        const DraftGeometry geometry(p);
        near(geometry.radius(0, 0), .0355 / 2, "The generated horn must meet the driver diameter");
        for (const auto invalid : {0.0, -2.0, std::numeric_limits<double>::quiet_NaN()}) {
            auto bad = request;
            bad.throatDiameterMm = invalid;
            bool rejected = false;
            try { (void)makeBasicDraft(bad); } catch (const std::invalid_argument &) { rejected = true; }
            require(rejected, "Invalid physical dimensions must be rejected");
        }
        for (const auto invalid : {-1.0, 180.0, std::numeric_limits<double>::infinity()}) {
            auto bad = request;
            bad.horizontalOpeningDegrees = invalid;
            bool rejected = false;
            try { (void)makeBasicDraft(bad); } catch (const std::invalid_argument &) { rejected = true; }
            require(rejected, "Invalid full angles must be rejected without correction");
        }
        DraftWizard wizard;
        wizard.setLocale(QLocale(QLocale::German, QLocale::Germany));
        wizard.show();
        app.processEvents();
        wizard.findChild<QDoubleSpinBox *>("wizardThroatDiameter")->setValue(35.5);
        wizard.findChild<QDoubleSpinBox *>("wizardLength")->setValue(240);
        wizard.findChild<QDoubleSpinBox *>("wizardWallThickness")->setValue(4.5);
        wizard.button(QWizard::NextButton)->click();
        wizard.findChild<QDoubleSpinBox *>("wizardHorizontalOpening")->setValue(90);
        wizard.findChild<QDoubleSpinBox *>("wizardVerticalOpening")->setValue(60);
        wizard.findChild<QComboBox *>("wizardContour")->setCurrentIndex(1);
        wizard.button(QWizard::NextButton)->click();
        require(wizard.currentId() == 2, "Native wizard must reach the review page");
        require(wizard.findChild<QLabel *>("wizardSummary")->text().contains("35,5"),
                "Summary must contain the actual edited dimensions");
        wizard.button(QWizard::BackButton)->click();
        require(wizard.findChild<QComboBox *>("wizardContour")->currentIndex() == 1,
                "Navigating back must preserve the user's contour choice");
        wizard.button(QWizard::NextButton)->click();
        wizard.button(QWizard::NextButton)->click();
        require(wizard.currentId() == 3 && wizard.startCalculation(),
                "Wizard must offer a real calculation after review");
        require(wizard.findChild<QWidget *>("wizardGeometryPreview"),
                "Wizard must expose its geometric preview");
        wizard.findChild<QDoubleSpinBox *>("wizardMinimumFrequency")->setValue(1000);
        wizard.findChild<QDoubleSpinBox *>("wizardMaximumFrequency")->setValue(500);
        wizard.button(QWizard::FinishButton)->click();
        require(wizard.result() != QDialog::Accepted, "Invalid frequency raster must be rejected");
        wizard.findChild<QDoubleSpinBox *>("wizardMaximumFrequency")->setValue(4000);
        wizard.button(QWizard::FinishButton)->click();
        require(wizard.result() == QDialog::Accepted, "Validated wizard must accept explicitly");
        near(wizard.parameters().horizontalAngle, p.horizontalAngle, "Wizard must generate the reviewed geometry");
        near(wizard.parameters().crossSectionExponent, 4, "Wizard must preserve the selected contour");
        DraftWizard cancelled;
        cancelled.show();
        cancelled.button(QWizard::CancelButton)->click();
        require(cancelled.result() == QDialog::Rejected, "Cancel must reject the draft proposal");
        std::puts("Draft wizard units, geometry, validation, review, back navigation and cancellation passed.");
        return 0;
    } catch (const std::exception &error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
