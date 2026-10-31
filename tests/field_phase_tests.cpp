#include "acoustic_field_widget.hpp"
#include "fixture_geometry.hpp"
#include "phase_response_widget.hpp"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QDir>
#include <QDialog>
#include <QElapsedTimer>
#include <QLabel>
#include <QPixmap>
#include <QPushButton>
#include <QSettings>
#include <QSlider>
#include <QTemporaryDir>
#include <QToolButton>
#include <algorithm>
#include <cmath>
#include <hornsculpt/plot_widget.hpp>
#include <hornsculpt/mesh_viewport.hpp>
#include <hornsculpt/pressure_field_plane.hpp>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>

using namespace hornsculpt;
namespace {
QString screenshotDirectory;
void require(bool value, const char *message);
void screenshot(QWidget &widget, const QString &name) {
    if (screenshotDirectory.isEmpty())
        return;
    widget.resize(1200, 650);
    widget.show();
    QApplication::processEvents();
    require(widget.grab().save(QDir(screenshotDirectory).filePath(name)),
            "Cannot save native widget QA image.");
}
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
void near(double value, double expected) {
    require(std::isfinite(value) && std::abs(value - expected) < 1e-9,
            "Unexpected phase from complex pressure.");
}
void phaseControls() {
    PhaseResponseWidget widget;
    auto *h = dynamic_cast<PlotWidget *>(widget.findChild<QWidget *>("phaseResponseH"));
    auto *v = dynamic_cast<PlotWidget *>(widget.findChild<QWidget *>("phaseResponseV"));
    require(h && v, "Missing phase plots.");
    std::vector<PhaseResponseRow> rows;
    for (const auto [frequency, angle] : {std::pair{500., 170.}, {1000., -170.}, {2000., -150.}}) {
        const auto p = std::polar(1., angle * std::numbers::pi / 180);
        const auto shifted = p * std::polar(1., 20 * std::numbers::pi / 180);
        rows.push_back(
            {frequency, {p, shifted}, {p, p * std::polar(1., -30 * std::numbers::pi / 180)}});
    }
    widget.setRows({0, 45}, rows, "Rev. 9");
    near(h->series()[0].points[0].y(), 170);
    require(h->series()[0].points.size() == 4 && std::isnan(h->series()[0].points[1].y()),
            "Wrapped phase discontinuity was drawn as a diagonal ramp.");
    auto s = widget.settings();
    s.unwrap = true;
    widget.setSettings(s);
    require(h->series()[0].points.size() == 3, "Unwrapped phase has artificial points.");
    near(h->series()[0].points[1].y(), 190);
    near(h->series()[0].points[2].y(), 210);
    s.angleDegrees = 45;
    s.relative = true;
    s.referenceAngleDegrees = 0;
    widget.setSettings(s);
    for (const auto &p : h->series()[0].points)
        near(p.y(), 20);
    for (const auto &p : v->series()[0].points)
        near(p.y(), -30);
    s.delaySeconds = .002;
    widget.setSettings(s);
    near(h->series()[0].points[2].y(), 1460);
    s.delaySeconds = .001;
    s.unwrap = false;
    widget.setSettings(s);
    near(h->series()[0].points.front().y(), -160);
    rows[1].horizontal[0] = {};
    widget.setRows({0, 45}, rows, "Rev. 9 · veraltet");
    require(std::isnan(h->series()[0].points[1].y()), "Zero reference pressure acquired a phase.");
    screenshot(widget, "Phase-Referenz-QA.png");
    for (auto &row : rows) {
        std::fill(row.horizontal.begin(), row.horizontal.end(), std::complex<double>{});
        std::fill(row.vertical.begin(), row.vertical.end(), std::complex<double>{});
    }
    s.relative = false;
    widget.setSettings(s);
    widget.setRows({0, 45}, rows, "Rev. 9 · Nullanregung");
    for (const auto *plot : {h, v})
        for (const auto &point : plot->series()[0].points)
            require(std::isnan(point.y()), "A zero excitation fabricated an absolute phase.");
}
void planeContract() {
    PressureFieldPlane p;
    p.columns = p.rows = 3;
    p.minimumLateral = p.minimumAxial = -1;
    p.maximumLateral = p.maximumAxial = 1;
    p.frequency = 1000;
    for (const double z : {-1., 0., 1.})
        for (const double x : {-1., 0., 1.}) {
            p.positions.push_back({x, 0, z});
            p.pressure.push_back({1, 0});
            p.valid.push_back(true);
        }
    require(pressureFieldTriangleIndices(p).size() == 24, "Field grid lost triangles.");
    MeshViewport viewport;
    viewport.resize(640, 480);
    viewport.setPressureField(std::make_shared<PressureFieldPlane>(p));
    viewport.viewFront({0, 1, 0});
    viewport.focusPressureField();
    const auto centerSample = viewport.pickFieldSample(QPointF(320, 240));
    require(centerSample && *centerSample == 4, "3D probe did not select the actual center sample.");
    require(!viewport.pickFieldSample(QPointF(-1, -1)), "3D probe accepted an outside click.");
    auto screen = std::make_shared<SurfaceMesh>();
    screen->vertices = {{-1, .5, -1}, {1, .5, -1}, {0, .5, 1}};
    screen->triangles.push_back({{0, 1, 2}, 1, BoundaryCondition::wall});
    viewport.setMesh(screen);
    viewport.viewFront({0, 1, 0});
    viewport.focusPressureField();
    require(!viewport.pickFieldSample(QPointF(320, 240)), "3D probe selected a point behind CAD.");
    viewport.setMesh({});
    viewport.focusPressureField();
    auto zero = p;
    std::fill(zero.pressure.begin(), zero.pressure.end(), std::complex<double>{});
    require(pressureFieldTriangleIndices(zero).size() == 24,
            "A valid zero-pressure field was misclassified as a spatial mask.");
    p.valid[4] = false;
    viewport.setPressureField(std::make_shared<PressureFieldPlane>(p));
    require(!viewport.pickFieldSample(QPointF(320, 240)), "3D probe accepted a masked sample.");
    const auto masked = pressureFieldTriangleIndices(p);
    require(!masked.empty(), "A single invalid point discarded all valid field cells.");
    for (const auto index : masked)
        require(index != 4, "An invalid sample was included in field triangles.");
    const auto coherent = p;
    const auto rejects = [&](PressureFieldPlane invalid) {
        bool rejected = false;
        try {
            validatePressureFieldPlane(invalid);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "World field points incompatible with their frame/raster were accepted.");
    };
    auto misplaced = coherent;
    misplaced.frame.origin.x += .001;
    rejects(misplaced);
    misplaced = coherent;
    misplaced.maximumAxial += .01;
    rejects(misplaced);
    misplaced = coherent;
    misplaced.positions[4].y += .001; // Invalid samples must retain their world identity too.
    rejects(misplaced);
    auto rotated = coherent;
    rotated.frame = normalizedAcousticFrame({{1e6, -2e6, 3e6}, {.3, .4, .5}, {.8, .2, -.3}});
    rotated.vertical = true;
    const auto lateral = verticalDirection(rotated.frame);
    for (unsigned row = 0; row < rotated.rows; ++row)
        for (unsigned col = 0; col < rotated.columns; ++col) {
            const double axial = -1. + row, across = -1. + col;
            rotated.positions[row * rotated.columns + col] = {
                rotated.frame.origin.x + axial * rotated.frame.forward.x + across * lateral.x,
                rotated.frame.origin.y + axial * rotated.frame.forward.y + across * lateral.y,
                rotated.frame.origin.z + axial * rotated.frame.forward.z + across * lateral.z};
        }
    validatePressureFieldPlane(rotated);
    rotated.positions[0].z += .001;
    rejects(rotated);
    bool rejected = false;
    p.positions[0].x = std::numeric_limits<double>::quiet_NaN();
    try {
        validatePressureFieldPlane(p);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "The 3D field accepted a non-finite world coordinate.");
}
void asynchronousField() {
    AcousticFieldWidget widget;
    std::shared_ptr<const PressureFieldPlane> published;
    PressureFieldAppearance publishedAppearance;
    widget.onFieldChanged = [&](auto plane, auto appearance) {
        published = std::move(plane);
        publishedAppearance = appearance;
    };
    auto *calculate = widget.findChild<QPushButton *>("evaluateField");
    require(calculate && !calculate->isEnabled(), "Field evaluation is enabled without traces.");
    auto *probeControl = widget.findChild<QToolButton *>("fieldProbeMode");
    require(probeControl && !probeControl->isEnabled(), "Probe mode is enabled without a field.");
    widget.requestEvaluation();
    require(!widget.evaluationInProgress() && !widget.fieldReady() && !widget.lastError().isEmpty(),
            "Explicit field request without a result did not report its unavailable input.");
    auto mesh = std::make_shared<SurfaceMesh>(testing::sphere(1));
    BemRequest request;
    auto result = std::make_shared<BemResult>();
    result->frequencies.push_back({350});
    auto &f = result->frequencies.front();
    f.excitations.emplace_back();
    auto &e = f.excitations.front();
    e.name = "Analytic Cauchy trace";
    BemBoundaryData boundary;
    const double k = 2 * std::numbers::pi * f.frequency / request.medium.soundSpeed;
    for (const auto p : mesh->vertices)
        boundary.pressureTrace.push_back(std::polar(10., -k * .1));
    for (const auto &t : mesh->triangles) {
        const auto a = mesh->vertices[t.vertices[0]], b = mesh->vertices[t.vertices[1]],
                   c = mesh->vertices[t.vertices[2]];
        const Vec3 p{(a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3, (a.z + b.z + c.z) / 3};
        const double r = std::hypot(p.x, p.y, p.z);
        const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
        const Vec3 n{u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
        const double cosine = (p.x * n.x + p.y * n.y + p.z * n.z) / (r * std::hypot(n.x, n.y, n.z));
        boundary.normalDerivative.push_back(std::polar(1 / r, -k * r) * std::complex<double>{-1 / r, -k} *
                                     cosine);
    }
    e.setBoundary(std::move(boundary));
    f.excitations.push_back(f.excitations.front());
    f.excitations.back().name = "Second source";
    result->frequencies.push_back(f);
    result->frequencies.back().frequency = 700;
    auto settings = widget.settings();
    settings.frequencyHz = 350;
    settings.excitationIndex = 0;
    settings.lateralCount = settings.axialCount = 3;
    settings.minimumLateral = -.03;
    settings.maximumLateral = .03;
    settings.minimumAxial = .2;
    settings.maximumAxial = .3;
    widget.setSettings(settings);
    widget.setSnapshot(mesh, request, result, {}, 9, true);
    require(calculate->isEnabled(), "Actual Cauchy traces did not enable field evaluation.");
    calculate->click();
    QElapsedTimer timer;
    timer.start();
    while (!widget.fieldPlane() && timer.elapsed() < 5000)
        QApplication::processEvents();
    require(widget.fieldPlane() && widget.fieldPlane()->pressure.size() == 9 &&
                widget.fieldPlane()->valid[4] && widget.fieldReady() &&
                !widget.evaluationInProgress() && widget.lastError().isEmpty(),
            "Asynchronous field did not publish calculated samples.");
    require(published == widget.fieldPlane() && published->positions.size() == 9 &&
                published->positions[0] == Vec3{-.03, 0, .2} && published->frequency == 350 &&
                published->revision == 9 && published->current,
            "The 3D field event lost actual world positions or snapshot identity.");
    const auto pressure = widget.fieldPlane()->pressure;
    auto *quantity = widget.findChild<QComboBox *>("fieldQuantity");
    auto *minimumDb = widget.findChild<QDoubleSpinBox *>("fieldMinimumDb");
    auto *maximumDb = widget.findChild<QDoubleSpinBox *>("fieldMaximumDb");
    auto *phase = widget.findChild<QDoubleSpinBox *>("fieldPhase");
    auto *exclusion = widget.findChild<QDoubleSpinBox *>("fieldExclusion");
    require(quantity && minimumDb && maximumDb && phase && exclusion,
            "Field display controls are missing.");
    require(minimumDb->isEnabled() && maximumDb->isEnabled() && !phase->isEnabled() &&
                exclusion->isEnabled(), "Level display enabled inappropriate field inputs.");
    quantity->setCurrentIndex(1);
    require(!minimumDb->isEnabled() && !maximumDb->isEnabled() && !phase->isEnabled() &&
                exclusion->isEnabled(), "Phase display enabled inappropriate field inputs.");
    require(widget.fieldPlane()->pressure == pressure,
            "Phase display unexpectedly recalculated or changed pressure samples.");
    require(published == widget.fieldPlane() && publishedAppearance.quantity == 1,
            "A phase-only display update did not reach the 3D field appearance.");
    quantity->setCurrentIndex(2);
    require(!minimumDb->isEnabled() && maximumDb->isEnabled() && phase->isEnabled() &&
                exclusion->isEnabled() && widget.fieldPlane()->pressure == pressure,
            "Real-part display enabled inappropriate inputs or recalculated the field.");
    auto *frequency = widget.findChild<QComboBox *>("fieldFrequency");
    auto *excitation = widget.findChild<QComboBox *>("fieldExcitation");
    auto *slider = widget.findChild<QSlider *>("fieldFrequencySlider");
    auto *animate = widget.findChild<QPushButton *>("fieldAnimate");
    require(frequency && excitation && slider && animate, "Field selectors missing");
    const auto firstField = widget.fieldPlane();
    slider->setValue(1);
    require(!widget.fieldPlane() && frequency->currentIndex() == 1,
            "Uncached frequency did not clear the previous field.");
    slider->setValue(0);
    require(widget.fieldPlane() == firstField && !widget.evaluationInProgress(),
            "Returning to a cached frequency did not restore its field.");
    quantity->setCurrentIndex(2);
    int edits = 0;
    widget.onChanged = [&] { ++edits; };
    const auto savedPhase = widget.settings().phaseDegrees;
    animate->click();
    timer.restart();
    while (publishedAppearance.phaseDegrees == savedPhase && timer.elapsed() < 5000)
        QApplication::processEvents();
    require(publishedAppearance.phaseDegrees != savedPhase && edits == 0 &&
                widget.settings().phaseDegrees == savedPhase && widget.fieldPlane() == firstField,
            "Animation changed stored settings or recomputed the field.");
    animate->click();
    require(publishedAppearance.phaseDegrees == savedPhase && edits == 0,
            "Pausing animation did not restore the saved phase.");
    widget.show();
    QApplication::processEvents();
    animate->click();
    auto *details = widget.findChild<QToolButton *>("fieldDetailsToggle");
    auto *dialog = widget.findChild<QDialog *>("fieldDetailsDialog");
    require(details && dialog, "Field details controls missing.");
    details->click();
    QApplication::processEvents();
    require(dialog->isVisible() && animate->text().contains("pausieren"),
            "Opening field settings stopped animation.");
    auto *probeMode = widget.findChild<QToolButton *>("fieldProbeMode");
    bool probing = false;
    widget.onProbeModeChanged = [&](bool enabled) { probing = enabled; };
    require(probeMode, "Field probe control missing.");
    probeMode->setChecked(true);
    require(probing, "Explicit field probe mode did not activate.");
    dialog->close();
    widget.hide();
    QApplication::processEvents();
    require(animate->text().contains("starten") && publishedAppearance.phaseDegrees == savedPhase &&
                edits == 0 && !probing && !probeMode->isChecked(),
            "Hidden field continued animating or retained probe mode.");
    calculate->click();
    frequency->setCurrentIndex(1);
    excitation->setCurrentIndex(1);
    calculate->click();
    timer.restart();
    while (!widget.fieldPlane() && timer.elapsed() < 5000) QApplication::processEvents();
    require(widget.fieldPlane() && widget.fieldPlane()->frequency == 700 &&
            widget.fieldPlane()->excitation == "Second source" &&
            widget.settings().frequencyHz == 700 && widget.settings().excitationIndex == 1,
            "Canceled field job replaced the latest selected result");
    auto movedFrame = AcousticFrame{};
    movedFrame.origin.x = .02;
    widget.setSnapshot(mesh, request, result, movedFrame, 9, true);
    require(!widget.fieldPlane(), "Changed frame retained the old field plane");
    widget.requestEvaluation(); timer.restart();
    while (!widget.fieldPlane() && timer.elapsed() < 5000) QApplication::processEvents();
    require(widget.fieldPlane() && widget.fieldPlane()->frame.origin.x == .02,
            "Changed frame did not reach the next field");
    auto changedRequest = request; changedRequest.medium.soundSpeed += 1;
    widget.setSnapshot(mesh, changedRequest, result, movedFrame, 9, true);
    require(!widget.fieldPlane(), "Changed request retained the old field plane");
    widget.setSnapshot(mesh, request, result, {}, 9, true);
    frequency->setCurrentIndex(0); excitation->setCurrentIndex(0);
    if (!screenshotDirectory.isEmpty()) {
        settings.lateralCount = 41;
        settings.axialCount = 51;
        settings.minimumLateral = -.2;
        settings.maximumLateral = .2;
        settings.minimumAxial = -.11;
        settings.maximumAxial = .35;
        widget.setSettings(settings);
        calculate->click();
        timer.restart();
        auto *cancel = widget.findChild<QPushButton *>("cancelField");
        while (cancel->isEnabled() && timer.elapsed() < 10000)
            QApplication::processEvents();
        require(widget.fieldPlane() && widget.fieldPlane()->pressure.size() == 41 * 51,
                "QA field evaluation did not finish.");
        screenshot(widget, "Feld-Steuerung-Pegel-QA.png");
        quantity->setCurrentIndex(1);
        screenshot(widget, "Feld-Steuerung-Phase-QA.png");
        quantity->setCurrentIndex(2);
        screenshot(widget, "Feld-Steuerung-Realteil-QA.png");
    }
    widget.requestEvaluation();
    require(widget.evaluationInProgress() && !widget.fieldReady(),
            "Explicit field start reused a previous displayed result as completion.");
    widget.cancel();
    QApplication::processEvents();
    require(
        !widget.evaluationInProgress() && !widget.fieldReady() && !widget.lastError().isEmpty(),
        "A cancelled field request was reported as successful or accepted a late worker result.");
    require(!published, "Cancellation left a current-looking field in the 3D viewport contract.");
    auto failed = std::make_shared<BemResult>(*result);
    auto failedBoundary = *failed->frequencies[0].excitations[0].boundary();
    failedBoundary.normalDerivative[0] = {std::numeric_limits<double>::quiet_NaN(), 0};
    failed->frequencies[0].excitations[0].setBoundary(std::move(failedBoundary));
    widget.setSnapshot(mesh, request, failed, {}, 9, true);
    widget.requestEvaluation();
    timer.restart();
    while (widget.evaluationInProgress() && timer.elapsed() < 5000)
        QApplication::processEvents();
    require(!widget.evaluationInProgress() && !widget.fieldReady() && !widget.lastError().isEmpty(),
            "A failed background field evaluation did not publish its error state.");
    widget.setSnapshot(mesh, request, result, {}, 9, false);
    require(!widget.fieldPlane(), "An old field stayed unmarked after result invalidation.");
    auto incomplete = std::make_shared<BemResult>(*result);
    auto incompleteBoundary = *incomplete->frequencies[0].excitations[0].boundary();
    incompleteBoundary.pressureTrace.clear();
    incomplete->frequencies[0].excitations[0].setBoundary(std::move(incompleteBoundary));
    widget.setSnapshot(mesh, request, incomplete, {}, 9, false);
    require(!calculate->isEnabled(),
            "An old archive without complete traces can fabricate a pressure field.");
    widget.requestEvaluation();
    require(!widget.evaluationInProgress() && !widget.fieldReady() && !widget.lastError().isEmpty(),
            "Explicit field request fabricated missing saved traces.");
}
} // namespace
int main(int argc, char **argv) {
    QApplication app(argc, argv);
    QTemporaryDir settings;
    QCoreApplication::setOrganizationName("HornSculptTests");
    QCoreApplication::setApplicationName("FieldAndPhase");
    QSettings::setDefaultFormat(QSettings::IniFormat);
    QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, settings.path());
    const auto arguments = QCoreApplication::arguments();
    const int screenshotOption = arguments.indexOf("--screenshots");
    if (screenshotOption >= 0 && screenshotOption + 1 < arguments.size()) {
        screenshotDirectory = arguments[screenshotOption + 1];
        QDir().mkpath(screenshotDirectory);
    }
    try {
        phaseControls();
        planeContract();
        asynchronousField();
        std::cout << "Phase controls and asynchronous pressure field passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
