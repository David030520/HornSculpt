#include "acoustic_partition_fixture.hpp"
#include "angular_response_widget.hpp"
#include "baffle_dialog.hpp"
#include "cad_fixture.hpp"
#include "fixture_geometry.hpp"
#include "freecad_dock.hpp"
#include "frequency_schedule.hpp"
#include "gmsh_source.hpp"
#include "latest_publication.hpp"
#include "lem_editor.hpp"
#include "main_window.hpp"
#include "mesh_panel.hpp"
#include "parameter_panel.hpp"
#include "project_assembly.hpp"
#include "simulation_dock.hpp"
#include "simulation_layout_checks.hpp"
#include "simulation_request.hpp"
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QCloseEvent>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QElapsedTimer>
#include <QEvent>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QHostAddress>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLabel>
#include <QLineEdit>
#include <QMainWindow>
#include <QMessageBox>
#include <QMouseEvent>
#include <QPixmap>
#include <QProgressBar>
#include <QProgressDialog>
#include <QPushButton>
#include <QSemaphore>
#include <QSettings>
#include <QSpinBox>
#include <QStatusBar>
#include <QTabBar>
#include <QTabWidget>
#include <QTcpSocket>
#include <QTemporaryDir>
#include <QTest>
#include <QTimer>
#include <QToolButton>
#include <QTreeWidget>
#include <QTreeWidgetItemIterator>
#include <QUndoStack>
#include <QWheelEvent>
#include <QtEndian>
#include <algorithm>
#include <array>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <hornsculpt/axisymmetric_geometry.hpp>
#include <hornsculpt/axisymmetric_mesh.hpp>
#include <hornsculpt/cad_geometry.hpp>
#include <hornsculpt/mesh_viewport.hpp>
#include <hornsculpt/plot_widget.hpp>
#include <hornsculpt/polar_plot_widget.hpp>
#include <iostream>
#include <limits>
#include <numbers>
#include <set>
#include <sstream>
#include <stdexcept>
#include <system_error>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class Predicate> void waitUntil(Predicate predicate) {
    QElapsedTimer timer;
    timer.start();
    while (!predicate()) {
        if (timer.elapsed() > 30000)
            throw std::runtime_error("Realtime worker timed out.");
        QTest::qWait(5);
    }
}
void saveAndWait(MainWindow &window, QAction *save) {
    require(save && save->isEnabled(), "Save action was unavailable.");
    save->trigger();
    require(window.saveInProgress(), "Save worker did not start.");
    waitUntil([&] { return !window.saveInProgress(); });
}
std::uint32_t clickViewportFace(MeshViewport &viewport) {
    for (int y = 1; y < 20; ++y)
        for (int x = 1; x < 20; ++x) {
            const QPointF position(viewport.width() * x / 20.0, viewport.height() * y / 20.0);
            if (const auto face = viewport.pickFace(position)) {
                QMouseEvent press(QEvent::MouseButtonPress, position,
                                  viewport.mapToGlobal(position.toPoint()), Qt::LeftButton,
                                  Qt::LeftButton, Qt::NoModifier);
                QMouseEvent release(QEvent::MouseButtonRelease, position,
                                    viewport.mapToGlobal(position.toPoint()), Qt::LeftButton,
                                    Qt::NoButton, Qt::NoModifier);
                QApplication::sendEvent(&viewport, &press);
                QApplication::sendEvent(&viewport, &release);
                return *face;
            }
        }
    throw std::runtime_error("No physical face can be picked in the CAD viewport.");
}
void testBaffleControls() {
    QMainWindow window;
    SimulationDock dock(&window);
    const auto defaults = dock.archive();
    const InfiniteBaffle baffle{{.00125, -.002, .003}, {0, 0, 1}, {3, 4}};
    auto *fmm = dock.findChild<QCheckBox *>("fastMultipoleEnabled");
    auto *button = dock.findChild<QPushButton *>("infiniteBaffle");
    require(fmm && button && !defaults.settings.infiniteBaffle,
            "Baffle mode is missing or changed the free-field default.");
    unsigned changes = 0;
    dock.onChanged = [&] { ++changes; };
    dock.setInfiniteBaffle(baffle);
    require(changes == 1 && dock.archive().settings.infiniteBaffle == baffle && !fmm->isEnabled(),
            "Baffle edits did not notify revisions, persist or exclude FMM.");
    const auto saved = dock.archive();
    dock.restoreArchive(defaults);
    require(fmm->isEnabled() && !dock.archive().settings.infiniteBaffle,
            "Restoring free field did not restore the supported solver controls.");
    fmm->setChecked(true);
    dock.setInfiniteBaffle(baffle);
    require(fmm->isChecked() && fmm->isEnabled() &&
                dock.archive().settings.acceleration == BemAcceleration::fastMultipole,
            "Enabling the baffle silently selected another solver.");
    fmm->setChecked(false);
    require(!fmm->isEnabled(), "FMM can be enabled for an unsupported baffle solve.");
    dock.restoreArchive(saved);
    require(dock.archive().settings.infiniteBaffle == baffle,
            "UI restoration lost the baffle domain.");

    BaffleDialog dialog(baffle);
    auto *mode = dialog.findChild<QComboBox *>("baffleMode");
    auto *faces = dialog.findChild<QLineEdit *>("baffleApertureFaces");
    auto *normal = dialog.findChild<QDoubleSpinBox *>("baffleNormalZ");
    auto *origin = dialog.findChild<QDoubleSpinBox *>("baffleOriginX");
    auto *buttons = dialog.findChild<QDialogButtonBox *>();
    require(mode && faces && normal && origin && buttons && mode->currentIndex() == 1 &&
                origin->value() == 1.25,
            "Baffle dialog did not restore its plane in millimetres.");
    faces->setText("3, 3");
    dialog.accept();
    require(dialog.result() != QDialog::Accepted &&
                !buttons->button(QDialogButtonBox::Ok)->isEnabled(),
            "The baffle dialog accepted duplicate aperture IDs.");
    faces->setText("3, 4");
    normal->setValue(0);
    dialog.accept();
    require(dialog.result() != QDialog::Accepted,
            "The baffle dialog accepted a zero plane normal.");
    normal->setValue(2);
    dialog.accept();
    require(dialog.result() == QDialog::Accepted && dialog.infiniteBaffle() == baffle,
            "Explicit normal normalization or millimetre conversion changed the domain.");
    BaffleDialog disabled(baffle);
    disabled.findChild<QComboBox *>("baffleMode")->setCurrentIndex(0);
    disabled.accept();
    require(disabled.result() == QDialog::Accepted && !disabled.infiniteBaffle(),
            "Freifeld did not remove the infinite-baffle domain.");
    std::cout << "CPU GUI: baffle plane, aperture validation, explicit solver mode and restoration "
                 "passed.\n";
}

std::shared_ptr<const DraftGeometry> rotationalDraft() {
    DraftParameters p;
    p.length = .02;
    p.entry.width = p.entry.height = .01;
    p.wallThickness = .002;
    p.horizontalAngle = p.verticalAngle = .3;
    p.verticalProfile = p.profile;
    return std::make_shared<DraftGeometry>(p);
}
struct RecessFixture {
    std::shared_ptr<const CadSnapshot> snapshot;
    std::uint32_t drivenFace{}, apertureFace{};
};
RecessFixture rotationalRecess() {
    // An actual exact air solid: bottom piston, cylindrical wall and planar
    // aperture. No quality mesh or customer artifact enters this fixture.
    const auto shape =
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -.005), gp_Dir(0, 0, 1)), .01, .005).Shape();
    std::ostringstream bytes;
    BRepTools::Write(shape, bytes, Standard_False, Standard_False, TopTools_FormatVersion_CURRENT);
    CadBrepPart part{bytes.str(), {1, 2, 3}};
    const auto geometry = CadGeometry::fromBrepParts(std::span{&part, 1}, 1);
    auto cad = std::make_shared<CadSnapshot>();
    cad->documentId = "axisymmetric-air-recess";
    cad->documentName = "Luftkavität";
    cad->sessionId = "offline-fixture";
    cad->revision = 1;
    cad->brepMetresPerUnit = 1;
    cad->brepObjects.push_back({"air", QByteArray::fromStdString(part.data), part.faceIds});
    cad->objects.push_back({"air",
                            "Air",
                            "Air cylinder",
                            "three exact surfaces",
                            {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}});
    RecessFixture fixture;
    for (const auto &measure : geometry->faces()) {
        CadFace face;
        face.id = measure.faceId;
        face.objectId = "air";
        face.sourceKey = QStringLiteral("native-%1").arg(face.id);
        face.label = face.sourceKey;
        face.state = CadFaceState::confirmed;
        face.signature = {QString::fromStdString(measure.surfaceType),
                          face.sourceKey,
                          face.sourceKey,
                          measure.area,
                          measure.center,
                          measure.normal,
                          measure.edgeCount,
                          measure.wireCount};
        if (measure.surfaceType == "Plane") {
            if (measure.center.z < -.0025) {
                fixture.drivenFace = face.id;
                face.boundary = BoundaryCondition::driven;
                face.label = "Kolben";
            } else {
                fixture.apertureFace = face.id;
                face.label = "Apertur";
            }
        }
        cad->faces.push_back(std::move(face));
    }
    require(fixture.drivenFace && fixture.apertureFace &&
                fixture.drivenFace != fixture.apertureFace,
            "Exact air cavity has no distinct piston and planar aperture.");
    prepareCadGeometry(*cad, true);
    fixture.snapshot = std::move(cad);
    return fixture;
}
void testAxisymmetricControls() {
    QMainWindow window;
    SimulationDock dock(&window);
    auto *mode = dock.findChild<QComboBox *>("bemGeometryMode");
    auto *maximum = dock.findChild<QDoubleSpinBox *>("axisymmetricMaximumLength");
    auto *tolerance = dock.findChild<QDoubleSpinBox *>("axisymmetricChordTolerance");
    auto *symmetry = dock.findChild<QDoubleSpinBox *>("axisymmetricRotationalTolerance");
    auto *sectors = dock.findChild<QSpinBox *>("axisymmetricSectors");
    auto *fmm = dock.findChild<QCheckBox *>("fastMultipoleEnabled");
    QPushButton *solve = nullptr;
    for (auto *button : dock.findChildren<QPushButton *>())
        if (button->text() == "Berechnen")
            solve = button;
    require(mode && maximum && tolerance && symmetry && sectors && fmm && solve &&
                mode->currentIndex() == 0 && !maximum->isEnabled(),
            "2D controls changed the default 3D mode or are missing.");
    const auto source = rotationalRecess();
    const std::vector<FaceVelocity> velocities{{source.drivenFace, {1, 0}}};
    dock.setMesh(nullptr, 9, velocities);
    dock.setGeometrySource(nullptr, source.snapshot, 9);
    require(!solve->isEnabled(), "3D can solve without its quality mesh.");
    mode->setCurrentIndex(1);
    require(solve->isEnabled() && maximum->isEnabled() && !fmm->isEnabled(),
            "2D still requires a 3D quality mesh or offers unsupported FMM.");
    dock.requestSolve();
    require(!dock.isBusy() && dock.lastError().contains("Schallwand"),
            "Reduced free-field mode was dispatched instead of requiring the air-baffle domain.");
    dock.setInfiniteBaffle(InfiniteBaffle{{}, {0, 0, 1}, {source.apertureFace}});
    maximum->setValue(4);
    tolerance->setValue(.075);
    symmetry->setValue(.02);
    sectors->setValue(16);
    dock.setFrequencySweep(300, 300, 1);
    auto saved = dock.archive();
    require(saved.settings.geometryMode == BemGeometryMode::axisymmetric &&
                saved.settings.axisymmetricSampling ==
                    AxisymmetricSamplingSettings{.004, .000075, 2e-8, 16},
            "2D GUI units do not round-trip to SI.");
    sectors->setValue(10);
    dock.requestSolve();
    require(!dock.isBusy() && !dock.lastError().isEmpty(),
            "Invalid angular sectors reached the numerical worker.");
    dock.restoreArchive(saved);
    MeridianSampling sampling;
    sampling.maximumLength = saved.settings.axisymmetricSampling.maximumLength;
    sampling.chordTolerance = saved.settings.axisymmetricSampling.chordTolerance;
    sampling.rotationalTolerance = saved.settings.axisymmetricSampling.rotationalTolerance;
    for (const auto &face : source.snapshot->faces)
        sampling.faceConditions.emplace(face.id, face.boundary);
    const auto meridian =
        makeAxisymmetricMeridian(*source.snapshot->geometry, {}, {0, 0, 1}, sampling);
    auto generated = buildAxisymmetricMesh(meridian, 16);
    saved.mesh = std::make_shared<SurfaceMesh>(generated.mesh);
    saved.request = makeSimulationRequest(saved.settings, velocities);
    saved.request->axisymmetric =
        AxisymmetricBem{meridian.origin,        meridian.axis,      16,
                        generated.pressureDofs, generated.fluxDofs, generated.testPanels};
    saved.axisymmetricSampling = saved.settings.axisymmetricSampling;
    auto result = std::make_shared<BemResult>();
    result->device = "SYNTHETIC GUI SNAPSHOT - no GPU computation";
    result->geometryMode = BemGeometryMode::axisymmetric;
    result->pressureDofCount = meridian.points.size();
    result->triangleCount = generated.mesh.triangles.size();
    BemExcitationResult output;
    output.name = saved.request->excitations.front().name;
    BemBoundaryData boundary;
    boundary.surfacePressure.assign(result->triangleCount, {.01, 0});
    output.setBoundary(std::move(boundary));
    output.observerPressure.assign(saved.request->observers.size(), {.01, 0});
    result->frequencies.push_back({300, {std::move(output)}});
    saved.result = result;
    saved.current = true;
    saved.revision = 9;
    dock.restoreArchive(saved);
    require(dock.hasCurrentResult() &&
                dock.archive().axisymmetricSampling == saved.axisymmetricSampling,
            "A reduced snapshot without a 3D quality mesh cannot be restored as current.");
    auto legacy = saved;
    legacy.request->interactionPlanFrequencies.clear();
    dock.restoreArchive(legacy);
    require(dock.hasCurrentResult(),
            "An absent legacy planning raster invalidated the same full sweep.");
    dock.restoreArchive(saved);
    maximum->setValue(5);
    require(!dock.hasCurrentResult() && dock.result() == result &&
                dock.archive().axisymmetricSampling == saved.axisymmetricSampling,
            "Changing 2D discretization failed to invalidate or rewrote the solved snapshot.");
    dock.restoreArchive(saved);
    fmm->setChecked(true);
    require(fmm->isChecked() && fmm->isEnabled(), "An unsupported solver was silently replaced.");
    dock.requestSolve();
    require(!dock.isBusy() && !dock.lastError().isEmpty(), "2D FMM was silently dispatched.");
    dock.restoreArchive(saved);
    mode->setCurrentIndex(0);
    require(!dock.hasCurrentResult() && !solve->isEnabled() && !maximum->isEnabled(),
            "The 3D/2D switch reused the other geometry mode's current result.");
    dock.restoreArchive(saved);
    auto p = rotationalDraft()->parameters();
    p.verticalAngle *= .8;
    dock.setGeometrySource(std::make_shared<DraftGeometry>(p), nullptr, 9);
    require(!dock.hasCurrentResult(), "Changing the immutable 2D source kept a current result.");
    auto *workspace = dock.findChild<QTabWidget *>("simulationWorkspace");
    auto *resultTabs = dock.findChild<QTabWidget *>("acousticResults");
    auto *field = dock.findChild<QWidget *>("acousticField");
    auto *commands = dock.findChild<QWidget *>("solverCommands");
    auto *stop = dock.findChild<QPushButton *>("stopSimulation");
    require(workspace && resultTabs && field && commands && stop,
            "Field and solver commands are missing.");
    int fieldIndex = -1;
    for (int i = 0; i < resultTabs->count(); ++i)
        if (resultTabs->widget(i)->isAncestorOf(field))
            fieldIndex = i;
    require(fieldIndex >= 0, "Field result page is missing.");
    resultTabs->setCurrentIndex(fieldIndex);
    dock.revealResultView();
    require(commands->isHidden(), "Idle field view retained the redundant solver commands.");
    workspace->setCurrentIndex(0);
    require(!commands->isHidden(), "Solver commands disappeared from the calculation page.");
    dock.revealResultView();
    resultTabs->setCurrentIndex(0);
    require(!commands->isHidden(), "Solver commands disappeared from other result pages.");
    resultTabs->setCurrentIndex(fieldIndex);
    dock.setRealtimeEnabled(true);
    require(!commands->isHidden(), "Enabled realtime cannot be switched off in the field view.");
    dock.setRealtimeEnabled(false);
    require(commands->isHidden(), "Stopped realtime retained the redundant command row.");
    dock.requestSolve();
    require(dock.isBusy() && !commands->isHidden() && stop->isEnabled(),
            "Active solver cannot be cancelled from the field view.");
    waitUntil([&] { return !dock.isBusy(); });
    require(commands->isHidden(), "Completed solver retained the redundant command row.");
    require(!dock.lastError().isEmpty() && dock.result() == result && !dock.hasCurrentResult(),
            "An asymmetric source was accepted or replaced the old visible result.");
    std::cout << "CPU GUI: independent 2D mode, SI sampling, mesh-free readiness, immutable "
                 "snapshot restore, invalidation and asymmetric-source rejection passed.\n";
}
void testAxisymmetricSimulation() {
    QTemporaryDir directory;
    require(directory.isValid(), "2D project directory is unavailable.");
    const auto source = rotationalRecess();
    Project project;
    project.source = ProjectSource::freecad;
    project.revision = project.geometryRevision = 1;
    project.cad = source.snapshot;
    project.preview = std::make_shared<SurfaceMesh>(source.snapshot->mesh);
    auto &settings = project.simulation.settings;
    settings.geometryMode = BemGeometryMode::axisymmetric;
    settings.axisymmetricSampling = {.01, .001, 1e-8, 8};
    settings.minimumFrequency = settings.maximumFrequency = 300;
    settings.frequencyCount = 1;
    settings.infiniteBaffle = InfiniteBaffle{{}, {0, 0, 1}, {source.apertureFace}};
    settings.observationFrame = AcousticFrame{{.01, -.02, .03}, {1, 0, 0}, {0, 1, 0}};
    const auto path = directory.filePath("rotational.hsproj");
    writeProject(path, project);
    MainWindow main;
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() && main.geometryReady() && !main.meshReady(),
            "The 2D project required a 3D quality mesh on import.");
    auto *dock = dynamic_cast<SimulationDock *>(main.findChild<QDockWidget *>("simulationDock"));
    require(dock, "2D simulation dock is missing.");
    main.requestSolve();
    waitUntil([&] { return !main.simulationInProgress(); });
    if (!main.lastSimulationError().isEmpty())
        throw std::runtime_error(main.lastSimulationError().toStdString());
    const auto computed = dock->archive();
    require(main.simulationReady() && !main.meshReady() && computed.request->axisymmetric &&
                computed.result->geometryMode == BemGeometryMode::axisymmetric &&
                computed.request->axisymmetric->origin == Vec3{} &&
                computed.request->axisymmetric->axis == Vec3{0, 0, 1} &&
                computed.request->infiniteBaffle->normal == Vec3{0, 0, 1} &&
                computed.request->infiniteBaffle->origin == Vec3{} &&
                computed.request->farField->origin == settings.observationFrame->origin &&
                computed.result->pressureDofCount < computed.mesh->vertices.size() &&
                computed.result->frequencies[0].excitations[0].boundary()->surfacePressure.size() ==
                    computed.mesh->triangles.size(),
            "The GUI did not run a real reduced system independently of observer rotation.");
    project.revision = project.geometryRevision = main.revision();
    project.simulation = computed;
    writeProject(path, project);
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() && main.simulationReady() && !main.meshReady(),
            "A computed 2D project lost current state after save/open without a 3D mesh.");
    dock->setFrequencySweep(300, 1200, 3);
    const std::vector<double> raster{300, 600, 1200};
    bool validPartials = true;
    unsigned partials = 0;
    const auto previousChanged = dock->onChanged;
    dock->onChanged = [&] {
        if (previousChanged)
            previousChanged();
        const auto partial = dock->archive();
        if (dock->isBusy() && partial.request &&
            partial.request->interactionPlanFrequencies == raster && partial.result &&
            partial.result->frequencies.size() < raster.size()) {
            ++partials;
            validPartials &= !partial.current && partial.request->frequencies.size() ==
                                                     partial.result->frequencies.size();
        }
    };
    dock->setRealtimeEnabled(true);
    waitUntil([&] {
        return dock->hasCurrentResult() && dock->archive().result->frequencies.size() == 3;
    });
    dock->onChanged = previousChanged;
    const auto realtime = dock->archive();
    unsigned hits = 0, builds = 0;
    for (const auto &f : realtime.result->frequencies) {
        hits += f.interactionPlanReused ? 1 : 0;
        builds += f.interactionPlanMilliseconds > 0 ? 1 : 0;
    }
    require(
        partials > 0 && validPartials && realtime.request->interactionPlanFrequencies == raster &&
            hits == 2 && builds == 1,
        "Realtime subsets lost the complete raster or rebuilt the cold host plan per frequency.");
    const auto restoredResult = dock->result();
    dock->invalidate(main.revision() + 1);
    require(!dock->hasCurrentResult() && dock->result() == restoredResult,
            "A newer geometry failed to stale the reduced result.");
    std::cout << "GPU GUI: actual mesh-free m=0 generation/solve, expanded fields and "
                 "independent observation frame passed.\n";
}

void testInterfaceControls() {
    QMainWindow host;
    SimulationDock dock(&host);
    const lem::Terminals pair{{1}, {2}};
    lem::Network network;
    network.nodes = {{{1}, lem::Domain::acoustic, false, "Chamber", {}},
                     {{2}, lem::Domain::acoustic, true, "Reference", {}}};
    network.elements = {
        {{1}, "Interface", {}, lem::BemPort{pair, {1}, 1, lem::BemPortKind::interfacePort}},
        {{2}, "Unspecified chamber", {}, lem::Cavity{pair, 0}}};
    auto knownMesh = std::make_shared<SurfaceMesh>(testing::sphere(0));
    for (std::size_t i = 0; i < knownMesh->triangles.size(); ++i)
        knownMesh->triangles[i].faceId = i % 2 ? 3 : 2;
    dock.setMesh(knownMesh, 7, {{2, {1, 0}}, {3, {.5, 0}}});
    dock.setInfiniteBaffle(InfiniteBaffle{{}, {0, 0, 1}, {4}});
    dock.setCoupling(network, true);
    auto *scale = dock.findChild<QDoubleSpinBox *>("prescribedVelocityScale");
    auto *phase = dock.findChild<QDoubleSpinBox *>("prescribedVelocityPhase");
    require(scale && phase && scale->isEnabled() && phase->isEnabled(),
            "Interface coupling disabled the independent dome/surround excitation controls");
    scale->setValue(.5);
    phase->setValue(31);
    const auto archive = dock.archive();
    dock.setInfiniteBaffle(std::nullopt);
    dock.setCoupling(
        lem::makeDriverNetwork({6.4, .0007, 4.98, .004, .0016, .262, .00544}, {{1}}, {1, 0}, .0004),
        true);
    require(!scale->isEnabled() && !phase->isEnabled(),
            "Driven-port coupling exposes an unused prescribed excitation");
    dock.restoreArchive(archive);
    require(scale->isEnabled() && phase->isEnabled() && scale->value() == .5 &&
                phase->value() == 31,
            "Restoring interface settings lost the known excitation controls");
    auto portBound = network;
    std::get<lem::Cavity>(portBound.elements[1].component).volume = 5e-6;
    portBound.elements.push_back({{3}, "Network drive", {}, lem::BemPort{pair, {2, 3}, 1}});
    std::get<lem::BemPort>(portBound.elements.back().component).velocityWeights = {{3, .5}};
    dock.setCoupling(portBound, true);
    require(!scale->isEnabled() && !phase->isEnabled(),
            "A fully port-driven model exposes unused known-velocity controls");
    auto saved = dock.archive();
    saved.request =
        makeSimulationRequest(saved.settings, std::vector<FaceVelocity>{{2, {1, 0}}, {3, {.5, 0}}});
    saved.mesh = knownMesh;
    saved.revision = 7;
    saved.current = true;
    // This fixture only checks immutable GUI snapshot state, not acoustic physics.
    auto synthetic = std::make_shared<BemResult>();
    for (const auto frequency : saved.request->frequencies) {
        BemExcitationResult output;
        output.name = "LEM";
        output.observerPressure.resize(saved.request->observers.size());
        synthetic->frequencies.push_back({frequency, {std::move(output)}});
    }
    saved.result = synthetic;
    dock.restoreArchive(saved);
    require(dock.hasCurrentResult(),
            "Weighted network fixture did not restore a current GUI snapshot");
    auto archived = saved;
    archived.current = false;
    archived.resultComplete = false;
    dock.restoreArchive(archived);
    require(!dock.hasCurrentResult() && dock.windowTitle().contains("Teilraster") &&
                dock.archive().resultComplete == std::optional{false},
            "Restoring a partial raster lost its completeness state.");
    archived.resultComplete.reset();
    dock.restoreArchive(archived);
    require(!dock.hasCurrentResult() && dock.windowTitle().contains("Vollständigkeit unbekannt") &&
                !dock.archive().resultComplete,
            "Legacy stale data was incorrectly labeled as complete or partial.");
    archived.resultComplete = true;
    dock.restoreArchive(archived);
    require(!dock.hasCurrentResult() && !dock.windowTitle().contains("Teilraster") &&
                !dock.windowTitle().contains("Vollständigkeit unbekannt"),
            "A known complete stale raster was mislabeled partial or unknown.");
    dock.restoreArchive(saved);
    auto changedWeights = portBound;
    std::get<lem::BemPort>(changedWeights.elements.back().component).velocityWeights[3] = .75;
    dock.setCoupling(changedWeights, true);
    require(!dock.hasCurrentResult() && dock.result() == synthetic &&
                dock.archive().request->coupling == portBound,
            "A weight change kept results current or mutated the immutable solve request");

    QTemporaryDir directory;
    require(directory.isValid(), "Interface view fixture directory unavailable");
    auto cad = std::make_shared<CadSnapshot>(testing::exactCadFixture());
    Project project;
    project.source = ProjectSource::freecad;
    project.cad = cad;
    project.preview = std::shared_ptr<const SurfaceMesh>(cad, &cad->mesh);
    project.geometryRevision = project.revision;
    const auto path = directory.filePath("interface-view.hsproj");
    writeProject(path, project);
    MainWindow main;
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    main.show();
    QApplication::processEvents();
    auto *editor = dynamic_cast<LemEditor *>(main.findChild<QWidget *>("lemEditor"));
    auto *visibility = main.findChild<QCheckBox *>("showBemLemInterfaces");
    auto *viewport = dynamic_cast<MeshViewport *>(main.findChild<QOpenGLWidget *>("meshViewport"));
    auto *cadDock = main.findChild<QDockWidget *>("freecadDock");
    auto *faces = cadDock ? cadDock->findChild<QTreeWidget *>("cadFaces") : nullptr;
    require(main.sourceImportError().isEmpty() && editor && visibility && viewport && faces &&
                !visibility->isEnabled() && !visibility->isChecked(),
            "Projects without virtual interfaces expose an active interface view switch");
    const auto visibleFace = clickViewportFace(*viewport);
    std::get<lem::BemPort>(network.elements[0].component).faceIds = {visibleFace};
    editor->setNetwork(network);
    editor->onChanged();
    require(visibility->isEnabled() && !visibility->isChecked(),
            "Assigned mathematical interfaces are shown as material walls by default");
    const auto canPick = [&] {
        for (int y = 1; y < 20; ++y)
            for (int x = 1; x < 20; ++x)
                if (viewport->pickFace(
                        {viewport->width() * x / 20., viewport->height() * y / 20.}) == visibleFace)
                    return true;
        return false;
    };
    require(!canPick(), "A hidden mathematical interface remains pickable as a material surface");
    const auto revision = main.revision();
    visibility->setChecked(true);
    require(canPick() && main.revision() == revision && editor->network() == network &&
                !main.simulationInProgress(),
            "Interface visibility changed physics or failed to restore 3D selection");
    bool labelled = false;
    for (int row = 0; row < faces->topLevelItemCount(); ++row)
        if (faces->topLevelItem(row)->data(0, Qt::UserRole).toUInt() == visibleFace)
            labelled = faces->topLevelItem(row)->text(1).contains("Interface");
    require(labelled, "The CAD face table presents a mathematical interface as ordinary Wall");
    std::cout << "CPU GUI: interface role, known-drive controls, archive restore and "
                 "visibility-only picking passed.\n";
}
void testAngularResponseControls() {
    QMainWindow window;
    SimulationDock dock(&window);
    auto mesh = std::make_shared<SurfaceMesh>(testing::sphere(0));
    dock.setMesh(mesh, 7);
    dock.setFrequencySweep(100, 1000, 3);
    auto saved = dock.archive();
    const std::vector<FaceVelocity> velocities{{1, {1, 0}}};
    saved.request = makeSimulationRequest(saved.settings, velocities);
    saved.mesh = mesh;
    auto result = std::make_shared<BemResult>();
    for (const auto frequency : saved.request->frequencies) {
        BemExcitationResult excitation;
        excitation.observerPressure.assign(147, {.01, 0});
        excitation.observerPressure[0] = {.02, 0};
        // Distinct duplicated axis samples exercise the dedicated on-axis index.
        excitation.observerPressure[37] = {.005, 0};
        excitation.observerPressure[110] = {.0025, 0};
        excitation.name = "Source A";
        auto second = excitation;
        second.name = "Source B";
        for (auto &pressure : second.observerPressure)
            pressure *= 2.;
        result->frequencies.push_back({frequency, {std::move(excitation), std::move(second)}});
    }
    saved.result = result;
    saved.revision = 7;
    saved.current = true;
    dock.restoreArchive(saved);
    require(dock.hasCurrentResult(), "Angular test archive failed current-result restoration.");
    auto *progress = dock.findChild<QProgressBar *>("simulationProgress");
    require(progress && progress->isHidden(),
            "Restored result displayed an inactive 0% progress bar.");
    auto *responses =
        dynamic_cast<AngularResponseWidget *>(dock.findChild<QWidget *>("angularResponses"));
    auto *mode = dock.findChild<QComboBox *>("angularResponseNormalization");
    auto *reference = dock.findChild<QDoubleSpinBox *>("angularResponseReferenceAngle");
    auto *h = dynamic_cast<PlotWidget *>(dock.findChild<QWidget *>("angularFrequencyH"));
    auto *v = dynamic_cast<PlotWidget *>(dock.findChild<QWidget *>("angularFrequencyV"));
    require(responses && mode && reference && h && v, "Angular results are not in the dock.");
    unsigned edits = 0;
    dock.onChanged = [&] { ++edits; };
    mode->setCurrentIndex(1);
    require(std::abs(h->series()[2].points[0].y() + 6.020599913279624) < 1e-10 &&
                h->series()[0].points[0].y() == 0 && v->series()[0].points[0].y() == 0,
            "Angular normalization did not use the actual dedicated on-axis pressure.");
    reference->setValue(30);
    require(std::abs(h->series()[0].points[0].y() - 6.020599913279624) < 1e-10 &&
                h->series()[2].points[0].y() == 0,
            "Editable angular reference did not update the dock curves.");
    const auto changed = dock.archive();
    require(edits == 2 && dock.hasCurrentResult() && !dock.isBusy() && dock.result() == result &&
                changed.request->observers == saved.request->observers &&
                changed.settings.angularResponses == responses->settings(),
            "Display-only edits invalidated results, started a solve or lost settings.");
    dock.restoreArchive(changed);
    require(dock.hasCurrentResult() && responses->settings().referenceAngleDegrees == 30,
            "Angular display settings did not survive dock restoration.");
    auto *tabs = dock.findChild<QTabWidget *>("acousticResults");
    auto *frequencyControls = dock.findChild<QWidget *>("polarFrequencyControls");
    auto *normalizationControls = dock.findChild<QWidget *>("polarNormalizationControls");
    auto *scaleControls = dock.findChild<QWidget *>("directivityScaleControls");
    require(tabs && frequencyControls && normalizationControls && scaleControls,
            "Result controls missing");
    auto *family = dock.findChild<QComboBox *>("resultFamily");
    auto *view = dock.findChild<QComboBox *>("resultView");
    require(family && view && family->count() == 3, "Result family navigation is missing.");
    int navigablePages = 0;
    for (int group = 0; group < family->count(); ++group) {
        family->setCurrentIndex(group);
        for (int i = 0; i < view->count(); ++i) {
            view->setCurrentIndex(i);
            require(tabs->currentIndex() == view->currentData().toInt(),
                    "Result navigation selected a different page.");
            ++navigablePages;
        }
    }
    require(navigablePages == tabs->count(), "Not every result page is reachable.");
    const auto showPageContaining = [&](const char *name) {
        auto *content = dock.findChild<QWidget *>(name);
        require(content, "Result page content missing");
        for (int i = 0; i < tabs->count(); ++i)
            if (tabs->widget(i)->isAncestorOf(content)) {
                tabs->setCurrentIndex(i);
                return;
            }
        throw std::runtime_error("Result page not found");
    };
    showPageContaining("cta2034Levels");
    require(family->currentIndex() == 0 && view->currentData().toInt() == tabs->currentIndex(),
            "Programmatic result navigation did not synchronize the selectors.");
    require(frequencyControls->isHidden() && normalizationControls->isHidden() &&
                scaleControls->isHidden(),
            "CTA page displayed polar controls");
    showPageContaining("acousticPolarH");
    require(!frequencyControls->isHidden() && !normalizationControls->isHidden() &&
                !scaleControls->isHidden(),
            "Polar page lost frequency or scale controls");
    showPageContaining("acousticIsobarH");
    require(frequencyControls->isHidden() && !normalizationControls->isHidden() &&
                !scaleControls->isHidden(),
            "Isobar page lost normalization or scale controls");
    showPageContaining("acousticRadialIsobarH");
    require(!normalizationControls->isHidden() && !scaleControls->isHidden(),
            "Radial isobar page lost normalization or scale controls");
    auto *ctaSelection = dock.findChild<QComboBox *>("ctaExcitation");
    auto *cta = dynamic_cast<PlotWidget *>(dock.findChild<QWidget *>("cta2034Levels"));
    auto *quality = dock.findChild<QToolButton *>("resultQuality");
    require(quality && quality->toolTip().contains("Konvergenz", Qt::CaseInsensitive) &&
                quality->toolTip().contains("70-Richtungs-Näherung"),
            "Result quality did not disclose its convergence and CTA limits.");
    require(ctaSelection && cta && ctaSelection->count() == 2 &&
                dock.findChild<QPushButton *>("exportCta2034"),
            "CTA controls missing");
    const double firstLevel = cta->series().front().points.front().y();
    ctaSelection->setCurrentIndex(1);
    require(std::abs(cta->series().front().points.front().y() - firstLevel - 20 * std::log10(2.)) <
                    1e-9 &&
                dock.hasCurrentResult() && !dock.isBusy(),
            "CTA selection changed physics or wrong source");
    const auto ctaSaved = dock.archive();
    require(ctaSaved.settings.ctaExcitationIndex == 1, "CTA selection not saved");
    dock.restoreArchive(ctaSaved);
    require(ctaSelection->currentIndex() == 1 && dock.hasCurrentResult(),
            "CTA selection not restored with current result");
    const auto previousResult = dock.result();
    dock.setFrequencySweep(20, 20000, static_cast<std::size_t>(std::numeric_limits<int>::max()));
    dock.requestSolve();
    require(!dock.isBusy() && dock.result() == previousResult &&
                dock.lastError().contains("CPU-Ergebnisdaten"),
            "Oversized frequency request allocated or discarded existing results instead of "
            "preflight rejection.");
    dock.invalidate(8);
    require(!dock.hasCurrentResult() && !h->series().empty(),
            "Geometry invalidation discarded the visible old angular traces.");
    require(quality->text().contains("veraltet") && quality->toolTip().contains("veraltet"),
            "Result quality did not follow geometry invalidation.");
    std::cout << "CPU GUI: actual pressure extraction, angular normalization and display-only "
                 "archive restoration passed.\n";
}
void testLemIntegration() {
    QTemporaryDir directory;
    require(directory.isValid(), "LEM integration directory is unavailable.");
    Project project;
    project.parameters.entry.width = .004;
    project.parameters.entry.height = .004;
    project.parameters.length = .02;
    project.parameters.wallThickness = .001;
    project.parameters.horizontalAngle = 20 * std::numbers::pi / 180;
    project.parameters.verticalAngle = 15 * std::numbers::pi / 180;
    project.simulation.settings.minimumFrequency = 300;
    project.simulation.settings.maximumFrequency = 300;
    project.simulation.settings.frequencyCount = 1;
    project.simulation.settings.farField = false;
    project.simulation.settings.observationDistance = 2;
    project.simulation.settings.couplingEnabled = true;
    project.simulation.settings.couplingNetwork =
        lem::makeDriverNetwork({6.4, .0007, 4.98, .004, .001647, .262, .00544},
                               {{static_cast<std::uint32_t>(DraftFace::driven)}}, {1, 0}, .000344);
    const auto path = directory.filePath("coupled.hsproj");
    writeProject(path, project);
    MainWindow main;
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() && main.geometryReady(),
            "Coupled project could not be opened.");
    auto *editor = dynamic_cast<LemEditor *>(main.findChild<QWidget *>("lemEditor"));
    auto *dock = dynamic_cast<SimulationDock *>(main.findChild<QDockWidget *>("simulationDock"));
    auto *undo = main.findChild<QAction *>("undoEdit");
    auto *redo = main.findChild<QAction *>("redoEdit");
    require(editor && dock && undo && redo && editor->couplingEnabled() &&
                editor->network() == project.simulation.settings.couplingNetwork,
            "Project LEM settings were not restored into the editor.");
    main.requestMesh();
    waitUntil([&] { return !main.meshInProgress(); });
    require(main.meshReady(), "Coupled Draft mesh failed.");
    main.requestSolve();
    waitUntil([&] { return !main.simulationInProgress(); });
    require(main.simulationReady(), "Coupled GPU solve failed.");
    const auto first = dock->archive();
    require(first.request && first.request->coupling && !first.request->farField &&
                first.request->excitations.empty() &&
                first.result->frequencies.front().excitations.front().lemSolution.has_value(),
            "Coupled request was replaced by prescribed velocities or far field.");
    const auto pressure = first.result->frequencies.front().excitations.front().observerPressure[0];
    const auto revision = main.revision();
    auto layoutOnly = editor->network();
    layoutOnly.elements.front().position.x += 100;
    layoutOnly.elements.front().name = "Renamed source";
    editor->setNetwork(layoutOnly);
    editor->onChanged();
    require(main.simulationReady() && dock->result() == first.result && main.revision() == revision,
            "A diagram layout/name edit invalidated unchanged physical results.");
    main.exerciseLoadedProjectRoundTrip();
    require(main.simulationReady(), "Layout-only save/open invalidated the physical snapshot.");
    auto doubled = editor->network();
    std::get<lem::EffortSource>(doubled.elements.front().component).amplitude *= 2;
    editor->setNetwork(doubled);
    editor->onChanged();
    require(!main.simulationReady() && main.meshReady() && main.revision() == revision,
            "Network edit did not invalidate results or unnecessarily invalidated the mesh.");
    undo->trigger();
    require(editor->network() == layoutOnly && main.meshReady(),
            "Undo failed to restore the circuit without remeshing.");
    redo->trigger();
    require(editor->network() == doubled && dock->archive().settings.couplingNetwork == doubled,
            "Redo failed to restore the actual solver circuit.");
    main.requestSolve();
    waitUntil([&] { return !main.simulationInProgress(); });
    require(main.simulationReady(), "Edited coupled GPU solve failed.");
    const auto second = dock->archive();
    const auto secondPressure =
        second.result->frequencies.front().excitations.front().observerPressure[0];
    require(std::abs(secondPressure - 2.0 * pressure) <= 1e-4 * std::abs(pressure),
            "Source voltage edit did not affect the actual BEM field.");
    main.exerciseLoadedProjectRoundTrip();
    require(main.simulationReady() && editor->network() == doubled &&
                dock->archive().request->coupling == doubled,
            "Save/open lost the coupled snapshot or current-result state.");
    editor->setCouplingEnabled(false);
    editor->onEnabledChanged(false);
    main.requestSolve();
    waitUntil([&] { return !main.simulationInProgress(); });
    require(main.simulationReady(), "Fixed-velocity mode was not restored after disabling LEM.");
    const auto fixedResult = dock->result();
    auto inactive = editor->network();
    std::get<lem::EffortSource>(inactive.elements.front().component).amplitude *= 3;
    editor->setNetwork(inactive);
    editor->onChanged();
    require(main.simulationReady() && dock->result() == fixedResult,
            "Editing an inactive circuit invalidated fixed-velocity results.");
    std::cout << "Real CUDA LEM GUI: Draft, finite observer, circuit edit, stale state, undo/redo "
                 "without remeshing, voltage response and project/result roundtrip passed.\n";
}

void testOrientationEdits() {
    MainWindow main;
    waitUntil([&] { return main.geometryReady(); });
    auto *faces = main.findChild<QTreeWidget *>("meshFaces");
    auto *panel = faces ? dynamic_cast<MeshPanel *>(faces->parentWidget()) : nullptr;
    require(panel, "Shared mesh panel is missing.");
    auto *reverse = panel->findChild<QPushButton *>("meshReverseFaces");
    auto *undo = main.findChild<QAction *>("undoEdit");
    auto *redo = main.findChild<QAction *>("redoEdit");
    require(faces && reverse && undo && redo, "Orientation edit controls are missing.");
    const auto reversed = [&] {
        int count = 0;
        for (int i = 0; i < faces->topLevelItemCount(); ++i)
            count += faces->topLevelItem(i)->text(2) == QString::fromUtf8("−");
        return count;
    };
    panel->setSelectedFaces({1, 5});
    const auto initialRevision = main.revision();
    reverse->click();
    require(reversed() == 2 && main.revision() > initialRevision && !main.simulationReady(),
            "Selected normal edits do not invalidate the current solver state.");
    const auto editRevision = main.revision();
    undo->trigger();
    waitUntil([&] { return main.geometryReady(); });
    require(reversed() == 0 && main.revision() > editRevision,
            "Undo did not restore the original normal orientations.");
    const auto undoRevision = main.revision();
    redo->trigger();
    waitUntil([&] { return main.geometryReady(); });
    require(reversed() == 2 && main.revision() > undoRevision,
            "Redo did not restore the normal edits with a new revision.");
    main.requestMesh();
    waitUntil([&] { return !main.meshInProgress(); });
    require(!main.meshReady() && !main.simulationReady() &&
                main.lastMeshError().startsWith(QString::fromUtf8("Flächennormalen:")),
            "An inconsistently oriented boundary was released for BEM.");
    main.exerciseLoadedProjectRoundTrip();
    require(reversed() == 2 && !main.simulationInProgress(),
            "Save/open lost normal edits or started a solver job.");
    std::cout << "CPU GUI: normal multi-selection, revision invalidation, undo/redo and save/open "
                 "passed.\n";
}
void testDraftAssemblyEditing() {
    QTemporaryDir directory;
    const auto path = directory.filePath("assembly.hsproj");
    Project p;
    p.parameters.length = .04;
    p.parameters.asymmetricEnabled = true;
    p.parameters.leftAngle = .3;
    p.meshing.adaptive.geometricTolerance = .001;
    p.assembly.instances.front().position = {-.3, 0, 0};
    auto right = p.assembly.instances.front();
    right.id = 2;
    right.position.x = .3;
    right.mirrorHorizontal = true;
    right.mouthAnchored = true;
    p.assembly.instances.push_back(right);
    p.assembly.includeCad = true;
    p.cad = std::make_shared<CadSnapshot>(testing::exactCadFixture());
    const auto local = DraftGeometry(p.parameters).previewMesh(.001);
    p.preview = std::make_shared<SurfaceMesh>(assemblyPreview(local, p));
    p.geometryRevision = p.revision;
    writeProject(path, p);
    MainWindow main;
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() && main.geometryReady(),
            "Mixed Draft/CAD project did not open");
    QAction *save = nullptr, *undo = nullptr;
    for (auto *action : main.findChildren<QAction *>()) {
        if (action->shortcut() == QKeySequence::Save)
            save = action;
        if (action->objectName() == "undoEdit")
            undo = action;
    }
    auto *length = main.findChild<QDoubleSpinBox *>("length");
    auto *dock = dynamic_cast<SimulationDock *>(main.findChild<QDockWidget *>("simulationDock"));
    require(save && undo && length && dock, "Assembly editing controls are missing");
    const auto before = main.revision();
    length->setValue(60);
    waitUntil([&] { return main.geometryReady() && main.revision() > before; });
    saveAndWait(main, save);
    const auto changed = readProject(path);
    require(changed.parameters.length == .06 && changed.assembly == p.assembly && changed.cad &&
                changed.cad->mesh == p.cad->mesh && !main.meshReady() && !dock->hasCurrentResult(),
            "Shared edit lost placement/CAD or retained a current mesh/result");
    const auto changedLocal = DraftGeometry(changed.parameters).previewMesh(.001);
    require(*changed.preview == assemblyPreview(changedLocal, changed),
            "Shared edit did not update both live horn instances");
    std::set<std::uint32_t> driven;
    for (const auto &triangle : changed.preview->triangles)
        if (triangle.boundary == BoundaryCondition::driven &&
            triangle.faceId < assemblyCadFaceOffset)
            driven.insert(triangle.faceId);
    require(driven == std::set<std::uint32_t>{5, 13}, "Shared edit lost either horn source");
    require(undo->isEnabled(), "Shared parameter edit cannot be undone");
    undo->trigger();
    waitUntil([&] { return main.geometryReady() && length->value() == 40; });
    saveAndWait(main, save);
    const auto restored = readProject(path);
    require(restored.assembly == p.assembly && *restored.preview == *p.preview && restored.cad,
            "Undo did not restore both horns and the CAD backdrop");
    main.exerciseLoadedProjectRoundTrip();
    std::cout << "CPU GUI: linked Draft edits, CAD coexistence, source IDs, invalidation, undo and "
                 "persistence passed.\n";
}
void testDiffractionSlotEditing(const QString &examplePath = {}) {
    QTemporaryDir directory;
    require(directory.isValid(), "Slot test directory is unavailable.");
    const auto path = directory.filePath("slot-editing.hsproj");
    Project baseline;
    baseline.preview = std::make_shared<SurfaceMesh>(
        DraftGeometry(baseline.parameters)
            .previewMesh(baseline.meshing.adaptive.geometricTolerance));
    baseline.geometryRevision = baseline.revision;
    writeProject(path, baseline);
    MainWindow main;
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() && main.geometryReady(),
            "Baseline Draft project did not open.");
    auto *dock = dynamic_cast<SimulationDock *>(main.findChild<QDockWidget *>("simulationDock"));
    auto *enabled = main.findChild<QCheckBox *>("slotEnabled");
    auto *repair = main.findChild<QPushButton *>("slotFormStart");
    auto *width = main.findChild<QDoubleSpinBox *>("slotWidth");
    auto *formStart = main.findChild<QDoubleSpinBox *>("transitionStart");
    auto *tabs = main.findChild<QTabWidget *>("draftSections");
    auto *panel = tabs ? dynamic_cast<ParameterPanel *>(tabs->parentWidget()) : nullptr;
    auto *undo = main.findChild<QAction *>("undoEdit");
    auto *redo = main.findChild<QAction *>("redoEdit");
    auto *history = main.findChild<QUndoStack *>();
    QAction *save = nullptr;
    for (auto *action : main.findChildren<QAction *>())
        if (action->shortcut() == QKeySequence::Save)
            save = action;
    require(dock && enabled && repair && width && formStart && panel && undo && redo && history &&
                save,
            "Slot editing, save or history controls are missing.");

    // This fixture exercises stale-result presentation only. It is never saved or exported.
    const auto cleanArchive = dock->archive();
    const std::vector<FaceVelocity> velocities{{5, {1, 0}}};
    dock->setMesh(baseline.preview, main.revision(), velocities);
    auto synthetic = cleanArchive;
    synthetic.settings.minimumFrequency = synthetic.settings.maximumFrequency = 1000;
    synthetic.settings.frequencyCount = 1;
    synthetic.mesh = baseline.preview;
    const auto frame = dock->observationFrame();
    synthetic.request =
        makeSimulationRequest(synthetic.settings, velocities, frame.origin, frame.forward);
    auto oldResult = std::make_shared<BemResult>();
    oldResult->device = "SYNTHETIC GUI SNAPSHOT - no GPU computation";
    oldResult->triangleCount = baseline.preview->triangles.size();
    oldResult->pressureDofCount = baseline.preview->vertices.size();
    BemExcitationResult output;
    output.name = synthetic.request->excitations.front().name;
    BemBoundaryData boundary;
    boundary.surfacePressure.assign(oldResult->triangleCount, {.01, 0});
    output.setBoundary(std::move(boundary));
    output.observerPressure.assign(synthetic.request->observers.size(), {.01, 0});
    oldResult->frequencies.push_back({1000, {std::move(output)}});
    synthetic.result = oldResult;
    synthetic.current = true;
    synthetic.revision = main.revision();
    dock->restoreArchive(synthetic);
    require(dock->hasCurrentResult(), "The synthetic stale-state fixture is not current.");

    const auto oldRevision = main.revision();
    enabled->setChecked(true);
    require(main.revision() > oldRevision && !main.geometryReady() && !dock->hasCurrentResult() &&
                dock->result() == oldResult && formStart->value() == 0 && repair->isEnabled(),
            "Slot activation failed to invalidate the old snapshot or silently repaired inputs.");
    repair->click();
    waitUntil([&] { return main.geometryReady(); });
    require(panel->parameters().slot.enabled && formStart->value() > 0 &&
                dock->result() == oldResult && !dock->hasCurrentResult() &&
                !main.simulationInProgress() && !main.meshReady(),
            "Slot preview overwrote old results or incorrectly claimed a current BEM mesh.");
    dock->restoreArchive(cleanArchive);
    saveAndWait(main, save);
    const auto slotted = readProject(path);
    require(slotted.parameters.slot.enabled && slotted.preview &&
                slotted.geometryRevision == slotted.revision && !slotted.simulation.result &&
                slotted.parameters.length == baseline.parameters.length,
            "Live slot project did not persist its current geometry or retained the GUI fixture.");
    std::set<std::uint32_t> tags, sources;
    for (const auto &triangle : slotted.preview->triangles) {
        tags.insert(triangle.faceId);
        if (triangle.boundary == BoundaryCondition::driven) {
            sources.insert(triangle.faceId);
            for (const auto vertex : triangle.vertices)
                require(slotted.preview->vertices[vertex].z == 0,
                        "The source moved from the driver entry onto the slot exit.");
        }
    }
    require(tags.contains(6) && tags.contains(7) && tags.contains(8) &&
                sources == std::set<std::uint32_t>{5},
            "Live slot geometry lost its material tags or original source surface.");
    const auto direct =
        DraftGeometry(slotted.parameters).previewMesh(slotted.meshing.adaptive.geometricTolerance);
    require(*slotted.preview == direct,
            "Saved live slot preview differs from the neutral geometry.");
    waitUntil([&] { return history->count() > 0; });
    const auto beforeWidthEdit = main.revision();
    width->setValue(14);
    waitUntil([&] { return main.geometryReady() && main.revision() > beforeWidthEdit; });
    require(undo->isEnabled(), "The slot width edit cannot be undone.");
    undo->trigger();
    waitUntil([&] { return main.geometryReady(); });
    require(panel->parameters().slot.enabled &&
                panel->parameters().slot.aperture.width == slotted.parameters.slot.aperture.width,
            "Undo lost slot activation or failed to restore the aperture.");
    redo->trigger();
    waitUntil([&] { return main.geometryReady(); });
    require(panel->parameters().slot.aperture.width == .014,
            "Redo failed to restore the slot aperture.");
    saveAndWait(main, save);
    const auto edited = readProject(path);
    main.exerciseLoadedProjectRoundTrip();
    require(panel->parameters().slot.aperture.width == .014 && main.geometryReady() &&
                !main.simulationInProgress(),
            "Slot project reload lost parameters, its current preview or started a solver.");
    if (!examplePath.isEmpty()) {
        writeProject(examplePath, edited);
        const auto restored = readProject(examplePath);
        require(restored.parameters.slot.enabled &&
                    restored.parameters.slot.aperture.width == .014 && restored.preview &&
                    *restored.preview == *edited.preview && !restored.simulation.result,
                "The exported slot starting project did not round-trip.");
    }
    std::cout
        << "CPU GUI: slot activation, explicit form start, live surfaces/source, stale result "
           "retention, undo/redo and project round-trip passed. No GPU calculation.\n";
}

void testDiffractionSlotRadii(const QString &examplePath = {}) {
    QTemporaryDir directory;
    require(directory.isValid(), "Slot-radius test directory is unavailable.");
    const auto path = directory.filePath("slot-radii.hsproj");
    Project baseline;
    auto &p = baseline.parameters;
    p.slot.enabled = true;
    p.slot.contour = SlotContour::roundedRectangle;
    p.slot.cornerRadius = .002;
    p.slot.exit = SlotExit::circularFillet;
    p.slot.exitRadius = .002;
    p.slot.aperture = {.012, .035, 4};
    p.slot.adapterLength = .020;
    p.slot.straightLength = .010;
    p.profile.k = p.verticalProfile.k = 1;
    p.profile.throatAngle = .10;
    p.verticalProfile.throatAngle = .15;
    p.profile.terminationScale = p.verticalProfile.terminationScale = .12;
    p.horizontalAngle = .5;
    p.verticalAngle = .3;
    p.crossSectionExponent = 4;
    p.transition.startFraction = .3;
    p.transition.endFraction = .7;
    const DraftGeometry baselineGeometry(p);
    baseline.preview = std::make_shared<SurfaceMesh>(
        baselineGeometry.previewMesh(baseline.meshing.adaptive.geometricTolerance));
    baseline.geometryRevision = baseline.revision;
    writeProject(path, baseline);

    MainWindow main;
    main.resize(1480, 1000);
    main.show();
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() && main.geometryReady(),
            "The OS-SE project with two physical slot radii did not open.");
    auto *tabs = main.findChild<QTabWidget *>("draftSections");
    auto *panel = tabs ? dynamic_cast<ParameterPanel *>(tabs->parentWidget()) : nullptr;
    auto *corner = main.findChild<QDoubleSpinBox *>("slotCornerRadius");
    auto *exit = main.findChild<QDoubleSpinBox *>("slotExitRadius");
    auto *contourMode = main.findChild<QComboBox *>("slotContour");
    auto *exitMode = main.findChild<QComboBox *>("slotExit");
    auto *status = main.findChild<QLabel *>("draftGeometryStatus");
    auto *repair = main.findChild<QPushButton *>("slotFormStart");
    auto *parameterDock = main.findChild<QDockWidget *>("parameterDock");
    QAction *save = nullptr;
    for (auto *action : main.findChildren<QAction *>())
        if (action->shortcut() == QKeySequence::Save)
            save = action;
    require(tabs && panel && corner && exit && contourMode && exitMode && status && repair &&
                parameterDock && save,
            "Independent slot-radius inputs, status or project actions are missing.");
    tabs->setCurrentIndex(4);
    parameterDock->show();
    parameterDock->raise();
    main.resizeDocks({parameterDock}, {440}, Qt::Horizontal);
    QCoreApplication::processEvents();
    require(corner->isEnabled() && exit->isEnabled() && corner->isVisible() && exit->isVisible() &&
                contourMode->currentIndex() == static_cast<int>(SlotContour::roundedRectangle) &&
                exitMode->currentIndex() == static_cast<int>(SlotExit::circularFillet) &&
                corner->value() == 2 && exit->value() == 2,
            "Opening the radius project did not expose both independent millimetre inputs.");

    const auto savedPreview = [&] {
        saveAndWait(main, save);
        auto saved = readProject(path);
        require(saved.preview && saved.geometryRevision == saved.revision && !saved.mesh &&
                    !saved.simulation.result && !main.simulationInProgress(),
                "Radius editing saved stale geometry or started an unsolicited simulation.");
        const auto direct =
            DraftGeometry(saved.parameters).previewMesh(saved.meshing.adaptive.geometricTolerance);
        require(*saved.preview == direct,
                "The live radius preview differs from its saved neutral geometry.");
        return saved;
    };
    const auto awaitGeometry = [&] {
        waitUntil([&] {
            return main.geometryReady() || status->text().contains("Vorschau nicht aktuell");
        });
        if (!main.geometryReady())
            throw std::runtime_error(
                ("The live radius edit failed: " + status->text()).toStdString());
        require(status->isHidden(), "A successful radius edit retained an error banner.");
    };
    const auto beforeCorner = main.revision();
    corner->setValue(3);
    require(main.revision() > beforeCorner && !main.geometryReady() &&
                panel->parameters().slot.exitRadius == .002,
            "A corner-radius edit did not invalidate the preview or changed the exit radius.");
    awaitGeometry();
    const auto changedCorner = savedPreview();
    require(changedCorner.parameters.slot.cornerRadius == .003 &&
                changedCorner.parameters.slot.exitRadius == .002 &&
                *changedCorner.preview != *baseline.preview,
            "The corner-radius edit did not change its live geometry independently.");
    const auto beforeExit = main.revision();
    exit->setValue(4);
    require(main.revision() > beforeExit && !main.geometryReady() &&
                panel->parameters().slot.cornerRadius == .003,
            "An exit-radius edit changed the independent corner radius.");
    awaitGeometry();
    const auto changedExit = savedPreview();
    require(changedExit.parameters.slot.cornerRadius == .003 &&
                changedExit.parameters.slot.exitRadius == .004 &&
                *changedExit.preview != *changedCorner.preview &&
                changedExit.parameters.profile.k == 1 &&
                changedExit.parameters.verticalProfile.k == 1 &&
                changedExit.parameters.profile.throatAngle == p.profile.throatAngle &&
                changedExit.parameters.verticalProfile.throatAngle == p.verticalProfile.throatAngle,
            "The exit-radius edit did not update the true OS-SE geometry independently.");

    corner->setValue(18);
    waitUntil([&] { return status->text().contains("Vorschau nicht aktuell"); });
    require(status->isVisible() && !main.geometryReady() && !repair->isEnabled() &&
                panel->parameters().slot.cornerRadius == .018 &&
                panel->parameters().slot.exitRadius == .004 &&
                status->text().contains("Letzte gültige Vorschau"),
            "An impossible corner radius was silently corrected or kept a current preview.");
    corner->setValue(3);
    awaitGeometry();
    exit->setValue(10000);
    waitUntil([&] { return status->text().contains("Vorschau nicht aktuell"); });
    require(status->isVisible() && !main.geometryReady() && !repair->isEnabled() &&
                panel->parameters().slot.exitRadius == 10 &&
                panel->parameters().slot.cornerRadius == .003,
            "An impossible exit radius was silently corrected or changed the corner radius.");
    exit->setValue(4);
    awaitGeometry();
    const auto saved = savedPreview();
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() && main.geometryReady() &&
                panel->parameters().slot.cornerRadius == .003 &&
                panel->parameters().slot.exitRadius == .004 &&
                panel->parameters().slot.contour == SlotContour::roundedRectangle &&
                panel->parameters().slot.exit == SlotExit::circularFillet && corner->value() == 3 &&
                exit->value() == 4 && *readProject(path).preview == *saved.preview,
            "Reopening the edited project lost a radius, mode or exact current preview.");

    corner->setValue(2);
    exit->setValue(2);
    awaitGeometry();
    const auto example = savedPreview();
    if (!examplePath.isEmpty()) {
        writeProject(examplePath, example);
        main.openProject(examplePath);
        waitUntil([&] { return !main.sourceImportInProgress(); });
        const auto restored = readProject(examplePath);
        require(main.sourceImportError().isEmpty() && main.geometryReady() &&
                    restored.parameters.slot.cornerRadius == .002 &&
                    restored.parameters.slot.exitRadius == .002 && restored.preview &&
                    *restored.preview == *example.preview && !restored.simulation.result,
                "The exported OS-SE example with independent 2-mm radii did not round-trip.");
        tabs->setCurrentIndex(4);
        parameterDock->show();
        parameterDock->raise();
        main.resizeDocks({parameterDock}, {440}, Qt::Horizontal);
        QCoreApplication::processEvents();
        require(main.grab().save(examplePath + ".png"),
                "The slot-radius example screenshot could not be saved.");
    }
    std::cout << "CPU GUI: independent physical slot radii on OS-SE, live geometry, visible "
                 "invalid-radius diagnostics and saved/reopened preview passed. "
                 "No GPU calculation.\n";
}

void testSavePendingEdits() {
    QTemporaryDir directory;
    require(directory.isValid(), "Pending edit project directory is unavailable.");
    const auto path = directory.filePath("pending-edits.hsproj");
    writeProject(path, Project{});
    MainWindow main;
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty(), "Pending edit project import failed.");
    main.show();
    main.activateWindow();
    auto *length = main.findChild<QDoubleSpinBox *>("length");
    auto *editor = length ? length->findChild<QLineEdit *>() : nullptr;
    QAction *save = nullptr;
    for (auto *action : main.findChildren<QAction *>())
        if (action->shortcut() == QKeySequence::Save)
            save = action;
    require(length && editor && save, "Pending edit save controls are missing.");
    const auto typeLength = [&](const char *text) {
        editor->setFocus();
        waitUntil([&] { return editor->hasFocus(); });
        editor->selectAll();
        QTest::keyClicks(editor, text);
    };
    const double before = length->value();
    typeLength("211");
    require(length->value() == before, "The test did not leave an uncommitted numeric edit.");
    saveAndWait(main, save);
    const auto saved = readProject(path);
    const bool savedInput = saved.parameters.length == .211 && length->value() == 211;

    // A new uncommitted first edit must reach the unsaved-changes decision too.
    typeLength("212");
    bool asked = false;
    QTimer cancelQuestion;
    QObject::connect(&cancelQuestion, &QTimer::timeout, &main, [&] {
        if (auto *question = qobject_cast<QMessageBox *>(QApplication::activeModalWidget())) {
            asked = true;
            question->done(QMessageBox::Cancel);
        }
    });
    cancelQuestion.start(5);
    QCloseEvent close;
    QApplication::sendEvent(&main, &close);
    cancelQuestion.stop();
    const bool guardedClose = asked && !close.isAccepted() && length->value() == 212;

    auto source = std::make_shared<CadSnapshot>(testing::exactCadFixture());
    Project cadProject;
    cadProject.source = ProjectSource::freecad;
    cadProject.cad = source;
    cadProject.preview = std::shared_ptr<const SurfaceMesh>(source, &source->mesh);
    cadProject.geometryRevision = cadProject.revision;
    const auto cadPath = directory.filePath("pending-cad-target.hsproj");
    writeProject(cadPath, cadProject);
    main.hide();
    MainWindow cadMain;
    cadMain.openProject(cadPath);
    waitUntil([&] { return !cadMain.sourceImportInProgress(); });
    require(cadMain.sourceImportError().isEmpty(), "Pending CAD target project import failed.");
    cadMain.show();
    cadMain.activateWindow();
    auto *cad = cadMain.findChild<QDockWidget *>("freecadDock");
    auto *faces = cad ? cad->findChild<QTreeWidget *>("cadFaces") : nullptr;
    require(faces && faces->topLevelItemCount(), "CAD target table is unavailable.");
    auto *item = faces->topLevelItem(0);
    const auto faceId = item->data(0, Qt::UserRole).toUInt();
    faces->setCurrentItem(item, 3);
    faces->editItem(item, 3);
    auto *target = faces->findChild<QDoubleSpinBox *>("cadFaceTargetEditor");
    auto *targetText = target ? target->findChild<QLineEdit *>() : nullptr;
    require(target && targetText, "The real CAD target delegate editor did not open.");
    targetText->setFocus();
    waitUntil([&] { return targetText->hasFocus(); });
    targetText->selectAll();
    QTest::keyClicks(targetText, "2");
    QAction *cadSave = nullptr;
    for (auto *action : cadMain.findChildren<QAction *>())
        if (action->shortcut() == QKeySequence::Save)
            cadSave = action;
    saveAndWait(cadMain, cadSave);
    const auto savedCad = readProject(cadPath);
    const auto found = savedCad.cadFaceTargetEdgeLengths.find(faceId);
    const bool savedDelegate =
        found != savedCad.cadFaceTargetEdgeLengths.end() && found->second == .002;
    std::cout << "Pending edits: save=" << savedInput << ", close guard=" << guardedClose
              << ", CAD delegate=" << savedDelegate << '\n';
    require(savedInput, "Save lost a focused numeric edit that had not received Enter.");
    require(guardedClose, "Close discarded the first focused edit without an unsaved warning.");
    require(savedDelegate, "Save lost an open CAD face target delegate edit.");
}
void testAsyncSaveSnapshotAndClose() {
    QTemporaryDir directory;
    require(directory.isValid(), "Async save test directory is unavailable.");
    const auto path = directory.filePath("async-save.hsproj");
    writeProject(path, Project{});
    QSemaphore entered, proceed;
    MainWindow main(nullptr, [&](const QString &target, const Project &snapshot) {
        entered.release();
        if (!proceed.tryAcquire(1, 30000))
            throw std::runtime_error("Save test gate timed out.");
        writeProject(target, snapshot);
    });
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty(), "Async save source import failed.");
    auto *length = main.findChild<QDoubleSpinBox *>("length");
    QAction *save = nullptr;
    QAction *newProject = nullptr;
    for (auto *action : main.findChildren<QAction *>()) {
        if (action->shortcut() == QKeySequence::Save)
            save = action;
        if (action->shortcut() == QKeySequence::New)
            newProject = action;
    }
    require(length && save && newProject, "Async save controls are missing.");
    length->setValue(211);
    save->trigger();
    waitUntil([&] { return entered.tryAcquire(); });
    require(main.saveInProgress() && !save->isEnabled() && !newProject->isEnabled(),
            "A blocked save did not disable overlapping file actions.");
    length->setValue(212);
    QCloseEvent close;
    QApplication::sendEvent(&main, &close);
    require(!close.isAccepted(), "Close bypassed an active save.");
    bool asked = false;
    QTimer cancelQuestion;
    QObject::connect(&cancelQuestion, &QTimer::timeout, &main, [&] {
        if (auto *question = qobject_cast<QMessageBox *>(QApplication::activeModalWidget())) {
            asked = true;
            question->done(QMessageBox::Cancel);
        }
    });
    cancelQuestion.start(5);
    proceed.release();
    waitUntil([&] { return !main.saveInProgress(); });
    cancelQuestion.stop();
    require(asked && main.windowTitle().contains('*') && length->value() == 212,
            "A newer edit was lost or the deferred close skipped its unsaved decision.");
    require(readProject(path).parameters.length == .211,
            "Worker did not persist its original immutable snapshot.");
    require(save->isEnabled() && newProject->isEnabled(), "File actions were not restored.");
}

void testAsyncSaveFailureKeepsProject() {
    QTemporaryDir directory;
    require(directory.isValid(), "Save failure test directory is unavailable.");
    const auto path = directory.filePath("save-failure.hsproj");
    writeProject(path, Project{});
    const auto original = [&] {
        QFile file(path);
        require(file.open(QIODevice::ReadOnly), "Cannot read save failure sentinel.");
        return file.readAll();
    }();
    QSemaphore entered, proceed;
    MainWindow main(nullptr, [&](const QString &, const Project &) {
        entered.release();
        if (!proceed.tryAcquire(1, 30000))
            throw std::runtime_error("Save test gate timed out.");
        throw std::runtime_error("Injected save failure");
    });
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty(), "Save failure source import failed.");
    auto *length = main.findChild<QDoubleSpinBox *>("length");
    QAction *save = nullptr;
    for (auto *action : main.findChildren<QAction *>())
        if (action->shortcut() == QKeySequence::Save)
            save = action;
    require(length && save, "Save failure controls are missing.");
    length->setValue(213);
    save->trigger();
    waitUntil([&] { return entered.tryAcquire(); });
    QCloseEvent close;
    QApplication::sendEvent(&main, &close);
    require(!close.isAccepted(), "Close bypassed a pending failing save.");
    bool errorSeen = false;
    QTimer dismissError;
    QObject::connect(&dismissError, &QTimer::timeout, &main, [&] {
        if (auto *box = qobject_cast<QMessageBox *>(QApplication::activeModalWidget())) {
            errorSeen = box->text().contains("Injected save failure");
            box->accept();
        }
    });
    dismissError.start(5);
    proceed.release();
    waitUntil([&] { return !main.saveInProgress(); });
    dismissError.stop();
    QFile file(path);
    require(file.open(QIODevice::ReadOnly) && file.readAll() == original,
            "Failed save changed the existing target.");
    require(errorSeen && main.windowTitle().contains('*') && length->value() == 213 &&
                save->isEnabled(),
            "Failed save lost edits, path, actions, or its error notice.");
}

void testAsyncSaveIncomingResult() {
    QTemporaryDir directory;
    require(directory.isValid(), "Incoming result test directory is unavailable.");
    const auto path = directory.filePath("incoming-result.hsproj");
    writeProject(path, Project{});
    QSemaphore entered, proceed;
    MainWindow main(nullptr, [&](const QString &target, const Project &snapshot) {
        entered.release();
        if (!proceed.tryAcquire(1, 30000))
            throw std::runtime_error("Save test gate timed out.");
        writeProject(target, snapshot);
    });
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty(), "Incoming result source import failed.");
    SimulationDock *dock = nullptr;
    for (auto *candidate : main.findChildren<QDockWidget *>())
        if (auto *simulation = dynamic_cast<SimulationDock *>(candidate))
            dock = simulation;
    QAction *save = nullptr;
    for (auto *action : main.findChildren<QAction *>())
        if (action->shortcut() == QKeySequence::Save)
            save = action;
    require(dock && save && !main.windowTitle().contains('*'),
            "Incoming result test must start with a clean loaded project.");
    save->trigger();
    waitUntil([&] { return entered.tryAcquire(); });
    auto archive = dock->archive();
    archive.request =
        makeSimulationRequest(archive.settings, std::vector<FaceVelocity>{{1, {1, 0}}});
    archive.mesh = std::make_shared<SurfaceMesh>();
    auto result = std::make_shared<BemResult>();
    result->device = "SYNTHETIC SAVE LIFECYCLE FIXTURE";
    archive.result = result;
    archive.resultComplete = false;
    archive.current = false;
    dock->restoreArchive(archive);
    QCloseEvent close;
    QApplication::sendEvent(&main, &close);
    require(!close.isAccepted(), "Close bypassed the pending result save.");
    bool asked = false;
    QTimer cancelQuestion;
    QObject::connect(&cancelQuestion, &QTimer::timeout, &main, [&] {
        if (auto *question = qobject_cast<QMessageBox *>(QApplication::activeModalWidget())) {
            asked = true;
            question->done(QMessageBox::Cancel);
        }
    });
    cancelQuestion.start(5);
    proceed.release();
    waitUntil([&] { return !main.saveInProgress(); });
    cancelQuestion.stop();
    require(asked && main.windowTitle().contains('*') && dock->archive().result == result,
            "Incoming result did not keep a clean project dirty and guard deferred close.");
    require(!readProject(path).simulation.result,
            "Save incorporated a later result into its snapshot.");
}

void testCloseSupersedesPendingNew() {
    QTemporaryDir directory;
    require(directory.isValid(), "Pending new test directory is unavailable.");
    const auto path = directory.filePath("pending-new.hsproj");
    writeProject(path, Project{});
    QSemaphore entered, proceed;
    MainWindow main(nullptr, [&](const QString &target, const Project &snapshot) {
        entered.release();
        if (!proceed.tryAcquire(1, 30000))
            throw std::runtime_error("Save test gate timed out.");
        writeProject(target, snapshot);
    });
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty(), "Pending new source import failed.");
    main.show();
    auto *length = main.findChild<QDoubleSpinBox *>("length");
    QAction *newProject = nullptr;
    for (auto *action : main.findChildren<QAction *>())
        if (action->shortcut() == QKeySequence::New)
            newProject = action;
    require(length && newProject, "Pending new controls are missing.");
    length->setValue(214);
    bool asked = false;
    QTimer chooseSave;
    QObject::connect(&chooseSave, &QTimer::timeout, &main, [&] {
        if (auto *question = qobject_cast<QMessageBox *>(QApplication::activeModalWidget())) {
            asked = true;
            if (auto *button = question->button(QMessageBox::Save))
                button->click();
        }
    });
    chooseSave.start(5);
    newProject->trigger();
    chooseSave.stop();
    waitUntil([&] { return entered.tryAcquire(); });
    require(asked && main.saveInProgress(), "New did not defer until the requested save.");
    QCloseEvent close;
    QApplication::sendEvent(&main, &close);
    require(!close.isAccepted(), "Close bypassed a pending save for New.");
    proceed.release();
    waitUntil([&] { return !main.saveInProgress(); });
    require(!main.isVisible() && length->value() == 214 &&
                readProject(path).parameters.length == .214,
            "Close did not replace the pending New continuation after saving.");
}

void testAsyncSaveAsAndRepeat() {
    QTemporaryDir directory;
    require(directory.isValid(), "Save As directory is unavailable.");
    const auto originalPath = directory.filePath("original.hsproj");
    const auto newPath = directory.filePath("copy.hsproj");
    writeProject(originalPath, Project{});
    QSemaphore entered, proceed;
    MainWindow main(nullptr, [&](const QString &target, const Project &snapshot) {
        entered.release();
        if (!proceed.tryAcquire(1, 30000))
            throw std::runtime_error("Save test gate timed out.");
        writeProject(target, snapshot);
    });
    main.openProject(originalPath);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty(), "Save As source import failed.");
    auto *length = main.findChild<QDoubleSpinBox *>("length");
    QAction *save = nullptr, *saveAs = nullptr;
    for (auto *action : main.findChildren<QAction *>()) {
        if (action->shortcut() == QKeySequence::Save)
            save = action;
        if (action->shortcut() == QKeySequence::SaveAs)
            saveAs = action;
    }
    require(length && save && saveAs, "Save As controls are missing.");
    length->setValue(215);
    bool selected = false;
    QTimer choosePath;
    QObject::connect(&choosePath, &QTimer::timeout, &main, [&] {
        if (auto *dialog = qobject_cast<QFileDialog *>(QApplication::activeModalWidget())) {
            selected = true;
            dialog->selectFile(newPath);
            QMetaObject::invokeMethod(dialog, "accept", Qt::DirectConnection);
        }
    });
    const bool nativeDisabled = QCoreApplication::testAttribute(Qt::AA_DontUseNativeDialogs);
    QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs, true);
    choosePath.start(5);
    saveAs->trigger();
    choosePath.stop();
    QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs, nativeDisabled);
    waitUntil([&] { return entered.tryAcquire(); });
    require(selected && !save->isEnabled() && !saveAs->isEnabled(),
            "Save As did not block overlapping save actions.");
    proceed.release();
    waitUntil([&] { return !main.saveInProgress(); });
    require(main.windowTitle().contains("copy.hsproj") && !main.windowTitle().contains('*'),
            "Successful Save As did not adopt its new path.");
    length->setValue(216);
    save->trigger();
    waitUntil([&] { return entered.tryAcquire(); });
    proceed.release();
    waitUntil([&] { return !main.saveInProgress(); });
    require(readProject(newPath).parameters.length == .216 &&
                readProject(originalPath).parameters.length != .216,
            "Repeated Save did not use the successful Save As path.");
}

void testSaveWorkerDestruction() {
    QSemaphore entered, proceed, destructionStarted, writerFinished;
    bool completionCalled = false;
    auto worker =
        std::make_unique<ProjectSaveWorker>(nullptr, [&](const QString &, const Project &) {
            entered.release();
            if (!proceed.tryAcquire(1, 30000))
                throw std::runtime_error("Save test gate timed out.");
            writerFinished.release();
        });
    require(worker->start("unused.hsproj", Project{}, [&](QString) { completionCalled = true; }),
            "Save worker did not start.");
    waitUntil([&] { return entered.tryAcquire(); });
    std::jthread release([&] {
        destructionStarted.acquire();
        proceed.release();
    });
    destructionStarted.release();
    worker.reset();
    release.join();
    require(writerFinished.tryAcquire(), "Save worker was destroyed before its writer finished.");
    QCoreApplication::processEvents();
    require(!completionCalled, "Destroyed save worker delivered a stale completion.");
}

void testSaveWorkerRealIoFailureAndRecovery() {
    QTemporaryDir directory;
    require(directory.isValid(), "Save I/O test directory is unavailable.");
    const auto sentinel = directory.filePath("not-a-directory");
    QFile file(sentinel);
    require(file.open(QIODevice::WriteOnly), "Cannot create save I/O sentinel.");
    const QByteArray original("save-io-sentinel");
    require(file.write(original) == original.size(), "Cannot write save I/O sentinel.");
    file.close();
    ProjectSaveWorker worker;
    bool completed = false;
    QString error;
    require(worker.start(sentinel + "/child.hsproj", Project{},
                         [&](QString message) {
                             error = std::move(message);
                             completed = true;
                         }),
            "Save I/O failure test did not start.");
    waitUntil([&] { return completed; });
    require(!worker.busy() && !error.isEmpty(), "Real I/O failure was not reported.");
    require(file.open(QIODevice::ReadOnly) && file.readAll() == original,
            "Real I/O failure altered the sentinel.");
    file.close();
    const auto valid = directory.filePath("recovered.hsproj");
    completed = false;
    error.clear();
    require(worker.start(valid, Project{},
                         [&](QString message) {
                             error = std::move(message);
                             completed = true;
                         }),
            "Save worker did not restart after I/O failure.");
    waitUntil([&] { return completed; });
    require(!worker.busy() && error.isEmpty(), "Save worker did not recover from I/O failure.");
    (void)readProject(valid);
}

void testSubdomainGuiArchive(bool coupled = false) {
    std::cerr << "GUI stage window\n" << std::flush;
    QMainWindow window;
    std::cerr << "GUI stage dock\n" << std::flush;
    SimulationDock dock(&window);
    auto mesh = std::make_shared<SurfaceMesh>(testing::acousticCup());
    std::cerr << "GUI stage setMesh\n" << std::flush;
    dock.setMesh(mesh, 1, {{3, {1, 0}}});
    MeshingOptions options;
    options.targetEdgeLength = .02;
    dock.setPartitionMeshingOptions(options);
    dock.setAcousticRegions({{7, "Hornluft", {2, 3}}});
    if (coupled) {
        lem::Network network;
        network.nodes = {{{1}, lem::Domain::acoustic, true, "reference", {}},
                         {{2}, lem::Domain::acoustic, false, "driver", {}}};
        network.elements.push_back({{1}, "flow", {}, lem::FlowSource{{{1}, {2}}, {.001, 0}}});
        network.elements.push_back({{2}, "port", {}, lem::BemPort{{{2}, {1}}, {3}, 1}});
        dock.setCoupling(std::move(network), true);
    }
    dock.setFrequencySweep(500, 500, 1);
    std::cerr << "GUI stage requestSolve\n" << std::flush;
    dock.requestSolve();
    std::cerr << "GUI stage wait busy=" << dock.isBusy()
              << " error=" << dock.lastError().toStdString() << "\n"
              << std::flush;
    waitUntil([&] { return !dock.isBusy(); });
    std::cerr << "GUI stage completed error=" << dock.lastError().toStdString() << "\n"
              << std::flush;
    require(dock.lastError().isEmpty() && dock.hasCurrentResult(), "GUI subdomain solve failed.");
    auto archive = dock.archive();
    require(archive.request && archive.request->subdomains && archive.sourceMesh &&
                archive.result &&
                !archive.result->frequencies.front().excitations.front().subdomainTraces().empty(),
            "GUI subdomain snapshot is incomplete.");
    // A stale result remains a fully readable numerical snapshot even if its
    // source model has subsequently been replaced.
    Project project;
    project.revision = 2;
    project.simulation = archive;
    project.simulation.current = false;
    QTemporaryDir directory;
    require(directory.isValid(), "Temporary directory unavailable.");
    const auto path = directory.filePath("subdomain.hsproj");
    std::cerr << "GUI stage write\n" << std::flush;
    writeProject(path, project);
    std::cerr << "GUI stage read\n" << std::flush;
    const auto restored = readProject(path).simulation;
    require(restored.request && restored.request->subdomains == archive.request->subdomains &&
                restored.resultRegions == archive.resultRegions && restored.sourceMesh &&
                *restored.sourceMesh == *mesh &&
                !restored.result->frequencies.front().excitations.front().subdomainTraces().empty(),
            "GUI subdomain archive did not preserve its source and local traces.");
    for (const bool damagePressure : {true, false}) {
        auto altered = std::make_shared<BemResult>(*archive.result);
        auto &excitation = altered->frequencies.front().excitations.front();
        auto traces = excitation.subdomainTraces();
        auto &trace = traces.front();
        require(!trace.interfaceFluxes.empty() &&
                    !trace.interfaceFluxes.front().vertexIndices.empty(),
                "Missing interface flux fixture.");
        if (damagePressure) {
            const auto vertex = trace.interfaceFluxes.front().vertexIndices.front();
            const auto found =
                std::find(trace.vertexIndices.begin(), trace.vertexIndices.end(), vertex);
            require(found != trace.vertexIndices.end(), "Interface pressure fixture missing.");
            trace.pressureTrace[static_cast<std::size_t>(found - trace.vertexIndices.begin())] +=
                std::complex<double>{1, 0};
        } else
            trace.interfaceFluxes.front().normalDerivative.front() += std::complex<double>{1, 0};
        excitation.setSubdomainTraces(std::move(traces));
        project.simulation.result = altered;
        bool rejected = false;
        try {
            writeProject(
                directory.filePath(damagePressure ? "bad-pressure.hsproj" : "bad-flux.hsproj"),
                project);
        } catch (const std::exception &) {
            rejected = true;
        }
        require(rejected, "Inconsistent interface trace was accepted into project archive.");
    }
    std::cerr << "GUI stage test exit\n" << std::flush;
}
void testGmshSubdomainEndToEnd() {
    auto importedMesh = testing::acousticCup();
    for (auto &triangle : importedMesh.triangles)
        triangle.boundary = BoundaryCondition::wall;
    const auto center = static_cast<std::uint32_t>(importedMesh.vertices.size());
    importedMesh.vertices.push_back({0, 0, .08});
    for (std::uint32_t i = 0; i < 16; ++i)
        importedMesh.triangles.push_back(
            {{center, 64 + (i + 1) % 16, 64 + i}, 6, BoundaryCondition::wall});
    require(importedMesh.triangles.size() < 10000, "Gmsh test fixture exceeds triangle budget.");
    MeshImportSelection selection;
    selection.lengthScaleToMetres = .001;
    selection.imported.mesh = importedMesh;
    for (std::uint64_t i = 0; i < importedMesh.vertices.size(); ++i)
        selection.imported.nodeTags.push_back(i + 1);
    for (std::size_t i = 0; i < importedMesh.triangles.size(); ++i)
        selection.imported.triangleMetadata.push_back(
            {i + 1, importedMesh.triangles[i].faceId, {}});
    auto cad = std::make_shared<CadSnapshot>(*makeGmshCadSnapshot(selection, "cup-interface.msh"));
    for (auto &face : cad->faces)
        if (face.id == 3)
            face.boundary = BoundaryCondition::driven;
    for (auto &triangle : cad->mesh.triangles)
        if (triangle.faceId == 3)
            triangle.boundary = BoundaryCondition::driven;
    Project project;
    project.source = ProjectSource::freecad;
    project.gmsh = selection;
    project.cad = cad;
    project.preview = std::shared_ptr<const SurfaceMesh>(cad, &cad->mesh);
    project.geometryRevision = project.revision;
    project.meshing.targetEdgeLength = .02;
    project.meshing.quality.enabled = false;
    project.simulation.settings.acousticRegions = {{7, "Hornluft", {2, 3}, {6}}};
    project.simulation.settings.minimumFrequency = project.simulation.settings.maximumFrequency =
        500;
    project.simulation.settings.frequencyCount = 1;
    QTemporaryDir directory;
    require(directory.isValid(), "Gmsh integration directory unavailable.");
    const auto path = directory.filePath("gmsh-cup.hsproj");
    writeProject(path, project);
    std::cerr << "Gmsh GUI stage open\n" << std::flush;
    MainWindow main;
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() && main.geometryReady(),
            "Gmsh source project did not open.");
    main.show();
    std::cerr << "Gmsh GUI stage mesh\n" << std::flush;
    main.requestMesh();
    waitUntil([&] { return !main.meshInProgress(); });
    require(main.lastMeshError().isEmpty() && main.meshReady(),
            "Gmsh source mesh did not validate.");
    auto *dock = dynamic_cast<SimulationDock *>(main.findChild<QDockWidget *>("simulationDock"));
    require(dock && dock->archive().settings.acousticRegions ==
                        project.simulation.settings.acousticRegions,
            "Gmsh region definition changed.");
    std::cerr << "Gmsh GUI stage solve\n" << std::flush;
    main.requestSolve();
    waitUntil([&] { return !main.simulationInProgress(); });
    require(main.lastSimulationError().isEmpty() && main.simulationReady(),
            "Gmsh subdomain solve failed.");
    require(dock->archive().sourceMesh && *dock->archive().sourceMesh == cad->mesh,
            "Gmsh physical mesh was changed during validation.");
    auto *tree = main.findChild<QTreeWidget *>("projectTree");
    require(tree, "Project tree unavailable.");
    QTreeWidgetItemIterator item(tree);
    bool found = false;
    for (; *item; ++item)
        if ((*item)->text(0).contains("Interface"))
            found = true;
    require(found, "Solved BEM interface absent from model tree.");
    auto *interfaces = main.findChild<QCheckBox *>("showBemLemInterfaces");
    auto *workflow = main.findChild<QLabel *>("workflowStatus");
    require(interfaces && interfaces->isEnabled() && workflow,
            "Solved BEM interface presentation is unavailable.");
    const auto solvedRevision = main.revision();
    interfaces->setChecked(true);
    waitUntil([&] { return workflow->text().startsWith("Ergebnis aktuell"); });
    require(main.revision() == solvedRevision && main.simulationReady(),
            "Interface visibility changed the acoustic snapshot.");
    if (const auto screenshot = qEnvironmentVariable("HS_SUBDOMAIN_SCREENSHOT");
        !screenshot.isEmpty())
        require(main.grab().save(screenshot), "Gmsh GUI screenshot failed.");
    QAction *save = nullptr;
    for (auto *action : main.findChildren<QAction *>())
        if (action->shortcut() == QKeySequence::Save)
            save = action;
    require(save, "Save action unavailable.");
    std::cerr << "Gmsh GUI stage save\n" << std::flush;
    QString saveError;
    QTimer saveErrorWatcher;
    QObject::connect(&saveErrorWatcher, &QTimer::timeout, &main, [&] {
        auto *box = qobject_cast<QMessageBox *>(QApplication::activeModalWidget());
        if (box && box->icon() == QMessageBox::Critical) {
            saveError = box->text();
            box->reject();
        }
    });
    saveErrorWatcher.start(10);
    saveAndWait(main, save);
    saveErrorWatcher.stop();
    if (!saveError.isEmpty())
        throw std::runtime_error(
            (QStringLiteral("Gmsh project save failed: ") + saveError).toStdString());
    const auto saved = readProject(path);
    require(saved.gmsh && saved.simulation.current && saved.simulation.request &&
                saved.simulation.request->subdomains,
            "Saved Gmsh subdomain result is incomplete.");
    MainWindow reopened;
    reopened.openProject(path);
    waitUntil([&] { return !reopened.sourceImportInProgress(); });
    require(reopened.sourceImportError().isEmpty() && reopened.simulationReady(),
            "Saved Gmsh project did not reopen current.");
    std::cerr << "Gmsh GUI stage done\n" << std::flush;
}
void testGmshImportPreservesLiveLink() {
    MainWindow main;
    auto *cad = dynamic_cast<FreeCadDock *>(main.findChild<QDockWidget *>("freecadDock"));
    require(cad, "CAD dock is missing.");
    auto prior = std::make_shared<CadSnapshot>(testing::exactCadFixture());
    cad->restoreSnapshot(prior);
    cad->ensureListening();
    auto *link = cad->findChild<FreeCadLink *>();
    require(link && cad->hasLiveLink(), "CAD listener did not start.");
    QTcpSocket client;
    client.connectToHost(QHostAddress::LocalHost, link->port());
    waitUntil([&] { return client.state() == QAbstractSocket::ConnectedState; });
    const auto hello = QJsonDocument(QJsonObject{{"type", "hello"},
                                                 {"protocol", 2},
                                                 {"sessionId", "gmsh-preservation-test"},
                                                 {"token", link->token()}})
                           .toJson(QJsonDocument::Compact);
    QByteArray frame(4, Qt::Uninitialized);
    qToBigEndian(static_cast<quint32>(hello.size()), frame.data());
    frame += hello;
    require(client.write(frame) == frame.size(), "CAD hello was not queued.");
    waitUntil([&] { return link->isConnected(); });
    const auto old = cad->snapshot();
    std::istringstream input("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n3\n1 0 0 0\n2 1 0 0\n3 "
                             "0 1 0\n$EndNodes\n$Elements\n1\n1 2 2 7 5 1 2 3\n$EndElements\n");
    MeshImportSelection selection{importGmsh(input, {1.0}), 1.0};
    Project project;
    project.source = ProjectSource::freecad;
    project.gmsh = selection;
    project.cad = makeGmshCadSnapshot(selection, "fixture.msh");
    project.preview = std::shared_ptr<const SurfaceMesh>(project.cad, &project.cad->mesh);
    project.geometryRevision = project.revision;
    QTemporaryDir directory;
    require(directory.isValid(), "Temporary directory unavailable.");
    const auto path = directory.filePath("gmsh.hsproj");
    writeProject(path, project);
    const auto roundtrip = readProject(path);
    require(roundtrip.gmsh && roundtrip.cad->sessionId == "gmsh-import-offline",
            "Gmsh project metadata was not retained.");
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(!main.sourceImportError().isEmpty() && cad->hasLiveLink() && link->isConnected() &&
                cad->snapshot() == old,
            "Gmsh project import disturbed a live CAD link or its snapshot.");
}
void testCadEditIdentity() {
    MainWindow main;
    waitUntil([&] { return main.geometryReady(); });
    auto source = std::make_shared<CadSnapshot>(testing::exactCadFixture());
    const auto editedId = source->faces.front().id;
    const auto topologyKey = source->faces.front().signature.topologyKey;
    source->faces.front().boundary = BoundaryCondition::driven;
    for (auto &triangle : source->mesh.triangles)
        if (triangle.faceId == editedId)
            triangle.boundary = BoundaryCondition::driven;
    Project project;
    project.source = ProjectSource::freecad;
    project.cad = source;
    project.preview = std::shared_ptr<const SurfaceMesh>(source, &source->mesh);
    auto mesh = std::make_shared<MeshingResult>();
    mesh->mesh = source->mesh;
    mesh->quality = inspectMesh(mesh->mesh);
    project.mesh = mesh;
    project.geometryRevision = project.meshRevision = project.revision;
    QTemporaryDir directory;
    require(directory.isValid(), "Temporary CAD identity project directory is unavailable.");
    const auto path = directory.filePath("cad-identity.hsproj");
    writeProject(path, project);
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty(), "CAD identity archive import failed.");
    require(main.meshReady(), "The initial CAD identity mesh was not restored.");

    auto *cad = dynamic_cast<FreeCadDock *>(main.findChild<QDockWidget *>("freecadDock"));
    auto *link = cad ? cad->findChild<FreeCadLink *>() : nullptr;
    auto *faces = main.findChild<QTreeWidget *>("meshFaces");
    auto *panel = faces ? dynamic_cast<MeshPanel *>(faces->parentWidget()) : nullptr;
    auto *cadFaces = cad ? cad->findChild<QTreeWidget *>("cadFaces") : nullptr;
    MeshViewport *viewport = nullptr;
    for (auto *widget : main.findChildren<QOpenGLWidget *>())
        if (auto *candidate = dynamic_cast<MeshViewport *>(widget))
            viewport = candidate;
    auto *undo = main.findChild<QAction *>("undoEdit");
    require(cad && link && panel && cadFaces && viewport && undo,
            "CAD identity controls are missing.");
    auto *representation = main.findChild<QComboBox *>("geometryRepresentation");
    require(representation, "Geometry representation control is missing.");
    const auto currentMeshPick = clickViewportFace(*viewport);
    require(cad->selectedFaces() == std::vector<std::uint32_t>{currentMeshPick} &&
                panel->selectedFaces() == cad->selectedFaces() &&
                viewport->selectedFaces() == cad->selectedFaces(),
            "A current CAD mesh pick did not select its physical face consistently.");
    auto *simulation = main.findChild<QDockWidget *>("simulationDock");
    QPushButton *solve = nullptr;
    if (simulation)
        for (auto *button : simulation->findChildren<QPushButton *>())
            if (button->text() == QStringLiteral("Berechnen"))
                solve = button;
    require(solve && solve->isEnabled(), "The initial CAD mesh was not released for simulation.");
    cad->ensureListening();
    QTcpSocket client;
    client.connectToHost(QHostAddress::LocalHost, link->port());
    waitUntil([&] { return client.state() == QAbstractSocket::ConnectedState; });
    const auto send = [&](const QJsonObject &object) {
        const auto payload = QJsonDocument(object).toJson(QJsonDocument::Compact);
        QByteArray frame(4, Qt::Uninitialized);
        qToBigEndian(static_cast<quint32>(payload.size()), frame.data());
        frame += payload;
        require(client.write(frame) == frame.size(), "CAD identity frame was not queued.");
    };
    const auto session = QStringLiteral("gui-identity-session");
    send({{"type", "hello"}, {"protocol", 2}, {"sessionId", session}, {"token", link->token()}});
    waitUntil([&] { return link->isConnected(); });
    const auto publish = [&](bool enlarged, const QString &document, unsigned revision) {
        auto frame = testing::cadFixture(enlarged);
        frame.insert("sessionId", session);
        frame.insert("documentId", document);
        frame.insert("revision", QString::number(revision));
        // Reorder real source faces and their BRep binding together. FaceN remains
        // a traversal label; geometric/topological signatures identify the surface.
        const auto original = frame["faces"].toArray();
        QJsonArray reordered;
        for (qsizetype i = original.size(); i-- > 0;) {
            auto face = original[i].toObject();
            face.insert("mappedName", QStringLiteral("Face%1").arg(reordered.size() + 1));
            reordered.append(face);
        }
        auto parts = frame["brepParts"].toArray();
        for (qsizetype i = 0; i < parts.size(); ++i) {
            auto part = parts[i].toObject();
            QJsonArray indices;
            for (const auto &index : part["faceIndices"].toArray())
                indices.append(static_cast<int>(original.size() - 1 - index.toInt()));
            part.insert("faceIndices", indices);
            parts[i] = part;
        }
        frame.insert("faces", reordered);
        frame.insert("brepParts", parts);
        send(frame);
        waitUntil([&] {
            const auto current = cad->snapshot();
            return current && current->sessionId == session && current->documentId == document &&
                   current->revision == revision && cad->isSourceAvailable();
        });
    };
    publish(false, source->documentId, 1);
    cad->selectFace(editedId);
    auto *target = panel->findChild<QDoubleSpinBox *>("meshLocalTarget");
    auto *assign = panel->findChild<QPushButton *>("meshAssignLocalTarget");
    auto *reverse = panel->findChild<QPushButton *>("meshReverseFaces");
    require(target && assign && reverse, "CAD face edit controls are missing.");
    target->setValue(2);
    assign->click();
    reverse->click();
    const auto editRevision = main.revision();
    const auto verifyEditedFace = [&] {
        const auto current = cad->snapshot();
        const auto found =
            std::find_if(current->faces.begin(), current->faces.end(), [&](const auto &face) {
                return face.signature.topologyKey == topologyKey;
            });
        require(found != current->faces.end() && found->id == editedId &&
                    found->boundary == BoundaryCondition::driven,
                "A uniquely matched source face lost its identity or boundary.");
        const std::vector<std::uint32_t> selected{editedId};
        require(cad->selectedFaces() == selected && panel->selectedFaces() == selected &&
                    viewport->selectedFaces() == selected,
                "Same-source CAD edits lost synchronized face selection.");
        for (int row = 0; row < cadFaces->topLevelItemCount(); ++row) {
            const auto *item = cadFaces->topLevelItem(row);
            if (item->data(0, Qt::UserRole).toUInt() == editedId)
                require(std::abs(item->data(3, Qt::UserRole).toDouble() - .002) < 1e-12 &&
                            item->text(4) == QString::fromUtf8("−"),
                        "Same-source CAD edit lost its local target or normal reversal.");
        }
        for (int row = 0; row < faces->topLevelItemCount(); ++row) {
            const auto *item = faces->topLevelItem(row);
            if (item->data(0, Qt::UserRole).toUInt() == editedId)
                require(item->text(1) != QStringLiteral("global") &&
                            item->text(2) == QString::fromUtf8("−"),
                        "Shared mesh panel disagrees with the preserved CAD edits.");
        }
    };
    verifyEditedFace();
    require(!main.meshReady() && !main.simulationReady() && !solve->isEnabled(),
            "A CAD orientation edit retained the old solver-ready mesh.");
    publish(true, source->documentId, 2);
    verifyEditedFace();
    require(main.revision() > editRevision && !main.meshReady() && !main.simulationReady() &&
                !solve->isEnabled() && !undo->isEnabled(),
            "A new CAD geometry retained an old mesh, result or incompatible undo history.");

    const auto refuseStalePick = [&] {
        representation->setCurrentIndex(1);
        cad->selectFace(cad->snapshot()->faces.front().id);
        const auto before = main.revision();
        (void)clickViewportFace(*viewport);
        require(cad->selectedFaces().empty() && panel->selectedFaces().empty() &&
                    viewport->selectedFaces().empty() && !reverse->isEnabled() &&
                    !assign->isEnabled(),
                "A stale CAD mesh pick retained an editable selection from another revision.");
        reverse->click();
        assign->click();
        require(main.revision() == before,
                "A refused stale pick changed normals or local sizing on a previous selection.");
        representation->setCurrentIndex(0);
        const auto currentPreviewPick = clickViewportFace(*viewport);
        require(cad->selectedFaces() == std::vector<std::uint32_t>{currentPreviewPick} &&
                    panel->selectedFaces() == cad->selectedFaces() &&
                    viewport->selectedFaces() == cad->selectedFaces() && reverse->isEnabled(),
                "Refusing a stale mesh disabled valid current CAD preview picking.");
    };
    // The new revision retains the same physical IDs. Checking contains(faceId)
    // alone must therefore never authorize a pick from the old displayed mesh.
    refuseStalePick();
    cad->selectFace(editedId);
    verifyEditedFace();

    publish(true, QStringLiteral("another-gui-document"), 1);
    require(cad->selectedFaces().empty() && panel->selectedFaces().empty() &&
                viewport->selectedFaces().empty(),
            "Switching CAD documents retained another document's selection.");
    for (const auto &face : cad->snapshot()->faces)
        require(std::none_of(source->faces.begin(), source->faces.end(),
                             [&](const auto &old) { return old.id == face.id; }) &&
                    face.boundary == BoundaryCondition::wall,
                "Different documents reused physical face identities or boundary assignments.");
    for (int row = 0; row < cadFaces->topLevelItemCount(); ++row) {
        const auto *item = cadFaces->topLevelItem(row);
        require(item->data(3, Qt::UserRole).toDouble() == 0 && item->text(4) == "+",
                "A different CAD document inherited local sizing or normal reversals.");
    }
    for (int row = 0; row < faces->topLevelItemCount(); ++row) {
        const auto *item = faces->topLevelItem(row);
        require(item->text(1) == QStringLiteral("global") && item->text(2) == "+",
                "The shared mesh panel retained another CAD document's edits.");
    }
    require(!main.meshReady() && !main.simulationReady() && !solve->isEnabled() &&
                !undo->isEnabled(),
            "A different CAD document inherited an old solver state or edit history.");
    refuseStalePick();

    // Archives may retain an old mesh while a separately imported document uses
    // the same export-local IDs. Its revision is still stale and must win over IDs.
    auto independent = std::make_shared<CadSnapshot>(testing::exactCadFixture());
    independent->documentId = QStringLiteral("independent-archive-document");
    Project archive = project;
    archive.cad = independent;
    archive.preview = std::shared_ptr<const SurfaceMesh>(independent, &independent->mesh);
    archive.revision = archive.geometryRevision = 2;
    archive.meshRevision = 1;
    const auto reusedPath = directory.filePath("reused-face-ids.hsproj");
    writeProject(reusedPath, archive);
    main.openProject(reusedPath);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() &&
                cad->snapshot()->documentId == independent->documentId && !main.meshReady(),
            "The independent CAD archive did not retain its old mesh as stale.");
    refuseStalePick();
    const auto selected = cad->selectedFaces();
    reverse->click();
    require(main.revision() > archive.revision && undo->isEnabled(),
            "Current CAD preview normal editing did not remain undoable.");
    undo->trigger();
    require(cad->selectedFaces() == selected,
            "Undo after current preview picking lost the CAD selection.");
    for (int row = 0; row < cadFaces->topLevelItemCount(); ++row)
        require(cadFaces->topLevelItem(row)->text(4) == "+",
                "Undo did not restore current CAD face orientations after a refused stale pick.");
    std::cout << "CPU GUI: authenticated CAD updates preserve uniquely matched face edits; "
                 "document changes clear selection, normal reversals, local sizing and old "
                 "solver readiness; stale picks are refused before IDs and cannot edit faces.\n";
}
class CompletionFilter final : public QObject {
  public:
    std::function<bool()> isCompletion;
    std::function<void()> beforeCompletion;

  protected:
    bool eventFilter(QObject *, QEvent *event) override {
        if (event->type() == QEvent::MetaCall && beforeCompletion &&
            (!isCompletion || isCompletion())) {
            auto action = std::move(beforeCompletion);
            action();
        }
        return false;
    }
};

void testSourceImportStartFailure() {
    QTemporaryDir directory;
    require(directory.isValid(), "Import start test directory is unavailable.");
    const auto path = directory.filePath("import.hsproj");
    writeProject(path, Project{});

    for (const int failCall : {1, 2}) {
        int calls = 0;
        MainWindow main(nullptr, {}, [&](std::function<void(std::stop_token)> job) {
            if (++calls == failCall)
                throw std::system_error(
                    std::make_error_code(std::errc::resource_unavailable_try_again));
            return std::jthread(std::move(job));
        });
        waitUntil([&] { return main.geometryReady(); });
        const auto revision = main.revision();
        const auto title = main.windowTitle();
        main.openProject(path);
        if (failCall == 2)
            main.openProject(path); // Completion is queued; the second request is pending.
        waitUntil([&] { return !main.sourceImportInProgress(); });
        require(calls == failCall && !main.sourceImportError().isEmpty() &&
                    main.revision() == revision && main.windowTitle() == title,
                "A failed import thread start damaged the current project or missed its error.");
        for (auto *progress : main.findChildren<QProgressDialog *>("sourceImportProgress"))
            require(!progress->isVisible(),
                    "Import progress remained visible after start failure.");
        main.openProject(path);
        waitUntil([&] { return !main.sourceImportInProgress(); });
        require(main.sourceImportError().isEmpty() && calls == failCall + 1,
                "Import did not recover after a thread-start failure.");
    }

    int calls = 0;
    MainWindow main(nullptr, {}, [&](std::function<void(std::stop_token)> job) {
        if (++calls == 1)
            return std::jthread{};
        return std::jthread(std::move(job));
    });
    waitUntil([&] { return main.geometryReady(); });
    const auto revision = main.revision();
    main.openProject(path);
    require(!main.sourceImportInProgress() && !main.sourceImportError().isEmpty() &&
                main.revision() == revision,
            "A non-joinable import worker was treated as a running job.");
    for (auto *progress : main.findChildren<QProgressDialog *>("sourceImportProgress"))
        require(!progress->isVisible(), "Import progress remained visible without a worker.");
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() && calls == 2,
            "Import did not recover after a non-joinable worker.");
}

void testAsyncProjectOpen() {
    MainWindow main;
    waitUntil([&] { return main.geometryReady(); });
    QTemporaryDir directory;
    require(directory.isValid(), "Async import directory is unavailable.");
    auto cadSource = std::make_shared<CadSnapshot>(testing::exactCadFixture());
    cadSource->faces.front().boundary = BoundaryCondition::driven;
    for (auto &triangle : cadSource->mesh.triangles)
        if (triangle.faceId == cadSource->faces.front().id)
            triangle.boundary = BoundaryCondition::driven;
    Project project;
    project.source = ProjectSource::freecad;
    project.cad = cadSource;
    project.preview = std::shared_ptr<const SurfaceMesh>(cadSource, &cadSource->mesh);
    auto mesh = std::make_shared<MeshingResult>();
    mesh->mesh = cadSource->mesh;
    mesh->quality = inspectMesh(mesh->mesh);
    project.mesh = mesh;
    project.geometryRevision = project.meshRevision = project.revision;
    const auto path = directory.filePath("async-cad.hsproj");
    const auto missing = directory.filePath("missing.hsproj");
    writeProject(path, project);
    auto *mode = main.findChild<QComboBox *>("sourceMode");
    auto *cad = dynamic_cast<FreeCadDock *>(main.findChild<QDockWidget *>("freecadDock"));
    auto *length = main.findChild<QDoubleSpinBox *>("length");
    require(mode && cad && length, "Async import controls are missing.");
    int heartbeats = 0;
    QTimer heartbeat;
    QObject::connect(&heartbeat, &QTimer::timeout, [&] {
        if (main.sourceImportInProgress())
            ++heartbeats;
    });
    heartbeat.start(0);
    main.openProject(path);
    require(main.sourceImportInProgress() && mode->currentIndex() == 0 && !main.meshReady(),
            "Opening an archive did not return before committing it.");
    waitUntil([&] { return !main.sourceImportInProgress(); });
    heartbeat.stop();
    require(heartbeats > 0 && main.sourceImportError().isEmpty() && main.meshReady() &&
                cad->snapshot()->documentId == cadSource->documentId,
            "Asynchronous CAD opening blocked the event loop or failed to restore its mesh.");
    const auto retained = cad->snapshot();
    const auto retainedRevision = main.revision();
    const auto retainedTitle = main.windowTitle();
    main.openProject(missing);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(!main.sourceImportError().isEmpty() && main.meshReady() &&
                main.revision() == retainedRevision && cad->snapshot() == retained &&
                main.windowTitle() == retainedTitle,
            "An invalid archive replaced the current project or its clean state.");

    main.openProject(missing);
    main.openProject(missing);
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() && main.meshReady(),
            "A superseded import published an error instead of the latest pending archive.");
    QAction *newProject = nullptr;
    for (auto *action : main.findChildren<QAction *>())
        if (action->text() == QStringLiteral("Neu"))
            newProject = action;
    require(newProject, "New-project action is missing.");
    main.openProject(missing);
    newProject->trigger();
    waitUntil([&] { return !main.sourceImportInProgress() && main.geometryReady(); });
    require(mode->currentIndex() == 0 && main.sourceImportError().isEmpty() && !cad->snapshot(),
            "An old project's archive error overwrote a new document.");
    main.openProject(path);
    length->setValue(length->value() + 1);
    waitUntil([&] { return !main.sourceImportInProgress() && main.geometryReady(); });
    require(mode->currentIndex() == 0 && main.sourceImportError().isEmpty() && !main.meshReady(),
            "An archive response replaced a more recent parameter edit.");
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() && main.meshReady(),
            "A current archive did not load after obsolete imports retired.");
    auto *normalization = main.findChild<QComboBox *>("polarNormalization");
    require(normalization, "Display normalization control is missing.");
    const auto beforeDisplayEdit = cad->snapshot();
    const auto displayRevision = main.revision();
    main.openProject(path);
    normalization->setCurrentIndex(1 - normalization->currentIndex());
    const auto editedNormalization = normalization->currentIndex();
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.revision() == displayRevision && cad->snapshot() == beforeDisplayEdit &&
                normalization->currentIndex() == editedNormalization &&
                main.sourceImportError().contains("Import verworfen") &&
                main.windowTitle().contains("*"),
            "A background import overwrote a newer non-geometry project edit.");
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() && main.meshReady(),
            "A fresh import after a refused overwrite failed.");
    const auto beforeCancel = cad->snapshot();
    const auto cancelRevision = main.revision();
    CompletionFilter filter;
    bool cancelledQueued = false;
    filter.beforeCompletion = [&] {
        QProgressDialog *progress = nullptr;
        for (auto *candidate : main.findChildren<QProgressDialog *>("sourceImportProgress"))
            if (candidate->isVisible()) {
                require(!progress, "More than one import progress dialog is visible.");
                progress = candidate;
            }
        require(progress && main.sourceImportInProgress(), "Queued import progress is missing.");
        for (auto *button : progress->findChildren<QPushButton *>())
            if (button->text() == QStringLiteral("Abbrechen")) {
                button->click();
                cancelledQueued = true;
            }
        require(cancelledQueued, "The queued import could not be cancelled.");
    };
    main.installEventFilter(&filter);
    main.openProject(missing);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    main.removeEventFilter(&filter);
    require(cancelledQueued && main.sourceImportError().isEmpty() && main.meshReady() &&
                cad->snapshot() == beforeCancel && main.revision() == cancelRevision,
            "Cancellation after queued completion did not preserve the current project.");
    std::cout << "CPU GUI: asynchronous archive heartbeat, valid/invalid adoption, pending "
                 "latest-wins, epoch/revision rejection and queued cancellation passed.\n";
}

void testLatestPublication() {
    LatestPublication<std::pair<int, int>> slot;
    QObject receiver;
    int scheduled = 0, delivered = -1, generation = 1;
    const auto enqueue = [&](int job, int sample) {
        if (slot.offer({job, sample})) {
            ++scheduled;
            QMetaObject::invokeMethod(
                &receiver,
                [&] {
                    const auto value = slot.take();
                    if (value && value->first == generation)
                        delivered = value->second;
                },
                Qt::QueuedConnection);
        }
    };
    for (int i = 0; i < 1000; ++i)
        enqueue(1, i);
    require(scheduled == 1 && delivered == -1, "Unbounded queued partial snapshots");
    QCoreApplication::processEvents();
    require(delivered == 999, "The newest completed partial snapshot was lost");
    enqueue(1, 1000);
    generation = 2; // An obsolete queued callback cannot publish into the next job.
    QCoreApplication::processEvents();
    require(delivered == 999, "Canceled job published into a new generation");
    for (int i = 0; i < 1000; ++i)
        enqueue(2, 2000 + i);
    require(scheduled == 3, "Restart lost its bounded delivery scheduling");
    QCoreApplication::processEvents();
    require(delivered == 2999, "Restart failed to publish its latest partial snapshot");
}

void testJobLifecycle() {
    testLatestPublication();
    MainWindow main;
    waitUntil([&] { return main.geometryReady(); });
    main.requestMesh();
    require(main.meshInProgress(), "Draft meshing did not start.");
    const auto meshingRevision = main.revision();
    auto *length = main.findChild<QDoubleSpinBox *>("length");
    require(length, "Draft length input is missing.");
    length->setValue(length->value() + 1);
    waitUntil([&] { return !main.meshInProgress() && main.geometryReady(); });
    require(main.revision() > meshingRevision && main.lastMeshError().isEmpty() &&
                !main.meshReady() && !main.simulationReady(),
            "An obsolete Draft cancellation overwrote the current edit state.");

    QAction *cancelMeshing = nullptr;
    for (auto *action : main.findChildren<QAction *>())
        if (action->objectName() == QStringLiteral("cancelMeshing"))
            cancelMeshing = action;
    require(cancelMeshing, "The mesh cancellation action is missing.");
    CompletionFilter meshCompletion;
    bool cancelledQueued = false;
    const auto queuedRevision = main.revision();
    meshCompletion.beforeCompletion = [&] {
        require(main.meshInProgress() && cancelMeshing->isEnabled(),
                "The queued mesh completion no longer owns the active job.");
        cancelMeshing->trigger();
        cancelledQueued = true;
    };
    main.installEventFilter(&meshCompletion);
    main.requestMesh();
    waitUntil([&] { return !main.meshInProgress(); });
    main.removeEventFilter(&meshCompletion);
    require(cancelledQueued && main.revision() == queuedRevision && main.geometryReady() &&
                !main.meshReady() && !main.simulationReady() &&
                main.lastMeshError() == QStringLiteral("Vernetzung abgebrochen"),
            "Cancellation after queued mesh completion committed the cancelled result.");

    auto source = std::make_shared<CadSnapshot>(testing::exactCadFixture());
    source->faces.front().boundary = BoundaryCondition::driven;
    const auto drivenId = source->faces.front().id;
    for (auto &triangle : source->mesh.triangles)
        if (triangle.faceId == drivenId)
            triangle.boundary = BoundaryCondition::driven;
    Project project;
    project.source = ProjectSource::freecad;
    project.cad = source;
    project.preview = std::shared_ptr<const SurfaceMesh>(source, &source->mesh);
    auto mesh = std::make_shared<MeshingResult>();
    mesh->mesh = source->mesh;
    mesh->quality = inspectMesh(mesh->mesh);
    project.mesh = mesh;
    project.geometryRevision = project.meshRevision = project.revision;
    project.meshing.targetEdgeLength = .00025;
    QTemporaryDir directory;
    require(directory.isValid(), "Temporary CAD project directory is unavailable.");
    const auto path = directory.filePath("cad-lifecycle.hsproj");
    writeProject(path, project);
    main.openProject(path);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty(), "CAD lifecycle archive import failed.");
    require(main.geometryReady() && main.meshReady(), "Offline CAD mesh was not restored.");
    auto *cad = dynamic_cast<FreeCadDock *>(main.findChild<QDockWidget *>("freecadDock"));
    require(cad, "CAD assignment controls are missing.");
    const auto original = cad->snapshot();
    std::vector<CadFaceAssignment> assignments;
    for (const auto &face : original->faces)
        assignments.push_back({face.id, face.boundary, face.normalVelocity, face.state});
    assignments.front().normalVelocity = {.25, -.1};
    main.requestMesh();
    require(main.meshInProgress(), "CAD remeshing did not start.");
    const auto cadMeshingRevision = main.revision();
    require(cad->setAssignments(original->documentId, original->sessionId, original->revision,
                                assignments, original->unresolvedPreviousFaces),
            "CAD assignment edit was rejected.");
    waitUntil([&] { return !main.meshInProgress(); });
    require(main.revision() > cadMeshingRevision && main.meshReady() &&
                main.lastMeshError().isEmpty() &&
                cad->snapshot()->faces.front().normalVelocity == std::complex<double>(.25, -.1),
            "An obsolete CAD cancellation replaced a retained mesh or its new assignment.");
    auto *simulation = main.findChild<QDockWidget *>("simulationDock");
    require(simulation, "Simulation controls are missing.");
    bool solveEnabled = false;
    for (const auto *button : simulation->findChildren<QPushButton *>())
        if (button->text() == QStringLiteral("Berechnen"))
            solveEnabled = button->isEnabled();
    require(solveEnabled && !main.simulationInProgress(),
            "The retained current CAD mesh was not released for simulation.");
    const QDir testDirectory(QFileInfo(QString::fromUtf8(__FILE__)).absolutePath());
    const auto stepPath = testDirectory.filePath("fixtures/freecad-box-20x30x40mm.step");
    require(QFileInfo::exists(stepPath), "The exact STEP import fixture is missing.");
    main.openStep(stepPath);
    require(main.sourceImportInProgress(), "STEP import did not start.");
    main.openProject(path);
    const auto loadedRevision = project.revision;
    waitUntil([&] { return !main.sourceImportInProgress(); });
    const bool epochPreserved = main.sourceImportError().isEmpty() && main.meshReady() &&
                                main.revision() == loadedRevision &&
                                cad->snapshot()->documentId == source->documentId;

    auto *mode = main.findChild<QComboBox *>("sourceMode");
    require(mode, "The geometry source selector is missing.");
    mode->setCurrentIndex(0);
    waitUntil([&] { return main.geometryReady(); });
    const auto missingStep = directory.filePath("missing.step");
    main.openStep(missingStep);
    require(main.sourceImportInProgress(), "Failing STEP import did not start.");
    const auto oldRevision = main.revision();
    length->setValue(length->value() + 1);
    waitUntil([&] { return !main.sourceImportInProgress() && main.geometryReady(); });
    const bool revisionPreserved =
        main.sourceImportError().isEmpty() && main.revision() > oldRevision;

    const auto currentRevision = main.revision();
    main.openStep(missingStep);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(!main.sourceImportError().isEmpty() && main.revision() == currentRevision &&
                main.geometryReady(),
            "A current STEP failure was hidden or replaced the current geometry.");
    main.openStep(stepPath);
    waitUntil([&] { return !main.sourceImportInProgress(); });
    require(main.sourceImportError().isEmpty() && main.geometryReady() &&
                mode->currentIndex() == 1 && cad->snapshot() && cad->snapshot()->faces.size() == 6,
            "A current successful STEP import was rejected or failed to replace the source.");
    std::cout << "STEP lifecycle: old-project reply preserved=" << epochPreserved
              << ", old-revision reply preserved=" << revisionPreserved
              << "; current failure/success handled.\n";
    require(epochPreserved && revisionPreserved,
            "An obsolete STEP response overwrote a newer project or parameter state.");
    std::cout << "CPU GUI: obsolete Draft/CAD cancellations preserve current revisions, "
                 "assignments and simulation readiness.\n";
}
} // namespace
int main(int argc, char **argv) {
    QApplication app(argc, argv);
    QCoreApplication::setOrganizationName("HornSculptTests");
    QCoreApplication::setApplicationName("Realtime");
    QTemporaryDir settings;
    QSettings::setDefaultFormat(QSettings::IniFormat);
    QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, settings.path());
    try {
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == "--gmsh-source-only") {
            testGmshImportPreservesLiveLink();
            testGmshSubdomainEndToEnd();
            return 0;
        }
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == "--subdomain-lem-gui-only") {
            testSubdomainGuiArchive(true);
            return 0;
        }
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == "--subdomain-gui-only") {
            testSubdomainGuiArchive();
            return 0;
        }
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == "--draft-preview-only") {
            MainWindow main;
            main.show();
            waitUntil([&] { return main.geometryReady(); });
            auto *length = main.findChild<QDoubleSpinBox *>("length");
            auto *dock = main.findChild<QDockWidget *>("geometryAnalysisDock");
            require(length && dock, "Draft preview controls missing.");
            const auto initial = main.revision();
            waitUntil([&] {
                return dock->windowTitle() ==
                       QStringLiteral("Geometrieanalyse · Rev. %1").arg(initial);
            });
            bool sawStaleAfterPreview = false;
            QObject::connect(dock, &QDockWidget::windowTitleChanged, &main,
                             [&](const QString &title) {
                                 if (main.revision() > initial && main.geometryReady() &&
                                     title.contains(QStringLiteral("veraltet")) &&
                                     title.contains(QStringLiteral("wird aktualisiert")))
                                     sawStaleAfterPreview = true;
                             });
            for (int i = 0; i < 8; ++i) {
                length->setValue(length->value() + 1);
                QTest::qWait(50);
            }
            const auto revision = main.revision();
            require(revision > initial && !main.geometryReady(),
                    "An edit retained a current preview.");
            waitUntil([&] { return main.geometryReady() && main.revision() == revision; });
            waitUntil([&] {
                return dock->windowTitle() ==
                       QStringLiteral("Geometrieanalyse · Rev. %1").arg(revision);
            });
            require(sawStaleAfterPreview,
                    "Old analysis was presented as current while the new preview was ready.");
            QTest::qWait(100);
            require(main.geometryReady() && main.revision() == revision &&
                        dock->windowTitle() ==
                            QStringLiteral("Geometrieanalyse · Rev. %1").arg(revision),
                    "A cancelled analysis overwrote the current revision.");
            QMainWindow host;
            SimulationDock dockCheck(&host);
            auto testMesh = std::make_shared<SurfaceMesh>(testing::sphere(0));
            dockCheck.setMesh(testMesh, 1);
            int changes = 0, frameNotifications = 0;
            dockCheck.onChanged = [&] { ++changes; };
            dockCheck.onObservationFrameChanged = [&] { ++frameNotifications; };
            dockCheck.invalidate(2);
            dockCheck.setMesh(testMesh, 2);
            dockCheck.setAcousticAxis({0, 0, .2}, {0, 0, 1}, true);
            require(changes == 0 && frameNotifications == 1 &&
                        dockCheck.observationFrame().origin.z == .2 &&
                        !dockCheck.hasCurrentResult(),
                    "An already-invalidated source frame lost its notification or dirtied the "
                    "project.");
            dockCheck.setAcousticAxis({0, 0, .3}, {0, 0, 1});
            require(changes == 1 && frameNotifications == 2,
                    "An independent source-frame change failed to invalidate its settings.");
            std::cout << "Draft preview and rapid-edit analysis revisions passed.\n";
            return 0;
        }
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == "--draft-compute-only") {
            MainWindow main;
            main.show();
            waitUntil([&] { return main.geometryReady(); });
            auto *length = main.findChild<QDoubleSpinBox *>("length");
            require(length, "Draft performance control missing.");
            std::vector<qint64> samples;
            for (int i = 0; i < 12; ++i) {
                QElapsedTimer timer;
                timer.start();
                length->setValue(length->value() + (i % 2 ? -1 : 1));
                waitUntil([&] { return main.geometryReady(); });
                samples.push_back(timer.elapsed());
            }
            auto *tabs = main.findChild<QTabWidget *>("draftSections");
            auto *panel = tabs ? dynamic_cast<ParameterPanel *>(tabs->parentWidget()) : nullptr;
            require(panel, "Draft parameter panel missing.");
            const auto snapshot = panel->parameters();
            QElapsedTimer stage;
            stage.start();
            const DraftGeometry geometry(snapshot);
            const auto validationMs = stage.elapsed();
            stage.restart();
            const auto mesh = geometry.previewMesh(MeshSettings{}.adaptive.geometricTolerance);
            const auto tessellationMs = stage.elapsed();
            stage.restart();
            const auto analysis = analyzeGeometry(geometry);
            const auto analysisMs = stage.elapsed();
            require(!mesh.triangles.empty() && !analysis.samples.empty(),
                    "Stage measurements failed.");
            std::cout << "Single-stage CPU (ms): geometry validation=" << validationMs
                      << " tessellation=" << tessellationMs << " analysis=" << analysisMs
                      << " triangles=" << mesh.triangles.size() << '\n';
            std::sort(samples.begin(), samples.end());
            std::cout << "Draft edit to accepted preview, offscreen event loop (ms): p50="
                      << samples[samples.size() / 2] << " p95=" << samples.back()
                      << " n=" << samples.size() << '\n';
            return 0;
        }
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == "--draft-performance-only") {
            MainWindow main;
            main.resize(1280, 800);
            main.show();
            waitUntil([&] { return main.geometryReady(); });
            auto *length = main.findChild<QDoubleSpinBox *>("length");
            auto *viewport = main.findChild<QOpenGLWidget *>("meshViewport");
            require(length && viewport && viewport->isVisible(),
                    "Draft performance controls missing.");
            std::vector<qint64> samples;
            for (int i = 0; i < 12; ++i) {
                bool presented = false;
                QElapsedTimer timer;
                const auto connection =
                    QObject::connect(viewport, &QOpenGLWidget::frameSwapped, &main, [&] {
                        if (main.geometryReady())
                            presented = true;
                    });
                timer.start();
                length->setValue(length->value() + (i % 2 ? -1 : 1));
                waitUntil([&] { return presented; });
                samples.push_back(timer.elapsed());
                QObject::disconnect(connection);
            }
            std::sort(samples.begin(), samples.end());
            std::cout << "Draft edit to visible frame (ms): p50=" << samples[samples.size() / 2]
                      << " p95=" << samples.back() << " n=" << samples.size() << '\n';
            return 0;
        }
        if ((argc == 2 || argc == 3) && QString::fromLocal8Bit(argv[1]) == "--slot-radii-only") {
            testDiffractionSlotRadii(argc == 3 ? QString::fromLocal8Bit(argv[2]) : QString{});
            return 0;
        }
        if ((argc == 2 || argc == 3) && QString::fromLocal8Bit(argv[1]) == "--slot-only") {
            testDiffractionSlotEditing(argc == 3 ? QString::fromLocal8Bit(argv[2]) : QString{});
            return 0;
        }
        if ((argc == 2 || argc == 3) && QString::fromLocal8Bit(argv[1]) == "--usability-only") {
            const auto capture = [&](MainWindow &window, const QString &suffix) {
                if (argc != 3)
                    return;
                QCoreApplication::processEvents();
                require(window.grab().save(QString::fromLocal8Bit(argv[2]) + suffix + ".png"),
                        "Usability screenshot could not be saved.");
            };
            hornsculpt::testing::checkSimulationDockLayout();
            MainWindow main;
            main.resize(1280, 800);
            main.show();
            waitUntil([&] { return main.geometryReady(); });
            const auto revision = main.revision();
            auto *parameters = main.findChild<QDockWidget *>("parameterDock");
            auto *reset = main.findChild<QAction *>("resetDockLayout");
            auto *representation = main.findChild<QComboBox *>("geometryRepresentation");
            require(parameters && reset && representation, "Usability controls are missing.");
            parameters->setFloating(true);
            parameters->hide();
            reset->trigger();
            QCoreApplication::processEvents();
            require(parameters->isVisible() && !parameters->isFloating() &&
                        main.dockWidgetArea(parameters) == Qt::RightDockWidgetArea &&
                        main.revision() == revision,
                    "Resetting the dock layout failed or modified the project.");
            const int display = representation->currentIndex();
            QWheelEvent wheel(QPointF(10, 10), representation->mapToGlobal(QPoint(10, 10)), {},
                              QPoint(0, -120), Qt::NoButton, Qt::NoModifier, Qt::NoScrollPhase,
                              false);
            QApplication::sendEvent(representation, &wheel);
            require(representation->currentIndex() == display + 1 && wheel.isAccepted(),
                    "Wheel cannot change a closed selector.");
            auto *guideEnabled = main.findChild<QCheckBox *>("guideEnabled");
            auto *guidePosition = main.findChild<QDoubleSpinBox *>("guidePosition");
            auto *draftStatus = main.findChild<QLabel *>("draftGeometryStatus");
            auto *draftTabs = main.findChild<QTabWidget *>("draftSections");
            require(guideEnabled && guidePosition && draftStatus && draftTabs,
                    "Guide validation controls are missing.");
            auto *panel = dynamic_cast<ParameterPanel *>(draftTabs->parentWidget());
            require(panel != nullptr, "Draft parameter panel is missing.");
            auto guided = panel->parameters();
            guided.guideEnabled = true;
            guided.guide.position = .05;
            guided.transition.startFraction = .2;
            panel->setParameters(guided);
            panel->onChanged();
            waitUntil([&] { return draftStatus->text().contains("Letzte gültige Vorschau"); });
            require(draftStatus->isVisible() && draftStatus->wordWrap() &&
                        draftStatus->text().contains("Leitkurve") &&
                        draftStatus->text().contains("Vorschau nicht aktuell") &&
                        !main.geometryReady(),
                    "Invalid guide leaves a seemingly current preview or hides its diagnosis.");
            draftTabs->setCurrentIndex(2);
            capture(main, "-leitkurve-ungueltig");
            guidePosition->setValue(100);
            waitUntil([&] { return main.geometryReady(); });
            require(draftStatus->isHidden() && draftStatus->text().isEmpty(),
                    "A corrected guide retains the stale-preview warning.");
            auto *simulation = main.findChild<QDockWidget *>("simulationDock");
            auto *viewport = main.findChild<QWidget *>("meshViewport");
            require(simulation && viewport, "Field layout widgets are missing.");
            capture(main, "-leitkurve");
            auto *results = main.findChild<QTabWidget *>("acousticResults");
            require(results, "Acoustic result tabs are missing.");
            for (const auto &[view, widget] : std::array{
                     std::pair{"frequency", "acousticOnAxis"}, std::pair{"polar", "acousticPolarH"},
                     std::pair{"isobar", "acousticIsobarH"},
                     std::pair{"angular", "angularResponses"}, std::pair{"phase", "phaseResponses"},
                     std::pair{"field", "acousticField"}}) {
                main.setSimulationView(QString::fromLatin1(view));
                require(results->currentWidget()->findChild<QWidget *>(QString::fromLatin1(widget)),
                        "Result navigation selected a different page.");
            }
            main.resize(853, 533);
            reset->trigger();
            QTest::qWait(100);
            require(main.width() <= 853 && main.height() <= 533,
                    "The main window cannot fit a 1280x800 display at 150% scaling.");
            auto *treeDock = main.findChild<QDockWidget *>("projectDock");
            auto *meshDock = main.findChild<QDockWidget *>("meshSettingsDock");
            auto *compactParameters = main.findChild<QDockWidget *>("parameterDock");
            require(treeDock && meshDock && compactParameters &&
                        main.tabifiedDockWidgets(treeDock).contains(meshDock) &&
                        main.tabifiedDockWidgets(treeDock).contains(compactParameters) &&
                        compactParameters->width() >= 300,
                    "The compact layout did not give the side controls a shared readable dock.");
            require(main.findChild<QScrollArea *>("parameterDockScroll") &&
                        main.findChild<QScrollArea *>("meshSettingsDockScroll"),
                    "Narrow side docks cannot scroll their full controls.");
            auto *enlargeResults = main.findChild<QAction *>("enlargeResults");
            require(enlargeResults, "Result enlargement action is missing.");
            main.setSimulationView("field");
            QTest::qWait(100);
            auto *field = main.findChild<QWidget *>("acousticField");
            require(field, "Field controls are missing.");
            for (const auto *name :
                 {"fieldPlane", "fieldFrequency", "fieldExcitation", "fieldQuantity"}) {
                auto *control = field->findChild<QWidget *>(name);
                const auto labels = field->findChildren<QLabel *>();
                require(control && std::any_of(labels.begin(), labels.end(),
                                               [control](QLabel *label) {
                                                   return label->buddy() == control &&
                                                          !label->text().isEmpty();
                                               }),
                        "A field selector lacks its visible label and keyboard association.");
            }
            const auto beforeEnlargement = main.revision();
            const auto previousDockHeight = simulation->height();
            enlargeResults->trigger();
            QTest::qWait(100);
            require(simulation->isVisible() && simulation->height() > previousDockHeight &&
                        main.revision() == beforeEnlargement && !main.simulationInProgress(),
                    "Result enlargement failed or changed the project state.");
            main.setSimulationView("field");
            QTest::qWait(100);
            require(viewport->isVisible() && !enlargeResults->isChecked() &&
                        main.revision() == beforeEnlargement,
                    "Selecting the field did not restore the 3D view.");
            enlargeResults->trigger();
            auto *simulationControls = dynamic_cast<SimulationDock *>(simulation);
            require(simulationControls && simulationControls->onFocusField,
                    "Field camera action is missing.");
            simulationControls->onFocusField();
            QTest::qWait(100);
            require(viewport->isVisible() && !enlargeResults->isChecked(),
                    "Focusing the field did not restore the 3D view.");
            enlargeResults->trigger();
            auto plane = std::make_shared<PressureFieldPlane>();
            plane->columns = plane->rows = 2;
            plane->minimumLateral = plane->minimumAxial = -1;
            plane->maximumLateral = plane->maximumAxial = 1;
            plane->frequency = 1000;
            plane->positions = {{-1, 0, -1}, {1, 0, -1}, {-1, 0, 1}, {1, 0, 1}};
            plane->pressure.assign(4, {1, 0});
            plane->valid.assign(4, true);
            require(bool(simulationControls->onFieldChanged),
                    "Field publication callback missing.");
            simulationControls->onFieldChanged(plane, {});
            QTest::qWait(100);
            require(viewport->isVisible() && !enlargeResults->isChecked(),
                    "A completed field stayed hidden behind enlarged results.");
            enlargeResults->trigger();
            QTest::qWait(100);
            reset->trigger();
            QTest::qWait(100);
            require(viewport->isVisible() && !enlargeResults->isChecked() &&
                        main.revision() == beforeEnlargement,
                    "Resetting the layout did not restore the enlarged viewport.");
            const auto displayRevision = main.revision();
            for (const auto resolution : {QSize(1480, 980), QSize(1050, 720)}) {
                main.resize(resolution);
                main.setSimulationView("field");
                try {
                    waitUntil([&] {
                        return simulation->height() <= main.height() * .42 &&
                               viewport->height() >= main.height() * .40;
                    });
                } catch (const std::exception &) {
                    std::cerr << "Field layout " << main.width() << 'x' << main.height()
                              << ": viewport=" << viewport->height()
                              << ", simulation=" << simulation->height() << '\n';
                    throw;
                }
                require(
                    viewport->height() >= main.height() * .40 &&
                        main.revision() == displayRevision && !main.simulationInProgress(),
                    "Inactive geometry/LEM tab minimums still squeeze the upper 3D field view.");
            }
            // The synthetic snapshot exercises layout only; no simulated result is asserted.
            auto *dock = dynamic_cast<SimulationDock *>(simulation);
            require(dock, "Simulation dock is missing.");
            const auto surface = std::make_shared<SurfaceMesh>(
                DraftGeometry(panel->parameters()).previewMesh(.00025));
            const std::vector<FaceVelocity> velocities{{5, {1, 0}}};
            dock->setMesh(surface, main.revision(), velocities);
            auto archive = dock->archive();
            archive.settings.minimumFrequency = archive.settings.maximumFrequency = 1000;
            archive.settings.frequencyCount = 1;
            archive.mesh = surface;
            const auto frame = dock->observationFrame();
            archive.request =
                makeSimulationRequest(archive.settings, velocities, frame.origin, frame.forward);
            auto result = std::make_shared<BemResult>();
            result->device = "SYNTHETIC GUI LAYOUT - no GPU computation";
            result->triangleCount = surface->triangles.size();
            result->pressureDofCount = surface->vertices.size();
            BemExcitationResult excitation;
            excitation.name = archive.request->excitations.front().name;
            BemBoundaryData boundary;
            boundary.surfacePressure.assign(result->triangleCount, {.01, 0});
            excitation.setBoundary(std::move(boundary));
            excitation.observerPressure.assign(archive.request->observers.size(), {.01, 0});
            result->frequencies.push_back({1000, {std::move(excitation)}});
            archive.result = result;
            archive.current = true;
            archive.revision = main.revision();
            dock->restoreArchive(archive);
            auto *detailsToggle = main.findChild<QToolButton *>("fieldDetailsToggle");
            require(detailsToggle, "Field details toggle missing.");
            for (const auto resolution :
                 {QSize(853, 533), QSize(1050, 720), QSize(1280, 800), QSize(1480, 980)}) {
                main.resize(resolution);
                reset->trigger();
                main.setSimulationView("field");
                QTest::qWait(100);
                if (main.size() != resolution) {
                    std::cerr << "Requested " << resolution.width() << 'x' << resolution.height()
                              << ", window " << main.width() << 'x' << main.height() << '\n';
                    throw std::runtime_error("Field layout exceeds the requested window size.");
                }
                if (resolution == QSize(853, 533)) {
                    auto *compactTools = main.findChild<QToolButton *>("compactViewTools");
                    auto *viewTools = main.findChild<QScrollArea *>("viewControlsScroll");
                    require(compactTools && compactTools->isVisible() &&
                                !compactTools->isChecked() && !compactTools->text().isEmpty() &&
                                viewTools && !viewTools->isVisible(),
                            "Compact toolbars are not available in the small window.");
                    compactTools->click();
                    require(viewTools->isVisible(), "Compact view controls cannot be reopened.");
                    auto *toolsDialog = main.findChild<QDialog *>("compactToolsDialog");
                    require(toolsDialog && toolsDialog->isVisible(),
                            "Compact tools dialog is missing.");
                    toolsDialog->close();
                    require(!compactTools->isChecked() && !viewTools->isVisible(),
                            "Closing compact tools left their state open.");
                    compactTools->click();
                    QTest::qWait(100);
                    for (const auto *name :
                         {"fieldPlane", "fieldFrequency", "fieldExcitation", "fieldQuantity",
                          "evaluateField", "cancelField", "exportField", "focusField3D",
                          "fieldStatus", "simulationStatus"}) {
                        auto *control = main.findChild<QWidget *>(name);
                        if (control &&
                            (!control->isVisible() ||
                             control->height() < control->minimumSizeHint().height() ||
                             !simulation->rect().contains(
                                 QRect(control->mapTo(simulation, QPoint()), control->size()))))
                            std::cerr << "Expanded tools: " << name << " at "
                                      << control->mapTo(simulation, QPoint()).y() << " height "
                                      << control->height() << " min "
                                      << control->minimumSizeHint().height() << " dock "
                                      << simulation->height() << '\n';
                        require(control && control->isVisible() &&
                                    control->height() >= control->minimumSizeHint().height() &&
                                    control->visibleRegion().contains(control->rect()) &&
                                    simulation->rect().contains(QRect(
                                        control->mapTo(simulation, QPoint()), control->size())),
                                "Expanded compact tools squeeze field controls or status.");
                        const QRect controlRect(control->mapTo(&main, QPoint()), control->size());
                        for (auto *bar : main.findChildren<QTabBar *>()) {
                            bool dockTabs = false;
                            for (int i = 0; i < bar->count(); ++i)
                                dockTabs = dockTabs || bar->tabText(i).contains("Akustik");
                            if (dockTabs && bar->isVisible())
                                require(!controlRect.intersects(
                                            QRect(bar->mapTo(&main, QPoint()), bar->size())),
                                        "Expanded compact tools overlap the dock tabs.");
                        }
                    }
                    compactTools->click();
                    QTest::qWait(100);
                }
                for (const bool expanded : {false, true}) {
                    auto *detailsDialog = main.findChild<QDialog *>("fieldDetailsDialog");
                    require(detailsDialog, "Field settings dialog missing.");
                    if (expanded)
                        detailsToggle->click();
                    else
                        detailsDialog->close();
                    QTest::qWait(100);
                    auto *fieldPage = main.findChild<QWidget *>("acousticField");
                    require(fieldPage, "Field page missing.");
                    for (const auto *name :
                         {"fieldPlane", "fieldFrequency", "fieldExcitation", "fieldQuantity",
                          "evaluateField", "cancelField", "exportField", "focusField3D",
                          "fieldDetailsToggle", "fieldProbeMode", "fieldFrequencySlider",
                          "fieldPreviousFrequency", "fieldNextFrequency", "fieldStatus"}) {
                        auto *control = main.findChild<QWidget *>(name);
                        if (!control || !control->isVisible() ||
                            !control->visibleRegion().contains(control->rect()) ||
                            control->height() < control->minimumSizeHint().height() ||
                            !simulation->rect().contains(
                                QRect(control->mapTo(simulation, QPoint()), control->size())) ||
                            !fieldPage->rect().contains(
                                QRect(control->mapTo(fieldPage, QPoint()), control->size()))) {
                            std::cerr
                                << resolution.width() << 'x' << resolution.height()
                                << ": clipped field control " << name << ", details=" << expanded
                                << ", height=" << (control ? control->height() : -1)
                                << ", minimum hint="
                                << (control ? control->minimumSizeHint().height() : -1)
                                << ", field width=" << fieldPage->width() << ", x="
                                << (control ? control->mapTo(fieldPage, QPoint()).x() : -1)
                                << ", width=" << (control ? control->width() : -1)
                                << ", widthHint=" << (control ? control->sizeHint().width() : -1)
                                << ", widthMin="
                                << (control ? control->minimumSizeHint().width() : -1) << '\n';
                            throw std::runtime_error("Field controls are not fully visible.");
                        }
                    }
                    auto *details = main.findChild<QScrollArea *>("fieldDetailsScroll");
                    auto *solverStatus = main.findChild<QLabel *>("simulationStatus");
                    require(solverStatus && (!solverStatus->isVisible() ||
                                             simulation->rect().contains(
                                                 QRect(solverStatus->mapTo(simulation, QPoint()),
                                                       solverStatus->size()))),
                            "Solver status extends outside its dock.");
                    require(details && (!expanded || (details->isVisible() &&
                                                      details->viewport()->height() >= 40)),
                            "Expanded field settings cannot be scrolled.");
                    if (expanded)
                        detailsDialog->close();
                }
                require(viewport->isVisible() && viewport->height() >= 100,
                        "Field controls leave no useful 3D view.");
            }
            auto *compactTools = main.findChild<QToolButton *>("compactViewTools");
            auto *viewTools = main.findChild<QScrollArea *>("viewControlsScroll");
            main.setSimulationView("frequency");
            for (int i = 0; i < results->count(); ++i)
                if (results->widget(i)->findChild<QWidget *>("cta2034Levels"))
                    results->setCurrentIndex(i);
            enlargeResults->trigger();
            main.setSimulationView("field");
            main.resize(853, 533);
            QTest::qWait(100);
            require(compactTools && compactTools->isVisible() && !compactTools->isChecked(),
                    "Shrinking the window did not enter compact toolbar mode.");
            compactTools->click();
            require(viewTools && viewTools->isVisible(), "Compact toolbar cannot be expanded.");
            main.resize(1050, 720);
            QTest::qWait(100);
            require(!compactTools->isVisible() && viewTools->isVisible(),
                    "Growing the window did not restore the regular toolbars.");
            main.resize(853, 533);
            QTest::qWait(100);
            require(compactTools->isVisible() && compactTools->isChecked() &&
                        viewTools->isVisible(),
                    "Resizing lost the user's expanded toolbar choice.");
            for (const auto *name : {"fieldPlane", "fieldFrequency", "fieldExcitation",
                                     "fieldQuantity", "evaluateField", "cancelField", "exportField",
                                     "focusField3D", "fieldStatus", "simulationStatus"}) {
                auto *control = main.findChild<QWidget *>(name);
                if (control && (!control->isVisible() ||
                                control->height() < control->minimumSizeHint().height() ||
                                !simulation->rect().contains(
                                    QRect(control->mapTo(simulation, QPoint()), control->size()))))
                    std::cerr << "After resize: " << name << " at "
                              << control->mapTo(simulation, QPoint()).y() << " height "
                              << control->height() << " dock " << simulation->height() << '\n';
                require(control && control->isVisible() &&
                            control->height() >= control->minimumSizeHint().height() &&
                            simulation->rect().contains(
                                QRect(control->mapTo(simulation, QPoint()), control->size())),
                        "Resizing with open compact tools squeezed field controls or status.");
                const QRect controlRect(control->mapTo(&main, QPoint()), control->size());
                require(!controlRect.intersects(QRect(main.statusBar()->mapTo(&main, QPoint()),
                                                      main.statusBar()->size())),
                        "Resizing with open compact tools overlaps the main status bar.");
                for (auto *bar : main.findChildren<QTabBar *>()) {
                    bool dockTabs = false;
                    for (int i = 0; i < bar->count(); ++i)
                        dockTabs = dockTabs || bar->tabText(i).contains("Akustik");
                    if (dockTabs && bar->isVisible())
                        require(!controlRect.intersects(
                                    QRect(bar->mapTo(&main, QPoint()), bar->size())),
                                "Resizing with open compact tools overlaps the dock tabs.");
                }
            }
            compactTools->click();
            std::cout
                << "Dock layout, wheel editing, guide diagnostics and result navigation passed.\n";
            return 0;
        }
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == "--assembly-only") {
            testDraftAssemblyEditing();
            return 0;
        }
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == "--axisymmetric-controls-only") {
            testAxisymmetricControls();
            return 0;
        }
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == "--axisymmetric-only") {
            testAxisymmetricSimulation();
            return 0;
        }
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == "--angular-controls-only") {
            testAngularResponseControls();
            return 0;
        }
        if (argc == 2 &&
            QString::fromLocal8Bit(argv[1]) == QStringLiteral("--baffle-controls-only")) {
            testBaffleControls();
            return 0;
        }
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == "--interface-controls-only") {
            testInterfaceControls();
            return 0;
        }
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == QStringLiteral("--lem-only")) {
            testLemIntegration();
            return 0;
        }
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == QStringLiteral("--orientation-only")) {
            testSavePendingEdits();
            testOrientationEdits();
            return 0;
        }
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == QStringLiteral("--lifecycle-only")) {
            testJobLifecycle();
            testAsyncSaveSnapshotAndClose();
            testAsyncSaveFailureKeepsProject();
            testAsyncSaveIncomingResult();
            testCloseSupersedesPendingNew();
            testAsyncSaveAsAndRepeat();
            testSaveWorkerDestruction();
            testSaveWorkerRealIoFailureAndRecovery();
            testCadEditIdentity();
            testGmshImportPreservesLiveLink();
            testSourceImportStartFailure();
            testAsyncProjectOpen();
            return 0;
        }
        QMainWindow window;
        SimulationDock dock(&window);
        const auto defaults = dock.archive();
        auto *fastMultipole = dock.findChild<QCheckBox *>("fastMultipoleEnabled");
        auto *multipoleLeafCapacity = dock.findChild<QSpinBox *>("fastMultipoleLeafCapacity");
        require(
            fastMultipole && multipoleLeafCapacity && !fastMultipole->isChecked() &&
                !multipoleLeafCapacity->isEnabled() &&
                defaults.settings.acceleration == BemAcceleration::dense,
            "Fast Multipole is not explicitly selectable or defaults to an unexpected backend.");
        unsigned solverChanges = 0;
        dock.onChanged = [&] { ++solverChanges; };
        fastMultipole->setChecked(true);
        multipoleLeafCapacity->setValue(64);
        const auto multipoleArchive = dock.archive();
        require(solverChanges == 2 && multipoleLeafCapacity->isEnabled() &&
                    multipoleArchive.settings.acceleration == BemAcceleration::fastMultipole &&
                    multipoleArchive.settings.fastMultipole.leafSupportCapacity == 64,
                "Solver controls did not notify the revision path or archive their state.");
        dock.restoreArchive(defaults);
        require(!fastMultipole->isChecked() && !multipoleLeafCapacity->isEnabled(),
                "Restoring the dense backend did not restore its controls.");
        dock.restoreArchive(multipoleArchive);
        require(fastMultipole->isChecked() && multipoleLeafCapacity->isEnabled() &&
                    multipoleLeafCapacity->value() == 64,
                "Fast Multipole configuration did not survive UI restoration.");
        dock.restoreArchive(defaults);
        auto *nearDepth = dock.findChild<QSpinBox *>("nearMaximumDepth");
        require(nearDepth && nearDepth->isEnabled() && nearDepth->minimum() == 0 &&
                    nearDepth->maximum() == 10 && nearDepth->value() == 6 &&
                    defaults.settings.nearMaximumDepth == 6,
                "Near quadrature depth is unavailable or has an unexpected default or range.");
        const auto changesBeforeDepth = solverChanges;
        nearDepth->setValue(8);
        const auto deeperQuadrature = dock.archive();
        require(solverChanges == changesBeforeDepth + 1 &&
                    deeperQuadrature.settings.nearMaximumDepth == 8,
                "The near quadrature control did not notify the revision path or save its value.");
        dock.restoreArchive(defaults);
        require(nearDepth->value() == 6, "The default near quadrature depth was not restored.");
        dock.restoreArchive(deeperQuadrature);
        require(nearDepth->value() == 8 && dock.archive().settings.nearMaximumDepth == 8,
                "The selected near quadrature depth did not survive UI restoration.");
        auto invalidDepth = deeperQuadrature;
        invalidDepth.settings.nearMaximumDepth = 11;
        bool rejectedDepth = false;
        try {
            dock.restoreArchive(invalidDepth);
        } catch (const std::invalid_argument &) {
            rejectedDepth = true;
        }
        require(rejectedDepth && nearDepth->value() == 8,
                "An invalid near quadrature depth was silently clamped during restoration.");
        dock.restoreArchive(defaults);
        dock.onChanged = {};
        if (argc == 2 && QString::fromLocal8Bit(argv[1]) == QStringLiteral("--fmm-controls-only")) {
            std::cout << "Solver controls, change notification and archive restoration passed.\n";
            return 0;
        }
        require(defaults.settings.directivityMaximumDb == 3,
                "The directivity scale clips off-axis gain by default.");
        auto *minimum = dock.findChild<QLineEdit *>("directivityMinimumDb");
        auto *step = dock.findChild<QLineEdit *>("directivityStepDb");
        auto *referenceAngle = dock.findChild<QDoubleSpinBox *>("polarNormalizationAngle");
        auto *polar = dynamic_cast<PolarPlotWidget *>(dock.findChild<QWidget *>("acousticPolarH"));
        require(minimum && step && referenceAngle && polar,
                "Directivity scaling and angle controls are unavailable.");
        minimum->setText("-42");
        step->setText("0,125");
        QMetaObject::invokeMethod(step, "editingFinished", Qt::DirectConnection);
        referenceAngle->setValue(23.5);
        auto display = dock.archive();
        require(display.settings.directivityMinimumDb == -42 &&
                    display.settings.directivityStepDb == .125 &&
                    display.settings.normalizationAngleDegrees == 23.5 &&
                    polar->decibelScale() == PlotDbScale{-42, 3, .125},
                "Shared directivity controls do not reach plots and project state.");
        polar->onDecibelScaleChanged({-.000002, .000003, .0000001});
        display = dock.archive();
        require(display.settings.directivityStepDb == .0000001 &&
                    display.settings.directivityMaximumDb == .000003,
                "A confirmed plot scale was silently rounded by shared controls.");
        dock.restoreArchive(display);
        require(dock.archive().settings.directivityStepDb == .0000001,
                "A fine dB scale cannot be restored.");
        dock.restoreArchive(defaults);
        auto mesh = std::make_shared<SurfaceMesh>(testing::sphere(2));
        dock.setMesh(mesh, 1);
        dock.setAcousticAxis({0, 0, .1}, {0, 0, 1});
        dock.setFrequencySweep(200, 800, 3);
        std::vector<double> arrival;
        std::vector<std::uint64_t> revisions;
        std::set<double> arrived;
        bool retainedDisplay = false;
        bool partialSaved = false;
        dock.onChanged = [&] {
            if (!dock.result() || !dock.isBusy())
                return;
            const auto archive = dock.archive();
            if (archive.resultComplete != std::optional{false})
                return;
            require(!archive.current, "A partial raster was marked complete.");
            require(archive.result->frequencies.size() == archive.request->frequencies.size(),
                    "Partial results and request do not match.");
            require(std::is_sorted(archive.request->frequencies.begin(),
                                   archive.request->frequencies.end()),
                    "Partial frequencies are not sorted for plots and persistence.");
            if (archive.result->frequencies.size() == 1)
                arrived.clear();
            for (const auto &computed : archive.result->frequencies)
                if (arrived.insert(computed.frequency).second)
                    arrival.push_back(computed.frequency);
            revisions.push_back(archive.revision);
            if (archive.request->frequencies.back() == 1200) {
                const auto *selector = dock.findChild<QComboBox *>("resultFrequency");
                require(selector &&
                            selector->count() ==
                                static_cast<int>(3 + archive.result->frequencies.size()) &&
                            selector->findData(400.0) >= 0 &&
                            selector->currentData().toDouble() == 200,
                        "New frequencies erased the old/interpolated display or selection.");
                retainedDisplay = true;
            }
            if (!partialSaved) {
                QTemporaryDir directory;
                Project project;
                project.revision = 1;
                project.simulation = archive;
                const auto path = directory.filePath("partial.hsproj");
                writeProject(path, project);
                require(readProject(path).simulation.result->frequencies.size() ==
                            archive.result->frequencies.size(),
                        "A streamed partial raster did not survive project save/open.");
                partialSaved = true;
            }
        };
        auto *button = dock.findChild<QPushButton *>("realtimeSolve");
        require(button && button->isCheckable(), "Realtime toggle is missing.");
        button->click();
        waitUntil([&] { return dock.hasCurrentResult() || !dock.lastError().isEmpty(); });
        require(dock.hasCurrentResult(), "Realtime solve failed.");
        require(arrived == std::set<double>({200, 400, 800}) && partialSaved,
                "The final coalesced partial raster lost a computed frequency or persistence.");
        require(frequencySchedule(9) == std::vector<std::size_t>({8, 0, 4, 6, 2, 7, 5, 3, 1}),
                "Frequency gaps were not refined from coarse to fine with high-gap priority.");
        for (std::size_t count = 0; count <= 82; ++count) {
            auto order = frequencySchedule(count);
            std::sort(order.begin(), order.end());
            require(order.size() == count, "The progressive sweep lost frequencies.");
            for (std::size_t i = 0; i < count; ++i)
                require(order[i] == i, "The progressive sweep duplicates or omits a frequency.");
        }
        const auto completed = dock.result();
        QTest::qWait(650);
        require(dock.result() == completed && !dock.isBusy() && dock.realtimeEnabled(),
                "An unchanged completed raster was needlessly repeated.");
        dock.beginExternalJob();
        waitUntil([&] { return dock.externalJobReady(); });
        require(dock.hasCurrentResult() && dock.result() == completed && !dock.realtimeEnabled(),
                "Preparing a study discarded an unchanged completed result.");
        dock.requestSolve();
        require(!dock.isBusy(), "A desktop solve competed with the external study GPU owner.");
        dock.endExternalJob();
        dock.setRealtimeEnabled(true);

        auto *selector = dock.findChild<QComboBox *>("resultFrequency");
        selector->setCurrentIndex(selector->findData(200.0));
        dock.setAcousticAxis({0, 0, .12}, {0, 0, 1});
        dock.setFrequencySweep(250, 1000, 3);
        dock.setFrequencySweep(300, 1200, 3);
        waitUntil([&] { return dock.hasCurrentResult() || !dock.lastError().isEmpty(); });
        require(dock.hasCurrentResult() && dock.result()->frequencies.back().frequency == 1200 &&
                    retainedDisplay,
                "Edits did not automatically replace the previous sweep.");

        dock.setFrequencySweep(350, 1400, 4);
        waitUntil([&] { return dock.isBusy(); });
        const auto previousSolveId = dock.activeSolveJobId();
        require(previousSolveId != 0, "Active solver job has no identity.");
        auto *activeProgress = dock.findChild<QProgressBar *>("simulationProgress");
        require(activeProgress && !activeProgress->isHidden(),
                "An active simulation hid its progress indicator.");
        dock.invalidate(2);
        dock.setMesh(mesh, 2);
        dock.setAcousticAxis({0, 0, .1}, {0, 0, 1});
        waitUntil([&] { return dock.hasCurrentResult() || !dock.lastError().isEmpty(); });
        require(dock.hasCurrentResult() && dock.resultRevision() == 2,
                "Revision changes did not cancel and replace the obsolete job.");
        const auto current = dock.result();
        dock.setRealtimeEnabled(false);
        dock.setFrequencySweep(400, 1600, 4);
        QTest::qWait(650);
        require(!dock.isBusy() && dock.result() == current,
                "Disabled realtime mode continued scheduling jobs.");
        CompletionFilter solveCompletion;
        bool stoppedQueued = false;
        std::shared_ptr<const BemResult> retainedPartial;
        solveCompletion.isCompletion = [&] {
            const auto result = dock.result();
            return dock.isBusy() && result && result->frequencies.size() == 4 &&
                   result->frequencies.front().frequency == 400 &&
                   result->frequencies.back().frequency == 1600;
        };
        solveCompletion.beforeCompletion = [&] {
            retainedPartial = dock.result();
            for (auto *stop : dock.findChildren<QPushButton *>())
                if (stop->text() == QStringLiteral("Abbrechen") && stop->isEnabled()) {
                    stop->click();
                    stoppedQueued = true;
                }
        };
        dock.installEventFilter(&solveCompletion);
        dock.setRealtimeEnabled(true);
        waitUntil([&] { return dock.isBusy(); });
        require(dock.activeSolveJobId() != 0 && dock.activeSolveJobId() != previousSolveId,
                "A replacement solve reused the previous solver job identity.");
        waitUntil([&] { return !dock.isBusy(); });
        require(dock.activeSolveJobId() == 0, "An idle solver retained an active job identity.");
        dock.removeEventFilter(&solveCompletion);
        const auto stopped = dock.result();
        require(stoppedQueued && dock.lastError() == QStringLiteral("Berechnung abgebrochen") &&
                    stopped == retainedPartial && !dock.hasCurrentResult(),
                "Explicit stop after queued CUDA completion committed the cancelled result.");
        QTest::qWait(400);
        require(activeProgress->isHidden(),
                "Canceled simulation left a progress indicator visible.");
        require(!dock.isBusy() && dock.result() == stopped && !dock.realtimeEnabled(),
                "Toggling off during a solve did not stop and stay stopped.");
        dock.beginExternalJob();
        waitUntil([&] { return dock.externalJobReady(); });
        dock.endExternalJob();

        dock.onChanged = {};
        const AcousticFrame rotated{{.2, -.1, .3}, {0, -1, 0}, {1, 0, 0}};
        const auto previousResult = dock.result();
        dock.setObservationFrame(rotated);
        require(!dock.hasCurrentResult() && dock.result() == previousResult,
                "Changing the observation frame erased the retained result or left it current.");
        dock.setAcousticAxis({0, 0, .5}, {1, 0, 0});
        require(dock.observationFrame() == rotated,
                "A source revision replaced the user's observation frame.");
        dock.setFrequencySweep(500, 1000, 2);
        bool cleanFrameDisplay = false;
        dock.onChanged = [&] {
            const auto archive = dock.archive();
            if (!dock.isBusy() || !archive.request || !archive.request->farField ||
                archive.request->farField->origin != rotated.origin ||
                archive.resultComplete != std::optional{false})
                return;
            const auto *frequencies = dock.findChild<QComboBox *>("resultFrequency");
            cleanFrameDisplay =
                frequencies &&
                frequencies->count() == static_cast<int>(archive.request->frequencies.size());
        };
        dock.setRealtimeEnabled(true);
        dock.requestSolve();
        waitUntil([&] { return dock.hasCurrentResult() || !dock.lastError().isEmpty(); });
        require(dock.hasCurrentResult() && cleanFrameDisplay,
                "A changed frame failed to solve or mixed incompatible H/V display samples.");
        dock.setRealtimeEnabled(false);
        dock.onChanged = {};
        const auto rotatedArchive = dock.archive();
        require(rotatedArchive.mesh == mesh &&
                    rotatedArchive.settings.observationFrame == rotated &&
                    rotatedArchive.request->farField->origin == rotated.origin &&
                    std::abs(rotatedArchive.request->observers.front().y + 1.1) < 1e-12,
                "The CUDA request did not use the custom frame on the existing mesh.");
        QTemporaryDir frameDirectory;
        Project frameProject;
        frameProject.revision = rotatedArchive.revision;
        frameProject.geometryRevision = frameProject.meshRevision = frameProject.revision;
        frameProject.preview = mesh;
        frameProject.mesh =
            std::make_shared<MeshingResult>(MeshingResult{*mesh, inspectMesh(*mesh)});
        frameProject.simulation = rotatedArchive;
        const auto framePath = frameDirectory.filePath("rotated.hsproj");
        writeProject(framePath, frameProject);
        dock.restoreArchive(readProject(framePath).simulation);
        require(dock.hasCurrentResult() && dock.observationFrame() == rotated,
                "Restored custom-frame results were incorrectly marked stale.");
        auto changedDomain = rotatedArchive;
        changedDomain.settings.infiniteBaffle = InfiniteBaffle{{}, {0, 0, 1}, {1}};
        dock.restoreArchive(changedDomain);
        require(!dock.hasCurrentResult() && dock.result() == rotatedArchive.result &&
                    !dock.archive().request->infiniteBaffle,
                "A changed baffle domain kept old free-field results current or rewrote them.");
        dock.restoreArchive(rotatedArchive);
        require(dock.hasCurrentResult(),
                "Restoring the actual computed domain failed to restore result validity.");
        dock.setObservationFrame({});
        require(dock.observationFrame() == defaultAcousticFrame({0, 0, .5}, {1, 0, 0}) &&
                    !dock.hasCurrentResult(),
                "Resetting observation settings did not restore the current source frame.");
        dock.beginExternalJob();
        waitUntil([&] { return dock.externalJobReady(); });
        dock.endExternalJob();

        MainWindow main;
        for (const int importVersion : {5, 6}) {
            QTemporaryDir imported;
            Project project;
            auto oldSource = std::make_shared<SurfaceMesh>(
                DraftGeometry({}).previewMesh(importVersion == 5 ? .00025 : .000125));
            if (importVersion == 5)
                for (auto &triangle : oldSource->triangles)
                    if (triangle.faceId == 1)
                        triangle.faceId = 6;
            project.preview = oldSource;
            project.geometryRevision = project.revision;
            auto oldMesh = std::make_shared<MeshingResult>();
            oldMesh->mesh = *oldSource;
            oldMesh->quality = inspectMesh(*oldSource);
            project.mesh = oldMesh;
            project.meshRevision = project.revision;
            const auto path = imported.filePath("retired-draft.hsproj");
            writeProject(path, project);
            QFile file(path);
            require(file.open(QIODevice::ReadOnly), "Cannot read the migration fixture.");
            auto json = QJsonDocument::fromJson(file.readAll()).object();
            file.close();
            json.insert("version", importVersion);
            json.remove("slot");
            json.remove("draftReversedFaces");
            json.remove("cadReversedFaces");
            if (importVersion == 5) {
                auto sizing = json["meshing"].toObject();
                sizing.insert("entryAngleToleranceDegrees", .5);
                json.insert("meshing", sizing);
            }
            require(file.open(QIODevice::WriteOnly | QIODevice::Truncate),
                    "Cannot write the migration fixture.");
            require(file.write(QJsonDocument(json).toJson()) > 0,
                    "Cannot store the migration fixture.");
            file.close();
            // Verify the migration fixture before exercising asynchronous GUI adoption.
            const auto legacy = readProject(path);
            require(legacy.draftReversedFaceIds.empty() && legacy.cadReversedFaceIds.empty(),
                    "Older projects do not default to their original face orientations.");
            main.openProject(path);
            waitUntil([&] { return !main.sourceImportInProgress(); });
            require(main.sourceImportError().isEmpty(), "Draft migration import failed.");
            require(main.geometryReady() && main.revision() == project.revision + 1 &&
                        !main.meshReady(),
                    "Imported Draft was not retessellated with its old mesh marked stale.");
            main.exerciseLoadedProjectRoundTrip();
        }
        const auto set = [&](const char *name, double value) {
            auto *input = main.findChild<QDoubleSpinBox *>(name);
            require(input, "A Draft input is missing.");
            input->setValue(value);
        };
        set("entryWidth", 4);
        set("entryHeight", 4);
        set("length", 20);
        set("thickness", 1);
        set("angleH", 20);
        set("angleV", 15);
        main.setSimulationFrequencySweep(400, 800, 3);
        auto *live =
            dynamic_cast<SimulationDock *>(main.findChild<QDockWidget *>("simulationDock"));
        require(live, "Main-window realtime controls are missing.");
        live->setRealtimeEnabled(true);
        waitUntil([&] { return main.simulationReady(); });
        require(main.meshReady() && live->resultRevision() == main.revision(),
                "Enabling realtime did not automatically mesh and solve Draft.");
        const auto previousRevision = main.revision();
        set("length", 21);
        waitUntil([&] { return main.simulationReady(); });
        require(main.revision() > previousRevision && main.meshReady() &&
                    live->resultRevision() == main.revision(),
                "A geometry edit did not automatically remesh and refresh the solve.");
        live->setRealtimeEnabled(false);
        const auto retainedResult = live->result();
        const auto retainedRevision = main.revision();
        QAction *cancelMeshing = nullptr;
        for (auto *action : main.findChildren<QAction *>())
            if (action->objectName() == QStringLiteral("cancelMeshing"))
                cancelMeshing = action;
        require(cancelMeshing, "The mesh cancellation action is missing.");
        main.requestMesh();
        require(main.meshInProgress(), "Repeated meshing did not start.");
        cancelMeshing->trigger();
        waitUntil([&] { return !main.meshInProgress(); });
        require(main.lastMeshError() == QStringLiteral("Vernetzung abgebrochen") &&
                    main.meshReady() && main.revision() == retainedRevision &&
                    main.simulationReady() && live->result() == retainedResult,
                "Cancelling an unchanged mesh rebuild invalidated its current CUDA result.");
        CompletionFilter meshCompletion;
        bool cancelledQueued = false;
        meshCompletion.beforeCompletion = [&] {
            require(main.meshInProgress() && cancelMeshing->isEnabled(),
                    "The queued mesh reply retired before cancellation.");
            cancelMeshing->trigger();
            cancelledQueued = true;
        };
        main.installEventFilter(&meshCompletion);
        main.requestMesh();
        waitUntil([&] { return !main.meshInProgress(); });
        main.removeEventFilter(&meshCompletion);
        require(cancelledQueued &&
                    main.lastMeshError() == QStringLiteral("Vernetzung abgebrochen") &&
                    main.meshReady() && main.revision() == retainedRevision &&
                    main.simulationReady() && live->result() == retainedResult,
                "A queued cancelled mesh reply replaced the current CUDA result.");
        std::cout << "Real CUDA: progressive streaming, retained/interpolated display, partial "
                     "save/open, edit coalescing, "
                     "revision replacement, idle completion, toggle stop and automatic Draft "
                     "remeshing with result-preserving cancellation passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
