#include "../apps/desktop/parameter_sweep.hpp"
#include "../apps/desktop/sweep_model.hpp"
#include "cad_fixture.hpp"
#include "fixture_geometry.hpp"

#include <QApplication>
#include <QComboBox>
#include <QDir>
#include <QDirIterator>
#include <QDoubleSpinBox>
#include <QElapsedTimer>
#include <QEvent>
#include <QEventLoop>
#include <QFileInfo>
#include <QItemSelectionModel>
#include <QLineEdit>
#include <QLabel>
#include <QCheckBox>
#include <QMessageBox>
#include <QPointer>
#include <QProgressBar>
#include <QScrollArea>
#include <QScrollBar>
#include <QScreen>
#include <QPushButton>
#include <QSpinBox>
#include <QStringList>
#include <QTableView>
#include <QTemporaryDir>
#include <QTimer>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <hornsculpt/acoustic_analysis.hpp>
#include <hornsculpt/cad_geometry.hpp>
#include <hornsculpt/heatmap_widget.hpp>
#include <hornsculpt/mesh_orientation.hpp>
#include <hornsculpt/meshing.hpp>
#include <hornsculpt/plot_widget.hpp>
#include <hornsculpt/polar_plot_widget.hpp>
#include <iostream>
#include <memory>
#include <mutex>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace hornsculpt;
namespace {
std::mutex messageMutex;
QStringList lifetimeWarnings;
void qtMessages(QtMsgType, const QMessageLogContext &, const QString &message) {
    if (message.contains("QObject") || message.contains("QBasicTimer") ||
        message.contains("Timers cannot")) {
        const std::lock_guard lock(messageMutex);
        lifetimeWarnings.push_back(message);
    }
}
void require(bool condition, const QString &message) {
    if (!condition)
        throw std::runtime_error(message.toStdString());
}
void near(double actual, double expected, const char *message) {
    require(std::isfinite(actual) &&
                std::abs(actual - expected) <= 1e-10 * std::max(1.0, std::abs(expected)),
            message);
}
template <class T> T *control(QObject &owner, const char *name) {
    auto *result = dynamic_cast<T *>(owner.findChild<QObject *>(QString::fromLatin1(name)));
    require(result != nullptr, QString("Missing sweep control: %1").arg(name));
    return result;
}
class Events {
  public:
    Events() {
        QObject::connect(&guard_, &QTimer::timeout, [&] {
            for (auto *widget : QApplication::topLevelWidgets())
                if (auto *message = dynamic_cast<QMessageBox *>(widget);
                    message && message->isVisible()) {
                    errors_.push_back(message->windowTitle() + ": " + message->text());
                    message->done(QMessageBox::Ok);
                }
        });
        guard_.start(10);
    }
    template <class Predicate> void until(Predicate predicate, int timeoutMs, const char *message) {
        QElapsedTimer timer;
        timer.start();
        while (!predicate()) {
            QCoreApplication::processEvents(QEventLoop::AllEvents, 20);
            require(errors_.empty(), errors_.join("; "));
            require(timer.elapsed() < timeoutMs, message);
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
        require(errors_.empty(), errors_.join("; "));
    }
    void drain() {
        QElapsedTimer timer;
        timer.start();
        until([&] { return timer.elapsed() >= 100; }, 1000, "Qt event cleanup timed out.");
        QCoreApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
    }

  private:
    QTimer guard_;
    QStringList errors_;
};
Project smallBaseline() {
    Project project;
    project.parameters.entry.width = project.parameters.entry.height = .004;
    project.parameters.length = .020;
    project.parameters.wallThickness = .001;
    project.parameters.horizontalAngle = 20.0 * std::numbers::pi / 180.0;
    project.parameters.verticalAngle = 15.0 * std::numbers::pi / 180.0;
    project.meshing.targetEdgeLength = .008;
    project.meshing.adaptive.geometricTolerance = .00025;
    project.simulation.settings.minimumFrequency = 100;
    project.simulation.settings.maximumFrequency = 200;
    project.simulation.settings.frequencyCount = 1;
    return project;
}
void catalogAndConvergence(const QString &directory) {
    auto symmetric = smallBaseline();
    auto asymmetric = symmetric;
    asymmetric.parameters.asymmetricEnabled = true;
    const auto s = createSweepStudy(directory, symmetric, {{"horizontalAngle", .25, .35, 2}});
    const auto a = createSweepStudy(directory, asymmetric, {{"leftAngle", .25, .35, 4}});
    const auto c = createSweepStudy(directory, symmetric, {{"meshTargetEdgeLength", .008, .004, 3}},
                                    {100, 200});
    ParameterSweepDialog dialog(symmetric);
    const auto available = QGuiApplication::primaryScreen()->availableGeometry().size();
    require(dialog.width() <= available.width() - 32 &&
                dialog.height() <= available.height() - 80,
            "Initial sweep dialog size exceeds the available screen with frame reserve.");
    int gpuRequests = 0;
    dialog.onAcquireGpu = [&](auto) { ++gpuRequests; };
    dialog.show();
    dialog.resize(529, 320);
    QCoreApplication::processEvents();
    require(dialog.width() <= 529 && dialog.height() <= 320,
            "Sweep dialog cannot fit a compact screen.");
    auto *scroll = control<QScrollArea>(dialog, "sweepDialogScroll");
    require(scroll->verticalScrollBar()->maximum() > 0,
            "Compact sweep controls cannot be reached by scrolling.");
    auto *analysisAction = control<QPushButton>(dialog, "sweepAnalyzeVariant");
    scroll->ensureWidgetVisible(analysisAction);
    require(scroll->verticalScrollBar()->value() > 0,
            "The analysis action cannot be reached in the compact dialog.");
    auto *first = control<QComboBox>(dialog, "sweepParameter1");
    auto *second = control<QComboBox>(dialog, "sweepParameter2");
    require(first->findData("leftAngle") < 0 && second->findData("bottomAngle") < 0,
            "Symmetric catalogue contains independent left/bottom axes.");
    dialog.openStudy(sweepManifestPath(a));
    require(first->currentData().toString() == "leftAngle" && second->findData("bottomAngle") >= 0,
            "Asymmetric saved study did not rebuild its catalogue.");
    auto *table = control<QTableView>(dialog, "sweepVariants");
    table->selectionModel()->select(table->model()->index(3, 0),
                                    QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
    dialog.openStudy(sweepManifestPath(c));
    require(first->isHidden() && control<QCheckBox>(dialog, "sweepSecondAxis")->isHidden() &&
                !control<QLabel>(dialog, "convergenceMeshAxis")->isHidden() &&
                control<QSpinBox>(dialog, "sweepPoints1")->value() == 3 &&
                control<QSpinBox>(dialog, "sweepPoints1")->isEnabled(),
            "Convergence controls are not editable and exclusive.");
    dialog.openStudy(sweepManifestPath(s));
    require(first->currentData().toString() == "horizontalAngle" &&
                first->findData("leftAngle") < 0 && !first->isHidden() &&
                !control<QCheckBox>(dialog, "sweepSecondAxis")->isHidden(),
            "Switching back to symmetric study left stale controls.");
    require(gpuRequests == 0, "Opening or editing a study unexpectedly acquired the GPU.");
    PolarPlotWidget polar("Test", Qt::black);
    polar.setSeries({{"A", Qt::red, {{-1., -3.}, {0., 0.}, {1., -3.}}},
                     {"B", Qt::blue, {{-1., -6.}, {0., 0.}, {1., -6.}}}});
    require(polar.series().size() == 2 && polar.series()[1].color == QColor(Qt::blue) &&
                polar.series()[1].levels.front().levelDecibels == -6.,
            "Polar renderer lost a selected curve or its color.");
    polar.setLevels({{0., 0.}});
    require(polar.series().size() == 1 && polar.series().front().levels.size() == 1,
            "The existing single-curve polar contract was not retained.");
}
Project polygonCadBaseline(bool correctSecondComponent = false) {
    // Two separate synthetic octahedra exercise source-relative face masks.
    // The optional inward smaller component requires a nonempty valid correction.
    auto project = smallBaseline();
    auto cad = std::make_shared<CadSnapshot>();
    cad->documentId = "sweep-orientation-document";
    cad->sessionId = "sweep-orientation-session";
    cad->revision = 1;
    cad->objects.push_back({"solid",
                            "Solid",
                            "Orientation test",
                            "two-octahedra",
                            {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}});
    for (std::uint32_t component = 0; component < 2; ++component) {
        auto part = testing::sphere(0);
        const auto offset = static_cast<std::uint32_t>(cad->mesh.vertices.size());
        const double scale = component == 0 ? .06 : .03;
        for (auto point : part.vertices)
            cad->mesh.vertices.push_back(
                {point.x * scale + .03 * component, point.y * scale, point.z * scale});
        for (std::uint32_t i = 0; i < part.triangles.size(); ++i) {
            auto triangle = part.triangles[i];
            triangle.faceId = 1 + component * 8 + i;
            triangle.boundary = i == 0 ? BoundaryCondition::driven : BoundaryCondition::wall;
            for (auto &vertex : triangle.vertices)
                vertex += offset;
            if (component == 1 && correctSecondComponent) {
                std::swap(triangle.vertices[1], triangle.vertices[2]);
                project.cadReversedFaceIds.insert(triangle.faceId);
            }
            cad->mesh.triangles.push_back(triangle);
            CadFace face;
            face.id = triangle.faceId;
            face.objectId = "solid";
            face.sourceKey = QString("triangle-%1").arg(face.id);
            face.signature = {"mesh", face.sourceKey, "triangle", .00001, {}, {}, 3, 1};
            face.boundary = triangle.boundary;
            face.normalVelocity =
                component == 0 ? std::complex<double>{.4, -.7} : std::complex<double>{-.2, .3};
            face.state = CadFaceState::confirmed;
            cad->faces.push_back(std::move(face));
        }
    }
    (void)validateCadSnapshot(*cad);
    project.source = ProjectSource::freecad;
    project.cad = std::move(cad);
    return project;
}
void configureLength(ParameterSweepDialog &dialog, const QString &directory, int count = 2) {
    auto *parameter = control<QComboBox>(dialog, "sweepParameter1");
    const int length = parameter->findData("length");
    require(length >= 0, "The sweep has no length parameter.");
    parameter->setCurrentIndex(length);
    control<QDoubleSpinBox>(dialog, "sweepStart1")->setValue(20.0);
    control<QDoubleSpinBox>(dialog, "sweepEnd1")->setValue(21.0);
    control<QSpinBox>(dialog, "sweepPoints1")->setValue(count);
    control<QLineEdit>(dialog, "sweepOutputDirectory")->setText(directory);
}
QString manifest(const QString &parent) {
    QDirIterator files(parent, {"study.hsstudy"}, QDir::Files, QDirIterator::Subdirectories);
    require(files.hasNext(), "The sweep did not write its study manifest.");
    const auto result = files.next();
    require(!files.hasNext(), "A single sweep created more than one study.");
    return result;
}
SweepStudy run(ParameterSweepDialog &dialog, Events &events, const QString &directory,
               int count = 2, bool expectSuccess = true) {
    auto *start = control<QPushButton>(dialog, "sweepStart");
    auto *progress = control<QProgressBar>(dialog, "sweepProgress");
    auto *table = control<QTableView>(dialog, "sweepVariants");
    dialog.show();
    require(start->isEnabled(), "The prepared sweep cannot be started.");
    start->click();
    events.until(
        [&] {
            return start->isEnabled() && progress->maximum() == count && progress->value() == count;
        },
        180000, "The sweep did not finish its variants.");
    require(table->model()->rowCount() == count, "The GUI lost a sweep variant.");
    auto study = readSweepStudy(manifest(directory));
    require(study.variants.size() == static_cast<std::size_t>(count),
            "The manifest lost a sweep variant.");
    if (expectSuccess)
        for (const auto &variant : study.variants)
            require(variant.state == SweepState::completed,
                    "A genuine sweep variant failed: " + variant.error);
    return study;
}
void configureMeshSweep(ParameterSweepDialog &dialog, const QString &directory,
                        bool convergence = false) {
    if (convergence) {
        dialog.configureConvergence({100});
        control<QLineEdit>(dialog, "sweepOutputDirectory")->setText(directory);
        return;
    }
    auto *parameter = control<QComboBox>(dialog, "sweepParameter1");
    parameter->setCurrentIndex(parameter->findData("meshTargetEdgeLength"));
    control<QDoubleSpinBox>(dialog, "sweepStart1")->setValue(8);
    control<QDoubleSpinBox>(dialog, "sweepEnd1")->setValue(8);
    control<QSpinBox>(dialog, "sweepPoints1")->setValue(1);
    control<QLineEdit>(dialog, "sweepOutputDirectory")->setText(directory);
}
void orientationConflicts(Events &events, const QString &directory) {
    for (const bool cadMode : {false, true})
        for (const bool complete : {false, true}) {
            auto baseline = cadMode ? polygonCadBaseline() : smallBaseline();
            auto &mask = cadMode ? baseline.cadReversedFaceIds : baseline.draftReversedFaceIds;
            mask.insert(1);
            if (complete) {
                const auto mesh =
                    cadMode ? baseline.cad->mesh : DraftGeometry(baseline.parameters).previewMesh();
                for (const auto &triangle : mesh.triangles)
                    mask.insert(triangle.faceId);
            }
            const auto raw = cadMode ? baseline.cad->mesh : SurfaceMesh{};
            const auto path = QDir(directory).filePath(
                QString("%1-%2").arg(cadMode ? "cad" : "draft", complete ? "inward" : "partial"));
            require(QDir().mkpath(path), "Cannot create the orientation test directory.");
            ParameterSweepDialog dialog(baseline);
            configureMeshSweep(dialog, path, cadMode);
            const auto study = run(dialog, events, path, cadMode ? 2 : 1, false);
            for (const auto &variant : study.variants) {
                require(variant.state == SweepState::failed &&
                            variant.error.startsWith("Flächennormalen: "),
                        "The sweep ignored an invalid orientation mask: " + variant.error);
                require(!QFileInfo::exists(QDir(study.directory).filePath(variant.projectFile)),
                        "An invalid orientation produced a computed project archive.");
            }
            require((cadMode ? study.baseline.cadReversedFaceIds
                             : study.baseline.draftReversedFaceIds) == mask,
                    "The study archive lost the rejected orientation edit.");
            if (cadMode)
                require(baseline.cad->mesh == raw && study.baseline.cad->mesh == raw,
                        "The orientation edit mutated the immutable CAD source.");
        }
}
void correctedCadArchive(Events &events, const QString &directory, const std::string &device) {
    const auto baseline = polygonCadBaseline(true);
    const auto raw = baseline.cad->mesh;
    ParameterSweepDialog dialog(baseline);
    configureMeshSweep(dialog, directory, true);
    const auto study = run(dialog, events, directory);
    for (const auto &variant : study.variants) {
        const auto archived = readProject(QDir(study.directory).filePath(variant.projectFile));
        const auto expected = applyFaceOrientation(raw, baseline.cadReversedFaceIds);
        require(archived.cad && archived.cad->mesh == raw && archived.preview &&
                    *archived.preview == expected &&
                    archived.cadReversedFaceIds == baseline.cadReversedFaceIds,
                "The corrected CAD preview was not derived from the unchanged archived source.");
        require(archived.mesh && archived.mesh->quality.watertight() &&
                    archived.simulation.current && archived.simulation.request &&
                    archived.simulation.result && archived.simulation.result->device == device &&
                    archived.simulation.result->frequencies.size() == 1,
                "The valid CAD orientation did not produce a current real CUDA archive.");
        const auto &velocities = archived.simulation.request->excitations.front().faces;
        for (const auto &face : baseline.cad->faces)
            if (face.boundary == BoundaryCondition::driven) {
                const auto found = std::find_if(velocities.begin(), velocities.end(),
                                                [&](const auto &v) { return v.faceId == face.id; });
                require(found != velocities.end() && found->normalVelocity == face.normalVelocity,
                        "Orientation silently changed the complex driven velocity.");
            }
        require(baseline.cad->mesh == raw, "A computed sweep mutated its CAD input snapshot.");
    }
}
void verifiedCadProvider(Events &events, const QString &directory) {
    for (const bool stale : {false, true}) {
        auto baseline = smallBaseline();
        auto cad = std::make_shared<CadSnapshot>(testing::exactCadFixture());
        baseline.source = ProjectSource::freecad;
        cad->faces.front().boundary = BoundaryCondition::driven;
        for (auto &triangle : cad->mesh.triangles)
            if (triangle.faceId == cad->faces.front().id)
                triangle.boundary = BoundaryCondition::driven;
        for (const auto &face : cad->faces)
            baseline.cadReversedFaceIds.insert(face.id);
        cad->geometry.reset();
        if (stale) {
            std::vector<CadBrepPart> parts;
            for (const auto &part : cad->brepObjects)
                parts.push_back({part.data.toStdString(), part.faceIds});
            cad->geometry = CadGeometry::fromBrepParts(parts, .002);
        }
        baseline.cad = cad;
        const auto raw = cad->mesh;
        const auto originalProvider = cad->geometry;
        const auto path = QDir(directory).filePath(stale ? "stale-provider" : "missing-provider");
        require(QDir().mkpath(path), "Cannot create the CAD provider test directory.");
        ParameterSweepDialog dialog(baseline);
        configureMeshSweep(dialog, path, true);
        const auto study = run(dialog, events, path, 2, false);
        for (const auto &variant : study.variants) {
            require(variant.state == SweepState::failed &&
                        variant.error.startsWith("Flächennormalen: "),
                    "The exact sweep did not use its verified BRep provider before the orientation "
                    "gate: " +
                        variant.error);
            require(!QFileInfo::exists(QDir(study.directory).filePath(variant.projectFile)),
                    "The inward exact CAD sweep wrote a computed archive.");
        }
        require(cad->geometry == originalProvider && cad->mesh == raw,
                "Preparing a sweep provider mutated the caller's CAD snapshot.");
    }
}
Project checkedVariant(const SweepStudy &study, std::size_t index, const std::string &device) {
    auto project =
        readProject(QDir(study.directory).filePath(study.variants.at(index).projectFile));
    const auto &archive = project.simulation;
    require(project.preview && project.mesh && archive.mesh && archive.request && archive.result &&
                archive.current && archive.revision == project.revision &&
                project.geometryRevision == project.revision &&
                project.meshRevision == project.revision,
            "A saved sweep variant has no complete current geometry/mesh/solver snapshot.");
    const auto &result = *archive.result;
    require(result.device == device && result.peakDeviceBytes > 0 &&
                result.triangleCount == archive.mesh->triangles.size() &&
                result.pressureDofCount == archive.mesh->vertices.size() &&
                result.frequencies.size() == 1,
            "The saved result does not describe a real CUDA Galerkin solve.");
    require(archive.request->frequencies == std::vector<double>{100.0},
            "The single-frequency settings did not produce the actual 100-Hz request.");
    const auto &frequency = result.frequencies.front();
    near(frequency.frequency, 100.0, "The stored result frequency is wrong.");
    require(frequency.assemblyMilliseconds > 0.0 && frequency.factorizationMilliseconds > 0.0 &&
                frequency.excitations.size() == 1,
            "The archive contains no measured CUDA assembly/factorization.");
    const auto &excitation = frequency.excitations.front();
    require(excitation.boundary()->surfacePressure.size() == result.triangleCount &&
                excitation.observerPressure.size() == archive.request->observers.size() &&
                !excitation.observerPressure.empty() && !excitation.impedances.empty(),
            "The sweep result is missing pressure or impedance data.");
    require(std::isfinite(std::abs(excitation.observerPressure.front())) &&
                std::abs(excitation.observerPressure.front()) > 0.0 &&
                excitation.boundaryResidual <= archive.request->numerics.relativeResidualTolerance,
            "The actual CUDA result is zero, non-finite or outside the residual gate.");
    near(archive.request->observers.front().z, 1.0 + project.parameters.length,
         "Draft sweep observers are not referenced to the same mouth origin as the app.");
    return project;
}
void comparePlots(ParameterSweepDialog &dialog, const std::vector<Project> &projects) {
    auto *table = control<QTableView>(dialog, "sweepVariants");
    QItemSelection selection(table->model()->index(0, 0), table->model()->index(1, 3));
    table->selectionModel()->select(selection, QItemSelectionModel::ClearAndSelect |
                                                   QItemSelectionModel::Rows);
    require(table->selectionModel()->selectedRows().size() == 2,
            "The variant table did not retain its multiple selection.");
    const auto *onAxis = control<PlotWidget>(dialog, "sweepOnAxis");
    const auto *impedance = control<PlotWidget>(dialog, "sweepImpedance");
    require(onAxis->series().size() == 2 && impedance->series().size() == 4,
            "Multiple variant selection did not populate the comparison curves.");
    const auto *polarH = control<PolarPlotWidget>(dialog, "sweepPolarH");
    const auto *polarV = control<PolarPlotWidget>(dialog, "sweepPolarV");
    require(polarH->series().size() == 2 && polarV->series().size() == 2,
            "Selected saved CUDA variants did not populate both polar cuts.");
    auto *mapChoice = control<QComboBox>(dialog, "sweepMapVariant");
    require(mapChoice->count() == 3 && mapChoice->currentIndex() == 0,
            "Multiple selected variants must require an explicit map selection.");
    mapChoice->setCurrentIndex(mapChoice->findData(QVariant::fromValue(qulonglong(1))));
    const auto *mapH = control<HeatmapWidget>(dialog, "sweepMapH");
    const auto *mapV = control<HeatmapWidget>(dialog, "sweepMapV");
    const auto &mapPressure =
        projects[1].simulation.result->frequencies.front().excitations.front().observerPressure;
    require(!mapH->anglesDegrees().empty() && mapH->frequencies() == std::vector<double>{100.0} &&
                mapH->anglesDegrees() == mapV->anglesDegrees(),
            "Saved frequency and H/V angle grids were not preserved in polar maps.");
    near(mapH->levelsDecibels().front(),
         20. * std::log10(std::abs(mapPressure[1]) / std::abs(mapPressure.front())),
         "Polar map does not show the selected variant's saved complex pressure.");
    for (std::size_t i = 0; i < projects.size(); ++i) {
        const auto &project = projects[i];
        const auto &excitation = project.simulation.result->frequencies.front().excitations.front();
        require(!polarH->series()[i].levels.empty() &&
                    polarH->series()[i].levels.size() == polarV->series()[i].levels.size(),
                "A saved polar cut is missing.");
        const auto expected = 20. * std::log10(std::abs(excitation.observerPressure[1]) /
                                                std::abs(excitation.observerPressure.front()));
        near(polarH->series()[i].levels.front().levelDecibels, expected,
             "Polar data differs from the saved complex pressure.");
        const auto &on = onAxis->series()[i];
        require(on.points.size() == 1, "An on-axis comparison curve has no actual data point.");
        near(on.points.front().x(), 100.0,
             "The plotted frequency differs from the solved frequency.");
        near(on.points.front().y(), soundPressureLevelDecibels(excitation.observerPressure.front()),
             "The on-axis comparison does not display the saved CUDA pressure.");
        const auto &z = excitation.impedances.front();
        const auto medium = project.simulation.request->medium;
        const auto normalized =
            normalizedVolumeImpedance(z.volumeImpedance, z.area, medium.density, medium.soundSpeed);
        require(impedance->series()[2 * i].points.size() == 1 &&
                    impedance->series()[2 * i + 1].points.size() == 1,
                "An impedance comparison curve has no actual data point.");
        near(impedance->series()[2 * i].points.front().y(), normalized.real(),
             "The comparison real impedance uses the wrong original normalization.");
        near(impedance->series()[2 * i + 1].points.front().y(), normalized.imag(),
             "The comparison imaginary impedance uses the wrong original normalization.");
    }
    // At 100 Hz this small radiator may have no -6 dB crossing. Keep those
    // physical gaps while still checking that both variants bind to each plot.
    for (const auto *name : {"sweepBeamwidthH", "sweepBeamwidthV"}) {
        const auto *plot = control<PlotWidget>(dialog, name);
        require(plot->series().size() == 2 && plot->series()[0].points.size() == 1 &&
                    plot->series()[1].points.size() == 1,
                "The beamwidth comparison lost a selected variant.");
    }
}
void geometrySweepAndReopen(const Project &baseline, Events &events, const QString &directory,
                            const std::string &device) {
    ParameterSweepDialog dialog(baseline);
    configureLength(dialog, directory);
    const auto study = run(dialog, events, directory);
    require(study.axes.size() == 1 && study.axes.front().parameter == "length",
            "The GUI started the wrong parameter sweep.");
    std::vector<Project> projects;
    for (std::size_t i = 0; i < 2; ++i) {
        projects.push_back(checkedVariant(study, i, device));
        near(projects.back().parameters.length, .020 + .001 * i,
             "The saved variant contains the wrong length input.");
        near(projects.back().parameters.entry.width, .004, "The sweep altered the baseline entry.");
    }
    comparePlots(dialog, projects);
    ParameterSweepDialog reopened(baseline);
    std::optional<Project> chosen;
    reopened.onVariantChosen = [&](Project project) { chosen = std::move(project); };
    reopened.openStudy(sweepManifestPath(study));
    reopened.show();
    comparePlots(reopened, projects);
    auto *table = control<QTableView>(reopened, "sweepVariants");
    table->selectionModel()->select(table->model()->index(1, 0),
                                    QItemSelectionModel::ClearAndSelect |
                                        QItemSelectionModel::Rows);
    auto *adopt = control<QPushButton>(reopened, "sweepAdopt");
    require(adopt->isEnabled(), "A reopened computed variant cannot be adopted.");
    adopt->click();
    require(chosen && chosen->simulation.current,
            "The variant callback lost its computed snapshot.");
    near(chosen->parameters.length, .021, "The callback adopted the wrong selected variant.");
}
void convergence(const Project &baseline, Events &events, const QString &directory,
                 const std::string &device) {
    ParameterSweepDialog dialog(baseline);
    dialog.configureConvergence({100.0});
    control<QSpinBox>(dialog, "sweepPoints1")->setValue(3);
    control<QLineEdit>(dialog, "sweepOutputDirectory")->setText(directory);
    require(control<QComboBox>(dialog, "sweepParameter1")->currentData().toString() ==
                "meshTargetEdgeLength",
            "Convergence did not switch to the shared mesh parameter.");
    const auto study = run(dialog, events, directory, 3);
    require(study.frequencyOverride == std::vector<double>{100.0},
            "Convergence did not persist its selected frequency.");
    std::vector<Project> projects;
    for (std::size_t i = 0; i < 3; ++i) {
        projects.push_back(checkedVariant(study, i, device));
        near(projects.back().meshing.targetEdgeLength, .008 - .002 * i,
             "Convergence did not solve the actual three linear mesh inputs.");
        near(projects.back().parameters.length, baseline.parameters.length,
             "Mesh convergence altered the model geometry.");
    }
    projects.resize(2);
    comparePlots(dialog, projects);
}
void cancelOnClose(const Project &baseline, Events &events, const QString &directory) {
    auto dialog = std::make_unique<ParameterSweepDialog>(baseline);
    configureLength(*dialog, directory, 16);
    dialog->show();
    auto *table = control<QTableView>(*dialog, "sweepVariants");
    control<QPushButton>(*dialog, "sweepStart")->click();
    events.until(
        [&] {
            for (int row = 0; row < table->model()->rowCount(); ++row)
                if (table->model()->data(table->model()->index(row, 2)).toString() == "In Arbeit")
                    return true;
            return false;
        },
        10000, "The cancellation test never started its real background worker.");
    const QPointer<ParameterSweepDialog> lifetime(dialog.get());
    QElapsedTimer closeTime;
    closeTime.start();
    dialog->close();
    events.until([&] { return !dialog->isVisible(); }, 15000,
                 "Closing the study did not finish its pending worker cancellation.");
    dialog.reset();
    require(lifetime.isNull() && closeTime.elapsed() < 15000,
            "Closing a busy study did not promptly stop and join its worker.");
    events.drain();
    const auto study = readSweepStudy(manifest(directory));
    require(
        std::none_of(study.variants.begin(), study.variants.end(),
                     [](const auto &variant) { return variant.state == SweepState::running; }) &&
            std::any_of(study.variants.begin(), study.variants.end(),
                        [](const auto &variant) { return variant.state == SweepState::cancelled; }),
        "The closed study left a running job or did not persist cancellation.");
}
} // namespace

int main(int argc, char **argv) {
    qputenv("QT_QPA_PLATFORM", "offscreen");
    QApplication application(argc, argv);
    application.setQuitOnLastWindowClosed(false);
    const auto previousHandler = qInstallMessageHandler(qtMessages);
    try {
        const bool orientationOnly = application.arguments().contains("--orientation-only");
        QTemporaryDir output(QDir::tempPath() + "/HornSculpt-sweep-XXXXXX");
        require(output.isValid(), "Cannot create the temporary study directory.");
        for (const auto *name : {"grid", "convergence", "cancel", "orientation", "corrected-cad"})
            require(QDir(output.path()).mkdir(name), "Cannot create a study output directory.");
        Events events;
        const bool uiOnly = application.arguments().contains("--sweep-ui-only");
        if (uiOnly) {
            catalogAndConvergence(output.path());
            qInstallMessageHandler(previousHandler);
            std::cout << "Sweep catalogue and convergence UI passed without CUDA.\n";
            return 0;
        }
        orientationConflicts(events, output.filePath("orientation"));
        verifiedCadProvider(events, output.filePath("orientation"));
        std::string device;
        if (!orientationOnly) {
            device = CudaBemSolver::deviceDescription();
            require(!device.empty(),
                    "No genuine CUDA device is available for the integration test.");
            const auto baseline = smallBaseline();
            geometrySweepAndReopen(baseline, events, output.filePath("grid"), device);
            convergence(baseline, events, output.filePath("convergence"), device);
            correctedCadArchive(events, output.filePath("corrected-cad"), device);
            cancelOnClose(baseline, events, output.filePath("cancel"));
        }
        events.drain();
        {
            const std::lock_guard lock(messageMutex);
            require(lifetimeWarnings.empty(),
                    "Qt worker/lifetime warning: " + lifetimeWarnings.join("; "));
        }
        qInstallMessageHandler(previousHandler);
        if (orientationOnly)
            std::cout << "Sweep orientation CPU gates passed: inconsistent and inward Draft/CAD "
                         "masks rejected before BEM; input masks and sources preserved.\n";
        else
            std::cout << "Parameter sweep GUI passed on " << device
                      << ": genuine CUDA variants, archived plots/reopen/adopt, h/h2 convergence, "
                         "corrected CAD orientation/archive/complex velocity, close/cancel.\n";
        return 0;
    } catch (const std::exception &error) {
        qInstallMessageHandler(previousHandler);
        std::cerr << error.what() << '\n';
        return 1;
    }
}
