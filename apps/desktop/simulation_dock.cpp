#include "simulation_dock.hpp"
#include "acoustic_field_widget.hpp"
#include "acoustic_frame_dialog.hpp"
#include "angular_response_widget.hpp"
#include "baffle_dialog.hpp"
#include "frequency_schedule.hpp"
#include "galerkin/exterior_domain.hpp"
#include "live_trace.hpp"
#include "latest_publication.hpp"
#include "numeric_input.hpp"
#include "phase_response_widget.hpp"
#include "planning_frequencies.hpp"
#include "result_export.hpp"
#include "simulation_request.hpp"
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QLabel>
#include <QLineEdit>
#include <QMainWindow>
#include <QMenu>
#include <QMessageBox>
#include <QProgressBar>
#include <QPushButton>
#include <QSaveFile>
#include <QScrollArea>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QTabWidget>
#include <QTabBar>
#include <QToolButton>
#include <QStyle>
#include <QTimer>
#include <QVBoxLayout>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <hornsculpt/acoustic_analysis.hpp>
#include <hornsculpt/axisymmetric_geometry.hpp>
#include <hornsculpt/axisymmetric_mesh.hpp>
#include <hornsculpt/heatmap_widget.hpp>
#include <hornsculpt/mesh_orientation.hpp>
#include <hornsculpt/plot_axes.hpp>
#include <hornsculpt/plot_scale.hpp>
#include <hornsculpt/plot_widget.hpp>
#include <hornsculpt/polar_plot_widget.hpp>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>
#include <utility>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace hornsculpt {
namespace {
constexpr std::size_t angularCount = 73;
constexpr double pressureReference = 20e-6;
std::optional<std::uint64_t> availablePhysicalMemory() {
#ifdef _WIN32
    MEMORYSTATUSEX status{};
    status.dwLength = sizeof(status);
    if (GlobalMemoryStatusEx(&status))
        return status.ullAvailPhys;
#elif defined(_SC_AVPHYS_PAGES)
    const long pages = sysconf(_SC_AVPHYS_PAGES), pageSize = sysconf(_SC_PAGESIZE);
    if (pages > 0 && pageSize > 0)
        return static_cast<std::uint64_t>(pages) * static_cast<std::uint64_t>(pageSize);
#endif
    return {};
}
void requireResultMemory(std::size_t panels, std::size_t vertices, std::size_t exteriorPanels,
                         std::size_t faces, std::size_t frequencies, std::size_t excitations,
                         bool realtime) {
    // The common estimator counts shared boundary arrays once and independent
    // snapshot data by simultaneous owners. Plot points are GUI-only overhead.
    const std::size_t snapshots = realtime ? 4 : 2;
    const long double estimate = estimateResultBytes(
        {panels, vertices, exteriorPanels, 1 + 2 * angularCount, faces, frequencies,
         excitations, snapshots}) +
        static_cast<long double>(frequencies) * snapshots *
            (3 + 2 * angularCount) * sizeof(QPointF);
    const auto available = availablePhysicalMemory();
    if (frequencies > std::vector<double>{}.max_size() ||
        estimate > static_cast<long double>(std::numeric_limits<std::size_t>::max()) ||
        (available && estimate > static_cast<long double>(*available))) {
        const auto message =
            QStringLiteral("Die CPU-Ergebnisdaten benötigen geschätzt %1 MiB%2. "
                           "Frequenzraster verkleinern oder mehr Arbeitsspeicher bereitstellen.")
                .arg(static_cast<double>(estimate / (1024 * 1024)), 0, 'f', 1)
                .arg(available
                         ? QStringLiteral(", verfügbar sind %1 MiB")
                               .arg(static_cast<double>(*available) / (1024 * 1024), 0, 'f', 1)
                         : QString{});
        throw std::invalid_argument(message.toStdString());
    }
}
QDoubleSpinBox *number(double value, double minimum, double maximum, int decimals,
                       const QString &suffix, const QString &explanation, QWidget *parent) {
    auto *spin = new DecimalInput(parent);
    spin->setRange(minimum, maximum);
    spin->setDecimals(decimals);
    spin->setValue(value);
    spin->setSuffix(suffix);
    spin->setToolTip(explanation);
    spin->setWhatsThis(explanation);
    return spin;
}
Vec3 scaled(Vec3 a, double scale) {
    return {a.x * scale, a.y * scale, a.z * scale};
}
Vec3 added(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 normalized(Vec3 a) {
    const double norm = std::hypot(a.x, a.y, a.z);
    if (!std::isfinite(norm) || norm <= 0)
        throw std::invalid_argument("Akustische Achse ist ungültig.");
    return scaled(a, 1 / norm);
}
double decibels(double magnitude, double reference) {
    // The display floor prevents logarithmic infinities; complex result data stays unchanged.
    return std::max(-300.0, soundPressureLevelDecibels({magnitude, 0}, reference));
}
std::vector<double> anglesDegrees(std::size_t count = angularCount) {
    std::vector<double> angles;
    angles.reserve(count);
    const double first = count == 73 ? -180 : -90;
    for (std::size_t i = 0; i < count; ++i)
        angles.push_back(first + 5 * static_cast<double>(i));
    return angles;
}
std::vector<AngularLevel> displayedCut(const std::vector<double> &levelsDb, double onAxisDb,
                                       bool atAngle, double referenceAngle) {
    double reference = onAxisDb;
    const auto angles = anglesDegrees(levelsDb.size());
    if (atAngle && referenceAngle != 0) {
        std::vector<AcousticsAngularSample> samples;
        for (std::size_t i = 0; i < levelsDb.size(); ++i)
            samples.push_back({angles[i], levelsDb[i]});
        const auto level = angularLevelDecibels(samples, referenceAngle);
        if (!level)
            return {};
        reference = *level;
    } else if (!atAngle) {
        reference = -std::numeric_limits<double>::infinity();
        for (const auto value : levelsDb)
            if (std::isfinite(value))
                reference = std::max(reference, value);
    }
    if (!std::isfinite(reference) || reference <= -300)
        return {};
    std::vector<AngularLevel> levels;
    for (std::size_t i = 0; i < levelsDb.size(); ++i)
        levels.push_back({angles[i] * std::numbers::pi / 180, levelsDb[i] - reference});
    return levels;
}
bool sameDisplayDirections(const BemRequest &a, const BemRequest &b) {
    if (a.infiniteBaffle != b.infiniteBaffle)
        return false;
    if (a.observers.size() != b.observers.size() ||
        a.farField.has_value() != b.farField.has_value())
        return false;
    if (!a.farField)
        return a.observers == b.observers;
    if (a.farField->referenceDistance != b.farField->referenceDistance)
        return false;
    for (std::size_t i = 0; i < a.observers.size(); ++i) {
        const auto u = normalized(added(a.observers[i], scaled(a.farField->origin, -1)));
        const auto v = normalized(added(b.observers[i], scaled(b.farField->origin, -1)));
        if (std::hypot(u.x - v.x, u.y - v.y, u.z - v.z) > 1e-12)
            return false;
    }
    return true;
}

bool sameRequest(const BemRequest &a, const BemRequest &b) {
    if (a.infiniteBaffle != b.infiniteBaffle || a.axisymmetric != b.axisymmetric)
        return false;
    // An absent raster is the complete solve raster in earlier archives.
    const auto &aPlan =
        a.interactionPlanFrequencies.empty() ? a.frequencies : a.interactionPlanFrequencies;
    const auto &bPlan =
        b.interactionPlanFrequencies.empty() ? b.frequencies : b.interactionPlanFrequencies;
    if (aPlan != bPlan)
        return false;
    if (a.coupling.has_value() != b.coupling.has_value() ||
        (a.coupling && !lem::samePhysics(*a.coupling, *b.coupling)))
        return false;
    if (a.farField.has_value() != b.farField.has_value() ||
        (a.farField && (a.farField->origin != b.farField->origin ||
                        a.farField->referenceDistance != b.farField->referenceDistance)))
        return false;
    if (a.frequencies != b.frequencies || a.medium.density != b.medium.density ||
        a.medium.soundSpeed != b.medium.soundSpeed || a.observers.size() != b.observers.size() ||
        a.excitations.size() != b.excitations.size())
        return false;
    const auto &x = a.numerics;
    const auto &y = b.numerics;
    if (x.quadratureOrder != y.quadratureOrder || x.singularOrder != y.singularOrder ||
        x.nearDistanceRatio != y.nearDistanceRatio ||
        x.nearRelativeTolerance != y.nearRelativeTolerance ||
        x.nearMaximumDepth != y.nearMaximumDepth ||
        x.relativeResidualTolerance != y.relativeResidualTolerance ||
        x.maximumRefinementSteps != y.maximumRefinementSteps || x.acceleration != y.acceleration ||
        x.fastMultipole != y.fastMultipole)
        return false;
    for (std::size_t i = 0; i < a.observers.size(); ++i) {
        const auto u = a.observers[i], v = b.observers[i];
        if (std::hypot(u.x - v.x, u.y - v.y, u.z - v.z) > 1e-12)
            return false;
    }
    for (std::size_t i = 0; i < a.excitations.size(); ++i) {
        auto first = a.excitations[i].faces, second = b.excitations[i].faces;
        const auto byFace = [](const FaceVelocity &u, const FaceVelocity &v) {
            return u.faceId < v.faceId;
        };
        std::sort(first.begin(), first.end(), byFace);
        std::sort(second.begin(), second.end(), byFace);
        if (first.size() != second.size())
            return false;
        for (std::size_t face = 0; face < first.size(); ++face)
            if (first[face].faceId != second[face].faceId ||
                first[face].normalVelocity != second[face].normalVelocity)
                return false;
    }
    return true;
}
} // namespace

struct SimulationDock::GeometrySource {
    std::shared_ptr<const DraftGeometry> draft;
    std::shared_ptr<const CadSnapshot> cad;
    std::uint64_t revision{};
    std::map<std::uint32_t, double> faceTargets;
    std::set<std::uint32_t> reversedFaces;
};
struct SimulationDock::SolveJob {
    std::shared_ptr<const SurfaceMesh> mesh;
    std::vector<FaceVelocity> boundaryDefinitions;
    BemRequest request;
    std::shared_ptr<const GeometrySource> source;
    std::optional<AxisymmetricSamplingSettings> axisymmetricSampling;
    std::uint64_t revision{};
    std::uint64_t generation{};
    std::uint64_t id{};
    bool realtime{};
    std::stop_token stop;
    struct Partial {
        std::shared_ptr<const BemResult> result;
        BemRequest request;
        double frequency{};
    };
    mutable LatestPublication<Partial> partial;
    mutable LatestPublication<BemProgress> progress;
};

void SimulationDock::requestFieldEvaluation() {
    acousticField_->requestEvaluation();
}
bool SimulationDock::fieldEvaluationInProgress() const noexcept {
    return acousticField_->evaluationInProgress();
}
bool SimulationDock::fieldReady() const noexcept {
    return acousticField_->fieldReady();
}
QString SimulationDock::lastFieldError() const {
    return acousticField_->lastError();
}

SimulationDock::SimulationDock(QMainWindow *parent) : QDockWidget(tr("Akustik"), parent) {
    setObjectName("simulationDock");
    auto *root = new QWidget;
    auto *rootLayout = new QVBoxLayout(root);
    rootLayout->setContentsMargins(7, 6, 7, 6);
    rootLayout->setSpacing(5);
    solverCommands_ = new QWidget;
    solverCommands_->setObjectName("solverCommands");
    auto *commands = new QHBoxLayout(solverCommands_);
    commands->setContentsMargins(0, 0, 0, 0);
    workspaceTabs_ = new QTabWidget;
    workspaceTabs_->setObjectName("simulationWorkspace");
    workspaceTabs_->setMinimumSize(300, 100);
    auto *content = new QWidget;
    auto *layout = new QVBoxLayout(content);
    layout->setContentsMargins(7, 6, 7, 6);
    layout->setSpacing(5);
    auto *controls = new QHBoxLayout;
    minimumFrequency_ =
        number(SimulationSettings{}.minimumFrequency, 1, 100000, 0, " Hz",
               tr("Untere Frequenz des logarithmischen Frequenzrasters. Einheit: Hz."), content);
    maximumFrequency_ = number(SimulationSettings{}.maximumFrequency, 1, 100000, 0, " Hz",
                               tr("Obere Frequenz. Die benötigte Netzauflösung hängt von der "
                                  "Wellenlänge ab. Einheit: Hz."),
                               content);
    frequencyCount_ = new IntegerInput;
    minimumFrequency_->setObjectName("simulationMinimumFrequency");
    maximumFrequency_->setObjectName("simulationMaximumFrequency");
    frequencyCount_->setObjectName("simulationFrequencyCount");
    minimumFrequency_->setAccessibleName(tr("Untere Frequenz"));
    maximumFrequency_->setAccessibleName(tr("Obere Frequenz"));
    frequencyCount_->setAccessibleName(tr("Frequenzpunkte"));
    frequencyCount_->setRange(1, std::numeric_limits<int>::max());
    frequencyCount_->setValue(17);
    frequencyCount_->setToolTip(tr("Anzahl logarithmisch verteilter Frequenzen. Bei einer Frequenz "
                                   "wird der untere Wert gerechnet."));
    controls->addWidget(new QLabel(tr("Frequenz")));
    controls->addWidget(minimumFrequency_);
    controls->addWidget(new QLabel("–"));
    controls->addWidget(maximumFrequency_);
    controls->addWidget(new QLabel(tr("Punkte")));
    controls->addWidget(frequencyCount_);
    solve_ = new QPushButton(tr("Berechnen"));
    solve_->setObjectName("startSimulation");
    solve_->setProperty("primaryAction", true);
    solve_->setToolTip(tr("Berechnet den eingestellten Frequenzbereich mit dem aktuellen Netz "
                          "und den zugeordneten Antrieben."));
    stop_ = new QPushButton(tr("Abbrechen"));
    stop_->setObjectName("stopSimulation");
    realtime_ = new QPushButton(tr("Echtzeit"));
    realtime_->setObjectName(QStringLiteral("realtimeSolve"));
    realtime_->setCheckable(true);
    realtime_->setToolTip(
        tr("Nach Änderungen automatisch neu berechnen: erst beide Frequenzgrenzen, "
           "dann große Lücken, hohe Bereiche bevorzugt. Neue Stützpunkte treiben "
           "die bisherige Anzeige; echte Zwischenwerte ersetzen die Schätzung. "
           "Nach dem Raster auf Änderungen warten."));
    solve_->setEnabled(false);
    stop_->setEnabled(false);
    commands->addWidget(solve_);
    commands->addWidget(stop_);
    commands->addWidget(realtime_);
    commands->addStretch();
    rootLayout->addWidget(solverCommands_);
    layout->addLayout(controls);
    auto *solverControls = new QHBoxLayout;
    fastMultipole_ = new QCheckBox(tr("Fast Multipole"));
    fastMultipole_->setObjectName("fastMultipoleEnabled");
    fastMultipole_->setToolTip(
        tr("Bündelt entfernte Wechselwirkungen durch eine hierarchische Approximation "
           "und löst Burton–Miller iterativ. Nah- und Singularintegrale bleiben direkt. "
           "Der ursprüngliche Gleichungsrest wird unabhängig geprüft. "
           "Ausgeschaltet: dichte Matrix mit direkter Lösung. Bei kleinen Netzen kann diese "
           "schneller sein."));
    multipoleOrder_ = new IntegerInput;
    multipoleOrder_->setObjectName("fastMultipoleOrder");
    multipoleOrder_->setRange(3, 8);
    multipoleOrder_->setValue(static_cast<int>(fastMultipoleSettings_.interpolationOrder));
    multipoleOrder_->setEnabled(false);
    multipoleOrder_->setToolTip(
        tr("Anzahl der Interpolationspunkte je Raumrichtung. Höhere Ordnung verbessert die "
           "Fernwechselwirkungen und benötigt mehr Rechenzeit und Speicher. Dimensionslos; "
           "Konvergenz durch Vergleich mehrerer Ordnungen prüfen."));
    solverControls->addWidget(fastMultipole_);
    solverControls->addWidget(new QLabel(tr("Ordnung")));
    solverControls->addWidget(multipoleOrder_);
    baffle_ = new QPushButton(tr("Schallwand …"));
    baffle_->setObjectName("infiniteBaffle");
    baffle_->setToolTip(tr("Freifeld oder unendliche Schallwand mit Luftkavität und Apertur. "
                           "Ebene und Flächen-IDs in Modellkoordinaten einstellen."));
    baffle_->setWhatsThis(baffle_->toolTip());
    solverControls->addWidget(baffle_);
    connect(baffle_, &QPushButton::clicked, this, [this] {
        BaffleDialog dialog(infiniteBaffle_, this);
        if (dialog.exec() == QDialog::Accepted)
            setInfiniteBaffle(dialog.infiniteBaffle());
    });
    solverControls->addStretch();
    layout->addLayout(solverControls);
    auto *geometryControls = new QHBoxLayout;
    geometryMode_ = new QComboBox;
    geometryMode_->setObjectName("bemGeometryMode");
    geometryMode_->addItems({tr("3D"), tr("2D-Rotation (m=0)")});
    geometryMode_->setToolTip(tr(
        "3D: vollständiger Burton–Miller-Randraum. 2D-Rotation: echte m=0-Galerkin-Reduktion "
        "auf dem Meridian einer geprüften Rotationsgeometrie; nur rotationssymmetrische "
        "Anregung. Benötigt eine Luftkavität mit ebener Apertur in einer unendlichen "
        "Schallwand. Die Quellachse kommt aus CAD bzw. Draft +Z, unabhängig vom Beobachterrahmen. "
        "Kein 3D-Qualitätsnetz erforderlich."));
    geometryMode_->setWhatsThis(geometryMode_->toolTip());
    geometryControls->addWidget(new QLabel(tr("Geometrie")));
    geometryControls->addWidget(geometryMode_);
    meridianLength_ =
        number(2, .000001, 1000000, 6, " mm",
               tr("Maximale Segmentlänge auf dem Meridian [mm], intern m. Diskretisiert die echten "
                  "meridionalen Unbekannten, nicht ein isotropes 3D-Netz."),
               content);
    meridianLength_->setObjectName("axisymmetricMaximumLength");
    meridianTolerance_ = number(.05, .000001, 1000000, 6, " mm",
                                tr("Maximale geometrische Sehnenabweichung des Meridians [mm], "
                                   "intern m. Die Quelle bleibt die exakte Geometrie."),
                                content);
    meridianTolerance_->setObjectName("axisymmetricChordTolerance");
    rotationalTolerance_ = number(
        .01, .000001, 1000000, 6, " µm",
        tr("Geometrische Toleranz für die nachgewiesene Rotation um die explizite Quellachse [µm], "
           "intern m. Keine automatische Achsschätzung oder Symmetrisierung."),
        content);
    rotationalTolerance_->setObjectName("axisymmetricRotationalTolerance");
    angularSectors_ = new IntegerInput;
    angularSectors_->setObjectName("axisymmetricSectors");
    angularSectors_->setRange(8, 4096);
    angularSectors_->setSingleStep(4);
    angularSectors_->setValue(64);
    angularSectors_->setToolTip(
        tr("Sektoren des gespeicherten 2D-Ergebnisnetzes, 8 bis 4096 und durch 4 teilbar. "
           "Mehr Sektoren verfeinern die Darstellung und vergrößern Projekt- und Exportdaten. "
           "Der kontinuierliche Rotationsoperator und seine Meridianunbekannten bleiben gleich; "
           "seine Winkelintegration wird unabhängig über das Fehlerbudget geprüft."));
    angularSectors_->setWhatsThis(angularSectors_->toolTip());
    axisymmetricControls_ = new QGroupBox(tr("Meridian und Winkelintegration"));
    axisymmetricControls_->setObjectName("axisymmetricControls");
    auto *meridianForm = new QFormLayout(axisymmetricControls_);
    meridianForm->setRowWrapPolicy(QFormLayout::WrapLongRows);
    for (const auto &[label, widget] :
         std::array<std::pair<QString, QWidget *>, 4>{{{tr("Kante"), meridianLength_},
                                                       {tr("Sehne"), meridianTolerance_},
                                                       {tr("Symmetrie"), rotationalTolerance_},
                                                       {tr("Netzsektoren"), angularSectors_}}}) {
        auto *caption = new QLabel(label);
        caption->setBuddy(widget);
        meridianForm->addRow(caption, widget);
    }
    geometryControls->addStretch();
    layout->addLayout(geometryControls);
    layout->addWidget(axisymmetricControls_);
    connect(geometryMode_, &QComboBox::currentIndexChanged, this, [this] {
        updateDomainControls();
        controlsChanged();
    });
    for (auto *input : {meridianLength_, meridianTolerance_, rotationalTolerance_})
        connect(input, &QDoubleSpinBox::valueChanged, this, [this] { controlsChanged(); });
    connect(angularSectors_, &QSpinBox::valueChanged, this, [this] { controlsChanged(); });
    nearMaximumDepth_ = new IntegerInput;
    nearMaximumDepth_->setObjectName("nearMaximumDepth");
    nearMaximumDepth_->setRange(0, 10);
    nearMaximumDepth_->setValue(6);
    nearMaximumDepth_->setToolTip(
        tr("Maximale Anzahl adaptiver Unterteilungen für nahe Dreieckspaare (0 bis 10). "
           "Höhere Werte ermöglichen die Integration sehr kleiner Flächenabstände und können "
           "mehr Rechenzeit benötigen. Das relative Nahfehlerziel bleibt 1e-4. "
           "Gilt für dichte Lösung und Fast Multipole."));
    nearMaximumDepth_->setWhatsThis(nearMaximumDepth_->toolTip());
    auto *nearControls = new QHBoxLayout;
    auto *nearLabel = new QLabel(tr("Nahquadratur: max. Tiefe"));
    nearLabel->setBuddy(nearMaximumDepth_);
    nearControls->addWidget(nearLabel);
    nearControls->addWidget(nearMaximumDepth_);
    nearControls->addStretch();
    layout->addLayout(nearControls);
    auto *mediumGroup = new QGroupBox(tr("Medium / Antrieb"));
    auto *mediumLayout = new QGridLayout(mediumGroup);
    mediumLayout->setContentsMargins(7, 10, 7, 5);
    density_ = number(1.21, .01, 100, 3, " kg/m³", tr("Dichte des akustischen Mediums in kg/m³."),
                      content);
    soundSpeed_ =
        number(343, 1, 10000, 1, " m/s", tr("Schallgeschwindigkeit des Mediums in m/s."), content);
    velocity_ =
        number(1, 0, 100, 4, " ×",
               tr("Dimensionsloser Faktor für die zugeordneten komplexen "
                  "RMS-Normalgeschwindigkeiten. Positive Normalrichtung: aus dem festen Körper in "
                  "das Schallfeld. Die Flächenzuordnung enthält die Basisgeschwindigkeit in m/s."),
               content);
    phase_ = number(0, -180, 180, 1, " °",
                    tr("Zusätzliche Phase der Anregung. RMS-Zeitkonvention exp(+iωt), positive "
                       "Normalrichtung aus dem Körper. Einheit: Grad."),
                    content);
    velocity_->setObjectName("prescribedVelocityScale");
    phase_->setObjectName("prescribedVelocityPhase");
    mediumLayout->addWidget(new QLabel(tr("Dichte ρ")), 0, 0);
    mediumLayout->addWidget(density_, 0, 1);
    mediumLayout->addWidget(new QLabel(tr("Schallgeschw. c")), 0, 2);
    mediumLayout->addWidget(soundSpeed_, 0, 3);
    mediumLayout->addWidget(new QLabel(tr("Antriebsfaktor")), 1, 0);
    mediumLayout->addWidget(velocity_, 1, 1);
    mediumLayout->addWidget(new QLabel(tr("Phase")), 1, 2);
    mediumLayout->addWidget(phase_, 1, 3);
    density_->setAccessibleName(tr("Dichte"));
    soundSpeed_->setAccessibleName(tr("Schallgeschwindigkeit"));
    velocity_->setAccessibleName(tr("Antriebsfaktor"));
    phase_->setAccessibleName(tr("Antriebsphase"));
    excitation_ = new QLabel(tr("Kein Antrieb"));
    mediumLayout->addWidget(excitation_, 2, 0, 1, 2);
    auto *frameButton = new QPushButton(tr("Abstrahlung …"));
    frameButton->setObjectName("acousticFrame");
    frameButton->setToolTip(tr("Ursprung, 0°-Richtung und H/V-Schnitte in Modellkoordinaten. "
                               "Unabhängig von der Kamera; Abstand ab diesem Ursprung."));
    mediumLayout->addWidget(frameButton, 2, 2, 1, 2);
    connect(frameButton, &QPushButton::clicked, this, [this] {
        AcousticFrameDialog dialog(defaultAcousticFrame(origin_, axis_), observationFrame_, this);
        if (dialog.exec() == QDialog::Accepted)
            setObservationFrame(dialog.observationFrame());
    });
    layout->addWidget(mediumGroup);
    auto *fieldControls = new QHBoxLayout;
    fieldMode_ = new QComboBox;
    fieldMode_->setObjectName("observationFieldMode");
    fieldMode_->addItems({tr("Fernfeld · auf Abstand normiert"), tr("Endlicher Messabstand")});
    fieldMode_->setToolTip(tr("Fernfeld nutzt die asymptotische Richtwirkung. Endlicher Abstand "
                              "berechnet den Druck an den tatsächlichen Messpunkten."));
    observationDistance_ = number(1, .001, 100000, 3, " m",
                                  tr("Abstand ab dem Abstrahlungsursprung. Bei endlichem Abstand "
                                     "müssen alle Messpunkte im äußeren Schallfeld liegen."),
                                  content);
    observationDistance_->setObjectName("observationDistance");
    fieldControls->addWidget(new QLabel(tr("Auswertung")));
    fieldControls->addWidget(fieldMode_);
    fieldControls->addWidget(observationDistance_);
    fieldControls->addStretch();
    layout->addLayout(fieldControls);
    connect(fieldMode_, &QComboBox::currentIndexChanged, this, [this] { controlsChanged(); });
    connect(observationDistance_, &QDoubleSpinBox::valueChanged, this,
            [this] { controlsChanged(); });
    layout->addStretch();
    auto *settingsScroll = new QScrollArea;
    settingsScroll->setObjectName("simulationSettingsScroll");
    settingsScroll->setFrameShape(QFrame::NoFrame);
    settingsScroll->setWidgetResizable(true);
    settingsScroll->setWidget(content);
    workspaceTabs_->addTab(settingsScroll, tr("Berechnung"));
    resultPage_ = new QWidget;
    auto *resultLayout = new QVBoxLayout(resultPage_);
    resultLayout->setContentsMargins(4, 4, 4, 4);
    resultLayout->setSpacing(4);
    resultSelectorBar_ = new QWidget;
    auto *selectorBar = new QHBoxLayout(resultSelectorBar_);
    selectorBar->setContentsMargins(0, 0, 0, 0);
    frequencyControls_ = new QWidget;
    frequencyControls_->setObjectName("polarFrequencyControls");
    auto *frequencyControls = new QHBoxLayout(frequencyControls_);
    frequencyControls->setContentsMargins(0, 0, 0, 0);
    frequencyControls->addWidget(new QLabel(tr("Schnitt bei")));
    frequencySelector_ = new QComboBox;
    frequencySelector_->setObjectName("resultFrequency");
    frequencySelector_->setMinimumWidth(105);
    frequencySelector_->setAccessibleName(tr("Polarfrequenz"));
    frequencyControls->addWidget(frequencySelector_);
    selectorBar->addWidget(frequencyControls_);
    impedanceControls_ = new QWidget;
    impedanceControls_->setObjectName("impedanceControls");
    auto *impedanceControls = new QHBoxLayout(impedanceControls_);
    impedanceControls->setContentsMargins(0, 0, 0, 0);
    impedanceControls->addWidget(new QLabel(tr("Impedanzfläche")));
    impedanceSelector_ = new QComboBox;
    impedanceSelector_->setMinimumWidth(100);
    impedanceSelector_->setObjectName("impedanceFace");
    impedanceSelector_->setAccessibleName(tr("Impedanzfläche"));
    impedanceControls->addWidget(impedanceSelector_);
    selectorBar->addWidget(impedanceControls_);
    normalizationControls_ = new QWidget;
    normalizationControls_->setObjectName("polarNormalizationControls");
    auto *normalizationControls = new QHBoxLayout(normalizationControls_);
    normalizationControls->setContentsMargins(0, 0, 0, 0);
    normalizationControls->addWidget(new QLabel(tr("Normierung")));
    normalization_ = new QComboBox;
    normalization_->setObjectName("polarNormalization");
    normalization_->addItems({tr("Winkel"), tr("Maximum")});
    normalization_->setToolTip(
        tr("Winkel: Pegel relativ zur gewählten Richtung derselben Frequenz, je H/V-Schnitt. "
           "Maximum: relativ zum Maximum des jeweiligen Schnitts. "
           "Gilt für Polars und Isobaren; Rohdaten bleiben erhalten. "
           "Bei null Referenzdruck wird kein normierter Schnitt angezeigt."));
    normalizationControls->addWidget(normalization_);
    normalizationAngle_ = number(0, -180, 180, 2, "°",
                                 tr("Referenzwinkel im jeweiligen H/V-Schnitt. Zwischen "
                                    "berechneten Winkeln wird der Pegel linear in dB "
                                    "interpoliert. 0° verwendet den Achsdruck."),
                                 resultPage_);
    normalizationAngle_->setObjectName("polarNormalizationAngle");
    normalizationControls->addWidget(normalizationAngle_);
    selectorBar->addWidget(normalizationControls_);
    selectorBar->addStretch();
    auto *measurementButton = new QPushButton(tr("Messung …"));
    measurementButton->setObjectName("measurementMenu");
    auto *measurementMenu = new QMenu(measurementButton);
    measurementMenu->addAction(tr("CSV / TXT importieren …"), this,
                               &SimulationDock::importMeasurement);
    measurementMenu->addAction(tr("Pegelversatz …"), this, [this] {
        if (measurements_.empty())
            return;
        QStringList names;
        for (const auto &curve : measurements_)
            names.append(curve.name);
        bool ok = false;
        const auto name =
            QInputDialog::getItem(this, tr("Messung"), tr("Kurve"), names, 0, false, &ok);
        if (!ok)
            return;
        auto &curve = measurements_[static_cast<std::size_t>(names.indexOf(name))];
        const auto offset = QInputDialog::getDouble(this, tr("Messung"), tr("Pegelversatz [dB]"),
                                                    curve.offsetDb, -300, 300, 2, &ok);
        if (!ok)
            return;
        curve.offsetDb = offset;
        updateFrequencyPlots();
        if (onChanged)
            onChanged();
    });
    measurementMenu->addAction(tr("Messkurven entfernen"), this, [this] {
        measurements_.clear();
        updateFrequencyPlots();
        if (onChanged)
            onChanged();
    });
    measurementButton->setMenu(measurementMenu);
    measurementButton->setToolTip(tr(
        "Messdatei: Hz und dB, optional Phase; CSV mit Komma/Strichpunkt oder TXT mit Leerzeichen. "
        "Messabstand und Anregung müssen für einen absoluten Vergleich übereinstimmen. "
        "Ein Pegelversatz verändert nur die Anzeige. Phase wird nicht importiert."));
    selectorBar->addWidget(measurementButton);
    auto *exportButton = new QPushButton(tr("CSV …"));
    exportButton->setToolTip(tr(
        "Echte komplexe Ergebnisdaten exportieren. Anzeigeabschätzungen werden nicht exportiert. "
        "Diagramme: Rechtsklick → PNG exportieren."));
    connect(exportButton, &QPushButton::clicked, this, &SimulationDock::exportResults);
    selectorBar->addWidget(exportButton);
    resultLayout->addWidget(resultSelectorBar_);
    directivityControls_ = new QWidget;
    directivityControls_->setObjectName("directivityScaleControls");
    auto *scaleBar = new QHBoxLayout(directivityControls_);
    scaleBar->setContentsMargins(0, 0, 0, 0);
    directivityMinimum_ = new QLineEdit("-30", resultPage_);
    directivityMaximum_ = new QLineEdit("3", resultPage_);
    directivityStep_ = new QLineEdit("3", resultPage_);
    directivityMinimum_->setToolTip(tr("Untere Grenze für Polars und Isobaren in dB."));
    directivityMaximum_->setToolTip(tr("Obere Grenze für Polars und Isobaren in dB. +3 dB zeigt "
                                       "Überhöhungen neben der Referenzrichtung."));
    directivityStep_->setToolTip(tr("Pegelabstand der Polarringe und Farbstufen in dB."));
    for (auto *input : {directivityMinimum_, directivityMaximum_, directivityStep_}) {
        input->setMaximumWidth(90);
        input->setAlignment(Qt::AlignRight);
    }
    directivityMinimum_->setObjectName("directivityMinimumDb");
    directivityMaximum_->setObjectName("directivityMaximumDb");
    directivityStep_->setObjectName("directivityStepDb");
    directivityMinimum_->setAccessibleName(tr("Polars / Isobaren: Minimum in dB"));
    directivityMaximum_->setAccessibleName(tr("Polars / Isobaren: Maximum in dB"));
    directivityStep_->setAccessibleName(tr("Polars / Isobaren: Schritt in dB"));
    scaleBar->addWidget(new QLabel(tr("Polars / Isobaren")));
    scaleBar->addWidget(new QLabel(tr("Min. [dB]")));
    scaleBar->addWidget(directivityMinimum_);
    scaleBar->addWidget(new QLabel(tr("Max. [dB]")));
    scaleBar->addWidget(directivityMaximum_);
    scaleBar->addWidget(new QLabel(tr("Schritt [dB]")));
    scaleBar->addWidget(directivityStep_);
    scaleBar->addStretch();
    resultLayout->addWidget(directivityControls_);
    auto *tabs = new QTabWidget;
    resultTabs_ = tabs;
    tabs->setObjectName("acousticResults");
    // A page owns its own overflow. A wide inactive Schallfeld form must not
    // enlarge the H/V plots on every other result page.
    const auto addResultPage = [tabs](QWidget *page, const QString &title) {
        auto *scroll = new QScrollArea;
        scroll->setFrameShape(QFrame::NoFrame);
        scroll->setWidgetResizable(true);
        scroll->setWidget(page);
        tabs->addTab(scroll, title);
    };
    auto *responses = new QWidget;
    auto *responseLayout = new QGridLayout(responses);
    responseLayout->setContentsMargins(0, 0, 0, 0);
    onAxis_ = new PlotWidget(tr("On-Axis"), tr("Frequenz [Hz]"), tr("dB SPL"));
    onAxis_->setFrequencyAxis();
    onAxis_->setObjectName("acousticOnAxis");
    onAxis_->setToolTip(
        tr("RMS-Schalldruckpegel auf der akustischen Achse, Referenz 20 µPa. Neue Rechnungen "
           "nutzen den gewählten Feldmodus und Abstand ab dem Abstrahlungsursprung."));
    impedance_ =
        new PlotWidget(tr("Normierte Strahlungsimpedanz"), tr("Frequenz [Hz]"), tr("Z S/(ρc) [1]"));
    impedance_->setFrequencyAxis();
    impedance_->setObjectName("acousticImpedance");
    impedance_->setToolTip(
        tr("Dimensionslose Impedanz Z·S/(ρc). Z ist der flächengemittelte Druck geteilt durch "
           "den Volumenstrom, S die zugehörige Driven-Fläche. Dichte und Schallgeschwindigkeit "
           "stammen aus dem ursprünglichen berechneten Request."));
    responseLayout->addWidget(onAxis_, 0, 0);
    responseLayout->addWidget(impedance_, 0, 1);
    addResultPage(responses, tr("Frequenzgänge"));
    auto *cta = new QWidget;
    auto *ctaLayout = new QGridLayout(cta);
    ctaLayout->setContentsMargins(0, 0, 0, 0);
    auto *ctaControls = new QHBoxLayout;
    ctaControls->addWidget(new QLabel(tr("Anregung")));
    ctaExcitation_ = new QComboBox;
    ctaExcitation_->setObjectName("ctaExcitation");
    ctaControls->addWidget(ctaExcitation_);
    auto *ctaExport = new QPushButton(tr("CTA-CSV …"));
    ctaExport->setObjectName("exportCta2034");
    ctaControls->addWidget(ctaExport);
    ctaControls->addStretch();
    ctaLayout->addLayout(ctaControls, 0, 0, 1, 2);
    connect(ctaExcitation_, &QComboBox::currentIndexChanged, this, [this](int index) {
        ctaSelection_ = index;
        updateCta2034(hasCurrentResult());
        if (onChanged) onChanged();
    });
    connect(ctaExport, &QPushButton::clicked, this, [this] {
        if (ctaExcitation_->currentIndex() < 0 || !angularSampleCount(resultRequest_)) return;
        const auto path = QFileDialog::getSaveFileName(this, tr("CTA-2034 exportieren"),
                                                        "HornSculpt-CTA2034.csv", tr("CSV (*.csv)"));
        if (path.isEmpty()) return;
        try {
            const auto bytes = ctaCsv(archive(), static_cast<std::size_t>(ctaExcitation_->currentIndex()));
            QSaveFile file(path);
            if (!file.open(QIODevice::WriteOnly) || file.write(bytes) != bytes.size() || !file.commit())
                throw std::runtime_error(file.errorString().toStdString());
        } catch (const std::exception &e) {
            QMessageBox::warning(this, tr("CTA-CSV"), QString::fromUtf8(e.what()));
        }
    });
    ctaStatus_ = new QLabel;
    ctaStatus_->setWordWrap(true);
    ctaStatus_->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard);
    ctaLayout->addWidget(ctaStatus_, 1, 0, 1, 2);
    ctaLevels_ = new PlotWidget(tr("CTA-2034 · Richtungspegel"), tr("Frequenz [Hz]"), tr("dB SPL"));
    ctaDi_ = new PlotWidget(tr("Directivity Index"), tr("Frequenz [Hz]"), tr("dB"));
    ctaLevels_->setObjectName("cta2034Levels");
    ctaDi_->setObjectName("cta2034EarlyReflectionsDi");
    ctaLevels_->setFrequencyAxis();
    ctaDi_->setFrequencyAxis();
    ctaLayout->addWidget(ctaLevels_, 2, 0);
    ctaLayout->addWidget(ctaDi_, 2, 1);
    addResultPage(cta, tr("CTA-2034"));
    auto *beamwidths = new QWidget;
    auto *beamwidthLayout = new QGridLayout(beamwidths);
    beamwidthThreshold_ = new QComboBox;
    beamwidthThreshold_->addItem(tr("−6 dB · Achse"), 6.0);
    beamwidthThreshold_->addItem(tr("−3 dB · Achse"), 3.0);
    beamwidthThreshold_->addItem(tr("−10 dB · Achse"), 10.0);
    beamwidthThreshold_->setToolTip(tr("Volle Breite der zusammenhängenden Keule um 0°. "
                                       "Erste Schwellenübertritte relativ zum Achspegel, linear "
                                       "zwischen 5°-Messpunkten interpoliert. "
                                       "Ohne beidseitigen Übertritt bleibt eine Lücke."));
    auto *widthControls = new QHBoxLayout;
    widthControls->addWidget(beamwidthThreshold_);
    beamwidthTarget_ = new QCheckBox(tr("Ziel"));
    targetH_ = number(
        90, 0, 360, 1, " °",
        tr("Horizontaler voller Zielwinkel. Nur Vergleichsreferenz, verändert keine Geometrie."),
        beamwidths);
    targetV_ = number(60, 0, 360, 1, " °",
                      tr("Vertikaler voller Zielwinkel. Nur Vergleichsreferenz."), beamwidths);
    targetTolerance_ =
        number(10, 0, 180, 1, " °",
               tr("Toleranz um den Zielwinkel; sichtbarer Korridor auf 0° bis 360° begrenzt."),
               beamwidths);
    widthControls->addWidget(beamwidthTarget_);
    widthControls->addWidget(new QLabel("H"));
    widthControls->addWidget(targetH_);
    widthControls->addWidget(new QLabel("V"));
    widthControls->addWidget(targetV_);
    widthControls->addWidget(new QLabel("±"));
    widthControls->addWidget(targetTolerance_);
    widthControls->addStretch();
    beamwidthLayout->addLayout(widthControls, 0, 0, 1, 2);
    beamwidthH_ =
        new PlotWidget(tr("Abstrahlwinkel H"), tr("Frequenz [Hz]"), tr("Voller Winkel [°]"));
    beamwidthV_ =
        new PlotWidget(tr("Abstrahlwinkel V"), tr("Frequenz [Hz]"), tr("Voller Winkel [°]"));
    beamwidthH_->setFrequencyAxis();
    beamwidthV_->setFrequencyAxis();
    beamwidthH_->setAxisScales({20, 20000, {}}, {0, 180, 30});
    beamwidthV_->setAxisScales({20, 20000, {}}, {0, 180, 30});
    beamwidthH_->setObjectName("acousticBeamwidthH");
    beamwidthV_->setObjectName("acousticBeamwidthV");
    beamwidthLayout->addWidget(beamwidthH_, 1, 0);
    beamwidthLayout->addWidget(beamwidthV_, 1, 1);
    connect(beamwidthThreshold_, &QComboBox::currentIndexChanged, this,
            [this] { updateFrequencyPlots(); });
    const auto targetsChanged = [this] {
        for (auto *control : {targetH_, targetV_, targetTolerance_})
            control->setEnabled(beamwidthTarget_->isChecked());
        updateFrequencyPlots();
        if (onChanged)
            onChanged();
    };
    connect(beamwidthTarget_, &QCheckBox::toggled, this, targetsChanged);
    for (auto *control : {targetH_, targetV_, targetTolerance_}) {
        control->setEnabled(false);
        connect(control, &QDoubleSpinBox::valueChanged, this, targetsChanged);
    }
    auto *cuts = new QWidget;
    auto *cutLayout = new QGridLayout(cuts);
    cutLayout->setContentsMargins(0, 0, 0, 0);
    polarH_ = new PolarPlotWidget(tr("Polar H"), QColor(28, 102, 180));
    polarV_ = new PolarPlotWidget(tr("Polar V"), QColor(205, 112, 33));
    polarH_->setObjectName("acousticPolarH");
    polarV_->setObjectName("acousticPolarV");
    cutLayout->addWidget(polarH_, 0, 0);
    cutLayout->addWidget(polarV_, 0, 1);
    addResultPage(cuts, tr("Polars"));
    auto *maps = new QWidget;
    auto *mapLayout = new QGridLayout(maps);
    mapLayout->setContentsMargins(0, 0, 0, 0);
    isobarH_ = new HeatmapWidget(tr("Isobaren H"));
    isobarV_ = new HeatmapWidget(tr("Isobaren V"));
    isobarH_->setObjectName("acousticIsobarH");

    isobarV_->setObjectName("acousticIsobarV");

    mapLayout->addWidget(isobarH_, 0, 0);
    mapLayout->addWidget(isobarV_, 0, 1);
    addResultPage(maps, tr("Isobaren"));
    auto *radialMaps = new QWidget;
    auto *radialLayout = new QGridLayout(radialMaps);
    radialLayout->setContentsMargins(0, 0, 0, 0);
    radialIsobarH_ = new HeatmapWidget(tr("Isobaren rund H"));
    radialIsobarV_ = new HeatmapWidget(tr("Isobaren rund V"));
    radialIsobarH_->setObjectName("acousticRadialIsobarH");
    radialIsobarV_->setObjectName("acousticRadialIsobarV");
    radialIsobarH_->setRadialLayout(true);
    radialIsobarV_->setRadialLayout(true);
    radialLayout->addWidget(radialIsobarH_, 0, 0);
    radialLayout->addWidget(radialIsobarV_, 0, 1);
    angularResponses_ = new AngularResponseWidget;
    angularResponses_->onChanged = [this] {
        if (onChanged)
            onChanged();
    };
    addResultPage(angularResponses_, tr("Einzelfrequenzgänge"));
    phaseResponses_ = new PhaseResponseWidget;
    acousticField_ = new AcousticFieldWidget;
    addResultPage(phaseResponses_, tr("Phase"));
    // The field's own detail area scrolls; keeping its core outside a second
    // scroll viewport leaves selection and actions visible while editing details.
    auto *fieldPage = new QWidget;
    auto *fieldLayout = new QVBoxLayout(fieldPage);
    fieldLayout->setContentsMargins(0, 0, 0, 0);
    fieldLayout->addWidget(acousticField_);
    tabs->addTab(fieldPage, tr("Schallfeld"));
    for (auto *widget :
         {static_cast<QWidget *>(phaseResponses_), static_cast<QWidget *>(acousticField_)})
        widget->setToolTip(tr("Darstellung aus einem gespeicherten, komplexen BEM-Ergebnis. "
                              "Anzeigeänderungen starten keine neue Simulation."));
    phaseResponses_->onChanged = [this] {
        if (onChanged)
            onChanged();
    };
    acousticField_->onChanged = [this] {
        if (onChanged)
            onChanged();
    };
    acousticField_->onFieldChanged = [this](std::shared_ptr<const PressureFieldPlane> plane,
                                            PressureFieldAppearance appearance) {
        if (onFieldChanged)
            onFieldChanged(std::move(plane), appearance);
    };
    acousticField_->onFocusField = [this] {
        if (onFocusField)
            onFocusField();
    };
    acousticField_->onProbeModeChanged = [this](bool enabled) {
        if (onFieldProbeModeChanged)
            onFieldProbeModeChanged(enabled);
    };
    addResultPage(beamwidths, tr("Abstrahlwinkel"));
    addResultPage(radialMaps, tr("Isobaren rund"));
    // Keep the existing pages and their identities, but avoid a row of tabs
    // that overflows before users can reach the later result views.
    auto *navigation = new QWidget;
    navigation->setObjectName("resultNavigation");
    auto *navigationLayout = new QHBoxLayout(navigation);
    navigationLayout->setContentsMargins(0, 0, 0, 0);
    auto *family = new QComboBox;
    family->setObjectName("resultFamily");
    family->setAccessibleName(tr("Ergebnisbereich"));
    family->addItems({tr("Übersicht"), tr("Abstrahlung"), tr("Feld und Phase")});
    auto *view = new QComboBox;
    view->setObjectName("resultView");
    view->setAccessibleName(tr("Ergebnisansicht"));
    navigationLayout->addWidget(new QLabel(tr("Bereich")));
    navigationLayout->addWidget(family);
    navigationLayout->addWidget(new QLabel(tr("Ansicht")));
    navigationLayout->addWidget(view, 1);
    quality_ = new QToolButton;
    quality_->setObjectName("resultQuality");
    quality_->setText(tr("Qualität"));
    navigationLayout->addWidget(quality_);
    connect(quality_, &QToolButton::clicked, this, [this] {
        QMessageBox dialog(QMessageBox::Information, tr("Ergebnisqualität"),
                           quality_->toolTip(), QMessageBox::Close, this);
        auto *convergence = dialog.addButton(tr("Netzkonvergenz …"), QMessageBox::ActionRole);
        dialog.exec();
        if (dialog.clickedButton() == convergence && onRequestConvergence)
            onRequestConvergence();
    });
    resultLayout->insertWidget(0, navigation);
    tabs->tabBar()->hide();
    const auto familyForPage = [tabs, responses, cta, phase = phaseResponses_,
                                field = acousticField_](int index) {
        const auto *page = tabs->widget(index);
        if (!page)
            return 0;
        if (page == responses || page->isAncestorOf(responses) || page == cta ||
            page->isAncestorOf(cta))
            return 0;
        return page == phase || page->isAncestorOf(phase) || page == field ||
                       page->isAncestorOf(field) ? 2 : 1;
    };
    const auto populateViews = [tabs, family, view, familyForPage](int selected) {
        const QSignalBlocker blocker(view);
        view->clear();
        for (int i = 0; i < tabs->count(); ++i)
            if (familyForPage(i) == family->currentIndex())
                view->addItem(tabs->tabText(i), i);
        const int index = view->findData(selected);
        view->setCurrentIndex(index >= 0 ? index : 0);
    };
    connect(family, &QComboBox::currentIndexChanged, this, [tabs, view, populateViews](int) {
        populateViews(-1);
        if (view->currentIndex() >= 0)
            tabs->setCurrentIndex(view->currentData().toInt());
    });
    connect(view, &QComboBox::currentIndexChanged, this, [tabs, view](int index) {
        if (index >= 0)
            tabs->setCurrentIndex(view->currentData().toInt());
    });
    connect(tabs, &QTabWidget::currentChanged, this,
            [family, populateViews, familyForPage](int index) {
                const QSignalBlocker blocker(family);
                family->setCurrentIndex(familyForPage(index));
                populateViews(index);
            });
    populateViews(tabs->currentIndex());
    emptyResults_ = new QLabel;
    emptyResults_->setObjectName("simulationEmptyResults");
    emptyResults_->setWordWrap(true);
    emptyResults_->setAlignment(Qt::AlignCenter);
    resultLayout->addWidget(emptyResults_, 1);
    resultLayout->addWidget(tabs, 1);
    workspaceTabs_->addTab(resultPage_, tr("Ergebnisse"));
    rootLayout->addWidget(workspaceTabs_, 1);
    auto *statusBar = new QHBoxLayout;
    status_ = new QLabel(tr("Netz erforderlich"));
    status_->setObjectName("simulationStatus");
    status_->setWordWrap(true);
    status_->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Preferred);
    status_->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard);
    status_->setFocusPolicy(Qt::StrongFocus);
    statusBar->addWidget(status_, 1);
    progress_ = new QProgressBar;
    progress_->setObjectName("simulationProgress");
    progress_->setRange(0, 17);
    progress_->setValue(0);
    progress_->setVisible(false);
    progress_->setMaximumWidth(150);
    statusBar->addWidget(progress_);
    rootLayout->addLayout(statusBar);
    setWidget(root);
    connect(tabs, &QTabWidget::currentChanged, this, [this] {
        updateResultControls();
        QTimer::singleShot(0, this, [this] { updateResultControls(); });
    });
    connect(workspaceTabs_, &QTabWidget::currentChanged, this, [this] { updateResultControls(); });
    for (auto *spin :
         {minimumFrequency_, maximumFrequency_, density_, soundSpeed_, velocity_, phase_})
        connect(spin, &QDoubleSpinBox::valueChanged, this, [this] { controlsChanged(); });
    connect(frequencyCount_, &QSpinBox::valueChanged, this, [this](int count) {
        maximumFrequency_->setEnabled(count > 1);
        controlsChanged();
    });
    connect(fastMultipole_, &QCheckBox::toggled, this, [this](bool enabled) {
        multipoleOrder_->setEnabled(enabled && !infiniteBaffle_);
        updateDomainControls();
        controlsChanged();
    });
    connect(multipoleOrder_, &QSpinBox::valueChanged, this, [this](int order) {
        fastMultipoleSettings_.interpolationOrder = static_cast<unsigned>(order);
        controlsChanged();
    });
    connect(nearMaximumDepth_, &QSpinBox::valueChanged, this, [this] { controlsChanged(); });
    connect(solve_, &QPushButton::clicked, this, &SimulationDock::requestSolve);
    realtimeTimer_ = new QTimer(this);
    realtimeTimer_->setSingleShot(true);
    realtimeTimer_->setInterval(200);
    connect(realtime_, &QPushButton::toggled, this, &SimulationDock::setRealtimeEnabled);
    connect(realtimeTimer_, &QTimer::timeout, this, [this] {
        if (realtimeEnabled_ && realtimePending_ && !busy_ && hasCurrentGeometry() &&
            hasExcitation())
            requestSolve();
    });
    connect(stop_, &QPushButton::clicked, this, &SimulationDock::cancelSolve);
    connect(frequencySelector_, &QComboBox::currentIndexChanged, this,
            [this] { updateAngularPlots(); });
    connect(impedanceSelector_, &QComboBox::currentIndexChanged, this,
            [this] { updateFrequencyPlots(); });
    const auto normalizationChanged = [this] {
        normalizationAngle_->setEnabled(normalization_->currentIndex() == 0);
        updateAngularPlots();
        updateHeatmaps();
        if (onChanged)
            onChanged();
    };
    connect(normalization_, &QComboBox::currentIndexChanged, this, normalizationChanged);
    connect(normalizationAngle_, &QDoubleSpinBox::valueChanged, this, normalizationChanged);
    const auto scaleChanged = [this] {
        updateDirectivityScale();
        if (onChanged)
            onChanged();
    };
    for (auto *control : {directivityMinimum_, directivityMaximum_, directivityStep_})
        connect(control, &QLineEdit::editingFinished, this, scaleChanged);
    const auto plotScaleChanged = [this](PlotDbScale scale) {
        applyDirectivityScale(scale);
        if (onChanged)
            onChanged();
    };
    polarH_->onDecibelScaleChanged = plotScaleChanged;
    polarV_->onDecibelScaleChanged = plotScaleChanged;
    isobarH_->onDecibelScaleChanged = plotScaleChanged;
    isobarV_->onDecibelScaleChanged = plotScaleChanged;
    radialIsobarH_->onDecibelScaleChanged = plotScaleChanged;
    radialIsobarV_->onDecibelScaleChanged = plotScaleChanged;
    updateDirectivityScale();
    updateDomainControls();
    updateResultControls();
    worker_ = std::jthread([this](std::stop_token shutdown) { runWorker(shutdown); });
}

SimulationDock::~SimulationDock() {
    acousticField_->cancel();
    // Child controls can emit while QWidget tears them down, after value members have died.
    for (auto *child : findChildren<QObject *>())
        QObject::disconnect(child, nullptr, this, nullptr);
    onChanged = {};
    onRealtimeChanged = {};
    onLemResultChanged = {};
    onFieldProbeModeChanged = {};
    onRequestConvergence = {};
    onFieldChanged = {};
    onFocusField = {};
    jobStop_.request_stop();
    worker_.request_stop();
    workerWake_.notify_all();
    if (worker_.joinable())
        worker_.join();
}
void SimulationDock::revealResultView() {
    workspaceTabs_->setCurrentWidget(resultPage_);
}

void SimulationDock::updateResultControls() {
    updateResultQuality();
    const auto *page = resultTabs_->currentWidget();
    const bool polarPage = page && page->isAncestorOf(polarH_);
    const bool isobarPage = page && (page->isAncestorOf(isobarH_) || page->isAncestorOf(radialIsobarH_));
    const bool responsePage = page && page->isAncestorOf(impedance_);
    const bool fieldPage = page && page->isAncestorOf(acousticField_);
    const bool empty = !result_ && measurements_.empty();
    const bool idleField = !empty && fieldPage && workspaceTabs_->currentWidget() == resultPage_ &&
                           !busy_ && !realtimeEnabled_;
    solverCommands_->setVisible(!idleField);
    const bool showingField = !empty && fieldPage && workspaceTabs_->currentWidget() == resultPage_;
    emptyResults_->setText(
        busy_ ? tr("Berechnung läuft. Ergebnisse erscheinen nach dem ersten Frequenzpunkt.")
              : tr("Noch keine berechneten Ergebnisse.\n"
                   "Unter „Berechnung“ die Einstellungen prüfen und mit „Berechnen“ starten.\n"
                   "Eine vorhandene Messkurve lässt sich über „Messung …“ importieren."));
    emptyResults_->setVisible(empty);
    resultTabs_->setVisible(!empty);
    resultSelectorBar_->setVisible(!fieldPage || empty);
    frequencyControls_->setVisible(!empty && polarPage);
    impedanceControls_->setVisible(!empty && responsePage);
    normalizationControls_->setVisible(!empty && (polarPage || isobarPage));
    directivityControls_->setVisible(!empty && (polarPage || isobarPage));
    const auto margins = resultPage_->layout()->contentsMargins();
    const auto *navigation = resultPage_->findChild<QWidget *>("resultNavigation");
    const int fieldMinimum = acousticField_->minimumSizeHint().height() +
                             workspaceTabs_->tabBar()->sizeHint().height() +
                             (navigation ? navigation->minimumSizeHint().height() : 0) + margins.top() +
                             margins.bottom() + resultPage_->layout()->spacing() +
                             4 * style()->pixelMetric(QStyle::PM_DefaultFrameWidth);
    workspaceTabs_->setMinimumHeight(showingField ? fieldMinimum : 100);
    setMinimumHeight(0);
    widget()->setMinimumHeight(0);
    if (showingField)
        widget()->setMinimumHeight(widget()->minimumSizeHint().height());
    widget()->updateGeometry();
    updateGeometry();
    if (showingField)
        setMinimumHeight(minimumSizeHint().height());
}

void SimulationDock::updateResultQuality() {
    if (!quality_)
        return;
    if (!result_) {
        quality_->setText(tr("Qualität · kein Ergebnis"));
        quality_->setToolTip(tr("Noch kein BEM-Ergebnis vorhanden."));
        return;
    }
    const bool current = hasCurrentGeometry() && resultRevision_ == currentRevision_ &&
                         resultGeneration_ == generation_;
    const QString freshness = current ? tr("aktuell") : tr("veraltet");
    const QString completeness = !resultCompletenessKnown_ ? tr("Vollständigkeit unbekannt")
                                 : resultComplete_ ? tr("vollständig") : tr("Teilraster");
    double maximumResidual = 0;
    bool haveResidual = false;
    for (const auto &frequency : result_->frequencies)
        for (const auto &excitation : frequency.excitations)
            if (std::isfinite(excitation.boundaryResidual) && excitation.boundaryResidual >= 0) {
                maximumResidual = std::max(maximumResidual, excitation.boundaryResidual);
                haveResidual = true;
            }
    quality_->setText(!current ? tr("Qualität: veraltet")
                      : resultCompletenessKnown_ && !resultComplete_ ? tr("Qualität: Teil")
                      : tr("Qualität: aktuell"));
    quality_->setToolTip(tr("Ergebnis Rev. %1 · %2 · %3\n%4\nNetz-/Quadraturkonvergenz "
                            "dieses Ergebnisses nicht geprüft. CTA-Schallleistung ist eine "
                            "70-Richtungs-Näherung, keine Zertifizierung.")
                             .arg(resultRevision_).arg(freshness, completeness)
                             .arg(haveResidual ? tr("Max. relatives Randresiduum: %1").arg(maximumResidual, 0, 'g', 4)
                                               : tr("Randresiduum nicht verfügbar")));
}

void SimulationDock::invalidate(std::uint64_t revision) {
    currentRevision_ = revision;
    meshCurrent_ = false;
    geometrySourceCurrent_ = false;
    ++generation_;
    updateResultQuality();
    jobStop_.request_stop();
    solve_->setEnabled(false);
    if (result_) {
        const auto completeness = !resultCompletenessKnown_ ? tr(" · Vollständigkeit unbekannt")
                                  : !resultComplete_ ? tr(" · Teilraster") : QString{};
        setWindowTitle(tr("Akustik · Rev. %1 · veraltet%2").arg(resultRevision_).arg(completeness));
        setStatus(tr("Geometrie geändert · Ergebnisse veraltet"));
    } else
        setStatus(tr("Netz erforderlich"));
    scheduleRealtime();
    updateAngularPlots();
    updateAngularFrequencyPlots(false);
    updateFieldAndPhase(false);
}

void SimulationDock::setMesh(std::shared_ptr<const SurfaceMesh> mesh, std::uint64_t revision,
                             std::vector<FaceVelocity> velocities) {
    if (currentRevision_ != revision || meshCurrent_ || geometrySourceCurrent_)
        invalidate(revision);
    mesh_ = std::move(mesh);
    meshRevision_ = revision;
    meshCurrent_ = static_cast<bool>(mesh_);
    velocities_ = std::move(velocities);
    if (mesh_ && velocities_.empty()) {
        std::set<std::uint32_t> driven;
        for (const auto &triangle : mesh_->triangles)
            if (triangle.boundary == BoundaryCondition::driven)
                driven.insert(triangle.faceId);
        for (const auto face : driven)
            velocities_.push_back({face, {1, 0}});
    }
    excitation_->setText(tr("%1 Driven").arg(velocities_.size()));
    QStringList definitions;
    for (const auto &v : velocities_)
        definitions.push_back(tr("Fläche %1: (%2 + i·%3) m/s RMS")
                                  .arg(v.faceId)
                                  .arg(v.normalVelocity.real(), 0, 'g', 6)
                                  .arg(v.normalVelocity.imag(), 0, 'g', 6));
    excitation_->setToolTip(definitions.join('\n'));
    updateDomainControls();
    solve_->setEnabled(!busy_ && mesh_ && hasExcitation());
    if (!result_)
        setStatus(mesh_
                      ? tr("Netz Rev. %1 · %2 Dreiecke").arg(revision).arg(mesh_->triangles.size())
                      : tr("Netz erforderlich"));
    scheduleRealtime();
}
void SimulationDock::setGeometrySource(std::shared_ptr<const DraftGeometry> draft,
                                       std::shared_ptr<const CadSnapshot> cad,
                                       std::uint64_t revision,
                                       std::map<std::uint32_t, double> faceTargets,
                                       std::set<std::uint32_t> reversedFaces) {
    const bool same = geometrySource_ && geometrySource_->draft == draft &&
                      geometrySource_->cad == cad && geometrySource_->revision == revision &&
                      geometrySource_->faceTargets == faceTargets &&
                      geometrySource_->reversedFaces == reversedFaces;
    if (!same) {
        jobStop_.request_stop();
        ++generation_;
        geometrySource_ = draft || cad ? std::make_shared<GeometrySource>(GeometrySource{
                                             std::move(draft), std::move(cad), revision,
                                             std::move(faceTargets), std::move(reversedFaces)})
                                       : nullptr;
    }
    geometrySourceCurrent_ = geometrySource_ && revision == currentRevision_;
    updateResultQuality();
    solve_->setEnabled(!busy_ && hasCurrentGeometry() && hasExcitation());
    if (usesAxisymmetricGeometry() && geometrySourceCurrent_ && !busy_ && !result_)
        setStatus(tr("2D-Quelle bereit · Meridianprüfung beim Start"));
    scheduleRealtime();
}
bool SimulationDock::usesAxisymmetricGeometry() const noexcept {
    return geometryMode_->currentIndex() == 1;
}
bool SimulationDock::hasCurrentGeometry() const {
    return usesAxisymmetricGeometry() ? geometrySourceCurrent_ && geometrySource_ &&
                                            geometrySource_->revision == currentRevision_
                                      : meshCurrent_ && mesh_ && meshRevision_ == currentRevision_;
}

void SimulationDock::setAcousticAxis(Vec3 origin, Vec3 axis, bool geometryAlreadyInvalidated) {
    const auto frame = defaultAcousticFrame(origin, axis);
    if (origin_ == frame.origin && axis_ == frame.forward)
        return;
    origin_ = frame.origin;
    axis_ = frame.forward;
    if (!observationFrame_) {
        if (!geometryAlreadyInvalidated)
            controlsChanged();
        if (onObservationFrameChanged)
            onObservationFrameChanged();
    }
}

void SimulationDock::setObservationFrame(std::optional<AcousticFrame> frame) {
    if (frame)
        frame = normalizedAcousticFrame(*frame);
    if (frame == observationFrame_)
        return;
    observationFrame_ = frame;
    controlsChanged();
    if (onObservationFrameChanged)
        onObservationFrameChanged();
}

AcousticFrame SimulationDock::observationFrame() const {
    return observationFrame_.value_or(defaultAcousticFrame(origin_, axis_));
}

void SimulationDock::setCoupling(lem::Network network, bool enabled) {
    if (network == couplingNetwork_ && enabled == couplingEnabled_)
        return;
    const bool physicsChanged =
        enabled != couplingEnabled_ ||
        ((enabled || hasBoundInterfacePorts(network) || hasBoundInterfacePorts(couplingNetwork_)) &&
         !lem::samePhysics(network, couplingNetwork_));
    couplingNetwork_ = std::move(network);
    couplingEnabled_ = enabled;
    updateDomainControls();
    if (physicsChanged)
        controlsChanged();
    else {
        // Layout and disabled circuits are still project edits, but cannot change this solve.
        if (onChanged)
            onChanged();
        updateAngularPlots();
    }
}

void SimulationDock::setInfiniteBaffle(std::optional<InfiniteBaffle> baffle) {
    if (baffle)
        validateInfiniteBaffleSettings(*baffle);
    if (baffle == infiniteBaffle_)
        return;
    infiniteBaffle_ = std::move(baffle);
    updateDomainControls();
    controlsChanged();
}

void SimulationDock::cancelSolve() {
    setRealtimeEnabled(false);
    jobStop_.request_stop();
    setStatus(tr("Berechnung wird abgebrochen …"));
}

bool SimulationDock::hasConfiguredExcitation() const { return hasExcitation(); }

bool SimulationDock::hasExcitation() const {
    return couplingEnabled_ || !velocities_.empty();
}

void SimulationDock::updateDomainControls() {
    const bool reduced = usesAxisymmetricGeometry();
    axisymmetricControls_->setVisible(reduced);
    for (auto *control : {meridianLength_, meridianTolerance_, rotationalTolerance_})
        control->setEnabled(reduced);
    angularSectors_->setEnabled(reduced);
    std::set<std::uint32_t> portFaces;
    for (const auto *element : lem::bemPorts(couplingNetwork_)) {
        const auto &port = std::get<lem::BemPort>(element->component);
        portFaces.insert(port.faceIds.begin(), port.faceIds.end());
    }
    const bool knownDrive =
        !couplingEnabled_ ||
        (infiniteBaffle_ &&
         std::any_of(velocities_.begin(), velocities_.end(),
                     [&](const auto &velocity) { return !portFaces.contains(velocity.faceId); }));
    velocity_->setEnabled(knownDrive);
    phase_->setEnabled(knownDrive);
    // A checked unsupported mode stays visible and can be explicitly switched off.
    fastMultipole_->setEnabled((!infiniteBaffle_ && !reduced) || fastMultipole_->isChecked());
    multipoleOrder_->setEnabled(!infiniteBaffle_ && !reduced && fastMultipole_->isChecked());
    baffle_->setText(infiniteBaffle_ ? tr("Schallwand · unendlich …") : tr("Schallwand …"));
    const auto orientation =
        infiniteBaffle_
            ? tr("RMS-Anregung in die Luftkavität, entgegen der aus der Luft gerichteten "
                 "Driven-Netznormale. Die Flächenzuordnung enthält die Basisgeschwindigkeit "
                 "in m/s.")
            : tr("RMS-Normalgeschwindigkeit aus dem festen Körper in das Schallfeld. "
                 "Die Flächenzuordnung enthält die Basisgeschwindigkeit in m/s.");
    velocity_->setToolTip(tr("Dimensionsloser Faktor für die komplexe ") + orientation +
                          tr(" Bei BEM↔LEM nur für Driven-Flächen ohne Portzuordnung; "
                             "Portgeschwindigkeiten bestimmt das Netzwerk."));
    velocity_->setWhatsThis(velocity_->toolTip());
    phase_->setToolTip(tr("Zusätzliche Phase [Grad], RMS-Zeitkonvention exp(+iωt). ") +
                       orientation);
    phase_->setWhatsThis(phase_->toolTip());
    if (reduced)
        fastMultipole_->setToolTip(
            tr("2D-Rotation (m=0) verwendet den dichten reduzierten Burton–Miller-Solver. "
               "Fast Multipole ausdrücklich ausschalten; es erfolgt kein Solverwechsel."));
    else if (infiniteBaffle_)
        fastMultipole_->setToolTip(
            tr("Fast Multipole ist für die unendliche Schallwand nicht unterstützt. "
               "Vor einer Berechnung ausdrücklich ausschalten; es erfolgt kein Solverwechsel."));
    else
        fastMultipole_->setToolTip(
            tr("Hierarchische Fernwechselwirkungen und iterative Burton–Miller-Lösung. "
               "Ausgeschaltet: dichte Matrix mit direkter Lösung."));
    fastMultipole_->setWhatsThis(fastMultipole_->toolTip());
    fieldMode_->setToolTip(
        infiniteBaffle_
            ? tr("Die H/V-Auswertung der unendlichen Schallwand unterstützt nur Fernfeld. "
                 "Der rückwärtige Halbraum ist abgeschirmt (Druck null). Endliche "
                 "H/V-Messkreise werden mit einer Fehlermeldung abgewiesen.")
            : tr("Fernfeld nutzt die asymptotische Richtwirkung. Endlicher Abstand berechnet "
                 "den Druck an den tatsächlichen Messpunkten."));
    fieldMode_->setWhatsThis(fieldMode_->toolTip());
}

void SimulationDock::setFrequencySweep(double minimum, double maximum, std::size_t count) {
    if (!std::isfinite(minimum) || !std::isfinite(maximum) || minimum < 1 || maximum > 100000 ||
        maximum < minimum || count < 1 ||
        count > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        (count > 1 && maximum == minimum))
        throw std::invalid_argument("Ungültiges Frequenzraster.");
    const QSignalBlocker first(minimumFrequency_), second(maximumFrequency_),
        third(frequencyCount_);
    minimumFrequency_->setValue(minimum);
    maximumFrequency_->setValue(maximum);
    frequencyCount_->setValue(static_cast<int>(count));
    maximumFrequency_->setEnabled(count > 1);
    controlsChanged();
}

void SimulationDock::controlsChanged() {
    ++generation_;
    updateResultQuality();
    jobStop_.request_stop();
    if (result_)
        setWindowTitle(tr("Akustik · Rev. %1 · veraltet").arg(resultRevision_));
    solve_->setEnabled(!busy_ && hasCurrentGeometry() && hasExcitation());
    if (hasCurrentGeometry())
        setStatus(result_ ? tr("Einstellungen geändert · Ergebnisse veraltet") : tr("Bereit"));
    if (onChanged)
        onChanged();
    scheduleRealtime();
    updateAngularPlots();
    updateAngularFrequencyPlots(false);
    updateFieldAndPhase(false);
}
void SimulationDock::setRealtimeEnabled(bool enabled) {
    if (externalJob_ && enabled) {
        const QSignalBlocker blocker(realtime_);
        realtime_->setChecked(false);
        return;
    }
    if (realtimeEnabled_ == enabled)
        return;
    realtimeEnabled_ = enabled;
    const QSignalBlocker blocker(realtime_);
    realtime_->setChecked(enabled);
    realtimeTimer_->stop();
    if (!enabled)
        jobStop_.request_stop();
    realtimePending_ = enabled && !hasCurrentResult();
    if (realtimePending_)
        realtimeTimer_->start();
    if (onRealtimeChanged)
        onRealtimeChanged(enabled);
    updateResultControls();
}
bool SimulationDock::realtimeEnabled() const noexcept {
    return realtimeEnabled_;
}
void SimulationDock::scheduleRealtime() {
    if (!realtimeEnabled_)
        return;
    realtimePending_ = true;
    realtimeTimer_->start();
}
bool SimulationDock::isBusy() const noexcept {
    return busy_;
}
bool SimulationDock::hasCurrentResult() const noexcept {
    return hasCurrentGeometry() && result_ && resultRevision_ == currentRevision_ &&
           resultGeneration_ == generation_ && resultComplete_ && !busy_;
}
QString SimulationDock::lastError() const {
    return error_;
}
void SimulationDock::beginExternalJob() {
    const bool current = hasCurrentResult();
    setRealtimeEnabled(false);
    cancelPendingJobs();
    if (current)
        resultGeneration_ = generation_;
    updateResultQuality();
    externalJob_ = true;
    solve_->setEnabled(false);
    {
        std::lock_guard lock(workerMutex_);
        releaseGpuRequested_ = true;
        gpuReleased_.store(false);
    }
    workerWake_.notify_one();
}
bool SimulationDock::externalJobReady() const noexcept {
    return gpuReleased_.load();
}
void SimulationDock::endExternalJob() {
    externalJob_ = false;
    solve_->setEnabled(hasCurrentGeometry() && hasExcitation());
}
std::shared_ptr<const BemResult> SimulationDock::result() const noexcept {
    return result_;
}
std::uint64_t SimulationDock::resultRevision() const noexcept {
    return resultRevision_;
}

SimulationArchive SimulationDock::archive() const {
    return {
        {minimumFrequency_->value(),
         maximumFrequency_->value(),
         static_cast<std::size_t>(frequencyCount_->value()),
         {density_->value(), soundSpeed_->value()},
         velocity_->value(),
         phase_->value(),
         normalization_->currentIndex() == 0,
         beamwidthTarget_->isChecked(),
         targetH_->value(),
         targetV_->value(),
         targetTolerance_->value(),
         normalizationAngle_->value(),
         directivityScale_.minimum,
         directivityScale_.maximum,
         directivityScale_.tickStep,
         observationFrame_,
         fieldMode_->currentIndex() == 0,
         observationDistance_->value(),
         couplingEnabled_,
         couplingNetwork_,
         fastMultipole_->isChecked() ? BemAcceleration::fastMultipole : BemAcceleration::dense,
         fastMultipoleSettings_,
         static_cast<unsigned>(nearMaximumDepth_->value()),
         infiniteBaffle_,
         angularResponses_->settings(),
         acousticField_->settings(),
         phaseResponses_->settings(),
         usesAxisymmetricGeometry() ? BemGeometryMode::axisymmetric
                                    : BemGeometryMode::threeDimensional,
         {meridianLength_->value() * .001, meridianTolerance_->value() * .001,
          rotationalTolerance_->value() * 1e-6, static_cast<unsigned>(angularSectors_->value())},
         ctaSelection_},
        resultMesh_,
        resultRequest_,
        result_,
        resultRevision_,
        hasCurrentResult(),
        measurements_,
        resultAxisymmetricSampling_,
        result_ && resultCompletenessKnown_ ? std::optional{resultComplete_} : std::nullopt};
}

void SimulationDock::cancelPendingJobs() {
    jobStop_.request_stop();
    ++generation_;
    ++activeJob_;
    {
        std::lock_guard lock(workerMutex_);
        pendingJob_.reset();
    }
    busy_ = false;
    stop_->setEnabled(false);
    progress_->setValue(0);
    progress_->setVisible(false);
    updateResultControls();
}

void SimulationDock::validateArchive(const SimulationArchive &saved) const {
    if ((!saved.result && saved.resultComplete.has_value()) ||
        (saved.current && saved.resultComplete.has_value() && !*saved.resultComplete))
        throw std::invalid_argument("Ergebnisvollständigkeit passt nicht zum Archivzustand.");
    for (const auto &curve : saved.measurements)
        validateMeasuredResponse(curve);
    const auto &s = saved.settings;
    validateAxisymmetricSamplingSettings(s.axisymmetricSampling);
    if (s.geometryMode != BemGeometryMode::threeDimensional &&
        s.geometryMode != BemGeometryMode::axisymmetric)
        throw std::invalid_argument("Unbekannter Geometriemodus im Archiv.");
    validateAngularResponseSettings(s.angularResponses);
    validateFieldDisplaySettings(s.fieldDisplay);
    validatePhaseDisplaySettings(s.phaseDisplay);
    if (s.nearMaximumDepth > 10)
        throw std::invalid_argument(
            "Die maximale Nahquadraturtiefe muss zwischen 0 und 10 liegen.");
    if (s.observationFrame)
        normalizedAcousticFrame(*s.observationFrame);
    if (s.infiniteBaffle)
        validateInfiniteBaffleSettings(*s.infiniteBaffle);
    const auto &multipole = s.fastMultipole;
    if ((s.acceleration != BemAcceleration::dense &&
         s.acceleration != BemAcceleration::fastMultipole) ||
        multipole.interpolationOrder < 3 || multipole.interpolationOrder > 8 ||
        multipole.leafPanelCapacity < 4 || multipole.leafPanelCapacity > 256 ||
        !std::isfinite(multipole.separationRatio) || multipole.separationRatio < 1.5 ||
        multipole.separationRatio > 10 || multipole.maximumIterations < 1 ||
        multipole.maximumIterations > 10000 || multipole.restartLength < 2 ||
        multipole.restartLength > 300)
        throw std::invalid_argument("Ungültige Fast-Multipole-Einstellungen.");
    const auto valid = [](double value, QDoubleSpinBox *control) {
        if (!std::isfinite(value) || value < control->minimum() || value > control->maximum())
            return false;
        QDoubleSpinBox checked;
        checked.setDecimals(control->decimals());
        checked.setRange(control->minimum(), control->maximum());
        checked.setValue(value);
        return checked.value() == value;
    };
    if (!valid(s.minimumFrequency, minimumFrequency_) ||
        !valid(s.maximumFrequency, maximumFrequency_) || s.maximumFrequency < s.minimumFrequency ||
        s.frequencyCount < 1 ||
        s.frequencyCount > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        (s.frequencyCount > 1 && s.maximumFrequency == s.minimumFrequency) ||
        !valid(s.medium.density, density_) || !valid(s.medium.soundSpeed, soundSpeed_) ||
        !valid(s.velocityScale, velocity_) || !valid(s.phaseDegrees, phase_) ||
        !valid(s.targetHorizontalDegrees, targetH_) || !valid(s.targetVerticalDegrees, targetV_) ||
        !valid(s.targetToleranceDegrees, targetTolerance_) ||
        !valid(s.normalizationAngleDegrees, normalizationAngle_) ||
        !valid(s.observationDistance, observationDistance_) ||
        !valid(s.axisymmetricSampling.maximumLength * 1000, meridianLength_) ||
        !valid(s.axisymmetricSampling.chordTolerance * 1000, meridianTolerance_) ||
        !valid(s.axisymmetricSampling.rotationalTolerance * 1e6, rotationalTolerance_) ||
        !std::isfinite(s.directivityMinimumDb) || !std::isfinite(s.directivityMaximumDb) ||
        !std::isfinite(s.directivityStepDb) || s.directivityStepDb <= 0 ||
        s.directivityMaximumDb <= s.directivityMinimumDb ||
        std::ceil((s.directivityMaximumDb - s.directivityMinimumDb) / s.directivityStepDb) + 1 >
            1000 ||
        (static_cast<bool>(saved.result) != static_cast<bool>(saved.mesh)) ||
        (static_cast<bool>(saved.result) != saved.request.has_value()))
        throw std::invalid_argument("Ungültiger gespeicherter Simulationszustand.");
    if (saved.result) {
        bem_detail::planningWavenumber(saved.request->frequencies,
                                       saved.request->interactionPlanFrequencies,
                                       saved.request->medium.soundSpeed);
        const bool reduced = saved.request->axisymmetric.has_value();
        if (saved.axisymmetricSampling.has_value() != reduced ||
            (saved.result->geometryMode == BemGeometryMode::axisymmetric) != reduced ||
            (saved.result->geometryMode != BemGeometryMode::axisymmetric &&
             saved.result->geometryMode != BemGeometryMode::threeDimensional))
            throw std::invalid_argument(
                "Geometriemodus und gespeicherter 2D-Snapshot passen nicht zusammen.");
        if (reduced) {
            if (!saved.request->infiniteBaffle)
                throw std::invalid_argument(
                    "2D-Ergebnis benötigt den gespeicherten Luftkavitäts-/Schallwandvertrag.");
            validateAxisymmetricSamplingSettings(*saved.axisymmetricSampling);
            if (saved.axisymmetricSampling->sectors != saved.request->axisymmetric->sectors)
                throw std::invalid_argument(
                    "2D-Sektoren und Integrationsnetz passen nicht zusammen.");
        }
        if (saved.current && ((s.geometryMode == BemGeometryMode::axisymmetric) != reduced ||
                              (reduced && s.axisymmetricSampling != *saved.axisymmetricSampling)))
            throw std::invalid_argument(
                "Aktuelles 2D-Ergebnis gehört zu anderen Meridian-Einstellungen.");
        double previous = 0;
        for (const auto &frequency : saved.result->frequencies) {
            if (!std::isfinite(frequency.frequency) || frequency.frequency <= previous)
                throw std::invalid_argument(
                    "Gespeicherte Ergebnisfrequenzen müssen positiv und streng aufsteigend sein.");
            previous = frequency.frequency;
        }
    } else if (saved.axisymmetricSampling)
        throw std::invalid_argument("2D-Netzsnapshot ohne Ergebnis.");
}

void SimulationDock::restoreArchive(const SimulationArchive &saved) {
    validateArchive(saved);
    setRealtimeEnabled(false);
    const auto &s = saved.settings;
    cancelPendingJobs();
    {
        const QSignalBlocker mode(geometryMode_), length(meridianLength_),
            chord(meridianTolerance_), symmetry(rotationalTolerance_), sectors(angularSectors_);
        geometryMode_->setCurrentIndex(s.geometryMode == BemGeometryMode::axisymmetric ? 1 : 0);
        meridianLength_->setValue(s.axisymmetricSampling.maximumLength * 1000);
        meridianTolerance_->setValue(s.axisymmetricSampling.chordTolerance * 1000);
        rotationalTolerance_->setValue(s.axisymmetricSampling.rotationalTolerance * 1e6);
        angularSectors_->setValue(static_cast<int>(s.axisymmetricSampling.sectors));
    }
    observationFrame_ = s.observationFrame;
    if (observationFrame_)
        observationFrame_ = normalizedAcousticFrame(*observationFrame_);
    if (onObservationFrameChanged)
        onObservationFrameChanged();
    const QSignalBlocker a(minimumFrequency_), b(maximumFrequency_), c(frequencyCount_),
        d(density_), e(soundSpeed_), f(velocity_), g(phase_), h(normalization_),
        i(beamwidthTarget_), j(targetH_), k(targetV_), l(targetTolerance_), m(normalizationAngle_),
        n(directivityMinimum_), o(directivityMaximum_), p(directivityStep_), q(fieldMode_),
        r(observationDistance_), multipoleBlock(fastMultipole_), orderBlock(multipoleOrder_),
        nearDepthBlock(nearMaximumDepth_);
    minimumFrequency_->setValue(s.minimumFrequency);
    maximumFrequency_->setValue(s.maximumFrequency);
    frequencyCount_->setValue(static_cast<int>(s.frequencyCount));
    maximumFrequency_->setEnabled(s.frequencyCount > 1);
    density_->setValue(s.medium.density);
    soundSpeed_->setValue(s.medium.soundSpeed);
    velocity_->setValue(s.velocityScale);
    phase_->setValue(s.phaseDegrees);
    fieldMode_->setCurrentIndex(s.farField ? 0 : 1);
    observationDistance_->setValue(s.observationDistance);
    fastMultipoleSettings_ = s.fastMultipole;
    fastMultipole_->setChecked(s.acceleration == BemAcceleration::fastMultipole);
    multipoleOrder_->setValue(static_cast<int>(s.fastMultipole.interpolationOrder));
    multipoleOrder_->setEnabled(fastMultipole_->isChecked());
    nearMaximumDepth_->setValue(static_cast<int>(s.nearMaximumDepth));
    couplingNetwork_ = s.couplingNetwork;
    couplingEnabled_ = s.couplingEnabled;
    infiniteBaffle_ = s.infiniteBaffle;
    angularResponses_->setSettings(s.angularResponses);
    acousticField_->setSettings(s.fieldDisplay);
    phaseResponses_->setSettings(s.phaseDisplay);
    updateDomainControls();
    normalization_->setCurrentIndex(s.normalizeToOnAxis ? 0 : 1);
    normalizationAngle_->setValue(s.normalizationAngleDegrees);
    normalizationAngle_->setEnabled(s.normalizeToOnAxis);
    applyDirectivityScale({s.directivityMinimumDb, s.directivityMaximumDb, s.directivityStepDb});
    beamwidthTarget_->setChecked(s.beamwidthTargetEnabled);
    targetH_->setValue(s.targetHorizontalDegrees);
    targetV_->setValue(s.targetVerticalDegrees);
    targetTolerance_->setValue(s.targetToleranceDegrees);
    for (auto *control : {targetH_, targetV_, targetTolerance_})
        control->setEnabled(s.beamwidthTargetEnabled);
    ctaSelection_ = s.ctaExcitationIndex;
    result_ = saved.result;
    measurements_ = saved.measurements;
    displayFrequencies_.clear();
    displayBaseline_.clear();
    resultMesh_ = saved.mesh;
    resultRequest_ = saved.request;
    resultAxisymmetricSampling_ = saved.axisymmetricSampling;
    resultRevision_ = saved.revision;
    resultCompletenessKnown_ = saved.result && (saved.resultComplete.has_value() || saved.current);
    resultComplete_ = saved.result && saved.resultComplete.value_or(saved.current);
    bool current = saved.current && hasCurrentGeometry() && result_ && resultRequest_ &&
                   resultRevision_ == currentRevision_;
    if (current) {
        try {
            current = sameRequest(*resultRequest_, makeRequest());
        } catch (const std::exception &) {
            current = false;
        }
    }
    resultGeneration_ = current ? generation_ : 0;
    error_.clear();
    presentResult(current);
    workspaceTabs_->setCurrentIndex(result_ || !measurements_.empty() ? 1 : 0);
    solve_->setEnabled(hasCurrentGeometry() && hasExcitation());
    setStatus(
        result_
            ? (current ? tr("Gespeicherte Ergebnisse") : tr("Gespeicherte Ergebnisse · veraltet"))
            : (hasCurrentGeometry() ? tr("Bereit") : tr("Geometrie / Netz erforderlich")));
}

void SimulationDock::updateDisplayEstimate() {
    std::map<double, DisplayFrequency> computed;
    const auto missing = std::numeric_limits<double>::quiet_NaN();
    if (result_ && resultRequest_) {
        const auto angularSamples = angularSampleCount(resultRequest_);
        for (const auto &frequency : result_->frequencies) {
            if (frequency.excitations.empty())
                continue;
            const auto &excitation = frequency.excitations.front();
            auto &display = computed[frequency.frequency];
            display.onAxisDb =
                excitation.observerPressure.empty()
                    ? missing
                    : decibels(std::abs(excitation.observerPressure[0]), pressureReference);
            for (const auto &z : excitation.impedances)
                display.impedances[z.faceId] = normalizedVolumeImpedance(
                    z.volumeImpedance, z.area, resultRequest_->medium.density,
                    resultRequest_->medium.soundSpeed);
            display.horizontalDb.assign(angularCount, missing);
            display.verticalDb.assign(angularCount, missing);
            if (angularSamples && excitation.observerPressure.size() == 1 + 2 * angularSamples) {
                const auto offset = (angularCount - angularSamples) / 2;
                for (std::size_t i = 0; i < angularSamples; ++i) {
                    display.horizontalDb[offset + i] =
                        decibels(std::abs(excitation.observerPressure[1 + i]), pressureReference);
                    display.verticalDb[offset + i] =
                        decibels(std::abs(excitation.observerPressure[1 + angularSamples + i]),
                                 pressureReference);
                }
            }
        }
    }
    if (resultComplete_ || !resultCompletenessKnown_ || !result_) {
        displayFrequencies_ = std::move(computed);
        displayBaseline_.clear();
        return;
    }
    displayFrequencies_ = displayBaseline_;
    for (const auto &[frequency, sample] : computed)
        displayFrequencies_[frequency] = sample;
    const auto drive = [&](auto get, auto set) {
        live::Trace baseline, fresh;
        for (const auto &[frequency, sample] : displayBaseline_)
            baseline[frequency] = get(sample);
        for (const auto &[frequency, sample] : computed)
            fresh[frequency] = get(sample);
        for (const auto &[frequency, value] : live::estimate(baseline, fresh))
            set(displayFrequencies_.at(frequency), value);
    };
    drive([](const auto &s) { return s.onAxisDb; },
          [](auto &s, double value) { s.onAxisDb = value; });
    std::set<std::uint32_t> faces;
    for (const auto *samples : {&displayBaseline_, &computed})
        for (const auto &[frequency, sample] : *samples)
            for (const auto &[face, impedance] : sample.impedances)
                faces.insert(face);
    for (const auto face : faces)
        for (bool imaginary : {false, true})
            drive(
                [&](const auto &s) {
                    const auto found = s.impedances.find(face);
                    return found == s.impedances.end() ? missing
                           : imaginary                 ? found->second.imag()
                                                       : found->second.real();
                },
                [&](auto &s, double value) {
                    if (imaginary)
                        s.impedances[face].imag(value);
                    else
                        s.impedances[face].real(value);
                });
    for (std::size_t i = 0; i < angularCount; ++i) {
        drive([&](const auto &s) { return s.horizontalDb[i]; },
              [&](auto &s, double value) { s.horizontalDb[i] = value; });
        drive([&](const auto &s) { return s.verticalDb[i]; },
              [&](auto &s, double value) { s.verticalDb[i] = value; });
    }
}

void SimulationDock::presentResult(bool current) {
    updateResultQuality();
    updateFieldAndPhase(current);
    updateCta2034(current);
    updateAngularFrequencyPlots(current);
    const bool farField = !resultRequest_ || resultRequest_->farField.has_value();
    const double reference =
        resultRequest_ && resultRequest_->farField
            ? resultRequest_->farField->referenceDistance
            : (resultRequest_ && resultRequest_->observers.size() == 147
                   ? .5 *
                         std::hypot(resultRequest_->observers[0].x - resultRequest_->observers[1].x,
                                    resultRequest_->observers[0].y - resultRequest_->observers[1].y,
                                    resultRequest_->observers[0].z - resultRequest_->observers[1].z)
                   : observationDistance_->value());
    const QString location = farField ? tr("Fernfeld · Bezug %1 m").arg(reference)
                                      : tr("Messabstand %1 m").arg(reference);
    onAxis_->setTitle(tr("On-Axis · %1").arg(location));
    isobarH_->setTitle(tr("Isobaren H · %1").arg(location));
    radialIsobarH_->setTitle(tr("Isobaren rund H · %1").arg(location));
    isobarV_->setTitle(tr("Isobaren V · %1").arg(location));
    radialIsobarV_->setTitle(tr("Isobaren rund V · %1").arg(location));
    const auto selectedFrequency = frequencySelector_->currentData();
    const auto selectedFace = impedanceSelector_->currentData();
    const QSignalBlocker frequencies(frequencySelector_), impedances(impedanceSelector_);
    frequencySelector_->clear();
    impedanceSelector_->clear();
    updateDisplayEstimate();
    if (displayFrequencies_.empty()) {
        setWindowTitle(tr("Akustik"));
        onAxis_->setSeries({});
        impedance_->setSeries({});
        polarH_->setLevels({});
        polarV_->setLevels({});
        isobarH_->setLevels({}, {}, {});
        radialIsobarH_->setLevels({}, {}, {});
        isobarV_->setLevels({}, {}, {});
        radialIsobarV_->setLevels({}, {}, {});
        updateFrequencyPlots();
        updateAngularPlots();
        return;
    }
    const auto completeness = !resultCompletenessKnown_ ? tr(" · Vollständigkeit unbekannt")
                              : !resultComplete_ ? tr(" · Teilraster") : QString{};
    setWindowTitle(tr("Akustik · Rev. %1%2%3")
                       .arg(resultRevision_)
                       .arg(current ? QString{} : tr(" · veraltet"))
                       .arg(completeness));
    for (const auto &[frequency, display] : displayFrequencies_)
        frequencySelector_->addItem(frequencyLabel(frequency) + " Hz", frequency);
    const int frequencyIndex =
        selectedFrequency.isValid() ? frequencySelector_->findData(selectedFrequency) : -1;
    frequencySelector_->setCurrentIndex(
        frequencyIndex >= 0 ? frequencyIndex : static_cast<int>(displayFrequencies_.size() / 2));
    std::set<std::uint32_t> faces;
    for (const auto &[frequency, display] : displayFrequencies_)
        for (const auto &[face, impedance] : display.impedances)
            faces.insert(face);
    for (const auto face : faces)
        impedanceSelector_->addItem(tr("Fläche %1").arg(face), face);
    const int faceIndex = selectedFace.isValid() ? impedanceSelector_->findData(selectedFace) : -1;
    if (faceIndex >= 0)
        impedanceSelector_->setCurrentIndex(faceIndex);
    updateFrequencyPlots();
    updateAngularPlots();
    updateHeatmaps();
}
void SimulationDock::updateCta2034(bool current) {
    PlotSeries axis{tr("On-Axis"), QColor(28, 102, 180), {}},
        listening{tr("Listening Window"), QColor(30, 145, 85), {}},
        reflections{tr("Early Reflections"), QColor(205, 112, 33), {}},
        power{tr("Sound Power · Näherung"), QColor(175, 55, 55), {}},
        predicted{tr("Raumprognose"), QColor(60, 120, 120), {}},
        di{tr("Early Reflections DI"), QColor(125, 65, 160), {}},
        powerDi{tr("Sound Power DI"), QColor(175, 55, 55), {}};
    const auto count = angularSampleCount(resultRequest_);
    {
        const QSignalBlocker blocker(ctaExcitation_);
        ctaExcitation_->clear();
        if (result_ && !result_->frequencies.empty())
            for (const auto &e : result_->frequencies.front().excitations)
                ctaExcitation_->addItem(QString::fromStdString(e.name));
        if (ctaExcitation_->count()) {
            const int selected = ctaSelection_ >= 0 && ctaSelection_ < ctaExcitation_->count()
                                     ? ctaSelection_ : 0;
            ctaExcitation_->setCurrentIndex(selected);
            ctaSelection_ = selected;
        }
    }
    if (result_ && count && resultRequest_ && resultRequest_->farField &&
        ctaExcitation_->currentIndex() >= 0) {
        for (const auto &row : ctaResults(*resultRequest_, *result_, static_cast<std::size_t>(ctaExcitation_->currentIndex()))) {
            const auto &curves = row.curves;
            const auto append = [&](PlotSeries &series, std::optional<double> value) {
                series.points.emplace_back(row.frequency,
                    value.value_or(std::numeric_limits<double>::quiet_NaN()));
            };
            append(axis, curves.onAxisDb);
            append(listening, curves.listeningWindowDb);
            append(reflections, curves.earlyReflectionsDb);
            append(di, curves.earlyReflectionsDiDb);
            append(power, curves.soundPowerDb);
            append(powerDi, curves.soundPowerDiDb);
            append(predicted, curves.predictedInRoomDb);
        }
    }
    ctaLevels_->setSeries({std::move(axis), std::move(listening), std::move(reflections),
                           std::move(power), std::move(predicted)});
    ctaDi_->setSeries({std::move(di), std::move(powerDi)});
    ctaStatus_->setText(result_ && resultRequest_ && resultRequest_->farField && count
        ? tr("Anregung: %5 · Simulation · Fernfeld %1 m · Rev. %2%3%4. CTA-Schallleistung ist eine "
             "70-Richtungs-Näherung; die Raumkurve ist eine Schätzung. Fehlende Winkel ergeben Lücken.")
              .arg(resultRequest_->farField->referenceDistance).arg(resultRevision_)
              .arg(current ? QString{} : tr(" · veraltet"))
              .arg(!resultCompletenessKnown_ ? tr(" · Vollständigkeit unbekannt")
                                            : resultComplete_ ? QString{} : tr(" · Teilraster"))
              .arg(ctaExcitation_->currentText())
        : tr("CTA-Auswertung benötigt ein Fernfeldergebnis mit H/V-Winkelstützstellen. "
             "Für Schallleistung und Raumprognose sind vollständige 360°-Schnitte erforderlich."));
}
void SimulationDock::updateFieldAndPhase(bool current) {
    const auto count = angularSampleCount(resultRequest_);
    std::vector<PhaseResponseRow> rows;
    std::vector<double> angles;
    auto frame = observationFrame();
    if (result_ && count && resultRequest_) {
        angles = anglesDegrees(count);
        const auto &points = resultRequest_->observers;
        const auto center = count == 73 ? scaled(added(points[0], points[1]), .5)
                                        : scaled(added(points[1], points[count]), .5);
        const auto horizontalIndex = count == 73 ? 54u : 36u;
        frame = normalizedAcousticFrame({center, added(points[0], scaled(center, -1)),
                                         added(points[1 + horizontalIndex], scaled(center, -1))});
        for (const auto &frequency : result_->frequencies) {
            PhaseResponseRow row{frequency.frequency, {}, {}};
            if (!frequency.excitations.empty()) {
                const auto &pressure = frequency.excitations.front().observerPressure;
                if (pressure.size() == 1 + 2 * count)
                    for (std::size_t i = 0; i < count; ++i) {
                        row.horizontal.push_back(angles[i] == 0 ? pressure[0] : pressure[1 + i]);
                        row.vertical.push_back(angles[i] == 0 ? pressure[0]
                                                              : pressure[1 + count + i]);
                    }
            }
            rows.push_back(std::move(row));
        }
    }
    phaseResponses_->setRows(std::move(angles), std::move(rows),
                             tr("Komplexe Ergebnisdaten · Rev. %1%2 · exp(+iωt)")
                                 .arg(resultRevision_)
                                 .arg(current ? QString{} : tr(" · veraltet")));
    acousticField_->setSnapshot(resultMesh_, resultRequest_, result_, frame, resultRevision_,
                                current);
}
void SimulationDock::updateAngularFrequencyPlots(bool current) {
    std::vector<AngularResponseRow> rows;
    const auto count = angularSampleCount(resultRequest_);
    if (result_ && count) {
        const auto angles = anglesDegrees(count);
        for (const auto &frequency : result_->frequencies) {
            AngularResponseRow row{frequency.frequency, {}, {}};
            if (!frequency.excitations.empty()) {
                const auto &pressure = frequency.excitations.front().observerPressure;
                if (pressure.size() == 1 + 2 * count) {
                    for (std::size_t i = 0; i < count; ++i) {
                        // The dedicated on-axis sample defines zero in both planes.
                        const auto h = angles[i] == 0 ? pressure[0] : pressure[1 + i];
                        const auto v = angles[i] == 0 ? pressure[0] : pressure[1 + count + i];
                        row.horizontal.push_back({angles[i], soundPressureLevelDecibels(h)});
                        row.vertical.push_back({angles[i], soundPressureLevelDecibels(v)});
                    }
                }
            }
            rows.push_back(std::move(row));
        }
    }
    QString context;
    if (result_ && resultRequest_) {
        context = tr("Rev. %1").arg(resultRevision_);
        if (!current)
            context += tr(" · veraltet");
        const bool completeRaster = std::equal(
            result_->frequencies.begin(), result_->frequencies.end(),
            resultRequest_->frequencies.begin(), resultRequest_->frequencies.end(),
            [](const auto &row, double frequency) { return row.frequency == frequency; });
        if (!resultCompletenessKnown_)
            context += tr(" · Vollständigkeit unbekannt");
        else if (!resultComplete_ || !completeRaster)
            context += tr(" · Teilraster");
        if (resultRequest_->farField)
            context += tr(" · Fernfeld %1 m").arg(resultRequest_->farField->referenceDistance);
        else
            context += tr(" · Messabstand");
    }
    angularResponses_->setRows(std::move(rows), context);
}
void SimulationDock::setStatus(QString text, bool error) {
    progress_->setVisible(busy_);
    status_->setText(text);
    status_->setToolTip(text);
    status_->setStyleSheet(error ? "color: #a72e2b" : "");
    updateResultControls();
}

BemRequest SimulationDock::makeRequest() const {
    if (!hasCurrentGeometry() || !hasExcitation())
        throw std::invalid_argument("Aktuelles Netz mit Antriebsfläche erforderlich.");
    const double minimum = minimumFrequency_->value(), maximum = maximumFrequency_->value();
    const int count = frequencyCount_->value();
    if (count > 1 && maximum <= minimum)
        throw std::invalid_argument("Obere Frequenz muss über der unteren Frequenz liegen.");
    requireResultMemory(!usesAxisymmetricGeometry() && mesh_ ? mesh_->triangles.size() : 0,
                        !usesAxisymmetricGeometry() && mesh_ ? mesh_->vertices.size() : 0,
                        infiniteBaffle_ ? infiniteBaffle_->exteriorDrivenMesh.triangles.size() : 0,
                        velocities_.size(), static_cast<std::size_t>(count),
                        couplingEnabled_ && infiniteBaffle_ && !velocities_.empty() ? 2 : 1,
                        realtimeEnabled_);
    const auto settings = archive().settings;
    auto built = makeSimulationRequest(settings, velocities_, origin_, axis_);
    if (usesAxisymmetricGeometry() && resultRequest_ && resultRequest_->axisymmetric &&
        resultAxisymmetricSampling_ == settings.axisymmetricSampling &&
        resultRevision_ == currentRevision_) {
        built.axisymmetric = resultRequest_->axisymmetric;
        built.axisymmetric->surfaceContract = AxisymmetricSurfaceContract::continuousRevolution;
    }
    return built;
}

void SimulationDock::requestSolve() {
    if (busy_ || externalJob_)
        return;
    std::shared_ptr<SolveJob> job;
    std::map<double, DisplayFrequency> baseline;
    std::stop_source newStop{std::nostopstate};
    try {
        auto request = makeRequest();
        if (resultRequest_ && sameDisplayDirections(*resultRequest_, request))
            baseline = displayFrequencies_;
        newStop = std::stop_source{};
        job = std::make_shared<SolveJob>();
        job->mesh = mesh_;
        job->boundaryDefinitions = velocities_;
        std::sort(job->boundaryDefinitions.begin(), job->boundaryDefinitions.end(),
                  [](const FaceVelocity &a, const FaceVelocity &b) { return a.faceId < b.faceId; });
        job->request = std::move(request);
        if (usesAxisymmetricGeometry()) {
            job->source = geometrySource_;
            job->axisymmetricSampling = archive().settings.axisymmetricSampling;
        }
        job->revision = currentRevision_;
        job->generation = generation_;
        job->id = activeJob_ + 1;
        job->realtime = realtimeEnabled_;
        job->stop = newStop.get_token();
    } catch (const std::exception &e) {
        error_ = QString::fromUtf8(e.what());
        setStatus(error_, true);
        setRealtimeEnabled(false);
        return;
    }
    displayBaseline_ = std::move(baseline);
    jobStop_ = std::move(newStop);
    activeJob_ = job->id;
    busy_ = true;
    realtimePending_ = false;
    realtimeTimer_->stop();
    error_.clear();
    solve_->setEnabled(false);
    stop_->setEnabled(true);
    progress_->setRange(0, static_cast<int>(job->request.frequencies.size()));
    progress_->setValue(0);
    setStatus(usesAxisymmetricGeometry() ? tr("2D · Meridianprüfung und Integrationsnetz …")
                                         : tr("CUDA · Operatoraufbau …"));
    revealResultView();
    {
        std::lock_guard lock(workerMutex_);
        pendingJob_ = std::move(job);
    }
    workerWake_.notify_one();
}

void SimulationDock::runWorker(std::stop_token shutdown) {
    // CUDA allocation, cache replacement and destruction all stay on this worker.
    std::shared_ptr<const SurfaceMesh> cachedMesh;
    std::vector<FaceVelocity> cachedDefinitions;
    std::unique_ptr<CudaBemSolver> solver;
    std::shared_ptr<const GeometrySource> cachedSource;
    std::optional<AxisymmetricSamplingSettings> cachedSampling;
    std::shared_ptr<const SurfaceMesh> reducedMesh;
    std::optional<AxisymmetricBem> reducedLayout;
    while (!shutdown.stop_requested()) {
        std::shared_ptr<const SolveJob> job;
        {
            std::unique_lock lock(workerMutex_);
            if (!workerWake_.wait(lock, shutdown,
                                  [this] { return pendingJob_ || releaseGpuRequested_; }))
                break;
            if (releaseGpuRequested_) {
                releaseGpuRequested_ = false;
                lock.unlock();
                solver.reset();
                cachedMesh.reset();
                cachedDefinitions.clear();
                gpuReleased_.store(true);
                continue;
            }
            job = std::exchange(pendingJob_, {});
            gpuReleased_.store(false);
        }
        if (shutdown.stop_requested())
            break;
        const auto revision = job->revision;
        const auto generation = job->generation;
        const auto id = job->id;
        const auto stop = job->stop;
        std::shared_ptr<const BemResult> result;
        QString error;
        bool reusedGeometry = false;
        auto actualMesh = job->mesh;
        BemRequest actualRequest;
        const auto start = std::chrono::steady_clock::now();
        try {
            actualRequest = job->request;
            if (job->axisymmetricSampling) {
                if (!job->source)
                    throw std::invalid_argument(
                        "2D-Rotation benötigt einen aktuellen exakten Geometriesnapshot.");
                if (cachedSource != job->source || cachedSampling != job->axisymmetricSampling ||
                    !reducedMesh) {
                    const auto &s = *job->axisymmetricSampling;
                    MeridianSampling sampling;
                    sampling.maximumLength = s.maximumLength;
                    sampling.chordTolerance = s.chordTolerance;
                    sampling.rotationalTolerance = s.rotationalTolerance;
                    sampling.faceTargetLengths = job->source->faceTargets;
                    AxisymmetricMeridian meridian;
                    if (job->source->cad) {
                        const auto &cad = *job->source->cad;
                        if (!cad.geometry)
                            throw std::invalid_argument("2D-Rotation: Exakte CAD-Geometrie fehlt.");
                        for (const auto &face : cad.faces)
                            sampling.faceConditions.emplace(face.id, face.boundary);
                        meridian = makeAxisymmetricMeridian(*cad.geometry, cad.acousticOrigin,
                                                            cad.acousticAxis, sampling, stop);
                    } else if (job->source->draft) {
                        meridian = makeAxisymmetricMeridian(*job->source->draft, sampling, stop);
                    } else
                        throw std::invalid_argument("2D-Rotation: Geometriequelle fehlt.");
                    auto generated = buildAxisymmetricMesh(meridian, s.sectors, stop);
                    if (!job->source->reversedFaces.empty()) {
                        generated.mesh =
                            applyFaceOrientation(generated.mesh, job->source->reversedFaces);
                        if (!inspectMesh(generated.mesh).watertight())
                            throw std::invalid_argument(
                                "2D-Rotation: Normalenumkehr ergibt keine geschlossene, konsistent "
                                "orientierte Fläche.");
                        galerkin::validateExteriorComponents(generated.mesh);
                    }
                    reducedLayout = AxisymmetricBem{meridian.origin,
                                                    meridian.axis,
                                                    s.sectors,
                                                    std::move(generated.pressureDofs),
                                                    std::move(generated.fluxDofs),
                                                    std::move(generated.testPanels)};
                    reducedMesh = std::make_shared<SurfaceMesh>(std::move(generated.mesh));
                    cachedSource = job->source;
                    cachedSampling = job->axisymmetricSampling;
                }
                actualMesh = reducedMesh;
                actualRequest.axisymmetric = reducedLayout;
                requireResultMemory(actualMesh->triangles.size(), actualMesh->vertices.size(),
                                    actualRequest.infiniteBaffle
                                        ? actualRequest.infiniteBaffle->exteriorDrivenMesh.triangles.size()
                                        : 0,
                                    job->boundaryDefinitions.size(),
                                    actualRequest.frequencies.size(),
                                    actualRequest.excitations.size() + (actualRequest.coupling ? 1 : 0),
                                    job->realtime);
            }
            reusedGeometry =
                solver && cachedMesh == actualMesh &&
                std::equal(cachedDefinitions.begin(), cachedDefinitions.end(),
                           job->boundaryDefinitions.begin(), job->boundaryDefinitions.end(),
                           [](const FaceVelocity &a, const FaceVelocity &b) {
                               return a.faceId == b.faceId && a.normalVelocity == b.normalVelocity;
                           });
            if (!reusedGeometry) {
                solver.reset();
                cachedMesh = actualMesh;
                cachedDefinitions = job->boundaryDefinitions;
                if (!stop.stop_requested())
                    solver = std::make_unique<CudaBemSolver>(*cachedMesh);
            }
            if (stop.stop_requested())
                throw std::runtime_error("Berechnung abgebrochen");
            if (job->realtime) {
                BemResult accumulated;
                auto singleRequest = actualRequest;
                for (const auto index : frequencySchedule(actualRequest.frequencies.size())) {
                    const auto frequency = actualRequest.frequencies[index];
                    if (stop.stop_requested())
                        throw std::runtime_error("Berechnung abgebrochen");
                    singleRequest.frequencies = {frequency};
                    auto single = solver->solve(singleRequest, stop);
                    if (stop.stop_requested())
                        throw std::runtime_error("Berechnung abgebrochen");
                    accumulated.device = single.device;
                    accumulated.acceleration = single.acceleration;
                    accumulated.geometryMode = single.geometryMode;
                    accumulated.triangleCount = single.triangleCount;
                    accumulated.pressureDofCount = single.pressureDofCount;
                    accumulated.peakDeviceBytes =
                        std::max(accumulated.peakDeviceBytes, single.peakDeviceBytes);
                    const auto position = std::lower_bound(
                        accumulated.frequencies.begin(), accumulated.frequencies.end(), frequency,
                        [](const auto &entry, double f) { return entry.frequency < f; });
                    accumulated.frequencies.insert(position, std::move(single.frequencies.at(0)));
                    // Every queued GUI publication owns a completed immutable frequency set.
                    std::shared_ptr<const BemResult> partial =
                        std::make_shared<BemResult>(accumulated);
                    auto partialRequest = actualRequest;
                    partialRequest.frequencies.clear();
                    for (const auto &computed : accumulated.frequencies)
                        partialRequest.frequencies.push_back(computed.frequency);
                    if (job->partial.offer({std::move(partial), std::move(partialRequest), frequency})) {
                        const bool queued = QMetaObject::invokeMethod(
                            this,
                            [this, job, actualMesh] {
                                auto latest = job->partial.take();
                                if (!latest || job->id != activeJob_ || job->generation != generation_ ||
                                    job->revision != currentRevision_ || job->stop.stop_requested())
                                    return;
                                result_ = std::move(latest->result);
                                resultMesh_ = actualMesh;
                                resultRequest_ = std::move(latest->request);
                                resultAxisymmetricSampling_ = job->axisymmetricSampling;
                                resultRevision_ = job->revision;
                                resultGeneration_ = job->generation;
                                resultComplete_ = false;
                                resultCompletenessKnown_ = true;
                                presentResult(true);
                                setWindowTitle(
                                    tr("Akustik · Rev. %1 · Teilraster %2/%3 · Zwischenwerte geschätzt")
                                        .arg(resultRevision_)
                                        .arg(result_->frequencies.size())
                                        .arg(job->request.frequencies.size()));
                                progress_->setValue(static_cast<int>(result_->frequencies.size()));
                                setStatus(tr("CUDA · %1 / %2 Frequenzen · zuletzt %3 Hz")
                                              .arg(result_->frequencies.size())
                                              .arg(job->request.frequencies.size())
                                              .arg(frequencyLabel(latest->frequency)));
                                if (onChanged)
                                    onChanged();
                            },
                            Qt::QueuedConnection);
                        if (!queued) {
                            (void)job->partial.take();
                            throw std::runtime_error("Teilstand konnte nicht an die Oberfläche übergeben werden.");
                        }
                    }
                }
                result = std::make_shared<BemResult>(std::move(accumulated));
            } else {
                result = std::make_shared<BemResult>(solver->solve(
                    actualRequest, stop, [this, job](const BemProgress &progress) {
                        if (!job->progress.offer(progress)) return;
                        const bool queued = QMetaObject::invokeMethod(
                            this,
                            [this, job] {
                                auto latest = job->progress.take();
                                if (!latest || job->generation != generation_ || job->id != activeJob_)
                                    return;
                                progress_->setValue(static_cast<int>(latest->completedFrequencies));
                                setStatus(tr("CUDA · %1 / %2 · %3 Hz")
                                              .arg(latest->completedFrequencies)
                                              .arg(latest->totalFrequencies)
                                              .arg(frequencyLabel(latest->frequency)));
                            }, Qt::QueuedConnection);
                        if (!queued) {
                            (void)job->progress.take();
                            throw std::runtime_error("Fortschritt konnte nicht an die Oberfläche übergeben werden.");
                        }
                    }));
            }
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        }
        const double elapsedMilliseconds =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                .count();
        QMetaObject::invokeMethod(
            this,
            [this, job, actualMesh, actualRequest = std::move(actualRequest), id, revision,
             generation, result = std::move(result), error = std::move(error), elapsedMilliseconds,
             reusedGeometry] {
                if (id != activeJob_)
                    return;
                busy_ = false;
                stop_->setEnabled(false);
                solve_->setEnabled(hasCurrentGeometry() && hasExcitation());
                if (generation != generation_ || revision != currentRevision_) {
                    setStatus(tr("Eingabe geändert · bisherige Anzeige bleibt"));
                    scheduleRealtime();
                    return;
                }
                if (job->stop.stop_requested()) {
                    error_ = tr("Berechnung abgebrochen");
                    setStatus(error_);
                    scheduleRealtime();
                    return;
                }
                if (!error.isEmpty()) {
                    error_ = error;
                    setStatus(error_, true);
                    setRealtimeEnabled(false);
                    return;
                }
                result_ = result;
                resultMesh_ = actualMesh;
                resultRequest_ = actualRequest;
                resultAxisymmetricSampling_ = job->axisymmetricSampling;
                resultRevision_ = revision;
                resultGeneration_ = generation;
                resultComplete_ = true;
                resultCompletenessKnown_ = true;
                realtimePending_ = false;
                realtimeTimer_->stop();
                presentResult(true);
                if (onChanged)
                    onChanged();
                double residual = 0, stageMilliseconds = 0, auditMilliseconds = 0,
                       planMilliseconds = 0;
                unsigned maximumIterations = 0, planHits = 0;
                for (const auto &f : result_->frequencies) {
                    stageMilliseconds +=
                        f.assemblyMilliseconds + f.solveMilliseconds + f.fieldMilliseconds;
                    auditMilliseconds += f.residualCheckMilliseconds;
                    planMilliseconds += f.interactionPlanMilliseconds;
                    planHits += f.interactionPlanReused ? 1 : 0;
                    for (const auto &e : f.excitations) {
                        residual = std::max(residual, e.boundaryResidual);
                        maximumIterations = std::max(maximumIterations, e.iterations);
                    }
                }
                const QString method =
                    result_->geometryMode == BemGeometryMode::axisymmetric ? tr("2D-Rotation (m=0)")
                    : result_->acceleration == BemAcceleration::fastMultipole ? tr("Fast Multipole")
                                                                              : tr("Dicht");
                setStatus(tr("%1 · %2 Frequenzen · %3 s · Residuum %4")
                              .arg(method)
                              .arg(result_->frequencies.size())
                              .arg(elapsedMilliseconds * .001, 0, 'f', 2)
                              .arg(residual, 0, 'g', 3));
                status_->setToolTip(
                    tr("Gesamtzeit: %1 ms\nNumerische Phasen: %2 ms\nGeometriedaten: %3\n"
                       "GPU: %4\nDirekte Residuumprüfung: %5 ms\nMaximale Iterationen: %6\n"
                       "CPU-Interaktionsplanung: %7 ms\nHost-Plan-Cache: %8 Treffer")
                        .arg(elapsedMilliseconds, 0, 'f', 1)
                        .arg(stageMilliseconds, 0, 'f', 1)
                        .arg(reusedGeometry ? tr("wiederverwendet") : tr("neu vorbereitet"))
                        .arg(QString::fromStdString(result_->device))
                        .arg(auditMilliseconds, 0, 'f', 1)
                        .arg(maximumIterations)
                        .arg(planMilliseconds, 0, 'f', 1)
                        .arg(planHits) +
                    (result_->geometryMode == BemGeometryMode::axisymmetric
                         ? tr("\nResiduum des reduzierten m=0-Burton–Miller-Systems.")
                         : QString{}));
                std::fprintf(stdout, "CUDA job completed: elapsed_ms=%.3f geometry_cache=%s\n",
                             elapsedMilliseconds, reusedGeometry ? "reused" : "new");
                std::fflush(stdout);
            },
            Qt::QueuedConnection);
    }
}

void SimulationDock::importMeasurement() {
    const auto path =
        QFileDialog::getOpenFileName(this, tr("Messkurve importieren"), {},
                                     tr("Frequenzgänge (*.csv *.txt *.frd);;Alle Dateien (*)"));
    if (path.isEmpty())
        return;
    try {
        auto response = readMeasuredResponse(path);
        bool ok = false;
        const auto reference =
            QInputDialog::getItem(this, tr("Messbezug"), tr("Pegel"),
                                  {tr("Absolut · dB SPL"), tr("Relativ · dB")}, 0, false, &ok);
        if (!ok)
            return;
        response.absoluteSpl = reference == tr("Absolut · dB SPL");
        if (!response.absoluteSpl) {
            response.offsetDb = QInputDialog::getDouble(
                this, tr("Messbezug"), tr("Pegelversatz zur SPL-Achse [dB]"), 0, -300, 300, 2, &ok);
            if (!ok)
                return;
        }
        measurements_.push_back(std::move(response));
        updateFrequencyPlots();
        resultTabs_->setCurrentIndex(0);
        revealResultView();
        if (onChanged)
            onChanged();
    } catch (const std::exception &e) {
        QMessageBox::warning(this, tr("Messkurve"), QString::fromUtf8(e.what()));
    }
}

void SimulationDock::exportResults() {
    if (!result_ || !resultRequest_) {
        setStatus(tr("Keine berechneten Ergebnisse vorhanden."), true);
        return;
    }
    const auto path = QFileDialog::getSaveFileName(this, tr("Ergebnisse exportieren"),
                                                   "HornSculpt-Ergebnisse.csv", tr("CSV (*.csv)"));
    if (path.isEmpty())
        return;
    try {
        const auto bytes = resultCsv(archive());
        QSaveFile file(path);
        if (!file.open(QIODevice::WriteOnly) || file.write(bytes) != bytes.size() || !file.commit())
            throw std::runtime_error(file.errorString().toStdString());
    } catch (const std::exception &e) {
        QMessageBox::warning(this, tr("CSV-Export"), QString::fromUtf8(e.what()));
    }
}

void SimulationDock::updateFrequencyPlots() {
    updateResultControls();
    PlotSeries onaxis{{}, QColor(28, 102, 180), {}}, real{tr("Re"), QColor(28, 102, 180), {}},
        imaginary{tr("Im"), QColor(205, 112, 33), {}};
    const auto face = impedanceSelector_->currentData().toUInt();
    for (const auto &[frequency, display] : displayFrequencies_) {
        onaxis.points.emplace_back(frequency, display.onAxisDb);
        if (const auto found = display.impedances.find(face); found != display.impedances.end()) {
            real.points.emplace_back(frequency, found->second.real());
            imaginary.points.emplace_back(frequency, found->second.imag());
        }
    }
    std::vector<PlotSeries> responses;
    if (!onaxis.points.empty())
        responses.push_back(std::move(onaxis));
    for (const auto &curve : measurements_) {
        PlotSeries measured{
            tr("M: %1%2").arg(curve.name.left(18)).arg(curve.absoluteSpl ? "" : tr(" (rel.)")),
            QColor(105, 105, 105),
            {}};
        for (const auto &sample : curve.samples)
            measured.points.emplace_back(sample.frequencyHz, sample.decibels + curve.offsetDb);
        responses.push_back(std::move(measured));
    }
    onAxis_->setSeries(std::move(responses));
    impedance_->setSeries({std::move(real), std::move(imaginary)});
    const auto widthSeries = [&](bool horizontal) {
        PlotSeries series{{}, horizontal ? QColor(28, 102, 180) : QColor(205, 112, 33), {}};
        for (const auto &[frequency, display] : displayFrequencies_) {
            const auto &levels = horizontal ? display.horizontalDb : display.verticalDb;
            std::vector<AcousticsAngularSample> samples;
            for (std::size_t i = 0; i < levels.size(); ++i)
                samples.push_back({-180 + 5 * static_cast<double>(i), levels[i]});
            const auto width =
                beamwidthDegrees(samples, beamwidthThreshold_->currentData().toDouble());
            series.points.emplace_back(frequency, width ? width->widthDegrees
                                                        : std::numeric_limits<double>::quiet_NaN());
        }
        return series;
    };
    beamwidthH_->setSeries({widthSeries(true)});
    beamwidthV_->setSeries({widthSeries(false)});
    const auto targetBand = [&](double value) -> std::optional<PlotReferenceBand> {
        if (!beamwidthTarget_->isChecked())
            return {};
        return PlotReferenceBand{std::max(0.0, value - targetTolerance_->value()),
                                 std::min(360.0, value + targetTolerance_->value()),
                                 QColor(51, 129, 75),
                                 tr("Ziel %1° ± %2°").arg(value).arg(targetTolerance_->value())};
    };
    beamwidthH_->setReferenceBand(targetBand(targetH_->value()));
    beamwidthV_->setReferenceBand(targetBand(targetV_->value()));
}

void SimulationDock::updateAngularPlots() {
    if (onLemResultChanged) {
        std::optional<lem::Solution> solution;
        const double selected = frequencySelector_->currentData().toDouble();
        if (result_)
            for (const auto &frequency : result_->frequencies)
                if (frequency.frequency == selected && !frequency.excitations.empty())
                    solution = frequency.excitations.front().lemSolution;
        onLemResultChanged(std::move(solution),
                           !hasCurrentGeometry() || resultGeneration_ != generation_);
    }
    const auto found = displayFrequencies_.find(frequencySelector_->currentData().toDouble());
    if (found == displayFrequencies_.end())
        return;
    const auto &display = found->second;
    polarH_->setLevels(displayedCut(display.horizontalDb, display.onAxisDb,
                                    normalization_->currentIndex() == 0,
                                    normalizationAngle_->value()));
    polarV_->setLevels(displayedCut(display.verticalDb, display.onAxisDb,
                                    normalization_->currentIndex() == 0,
                                    normalizationAngle_->value()));
}

void SimulationDock::updateHeatmaps() {
    std::vector<double> frequencies, horizontal, vertical;
    const auto validCount = angularCount;
    const auto appendRow = [validCount](std::vector<double> &target,
                                        const std::vector<AngularLevel> &cut) {
        const auto begin = target.size();
        target.resize(begin + validCount, std::numeric_limits<double>::quiet_NaN());
        if (!cut.empty() && cut.size() <= validCount) {
            const auto offset = (validCount - cut.size()) / 2;
            for (std::size_t i = 0; i < cut.size(); ++i)
                target[begin + offset + i] = cut[i].levelDecibels;
        }
    };
    for (const auto &[frequency, display] : displayFrequencies_) {
        frequencies.push_back(frequency);
        const auto h = displayedCut(display.horizontalDb, display.onAxisDb,
                                    normalization_->currentIndex() == 0,
                                    normalizationAngle_->value()),
                   v = displayedCut(display.verticalDb, display.onAxisDb,
                                    normalization_->currentIndex() == 0,
                                    normalizationAngle_->value());
        appendRow(horizontal, h);
        appendRow(vertical, v);
    }
    if (!frequencies.empty()) {
        isobarH_->setLevels(anglesDegrees(validCount), frequencies, horizontal);
        radialIsobarH_->setLevels(anglesDegrees(validCount), frequencies, horizontal);
        isobarV_->setLevels(anglesDegrees(validCount), frequencies, vertical);
        radialIsobarV_->setLevels(anglesDegrees(validCount), frequencies, vertical);
    } else {
        isobarH_->setLevels({}, {}, {});
        radialIsobarH_->setLevels({}, {}, {});
        isobarV_->setLevels({}, {}, {});
        radialIsobarV_->setLevels({}, {}, {});
    }
}

void SimulationDock::updateDirectivityScale() {
    try {
        const auto number = [](const QLineEdit *input) {
            bool valid = false;
            const double value = QString(input->text()).replace(',', '.').toDouble(&valid);
            if (!valid || !std::isfinite(value))
                throw std::invalid_argument(
                    "Pegelgrenzen und Schritt müssen endliche Zahlen sein.");
            return value;
        };
        applyDirectivityScale(
            {number(directivityMinimum_), number(directivityMaximum_), number(directivityStep_)});
    } catch (const std::invalid_argument &error) {
        setStatus(QString::fromUtf8(error.what()), true);
    }
}
void SimulationDock::applyDirectivityScale(PlotDbScale scale) {
    polarH_->setDecibelScale(scale);
    polarV_->setDecibelScale(scale);
    isobarH_->setDecibelScale(scale);
    radialIsobarH_->setDecibelScale(scale);
    isobarV_->setDecibelScale(scale);
    radialIsobarV_->setDecibelScale(scale);
    directivityScale_ = scale;
    directivityMinimum_->setText(QString::number(scale.minimum, 'g', 17));
    directivityMaximum_->setText(QString::number(scale.maximum, 'g', 17));
    directivityStep_->setText(QString::number(scale.tickStep, 'g', 17));
}
} // namespace hornsculpt
