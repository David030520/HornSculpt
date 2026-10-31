#include "parameter_sweep.hpp"
#include "galerkin/exterior_domain.hpp"
#include "numeric_input.hpp"
#include "result_export.hpp"
#include "simulation_request.hpp"
#include "sweep_model.hpp"
#include <QAbstractTableModel>
#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QGridLayout>
#include <QGuiApplication>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QMetaObject>
#include <QPointer>
#include <QProgressBar>
#include <QPushButton>
#include <QScrollArea>
#include <QScreen>
#include <QTabWidget>
#include <QSignalBlocker>
#include <QSplitter>
#include <QTableView>
#include <QVBoxLayout>
#include <QWindow>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <hornsculpt/acoustic_analysis.hpp>
#include <hornsculpt/cad_geometry.hpp>
#include <hornsculpt/heatmap_widget.hpp>
#include <hornsculpt/mesh_orientation.hpp>
#include <hornsculpt/plot_widget.hpp>
#include <limits>
#include <numbers>
#include <hornsculpt/polar_plot_widget.hpp>
#include <set>
#include <stdexcept>
#include <thread>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace hornsculpt {
namespace {
const SweepParameter &parameter(const QString &key) {
    const auto all = sweepParameters();
    const auto found =
        std::find_if(all.begin(), all.end(), [&](const auto &p) { return p.key == key; });
    if (found == all.end())
        throw std::invalid_argument("Unbekannter Studienparameter.");
    return *found;
}
QString stateText(SweepState state) {
    switch (state) {
    case SweepState::pending:
        return "Offen";
    case SweepState::running:
        return "In Arbeit";
    case SweepState::completed:
        return "Berechnet";
    case SweepState::failed:
        return "Fehler";
    case SweepState::cancelled:
        return "Abgebrochen";
    }
    return {};
}
class StudyTableModel final : public QAbstractTableModel {
  public:
    explicit StudyTableModel(QObject *parent) : QAbstractTableModel(parent) {}
    void setStudy(const SweepStudy *study) {
        beginResetModel();
        study_ = study;
        endResetModel();
    }
    void updateRow(std::size_t row) {
        emit dataChanged(index(static_cast<int>(row), 0), index(static_cast<int>(row), 3));
    }
    int rowCount(const QModelIndex &parent = {}) const override {
        return parent.isValid() || !study_ ? 0 : static_cast<int>(study_->variants.size());
    }
    int columnCount(const QModelIndex &parent = {}) const override {
        return parent.isValid() ? 0 : 4;
    }
    QVariant headerData(int section, Qt::Orientation orientation, int role) const override {
        if (role != Qt::DisplayRole)
            return {};
        if (orientation == Qt::Vertical)
            return section + 1;
        const QStringList headers{"Nr.", "Parameter", "Status", "Datei / Fehler"};
        return section >= 0 && section < headers.size() ? headers[section] : QVariant{};
    }
    QVariant data(const QModelIndex &index, int role) const override {
        if (!study_ || !index.isValid() || index.row() >= rowCount())
            return {};
        const auto &v = study_->variants[static_cast<std::size_t>(index.row())];
        if (role == Qt::ToolTipRole)
            return v.error.isEmpty() ? v.projectFile : v.error;
        if (role != Qt::DisplayRole)
            return {};
        switch (index.column()) {
        case 0:
            return index.row() + 1;
        case 1:
            return sweepVariantLabel(*study_, static_cast<std::size_t>(index.row()));
        case 2:
            return stateText(v.state);
        case 3:
            return v.error.isEmpty() ? v.projectFile : v.error.left(120);
        }
        return {};
    }

  private:
    const SweepStudy *study_{};
};
void resultMemoryCheck(const SurfaceMesh &mesh, std::size_t exteriorPanels, std::size_t count,
                       std::size_t drivenCount, std::size_t excitationCount) {
    // Variants are evaluated sequentially; one result owns its own shared blocks.
    const long double bytes = estimateResultBytes(
        {mesh.triangles.size(), mesh.vertices.size(), exteriorPanels, 147, drivenCount,
         count, excitationCount, 1});
    if (count > std::vector<BemFrequencyResult>{}.max_size() ||
        bytes > std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("Das Frequenzraster überschreitet den verfügbaren Adressraum.");
#ifdef _WIN32
    MEMORYSTATUSEX memory{};
    memory.dwLength = sizeof(memory);
    if (GlobalMemoryStatusEx(&memory) && bytes > memory.ullAvailPhys)
        throw std::invalid_argument(
            "Das Frequenzraster benötigt mehr verfügbaren Arbeitsspeicher.");
#endif
}
void checkStop(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Studie abgebrochen.");
}
Project calculateVariant(const SweepStudy &study, std::size_t index, std::stop_token stop,
                         const std::function<void(const BemProgress &)> &progress) {
    checkStop(stop);
    auto project = sweepVariantInputs(study, index);
    if (project.simulation.settings.geometryMode == BemGeometryMode::axisymmetric)
        throw std::invalid_argument(
            "Diese Parameterstudie benötigt den 3D-Geometriemodus. 2D-Rotation direkt im "
            "Akustikdock berechnen; es erfolgt kein Wechsel auf 3D.");
    validateMeshSettings(project.meshing);
    std::vector<FaceVelocity> velocities;
    Vec3 origin{}, axis{0, 0, 1};
    std::shared_ptr<const SurfaceMesh> source;
    std::shared_ptr<const SurfaceParameterization> parameterization;
    if (project.source == ProjectSource::draft) {
        origin = {0, 0, project.parameters.length};
        source = std::make_shared<SurfaceMesh>(
            DraftGeometry(project.parameters)
                .previewMesh(project.meshing.adaptive.geometricTolerance, stop));
        std::set<std::uint32_t> driven;
        for (const auto &t : source->triangles)
            if (t.boundary == BoundaryCondition::driven)
                driven.insert(t.faceId);
        for (const auto face : driven)
            velocities.push_back({face, {1, 0}});
    } else {
        if (!project.cad || project.cad->hasUnresolvedAssignments())
            throw std::invalid_argument("Bestätigter CAD-Snapshot erforderlich.");
        auto cad = std::make_shared<CadSnapshot>(*project.cad);
        cad->geometry = validateCadSnapshot(*cad);
        project.cad = std::move(cad);
        auto cadSource =
            makeCadMeshingSource(*project.cad, project.meshing.adaptive.geometricTolerance, stop);
        parameterization = std::move(cadSource.parameters);
        source = std::make_shared<SurfaceMesh>(std::move(cadSource.mesh));
        origin = project.cad->acousticOrigin;
        axis = project.cad->acousticAxis;
        for (const auto &face : project.cad->faces)
            if (face.boundary == BoundaryCondition::driven)
                velocities.push_back({face.id, face.normalVelocity});
    }
    if (velocities.empty())
        throw std::invalid_argument("Keine Antriebsfläche definiert.");
    MeshingOptions options;
    options.targetEdgeLength = project.meshing.targetEdgeLength;
    options.adaptive = project.meshing.adaptive;
    options.quality = project.meshing.quality;
    options.faceTargetEdgeLengths = project.source == ProjectSource::draft
                                        ? project.draftFaceTargetEdgeLengths
                                        : project.cadFaceTargetEdgeLengths;
    if (project.source == ProjectSource::freecad && project.cad)
        options.geometry = project.cad->geometry;
    options.parameterization = std::move(parameterization);
    auto mesh = std::make_shared<MeshingResult>(remeshIsotropic(*source, options, stop));
    checkStop(stop);
    const auto &reversed = project.source == ProjectSource::draft ? project.draftReversedFaceIds
                                                                  : project.cadReversedFaceIds;
    project.preview = project.source == ProjectSource::freecad
                          ? std::shared_ptr<const SurfaceMesh>(project.cad, &project.cad->mesh)
                          : std::move(source);
    if (!reversed.empty()) {
        try {
            mesh->mesh = applyFaceOrientation(mesh->mesh, reversed);
            mesh->quality = inspectMesh(mesh->mesh);
            if (!mesh->quality.watertight())
                throw std::invalid_argument(
                    "Flächenauswahl ergibt keine konsistent orientierte geschlossene Hülle.");
            galerkin::validateExteriorComponents(mesh->mesh);
            project.preview =
                std::make_shared<SurfaceMesh>(applyFaceOrientation(*project.preview, reversed));
        } catch (const std::exception &orientationError) {
            throw std::invalid_argument(std::string("Flächennormalen: ") + orientationError.what());
        }
    }
    checkStop(stop);
    auto settings = project.simulation.settings;
    if (!study.frequencyOverride.empty()) {
        settings.minimumFrequency = study.frequencyOverride.front();
        settings.maximumFrequency = study.frequencyOverride.back();
        settings.frequencyCount = study.frequencyOverride.size();
    }
    resultMemoryCheck(mesh->mesh,
                      settings.infiniteBaffle
                          ? settings.infiniteBaffle->exteriorDrivenMesh.triangles.size()
                          : 0,
                      settings.frequencyCount, velocities.size(),
                      settings.couplingEnabled && settings.infiniteBaffle && !velocities.empty()
                          ? 2 : 1);
    auto request = makeSimulationRequest(settings, velocities, origin, axis);
    if (!study.frequencyOverride.empty())
        request.frequencies = study.frequencyOverride;
    CudaBemSolver solver(mesh->mesh);
    auto result = std::make_shared<BemResult>(solver.solve(request, stop, progress));
    checkStop(stop);
    if (result->frequencies.size() != request.frequencies.size())
        throw std::runtime_error("Berechnung lieferte ein unvollständiges Frequenzraster.");
    for (std::size_t i = 0; i < request.frequencies.size(); ++i)
        if (result->frequencies[i].frequency != request.frequencies[i])
            throw std::runtime_error("Berechnung lieferte ein abweichendes Frequenzraster.");
    project.geometryRevision = project.revision;
    project.mesh = std::move(mesh);
    project.meshRevision = project.revision;
    project.simulation.settings = settings;
    project.simulation.mesh = std::make_shared<SurfaceMesh>(project.mesh->mesh);
    project.simulation.request = std::move(request);
    project.simulation.result = std::move(result);
    project.simulation.revision = project.revision;
    project.simulation.current = true;
    project.simulation.resultComplete = true;
    return project;
}
constexpr std::array<QRgb, 8> colors{0xff2166ac, 0xffc45b19, 0xff238b45, 0xff8e4a9f,
                                     0xffa83246, 0xff158a95, 0xff786a14, 0xff555555};
std::vector<double> normalizedCut(const std::vector<std::complex<double>> &pressure,
                                  std::size_t count, std::size_t cut, bool maximum) {
    const double missing = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> levels(count, missing);
    if (pressure.size() != 1 + 2 * count)
        return levels;
    double reference = std::abs(pressure.front());
    if (maximum) {
        reference = 0.;
        for (std::size_t i = 0; i < count; ++i) {
            const auto magnitude = std::abs(pressure[1 + cut * count + i]);
            if (std::isfinite(magnitude))
                reference = std::max(reference, magnitude);
        }
    }
    if (!std::isfinite(reference) || reference <= 0.)
        return levels;
    for (std::size_t i = 0; i < count; ++i) {
        const auto magnitude = std::abs(pressure[1 + cut * count + i]);
        if (std::isfinite(magnitude))
            levels[i] = magnitude > 0. ? 20. * std::log10(magnitude / reference)
                                       : -std::numeric_limits<double>::infinity();
    }
    return levels;
}
bool completeVariantResult(const SweepStudy &study, const Project &project) {
    const auto &simulation = project.simulation;
    if (!simulation.current || !simulation.resultComplete.value_or(true) || !simulation.mesh ||
        !simulation.request || !simulation.result)
        return false;
    const auto &requested = simulation.request->frequencies;
    const auto &actual = simulation.result->frequencies;
    if (requested.empty() || requested.size() != actual.size() ||
        (!study.frequencyOverride.empty() && requested != study.frequencyOverride))
        return false;
    for (std::size_t i = 0; i < requested.size(); ++i)
        if (requested[i] != actual[i].frequency)
            return false;
    return true;
}
} // namespace

struct ParameterSweepDialog::Implementation {
    ParameterSweepDialog *dialog;
    Project baseline;
    std::optional<SweepStudy> study;
    std::vector<double> frequencyOverride;
    std::jthread worker;
    bool busy{}, convergence{}, closeWhenIdle{}, waitingForGpu{}, gpuLease{};
    std::uint64_t acquisitionId{};
    std::array<QComboBox *, 2> parameters{};
    std::array<DecimalInput *, 2> starts{}, ends{};
    std::array<IntegerInput *, 2> counts{};
    QCheckBox *second{};
    QLineEdit *directory{};
    QPushButton *browse{}, *start{}, *cancel{}, *open{}, *adopt{}, *analyze{};
    QLabel *status{}, *difference{};
    QLabel *parameterMode{}, *convergenceAxis{};
    QProgressBar *progress{};
    QTableView *table{};
    StudyTableModel *model{};
    QComboBox *impedanceFace{}, *beamwidthThreshold{};
    QCheckBox *targetEnabled{};
    std::array<PlotWidget *, 4> plots{};
    QComboBox *polarFrequency{}, *polarNormalization{};
    QComboBox *mapVariant{}, *mapNormalization{};
    QLabel *polarStatus{}, *polarLegend{};
    QLabel *mapStatus{};
    PolarPlotWidget *polarH{}, *polarV{};
    HeatmapWidget *mapH{}, *mapV{};
    struct PolarInput { std::size_t row; Project project; };
    std::vector<PolarInput> polarInputs;

    Implementation(ParameterSweepDialog *owner, Project input);
    ~Implementation() {
        worker.request_stop();
        if (worker.joinable())
            worker.join();
        if (waitingForGpu && dialog->onCancelGpuAcquire) dialog->onCancelGpuAcquire();
        if (gpuLease && dialog->onReleaseGpu) dialog->onReleaseGpu();
    }
    void axisChanged(std::size_t axis);
    void rebuildParameterCatalog(const std::array<QString, 2> &keys);
    void updateControls();
    void begin();
    void beginWorker();
    void cancelWaiting();
    void load(const QString &path);
    void compare();
    void configure(std::vector<double> frequencies);
    void acceptVariant();
    void analyzeVariant();
    void updatePolarFrequencies();
    void updatePolars();
    void updateMapVariants();
    void updateMaps();
    void applyStudyScale();
    void run(SweepStudy task, std::stop_token stop);
};

ParameterSweepDialog::Implementation::Implementation(ParameterSweepDialog *owner, Project input)
    : dialog(owner), baseline(std::move(input)) {
    dialog->setWindowTitle("Parameterstudie");
    dialog->setObjectName("parameterSweepDialog");
    QScreen *screen = nullptr;
    if (auto *parent = dialog->parentWidget()) {
        if (auto *window = parent->window()->windowHandle())
            screen = window->screen();
        if (!screen)
            screen = QGuiApplication::screenAt(parent->mapToGlobal(parent->rect().center()));
    }
    if (!screen)
        screen = QGuiApplication::primaryScreen();
    const QSize available = screen ? screen->availableGeometry().size() : QSize(1180, 860);
    // Leave space for the native frame and title bar; screen coordinates are Qt logical pixels.
    dialog->resize(std::min(1180, std::max(1, available.width() - 32)),
                   std::min(860, std::max(1, available.height() - 80)));
    auto *outerLayout = new QVBoxLayout(dialog);
    outerLayout->setContentsMargins(0, 0, 0, 0);
    auto *scroll = new QScrollArea;
    scroll->setObjectName("sweepDialogScroll");
    scroll->setFrameShape(QFrame::NoFrame);
    scroll->setWidgetResizable(true);
    outerLayout->addWidget(scroll);
    auto *content = new QWidget;
    auto *layout = new QVBoxLayout(content);
    scroll->setWidget(content);
    parameterMode = new QLabel;
    parameterMode->setObjectName("sweepParameterMode");
    parameterMode->setWordWrap(true);
    layout->addWidget(parameterMode);
    auto *axes = new QGridLayout;
    for (int column = 0; column < 4; ++column)
        axes->addWidget(new QLabel(QStringList{"Parameter", "Start", "Ende", "Punkte"}[column]), 0,
                        column);
    for (std::size_t a = 0; a < 2; ++a) {
        parameters[a] = new QComboBox;
        parameters[a]->setObjectName(QString("sweepParameter%1").arg(a + 1));
        starts[a] = new DecimalInput;
        ends[a] = new DecimalInput;
        counts[a] = new IntegerInput;
        starts[a]->setObjectName(QString("sweepStart%1").arg(a + 1));
        ends[a]->setObjectName(QString("sweepEnd%1").arg(a + 1));
        counts[a]->setObjectName(QString("sweepPoints%1").arg(a + 1));
        for (auto *value : {starts[a], ends[a]}) {
            value->setDecimals(8);
            value->setRange(-1e9, 1e9);
        }
        counts[a]->setRange(1, std::numeric_limits<int>::max());
        counts[a]->setValue(1);
        axes->addWidget(parameters[a], static_cast<int>(a) + 1, 0);
        axes->addWidget(starts[a], static_cast<int>(a) + 1, 1);
        axes->addWidget(ends[a], static_cast<int>(a) + 1, 2);
        axes->addWidget(counts[a], static_cast<int>(a) + 1, 3);
        QObject::connect(parameters[a], &QComboBox::currentIndexChanged, dialog,
                         [this, a] { axisChanged(a); });
    }
    convergenceAxis = new QLabel(QObject::tr("Zielkante [mm]"));
    convergenceAxis->setObjectName("convergenceMeshAxis");
    convergenceAxis->setToolTip(QObject::tr("Grobe und feine Zielkante für dasselbe Modell. "
                                           "Die Zwischenstufen werden linear erzeugt."));
    axes->addWidget(convergenceAxis, 1, 0);
    convergenceAxis->hide();
    rebuildParameterCatalog({QStringLiteral("length"), QStringLiteral("horizontalAngle")});
    second = new QCheckBox("Zweite Achse");
    second->setObjectName("sweepSecondAxis");
    axes->addWidget(second, 3, 0);
    layout->addLayout(axes);
    QObject::connect(second, &QCheckBox::toggled, dialog, [this] { updateControls(); });
    auto *paths = new QHBoxLayout;
    paths->addWidget(new QLabel("Ausgabe"));
    directory = new QLineEdit;
    directory->setObjectName("sweepOutputDirectory");
    paths->addWidget(directory, 1);
    browse = new QPushButton("Verzeichnis …");
    paths->addWidget(browse);
    open = new QPushButton("Studie öffnen …");
    paths->addWidget(open);
    layout->addLayout(paths);
    QObject::connect(browse, &QPushButton::clicked, dialog, [this] {
        const auto path =
            QFileDialog::getExistingDirectory(dialog, "Studienverzeichnis", directory->text());
        if (!path.isEmpty())
            directory->setText(path);
    });
    QObject::connect(open, &QPushButton::clicked, dialog, [this] {
        const auto path = QFileDialog::getOpenFileName(dialog, "Studie öffnen", {},
                                                       "HornSculpt-Studie (*.hsstudy)");
        if (!path.isEmpty())
            try {
                load(path);
            } catch (const std::exception &e) {
                QMessageBox::critical(dialog, "Studie öffnen", QString::fromUtf8(e.what()));
            }
    });
    auto *split = new QSplitter(Qt::Vertical);
    layout->addWidget(split, 1);
    table = new QTableView;
    table->setObjectName("sweepVariants");
    model = new StudyTableModel(table);
    table->setModel(model);
    table->setSelectionBehavior(QAbstractItemView::SelectRows);
    table->setSelectionMode(QAbstractItemView::ExtendedSelection);
    table->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    table->horizontalHeader()->setSectionResizeMode(3, QHeaderView::Stretch);
    table->setColumnWidth(0, 45);
    table->setColumnWidth(2, 100);
    split->addWidget(table);
    QObject::connect(table->selectionModel(), &QItemSelectionModel::selectionChanged, dialog,
                     [this] { compare(); });
    auto *chartArea = new QWidget;
    auto *chartLayout = new QGridLayout(chartArea);
    impedanceFace = new QComboBox;
    auto *impedanceLabel = new QLabel("Impedanzfläche");
    auto *selector = new QHBoxLayout;
    selector->addWidget(impedanceLabel);
    selector->addWidget(impedanceFace);
    selector->addStretch();
    targetEnabled = new QCheckBox("Zielkorridor");
    targetEnabled->setObjectName("sweepTargetBand");
    targetEnabled->setChecked(baseline.simulation.settings.beamwidthTargetEnabled);
    targetEnabled->setToolTip(
        "H/V-Ziel und Toleranz aus den gemeinsamen Simulationseinstellungen.");
    selector->addWidget(targetEnabled);
    beamwidthThreshold = new QComboBox;
    beamwidthThreshold->setObjectName("sweepBeamwidthThreshold");
    for (const int attenuation : {3, 6, 10})
        beamwidthThreshold->addItem(QString("−%1 dB").arg(attenuation), attenuation);
    beamwidthThreshold->setCurrentIndex(1);
    beamwidthThreshold->setToolTip(
        "Schwelle gegenüber 0°; erste Übertritte der verbundenen Hauptkeule.");
    selector->addWidget(beamwidthThreshold);
    QObject::connect(beamwidthThreshold, &QComboBox::currentIndexChanged, dialog,
                     [this] { compare(); });
    QObject::connect(targetEnabled, &QCheckBox::toggled, dialog, [this] { compare(); });
    chartLayout->addLayout(selector, 0, 0, 1, 2);
    plots = {new PlotWidget("On-Axis · 1 m", "f [Hz]", "SPL [dB]"),
             new PlotWidget("Öffnungswinkel H", "f [Hz]", "Breite [°]"),
             new PlotWidget("Öffnungswinkel V", "f [Hz]", "Breite [°]"),
             new PlotWidget("Strahlungsimpedanz · Re / Im", "f [Hz]", "Z / (ρc/S)")};
    for (std::size_t i = 0; i < plots.size(); ++i) {
        if (i == 1 || i == 2)
            plots[i]->setAxisScales({20, 20000, {}}, {0, 180, 30});
        plots[i]->setObjectName(QStringList{"sweepOnAxis", "sweepBeamwidthH", "sweepBeamwidthV",
                                            "sweepImpedance"}[static_cast<qsizetype>(i)]);
        plots[i]->setFrequencyAxis();
        chartLayout->addWidget(plots[i], 1 + static_cast<int>(i / 2), static_cast<int>(i % 2));
    }
    auto *chartTabs = new QTabWidget;
    chartTabs->setObjectName("sweepChartViews");
    chartTabs->addTab(chartArea, "Übersicht");
    auto *polarPage = new QWidget;
    auto *polarLayout = new QVBoxLayout(polarPage);
    auto *polarControls = new QHBoxLayout;
    polarControls->addWidget(new QLabel("Frequenz"));
    polarFrequency = new QComboBox;
    polarFrequency->setObjectName("sweepPolarFrequency");
    polarControls->addWidget(polarFrequency);
    polarControls->addWidget(new QLabel("Bezug"));
    polarNormalization = new QComboBox;
    polarNormalization->setObjectName("sweepPolarNormalization");
    polarNormalization->addItem("0° je Variante", 0);
    polarNormalization->addItem("Maximum je Schnitt", 1);
    polarControls->addWidget(polarNormalization);
    polarControls->addStretch();
    polarLayout->addLayout(polarControls);
    auto *polarScroll = new QScrollArea;
    polarScroll->setWidgetResizable(true);
    auto *polarContent = new QWidget;
    auto *polarPlots = new QHBoxLayout(polarContent);
    polarH = new PolarPlotWidget("Horizontal", QColor(33, 102, 172));
    polarV = new PolarPlotWidget("Vertikal", QColor(33, 102, 172));
    polarH->setObjectName("sweepPolarH");
    polarV->setObjectName("sweepPolarV");
    polarPlots->addWidget(polarH);
    polarPlots->addWidget(polarV);
    polarScroll->setWidget(polarContent);
    polarLayout->addWidget(polarScroll, 1);
    polarLegend = new QLabel;
    polarLegend->setWordWrap(true);
    polarLayout->addWidget(polarLegend);
    polarStatus = new QLabel("Varianten und Frequenz auswählen.");
    polarStatus->setWordWrap(true);
    polarLayout->addWidget(polarStatus);
    chartTabs->addTab(polarPage, "Polars H/V");
    auto *mapPage = new QWidget;
    auto *mapLayout = new QVBoxLayout(mapPage);
    auto *mapControls = new QHBoxLayout;
    mapControls->addWidget(new QLabel("Variante"));
    mapVariant = new QComboBox;
    mapVariant->setObjectName("sweepMapVariant");
    mapControls->addWidget(mapVariant, 1);
    mapControls->addWidget(new QLabel("Bezug"));
    mapNormalization = new QComboBox;
    mapNormalization->setObjectName("sweepMapNormalization");
    mapNormalization->addItem("0° je Variante", 0);
    mapNormalization->addItem("Maximum je Schnitt", 1);
    mapControls->addWidget(mapNormalization);
    mapLayout->addLayout(mapControls);
    auto *mapScroll = new QScrollArea;
    mapScroll->setWidgetResizable(true);
    auto *mapContent = new QWidget;
    auto *mapPlots = new QHBoxLayout(mapContent);
    mapH = new HeatmapWidget("Isobaren H");
    mapV = new HeatmapWidget("Isobaren V");
    mapH->setObjectName("sweepMapH");
    mapV->setObjectName("sweepMapV");
    mapPlots->addWidget(mapH);
    mapPlots->addWidget(mapV);
    mapScroll->setWidget(mapContent);
    mapLayout->addWidget(mapScroll, 1);
    mapStatus = new QLabel("Variante auswählen.");
    mapStatus->setWordWrap(true);
    mapLayout->addWidget(mapStatus);
    chartTabs->addTab(mapPage, "Polarmaps H/V");
    split->addWidget(chartTabs);
    QObject::connect(polarFrequency, &QComboBox::currentIndexChanged, dialog,
                     [this] { updatePolars(); });
    QObject::connect(polarNormalization, &QComboBox::currentIndexChanged, dialog,
                     [this] {
                         mapNormalization->setCurrentIndex(polarNormalization->currentIndex());
                         updatePolars();
                         updateMaps();
                     });
    QObject::connect(mapNormalization, &QComboBox::currentIndexChanged, dialog,
                     [this] { polarNormalization->setCurrentIndex(mapNormalization->currentIndex()); });
    QObject::connect(mapVariant, &QComboBox::currentIndexChanged, dialog,
                     [this] { updateMaps(); });
    const auto setAllScales = [this](PlotDbScale scale) {
        polarH->setDecibelScale(scale);
        polarV->setDecibelScale(scale);
        mapH->setDecibelScale(scale);
        mapV->setDecibelScale(scale);
    };
    polarH->onDecibelScaleChanged = setAllScales;
    polarV->onDecibelScaleChanged = setAllScales;
    mapH->onDecibelScaleChanged = setAllScales;
    mapV->onDecibelScaleChanged = setAllScales;
    applyStudyScale();
    split->setStretchFactor(0, 1);
    split->setStretchFactor(1, 3);
    QObject::connect(impedanceFace, &QComboBox::currentIndexChanged, dialog, [this] { compare(); });
    difference = new QLabel;
    difference->setWordWrap(true);
    difference->setToolTip(
        "Differenzen beziehen sich auf die erste markierte Variante bei identischen Frequenzen.");
    layout->addWidget(difference);
    auto *bottom = new QHBoxLayout;
    start = new QPushButton("Start");
    cancel = new QPushButton("Abbrechen");
    adopt = new QPushButton("Variante übernehmen");
    analyze = new QPushButton("Im Hauptfenster analysieren");
    analyze->setObjectName("sweepAnalyzeVariant");
    start->setObjectName("sweepStart");
    cancel->setObjectName("sweepCancel");
    adopt->setObjectName("sweepAdopt");
    start->setToolTip(
        "Jeder Start legt ein neues Studienverzeichnis an. Vorhandene Dateien bleiben erhalten.");
    bottom->addWidget(start);
    bottom->addWidget(cancel);
    bottom->addWidget(adopt);
    bottom->addWidget(analyze);
    bottom->addStretch();
    auto *close = new QPushButton("Schließen");
    bottom->addWidget(close);
    layout->addLayout(bottom);
    progress = new QProgressBar;
    progress->setObjectName("sweepProgress");
    progress->setRange(0, 1);
    progress->setValue(0);
    layout->addWidget(progress);
    status = new QLabel("Bereit");
    status->setWordWrap(true);
    status->setObjectName("sweepStatus");
    status->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(status);
    QObject::connect(start, &QPushButton::clicked, dialog, [this] { begin(); });
    QObject::connect(cancel, &QPushButton::clicked, dialog, [this] {
        if (waitingForGpu) { cancelWaiting(); return; }
        worker.request_stop();
        cancel->setEnabled(false);
        status->setText("Abbruch angefordert …");
    });
    QObject::connect(adopt, &QPushButton::clicked, dialog, [this] { acceptVariant(); });
    QObject::connect(analyze, &QPushButton::clicked, dialog, [this] { analyzeVariant(); });
    QObject::connect(close, &QPushButton::clicked, dialog, &QDialog::reject);
    // Enter in a numeric editor confirms the edit; it must never launch a study.
    for (auto *button : dialog->findChildren<QPushButton *>()) {
        button->setAutoDefault(false);
        button->setDefault(false);
    }
    updateControls();
    updateMapVariants();
}
void ParameterSweepDialog::Implementation::rebuildParameterCatalog(const std::array<QString, 2> &keys) {
    const bool asymmetric = baseline.parameters.asymmetricEnabled;
    parameterMode->setText(baseline.source == ProjectSource::freecad
                               ? QObject::tr("CAD-Basis · Netzkonvergenz")
                               : asymmetric ? QObject::tr("Draft · asymmetrisch: rechts, links, oben, unten")
                                            : QObject::tr("Draft · symmetrisch: H und V gelten für beide Seiten"));
    const auto catalogue = sweepParameters(baseline.parameters);
    for (std::size_t a = 0; a < parameters.size(); ++a) {
        const QSignalBlocker blocker(parameters[a]);
        parameters[a]->clear();
        parameters[a]->setToolTip(asymmetric
            ? QObject::tr("H rechts/links und V oben/unten werden getrennt variiert. Winkel sind geometrische Halbwinkel in Grad.")
            : QObject::tr("H und V gelten jeweils für beide Seiten. Winkel sind geometrische Halbwinkel in Grad."));
        for (const auto &p : catalogue)
            parameters[a]->addItem(p.label + (p.unit.isEmpty() ? QString{} : " [" + p.unit + "]"), p.key);
        const int selected = parameters[a]->findData(keys[a]);
        if (selected < 0)
            throw std::invalid_argument("Studienparameter passt nicht zum Symmetriemodus der Basis.");
        parameters[a]->setCurrentIndex(selected);
    }
    axisChanged(0);
    axisChanged(1);
}
void ParameterSweepDialog::Implementation::axisChanged(std::size_t a) {
    const auto key = parameters[a]->currentData().toString();
    if (key.isEmpty())
        return;
    const auto &p = parameter(key);
    const double value = key == "meshTargetEdgeLength"
                             ? baseline.meshing.targetEdgeLength
                             : sweepParameterValue(baseline.parameters, key);
    starts[a]->setValue(value * p.displayScale);
    ends[a]->setValue(value * p.displayScale);
}
void ParameterSweepDialog::Implementation::updateControls() {
    if (convergence) {
        QStringList frequencies;
        for (const auto frequency : frequencyOverride)
            frequencies << QString::number(frequency, 'g', 8);
        parameterMode->setText(QObject::tr("Netzkonvergenz · Prüffrequenzen: %1 Hz · Geometrie bleibt gleich")
                                   .arg(frequencies.join(", ")));
    } else {
        parameterMode->setText(baseline.source == ProjectSource::freecad
                                   ? QObject::tr("CAD-Basis · Netzkonvergenz")
                                   : baseline.parameters.asymmetricEnabled
                                       ? QObject::tr("Draft · asymmetrisch: rechts, links, oben, unten")
                                       : QObject::tr("Draft · symmetrisch: H und V gelten für beide Seiten"));
    }
    for (std::size_t a = 0; a < 2; ++a) {
        const bool enabled = !busy && (a == 0 || second->isChecked());
        parameters[a]->setVisible(!convergence || a != 0);
        if (a == 1) {
            parameters[a]->setVisible(!convergence);
            starts[a]->setVisible(!convergence);
            ends[a]->setVisible(!convergence);
            counts[a]->setVisible(!convergence);
        }
        parameters[a]->setEnabled(enabled && !convergence);
        starts[a]->setEnabled(enabled);
        ends[a]->setEnabled(enabled);
        counts[a]->setEnabled(enabled);
    }
    convergenceAxis->setVisible(convergence);
    second->setVisible(!convergence);
    second->setEnabled(!busy && !convergence);
    counts[0]->setMinimum(convergence ? 2 : 1);
    starts[0]->setToolTip(convergence ? QObject::tr("Grobe Zielkante in mm; muss größer als die feine Zielkante sein.") : QString{});
    ends[0]->setToolTip(convergence ? QObject::tr("Feine Zielkante in mm; größer als null und kleiner als die grobe Zielkante.") : QString{});
    counts[0]->setToolTip(convergence ? QObject::tr("Mindestens zwei linear verteilte Netzstufen; für alle gelten dieselben Prüffrequenzen.") : QString{});
    directory->setEnabled(!busy);
    browse->setEnabled(!busy);
    start->setEnabled(!busy);
    open->setEnabled(!busy);
    cancel->setEnabled(busy);
    const auto selected = table->selectionModel()->selectedRows();
    const bool completed =
        study && selected.size() == 1 &&
        study->variants[static_cast<std::size_t>(selected.front().row())].state ==
            SweepState::completed;
    const bool validSelection = completed && polarInputs.size() == 1 &&
        polarInputs.front().row == static_cast<std::size_t>(selected.front().row());
    adopt->setEnabled(!busy && validSelection && static_cast<bool>(dialog->onVariantChosen));
    analyze->setEnabled(!busy && validSelection &&
                        static_cast<bool>(dialog->onVariantAnalysisRequested));
}
void ParameterSweepDialog::Implementation::begin() {
    try {
        std::vector<SweepAxis> axes;
        for (std::size_t a = 0; a < (second->isChecked() ? 2u : 1u); ++a) {
            const auto key = parameters[a]->currentData().toString();
            const auto &p = parameter(key);
            axes.push_back({key, starts[a]->value() / p.displayScale,
                            ends[a]->value() / p.displayScale,
                            static_cast<std::size_t>(counts[a]->value())});
        }
        if (baseline.source != ProjectSource::draft && !convergence)
            throw std::invalid_argument("Geometrieparameterstudien sind im Draft-Modus verfügbar.");
        auto common = baseline;
        common.simulation.settings.beamwidthTargetEnabled = targetEnabled->isChecked();
        study = createSweepStudy(directory->text(), common, std::move(axes), frequencyOverride);
        model->setStudy(&*study);
        progress->setRange(0, static_cast<int>(study->variants.size()));
        progress->setValue(0);
        busy = true;
        updateControls();
        if (!dialog->onAcquireGpu) { beginWorker(); return; }
        waitingForGpu = true;
        const auto id = ++acquisitionId;
        const QPointer<ParameterSweepDialog> receiver(dialog);
        status->setText(dialog->gpuWaitMessage ? dialog->gpuWaitMessage()
                                               : QStringLiteral("GPU-Speicher wird freigegeben …"));
        dialog->onAcquireGpu([this, receiver, id](bool acquired) {
            if (!receiver || !waitingForGpu || id != acquisitionId) return;
            waitingForGpu = false;
            if (!acquired) {
                busy = false;
                status->setText("GPU-Übergabe abgebrochen");
                updateControls();
                return;
            }
            gpuLease = true;
            beginWorker();
        });
    } catch (const std::exception &e) {
        if (waitingForGpu && dialog->onCancelGpuAcquire) dialog->onCancelGpuAcquire();
        waitingForGpu = false;
        if (gpuLease && dialog->onReleaseGpu) dialog->onReleaseGpu();
        gpuLease = false;
        busy = false;
        updateControls();
        QMessageBox::critical(dialog, "Studie starten", QString::fromUtf8(e.what()));
    }
}
void ParameterSweepDialog::Implementation::beginWorker() {
    try {
        if (worker.joinable()) worker.join();
        status->setText(sweepManifestPath(*study));
        worker = std::jthread([this, task = *study](std::stop_token stop) { run(task, stop); });
    } catch (const std::exception &error) {
        if (gpuLease && dialog->onReleaseGpu) dialog->onReleaseGpu();
        gpuLease = false;
        busy = false;
        updateControls();
        QMessageBox::critical(dialog, "Studie starten", QString::fromUtf8(error.what()));
    }
}
void ParameterSweepDialog::Implementation::cancelWaiting() {
    if (!waitingForGpu) return;
    waitingForGpu = false;
    ++acquisitionId;
    if (dialog->onCancelGpuAcquire) dialog->onCancelGpuAcquire();
    busy = false;
    status->setText("GPU-Übergabe abgebrochen");
    updateControls();
}
void ParameterSweepDialog::Implementation::run(SweepStudy task, std::stop_token stop) {
    QString error;
    try {
        for (std::size_t i = 0; i < task.variants.size(); ++i) {
            if (stop.stop_requested())
                break;
            auto &variant = task.variants[i];
            variant.state = SweepState::running;
            writeSweepStudy(task);
            const auto publish = [this, i](SweepVariant value) {
                QMetaObject::invokeMethod(
                    dialog,
                    [this, i, value = std::move(value)] {
                        if (!study)
                            return;
                        study->variants[i] = value;
                        model->updateRow(i);
                        if (value.state != SweepState::running) {
                            progress->setValue(static_cast<int>(i + 1));
                            if (table->selectionModel()->selectedRows().empty() &&
                                value.state == SweepState::completed)
                                table->selectRow(static_cast<int>(i));
                            compare();
                        }
                    },
                    Qt::QueuedConnection);
            };
            publish(variant);
            try {
                auto last = std::chrono::steady_clock::now();
                auto result = calculateVariant(task, i, stop, [&, i](const BemProgress &p) {
                    const auto now = std::chrono::steady_clock::now();
                    if (now - last < std::chrono::milliseconds(100) &&
                        p.completedFrequencies != p.totalFrequencies)
                        return;
                    last = now;
                    const auto text = QString("Variante %1/%2 · %3/%4 Frequenzen")
                                          .arg(i + 1)
                                          .arg(task.variants.size())
                                          .arg(p.completedFrequencies)
                                          .arg(p.totalFrequencies);
                    QMetaObject::invokeMethod(
                        dialog, [this, text] { status->setText(text); }, Qt::QueuedConnection);
                });
                const auto path = QDir(task.directory).filePath(variant.projectFile);
                if (QFileInfo::exists(path))
                    throw std::runtime_error(
                        "Die Variantendatei existiert bereits; sie wird nicht überschrieben.");
                writeProject(path, result);
                variant.state = SweepState::completed;
            } catch (const std::exception &e) {
                variant.state = stop.stop_requested() ? SweepState::cancelled : SweepState::failed;
                variant.error = QString::fromUtf8(e.what()).left(4096);
            }
            writeSweepStudy(task);
            publish(variant);
            if (stop.stop_requested())
                break;
        }
    } catch (const std::exception &e) {
        error = QString::fromUtf8(e.what());
    }
    QMetaObject::invokeMethod(
        dialog,
        [this, error, stopped = stop.stop_requested()] {
            if (worker.joinable()) worker.join();
            if (gpuLease && dialog->onReleaseGpu) dialog->onReleaseGpu();
            gpuLease = false;
            busy = false;
            updateControls();
            status->setText(!error.isEmpty() ? error
                            : stopped        ? "Abgebrochen"
                                             : "Studie abgeschlossen");
            if (closeWhenIdle)
                dialog->reject();
        },
        Qt::QueuedConnection);
}
void ParameterSweepDialog::Implementation::load(const QString &path) {
    if (busy)
        throw std::runtime_error("Laufende Studie zuerst abbrechen.");
    auto loaded = readSweepStudy(path);
    const QSignalBlocker selectionBlocker(table->selectionModel());
    const QSignalBlocker targetBlocker(targetEnabled);
    const QSignalBlocker secondBlocker(second);
    model->setStudy(nullptr);
    polarInputs.clear();
    baseline = loaded.baseline;
    applyStudyScale();
    study = std::move(loaded);
    const std::array<QString, 2> keys{
        study->axes.empty() ? QStringLiteral("length") : study->axes[0].parameter,
        study->axes.size() > 1 ? study->axes[1].parameter : QStringLiteral("horizontalAngle")};
    rebuildParameterCatalog(keys);
    for (auto &variant : study->variants)
        if (variant.state == SweepState::running) {
            variant.state = SweepState::cancelled;
            variant.error = "Unterbrochener Auftrag";
        }
    targetEnabled->setChecked(baseline.simulation.settings.beamwidthTargetEnabled);
    frequencyOverride = study->frequencyOverride;
    convergence = !frequencyOverride.empty();
    counts[0]->setMinimum(convergence ? 2 : 1);
    dialog->setWindowTitle(convergence ? "Netzkonvergenz" : "Parameterstudie");
    second->setChecked(study->axes.size() == 2);
    for (std::size_t a = 0; a < study->axes.size(); ++a) {
        const auto &axis = study->axes[a];
        parameters[a]->setCurrentIndex(parameters[a]->findData(axis.parameter));
        starts[a]->setValue(axis.start * parameter(axis.parameter).displayScale);
        ends[a]->setValue(axis.end * parameter(axis.parameter).displayScale);
        counts[a]->setValue(static_cast<int>(axis.count));
    }
    directory->setText(QFileInfo(study->directory).absolutePath());
    model->setStudy(&*study);
    status->setText(path);
    progress->setRange(0, static_cast<int>(study->variants.size()));
    progress->setValue(static_cast<int>(
        std::count_if(study->variants.begin(), study->variants.end(),
                      [](const auto &v) { return v.state == SweepState::completed; })));
    for (std::size_t i = 0; i < study->variants.size(); ++i)
        if (study->variants[i].state == SweepState::completed) {
            table->selectRow(static_cast<int>(i));
            break;
        }
    updateControls();
    compare();
}
void ParameterSweepDialog::Implementation::configure(std::vector<double> frequencies) {
    if (busy || frequencies.empty())
        throw std::invalid_argument("Ausgewählte Frequenzen erforderlich.");
    double previous = 0;
    for (const double f : frequencies) {
        if (!std::isfinite(f) || f <= previous)
            throw std::invalid_argument("Frequenzen müssen positiv und streng aufsteigend sein.");
        previous = f;
    }
    const double fine = baseline.meshing.targetEdgeLength / 2;
    if (baseline.source == ProjectSource::freecad &&
        (!baseline.cad || baseline.cad->hasUnresolvedAssignments()))
        throw std::invalid_argument("Bestätigter CAD-Snapshot erforderlich.");
    convergence = true;
    frequencyOverride = std::move(frequencies);
    second->setChecked(false);
    parameters[0]->setCurrentIndex(parameters[0]->findData("meshTargetEdgeLength"));
    starts[0]->setValue(baseline.meshing.targetEdgeLength * 1000);
    ends[0]->setValue(fine * 1000);
    counts[0]->setValue(2);
    dialog->setWindowTitle("Netzkonvergenz");
    updateControls();
}
void ParameterSweepDialog::Implementation::acceptVariant() {
    if (!study || busy)
        return;
    const auto rows = table->selectionModel()->selectedRows();
    if (rows.size() != 1)
        return;
    const auto &v = study->variants[static_cast<std::size_t>(rows.front().row())];
    if (v.state != SweepState::completed)
        return;
    try {
        auto project = readProject(QDir(study->directory).filePath(v.projectFile));
        if (!sweepVariantMatches(*study, static_cast<std::size_t>(rows.front().row()), project) ||
            !completeVariantResult(*study, project))
            throw std::invalid_argument("Gespeicherte, gültige Variantensimulation fehlt.");
        if (dialog->onVariantChosen)
            dialog->onVariantChosen(std::move(project));
    } catch (const std::exception &e) {
        QMessageBox::critical(dialog, "Variante übernehmen", QString::fromUtf8(e.what()));
    }
}
void ParameterSweepDialog::Implementation::compare() {
    std::array<std::vector<PlotSeries>, 4> series;
    polarInputs.clear();
    struct Values {
        std::size_t row{};
        std::vector<double> frequencies, levels;
        std::vector<std::complex<double>> pressure, impedance;
    };
    std::vector<Values> compared;
    QStringList errors;
    if (study) {
        auto rows = table->selectionModel()->selectedRows();
        std::sort(rows.begin(), rows.end(),
                  [](const auto &a, const auto &b) { return a.row() < b.row(); });
        bool faceListSet = false;
        for (const auto &row : rows) {
            const auto index = static_cast<std::size_t>(row.row());
            const auto &variant = study->variants[index];
            if (variant.state != SweepState::completed)
                continue;
            try {
                auto project =
                    readProject(QDir(study->directory).filePath(variant.projectFile));
                if (!completeVariantResult(*study, project))
                    throw std::invalid_argument("Gespeicherte Simulation fehlt.");
                if (!sweepVariantMatches(*study, index, project))
                    throw std::invalid_argument("Variantendatei passt nicht zur Studienbasis.");
                const auto &result = *project.simulation.result;
                if (result.frequencies.empty())
                    continue;
                if (!faceListSet) {
                    const auto selected = impedanceFace->currentData();
                    const QSignalBlocker block(impedanceFace);
                    impedanceFace->clear();
                    if (!result.frequencies.front().excitations.empty())
                        for (const auto &z :
                             result.frequencies.front().excitations.front().impedances)
                            impedanceFace->addItem(QString("Fläche %1").arg(z.faceId), z.faceId);
                    const int previous = impedanceFace->findData(selected);
                    if (previous >= 0)
                        impedanceFace->setCurrentIndex(previous);
                    faceListSet = true;
                }
                const auto color = QColor::fromRgba(colors[index % colors.size()]);
                const auto &field = project.simulation.request->farField;
                const auto label =
                    QString("#%1 · %2 · %3")
                        .arg(index + 1)
                        .arg(sweepVariantLabel(*study, index))
                        .arg(field ? QString("Fernfeld / %1 m").arg(field->referenceDistance)
                                   : QString("endlicher Abstand"));
                PlotSeries on{label, color, {}}, h{label, color, {}}, v{label, color, {}},
                    real{QString("#%1 Re").arg(index + 1), color, {}},
                    imaginary{QString("#%1 Im").arg(index + 1), color.lighter(150), {}};
                Values values;
                values.row = index;
                const double missing = std::numeric_limits<double>::quiet_NaN();
                for (const auto &frequency : result.frequencies) {
                    if (frequency.excitations.empty())
                        continue;
                    const auto &excitation = frequency.excitations.front();
                    const auto &pressure = excitation.observerPressure;
                    const auto axisPressure = pressure.empty()
                                                  ? std::complex<double>{missing, missing}
                                                  : pressure.front();
                    const double level = soundPressureLevelDecibels(axisPressure);
                    on.points.emplace_back(frequency.frequency, level);
                    values.frequencies.push_back(frequency.frequency);
                    values.levels.push_back(level);
                    values.pressure.push_back(axisPressure);
                    const auto count = angularSampleCount(*project.simulation.request);
                    std::array<double, 2> width{missing, missing};
                    if ((count == 37 || count == 73) && pressure.size() == 1 + 2 * count) {
                        for (std::size_t cut = 0; cut < 2; ++cut) {
                            std::vector<AcousticsAngularSample> samples;
                            for (std::size_t a = 0; a < count; ++a)
                                samples.push_back(
                                    {(count == 73 ? -180. : -90.) + 5. * a,
                                     soundPressureLevelDecibels(pressure[1 + cut * count + a])});
                            if (const auto measured = beamwidthDegrees(
                                    samples, beamwidthThreshold->currentData().toDouble()))
                                width[cut] = measured->widthDegrees;
                        }
                    }
                    h.points.emplace_back(frequency.frequency, width[0]);
                    v.points.emplace_back(frequency.frequency, width[1]);
                    std::complex<double> normalized{missing, missing};
                    for (const auto &z : excitation.impedances)
                        if (z.faceId == impedanceFace->currentData().toUInt()) {
                            const auto medium = project.simulation.request->medium;
                            normalized = normalizedVolumeImpedance(
                                z.volumeImpedance, z.area, medium.density, medium.soundSpeed);
                        }
                    real.points.emplace_back(frequency.frequency, normalized.real());
                    imaginary.points.emplace_back(frequency.frequency, normalized.imag());
                    values.impedance.push_back(normalized);
                }
                series[0].push_back(std::move(on));
                series[1].push_back(std::move(h));
                series[2].push_back(std::move(v));
                series[3].push_back(std::move(real));
                series[3].push_back(std::move(imaginary));
                compared.push_back(std::move(values));
                polarInputs.push_back({index, std::move(project)});
            } catch (const std::exception &e) {
                errors << QString("#%1: %2").arg(index + 1).arg(QString::fromUtf8(e.what()));
            }
        }
    }
    for (std::size_t i = 0; i < plots.size(); ++i)
        plots[i]->setSeries(std::move(series[i]));
    const auto &settings = baseline.simulation.settings;
    targetEnabled->setToolTip(QString("H %1° · V %2° · ±%3°; gemeinsame Simulationseinstellungen.")
                                  .arg(settings.targetHorizontalDegrees)
                                  .arg(settings.targetVerticalDegrees)
                                  .arg(settings.targetToleranceDegrees));
    for (std::size_t i = 1; i <= 2; ++i) {
        if (targetEnabled->isChecked()) {
            const double target =
                i == 1 ? settings.targetHorizontalDegrees : settings.targetVerticalDegrees;
            plots[i]->setReferenceBand(
                PlotReferenceBand{std::max(0., target - settings.targetToleranceDegrees),
                                  std::min(360., target + settings.targetToleranceDegrees),
                                  QColor(90, 135, 90, 55), "Ziel"});
        } else
            plots[i]->setReferenceBand(std::nullopt);
    }
    QStringList differences;
    if (compared.size() > 1) {
        const auto &reference = compared.front();
        for (std::size_t i = 1; i < compared.size(); ++i) {
            const auto &value = compared[i];
            if (value.frequencies != reference.frequencies)
                continue;
            double delta = 0;
            bool haveLevel = false;
            long double pNumerator = 0, pDenominator = 0, zNumerator = 0, zDenominator = 0;
            for (std::size_t f = 0; f < value.frequencies.size(); ++f) {
                if (std::isfinite(value.levels[f]) && std::isfinite(reference.levels[f])) {
                    delta = std::max(delta, std::abs(value.levels[f] - reference.levels[f]));
                    haveLevel = true;
                }
                if (std::isfinite(std::abs(value.pressure[f])) &&
                    std::isfinite(std::abs(reference.pressure[f]))) {
                    pNumerator += std::norm(value.pressure[f] - reference.pressure[f]);
                    pDenominator += std::norm(reference.pressure[f]);
                }
                if (std::isfinite(std::abs(value.impedance[f])) &&
                    std::isfinite(std::abs(reference.impedance[f]))) {
                    zNumerator += std::norm(value.impedance[f] - reference.impedance[f]);
                    zDenominator += std::norm(reference.impedance[f]);
                }
            }
            const auto relative = [](long double a, long double b) {
                return b > 0
                           ? QString::number(static_cast<double>(std::sqrt(a / b)) * 100, 'g', 4) +
                                 " %"
                           : QString("—");
            };
            differences << QString("#%1 / #%2: |ΔSPL|max %3 dB · Δp %4 · ΔZ %5")
                               .arg(value.row + 1)
                               .arg(reference.row + 1)
                               .arg(haveLevel ? QString::number(delta, 'g', 4) : QString("—"))
                               .arg(relative(pNumerator, pDenominator))
                               .arg(relative(zNumerator, zDenominator));
        }
    }
    difference->setText(differences.join("   "));
    if (!errors.empty()) {
        status->setText(errors.join("; ").left(180));
        status->setToolTip(errors.join("\n"));
    }
    updatePolarFrequencies();
    updateMapVariants();
    updateControls();
}

void ParameterSweepDialog::Implementation::analyzeVariant() {
    if (busy || !study || polarInputs.size() != 1 || !dialog->onVariantAnalysisRequested)
        return;
    const auto row = polarInputs.front().row;
    const QFileInfo file(QDir(study->directory).filePath(study->variants[row].projectFile));
    const QString directoryPath = QFileInfo(study->directory).canonicalFilePath();
    if (file.isSymLink() || file.canonicalFilePath().isEmpty() ||
        file.canonicalPath() != directoryPath) {
        QMessageBox::critical(dialog, "Variante analysieren", "Die Variantendatei liegt nicht im Studienverzeichnis.");
        return;
    }
    try {
        const auto project = readProject(file.canonicalFilePath());
        if (!sweepVariantMatches(*study, row, project) || !completeVariantResult(*study, project))
            throw std::invalid_argument("Ein vollständiges gültiges Variantenergebnis fehlt.");
        dialog->onVariantAnalysisRequested(file.canonicalFilePath());
        dialog->accept();
    } catch (const std::exception &e) {
        QMessageBox::critical(dialog, "Variante analysieren", QString::fromUtf8(e.what()));
    }
}

void ParameterSweepDialog::Implementation::updatePolarFrequencies() {
    const auto previous = polarFrequency->currentData();
    std::vector<double> frequencies;
    for (const auto &input : polarInputs)
        for (const auto &frequency : input.project.simulation.result->frequencies)
            frequencies.push_back(frequency.frequency);
    std::sort(frequencies.begin(), frequencies.end());
    frequencies.erase(std::unique(frequencies.begin(), frequencies.end()), frequencies.end());
    {
        const QSignalBlocker blocker(polarFrequency);
        polarFrequency->clear();
        for (const double frequency : frequencies)
            polarFrequency->addItem(QString("%1 Hz").arg(frequency, 0, 'g', 8), frequency);
        const int old = polarFrequency->findData(previous);
        if (old >= 0)
            polarFrequency->setCurrentIndex(old);
    }
    updatePolars();
}

void ParameterSweepDialog::Implementation::updatePolars() {
    std::array<std::vector<PolarSeries>, 2> polarSeries;
    QStringList missing, legend;
    const double selected = polarFrequency->currentData().toDouble();
    for (const auto &input : polarInputs) {
        const auto &project = input.project;
        const auto count = angularSampleCount(*project.simulation.request);
        const auto &frequencies = project.simulation.result->frequencies;
        const auto found = std::find_if(frequencies.begin(), frequencies.end(),
                                        [selected](const auto &f) { return f.frequency == selected; });
        if ((count != 37 && count != 73) || found == frequencies.end() ||
            found->excitations.empty() ||
            found->excitations.front().observerPressure.size() != 1 + 2 * count) {
            missing << QString("#%1").arg(input.row + 1);
            continue;
        }
        const auto &pressure = found->excitations.front().observerPressure;
        const QColor color = QColor::fromRgba(colors[input.row % colors.size()]);
        const QString label = QString("#%1 · %2").arg(input.row + 1)
                                  .arg(sweepVariantLabel(*study, input.row));
        bool valid = true;
        std::array<PolarSeries, 2> pair{PolarSeries{label, color, {}}, PolarSeries{label, color, {}}};
        for (std::size_t cut = 0; cut < 2; ++cut) {
            const auto levels = normalizedCut(pressure, count, cut,
                                              polarNormalization->currentData().toInt() == 1);
            if (std::all_of(levels.begin(), levels.end(), [](double x) { return std::isnan(x); })) {
                valid = false;
                break;
            }
            for (std::size_t i = 0; i < count; ++i) {
                pair[cut].levels.push_back({((count == 73 ? -180. : -90.) + 5. * i) *
                                                std::numbers::pi / 180., levels[i]});
            }
        }
        if (!valid) {
            missing << QString("#%1").arg(input.row + 1);
            continue;
        }
        for (std::size_t cut = 0; cut < 2; ++cut)
            polarSeries[cut].push_back(std::move(pair[cut]));
        legend << QString("<span style=\"color:%1\">●</span> %2")
                      .arg(color.name(), label.toHtmlEscaped());
    }
    polarH->setSeries(std::move(polarSeries[0]));
    polarV->setSeries(std::move(polarSeries[1]));
    polarLegend->setText(legend.join(" &nbsp; "));
    polarStatus->setText(missing.isEmpty()
        ? QString("H/V aus gespeicherten Beobachtern; keine Winkel ergänzt.")
        : QString("Für %1 fehlen gültige H/V-Daten bei dieser Frequenz oder ein gültiger Bezug.")
              .arg(missing.join(", ")));
}

void ParameterSweepDialog::Implementation::applyStudyScale() {
    const auto &settings = baseline.simulation.settings;
    const PlotDbScale scale{settings.directivityMinimumDb, settings.directivityMaximumDb,
                            settings.directivityStepDb};
    polarH->setDecibelScale(scale);
    polarV->setDecibelScale(scale);
    mapH->setDecibelScale(scale);
    mapV->setDecibelScale(scale);
}

void ParameterSweepDialog::Implementation::updateMapVariants() {
    const auto previous = mapVariant->currentData();
    {
        const QSignalBlocker blocker(mapVariant);
        mapVariant->clear();
        mapVariant->addItem("Variante auswählen …", QVariant{});
        for (const auto &input : polarInputs)
            mapVariant->addItem(QString("#%1 · %2").arg(input.row + 1)
                                    .arg(sweepVariantLabel(*study, input.row)),
                                QVariant::fromValue(static_cast<qulonglong>(input.row)));
        const int selected = previous.isValid() ? mapVariant->findData(previous) : 0;
        mapVariant->setCurrentIndex(std::max(0, selected));
    }
    updateMaps();
}

void ParameterSweepDialog::Implementation::updateMaps() {
    mapH->setLevels({}, {}, {});
    mapV->setLevels({}, {}, {});
    if (!mapVariant->currentData().isValid()) {
        mapStatus->setText("Variante für die Polarmaps auswählen.");
        return;
    }
    const auto row = static_cast<std::size_t>(mapVariant->currentData().toULongLong());
    const auto input = std::find_if(polarInputs.begin(), polarInputs.end(),
                                    [row](const auto &value) { return value.row == row; });
    if (input == polarInputs.end()) {
        mapStatus->setText("Die ausgewählte Variante ist nicht verfügbar.");
        return;
    }
    const auto &project = input->project;
    const auto count = angularSampleCount(*project.simulation.request);
    if (count != 37 && count != 73) {
        mapStatus->setText("Kein gültiges gespeichertes H/V-Winkelraster.");
        return;
    }
    std::vector<double> angles, frequencies;
    std::array<std::vector<double>, 2> levels;
    for (std::size_t i = 0; i < count; ++i)
        angles.push_back((count == 73 ? -180. : -90.) + 5. * i);
    const double missing = std::numeric_limits<double>::quiet_NaN();
    bool gaps = false;
    for (const auto &frequency : project.simulation.result->frequencies) {
        frequencies.push_back(frequency.frequency);
        for (std::size_t cut = 0; cut < 2; ++cut) {
            const auto rowOffset = levels[cut].size();
            levels[cut].resize(rowOffset + count, missing);
            if (frequency.excitations.empty() ||
                frequency.excitations.front().observerPressure.size() != 1 + 2 * count) {
                gaps = true;
                continue;
            }
            const auto values = normalizedCut(frequency.excitations.front().observerPressure,
                                              count, cut,
                                              mapNormalization->currentData().toInt() == 1);
            gaps = gaps || std::any_of(values.begin(), values.end(),
                                       [](double value) { return std::isnan(value); });
            std::copy(values.begin(), values.end(), levels[cut].begin() + rowOffset);
        }
    }
    if (frequencies.empty()) {
        mapStatus->setText("Keine gespeicherten Frequenzdaten.");
        return;
    }
    mapH->setLevels(angles, frequencies, std::move(levels[0]));
    mapV->setLevels(std::move(angles), std::move(frequencies), std::move(levels[1]));
    mapStatus->setText(gaps ? "Fehlende oder ungültige Werte bleiben als Lücken erhalten."
                            : "Gespeicherte H/V-Schnitte; Farbskala wie bei den Polars.");
}

ParameterSweepDialog::ParameterSweepDialog(Project baseline, QWidget *parent)
    : QDialog(parent),
      implementation_(std::make_unique<Implementation>(this, std::move(baseline))) {}
ParameterSweepDialog::~ParameterSweepDialog() {
    // The model and control connections must stop referencing the implementation before it dies.
    for (auto *child : findChildren<QObject *>())
        QObject::disconnect(child, nullptr, this, nullptr);
    implementation_->model->setStudy(nullptr);
}
void ParameterSweepDialog::reject() {
    if (implementation_->busy) {
        if (implementation_->waitingForGpu) {
            implementation_->cancelWaiting();
            QDialog::reject();
            return;
        }
        implementation_->closeWhenIdle = true;
        implementation_->worker.request_stop();
        implementation_->cancel->setEnabled(false);
        implementation_->status->setText("Abbruch angefordert …");
        return;
    }
    QDialog::reject();
}
void ParameterSweepDialog::openStudy(const QString &path) {
    implementation_->load(path);
}
void ParameterSweepDialog::configureConvergence(std::vector<double> frequencies) {
    implementation_->configure(std::move(frequencies));
}
} // namespace hornsculpt
