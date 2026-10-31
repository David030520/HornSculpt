#pragma once

#include "project_file.hpp"
#include <QDockWidget>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <hornsculpt/bem.hpp>
#include <hornsculpt/plot_scale.hpp>
#include <hornsculpt/pressure_field_plane.hpp>
#include <map>
#include <memory>
#include <mutex>
#include <thread>

class QComboBox;
class QCheckBox;
class QDoubleSpinBox;
class QLabel;
class QLineEdit;
class QMainWindow;
class QProgressBar;
class QPushButton;
class QSpinBox;
class QTabWidget;
class QTimer;
class QToolButton;

namespace hornsculpt {
class HeatmapWidget;
class PlotWidget;
class PolarPlotWidget;
class AngularResponseWidget;
class PhaseResponseWidget;
class AcousticFieldWidget;

class SimulationDock final : public QDockWidget {
  public:
    explicit SimulationDock(QMainWindow *parent);
    ~SimulationDock() override;
    void setMesh(std::shared_ptr<const SurfaceMesh> mesh, std::uint64_t revision,
                 std::vector<FaceVelocity> velocities = {});
    void setGeometrySource(std::shared_ptr<const DraftGeometry> draft,
                           std::shared_ptr<const CadSnapshot> cad, std::uint64_t revision,
                           std::map<std::uint32_t, double> faceTargets = {},
                           std::set<std::uint32_t> reversedFaces = {});
    bool usesAxisymmetricGeometry() const noexcept;
    void invalidate(std::uint64_t revision);
    // A new geometry source may already have invalidated the result in setMesh.
    void setAcousticAxis(Vec3 origin, Vec3 axis, bool geometryAlreadyInvalidated = false);
    void setObservationFrame(std::optional<AcousticFrame> frame);
    AcousticFrame observationFrame() const;
    void setCoupling(lem::Network network, bool enabled);
    void setInfiniteBaffle(std::optional<InfiniteBaffle> baffle);
    void setAcousticRegions(std::vector<AcousticRegionDefinition> regions);
    void setPartitionMeshingOptions(MeshingOptions options);
    std::function<void()> onAcousticRegionsChanged;
    void setFrequencySweep(double minimum, double maximum, std::size_t count);
    void requestSolve();
    void cancelSolve();
    void setRealtimeEnabled(bool enabled);
    bool realtimeEnabled() const noexcept;
    bool isBusy() const noexcept;
    // GUI-thread identity of the currently active solve; zero means no active job.
    std::uint64_t activeSolveJobId() const noexcept {
        return busy_ ? activeJob_ : 0;
    }
    bool hasCurrentResult() const noexcept;
    bool hasConfiguredExcitation() const;
    std::uint64_t inputGeneration() const noexcept {
        return generation_;
    }
    QString lastError() const;
    void requestFieldEvaluation();
    bool fieldEvaluationInProgress() const noexcept;
    bool fieldReady() const noexcept;
    QString lastFieldError() const;
    void revealResultView();
    void beginExternalJob();
    bool externalJobReady() const noexcept;
    void endExternalJob();
    std::shared_ptr<const BemResult> result() const noexcept;
    std::uint64_t resultRevision() const noexcept;
    SimulationArchive archive() const;
    void validateArchive(const SimulationArchive &archive) const;
    void restoreArchive(const SimulationArchive &archive);
    std::function<void()> onChanged;
    std::function<void(bool)> onRealtimeChanged;
    std::function<void()> onObservationFrameChanged;
    std::function<void(std::optional<lem::Solution>, bool)> onLemResultChanged;
    std::function<void(std::shared_ptr<const PressureFieldPlane>, PressureFieldAppearance)>
        onFieldChanged;
    std::function<void()> onFocusField;
    std::function<void(bool)> onFieldProbeModeChanged;
    std::function<void()> onRequestConvergence;

  private:
    struct SolveJob;
    struct GeometrySource;
    void runWorker(std::stop_token shutdown);
    void controlsChanged();
    void updateDomainControls();
    void updateResultControls();
    void updateResultQuality();
    bool hasExcitation() const;
    bool hasCurrentGeometry() const;
    void scheduleRealtime();
    void cancelPendingJobs();
    void presentResult(bool current);
    struct DisplayFrequency {
        double onAxisDb{};
        std::map<std::uint32_t, std::complex<double>> impedances;
        std::vector<double> horizontalDb;
        std::vector<double> verticalDb;
    };
    // Presentation only: new anchors drive old rows until exact replacements arrive.
    // Project archives always use result_/resultRequest_ from one solver snapshot.
    std::map<double, DisplayFrequency> displayFrequencies_;
    std::map<double, DisplayFrequency> displayBaseline_;
    std::vector<MeasuredResponse> measurements_;
    void updateDisplayEstimate();
    void updateFrequencyPlots();
    void updateAngularPlots();
    void updateAngularFrequencyPlots(bool current);
    void updateFieldAndPhase(bool current);
    void updateCta2034(bool current);
    void updateHeatmaps();
    void updateDirectivityScale();
    void applyDirectivityScale(PlotDbScale scale);
    void importMeasurement();
    void exportResults();
    void setStatus(QString text, bool error = false);
    BemRequest makeRequest() const;
    std::shared_ptr<const SurfaceMesh> mesh_;
    std::vector<FaceVelocity> velocities_;
    std::shared_ptr<const BemResult> result_;
    std::shared_ptr<const SurfaceMesh> resultMesh_;
    std::shared_ptr<const SurfaceMesh> resultSourceMesh_;
    std::vector<AcousticRegionDefinition> resultRegions_;
    std::optional<BemRequest> resultRequest_;
    std::optional<AxisymmetricSamplingSettings> resultAxisymmetricSampling_;
    std::shared_ptr<const GeometrySource> geometrySource_;
    bool geometrySourceCurrent_{};
    Vec3 origin_{};
    Vec3 axis_{0, 0, 1};
    std::optional<AcousticFrame> observationFrame_;
    lem::Network couplingNetwork_;
    bool couplingEnabled_{};
    std::optional<InfiniteBaffle> infiniteBaffle_;
    std::vector<AcousticRegionDefinition> acousticRegions_;
    MeshingOptions partitionOptions_;
    FastMultipoleSettings fastMultipoleSettings_;
    std::uint64_t currentRevision_{};
    std::uint64_t meshRevision_{};
    std::uint64_t resultRevision_{};
    std::uint64_t generation_{};
    std::uint64_t resultGeneration_{};
    std::uint64_t activeJob_{};
    std::jthread worker_;
    std::stop_source jobStop_;
    std::mutex workerMutex_;
    std::condition_variable_any workerWake_;
    std::shared_ptr<const SolveJob> pendingJob_;
    bool busy_{};
    bool externalJob_{};
    bool releaseGpuRequested_{}; // Protected by workerMutex_.
    std::atomic_bool gpuReleased_{true};
    bool meshCurrent_{};
    bool resultComplete_{};
    bool resultCompletenessKnown_{};
    bool realtimeEnabled_{};
    bool realtimePending_{};
    QString error_;
    QDoubleSpinBox *minimumFrequency_{};
    QDoubleSpinBox *maximumFrequency_{};
    QSpinBox *frequencyCount_{};
    QCheckBox *fastMultipole_{};
    QSpinBox *multipoleLeafCapacity_{};
    QComboBox *multipoleJacobiBlock_{};
    QSpinBox *nearMaximumDepth_{};
    QComboBox *geometryMode_{};
    QDoubleSpinBox *meridianLength_{};
    QDoubleSpinBox *meridianTolerance_{};
    QDoubleSpinBox *rotationalTolerance_{};
    QSpinBox *angularSectors_{};
    QSpinBox *singularOrder_{};
    QDoubleSpinBox *density_{};
    QDoubleSpinBox *soundSpeed_{};
    QDoubleSpinBox *velocity_{};
    QDoubleSpinBox *phase_{};
    QComboBox *fieldMode_{};
    QToolButton *quality_{};
    QDoubleSpinBox *observationDistance_{};
    QComboBox *frequencySelector_{};
    QComboBox *impedanceSelector_{};
    QComboBox *normalization_{};
    QDoubleSpinBox *normalizationAngle_{};
    QLineEdit *directivityMinimum_{};
    QLineEdit *directivityMaximum_{};
    QLineEdit *directivityStep_{};
    PlotDbScale directivityScale_;
    QComboBox *beamwidthThreshold_{};
    QCheckBox *beamwidthTarget_{};
    QDoubleSpinBox *targetH_{};
    QDoubleSpinBox *targetV_{};
    QDoubleSpinBox *targetTolerance_{};
    QPushButton *solve_{};
    QPushButton *baffle_{};
    QPushButton *stop_{};
    QPushButton *realtime_{};
    QTimer *realtimeTimer_{};
    QLabel *status_{};
    QLabel *excitation_{};
    QProgressBar *progress_{};
    PlotWidget *onAxis_{};
    PlotWidget *impedance_{};
    PlotWidget *ctaLevels_{};
    PlotWidget *ctaDi_{};
    QLabel *ctaStatus_{};
    QComboBox *ctaExcitation_{};
    int ctaSelection_{-1};
    PlotWidget *beamwidthH_{};
    PlotWidget *beamwidthV_{};
    PolarPlotWidget *polarH_{};
    PolarPlotWidget *polarV_{};
    HeatmapWidget *isobarH_{};
    HeatmapWidget *isobarV_{};
    HeatmapWidget *radialIsobarH_{};
    HeatmapWidget *radialIsobarV_{};
    AngularResponseWidget *angularResponses_{};
    PhaseResponseWidget *phaseResponses_{};
    AcousticFieldWidget *acousticField_{};
    QTabWidget *workspaceTabs_{};
    QWidget *resultPage_{};
    QWidget *solverCommands_{};
    QWidget *resultSelectorBar_{};
    QLabel *emptyResults_{};
    QWidget *axisymmetricControls_{};
    QWidget *frequencyControls_{};
    QWidget *impedanceControls_{};
    QWidget *normalizationControls_{};
    QWidget *directivityControls_{};
    QTabWidget *resultTabs_{};
};

} // namespace hornsculpt
