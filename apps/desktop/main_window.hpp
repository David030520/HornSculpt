#pragma once

#include "edit_history.hpp"
#include "geometry_mesh_worker.hpp"
#include "mesh_settings.hpp"
#include "project_file.hpp"
#include "project_save_worker.hpp"
#include "recovery_store.hpp"
#include <QMainWindow>
#include <cstdint>
#include <functional>
#include <hornsculpt/freecad_link.hpp>
#include <hornsculpt/geometry_analysis.hpp>
#include <hornsculpt/meshing.hpp>
#include <map>
#include <memory>
#include <set>
#include <thread>

class QAction;
class QCheckBox;
class QComboBox;
class QDockWidget;
class QDialog;
class QLabel;
class QTimer;
class QTabWidget;
class QScrollArea;
class QVBoxLayout;
class QTreeWidgetItem;
class QToolButton;
class QUndoStack;

namespace hornsculpt {
class MeshViewport;
class MeshSectionView;
class ParameterPanel;
class MeshPanel;
class PlotWidget;
class FreeCadDock;
class SimulationDock;
class LemEditor;

class MainWindow final : public QMainWindow {
  public:
    // The starter must return a joinable worker, or throw without starting the job.
    using SourceImportStarter = std::function<std::jthread(std::function<void(std::stop_token)>)>;
    explicit MainWindow(QWidget *parent = nullptr, ProjectSaveWorker::Writer writer = {},
                        SourceImportStarter sourceImportStarter = {}, QString recoveryRoot = {});
    bool saveInProgress() const noexcept;
    ~MainWindow() override;
    void openProject(const QString &path);
    void openStep(const QString &path);
    bool sourceImportInProgress() const noexcept;
    QString sourceImportError() const;
    bool geometryReady() const noexcept;
    std::uint64_t revision() const noexcept;
    void requestMesh();
    bool meshReady() const noexcept;
    bool meshInProgress() const noexcept;
    QString lastMeshError() const;
    void exerciseProjectRoundTrip();
    void exerciseResultProjectRoundTrip();
    void exerciseLoadedProjectRoundTrip();
    void requestSolve();
    bool simulationReady() const;
    bool simulationInProgress() const;
    QString lastSimulationError() const;
    void setSimulationFrequencySweep(double minimum, double maximum, int count);
    void setSimulationView(const QString &view);
    void requestFieldEvaluation();
    bool fieldEvaluationInProgress() const noexcept;
    bool fieldReady() const noexcept;
    QString lastFieldError() const;

  protected:
    void closeEvent(QCloseEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;
    void setCompactToolsMode(bool compact);
    bool eventFilter(QObject *object, QEvent *event) override;

  private:
    struct PreparedProject {
        Project project;
        std::shared_ptr<const DraftGeometry> geometry;
        GeometryAnalysis analysis;
        std::shared_ptr<const SurfaceMesh> rawPreview;
        PreparedCadSnapshot cad{PreparedCadSnapshot::prepare({})};
        bool rebuiltDraft{};
    };
    struct SourceImport {
        QString path;
        bool step{};
        std::uint64_t id{}, epoch{}, revision{}, changeSerial{};
        MeshSettings meshing;
        std::shared_ptr<const Project> appendToDraft;
        bool openAsCopy{};
    };
    static PreparedProject prepareProject(Project project, bool verifyDraft,
                                          std::stop_token stop = {});
    void applyProject(const PreparedProject &prepared);
    void adoptProject(const PreparedProject &prepared, const QString &path);
    void queueSourceImport(const QString &path, bool step, bool openAsCopy = false);
    void startSourceImport(SourceImport request);
    std::jthread sourceImportWorker_;
    SourceImportStarter sourceImportStarter_;
    bool sourceImportBusy_{};
    QString sourceImportError_;
    std::uint64_t sourceImportId_{};
    std::optional<SourceImport> pendingSourceImport_;
    void parametersChanged();
    void startGeometry();
    void acceptGeometry(std::uint64_t revision, std::shared_ptr<const DraftGeometry> geometry,
                        std::shared_ptr<const SurfaceMesh> preview, QString error);
    void finishGeometry(std::uint64_t revision, GeometryAnalysis analysis, QString error);
    void updateDisplay();
    std::vector<std::uint32_t> faceColorSelection() const;
    FaceColors &activeFaceColors();
    void updateFaceColorActions();
    void chooseFaceColor();
    void setSelectedFaceColor(std::optional<FaceColor> color);
    void exerciseFaceColorEditing();
    void updatePlots(const GeometryAnalysis &analysis);
    void updateTitle();
    void markDirty();
    void updateMeshStatus();
    void updateSourcePresentation();
    Project projectState() const;
    void applyProject(const Project &project);
    void loadProject(const QString &path);
    void setDiagnostic(const QString &message, bool error = false);
    void saveProject(bool choosePath, std::function<void()> continuation = {});
    void confirmUnsaved(std::function<void()> continuation);
    void setSaveActionsEnabled(bool enabled);
    void newProject();
    void openDraftWizard();
    void startGuidedCalculation();
    void startRecovery();
    void offerRecovery();
    void saveRecovery();
    void clearRecovery(bool newSession = true);
    void updateWorkflow();
    void exportDraftStep();
    void editDraftAssembly();
    void addCadToDraft();
    std::vector<std::uint32_t> placedSelection(const std::vector<std::uint32_t> &ids) const;
    void sourceModeChanged();
    void acceptCadSnapshot(std::shared_ptr<const CadSnapshot> snapshot);
    void cadAvailabilityChanged(bool available);
    void meshingChanged();
    void updateMeshFaces();
    void setFaceTargets(const std::vector<std::uint32_t> &ids, std::optional<double> target);
    void reverseFaces(const std::vector<std::uint32_t> &ids);
    void updateOrientedPreview();
    std::set<std::uint32_t> &activeReversedFaces();
    void validateFaceTargets(const MeshSettings &settings) const;
    std::map<std::uint32_t, double> &activeFaceTargets();
    void configureSimulation();
    void configureLemFaces();
    void updateInterfacePresentation();
    void updateAcousticDomainDisplay(const SimulationArchive &archive);
    void couplingChanged();
    bool cadMode() const;
    EditState editState() const;
    void rememberEdit(bool immediate = false);
    void flushEdit();
    void resetEditHistory();
    void applyEditState(const EditState &state);
    void undoEdit();
    void redoEdit();
    void updateEditActions();
    void openParameterStudy(bool convergence = false, const QString &manifestPath = {});
    void adoptStudyVariant(Project variant);
    ParameterPanel *parameters_{};
    MeshPanel *meshPanel_{};
    MeshSettings meshing_;
    DraftAssembly assembly_;
    std::map<std::uint32_t, double> draftFaceTargetEdgeLengths_;
    std::map<std::uint32_t, double> cadFaceTargetEdgeLengths_;
    FaceColors draftFaceColors_;
    FaceColors cadFaceColors_;
    std::set<std::uint32_t> draftReversedFaceIds_;
    std::set<std::uint32_t> cadReversedFaceIds_;
    MeshViewport *viewport_{};
    PlotWidget *volumePlot_{};
    PlotWidget *areaPlot_{};
    PlotWidget *derivativePlot_{};
    PlotWidget *anglePlot_{};
    MeshSectionView *horizontalSection_{};
    MeshSectionView *verticalSection_{};
    QTabWidget *geometryTabs_{};
    QCheckBox *sectionMesh_{};
    QCheckBox *showAperture_{};
    QCheckBox *showInterfaces_{};
    QComboBox *representation_{};
    QComboBox *meshColor_{};
    QComboBox *sourceMode_{};
    QLabel *revisionLabel_{};
    QLabel *qualityLabel_{};
    QLabel *diagnostic_{};
    QTreeWidgetItem *meshItem_{};
    QTreeWidgetItem *sourceItem_{};
    QTreeWidgetItem *facesItem_{};
    QTreeWidgetItem *projectItem_{};
    QTreeWidgetItem *domainsItem_{};
    std::shared_ptr<const SurfaceMesh> interfaceOverlaySource_;
    std::shared_ptr<const SurfaceMesh> interfaceOverlay_;
    std::set<std::uint32_t> interfaceOverlayFaceIds_;
    std::vector<AcousticRegionDefinition> displayedAcousticRegions_;
    bool interfaceOverlayFromResult_{};
    QDockWidget *analysisDock_{};
    QDockWidget *parameterDock_{};
    QDockWidget *meshSettingsDock_{};
    FreeCadDock *cadDock_{};
    SimulationDock *simulationDock_{};
    LemEditor *lemEditor_{};
    QTimer *debounce_{};
    QTimer *realtimeMeshTimer_{};
    std::uint64_t realtimeMeshAttempt_{};
    QAction *meshAction_{};
    QAction *draftWizardAction_{};
    QToolButton *workflowMesh_{};
    QToolButton *workflowSimulation_{};
    QToolButton *compactTools_{};
    QDialog *compactToolsDialog_{};
    QScrollArea *workflowScroll_{};
    QScrollArea *viewScroll_{};
    QVBoxLayout *centerLayout_{};
    bool compactHeight_{};
    std::uint64_t dockResizeGeneration_{};
    QLabel *workflowStatus_{};
    bool initialModelFitPending_{true};
    QAction *stepExportAction_{};
    QAction *cancelAction_{};
    QAction *assignFaceColor_{};
    QAction *resetFaceColors_{};
    QAction *fieldVisibility_{};
    bool fieldCameraInitialized_{};
    GeometryMeshWorker geometryWorker_;
    GeometryMeshWorker meshWorker_;
    std::shared_ptr<const DraftGeometry> geometry_;
    std::shared_ptr<const SurfaceMesh> preview_;
    std::shared_ptr<const SurfaceMesh> rawPreview_;
    std::shared_ptr<const MeshingResult> mesh_;
    std::shared_ptr<const CadSnapshot> cadSnapshot_;
    std::optional<MeshImportSelection> gmsh_;
    std::uint64_t revision_{1};
    std::uint64_t geometryRevision_{};
    std::uint64_t analysisRevision_{};
    std::uint64_t geometryJobRevision_{};
    std::uint64_t meshRevision_{};
    std::uint64_t meshJobId_{};
    std::uint64_t projectEpoch_{};
    std::uint64_t changeSerial_{};
    bool geometryBusy_{};
    bool geometryPending_{};
    bool meshBusy_{};
    bool dirty_{};
    bool closing_{};
    bool closeApproved_{};
    bool loading_{};
    bool meshSettingsValid_{true};
    QUndoStack *undoStack_{};
    QTimer *editTimer_{};
    QAction *undoAction_{};
    QAction *redoAction_{};
    EditState historyState_;
    std::optional<EditState> pendingEdit_;
    bool historyReady_{};
    bool historyReplaying_{};
    std::unique_ptr<ProjectSaveWorker> saveWorker_;
    std::unique_ptr<RecoveryStore> recoveryStore_;
    std::shared_ptr<RecoveryStore::Entry> recoveryOpening_;
    std::shared_ptr<RecoveryStore::Entry> recoveredEntry_;
    QTimer *recoveryTimer_{};
    QString recoveryRoot_;
    std::uint64_t recoveryGeneration_{}, recoverySavedSerial_{};
    std::weak_ptr<const BemResult> recoverySavedResult_;
    std::optional<bool> recoverySavedComplete_;
    bool recoveryHasSaved_{};
    std::vector<QAction *> saveBlockedActions_;
    std::function<void()> afterSave_;
    QString projectPath_;
    std::vector<NamedProjectSnapshot> comparisonSnapshots_;
    QString meshError_;
};
} // namespace hornsculpt
