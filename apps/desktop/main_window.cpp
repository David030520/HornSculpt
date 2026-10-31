#include "main_window.hpp"
#include "draft_assembly_dialog.hpp"
#include "draft_wizard.hpp"
#include "freecad_dock.hpp"
#include "galerkin/exterior_domain.hpp"
#include "gmsh_source.hpp"
#include "lem_editor.hpp"
#include "mesh_import_dialog.hpp"
#include "mesh_panel.hpp"
#include "parameter_panel.hpp"
#include "parameter_sweep.hpp"
#include "project_assembly.hpp"
#include "result_comparison.hpp"
#include "simulation_dock.hpp"
#include "simulation_request.hpp"
#include "step_export_dialog.hpp"
#include "step_import.hpp"
#include <QAbstractItemDelegate>
#include <QAbstractItemView>
#include <QAbstractSpinBox>
#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QCloseEvent>
#include <QColorDialog>
#include <QComboBox>
#include <QCoreApplication>
#include <QDialog>
#include <QDir>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QGridLayout>
#include <QHash>
#include <QInputDialog>
#include <QKeyEvent>
#include <QLabel>
#include <QLineEdit>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QPointer>
#include <QProgressDialog>
#include <QPushButton>
#include <QRegularExpression>
#include <QResizeEvent>
#include <QScopedValueRollback>
#include <QScreen>
#include <QScrollArea>
#include <QScrollBar>
#include <QSettings>
#include <QSignalBlocker>
#include <QSplitter>
#include <QStandardPaths>
#include <QStatusBar>
#include <QStyle>
#include <QStyleOptionViewItem>
#include <QTabWidget>
#include <QTemporaryDir>
#include <QTimer>
#include <QToolBar>
#include <QToolButton>
#include <QTreeWidget>
#include <QUndoStack>
#include <QVBoxLayout>
#include <QWhatsThis>
#include <algorithm>
#include <cmath>
#include <complex>
#include <hornsculpt/cad_geometry.hpp>
#include <hornsculpt/mesh_orientation.hpp>
#include <hornsculpt/mesh_section_view.hpp>
#include <hornsculpt/mesh_viewport.hpp>
#include <hornsculpt/plot_scale.hpp>
#include <hornsculpt/plot_widget.hpp>
#include <iterator>
#include <limits>
#include <numbers>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>

namespace hornsculpt {
namespace {
class GuidedDraftCommand final : public QUndoCommand {
  public:
    explicit GuidedDraftCommand(std::function<void(bool)> apply) : apply_(std::move(apply)) {
        setText(QObject::tr("Horn-Assistent"));
    }
    void undo() override {
        apply_(false);
    }
    void redo() override {
        if (initial_)
            initial_ = false;
        else
            apply_(true);
    }

  private:
    std::function<void(bool)> apply_;
    bool initial_{true};
};
QDockWidget *dock(QMainWindow *window, const QString &title, const QString &id, QWidget *content,
                  Qt::DockWidgetArea area) {
    auto *result = new QDockWidget(title, window);
    result->setObjectName(id);
    result->setWidget(content);
    window->addDockWidget(area, result);
    return result;
}
QScrollArea *scrollableDockContent(QWidget *content, const QString &name, int minimumWidth = 300) {
    // A tabified QMainWindow dock stack combines even inactive minimum hints.
    // Preserve the actual content constraints inside a viewport, so an inactive
    // graph/editor cannot reserve half the window while a short field form is active.
    auto *scroll = new QScrollArea;
    scroll->setObjectName(name);
    scroll->setWidgetResizable(true);
    scroll->setFrameShape(QFrame::NoFrame);
    scroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    scroll->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    scroll->setMinimumSize(minimumWidth, 160);
    scroll->setWidget(content);
    return scroll;
}
void commitFocusedEdit() {
    QPointer<QWidget> focus = QApplication::focusWidget();
    if (focus) {
        auto *spin = qobject_cast<QAbstractSpinBox *>(focus);
        if (!spin)
            spin = qobject_cast<QAbstractSpinBox *>(focus->parentWidget());
        if (spin)
            spin->interpretText();
    }
    // A delegate owns the model commit: interpreting its spin box alone only
    // changes the editor. The model callback may rebuild and delete that editor.
    for (auto *editor = focus.data(); editor; editor = editor->parentWidget()) {
        auto *viewport = editor->parentWidget();
        auto *view =
            viewport ? qobject_cast<QAbstractItemView *>(viewport->parentWidget()) : nullptr;
        if (!view || view->viewport() != viewport)
            continue;
        QPointer<QWidget> guardedEditor = editor;
        QPointer<QAbstractItemDelegate> delegate = view->itemDelegateForIndex(view->currentIndex());
        if (delegate) {
            delegate->commitData(editor);
            if (guardedEditor && delegate)
                delegate->closeEditor(guardedEditor, QAbstractItemDelegate::NoHint);
        }
        return;
    }
}
} // namespace

MainWindow::MainWindow(QWidget *parent, ProjectSaveWorker::Writer writer,
                       SourceImportStarter sourceImportStarter, QString recoveryRoot)
    : QMainWindow(parent), sourceImportStarter_(std::move(sourceImportStarter)),
      saveWorker_(std::make_unique<ProjectSaveWorker>(this, std::move(writer))),
      recoveryRoot_(std::move(recoveryRoot)) {
    if (!sourceImportStarter_)
        sourceImportStarter_ = [](std::function<void(std::stop_token)> job) {
            return std::jthread(std::move(job));
        };
    setObjectName("HornSculptMainWindow");
    // Keep the same presentation for the application and embedded desktop tests.
    setStyleSheet(QStringLiteral(R"(
        QMainWindow { background: #edf0f2; color: #202b33; }
        QDockWidget::title { background: #e1e7eb; color: #263743;
                             font-weight: 600; padding: 3px 7px; }
        QToolBar { background: #edf0f2; spacing: 3px; }
        QToolBar QToolButton { padding: 2px; }
        QToolBar QToolButton:focus { border: 2px solid #256b91; padding: 0px; }
        QToolBar QToolButton[primaryAction="true"] { color: #155777; font-weight: 600; }
        #draftWorkflow QToolButton { padding: 3px 7px; color: #21485e; }
        #draftWorkflow QToolButton:hover { background: #e1edf3; }
        #draftWorkflow QToolButton:focus { border: 2px solid #256b91; padding: 1px 5px; }
        QGroupBox[fieldSection="true"] { border: 1px solid #c4d0d6; border-radius: 3px;
            margin-top: 10px; padding: 8px 5px 5px; }
        QGroupBox[fieldSection="true"]::title { subcontrol-origin: margin;
            left: 6px; padding: 0 3px; color: #304c5b; font-weight: 600; }
        QLineEdit, QAbstractSpinBox, QComboBox { background: white; color: #202b33; }
        QLineEdit:focus, QAbstractSpinBox:focus, QComboBox:focus { border: 2px solid #256b91; }
        QLineEdit:disabled, QAbstractSpinBox:disabled, QComboBox:disabled {
            background: #e8ecee; color: #697880; }
        QPushButton[primaryAction="true"] { background: #256b91; color: white;
            border: 1px solid #205c7c; border-radius: 3px; padding: 5px 10px; }
        QPushButton[primaryAction="true"]:hover { background: #1d5b7d; }
        QPushButton[primaryAction="true"]:disabled { background: #dce3e7;
            color: #68777f; border-color: #c7d0d5; }
        QPushButton[primaryAction="true"]:focus { border: 2px solid #153f58; padding: 4px 9px; }
    )"));
    resize(1480, 980);
    setMinimumSize(800, 450);
    setDockOptions(AnimatedDocks | AllowNestedDocks | AllowTabbedDocks);
    setCorner(Qt::BottomLeftCorner, Qt::LeftDockWidgetArea);
    setCorner(Qt::BottomRightCorner, Qt::RightDockWidgetArea);
    auto *file = menuBar()->addMenu(tr("&Datei"));
    auto *edit = menuBar()->addMenu(tr("&Bearbeiten"));
    auto *view = menuBar()->addMenu(tr("&Ansicht"));
    auto *geometryMenu = menuBar()->addMenu(tr("&Geometrie"));
    geometryMenu->addAction(tr("Draft-Baugruppe …"), this, [this] { editDraftAssembly(); });
    geometryMenu->addAction(tr("CAD zur Draft-Baugruppe hinzufügen …"), this,
                            [this] { addCadToDraft(); });
    auto *studyMenu = menuBar()->addMenu(tr("&Studie"));
    auto *comparisonAction = studyMenu->addAction(tr("Ergebnisstände und Bericht …"), this, [this] {
        commitFocusedEdit();
        QImage image;
        if (viewport_ && viewport_->isValid())
            image = viewport_->grabFramebuffer();
        const auto dialogEpoch = projectEpoch_;
        ResultComparisonDialog dialog(projectState(), std::move(image), this);
        dialog.onSnapshotsChanged = [this,
                                     dialogEpoch](std::vector<NamedProjectSnapshot> snapshots) {
            if (projectEpoch_ != dialogEpoch) {
                setDiagnostic(tr("Projekt wurde gewechselt · Ergebnisvergleich erneut öffnen."),
                              true);
                return;
            }
            comparisonSnapshots_ = std::move(snapshots);
            markDirty();
        };
        dialog.exec();
    });
    comparisonAction->setObjectName("openResultComparison");
    auto *sweepAction =
        studyMenu->addAction(tr("Parameter-Sweep …"), this, [this] { openParameterStudy(); });
    sweepAction->setObjectName("openParameterSweep");
    studyMenu->addAction(tr("Netzkonvergenz …"), this, [this] { openParameterStudy(true); });
    studyMenu->addAction(tr("Studie öffnen …"), this, [this] {
        const auto path = QFileDialog::getOpenFileName(this, tr("Studie öffnen"), {},
                                                       tr("HornSculpt-Studie (*.hsstudy)"));
        if (!path.isEmpty())
            openParameterStudy(false, path);
    });
    auto *toolbar = addToolBar(tr("Projekt"));
    toolbar->setObjectName("projectToolbar");
    toolbar->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    undoStack_ = new QUndoStack(this);
    undoStack_->setUndoLimit(100);
    editTimer_ = new QTimer(this);
    editTimer_->setSingleShot(true);
    editTimer_->setInterval(250);
    undoAction_ = edit->addAction(tr("Rückgängig"));
    redoAction_ = edit->addAction(tr("Wiederholen"));
    undoAction_->setObjectName(QStringLiteral("undoEdit"));
    redoAction_->setObjectName(QStringLiteral("redoEdit"));
    undoAction_->setShortcut(QKeySequence::Undo);
    redoAction_->setShortcut(QKeySequence(QStringLiteral("Ctrl+Y")));
    connect(undoAction_, &QAction::triggered, this, [this] { undoEdit(); });
    connect(redoAction_, &QAction::triggered, this, [this] { redoEdit(); });
    connect(editTimer_, &QTimer::timeout, this, [this] { flushEdit(); });
    connect(undoStack_, &QUndoStack::canUndoChanged, this, [this] { updateEditActions(); });
    connect(undoStack_, &QUndoStack::canRedoChanged, this, [this] { updateEditActions(); });
    updateEditActions();
    auto *newAction = new QAction(style()->standardIcon(QStyle::SP_FileIcon), tr("Neu"), this);
    newAction->setShortcut(QKeySequence::New);
    auto *openAction =
        new QAction(style()->standardIcon(QStyle::SP_DirOpenIcon), tr("Öffnen"), this);
    openAction->setShortcut(QKeySequence::Open);
    auto *saveAction =
        new QAction(style()->standardIcon(QStyle::SP_DialogSaveButton), tr("Speichern"), this);
    saveAction->setShortcut(QKeySequence::Save);
    file->addAction(newAction);
    file->addAction(openAction);
    auto *importStepAction = file->addAction(tr("STEP importieren …"), this, [this] {
        confirmUnsaved([this] {
            const auto path = QFileDialog::getOpenFileName(this, tr("STEP importieren"), {},
                                                           tr("STEP (*.step *.stp)"));
            if (!path.isEmpty())
                openStep(path);
        });
    });
    auto *importMeshAction = file->addAction(tr("Gmsh-Netz importieren …"), this, [this] {
        confirmUnsaved([this] {
            const auto path = QFileDialog::getOpenFileName(this, tr("Gmsh-Netz importieren"), {},
                                                           tr("Gmsh (*.msh)"));
            if (path.isEmpty())
                return;
            MeshImportDialog dialog(path, this);
            if (dialog.exec() != QDialog::Accepted || !dialog.selection())
                return;
            try {
                Project project;
                project.source = ProjectSource::freecad;
                project.meshing = meshPanel_->inputSettings();
                project.gmsh = *dialog.selection();
                project.cad = makeGmshCadSnapshot(*project.gmsh, QFileInfo(path).fileName());
                project.preview =
                    std::shared_ptr<const SurfaceMesh>(project.cad, &project.cad->mesh);
                project.geometryRevision = project.revision;
                adoptProject(prepareProject(std::move(project), false), {});
                markDirty();
                setDiagnostic(tr("Gmsh-Polygonquelle importiert · Randzuordnung und Netz prüfen."));
            } catch (const std::exception &e) {
                QMessageBox::warning(this, tr("Gmsh-Netz importieren"),
                                     QString::fromUtf8(e.what()));
            }
        });
    });
    file->addAction(saveAction);
    auto *saveAs = file->addAction(tr("Speichern unter …"));
    saveAs->setShortcut(QKeySequence::SaveAs);
    stepExportAction_ = file->addAction(tr("Draft als STEP exportieren …"));
    stepExportAction_->setObjectName("exportDraftStep");
    stepExportAction_->setToolTip(tr("Aktuellen Draft als Hornwand, Innenfläche oder Konturkurven "
                                     "mit eigener CAD-Genauigkeit exportieren."));
    stepExportAction_->setWhatsThis(stepExportAction_->toolTip());
    stepExportAction_->setEnabled(false);
    connect(stepExportAction_, &QAction::triggered, this, [this] { exportDraftStep(); });
    file->addSeparator();
    file->addAction(tr("Beenden"), QKeySequence::Quit, this, &QWidget::close);
    toolbar->addAction(newAction);
    toolbar->addAction(openAction);
    toolbar->addAction(saveAction);
    toolbar->addSeparator();
    draftWizardAction_ = new QAction(tr("Horn-Assistent …"), this);
    draftWizardAction_->setObjectName("openDraftWizard");
    draftWizardAction_->setToolTip(tr("Ein OS-SE-Grundhorn schrittweise aus Anschluss, Bauraum, "
                                      "geometrischen Öffnungswinkeln und Querschnitt erzeugen."));
    draftWizardAction_->setWhatsThis(draftWizardAction_->toolTip());
    toolbar->addAction(draftWizardAction_);
    toolbar->widgetForAction(draftWizardAction_)->setProperty("primaryAction", true);
    geometryMenu->insertAction(geometryMenu->actions().front(), draftWizardAction_);
    connect(draftWizardAction_, &QAction::triggered, this, [this] { openDraftWizard(); });
    toolbar->addSeparator();
    meshAction_ = new QAction(style()->standardIcon(QStyle::SP_FileDialogDetailedView),
                              tr("Netz erzeugen"), this);
    meshAction_->setToolTip(tr("Isotropes Netz erzeugen und auf geschlossenen Rand, Orientierung "
                               "und Elementqualität prüfen."));
    cancelAction_ =
        new QAction(style()->standardIcon(QStyle::SP_BrowserStop), tr("Netz abbrechen"), this);
    cancelAction_->setObjectName("cancelMeshing");
    cancelAction_->setEnabled(false);
    toolbar->addAction(meshAction_);
    toolbar->addAction(cancelAction_);
    geometryMenu->addAction(meshAction_);
    connect(newAction, &QAction::triggered, this, [this] { newProject(); });
    connect(openAction, &QAction::triggered, this, [this] {
        confirmUnsaved([this] {
            const auto path = QFileDialog::getOpenFileName(this, tr("Projekt öffnen"), {},
                                                           tr("HornSculpt (*.hsproj)"));
            if (!path.isEmpty())
                openProject(path);
        });
    });
    connect(saveAction, &QAction::triggered, this, [this] { saveProject(false); });
    connect(saveAs, &QAction::triggered, this, [this] { saveProject(true); });
    saveBlockedActions_ = {newAction,        openAction, importStepAction,
                           importMeshAction, saveAction, saveAs};
    connect(meshAction_, &QAction::triggered, this, &MainWindow::requestMesh);
    connect(cancelAction_, &QAction::triggered, this, [this] {
        simulationDock_->setRealtimeEnabled(false);
        meshWorker_.request_stop();
        setDiagnostic(tr("Vernetzung wird abgebrochen …"));
    });

    auto *center = new QWidget;
    auto *centerLayout = new QVBoxLayout(center);
    centerLayout_ = centerLayout;
    centerLayout->setContentsMargins(0, 0, 0, 0);
    centerLayout->setSpacing(0);
    auto *workflow = new QWidget;
    workflow->setObjectName("draftWorkflow");
    auto *workflowLayout = new QHBoxLayout(workflow);
    workflowLayout->setContentsMargins(10, 5, 7, 5);
    auto *draftStep = new QToolButton;
    draftStep->setText(tr("1 · Entwurf"));
    draftStep->setToolTip(
        tr("Draft-Parameter einblenden. Für ein neues Grundhorn den Horn-Assistenten öffnen."));
    connect(draftStep, &QToolButton::clicked, this, [this] {
        auto *panel = cadMode() ? static_cast<QDockWidget *>(cadDock_) : parameterDock_;
        panel->show();
        panel->raise();
    });
    workflowMesh_ = new QToolButton;
    workflowMesh_->setObjectName("workflowMesh");
    workflowMesh_->setText(tr("2 · Netz erzeugen"));
    workflowMesh_->setToolTip(meshAction_->toolTip());
    connect(workflowMesh_, &QToolButton::clicked, this, [this] {
        meshSettingsDock_->show();
        meshSettingsDock_->raise();
        requestMesh();
    });
    workflowSimulation_ = new QToolButton;
    workflowSimulation_->setObjectName("workflowSimulation");
    workflowSimulation_->setText(tr("3 · Berechnung …"));
    workflowSimulation_->setToolTip(tr("Berechnungseinstellungen öffnen, Frequenzen und Anregung "
                                       "prüfen, dann Berechnen starten."));
    connect(workflowSimulation_, &QToolButton::clicked, this, [this] {
        simulationDock_->show();
        simulationDock_->raise();
    });
    workflowStatus_ = new QLabel;
    workflowStatus_->setObjectName("workflowStatus");
    workflowStatus_->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Preferred);
    for (auto *button : {draftStep, workflowMesh_, workflowSimulation_}) {
        button->setFocusPolicy(Qt::StrongFocus);
        workflowLayout->addWidget(button);
    }
    workflowLayout->addWidget(workflowStatus_, 1);
    auto *workflowScroll = new QScrollArea;
    workflowScroll_ = workflowScroll;
    workflowScroll->setObjectName("workflowScroll");
    workflowScroll->setFrameShape(QFrame::NoFrame);
    workflowScroll->setWidgetResizable(true);
    workflowScroll->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    workflowScroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    workflowScroll->setMinimumWidth(0);
    workflowScroll->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Fixed);
    workflowScroll->setFixedHeight(workflow->sizeHint().height());
    workflowScroll->setWidget(workflow);
    connect(workflowScroll->horizontalScrollBar(), &QScrollBar::rangeChanged, this,
            [this, workflowScroll, workflow](int, int maximum) {
                workflowScroll->setFixedHeight(
                    workflow->sizeHint().height() +
                    (maximum > 0 ? style()->pixelMetric(QStyle::PM_ScrollBarExtent) : 0));
            });
    centerLayout->addWidget(workflowScroll);
    auto *viewBar = new QWidget;
    auto *viewLayout = new QGridLayout(viewBar);
    viewLayout->setSpacing(4);
    viewLayout->setContentsMargins(10, 4, 7, 4);
    sourceMode_ = new QComboBox;
    sourceMode_->setObjectName("sourceMode");
    sourceMode_->setAccessibleName(tr("Geometriequelle"));
    sourceMode_->setToolTip(tr("Draft bearbeitet das parametrische Horn. FreeCAD-Link verwendet "
                               "die importierte oder live verbundene CAD-Geometrie."));
    sourceMode_->addItems({tr("Draft"), tr("FreeCAD-Link")});
    viewLayout->addWidget(sourceMode_, 0, 0);
    representation_ = new QComboBox;
    representation_->setObjectName(QStringLiteral("geometryRepresentation"));
    representation_->setAccessibleName(tr("Geometriedarstellung"));
    representation_->setToolTip(tr("Vorschau folgt dem Entwurf. Netz zeigt die zuletzt erzeugte "
                                   "Vernetzung mit ihrer Revision."));
    representation_->addItems({tr("Vorschau"), tr("Netz")});
    representation_->setMinimumWidth(110);
    viewLayout->addWidget(representation_, 0, 1);
    auto *meshColor = new QComboBox;
    meshColor_ = meshColor;
    meshColor->setObjectName(QStringLiteral("meshColorMode"));
    meshColor->setAccessibleName(tr("Flächendarstellung"));
    meshColor->addItems({tr("Randbedingungen"), tr("Elementgröße"), tr("Min. Winkel"),
                         tr("Seitenverhältnis"), tr("Flächennormalen"), tr("Flächenfarben")});
    meshColor->setToolTip(
        tr("Färbung des angezeigten Netzes. Flächennormalen: Vorderseite (+n) grün, Rückseite "
           "(−n) rot, entsprechend der tatsächlichen Dreiecksorientierung. Die Auswahl bleibt "
           "hell markiert. Flächenfarben: frei zugewiesene RGB-Farben, sonst Neutralgrau; "
           "ohne Einfluss auf Geometrie oder Akustik."));
    meshColor->setWhatsThis(meshColor->toolTip());
    viewLayout->addWidget(meshColor, 0, 2);
    auto *faceColorsMenu = view->addMenu(tr("Flächen&farben"));
    assignFaceColor_ = faceColorsMenu->addAction(tr("Farbe zuweisen …"));
    assignFaceColor_->setObjectName("assignFaceColor");
    assignFaceColor_->setShortcut(QKeySequence("Ctrl+Shift+C"));
    assignFaceColor_->setToolTip(tr("Ausgewählten Flächen eine frei wählbare Farbe zuweisen. "
                                    "Flächen zuerst im 3D-Modell oder der Flächentabelle wählen; "
                                    "Strg erweitert die Auswahl. Reine Anzeigeeinstellung."));
    assignFaceColor_->setWhatsThis(assignFaceColor_->toolTip());
    resetFaceColors_ = faceColorsMenu->addAction(tr("Farbe zurücksetzen"));
    resetFaceColors_->setObjectName("resetFaceColors");
    resetFaceColors_->setToolTip(tr("Farbzuweisung der ausgewählten Flächen entfernen. "
                                    "Im Modus Flächenfarben erscheinen diese wieder neutralgrau."));
    resetFaceColors_->setWhatsThis(resetFaceColors_->toolTip());
    assignFaceColor_->setEnabled(false);
    resetFaceColors_->setEnabled(false);
    auto *colorTool = new QToolButton;
    colorTool->setObjectName("faceColorTool");
    colorTool->setDefaultAction(assignFaceColor_);
    colorTool->setText(tr("Farbe …"));
    colorTool->setMenu(faceColorsMenu);
    colorTool->setPopupMode(QToolButton::MenuButtonPopup);
    colorTool->setFocusPolicy(Qt::StrongFocus);
    viewLayout->addWidget(colorTool, 0, 3);
    connect(assignFaceColor_, &QAction::triggered, this, [this] { chooseFaceColor(); });
    connect(resetFaceColors_, &QAction::triggered, this,
            [this] { setSelectedFaceColor(std::nullopt); });
    auto *wireframe = new QCheckBox(tr("Kanten"));
    wireframe->setObjectName("showMeshEdges");
    wireframe->setChecked(QSettings().value("view/wireframe", false).toBool());
    viewLayout->addWidget(wireframe, 1, 0);
    showAperture_ = new QCheckBox(tr("Apertur"));
    showAperture_->setObjectName("showRadiationAperture");
    showAperture_->setToolTip(tr("Rechenöffnung zur unendlichen Schallwand anzeigen. "
                                 "Ausblenden verändert nur die Ansicht, nicht das Netz."));
    showAperture_->setEnabled(false);
    viewLayout->addWidget(showAperture_, 1, 1);
    showInterfaces_ = new QCheckBox(tr("Interfaces"));
    showInterfaces_->setObjectName("showBemLemInterfaces");
    showInterfaces_->setToolTip(tr(
        "BE-Interfaces als cyanfarbenes Drahtnetz; Auswahl importierter Flächen in der "
        "Flächentabelle. "
        "BEM↔LEM-Flächen bleiben getrennt zugeordnet. Die Sichtbarkeit verändert keine Rechnung."));
    showInterfaces_->setWhatsThis(showInterfaces_->toolTip());
    showInterfaces_->setEnabled(false);
    viewLayout->addWidget(showInterfaces_, 1, 2);
    auto *showGuides = new QCheckBox(tr("Modellebenen"));
    showGuides->setObjectName("showAcousticGuides");
    showGuides->setChecked(true);
    showGuides->setToolTip(tr("Symmetrie- und Schallwandebenen anzeigen. Ihre dargestellte "
                              "Ausdehnung ist symbolisch; die Ebenen sind unendlich."));
    showGuides->setWhatsThis(showGuides->toolTip());
    viewLayout->addWidget(showGuides, 2, 0, 1, 2);
    auto *showInterfaceNormals = new QCheckBox(tr("Interface-Normalen"));
    showInterfaceNormals->setObjectName("showInterfaceNormals");
    showInterfaceNormals->setChecked(true);
    showInterfaceNormals->setToolTip(tr("+n zeigt in Subdomain 1, −n in Subdomain 2. "
                                        "Die Pfeile ändern weder Netz noch Rechnung."));
    showInterfaceNormals->setWhatsThis(showInterfaceNormals->toolTip());
    viewLayout->addWidget(showInterfaceNormals, 2, 2, 1, 2);
    showFullSymmetry_ = new QCheckBox(tr("Ganzes Modell"));
    showFullSymmetry_->setObjectName("showFullSymmetry");
    showFullSymmetry_->setChecked(QSettings().value("view/fullSymmetry", true).toBool());
    showFullSymmetry_->setEnabled(false);
    showFullSymmetry_->setToolTip(tr("Das Viertelmodell an beiden Symmetrieebenen vollständig "
                                     "spiegeln. Auswahl und Flächenfarben gelten auch für die "
                                     "Spiegelbilder; das Rechennetz bleibt unverändert."));
    showFullSymmetry_->setWhatsThis(showFullSymmetry_->toolTip());
    viewLayout->addWidget(showFullSymmetry_, 3, 0, 1, 2);
    // Keep the native controls reachable on narrow displays without forcing
    // the entire center column to their combined minimum width.
    auto *viewScroll = new QScrollArea;
    viewScroll_ = viewScroll;
    viewScroll->setObjectName("viewControlsScroll");
    viewScroll->setFrameShape(QFrame::NoFrame);
    viewScroll->setWidgetResizable(true);
    viewScroll->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    viewScroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    viewScroll->setMinimumWidth(0);
    viewScroll->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Fixed);
    viewScroll->setFixedHeight(viewBar->sizeHint().height());
    viewScroll->setWidget(viewBar);
    connect(viewScroll->horizontalScrollBar(), &QScrollBar::rangeChanged, this,
            [this, viewScroll, viewBar](int, int maximum) {
                viewScroll->setFixedHeight(
                    viewBar->sizeHint().height() +
                    (maximum > 0 ? style()->pixelMetric(QStyle::PM_ScrollBarExtent) : 0));
            });
    centerLayout->addWidget(viewScroll);
    auto *compactTools = new QToolButton;
    compactTools_ = compactTools;
    compactTools->setObjectName("compactViewTools");
    compactTools->setCheckable(true);
    compactTools->setText(tr("Werkzeuge öffnen"));
    compactTools->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    compactTools->setToolTip(tr("Arbeitsablauf und 3D-Ansichtssteuerung ein- oder ausblenden."));
    centerLayout->insertWidget(0, compactTools, 0, Qt::AlignLeft);
    compactToolsDialog_ = new QDialog(this);
    compactToolsDialog_->setObjectName("compactToolsDialog");
    compactToolsDialog_->setWindowTitle(tr("Arbeitsablauf und 3D-Ansicht"));
    compactToolsDialog_->resize(700, 170);
    new QVBoxLayout(compactToolsDialog_);
    connect(compactTools, &QToolButton::toggled, this, [this, compactTools](bool expanded) {
        if (compactHeight_ && expanded)
            compactToolsDialog_->show();
        else
            compactToolsDialog_->hide();
        compactTools->setText(expanded ? tr("Werkzeuge geöffnet") : tr("Werkzeuge öffnen"));
    });
    connect(compactToolsDialog_, &QDialog::finished, compactTools, [compactTools] {
        compactTools->setChecked(false);
        QSettings().setValue("view/compactToolsExpanded", false);
    });
    viewport_ = new MeshViewport;
    viewport_->setObjectName("meshViewport");
    viewport_->setMinimumHeight(100);
    viewport_->setWireframe(wireframe->isChecked());
    centerLayout->addWidget(viewport_, 1);
    setCentralWidget(center);
    connect(wireframe, &QCheckBox::toggled, viewport_, &MeshViewport::setWireframe);
    connect(showGuides, &QCheckBox::toggled, viewport_, &MeshViewport::setAcousticGuidesVisible);
    connect(showInterfaceNormals, &QCheckBox::toggled, viewport_,
            &MeshViewport::setInterfaceNormalsVisible);
    viewport_->setFullSymmetryVisible(showFullSymmetry_->isChecked());
    connect(showFullSymmetry_, &QCheckBox::toggled, viewport_,
            &MeshViewport::setFullSymmetryVisible);
    connect(showFullSymmetry_, &QCheckBox::toggled, this,
            [](bool enabled) { QSettings().setValue("view/fullSymmetry", enabled); });
    connect(wireframe, &QCheckBox::toggled, this,
            [](bool enabled) { QSettings().setValue("view/wireframe", enabled); });
    auto *backgroundGrid = view->addAction(tr("Hintergrundraster"));
    backgroundGrid->setObjectName("showBackgroundGrid");
    backgroundGrid->setCheckable(true);
    backgroundGrid->setToolTip(tr("Orientierungsraster im 3D-Hintergrund anzeigen. "
                                  "Dreieckskanten werden separat über Kanten geschaltet."));
    backgroundGrid->setWhatsThis(backgroundGrid->toolTip());
    backgroundGrid->setChecked(QSettings().value("view/backgroundGrid", false).toBool());
    viewport_->setBackgroundGrid(backgroundGrid->isChecked());
    connect(backgroundGrid, &QAction::toggled, this, [this](bool enabled) {
        viewport_->setBackgroundGrid(enabled);
        QSettings().setValue("view/backgroundGrid", enabled);
    });
    connect(meshColor, &QComboBox::currentIndexChanged, this,
            [this](int index) { viewport_->setColorMode(static_cast<MeshColorMode>(index)); });
    connect(representation_, &QComboBox::currentIndexChanged, this, [this] { updateDisplay(); });
    connect(showAperture_, &QCheckBox::toggled, this, [this] { updateDisplay(); });
    connect(showInterfaces_, &QCheckBox::toggled, this, [this] { updateDisplay(); });
    auto *fitAction = view->addAction(tr("Einpassen"), viewport_, &MeshViewport::fitView);
    fitAction->setShortcut(QKeySequence("F"));
    fitAction->setToolTip(tr("Gesamtes Modell in die 3D-Ansicht einpassen (F)."));
    auto *fitTool = new QToolButton;
    fitTool->setObjectName("fitViewTool");
    fitTool->setDefaultAction(fitAction);
    fitTool->setFocusPolicy(Qt::StrongFocus);
    viewLayout->addWidget(fitTool, 1, 3);
    auto *frontAction = view->addAction(tr("Vorderansicht"));
    frontAction->setObjectName("frontView");
    frontAction->setShortcut(QKeySequence("Ctrl+1"));
    frontAction->setToolTip(
        tr("Blick von der positiven akustischen Vorwärtsachse auf das Modell. "
           "Verwendet den eingestellten Beobachter-/Akustikrahmen und passt die Ansicht ein. "
           "Geometrie und Simulation bleiben unverändert."));
    frontAction->setWhatsThis(frontAction->toolTip());
    connect(frontAction, &QAction::triggered, this,
            [this] { viewport_->viewFront(simulationDock_->observationFrame().forward); });
    auto *frontTool = new QToolButton;
    frontTool->setObjectName("frontViewTool");
    frontTool->setDefaultAction(frontAction);
    frontTool->setText(tr("Vorne"));
    frontTool->setFocusPolicy(Qt::StrongFocus);
    viewLayout->addWidget(frontTool, 1, 4);
    view->addSeparator();

    fieldVisibility_ = view->addAction(tr("Schallfeldebene"));
    fieldVisibility_->setObjectName("showPressureField");
    fieldVisibility_->setCheckable(true);
    fieldVisibility_->setChecked(true);
    fieldVisibility_->setEnabled(false);
    fieldVisibility_->setToolTip(tr("Berechnete H/V-Ebene im oberen 3D-Modell zeigen. "
                                    "CAD und Feld verwenden dieselbe Kamera und Tiefenprüfung."));
    connect(fieldVisibility_, &QAction::toggled, this, [this](bool visible) {
        viewport_->setPressureFieldVisible(visible);
        if (visible && viewport_->pressureField() && !fieldCameraInitialized_) {
            viewport_->focusPressureField();
            fieldCameraInitialized_ = true;
        }
    });
    view->addAction(tr("3D-Bild exportieren …"), this, [this] {
        const auto path = QFileDialog::getSaveFileName(this, tr("3D-Ansicht exportieren"),
                                                       "HornSculpt-3D.png", tr("PNG (*.png)"));
        if (!path.isEmpty() && !viewport_->grabFramebuffer().save(path, "PNG"))
            QMessageBox::warning(this, tr("3D-Bild exportieren"),
                                 tr("Die Bilddatei konnte nicht gespeichert werden."));
    });

    auto *tree = new QTreeWidget;
    tree->setObjectName("projectTree");
    tree->setHeaderHidden(true);
    tree->setMinimumWidth(170);
    projectItem_ = new QTreeWidgetItem(tree, {tr("Entwurf")});
    projectItem_->setIcon(0, style()->standardIcon(QStyle::SP_DirIcon));
    sourceItem_ = new QTreeWidgetItem(projectItem_, {tr("Geometrie · OS-SE")});
    facesItem_ = new QTreeWidgetItem(projectItem_, {tr("Flächen")});
    new QTreeWidgetItem(facesItem_, {tr("Hornwand · Wall")});
    new QTreeWidgetItem(facesItem_, {tr("Hals · Driven")});
    meshItem_ = new QTreeWidgetItem(projectItem_, {tr("Netz · nicht erzeugt")});
    domainsItem_ = new QTreeWidgetItem(projectItem_, {tr("Akustische Gebiete")});
    domainsItem_->setHidden(true);
    symmetryItem_ = new QTreeWidgetItem(projectItem_, {tr("Symmetrie · x = 0, y = 0")});
    symmetryItem_->setHidden(true);
    tree->expandAll();
    auto *treeDock = dock(this, tr("Projekt"), "projectDock", tree, Qt::LeftDockWidgetArea);
    view->addAction(treeDock->toggleViewAction());

    parameters_ = new ParameterPanel;
    parameters_->setMinimumWidth(310);
    parameterDock_ = dock(this, tr("Parameter"), "parameterDock",
                          scrollableDockContent(parameters_, "parameterDockScroll", 160),
                          Qt::RightDockWidgetArea);
    view->addAction(parameterDock_->toggleViewAction());
    cadDock_ = new FreeCadDock(this);
    addDockWidget(Qt::RightDockWidgetArea, cadDock_);
    view->addAction(cadDock_->toggleViewAction());
    meshPanel_ = new MeshPanel;
    meshSettingsDock_ = dock(this, tr("Netzeinstellungen"), "meshSettingsDock",
                             scrollableDockContent(meshPanel_, "meshSettingsDockScroll", 160),
                             Qt::LeftDockWidgetArea);
    splitDockWidget(treeDock, meshSettingsDock_, Qt::Vertical);
    view->addAction(meshSettingsDock_->toggleViewAction());

    auto *plots = new QWidget;
    auto *plotLayout = new QGridLayout(plots);
    plotLayout->setContentsMargins(1, 1, 1, 1);
    plotLayout->setSpacing(1);
    volumePlot_ = new PlotWidget(tr("Volumen V(z)"), tr("z [mm]"), tr("cm³"));
    areaPlot_ = new PlotWidget(tr("Querschnitt dV/dz"), tr("z [mm]"), tr("mm²"));
    derivativePlot_ = new PlotWidget(tr("Flächenänderung d²V/dz²"), tr("z [mm]"), tr("mm"));
    anglePlot_ = new PlotWidget(tr("Wandwinkel"), tr("z [mm]"), tr("°"));
    volumePlot_->setObjectName("geometryVolume");
    areaPlot_->setObjectName("geometryArea");
    derivativePlot_->setObjectName("geometryAreaDerivative");
    anglePlot_->setAxisScales({}, {0, 90, 15});
    anglePlot_->setObjectName("geometryWallAngles");
    plotLayout->addWidget(volumePlot_, 0, 0);
    plotLayout->addWidget(areaPlot_, 0, 1);
    plotLayout->addWidget(derivativePlot_, 1, 0);
    plotLayout->addWidget(anglePlot_, 1, 1);
    geometryTabs_ = new QTabWidget;
    geometryTabs_->setObjectName(QStringLiteral("geometryTabs"));
    auto *contours = new QWidget;
    auto *contourLayout = new QVBoxLayout(contours);
    contourLayout->setContentsMargins(3, 3, 3, 3);
    contourLayout->setSpacing(2);
    auto *contourTools = new QHBoxLayout;
    auto *fitContours = new QPushButton(tr("Einpassen"));
    auto *focusThroat = new QPushButton(tr("Hals"));
    fitContours->setObjectName(QStringLiteral("fitContours"));
    focusThroat->setObjectName(QStringLiteral("focusContourThroat"));
    focusThroat->setToolTip(tr("Auf den Schnitt der angetriebenen Halsfläche vergrößern."));
    sectionMesh_ = new QCheckBox(tr("Netz"));
    sectionMesh_->setObjectName(QStringLiteral("sectionMeshOverlay"));
    sectionMesh_->setToolTip(tr("BEM-Netz gestrichelt in Grün mit der Vorschau vergleichen."));
    auto *drivenLegend = new QLabel(tr("Driven"));
    auto legendPalette = drivenLegend->palette();
    legendPalette.setColor(QPalette::WindowText, QColor(177, 49, 43));
    drivenLegend->setPalette(legendPalette);
    contourTools->addWidget(fitContours);
    contourTools->addWidget(focusThroat);
    contourTools->addStretch();
    contourTools->addWidget(drivenLegend);
    contourTools->addWidget(sectionMesh_);
    contourLayout->addLayout(contourTools);
    auto *sections = new QSplitter(Qt::Horizontal);
    horizontalSection_ = new MeshSectionView(tr("H · Horizontal"), SectionPlane::xz);
    verticalSection_ = new MeshSectionView(tr("V · Vertikal"), SectionPlane::yz);
    horizontalSection_->setObjectName(QStringLiteral("horizontalContour"));
    verticalSection_->setObjectName(QStringLiteral("verticalContour"));
    sections->addWidget(horizontalSection_);
    sections->addWidget(verticalSection_);
    contourLayout->addWidget(sections, 1);
    geometryTabs_->addTab(contours, tr("Kontur H/V"));
    geometryTabs_->addTab(scrollableDockContent(plots, "plotsScroll"), tr("Analyse"));
    connect(fitContours, &QPushButton::clicked, this, [this] {
        horizontalSection_->fitView();
        verticalSection_->fitView();
    });
    connect(focusThroat, &QPushButton::clicked, this, [this] {
        horizontalSection_->focusDriven();
        verticalSection_->focusDriven();
    });
    connect(sectionMesh_, &QCheckBox::toggled, this, [this] { updateDisplay(); });
    analysisDock_ =
        dock(this, tr("Geometrieanalyse"), "geometryAnalysisDock",
             scrollableDockContent(geometryTabs_, "geometryDockScroll"), Qt::BottomDockWidgetArea);
    view->addAction(analysisDock_->toggleViewAction());
    simulationDock_ = new SimulationDock(this);
    addDockWidget(Qt::BottomDockWidgetArea, simulationDock_);
    view->addAction(simulationDock_->toggleViewAction());
    auto *enlargeResults = view->addAction(tr("Ergebnisse vergrößern"));
    enlargeResults->setObjectName("enlargeResults");
    enlargeResults->setShortcut(QKeySequence("Ctrl+Shift+E"));
    enlargeResults->setCheckable(true);
    enlargeResults->setToolTip(
        tr("Ergebnisdock vergrößern; erneut wählen zeigt die 3D-Ansicht wieder."));
    connect(enlargeResults, &QAction::toggled, this, [this](bool enlarged) {
        viewport_->setVisible(!enlarged);
        simulationDock_->show();
        simulationDock_->raise();
        simulationDock_->revealResultView();
        const auto generation = ++dockResizeGeneration_;
        QTimer::singleShot(0, this, [this, enlarged, generation] {
            if (generation != dockResizeGeneration_)
                return;
            resizeDocks({simulationDock_}, {height() * (enlarged ? 7 : 3) / 10}, Qt::Vertical);
        });
    });
    tabifyDockWidget(analysisDock_, simulationDock_);
    lemEditor_ = new LemEditor(this);
    lemEditor_->setObjectName("lemEditor");
    auto *lemDock =
        dock(this, tr("Treiber / LEM"), "lemDock",
             scrollableDockContent(lemEditor_, "lemDockScroll"), Qt::BottomDockWidgetArea);
    view->addAction(lemDock->toggleViewAction());
    tabifyDockWidget(simulationDock_, lemDock);
    lemEditor_->onChanged = [this] { couplingChanged(); };
    lemEditor_->onEnabledChanged = [this](bool) { couplingChanged(); };
    lemEditor_->onFacesSelected = [this](std::vector<std::uint32_t> faces) {
        viewport_->setSelectedFaces(faces);
        meshPanel_->setSelectedFaces(faces);
        if (cadMode())
            cadDock_->setSelectedFaces(faces);
        updateFaceColorActions();
    };
    auto *resultTabs = simulationDock_->findChild<QTabWidget *>("acousticResults");
    connect(resultTabs, &QTabWidget::currentChanged, this, [resultTabs, enlargeResults](int index) {
        auto *page = resultTabs->widget(index);
        if (page && page->findChild<QWidget *>("acousticField") && enlargeResults->isChecked())
            enlargeResults->setChecked(false);
    });
    simulationDock_->onLemResultChanged = [this](std::optional<lem::Solution> solution,
                                                 bool stale) {
        lemEditor_->setSolution(std::move(solution), stale);
    };
    simulationDock_->onFieldChanged = [this, enlargeResults,
                                       resultTabs](std::shared_ptr<const PressureFieldPlane> plane,
                                                   PressureFieldAppearance appearance) {
        auto *page = resultTabs->currentWidget();
        if (plane && page && page->findChild<QWidget *>("acousticField") &&
            enlargeResults->isChecked())
            enlargeResults->setChecked(false);
        fieldVisibility_->setEnabled(bool(plane));
        viewport_->setPressureField(std::move(plane), appearance);
        viewport_->setPressureFieldVisible(fieldVisibility_->isChecked());
        if (viewport_->pressureFieldVisible() && !fieldCameraInitialized_) {
            viewport_->focusPressureField();
            fieldCameraInitialized_ = true;
        }
    };
    simulationDock_->onFieldProbeModeChanged = [this](bool enabled) {
        viewport_->setFieldProbeEnabled(enabled);
    };
    simulationDock_->onRequestConvergence = [this] { openParameterStudy(true); };
    simulationDock_->onFocusField = [this, enlargeResults] {
        if (enlargeResults->isChecked())
            enlargeResults->setChecked(false);
        fieldVisibility_->setChecked(true);
        viewport_->focusPressureField();
        fieldCameraInitialized_ = true;
    };
    analysisDock_->raise();
    resizeDocks({treeDock, parameterDock_}, {200, 340}, Qt::Horizontal);
    resizeDocks({analysisDock_}, {260}, Qt::Vertical);

    revisionLabel_ = new QLabel;
    qualityLabel_ = new QLabel(tr("Netz: nicht erzeugt"));
    diagnostic_ = new QLabel;
    statusBar()->addWidget(revisionLabel_);
    statusBar()->addWidget(diagnostic_, 1);
    statusBar()->addPermanentWidget(qualityLabel_);
    debounce_ = new QTimer(this);
    debounce_->setSingleShot(true);
    debounce_->setInterval(40);
    connect(debounce_, &QTimer::timeout, this, &MainWindow::startGeometry);
    parameters_->onChanged = [this] { parametersChanged(); };
    meshPanel_->onChanged = [this] { meshingChanged(); };
    simulationDock_->onChanged = [this] {
        markDirty();
        configureLemFaces();
        updateDisplay();
    };
    simulationDock_->onAcousticRegionsChanged = [this] {
        rememberEdit(true);
        if (gmsh_)
            meshingChanged();
    };
    simulationDock_->onSymmetryChanged = [this] {
        rememberEdit(true);
        if (gmsh_)
            meshingChanged();
        updateDisplay();
    };
    simulationDock_->onObservationFrameChanged = [this] {
        viewport_->setObservationFrame(simulationDock_->observationFrame());
    };
    realtimeMeshTimer_ = new QTimer(this);
    realtimeMeshTimer_->setInterval(200);
    connect(realtimeMeshTimer_, &QTimer::timeout, this, [this] {
        if (closing_ || loading_ || !simulationDock_->realtimeEnabled() || geometryBusy_ ||
            meshBusy_ || !geometryReady() || meshReady() ||
            simulationDock_->usesAxisymmetricGeometry() || realtimeMeshAttempt_ == revision_)
            return;
        realtimeMeshAttempt_ = revision_;
        requestMesh();
    });
    simulationDock_->onRealtimeChanged = [this](bool enabled) {
        realtimeMeshAttempt_ = 0;
        if (enabled)
            realtimeMeshTimer_->start();
        else
            realtimeMeshTimer_->stop();
    };
    cadDock_->onSnapshotChanged = [this](std::shared_ptr<const CadSnapshot> snapshot) {
        acceptCadSnapshot(std::move(snapshot));
    };
    cadDock_->onAssignmentEditStarted = [this] { flushEdit(); };
    cadDock_->onEditBemRegions = [this] {
        simulationDock_->editAcousticRegions(cadDock_->selectedFaces());
    };
    cadDock_->onAvailabilityChanged = [this](bool available) { cadAvailabilityChanged(available); };
    cadDock_->onError = [this](QString error) {
        if (cadMode())
            setDiagnostic(error, true);
    };
    cadDock_->onFacesSelected = [this](std::vector<std::uint32_t> faces) {
        if (cadMode()) {
            meshPanel_->setSelectedFaces(faces);
            viewport_->setSelectedFaces(std::move(faces));
            updateFaceColorActions();
        }
    };
    meshPanel_->onFacesSelected = [this](std::vector<std::uint32_t> faces) {
        viewport_->setSelectedFaces(placedSelection(faces));
        if (cadMode())
            cadDock_->setSelectedFaces(faces);
        updateFaceColorActions();
    };
    meshPanel_->onFaceTargetChanged = [this](std::vector<std::uint32_t> ids,
                                             std::optional<double> target) {
        setFaceTargets(ids, target);
    };
    meshPanel_->onFaceOrientationChanged = [this](std::vector<std::uint32_t> ids) {
        reverseFaces(ids);
    };
    cadDock_->onFaceTargetChanged = [this](std::vector<std::uint32_t> ids,
                                           std::optional<double> target) {
        if (cadMode())
            setFaceTargets(ids, target);
    };
    viewport_->onFacePicked = [this](std::optional<std::uint32_t> face,
                                     Qt::KeyboardModifiers modifiers) {
        if (cadMode()) {
            const auto current = cadDock_->snapshot();
            const bool currentSource = current && cadSnapshot_ &&
                                       current->documentId == cadSnapshot_->documentId &&
                                       current->sessionId == cadSnapshot_->sessionId &&
                                       current->revision == cadSnapshot_->revision;
            const bool currentView =
                geometryRevision_ == revision_ &&
                (representation_->currentIndex() != 1 || !mesh_ || meshRevision_ == revision_);
            // Historical meshes remain visible, but their IDs cannot authorize
            // edits of the current CAD source, even when the numeric IDs match.
            if (!currentSource || !currentView) {
                cadDock_->clearSelection();
                meshPanel_->setSelectedFaces({});
                viewport_->setSelectedFaces({});
                updateFaceColorActions();
                setDiagnostic(tr("Veraltete CAD-Ansicht · Flächenauswahl verworfen."));
                return;
            }
            if (face)
                cadDock_->selectFace(*face, modifiers);
            else if (modifiers == Qt::NoModifier)
                cadDock_->clearSelection();
        } else if (face) {
            if (*face < assemblyCadFaceOffset)
                meshPanel_->selectFace((*face - 1) % 8 + 1, modifiers);
        } else if (modifiers == Qt::NoModifier) {
            meshPanel_->setSelectedFaces({});
            viewport_->setSelectedFaces({});
        }
        updateFaceColorActions();
    };
    connect(sourceMode_, &QComboBox::currentIndexChanged, this, [this] { sourceModeChanged(); });
    updateSourcePresentation();
    const auto defaultLayout = saveState();
    view->addSeparator();
    auto *resetLayout = view->addAction(tr("Docklayout zurücksetzen"));
    resetLayout->setObjectName("resetDockLayout");
    resetLayout->setToolTip(tr("Projekt links, Parameter rechts und Ergebnisse unten anordnen. "
                               "Modell und Ergebnisse bleiben erhalten."));
    const auto applyDefaultLayout = [this, defaultLayout, treeDock, compactTools] {
        restoreState(defaultLayout);
        updateSourcePresentation();
        const bool compactHeight = height() < 650;
        setCompactToolsMode(compactHeight);
        compactTools->setChecked(false);
        const bool compact = width() < 1100;
        horizontalSection_->setMinimumWidth(compact ? 220 : 280);
        verticalSection_->setMinimumWidth(compact ? 220 : 280);
        if (compact) {
            tabifyDockWidget(treeDock, meshSettingsDock_);
            tabifyDockWidget(treeDock, parameterDock_);
            tabifyDockWidget(treeDock, cadDock_);
            (cadMode() ? static_cast<QDockWidget *>(cadDock_) : parameterDock_)->raise();
            resizeDocks({cadMode() ? static_cast<QDockWidget *>(cadDock_) : parameterDock_}, {320},
                        Qt::Horizontal);
        } else {
            resizeDocks({treeDock, parameterDock_}, {250, 340}, Qt::Horizontal);
        }
        resizeDocks({analysisDock_}, {std::clamp(height() * 35 / 100, 200, 380)}, Qt::Vertical);
        analysisDock_->show();
        analysisDock_->raise();
    };
    connect(resetLayout, &QAction::triggered, this, [this, enlargeResults, applyDefaultLayout] {
        ++dockResizeGeneration_;
        if (enlargeResults->isChecked()) {
            const QSignalBlocker blocker(enlargeResults);
            enlargeResults->setChecked(false);
            viewport_->show();
        }
        applyDefaultLayout();
    });
    QSettings settings;
    const bool restoredGeometry = restoreGeometry(settings.value("window/geometry").toByteArray());
    if (!restoredGeometry && screen()) {
        const auto available = screen()->availableGeometry();
        resize(std::min(width(), available.width()), std::min(height(), available.height()));
    }
    const bool restoredLayout = restoreState(settings.value("window/state").toByteArray());
    const bool compactHeight = height() < 650;
    setCompactToolsMode(compactHeight);
    compactTools->setChecked(compactHeight &&
                             settings.value("view/compactToolsExpanded", false).toBool());
    connect(compactTools, &QToolButton::clicked, this,
            [](bool expanded) { QSettings().setValue("view/compactToolsExpanded", expanded); });
    updateSourcePresentation();
    analysisDock_->show();
    analysisDock_->raise();
    updateTitle();
    resetEditHistory();
    qApp->installEventFilter(this);
    startGeometry();
    auto *workflowTimer = new QTimer(this);
    workflowTimer->setInterval(300);
    connect(workflowTimer, &QTimer::timeout, this, [this] { updateWorkflow(); });
    workflowTimer->start();
    updateWorkflow();
    if (!restoredLayout)
        QTimer::singleShot(0, this, applyDefaultLayout);
    const auto arguments = QCoreApplication::arguments();
    const bool automated =
        arguments.contains("--smoke-test") || arguments.contains("--mesh-test") ||
        arguments.contains("--solver-test") || arguments.contains("--project-test") ||
        arguments.contains("--screenshot") || arguments.contains("--solve-project");
    if (!recoveryRoot_.isEmpty() ||
        (QCoreApplication::applicationName() == QStringLiteral("HornSculpt") && !automated))
        QTimer::singleShot(0, this, &MainWindow::startRecovery);
}

MainWindow::~MainWindow() {
    qApp->removeEventFilter(this);
    closing_ = true;
    editTimer_->stop();
    disconnect(editTimer_, nullptr, this, nullptr);
    // QWidget deletes its children after our edit state and action members are gone.
    // QUndoStack destruction can emit availability signals while clearing commands.
    disconnect(undoStack_, nullptr, this, nullptr);
    undoStack_->clear();
    for (auto *child : findChildren<QObject *>())
        disconnect(child, nullptr, this, nullptr);
    parameters_->onChanged = {};
    meshPanel_->onChanged = {};
    meshPanel_->onFacesSelected = {};
    meshPanel_->onFaceTargetChanged = {};
    meshPanel_->onFaceOrientationChanged = {};
    cadDock_->onSnapshotChanged = {};
    cadDock_->onAssignmentEditStarted = {};
    cadDock_->onEditBemRegions = {};
    cadDock_->onAvailabilityChanged = {};
    cadDock_->onError = {};
    cadDock_->onFacesSelected = {};
    cadDock_->onFaceTargetChanged = {};
    simulationDock_->onChanged = {};
    simulationDock_->onAcousticRegionsChanged = {};
    simulationDock_->onSymmetryChanged = {};
    simulationDock_->onRealtimeChanged = {};
    simulationDock_->onObservationFrameChanged = {};
    simulationDock_->onLemResultChanged = {};
    simulationDock_->onFieldChanged = {};
    simulationDock_->onFieldProbeModeChanged = {};
    simulationDock_->onRequestConvergence = {};
    simulationDock_->onFocusField = {};
    lemEditor_->onChanged = {};
    lemEditor_->onEnabledChanged = {};
    lemEditor_->onFacesSelected = {};
    viewport_->onFacePicked = {};
    geometryWorker_.request_stop();
    meshWorker_.request_stop();
    sourceImportWorker_.request_stop();
    geometryWorker_.shutdown();
    meshWorker_.shutdown();
    if (sourceImportWorker_.joinable())
        sourceImportWorker_.join();
}

void MainWindow::parametersChanged() {
    if (loading_ || cadMode())
        return;
    rememberEdit();
    ++revision_;
    dirty_ = true;
    geometryPending_ = true;
    parameters_->setSlotTransitionEnd(std::nullopt);
    parameters_->setGeometryStatus(tr("Vorschau wird aktualisiert …"));
    // Edits belong to the live geometry. A previously selected mesh view must
    // not hide them behind the last, now stale BEM mesh.
    representation_->setCurrentIndex(0);
    geometryWorker_.request_stop();
    meshWorker_.request_stop();
    meshAction_->setEnabled(false);
    simulationDock_->invalidate(revision_);
    revisionLabel_->setText(tr("Geometrie Rev. %1 · Änderung").arg(revision_));
    analysisDock_->setWindowTitle(
        tr("Geometrieanalyse · Rev. %1 · veraltet · wird aktualisiert").arg(analysisRevision_));
    if (mesh_)
        qualityLabel_->setText(tr("Netz Rev. %1 · veraltet").arg(meshRevision_));
    updateTitle();
    updateDisplay();
    debounce_->start();
}

void MainWindow::startGeometry() {
    if (cadMode())
        return;
    if (geometryBusy_) {
        geometryPending_ = geometryJobRevision_ != revision_;
        return;
    }
    geometryPending_ = false;
    if (!meshSettingsValid_) {
        parameters_->setGeometryStatus(tr("Vorschau nicht aktuell\nNetzeinstellungen korrigieren."),
                                       true);
        return;
    }
    DraftParameters snapshot;
    try {
        snapshot = parameters_->parameters();
    } catch (const std::exception &e) {
        parameters_->setGeometryStatus(tr("Vorschau nicht aktuell\n") + QString::fromUtf8(e.what()),
                                       true);
        setDiagnostic(QString::fromUtf8(e.what()), true);
        return;
    }
    geometryBusy_ = true;
    const auto revision = revision_;
    geometryJobRevision_ = revision;
    const auto epoch = projectEpoch_;
    const double geometricTolerance = meshing_.adaptive.geometricTolerance;
    try {
        geometryWorker_.submit(
            [this, snapshot, revision, epoch, geometricTolerance](std::stop_token stop) {
                std::shared_ptr<const DraftGeometry> geometry;
                std::shared_ptr<const SurfaceMesh> preview;
                QString error;
                try {
                    if (!stop.stop_requested()) {
                        geometry = std::make_shared<DraftGeometry>(snapshot);
                        preview = std::make_shared<SurfaceMesh>(
                            geometry->previewMesh(geometricTolerance, stop));
                    }
                } catch (const std::exception &e) {
                    error = QString::fromUtf8(e.what());
                }
                if (!stop.stop_requested()) {
                    QMetaObject::invokeMethod(
                        this,
                        [this, revision, epoch, geometry, preview, error]() mutable {
                            if (epoch == projectEpoch_)
                                acceptGeometry(revision, std::move(geometry), std::move(preview),
                                               error);
                        },
                        Qt::QueuedConnection);
                }
                GeometryAnalysis analysis;
                if (error.isEmpty() && !stop.stop_requested()) {
                    try {
                        analysis = analyzeGeometry(*geometry, 101, stop);
                    } catch (const std::exception &e) {
                        error = QString::fromUtf8(e.what());
                    }
                }
                QMetaObject::invokeMethod(
                    this,
                    [this, revision, epoch, analysis = std::move(analysis), error,
                     cancelled = stop.stop_requested()]() mutable {
                        if (epoch == projectEpoch_)
                            finishGeometry(revision, std::move(analysis),
                                           cancelled ? QString{} : error);
                    },
                    Qt::QueuedConnection);
            });
    } catch (const std::exception &e) {
        geometryBusy_ = false;
        geometryPending_ = false;
        parameters_->setGeometryStatus(
            tr("Vorschau konnte nicht gestartet werden: ") + QString::fromUtf8(e.what()), true);
        setDiagnostic(QString::fromUtf8(e.what()), true);
    }
}

void MainWindow::acceptGeometry(std::uint64_t revision,
                                std::shared_ptr<const DraftGeometry> geometry,
                                std::shared_ptr<const SurfaceMesh> preview, QString error) {
    if (closing_)
        return;
    if (revision == revision_ && !cadMode()) {
        if (error.isEmpty()) {
            parameters_->setGeometryStatus({});
            parameters_->setSlotTransitionEnd(geometry->slotTransitionEnd());
            geometry_ = std::move(geometry);
            rawPreview_ = std::move(preview);
            geometryRevision_ = revision;
            std::set<std::uint32_t> available;
            for (const auto &triangle : rawPreview_->triangles)
                available.insert(triangle.faceId);
            const auto before = draftFaceTargetEdgeLengths_.size();
            std::erase_if(draftFaceTargetEdgeLengths_,
                          [&](const auto &entry) { return !available.contains(entry.first); });
            std::erase_if(draftReversedFaceIds_,
                          [&](const auto id) { return !available.contains(id); });
            std::erase_if(draftFaceColors_,
                          [&](const auto &entry) { return !available.contains(entry.first); });
            updateOrientedPreview();
            if (before != draftFaceTargetEdgeLengths_.size())
                meshPanel_->setFaceNotice(tr("Draft · %1 entfallene lokale Zielkanten entfernt")
                                              .arg(before - draftFaceTargetEdgeLengths_.size()));
            updateMeshFaces();
            configureSimulation();
            updateDisplay();
            if (initialModelFitPending_) {
                viewport_->fitView();
                initialModelFitPending_ = false;
            }
            setDiagnostic({});
            analysisDock_->setWindowTitle(
                tr("Geometrieanalyse · Rev. %1 · veraltet · wird aktualisiert für Rev. %2")
                    .arg(analysisRevision_)
                    .arg(revision));
            revisionLabel_->setText(tr("Geometrie Rev. %1").arg(revision_));
            meshAction_->setEnabled(!meshBusy_);
        } else {
            parameters_->setGeometryStatus(tr("Vorschau nicht aktuell\n") + error +
                                               tr("\nLetzte gültige Vorschau bleibt sichtbar."),
                                           true);
            setDiagnostic(error + tr(" Letzte gültige Vorschau bleibt sichtbar."), true);
            analysisDock_->setWindowTitle(
                tr("Geometrieanalyse · Rev. %1 · veraltet · Geometrie ungültig")
                    .arg(analysisRevision_));
            revisionLabel_->setText(tr("Geometrie Rev. %1 · ungültig").arg(revision_));
        }
    }
    meshAction_->setEnabled(geometryReady() && !meshBusy_);
}

void MainWindow::finishGeometry(std::uint64_t revision, GeometryAnalysis analysis, QString error) {
    geometryBusy_ = false;
    if (closing_)
        return;
    if (revision == revision_ && !cadMode() && geometryRevision_ == revision) {
        if (error.isEmpty() && !analysis.samples.empty()) {
            updatePlots(analysis);
            analysisRevision_ = revision;
            analysisDock_->setWindowTitle(tr("Geometrieanalyse · Rev. %1").arg(revision));
        } else if (!error.isEmpty()) {
            analysisDock_->setWindowTitle(
                tr("Geometrieanalyse · Rev. %1 · veraltet · Analysefehler").arg(analysisRevision_));
            setDiagnostic(error, true);
        }
    }
    if (geometryPending_ && !cadMode())
        startGeometry();
    updateDisplay();
    meshAction_->setEnabled(geometryReady() && !meshBusy_);
}

void MainWindow::requestMesh() {
    if (!geometryReady() || meshBusy_)
        return;
    try {
        const auto snapshot = rawPreview_;
        const auto assemblyProject =
            !cadMode() ? std::make_shared<Project>(projectState()) : nullptr;
        const auto reversed = activeReversedFaces();
        const auto revision = revision_;
        const auto epoch = projectEpoch_;
        const auto jobId = ++meshJobId_;
        MeshingOptions options;
        options.targetEdgeLength = meshing_.targetEdgeLength;
        options.adaptive = meshing_.adaptive;
        options.quality = meshing_.quality;
        options.faceTargetEdgeLengths = activeFaceTargets();
        const auto cad = cadMode() ? cadSnapshot_ : nullptr;
        const bool importedGmsh =
            gmsh_.has_value() && cad && cad->sessionId == "gmsh-import-offline";
        const auto simulationSettings = simulationDock_->archive().settings;
        const auto regions = simulationSettings.acousticModel;
        if (cad)
            options.geometry = cad->geometry;
        meshBusy_ = true;
        meshError_.clear();
        meshAction_->setEnabled(false);
        cancelAction_->setEnabled(true);
        setDiagnostic(importedGmsh ? tr("Importnetz prüfen …") : tr("Isotrope Vernetzung …"));
        meshWorker_.submit([this, snapshot, revision, epoch, jobId, options, reversed, cad,
                            assemblyProject, importedGmsh, regions,
                            simulationSettings](std::stop_token stop) {
            std::shared_ptr<const MeshingResult> result;
            QString error;
            try {
                if (!stop.stop_requested()) {
                    auto jobOptions = options;
                    std::optional<MeshingSource> source;
                    if (cad) {
                        source =
                            makeCadMeshingSource(*cad, options.adaptive.geometricTolerance, stop);
                        jobOptions.parameterization = source->parameters;
                    }
                    const auto importedMesh =
                        importedGmsh
                            ? (reversed.empty() ? source->mesh
                                                : applyFaceOrientation(source->mesh, reversed))
                            : SurfaceMesh{};
                    const auto symmetry =
                        importedGmsh && simulationSettings.symmetryXYEnabled
                            ? std::optional{inferXYSymmetry(
                                  importedMesh, simulationSettings.symmetryPlaneTolerance)}
                            : std::nullopt;
                    auto generated =
                        importedGmsh ? validateGmshMesh(importedMesh, regions, jobOptions, stop,
                                                        symmetry ? &*symmetry : nullptr)
                                     : (assemblyProject
                                            ? meshDraftAssembly(*assemblyProject, *snapshot, stop)
                                            : remeshIsotropic(source ? source->mesh : *snapshot,
                                                              jobOptions, stop));
                    if (!importedGmsh && !assemblyProject && !reversed.empty()) {
                        try {
                            generated.mesh = applyFaceOrientation(generated.mesh, reversed);
                            generated.quality = inspectMesh(generated.mesh);
                            if (!generated.quality.watertight())
                                throw std::invalid_argument(
                                    "Flächenauswahl ergibt keine konsistent "
                                    "orientierte geschlossene Hülle.");
                            galerkin::validateExteriorComponents(generated.mesh);
                        } catch (const std::exception &orientationError) {
                            throw std::invalid_argument(std::string("Flächennormalen: ") +
                                                        orientationError.what());
                        }
                    }
                    result = std::make_shared<MeshingResult>(std::move(generated));
                }
            } catch (const std::exception &e) {
                error = QString::fromUtf8(e.what());
            }
            QMetaObject::invokeMethod(
                this,
                [this, revision, epoch, jobId, result = std::move(result), error = std::move(error),
                 stop] {
                    if (epoch != projectEpoch_ || jobId != meshJobId_)
                        return;
                    meshBusy_ = false;
                    cancelAction_->setEnabled(false);
                    meshAction_->setEnabled(geometryReady());
                    if (closing_)
                        return;
                    if (revision != revision_) {
                        // A CAD assignment edit can retain a current mesh while this
                        // obsolete job is retiring. Release that mesh, not its old status.
                        if (meshReady() && geometryReady())
                            configureSimulation();
                        return;
                    }
                    if (stop.stop_requested()) {
                        meshError_ = tr("Vernetzung abgebrochen");
                        setDiagnostic(meshError_);
                        return;
                    }
                    if (!error.isEmpty()) {
                        meshError_ = error;
                        setDiagnostic(error, true);
                        simulationDock_->setRealtimeEnabled(false);
                        if (error.startsWith(QStringLiteral("Flächennormalen:")))
                            meshPanel_->setOrientationStatus(
                                tr("Normalen inkonsistent · BEM gesperrt"), true);
                        return;
                    }
                    mesh_ = result;
                    meshRevision_ = revision;
                    updateMeshStatus();
                    updateMeshFaces();
                    markDirty();
                    configureSimulation();
                    representation_->setCurrentIndex(1);
                    updateDisplay();
                    setDiagnostic({});
                },
                Qt::QueuedConnection);
        });
    } catch (const std::exception &e) {
        meshBusy_ = false;
        cancelAction_->setEnabled(false);
        meshAction_->setEnabled(geometryReady());
        meshError_ = tr("Vernetzung konnte nicht gestartet werden: ") + QString::fromUtf8(e.what());
        setDiagnostic(meshError_, true);
    }
}

void MainWindow::updateDisplay() {
    stepExportAction_->setEnabled(!cadMode() && geometry_ && !geometryBusy_ && !geometryPending_ &&
                                  geometryRevision_ == revision_);
    const bool currentView =
        geometryRevision_ == revision_ &&
        (representation_->currentIndex() != 1 || !mesh_ || meshRevision_ == revision_);
    // A historical mesh can reuse numeric IDs from a different CAD source.
    // It must not inherit the current document's face colors.
    auto colors = activeFaceColors();
    if (!cadMode()) {
        colors.clear();
        for (const auto &instance : assembly_.instances)
            for (const auto &[id, color] : draftFaceColors_)
                colors.emplace(instanceFaceId(instance.id, id), color);
        if (assembly_.includeCad)
            for (const auto &[id, color] : cadFaceColors_)
                colors.emplace(id + assemblyCadFaceOffset, color);
    }
    viewport_->setFaceColors(currentView ? colors : FaceColors{});
    const auto simulationState = simulationDock_->archive();
    const auto baffle = simulationState.settings.infiniteBaffle;
    showAperture_->setEnabled(baffle.has_value());
    std::set<std::uint32_t> interfaceFaces;
    for (const auto *element : lem::bemPorts(lemEditor_->network())) {
        const auto &port = std::get<lem::BemPort>(element->component);
        if (port.kind == lem::BemPortKind::interfacePort)
            interfaceFaces.insert(port.faceIds.begin(), port.faceIds.end());
    }
    std::set<std::uint32_t> beInterfaces;
    for (const auto &interface : simulationState.settings.acousticModel.interfaces)
        beInterfaces.insert(interface.faceIds.begin(), interface.faceIds.end());
    const bool computedInterfaces = simulationState.current && simulationState.request &&
                                    simulationState.request->subdomains.has_value();
    showInterfaces_->setEnabled(!interfaceFaces.empty() || !beInterfaces.empty() ||
                                computedInterfaces);
    std::set<std::uint32_t> hidden;
    if (baffle && !showAperture_->isChecked())
        hidden.insert(baffle->apertureFaceIds.begin(), baffle->apertureFaceIds.end());
    if (!showInterfaces_->isChecked())
        hidden.insert(interfaceFaces.begin(), interfaceFaces.end());
    // Transparent BE interfaces are drawn separately, never as opaque material.
    hidden.insert(beInterfaces.begin(), beInterfaces.end());
    viewport_->setHiddenFaces({hidden.begin(), hidden.end()});
    if (cadMode() && cadSnapshot_)
        simulationDock_->setAcousticAxis(cadSnapshot_->acousticOrigin, cadSnapshot_->acousticAxis);
    else if (!cadMode() && geometry_) {
        const auto frame = assemblyAcousticFrame(assembly_, geometry_->length());
        simulationDock_->setAcousticAxis(frame.origin, frame.forward);
    }
    viewport_->setObservationFrame(simulationDock_->observationFrame());
    const auto draftPreview = !cadMode() ? rawPreview_ : nullptr;
    const auto overlay =
        !cadMode() && !isPlacedAssembly(assembly_) && sectionMesh_->isChecked() && mesh_
            ? std::shared_ptr<const SurfaceMesh>(mesh_, &mesh_->mesh)
            : nullptr;
    horizontalSection_->setMeshes(draftPreview, overlay);
    verticalSection_->setMeshes(draftPreview, overlay);
    sectionMesh_->setEnabled(!cadMode() && static_cast<bool>(mesh_));
    sectionMesh_->setText(mesh_
                              ? tr("Netz · Rev. %1%2")
                                    .arg(meshRevision_)
                                    .arg(meshRevision_ == revision_ ? QString{} : tr(" · veraltet"))
                              : tr("Netz"));
    if (representation_->currentIndex() == 1 && mesh_) {
        viewport_->setMesh(std::shared_ptr<const SurfaceMesh>(mesh_, &mesh_->mesh));
        viewport_->setLabel(tr("Netz · Rev. %1%2")
                                .arg(meshRevision_)
                                .arg(meshRevision_ == revision_ ? "" : tr(" · veraltet")));
    } else {
        viewport_->setMesh(preview_);
        viewport_->setLabel(tr("Vorschau · Rev. %1%2")
                                .arg(geometryRevision_)
                                .arg(geometryRevision_ == revision_ ? "" : tr(" · veraltet")));
    }
    viewport_->setSelectedFaces(
        placedSelection(cadMode() ? cadDock_->selectedFaces() : meshPanel_->selectedFaces()));
    updateAcousticDomainDisplay(simulationState);
    updateInterfacePresentation();
    updateFaceColorActions();
}

FaceColors &MainWindow::activeFaceColors() {
    return cadMode() ? cadFaceColors_ : draftFaceColors_;
}

std::vector<std::uint32_t> MainWindow::faceColorSelection() const {
    if (!preview_ || geometryRevision_ != revision_ ||
        (representation_->currentIndex() == 1 && mesh_ && meshRevision_ != revision_))
        return {};
    if (cadMode()) {
        const auto current = cadDock_->snapshot();
        if (!current || !cadSnapshot_ || cadDock_->isGeometryPending() ||
            current->documentId != cadSnapshot_->documentId ||
            current->sessionId != cadSnapshot_->sessionId ||
            current->revision != cadSnapshot_->revision)
            return {};
    }
    auto selected = viewport_->selectedFaces();
    std::set<std::uint32_t> known;
    for (const auto &triangle : preview_->triangles)
        known.insert(triangle.faceId);
    if (std::any_of(selected.begin(), selected.end(), [&](auto id) { return !known.contains(id); }))
        return {};
    return selected;
}

void MainWindow::updateFaceColorActions() {
    const auto selected = faceColorSelection();
    lemEditor_->setSelectedFaces(selected);
    assignFaceColor_->setEnabled(!selected.empty());
    const auto &colors = activeFaceColors();
    resetFaceColors_->setEnabled(std::any_of(selected.begin(), selected.end(),
                                             [&](auto id) { return colors.contains(id); }));
}

void MainWindow::chooseFaceColor() {
    const auto selected = faceColorSelection();
    if (selected.empty()) {
        setDiagnostic(tr("Zuerst Flächen der aktuellen Ansicht auswählen."));
        return;
    }
    const auto currentRevision = revision_;
    const auto currentEpoch = projectEpoch_;
    const auto &colors = activeFaceColors();
    const auto first = colors.find(selected.front());
    const auto color = first == colors.end() ? FaceColor{255, 145, 0} : first->second;
    const auto chosen = QColorDialog::getColor(QColor(color.red, color.green, color.blue), this,
                                               tr("Flächenfarbe wählen"));
    if (!chosen.isValid())
        return;
    if (currentEpoch != projectEpoch_ || currentRevision != revision_ ||
        selected != faceColorSelection()) {
        setDiagnostic(tr("Flächenauswahl oder Geometrie geändert · Farbe nicht übernommen."));
        return;
    }
    setSelectedFaceColor(FaceColor{static_cast<std::uint8_t>(chosen.red()),
                                   static_cast<std::uint8_t>(chosen.green()),
                                   static_cast<std::uint8_t>(chosen.blue())});
}

void MainWindow::setSelectedFaceColor(std::optional<FaceColor> color) {
    const auto selected = faceColorSelection();
    if (selected.empty()) {
        setDiagnostic(tr("Zuerst Flächen der aktuellen Ansicht auswählen."));
        return;
    }
    auto &colors = activeFaceColors();
    bool changed = false;
    for (const auto id : selected) {
        if (color) {
            const auto found = colors.find(id);
            changed |= found == colors.end() || found->second != *color;
            colors[id] = *color;
        } else
            changed |= colors.erase(id) != 0;
    }
    meshColor_->setCurrentIndex(static_cast<int>(MeshColorMode::faceColors));
    if (changed)
        markDirty();
    updateDisplay();
}

void MainWindow::updatePlots(const GeometryAnalysis &analysis) {
    PlotSeries volume{{}, QColor(28, 102, 180), {}}, area = volume, derivative = volume;
    const bool asymmetric = geometry_ && geometry_->parameters().asymmetricEnabled;
    PlotSeries h{asymmetric ? tr("Rechts") : tr("H"), QColor(28, 102, 180), {}},
        v{asymmetric ? tr("Oben") : tr("V"), QColor(205, 112, 33), {}},
        left{tr("Links"), QColor(32, 140, 104), {}}, bottom{tr("Unten"), QColor(146, 77, 166), {}};
    for (const auto &s : analysis.samples) {
        const double z = s.z * 1000;
        volume.points.emplace_back(z, s.volume * 1e6);
        area.points.emplace_back(z, s.area * 1e6);
        derivative.points.emplace_back(z, s.areaDerivative * 1000);
        h.points.emplace_back(z, s.horizontalWallAngle * 180 / std::numbers::pi);
        v.points.emplace_back(z, s.verticalWallAngle * 180 / std::numbers::pi);
        if (asymmetric) {
            left.points.emplace_back(z, s.leftWallAngle * 180 / std::numbers::pi);
            bottom.points.emplace_back(z, s.bottomWallAngle * 180 / std::numbers::pi);
        }
    }
    volumePlot_->setSeries({std::move(volume)});
    areaPlot_->setSeries({std::move(area)});
    derivativePlot_->setSeries({std::move(derivative)});
    if (asymmetric)
        anglePlot_->setSeries({std::move(h), std::move(left), std::move(v), std::move(bottom)});
    else
        anglePlot_->setSeries({std::move(h), std::move(v)});
}

bool MainWindow::geometryReady() const noexcept {
    return meshSettingsValid_ && preview_ && rawPreview_ && geometryRevision_ == revision_ &&
           !geometryPending_ &&
           (cadMode() ? cadDock_->isSourceAvailable() && cadSnapshot_ &&
                            !cadSnapshot_->hasUnresolvedAssignments()
                      : geometry_ && (!assembly_.includeCad ||
                                      (cadSnapshot_ && !cadSnapshot_->hasUnresolvedAssignments())));
}
bool MainWindow::meshReady() const noexcept {
    return mesh_ && meshRevision_ == revision_ && !meshBusy_;
}
bool MainWindow::meshInProgress() const noexcept {
    return meshBusy_;
}
QString MainWindow::lastMeshError() const {
    return meshError_;
}
std::uint64_t MainWindow::revision() const noexcept {
    return revision_;
}
void MainWindow::setDiagnostic(const QString &message, bool error) {
    diagnostic_->setText(message.left(95));
    diagnostic_->setToolTip(message);
    diagnostic_->setStyleSheet(error ? "color: #a72e2b" : "");
}
void MainWindow::openParameterStudy(bool convergence, const QString &manifestPath) {
    if (isPlacedAssembly(assembly_)) {
        setDiagnostic(tr(
            "Baugruppen-Sweeps sind nicht verfügbar; die gemeinsame Simulation bleibt nutzbar."));
        return;
    }
    commitFocusedEdit();
    flushEdit();
    try {
        auto baseline = projectState();
        const auto studyEpoch = projectEpoch_;
        if (manifestPath.isEmpty()) {
            if (!convergence && cadMode())
                throw std::invalid_argument("Parameter-Sweeps sind im Draft-Modus verfügbar. Für "
                                            "CAD Netzkonvergenz wählen.");
            if (cadMode() && (cadDock_->isGeometryPending() || !baseline.cad ||
                              baseline.cad->hasUnresolvedAssignments()))
                throw std::invalid_argument("Bestätigter CAD-Geometriestand erforderlich.");
            if (!cadMode())
                (void)DraftGeometry(baseline.parameters);
        }
        std::vector<double> frequencies;
        if (convergence) {
            const auto &s = baseline.simulation.settings;
            const auto suggestion = QString("%1; %2; %3")
                                        .arg(s.minimumFrequency)
                                        .arg(std::sqrt(s.minimumFrequency * s.maximumFrequency))
                                        .arg(s.maximumFrequency);
            bool ok = false;
            const auto text = QInputDialog::getText(this, tr("Netzkonvergenz"),
                                                    tr("Frequenzen [Hz], mit Strichpunkt getrennt"),
                                                    QLineEdit::Normal, suggestion, &ok);
            if (!ok)
                return;
            for (auto token : text.split(QRegularExpression("[;\\s]+"), Qt::SkipEmptyParts)) {
                token.replace(',', '.');
                bool valid = false;
                const auto value = token.toDouble(&valid);
                if (!valid || !std::isfinite(value) || value < 1 || value > 100000)
                    throw std::invalid_argument(
                        "Frequenzen müssen zwischen 1 und 100000 Hz liegen.");
                frequencies.push_back(value);
            }
            std::sort(frequencies.begin(), frequencies.end());
            frequencies.erase(std::unique(frequencies.begin(), frequencies.end()),
                              frequencies.end());
            if (frequencies.empty())
                throw std::invalid_argument("Mindestens eine Prüffrequenz erforderlich.");
        }
        ParameterSweepDialog study(baseline, this);
        QString analysisPath;
        study.onVariantAnalysisRequested = [&analysisPath](const QString &path) {
            analysisPath = path;
        };
        study.onVariantChosen = [this, studyEpoch](Project variant) {
            if (projectEpoch_ != studyEpoch)
                throw std::invalid_argument("Projekt wurde gewechselt · Studie erneut öffnen.");
            adoptStudyVariant(std::move(variant));
        };
        if (convergence)
            study.configureConvergence(std::move(frequencies));
        if (!manifestPath.isEmpty())
            study.openStudy(manifestPath);
        struct GpuHandoff {
            QTimer *timer{};
            bool reserved{};
            std::function<void(bool)> ready;
        };
        auto handoff = std::make_shared<GpuHandoff>();
        handoff->timer = new QTimer(&study);
        handoff->timer->setObjectName("studyGpuHandoffTimer");
        handoff->timer->setInterval(30);
        const auto release = [this, handoff] {
            handoff->timer->stop();
            handoff->ready = {};
            if (handoff->reserved) {
                handoff->reserved = false;
                simulationDock_->endExternalJob();
            }
        };
        connect(handoff->timer, &QTimer::timeout, &study, [this, handoff, release, studyEpoch] {
            if (projectEpoch_ != studyEpoch) {
                auto ready = std::move(handoff->ready);
                release();
                setDiagnostic(tr("Projekt wurde gewechselt · Studie erneut öffnen."), true);
                if (ready)
                    ready(false);
                return;
            }
            if (!handoff->reserved || !simulationDock_->externalJobReady())
                return;
            handoff->timer->stop();
            auto ready = std::move(handoff->ready);
            if (ready)
                ready(true);
        });
        study.gpuWaitMessage = [this] {
            return simulationDock_->isBusy() ? tr("Laufende Berechnung wird beendet …")
                                             : tr("GPU-Speicher wird freigegeben …");
        };
        study.onAcquireGpu = [this, handoff, release, studyEpoch](std::function<void(bool)> ready) {
            if (projectEpoch_ != studyEpoch) {
                setDiagnostic(tr("Projekt wurde gewechselt · Studie erneut öffnen."), true);
                ready(false);
                return;
            }
            handoff->ready = std::move(ready);
            handoff->reserved = true;
            try {
                simulationDock_->beginExternalJob();
                handoff->timer->start();
            } catch (...) {
                release();
                throw;
            }
        };
        study.onCancelGpuAcquire = release;
        study.onReleaseGpu = release;
        study.exec();
        if (!analysisPath.isEmpty()) {
            if (projectEpoch_ != studyEpoch) {
                setDiagnostic(tr("Projekt wurde gewechselt · Studie erneut öffnen."), true);
                return;
            }
            confirmUnsaved([this, analysisPath, studyEpoch] {
                if (projectEpoch_ != studyEpoch) {
                    setDiagnostic(tr("Projekt wurde gewechselt · Studie erneut öffnen."), true);
                    return;
                }
                fieldCameraInitialized_ = false;
                queueSourceImport(analysisPath, false, true);
            });
        }
    } catch (const std::exception &error) {
        QMessageBox::warning(this, tr("Studie"), QString::fromUtf8(error.what()));
    }
}

void MainWindow::adoptStudyVariant(Project variant) {
    if ((variant.source == ProjectSource::freecad) != cadMode())
        throw std::invalid_argument("Zum Übernehmen in den passenden Arbeitsbereich wechseln.");
    const auto before = editState();
    if (cadMode()) {
        const auto current = cadDock_->snapshot();
        if (!current || !variant.cad || current->documentId != variant.cad->documentId ||
            current->sessionId != variant.cad->sessionId ||
            current->revision != variant.cad->revision ||
            current->faces.size() != variant.cad->faces.size())
            throw std::invalid_argument(
                "Der CAD-Geometriestand der Variante ist nicht mehr aktuell.");
        for (const auto &face : variant.cad->faces) {
            const auto found = std::find_if(current->faces.begin(), current->faces.end(),
                                            [&](const auto &c) { return c.id == face.id; });
            if (found == current->faces.end() || found->objectId != face.objectId ||
                found->boundary != face.boundary || found->normalVelocity != face.normalVelocity ||
                found->state != face.state)
                throw std::invalid_argument("CAD-Zuordnungen wurden seit der Studie geändert.");
        }
    }
    auto after = before;
    if (!cadMode())
        after.parameters = variant.parameters;
    after.meshing = variant.meshing;
    if (cadMode())
        after.cadFaceTargets = variant.cadFaceTargetEdgeLengths;
    else
        after.draftFaceTargets = variant.draftFaceTargetEdgeLengths;
    applyEditState(after);
    after = editState();
    historyState_ = after;
    if (!sameEditState(before, after))
        undoStack_->push(
            new EditCommand(before, after, [this](const EditState &s) { applyEditState(s); }));
    // Adoption edits the current document. Archived study data stays visibly stale until a new
    // solve.
    variant.simulation.current = false;
    simulationDock_->restoreArchive(variant.simulation);
    markDirty();
    updateEditActions();
    setDiagnostic(tr("Variante übernommen · Netz und Akustik neu berechnen."));
}

void MainWindow::updateTitle() {
    const QString name =
        projectPath_.isEmpty() ? tr("Unbenannt") : QFileInfo(projectPath_).fileName();
    setWindowTitle(name + (dirty_ ? " *" : "") + " — HornSculpt " +
                   QCoreApplication::applicationVersion());
}

void MainWindow::markDirty() {
    if (!loading_) {
        ++changeSerial_;
        dirty_ = true;
        updateTitle();
    }
}

void MainWindow::updateMeshStatus() {
    if (!mesh_) {
        qualityLabel_->setText(tr("Netz: nicht erzeugt"));
        qualityLabel_->setToolTip({});
        meshItem_->setText(0, tr("Netz · nicht erzeugt"));
        return;
    }
    qualityLabel_->setText(
        tr("Netz Rev. %1 · %2 Dreiecke%3 · kleinste Kante %4 mm")
            .arg(meshRevision_)
            .arg(mesh_->mesh.triangles.size())
            .arg(meshRevision_ != revision_
                     ? tr(" · veraltet")
                     : (gmsh_ && !simulationDock_->archive().settings.acousticModel.domains.empty()
                            ? tr(" · Gebiete geprüft")
                            : tr(" · geschlossen")))
            .arg(locale().toString(mesh_->quality.minimumEdgeLength * 1000, 'f', 3)));
    const auto &q = mesh_->quality;
    qualityLabel_->setToolTip(
        tr("Randkanten: %1\nNicht-manifold Kanten: %2\nMin. Winkel: %3°\nMax. "
           "Seitenverhältnis: %4\nKanten min / median / max: %5 / %6 / %7 "
           "mm\nNicht-manifold Knoten: %8\nInkonsistente Kanten: %9\nSelbstschnittpaare: %10")
            .arg(q.boundaryEdges)
            .arg(q.nonManifoldEdges)
            .arg(q.minimumAngleDegrees, 0, 'f', 1)
            .arg(q.maximumAspectRatio, 0, 'f', 2)
            .arg(q.minimumEdgeLength * 1000, 0, 'f', 2)
            .arg(q.medianEdgeLength * 1000, 0, 'f', 2)
            .arg(q.maximumEdgeLength * 1000, 0, 'f', 2)
            .arg(q.nonManifoldVertices)
            .arg(q.inconsistentEdges)
            .arg(q.selfIntersectionPairs));
    meshItem_->setText(0, tr("Netz · %1 Dreiecke").arg(mesh_->mesh.triangles.size()));
}

Project MainWindow::projectState() const {
    Project project;
    project.parameters = parameters_->parameters();
    project.assembly = assembly_;
    project.meshing = meshPanel_->settings();
    project.draftFaceTargetEdgeLengths = draftFaceTargetEdgeLengths_;
    project.cadFaceTargetEdgeLengths = cadFaceTargetEdgeLengths_;
    project.draftFaceColors = draftFaceColors_;
    project.cadFaceColors = cadFaceColors_;
    project.draftReversedFaceIds = draftReversedFaceIds_;
    project.cadReversedFaceIds = cadReversedFaceIds_;
    project.source = cadMode() ? ProjectSource::freecad : ProjectSource::draft;
    project.revision = revision_;
    project.geometryRevision = geometryRevision_;
    project.meshRevision = meshRevision_;
    project.preview = preview_;
    project.mesh = mesh_;
    project.cad = cadDock_->snapshot();
    if (project.cad && project.cad->sessionId == "gmsh-import-offline")
        project.gmsh = gmsh_;
    project.simulation = simulationDock_->archive();
    project.comparisonSnapshots = comparisonSnapshots_;
    return project;
}

MainWindow::PreparedProject MainWindow::prepareProject(Project project, bool verifyDraft,
                                                       std::stop_token stop) {
    const auto checkStop = [&] {
        if (stop.stop_requested())
            throw std::runtime_error("Import abgebrochen.");
    };
    checkStop();
    PreparedProject prepared;
    if (project.source == ProjectSource::draft && isPlacedAssembly(project.assembly)) {
        prepared.cad = PreparedCadSnapshot::prepare(project.cad);
        project.cad = prepared.cad.snapshot();
        prepared.geometry = std::make_shared<DraftGeometry>(project.parameters);
        prepared.rawPreview = std::make_shared<SurfaceMesh>(
            prepared.geometry->previewMesh(project.meshing.adaptive.geometricTolerance, stop));
        auto composed =
            std::make_shared<SurfaceMesh>(assemblyPreview(*prepared.rawPreview, project));
        if (project.preview && project.geometryRevision == project.revision &&
            *project.preview != *composed) {
            ++project.revision;
            project.simulation.current = false;
            prepared.rebuiltDraft = true;
        }
        project.preview = std::move(composed);
        project.geometryRevision = project.revision;
        prepared.analysis = analyzeGeometry(*prepared.geometry);
        prepared.project = std::move(project);
        return prepared;
    }
    if (project.source == ProjectSource::draft) {
        prepared.geometry = std::make_shared<DraftGeometry>(project.parameters);
        if (verifyDraft && project.preview && project.geometryRevision == project.revision) {
            prepared.rawPreview = std::make_shared<SurfaceMesh>(
                prepared.geometry->previewMesh(project.meshing.adaptive.geometricTolerance));
            auto current = std::make_shared<SurfaceMesh>(
                applyFaceOrientation(*prepared.rawPreview, project.draftReversedFaceIds));
            if (*current != *project.preview) {
                if (project.revision >= std::numeric_limits<std::uint64_t>::max() - 1)
                    throw std::runtime_error("Projektrevision kann nicht erneuert werden.");
                ++project.revision;
                project.geometryRevision = project.revision;
                project.preview = std::move(current);
                project.simulation.current = false;
                prepared.rebuiltDraft = true;
            }
        }
        checkStop();
        prepared.analysis = analyzeGeometry(*prepared.geometry);
        if (!project.preview || project.geometryRevision != project.revision) {
            prepared.rawPreview = std::make_shared<SurfaceMesh>(
                prepared.geometry->previewMesh(project.meshing.adaptive.geometricTolerance));
            project.preview = std::make_shared<SurfaceMesh>(
                applyFaceOrientation(*prepared.rawPreview, project.draftReversedFaceIds));
            project.geometryRevision = project.revision;
        } else if (!prepared.rawPreview)
            prepared.rawPreview = std::make_shared<SurfaceMesh>(
                applyFaceOrientation(*project.preview, project.draftReversedFaceIds));
    }
    checkStop();
    prepared.cad = PreparedCadSnapshot::prepare(project.cad);
    project.cad = prepared.cad.snapshot();
    if (project.source == ProjectSource::freecad && project.cad)
        prepared.rawPreview = std::shared_ptr<const SurfaceMesh>(project.cad, &project.cad->mesh);
    checkStop();
    prepared.project = std::move(project);
    return prepared;
}

void MainWindow::applyProject(const Project &project) {
    applyProject(prepareProject(project, false));
}

void MainWindow::applyProject(const PreparedProject &prepared) {
    const auto &project = prepared.project;
    if (project.gmsh && cadDock_->hasLiveLink())
        throw std::invalid_argument("Gmsh-Quelle kann bei aktiver FreeCAD-Verbindung nicht "
                                    "geöffnet werden. Verbindung zuerst trennen.");
    // Widgets remain on their Qt thread. Validate all control domains before mutation;
    // the owned CAD snapshot and numerical geometry were already prepared separately.
    ParameterPanel checkedParameters;
    checkedParameters.setParameters(project.parameters);
    MeshPanel checkedMesh;
    checkedMesh.setSettings(project.meshing);
    simulationDock_->validateArchive(project.simulation);
    const QSignalBlocker mode(sourceMode_);
    loading_ = true;
    ++projectEpoch_;
    comparisonSnapshots_ = project.comparisonSnapshots;
    debounce_->stop();
    geometryWorker_.request_stop();
    meshWorker_.request_stop();
    geometryBusy_ = geometryPending_ = meshBusy_ = false;
    parameters_->setParameters(project.parameters);
    assembly_ = project.assembly;
    meshPanel_->setSettings(project.meshing);
    meshing_ = project.meshing;
    meshSettingsValid_ = true;
    draftFaceTargetEdgeLengths_ = project.draftFaceTargetEdgeLengths;
    cadFaceTargetEdgeLengths_ = project.cadFaceTargetEdgeLengths;
    draftFaceColors_ = project.draftFaceColors;
    cadFaceColors_ = project.cadFaceColors;
    draftReversedFaceIds_ = project.draftReversedFaceIds;
    cadReversedFaceIds_ = project.cadReversedFaceIds;
    cadDock_->restoreSnapshot(prepared.cad);
    sourceMode_->setCurrentIndex(project.source == ProjectSource::freecad ? 1 : 0);
    revision_ = project.revision;
    geometryRevision_ = project.geometryRevision;
    analysisRevision_ = prepared.analysis.samples.empty() ? 0 : geometryRevision_;
    meshRevision_ = project.meshRevision;
    geometry_ = prepared.geometry;
    parameters_->setSlotTransitionEnd(geometry_ ? std::optional{geometry_->slotTransitionEnd()}
                                                : std::nullopt);
    preview_ = project.preview;
    rawPreview_ = prepared.rawPreview;
    mesh_ = project.mesh;
    cadSnapshot_ = project.cad;
    gmsh_ = project.gmsh;
    meshError_.clear();
    cancelAction_->setEnabled(false);
    viewport_->setSelectedFaces({});
    meshPanel_->setSelectedFaces({});
    meshPanel_->setFaceNotice({});
    {
        const QSignalBlocker visibility(showInterfaces_);
        showInterfaces_->setChecked(!project.simulation.settings.acousticModel.interfaces.empty());
    }
    updateSourcePresentation();
    updatePlots(prepared.analysis);
    updateMeshStatus();
    lemEditor_->setNetwork(project.simulation.settings.couplingNetwork);
    lemEditor_->setCouplingEnabled(project.simulation.settings.couplingEnabled);
    configureSimulation();
    simulationDock_->restoreArchive(project.simulation);
    configureLemFaces();
    updateInterfacePresentation();
    representation_->setCurrentIndex(meshReady() ? 1 : 0);
    meshColor_->setCurrentIndex(activeFaceColors().empty()
                                    ? static_cast<int>(MeshColorMode::boundary)
                                    : static_cast<int>(MeshColorMode::faceColors));
    revisionLabel_->setText(tr("Geometrie Rev. %1").arg(revision_));
    meshAction_->setEnabled(geometryReady());
    updateDisplay();
    if ((cadMode() && cadSnapshot_) || isPlacedAssembly(assembly_))
        viewport_->viewFront(simulationDock_->observationFrame().forward);
    else
        viewport_->fitView();
    initialModelFitPending_ = false;
    setDiagnostic({});
    dirty_ = false;
    loading_ = false;
    if (cadMode()) {
        if (!project.gmsh)
            cadDock_->ensureListening();
        simulationDock_->show();
        simulationDock_->raise();
    }
    updateTitle();
    resetEditHistory();
}

void MainWindow::adoptProject(const PreparedProject &prepared, const QString &path) {
    applyProject(prepared);
    if (recoveryStore_) {
        ++recoveryGeneration_;
        recoveryHasSaved_ = false;
        try {
            recoveryStore_->newSession();
        } catch (const std::exception &e) {
            setDiagnostic(QString::fromUtf8(e.what()), true);
        }
    }
    projectPath_ = path;
    if (prepared.rebuiltDraft) {
        markDirty();
        setDiagnostic(tr("Geometrie neu tesselliert; vorhandenes Netz ist veraltet."));
    }
    updateTitle();
}

void MainWindow::loadProject(const QString &path) {
    adoptProject(prepareProject(readProject(path), true), path);
}

bool MainWindow::sourceImportInProgress() const noexcept {
    return sourceImportBusy_;
}

QString MainWindow::sourceImportError() const {
    return sourceImportError_;
}

void MainWindow::openStep(const QString &path) {
    fieldCameraInitialized_ = false;
    queueSourceImport(path, true);
}

void MainWindow::openProject(const QString &path) {
    fieldCameraInitialized_ = false;
    queueSourceImport(path, false);
}

void MainWindow::queueSourceImport(const QString &path, bool step, bool openAsCopy) {
    if (closing_)
        return;
    if (recoveryOpening_ && path != recoveryOpening_->path)
        recoveryOpening_.reset();
    SourceImport request{path,      step,          ++sourceImportId_, projectEpoch_,
                         revision_, changeSerial_, meshing_};
    request.openAsCopy = openAsCopy;
    sourceImportError_.clear();
    if (sourceImportBusy_) {
        pendingSourceImport_ = std::move(request);
        sourceImportWorker_.request_stop();
        return;
    }
    startSourceImport(std::move(request));
}

void MainWindow::startSourceImport(SourceImport request) {
    try {
        auto progressOwner = std::make_unique<QProgressDialog>(
            request.step ? tr("STEP importieren …") : tr("Projekt öffnen …"), tr("Abbrechen"), 0, 0,
            this);
        auto *progress = progressOwner.get();
        progress->setObjectName(QStringLiteral("sourceImportProgress"));
        progress->setWindowModality(Qt::WindowModal);
        progress->setMinimumDuration(0);
        progress->setToolTip(tr("Abbrechen verwirft den Import. Ein laufender Projekt-Lesevorgang "
                                "muss zunächst zurückkehren."));
        progress->show();
        connect(progress, &QProgressDialog::canceled, this, [this] {
            ++sourceImportId_;
            pendingSourceImport_.reset();
            sourceImportWorker_.request_stop();
        });
        auto worker = sourceImportStarter_([this, request = std::move(request),
                                            progress](std::stop_token stop) {
            std::shared_ptr<const PreparedProject> prepared;
            QString error;
            try {
                Project project;
                if (request.step) {
                    auto snapshot = readStepSnapshot(
                        request.path, request.meshing.adaptive.geometricTolerance, stop);
                    project.source = ProjectSource::freecad;
                    project.meshing = request.meshing;
                    project.cad = snapshot;
                    project.preview = std::shared_ptr<const SurfaceMesh>(snapshot, &snapshot->mesh);
                    project.geometryRevision = project.revision;
                    if (request.appendToDraft) {
                        project = *request.appendToDraft;
                        project.cad = project.assembly.includeCad && project.cad
                                          ? appendAssemblyCad(*project.cad, *snapshot)
                                          : snapshot;
                        if (!project.assembly.includeCad) {
                            project.cadFaceTargetEdgeLengths.clear();
                            project.cadReversedFaceIds.clear();
                            project.cadFaceColors.clear();
                        }
                        project.assembly.includeCad = true;
                        project.source = ProjectSource::draft;
                        ++project.revision;
                        project.preview.reset();
                        project.mesh.reset();
                        project.geometryRevision = project.meshRevision = 0;
                        project.simulation.current = false;
                    }
                } else
                    project = readProject(request.path);
                prepared = std::make_shared<PreparedProject>(
                    prepareProject(std::move(project), !request.step, stop));
            } catch (const std::exception &e) {
                error = QString::fromUtf8(e.what());
            }
            QMetaObject::invokeMethod(
                this,
                [this, request, progress, stop, prepared = std::move(prepared),
                 error = std::move(error)] {
                    if (sourceImportWorker_.joinable())
                        sourceImportWorker_.join();
                    progress->hide();
                    progress->deleteLater();
                    sourceImportBusy_ = false;
                    // No next worker is assigned while the preceding read is still
                    // running. A pending request also retains the state against which
                    // it was issued.
                    if (pendingSourceImport_) {
                        auto next = std::move(*pendingSourceImport_);
                        pendingSourceImport_.reset();
                        if (!closing_ && next.epoch == projectEpoch_ &&
                            next.revision == revision_ && next.changeSerial == changeSerial_)
                            startSourceImport(std::move(next));
                        else if (!closing_ && next.epoch == projectEpoch_ &&
                                 next.revision == revision_ && next.changeSerial != changeSerial_) {
                            sourceImportError_ =
                                tr("Import verworfen: Das aktuelle Projekt wurde während "
                                   "des Lesens geändert.");
                            setDiagnostic(sourceImportError_, true);
                        }
                        return;
                    }
                    if (closing_ || stop.stop_requested() || request.id != sourceImportId_ ||
                        request.epoch != projectEpoch_ || request.revision != revision_) {
                        if (recoveryOpening_ && request.path == recoveryOpening_->path)
                            recoveryOpening_.reset();
                        return;
                    }
                    if (request.changeSerial != changeSerial_) {
                        if (recoveryOpening_ && request.path == recoveryOpening_->path)
                            recoveryOpening_.reset();
                        sourceImportError_ = tr("Import verworfen: Das aktuelle Projekt "
                                                "wurde während des Lesens geändert.");
                        setDiagnostic(sourceImportError_, true);
                        return;
                    }
                    sourceImportError_ = error;
                    if (sourceImportError_.isEmpty()) {
                        try {
                            const bool restoring =
                                recoveryOpening_ && request.path == recoveryOpening_->path;
                            adoptProject(*prepared, restoring || request.openAsCopy ? QString{}
                                                    : request.appendToDraft         ? projectPath_
                                                    : request.step                  ? QString{}
                                                                                    : request.path);
                            if (restoring) {
                                recoveredEntry_ = std::move(recoveryOpening_);
                                markDirty();
                                setDiagnostic(
                                    tr("Wiederhergestellt · zum Erhalt als neue Datei speichern."));
                            } else if (request.step || request.openAsCopy)
                                markDirty();
                        } catch (const std::exception &e) {
                            sourceImportError_ = QString::fromUtf8(e.what());
                        }
                    }
                    if (!sourceImportError_.isEmpty())
                        setDiagnostic(sourceImportError_, true);
                    if (recoveryOpening_ && request.path == recoveryOpening_->path &&
                        !sourceImportError_.isEmpty())
                        recoveryOpening_.reset(); // Keep the failed recovery file for diagnosis.
                },
                Qt::QueuedConnection);
        });
        if (!worker.joinable())
            throw std::runtime_error("Import-Thread konnte nicht gestartet werden");
        sourceImportWorker_ = std::move(worker);
        sourceImportBusy_ = true;
        progressOwner.release();
    } catch (const std::exception &e) {
        sourceImportBusy_ = false;
        sourceImportError_ =
            tr("Import konnte nicht gestartet werden: ") + QString::fromUtf8(e.what());
        setDiagnostic(sourceImportError_, true);
    }
}

bool MainWindow::saveInProgress() const noexcept {
    return saveWorker_->busy();
}

void MainWindow::setSaveActionsEnabled(bool enabled) {
    for (auto *action : saveBlockedActions_)
        action->setEnabled(enabled);
}

void MainWindow::startRecovery() {
    if (closing_ || recoveryStore_)
        return;
    try {
        if (recoveryRoot_.isEmpty()) {
            const QString base =
                QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation);
            if (base.isEmpty())
                throw std::runtime_error("Kein Wiederherstellungsordner verfügbar.");
            recoveryRoot_ = QDir(base).filePath("recovery");
        }
        recoveryStore_ = std::make_unique<RecoveryStore>(recoveryRoot_);
        recoveryTimer_ = new QTimer(this);
        recoveryTimer_->setObjectName("recoveryTimer");
        recoveryTimer_->setInterval(60000);
        connect(recoveryTimer_, &QTimer::timeout, this, &MainWindow::saveRecovery);
        recoveryTimer_->start();
        offerRecovery();
    } catch (const std::exception &e) {
        setDiagnostic(
            tr("Automatische Sicherung nicht verfügbar: %1").arg(QString::fromUtf8(e.what())),
            true);
    }
}

void MainWindow::offerRecovery() {
    if (closing_ || !recoveryStore_)
        return;
    if (sourceImportBusy_) {
        QTimer::singleShot(250, this, &MainWindow::offerRecovery);
        return;
    }
    auto entries = RecoveryStore::available(recoveryRoot_);
    for (auto &entry : entries) {
        QMessageBox question(
            QMessageBox::Question, tr("Wiederherstellung"),
            tr("Eine nicht regulär beendete Sicherung vom %1 wurde gefunden.\n"
               "Ursprüngliches Projekt: %2\n\n"
               "Wiederhergestellte Daten werden zunächst ohne Dateinamen geöffnet; "
               "das Original wird nicht überschrieben.")
                .arg(locale().toString(entry.modified, QLocale::ShortFormat),
                     entry.originalPath.isEmpty() ? tr("Unbenannt") : entry.originalPath),
            QMessageBox::NoButton, this);
        auto *restore = question.addButton(tr("Wiederherstellen"), QMessageBox::AcceptRole);
        auto *discard = question.addButton(tr("Verwerfen"), QMessageBox::DestructiveRole);
        question.addButton(tr("Später"), QMessageBox::RejectRole);
        question.exec();
        if (question.clickedButton() == discard) {
            RecoveryStore::discard(entry);
            continue;
        }
        if (question.clickedButton() == restore) {
            auto selected = std::make_shared<RecoveryStore::Entry>(std::move(entry));
            confirmUnsaved([this, selected] {
                recoveryOpening_ = selected;
                openProject(selected->path);
            });
        }
        return;
    }
}

void MainWindow::saveRecovery() {
    if (closing_ || !dirty_ || !recoveryStore_ || recoveryStore_->busy())
        return;
    try {
        auto snapshot = projectState();
        const auto serial = changeSerial_, generation = recoveryGeneration_;
        const auto result = snapshot.simulation.result;
        const auto complete = snapshot.simulation.resultComplete;
        if (recoveryHasSaved_ && serial == recoverySavedSerial_ &&
            recoverySavedResult_.lock() == result && recoverySavedComplete_ == complete)
            return;
        recoveryStore_->save(
            std::move(snapshot), projectPath_,
            [this, serial, generation, result, complete](QString error) {
                if (!closing_ && error.isEmpty() && generation == recoveryGeneration_) {
                    recoverySavedSerial_ = serial;
                    recoverySavedResult_ = result;
                    recoverySavedComplete_ = complete;
                    recoveryHasSaved_ = true;
                } else if (!closing_ && !error.isEmpty())
                    setDiagnostic(tr("Automatische Sicherung fehlgeschlagen: %1").arg(error), true);
            });
    } catch (const std::exception &e) {
        setDiagnostic(
            tr("Automatische Sicherung fehlgeschlagen: %1").arg(QString::fromUtf8(e.what())), true);
    }
}

void MainWindow::clearRecovery(bool newSession) {
    ++recoveryGeneration_;
    recoveryHasSaved_ = false;
    if (recoveredEntry_) {
        RecoveryStore::discard(*recoveredEntry_);
        recoveredEntry_.reset();
    }
    if (!recoveryStore_)
        return;
    try {
        if (newSession)
            recoveryStore_->newSession();
        else
            recoveryStore_->retire();
    } catch (const std::exception &e) {
        setDiagnostic(
            tr("Automatische Sicherung nicht verfügbar: %1").arg(QString::fromUtf8(e.what())),
            true);
    }
}

void MainWindow::saveProject(bool choosePath, std::function<void()> continuation) {
    if (saveWorker_->busy()) {
        setDiagnostic(tr("Speichern läuft bereits · bitte Abschluss abwarten."));
        return;
    }
    commitFocusedEdit();
    QString path = projectPath_;
    if (path.isEmpty() || choosePath)
        path = QFileDialog::getSaveFileName(this, tr("Projekt speichern"), path,
                                            tr("HornSculpt (*.hsproj)"));
    if (path.isEmpty())
        return;
    if (!path.endsWith(".hsproj", Qt::CaseInsensitive))
        path += ".hsproj";
    try {
        auto snapshot = projectState();
        const auto epoch = projectEpoch_, revision = revision_, serial = changeSerial_;
        const auto result = snapshot.simulation.result;
        const auto complete = snapshot.simulation.resultComplete;
        afterSave_ = std::move(continuation);
        setSaveActionsEnabled(false);
        setDiagnostic(tr("Projekt wird gespeichert …"));
        const bool started = saveWorker_->start(
            path, std::move(snapshot),
            [this, path, epoch, revision, serial, result, complete](QString error) {
                setSaveActionsEnabled(true);
                auto continuation = std::move(afterSave_);
                afterSave_ = {};
                if (!error.isEmpty()) {
                    setDiagnostic(tr("Speichern fehlgeschlagen"), true);
                    QMessageBox::critical(this, tr("Projekt speichern"), error);
                    return;
                }
                if (epoch != projectEpoch_) {
                    setDiagnostic(tr("Gespeichert · Projekt inzwischen gewechselt"));
                    return;
                }
                projectPath_ = path;
                const auto current = simulationDock_->archive();
                const bool unchanged = revision == revision_ && serial == changeSerial_ &&
                                       current.result == result &&
                                       current.resultComplete == complete;
                dirty_ = !unchanged;
                if (unchanged)
                    clearRecovery();
                updateTitle();
                setDiagnostic(unchanged ? tr("Projekt gespeichert")
                                        : tr("Projekt gespeichert · neuere Änderungen offen"));
                if (continuation)
                    confirmUnsaved(std::move(continuation));
            });
        if (!started) {
            afterSave_ = {};
            setSaveActionsEnabled(true);
            setDiagnostic(tr("Speichern läuft bereits"), true);
        }
    } catch (const std::exception &e) {
        afterSave_ = {};
        setSaveActionsEnabled(true);
        QMessageBox::critical(this, tr("Projekt speichern"), QString::fromUtf8(e.what()));
    }
}

void MainWindow::confirmUnsaved(std::function<void()> continuation) {
    commitFocusedEdit();
    if (saveWorker_->busy()) {
        if (!afterSave_)
            afterSave_ = std::move(continuation);
        else
            setDiagnostic(tr("Speichern läuft · eine Aktion wartet bereits."));
        return;
    }
    if (!dirty_) {
        continuation();
        return;
    }
    const auto answer = QMessageBox::question(
        this, tr("Projekt speichern"), tr("Änderungen speichern?"),
        QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel, QMessageBox::Save);
    if (answer == QMessageBox::Discard) {
        clearRecovery();
        continuation();
    } else if (answer == QMessageBox::Save)
        saveProject(false, std::move(continuation));
}
void MainWindow::newProject() {
    confirmUnsaved([this] {
        fieldCameraInitialized_ = false;
        applyProject(Project{});
        projectPath_.clear();
        updateTitle();
    });
}
void MainWindow::openDraftWizard() {
    if (cadMode() || sourceImportBusy_ || loading_)
        return;
    if (meshBusy_ || simulationDock_->isBusy() || simulationDock_->realtimeEnabled()) {
        setDiagnostic(tr("Bitte laufende Netz-/Berechnungsaufträge und Echtzeit zuerst beenden."),
                      true);
        return;
    }
    commitFocusedEdit();
    flushEdit();
    DraftWizard wizard(this);
    if (wizard.exec() != QDialog::Accepted)
        return;
    if (cadMode() || sourceImportBusy_) {
        setDiagnostic(
            tr("Geometriequelle inzwischen geändert · Assistent erneut im Draft öffnen."));
        return;
    }
    try {
        const auto before = editState();
        auto after = before;
        after.parameters = wizard.parameters();
        const bool calculate = wizard.startCalculation();
        const auto previousSettings = simulationDock_->archive().settings;
        auto *mode = simulationDock_->findChild<QComboBox *>("bemGeometryMode");
        const int previousMode = mode->currentIndex();
        if (calculate) {
            const QScopedValueRollback<bool> guard(loading_, true);
            auto *preset = meshPanel_->findChild<QComboBox *>("meshPreset");
            preset->setCurrentIndex(wizard.meshPreset());
            after.meshing = meshPanel_->settings();
            meshPanel_->setSettings(before.meshing);
        }
        const auto minimum = wizard.minimumFrequency(), maximum = wizard.maximumFrequency();
        const auto count = static_cast<std::size_t>(wizard.frequencyCount());
        const bool stateChanged = !sameEditState(before, after);
        const bool solverChanged =
            calculate && (previousMode != 0 || previousSettings.minimumFrequency != minimum ||
                          previousSettings.maximumFrequency != maximum ||
                          previousSettings.frequencyCount != count);
        if (stateChanged)
            applyEditState(after);
        if (calculate) {
            mode->setCurrentIndex(0);
            simulationDock_->setFrequencySweep(minimum, maximum, count);
        }
        if (stateChanged || solverChanged) {
            undoStack_->push(
                new GuidedDraftCommand([this, before, after, previousSettings, previousMode,
                                        minimum, maximum, count, calculate](bool redo) {
                    applyEditState(redo ? after : before);
                    if (calculate) {
                        simulationDock_->findChild<QComboBox *>("bemGeometryMode")
                            ->setCurrentIndex(redo ? 0 : previousMode);
                        simulationDock_->setFrequencySweep(
                            redo ? minimum : previousSettings.minimumFrequency,
                            redo ? maximum : previousSettings.maximumFrequency,
                            redo ? count : previousSettings.frequencyCount);
                    }
                }));
        }
        initialModelFitPending_ = true;
        parameterDock_->show();
        parameterDock_->raise();
        setDiagnostic(tr("Grundhorn übernommen · Parameter rechts bearbeiten oder Netz erzeugen."));
        updateWorkflow();
        if (calculate)
            startGuidedCalculation();
    } catch (const std::invalid_argument &error) {
        setDiagnostic(QString::fromUtf8(error.what()), true);
    }
}

void MainWindow::startGuidedCalculation() {
    if (findChild<QProgressDialog *>("guidedCalculationProgress")) {
        setDiagnostic(tr("Ein geführter Ablauf läuft bereits."), true);
        return;
    }
    struct State {
        enum class Phase { geometry, mesh, solve } phase{Phase::geometry};
        std::uint64_t epoch{}, revision{}, serial{}, generation{}, meshJobId{}, solveJobId{};
        bool ownsMesh{}, ownsSolve{}, finished{};
    };
    auto state = std::make_shared<State>();
    state->epoch = projectEpoch_;
    state->revision = revision_;
    state->serial = changeSerial_;
    auto *progress = new QProgressDialog(tr("Entwurf wird geprüft …"), tr("Abbrechen"), 0, 0, this);
    progress->setObjectName("guidedCalculationProgress");
    progress->setWindowModality(Qt::NonModal);
    progress->setMinimumDuration(0);
    progress->setAutoClose(false);
    auto *timer = new QTimer(progress);
    timer->setObjectName("guidedCalculationTimer");
    timer->setInterval(50);
    const auto finish = [this, progress, timer, state](const QString &message, bool error) {
        if (state->finished)
            return;
        state->finished = true;
        timer->stop();
        progress->close();
        progress->deleteLater();
        setDiagnostic(message, error);
    };
    const auto cancelOwnedJobs = [this, state] {
        if (state->ownsMesh && state->epoch == projectEpoch_ && meshBusy_ &&
            state->meshJobId == meshJobId_)
            meshWorker_.request_stop();
        if (state->ownsSolve && state->epoch == projectEpoch_ && state->solveJobId != 0 &&
            state->solveJobId == simulationDock_->activeSolveJobId())
            simulationDock_->cancelSolve();
    };
    connect(progress, &QProgressDialog::canceled, this, [state, finish, cancelOwnedJobs] {
        if (state->finished)
            return;
        cancelOwnedJobs();
        finish(tr("Geführter Ablauf abgebrochen."), false);
    });
    connect(timer, &QTimer::timeout, this, [this, state, progress, finish, cancelOwnedJobs] {
        if (closing_ || state->epoch != projectEpoch_ || state->revision != revision_) {
            cancelOwnedJobs();
            finish(tr("Geführter Ablauf beendet: Der Entwurf wurde inzwischen geändert."), true);
            return;
        }
        if (state->phase == State::Phase::geometry) {
            if (changeSerial_ != state->serial) {
                finish(tr("Geführter Ablauf beendet: Einstellungen wurden geändert."), true);
                return;
            }
            if (!geometryReady()) {
                if (!geometryBusy_ && !geometryPending_ && !debounce_->isActive())
                    finish(tr("Entwurf ungültig · Geometriediagnose prüfen."), true);
                return;
            }
            progress->setLabelText(tr("Isotropes Qualitätsnetz wird erzeugt …"));
            state->phase = State::Phase::mesh;
            if (!meshReady()) {
                requestMesh();
                state->ownsMesh = meshBusy_;
                if (state->ownsMesh)
                    state->meshJobId = meshJobId_;
                if (!state->ownsMesh)
                    finish(meshError_.isEmpty() ? tr("Netz konnte nicht gestartet werden.")
                                                : meshError_,
                           true);
            }
            return;
        }
        if (state->phase == State::Phase::mesh) {
            if (meshBusy_) {
                if (changeSerial_ != state->serial) {
                    cancelOwnedJobs();
                    finish(tr("Geführter Ablauf beendet: Einstellungen wurden geändert."), true);
                }
                return;
            }
            const auto expectedSerial = state->serial + (state->ownsMesh && meshReady() ? 1 : 0);
            state->ownsMesh = false;
            if (changeSerial_ != expectedSerial) {
                finish(tr("Geführter Ablauf beendet: Einstellungen wurden geändert."), true);
                return;
            }
            if (!meshReady()) {
                finish(meshError_.isEmpty() ? tr("Das Qualitätsnetz ist nicht verfügbar.")
                                            : meshError_,
                       true);
                return;
            }
            if (!simulationDock_->hasConfiguredExcitation()) {
                finish(tr("Keine gültige Driven-Anregung · Randbedingungen prüfen."), true);
                return;
            }
            const double dofs = static_cast<double>(mesh_->mesh.vertices.size());
            const double denseMiB = dofs * dofs * sizeof(std::complex<double>) / (1024.0 * 1024.0);
            progress->setLabelText(
                tr("Netz: %1 Dreiecke, %2 P1-Knoten. Eine dichte komplexe "
                   "Knotenmatrix allein benötigt rechnerisch etwa %3 MiB; der tatsächliche Bedarf "
                   "hängt vom Verfahren ab. Die Berechnung läuft …")
                    .arg(mesh_->mesh.triangles.size())
                    .arg(mesh_->mesh.vertices.size())
                    .arg(locale().toString(denseMiB, 'f', 1)));
            state->phase = State::Phase::solve;
            state->generation = simulationDock_->inputGeneration();
            requestSolve();
            state->ownsSolve = simulationDock_->isBusy();
            if (state->ownsSolve)
                state->solveJobId = simulationDock_->activeSolveJobId();
            if (!state->ownsSolve)
                finish(simulationDock_->lastError().isEmpty()
                           ? tr("Berechnung konnte nicht gestartet werden.")
                           : simulationDock_->lastError(),
                       true);
            return;
        }
        if (state->generation != simulationDock_->inputGeneration()) {
            cancelOwnedJobs();
            state->ownsSolve = false;
            finish(tr("Geführter Ablauf beendet: Berechnungseingaben wurden geändert."), true);
            return;
        }
        if (simulationDock_->isBusy())
            return;
        state->ownsSolve = false;
        if (simulationDock_->hasCurrentResult())
            finish(tr("Netz und Berechnung abgeschlossen · Ergebnisse unten wählen."), false);
        else
            finish(simulationDock_->lastError().isEmpty()
                       ? tr("Berechnung lieferte kein aktuelles Ergebnis.")
                       : simulationDock_->lastError(),
                   true);
    });
    progress->show();
    timer->start();
}

void MainWindow::updateWorkflow() {
    if (!workflowStatus_)
        return;
    draftWizardAction_->setEnabled(!cadMode() && !sourceImportBusy_ && !loading_);
    workflowMesh_->setEnabled(geometryReady() && !meshBusy_ && !simulationInProgress());
    workflowSimulation_->setEnabled(!sourceImportBusy_);
    QString text;
    QString detail;
    if (sourceImportBusy_)
        text = tr("Geometrie wird geladen …");
    else if (geometryBusy_ || geometryPending_)
        text = tr("Entwurf wird geprüft …");
    else if (meshBusy_)
        text = tr("Netz wird geprüft · Abbruch in der Werkzeugleiste");
    else if (simulationInProgress())
        text = tr("Berechnung läuft · Fortschritt und Abbruch unten");
    else if (!geometryReady()) {
        text = cadMode() ? tr("CAD-Verbindung und Flächenzuordnung prüfen")
                         : tr("Entwurf prüfen · Hinweise rechts");
        detail = diagnostic_->toolTip();
    } else if (simulationReady())
        text = tr("Ergebnis aktuell · Ansichten unten wählen");
    else if (!lastSimulationError().isEmpty()) {
        text = tr("Berechnung fehlgeschlagen · Details unten");
        detail = lastSimulationError();
    } else if (meshReady() || simulationDock_->usesAxisymmetricGeometry())
        text = tr("Bereit · Frequenzen und Anregung unter Berechnung prüfen");
    else if (!meshError_.isEmpty()) {
        text = tr("Netz fehlgeschlagen · Einstellungen links prüfen");
        detail = meshError_;
    } else
        text = tr("Entwurf bereit · Netz erzeugen oder Horn-Assistent öffnen");
    workflowStatus_->setText(text);
    workflowStatus_->setToolTip(detail.isEmpty() ? text : detail);
}
void MainWindow::exportDraftStep() {
    commitFocusedEdit();
    flushEdit();
    if (cadMode() || !geometry_ || geometryBusy_ || geometryPending_ ||
        geometryRevision_ != revision_) {
        setDiagnostic(tr("STEP-Export benötigt eine aktuelle gültige Draft-Geometrie."), true);
        return;
    }
    auto *dialog = new StepExportDialog(geometry_, geometryRevision_, this);
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->show();
}

std::vector<std::uint32_t>
MainWindow::placedSelection(const std::vector<std::uint32_t> &ids) const {
    if (cadMode())
        return ids;
    std::vector<std::uint32_t> result;
    for (const auto &instance : assembly_.instances)
        for (const auto id : ids)
            if (id >= 1 && id <= 8)
                result.push_back(instanceFaceId(instance.id, id));
    return result;
}

void MainWindow::editDraftAssembly() {
    if (cadMode()) {
        setDiagnostic(tr("Zum Platzieren von Hörnern zuerst den Draft-Arbeitsbereich öffnen."));
        return;
    }
    commitFocusedEdit();
    DraftAssemblyDialog dialog(assembly_, static_cast<bool>(cadSnapshot_), this);
    if (dialog.exec() != QDialog::Accepted)
        return;
    flushEdit();
    assembly_ = dialog.assembly();
    parametersChanged();
    rememberEdit(true);
    updateSourcePresentation();
}

void MainWindow::addCadToDraft() {
    if (sourceImportBusy_)
        return;
    if (cadMode()) {
        setDiagnostic(tr("Zuerst den Draft öffnen, zu dem das CAD-Modell gehören soll."));
        return;
    }
    commitFocusedEdit();
    const auto path = QFileDialog::getOpenFileName(this, tr("CAD zur Draft-Baugruppe hinzufügen"),
                                                   {}, tr("STEP (*.step *.stp)"));
    if (path.isEmpty())
        return;
    SourceImport request{path,      true,          ++sourceImportId_, projectEpoch_,
                         revision_, changeSerial_, meshing_};
    request.appendToDraft = std::make_shared<Project>(projectState());
    startSourceImport(std::move(request));
}
void MainWindow::setCompactToolsMode(bool compact) {
    if (!compactTools_ || !compactToolsDialog_ || !workflowScroll_ || !viewScroll_)
        return;
    compactHeight_ = compact;
    compactTools_->setVisible(compact);
    if (compact) {
        compactToolsDialog_->layout()->addWidget(workflowScroll_);
        compactToolsDialog_->layout()->addWidget(viewScroll_);
        workflowScroll_->show();
        viewScroll_->show();
        if (compactTools_->isChecked())
            compactToolsDialog_->show();
    } else {
        compactToolsDialog_->hide();
        centerLayout_->insertWidget(1, workflowScroll_);
        centerLayout_->insertWidget(2, viewScroll_);
        workflowScroll_->show();
        viewScroll_->show();
    }
}

void MainWindow::resizeEvent(QResizeEvent *event) {
    QMainWindow::resizeEvent(event);
    const bool compact = event->size().height() < 650;
    if (!compactTools_ || compact == compactHeight_)
        return;
    setCompactToolsMode(compact);
    compactTools_->setChecked(compact &&
                              QSettings().value("view/compactToolsExpanded", false).toBool());
    QTimer::singleShot(0, this, [this] {
        layout()->invalidate();
        layout()->activate();
    });
}

void MainWindow::closeEvent(QCloseEvent *event) {
    commitFocusedEdit();
    if (!closeApproved_ && saveWorker_->busy() && afterSave_) {
        afterSave_ = [this] {
            closeApproved_ = true;
            close();
        };
        event->ignore();
        setDiagnostic(tr("Nach dem Speichern wird geschlossen."));
        return;
    }
    if (!closeApproved_ && (saveWorker_->busy() || dirty_)) {
        event->ignore();
        confirmUnsaved([this] {
            closeApproved_ = true;
            close();
        });
        return;
    }
    closeApproved_ = false;
    clearRecovery(false);
    QSettings settings;
    settings.setValue("window/geometry", saveGeometry());
    settings.setValue("window/state", saveState());
    event->accept();
}

bool MainWindow::cadMode() const {
    return sourceMode_->currentIndex() == 1;
}

void MainWindow::sourceModeChanged() {
    if (closing_ || loading_)
        return;
    if (isPlacedAssembly(assembly_)) {
        const QSignalBlocker blocker(sourceMode_);
        sourceMode_->setCurrentIndex(0);
        setDiagnostic(
            tr("Baugruppe aktiv · CAD über Geometrie → CAD zur Draft-Baugruppe hinzufügen."));
        return;
    }
    flushEdit();
    ++revision_;
    markDirty();
    debounce_->stop();
    geometryPending_ = false;
    geometryWorker_.request_stop();
    meshWorker_.request_stop();
    mesh_.reset();
    meshRevision_ = 0;
    preview_.reset();
    rawPreview_.reset();
    geometryRevision_ = 0;
    simulationDock_->setMesh({}, revision_);
    viewport_->fitView();
    viewport_->setSelectedFaces({});
    meshPanel_->setSelectedFaces({});
    cadDock_->clearSelection();
    representation_->setCurrentIndex(0);
    meshAction_->setEnabled(false);
    qualityLabel_->setText(tr("Netz: nicht erzeugt"));
    meshItem_->setText(0, tr("Netz · nicht erzeugt"));
    updateSourcePresentation();
    if (tabifiedDockWidgets(parameterDock_).contains(cadDock_))
        (cadMode() ? static_cast<QDockWidget *>(cadDock_) : parameterDock_)->raise();
    if (cadMode()) {
        if (!gmsh_)
            cadDock_->ensureListening();
        if (const auto current = cadDock_->snapshot())
            acceptCadSnapshot(current);
        else {
            revisionLabel_->setText(tr("CAD · wartet auf Geometrie"));
            updateDisplay();
        }
        simulationDock_->raise();
    } else {
        analysisDock_->show();
        analysisDock_->raise();
        geometryPending_ = true;
        startGeometry();
    }
}

void MainWindow::updateSourcePresentation() {
    meshAction_->setText(gmsh_ ? tr("Importnetz prüfen") : tr("Netz erzeugen"));
    meshAction_->setToolTip(
        gmsh_ ? tr("Importiertes Netz unverändert auf Qualität und Gebietsabschluss prüfen.")
              : tr("Isotropes Netz erzeugen und auf geschlossenen Rand, Orientierung und "
                   "Elementqualität prüfen."));
    updateMeshFaces();
    meshSettingsDock_->show();
    geometryTabs_->setTabVisible(0, !cadMode());
    while (facesItem_->childCount())
        delete facesItem_->takeChild(0);
    while (sourceItem_->childCount())
        delete sourceItem_->takeChild(0);
    if (cadMode()) {
        parameterDock_->hide();
        cadDock_->show();
        projectItem_->setText(0, cadSnapshot_ ? cadSnapshot_->documentName : tr("FreeCAD-Link"));
        sourceItem_->setText(0, cadSnapshot_ ? tr("CAD-Revision %1").arg(cadSnapshot_->revision)
                                             : tr("CAD-Geometrie"));
        updateInterfacePresentation();
        updatePlots({});
        analysisRevision_ = 0;
        analysisDock_->setWindowTitle(
            tr("Axiale Analyse · für diese CAD-Geometrie nicht verfügbar"));
    } else {
        cadDock_->hide();
        parameterDock_->show();
        projectItem_->setText(0,
                              isPlacedAssembly(assembly_) ? tr("Draft-Baugruppe") : tr("Entwurf"));
        sourceItem_->setText(0,
                             tr("%1 gekoppelte Draft-Instanz(en)").arg(assembly_.instances.size()));
        for (const auto &instance : assembly_.instances)
            new QTreeWidgetItem(sourceItem_,
                                {QString::fromStdString(instance.name) +
                                 (instance.mirrorHorizontal ? tr(" · gespiegelt") : QString{})});
        if (assembly_.includeCad && cadSnapshot_)
            new QTreeWidgetItem(sourceItem_, {tr("CAD · %1 · h = %2 mm")
                                                  .arg(cadSnapshot_->documentName)
                                                  .arg(assembly_.cadTargetEdgeLength * 1000)});
        sourceItem_->setExpanded(true);
        new QTreeWidgetItem(facesItem_, {tr("Hornwand · Wall")});
        new QTreeWidgetItem(facesItem_, {tr("Hals · Driven")});
        analysisDock_->setWindowTitle(
            analysisRevision_ == geometryRevision_ && analysisRevision_ != 0
                ? tr("Geometrieanalyse · Rev. %1").arg(analysisRevision_)
                : tr("Geometrieanalyse · Rev. %1 · veraltet · wird aktualisiert")
                      .arg(analysisRevision_));
    }
}

void MainWindow::acceptCadSnapshot(std::shared_ptr<const CadSnapshot> snapshot) {
    if (closing_ || loading_ || !snapshot)
        return;
    const bool newDocument = !cadSnapshot_ || cadSnapshot_->documentId != snapshot->documentId;
    const bool sameGeometry = cadSnapshot_ && cadSnapshot_->documentId == snapshot->documentId &&
                              cadSnapshot_->sessionId == snapshot->sessionId &&
                              cadSnapshot_->revision == snapshot->revision;
    const auto before = cadFaceTargetEdgeLengths_.size();
    const auto reversedBefore = cadReversedFaceIds_.size();
    const auto colorsBefore = cadFaceColors_.size();
    auto network =
        reconciledCadBemAssignments(lemEditor_->network(), cadSnapshot_.get(), snapshot.get());
    const bool removedPorts = !lem::samePhysics(network, lemEditor_->network());
    auto regions = reconciledCadAcousticModel(simulationDock_->archive().settings.acousticModel,
                                              cadSnapshot_.get(), snapshot.get());
    const bool changedRegions = regions != simulationDock_->archive().settings.acousticModel;
    cadFaceColors_ = retainedCadFaceColors(cadFaceColors_, cadSnapshot_.get(), snapshot.get());
    if (cadSnapshot_ && cadSnapshot_->documentId == snapshot->documentId) {
        // FreeCadLink carries an ID forward only after unique mapped-name,
        // geometry or topology reconciliation; ordinal FaceN names are not proof.
        const auto removedFace = [&](std::uint32_t id) {
            const auto find = [&](const auto &faces) {
                return std::find_if(faces.begin(), faces.end(),
                                    [&](const auto &face) { return face.id == id; });
            };
            const auto previous = find(cadSnapshot_->faces), current = find(snapshot->faces);
            return previous == cadSnapshot_->faces.end() || current == snapshot->faces.end() ||
                   previous->objectId != current->objectId ||
                   current->state == CadFaceState::needsReview;
        };
        std::erase_if(cadFaceTargetEdgeLengths_,
                      [&](const auto &entry) { return removedFace(entry.first); });
        std::erase_if(cadReversedFaceIds_, removedFace);
    } else {
        cadFaceTargetEdgeLengths_.clear();
        cadReversedFaceIds_.clear();
    }
    if (before != cadFaceTargetEdgeLengths_.size())
        meshPanel_->setFaceNotice(tr("CAD · %1 lokale Zielkanten nach Flächenänderung entfernt")
                                      .arg(before - cadFaceTargetEdgeLengths_.size()));
    if (reversedBefore != cadReversedFaceIds_.size())
        meshPanel_->setFaceNotice(tr("CAD · %1 nicht mehr zuordenbare Normalenumkehrungen entfernt")
                                      .arg(reversedBefore - cadReversedFaceIds_.size()));
    if (colorsBefore != cadFaceColors_.size())
        setDiagnostic(tr("CAD · %1 nicht mehr zuordenbare Flächenfarben entfernt")
                          .arg(colorsBefore - cadFaceColors_.size()));
    cadSnapshot_ = snapshot;
    if (changedRegions) {
        simulationDock_->setAcousticModel(std::move(regions));
        setDiagnostic(tr("Gebietszuordnungen nach CAD-Änderung prüfen."), true);
    }
    if (removedPorts) {
        lemEditor_->setNetwork(std::move(network));
        simulationDock_->setCoupling(lemEditor_->network(), lemEditor_->couplingEnabled());
        setDiagnostic(tr("Nicht mehr zuordenbare BEM↔LEM-Flächen entfernt · Anschlüsse prüfen."),
                      true);
    }
    cadDock_->setReversedFaces(cadReversedFaceIds_);
    markDirty();
    if (!cadMode()) {
        if (assembly_.includeCad) {
            parametersChanged();
            updateSourcePresentation();
        }
        if (!meshSettingsValid_ && before != cadFaceTargetEdgeLengths_.size())
            meshingChanged();
        if (sameGeometry && !changedRegions)
            rememberEdit(true);
        else
            resetEditHistory();
        return;
    }
    QString sizingError;
    try {
        const auto settings = meshPanel_->settings();
        validateFaceTargets(settings);
        meshing_ = settings;
        meshSettingsValid_ = true;
    } catch (const std::exception &error) {
        meshSettingsValid_ = false;
        sizingError = QString::fromUtf8(error.what());
    }
    const bool reusableMesh =
        sameGeometry && !changedRegions && mesh_ && meshRevision_ == revision_ &&
        !(gmsh_ && !simulationDock_->archive().settings.acousticModel.domains.empty());
    ++revision_;
    geometryPending_ = false;
    meshWorker_.request_stop();
    simulationDock_->invalidate(revision_);
    cadSnapshot_ = std::move(snapshot);
    if (cadSnapshot_ && cadSnapshot_->sessionId != "gmsh-import-offline")
        gmsh_.reset();
    rawPreview_ = std::shared_ptr<const SurfaceMesh>(cadSnapshot_, &cadSnapshot_->mesh);
    updateOrientedPreview();
    geometryRevision_ = revision_;
    updateSourcePresentation();
    revisionLabel_->setText(
        tr("CAD-Revision %1 · Geometrie Rev. %2").arg(cadSnapshot_->revision).arg(revision_));
    configureLemFaces();
    if (reusableMesh) {
        auto updated = std::make_shared<MeshingResult>(*mesh_);
        QHash<std::uint32_t, BoundaryCondition> boundaries;
        for (const auto &face : cadSnapshot_->faces)
            boundaries.insert(face.id, face.boundary);
        for (auto &triangle : updated->mesh.triangles)
            triangle.boundary = boundaries.value(triangle.faceId, BoundaryCondition::wall);
        mesh_ = std::move(updated);
        meshRevision_ = revision_;
        updateMeshStatus();
    } else {
        if (mesh_)
            qualityLabel_->setText(tr("Netz Rev. %1 · veraltet").arg(meshRevision_));
        representation_->setCurrentIndex(0);
    }
    meshAction_->setEnabled(geometryReady() && !meshBusy_);
    if (geometryReady())
        configureSimulation();
    updateDisplay();
    if (newDocument)
        viewport_->fitView();
    const auto currentRegions = simulationDock_->archive().settings.acousticModel;
    if (std::any_of(currentRegions.domains.begin(), currentRegions.domains.end(),
                    [](const auto &region) { return region.needsAssignmentReview; }) ||
        std::any_of(currentRegions.interfaces.begin(), currentRegions.interfaces.end(),
                    [](const auto &interface) { return interface.needsAssignmentReview; }))
        setDiagnostic(tr("Gebietszuordnungen nach CAD-Änderung prüfen."), true);
    else if (cadSnapshot_->hasUnresolvedAssignments())
        setDiagnostic(tr("Flächenzuordnungen prüfen."), true);
    else if (!sizingError.isEmpty())
        setDiagnostic(sizingError, true);
    else
        setDiagnostic({});
    if (sameGeometry)
        rememberEdit(true);
    else
        resetEditHistory();
}

void MainWindow::cadAvailabilityChanged(bool available) {
    if (closing_ || loading_ || !cadMode())
        return;
    if (!available) {
        ++revision_;
        markDirty();
        meshWorker_.request_stop();
        simulationDock_->invalidate(revision_);
        meshAction_->setEnabled(false);
        revisionLabel_->setText(tr("CAD · Änderung ausstehend oder Verbindung getrennt"));
        if (mesh_)
            qualityLabel_->setText(tr("Netz Rev. %1 · veraltet").arg(meshRevision_));
        updateDisplay();
    } else
        meshAction_->setEnabled(geometryReady() && !meshBusy_);
}

void MainWindow::meshingChanged() {
    if (loading_)
        return;
    rememberEdit();
    QString error;
    bool geometryAccuracyChanged = false;
    try {
        const auto next = meshPanel_->settings();
        validateFaceTargets(next);
        geometryAccuracyChanged =
            next.adaptive.geometricTolerance != meshing_.adaptive.geometricTolerance;
        meshing_ = next;
        meshSettingsValid_ = true;
    } catch (const std::exception &e) {
        meshSettingsValid_ = false;
        error = QString::fromUtf8(e.what());
    }
    if (!cadMode() && (geometryAccuracyChanged || geometryBusy_ || geometryPending_ ||
                       geometryRevision_ != revision_ || !geometry_)) {
        parametersChanged();
        if (!error.isEmpty())
            setDiagnostic(error, true);
        return;
    }
    markDirty();
    const bool analysisCurrent = analysisRevision_ != 0 && analysisRevision_ == geometryRevision_;
    ++revision_;
    if (preview_)
        geometryRevision_ = revision_;
    if (analysisCurrent)
        analysisRevision_ = geometryRevision_;
    meshWorker_.request_stop();
    simulationDock_->invalidate(revision_);
    if (mesh_)
        qualityLabel_->setText(tr("Netz Rev. %1 · veraltet").arg(meshRevision_));
    meshAction_->setEnabled(geometryReady() && !meshBusy_);
    revisionLabel_->setText(tr("Geometrie Rev. %1 · Netzparameter geändert").arg(revision_));
    if (!cadMode())
        analysisDock_->setWindowTitle(
            analysisCurrent ? tr("Geometrieanalyse · Rev. %1").arg(analysisRevision_)
                            : tr("Geometrieanalyse · Rev. %1 · veraltet").arg(analysisRevision_));
    updateDisplay();
    if (geometryReady())
        configureSimulation();
    setDiagnostic(error, !error.isEmpty());
}

std::map<std::uint32_t, double> &MainWindow::activeFaceTargets() {
    return cadMode() ? cadFaceTargetEdgeLengths_ : draftFaceTargetEdgeLengths_;
}
std::set<std::uint32_t> &MainWindow::activeReversedFaces() {
    return cadMode() ? cadReversedFaceIds_ : draftReversedFaceIds_;
}
void MainWindow::updateOrientedPreview() {
    if (!rawPreview_) {
        preview_.reset();
        return;
    }
    if (!cadMode()) {
        // The default instance is the identity and retains its source face IDs.
        if (!isPlacedAssembly(assembly_) && draftReversedFaceIds_.empty()) {
            preview_ = rawPreview_;
            return;
        }
        Project p;
        p.parameters = geometry_ ? geometry_->parameters() : parameters_->parameters();
        p.assembly = assembly_;
        p.cad = cadSnapshot_;
        p.draftReversedFaceIds = draftReversedFaceIds_;
        p.cadReversedFaceIds = cadReversedFaceIds_;
        preview_ = std::make_shared<SurfaceMesh>(assemblyPreview(*rawPreview_, p));
        return;
    }
    const auto &reversed = activeReversedFaces();
    preview_ = reversed.empty()
                   ? rawPreview_
                   : std::make_shared<SurfaceMesh>(applyFaceOrientation(*rawPreview_, reversed));
}
void MainWindow::reverseFaces(const std::vector<std::uint32_t> &ids) {
    if (loading_ || ids.empty() || !rawPreview_)
        return;
    if (cadMode() && (!cadDock_->isSourceAvailable() || cadDock_->isGeometryPending())) {
        setDiagnostic(tr("Normalen bearbeiten: auf den vollständigen CAD-Stand warten."), true);
        return;
    }
    std::set<std::uint32_t> available;
    for (const auto &triangle : rawPreview_->triangles)
        available.insert(triangle.faceId);
    if (std::any_of(ids.begin(), ids.end(), [&](auto id) { return !available.contains(id); }))
        return;
    flushEdit();
    auto &reversed = activeReversedFaces();
    for (const auto id : ids)
        if (reversed.contains(id))
            reversed.erase(id);
        else
            reversed.insert(id);
    updateOrientedPreview();
    representation_->setCurrentIndex(0);
    meshColor_->setCurrentIndex(static_cast<int>(MeshColorMode::orientation));
    meshingChanged();
    updateMeshFaces();
    setDiagnostic(
        tr("Normalen geändert · BEM erfordert ein geprüftes, nach außen orientiertes Netz."));
    rememberEdit(true);
}
void MainWindow::setFaceTargets(const std::vector<std::uint32_t> &ids,
                                std::optional<double> target) {
    std::set<std::uint32_t> available;
    if (rawPreview_)
        for (const auto &triangle : rawPreview_->triangles)
            available.insert(triangle.faceId);
    if (ids.empty() ||
        std::any_of(ids.begin(), ids.end(), [&](auto id) { return !available.contains(id); }))
        return;
    if (target) {
        const auto settings = meshPanel_->inputSettings();
        if (!std::isfinite(*target) || *target < settings.adaptive.minimumEdgeLength ||
            *target > settings.targetEdgeLength) {
            setDiagnostic(tr("Lokale Zielkante muss zwischen Mindestkante und globaler "
                             "Zielkante im Netzpanel liegen."),
                          true);
            return;
        }
    }
    auto next = activeFaceTargets();
    for (const auto id : ids)
        if (target)
            next[id] = *target;
        else
            next.erase(id);
    if (next == activeFaceTargets())
        return;
    flushEdit();
    activeFaceTargets() = std::move(next);
    meshPanel_->setFaceNotice({});
    updateMeshFaces();
    meshingChanged();
    rememberEdit(true);
}
void MainWindow::validateFaceTargets(const MeshSettings &settings) const {
    const auto validate = [&](const auto &targets, const QString &source) {
        for (const auto &[id, target] : targets)
            if (!std::isfinite(target) || target < settings.adaptive.minimumEdgeLength ||
                target > settings.targetEdgeLength)
                throw std::invalid_argument(
                    tr("%1 · Fläche %2: lokale Zielkante außerhalb Min./Zielkante. Zuweisung "
                       "korrigieren oder zurücksetzen.")
                        .arg(source)
                        .arg(id)
                        .toStdString());
    };
    validate(draftFaceTargetEdgeLengths_, tr("Draft"));
    validate(cadFaceTargetEdgeLengths_, tr("CAD"));
}
void MainWindow::updateMeshFaces() {
    std::vector<MeshFaceItem> faces;
    if (cadMode()) {
        if (cadSnapshot_)
            for (const auto &face : cadSnapshot_->faces)
                faces.push_back({face.id, face.label});
    } else if (rawPreview_) {
        std::set<std::uint32_t> ids;
        for (const auto &triangle : rawPreview_->triangles)
            ids.insert(triangle.faceId);
        for (const auto id : ids) {
            QString label;
            switch (static_cast<DraftFace>(id)) {
            case DraftFace::innerWall:
                label = tr("Innenwand");
                break;
            case DraftFace::outerWall:
                label = tr("Außenwand");
                break;
            case DraftFace::mouthRim:
                label = tr("Mundrand");
                break;
            case DraftFace::rear:
                label = tr("Rückseite");
                break;
            case DraftFace::driven:
                label = tr("Eintritt · Driven");
                break;
            case DraftFace::adapterWall:
                label = tr("Formadapter");
                break;
            case DraftFace::slotWall:
                label = tr("Spaltwand");
                break;
            case DraftFace::upstreamOuterWall:
                label = tr("Halsaußenwand");
                break;
            default:
                label = tr("Fläche %1").arg(id);
                break;
            }
            faces.push_back({id, std::move(label)});
        }
    }
    meshPanel_->setFaces(std::move(faces), activeFaceTargets());
    meshPanel_->setReversedFaces(activeReversedFaces());
    meshPanel_->setOrientationStatus(activeReversedFaces().empty() ? QString{}
                                     : meshReady()
                                         ? tr("Normalen geprüft")
                                         : tr("%1 Flächen umgekehrt · Netzprüfung erforderlich")
                                               .arg(activeReversedFaces().size()));
    cadDock_->setFaceTargets(cadFaceTargetEdgeLengths_);
    cadDock_->setReversedFaces(cadReversedFaceIds_);
    if (cadMode())
        meshPanel_->setSelectedFaces(cadDock_->selectedFaces());
}

void MainWindow::configureLemFaces() {
    std::vector<LemFaceItem> faces;
    std::map<std::uint32_t, BoundaryCondition> groups;
    if (preview_)
        for (const auto &triangle : preview_->triangles)
            groups.emplace(triangle.faceId, triangle.boundary);
    const auto simulationState = simulationDock_->archive();
    const auto baffle = simulationState.settings.infiniteBaffle;
    std::set<std::uint32_t> bemInterfaces;
    for (const auto &interface : simulationState.settings.acousticModel.interfaces)
        bemInterfaces.insert(interface.faceIds.begin(), interface.faceIds.end());
    for (const auto &[id, boundary] : groups) {
        if (bemInterfaces.contains(id))
            continue;
        const bool aperture =
            baffle && std::find(baffle->apertureFaceIds.begin(), baffle->apertureFaceIds.end(),
                                id) != baffle->apertureFaceIds.end();
        QString label = tr("%1 · Fläche %2")
                            .arg(boundary == BoundaryCondition::driven ? "Driven"
                                 : aperture                            ? "Apertur"
                                                                       : "Wall")
                            .arg(id);
        if (cadMode() && cadSnapshot_) {
            const auto found = std::find_if(cadSnapshot_->faces.begin(), cadSnapshot_->faces.end(),
                                            [id](const auto &face) { return face.id == id; });
            if (found != cadSnapshot_->faces.end())
                label += " · " + found->label;
        }
        faces.push_back({id, std::move(label), boundary, !aperture});
    }
    lemEditor_->setFaces(std::move(faces));
}

void MainWindow::updateInterfacePresentation() {
    std::set<std::uint32_t> interfaces;
    bool needsReview = false;
    for (const auto *element : lem::bemPorts(lemEditor_->network())) {
        const auto &port = std::get<lem::BemPort>(element->component);
        if (port.kind == lem::BemPortKind::interfacePort)
            interfaces.insert(port.faceIds.begin(), port.faceIds.end());
        needsReview |= port.needsAssignmentReview;
    }
    cadDock_->setInterfaceFaces(interfaces, lemEditor_->couplingEnabled());
    std::set<std::uint32_t> beInterfaces;
    for (const auto &interface : simulationDock_->archive().settings.acousticModel.interfaces)
        beInterfaces.insert(interface.faceIds.begin(), interface.faceIds.end());
    cadDock_->setBemInterfaceFaces(beInterfaces);
    if (cadMode()) {
        while (facesItem_->childCount())
            delete facesItem_->takeChild(0);
        std::size_t driven = 0, virtualFaces = 0, beFaces = 0;
        if (cadSnapshot_)
            for (const auto &face : cadSnapshot_->faces) {
                beFaces += beInterfaces.contains(face.id);
                driven +=
                    face.boundary == BoundaryCondition::driven && !beInterfaces.contains(face.id);
                virtualFaces += interfaces.contains(face.id) && !beInterfaces.contains(face.id) &&
                                face.boundary == BoundaryCondition::wall;
            }
        const auto count = cadSnapshot_ ? cadSnapshot_->faces.size() : 0;
        new QTreeWidgetItem(facesItem_,
                            {tr("%1 Wall").arg(count - driven - virtualFaces - beFaces)});
        new QTreeWidgetItem(facesItem_, {tr("%1 Driven").arg(driven)});
        if (virtualFaces)
            new QTreeWidgetItem(
                facesItem_, {tr("%1 Interface ↔ LEM%2")
                                 .arg(virtualFaces)
                                 .arg(lemEditor_->couplingEnabled() ? QString{} : tr(" · aus"))});
        if (beFaces)
            new QTreeWidgetItem(facesItem_, {tr("%1 Interface ↔ BEM").arg(beFaces)});
        if (needsReview)
            new QTreeWidgetItem(facesItem_, {tr("BEM↔LEM-Zuordnung prüfen")});
    }
}

void MainWindow::updateAcousticDomainDisplay(const SimulationArchive &archive) {
    const auto &regions = archive.settings.acousticModel;
    const bool computed = archive.current && archive.revision == revision_ && archive.mesh &&
                          archive.request && archive.request->subdomains;
    std::set<std::uint32_t> faces;
    bool review = false;
    for (const auto &domain : regions.domains)
        review |= domain.needsAssignmentReview;
    for (const auto &interface : regions.interfaces) {
        faces.insert(interface.faceIds.begin(), interface.faceIds.end());
        review |= interface.needsAssignmentReview;
    }
    const auto source =
        computed
            ? archive.mesh
            : (geometryRevision_ == revision_ && !review && !faces.empty()
                   ? (representation_->currentIndex() == 1 && mesh_ && meshRevision_ == revision_
                          ? std::shared_ptr<const SurfaceMesh>(mesh_, &mesh_->mesh)
                          : preview_)
                   : nullptr);
    const bool changed = source != interfaceOverlaySource_ ||
                         computed != interfaceOverlayFromResult_ ||
                         faces != interfaceOverlayFaceIds_ || displayedAcousticRegions_ != regions;
    if (changed) {
        interfaceOverlaySource_ = source;
        interfaceOverlayFromResult_ = computed;
        interfaceOverlayFaceIds_ = faces;
        displayedAcousticRegions_ = regions;
        interfaceOverlay_.reset();
        if (source) {
            auto overlay = std::make_shared<SurfaceMesh>();
            overlay->vertices = source->vertices;
            if (computed) {
                for (const auto &interface : archive.request->subdomains->interfaces)
                    for (const auto triangle : interface.triangleIndices)
                        overlay->triangles.push_back(source->triangles.at(triangle));
            } else {
                for (const auto &triangle : source->triangles)
                    if (faces.contains(triangle.faceId))
                        overlay->triangles.push_back(triangle);
            }
            if (!overlay->triangles.empty())
                interfaceOverlay_ = std::move(overlay);
        }
        while (domainsItem_->childCount())
            delete domainsItem_->takeChild(0);
        domainsItem_->setHidden(regions.domains.empty());
        {
            for (const auto &domain : regions.domains) {
                auto *item = new QTreeWidgetItem(
                    domainsItem_,
                    {QString::fromStdString(domain.name) +
                     (domain.kind == BemDomainKind::interior ? tr(" · Innenraum")
                                                             : tr(" · Außenraum")) +
                     (domain.needsAssignmentReview ? tr(" · Zuordnung prüfen") : QString{})});
                QStringList physical, twoSided;
                for (auto id : domain.faceIds)
                    physical << QString::number(id);
                for (auto id : domain.twoSidedFaceIds)
                    twoSided << QString::number(id);
                if (domain.kind == BemDomainKind::exterior && regions.exteriorBaffle) {
                    const auto &b = *regions.exteriorBaffle;
                    new QTreeWidgetItem(item, {tr("Unendliche Schallwand · Ursprung %1 / %2 / %3 "
                                                  "mm · Normale %4 / %5 / %6 · Toleranz %7 mm")
                                                   .arg(b.origin.x * 1000)
                                                   .arg(b.origin.y * 1000)
                                                   .arg(b.origin.z * 1000)
                                                   .arg(b.normal.x)
                                                   .arg(b.normal.y)
                                                   .arg(b.normal.z)
                                                   .arg(b.planeDistanceTolerance * 1000)});
                }
                if (!physical.empty())
                    new QTreeWidgetItem(item, {tr("Wall / Driven: %1").arg(physical.join(", "))});
                if (!twoSided.empty())
                    new QTreeWidgetItem(item, {tr("Beidseitig: %1").arg(twoSided.join(", "))});
                for (const auto &group : domain.groups) {
                    QStringList members;
                    for (auto id : group.faceIds)
                        members << QString::number(id);
                    new QTreeWidgetItem(item, {tr("%1: %2").arg(QString::fromStdString(group.name),
                                                                members.join(", "))});
                }
            }
            for (const auto &interface : regions.interfaces) {
                const auto domainName = [&](std::uint32_t id) {
                    const auto found = std::find_if(regions.domains.begin(), regions.domains.end(),
                                                    [id](const auto &d) { return d.id == id; });
                    return found == regions.domains.end() ? QString::number(id)
                                                          : QString::fromStdString(found->name);
                };
                auto *item = new QTreeWidgetItem(
                    domainsItem_,
                    {tr("%1 · Subdomain 1: %2 · Subdomain 2: %3%4")
                         .arg(QString::fromStdString(interface.name))
                         .arg(domainName(interface.firstDomainId),
                              domainName(interface.secondDomainId))
                         .arg(interface.needsAssignmentReview ? tr(" · Zuordnung prüfen")
                                                              : QString{})});
                item->setToolTip(0, tr("Flächennormale zeigt in Subdomain 1; die zweite Seite "
                                       "verwendet die Gegenrichtung."));
                QStringList faceLabels;
                for (auto id : interface.faceIds)
                    faceLabels << QString::number(id);
                new QTreeWidgetItem(item, {tr("Flächen: %1").arg(faceLabels.join(", "))});
                for (const auto &group : interface.groups) {
                    QStringList members;
                    for (auto id : group.faceIds)
                        members << QString::number(id);
                    new QTreeWidgetItem(item, {tr("%1: %2").arg(QString::fromStdString(group.name),
                                                                members.join(", "))});
                }
            }
        }
        domainsItem_->setExpanded(true);
    }
    const auto domainName = [&](std::uint32_t id) {
        const auto found = std::find_if(regions.domains.begin(), regions.domains.end(),
                                        [id](const auto &d) { return d.id == id; });
        return found == regions.domains.end() ? QString::number(id)
                                              : QString::fromStdString(found->name);
    };
    std::vector<InterfaceOrientationLabels> labels;
    if (interfaceOverlay_) {
        if (computed) {
            std::size_t offset = 0;
            for (const auto &interface : archive.request->subdomains->interfaces) {
                InterfaceOrientationLabels group{{},
                                                 QString::fromStdString(interface.name),
                                                 domainName(interface.firstDomainId),
                                                 domainName(interface.secondDomainId)};
                for (std::size_t i = 0; i < interface.triangleIndices.size(); ++i)
                    group.triangleIndices.push_back(offset++);
                labels.push_back(std::move(group));
            }
        } else {
            for (const auto &interface : regions.interfaces) {
                InterfaceOrientationLabels group{{},
                                                 QString::fromStdString(interface.name),
                                                 domainName(interface.firstDomainId),
                                                 domainName(interface.secondDomainId)};
                for (std::size_t i = 0; i < interfaceOverlay_->triangles.size(); ++i)
                    if (std::find(interface.faceIds.begin(), interface.faceIds.end(),
                                  interfaceOverlay_->triangles[i].faceId) !=
                        interface.faceIds.end())
                        group.triangleIndices.push_back(i);
                labels.push_back(std::move(group));
            }
        }
    }
    viewport_->setInterfaceMesh(interfaceOverlay_, std::move(labels));
    symmetryItem_->setHidden(!archive.settings.symmetryXYEnabled);
    symmetryItem_->setToolTip(0, tr("Symmetrieebenen x = 0 und y = 0; Toleranz %1 mm")
                                     .arg(archive.settings.symmetryPlaneTolerance * 1000));
    std::optional<OrthogonalSymmetryPlanes> planes;
    if (archive.settings.symmetryXYEnabled) {
        const auto sector = mesh_ && meshRevision_ == revision_
                                ? std::shared_ptr<const SurfaceMesh>(mesh_, &mesh_->mesh)
                                : (geometryRevision_ == revision_ ? preview_ : nullptr);
        if (sector) {
            try {
                planes = inferXYSymmetry(*sector, archive.settings.symmetryPlaneTolerance);
            } catch (const std::invalid_argument &) {
                // Invalid sectors cannot supply a retained-side orientation.
            }
        }
    }
    viewport_->setSymmetryPlanes(planes);
    showFullSymmetry_->setEnabled(planes.has_value());
    auto baffle = regions.exteriorBaffle;
    if (!baffle && archive.settings.infiniteBaffle) {
        const auto &b = *archive.settings.infiniteBaffle;
        baffle = ExteriorBaffle{b.origin, b.normal, 0};
    }
    viewport_->setExteriorBaffle(baffle);
    viewport_->setInterfaceVisible(showInterfaces_->isChecked());
}
void MainWindow::configureSimulation() {
    // Changes of the derived source frame are not user edits. In particular,
    // the first generated model must not dirty an untouched new project.
    const QScopedValueRollback<bool> configuring(loading_, true);
    configureLemFaces();
    const bool currentGeometry = geometryReady();
    std::vector<FaceVelocity> velocities;
    if (currentGeometry && cadMode()) {
        for (const auto &face : cadSnapshot_->faces)
            if (face.boundary == BoundaryCondition::driven)
                velocities.push_back({face.id, face.normalVelocity});
    } else if (currentGeometry && preview_) {
        std::set<std::uint32_t> driven;
        for (const auto &triangle : preview_->triangles)
            if (triangle.boundary == BoundaryCondition::driven)
                driven.insert(triangle.faceId);
        for (const auto face : driven) {
            std::complex<double> velocity{1, 0};
            if (face >= assemblyCadFaceOffset && cadSnapshot_) {
                const auto local = face - assemblyCadFaceOffset;
                const auto found =
                    std::find_if(cadSnapshot_->faces.begin(), cadSnapshot_->faces.end(),
                                 [local](const auto &f) { return f.id == local; });
                if (found != cadSnapshot_->faces.end())
                    velocity = found->normalVelocity;
            }
            velocities.push_back({face, velocity});
        }
    }
    {
        const auto settings = meshPanel_->inputSettings();
        MeshingOptions options;
        options.targetEdgeLength = settings.targetEdgeLength;
        options.adaptive = settings.adaptive;
        options.quality = settings.quality;
        options.faceTargetEdgeLengths = activeFaceTargets();
        simulationDock_->setPartitionMeshingOptions(std::move(options));
    }
    simulationDock_->setMesh(currentGeometry && meshReady()
                                 ? std::shared_ptr<const SurfaceMesh>(mesh_, &mesh_->mesh)
                                 : nullptr,
                             revision_, std::move(velocities));
    simulationDock_->setGeometrySource(
        currentGeometry && !cadMode() && !isPlacedAssembly(assembly_) ? geometry_ : nullptr,
        currentGeometry && cadMode() ? cadSnapshot_ : nullptr, revision_, activeFaceTargets(),
        activeReversedFaces());
    if (!currentGeometry)
        return;
    // The source frame remains separate from a user-defined observation frame.
    if (cadMode())
        simulationDock_->setAcousticAxis(cadSnapshot_->acousticOrigin, cadSnapshot_->acousticAxis,
                                         true);
    else {
        const auto frame = assemblyAcousticFrame(assembly_, geometry_->length());
        simulationDock_->setAcousticAxis(frame.origin, frame.forward, true);
    }
}

void MainWindow::couplingChanged() {
    if (loading_)
        return;
    simulationDock_->setCoupling(lemEditor_->network(), lemEditor_->couplingEnabled());
    updateInterfacePresentation();
    markDirty();
    rememberEdit();
}

void MainWindow::requestSolve() {
    if (!geometryReady() || (!simulationDock_->usesAxisymmetricGeometry() && !meshReady())) {
        setDiagnostic(simulationDock_->usesAxisymmetricGeometry()
                          ? tr("Eine aktuelle exakte Rotationsgeometrie wird benötigt.")
                          : tr("Ein aktuelles geprüftes Netz wird benötigt."),
                      true);
        return;
    }
    simulationDock_->show();
    simulationDock_->raise();
    simulationDock_->requestSolve();
}
bool MainWindow::simulationReady() const {
    return simulationDock_->hasCurrentResult();
}
bool MainWindow::simulationInProgress() const {
    return simulationDock_->isBusy();
}
QString MainWindow::lastSimulationError() const {
    return simulationDock_->lastError();
}
void MainWindow::setSimulationFrequencySweep(double minimum, double maximum, int count) {
    simulationDock_->setFrequencySweep(minimum, maximum, count);
}
void MainWindow::setSimulationView(const QString &view) {
    if (view == "lem") {
        auto *panel = findChild<QDockWidget *>("lemDock");
        panel->show();
        panel->raise();
        resizeDocks({panel}, {height() * 3 / 5}, Qt::Vertical);
        QTimer::singleShot(0, lemEditor_,
                           [this] { lemEditor_->findChild<QPushButton *>("lemFit")->click(); });
        return;
    }
    const auto target = view == "frequency" ? "acousticOnAxis"
                        : view == "polar"   ? "acousticPolarH"
                        : view == "isobar"  ? "acousticIsobarH"
                        : view == "angular" ? "angularResponses"
                        : view == "phase"   ? "phaseResponses"
                        : view == "field"   ? "acousticField"
                                            : nullptr;
    if (!target)
        throw std::invalid_argument("Unbekannte Ergebnisansicht.");
    if (view == "field") {
        auto *enlargeResults = findChild<QAction *>("enlargeResults");
        if (enlargeResults && enlargeResults->isChecked())
            enlargeResults->setChecked(false);
        ++dockResizeGeneration_;
    }
    auto *tabs = simulationDock_->findChild<QTabWidget *>("acousticResults");
    int index = -1;
    for (int i = 0; i < tabs->count(); ++i)
        if (tabs->widget(i)->findChild<QWidget *>(QString::fromLatin1(target))) {
            index = i;
            break;
        }
    if (index < 0)
        throw std::logic_error("Ergebnisansicht fehlt.");
    tabs->setCurrentIndex(index);
    simulationDock_->show();
    simulationDock_->raise();
    simulationDock_->revealResultView();
    // Explicitly opening an acoustic view should give the result room; the
    // native dock layout still respects the viewport and side-dock limits.
    resizeDocks({simulationDock_},
                {view == "field" ? std::clamp(height() * 32 / 100, 200, 280) : height() * 3 / 5},
                Qt::Vertical);
    if (view == "field")
        fieldVisibility_->setChecked(true);
}
void MainWindow::requestFieldEvaluation() {
    simulationDock_->requestFieldEvaluation();
}
bool MainWindow::fieldEvaluationInProgress() const noexcept {
    return simulationDock_->fieldEvaluationInProgress();
}
bool MainWindow::fieldReady() const noexcept {
    return simulationDock_->fieldReady();
}
QString MainWindow::lastFieldError() const {
    return simulationDock_->lastFieldError();
}

void MainWindow::exerciseFaceColorEditing() {
    if (!preview_ || geometryRevision_ != revision_)
        return;
    std::set<std::uint32_t> known;
    for (const auto &triangle : preview_->triangles)
        known.insert(triangle.faceId);
    if (known.empty())
        throw std::runtime_error("Keine Farbflächen im aktuellen Modell.");
    const auto colors = activeFaceColors();
    const auto selection = viewport_->selectedFaces();
    const auto colorMode = meshColor_->currentIndex();
    const auto revision = revision_, geometryRevision = geometryRevision_,
               meshRevision = meshRevision_;
    const auto mesh = mesh_;
    const auto result = simulationDock_->result();
    const bool current = simulationReady(), busy = simulationInProgress(), dirty = dirty_;
    const std::vector<std::uint32_t> selected(
        known.begin(), std::next(known.begin(), std::min<std::size_t>(2, known.size())));
    const auto select = [&](const std::vector<std::uint32_t> &ids) {
        if (cadMode())
            cadDock_->setSelectedFaces(ids);
        else
            meshPanel_->setSelectedFaces(ids);
        viewport_->setSelectedFaces(ids);
        updateFaceColorActions();
    };
    select({});
    if (assignFaceColor_->isEnabled() || resetFaceColors_->isEnabled())
        throw std::runtime_error("Farbwerkzeug erlaubt eine Zuweisung ohne Flächenauswahl.");
    select(selected);
    if (!assignFaceColor_->isEnabled())
        throw std::runtime_error("Farbwerkzeug folgt der Flächenauswahl nicht.");
    const FaceColor chosen{23, 145, 231};
    setSelectedFaceColor(chosen);
    for (const auto id : selected)
        if (!activeFaceColors().contains(id) || activeFaceColors().at(id) != chosen)
            throw std::runtime_error("Farbwerkzeug übernimmt die Mehrfachauswahl nicht.");
    if (revision_ != revision || geometryRevision_ != geometryRevision ||
        meshRevision_ != meshRevision || mesh_ != mesh || simulationDock_->result() != result ||
        simulationReady() != current || simulationInProgress() != busy || !dirty_ ||
        meshColor_->currentIndex() != static_cast<int>(MeshColorMode::faceColors))
        throw std::runtime_error(
            "Eine reine Flächenfarbe verändert Geometrie, Netz oder Solverzustand.");
    QTemporaryDir directory;
    if (!directory.isValid())
        throw std::runtime_error("Temporäres Verzeichnis fehlt.");
    const auto path = directory.filePath("face-colors.hsproj");
    writeProject(path, projectState());
    const auto restored = readProject(path);
    if (restored.draftFaceColors != draftFaceColors_ || restored.cadFaceColors != cadFaceColors_)
        throw std::runtime_error("Farbwerkzeug verliert Farbzuweisungen beim Speichern.");
    resetFaceColors_->trigger();
    for (const auto id : selected)
        if (activeFaceColors().contains(id))
            throw std::runtime_error("Farbwerkzeug setzt ausgewählte Farben nicht zurück.");
    if (resetFaceColors_->isEnabled() || revision_ != revision || mesh_ != mesh ||
        simulationDock_->result() != result || simulationReady() != current)
        throw std::runtime_error("Zurücksetzen der Flächenfarbe verändert den Rechenzustand.");
    activeFaceColors() = colors;
    meshColor_->setCurrentIndex(colorMode);
    select(selection);
    updateDisplay();
    dirty_ = dirty;
    updateTitle();
}

void MainWindow::exerciseProjectRoundTrip() {
    exerciseFaceColorEditing();
    QTemporaryDir directory;
    if (!directory.isValid())
        throw std::runtime_error("Temporäres Verzeichnis fehlt.");
    const double edge = meshPanel_->settings().targetEdgeLength;
    const auto path = directory.filePath("roundtrip.hsproj");
    writeProject(path, projectState());
    const auto after = readProject(path);
    const DraftGeometry reloaded(after.parameters);
    if ((geometry_ && std::abs(geometry_->radius(geometry_->length() * .73, .91) -
                               reloaded.radius(reloaded.length() * .73, .91)) > 1e-12) ||
        after.meshing.targetEdgeLength != edge)
        throw std::runtime_error("Projekt-Roundtrip verändert Geometrie.");
    if (preview_) {
        std::set<std::uint32_t> ids;
        for (const auto &triangle : preview_->triangles)
            ids.insert(triangle.faceId);
        if (ids.size() >= 2) {
            const std::vector<std::uint32_t> selection{*ids.begin(), *std::next(ids.begin())};
            const auto originalTargets = activeFaceTargets();
            meshPanel_->selectFace(selection[0]);
            meshPanel_->selectFace(selection[1], Qt::ControlModifier);
            if (viewport_->selectedFaces() != selection ||
                (cadMode() && cadDock_->selectedFaces() != selection))
                throw std::runtime_error("Netz-Flächenauswahl und Modell sind nicht synchron.");
            auto *local = meshPanel_->findChild<QDoubleSpinBox *>("meshLocalTarget");
            auto *assign = meshPanel_->findChild<QPushButton *>("meshAssignLocalTarget");
            auto *reset = meshPanel_->findChild<QPushButton *>("meshResetLocalTarget");
            const auto settings = meshPanel_->settings();
            local->setValue((settings.adaptive.minimumEdgeLength + settings.targetEdgeLength) *
                            500);
            const double requested = local->value() * .001;
            assign->click();
            auto expected = originalTargets;
            for (const auto id : selection)
                expected[id] = requested;
            if (activeFaceTargets() != expected ||
                (originalTargets != expected && (!dirty_ || meshReady())))
                throw std::runtime_error("Lokale Zielkanten werden nicht gemeinsam übernommen.");
            const auto localPath = directory.filePath("local-sizing.hsproj");
            writeProject(localPath, projectState());
            const auto saved = readProject(localPath);
            if (saved.draftFaceTargetEdgeLengths != draftFaceTargetEdgeLengths_ ||
                saved.cadFaceTargetEdgeLengths != cadFaceTargetEdgeLengths_)
                throw std::runtime_error("Lokale Zielkanten gehen beim Speichern verloren.");
            reset->click();
            for (const auto id : selection)
                if (activeFaceTargets().contains(id))
                    throw std::runtime_error("Lokale Zielkanten werden nicht zurückgesetzt.");
            activeFaceTargets() = originalTargets;
            updateMeshFaces();
            meshingChanged();
            meshPanel_->setSelectedFaces({});
            viewport_->setSelectedFaces({});
            if (cadMode())
                cadDock_->clearSelection();
        }
    }
    if (!cadMode()) {
        flushEdit();
        auto *length = parameters_->findChild<QDoubleSpinBox *>("length");
        const double before = length->value();
        length->setValue(before + 1);
        const auto editedRevision = revision_;
        undoEdit();
        if (length->value() != before || revision_ <= editedRevision)
            throw std::runtime_error(
                "Draft-Rückgängig verändert keine monotone Geometrierevision.");
        const auto undoRevision = revision_;
        redoEdit();
        if (length->value() != before + 1 || revision_ <= undoRevision)
            throw std::runtime_error(
                "Draft-Wiederholen stellt Parameter nicht vollständig wieder her.");
    }
}

void MainWindow::exerciseResultProjectRoundTrip() {
    if (!simulationReady() || simulationInProgress())
        throw std::runtime_error("Aktuelles berechnetes Ergebnis für Projektprüfung erforderlich.");
    exerciseFaceColorEditing();
    const auto original = projectState();
    QTemporaryDir directory;
    if (!directory.isValid())
        throw std::runtime_error("Temporäres Verzeichnis fehlt.");
    const auto saved = directory.filePath(QString::fromUtf8("akustik-ü.hsproj"));
    const auto again = directory.filePath("again.hsproj");
    writeProject(saved, original);
    loadProject(saved);
    if (!simulationReady() || simulationInProgress() || !meshReady() || dirty_)
        throw std::runtime_error("Aktuelles Projektergebnis wurde nicht ohne Berechnung geladen.");
    auto *normalization = simulationDock_->findChild<QComboBox *>("polarNormalization");
    const int originalNormalization = normalization->currentIndex();
    const auto computed = simulationDock_->result();
    normalization->setCurrentIndex(1 - originalNormalization);
    if (!simulationReady() || simulationInProgress() || simulationDock_->result() != computed ||
        simulationDock_->archive().settings.normalizeToOnAxis ==
            original.simulation.settings.normalizeToOnAxis)
        throw std::runtime_error(
            "Normierungswechsel verändert oder invalidiert numerische Ergebnisse.");
    const auto normalizedPath = directory.filePath("normalized.hsproj");
    writeProject(normalizedPath, projectState());
    loadProject(normalizedPath);
    if (!simulationReady() || simulationDock_->archive().settings.normalizeToOnAxis ==
                                  original.simulation.settings.normalizeToOnAxis)
        throw std::runtime_error("Normierungswahl geht beim Speichern verloren.");
    loadProject(saved);
    writeProject(again, projectState());
    QFile first(saved), second(again);
    if (!first.open(QIODevice::ReadOnly) || !second.open(QIODevice::ReadOnly) ||
        first.readAll() != second.readAll())
        throw std::runtime_error("Projekt-Roundtrip verändert Netz, Request oder Ergebnisdaten.");
    const auto preserved = projectState();
    bool rejected = false;
    try {
        loadProject(directory.filePath("missing.hsproj"));
    } catch (const std::exception &) {
        rejected = true;
    }
    if (!rejected || simulationDock_->result() != preserved.simulation.result ||
        mesh_ != preserved.mesh || revision_ != preserved.revision || dirty_)
        throw std::runtime_error("Öffnungsfehler verändert das aktuelle Projekt.");

    simulationDock_->setFrequencySweep(400, 4000, 5000);
    if (simulationReady() || !dirty_)
        throw std::runtime_error("Geänderte Einstellungen markieren das Ergebnis nicht veraltet.");
    const auto stale = directory.filePath("stale.hsproj");
    writeProject(stale, projectState());
    loadProject(stale);
    const auto restored = simulationDock_->archive();
    if (simulationReady() || !restored.result || !restored.mesh || !restored.request ||
        restored.request->frequencies != original.simulation.request->frequencies ||
        restored.settings.frequencyCount != 5000 || dirty_)
        throw std::runtime_error("Veraltetes Ergebnis verliert seinen ursprünglichen Request.");
    applyProject(Project{});
    const auto cleared = simulationDock_->archive();
    if (cleared.result || cleared.mesh || cleared.request || mesh_ || cadDock_->snapshot())
        throw std::runtime_error("Neues Projekt behält alte Netz-/CAD-/Ergebnisdaten.");
    loadProject(saved);
    if (!simulationReady() || simulationInProgress())
        throw std::runtime_error(
            "Gespeichertes Ergebnis wurde bei Wiederherstellung neu gerechnet.");
    simulationDock_->show();
    simulationDock_->raise();
}

void MainWindow::exerciseLoadedProjectRoundTrip() {
    const auto original = projectState();
    if (simulationInProgress() || !geometryReady() || meshBusy_)
        throw std::runtime_error("Projektimport hat unerwartet Hintergrundberechnungen gestartet.");
    if (cadMode() && (!cadSnapshot_ || cadDock_->snapshot() != original.cad))
        throw std::runtime_error("Gespeicherter CAD-Stand fehlt im Offline-Dock.");
    if (cadMode()) {
        const auto *missing = cadDock_->findChild<QPushButton *>("cadConfirmMissingFaces");
        if (!missing || missing->isHidden() != cadSnapshot_->unresolvedPreviousFaces.empty())
            throw std::runtime_error(
                "Bestätigung entfallener Flächen passt nicht zum gespeicherten CAD-Stand.");
    }
    if (simulationReady() != original.simulation.current)
        throw std::runtime_error(
            "Gespeicherter Ergebnisstatus stimmt nicht mit der Anzeige überein.");
    if (!activeFaceColors().empty() &&
        meshColor_->currentIndex() != static_cast<int>(MeshColorMode::faceColors))
        throw std::runtime_error("Gespeicherte Flächenfarben sind nach dem Öffnen nicht sichtbar.");
    exerciseFaceColorEditing();
    QTemporaryDir directory;
    if (!directory.isValid())
        throw std::runtime_error("Temporäres Verzeichnis fehlt.");
    const auto firstPath = directory.filePath("import.hsproj");
    const auto secondPath = directory.filePath("reload.hsproj");
    writeProject(firstPath, original);
    loadProject(firstPath);
    writeProject(secondPath, projectState());
    QFile first(firstPath), second(secondPath);
    if (!first.open(QIODevice::ReadOnly) || !second.open(QIODevice::ReadOnly))
        throw std::runtime_error("Gespeicherte Offline-Prüfprojekte sind nicht lesbar.");
    const auto firstBytes = first.readAll(), secondBytes = second.readAll();
    if (firstBytes != secondBytes) {
        qsizetype mismatch = 0;
        while (mismatch < std::min(firstBytes.size(), secondBytes.size()) &&
               firstBytes[mismatch] == secondBytes[mismatch])
            ++mismatch;
        const auto start = std::max<qsizetype>(0, mismatch - 40);
        throw std::runtime_error(
            QString("Offline-Projektimport verändert Daten bei Byte %1: [%2] → [%3].")
                .arg(mismatch)
                .arg(QString::fromUtf8(firstBytes.mid(start, 120)),
                     QString::fromUtf8(secondBytes.mid(start, 120)))
                .toStdString());
    }
    if (simulationInProgress() || meshBusy_ || dirty_)
        throw std::runtime_error(
            "Offline-Projektimport startet Berechnungen oder markiert das Projekt geändert.");
    std::set<std::uint32_t> interfaceFaces;
    for (const auto &interface : original.simulation.settings.acousticModel.interfaces)
        interfaceFaces.insert(interface.faceIds.begin(), interface.faceIds.end());
    std::vector<std::uint32_t> selected;
    if (cadMode())
        for (const auto &face : cadSnapshot_->faces)
            if (!interfaceFaces.contains(face.id) && selected.size() < 2)
                selected.push_back(face.id);
    // The group-assignment exercise requires two independently assignable
    // material faces. The import and temporary round-trip above still run for
    // models with fewer such faces.
    if (selected.size() == 2) {
        std::sort(selected.begin(), selected.end());
        cadDock_->selectFace(selected[0]);
        cadDock_->selectFace(selected[1], Qt::ControlModifier);
        if (cadDock_->selectedFaces() != selected || viewport_->selectedFaces() != selected ||
            meshPanel_->selectedFaces() != selected)
            throw std::runtime_error("CAD-Mehrfachauswahl und Modell sind nicht synchron.");
        auto *table = cadDock_->findChild<QTreeWidget *>("cadFaces");
        if (!table || table->columnCount() != 5)
            throw std::runtime_error("Lokale CAD-Netzspalte fehlt.");
        const auto meshIndex = table->model()->index(0, 3);
        auto *delegate = table->itemDelegate();
        const auto editTarget = [&](double millimetres) {
            std::unique_ptr<QWidget> editor(
                delegate->createEditor(table, QStyleOptionViewItem{}, meshIndex));
            auto *input = qobject_cast<QDoubleSpinBox *>(editor.get());
            if (!input)
                throw std::runtime_error("Lokale CAD-Netzgröße ist nicht editierbar.");
            delegate->setEditorData(input, meshIndex);
            input->setValue(millimetres);
            delegate->setModelData(input, table->model(), meshIndex);
        };
        const auto meshSettings = meshPanel_->settings();
        const double target =
            (meshSettings.adaptive.minimumEdgeLength + meshSettings.targetEdgeLength) * .5;
        editTarget(target * 1000);
        // The table editor rounds to four decimal places in millimetres.
        const double committed = std::round(target * 1e7) / 1e7;
        for (const auto id : selected)
            if (!cadFaceTargetEdgeLengths_.contains(id) ||
                std::abs(cadFaceTargetEdgeLengths_.at(id) - committed) > 1e-12)
                throw std::runtime_error("CAD-Netzspalte verliert Mehrfachzuweisung.");
        const auto localPath = directory.filePath("local-cad.hsproj");
        writeProject(localPath, projectState());
        if (readProject(localPath).cadFaceTargetEdgeLengths != cadFaceTargetEdgeLengths_)
            throw std::runtime_error("Lokale CAD-Netzgrößen werden nicht gespeichert.");
        const auto assignedTargets = cadFaceTargetEdgeLengths_;
        editTarget((meshSettings.targetEdgeLength + 1) * 1000);
        if (cadFaceTargetEdgeLengths_ != assignedTargets)
            throw std::runtime_error("CAD-Netzspalte akzeptiert ungültige Zielkante.");
        editTarget(0);
        for (const auto id : selected)
            if (cadFaceTargetEdgeLengths_.contains(id))
                throw std::runtime_error("CAD-Netzspalte setzt globale Zielkante nicht zurück.");
        undoEdit();
        if (cadFaceTargetEdgeLengths_ != assignedTargets ||
            std::abs(table->topLevelItem(0)->data(3, Qt::UserRole).toDouble() - committed) > 1e-12)
            throw std::runtime_error(
                "CAD-Netzgröße wird nicht in Modell und Tabelle zurückgenommen.");
        redoEdit();
        for (const auto id : selected)
            if (cadFaceTargetEdgeLengths_.contains(id))
                throw std::runtime_error("CAD-Netzgröße wird nicht wiederholt.");
        loadProject(firstPath);
        cadDock_->setSelectedFaces(selected);
        QComboBox *boundary = nullptr;
        for (auto *candidate : cadDock_->findChildren<QComboBox *>())
            if (candidate->count() == 2 && candidate->itemText(0) == "Wall" &&
                candidate->itemText(1) == "Driven")
                boundary = candidate;
        auto *real = cadDock_->findChild<QDoubleSpinBox *>("cadVelocityReal");
        auto *imaginary = cadDock_->findChild<QDoubleSpinBox *>("cadVelocityImaginary");
        auto *assign = cadDock_->findChild<QPushButton *>("cadAssignBoundary");
        if (!boundary || !real || !imaginary || !assign)
            throw std::runtime_error("CAD-Gruppenzuweisung ist nicht bedienbar.");
        boundary->setCurrentIndex(1);
        real->setValue(.2);
        imaginary->setValue(-.3);
        assign->click();
        const auto updated = cadDock_->snapshot();
        for (const auto &face : updated->faces) {
            if (std::binary_search(selected.begin(), selected.end(), face.id) &&
                (face.boundary != BoundaryCondition::driven ||
                 face.normalVelocity != std::complex<double>(.2, -.3)))
                throw std::runtime_error(
                    "CAD-Gruppenzuweisung ändert nicht alle gewählten Flächen.");
        }
        if (cadDock_->selectedFaces() != selected || viewport_->selectedFaces() != selected ||
            meshPanel_->selectedFaces() != selected || !dirty_ || simulationInProgress())
            throw std::runtime_error(
                "CAD-Gruppenzuweisung verliert Auswahl oder startet einen Solve.");
        const auto assignmentRevision = revision_;
        undoEdit();
        const auto reverted = cadDock_->snapshot();
        for (const auto &face : reverted->faces)
            if (std::binary_search(selected.begin(), selected.end(), face.id)) {
                const auto previous =
                    std::find_if(original.cad->faces.begin(), original.cad->faces.end(),
                                 [&](const auto &candidate) { return candidate.id == face.id; });
                if (previous == original.cad->faces.end() || face.boundary != previous->boundary ||
                    face.normalVelocity != previous->normalVelocity ||
                    face.state != previous->state)
                    throw std::runtime_error("CAD-Rückgängig verliert ursprüngliche Zuordnungen.");
            }
        if (revision_ <= assignmentRevision || reverted->revision != original.cad->revision)
            throw std::runtime_error("CAD-Rückgängig darf keinen anderen Geometriestand laden.");
        redoEdit();
        for (const auto &face : cadDock_->snapshot()->faces)
            if (std::binary_search(selected.begin(), selected.end(), face.id) &&
                (face.boundary != BoundaryCondition::driven ||
                 face.normalVelocity != std::complex<double>(.2, -.3)))
                throw std::runtime_error(
                    "CAD-Wiederholen stellt Gruppenänderung nicht wieder her.");
        cadDock_->clearSelection();
        if (!viewport_->selectedFaces().empty() || !meshPanel_->selectedFaces().empty())
            throw std::runtime_error("CAD-Auswahl im Modell wird nicht geleert.");
        loadProject(firstPath);
        auto *growth = meshPanel_->findChild<QDoubleSpinBox *>("meshSizeGrowth");
        if (!growth)
            throw std::runtime_error("Adaptive CAD-Netzeinstellung fehlt.");
        const double requestedGrowth = growth->value() == .45 ? .35 : .45;
        growth->setValue(requestedGrowth);
        const auto changed = projectState();
        if (!dirty_ || changed.meshing.adaptive.growthRate != requestedGrowth ||
            changed.revision == original.revision || simulationInProgress() || meshBusy_)
            throw std::runtime_error("Adaptive CAD-Netzeinstellung wird nicht übernommen.");
        const auto adaptivePath = directory.filePath("adaptive.hsproj");
        writeProject(adaptivePath, changed);
        loadProject(adaptivePath);
        if (meshPanel_->settings().adaptive.growthRate != requestedGrowth || dirty_ ||
            simulationInProgress() || meshBusy_)
            throw std::runtime_error("Adaptive CAD-Netzeinstellung wird nicht wiederhergestellt.");
        loadProject(firstPath);
    }
    if (original.simulation.result) {
        simulationDock_->show();
        simulationDock_->raise();
    }
}

EditState MainWindow::editState() const {
    EditState state;
    state.parameters = parameters_->parameters();
    state.assembly = assembly_;
    state.meshing = meshPanel_->inputSettings();
    state.draftFaceTargets = draftFaceTargetEdgeLengths_;
    state.cadFaceTargets = cadFaceTargetEdgeLengths_;
    state.draftReversedFaceIds = draftReversedFaceIds_;
    state.cadReversedFaceIds = cadReversedFaceIds_;
    state.couplingNetwork = lemEditor_->network();
    state.couplingEnabled = lemEditor_->couplingEnabled();
    state.acousticModel = simulationDock_->archive().settings.acousticModel;
    state.symmetryXYEnabled = simulationDock_->archive().settings.symmetryXYEnabled;
    state.symmetryPlaneTolerance = simulationDock_->archive().settings.symmetryPlaneTolerance;
    if (const auto cad = cadDock_->snapshot()) {
        state.cad = CadAssignmentState{
            cad->documentId, cad->sessionId, cad->revision, {}, cad->unresolvedPreviousFaces};
        for (const auto &face : cad->faces)
            state.cad->faces.push_back({face.id, face.boundary, face.normalVelocity, face.state});
        std::sort(state.cad->faces.begin(), state.cad->faces.end(),
                  [](const auto &a, const auto &b) { return a.id < b.id; });
    }
    return state;
}
void MainWindow::rememberEdit(bool immediate) {
    if (loading_ || historyReplaying_ || !undoStack_)
        return;
    if (immediate)
        flushEdit();
    EditState state;
    try {
        state = editState();
    } catch (const std::invalid_argument &) {
        // An unfinished LUT cell is not a committed numeric model edit.
        return;
    }
    if (!historyReady_ || !sameCadGeometry(historyState_.cad, state.cad)) {
        resetEditHistory();
        return;
    }
    if (sameEditState(historyState_, state)) {
        editTimer_->stop();
        pendingEdit_.reset();
    } else {
        pendingEdit_ = std::move(state);
        if (immediate)
            flushEdit();
        else
            editTimer_->start();
    }
    updateEditActions();
}
void MainWindow::flushEdit() {
    if (!pendingEdit_ || historyReplaying_)
        return;
    editTimer_->stop();
    auto before = historyState_;
    auto after = std::move(*pendingEdit_);
    pendingEdit_.reset();
    historyState_ = after;
    undoStack_->push(new EditCommand(std::move(before), std::move(after),
                                     [this](const EditState &state) { applyEditState(state); }));
    updateEditActions();
}
void MainWindow::resetEditHistory() {
    if (!undoStack_)
        return;
    editTimer_->stop();
    pendingEdit_.reset();
    historyReady_ = false;
    undoStack_->clear();
    try {
        historyState_ = editState();
        historyReady_ = true;
    } catch (const std::invalid_argument &) {
        // No numeric history anchor is available until the LUT edit is valid.
    }
    updateEditActions();
}
void MainWindow::applyEditState(const EditState &state) {
    auto withoutCoupling = editState();
    const bool regionsChanged = withoutCoupling.acousticModel != state.acousticModel;
    const bool symmetryChanged =
        withoutCoupling.symmetryXYEnabled != state.symmetryXYEnabled ||
        withoutCoupling.symmetryPlaneTolerance != state.symmetryPlaneTolerance;
    withoutCoupling.couplingNetwork = state.couplingNetwork;
    withoutCoupling.couplingEnabled = state.couplingEnabled;
    withoutCoupling.acousticModel = state.acousticModel;
    withoutCoupling.symmetryXYEnabled = state.symmetryXYEnabled;
    withoutCoupling.symmetryPlaneTolerance = state.symmetryPlaneTolerance;
    const bool onlyCoupling = sameEditState(withoutCoupling, state);
    const auto current = cadDock_->snapshot();
    if (state.cad.has_value() != static_cast<bool>(current) ||
        (state.cad && (state.cad->documentId != current->documentId ||
                       state.cad->sessionId != current->sessionId ||
                       state.cad->revision != current->revision || cadDock_->isGeometryPending())))
        throw std::invalid_argument(
            "CAD-Geometrie wurde geändert; Bearbeitungshistorie ist nicht mehr anwendbar.");
    ParameterPanel checkedParameters;
    checkedParameters.setParameters(state.parameters);
    MeshPanel checkedMesh;
    checkedMesh.setInputSettings(state.meshing);
    const QScopedValueRollback<bool> replay(historyReplaying_, true);
    {
        const QScopedValueRollback<bool> loading(loading_, true);
        lemEditor_->setNetwork(state.couplingNetwork);
        lemEditor_->setCouplingEnabled(state.couplingEnabled);
    }
    simulationDock_->setCoupling(state.couplingNetwork, state.couplingEnabled);
    simulationDock_->setAcousticModel(state.acousticModel);
    simulationDock_->setSymmetryXY(state.symmetryXYEnabled, state.symmetryPlaneTolerance);
    updateInterfacePresentation();
    if (onlyCoupling) {
        if (gmsh_ && (regionsChanged || symmetryChanged))
            meshingChanged();
        markDirty();
        historyState_ = state;
        historyReady_ = true;
        pendingEdit_.reset();
        editTimer_->stop();
        updateEditActions();
        return;
    }
    {
        const QScopedValueRollback<bool> loading(loading_, true);
        if (state.cad &&
            !cadDock_->setAssignments(state.cad->documentId, state.cad->sessionId,
                                      state.cad->revision, state.cad->faces, state.cad->missing))
            throw std::invalid_argument("CAD-Zuordnungen konnten nicht zurückgesetzt werden.");
        parameters_->setParameters(state.parameters);
        assembly_ = state.assembly;
        meshPanel_->setInputSettings(state.meshing);
        draftFaceTargetEdgeLengths_ = state.draftFaceTargets;
        cadFaceTargetEdgeLengths_ = state.cadFaceTargets;
        draftReversedFaceIds_ = state.draftReversedFaceIds;
        cadReversedFaceIds_ = state.cadReversedFaceIds;
        cadSnapshot_ = cadDock_->snapshot();
        if (cadMode() && cadSnapshot_)
            rawPreview_ = std::shared_ptr<const SurfaceMesh>(cadSnapshot_, &cadSnapshot_->mesh);
        updateOrientedPreview();
        meshPanel_->setFaceNotice({});
        updateSourcePresentation();
    }
    if (!cadMode())
        geometryPending_ = true;
    meshingChanged();
    historyState_ = state;
    historyReady_ = true;
    pendingEdit_.reset();
    editTimer_->stop();
    updateEditActions();
}
void MainWindow::updateEditActions() {
    if (!undoAction_ || !redoAction_)
        return;
    undoAction_->setEnabled(pendingEdit_.has_value() || undoStack_->canUndo());
    redoAction_->setEnabled(!pendingEdit_ && undoStack_->canRedo());
}
void MainWindow::undoEdit() {
    if (loading_ || closing_)
        return;
    if (cadDock_->isGeometryPending()) {
        setDiagnostic(tr("CAD-Geometrie wird geändert · Rückgängig wartet auf den neuen Stand."));
        return;
    }
    commitFocusedEdit();
    flushEdit();
    try {
        undoStack_->undo();
    } catch (const std::exception &error) {
        resetEditHistory();
        setDiagnostic(QString::fromUtf8(error.what()), true);
    }
    updateEditActions();
}
void MainWindow::redoEdit() {
    if (loading_ || closing_)
        return;
    if (cadDock_->isGeometryPending()) {
        setDiagnostic(tr("CAD-Geometrie wird geändert · Wiederholen wartet auf den neuen Stand."));
        return;
    }
    commitFocusedEdit();
    flushEdit();
    try {
        undoStack_->redo();
    } catch (const std::exception &error) {
        resetEditHistory();
        setDiagnostic(QString::fromUtf8(error.what()), true);
    }
    updateEditActions();
}
bool MainWindow::eventFilter(QObject *object, QEvent *event) {
    if (event->type() == QEvent::KeyPress && !QApplication::activeModalWidget()) {
        auto *focus = QApplication::focusWidget();
        const auto *key = static_cast<QKeyEvent *>(event);
        if (focus && (focus == this || isAncestorOf(focus))) {
            if (key->key() == Qt::Key_F1 && key->modifiers() == Qt::NoModifier) {
                auto *field = focus;
                while (field && field != this && field->whatsThis().isEmpty() &&
                       field->toolTip().isEmpty())
                    field = field->parentWidget();
                if (field && field != this) {
                    QWhatsThis::showText(field->mapToGlobal(QPoint(0, field->height())),
                                         field->whatsThis().isEmpty() ? field->toolTip()
                                                                      : field->whatsThis(),
                                         field);
                    return true;
                }
            }
            if (key->matches(QKeySequence::Undo)) {
                undoEdit();
                return true;
            }
            if (key->matches(QKeySequence::Redo)) {
                redoEdit();
                return true;
            }
        }
    }
    return QMainWindow::eventFilter(object, event);
}

} // namespace hornsculpt
