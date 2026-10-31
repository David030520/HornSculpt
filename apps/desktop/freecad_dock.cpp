#include "freecad_dock.hpp"
#include "numeric_input.hpp"

#include <QColor>
#include <QComboBox>
#include <QCoreApplication>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QItemSelectionModel>
#include <QLabel>
#include <QMainWindow>
#include <QProcess>
#include <QProcessEnvironment>
#include <QPushButton>
#include <QSaveFile>
#include <QSignalBlocker>
#include <QStandardPaths>
#include <QStyledItemDelegate>
#include <QTreeWidget>
#include <QVBoxLayout>
#include <algorithm>
#include <array>
#include <stdexcept>

namespace hornsculpt {
namespace {
constexpr int meshColumn = 3;

class FaceTargetDelegate final : public QStyledItemDelegate {
  public:
    explicit FaceTargetDelegate(FreeCadDock *dock) : QStyledItemDelegate(dock), dock_(dock) {}

    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &,
                          const QModelIndex &index) const override {
        if (index.column() != meshColumn)
            return nullptr;
        auto *input = new DecimalInput(parent);
        input->setObjectName(QStringLiteral("cadFaceTargetEditor"));
        input->setDecimals(4);
        input->setRange(0, 1e6);
        input->setSingleStep(.1);
        input->setSpecialValueText(tr("global"));
        input->setSuffix(QStringLiteral(" mm"));
        input->setToolTip(tr("Lokale Zielkante in mm. 0 = globale Zielkante. "
                             "Gilt für alle ausgewählten Flächen; erlaubt zwischen "
                             "Mindestkante und globaler Zielkante im Netzpanel."));
        return input;
    }

    void setEditorData(QWidget *editor, const QModelIndex &index) const override {
        auto *input = static_cast<DecimalInput *>(editor);
        input->setValue(index.data(Qt::UserRole).toDouble() * 1000);
        input->selectAll();
    }

    void setModelData(QWidget *editor, QAbstractItemModel *,
                      const QModelIndex &index) const override {
        auto *input = static_cast<DecimalInput *>(editor);
        input->interpretText();
        const auto id = index.siblingAtColumn(0).data(Qt::UserRole).toUInt();
        auto ids = dock_->selectedFaces();
        if (!std::binary_search(ids.begin(), ids.end(), id))
            ids = {id};
        if (dock_->onFaceTargetChanged)
            dock_->onFaceTargetChanged(std::move(ids), input->value() > 0
                                                           ? std::optional(input->value() * .001)
                                                           : std::nullopt);
    }

  private:
    FreeCadDock *dock_;
};

QString freeCadUserDirectory() {
    QStringList installations;
    const auto executable = QStandardPaths::findExecutable(QStringLiteral("FreeCAD.exe"));
    if (!executable.isEmpty())
        installations.append(
            QDir::cleanPath(QFileInfo(executable).dir().absoluteFilePath(QStringLiteral(".."))));
    for (const auto *variable : {"ProgramW6432", "ProgramFiles", "ProgramFiles(x86)"}) {
        const QDir programs(qEnvironmentVariable(variable));
        if (programs.path().isEmpty() || programs.path() == QStringLiteral("."))
            continue;
        for (const auto &directory :
             programs.entryList({QStringLiteral("FreeCAD*")}, QDir::Dirs | QDir::NoDotAndDotDot))
            installations.append(programs.absoluteFilePath(directory));
    }
    installations.removeDuplicates();
    QStringList userDirectories;
    for (const auto &installation : installations) {
        const QDir root(installation);
        const auto python = root.absoluteFilePath(QStringLiteral("bin/python.exe"));
        const auto library = root.absoluteFilePath(QStringLiteral("lib"));
        if (!QFileInfo::exists(python) || !QDir(library).exists())
            continue;
        QProcess query;
        auto environment = QProcessEnvironment::systemEnvironment();
        environment.insert(QStringLiteral("PATH"), root.absoluteFilePath(QStringLiteral("bin")) +
                                                       QDir::listSeparator() +
                                                       environment.value(QStringLiteral("PATH")));
        query.setProcessEnvironment(environment);
        query.setWorkingDirectory(root.absolutePath());
        query.start(python,
                    {QStringLiteral("-c"),
                     QStringLiteral("import sys;sys.path.insert(0,sys.argv[1]);"
                                    "import FreeCAD as App;"
                                    "print('HORNSCULPT_USERAPPDATA='+App.getUserAppDataDir())"),
                     library});
        if (!query.waitForStarted(2000))
            continue;
        if (!query.waitForFinished(5000)) {
            query.kill();
            query.waitForFinished(1000);
            continue;
        }
        if (query.exitStatus() != QProcess::NormalExit || query.exitCode() != 0)
            continue;
        for (const auto &line : QString::fromUtf8(query.readAllStandardOutput()).split('\n')) {
            const QString prefix = QStringLiteral("HORNSCULPT_USERAPPDATA=");
            if (!line.startsWith(prefix))
                continue;
            const auto path = QDir::cleanPath(line.mid(prefix.size()).trimmed());
            if (QDir::isAbsolutePath(path))
                userDirectories.append(path);
        }
    }
    userDirectories.removeDuplicates();
    // Different installed FreeCAD versions may use different user module directories.
    return userDirectories.size() == 1 ? userDirectories.front() : QString{};
}

QDoubleSpinBox *velocityInput(const QString &name) {
    auto *result = new DecimalInput;
    result->setObjectName(name);
    result->setRange(-1e6, 1e6);
    result->setDecimals(9);
    result->setSingleStep(0.001);
    result->setSuffix(QStringLiteral(" m/s"));
    result->setToolTip(
        QObject::tr("Komplexe physische Normalgeschwindigkeit, RMS in m/s; exp(+iωt). "
                    "Freifeld: positiv mit der Netznormale aus dem Festkörper in die Luft. "
                    "Schallwand-Luftkavität: positiv in die Luft, entgegen der aus der Luft "
                    "gerichteten Driven-Netznormale."));
    result->setWhatsThis(result->toolTip());
    return result;
}
} // namespace

FreeCadDock::FreeCadDock(QMainWindow *parent)
    : QDockWidget(tr("FreeCAD-Link"), parent), link_(this) {
    setObjectName(QStringLiteral("freecadDock"));
    auto *content = new QWidget;
    auto *layout = new QVBoxLayout(content);
    layout->setContentsMargins(8, 8, 8, 8);
    auto *connection = new QHBoxLayout;
    connect_ = new QPushButton(tr("Verbinden"));
    install_ = new QPushButton(tr("Workbench installieren …"));
    install_->setObjectName(QStringLiteral("cadInstallWorkbench"));
    install_->setToolTip(tr("HornSculpt-Workbench im FreeCAD-Benutzerverzeichnis installieren. "
                            "Danach in FreeCAD die Workbench HornSculpt öffnen und verbinden."));
    connection->addWidget(connect_);
    connection->addWidget(install_);
    layout->addLayout(connection);
    status_ = new QLabel;
    status_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(status_);
    document_ = new QLabel;
    document_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(document_);
    faces_ = new QTreeWidget;
    faces_->setObjectName(QStringLiteral("cadFaces"));
    faces_->setHeaderLabels({tr("Fläche"), tr("Rand"), tr("Zuordnung"), tr("Netz [mm]"), tr("n")});
    faces_->setItemDelegate(new FaceTargetDelegate(this));
    faces_->setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::EditKeyPressed);
    faces_->setRootIsDecorated(false);
    faces_->setSelectionMode(QAbstractItemView::ExtendedSelection);
    faces_->setToolTip(tr("Strg: Auswahl erweitern oder umschalten · Umschalt: Bereich wählen"));
    faces_->setMinimumHeight(160);
    layout->addWidget(faces_, 1);
    selectionLabel_ = new QLabel;
    layout->addWidget(selectionLabel_);
    auto *fields = new QFormLayout;
    boundary_ = new QComboBox;
    boundary_->addItems({QStringLiteral("Wall"), QStringLiteral("Driven")});
    boundary_->setPlaceholderText(tr("Gemischt"));
    boundary_->setToolTip(tr("Wall: schallhart, Normalgeschwindigkeit null. Driven: "
                             "vorgeschriebene komplexe Normalgeschwindigkeit."));
    velocityReal_ = velocityInput(QStringLiteral("cadVelocityReal"));
    velocityImaginary_ = velocityInput(QStringLiteral("cadVelocityImaginary"));
    fields->addRow(tr("Randbedingung"), boundary_);
    fields->addRow(tr("Geschwindigkeit Re"), velocityReal_);
    fields->addRow(tr("Geschwindigkeit Im"), velocityImaginary_);
    layout->addLayout(fields);
    auto *assignments = new QHBoxLayout;
    apply_ = new QPushButton(tr("Zuweisen"));
    apply_->setObjectName(QStringLiteral("cadAssignBoundary"));
    confirm_ = new QPushButton(tr("Zuordnung bestätigen"));
    assignments->addWidget(apply_);
    assignments->addWidget(confirm_);
    layout->addLayout(assignments);
    bemRegions_ = new QPushButton(tr("BEM-Teilgebiete …"));
    bemRegions_->setObjectName(QStringLiteral("cadBemRegions"));
    bemRegions_->setToolTip(
        tr("Teilgebiete und transparente BEM-Interfaces für diese Flächen zuordnen."));
    bemRegions_->setWhatsThis(bemRegions_->toolTip());
    layout->addWidget(bemRegions_);
    connect(bemRegions_, &QPushButton::clicked, this, [this] {
        if (onEditBemRegions)
            onEditBemRegions();
    });
    mapping_ = new QLabel;
    layout->addWidget(mapping_);
    missing_ = new QPushButton(tr("Entfallene Flächen bestätigen"));
    missing_->setObjectName(QStringLiteral("cadConfirmMissingFaces"));
    missing_->setToolTip(tr("Entfallene Flächen und ihre alten Randbedingungen verwerfen."));
    layout->addWidget(missing_);
    setWidget(content);
    connect(connect_, &QPushButton::clicked, this, [this] {
        if (link_.isListening())
            link_.stop();
        else
            ensureListening();
        refresh();
        publishAvailability();
    });
    connect(install_, &QPushButton::clicked, this, [this] { installWorkbench(); });
    connect(faces_, &QTreeWidget::itemSelectionChanged, this, [this] { selectionChanged(); });
    connect(boundary_, &QComboBox::currentIndexChanged, this, [this] {
        const bool selected = !selectedFaces().empty();
        const bool driven = selected && boundary_->currentIndex() == 1;
        apply_->setEnabled(selected && boundary_->currentIndex() >= 0);
        velocityReal_->setEnabled(driven);
        velocityImaginary_->setEnabled(driven);
    });
    connect(apply_, &QPushButton::clicked, this, [this] {
        const auto ids = selectedFaces();
        if (!ids.empty() && boundary_->currentIndex() >= 0) {
            if (boundary_->currentIndex() == 1 &&
                std::any_of(ids.begin(), ids.end(),
                            [this](auto id) { return bemInterfaceFaces_.contains(id); })) {
                if (onError)
                    onError(tr("BEM-Interfaceflächen können nicht als Driven zugewiesen werden."));
                return;
            }
            if (onAssignmentEditStarted)
                onAssignmentEditStarted();
            link_.setBoundaries(ids,
                                boundary_->currentIndex() == 1 ? BoundaryCondition::driven
                                                               : BoundaryCondition::wall,
                                {velocityReal_->value(), velocityImaginary_->value()});
        }
    });
    connect(confirm_, &QPushButton::clicked, this, [this] {
        const auto ids = selectedFaces();
        if (!ids.empty()) {
            if (onAssignmentEditStarted)
                onAssignmentEditStarted();
            link_.confirmFaces(ids);
        }
    });
    connect(missing_, &QPushButton::clicked, this, [this] {
        if (onAssignmentEditStarted)
            onAssignmentEditStarted();
        link_.acknowledgeMissingFaces();
    });
    connect(&link_, &FreeCadLink::connectionChanged, this, [this] {
        if (link_.isConnected()) {
            offlineSnapshot_ = false;
            awaitingSnapshot_ = true;
        }
        refresh();
        publishAvailability();
    });
    connect(&link_, &FreeCadLink::sourceDirtyChanged, this, [this] {
        refresh();
        publishAvailability();
    });
    connect(&link_, &FreeCadLink::snapshotChanged, this, [this] {
        if (link_.isConnected())
            awaitingSnapshot_ = false;
        refresh();
        if (onSnapshotChanged)
            onSnapshotChanged(link_.snapshot());
        publishAvailability();
    });
    connect(&link_, &FreeCadLink::errorOccurred, this, [this](const QString &error) {
        status_->setToolTip(error);
        if (onError)
            onError(error);
    });
    refresh();
}

FreeCadDock::~FreeCadDock() {
    // Link shutdown and child selection teardown must precede our member destruction.
    disconnect(&link_, nullptr, this, nullptr);
    for (auto *child : findChildren<QObject *>())
        disconnect(child, nullptr, this, nullptr);
    link_.stop();
}

void FreeCadDock::ensureListening() {
    if (!link_.isListening())
        link_.listen();
    refresh();
    publishAvailability();
}
bool FreeCadDock::hasLiveLink() const {
    return link_.isListening() || link_.isConnected();
}

bool FreeCadDock::isSourceAvailable() const {
    const auto snapshot = link_.snapshot();
    return (offlineSnapshot_ || link_.isConnected()) && !awaitingSnapshot_ &&
           !link_.isGeometryPending() && snapshot && !snapshot->hasUnresolvedAssignments();
}
bool FreeCadDock::isGeometryPending() const {
    return link_.isGeometryPending();
}
bool FreeCadDock::setAssignments(const QString &documentId, const QString &sessionId,
                                 std::uint64_t revision,
                                 const std::vector<CadFaceAssignment> &faces,
                                 const std::vector<std::uint32_t> &missing) {
    return link_.setAssignments(documentId, sessionId, revision, faces, missing);
}
std::shared_ptr<const CadSnapshot> FreeCadDock::snapshot() const {
    return link_.snapshot();
}
void FreeCadDock::restoreSnapshot(std::shared_ptr<const CadSnapshot> snapshot) {
    restoreSnapshot(PreparedCadSnapshot::prepare(std::move(snapshot)));
}
void FreeCadDock::restoreSnapshot(const PreparedCadSnapshot &snapshot) {
    const QSignalBlocker linkSignals(link_);
    link_.restoreSnapshot(snapshot);
    {
        const QSignalBlocker selectionSignals(faces_);
        faces_->clearSelection();
    }
    selectionDocumentId_.clear();
    offlineSnapshot_ = static_cast<bool>(snapshot.snapshot());
    awaitingSnapshot_ = false;
    refresh();
    // The owning window installs the complete project transaction; avoid live callbacks.
    sourceAvailable_ = isSourceAvailable();
}
void FreeCadDock::resetProject() {
    reversedFaces_.clear();
    restoreSnapshot({});
}
std::vector<std::uint32_t> FreeCadDock::selectedFaces() const {
    std::vector<std::uint32_t> result;
    for (const auto *item : faces_->selectedItems())
        result.push_back(item->data(0, Qt::UserRole).toUInt());
    std::sort(result.begin(), result.end());
    return result;
}
void FreeCadDock::clearSelection() {
    setSelectedFaces({});
}
void FreeCadDock::setSelectedFaces(const std::vector<std::uint32_t> &selection) {
    if (selectedFaces() == selection)
        return;
    const QSignalBlocker blocker(faces_);
    faces_->clearSelection();
    faces_->setCurrentItem(nullptr);
    for (int row = 0; row < faces_->topLevelItemCount(); ++row) {
        auto *item = faces_->topLevelItem(row);
        item->setSelected(std::find(selection.begin(), selection.end(),
                                    item->data(0, Qt::UserRole).toUInt()) != selection.end());
    }
    selectionChanged();
}
void FreeCadDock::selectFace(std::uint32_t faceId, Qt::KeyboardModifiers modifiers) {
    for (int i = 0; i < faces_->topLevelItemCount(); ++i) {
        auto *item = faces_->topLevelItem(i);
        if (item->data(0, Qt::UserRole).toUInt() == faceId) {
            const QSignalBlocker blocker(faces_);
            const bool wasSelected = item->isSelected();
            if (!(modifiers & Qt::ControlModifier))
                faces_->clearSelection();
            faces_->setCurrentItem(item, 0, QItemSelectionModel::NoUpdate);
            item->setSelected((modifiers & Qt::ControlModifier) ? !wasSelected : true);
            faces_->scrollToItem(item);
            selectionChanged();
            return;
        }
    }
}
void FreeCadDock::publishAvailability() {
    const bool available = isSourceAvailable();
    if (available == sourceAvailable_)
        return;
    sourceAvailable_ = available;
    if (onAvailabilityChanged)
        onAvailabilityChanged(available);
}
void FreeCadDock::refresh() {
    connect_->setText(link_.isListening() ? tr("Trennen") : tr("Verbinden"));
    status_->setText(link_.isConnected()
                         ? (link_.isGeometryPending() ? tr("FreeCAD · Änderung ausstehend")
                                                      : tr("FreeCAD · verbunden"))
                     : offlineSnapshot_    ? tr("CAD · gespeicherter Stand · offline")
                     : link_.isListening() ? tr("Wartet auf FreeCAD · Port %1").arg(link_.port())
                                           : tr("FreeCAD · getrennt"));
    const auto current = link_.snapshot();
    bemRegions_->setEnabled(bool(current) && isSourceAvailable());
    if (current && current->sessionId == "gmsh-import-offline")
        status_->setText(tr("Gmsh · importierte Polygonquelle · offline"));
    const auto selection = current && current->documentId == selectionDocumentId_
                               ? selectedFaces()
                               : std::vector<std::uint32_t>{};
    std::uint32_t focused = 0;
    if (const auto *item = faces_->currentItem())
        focused = item->data(0, Qt::UserRole).toUInt();
    selectionDocumentId_ = current ? current->documentId : QString{};
    QSignalBlocker blocker(faces_);
    faces_->clear();
    if (current) {
        document_->setText(
            tr("%1 · CAD-Revision %2").arg(current->documentName).arg(current->revision));
        std::size_t reviewCount = 0;
        for (const auto &face : current->faces) {
            const QString state = face.state == CadFaceState::needsReview ? tr("prüfen")
                                  : face.state == CadFaceState::newFace   ? tr("neu")
                                                                          : tr("bestätigt");
            auto *item = new QTreeWidgetItem(faces_, {face.label,
                                                      face.boundary == BoundaryCondition::driven
                                                          ? QStringLiteral("Driven")
                                                          : QStringLiteral("Wall"),
                                                      state});
            item->setData(0, Qt::UserRole, face.id);
            item->setFlags(item->flags() | Qt::ItemIsEditable);
            item->setToolTip(0, tr("Flächen-ID %1 · %2 · %3 mm²")
                                    .arg(face.id)
                                    .arg(face.signature.surfaceType)
                                    .arg(face.signature.area * 1e6, 0, 'f', 2));
            if (face.state == CadFaceState::needsReview) {
                item->setForeground(2, QColor(151, 91, 23));
                ++reviewCount;
            }
            if (std::binary_search(selection.begin(), selection.end(), face.id)) {
                if (focused == face.id)
                    faces_->setCurrentItem(item, 0, QItemSelectionModel::NoUpdate);
                item->setSelected(true);
            }
        }
        const auto source =
            current->sessionId == "gmsh-import-offline"
                ? tr("Gmsh · Polygonquelle")
                : (current->geometry ? tr("CAD · BRep") : tr("CAD · Polygonquelle"));
        mapping_->setText(source + QStringLiteral(" · ") +
                          (reviewCount || !current->unresolvedPreviousFaces.empty()
                               ? tr("%1 zu prüfen · %2 entfallen")
                                     .arg(reviewCount)
                                     .arg(current->unresolvedPreviousFaces.size())
                               : tr("%1 Flächen").arg(current->faces.size())));
        missing_->setEnabled(!current->unresolvedPreviousFaces.empty());
        missing_->setVisible(!current->unresolvedPreviousFaces.empty());
    } else {
        document_->setText(tr("Kein Geometriestand"));
        mapping_->clear();
        missing_->setEnabled(false);
        missing_->setVisible(false);
    }
    faces_->resizeColumnToContents(1);
    faces_->resizeColumnToContents(2);
    setFaceTargets(faceTargets_);
    setReversedFaces(reversedFaces_);
    setInterfaceFaces(interfaceFaces_, interfacesActive_);
    selectionChanged();
}
void FreeCadDock::setBemInterfaceFaces(std::set<std::uint32_t> faces) {
    bemInterfaceFaces_ = std::move(faces);
    setInterfaceFaces(interfaceFaces_, interfacesActive_);
}
void FreeCadDock::setInterfaceFaces(std::set<std::uint32_t> faces, bool active) {
    interfaceFaces_ = std::move(faces);
    interfacesActive_ = active;
    const auto current = link_.snapshot();
    if (!current)
        return;
    const QSignalBlocker blocker(faces_);
    for (int row = 0; row < faces_->topLevelItemCount(); ++row) {
        auto *item = faces_->topLevelItem(row);
        const auto id = item->data(0, Qt::UserRole).toUInt();
        const auto found = std::find_if(current->faces.begin(), current->faces.end(),
                                        [id](const auto &face) { return face.id == id; });
        if (found == current->faces.end())
            continue;
        const bool interface = interfaceFaces_.contains(id);
        const bool bemInterface = bemInterfaceFaces_.contains(id);
        item->setText(
            1, bemInterface
                   ? (found->boundary == BoundaryCondition::wall ? tr("Interface ↔ BEM")
                                                                 : tr("Interface · Tagkonflikt"))
               : interface ? (found->boundary == BoundaryCondition::wall
                                  ? (active ? tr("Interface ↔ LEM") : tr("Interface ↔ LEM · aus"))
                                  : tr("Interface · Tagkonflikt"))
               : found->boundary == BoundaryCondition::driven ? QStringLiteral("Driven")
                                                              : QStringLiteral("Wall"));
        item->setToolTip(
            1,
            bemInterface ? tr("Transparente Grenzfläche zwischen zwei BEM-Gebieten; Basistag Wall, "
                              "keine Driven-Fläche.")
            : interface  ? tr("Virtuelle Rechenfläche zwischen BEM-Frontluft und LEM-Innenvolumen. "
                               "Das geometrische Basistag bleibt Wall; bei aktiver Kopplung ist "
                               "dies keine schallharte Materialwand. Flächenauswahl im "
                               "LEM-Anschluss hebt das Interface im 3D-Modell hervor.")
                         : QString{});
        item->setForeground(1, bemInterface
                                   ? QBrush(QColor(0, 128, 150))
                                   : (interface ? QBrush(QColor(126, 65, 160)) : QBrush()));
    }
    faces_->resizeColumnToContents(1);
}
void FreeCadDock::setReversedFaces(const std::set<std::uint32_t> &faces) {
    reversedFaces_ = faces;
    const QSignalBlocker blocker(faces_);
    for (int row = 0; row < faces_->topLevelItemCount(); ++row) {
        auto *item = faces_->topLevelItem(row);
        const bool reversed = reversedFaces_.contains(item->data(0, Qt::UserRole).toUInt());
        item->setText(4, reversed ? QStringLiteral("−") : QStringLiteral("+"));
        item->setToolTip(4, reversed ? tr("Gegenüber der CAD-Quelle umgekehrt")
                                     : tr("Orientierung der CAD-Quelle"));
    }
    faces_->resizeColumnToContents(4);
}
void FreeCadDock::setFaceTargets(const std::map<std::uint32_t, double> &targets) {
    faceTargets_ = targets;
    const QSignalBlocker blocker(faces_);
    for (int row = 0; row < faces_->topLevelItemCount(); ++row) {
        auto *item = faces_->topLevelItem(row);
        const auto found = faceTargets_.find(item->data(0, Qt::UserRole).toUInt());
        const double target = found == faceTargets_.end() ? 0 : found->second;
        item->setData(meshColumn, Qt::UserRole, target);
        item->setText(meshColumn,
                      target > 0 ? locale().toString(target * 1000, 'f', 4) : tr("global"));
        item->setToolTip(meshColumn, tr("Lokale Zielkante in mm · Doppelklick oder F2 zum Ändern. "
                                        "Mehrfachauswahl: Wert gilt für alle ausgewählten Flächen. "
                                        "0 = global. Grenzen im gemeinsamen Netzpanel."));
    }
    faces_->resizeColumnToContents(meshColumn);
}
void FreeCadDock::selectionChanged() {
    const auto ids = selectedFaces();
    const auto current = link_.snapshot();
    const CadFace *selected = nullptr;
    bool needsReview = false;
    bool sameBoundary = true;
    bool sameVelocity = true;
    if (!ids.empty() && current) {
        const auto *focused = faces_->currentItem();
        const auto focusedId = focused && focused->isSelected()
                                   ? focused->data(0, Qt::UserRole).toUInt()
                                   : ids.front();
        const auto found = std::find_if(current->faces.begin(), current->faces.end(),
                                        [&](const CadFace &face) { return face.id == focusedId; });
        if (found != current->faces.end())
            selected = &*found;
        for (const auto &face : current->faces)
            if (std::binary_search(ids.begin(), ids.end(), face.id)) {
                needsReview |= face.state == CadFaceState::needsReview;
                sameBoundary &= selected && face.boundary == selected->boundary;
                sameVelocity &= selected && face.normalVelocity == selected->normalVelocity;
            }
    }
    const bool mixed = !sameBoundary || !sameVelocity;
    selectionLabel_->setText(ids.empty() ? tr("Keine Auswahl")
                             : mixed     ? tr("%1 ausgewählt · Werte gemischt").arg(ids.size())
                                         : tr("%1 ausgewählt").arg(ids.size()));
    selectionLabel_->setToolTip(
        mixed ? tr("Die Eingabefelder zeigen Vorschlagswerte der fokussierten Fläche. Zuweisen "
                   "setzt die neuen Werte für die gesamte Auswahl.")
              : QString{});
    apply_->setEnabled(selected && sameBoundary);
    apply_->setText(ids.size() > 1 ? tr("%1 Flächen zuweisen").arg(ids.size()) : tr("Zuweisen"));
    confirm_->setEnabled(needsReview);
    boundary_->setEnabled(selected != nullptr);
    if (selected) {
        const QSignalBlocker blocker(boundary_);
        boundary_->setCurrentIndex(!sameBoundary                                     ? -1
                                   : selected->boundary == BoundaryCondition::driven ? 1
                                                                                     : 0);
        velocityReal_->setValue(selected->normalVelocity.real());
        velocityImaginary_->setValue(selected->normalVelocity.imag());
    }
    const bool driven = selected && boundary_->currentIndex() == 1;
    velocityReal_->setEnabled(driven);
    velocityImaginary_->setEnabled(driven);
    if (onFacesSelected)
        onFacesSelected(ids);
}
void FreeCadDock::installWorkbench() {
    const QDir source(QCoreApplication::applicationDirPath() +
                      QStringLiteral("/integrations/freecad/workbench"));
    constexpr std::array<const char *, 4> names{"Init.py", "InitGui.py", "hornsculpt_commands.py",
                                                "bridge.py"};
    std::array<QByteArray, names.size()> contents;
    try {
        for (std::size_t i = 0; i < names.size(); ++i) {
            QFile file(source.absoluteFilePath(QString::fromLatin1(names[i])));
            if (!file.open(QIODevice::ReadOnly))
                throw std::runtime_error(tr("Workbench-Datei fehlt oder ist nicht lesbar: %1")
                                             .arg(file.fileName())
                                             .toStdString());
            contents[i] = file.readAll();
            if (file.error() != QFileDevice::NoError)
                throw std::runtime_error(file.errorString().toStdString());
        }
        auto userDirectory = freeCadUserDirectory();
        if (userDirectory.isEmpty())
            userDirectory = QFileDialog::getExistingDirectory(
                this, tr("FreeCAD-Benutzerverzeichnis wählen (enthält Mod)"),
                QStandardPaths::writableLocation(QStandardPaths::GenericDataLocation) +
                    QStringLiteral("/FreeCAD"));
        if (userDirectory.isEmpty())
            return;
        const QDir target(QDir(userDirectory).absoluteFilePath(QStringLiteral("Mod/HornSculpt")));
        if (!QDir().mkpath(target.absolutePath()))
            throw std::runtime_error(tr("Workbench-Verzeichnis kann nicht angelegt werden: %1")
                                         .arg(target.absolutePath())
                                         .toStdString());
        for (std::size_t i = 0; i < names.size(); ++i) {
            QSaveFile file(target.absoluteFilePath(QString::fromLatin1(names[i])));
            if (!file.open(QIODevice::WriteOnly) || file.write(contents[i]) != contents[i].size() ||
                !file.commit())
                throw std::runtime_error(
                    tr("Workbench-Datei kann nicht installiert werden: %1 · %2")
                        .arg(file.fileName(), file.errorString())
                        .toStdString());
        }
        status_->setText(tr("Workbench installiert · FreeCAD neu starten"));
        status_->setToolTip(QDir::toNativeSeparators(target.absolutePath()));
    } catch (const std::runtime_error &error) {
        if (onError)
            onError(QString::fromUtf8(error.what()));
    }
}

} // namespace hornsculpt
