#include "mesh_import_dialog.hpp"
#include <QComboBox>
#include <QDialogButtonBox>
#include <QFileInfo>
#include <QFormLayout>
#include <QLabel>
#include <QListWidget>
#include <QMetaObject>
#include <QPushButton>
#include <QVBoxLayout>
#include <filesystem>
#include <memory>
#include <stdexcept>

namespace hornsculpt {
MeshImportDialog::MeshImportDialog(QString path, QWidget *parent)
    : QDialog(parent), path_(std::move(path)) {
    setObjectName("meshImportDialog");
    setWindowTitle(tr("Gmsh-Netz importieren"));
    auto *layout = new QVBoxLayout(this);
    auto *file = new QLabel(QFileInfo(path_).fileName());
    file->setObjectName("meshImportFile");
    file->setTextFormat(Qt::PlainText);
    layout->addWidget(file);
    auto *form = new QFormLayout;
    units_ = new QComboBox;
    units_->setObjectName("meshImportUnits");
    units_->addItem(tr("Einheit wählen"), 0.0);
    units_->addItem(tr("Meter"), 1.0);
    units_->addItem(tr("Millimeter"), .001);
    units_->setToolTip(tr(
        "Einheit der Koordinaten in der Datei. Die Geometrie wird intern in Meter umgerechnet."));
    units_->setWhatsThis(units_->toolTip());
    form->addRow(tr("Dateieinheit"), units_);
    layout->addLayout(form);
    groups_ = new QListWidget;
    groups_->setObjectName("meshImportGroups");
    groups_->setSelectionMode(QAbstractItemView::SingleSelection);
    groups_->setToolTip(tr("Physische Gruppen aus der Datei. Die Auswahl dient zur Ansicht; alle "
                           "Dreiecke werden importiert."));
    groups_->setWhatsThis(groups_->toolTip());
    layout->addWidget(new QLabel(tr("Physische Gruppen")));
    layout->addWidget(groups_);
    status_ = new QLabel;
    status_->setObjectName("meshImportStatus");
    status_->setWordWrap(true);
    status_->setTextFormat(Qt::PlainText);
    layout->addWidget(status_);
    buttons_ = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    buttons_->button(QDialogButtonBox::Ok)->setEnabled(false);
    connect(buttons_, &QDialogButtonBox::accepted, this, &MeshImportDialog::accept);
    connect(buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
    layout->addWidget(buttons_);
    connect(units_, &QComboBox::currentIndexChanged, this, [this] { reload(); });
    reload();
}
const std::optional<MeshImportSelection> &MeshImportDialog::selection() const noexcept {
    return committed_;
}
MeshImportDialog::~MeshImportDialog() {
    worker_.request_stop();
    if (worker_.joinable())
        worker_.join();
}
void MeshImportDialog::reload() {
    ++generation_;
    candidate_.reset();
    groups_->clear();
    buttons_->button(QDialogButtonBox::Ok)->setEnabled(false);
    if (worker_.joinable()) {
        worker_.request_stop();
        status_->setText(tr("Import wird abgebrochen …"));
        return;
    }
    startImport();
}
void MeshImportDialog::startImport() {
    const double scale = units_->currentData().toDouble();
    if (scale <= 0) {
        status_->setText(tr("Dateieinheit auswählen."));
        return;
    }
    status_->setText(tr("Netz wird eingelesen …"));
    const auto generation = generation_;
    const auto sourcePath = path_;
    worker_ = std::jthread([this, sourcePath, scale, generation](std::stop_token stop) {
        auto result = std::make_shared<std::optional<GmshSurfaceImport>>();
        QString error;
        try {
#ifdef _WIN32
            auto path = std::filesystem::path(sourcePath.toStdWString());
#else
            auto path = std::filesystem::path(sourcePath.toStdString());
#endif
            *result = importGmshFile(path, {scale, stop});
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        }
        QMetaObject::invokeMethod(
            this,
            [this, result, error, scale, generation] {
                if (worker_.joinable())
                    worker_.join();
                worker_ = std::jthread{};
                if (generation != generation_) {
                    startImport();
                    return;
                }
                if (!*result) {
                    status_->setText(tr("Import fehlgeschlagen: %1").arg(error));
                    return;
                }
                auto &data = **result;
                for (const auto &group : data.physicalGroups)
                    groups_->addItem(
                        tr("Dimension %1 · ID %2 · %3")
                            .arg(group.dimension)
                            .arg(group.tag)
                            .arg(QString::fromUtf8(group.name.data(),
                                                   static_cast<int>(group.name.size()))));
                const auto &d = data.diagnostics;
                status_->setText(
                    tr("%1 Dreiecke · %2 offene Kanten · %3 nichtmannigfaltige Kanten · "
                       "%4 widersprüchliche Kanten · %5 entartete Dreiecke")
                        .arg(static_cast<qulonglong>(data.mesh.triangles.size()))
                        .arg(static_cast<qulonglong>(d.boundaryEdges))
                        .arg(static_cast<qulonglong>(d.nonManifoldEdges))
                        .arg(static_cast<qulonglong>(d.inconsistentEdges))
                        .arg(static_cast<qulonglong>(d.degenerateTriangles)));
                candidate_ = MeshImportSelection{std::move(data), scale};
                buttons_->button(QDialogButtonBox::Ok)->setEnabled(true);
            },
            Qt::QueuedConnection);
    });
}
void MeshImportDialog::accept() {
    if (!candidate_)
        return;
    committed_ = std::move(candidate_);
    QDialog::accept();
}
} // namespace hornsculpt
