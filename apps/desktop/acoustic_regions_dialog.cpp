#include "acoustic_regions_dialog.hpp"
#include <QColor>
#include <QDialogButtonBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QRegularExpression>
#include <QStringList>
#include <QTableWidget>
#include <QVBoxLayout>
#include <algorithm>
#include <limits>
#include <set>

namespace hornsculpt {
namespace {
QString ids(const std::vector<std::uint32_t> &values) {
    QStringList parts;
    for (auto id : values)
        parts << QString::number(id);
    return parts.join(", ");
}
bool parseIds(const QString &text, const std::map<std::uint32_t, QString> &available,
              std::vector<std::uint32_t> &result) {
    if (text.trimmed().isEmpty())
        return true;
    std::set<std::uint32_t> seen;
    for (const auto &part : text.split(',')) {
        const auto token = part.trimmed();
        static const QRegularExpression pattern(QStringLiteral("^[1-9][0-9]*$"));
        bool ok = false;
        const auto value = token.toULongLong(&ok);
        if (!pattern.match(token).hasMatch() || !ok ||
            value > std::numeric_limits<std::uint32_t>::max())
            return false;
        auto id = static_cast<std::uint32_t>(value);
        if (!available.contains(id) || !seen.insert(id).second)
            return false;
        result.push_back(id);
    }
    return true;
}
} // namespace
AcousticRegionsDialog::AcousticRegionsDialog(std::vector<AcousticRegionDefinition> regions,
                                             std::map<std::uint32_t, QString> availableFaceLabels,
                                             QWidget *parent)
    : QDialog(parent), availableFaces_(std::move(availableFaceLabels)) {
    setWindowTitle(tr("Akustische Gebiete"));
    setObjectName("acousticRegionsDialog");
    auto *layout = new QVBoxLayout(this);
    table_ = new QTableWidget(0, 3);
    table_->setObjectName("acousticRegionsTable");
    table_->setHorizontalHeaderLabels({tr("Name"), tr("Wand und Antrieb"), tr("Interface")});
    table_->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    table_->setSelectionBehavior(QAbstractItemView::SelectItems);
    table_->setToolTip(
        tr("Flächen-IDs durch Komma trennen. Die Interface-Normale zeigt in das erste Gebiet."));
    table_->setWhatsThis(table_->toolTip());
    layout->addWidget(table_);
    for (const auto &region : regions) {
        const int row = table_->rowCount();
        table_->insertRow(row);
        auto *name = new QTableWidgetItem(QString::fromStdString(region.name));
        name->setData(Qt::UserRole, QVariant::fromValue<qulonglong>(region.id));
        if (region.needsAssignmentReview) {
            name->setForeground(QColor(Qt::red));
            name->setToolTip(tr("Flächenzuordnung prüfen"));
        }
        table_->setItem(row, 0, name);
        table_->setItem(row, 1, new QTableWidgetItem(ids(region.faceIds)));
        table_->setItem(row, 2, new QTableWidgetItem(ids(region.interfaceFaceIds)));
        nextId_ = std::max(nextId_, static_cast<std::uint64_t>(region.id) + 1);
    }
    auto *rowButtons = new QHBoxLayout;
    auto *add = new QPushButton(tr("Gebiet hinzufügen"));
    add->setObjectName("acousticRegionAdd");
    auto *remove = new QPushButton(tr("Gebiet entfernen"));
    remove->setObjectName("acousticRegionRemove");
    connect(add, &QPushButton::clicked, this, [this] { addRegion(); });
    connect(remove, &QPushButton::clicked, this, [this] { removeRegion(); });
    rowButtons->addWidget(add);
    rowButtons->addWidget(remove);
    rowButtons->addStretch();
    layout->addLayout(rowButtons);
    layout->addWidget(new QLabel(tr("Verfügbare Flächen")));
    faces_ = new QListWidget;
    faces_->setObjectName("acousticRegionFaces");
    faces_->setSelectionMode(QAbstractItemView::ExtendedSelection);
    for (const auto &[id, label] : availableFaces_) {
        auto *item = new QListWidgetItem(QStringLiteral("%1 · %2").arg(id).arg(label), faces_);
        item->setData(Qt::UserRole, QVariant::fromValue<qulonglong>(id));
    }
    layout->addWidget(faces_);
    auto *faceButtons = new QHBoxLayout;
    auto *physical = new QPushButton(tr("Als Wand / Antrieb"));
    physical->setObjectName("acousticRegionUsePhysical");
    auto *interface = new QPushButton(tr("Als Interface"));
    interface->setObjectName("acousticRegionUseInterface");
    connect(physical, &QPushButton::clicked, this, [this] { useSelectedFaces(1); });
    connect(interface, &QPushButton::clicked, this, [this] { useSelectedFaces(2); });
    faceButtons->addWidget(physical);
    faceButtons->addWidget(interface);
    faceButtons->addStretch();
    layout->addLayout(faceButtons);
    error_ = new QLabel;
    error_->setObjectName("acousticRegionsError");
    error_->setWordWrap(true);
    error_->setTextFormat(Qt::PlainText);
    layout->addWidget(error_);
    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(buttons, &QDialogButtonBox::accepted, this, &AcousticRegionsDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
    layout->addWidget(buttons);
}
const std::vector<AcousticRegionDefinition> &AcousticRegionsDialog::regions() const noexcept {
    return committed_;
}
void AcousticRegionsDialog::addRegion() {
    if (nextId_ > std::numeric_limits<std::uint32_t>::max()) {
        error_->setText(tr("Keine freie Gebiets-ID."));
        return;
    }
    const int row = table_->rowCount();
    table_->insertRow(row);
    auto *name = new QTableWidgetItem;
    name->setData(Qt::UserRole, QVariant::fromValue<qulonglong>(nextId_++));
    table_->setItem(row, 0, name);
    table_->setItem(row, 1, new QTableWidgetItem);
    table_->setItem(row, 2, new QTableWidgetItem);
    table_->setCurrentCell(row, 0);
    table_->editItem(name);
}
void AcousticRegionsDialog::removeRegion() {
    if (table_->currentRow() >= 0)
        table_->removeRow(table_->currentRow());
}
void AcousticRegionsDialog::useSelectedFaces(int column) {
    const int row = table_->currentRow();
    if (row < 0) {
        error_->setText(tr("Zuerst ein Gebiet auswählen."));
        return;
    }
    std::vector<std::uint32_t> selected;
    for (const auto *item : faces_->selectedItems())
        selected.push_back(static_cast<std::uint32_t>(item->data(Qt::UserRole).toULongLong()));
    std::sort(selected.begin(), selected.end());
    table_->item(row, column)->setText(ids(selected));
    table_->setCurrentCell(row, column);
    error_->clear();
}
void AcousticRegionsDialog::accept() {
    std::vector<AcousticRegionDefinition> candidate;
    std::set<std::uint32_t> idsSeen, physical;
    std::map<std::uint32_t, unsigned> interfaces;
    for (int row = 0; row < table_->rowCount(); ++row) {
        const auto *name = table_->item(row, 0);
        const auto *p = table_->item(row, 1);
        const auto *i = table_->item(row, 2);
        bool ok = false;
        const auto raw = name ? name->data(Qt::UserRole).toULongLong(&ok) : 0;
        if (!name || !p || !i || !ok || raw == 0 ||
            raw > std::numeric_limits<std::uint32_t>::max() || name->text().trimmed().isEmpty() ||
            !idsSeen.insert(static_cast<std::uint32_t>(raw)).second) {
            error_->setText(tr("Jedes Gebiet benötigt einen Namen und eine eindeutige ID."));
            return;
        }
        AcousticRegionDefinition region;
        region.id = static_cast<std::uint32_t>(raw);
        region.name = name->text().trimmed().toStdString();
        if (!parseIds(p->text(), availableFaces_, region.faceIds) ||
            !parseIds(i->text(), availableFaces_, region.interfaceFaceIds) ||
            (region.faceIds.empty() && region.interfaceFaceIds.empty())) {
            error_->setText(tr("Flächen-IDs prüfen: bekannt, positiv und je Liste eindeutig."));
            return;
        }
        for (auto id : region.faceIds)
            if (!physical.insert(id).second) {
                error_->setText(tr("Eine physische Fläche darf nur einem Gebiet gehören."));
                return;
            }
        for (auto id : region.interfaceFaceIds)
            if (++interfaces[id] > 2) {
                error_->setText(tr("Ein Interface darf höchstens zwei Gebiete verbinden."));
                return;
            }
        candidate.push_back(std::move(region));
    }
    for (const auto &[id, count] : interfaces) {
        (void)count;
        if (physical.contains(id)) {
            error_->setText(tr("Eine Fläche darf nicht zugleich physisch und Interface sein."));
            return;
        }
    }
    committed_ = std::move(candidate);
    error_->clear();
    QDialog::accept();
}
} // namespace hornsculpt
