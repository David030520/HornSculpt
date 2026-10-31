#include "acoustic_regions_dialog.hpp"
#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QInputDialog>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QPushButton>
#include <QTreeWidget>
#include <QVBoxLayout>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace hornsculpt {
namespace {
QString faceText(const std::vector<std::uint32_t> &faces) {
    QStringList result;
    for (auto id : faces)
        result << QString::number(id);
    return result.join(", ");
}
template <class T> std::uint32_t nextId(const std::vector<T> &items) {
    std::set<std::uint32_t> ids;
    for (const auto &item : items)
        ids.insert(item.id);
    std::uint32_t id = 1;
    while (ids.contains(id)) {
        if (id == std::numeric_limits<std::uint32_t>::max())
            return 0;
        ++id;
    }
    return id;
}
} // namespace
AcousticRegionsDialog::AcousticRegionsDialog(AcousticModelDefinition model,
                                             std::map<std::uint32_t, QString> availableFaceLabels,
                                             QWidget *parent)
    : QDialog(parent), working_(std::move(model)), availableFaces_(std::move(availableFaceLabels)) {
    setWindowTitle(tr("Akustische Gebiete"));
    setObjectName("acousticRegionsDialog");
    auto *layout = new QVBoxLayout(this);
    tree_ = new QTreeWidget;
    tree_->setObjectName("acousticModelTree");
    tree_->setHeaderLabels({tr("Gebiet / Interface"), tr("Zuordnung")});
    tree_->header()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    tree_->header()->setSectionResizeMode(1, QHeaderView::Stretch);
    connect(tree_, &QTreeWidget::itemDoubleClicked, this, [this](QTreeWidgetItem *item) {
        if (item->data(0, Qt::UserRole).toInt() == 5)
            editBaffle();
        else
            editCurrent();
    });
    layout->addWidget(tree_);
    auto *actions = new QHBoxLayout;
    auto addButton = [&](const QString &name, const char *object, auto fn) {
        auto *button = new QPushButton(name);
        button->setObjectName(object);
        connect(button, &QPushButton::clicked, this, fn);
        actions->addWidget(button);
    };
    addButton(tr("Innenraum"), "acousticAddInterior",
              [this] { addDomain(BemDomainKind::interior); });
    addButton(tr("Außenraum"), "acousticAddExterior",
              [this] { addDomain(BemDomainKind::exterior); });
    addButton(tr("Interface"), "acousticAddInterface", [this] { addInterface(); });
    addButton(tr("Gruppe"), "acousticAddGroup", [this] { addGroup(); });
    addButton(tr("Bearbeiten"), "acousticEdit", [this] { editCurrent(); });
    addButton(tr("Entfernen"), "acousticRemove", [this] { removeCurrent(); });
    addButton(tr("Zuordnung bestätigt"), "acousticConfirmAssignment", [this] {
        auto *item = tree_->currentItem();
        if (!item)
            return;
        if (item->parent())
            item = item->parent();
        const int type = item->data(0, Qt::UserRole).toInt();
        const auto id = item->data(0, Qt::UserRole + 1).toUInt();
        if (type == 1)
            for (auto &d : working_.domains)
                if (d.id == id)
                    d.needsAssignmentReview = false;
        if (type == 2)
            for (auto &i : working_.interfaces)
                if (i.id == id)
                    i.needsAssignmentReview = false;
        refresh();
    });
    addButton(tr("Schallwand"), "acousticBaffle", [this] { editBaffle(); });
    layout->addLayout(actions);
    connect(tree_, &QTreeWidget::itemClicked, this, [this](QTreeWidgetItem *item) {
        const int type = item->data(0, Qt::UserRole).toInt();
        if (type != 3 && type != 4)
            return;
        const auto groupId = item->data(0, Qt::UserRole + 1).toUInt();
        const auto owner = item->parent()->data(0, Qt::UserRole + 1).toUInt();
        faces_->clearSelection();
        if (type == 3)
            for (const auto &d : working_.domains)
                if (d.id == owner)
                    for (const auto &g : d.groups)
                        if (g.id == groupId)
                            selectFaces(g.faceIds);
        if (type == 4)
            for (const auto &i : working_.interfaces)
                if (i.id == owner)
                    for (const auto &g : i.groups)
                        if (g.id == groupId)
                            selectFaces(g.faceIds);
    });
    layout->addWidget(new QLabel(tr("Verfügbare Flächen")));
    faces_ = new QListWidget;
    faces_->setObjectName("acousticRegionFaces");
    faces_->setSelectionMode(QAbstractItemView::ExtendedSelection);
    for (const auto &[id, label] : availableFaces_) {
        auto *item = new QListWidgetItem(QStringLiteral("%1 · %2").arg(id).arg(label), faces_);
        item->setData(Qt::UserRole, QVariant::fromValue<qulonglong>(id));
    }
    layout->addWidget(faces_);
    auto *assignment = new QHBoxLayout;
    auto *assign = new QPushButton(tr("Auswahl zuordnen"));
    assign->setObjectName("acousticAssignFaces");
    connect(assign, &QPushButton::clicked, this, [this] { assignFaces(false); });
    auto *both = new QPushButton(tr("Beidseitig zuordnen"));
    both->setObjectName("acousticAssignTwoSided");
    connect(both, &QPushButton::clicked, this, [this] { assignFaces(true); });
    assignment->addWidget(assign);
    assignment->addWidget(both);
    layout->addLayout(assignment);
    error_ = new QLabel;
    error_->setObjectName("acousticRegionsError");
    error_->setWordWrap(true);
    layout->addWidget(error_);
    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(buttons, &QDialogButtonBox::accepted, this, &AcousticRegionsDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
    layout->addWidget(buttons);
    refresh();
}
const AcousticModelDefinition &AcousticRegionsDialog::model() const noexcept {
    return committed_;
}
void AcousticRegionsDialog::selectFaces(const std::vector<std::uint32_t> &ids) {
    for (int i = 0; i < faces_->count(); ++i) {
        auto *item = faces_->item(i);
        item->setSelected(std::find(ids.begin(), ids.end(), item->data(Qt::UserRole).toUInt()) !=
                          ids.end());
    }
}
void AcousticRegionsDialog::refresh() {
    tree_->clear();
    for (const auto &domain : working_.domains) {
        auto *item = new QTreeWidgetItem(
            tree_, {QString::fromStdString(domain.name),
                    domain.kind == BemDomainKind::interior ? tr("Innenraum") : tr("Außenraum")});
        item->setToolTip(0, QString::fromStdString(domain.name));
        item->setData(0, Qt::UserRole, 1);
        item->setData(0, Qt::UserRole + 1, domain.id);
        new QTreeWidgetItem(item, {tr("Wall / Driven"), faceText(domain.faceIds)});
        new QTreeWidgetItem(item, {tr("Beidseitig"), faceText(domain.twoSidedFaceIds)});
        for (const auto &group : domain.groups) {
            auto *child = new QTreeWidgetItem(
                item, {QString::fromStdString(group.name), faceText(group.faceIds)});
            child->setToolTip(0, QString::fromStdString(group.name));
            child->setData(0, Qt::UserRole, 3);
            child->setData(0, Qt::UserRole + 1, group.id);
        }
        if (domain.kind == BemDomainKind::exterior && working_.exteriorBaffle) {
            const auto &b = *working_.exteriorBaffle;
            auto *plane = new QTreeWidgetItem(
                item, {tr("Unendliche Schallwand"),
                       tr("Ursprung [%1, %2, %3] mm · Normale [%4, %5, %6] · Toleranz %7 mm")
                           .arg(b.origin.x * 1000)
                           .arg(b.origin.y * 1000)
                           .arg(b.origin.z * 1000)
                           .arg(b.normal.x)
                           .arg(b.normal.y)
                           .arg(b.normal.z)
                           .arg(b.planeDistanceTolerance * 1000)});
            plane->setData(0, Qt::UserRole, 5);
            plane->setToolTip(
                0, tr("Doppelklick bearbeitet die Ebene; die Normale zeigt in den Außenhalbraum."));
        }
        if (domain.needsAssignmentReview)
            item->setText(1, item->text(1) + tr(" · prüfen"));
    }
    for (const auto &interface : working_.interfaces) {
        const auto name = [this](std::uint32_t id) {
            const auto it = std::find_if(working_.domains.begin(), working_.domains.end(),
                                         [id](const auto &d) { return d.id == id; });
            return it == working_.domains.end() ? QString::number(id)
                                                : QString::fromStdString(it->name);
        };
        auto *item = new QTreeWidgetItem(
            tree_, {QString::fromStdString(interface.name),
                    tr("Subdomain 1: %1 · Subdomain 2: %2")
                        .arg(name(interface.firstDomainId), name(interface.secondDomainId))});
        item->setToolTip(0, QString::fromStdString(interface.name) + QStringLiteral("\n") +
                                tr("Flächennormale zeigt in Subdomain 1; die zweite Seite "
                                   "verwendet die Gegenrichtung."));
        item->setData(0, Qt::UserRole, 2);
        item->setData(0, Qt::UserRole + 1, interface.id);
        new QTreeWidgetItem(item, {tr("Flächen"), faceText(interface.faceIds)});
        for (const auto &group : interface.groups) {
            auto *child = new QTreeWidgetItem(
                item, {QString::fromStdString(group.name), faceText(group.faceIds)});
            child->setToolTip(0, QString::fromStdString(group.name));
            child->setData(0, Qt::UserRole, 4);
            child->setData(0, Qt::UserRole + 1, group.id);
        }
        if (interface.needsAssignmentReview)
            item->setText(1, item->text(1) + tr(" · prüfen"));
    }
    tree_->expandAll();
}
void AcousticRegionsDialog::addDomain(BemDomainKind kind) {
    if (kind == BemDomainKind::exterior &&
        std::any_of(working_.domains.begin(), working_.domains.end(),
                    [](const auto &d) { return d.kind == BemDomainKind::exterior; })) {
        error_->setText(tr("Nur ein Außenraum ist möglich."));
        return;
    }
    bool ok = false;
    const auto name =
        QInputDialog::getText(this, tr("Gebiet"), tr("Name"), QLineEdit::Normal, {}, &ok).trimmed();
    if (!ok || name.isEmpty())
        return;
    const auto id = nextId(working_.domains);
    if (!id) {
        error_->setText(tr("Keine freie Gebiets-ID."));
        return;
    }
    working_.domains.push_back({id, name.toStdString(), kind});
    refresh();
}
void AcousticRegionsDialog::addInterface() {
    if (working_.domains.size() < 2) {
        error_->setText(tr("Zwei Gebiete anlegen."));
        return;
    }
    bool ok = false;
    const auto name =
        QInputDialog::getText(this, tr("Interface"), tr("Name"), QLineEdit::Normal, {}, &ok)
            .trimmed();
    if (!ok || name.isEmpty())
        return;
    QStringList names;
    for (const auto &d : working_.domains)
        names << QStringLiteral("%1 · %2").arg(d.id).arg(QString::fromStdString(d.name));
    auto first =
        QInputDialog::getItem(this, tr("Interface"), tr("Subdomain 1"), names, 0, false, &ok);
    if (!ok)
        return;
    auto second =
        QInputDialog::getItem(this, tr("Interface"), tr("Subdomain 2"), names, 1, false, &ok);
    if (!ok)
        return;
    const auto a = names.indexOf(first), b = names.indexOf(second);
    if (a == b) {
        error_->setText(tr("Interface benötigt zwei verschiedene Gebiete."));
        return;
    }
    const auto id = nextId(working_.interfaces);
    if (!id) {
        error_->setText(tr("Keine freie Interface-ID."));
        return;
    }
    working_.interfaces.push_back(
        {id, name.toStdString(), working_.domains[a].id, working_.domains[b].id, {}});
    refresh();
}
void AcousticRegionsDialog::addGroup() {
    auto *item = tree_->currentItem();
    if (!item) {
        error_->setText(tr("Gebiet oder Interface auswählen."));
        return;
    }
    if (item->parent())
        item = item->parent();
    const int type = item->data(0, Qt::UserRole).toInt();
    const auto owner = item->data(0, Qt::UserRole + 1).toUInt();
    std::vector<std::uint32_t> selected;
    for (auto *f : faces_->selectedItems())
        selected.push_back(f->data(Qt::UserRole).toUInt());
    std::sort(selected.begin(), selected.end());
    if (selected.empty()) {
        error_->setText(tr("Flächen für die Gruppe auswählen."));
        return;
    }
    std::vector<AcousticElementGroup> *groups = nullptr;
    std::vector<std::uint32_t> allowed;
    if (type == 1)
        for (auto &d : working_.domains)
            if (d.id == owner) {
                groups = &d.groups;
                allowed = d.faceIds;
                allowed.insert(allowed.end(), d.twoSidedFaceIds.begin(), d.twoSidedFaceIds.end());
            }
    if (type == 2)
        for (auto &i : working_.interfaces)
            if (i.id == owner) {
                groups = &i.groups;
                allowed = i.faceIds;
            }
    if (!groups)
        return;
    for (auto id : selected) {
        if (std::find(allowed.begin(), allowed.end(), id) == allowed.end() ||
            std::any_of(groups->begin(), groups->end(), [id](const auto &g) {
                return std::find(g.faceIds.begin(), g.faceIds.end(), id) != g.faceIds.end();
            })) {
            error_->setText(
                tr("Gruppenflächen müssen zugeordnet und innerhalb des Eintrags eindeutig sein."));
            return;
        }
    }
    QString suggestion = selected.size() == 1 ? availableFaces_[selected.front()] : QString{};
    bool ok = false;
    const auto name = QInputDialog::getText(this, tr("Elementgruppe"), tr("Name"),
                                            QLineEdit::Normal, suggestion, &ok)
                          .trimmed();
    if (!ok || name.isEmpty())
        return;
    const auto id = nextId(*groups);
    if (!id) {
        error_->setText(tr("Keine freie Gruppen-ID."));
        return;
    }
    groups->push_back({id, name.toStdString(), std::move(selected)});
    refresh();
}
void AcousticRegionsDialog::editCurrent() {
    auto *item = tree_->currentItem();
    if (!item)
        return;
    const int childType = item->data(0, Qt::UserRole).toInt();
    if ((childType == 3 || childType == 4) && item->parent()) {
        const auto owner = item->parent()->data(0, Qt::UserRole + 1).toUInt();
        const auto id = item->data(0, Qt::UserRole + 1).toUInt();
        bool ok = false;
        const auto name = QInputDialog::getText(this, tr("Elementgruppe"), tr("Name"),
                                                QLineEdit::Normal, item->text(0), &ok)
                              .trimmed();
        if (!ok || name.isEmpty())
            return;
        if (childType == 3)
            for (auto &d : working_.domains)
                if (d.id == owner)
                    for (auto &g : d.groups)
                        if (g.id == id)
                            g.name = name.toStdString();
        if (childType == 4)
            for (auto &i : working_.interfaces)
                if (i.id == owner)
                    for (auto &g : i.groups)
                        if (g.id == id)
                            g.name = name.toStdString();
        refresh();
        return;
    }
    if (item->parent())
        item = item->parent();
    const int type = item->data(0, Qt::UserRole).toInt();
    const auto id = item->data(0, Qt::UserRole + 1).toUInt();
    bool ok = false;
    const auto name = QInputDialog::getText(this, tr("Bearbeiten"), tr("Name"), QLineEdit::Normal,
                                            item->text(0), &ok)
                          .trimmed();
    if (!ok || name.isEmpty())
        return;
    if (type == 1)
        for (auto &d : working_.domains)
            if (d.id == id)
                d.name = name.toStdString();
    if (type == 2)
        for (auto &i : working_.interfaces)
            if (i.id == id) {
                QStringList names;
                for (const auto &d : working_.domains)
                    names << QStringLiteral("%1 · %2").arg(d.id).arg(
                        QString::fromStdString(d.name));
                int a = 0, b = 0;
                for (int n = 0; n < static_cast<int>(working_.domains.size()); ++n) {
                    if (working_.domains[n].id == i.firstDomainId)
                        a = n;
                    if (working_.domains[n].id == i.secondDomainId)
                        b = n;
                }
                const auto first = QInputDialog::getItem(this, tr("Interface"), tr("Subdomain 1"),
                                                         names, a, false, &ok);
                if (!ok)
                    return;
                const auto second = QInputDialog::getItem(this, tr("Interface"), tr("Subdomain 2"),
                                                          names, b, false, &ok);
                if (!ok)
                    return;
                a = names.indexOf(first);
                b = names.indexOf(second);
                if (a == b) {
                    error_->setText(tr("Zwei verschiedene Gebiete wählen."));
                    return;
                }
                i.name = name.toStdString();
                i.firstDomainId = working_.domains[a].id;
                i.secondDomainId = working_.domains[b].id;
            }
    refresh();
}
void AcousticRegionsDialog::removeCurrent() {
    auto *item = tree_->currentItem();
    if (!item)
        return;
    const int childType = item->data(0, Qt::UserRole).toInt();
    if ((childType == 3 || childType == 4) && item->parent()) {
        const auto owner = item->parent()->data(0, Qt::UserRole + 1).toUInt();
        const auto id = item->data(0, Qt::UserRole + 1).toUInt();
        if (childType == 3)
            for (auto &d : working_.domains)
                if (d.id == owner)
                    std::erase_if(d.groups, [id](const auto &g) { return g.id == id; });
        if (childType == 4)
            for (auto &i : working_.interfaces)
                if (i.id == owner)
                    std::erase_if(i.groups, [id](const auto &g) { return g.id == id; });
        refresh();
        return;
    }
    if (item->parent())
        item = item->parent();
    const int type = item->data(0, Qt::UserRole).toInt();
    const auto id = item->data(0, Qt::UserRole + 1).toUInt();
    if (type == 1) {
        if (std::any_of(
                working_.interfaces.begin(), working_.interfaces.end(),
                [id](const auto &i) { return i.firstDomainId == id || i.secondDomainId == id; })) {
            error_->setText(tr("Zuerst zugehörige Interfaces entfernen."));
            return;
        }
        std::erase_if(working_.domains, [id](const auto &d) { return d.id == id; });
    } else if (type == 2)
        std::erase_if(working_.interfaces, [id](const auto &i) { return i.id == id; });
    refresh();
}
void AcousticRegionsDialog::assignFaces(bool twoSided) {
    auto *item = tree_->currentItem();
    if (!item) {
        error_->setText(tr("Gebiet oder Interface auswählen."));
        return;
    }
    std::vector<std::uint32_t> faces;
    for (auto *face : faces_->selectedItems())
        faces.push_back(face->data(Qt::UserRole).toUInt());
    std::sort(faces.begin(), faces.end());
    const int childType = item->data(0, Qt::UserRole).toInt();
    if ((childType == 3 || childType == 4) && item->parent()) {
        const auto groupId = item->data(0, Qt::UserRole + 1).toUInt();
        const auto owner = item->parent()->data(0, Qt::UserRole + 1).toUInt();
        std::vector<AcousticElementGroup> *groups = nullptr;
        std::vector<std::uint32_t> allowed;
        if (childType == 3)
            for (auto &d : working_.domains)
                if (d.id == owner) {
                    groups = &d.groups;
                    allowed = d.faceIds;
                    allowed.insert(allowed.end(), d.twoSidedFaceIds.begin(),
                                   d.twoSidedFaceIds.end());
                }
        if (childType == 4)
            for (auto &i : working_.interfaces)
                if (i.id == owner) {
                    groups = &i.groups;
                    allowed = i.faceIds;
                }
        if (!groups || faces.empty()) {
            error_->setText(tr("Gruppe benötigt zugeordnete Flächen."));
            return;
        }
        for (auto face : faces)
            if (std::find(allowed.begin(), allowed.end(), face) == allowed.end() ||
                std::any_of(groups->begin(), groups->end(), [&](const auto &g) {
                    return g.id != groupId &&
                           std::find(g.faceIds.begin(), g.faceIds.end(), face) != g.faceIds.end();
                })) {
                error_->setText(tr("Gruppenflächen müssen zugeordnet und eindeutig sein."));
                return;
            }
        for (auto &g : *groups)
            if (g.id == groupId)
                g.faceIds = faces;
        refresh();
        return;
    }
    if (item->parent())
        item = item->parent();
    const int type = item->data(0, Qt::UserRole).toInt();
    const auto id = item->data(0, Qt::UserRole + 1).toUInt();
    if (type == 1)
        for (auto &d : working_.domains)
            if (d.id == id) {
                (twoSided ? d.twoSidedFaceIds : d.faceIds) = faces;
            }
    if (type == 2 && !twoSided)
        for (auto &i : working_.interfaces)
            if (i.id == id) {
                i.faceIds = faces;
            }
    if (type == 2 && twoSided) {
        error_->setText(tr("Interface ist keine beidseitige Wand."));
        return;
    }
    auto prune = [](auto &groups, const auto &allowed) {
        for (auto &g : groups)
            std::erase_if(g.faceIds, [&](auto id) {
                return std::find(allowed.begin(), allowed.end(), id) == allowed.end();
            });
        std::erase_if(groups, [](const auto &g) { return g.faceIds.empty(); });
    };
    if (type == 1)
        for (auto &d : working_.domains)
            if (d.id == id) {
                auto allowed = d.faceIds;
                allowed.insert(allowed.end(), d.twoSidedFaceIds.begin(), d.twoSidedFaceIds.end());
                prune(d.groups, allowed);
            }
    if (type == 2)
        for (auto &i : working_.interfaces)
            if (i.id == id)
                prune(i.groups, i.faceIds);
    refresh();
}
void AcousticRegionsDialog::editBaffle() {
    QDialog dialog(this);
    dialog.setWindowTitle(tr("Außenraum-Schallwand"));
    auto *form = new QFormLayout(&dialog);
    QCheckBox enabled;
    enabled.setChecked(working_.exteriorBaffle.has_value());
    form->addRow(tr("Aktiv"), &enabled);
    const auto value = working_.exteriorBaffle.value_or(ExteriorBaffle{{}, {0, 0, 1}, 1e-6});
    QDoubleSpinBox ox, oy, oz, nx, ny, nz, tolerance;
    auto setup = [](QDoubleSpinBox &box, double v, double range) {
        box.setDecimals(9);
        box.setRange(-range, range);
        box.setValue(v);
    };
    setup(ox, value.origin.x * 1000, 1e9);
    setup(oy, value.origin.y * 1000, 1e9);
    setup(oz, value.origin.z * 1000, 1e9);
    setup(nx, value.normal.x, 1);
    setup(ny, value.normal.y, 1);
    setup(nz, value.normal.z, 1);
    tolerance.setDecimals(6);
    tolerance.setRange(.000001, 10);
    tolerance.setValue(value.planeDistanceTolerance * 1000);
    form->addRow(tr("Ursprung X [mm]"), &ox);
    form->addRow(tr("Ursprung Y [mm]"), &oy);
    form->addRow(tr("Ursprung Z [mm]"), &oz);
    form->addRow(tr("Normale X"), &nx);
    form->addRow(tr("Normale Y"), &ny);
    form->addRow(tr("Normale Z"), &nz);
    form->addRow(tr("Toleranz [mm]"), &tolerance);
    QDialogButtonBox buttons(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    form->addRow(&buttons);
    connect(&buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
    connect(&buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    if (dialog.exec() != QDialog::Accepted)
        return;
    if (!enabled.isChecked()) {
        working_.exteriorBaffle.reset();
        refresh();
        return;
    }
    const double length = std::hypot(nx.value(), ny.value(), nz.value());
    if (std::abs(length - 1) > 1e-8) {
        error_->setText(tr("Schallwandnormale muss Länge 1 haben."));
        return;
    }
    working_.exteriorBaffle =
        ExteriorBaffle{{ox.value() * .001, oy.value() * .001, oz.value() * .001},
                       {nx.value(), ny.value(), nz.value()},
                       tolerance.value() * .001};
    refresh();
}
void AcousticRegionsDialog::accept() {
    std::set<std::uint32_t> domains, physical, interfaces;
    std::map<std::uint32_t, unsigned> both;
    unsigned exterior = 0;
    for (const auto &d : working_.domains) {
        if (!d.id || d.name.empty() || !domains.insert(d.id).second) {
            error_->setText(tr("Gebietsname und ID prüfen."));
            return;
        }
        exterior += d.kind == BemDomainKind::exterior;
        for (auto id : d.faceIds)
            if (!availableFaces_.contains(id) || !physical.insert(id).second) {
                error_->setText(tr("Flächenzuordnung prüfen."));
                return;
            }
        for (auto id : d.twoSidedFaceIds)
            if (!availableFaces_.contains(id) || ++both[id] > 2) {
                error_->setText(tr("Beidseitige Flächen prüfen."));
                return;
            }
    }
    if ((!working_.domains.empty() && exterior != 1) ||
        (working_.exteriorBaffle && exterior != 1)) {
        error_->setText(tr("Genau ein Außenraum erforderlich."));
        return;
    }
    std::set<std::uint32_t> ids;
    for (const auto &i : working_.interfaces) {
        if (!i.id || i.name.empty() || !ids.insert(i.id).second ||
            !domains.contains(i.firstDomainId) || !domains.contains(i.secondDomainId) ||
            i.firstDomainId == i.secondDomainId || i.faceIds.empty()) {
            error_->setText(tr("Interface und Gebietsreferenzen prüfen."));
            return;
        }
        for (auto id : i.faceIds)
            if (!availableFaces_.contains(id) || !interfaces.insert(id).second ||
                physical.contains(id) || both.contains(id)) {
                error_->setText(tr("Interfaceflächen müssen eindeutig sein."));
                return;
            }
    }
    for (const auto &[id, count] : both) {
        (void)count;
        if (physical.contains(id)) {
            error_->setText(tr("Beidseitige Fläche ist bereits einseitig zugeordnet."));
            return;
        }
    }
    committed_ = working_;
    error_->clear();
    QDialog::accept();
}
} // namespace hornsculpt
