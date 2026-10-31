#include "draft_assembly_dialog.hpp"
#include "numeric_input.hpp"
#include <QCheckBox>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QHeaderView>
#include <QMessageBox>
#include <QPushButton>
#include <QTableWidget>
#include <QVBoxLayout>
#include <algorithm>
#include <cmath>
#include <numbers>

namespace hornsculpt {
namespace {
constexpr double degrees = 180 / std::numbers::pi;
}
DraftAssemblyDialog::DraftAssemblyDialog(DraftAssembly a, bool hasCad, QWidget *parent)
    : QDialog(parent) {
    setWindowTitle(tr("Draft-Baugruppe · gekoppelte Hörner"));
    setObjectName("draftAssemblyDialog");
    resize(1050, 420);
    auto *layout = new QVBoxLayout(this);
    table_ = new QTableWidget(0, 10);
    table_->setHorizontalHeaderLabels({tr("Horn"), tr("X [mm]"), tr("Y [mm]"), tr("Z [mm]"),
                                       tr("Rz [°]"), tr("Ry [°]"), tr("Rx [°]"), tr("Spiegel X"),
                                       tr("Anker Mund"), tr("ID")});
    table_->setColumnHidden(9, true);
    table_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    table_->setToolTip(
        tr("Alle Instanzen teilen die Draft-Parameter. Drehfolge: zuerst X, dann Y, dann Z. "
           "Der Anker liegt am Hals oder am Mund. Spiegelung erfolgt im lokalen X."));
    table_->setWhatsThis(table_->toolTip());
    layout->addWidget(table_);
    for (const auto &i : a.instances) {
        addRow(i);
        nextId_ = std::max(nextId_, i.id + 1);
    }
    auto *buttons = new QHBoxLayout;
    auto *add = new QPushButton(tr("Gekoppelte Kopie"));
    auto *remove = new QPushButton(tr("Entfernen"));
    buttons->addWidget(add);
    buttons->addWidget(remove);
    buttons->addStretch();
    layout->addLayout(buttons);
    connect(add, &QPushButton::clicked, this, [this] {
        if (table_->rowCount() >= 64)
            return;
        auto a = assembly();
        auto i = a.instances.at(std::max(0, table_->currentRow()));
        i.id = nextId_++;
        i.name = "Horn " + std::to_string(i.id);
        i.position.x += .5;
        addRow(i);
    });
    connect(remove, &QPushButton::clicked, this, [this] {
        if (table_->rowCount() > 1 && table_->currentRow() >= 0)
            table_->removeRow(table_->currentRow());
    });
    auto *form = new QFormLayout;
    includeCad_ = new QCheckBox(tr("CAD-Körper gemeinsam anzeigen und simulieren"));
    includeCad_->setChecked(a.includeCad);
    includeCad_->setEnabled(hasCad);
    form->addRow(includeCad_);
    cadSize_ = new DecimalInput;
    cadSize_->setRange(.05, 2000);
    cadSize_->setDecimals(3);
    cadSize_->setSuffix(" mm");
    cadSize_->setValue(a.cadTargetEdgeLength * 1000);
    cadSize_->setToolTip(
        tr("Eigenes skalares Größenfeld für CAD-Körper. Horn-Netze verwenden die "
           "Draft-Netzeinstellungen. Große Kanten begrenzen die akustische Genauigkeit."));
    cadSize_->setWhatsThis(cadSize_->toolTip());
    form->addRow(tr("CAD-Zielkante"), cadSize_);
    layout->addLayout(form);
    auto *ok = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(ok, &QDialogButtonBox::rejected, this, &QDialog::reject);
    connect(ok, &QDialogButtonBox::accepted, this, [this] {
        try {
            validateDraftAssembly(assembly());
            accept();
        } catch (const std::exception &e) {
            QMessageBox::warning(this, tr("Baugruppe"), QString::fromUtf8(e.what()));
        }
    });
    layout->addWidget(ok);
}
void DraftAssemblyDialog::addRow(const DraftInstance &i) {
    const int row = table_->rowCount();
    table_->insertRow(row);
    table_->setItem(row, 0, new QTableWidgetItem(QString::fromStdString(i.name)));
    table_->setItem(row, 9, new QTableWidgetItem(QString::number(i.id)));
    const double pitch = std::asin(std::clamp(-i.axes[0].z, -1., 1.));
    const double yaw = std::abs(std::cos(pitch)) > 1e-9 ? std::atan2(i.axes[0].y, i.axes[0].x)
                                                        : std::atan2(-i.axes[1].x, i.axes[1].y);
    const double roll = std::abs(std::cos(pitch)) > 1e-9 ? std::atan2(i.axes[1].z, i.axes[2].z) : 0;
    const double values[]{i.position.x * 1000, i.position.y * 1000, i.position.z * 1000,
                          yaw * degrees,       pitch * degrees,     roll * degrees};
    for (int k = 0; k < 6; k++) {
        auto *input = new DecimalInput;
        input->setDecimals(10);
        input->setRange(k < 3 ? -1e7 : -360, k < 3 ? 1e7 : 360);
        input->setValue(values[k]);
        input->setMinimumWidth(95);
        table_->setCellWidget(row, k + 1, input);
    }
    for (int k = 7; k < 9; k++) {
        auto *check = new QCheckBox;
        check->setChecked(k == 7 ? i.mirrorHorizontal : i.mouthAnchored);
        table_->setCellWidget(row, k, check);
    }
}
DraftAssembly DraftAssemblyDialog::assembly() const {
    DraftAssembly a;
    a.instances.clear();
    a.includeCad = includeCad_->isChecked();
    a.cadTargetEdgeLength = cadSize_->value() * .001;
    for (int r = 0; r < table_->rowCount(); r++) {
        DraftInstance i;
        i.id = table_->item(r, 9)->text().toUInt();
        i.name = table_->item(r, 0)->text().toStdString();
        const auto v = [&](int c) {
            return static_cast<DecimalInput *>(table_->cellWidget(r, c))->value();
        };
        i.position = {v(1) * .001, v(2) * .001, v(3) * .001};
        const double z = v(4) / degrees, y = v(5) / degrees, x = v(6) / degrees;
        const double cz = std::cos(z), sz = std::sin(z), cy = std::cos(y), sy = std::sin(y),
                     cx = std::cos(x), sx = std::sin(x);
        i.axes = {Vec3{cz * cy, sz * cy, -sy},
                  Vec3{cz * sy * sx - sz * cx, sz * sy * sx + cz * cx, cy * sx},
                  Vec3{cz * sy * cx + sz * sx, sz * sy * cx - cz * sx, cy * cx}};
        i.mirrorHorizontal = static_cast<QCheckBox *>(table_->cellWidget(r, 7))->isChecked();
        i.mouthAnchored = static_cast<QCheckBox *>(table_->cellWidget(r, 8))->isChecked();
        a.instances.push_back(std::move(i));
    }
    return a;
}
} // namespace hornsculpt
