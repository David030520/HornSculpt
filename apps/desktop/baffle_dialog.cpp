#include "baffle_dialog.hpp"
#include "numeric_input.hpp"
#include "simulation_request.hpp"
#include <QComboBox>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QGridLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QStringList>
#include <QVBoxLayout>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hornsculpt {
BaffleDialog::BaffleDialog(std::optional<InfiniteBaffle> baffle, QWidget *parent)
    : QDialog(parent), committed_(baffle) {
    setWindowTitle(tr("Schallwand"));
    setObjectName("baffleDialog");
    auto *layout = new QVBoxLayout(this);
    auto *form = new QFormLayout;
    mode_ = new QComboBox;
    mode_->setObjectName("baffleMode");
    mode_->addItems({tr("Freifeld"), tr("Unendliche Schallwand · Luftkavität")});
    mode_->setToolTip(tr("Unendliche schallharte Ebene mit ebener Apertur. Das geschlossene, "
                         "positiv orientierte Netz begrenzt die Luftkavität hinter der Ebene. "
                         "Aperturflächen schließen dieses Luftvolumen nur für die Kopplung zum "
                         "äußeren Halbraum. Kein äußerer Schallkörper. Nur dichter GPU-Solver "
                         "und Fernfeld. Treiber und Aperturen können gekoppelt in denselben "
                         "äußeren Halbraum abstrahlen."));
    mode_->setWhatsThis(mode_->toolTip());
    form->addRow(tr("Modell"), mode_);
    layout->addLayout(form);
    auto *grid = new QGridLayout;
    const std::array<const char *, 3> coordinates{"X", "Y", "Z"};
    const auto value = baffle.value_or(InfiniteBaffle{});
    exteriorDrivenMesh_ = value.exteriorDrivenMesh;
    const std::array<double, 3> origin{value.origin.x, value.origin.y, value.origin.z};
    const std::array<double, 3> normal{value.normal.x, value.normal.y, value.normal.z};
    grid->addWidget(new QLabel(tr("Ebenenursprung")), 1, 0);
    grid->addWidget(new QLabel(tr("Normale zum Außenraum")), 2, 0);
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        grid->addWidget(new QLabel(QString::fromLatin1(coordinates[i])), 0, int(i) + 1);
        origin_[i] = new DecimalInput;
        origin_[i]->setObjectName(QStringLiteral("baffleOrigin") + coordinates[i]);
        origin_[i]->setRange(-1e9, 1e9);
        origin_[i]->setDecimals(6);
        origin_[i]->setSuffix(" mm");
        origin_[i]->setValue(origin[i] * 1000);
        origin_[i]->setToolTip(tr("Punkt der unendlichen Ebene in Modellkoordinaten [mm]. "
                                  "Die gesamte Luftkavität liegt auf ihrer negativen Seite."));
        origin_[i]->setWhatsThis(origin_[i]->toolTip());
        normal_[i] = new DecimalInput;
        normal_[i]->setObjectName(QStringLiteral("baffleNormal") + coordinates[i]);
        normal_[i]->setRange(-1e6, 1e6);
        normal_[i]->setDecimals(12);
        normal_[i]->setSingleStep(.1);
        normal_[i]->setValue(normal[i]);
        normal_[i]->setToolTip(
            tr("Richtung von der Luftkavität in den äußeren Halbraum; "
               "dimensionslose XYZ-Komponenten, beim Übernehmen normiert. "
               "Nullvektor ist ungültig. Die Apertur-Netznormalen zeigen "
               "in diese Richtung. Auf Driven zeigt die Netznormale aus der "
               "Luft; positive physische Anregung in die Luft ist entgegengesetzt."));
        normal_[i]->setWhatsThis(normal_[i]->toolTip());
        grid->addWidget(origin_[i], 1, int(i) + 1);
        grid->addWidget(normal_[i], 2, int(i) + 1);
    }
    layout->addLayout(grid);
    auto *faces = new QFormLayout;
    aperture_ = new QLineEdit;
    aperture_->setObjectName("baffleApertureFaces");
    QStringList ids;
    for (const auto id : value.apertureFaceIds)
        ids.append(QString::number(id));
    aperture_->setText(ids.join(", "));
    aperture_->setToolTip(tr("IDs der ebenen Aperturflächen, durch Komma getrennt. "
                             "Positive, eindeutige IDs; Flächen müssen als Wall zugeordnet sein. "
                             "Alle übrigen Flächen begrenzen die Luftkavität. Die Apertur ist "
                             "eine Kopplungsfläche, keine physische Abschlusswand."));
    aperture_->setWhatsThis(aperture_->toolTip());
    faces->addRow(tr("Apertur-Flächen"), aperture_);
    layout->addLayout(faces);
    error_ = new QLabel;
    error_->setObjectName("baffleError");
    error_->setWordWrap(true);
    layout->addWidget(error_);
    buttons_ = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    buttons_->button(QDialogButtonBox::Cancel)->setText(tr("Abbrechen"));
    layout->addWidget(buttons_);
    mode_->setCurrentIndex(baffle ? 1 : 0);
    connect(mode_, &QComboBox::currentIndexChanged, this, [this] { validateInput(); });
    for (const auto &inputs : {origin_, normal_})
        for (auto *input : inputs)
            connect(input, &QDoubleSpinBox::valueChanged, this, [this] { validateInput(); });
    connect(aperture_, &QLineEdit::textChanged, this, [this] { validateInput(); });
    connect(buttons_, &QDialogButtonBox::accepted, this, &BaffleDialog::accept);
    connect(buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
    validateInput();
}

InfiniteBaffle BaffleDialog::inputBaffle() const {
    InfiniteBaffle value;
    value.exteriorDrivenMesh = exteriorDrivenMesh_;
    value.origin = {origin_[0]->value() * .001, origin_[1]->value() * .001,
                    origin_[2]->value() * .001};
    const double length = std::hypot(normal_[0]->value(), normal_[1]->value(), normal_[2]->value());
    if (!std::isfinite(length) || length == 0)
        throw std::invalid_argument("Die Ebenennormale darf nicht null sein.");
    value.normal = {normal_[0]->value() / length, normal_[1]->value() / length,
                    normal_[2]->value() / length};
    for (const auto &entry : aperture_->text().split(',')) {
        bool ok = false;
        const auto id = entry.trimmed().toULongLong(&ok);
        if (!ok || id > std::numeric_limits<std::uint32_t>::max())
            throw std::invalid_argument("Apertur-Flächen als positive IDs mit Komma trennen.");
        value.apertureFaceIds.push_back(static_cast<std::uint32_t>(id));
    }
    validateInfiniteBaffleSettings(value);
    return value;
}

void BaffleDialog::validateInput() {
    const bool enabled = mode_->currentIndex() == 1;
    for (const auto &inputs : {origin_, normal_})
        for (auto *input : inputs)
            input->setEnabled(enabled);
    aperture_->setEnabled(enabled);
    QString error;
    if (enabled)
        try {
            (void)inputBaffle();
        } catch (const std::invalid_argument &failure) {
            error = QString::fromUtf8(failure.what());
        }
    error_->setText(error);
    buttons_->button(QDialogButtonBox::Ok)->setEnabled(error.isEmpty());
}

void BaffleDialog::accept() {
    for (const auto &inputs : {origin_, normal_})
        for (auto *input : inputs)
            input->interpretText();
    validateInput();
    if (!buttons_->button(QDialogButtonBox::Ok)->isEnabled())
        return;
    committed_ = mode_->currentIndex() == 1 ? std::optional(inputBaffle()) : std::nullopt;
    QDialog::accept();
}

std::optional<InfiniteBaffle> BaffleDialog::infiniteBaffle() const {
    return committed_;
}
} // namespace hornsculpt
