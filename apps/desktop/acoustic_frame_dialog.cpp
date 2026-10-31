#include "acoustic_frame_dialog.hpp"
#include "numeric_input.hpp"

#include <QColor>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QGridLayout>
#include <QLabel>
#include <QLocale>
#include <QPalette>
#include <QPushButton>
#include <QSignalBlocker>
#include <QVBoxLayout>
#include <cmath>
#include <stdexcept>

namespace hornsculpt {
namespace {
constexpr std::array<Vec3, 6> axes{
    {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}}};
constexpr int freeDirection = 6;
constexpr std::array<const char *, 3> coordinateNames{"X", "Y", "Z"};

int presetIndex(Vec3 vector) {
    for (std::size_t i = 0; i < axes.size(); ++i)
        if (vector.x == axes[i].x && vector.y == axes[i].y && vector.z == axes[i].z)
            return static_cast<int>(i);
    return freeDirection;
}

Vec3 vectorValue(const std::array<DecimalInput *, 3> &inputs) {
    return {inputs[0]->value(), inputs[1]->value(), inputs[2]->value()};
}

void setVector(const std::array<DecimalInput *, 3> &inputs, Vec3 vector, double scale = 1) {
    inputs[0]->setValue(vector.x * scale);
    inputs[1]->setValue(vector.y * scale);
    inputs[2]->setValue(vector.z * scale);
}

bool finite(Vec3 vector) {
    return std::isfinite(vector.x) && std::isfinite(vector.y) && std::isfinite(vector.z);
}
} // namespace

AcousticFrameDialog::AcousticFrameDialog(const AcousticFrame &sourceFrame,
                                         std::optional<AcousticFrame> overrideFrame,
                                         QWidget *parent)
    : QDialog(parent), sourceFrame_(sourceFrame), customFrame_(overrideFrame.value_or(sourceFrame)),
      committedFrame_(overrideFrame), showingCustom_(overrideFrame.has_value()) {
    setWindowTitle(tr("Abstrahlung"));
    setObjectName(QStringLiteral("acousticFrameDialog"));
    auto *layout = new QVBoxLayout(this);
    auto *form = new QFormLayout;
    mode_ = new QComboBox;
    mode_->setObjectName(QStringLiteral("acousticFrameMode"));
    mode_->addItems({tr("Aus Geometrie"), tr("Benutzerdefiniert")});
    mode_->setToolTip(tr("Ursprung und Richtungen beziehen sich auf die Modellkoordinaten. "
                         "Eine Kameradrehung ändert die Messrichtung nicht."));
    mode_->setWhatsThis(mode_->toolTip());
    form->addRow(tr("Koordinatensystem"), mode_);
    layout->addLayout(form);

    auto *grid = new QGridLayout;
    for (std::size_t i = 0; i < coordinateNames.size(); ++i)
        grid->addWidget(new QLabel(QString::fromLatin1(coordinateNames[i])), 0,
                        static_cast<int>(i) + 2);
    const auto addCoordinates = [&](int row, const QString &prefix, bool position,
                                    const QString &tooltip) {
        std::array<DecimalInput *, 3> controls{};
        for (std::size_t i = 0; i < controls.size(); ++i) {
            auto *input = new DecimalInput;
            input->setObjectName(prefix + QString::fromLatin1(coordinateNames[i]));
            input->setDecimals(position ? 6 : 8);
            input->setRange(position ? -1e9 : -1e6, position ? 1e9 : 1e6);
            input->setSingleStep(position ? 1 : .1);
            if (position)
                input->setSuffix(QStringLiteral(" mm"));
            input->setToolTip(tooltip);
            input->setWhatsThis(tooltip);
            input->setMinimumWidth(112);
            input->setMaximumWidth(160);
            grid->addWidget(input, row, static_cast<int>(i) + 2);
            controls[i] = input;
        }
        return controls;
    };
    grid->addWidget(new QLabel(tr("Ursprung")), 1, 0);
    origin_ = addCoordinates(
        1, QStringLiteral("acousticOrigin"), true,
        tr("Ursprung in Modellkoordinaten [mm]. Bezugspunkt für Phase und Fernfeldrichtungen; "
           "der Fernfeldpegel wird auf 1 m bezogen angegeben."));
    const auto addDirection = [&](DirectionControls &controls, int row, const QString &label,
                                  const QString &prefix, const QString &tooltip) {
        grid->addWidget(new QLabel(label), row, 0);
        controls.preset = new QComboBox;
        controls.preset->setObjectName(prefix + QStringLiteral("Preset"));
        controls.preset->addItems({QStringLiteral("+X"), QStringLiteral("−X"), QStringLiteral("+Y"),
                                   QStringLiteral("−Y"), QStringLiteral("+Z"), QStringLiteral("−Z"),
                                   tr("Frei")});
        controls.preset->setToolTip(tooltip);
        controls.preset->setWhatsThis(tooltip);
        grid->addWidget(controls.preset, row, 1);
        controls.coordinates = addCoordinates(row, prefix, false, tooltip);
    };
    addDirection(forward_, 2, tr("Richtung (0°)"), QStringLiteral("acousticForward"),
                 tr("0°-Messrichtung in Modellkoordinaten, unabhängig von der Kamera. "
                    "Frei: XYZ-Komponenten eines Richtungsvektors; wird normiert."));
    addDirection(horizontal_, 3, tr("H-Richtung (+90°)"), QStringLiteral("acousticHorizontal"),
                 tr("Positive horizontale Messrichtung in Modellkoordinaten. Legt die Drehung "
                    "um die 0°-Achse fest. Wird senkrecht zur 0°-Richtung projiziert und normiert; "
                    "parallele Richtungen sind ungültig."));
    auto *verticalLabel = new QLabel(tr("V-Richtung (+90°)"));
    const auto verticalTooltip = tr("Positive vertikale Messrichtung: 0°-Richtung × H-Richtung. "
                                    "Rechtshändiges Koordinatensystem.");
    verticalLabel->setToolTip(verticalTooltip);
    verticalLabel->setWhatsThis(verticalTooltip);
    grid->addWidget(verticalLabel, 4, 0, 1, 2);
    for (std::size_t i = 0; i < vertical_.size(); ++i) {
        vertical_[i] = new QLabel;
        vertical_[i]->setObjectName(QStringLiteral("acousticVertical") +
                                    QString::fromLatin1(coordinateNames[i]));
        vertical_[i]->setToolTip(verticalTooltip);
        vertical_[i]->setWhatsThis(verticalTooltip);
        grid->addWidget(vertical_[i], 4, static_cast<int>(i) + 2);
    }
    layout->addLayout(grid);
    error_ = new QLabel;
    error_->setObjectName(QStringLiteral("acousticFrameError"));
    error_->setWordWrap(true);
    auto palette = error_->palette();
    palette.setColor(QPalette::WindowText, QColor(180, 40, 40));
    error_->setPalette(palette);
    layout->addWidget(error_);
    buttons_ = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    buttons_->button(QDialogButtonBox::Ok)->setObjectName(QStringLiteral("acousticFrameAccept"));
    buttons_->button(QDialogButtonBox::Cancel)
        ->setObjectName(QStringLiteral("acousticFrameCancel"));
    buttons_->button(QDialogButtonBox::Cancel)->setText(tr("Abbrechen"));
    layout->addWidget(buttons_);

    mode_->setCurrentIndex(showingCustom_ ? 1 : 0);
    connect(mode_, &QComboBox::currentIndexChanged, this, [this] { modeChanged(); });
    for (auto *input : origin_)
        connect(input, &QDoubleSpinBox::valueChanged, this, [this] { validateFrame(); });
    for (auto *controls : {&forward_, &horizontal_}) {
        connect(controls->preset, &QComboBox::currentIndexChanged, this,
                [this, controls] { presetChanged(*controls); });
        for (auto *input : controls->coordinates)
            connect(input, &QDoubleSpinBox::valueChanged, this,
                    [this, controls] { coordinatesChanged(*controls); });
    }
    connect(buttons_, &QDialogButtonBox::accepted, this, &AcousticFrameDialog::accept);
    connect(buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
    displayFrame(showingCustom_ ? customFrame_ : sourceFrame_);
}

std::optional<AcousticFrame> AcousticFrameDialog::observationFrame() const {
    return committedFrame_;
}

AcousticFrame AcousticFrameDialog::inputFrame() const {
    auto origin = vectorValue(origin_);
    origin.x *= .001;
    origin.y *= .001;
    origin.z *= .001;
    return {origin, vectorValue(forward_.coordinates), vectorValue(horizontal_.coordinates)};
}

void AcousticFrameDialog::displayFrame(const AcousticFrame &frame) {
    loading_ = true;
    setVector(origin_, frame.origin, 1000);
    setVector(forward_.coordinates, frame.forward);
    setVector(horizontal_.coordinates, frame.horizontal);
    forward_.preset->setCurrentIndex(presetIndex(frame.forward));
    horizontal_.preset->setCurrentIndex(presetIndex(frame.horizontal));
    loading_ = false;
    updateReadOnly();
    validateFrame();
}

void AcousticFrameDialog::modeChanged() {
    if (showingCustom_)
        customFrame_ = inputFrame();
    showingCustom_ = mode_->currentIndex() == 1;
    displayFrame(showingCustom_ ? customFrame_ : sourceFrame_);
}

void AcousticFrameDialog::presetChanged(DirectionControls &controls) {
    if (loading_)
        return;
    const int index = controls.preset->currentIndex();
    if (index >= 0 && index < freeDirection) {
        loading_ = true;
        setVector(controls.coordinates, axes[static_cast<std::size_t>(index)]);
        loading_ = false;
    }
    updateReadOnly();
    validateFrame();
}

void AcousticFrameDialog::coordinatesChanged(DirectionControls &controls) {
    if (loading_)
        return;
    // Free entry remains free while typing, even when one intermediate vector
    // happens to coincide with an axis preset.
    if (controls.preset->currentIndex() != freeDirection) {
        const QSignalBlocker blocked(controls.preset);
        controls.preset->setCurrentIndex(presetIndex(vectorValue(controls.coordinates)));
    }
    updateReadOnly();
    validateFrame();
}

void AcousticFrameDialog::updateReadOnly() {
    for (auto *input : origin_)
        input->setReadOnly(!showingCustom_);
    for (auto *controls : {&forward_, &horizontal_}) {
        controls->preset->setEnabled(showingCustom_);
        for (auto *input : controls->coordinates)
            input->setReadOnly(!showingCustom_ ||
                               controls->preset->currentIndex() != freeDirection);
    }
}

void AcousticFrameDialog::validateFrame() {
    if (loading_)
        return;
    validatedFrame_.reset();
    const auto frame = showingCustom_ ? inputFrame() : sourceFrame_;
    QString error;
    try {
        validatedFrame_ = normalizedAcousticFrame(frame);
    } catch (const std::invalid_argument &) {
        if (!finite(frame.origin) || !finite(frame.forward) || !finite(frame.horizontal))
            error = tr("Koordinaten und Richtungen müssen endlich sein.");
        else if (std::hypot(frame.forward.x, frame.forward.y, frame.forward.z) == 0)
            error = tr("Richtung (0°) darf nicht null sein.");
        else if (std::hypot(frame.horizontal.x, frame.horizontal.y, frame.horizontal.z) == 0)
            error = tr("H-Richtung darf nicht null sein.");
        else
            error = tr("Richtung (0°) und H-Richtung dürfen nicht parallel sein.");
    }
    error_->setText(error);
    buttons_->button(QDialogButtonBox::Ok)->setEnabled(validatedFrame_.has_value());
    if (validatedFrame_) {
        const auto vertical = verticalDirection(*validatedFrame_);
        const std::array<double, 3> values{vertical.x, vertical.y, vertical.z};
        for (std::size_t i = 0; i < vertical_.size(); ++i)
            vertical_[i]->setText(locale().toString(values[i], 'f', 6));
    } else {
        for (auto *label : vertical_)
            label->setText(QStringLiteral("—"));
    }
}

void AcousticFrameDialog::accept() {
    // Commit a focused numeric editor before validating, including Enter on
    // the dialog's default button while the editor still owns focus.
    for (auto *input : origin_)
        input->interpretText();
    for (auto *controls : {&forward_, &horizontal_})
        for (auto *input : controls->coordinates)
            input->interpretText();
    validateFrame();
    if (!validatedFrame_)
        return;
    committedFrame_ = showingCustom_ ? validatedFrame_ : std::nullopt;
    QDialog::accept();
}

} // namespace hornsculpt
