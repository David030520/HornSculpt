#include "angular_response_widget.hpp"
#include "numeric_input.hpp"
#include <QComboBox>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QLocale>
#include <QSignalBlocker>
#include <QStringList>
#include <array>
#include <cmath>
#include <hornsculpt/plot_widget.hpp>
#include <limits>
#include <stdexcept>

namespace hornsculpt {
namespace {
QString angleText(double angle) {
    return QLocale().toString(angle, 'g', 12);
}
} // namespace

AngularResponseWidget::AngularResponseWidget(QWidget *parent) : QWidget(parent) {
    setObjectName("angularResponses");
    auto *layout = new QGridLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    auto *controls = new QHBoxLayout;
    auto *anglesLabel = new QLabel(tr("Kurven [°]"));
    angles_ = new QLineEdit;
    angles_->setObjectName("angularResponseAngles");
    anglesLabel->setBuddy(angles_);
    const auto anglesHelp =
        tr("1 bis 12 verschiedene Winkel von −180° bis +180°, getrennt durch Semikolon. "
           "Beispiel: 0; 15; 22,5; 45. H und V beziehen sich auf den gespeicherten "
           "Beobachtungsrahmen. Zwischen berechneten Winkelpunkten wird linear in dB "
           "interpoliert; außerhalb des Rasters bleiben Lücken.");
    angles_->setToolTip(anglesHelp);
    angles_->setWhatsThis(anglesHelp);
    normalization_ = new QComboBox;
    normalization_->setObjectName("angularResponseNormalization");
    normalization_->addItems({tr("Absolut"), tr("Bezug auf Winkel")});
    const auto normalizationHelp =
        tr("Absolut: RMS-Schalldruckpegel mit Bezug 20 µPa. Bezug auf Winkel: "
           "Bei jeder Frequenz wird der Pegel am Bezugswinkel derselben Ebene abgezogen. "
           "Ein fehlender oder unendlicher Bezugspegel erzeugt eine Lücke.");
    normalization_->setToolTip(normalizationHelp);
    normalization_->setWhatsThis(normalizationHelp);
    auto *referenceLabel = new QLabel(tr("Bezug"));
    reference_ = new DecimalInput;
    reference_->setObjectName("angularResponseReferenceAngle");
    reference_->setDecimals(10);
    reference_->setRange(-180, 180);
    reference_->setSuffix(" °");
    referenceLabel->setBuddy(reference_);
    const auto referenceHelp =
        tr("Bezugswinkel in Grad, Standard 0°. Der Bezug wird für H und V getrennt "
           "aus den gespeicherten Daten bestimmt. Zwischen Winkelpunkten linear in dB; "
           "keine Extrapolation und kein Ersatz für fehlende Werte.");
    reference_->setToolTip(referenceHelp);
    reference_->setWhatsThis(referenceHelp);
    controls->addWidget(anglesLabel);
    controls->addWidget(angles_, 1);
    controls->addWidget(normalization_);
    controls->addWidget(referenceLabel);
    controls->addWidget(reference_);
    layout->addLayout(controls, 0, 0, 1, 2);
    error_ = new QLabel;
    error_->setObjectName("angularResponseError");
    error_->setStyleSheet("color: #a72e2b");
    error_->setWordWrap(true);
    layout->addWidget(error_, 1, 0, 1, 2);
    horizontal_ = new PlotWidget(tr("Frequenzgänge H"), tr("Frequenz [Hz]"), tr("dB SPL"));
    vertical_ = new PlotWidget(tr("Frequenzgänge V"), tr("Frequenz [Hz]"), tr("dB SPL"));
    horizontal_->setObjectName("angularFrequencyH");
    vertical_->setObjectName("angularFrequencyV");
    for (auto *plot : {horizontal_, vertical_}) {
        plot->setFrequencyAxis();
        plot->setToolTip(tr("Nur berechnete Frequenzpunkte. Linien verbinden diese Punkte; "
                            "fehlende Werte unterbrechen die Kurve. "
                            "Rechtsklick: Achsen und Diagrammexport."));
        plot->setWhatsThis(plot->toolTip());
    }
    layout->addWidget(horizontal_, 2, 0);
    layout->addWidget(vertical_, 2, 1);
    layout->setRowStretch(2, 1);
    coverage_ = new QLabel;
    coverage_->setObjectName("angularResponseCoverage");
    layout->addWidget(coverage_, 3, 0, 1, 2);
    connect(angles_, &QLineEdit::textChanged, this, [this] { selectionChanged(); });
    connect(normalization_, &QComboBox::currentIndexChanged, this, [this] { selectionChanged(); });
    connect(reference_, &QDoubleSpinBox::valueChanged, this, [this] { selectionChanged(); });
    setSettings(settings_);
}

void AngularResponseWidget::setSettings(const AngularResponseSettings &settings) {
    validateAngularResponseSettings(settings);
    const QSignalBlocker anglesBlock(angles_), normalizationBlock(normalization_),
        referenceBlock(reference_);
    settings_ = settings;
    QStringList values;
    for (const auto angle : settings_.anglesDegrees)
        values.push_back(angleText(angle));
    angles_->setText(values.join("; "));
    normalization_->setCurrentIndex(settings_.normalizeToAngle ? 1 : 0);
    reference_->setValue(settings_.referenceAngleDegrees);
    reference_->setEnabled(settings_.normalizeToAngle);
    error_->clear();
    error_->hide();
    refresh();
}

void AngularResponseWidget::selectionChanged() {
    auto selected = settings_;
    selected.anglesDegrees.clear();
    for (auto text : angles_->text().split(';')) {
        // Semicolons separate angles, leaving both common decimal separators usable.
        text = text.trimmed().replace(',', '.');
        bool valid = false;
        const double angle = QLocale::c().toDouble(text, &valid);
        if (!valid || !std::isfinite(angle)) {
            error_->setText(tr("Winkel mit Semikolon trennen. Anzeige: letzte gültige Auswahl."));
            error_->show();
            return;
        }
        selected.anglesDegrees.push_back(angle);
    }
    selected.normalizeToAngle = normalization_->currentIndex() == 1;
    selected.referenceAngleDegrees = reference_->value();
    try {
        validateAngularResponseSettings(selected);
    } catch (const std::invalid_argument &error) {
        error_->setText(QString::fromUtf8(error.what()) + tr(" Anzeige: letzte gültige Auswahl."));
        error_->show();
        return;
    }
    error_->clear();
    error_->hide();
    reference_->setEnabled(selected.normalizeToAngle);
    if (selected == settings_)
        return;
    settings_ = std::move(selected);
    refresh();
    if (onChanged)
        onChanged();
}

void AngularResponseWidget::setRows(std::vector<AngularResponseRow> rows, QString context) {
    double previous = 0;
    for (const auto &row : rows) {
        if (!std::isfinite(row.frequencyHz) || row.frequencyHz <= previous)
            throw std::invalid_argument("Frequenzen müssen positiv und streng aufsteigend sein.");
        previous = row.frequencyHz;
    }
    rows_ = std::move(rows);
    context_ = std::move(context);
    refresh();
}

void AngularResponseWidget::refresh() {
    const std::array<QColor, 12> colors{
        QColor(28, 102, 180), QColor(205, 112, 33),  QColor(44, 145, 81),  QColor(166, 68, 150),
        QColor(210, 60, 63),  QColor(35, 145, 155),  QColor(125, 106, 38), QColor(87, 91, 178),
        QColor(153, 80, 50),  QColor(125, 125, 125), QColor(53, 100, 58),  QColor(184, 68, 111)};
    const auto mode = settings_.normalizeToAngle ? AngularResponseNormalization::referenceAngle
                                                 : AngularResponseNormalization::absolute;
    const QString modeLabel = settings_.normalizeToAngle
                                  ? tr("Bezug %1°").arg(angleText(settings_.referenceAngleDegrees))
                                  : tr("absolut");
    bool hasMissing = false;
    for (const bool horizontal : {true, false}) {
        std::vector<AngularFrequencyRow> samples;
        for (const auto &row : rows_)
            samples.push_back({row.frequencyHz, horizontal ? row.horizontal : row.vertical});
        std::vector<PlotSeries> series;
        if (!rows_.empty())
            for (std::size_t i = 0; i < settings_.anglesDegrees.size(); ++i) {
                const auto angle = settings_.anglesDegrees[i];
                PlotSeries curve{angleText(angle) + "°", colors[i], {}};
                for (const auto &sample : angularFrequencyResponse(
                         samples, angle, mode, settings_.referenceAngleDegrees)) {
                    hasMissing |= !sample.levelDecibels.has_value();
                    // NaNs retain missing frequency rows as actual path breaks.
                    curve.points.emplace_back(
                        sample.frequencyHz,
                        sample.levelDecibels.value_or(std::numeric_limits<double>::quiet_NaN()));
                }
                series.push_back(std::move(curve));
            }
        auto *plot = horizontal ? horizontal_ : vertical_;
        plot->setTitle(tr("Frequenzgänge %1 · %2%3")
                           .arg(horizontal ? "H" : "V", modeLabel,
                                context_.isEmpty() ? QString{} : " · " + context_));
        plot->setAxisLabels(tr("Frequenz [Hz]"),
                            settings_.normalizeToAngle ? tr("dB relativ") : tr("dB SPL"));
        plot->setSeries(std::move(series));
    }
    coverage_->setText(
        rows_.empty()
            ? tr("Keine berechneten Winkelfrequenzgänge vorhanden.")
            : tr("%1 berechnete Frequenzen%2")
                  .arg(rows_.size())
                  .arg(hasMissing ? tr(" · Lücken: Winkel oder Bezug ohne Pegel") : QString{}));
}
} // namespace hornsculpt
