#include "phase_response_widget.hpp"
#include "numeric_input.hpp"
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QSignalBlocker>
#include <algorithm>
#include <cmath>
#include <hornsculpt/plot_widget.hpp>
#include <limits>
#include <numbers>

namespace hornsculpt {
namespace {
QDoubleSpinBox *angleInput(const char *name) {
    auto *input = new DecimalInput;
    input->setObjectName(name);
    input->setRange(-180, 180);
    input->setDecimals(3);
    input->setSuffix(" °");
    return input;
}
std::complex<double> sample(const std::vector<double> &angles,
                            const std::vector<std::complex<double>> &values, double angle) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (angles.size() != values.size() || angles.empty() || angle < angles.front() ||
        angle > angles.back())
        return {nan, nan};
    const auto found = std::lower_bound(angles.begin(), angles.end(), angle);
    const auto i = static_cast<std::size_t>(found - angles.begin());
    if (*found == angle || i == 0)
        return values[i];
    const double fraction = (angle - angles[i - 1]) / (angles[i] - angles[i - 1]);
    return (1 - fraction) * values[i - 1] + fraction * values[i];
}
} // namespace
PhaseResponseWidget::PhaseResponseWidget(QWidget *parent) : QWidget(parent) {
    setObjectName("phaseResponses");
    auto *layout = new QGridLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    auto *controls = new QHBoxLayout;
    angle_ = angleInput("phaseResponseAngle");
    reference_ = angleInput("phaseReferenceAngle");
    relative_ = new QCheckBox(tr("Bezug auf Winkel"));
    relative_->setObjectName("phaseRelative");
    unwrap_ = new QCheckBox(tr("Entwickeln"));
    unwrap_->setObjectName("phaseUnwrap");
    delay_ = new DecimalInput;
    delay_->setObjectName("phaseDelay");
    delay_->setRange(-1000, 1000);
    delay_->setDecimals(4);
    delay_->setSuffix(" ms");
    const auto phaseHelp = tr("Phase aus dem komplexen RMS-Schalldruck; exp(+iωt). Zwischen "
                              "Winkelstützpunkten wird komplex linear interpoliert. Fehlende "
                              "Werte und exakte Nullstellen bleiben Lücken. Bezug: Phase des "
                              "gewählten Winkels derselben Ebene abziehen.");
    for (auto *w : {angle_, reference_}) {
        w->setToolTip(phaseHelp);
        w->setWhatsThis(phaseHelp);
    }
    unwrap_->setToolTip(tr("Verbindet benachbarte berechnete Frequenzen über den nächsten "
                           "360°-Zweig. Bei grobem Frequenzraster sind zusätzliche Umläufe "
                           "nicht bestimmbar; kein Nachweis einer absoluten Laufzeit."));
    unwrap_->setWhatsThis(unwrap_->toolTip());
    delay_->setToolTip(tr("Laufzeitkompensation: +360°·f·t zur dargestellten Phase. Einheit ms. "
                          "Ändert ausschließlich die Anzeige."));
    delay_->setWhatsThis(delay_->toolTip());
    controls->addWidget(new QLabel(tr("Winkel")));
    controls->addWidget(angle_);
    controls->addWidget(relative_);
    controls->addWidget(reference_);
    controls->addWidget(unwrap_);
    controls->addWidget(new QLabel(tr("Laufzeit")));
    controls->addWidget(delay_);
    controls->addStretch();
    layout->addLayout(controls, 0, 0, 1, 2);
    horizontal_ = new PlotWidget(tr("Phase H"), tr("Frequenz [Hz]"), tr("Phase [°]"));
    vertical_ = new PlotWidget(tr("Phase V"), tr("Frequenz [Hz]"), tr("Phase [°]"));
    horizontal_->setObjectName("phaseResponseH");
    vertical_->setObjectName("phaseResponseV");
    for (auto *p : {horizontal_, vertical_}) {
        p->setFrequencyAxis();
        p->setToolTip(phaseHelp);
    }
    layout->addWidget(horizontal_, 1, 0);
    layout->addWidget(vertical_, 1, 1);
    layout->setRowStretch(1, 1);
    context_ = new QLabel;
    context_->setObjectName("phaseResponseContext");
    layout->addWidget(context_, 2, 0, 1, 2);
    for (auto *w : {angle_, reference_, delay_})
        connect(w, &QDoubleSpinBox::valueChanged, this, [this] { controlsChanged(); });
    for (auto *w : {relative_, unwrap_})
        connect(w, &QCheckBox::toggled, this, [this] { controlsChanged(); });
    setSettings({});
}
PhaseResponseWidget::~PhaseResponseWidget() {
    for (auto *child : findChildren<QObject *>())
        QObject::disconnect(child, nullptr, this, nullptr);
    onChanged = {};
}
void PhaseResponseWidget::setSettings(const PhaseDisplaySettings &s) {
    validatePhaseDisplaySettings(s);
    const QSignalBlocker a(angle_), b(reference_), c(delay_), d(relative_), e(unwrap_);
    settings_ = s;
    angle_->setValue(s.angleDegrees);
    reference_->setValue(s.referenceAngleDegrees);
    delay_->setValue(s.delaySeconds * 1000);
    relative_->setChecked(s.relative);
    unwrap_->setChecked(s.unwrap);
    reference_->setEnabled(s.relative);
    refresh();
}
void PhaseResponseWidget::setRows(std::vector<double> angles, std::vector<PhaseResponseRow> rows,
                                  QString context) {
    angles_ = std::move(angles);
    rows_ = std::move(rows);
    context_->setText(context);
    refresh();
}
void PhaseResponseWidget::controlsChanged() {
    settings_ = {angle_->value(), reference_->value(), delay_->value() / 1000,
                 relative_->isChecked(), unwrap_->isChecked()};
    reference_->setEnabled(settings_.relative);
    refresh();
    if (onChanged)
        onChanged();
}
void PhaseResponseWidget::refresh() {
    for (unsigned plane = 0; plane < 2; ++plane) {
        PlotSeries series{tr("%1°").arg(settings_.angleDegrees), QColor(36, 91, 155), {}};
        double previous = std::numeric_limits<double>::quiet_NaN();
        double previousRaw = std::numeric_limits<double>::quiet_NaN();
        for (const auto &row : rows_) {
            const auto &values = plane == 0 ? row.horizontal : row.vertical;
            auto pressure = sample(angles_, values, settings_.angleDegrees);
            if (settings_.relative) {
                const auto reference = sample(angles_, values, settings_.referenceAngleDegrees);
                if (std::abs(reference) > 0)
                    pressure /= reference;
                else
                    pressure = {std::numeric_limits<double>::quiet_NaN(), 0};
            }
            double phase = std::numeric_limits<double>::quiet_NaN();
            if (std::isfinite(pressure.real()) && std::isfinite(pressure.imag()) &&
                std::abs(pressure) > 0) {
                double raw = std::arg(pressure) * 180 / std::numbers::pi;
                if (settings_.unwrap && std::isfinite(previousRaw))
                    raw = previousRaw + std::remainder(raw - previousRaw, 360.);
                previousRaw = raw;
                // A known display delay must retain all of its turns; only
                // the sampled physical pressure phase requires unwrapping.
                phase = raw + 360 * row.frequency * settings_.delaySeconds;
                if (!settings_.unwrap)
                    phase = std::remainder(phase, 360.);
            } else
                previousRaw = std::numeric_limits<double>::quiet_NaN();
            // Wrapped discontinuities are gaps, not diagonal phase ramps.
            if (!settings_.unwrap && std::isfinite(previous) && std::isfinite(phase) &&
                std::abs(phase - previous) > 180)
                series.points.emplace_back(row.frequency, std::numeric_limits<double>::quiet_NaN());
            series.points.emplace_back(row.frequency, phase);
            previous = phase;
        }
        auto *plot = plane == 0 ? horizontal_ : vertical_;
        plot->setSeries({std::move(series)});
        plot->setAxisScales({},
                            settings_.unwrap ? PlotAxisScale{} : PlotAxisScale{-180., 180., 90.});
    }
}
} // namespace hornsculpt
