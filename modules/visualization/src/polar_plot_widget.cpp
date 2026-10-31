#include "chart_export.hpp"
#include "chart_scaling.hpp"
#include <hornsculpt/polar_plot_widget.hpp>

#include <QPainter>
#include <QPainterPath>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace hornsculpt {
PolarPlotWidget::PolarPlotWidget(QString title, QColor color, QWidget *parent)
    : QWidget(parent), title_(std::move(title)), color_(std::move(color)) {
    setMinimumSize(260, 210);
    enableChartExport(*this, [this] {
        if (auto s = editPolarScales(*this, {decibels_, angularStep_})) {
            setDecibelScale(s->decibels);
            setAngularTickStep(s->angularStep);
            if (onDecibelScaleChanged)
                onDecibelScaleChanged(decibels_);
        }
    });
    connect(this, &QObject::objectNameChanged, this, [this] {
        angularStep_ = readAngularStep(*this, angularStep_);
        update();
    });
    setToolTip(tr("Relativer Pegel über dem Winkel. 0° liegt auf der akustischen Achse oben; "
                  "positive Winkel laufen im Uhrzeigersinn. Skalierung über das Kontextmenü. "
                  "Die Normalisierung wird gemeinsam für Polars und Isobaren gewählt."));
}
void PolarPlotWidget::setTitle(QString title) {
    title_ = std::move(title);
    update();
}
void PolarPlotWidget::setDecibelScale(PlotDbScale s) {
    validatePlotDbScale(s);
    decibels_ = s;
    update();
}
void PolarPlotWidget::setAngularTickStep(double step) {
    if (!std::isfinite(step) || step <= 0)
        throw std::invalid_argument("Winkelschritt muss endlich und positiv sein.");
    angularStep_ = step;
    saveAngularStep(*this, step);
    update();
}

void PolarPlotWidget::setLevels(std::vector<AngularLevel> levels) {
    setSeries({PolarSeries{{}, color_, std::move(levels)}});
}
void PolarPlotWidget::setSeries(std::vector<PolarSeries> series) {
    bool any = false, full = false;
    for (const auto &curve : series) {
        const auto &levels = curve.levels;
        for (std::size_t index = 0; index < levels.size(); ++index)
            if (!std::isfinite(levels[index].angleRadians) ||
                std::abs(levels[index].angleRadians) > std::numbers::pi + 1e-10 ||
                (index && levels[index].angleRadians <= levels[index - 1].angleRadians))
                throw std::invalid_argument(
                    "Polarwinkel müssen zwischen −180° und 180° streng ansteigen.");
        if (!levels.empty()) {
            any = true;
            full = full || levels.front().angleRadians < -std::numbers::pi / 2 - 1e-10 ||
                   levels.back().angleRadians > std::numbers::pi / 2 + 1e-10;
        }
    }
    fullCircle_ = !any || full;
    series_ = std::move(series);
    update();
}

void PolarPlotWidget::paintEvent(QPaintEvent *) {
    QPainter painter(this);
    painter.fillRect(rect(), palette().base());
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setPen(palette().text().color());
    painter.drawText(QRect(4, 3, width() - 8, 21), Qt::AlignCenter, title_);
    const QRectF area(37, 50, std::max(1, width() - 74),
                      std::max(1, height() - (fullCircle_ ? 103 : 84)));
    const double radius = std::max(1.0, fullCircle_ ? std::min(area.width(), area.height()) / 2
                                                    : std::min(area.width() / 2, area.height()));
    const QPointF center = fullCircle_ ? area.center() : QPointF(area.center().x(), area.bottom());
    QFont ticks = font();
    ticks.setPointSizeF(std::max(7.0, font().pointSizeF() - 1));
    painter.setFont(ticks);
    double previousLabelRadius = -1000;
    for (const double db : plotTicks(decibels_.minimum, decibels_.maximum, decibels_.tickStep,
                                     static_cast<unsigned>(radius * 2))) {
        const double r =
            radius * (db - decibels_.minimum) / (decibels_.maximum - decibels_.minimum);
        const QRectF ring(center.x() - r, center.y() - r, 2 * r, 2 * r);
        painter.setPen(QPen(QColor(181, 190, 197, 145), 1));
        painter.setBrush(Qt::NoBrush);
        if (fullCircle_)
            painter.drawEllipse(ring);
        else
            painter.drawArc(ring, 0, 180 * 16);
        if (db > decibels_.minimum && r < radius - 15 && r - previousLabelRadius >= 14) {
            painter.setPen(palette().text().color());
            painter.drawText(QPointF(center.x() + 4, center.y() - r + 12), QString::number(db));
            previousLabelRadius = r;
        }
    }
    painter.setPen(QPen(QColor(181, 190, 197, 145), 1));
    const QRectF outer(center.x() - radius, center.y() - radius, 2 * radius, 2 * radius);
    if (fullCircle_)
        painter.drawEllipse(outer);
    else
        painter.drawArc(outer, 0, 180 * 16);
    painter.setPen(palette().text().color());
    painter.drawText(QPointF(center.x() + 4, center.y() - radius + 12),
                     QString::number(decibels_.maximum));
    painter.drawText(QPointF(center.x() + 4, center.y() + 12), QString::number(decibels_.minimum));
    const double start = fullCircle_ ? -180 : -90, end = fullCircle_ ? 180 : 90;
    for (const double angle :
         plotTicks(start, end, angularStep_, static_cast<unsigned>(radius * 6))) {
        if (fullCircle_ && angle == 180)
            continue;
        const double radians = angle * std::numbers::pi / 180;
        const QPointF ray(std::sin(radians), -std::cos(radians));
        painter.setPen(QPen(QColor(181, 190, 197, 145), 1));
        painter.drawLine(center, center + ray * radius);
        const auto label = center + ray * (radius + 17);
        painter.setPen(palette().text().color());
        painter.drawText(QRectF(label.x() - 22, label.y() - 8, 44, 16), Qt::AlignCenter,
                         angle == -180 ? QStringLiteral("±180°") : QString::number(angle) + "°");
    }
    for (const auto &series : series_) {
        QPainterPath path;
        bool started = false;
        for (const auto &level : series.levels) {
            if (std::isnan(level.levelDecibels)) {
                started = false;
                continue;
            }
            const double r = radius *
                             (std::clamp(level.levelDecibels, decibels_.minimum, decibels_.maximum) -
                              decibels_.minimum) /
                             (decibels_.maximum - decibels_.minimum);
            const QPointF point(center.x() + r * std::sin(level.angleRadians),
                                center.y() - r * std::cos(level.angleRadians));
            if (started)
                path.lineTo(point);
            else
                path.moveTo(point);
            started = true;
        }
        painter.setPen(QPen(series.color, 1.8));
        painter.drawPath(path);
    }
    painter.setPen(palette().text().color());
    painter.drawText(QRect(0, height() - 19, width(), 17), Qt::AlignCenter,
                     tr("Relativer Pegel [dB]"));
}
} // namespace hornsculpt
