#include <hornsculpt/heatmap_widget.hpp>
#include <hornsculpt/plot_axes.hpp>

#include <QPainter>
#include <algorithm>
#include <cmath>
#include <numbers>

namespace hornsculpt {
namespace {
// The central disk is a visual origin, not a physical distance or acoustic result.
constexpr double innerFraction = .15;
QPointF polarPoint(QPointF center, double radius, double degrees) {
    const double angle = degrees * std::numbers::pi / 180.;
    return center + QPointF(radius * std::sin(angle), -radius * std::cos(angle));
}
} // namespace

void HeatmapWidget::setRadialLayout(bool enabled) {
    radial_ = enabled;
    rasterDirty_ = true;
    cursor_.reset();
    if (enabled)
        setToolTip(tr("Frequenz radial logarithmisch: tief innen, hoch außen. "
                      "0° oben, positive Winkel im Uhrzeigersinn. Der Innenkreis ist ein "
                      "Bezugssymbol, keine Entfernungsskala. Fehlende Winkel bleiben leer. "
                      "Maus zum Ablesen; Skalierung und PNG-Export im Kontextmenü."));
    update();
}

std::optional<QPointF> HeatmapWidget::radialCoordinates(QPointF position, QRectF plot) const {
    const auto delta = position - plot.center();
    const double radius = std::hypot(delta.x(), delta.y()) / (plot.width() / 2);
    if (radius < innerFraction || radius > 1)
        return {};
    const double angle = std::atan2(delta.x(), -delta.y()) * 180 / std::numbers::pi;
    const auto bounds = visibleBounds();
    if (angle < bounds[2] || angle > bounds[3])
        return {};
    // A single frequency occupies a narrow ring instead of inventing a frequency band.
    if (frequencies_.size() == 1) {
        const double location =
            innerFraction + (1 - innerFraction) *
                                (std::log(frequencies_.front()) - std::log(bounds[0])) /
                                (std::log(bounds[1]) - std::log(bounds[0]));
        if (std::abs(radius - location) > 2. / (plot.width() / 2))
            return {};
        return QPointF(logarithmicFrequencies_.front(), angle);
    }
    const double logarithm = std::log(bounds[0]) + (radius - innerFraction) / (1 - innerFraction) *
                                                       (std::log(bounds[1]) - std::log(bounds[0]));
    return QPointF(logarithm, angle);
}

void HeatmapWidget::paintRadial(QPainter &painter) {
    const auto plot = plotRectangle();
    const auto bounds = visibleBounds();
    const auto center = plot.center();
    const double outer = plot.width() / 2;
    const double inner = outer * innerFraction;
    const double ratio = devicePixelRatioF();
    updateRaster({std::max(1, static_cast<int>(std::ceil(plot.width() * ratio))),
                  std::max(1, static_cast<int>(std::ceil(plot.height() * ratio)))},
                 ratio);
    painter.drawImage(plot, raster_);
    painter.setRenderHint(QPainter::Antialiasing);
    QFont ticks = font();
    ticks.setPointSizeF(std::max(7., font().pointSizeF() - 1));
    painter.setFont(ticks);
    const QPen grid(QColor(115, 125, 133, 125), 1, Qt::DotLine);
    const auto frequencies =
        scales_.x.tickStep
            ? plotTicks(bounds[0], bounds[1], scales_.x.tickStep, static_cast<unsigned>(outer))
            : frequencyTicks(bounds[0], bounds[1], outer);
    painter.setBrush(Qt::NoBrush);
    for (double frequency : frequencies) {
        const double radius = inner + (outer - inner) *
                                          (std::log(frequency) - std::log(bounds[0])) /
                                          (std::log(bounds[1]) - std::log(bounds[0]));
        painter.setPen(grid);
        painter.drawEllipse(center, radius, radius);
        painter.setPen(palette().text().color());
        const auto point = polarPoint(center, radius, 135);
        const QRectF label(point.x() - 23, point.y() - 8, 46, 16);
        painter.fillRect(label, palette().base());
        painter.drawText(label, Qt::AlignCenter, frequencyLabel(frequency));
    }
    for (double angle :
         plotTicks(-180, 180, scales_.y.tickStep.value_or(30), static_cast<unsigned>(outer * 4))) {
        if (angle == 180)
            continue;
        painter.setPen(grid);
        painter.drawLine(polarPoint(center, inner, angle), polarPoint(center, outer, angle));
        painter.setPen(palette().text().color());
        const auto point = polarPoint(center, outer + 17, angle);
        painter.drawText(QRectF(point.x() - 23, point.y() - 8, 46, 16), Qt::AlignCenter,
                         QString::number(angle, 'g', 4) + "°");
    }
    painter.setPen(palette().mid().color());
    painter.drawEllipse(center, outer, outer);
    painter.setBrush(palette().alternateBase());
    painter.drawEllipse(center, inner, inner);
    painter.setBrush(Qt::NoBrush);
    const QRectF bar(plot.right() + 45, plot.top(), 13, plot.height());
    for (int row = 0; row < static_cast<int>(std::ceil(bar.height())); ++row)
        painter.fillRect(QRectF(bar.left(), bar.top() + row, bar.width(), 1.1),
                         levelColor(decibels_.maximum - row / bar.height() *
                                                            (decibels_.maximum - decibels_.minimum),
                                    decibels_));
    painter.setPen(palette().text().color());
    painter.drawRect(bar);
    painter.drawText(QRectF(bar.left(), bar.top() - 20, 30, 16), "dB");
    for (double level : plotTicks(decibels_.minimum, decibels_.maximum, decibels_.tickStep,
                                  static_cast<unsigned>(bar.height()))) {
        const double y = bar.top() + (decibels_.maximum - level) /
                                         (decibels_.maximum - decibels_.minimum) * bar.height();
        painter.drawText(QRectF(bar.right() + 5, y - 8, 42, 16), Qt::AlignVCenter,
                         QString::number(level, 'g', 4));
    }
    painter.drawText(QRectF(4, height() - 43, width() - 8, 17), Qt::AlignCenter,
                     tr("Frequenz [Hz] · innen → außen · logarithmisch"));
    if (cursor_) {
        const auto coordinate = radialCoordinates(*cursor_, plot);
        if (coordinate) {
            const auto level = interpolatedLevel(coordinate->x(), coordinate->y());
            painter.drawEllipse(*cursor_, 3, 3);
            painter.drawText(
                QRectF(4, height() - 22, width() - 8, 18), Qt::AlignCenter,
                frequencyLabel(std::exp(coordinate->x())) + " Hz · " +
                    QString::number(coordinate->y(), 'f', 1) + "° · " +
                    (level ? QString::number(*level, 'f', 1) + " dB" : tr("keine Daten")));
        }
    }
}
} // namespace hornsculpt
