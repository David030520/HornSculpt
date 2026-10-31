#include "chart_export.hpp"
#include "chart_scaling.hpp"
#include <hornsculpt/heatmap_widget.hpp>
#include <hornsculpt/level_palette.hpp>
#include <hornsculpt/plot_axes.hpp>

#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hornsculpt {
QColor HeatmapWidget::levelColor(double level, PlotDbScale scale) {
    return steppedLevelColor(level, scale);
}
namespace {
struct Interval {
    std::size_t low{}, high{};
    double fraction{};
};
Interval interval(const std::vector<double> &axis, double value) {
    if (value <= axis.front())
        return {0, 0, 0};
    if (value >= axis.back())
        return {axis.size() - 1, axis.size() - 1, 0};
    const auto upper = std::upper_bound(axis.begin(), axis.end(), value);
    const auto high = static_cast<std::size_t>(upper - axis.begin()), low = high - 1;
    return {low, high, (value - axis[low]) / (axis[high] - axis[low])};
}
std::optional<double> sample(const std::vector<double> &levels, std::size_t angles,
                             Interval frequency, Interval angle) {
    const std::array<std::size_t, 2> rows{frequency.low, frequency.high},
        columns{angle.low, angle.high};
    const std::array<double, 2> x{1 - frequency.fraction, frequency.fraction},
        y{1 - angle.fraction, angle.fraction};
    double value = 0;
    for (unsigned row = 0; row < 2; ++row)
        for (unsigned column = 0; column < 2; ++column) {
            const double weight = x[row] * y[column];
            if (weight == 0)
                continue;
            const double level = levels[rows[row] * angles + columns[column]];
            if (std::isnan(level))
                return {};
            value += weight * level;
        }
    return std::isnan(value) ? std::nullopt : std::optional{value};
}
} // namespace

HeatmapWidget::HeatmapWidget(QString title, QWidget *parent)
    : QWidget(parent), title_(std::move(title)) {
    setMinimumSize(320, 210);
    enableChartExport(*this, [this] {
        if (auto s = editHeatmapScales(*this, {scales_, decibels_})) {
            setAxisScales(s->axes.x, s->axes.y);
            setDecibelScale(s->decibels);
            if (onDecibelScaleChanged)
                onDecibelScaleChanged(decibels_);
        }
    });
    connect(this, &QObject::objectNameChanged, this, [this] {
        auto x = readAxisScale(*this, "x", scales_.x);
        try {
            validatePlotAxisScale(x, true);
        } catch (const std::invalid_argument &) {
            x = scales_.x;
        }
        scales_ = {x, readAxisScale(*this, "y", scales_.y)};
        rasterDirty_ = true;
        update();
    });
    setMouseTracking(true);
    setToolTip(tr("Relative Pegel in einstellbaren dB-Farbstufen. Frequenz logarithmisch, "
                  "Winkel vertikal. Mauszeiger oder linke Maustaste zum Ablesen. "
                  "Fehlende Werte bleiben leer."));
}
void HeatmapWidget::setTitle(QString title) {
    title_ = std::move(title);
    update();
}
void HeatmapWidget::setAxisScales(PlotAxisScale x, PlotAxisScale y) {
    validatePlotAxisScale(x, true);
    validatePlotAxisScale(y);
    scales_ = {x, y};
    saveAxisScale(*this, "x", x);
    saveAxisScale(*this, "y", y);
    rasterDirty_ = true;
    update();
}
void HeatmapWidget::setDecibelScale(PlotDbScale s) {
    validatePlotDbScale(s);
    decibels_ = s;
    rasterDirty_ = true;
    update();
}
std::array<double, 4> HeatmapWidget::visibleBounds() const {
    const double f = frequencies_.front();
    return {scales_.x.minimum.value_or(frequencies_.size() == 1 ? f / 1.05 : f),
            scales_.x.maximum.value_or(frequencies_.size() == 1 ? f * 1.05 : frequencies_.back()),
            scales_.y.minimum.value_or(angles_.front()),
            scales_.y.maximum.value_or(angles_.back())};
}

void HeatmapWidget::setLevels(std::vector<double> angles, std::vector<double> frequencies,
                              std::vector<double> levels) {
    if (angles.empty() && frequencies.empty() && levels.empty()) {
        angles_.clear();
        frequencies_.clear();
        logarithmicFrequencies_.clear();
        levels_.clear();
        raster_ = {};
        cursor_.reset();
        dragging_ = false;
        rasterDirty_ = true;
        update();
        return;
    }
    if (angles.size() < 2 || frequencies.empty() ||
        frequencies.size() > std::numeric_limits<std::size_t>::max() / angles.size() ||
        levels.size() != angles.size() * frequencies.size())
        throw std::invalid_argument("Ungültige Isobardaten.");
    for (std::size_t index = 0; index < angles.size(); ++index)
        if (!std::isfinite(angles[index]) || (index && angles[index] <= angles[index - 1]))
            throw std::invalid_argument("Isobarwinkel müssen streng ansteigen.");
    std::vector<double> logarithms;
    logarithms.reserve(frequencies.size());
    for (std::size_t index = 0; index < frequencies.size(); ++index) {
        if (!std::isfinite(frequencies[index]) || frequencies[index] <= 0 ||
            (index && frequencies[index] <= frequencies[index - 1]))
            throw std::invalid_argument(
                "Isobarfrequenzen müssen positiv sein und streng ansteigen.");
        logarithms.push_back(std::log(frequencies[index]));
    }
    angles_ = std::move(angles);
    frequencies_ = std::move(frequencies);
    logarithmicFrequencies_ = std::move(logarithms);
    levels_ = std::move(levels);
    cursor_.reset();
    rasterDirty_ = true;
    update();
}

QRectF HeatmapWidget::plotRectangle() const {
    if (radial_) {
        const double side = std::max(1, std::min(width() - 160, height() - 120));
        return QRectF(45 + (width() - 160 - side) / 2., 45, side, side);
    }
    return QRectF(58, 29, std::max(1, width() - (frequencies_.size() == 1 ? 82 : 145)),
                  std::max(1, height() - 94));
}

std::optional<double> HeatmapWidget::interpolatedLevel(double logarithm, double angle) const {
    if (levels_.empty())
        return {};
    if (logarithm < logarithmicFrequencies_.front() || logarithm > logarithmicFrequencies_.back() ||
        angle < angles_.front() || angle > angles_.back())
        return {};
    return sample(levels_, angles_.size(), interval(logarithmicFrequencies_, logarithm),
                  interval(angles_, angle));
}

void HeatmapWidget::updateRaster(QSize size, double pixelRatio) {
    if (!rasterDirty_ && raster_.size() == size && raster_.devicePixelRatio() == pixelRatio)
        return;
    raster_ = QImage(size, QImage::Format_RGB32);
    raster_.setDevicePixelRatio(pixelRatio);
    if (raster_.isNull())
        return;
    const auto bounds = visibleBounds();
    const double low = std::log(bounds[0]), high = std::log(bounds[1]);
    if (radial_) {
        const QRectF plot(0, 0, size.width(), size.height());
        for (int row = 0; row < size.height(); ++row) {
            auto *pixels = reinterpret_cast<QRgb *>(raster_.scanLine(row));
            for (int column = 0; column < size.width(); ++column) {
                const auto coordinate = radialCoordinates(QPointF(column + .5, row + .5), plot);
                const auto value =
                    coordinate ? interpolatedLevel(coordinate->x(), coordinate->y()) : std::nullopt;
                pixels[column] =
                    value ? levelColor(*value, decibels_).rgb() : palette().base().color().rgb();
            }
        }
        rasterDirty_ = false;
        return;
    }
    std::vector<std::optional<Interval>> frequencyIntervals;
    frequencyIntervals.reserve(static_cast<std::size_t>(size.width()));
    for (int column = 0; column < size.width(); ++column) {
        const double fraction = (column + .5) / size.width();
        const double value = low + fraction * (high - low);
        frequencyIntervals.push_back(value < logarithmicFrequencies_.front() ||
                                             value > logarithmicFrequencies_.back()
                                         ? std::nullopt
                                         : std::optional{interval(logarithmicFrequencies_, value)});
    }
    const auto missing = palette().alternateBase().color().rgb();
    for (int row = 0; row < size.height(); ++row) {
        auto *pixels = reinterpret_cast<QRgb *>(raster_.scanLine(row));
        const double valueAngle = bounds[3] - (row + .5) / size.height() * (bounds[3] - bounds[2]);
        const auto angle = interval(angles_, valueAngle);
        for (int column = 0; column < size.width(); ++column) {
            const auto value =
                frequencyIntervals[column] && valueAngle >= angles_.front() &&
                        valueAngle <= angles_.back()
                    ? sample(levels_, angles_.size(), *frequencyIntervals[column], angle)
                    : std::nullopt;
            pixels[column] = value ? levelColor(*value, decibels_).rgb() : missing;
        }
    }
    rasterDirty_ = false;
}

void HeatmapWidget::paintEvent(QPaintEvent *) {
    QPainter painter(this);
    painter.fillRect(rect(), palette().base());
    painter.setPen(palette().text().color());
    const bool single = frequencies_.size() == 1;
    painter.drawText(QRect(4, 3, width() - 8, 20), Qt::AlignCenter,
                     single ? title_ + " · " + frequencyLabel(frequencies_.front()) + " Hz"
                            : title_);
    if (levels_.empty())
        return;
    if (radial_) {
        paintRadial(painter);
        return;
    }
    const auto plot = plotRectangle();
    const auto bounds = visibleBounds();
    const double logMinimum = std::log(bounds[0]), logMaximum = std::log(bounds[1]);
    const auto y = [&](double angle) {
        return plot.bottom() - (angle - bounds[2]) / (bounds[3] - bounds[2]) * plot.height();
    };
    const auto x = [&](double frequency) {
        return plot.left() +
               (std::log(frequency) - logMinimum) / (logMaximum - logMinimum) * plot.width();
    };
    const auto levelX = [&](double level) {
        return plot.left() +
               (std::clamp(level, decibels_.minimum, decibels_.maximum) - decibels_.minimum) /
                   (decibels_.maximum - decibels_.minimum) * plot.width();
    };
    if (!single) {
        const double ratio = devicePixelRatioF();
        updateRaster({std::max(1, static_cast<int>(std::ceil(plot.width() * ratio))),
                      std::max(1, static_cast<int>(std::ceil(plot.height() * ratio)))},
                     ratio);
        painter.drawImage(plot, raster_);
    }
    QFont ticks = font();
    ticks.setPointSizeF(std::max(7.0, font().pointSizeF() - 1));
    painter.setFont(ticks);
    const auto grid = QColor(115, 125, 133, single ? 65 : 100);
    for (const auto angle : plotTicks(bounds[2], bounds[3], scales_.y.tickStep,
                                      static_cast<unsigned>(plot.height()))) {
        painter.setPen(QPen(grid, 1, Qt::DotLine));
        painter.drawLine(QPointF(plot.left(), y(angle)), QPointF(plot.right(), y(angle)));
        painter.setPen(palette().text().color());
        painter.drawText(QRectF(13, y(angle) - 8, 39, 16), Qt::AlignRight | Qt::AlignVCenter,
                         QString::number(angle, 'g', 4));
    }
    if (single) {
        for (const double level :
             plotTicks(decibels_.minimum, decibels_.maximum, decibels_.tickStep,
                       static_cast<unsigned>(plot.width()))) {
            const double px = levelX(level);
            painter.setPen(QPen(grid, 1, Qt::DotLine));
            painter.drawLine(QPointF(px, plot.top()), QPointF(px, plot.bottom()));
            painter.setPen(palette().text().color());
            painter.drawText(QRectF(px - 22, plot.bottom() + 4, 44, 17), Qt::AlignCenter,
                             QString::number(level));
        }
        QPainterPath path;
        bool started = false;
        for (std::size_t index = 0; index < angles_.size(); ++index) {
            if (std::isnan(levels_[index])) {
                started = false;
                continue;
            }
            const QPointF point(levelX(levels_[index]), y(angles_[index]));
            if (started)
                path.lineTo(point);
            else
                path.moveTo(point);
            started = true;
        }
        painter.setRenderHint(QPainter::Antialiasing);
        painter.save();
        painter.setClipRect(plot);
        painter.setPen(QPen(QColor(34, 104, 160), 1.8));
        painter.drawPath(path);
        painter.restore();
    } else {
        const auto frequencyMarks = scales_.x.tickStep
                                        ? plotTicks(bounds[0], bounds[1], scales_.x.tickStep,
                                                    static_cast<unsigned>(plot.width()))
                                        : frequencyTicks(bounds[0], bounds[1], plot.width());
        for (const auto frequency : frequencyMarks) {
            const double px = x(frequency);
            painter.setPen(QPen(grid, 1, Qt::DotLine));
            painter.drawLine(QPointF(px, plot.top()), QPointF(px, plot.bottom()));
            painter.setPen(palette().text().color());
            painter.drawText(QRectF(px - 31, plot.bottom() + 4, 62, 17), Qt::AlignCenter,
                             frequencyLabel(frequency));
        }
        const QRectF bar(plot.right() + 13, plot.top(), 13, plot.height());
        for (int row = 0; row < static_cast<int>(std::ceil(bar.height())); ++row)
            painter.fillRect(
                QRectF(bar.left(), bar.top() + row, bar.width(), 1.1),
                levelColor(decibels_.maximum -
                               row / bar.height() * (decibels_.maximum - decibels_.minimum),
                           decibels_));
        painter.setPen(palette().text().color());
        painter.setBrush(Qt::NoBrush);
        painter.drawRect(bar);
        double previousLabel = std::numeric_limits<double>::infinity();
        for (const double level :
             plotTicks(decibels_.minimum, decibels_.maximum, decibels_.tickStep,
                       static_cast<unsigned>(plot.height()))) {
            const double py = bar.top() + (decibels_.maximum - level) /
                                              (decibels_.maximum - decibels_.minimum) *
                                              bar.height();
            painter.drawLine(QPointF(bar.right(), py), QPointF(bar.right() + 3, py));
            if (level == decibels_.minimum || level == decibels_.maximum ||
                (previousLabel - py >= 16 && py > bar.top() + 16)) {
                painter.drawText(QRectF(bar.right() + 5, py - 8, 45, 16),
                                 Qt::AlignLeft | Qt::AlignVCenter, QString::number(level, 'g', 4));
                previousLabel = py;
            }
        }
        painter.drawText(QRectF(bar.left() - 2, bar.top() - 18, 44, 16), Qt::AlignLeft, "dB");
    }
    painter.setPen(QColor(112, 120, 128));
    painter.setBrush(Qt::NoBrush);
    painter.drawRect(plot);
    painter.setPen(palette().text().color());
    painter.drawText(QRectF(plot.left(), height() - 43, plot.width(), 17), Qt::AlignCenter,
                     single ? tr("Relativer Pegel [dB]") : tr("Frequenz [Hz]"));
    painter.save();
    painter.translate(9, plot.center().y());
    painter.rotate(-90);
    painter.drawText(QRectF(-plot.height() / 2, -8, plot.height(), 16), Qt::AlignCenter,
                     tr("Winkel [°]"));
    painter.restore();
    if (cursor_) {
        const double angle =
            bounds[3] - (cursor_->y() - plot.top()) / plot.height() * (bounds[3] - bounds[2]);
        const double logarithm = single ? logarithmicFrequencies_.front()
                                        : logMinimum + (cursor_->x() - plot.left()) / plot.width() *
                                                           (logMaximum - logMinimum);
        const auto value = interpolatedLevel(logarithm, angle);
        painter.setPen(QPen(QColor(25, 25, 25, 190), 1, Qt::DashLine));
        painter.drawLine(QPointF(plot.left(), cursor_->y()), QPointF(plot.right(), cursor_->y()));
        if (!single)
            painter.drawLine(QPointF(cursor_->x(), plot.top()),
                             QPointF(cursor_->x(), plot.bottom()));
        else if (value)
            painter.drawEllipse(QPointF(levelX(*value), cursor_->y()), 3, 3);
        painter.setPen(palette().text().color());
        painter.drawText(QRectF(4, height() - 22, width() - 8, 18), Qt::AlignCenter,
                         frequencyLabel(std::exp(logarithm)) + " Hz  ·  " +
                             QString::number(angle, 'f', 1) + "°  ·  " +
                             (value ? QString::number(*value, 'f', 1) + " dB" : tr("keine Daten")));
    }
}

void HeatmapWidget::updateCursor(QPointF position) {
    if (levels_.empty())
        return;
    const auto plot = plotRectangle();
    if (dragging_)
        cursor_ = QPointF(std::clamp(position.x(), plot.left(), plot.right()),
                          std::clamp(position.y(), plot.top(), plot.bottom()));
    else if (plot.contains(position))
        cursor_ = position;
    else
        cursor_.reset();
    update();
}
void HeatmapWidget::mouseMoveEvent(QMouseEvent *event) {
    updateCursor(event->position());
}
void HeatmapWidget::mousePressEvent(QMouseEvent *event) {
    if (event->button() == Qt::LeftButton && !levels_.empty() &&
        plotRectangle().contains(event->position())) {
        dragging_ = true;
        updateCursor(event->position());
        event->accept();
    } else
        QWidget::mousePressEvent(event);
}
void HeatmapWidget::mouseReleaseEvent(QMouseEvent *event) {
    if (event->button() == Qt::LeftButton) {
        dragging_ = false;
        updateCursor(event->position());
    }
    QWidget::mouseReleaseEvent(event);
}
void HeatmapWidget::leaveEvent(QEvent *event) {
    if (!dragging_) {
        cursor_.reset();
        update();
    }
    QWidget::leaveEvent(event);
}
} // namespace hornsculpt
