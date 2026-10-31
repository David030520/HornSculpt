#include "chart_export.hpp"
#include "chart_scaling.hpp"
#include "plot_legend_layout.hpp"
#include <QEvent>
#include <QPainter>
#include <QPainterPath>
#include <QResizeEvent>
#include <algorithm>
#include <cmath>
#include <hornsculpt/plot_axes.hpp>
#include <hornsculpt/plot_widget.hpp>
#include <limits>
#include <stdexcept>

namespace hornsculpt {
namespace {
QString axisNumber(double value) {
    QString text = QString::number(value, 'g', 3);
    text.replace("e+0", "e");
    text.replace("e-0", "e−");
    text.replace("e+", "e");
    return text;
}
} // namespace

PlotWidget::PlotWidget(QString title, QString xLabel, QString yLabel, QWidget *parent)
    : QWidget(parent), title_(std::move(title)), xLabel_(std::move(xLabel)),
      yLabel_(std::move(yLabel)) {
    setMinimumSize(260, 155);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    enableChartExport(*this, [this] {
        if (auto scales =
                editAxisScales(*this, scales_, xLabel_, yLabel_, frequencyAxis_, decibelAxis_)) {
            setAxisScales(scales->x, scales->y);
            if (decibelAxis_ && scales_.y.minimum && scales_.y.tickStep && onDecibelScaleChanged)
                onDecibelScaleChanged(decibelScale());
        }
    });
    connect(this, &QObject::objectNameChanged, this, [this] {
        auto x = readAxisScale(*this, "x", scales_.x);
        try {
            validatePlotAxisScale(x, frequencyAxis_);
        } catch (const std::invalid_argument &) {
            x = scales_.x;
        }
        scales_ = {x, readAxisScale(*this, "y", scales_.y)};
        update();
    });
}
void PlotWidget::setTitle(QString title) {
    title_ = std::move(title);
    update();
}
void PlotWidget::setAxisLabels(QString xLabel, QString yLabel) {
    xLabel_ = std::move(xLabel);
    yLabel_ = std::move(yLabel);
    update();
}
void PlotWidget::setAxisScales(PlotAxisScale x, PlotAxisScale y) {
    validatePlotAxisScale(x, frequencyAxis_);
    validatePlotAxisScale(y);
    scales_ = {x, y};
    saveAxisScale(*this, "x", x);
    saveAxisScale(*this, "y", y);
    update();
}
void PlotWidget::setDecibelScale(PlotDbScale s) {
    validatePlotDbScale(s);
    decibelAxis_ = true;
    setAxisScales(scales_.x, {s.minimum, s.maximum, s.tickStep});
}
PlotDbScale PlotWidget::decibelScale() const {
    return {scales_.y.minimum.value_or(-30), scales_.y.maximum.value_or(3),
            scales_.y.tickStep.value_or(3)};
}
void PlotWidget::setSeries(std::vector<PlotSeries> series) {
    series_ = std::move(series);
    updateLegendMinimumHeight();
    update();
}
void PlotWidget::updateLegendMinimumHeight() {
    std::vector<QString> labels;
    labels.reserve(series_.size());
    for (const auto &series : series_)
        labels.push_back(series.label);
    QFont legendFont = font();
    legendFont.setPointSizeF(std::max(7.0, font().pointSizeF() - 1));
    const auto legend = chart::legendLayout(labels, QFontMetricsF(legendFont),
                                             std::max(260, width()) - 16.);
    // Reserve an actual data rectangle even at the minimum supported width.
    setMinimumHeight(std::max(155, int(std::ceil(28. + 80. + 40. + legend.height + 8.))));
    updateGeometry();
}
void PlotWidget::changeEvent(QEvent *event) {
    QWidget::changeEvent(event);
    if (event->type() == QEvent::FontChange)
        updateLegendMinimumHeight();
}
void PlotWidget::resizeEvent(QResizeEvent *event) {
    QWidget::resizeEvent(event);
    if (event->size().width() != event->oldSize().width())
        updateLegendMinimumHeight();
}
void PlotWidget::setFrequencyAxis() {
    validatePlotAxisScale(scales_.x, true);
    frequencyAxis_ = true;
    update();
}
void PlotWidget::setReferenceBand(std::optional<PlotReferenceBand> band) {
    if (band && (!std::isfinite(band->minimum) || !std::isfinite(band->maximum) ||
                 band->minimum > band->maximum))
        throw std::invalid_argument("Invalid plot reference interval.");
    referenceBand_ = std::move(band);
    update();
}

void PlotWidget::paintEvent(QPaintEvent *) {
    QPainter painter(this);
    painter.fillRect(rect(), palette().base());
    painter.setRenderHint(QPainter::Antialiasing);
    QFont tickFont = font();
    tickFont.setPointSizeF(std::max(7.0, font().pointSizeF() - 1));
    std::vector<QString> labels;
    labels.reserve(series_.size());
    for (const auto &series : series_)
        labels.push_back(series.label);
    const auto legend = chart::legendLayout(labels, QFontMetricsF(tickFont), width() - 16.);
    const double legendSpace = legend.height > 0 ? legend.height + 8. : 0.;
    const QRectF area(65, 28, width() - 81, std::max(1., height() - 67. - legendSpace));
    painter.setPen(palette().text().color());
    painter.drawText(QRectF(5, 3, width() - 10, 20), Qt::AlignCenter, title_);
    if (series_.empty())
        return;
    double xmin = std::numeric_limits<double>::infinity(), xmax = -xmin, ymin = xmin, ymax = xmax;
    for (const auto &s : series_)
        for (const auto &point : s.points) {
            if (!std::isfinite(point.x()) || std::isnan(point.y()) ||
                (!std::isfinite(point.y()) && !scales_.y.minimum) ||
                (frequencyAxis_ && point.x() <= 0))
                continue;
            const double coordinate = frequencyAxis_ ? std::log10(point.x()) : point.x();
            xmin = std::min(xmin, coordinate);
            xmax = std::max(xmax, coordinate);
            if (std::isfinite(point.y())) {
                ymin = std::min(ymin, point.y());
                ymax = std::max(ymax, point.y());
            }
        }
    if (!std::isfinite(xmin))
        return;
    if (scales_.y.minimum) {
        ymin = *scales_.y.minimum;
        ymax = *scales_.y.maximum;
    }
    if (referenceBand_) {
        ymin = std::min(ymin, referenceBand_->minimum);
        ymax = std::max(ymax, referenceBand_->maximum);
    }
    if (xmax <= xmin) {
        const double pad = frequencyAxis_ ? .05 : std::max(.5, std::abs(xmin) * .05);
        xmin -= pad;
        xmax += pad;
    }
    if (ymax <= ymin) {
        ymin -= .5;
        ymax += .5;
    } else {
        const double pad = (ymax - ymin) * .04;
        ymin -= pad;
        ymax += pad;
    }
    if (scales_.x.minimum) {
        xmin = frequencyAxis_ ? std::log10(*scales_.x.minimum) : *scales_.x.minimum;
        xmax = frequencyAxis_ ? std::log10(*scales_.x.maximum) : *scales_.x.maximum;
    }
    if (scales_.y.minimum) {
        ymin = *scales_.y.minimum;
        ymax = *scales_.y.maximum;
    }
    const auto map = [&](QPointF p) {
        return QPointF(area.left() + ((frequencyAxis_ ? std::log10(p.x()) : p.x()) - xmin) /
                                         (xmax - xmin) * area.width(),
                       area.bottom() -
                           ((std::isfinite(p.y()) ? p.y() : std::clamp(p.y(), ymin, ymax)) - ymin) /
                               (ymax - ymin) * area.height());
    };
    if (referenceBand_) {
        painter.save();
        painter.setClipRect(area);
        const double top =
            area.bottom() - (referenceBand_->maximum - ymin) / (ymax - ymin) * area.height();
        const double bottom =
            area.bottom() - (referenceBand_->minimum - ymin) / (ymax - ymin) * area.height();
        auto fill = referenceBand_->color;
        fill.setAlpha(30);
        painter.fillRect(QRectF(area.left(), top, area.width(), bottom - top), fill);
        painter.setPen(QPen(referenceBand_->color, 1, Qt::DashLine));
        painter.drawLine(QPointF(area.left(), top), QPointF(area.right(), top));
        painter.drawLine(QPointF(area.left(), bottom), QPointF(area.right(), bottom));
        painter.drawText(QRectF(area.left() + 5, area.bottom() - 20, area.width() - 10, 18),
                         Qt::AlignRight, referenceBand_->label);
        painter.restore();
    }
    painter.setFont(tickFont);
    double previousYLabel = std::numeric_limits<double>::infinity();
    for (const double tick :
         plotTicks(ymin, ymax, scales_.y.tickStep, static_cast<unsigned>(area.height()))) {
        const double f = (tick - ymin) / (ymax - ymin);
        const double y = area.bottom() - area.height() * f;
        painter.setPen(QPen(QColor(219, 223, 227), 1));
        painter.drawLine(QPointF(area.left(), y), QPointF(area.right(), y));
        if (tick == ymin || tick == ymax || (previousYLabel - y >= 16 && y > area.top() + 16)) {
            painter.setPen(palette().text().color());
            painter.drawText(QRectF(18, y - 8, 42, 16), Qt::AlignRight | Qt::AlignVCenter,
                             axisNumber(tick));
            previousYLabel = y;
        }
    }
    double previousXLabel = -std::numeric_limits<double>::infinity();
    const auto drawXTick = [&](double position, const QString &label) {
        const double x = area.left() + (position - xmin) / (xmax - xmin) * area.width();
        painter.setPen(QPen(QColor(219, 223, 227), 1));
        painter.drawLine(QPointF(x, area.top()), QPointF(x, area.bottom()));
        if (position == xmin || position == xmax ||
            (x - previousXLabel >= 65 && x < area.right() - 65)) {
            painter.setPen(palette().text().color());
            painter.drawText(QRectF(x - 35, area.bottom() + 4, 70, 16), Qt::AlignHCenter, label);
            previousXLabel = x;
        }
    };
    if (frequencyAxis_) {
        const auto ticks =
            scales_.x.tickStep
                ? plotTicks(std::pow(10., xmin), std::pow(10., xmax), scales_.x.tickStep,
                            static_cast<unsigned>(area.width()))
                : frequencyTicks(std::pow(10., xmin), std::pow(10., xmax), area.width());
        for (const double tick : ticks)
            drawXTick(std::log10(tick), frequencyLabel(tick));
    } else {
        for (const double tick :
             plotTicks(xmin, xmax, scales_.x.tickStep, static_cast<unsigned>(area.width()))) {
            drawXTick(tick, axisNumber(tick));
        }
    }
    painter.setPen(QPen(QColor(112, 120, 128), 1));
    painter.drawRect(area);
    painter.save();
    painter.setClipRect(area.adjusted(-1, -1, 1, 1));
    for (const auto &s : series_) {
        QPainterPath path;
        bool started = false;
        std::size_t segmentPointCount = 0;
        QPointF segmentPoint;
        std::vector<QPointF> isolatedPoints;
        for (const auto &point : s.points) {
            if (!std::isfinite(point.x()) || std::isnan(point.y()) ||
                (!std::isfinite(point.y()) && !scales_.y.minimum) ||
                (frequencyAxis_ && point.x() <= 0)) {
                if (segmentPointCount == 1)
                    isolatedPoints.push_back(segmentPoint);
                segmentPointCount = 0;
                started = false;
                continue;
            }
            ++segmentPointCount;
            segmentPoint = map(point);
            if (!started) {
                path.moveTo(map(point));
                started = true;
            } else
                path.lineTo(map(point));
        }
        painter.setPen(QPen(s.color, 1.8));
        painter.drawPath(path);
        if (segmentPointCount == 1)
            isolatedPoints.push_back(segmentPoint);
        if (!isolatedPoints.empty()) {
            painter.setBrush(s.color);
            for (const auto &point : isolatedPoints)
                painter.drawEllipse(point, 3.0, 3.0);
            painter.setBrush(Qt::NoBrush);
        }
    }
    painter.restore();
    painter.setPen(palette().text().color());
    painter.drawText(QRectF(area.left(), area.bottom() + 20, area.width(), 17), Qt::AlignCenter,
                     xLabel_);
    painter.save();
    painter.translate(8, area.center().y());
    painter.rotate(-90);
    painter.drawText(QRectF(-area.height() * .5, -8, area.height(), 16), Qt::AlignCenter, yLabel_);
    painter.restore();
    const QPointF origin(8., area.bottom() + 40.);
    for (const auto &item : legend.items) {
        const auto &series = series_[item.seriesIndex];
        const QRectF swatch = item.swatch.translated(origin);
        const QRectF text = item.text.translated(origin);
        painter.setPen(QPen(series.color, 2));
        painter.drawLine(QPointF(swatch.left(), swatch.center().y()),
                         QPointF(swatch.right(), swatch.center().y()));
        painter.setPen(palette().text().color());
        painter.drawText(text, Qt::AlignLeft | Qt::AlignTop | Qt::TextWordWrap |
                                   Qt::TextWrapAnywhere, series.label);
    }
}

} // namespace hornsculpt
