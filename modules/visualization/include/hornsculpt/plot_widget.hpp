#pragma once

#include <QColor>
#include <QPointF>
#include <QString>
#include <QWidget>
#include <functional>
#include <hornsculpt/plot_scale.hpp>
#include <optional>
#include <vector>

namespace hornsculpt {

struct PlotSeries {
    QString label;
    QColor color;
    std::vector<QPointF> points;
};
struct PlotReferenceBand {
    double minimum{};
    double maximum{};
    QColor color;
    QString label;
};

class PlotWidget final : public QWidget {
  public:
    PlotWidget(QString title, QString xLabel, QString yLabel, QWidget *parent = nullptr);
    void setSeries(std::vector<PlotSeries> series);
    const std::vector<PlotSeries> &series() const noexcept {
        return series_;
    }
    void setFrequencyAxis();
    void setReferenceBand(std::optional<PlotReferenceBand> band);
    void setTitle(QString title);
    void setAxisLabels(QString xLabel, QString yLabel);
    void setAxisScales(PlotAxisScale x, PlotAxisScale y);
    PlotAxisScales axisScales() const {
        return scales_;
    }
    void setDecibelScale(PlotDbScale scale);
    PlotDbScale decibelScale() const;
    std::function<void(PlotDbScale)> onDecibelScaleChanged;

  protected:
    void paintEvent(QPaintEvent *) override;
    void changeEvent(QEvent *) override;
    void resizeEvent(QResizeEvent *) override;

  private:
    void updateLegendMinimumHeight();
    QString title_;
    QString xLabel_;
    QString yLabel_;
    std::vector<PlotSeries> series_;
    bool frequencyAxis_{};
    std::optional<PlotReferenceBand> referenceBand_;
    PlotAxisScales scales_;
    bool decibelAxis_{};
};

} // namespace hornsculpt
