#pragma once

#include <QImage>
#include <QPointF>
#include <QString>
#include <QWidget>
#include <array>
#include <functional>
#include <hornsculpt/plot_scale.hpp>
#include <optional>
#include <vector>

namespace hornsculpt {
class HeatmapWidget final : public QWidget {
  public:
    explicit HeatmapWidget(QString title, QWidget *parent = nullptr);
    // Levels are row-major: each ascending frequency row contains all ascending angles.
    void setLevels(std::vector<double> anglesDegrees, std::vector<double> frequencies,
                   std::vector<double> levelsDecibels);
    const std::vector<double> &anglesDegrees() const noexcept { return angles_; }
    const std::vector<double> &frequencies() const noexcept { return frequencies_; }
    const std::vector<double> &levelsDecibels() const noexcept { return levels_; }
    void setTitle(QString title);
    void setRadialLayout(bool enabled);
    void setAxisScales(PlotAxisScale frequency, PlotAxisScale angle);
    PlotAxisScales axisScales() const {
        return scales_;
    }
    void setDecibelScale(PlotDbScale scale);
    PlotDbScale decibelScale() const {
        return decibels_;
    }
    std::function<void(PlotDbScale)> onDecibelScaleChanged;

  protected:
    void paintEvent(QPaintEvent *) override;
    void mouseMoveEvent(QMouseEvent *) override;
    void mousePressEvent(QMouseEvent *) override;
    void mouseReleaseEvent(QMouseEvent *) override;
    void leaveEvent(QEvent *) override;

  private:
    static QColor levelColor(double level, PlotDbScale scale);
    void paintRadial(QPainter &painter);
    std::optional<QPointF> radialCoordinates(QPointF position, QRectF plot) const;
    bool radial_{};
    QRectF plotRectangle() const;
    std::optional<double> interpolatedLevel(double logarithmicFrequency, double angle) const;
    void updateRaster(QSize size, double pixelRatio);
    void updateCursor(QPointF position);
    std::array<double, 4> visibleBounds() const;
    QString title_;
    std::vector<double> angles_;
    std::vector<double> frequencies_;
    std::vector<double> logarithmicFrequencies_;
    std::vector<double> levels_;
    QImage raster_;
    bool rasterDirty_{true};
    bool dragging_{};
    std::optional<QPointF> cursor_;
    PlotAxisScales scales_;
    PlotDbScale decibels_;
};
} // namespace hornsculpt
