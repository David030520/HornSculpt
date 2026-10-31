#pragma once

#include <QColor>
#include <QString>
#include <QWidget>
#include <functional>
#include <hornsculpt/plot_scale.hpp>
#include <vector>

namespace hornsculpt {
struct AngularLevel {
    double angleRadians{};
    double levelDecibels{};
};
struct PolarSeries {
    QString label;
    QColor color;
    std::vector<AngularLevel> levels;
};

class PolarPlotWidget final : public QWidget {
  public:
    explicit PolarPlotWidget(QString title, QColor color, QWidget *parent = nullptr);
    void setLevels(std::vector<AngularLevel> levels);
    void setSeries(std::vector<PolarSeries> series);
    const std::vector<PolarSeries> &series() const noexcept { return series_; }
    void setTitle(QString title);
    void setDecibelScale(PlotDbScale scale);
    PlotDbScale decibelScale() const {
        return decibels_;
    }
    void setAngularTickStep(double degrees);
    double angularTickStep() const {
        return angularStep_;
    }
    std::function<void(PlotDbScale)> onDecibelScaleChanged;

  protected:
    void paintEvent(QPaintEvent *) override;

  private:
    QString title_;
    QColor color_;
    std::vector<PolarSeries> series_;
    bool fullCircle_{true};
    PlotDbScale decibels_;
    double angularStep_{30};
};
} // namespace hornsculpt
