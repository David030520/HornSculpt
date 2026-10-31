#pragma once
#include <QString>
#include <hornsculpt/plot_scale.hpp>
#include <optional>
class QWidget;
namespace hornsculpt {
std::optional<PlotAxisScales> editAxisScales(QWidget &, PlotAxisScales, const QString &xLabel,
                                             const QString &yLabel, bool logarithmicX = false,
                                             bool decibelY = false);
struct PolarScales {
    PlotDbScale decibels;
    double angularStep;
};
std::optional<PolarScales> editPolarScales(QWidget &, PolarScales);
struct HeatmapScales {
    PlotAxisScales axes;
    PlotDbScale decibels;
};
std::optional<HeatmapScales> editHeatmapScales(QWidget &, HeatmapScales);
PlotAxisScale readAxisScale(const QWidget &, const QString &axis, PlotAxisScale fallback);
void saveAxisScale(const QWidget &, const QString &axis, PlotAxisScale);
double readAngularStep(const QWidget &, double fallback);
void saveAngularStep(const QWidget &, double);
} // namespace hornsculpt
