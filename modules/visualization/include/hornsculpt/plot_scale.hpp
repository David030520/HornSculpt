#pragma once

#include <optional>
#include <vector>

namespace hornsculpt {
struct PlotAxisScale {
    std::optional<double> minimum, maximum, tickStep;
    bool operator==(const PlotAxisScale &) const = default;
};
struct PlotAxisScales {
    PlotAxisScale x, y;
    bool operator==(const PlotAxisScales &) const = default;
};
struct PlotDbScale {
    double minimum{-30}, maximum{3}, tickStep{3};
    bool operator==(const PlotDbScale &) const = default;
};
void validatePlotAxisScale(const PlotAxisScale &scale, bool logarithmic = false);
void validatePlotDbScale(const PlotDbScale &scale);
// Dense tick sequences are sampled at integer multiples of their requested
// step to avoid subpixel grids; the stored scale and visible bounds are exact.
std::vector<double> plotTicks(double minimum, double maximum, std::optional<double> step = {},
                              unsigned pixelSpan = 600);
} // namespace hornsculpt
