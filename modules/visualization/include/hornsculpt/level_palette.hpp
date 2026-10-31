#pragma once

#include <QColor>
#include <algorithm>
#include <array>
#include <cmath>
#include <hornsculpt/plot_scale.hpp>

namespace hornsculpt {
// High levels are dark red; low levels are dark blue. Bands between these
// reference stops are interpolated by band index, never by sample position.
inline const std::array<QColor, 10> &levelPalette() {
    static const std::array<QColor, 10> colors{
        QColor(153, 0, 0),   QColor(224, 24, 0),   QColor(255, 105, 0), QColor(255, 211, 0),
        QColor(205, 255, 0), QColor(48, 255, 171), QColor(0, 216, 255), QColor(0, 132, 255),
        QColor(0, 48, 230),  QColor(0, 0, 132)};
    return colors;
}
inline QColor steppedLevelColor(double level, PlotDbScale scale) {
    // Infinite dB is meaningful for a zero level and clamps to the scale;
    // NaN is missing data and must never be converted to a colour index.
    if (std::isnan(level) || !std::isfinite(scale.minimum) || !std::isfinite(scale.maximum) ||
        !std::isfinite(scale.tickStep) || scale.maximum <= scale.minimum || scale.tickStep <= 0)
        return {};
    if (!std::isfinite(scale.maximum - scale.minimum))
        return {};
    const auto &colors = levelPalette();
    const double count = std::ceil((scale.maximum - scale.minimum) / scale.tickStep);
    const double band =
        std::clamp(std::floor((scale.maximum - std::clamp(level, scale.minimum, scale.maximum)) /
                              scale.tickStep),
                   0., std::max(0., count - 1));
    const double position =
        std::isfinite(count) ? (count > 1 ? band / (count - 1) * 9 : 0)
                             : (scale.maximum - std::clamp(level, scale.minimum, scale.maximum)) /
                                   (scale.maximum - scale.minimum) * 9;
    const auto index = static_cast<std::size_t>(std::clamp(std::floor(position), 0., 9.));
    const auto next = std::min<std::size_t>(9, index + 1);
    const double fraction = position - index;
    return QColor::fromRgbF(
        colors[index].redF() * (1 - fraction) + colors[next].redF() * fraction,
        colors[index].greenF() * (1 - fraction) + colors[next].greenF() * fraction,
        colors[index].blueF() * (1 - fraction) + colors[next].blueF() * fraction);
}
} // namespace hornsculpt
