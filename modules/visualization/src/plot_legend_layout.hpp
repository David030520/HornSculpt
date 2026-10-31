#pragma once

#include <QFontMetricsF>
#include <QRectF>
#include <QString>
#include <algorithm>
#include <vector>

namespace hornsculpt::chart {
struct LegendItem {
    std::size_t seriesIndex;
    QRectF swatch;
    QRectF text;
};
struct LegendLayout {
    std::vector<LegendItem> items;
    double height{};
};
// Positions are relative to the legend's top left. A long label wraps rather
// than disappearing or covering the data area.
inline LegendLayout legendLayout(const std::vector<QString> &labels, const QFontMetricsF &font,
                                 double width) {
    LegendLayout result;
    if (width <= 28.)
        return result;
    constexpr double swatchWidth = 14., textGap = 5., itemGap = 12., rowGap = 3.;
    const double maxText = width - swatchWidth - textGap;
    double x = 0., y = 0., rowHeight = 0.;
    for (std::size_t i = 0; i < labels.size(); ++i) {
        if (labels[i].isEmpty())
            continue;
        const double textWidth = std::min(maxText, std::max(1., font.horizontalAdvance(labels[i])));
        const auto bounds = font.boundingRect(QRectF(0, 0, textWidth, 1000000.),
                                              Qt::TextWordWrap | Qt::TextWrapAnywhere, labels[i]);
        const double textHeight = std::max(font.height(), bounds.height());
        const double itemWidth = swatchWidth + textGap + textWidth;
        if (x > 0. && x + itemWidth > width) {
            y += rowHeight + rowGap;
            x = 0.;
            rowHeight = 0.;
        }
        result.items.push_back({i, QRectF(x, y, swatchWidth, font.height()),
                                QRectF(x + swatchWidth + textGap, y, textWidth, textHeight)});
        rowHeight = std::max(rowHeight, textHeight);
        x += itemWidth + itemGap;
    }
    result.height = result.items.empty() ? 0. : y + rowHeight;
    return result;
}
} // namespace hornsculpt::chart
