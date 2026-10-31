#pragma once

#include <QFontMetrics>
#include <QPointF>
#include <QRectF>
#include <QSizeF>
#include <algorithm>
#include <array>
#include <optional>
#include <vector>

namespace hornsculpt::overlay {
struct GizmoLayout {
    QPointF origin;
    double arm{};
    QRectF reserved;
};
inline std::optional<GizmoLayout> gizmo(const QSizeF &viewport, const QFontMetrics &font,
                                        bool fieldContext, bool heading) {
    const double top = fieldContext ? 32. + font.height() + 6. : heading ? 31. : 4.;
    const double byHeight = (viewport.height() - top - font.height() - 12.) / 2.;
    const double byWidth = (viewport.width() - 24. - font.horizontalAdvance("W")) / 2.;
    const double arm = std::min({28., byHeight, byWidth});
    if (arm < 6.)
        return {};
    const QPointF origin(12. + arm, viewport.height() - 5. - arm);
    const double reservedTop = origin.y() - arm - font.height() - 3.;
    const QRectF reserved(0., reservedTop,
                          std::min(viewport.width(), 18. + 2. * arm + font.horizontalAdvance("W")),
                          viewport.height() - reservedTop);
    return GizmoLayout{origin, arm, reserved};
}
inline QRectF fieldLegend(const QSizeF &viewport, const QFontMetrics &font, const QString &title,
                          const QString &minimum, const QString &maximum) {
    const double top = 35.;
    const double textHeight = font.height();
    const double available = viewport.height() - top - textHeight - 18.;
    if (available < std::max(24., 2. * textHeight) || viewport.width() < 70.)
        return {};
    const double width = std::min(
        viewport.width() - 12.,
        std::max(double(font.horizontalAdvance(title)) + 12.,
                 29. + std::max(font.horizontalAdvance(minimum), font.horizontalAdvance(maximum))));
    const double barHeight = std::min(150., available);
    return {viewport.width() - width - 6., top, width, textHeight + barHeight + 12.};
}
inline std::optional<QRectF> label(const QPointF &anchor, const QSizeF &size,
                                   const QRectF &viewport, const std::vector<QRectF> &occupied) {
    if (size.width() > viewport.width() || size.height() > viewport.height())
        return {};
    const std::array<QPointF, 8> offsets{{{6, -size.height() - 4},
                                          {6, 4},
                                          {-size.width() - 6, -size.height() - 4},
                                          {-size.width() - 6, 4},
                                          {-size.width() / 2, -size.height() - 10},
                                          {-size.width() / 2, 10},
                                          {12, -size.height() / 2},
                                          {-size.width() - 12, -size.height() / 2}}};
    for (const auto &offset : offsets) {
        QRectF r(anchor + offset, size);
        r.moveLeft(std::clamp(r.left(), viewport.left(), viewport.right() - r.width()));
        r.moveTop(std::clamp(r.top(), viewport.top(), viewport.bottom() - r.height()));
        if (std::none_of(occupied.begin(), occupied.end(),
                         [&](const QRectF &other) { return r.intersects(other); }))
            return r;
    }
    return {};
}
} // namespace hornsculpt::overlay
