#include "../src/overlay_layout.hpp"
#include <QApplication>
#include <QFontMetrics>
#include <stdexcept>
#include <vector>
using namespace hornsculpt;
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
int main(int argc, char **argv) {
    QApplication app(argc, argv);
    const QFontMetrics font(app.font());
    const QString title = QStringLiteral("Pa RMS · 10 Stufen");
    for (const QSize size : {QSize(853, 130), QSize(853, 160), QSize(320, 100)}) {
        const auto r = overlay::fieldLegend(size, font, title, QStringLiteral("-123.45"),
                                            QStringLiteral("123.45"));
        require(!r.isEmpty() && QRectF(QPointF(0, 0), size).contains(r), "Legend exceeds viewport");
        require(r.width() >= font.horizontalAdvance(title) + 10, "Legend title does not fit");
        const auto gizmo = overlay::gizmo(size, font, true, true);
        require(gizmo.has_value() && gizmo->arm >= 6 &&
                    gizmo->reserved.top() >= 32 + font.height() + 6 &&
                    QRectF(QPointF(0, 0), size).contains(gizmo->reserved),
                "XYZ indicator overlaps the field context or leaves the viewport");
    }
    require(overlay::fieldLegend(QSize(150, 40), font, title, "-1", "1").isEmpty(),
            "Degenerate viewport must not draw a clipped legend");
    require(!overlay::gizmo(QSize(150, 40), font, true, true),
            "Degenerate viewport must not draw a clipped XYZ indicator");
    const QRectF viewport(0, 0, 220, 130);
    std::vector<QRectF> placed;
    for (int i = 0; i < 3; ++i) {
        const auto r = overlay::label(QPointF(110, 65), QSizeF(44, 20), viewport, placed);
        require(r.has_value() && viewport.contains(*r), "Axis label not contained");
        for (const auto &previous : placed)
            require(!r->intersects(previous), "Axis labels overlap");
        placed.push_back(*r);
    }
    require(!overlay::label(QPointF(2, 2), QSizeF(500, 20), viewport, placed).has_value(),
            "Oversized label must not be clipped");
    const std::vector<QRectF> reserved{{0, 0, 220, 24}, {0, 70, 88, 60}};
    const auto nearGizmo = overlay::label(QPointF(80, 80), QSizeF(44, 20), viewport, reserved);
    require(nearGizmo.has_value() && !nearGizmo->intersects(reserved[0]) &&
                !nearGizmo->intersects(reserved[1]),
            "Axis labels must avoid the heading and XYZ indicator");
}
