#include <hornsculpt/mesh_section_view.hpp>

#include <QLocale>
#include <QMouseEvent>
#include <QPainter>
#include <QWheelEvent>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace hornsculpt {
namespace {
double tickStep(double range) {
    const double base = std::pow(10.0, std::floor(std::log10(range / 5)));
    const double fraction = range / (5 * base);
    return base * (fraction <= 1 ? 1 : fraction <= 2 ? 2 : fraction <= 5 ? 5 : 10);
}
QString tickLabel(double value, double step) {
    const int decimals = std::clamp(-static_cast<int>(std::floor(std::log10(step))), 0, 6);
    return QLocale().toString(std::abs(value) < step * 1e-8 ? 0 : value, 'f', decimals);
}
} // namespace

std::vector<SectionSegment> meshSection(const SurfaceMesh &mesh, SectionPlane plane) {
    std::vector<SectionSegment> segments;
    if (mesh.vertices.empty())
        return segments;
    double scale = 0;
    for (const auto &p : mesh.vertices)
        scale = std::max({scale, std::abs(p.x), std::abs(p.y), std::abs(p.z)});
    const double tolerance = std::max(scale * 1e-10, 1e-14);
    const auto distance = [plane](Vec3 p) { return plane == SectionPlane::xz ? p.y : p.x; };
    const auto project = [plane](Vec3 p) {
        return QPointF(p.z, plane == SectionPlane::xz ? p.x : p.y);
    };
    struct CoplanarEdge {
        std::size_t count{};
        BoundaryCondition boundary{};
    };
    std::map<std::pair<std::uint32_t, std::uint32_t>, CoplanarEdge> coplanar;
    for (const auto &triangle : mesh.triangles) {
        std::array<Vec3, 3> p;
        std::array<double, 3> d;
        for (std::size_t i = 0; i < 3; ++i) {
            p[i] = mesh.vertices.at(triangle.vertices[i]);
            d[i] = distance(p[i]);
            if (std::abs(d[i]) <= tolerance)
                d[i] = 0;
        }
        if (d[0] == 0 && d[1] == 0 && d[2] == 0) {
            // Retain only the outline of a coplanar patch, not its tessellation diagonals.
            for (std::size_t i = 0; i < 3; ++i) {
                const auto a = triangle.vertices[i], b = triangle.vertices[(i + 1) % 3];
                auto &edge = coplanar[std::minmax(a, b)];
                ++edge.count;
                edge.boundary = triangle.boundary;
            }
            continue;
        }
        std::vector<QPointF> hits;
        const auto add = [&](QPointF point) {
            if (std::none_of(hits.begin(), hits.end(), [&](QPointF existing) {
                    return std::hypot(point.x() - existing.x(), point.y() - existing.y()) <=
                           tolerance;
                }))
                hits.push_back(point);
        };
        for (std::size_t i = 0; i < 3; ++i) {
            const auto j = (i + 1) % 3;
            if (d[i] == 0)
                add(project(p[i]));
            if ((d[i] < 0 && d[j] > 0) || (d[i] > 0 && d[j] < 0)) {
                const double t = d[i] / (d[i] - d[j]);
                add(project({p[i].x + t * (p[j].x - p[i].x), p[i].y + t * (p[j].y - p[i].y),
                             p[i].z + t * (p[j].z - p[i].z)}));
            }
        }
        if (hits.size() == 2)
            segments.push_back({hits[0], hits[1], triangle.boundary});
    }
    for (const auto &[vertices, edge] : coplanar)
        if (edge.count == 1)
            segments.push_back({project(mesh.vertices[vertices.first]),
                                project(mesh.vertices[vertices.second]), edge.boundary});
    return segments;
}

MeshSectionView::MeshSectionView(QString title, SectionPlane plane, QWidget *parent)
    : QWidget(parent), title_(std::move(title)), plane_(plane) {
    setMinimumSize(280, 180);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    setToolTip(tr("Gleicher Maßstab für beide Achsen · Mausrad: Zoom · Ziehen: Verschieben · "
                  "Doppelklick: Einpassen. "
                  "Rot: Driven. Gestrichelt: BEM-Netz, falls eingeblendet."));
}
void MeshSectionView::setMeshes(std::shared_ptr<const SurfaceMesh> preview,
                                std::shared_ptr<const SurfaceMesh> overlay) {
    if (preview == preview_ && overlay == overlay_)
        return;
    preview_ = std::move(preview);
    overlay_ = std::move(overlay);
    previewSegments_ = preview_ ? meshSection(*preview_, plane_) : std::vector<SectionSegment>{};
    overlaySegments_ = overlay_ ? meshSection(*overlay_, plane_) : std::vector<SectionSegment>{};
    double xmin = std::numeric_limits<double>::infinity(), xmax = -xmin, ymin = xmin, ymax = xmax;
    for (const auto *set : {&previewSegments_, &overlaySegments_})
        for (const auto &segment : *set)
            for (const auto p : {segment.first, segment.second}) {
                xmin = std::min(xmin, p.x());
                xmax = std::max(xmax, p.x());
                ymin = std::min(ymin, p.y());
                ymax = std::max(ymax, p.y());
            }
    dataBounds_ = std::isfinite(xmin) ? QRectF(QPointF(xmin, ymin), QPointF(xmax, ymax)) : QRectF{};
    if (autoFit_)
        fitView();
    update();
}
void MeshSectionView::fitView() {
    autoFit_ = true;
    const double pad = std::max(dataBounds_.width(), dataBounds_.height()) * .06;
    viewBounds_ = dataBounds_.adjusted(-pad, -pad, pad, pad);
    update();
}
void MeshSectionView::setViewBounds(QRectF bounds) {
    if (!std::isfinite(bounds.left()) || !std::isfinite(bounds.top()) ||
        !std::isfinite(bounds.width()) || !std::isfinite(bounds.height()) || bounds.width() <= 0 ||
        bounds.height() <= 0)
        throw std::invalid_argument("Invalid section view bounds.");
    autoFit_ = false;
    viewBounds_ = bounds;
    update();
}
void MeshSectionView::focusDriven() {
    double xmin = std::numeric_limits<double>::infinity(), xmax = -xmin, ymin = xmin, ymax = xmax;
    for (const auto &segment : previewSegments_)
        if (segment.boundary == BoundaryCondition::driven)
            for (const auto p : {segment.first, segment.second}) {
                xmin = std::min(xmin, p.x());
                xmax = std::max(xmax, p.x());
                ymin = std::min(ymin, p.y());
                ymax = std::max(ymax, p.y());
            }
    if (!std::isfinite(xmin))
        return;
    const double pad = std::max({xmax - xmin, ymax - ymin, 1e-6}) * .4;
    setViewBounds(QRectF(QPointF(xmin, ymin), QPointF(xmax, ymax)).adjusted(-pad, -pad, pad, pad));
}
QRectF MeshSectionView::plotArea() const {
    return {62, 27, static_cast<double>(width() - 78), static_cast<double>(height() - 69)};
}
QRectF MeshSectionView::visibleBounds() const {
    auto bounds = viewBounds_;
    if (bounds.width() <= 0 || bounds.height() <= 0)
        return {};
    const double ratio = plotArea().width() / plotArea().height();
    const auto center = bounds.center();
    if (bounds.width() / bounds.height() < ratio)
        bounds.setWidth(bounds.height() * ratio);
    else
        bounds.setHeight(bounds.width() / ratio);
    bounds.moveCenter(center);
    return bounds;
}
void MeshSectionView::paintEvent(QPaintEvent *) {
    QPainter painter(this);
    painter.fillRect(rect(), palette().base());
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setPen(palette().text().color());
    painter.drawText(QRectF(5, 3, width() - 10, 20), Qt::AlignCenter, title_);
    if (previewSegments_.empty())
        return;
    const auto area = plotArea(), bounds = visibleBounds();
    const auto map = [&](QPointF p) {
        return QPointF(area.left() + (p.x() - bounds.left()) / bounds.width() * area.width(),
                       area.bottom() - (p.y() - bounds.top()) / bounds.height() * area.height());
    };
    QFont ticks = font();
    ticks.setPointSizeF(std::max(7.0, font().pointSizeF() - 1));
    painter.setFont(ticks);
    const double xstep = tickStep(bounds.width() * 1000), ystep = tickStep(bounds.height() * 1000);
    for (double value = std::ceil(bounds.left() * 1000 / xstep) * xstep;
         value <= bounds.right() * 1000; value += xstep) {
        const auto position = map({value * .001, 0});
        painter.setPen(QPen(QColor(221, 225, 229), 1));
        painter.drawLine(QPointF(position.x(), area.top()), QPointF(position.x(), area.bottom()));
        painter.setPen(palette().text().color());
        painter.drawText(QRectF(position.x() - 30, area.bottom() + 4, 60, 16), Qt::AlignHCenter,
                         tickLabel(value, xstep));
    }
    for (double value = std::ceil(bounds.top() * 1000 / ystep) * ystep;
         value <= bounds.bottom() * 1000; value += ystep) {
        const auto position = map({0, value * .001});
        painter.setPen(QPen(QColor(221, 225, 229), 1));
        painter.drawLine(QPointF(area.left(), position.y()), QPointF(area.right(), position.y()));
        painter.setPen(palette().text().color());
        painter.drawText(QRectF(9, position.y() - 8, 47, 16), Qt::AlignRight | Qt::AlignVCenter,
                         tickLabel(value, ystep));
    }
    painter.save();
    painter.setClipRect(area);
    painter.setPen(QPen(QColor(146, 152, 158), 1, Qt::DashDotLine));
    painter.drawLine(map({bounds.left(), 0}), map({bounds.right(), 0}));
    const QColor contour = plane_ == SectionPlane::xz ? QColor(28, 102, 180) : QColor(197, 107, 30);
    const auto draw = [&](const auto &segments, bool overlay) {
        for (const auto &segment : segments) {
            const auto color = overlay                                         ? QColor(40, 126, 80)
                               : segment.boundary == BoundaryCondition::driven ? QColor(177, 49, 43)
                                                                               : contour;
            painter.setPen(
                QPen(color, overlay ? 1.3 : 1.8, overlay ? Qt::DashLine : Qt::SolidLine));
            painter.drawLine(map(segment.first), map(segment.second));
        }
    };
    draw(previewSegments_, false);
    draw(overlaySegments_, true);
    painter.restore();
    painter.setPen(QPen(QColor(112, 120, 128), 1));
    painter.drawRect(area);
    painter.setPen(palette().text().color());
    painter.drawText(QRectF(area.left(), height() - 19, area.width(), 17), Qt::AlignCenter,
                     QStringLiteral("z [mm]"));
    painter.save();
    painter.translate(8, area.center().y());
    painter.rotate(-90);
    painter.drawText(QRectF(-area.height() * .5, -8, area.height(), 16), Qt::AlignCenter,
                     plane_ == SectionPlane::xz ? QStringLiteral("x [mm]")
                                                : QStringLiteral("y [mm]"));
    painter.restore();
}
void MeshSectionView::wheelEvent(QWheelEvent *event) {
    if (previewSegments_.empty())
        return;
    const auto area = plotArea(), bounds = visibleBounds();
    if (!area.contains(event->position()))
        return;
    const double factor = std::exp(-event->angleDelta().y() * .001);
    const double x = (event->position().x() - area.left()) / area.width();
    const double y = (area.bottom() - event->position().y()) / area.height();
    const QPointF anchor(bounds.left() + bounds.width() * x, bounds.top() + bounds.height() * y);
    const double w = std::clamp(bounds.width() * factor, 1e-6, 1e4),
                 h = w * bounds.height() / bounds.width();
    setViewBounds({anchor.x() - x * w, anchor.y() - y * h, w, h});
    event->accept();
}
void MeshSectionView::mousePressEvent(QMouseEvent *event) {
    if ((event->button() == Qt::LeftButton || event->button() == Qt::MiddleButton) &&
        plotArea().contains(event->position())) {
        dragging_ = true;
        lastMouse_ = event->position();
    }
}
void MeshSectionView::mouseMoveEvent(QMouseEvent *event) {
    if (!dragging_ || previewSegments_.empty())
        return;
    const auto bounds = visibleBounds(), area = plotArea();
    const auto delta = event->position() - lastMouse_;
    setViewBounds(bounds.translated(-delta.x() * bounds.width() / area.width(),
                                    delta.y() * bounds.height() / area.height()));
    lastMouse_ = event->position();
}
void MeshSectionView::mouseReleaseEvent(QMouseEvent *) {
    dragging_ = false;
}
void MeshSectionView::mouseDoubleClickEvent(QMouseEvent *) {
    fitView();
}

} // namespace hornsculpt
