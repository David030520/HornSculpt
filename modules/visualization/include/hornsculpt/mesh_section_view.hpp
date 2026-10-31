#pragma once

#include <QPointF>
#include <QRectF>
#include <QWidget>
#include <hornsculpt/mesh.hpp>
#include <memory>
#include <vector>

namespace hornsculpt {

enum class SectionPlane { xz, yz };

struct SectionSegment {
    QPointF first;
    QPointF second;
    BoundaryCondition boundary;
};

// Intersections with the coordinate plane through the origin. Coordinates are
// (axial z, transverse x/y), in metres; no horn/profile assumptions are made.
std::vector<SectionSegment> meshSection(const SurfaceMesh &mesh, SectionPlane plane);

class MeshSectionView final : public QWidget {
  public:
    MeshSectionView(QString title, SectionPlane plane, QWidget *parent = nullptr);
    void setMeshes(std::shared_ptr<const SurfaceMesh> preview,
                   std::shared_ptr<const SurfaceMesh> overlay = {});
    void fitView();
    void focusDriven();
    void setViewBounds(QRectF bounds);

  protected:
    void paintEvent(QPaintEvent *) override;
    void wheelEvent(QWheelEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void mouseDoubleClickEvent(QMouseEvent *event) override;

  private:
    QRectF plotArea() const;
    QRectF visibleBounds() const;
    QString title_;
    SectionPlane plane_;
    std::shared_ptr<const SurfaceMesh> preview_, overlay_;
    std::vector<SectionSegment> previewSegments_, overlaySegments_;
    QRectF dataBounds_, viewBounds_;
    QPointF lastMouse_;
    bool autoFit_{true};
    bool dragging_{};
};

} // namespace hornsculpt
