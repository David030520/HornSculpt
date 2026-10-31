#pragma once

#include <QMatrix4x4>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <QPoint>
#include <cstddef>
#include <functional>
#include <hornsculpt/acoustic_frame.hpp>
#include <hornsculpt/face_color.hpp>
#include <hornsculpt/mesh.hpp>
#include <hornsculpt/pressure_field_plane.hpp>
#include <memory>
#include <optional>
#include <vector>

class QOpenGLShaderProgram;

namespace hornsculpt {

struct RenderBenchmark {
    QString renderer;
    QString driverVersion;
    std::size_t triangles{};
    std::size_t measuredFrames{};
    QSize framebufferSize;
    double medianMilliseconds{};
    double p95Milliseconds{};
};

enum class MeshColorMode {
    boundary,
    elementSize,
    minimumAngle,
    aspectRatio,
    orientation,
    faceColors
};

class MeshViewport final : public QOpenGLWidget, protected QOpenGLFunctions_3_3_Core {
  public:
    explicit MeshViewport(QWidget *parent = nullptr);
    ~MeshViewport() override;
    void setMesh(std::shared_ptr<const SurfaceMesh> mesh);
    void setWireframe(bool enabled);
    void setBackgroundGrid(bool enabled);
    void setColorMode(MeshColorMode mode);
    // Surface colors are presentation data; they never change mesh tags or geometry.
    void setFaceColors(FaceColors colors);
    void fitView();
    // View from the positive physical direction toward the model; geometry and
    // the acoustic/observer frame remain unchanged. An explicit command only.
    void viewFront(Vec3 forward);
    void setPressureField(std::shared_ptr<const PressureFieldPlane> plane,
                          PressureFieldAppearance appearance = {});
    void setPressureFieldVisible(bool visible);
    // Explicit probe mode leaves ordinary CAD face selection unchanged.
    void setFieldProbeEnabled(bool enabled);
    std::optional<std::size_t> pickFieldSample(QPointF position) const;
    bool pressureFieldVisible() const noexcept;
    std::shared_ptr<const PressureFieldPlane> pressureField() const noexcept;
    // Explicit camera command; later colour/frequency updates preserve camera.
    void focusPressureField();
    void setLabel(QString label);
    void setObservationFrame(std::optional<AcousticFrame> frame);
    void setSelectedFaces(std::vector<std::uint32_t> faces);
    // Visibility only: the immutable geometry and solver mesh retain these faces.
    void setHiddenFaces(std::vector<std::uint32_t> faces);
    const std::vector<std::uint32_t> &selectedFaces() const noexcept;
    std::optional<std::uint32_t> pickFace(QPointF position) const;
    std::function<void(std::optional<std::uint32_t>, Qt::KeyboardModifiers)> onFacePicked;
    // Unpaced draw-to-GPU-completion latency, including overlay and camera changes.
    // Window-system presentation and vertical sync are outside this measurement.
    RenderBenchmark benchmarkFrames(std::size_t measuredFrames = 240,
                                    std::size_t warmupFrames = 30);
    QString rendererError() const;

  protected:
    void initializeGL() override;
    void resizeGL(int width, int height) override;
    void paintGL() override;
    void mousePressEvent(QMouseEvent *) override;
    void mouseMoveEvent(QMouseEvent *) override;
    void mouseReleaseEvent(QMouseEvent *) override;
    void wheelEvent(QWheelEvent *) override;
    void mouseDoubleClickEvent(QMouseEvent *) override;

  private:
    void releaseResources();
    void uploadMesh();
    void uploadSelection();
    void updateFit();
    void updateSceneBounds();
    void uploadField();
    QMatrix4x4 viewMatrix() const;
    QMatrix4x4 projectionMatrix() const;
    void paintOverlay();
    std::shared_ptr<const SurfaceMesh> mesh_;
    std::optional<AcousticFrame> observationFrame_;
    std::unique_ptr<QOpenGLShaderProgram> program_;
    std::unique_ptr<QOpenGLShaderProgram> fieldProgram_;
    std::shared_ptr<const PressureFieldPlane> field_;
    PressureFieldAppearance fieldAppearance_;
    QOpenGLVertexArrayObject fieldVao_;
    QOpenGLBuffer fieldBuffer_{QOpenGLBuffer::VertexBuffer};
    QOpenGLBuffer fieldIndexBuffer_{QOpenGLBuffer::IndexBuffer};
    bool fieldVisible_{true};
    bool fieldProbeEnabled_{};
    bool fieldUploadPending_{};
    int fieldIndexCount_{};
    double fieldPressureScale_{1};
    QOpenGLVertexArrayObject vao_;
    QOpenGLBuffer vertexBuffer_{QOpenGLBuffer::VertexBuffer};
    QOpenGLBuffer indexBuffer_{QOpenGLBuffer::IndexBuffer};
    QOpenGLBuffer selectionBuffer_{QOpenGLBuffer::VertexBuffer};
    QOpenGLVertexArrayObject gridVao_;
    QOpenGLBuffer gridBuffer_{QOpenGLBuffer::VertexBuffer};
    QMetaObject::Connection contextConnection_;
    Vec3 center_{};
    Vec3 normalizedHalfSize_{};
    double extent_{1.0};
    double yaw_{0.95};
    double pitch_{-0.3};
    double zoom_{1.4};
    double fittedZoom_{1.4};
    QPoint lastMouse_;
    QPoint pressMouse_;
    QPointF pan_;
    std::vector<std::uint32_t> selectedFaces_;
    std::vector<std::uint32_t> hiddenFaces_;
    std::vector<std::uint32_t> renderFaceIds_;
    std::size_t paintedFrames_{};
    int indexCount_{};
    bool wireframe_{true};
    bool backgroundGrid_{};
    MeshColorMode colorMode_{MeshColorMode::boundary};
    FaceColors faceColors_;
    double colorMinimum_{};
    double colorMaximum_{};
    bool uploadPending_{};
    bool selectionUploadPending_{};
    bool dragStarted_{};
    bool automaticFit_{true};
    bool initialized_{};
    QString label_;
    QString renderer_;
    QString driverVersion_;
    QString rendererError_;
};

} // namespace hornsculpt
