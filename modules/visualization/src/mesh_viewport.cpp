#include "overlay_layout.hpp"
#include <QElapsedTimer>
#include <QMouseEvent>
#include <QOpenGLContext>
#include <QOpenGLPaintDevice>
#include <QOpenGLShaderProgram>
#include <QPainter>
#include <QStringList>
#include <QSurfaceFormat>
#include <QToolTip>
#include <QVector2D>
#include <QVector3D>
#include <QVector4D>
#include <QWheelEvent>
#include <algorithm>
#include <array>
#include <cmath>
#include <hornsculpt/level_palette.hpp>
#include <hornsculpt/mesh_viewport.hpp>
#include <hornsculpt/triangle_metrics.hpp>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace hornsculpt {
namespace {
constexpr float cameraDistance = 1.9f;
constexpr double halfFieldOfView = 20.0 * std::numbers::pi / 180.0;

struct RenderVertex {
    float position[3];
    float normal[3];
    float color[3];
};
struct FieldVertex {
    float position[3];
    float pressure[2];
};
constexpr auto fieldVertexShader = R"(
#version 330 core
layout(location=0) in vec3 position;
layout(location=1) in vec2 pressure;
uniform mat4 projection;
uniform mat4 view;
out vec2 phasor;
void main() {
    gl_Position = projection * view * vec4(position, 1.0);
    phasor = pressure;
}
)";
constexpr auto fieldFragmentShader = R"(
#version 330 core
in vec2 phasor;
uniform int quantity;
uniform vec2 limits;
uniform float phase;
uniform float levelOffsetDb;
uniform float realScale;
uniform float bandStep;
uniform vec3 levelPalette[10];
out vec4 fragmentColor;
void main() {
    float amplitude = length(phasor);
    // A zero phasor has no phase, but its level floor and zero real part are valid.
    bool zeroPressure = phasor.x == 0.0 && phasor.y == 0.0;
    if (quantity == 1 && zeroPressure) discard;
    float value = quantity == 1 ? degrees(atan(phasor.y, phasor.x)) :
        quantity == 2 ? (phasor.x * cos(phase) - phasor.y * sin(phase)) * realScale :
        zeroPressure ? limits.x :
        20.0 * log(max(amplitude, 1e-30)) / log(10.0) + levelOffsetDb;
    float count = ceil((limits.y - limits.x) / bandStep);
    float band = clamp(floor((limits.y - clamp(value, limits.x, limits.y)) / bandStep),
                       0.0, max(0.0, count - 1.0));
    if (quantity == 1) {
        // Ten 36° sectors: sector zero straddles the ±180° wrap.
        band = floor(mod(value + 198.0, 360.0) / 36.0);
        count = 10.0;
    }
    float position = count > 1.0 ? band / (count - 1.0) * 9.0 : 0.0;
    int index = int(clamp(floor(position), 0.0, 9.0));
    vec3 color = mix(levelPalette[index], levelPalette[min(index + 1, 9)], position - float(index));
    // Opaque quantitative colours preserve the colour legend. The normal
    // depth buffer handles CAD in front/behind; no through-solid overlay.
    fragmentColor = vec4(color, 1.0);
}
)";

constexpr auto vertexShader = R"(
#version 330 core
layout(location=0) in vec3 position;
layout(location=1) in vec3 normal;
layout(location=2) in vec3 color;
layout(location=3) in float selection;
uniform mat4 projection;
uniform mat4 view;
out vec3 viewNormal;
out vec3 baseColor;
flat out float selected;
void main() {
    gl_Position = projection * view * vec4(position, 1.0);
    viewNormal = mat3(view) * normal;
    baseColor = color;
    selected = selection;
}
)";

constexpr auto fragmentShader = R"(
#version 330 core
in vec3 viewNormal;
in vec3 baseColor;
flat in float selected;
uniform int wirePass;
uniform int unlit;
uniform int orientationPass;
out vec4 fragmentColor;
void main() {
    if (wirePass != 0) {
        fragmentColor = vec4(0.20, 0.24, 0.28, 0.54);
        return;
    }
    vec3 color = baseColor;
    if (unlit == 0 && orientationPass != 0)
        color = gl_FrontFacing ? vec3(0.20, 0.66, 0.39) : vec3(0.88, 0.29, 0.24);
    if (unlit == 0 && selected > 0.5)
        color = orientationPass != 0 ? mix(color, vec3(1.0, 0.87, 0.35), 0.35) :
            vec3(0.91, 0.65, 0.23);
    float light = unlit != 0 ? 1.0 :
        0.60 + 0.38 * abs(dot(normalize(viewNormal), normalize(vec3(0.3, 0.6, 0.74))));
    fragmentColor = vec4(color * light, 1.0);
}
)";

Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

bool usesMetricColors(MeshColorMode mode) {
    return mode == MeshColorMode::elementSize || mode == MeshColorMode::minimumAngle ||
           mode == MeshColorMode::aspectRatio;
}

} // namespace

MeshViewport::MeshViewport(QWidget *parent) : QOpenGLWidget(parent) {
    QSurfaceFormat requested;
    requested.setVersion(3, 3);
    requested.setProfile(QSurfaceFormat::CoreProfile);
    requested.setDepthBufferSize(24);
    requested.setStencilBufferSize(8);
    requested.setSamples(4);
    setFormat(requested);
    setMinimumSize(350, 260);
    setFocusPolicy(Qt::StrongFocus);
    setToolTip(tr("Klick: Fläche wählen · Strg + Klick: Auswahl umschalten · Ziehen: drehen · "
                  "Umschalt + Ziehen: verschieben · "
                  "Mausrad: zoomen · Doppelklick: einpassen"));
}

MeshViewport::~MeshViewport() {
    QObject::disconnect(contextConnection_);
    if (context()) {
        makeCurrent();
        releaseResources();
        doneCurrent();
    }
}

void MeshViewport::releaseResources() {
    interfaceBuffer_.destroy();
    interfaceVao_.destroy();
    fieldVao_.destroy();
    fieldBuffer_.destroy();
    fieldIndexBuffer_.destroy();
    fieldProgram_.reset();
    vertexBuffer_.destroy();
    indexBuffer_.destroy();
    selectionBuffer_.destroy();
    vao_.destroy();
    gridBuffer_.destroy();
    gridVao_.destroy();
    program_.reset();
    initialized_ = false;
}

void MeshViewport::initializeGL() {
    if (!initializeOpenGLFunctions()) {
        rendererError_ = tr("OpenGL 3.3 wird nicht unterstützt.");
        return;
    }
    renderer_ = QString::fromLatin1(reinterpret_cast<const char *>(glGetString(GL_RENDERER)));
    driverVersion_ = QString::fromLatin1(reinterpret_cast<const char *>(glGetString(GL_VERSION)));
    contextConnection_ = connect(
        context(), &QOpenGLContext::aboutToBeDestroyed, this,
        [this] {
            makeCurrent();
            releaseResources();
            doneCurrent();
        },
        Qt::DirectConnection);
    program_ = std::make_unique<QOpenGLShaderProgram>();
    if (!program_->addShaderFromSourceCode(QOpenGLShader::Vertex, vertexShader) ||
        !program_->addShaderFromSourceCode(QOpenGLShader::Fragment, fragmentShader) ||
        !program_->link()) {
        rendererError_ = program_->log();
        program_.reset();
        return;
    }
    vao_.create();
    fieldProgram_ = std::make_unique<QOpenGLShaderProgram>();
    if (!fieldProgram_->addShaderFromSourceCode(QOpenGLShader::Vertex, fieldVertexShader) ||
        !fieldProgram_->addShaderFromSourceCode(QOpenGLShader::Fragment, fieldFragmentShader) ||
        !fieldProgram_->link()) {
        rendererError_ = fieldProgram_->log();
        return;
    }
    fieldVao_.create();
    fieldBuffer_.create();
    fieldIndexBuffer_.create();
    fieldVao_.bind();
    fieldBuffer_.bind();
    fieldIndexBuffer_.bind();
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(FieldVertex),
                          reinterpret_cast<const void *>(offsetof(FieldVertex, position)));
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(FieldVertex),
                          reinterpret_cast<const void *>(offsetof(FieldVertex, pressure)));
    fieldVao_.release();
    fieldBuffer_.release();
    vertexBuffer_.create();
    indexBuffer_.create();
    selectionBuffer_.create();
    gridVao_.create();
    gridBuffer_.create();
    const auto attributes = [this] {
        for (GLuint i = 0; i < 3; ++i)
            glEnableVertexAttribArray(i);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(RenderVertex),
                              reinterpret_cast<const void *>(offsetof(RenderVertex, position)));
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(RenderVertex),
                              reinterpret_cast<const void *>(offsetof(RenderVertex, normal)));
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(RenderVertex),
                              reinterpret_cast<const void *>(offsetof(RenderVertex, color)));
    };
    vao_.bind();
    vertexBuffer_.bind();
    indexBuffer_.bind();
    attributes();
    selectionBuffer_.bind();
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(float), nullptr);
    vao_.release();
    std::vector<RenderVertex> grid;
    for (int i = -5; i <= 5; ++i) {
        const float t = static_cast<float>(i) * .2f;
        for (const auto &p : std::array<std::array<float, 3>, 4>{
                 {{t, -.62f, -1}, {t, -.62f, 1}, {-1, -.62f, t}, {1, -.62f, t}}})
            grid.push_back({{p[0], p[1], p[2]}, {0, 1, 0}, {.82f, .85f, .87f}});
    }
    interfaceVao_.create();
    interfaceBuffer_.create();
    interfaceVao_.bind();
    interfaceBuffer_.bind();
    attributes();
    glDisableVertexAttribArray(3);
    interfaceVao_.release();
    interfaceBuffer_.release();
    gridVao_.bind();
    gridBuffer_.bind();
    gridBuffer_.allocate(grid.data(), static_cast<int>(grid.size() * sizeof(RenderVertex)));
    attributes();
    glDisableVertexAttribArray(3);
    gridVao_.release();
    vertexBuffer_.release();
    gridBuffer_.release();
    selectionBuffer_.release();
    initialized_ = true;
    uploadPending_ = true;
    fieldUploadPending_ = true;
    interfaceUploadPending_ = true;
    rendererError_.clear();
}

void MeshViewport::setMesh(std::shared_ptr<const SurfaceMesh> mesh) {
    if (mesh_ == mesh)
        return;
    if (mesh) {
        for (const auto &p : mesh->vertices)
            if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
                throw std::invalid_argument("Viewer: nicht endliche Koordinate.");
        for (const auto &triangle : mesh->triangles)
            for (const auto id : triangle.vertices)
                if (id >= mesh->vertices.size())
                    throw std::invalid_argument("Viewer: ungültiger Dreiecksindex.");
        if (mesh->triangles.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max()) / (3 * sizeof(RenderVertex)))
            throw std::invalid_argument("Viewer: Indexpuffer zu groß.");
    }
    mesh_ = std::move(mesh);
    selectedFaces_.clear();
    updateSceneBounds();
    if (automaticFit_)
        updateFit();
    uploadPending_ = true;
    fieldUploadPending_ = true;
    update();
}
void MeshViewport::setInterfaceMesh(std::shared_ptr<const SurfaceMesh> mesh,
                                    std::vector<InterfaceOrientationLabels> groups) {
    if (interfaceMesh_ == mesh && interfaceGroups_ == groups)
        return;
    if (mesh) {
        for (const auto &p : mesh->vertices)
            if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
                throw std::invalid_argument("Viewer: nicht endliche Interface-Koordinate.");
        for (const auto &triangle : mesh->triangles)
            for (const auto id : triangle.vertices)
                if (id >= mesh->vertices.size())
                    throw std::invalid_argument("Viewer: ungültiger Interface-Dreiecksindex.");
        if (mesh->triangles.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max()) / (6 * sizeof(RenderVertex)))
            throw std::invalid_argument("Viewer: Interface-Puffer zu groß.");
    }
    for (const auto &group : groups)
        for (const auto index : group.triangleIndices)
            if (!mesh || index >= mesh->triangles.size())
                throw std::invalid_argument("Viewer: ungueltige Interface-Gruppenzuordnung.");
    // Account for the worst case: every sampled triangle has its own face.
    const auto maxVertices =
        static_cast<std::size_t>(std::numeric_limits<int>::max()) / sizeof(RenderVertex);
    std::size_t vertexBudget = mesh ? mesh->triangles.size() * 6 : 0;
    const auto addBudget = [&](std::size_t triangles) {
        if (vertexBudget > maxVertices - 78 || triangles > (maxVertices - 78 - vertexBudget) / 12)
            throw std::invalid_argument("Viewer: Interface-Puffer zu groß.");
        vertexBudget += triangles * 12;
    };
    if (groups.empty() && mesh)
        addBudget(mesh->triangles.size());
    else
        for (const auto &group : groups)
            addBudget(group.triangleIndices.size());
    interfaceMesh_ = std::move(mesh);
    interfaceGroups_ = std::move(groups);
    interfaceUploadPending_ = true;
    update();
}
void MeshViewport::setSymmetryPlanes(std::optional<OrthogonalSymmetryPlanes> planes) {
    if (symmetryPlanes_ == planes)
        return;
    const auto finite = [](Vec3 v) {
        return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
    };
    const auto nonzero = [](Vec3 v) { return v.x != 0 || v.y != 0 || v.z != 0; };
    if (planes &&
        (!finite(planes->origin) || !finite(planes->firstNormal) || !finite(planes->secondNormal) ||
         !nonzero(planes->firstNormal) || !nonzero(planes->secondNormal) ||
         !std::isfinite(planes->planeDistanceTolerance)))
        throw std::invalid_argument("Viewer: ungültige Symmetrieebene.");
    symmetryPlanes_ = std::move(planes);
    if (fullSymmetryVisible_) {
        updateSceneBounds();
        if (automaticFit_)
            updateFit();
    }
    interfaceUploadPending_ = true;
    update();
}
void MeshViewport::setFullSymmetryVisible(bool visible) {
    if (fullSymmetryVisible_ == visible)
        return;
    fullSymmetryVisible_ = visible;
    fitView();
}
void MeshViewport::setExteriorBaffle(std::optional<ExteriorBaffle> baffle) {
    if (exteriorBaffle_ == baffle)
        return;
    const auto finite = [](Vec3 v) {
        return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
    };
    if (baffle && (!finite(baffle->origin) || !finite(baffle->normal) ||
                   (baffle->normal.x == 0 && baffle->normal.y == 0 && baffle->normal.z == 0) ||
                   !std::isfinite(baffle->planeDistanceTolerance)))
        throw std::invalid_argument("Viewer: ungültige Schallwandebene.");
    exteriorBaffle_ = std::move(baffle);
    interfaceUploadPending_ = true;
    update();
}
void MeshViewport::setAcousticGuidesVisible(bool visible) {
    if (acousticGuidesVisible_ == visible)
        return;
    acousticGuidesVisible_ = visible;
    update();
}
void MeshViewport::setInterfaceNormalsVisible(bool visible) {
    if (interfaceNormalsVisible_ == visible)
        return;
    interfaceNormalsVisible_ = visible;
    update();
}
void MeshViewport::setInterfaceVisible(bool visible) {
    if (interfaceVisible_ == visible)
        return;
    interfaceVisible_ = visible;
    update();
}
void MeshViewport::updateSceneBounds() {
    bool hasPoints = false;
    Vec3 low{}, high{};
    const auto include = [&](Vec3 p) {
        if (!hasPoints) {
            low = high = p;
            hasPoints = true;
        }
        low.x = std::min(low.x, p.x);
        low.y = std::min(low.y, p.y);
        low.z = std::min(low.z, p.z);
        high.x = std::max(high.x, p.x);
        high.y = std::max(high.y, p.y);
        high.z = std::max(high.z, p.z);
    };
    if (mesh_)
        for (const auto p : mesh_->vertices) {
            include(p);
            if (fullSymmetryVisible_ && symmetryPlanes_) {
                const auto reflect = [&](Vec3 point, Vec3 normal) {
                    const double length = std::hypot(normal.x, normal.y, normal.z);
                    normal = {normal.x / length, normal.y / length, normal.z / length};
                    const double distance = dot(subtract(point, symmetryPlanes_->origin), normal);
                    return Vec3{point.x - 2 * distance * normal.x,
                                point.y - 2 * distance * normal.y,
                                point.z - 2 * distance * normal.z};
                };
                const auto first = reflect(p, symmetryPlanes_->firstNormal);
                include(first);
                include(reflect(p, symmetryPlanes_->secondNormal));
                include(reflect(first, symmetryPlanes_->secondNormal));
            }
        }
    if (field_)
        for (const auto p : field_->positions)
            include(p);
    if (hasPoints) {
        center_ = {low.x + (high.x - low.x) * .5, low.y + (high.y - low.y) * .5,
                   low.z + (high.z - low.z) * .5};
        extent_ = std::max({high.x - low.x, high.y - low.y, high.z - low.z, 1e-12});
        normalizedHalfSize_ = {(high.x - low.x) / (2 * extent_), (high.y - low.y) / (2 * extent_),
                               (high.z - low.z) / (2 * extent_)};
    } else
        normalizedHalfSize_ = {};
    uploadPending_ = true;
    fieldUploadPending_ = true;
    interfaceUploadPending_ = true;
}

void MeshViewport::setPressureField(std::shared_ptr<const PressureFieldPlane> plane,
                                    PressureFieldAppearance appearance) {
    validatePressureFieldAppearance(appearance);
    if (plane)
        validatePressureFieldPlane(*plane);
    const bool changed = plane != field_;
    field_ = std::move(plane);
    fieldAppearance_ = appearance;
    if (changed) {
        fieldUploadPending_ = true;
        // Clearing during a new frequency job must not move the user's camera
        // or briefly renormalize the CAD. Fit explicitly if a mesh-only view is wanted.
        if (field_) {
            updateSceneBounds();
            if (automaticFit_)
                updateFit();
        }
    }
    update();
}
void MeshViewport::setPressureFieldVisible(bool visible) {
    fieldVisible_ = visible;
    update();
}
bool MeshViewport::pressureFieldVisible() const noexcept {
    return fieldVisible_ && bool(field_);
}
std::shared_ptr<const PressureFieldPlane> MeshViewport::pressureField() const noexcept {
    return field_;
}
void MeshViewport::focusPressureField() {
    if (!field_)
        return;
    const auto frame = normalizedAcousticFrame(field_->frame);
    const auto lateral = field_->vertical ? verticalDirection(frame) : frame.horizontal;
    const auto normal = cross(lateral, frame.forward);
    // An oblique view of the actual field plane keeps its area and the CAD
    // visible together; looking down the acoustic axis would collapse it to a line.
    viewFront({normal.x + .65 * frame.forward.x + .25 * lateral.x,
               normal.y + .65 * frame.forward.y + .25 * lateral.y,
               normal.z + .65 * frame.forward.z + .25 * lateral.z});
}
void MeshViewport::uploadField() {
    fieldIndexCount_ = 0;
    if (field_) {
        auto indices = pressureFieldTriangleIndices(*field_);
        fieldPressureScale_ = 0;
        for (std::size_t i = 0; i < field_->pressure.size(); ++i)
            if (field_->valid[i])
                fieldPressureScale_ =
                    std::max({fieldPressureScale_, std::abs(field_->pressure[i].real()),
                              std::abs(field_->pressure[i].imag())});
        if (fieldPressureScale_ == 0)
            fieldPressureScale_ = 1;
        std::vector<FieldVertex> vertices;
        vertices.reserve(field_->positions.size());
        for (std::size_t i = 0; i < field_->positions.size(); ++i) {
            const auto p = field_->positions[i];
            const auto q = field_->valid[i] ? field_->pressure[i] / fieldPressureScale_
                                            : std::complex<double>{};
            vertices.push_back({{static_cast<float>((p.x - center_.x) / extent_),
                                 static_cast<float>((p.y - center_.y) / extent_),
                                 static_cast<float>((p.z - center_.z) / extent_)},
                                {static_cast<float>(q.real()), static_cast<float>(q.imag())}});
        }
        fieldVao_.bind();
        fieldBuffer_.bind();
        fieldBuffer_.allocate(vertices.data(),
                              static_cast<int>(vertices.size() * sizeof(FieldVertex)));
        fieldIndexBuffer_.bind();
        fieldIndexBuffer_.allocate(indices.data(),
                                   static_cast<int>(indices.size() * sizeof(std::uint32_t)));
        fieldVao_.release();
        fieldBuffer_.release();
        fieldIndexCount_ = static_cast<int>(indices.size());
    }
    fieldUploadPending_ = false;
}

void MeshViewport::uploadMesh() {
    std::vector<RenderVertex> vertices;
    renderFaceIds_.clear();
    std::vector<std::uint32_t> indices;
    const auto metricValue = [this](const TriangleMetrics &m) {
        switch (colorMode_) {
        case MeshColorMode::elementSize:
            return m.medianEdge * 1000;
        case MeshColorMode::minimumAngle:
            return m.minimumAngleDegrees;
        case MeshColorMode::aspectRatio:
            return m.maximumAspectRatio;
        default:
            return 0.0;
        }
    };
    const bool metricColor = usesMetricColors(colorMode_);
    if (mesh_ && metricColor) {
        colorMinimum_ = std::numeric_limits<double>::infinity();
        colorMaximum_ = -colorMinimum_;
        for (const auto &t : mesh_->triangles) {
            const auto value = metricValue(triangleMetrics(mesh_->vertices[t.vertices[0]],
                                                           mesh_->vertices[t.vertices[1]],
                                                           mesh_->vertices[t.vertices[2]]));
            if (std::isfinite(value)) {
                colorMinimum_ = std::min(colorMinimum_, value);
                colorMaximum_ = std::max(colorMaximum_, value);
            }
        }
        if (colorMode_ == MeshColorMode::minimumAngle) {
            colorMinimum_ = 0;
            colorMaximum_ = 60;
        } else if (colorMode_ == MeshColorMode::aspectRatio) {
            colorMinimum_ = 2 / std::sqrt(3.0);
            colorMaximum_ = 4;
        }
        if (!std::isfinite(colorMinimum_)) {
            colorMinimum_ = 0;
            colorMaximum_ = 1;
        }
    }
    if (mesh_) {
        vertices.reserve(mesh_->vertices.size());
        indices.reserve(mesh_->triangles.size() * 3);
        std::unordered_map<std::uint64_t, std::uint32_t> vertexMap;
        // A face boundary is a real shading crease; vertices shared inside a face keep smooth
        // normals.
        for (const auto &triangle : mesh_->triangles) {
            if (std::binary_search(hiddenFaces_.begin(), hiddenFaces_.end(), triangle.faceId))
                continue;
            const auto &p = mesh_->vertices[triangle.vertices[0]];
            const auto normal = cross(subtract(mesh_->vertices[triangle.vertices[1]], p),
                                      subtract(mesh_->vertices[triangle.vertices[2]], p));
            auto color = triangle.boundary == BoundaryCondition::driven
                             ? std::array<float, 3>{.24f, .55f, .83f}
                             : std::array<float, 3>{.80f, .82f, .84f};
            if (colorMode_ == MeshColorMode::faceColors) {
                color = {.75f, .75f, .75f};
                if (const auto found = faceColors_.find(triangle.faceId);
                    found != faceColors_.end()) {
                    const auto &rgb = found->second;
                    color = {rgb.red / 255.f, rgb.green / 255.f, rgb.blue / 255.f};
                }
            }
            if (metricColor) {
                const auto value =
                    metricValue(triangleMetrics(p, mesh_->vertices[triangle.vertices[1]],
                                                mesh_->vertices[triangle.vertices[2]]));
                double fraction =
                    std::isfinite(value)
                        ? std::clamp((value - colorMinimum_) /
                                         std::max(1e-12, colorMaximum_ - colorMinimum_),
                                     0.0, 1.0)
                        : 1;
                if (colorMode_ == MeshColorMode::minimumAngle)
                    fraction = 1 - fraction;
                const auto tint =
                    QColor::fromHsvF(static_cast<float>((1 - fraction) * .62), .75f, .9f);
                color = {static_cast<float>(tint.redF()), static_cast<float>(tint.greenF()),
                         static_cast<float>(tint.blueF())};
                const auto normalLength = std::hypot(normal.x, normal.y, normal.z);
                for (const auto original : triangle.vertices) {
                    const auto &v = mesh_->vertices[original];
                    indices.push_back(static_cast<std::uint32_t>(vertices.size()));
                    vertices.push_back(
                        {{static_cast<float>((v.x - center_.x) / extent_),
                          static_cast<float>((v.y - center_.y) / extent_),
                          static_cast<float>((v.z - center_.z) / extent_)},
                         {static_cast<float>(normal.x / std::max(1e-30, normalLength)),
                          static_cast<float>(normal.y / std::max(1e-30, normalLength)),
                          static_cast<float>(normal.z / std::max(1e-30, normalLength))},
                         {color[0], color[1], color[2]}});
                    renderFaceIds_.push_back(triangle.faceId);
                }
                continue;
            }
            for (const auto original : triangle.vertices) {
                const auto key = (static_cast<std::uint64_t>(triangle.faceId) << 32) | original;
                auto [found, inserted] =
                    vertexMap.try_emplace(key, static_cast<std::uint32_t>(vertices.size()));
                if (inserted) {
                    const auto &v = mesh_->vertices[original];
                    vertices.push_back({{static_cast<float>((v.x - center_.x) / extent_),
                                         static_cast<float>((v.y - center_.y) / extent_),
                                         static_cast<float>((v.z - center_.z) / extent_)},
                                        {0, 0, 0},
                                        {color[0], color[1], color[2]}});
                    renderFaceIds_.push_back(triangle.faceId);
                }
                auto &vertex = vertices[found->second];
                vertex.normal[0] += static_cast<float>(normal.x);
                vertex.normal[1] += static_cast<float>(normal.y);
                vertex.normal[2] += static_cast<float>(normal.z);
                indices.push_back(found->second);
            }
        }
        for (auto &vertex : vertices) {
            const float norm = std::hypot(vertex.normal[0], vertex.normal[1], vertex.normal[2]);
            if (norm > 0)
                for (auto &component : vertex.normal)
                    component /= norm;
        }
    }
    if (vertices.size() >
        static_cast<std::size_t>(std::numeric_limits<int>::max()) / sizeof(RenderVertex))
        throw std::length_error("Viewer: Vertexpuffer zu groß.");
    vao_.bind();
    vertexBuffer_.bind();
    vertexBuffer_.setUsagePattern(QOpenGLBuffer::StaticDraw);
    vertexBuffer_.allocate(vertices.data(),
                           static_cast<int>(vertices.size() * sizeof(RenderVertex)));
    indexBuffer_.bind();
    indexBuffer_.setUsagePattern(QOpenGLBuffer::StaticDraw);
    indexBuffer_.allocate(indices.data(), static_cast<int>(indices.size() * sizeof(std::uint32_t)));
    vao_.release();
    vertexBuffer_.release();
    indexCount_ = static_cast<int>(indices.size());
    uploadPending_ = false;
    selectionUploadPending_ = true;
}

void MeshViewport::uploadInterfaces() {
    std::vector<RenderVertex> vertices;
    if (interfaceMesh_) {
        vertices.reserve(interfaceMesh_->triangles.size() * 6);
        for (const auto &triangle : interfaceMesh_->triangles) {
            for (const auto local : std::array<unsigned, 6>{0, 1, 1, 2, 2, 0}) {
                const auto &v = interfaceMesh_->vertices[triangle.vertices[local]];
                vertices.push_back({{static_cast<float>((v.x - center_.x) / extent_),
                                     static_cast<float>((v.y - center_.y) / extent_),
                                     static_cast<float>((v.z - center_.z) / extent_)},
                                    {0, 0, 1},
                                    {0.0f, 0.65f, 0.72f}});
            }
        }
    }
    interfaceVertexCount_ = static_cast<int>(vertices.size());
    const auto add = [&](Vec3 p, std::array<float, 3> color) {
        vertices.push_back({{float((p.x - center_.x) / extent_), float((p.y - center_.y) / extent_),
                             float((p.z - center_.z) / extent_)},
                            {0, 0, 1},
                            {color[0], color[1], color[2]}});
    };
    const auto sum = [](Vec3 a, Vec3 b) { return Vec3{a.x + b.x, a.y + b.y, a.z + b.z}; };
    const auto scale = [](Vec3 a, double k) { return Vec3{a.x * k, a.y * k, a.z * k}; };
    const auto unit = [&](Vec3 a) {
        const double length = std::sqrt(dot(a, a));
        return length > 0 && std::isfinite(length) ? scale(a, 1 / length) : Vec3{};
    };
    const auto arrow = [&](Vec3 origin, Vec3 direction, double length, std::array<float, 3> color) {
        const Vec3 n = unit(direction);
        if (dot(n, n) == 0)
            return;
        const Vec3 side = unit(cross(n, std::abs(n.z) < .9 ? Vec3{0, 0, 1} : Vec3{0, 1, 0}));
        const Vec3 tip = sum(origin, scale(n, length));
        add(origin, color);
        add(tip, color);
        add(tip, color);
        add(sum(sum(tip, scale(n, -length * .25)), scale(side, length * .12)), color);
        add(tip, color);
        add(sum(sum(tip, scale(n, -length * .25)), scale(side, -length * .12)), color);
    };
    if (interfaceMesh_) {
        auto groups = interfaceGroups_;
        if (groups.empty()) {
            InterfaceOrientationLabels fallback;
            fallback.triangleIndices.resize(interfaceMesh_->triangles.size());
            std::iota(fallback.triangleIndices.begin(), fallback.triangleIndices.end(), 0);
            groups.push_back(std::move(fallback));
        }
        for (const auto &group : groups) {
            std::unordered_map<std::uint32_t, std::vector<std::size_t>> faces;
            for (const auto index : group.triangleIndices)
                faces[interfaceMesh_->triangles[index].faceId].push_back(index);
            for (const auto &[face, indices] : faces) {
                (void)face;
                const std::size_t step = std::max<std::size_t>(1, (indices.size() + 15) / 16);
                for (std::size_t i = 0; i < indices.size(); i += step) {
                    const auto &t = interfaceMesh_->triangles[indices[i]];
                    const auto a = interfaceMesh_->vertices[t.vertices[0]];
                    const auto b = interfaceMesh_->vertices[t.vertices[1]];
                    const auto c = interfaceMesh_->vertices[t.vertices[2]];
                    const Vec3 center = scale(sum(sum(a, b), c), 1. / 3);
                    const Vec3 n = unit(cross(subtract(b, a), subtract(c, a)));
                    arrow(center, n, extent_ * .10, {0.f, .55f, .7f});
                    arrow(center, scale(n, -1), extent_ * .10, {.88f, .39f, .08f});
                }
            }
        }
    }
    normalVertexCount_ = static_cast<int>(vertices.size()) - interfaceVertexCount_;
    const auto plane = [&](Vec3 origin, Vec3 normal, std::array<float, 3> color) {
        const Vec3 n = unit(normal);
        if (dot(n, n) == 0)
            return;
        const Vec3 u = unit(cross(n, std::abs(n.z) < .9 ? Vec3{0, 0, 1} : Vec3{0, 1, 0}));
        const Vec3 v = unit(cross(n, u));
        // A finite view patch represents the plane at its actual world position.
        const Vec3 base = subtract(center_, scale(n, dot(subtract(center_, origin), n)));
        const double radius = extent_ * .8;
        for (int i = -2; i <= 2; ++i) {
            const double offset = radius * i / 2.;
            const Vec3 a = sum(base, scale(u, offset)), b = sum(base, scale(v, offset));
            add(sum(a, scale(v, -radius)), color);
            add(sum(a, scale(v, radius)), color);
            add(sum(b, scale(u, -radius)), color);
            add(sum(b, scale(u, radius)), color);
        }
        arrow(base, n, extent_ * .22, color);
    };
    if (symmetryPlanes_) {
        plane(symmetryPlanes_->origin, symmetryPlanes_->firstNormal, {.42f, .34f, .74f});
        plane(symmetryPlanes_->origin, symmetryPlanes_->secondNormal, {.32f, .58f, .30f});
    }
    if (exteriorBaffle_)
        plane(exteriorBaffle_->origin, exteriorBaffle_->normal, {.75f, .36f, .13f});
    guideVertexCount_ =
        static_cast<int>(vertices.size()) - interfaceVertexCount_ - normalVertexCount_;
    interfaceVao_.bind();
    interfaceBuffer_.bind();
    interfaceBuffer_.setUsagePattern(QOpenGLBuffer::StaticDraw);
    interfaceBuffer_.allocate(vertices.data(),
                              static_cast<int>(vertices.size() * sizeof(RenderVertex)));
    interfaceVao_.release();
    interfaceBuffer_.release();
    interfaceUploadPending_ = false;
}

void MeshViewport::uploadSelection() {
    std::vector<float> selected;
    selected.reserve(renderFaceIds_.size());
    for (auto face : renderFaceIds_)
        selected.push_back(
            std::binary_search(selectedFaces_.begin(), selectedFaces_.end(), face) ? 1.f : 0.f);
    selectionBuffer_.bind();
    selectionBuffer_.setUsagePattern(QOpenGLBuffer::DynamicDraw);
    selectionBuffer_.allocate(selected.data(), static_cast<int>(selected.size() * sizeof(float)));
    selectionBuffer_.release();
    selectionUploadPending_ = false;
}

std::vector<QMatrix4x4> MeshViewport::symmetryTransforms() const {
    std::vector<QMatrix4x4> transforms(1);
    if (!fullSymmetryVisible_ || !symmetryPlanes_)
        return transforms;
    const auto reflection = [&](Vec3 normal) {
        const double length = std::hypot(normal.x, normal.y, normal.z);
        const double n[] = {normal.x / length, normal.y / length, normal.z / length};
        const Vec3 origin = subtract(symmetryPlanes_->origin, center_);
        const double distance = (origin.x * n[0] + origin.y * n[1] + origin.z * n[2]) / extent_;
        QMatrix4x4 matrix;
        for (int row = 0; row < 3; ++row) {
            for (int column = 0; column < 3; ++column)
                matrix(row, column) =
                    static_cast<float>((row == column ? 1.0 : 0.0) - 2 * n[row] * n[column]);
            matrix(row, 3) = static_cast<float>(2 * distance * n[row]);
        }
        return matrix;
    };
    const auto first = reflection(symmetryPlanes_->firstNormal);
    const auto second = reflection(symmetryPlanes_->secondNormal);
    transforms.push_back(first);
    transforms.push_back(second);
    transforms.push_back(second * first);
    return transforms;
}

QMatrix4x4 MeshViewport::viewMatrix() const {
    const double unitsPerPixel =
        2 * cameraDistance * std::tan(halfFieldOfView) / zoom_ / std::max(1, height());
    QMatrix4x4 view;
    view.translate(static_cast<float>(pan_.x() * unitsPerPixel),
                   static_cast<float>(-pan_.y() * unitsPerPixel), -cameraDistance);
    view.rotate(static_cast<float>(pitch_ * 180 / std::numbers::pi), 1, 0, 0);
    view.rotate(static_cast<float>(yaw_ * 180 / std::numbers::pi), 0, 1, 0);
    return view;
}

QMatrix4x4 MeshViewport::projectionMatrix() const {
    QMatrix4x4 projection;
    const float fov = static_cast<float>(2 * std::atan(std::tan(halfFieldOfView) / zoom_) * 180 /
                                         std::numbers::pi);
    projection.perspective(fov, static_cast<float>(width()) / std::max(1, height()), .01f, 100.f);
    return projection;
}

void MeshViewport::paintGL() {
    if (!initialized_) {
        paintOverlay();
        return;
    }
    glViewport(0, 0, static_cast<int>(width() * devicePixelRatioF()),
               static_cast<int>(height() * devicePixelRatioF()));
    glDisable(GL_SCISSOR_TEST);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glDepthMask(GL_TRUE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glClearColor(228.f / 255, 232.f / 255, 235.f / 255, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if (uploadPending_)
        uploadMesh();
    if (selectionUploadPending_)
        uploadSelection();
    if (fieldUploadPending_)
        uploadField();
    if (interfaceUploadPending_)
        uploadInterfaces();
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glDepthMask(GL_TRUE);
    glDisable(GL_CULL_FACE);
    glFrontFace(GL_CCW);
    glDisable(GL_BLEND);
    glEnable(GL_MULTISAMPLE);
    program_->bind();
    program_->setUniformValue("view", viewMatrix());
    program_->setUniformValue("projection", projectionMatrix());
    program_->setUniformValue("wirePass", 0);
    program_->setUniformValue("orientationPass", colorMode_ == MeshColorMode::orientation ? 1 : 0);
    program_->setUniformValue("unlit", 1);
    if (backgroundGrid_) {
        gridVao_.bind();
        glDrawArrays(GL_LINES, 0, 44);
        gridVao_.release();
    }
    if (indexCount_ > 0) {
        program_->setUniformValue("unlit", 0);
        vao_.bind();
        if (wireframe_) {
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.f, 1.f);
        }
        const auto transforms = symmetryTransforms();
        for (std::size_t i = 0; i < transforms.size(); ++i) {
            program_->setUniformValue("view", viewMatrix() * transforms[i]);
            glFrontFace(i == 1 || i == 2 ? GL_CW : GL_CCW);
            glDrawElements(GL_TRIANGLES, indexCount_, GL_UNSIGNED_INT, nullptr);
        }
        glDisable(GL_POLYGON_OFFSET_FILL);
        if (wireframe_) {
            program_->setUniformValue("wirePass", 1);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glDepthMask(GL_FALSE);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            for (std::size_t i = 0; i < transforms.size(); ++i) {
                program_->setUniformValue("view", viewMatrix() * transforms[i]);
                glFrontFace(i == 1 || i == 2 ? GL_CW : GL_CCW);
                glDrawElements(GL_TRIANGLES, indexCount_, GL_UNSIGNED_INT, nullptr);
            }
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glDepthMask(GL_TRUE);
            glDisable(GL_BLEND);
        }
        vao_.release();
        glFrontFace(GL_CCW);
        program_->setUniformValue("view", viewMatrix());
    }
    if ((interfaceVisible_ && (interfaceVertexCount_ > 0 || normalVertexCount_ > 0)) ||
        (acousticGuidesVisible_ && guideVertexCount_ > 0)) {
        program_->setUniformValue("wirePass", 0);
        program_->setUniformValue("orientationPass", 0);
        program_->setUniformValue("unlit", 1);
        glDepthMask(GL_FALSE);
        interfaceVao_.bind();
        if (interfaceVisible_) {
            for (const auto &transform : symmetryTransforms()) {
                program_->setUniformValue("view", viewMatrix() * transform);
                glDrawArrays(GL_LINES, 0, interfaceVertexCount_);
                if (interfaceNormalsVisible_)
                    glDrawArrays(GL_LINES, interfaceVertexCount_, normalVertexCount_);
            }
            program_->setUniformValue("view", viewMatrix());
        }
        if (acousticGuidesVisible_)
            glDrawArrays(GL_LINES, interfaceVertexCount_ + normalVertexCount_, guideVertexCount_);
        interfaceVao_.release();
        glDepthMask(GL_TRUE);
    }
    program_->release();
    if (fieldVisible_ && field_ && fieldIndexCount_ > 0 && fieldProgram_) {
        const auto limits = pressureFieldRange(fieldAppearance_);
        fieldProgram_->bind();
        fieldProgram_->setUniformValue("view", viewMatrix());
        fieldProgram_->setUniformValue("projection", projectionMatrix());
        fieldProgram_->setUniformValue("quantity", fieldAppearance_.quantity);
        const auto shaderLimits =
            fieldAppearance_.quantity == 2 ? std::array<double, 2>{-1, 1} : limits;
        fieldProgram_->setUniformValue("limits", QVector2D(static_cast<float>(shaderLimits[0]),
                                                           static_cast<float>(shaderLimits[1])));
        const float bandStep = fieldAppearance_.quantity == 0   ? 3.f
                               : fieldAppearance_.quantity == 1 ? 36.f
                                                                : .2f;
        fieldProgram_->setUniformValue("bandStep", bandStep);
        std::array<QVector3D, 10> colors;
        for (std::size_t i = 0; i < colors.size(); ++i) {
            const auto &color = levelPalette()[i];
            colors[i] = QVector3D(color.redF(), color.greenF(), color.blueF());
        }
        fieldProgram_->setUniformValueArray("levelPalette", colors.data(),
                                            static_cast<int>(colors.size()));
        fieldProgram_->setUniformValue(
            "levelOffsetDb",
            static_cast<float>(20 * (std::log10(fieldPressureScale_) - std::log10(20e-6))));
        const double realLimit = 20e-6 * std::pow(10., fieldAppearance_.maximumDb / 20);
        fieldProgram_->setUniformValue(
            "realScale", static_cast<float>(std::exp(std::clamp(
                             std::log(fieldPressureScale_) - std::log(realLimit), -69., 69.))));
        fieldProgram_->setUniformValue(
            "phase", static_cast<float>(fieldAppearance_.phaseDegrees * std::numbers::pi / 180));
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);
        glDisable(GL_BLEND);
        fieldVao_.bind();
        glDrawElements(GL_TRIANGLES, fieldIndexCount_, GL_UNSIGNED_INT, nullptr);
        fieldVao_.release();
        fieldProgram_->release();
    }
    glDisable(GL_DEPTH_TEST);
    paintOverlay();
    ++paintedFrames_;
}

void MeshViewport::paintOverlay() {
    QOpenGLPaintDevice device(QSize(static_cast<int>(width() * devicePixelRatioF()),
                                    static_cast<int>(height() * devicePixelRatioF())));
    device.setDevicePixelRatio(devicePixelRatioF());
    QPainter painter(&device);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setPen(QColor(55, 61, 66));
    painter.drawText(QRect(12, 8, width() - 24, 23), Qt::AlignLeft | Qt::AlignVCenter, label_);
    if (interfaceVisible_ && interfaceMesh_ && !interfaceMesh_->triangles.empty()) {
        painter.setPen(QColor(0, 140, 155));
        painter.drawText(QRect(12, 32, 100, 20), Qt::AlignLeft | Qt::AlignVCenter,
                         tr("Interfaces"));
        painter.setPen(QColor(55, 61, 66));
    }
    std::vector<QRectF> occupied;
    if (!label_.isEmpty()) {
        const auto textSize = painter.fontMetrics().size(Qt::TextSingleLine, label_);
        occupied.emplace_back(12, 8, std::max(0, std::min(width() - 24, textSize.width())), 23);
    }
    QStringList guideLines;
    if (acousticGuidesVisible_ && symmetryPlanes_) {
        const auto planeLabel = [&](Vec3 n) {
            if (std::abs(n.x) > .999 && std::abs(symmetryPlanes_->origin.x) < 1e-12)
                return n.x > 0 ? tr("x=0, x≥0") : tr("x=0, x≤0");
            if (std::abs(n.y) > .999 && std::abs(symmetryPlanes_->origin.y) < 1e-12)
                return n.y > 0 ? tr("y=0, y≥0") : tr("y=0, y≤0");
            return tr("beliebige Ebene");
        };
        guideLines << tr("Symmetrie violett: %1 / grün: %2; Pfeile zum Halbraum")
                          .arg(planeLabel(symmetryPlanes_->firstNormal),
                               planeLabel(symmetryPlanes_->secondNormal));
    }
    if (acousticGuidesVisible_ && exteriorBaffle_)
        guideLines << tr("Schallwand orange (unendlich, symbolisch): Pfeil zum Außenraum");
    if (interfaceVisible_ && interfaceNormalsVisible_)
        for (const auto &group : interfaceGroups_)
            guideLines << tr("%1: cyan +n → %2 / orange −n → %3")
                              .arg(group.name, group.firstDomain, group.secondDomain);
    if (interfaceVisible_ && interfaceNormalsVisible_ && interfaceMesh_ &&
        !interfaceMesh_->triangles.empty() && interfaceGroups_.empty())
        guideLines << tr("Interface: cyan +n / orange −n");
    const auto gizmo = overlay::gizmo(size(), painter.fontMetrics(), fieldVisible_ && bool(field_),
                                      !label_.isEmpty());
    if (gizmo)
        occupied.push_back(gizmo->reserved);
    int guideY = height() - 8 - int(guideLines.size()) * painter.fontMetrics().height();
    const int guideLeft = gizmo ? int(gizmo->reserved.right()) + 4 : 8;
    for (const auto &line : guideLines) {
        const QRectF lineRect(guideLeft, guideY, std::max(0, width() - guideLeft - 8),
                              painter.fontMetrics().height());
        painter.fillRect(lineRect, QColor(248, 250, 251, 220));
        painter.setPen(QColor(55, 61, 66));
        painter.drawText(
            lineRect, Qt::AlignLeft | Qt::AlignVCenter,
            painter.fontMetrics().elidedText(line, Qt::ElideRight, int(lineRect.width())));
        occupied.push_back(lineRect);
        guideY += painter.fontMetrics().height();
    }
    if (fieldVisible_ && field_) {
        const auto limits = pressureFieldRange(fieldAppearance_);
        const auto unit = fieldAppearance_.quantity == 0   ? tr("dB SPL · 3 dB")
                          : fieldAppearance_.quantity == 1 ? tr("Phase · 36°")
                                                           : tr("Pa RMS · 10 Stufen");
        const auto &font = painter.fontMetrics();
        const QRectF legend =
            overlay::fieldLegend(size(), font, unit, QString::number(limits[0], 'g', 5),
                                 QString::number(limits[1], 'g', 5));
        if (!legend.isEmpty()) {
            occupied.push_back(legend);
            const int barHeight = static_cast<int>(legend.height()) - font.height() - 12;
            const QRect bar(static_cast<int>(legend.left()) + 6,
                            static_cast<int>(legend.top()) + font.height() + 4, 16, barHeight);
            painter.fillRect(legend, QColor(248, 250, 251, 235));
            painter.drawText(
                QRectF(legend.left() + 5, legend.top(), legend.width() - 10, font.height()),
                Qt::AlignLeft | Qt::AlignVCenter,
                font.elidedText(unit, Qt::ElideRight, int(legend.width()) - 10));
            for (int i = 0; i < barHeight; ++i) {
                const double t = 1 - double(i) / std::max(1, barHeight - 1);
                painter.setPen(
                    pressureFieldColor(limits[0] + t * (limits[1] - limits[0]), fieldAppearance_));
                painter.drawLine(bar.left(), bar.top() + i, bar.right(), bar.top() + i);
            }
            painter.setPen(QColor(55, 61, 66));
            const double step = fieldAppearance_.quantity == 0   ? 3.
                                : fieldAppearance_.quantity == 1 ? 36.
                                                                 : (limits[1] - limits[0]) / 10.;
            // Mark actual band boundaries; phase's first boundary is -162°,
            // leaving the wrap in the middle of one colour band.
            const double first = fieldAppearance_.quantity == 1 ? -162. : limits[1] - step;
            for (double boundary = first;
                 fieldAppearance_.quantity == 1 ? boundary < limits[1] : boundary > limits[0];
                 boundary += fieldAppearance_.quantity == 1 ? step : -step) {
                const double t = (boundary - limits[0]) / (limits[1] - limits[0]);
                const int y = bar.bottom() - int(t * (barHeight - 1));
                painter.drawLine(bar.right() + 1, y, bar.right() + 3, y);
            }
            const int intervals = barHeight >= 5 * font.height()   ? 4
                                  : barHeight >= 3 * font.height() ? 2
                                                                   : 1;
            for (int i = 0; i <= intervals; ++i) {
                const double t = double(i) / intervals;
                const int y = std::clamp(
                    bar.bottom() - int(t * (barHeight - 1)) - font.height() / 2,
                    int(legend.top()) + font.height(), int(legend.bottom()) - font.height());
                painter.drawText(
                    QRectF(bar.right() + 5, y, legend.right() - bar.right() - 7, font.height()),
                    Qt::AlignLeft | Qt::AlignVCenter,
                    QString::number(limits[0] + t * (limits[1] - limits[0]), 'g', 5));
            }
        }
        const auto context = tr("Schallfeld %1 · %2 Hz · Rev. %3%4")
                                 .arg(field_->vertical ? "V" : "H")
                                 .arg(field_->frequency, 0, 'g', 8)
                                 .arg(field_->revision)
                                 .arg(field_->current ? QString{} : tr(" · veraltet"));
        const int contextWidth =
            std::max(0, int((legend.isEmpty() ? width() - 10 : legend.left()) - 14));
        if (contextWidth > 0) {
            const QRect contextRect(10, 32, contextWidth, font.height() + 6);
            occupied.push_back(contextRect);
            painter.fillRect(contextRect, QColor(248, 250, 251, 225));
            painter.drawText(
                contextRect.adjusted(4, 0, -4, 0), Qt::AlignLeft | Qt::AlignVCenter,
                font.elidedText(context, Qt::ElideRight, std::max(0, contextWidth - 8)));
        }
    } else if (colorMode_ == MeshColorMode::orientation && mesh_) {
        const int left = width() - 190;
        painter.drawText(QRect(left, 32, 175, 22), Qt::AlignLeft, tr("Flächennormalen"));
        painter.fillRect(QRect(left, 59, 12, 12), QColor::fromRgbF(.20f, .66f, .39f));
        painter.drawText(QRect(left + 20, 54, 155, 22), Qt::AlignLeft, tr("Vorderseite (+n)"));
        painter.fillRect(QRect(left, 84, 12, 12), QColor::fromRgbF(.88f, .29f, .24f));
        painter.drawText(QRect(left + 20, 79, 155, 22), Qt::AlignLeft, tr("Rückseite (−n)"));
    } else if (usesMetricColors(colorMode_) && mesh_) {
        const auto title =
            colorMode_ == MeshColorMode::elementSize    ? tr("Mittlere Dreieckskante [mm]")
            : colorMode_ == MeshColorMode::minimumAngle ? tr("Kleinster Dreieckswinkel [°]")
                                                        : tr("Längste Kante / kleinste Höhe [1]");
        const QRect bar(width() - 180, 58, 150, 12);
        painter.drawText(QRect(width() - 230, 32, 215, 22), Qt::AlignRight, title);
        for (int x = 0; x < bar.width(); ++x) {
            double fraction = static_cast<double>(x) / (bar.width() - 1);
            if (colorMode_ == MeshColorMode::minimumAngle)
                fraction = 1 - fraction;
            painter.setPen(QColor::fromHsvF(static_cast<float>((1 - fraction) * .62), .75f, .9f));
            painter.drawLine(bar.left() + x, bar.top(), bar.left() + x, bar.bottom());
        }
        painter.setPen(QColor(55, 61, 66));
        painter.drawText(QRect(bar.left(), bar.bottom() + 3, 75, 20), Qt::AlignLeft,
                         QString::number(colorMinimum_, 'f', 2));
        painter.drawText(QRect(bar.left() + 75, bar.bottom() + 3, 75, 20), Qt::AlignRight,
                         QString::number(colorMaximum_, 'f', 2));
    }
    if (!rendererError_.isEmpty())
        painter.drawText(rect().adjusted(15, 40, -15, -15), Qt::AlignCenter | Qt::TextWordWrap,
                         rendererError_);
    if (observationFrame_ && mesh_ && !mesh_->vertices.empty()) {
        const auto &frame = *observationFrame_;
        const auto transform = projectionMatrix() * viewMatrix();
        const auto project = [&](Vec3 point) -> std::optional<QPointF> {
            const auto clip =
                transform * QVector4D(static_cast<float>((point.x - center_.x) / extent_),
                                      static_cast<float>((point.y - center_.y) / extent_),
                                      static_cast<float>((point.z - center_.z) / extent_), 1);
            if (clip.w() <= 0 || clip.z() < -clip.w() || clip.z() > clip.w())
                return {};
            const auto ndc = clip.toVector3DAffine();
            return QPointF((ndc.x() + 1) * .5 * width(), (1 - ndc.y()) * .5 * height());
        };
        if (const auto start = project(frame.origin)) {
            const std::array<Vec3, 3> directions{frame.forward, frame.horizontal,
                                                 verticalDirection(frame)};
            const std::array<QColor, 3> colors{QColor(180, 75, 18), QColor(25, 83, 165),
                                               QColor(30, 119, 56)};
            const std::array<QString, 3> labels{tr("0°"), tr("H +90°"), tr("V +90°")};
            // Arrow length is a view aid; the frame's origin and directions stay in world
            // coordinates.
            const double length = extent_ * .22;
            for (std::size_t i = 0; i < directions.size(); ++i) {
                const auto &d = directions[i];
                const auto end =
                    project({frame.origin.x + length * d.x, frame.origin.y + length * d.y,
                             frame.origin.z + length * d.z});
                if (!end)
                    continue;
                painter.setPen(QPen(colors[i], 2));
                painter.drawLine(*start, *end);
                const auto delta = *end - *start;
                const double pixels = std::hypot(delta.x(), delta.y());
                if (pixels > 8) {
                    const auto unit = delta / pixels;
                    const QPointF side(-unit.y(), unit.x());
                    painter.drawLine(*end, *end - unit * 8 + side * 3);
                    painter.drawLine(*end, *end - unit * 8 - side * 3);
                }
                const auto textSize = painter.fontMetrics().size(Qt::TextSingleLine, labels[i]);
                const auto label =
                    overlay::label(*end, QSizeF(textSize.width() + 6, textSize.height() + 4),
                                   QRectF(rect()).adjusted(2, 2, -2, -2), occupied);
                if (label) {
                    occupied.push_back(*label);
                    painter.fillRect(*label, QColor(240, 243, 245, 220));
                    painter.drawText(*label, Qt::AlignCenter, labels[i]);
                }
            }
        }
    }
    QMatrix4x4 rotation;
    rotation.rotate(static_cast<float>(pitch_ * 180 / std::numbers::pi), 1, 0, 0);
    rotation.rotate(static_cast<float>(yaw_ * 180 / std::numbers::pi), 0, 1, 0);
    const std::array<QVector3D, 3> axes{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    const std::array<QColor, 3> colors{
        {QColor(178, 52, 48), QColor(48, 130, 60), QColor(35, 92, 180)}};
    const std::array<QString, 3> names{{"X", "Y", "Z"}};
    if (gizmo) {
        for (std::size_t i = 0; i < axes.size(); ++i) {
            const auto direction = rotation.mapVector(axes[i]);
            const QPointF end(gizmo->origin.x() + gizmo->arm * direction.x(),
                              gizmo->origin.y() - gizmo->arm * direction.y());
            painter.setPen(QPen(colors[i], 1.5));
            painter.drawLine(gizmo->origin, end);
            painter.drawText(end + QPointF(3, -3), names[i]);
        }
    }
}

void MeshViewport::setWireframe(bool enabled) {
    wireframe_ = enabled;
    update();
}
void MeshViewport::setBackgroundGrid(bool enabled) {
    backgroundGrid_ = enabled;
    update();
}
void MeshViewport::setObservationFrame(std::optional<AcousticFrame> frame) {
    if (frame)
        frame = normalizedAcousticFrame(*frame);
    if (observationFrame_ == frame)
        return;
    observationFrame_ = frame;
    update();
}
void MeshViewport::setColorMode(MeshColorMode mode) {
    if (colorMode_ == mode)
        return;
    colorMode_ = mode;
    uploadPending_ = true;
    update();
}
void MeshViewport::setFaceColors(FaceColors colors) {
    if (faceColors_ == colors)
        return;
    faceColors_ = std::move(colors);
    if (colorMode_ == MeshColorMode::faceColors) {
        uploadPending_ = true;
        update();
    }
}
void MeshViewport::fitView() {
    automaticFit_ = true;
    pan_ = {};
    updateSceneBounds();
    updateFit();
    update();
}
void MeshViewport::viewFront(Vec3 forward) {
    if (!std::isfinite(forward.x) || !std::isfinite(forward.y) || !std::isfinite(forward.z))
        throw std::invalid_argument("Vorderansicht: Die Blickachse muss endlich sein.");
    const double scale = std::max({std::abs(forward.x), std::abs(forward.y), std::abs(forward.z)});
    if (scale == 0)
        throw std::invalid_argument("Vorderansicht: Die Blickachse darf nicht null sein.");
    forward = {forward.x / scale, forward.y / scale, forward.z / scale};
    // The existing Rx(pitch)*Ry(yaw) camera maps this physical forward axis to
    // camera +Z. Positive forward therefore faces the observer, including poles.
    yaw_ = -std::atan2(forward.x, forward.z);
    pitch_ = std::atan2(forward.y, std::hypot(forward.x, forward.z));
    fitView();
}
void MeshViewport::resizeGL(int, int) {
    if (automaticFit_)
        updateFit();
}
void MeshViewport::updateFit() {
    if ((!mesh_ || mesh_->vertices.empty()) && !field_)
        return;
    QMatrix4x4 rotation;
    rotation.rotate(static_cast<float>(pitch_ * 180 / std::numbers::pi), 1, 0, 0);
    rotation.rotate(static_cast<float>(yaw_ * 180 / std::numbers::pi), 0, 1, 0);
    const double aspect = static_cast<double>(std::max(1, width())) / std::max(1, height());
    // A projected box encloses every triangle. Depth matters: an orthographic
    // extent alone cannot bound the perspective silhouette after rotation.
    constexpr double occupancy = .84;
    double fit = std::numeric_limits<double>::infinity();
    for (const int x : {-1, 1})
        for (const int y : {-1, 1})
            for (const int z : {-1, 1}) {
                const auto p =
                    rotation.map(QVector3D(static_cast<float>(x * normalizedHalfSize_.x),
                                           static_cast<float>(y * normalizedHalfSize_.y),
                                           static_cast<float>(z * normalizedHalfSize_.z)));
                const double available =
                    occupancy * std::tan(halfFieldOfView) * (cameraDistance - p.z());
                if (p.x() != 0)
                    fit = std::min(fit, available * aspect / std::abs(p.x()));
                if (p.y() != 0)
                    fit = std::min(fit, available / std::abs(p.y()));
            }
    if (std::isfinite(fit) && fit > 0)
        zoom_ = fittedZoom_ = fit;
}
void MeshViewport::setLabel(QString label) {
    label_ = std::move(label);
    update();
}
void MeshViewport::setHiddenFaces(std::vector<std::uint32_t> faces) {
    std::sort(faces.begin(), faces.end());
    faces.erase(std::unique(faces.begin(), faces.end()), faces.end());
    if (faces == hiddenFaces_)
        return;
    hiddenFaces_ = std::move(faces);
    uploadPending_ = true;
    update();
}

void MeshViewport::setSelectedFaces(std::vector<std::uint32_t> faces) {
    std::sort(faces.begin(), faces.end());
    faces.erase(std::unique(faces.begin(), faces.end()), faces.end());
    if (faces == selectedFaces_)
        return;
    selectedFaces_ = std::move(faces);
    selectionUploadPending_ = true;
    update();
}
const std::vector<std::uint32_t> &MeshViewport::selectedFaces() const noexcept {
    return selectedFaces_;
}
QString MeshViewport::rendererError() const {
    return rendererError_;
}

std::optional<std::uint32_t> MeshViewport::pickFace(QPointF position) const {
    if (!mesh_ || width() <= 0 || height() <= 0 || !QRectF(rect()).contains(position))
        return {};
    bool invertible = false;
    const auto inverse = (projectionMatrix() * viewMatrix()).inverted(&invertible);
    if (!invertible)
        return {};
    const float x = static_cast<float>(2 * position.x() / width() - 1),
                y = static_cast<float>(1 - 2 * position.y() / height());
    const auto nearPoint = (inverse * QVector4D(x, y, -1, 1)).toVector3DAffine();
    const auto farPoint = (inverse * QVector4D(x, y, 1, 1)).toVector3DAffine();
    const auto direction = (farPoint - nearPoint).normalized();
    const Vec3 origin{nearPoint.x(), nearPoint.y(), nearPoint.z()},
        ray{direction.x(), direction.y(), direction.z()};
    const auto normalize = [this](Vec3 p) {
        return Vec3{(p.x - center_.x) / extent_, (p.y - center_.y) / extent_,
                    (p.z - center_.z) / extent_};
    };
    double nearest = std::numeric_limits<double>::infinity();
    std::optional<std::uint32_t> selected;
    // Double-sided Moller–Trumbore intersection, in the same normalized coordinates as the VBO.
    for (const auto &transform : symmetryTransforms()) {
        const auto instancePoint = [&](Vec3 point) {
            const auto p = normalize(point);
            const auto mapped = transform * QVector4D(float(p.x), float(p.y), float(p.z), 1);
            return Vec3{mapped.x(), mapped.y(), mapped.z()};
        };
        for (const auto &triangle : mesh_->triangles) {
            if (std::binary_search(hiddenFaces_.begin(), hiddenFaces_.end(), triangle.faceId))
                continue;
            const auto a = instancePoint(mesh_->vertices[triangle.vertices[0]]);
            const auto e1 = subtract(instancePoint(mesh_->vertices[triangle.vertices[1]]), a);
            const auto e2 = subtract(instancePoint(mesh_->vertices[triangle.vertices[2]]), a);
            const auto p = cross(ray, e2);
            const double determinant = dot(e1, p);
            const double scale = std::sqrt(dot(e1, e1) * dot(e2, e2));
            if (std::abs(determinant) <= scale * 1e-12)
                continue;
            const auto t = subtract(origin, a);
            const double u = dot(t, p) / determinant;
            if (u < 0 || u > 1)
                continue;
            const auto q = cross(t, e1);
            const double v = dot(ray, q) / determinant;
            if (v < 0 || u + v > 1)
                continue;
            const double distance = dot(e2, q) / determinant;
            if (distance >= 0 && distance < nearest) {
                nearest = distance;
                selected = triangle.faceId;
            }
        }
    }
    return selected;
}

void MeshViewport::setFieldProbeEnabled(bool enabled) {
    fieldProbeEnabled_ = enabled;
    if (!enabled)
        QToolTip::hideText();
}

std::optional<std::size_t> MeshViewport::pickFieldSample(QPointF position) const {
    if (!fieldVisible_ || !field_ || width() <= 0 || height() <= 0 ||
        !QRectF(rect()).contains(position))
        return {};
    const auto indices = pressureFieldTriangleIndices(*field_);
    std::vector<bool> rendered(field_->positions.size(), false);
    for (const auto index : indices)
        rendered[index] = true;
    const auto transform = projectionMatrix() * viewMatrix();
    double closest = 12. * 12.;
    std::optional<std::size_t> selected;
    for (std::size_t i = 0; i < field_->positions.size(); ++i) {
        if (!field_->valid[i] || !rendered[i])
            continue;
        const auto &p = field_->positions[i];
        const auto clip = transform * QVector4D(static_cast<float>((p.x - center_.x) / extent_),
                                                static_cast<float>((p.y - center_.y) / extent_),
                                                static_cast<float>((p.z - center_.z) / extent_), 1);
        if (clip.w() <= 0 || clip.z() < -clip.w() || clip.z() > clip.w())
            continue;
        const auto ndc = clip.toVector3DAffine();
        const QPointF projected((ndc.x() + 1.) * .5 * width(), (1. - ndc.y()) * .5 * height());
        const auto delta = projected - position;
        const double distance = delta.x() * delta.x() + delta.y() * delta.y();
        if (distance <= closest) {
            closest = distance;
            selected = i;
        }
    }
    if (!selected || !mesh_)
        return selected;
    // A projected point hidden behind CAD must not be reported as a visible
    // measurement. Compare the camera ray against the displayed triangles.
    const auto &world = field_->positions[*selected];
    const Vec3 target{(world.x - center_.x) / extent_, (world.y - center_.y) / extent_,
                      (world.z - center_.z) / extent_};
    const auto targetClip =
        transform * QVector4D(float(target.x), float(target.y), float(target.z), 1);
    const auto ndc = targetClip.toVector3DAffine();
    bool invertible = false;
    const auto inverse = transform.inverted(&invertible);
    if (!invertible)
        return {};
    const auto near = (inverse * QVector4D(ndc.x(), ndc.y(), -1, 1)).toVector3DAffine();
    const Vec3 origin{near.x(), near.y(), near.z()};
    const auto difference = subtract(target, origin);
    const double targetDistance = std::sqrt(dot(difference, difference));
    if (targetDistance <= 0)
        return {};
    const Vec3 ray{difference.x / targetDistance, difference.y / targetDistance,
                   difference.z / targetDistance};
    const auto normalize = [this](Vec3 p) {
        return Vec3{(p.x - center_.x) / extent_, (p.y - center_.y) / extent_,
                    (p.z - center_.z) / extent_};
    };
    for (const auto &transform : symmetryTransforms()) {
        const auto instancePoint = [&](Vec3 point) {
            const auto p = normalize(point);
            const auto mapped = transform * QVector4D(float(p.x), float(p.y), float(p.z), 1);
            return Vec3{mapped.x(), mapped.y(), mapped.z()};
        };
        for (const auto &triangle : mesh_->triangles) {
            if (std::binary_search(hiddenFaces_.begin(), hiddenFaces_.end(), triangle.faceId))
                continue;
            const auto a = instancePoint(mesh_->vertices[triangle.vertices[0]]);
            const auto e1 = subtract(instancePoint(mesh_->vertices[triangle.vertices[1]]), a);
            const auto e2 = subtract(instancePoint(mesh_->vertices[triangle.vertices[2]]), a);
            const auto crossRay = cross(ray, e2);
            const double determinant = dot(e1, crossRay);
            const double scale = std::sqrt(dot(e1, e1) * dot(e2, e2));
            if (std::abs(determinant) <= scale * 1e-12)
                continue;
            const auto t = subtract(origin, a);
            const double u = dot(t, crossRay) / determinant;
            if (u < 0 || u > 1)
                continue;
            const auto q = cross(t, e1);
            const double v = dot(ray, q) / determinant;
            const double distance = dot(e2, q) / determinant;
            if (v >= 0 && u + v <= 1 && distance >= 0 && distance < targetDistance - 1e-6)
                return {};
        }
    }
    return selected;
}

void MeshViewport::mousePressEvent(QMouseEvent *event) {
    lastMouse_ = event->position().toPoint();
    pressMouse_ = lastMouse_;
    dragStarted_ = false;
}
void MeshViewport::mouseMoveEvent(QMouseEvent *event) {
    if (!(event->buttons() & (Qt::LeftButton | Qt::MiddleButton)))
        return;
    dragStarted_ |= (event->position().toPoint() - pressMouse_).manhattanLength() > 4;
    if (!dragStarted_)
        return;
    automaticFit_ = false;
    const QPoint delta = event->position().toPoint() - lastMouse_;
    if ((event->modifiers() & Qt::ShiftModifier) || (event->buttons() & Qt::MiddleButton))
        pan_ += delta;
    else {
        yaw_ += delta.x() * .008;
        pitch_ = std::clamp(pitch_ + delta.y() * .008, -1.5, 1.5);
    }
    lastMouse_ = event->position().toPoint();
    update();
}
void MeshViewport::mouseReleaseEvent(QMouseEvent *event) {
    if (event->button() != Qt::LeftButton || dragStarted_ ||
        (event->position().toPoint() - pressMouse_).manhattanLength() > 4 ||
        (event->modifiers() & (Qt::AltModifier | Qt::MetaModifier | Qt::ShiftModifier)))
        return;
    if (fieldProbeEnabled_ && !(event->modifiers() & Qt::ControlModifier)) {
        const auto sample = pickFieldSample(event->position());
        if (!sample) {
            QToolTip::showText(event->globalPosition().toPoint(),
                               tr("Kein gültiger Feldrasterpunkt in der Nähe."), this);
            return;
        }
        const auto &p = field_->positions[*sample];
        const auto clip = (projectionMatrix() * viewMatrix()) *
                          QVector4D(static_cast<float>((p.x - center_.x) / extent_),
                                    static_cast<float>((p.y - center_.y) / extent_),
                                    static_cast<float>((p.z - center_.z) / extent_), 1);
        const auto ndc = clip.toVector3DAffine();
        const QPointF projected((ndc.x() + 1.) * .5 * width(), (1. - ndc.y()) * .5 * height());
        const auto delta = projected - event->position();
        const double pixelDistance = std::hypot(delta.x(), delta.y());
        const auto pressure = field_->pressure[*sample];
        const double magnitude = std::abs(pressure);
        const QString spl = magnitude > 0
                                ? tr("%1 dB SPL").arg(20 * std::log10(magnitude / 20e-6), 0, 'g', 5)
                                : tr("−∞ dB SPL");
        const QString phase =
            magnitude > 0 ? tr("%1°").arg(std::arg(pressure) * 180 / std::numbers::pi, 0, 'g', 5)
                          : tr("nicht definiert");
        const QString message =
            tr("Nächster projizierter Rasterpunkt (%1 px; keine Interpolation)\n"
               "x=%2 m · y=%3 m · z=%4 m\n%5 · Phase %6")
                .arg(pixelDistance, 0, 'f', 1)
                .arg(p.x, 0, 'g', 5)
                .arg(p.y, 0, 'g', 5)
                .arg(p.z, 0, 'g', 5)
                .arg(spl, phase);
        QToolTip::showText(event->globalPosition().toPoint(), message, this);
        return;
    }
    const auto face = pickFace(event->position());
    auto selected = selectedFaces_;
    if (event->modifiers() & Qt::ControlModifier) {
        if (face) {
            const auto found = std::lower_bound(selected.begin(), selected.end(), *face);
            if (found != selected.end() && *found == *face)
                selected.erase(found);
            else
                selected.insert(found, *face);
        }
    } else
        selected = face ? std::vector<std::uint32_t>{*face} : std::vector<std::uint32_t>{};
    setSelectedFaces(std::move(selected));
    if (onFacePicked)
        onFacePicked(face, event->modifiers());
}
void MeshViewport::wheelEvent(QWheelEvent *event) {
    automaticFit_ = false;
    zoom_ = std::clamp(zoom_ * std::exp(event->angleDelta().y() * .001), fittedZoom_ * .2,
                       fittedZoom_ * 8.0);
    update();
    event->accept();
}
void MeshViewport::mouseDoubleClickEvent(QMouseEvent *) {
    fitView();
}

RenderBenchmark MeshViewport::benchmarkFrames(std::size_t measuredFrames,
                                              std::size_t warmupFrames) {
    if (!initialized_ || !isVisible() || !rendererError_.isEmpty())
        throw std::runtime_error("OpenGL-Viewport ist nicht messbereit.");
    if (measuredFrames < 2 || measuredFrames > 10000 || warmupFrames > 10000)
        throw std::invalid_argument("Ungültige Benchmark-Länge.");
    const double originalYaw = yaw_;
    std::vector<double> times;
    times.reserve(measuredFrames);
    for (std::size_t i = 0; i < warmupFrames + measuredFrames; ++i) {
        yaw_ = originalYaw + static_cast<double>(i) * .003;
        const auto before = paintedFrames_;
        QElapsedTimer timer;
        timer.start();
        makeCurrent();
        paintGL();
        glFinish();
        const auto error = glGetError();
        doneCurrent();
        if (error != GL_NO_ERROR) {
            yaw_ = originalYaw;
            throw std::runtime_error("OpenGL meldet einen Fehler beim Renderbenchmark: " +
                                     std::to_string(error));
        }
        const double milliseconds = static_cast<double>(timer.nsecsElapsed()) * 1e-6;
        if (paintedFrames_ == before) {
            yaw_ = originalYaw;
            throw std::runtime_error("Benchmark hat keinen Frame gezeichnet.");
        }
        if (i >= warmupFrames)
            times.push_back(milliseconds);
    }
    yaw_ = originalYaw;
    update();
    std::sort(times.begin(), times.end());
    const auto percentile = [&times](double fraction) {
        return times[static_cast<std::size_t>(
                         std::ceil(fraction * static_cast<double>(times.size()))) -
                     1];
    };
    return {renderer_,
            driverVersion_,
            mesh_ ? mesh_->triangles.size() : 0,
            measuredFrames,
            QSize(static_cast<int>(width() * devicePixelRatioF()),
                  static_cast<int>(height() * devicePixelRatioF())),
            percentile(.5),
            percentile(.95)};
}

} // namespace hornsculpt
