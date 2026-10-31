#include <QApplication>
#include <QCommandLineParser>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QMouseEvent>
#include <QTimer>
#include <QWheelEvent>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>
#include <hornsculpt/draft_geometry.hpp>
#include <hornsculpt/mesh_viewport.hpp>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>

namespace {
void click(hornsculpt::MeshViewport &viewport, QPointF position,
           Qt::KeyboardModifiers modifiers = Qt::NoModifier) {
    QMouseEvent press(QEvent::MouseButtonPress, position, viewport.mapToGlobal(position.toPoint()),
                      Qt::LeftButton, Qt::LeftButton, modifiers);
    QMouseEvent release(QEvent::MouseButtonRelease, position,
                        viewport.mapToGlobal(position.toPoint()), Qt::LeftButton, Qt::NoButton,
                        modifiers);
    QApplication::sendEvent(&viewport, &press);
    QApplication::sendEvent(&viewport, &release);
}
QJsonObject measure(hornsculpt::MeshViewport &viewport, std::size_t frames) {
    const auto result = viewport.benchmarkFrames(frames);
    return {{"renderer", result.renderer},
            {"driverVersion", result.driverVersion},
            {"triangles", static_cast<qint64>(result.triangles)},
            {"frames", static_cast<qint64>(result.measuredFrames)},
            {"framebufferWidth", result.framebufferSize.width()},
            {"framebufferHeight", result.framebufferSize.height()},
            {"medianMs", result.medianMilliseconds},
            {"p95Ms", result.p95Milliseconds}};
}
void drag(hornsculpt::MeshViewport &viewport, QPointF delta,
          Qt::KeyboardModifiers modifiers = Qt::NoModifier) {
    const QPointF start(viewport.width() * .5, viewport.height() * .5), end = start + delta;
    QMouseEvent press(QEvent::MouseButtonPress, start, viewport.mapToGlobal(start.toPoint()),
                      Qt::LeftButton, Qt::LeftButton, modifiers);
    QMouseEvent move(QEvent::MouseMove, end, viewport.mapToGlobal(end.toPoint()), Qt::NoButton,
                     Qt::LeftButton, modifiers);
    QMouseEvent release(QEvent::MouseButtonRelease, end, viewport.mapToGlobal(end.toPoint()),
                        Qt::LeftButton, Qt::NoButton, modifiers);
    QApplication::sendEvent(&viewport, &press);
    QApplication::sendEvent(&viewport, &move);
    QApplication::sendEvent(&viewport, &release);
}
std::shared_ptr<hornsculpt::SurfaceMesh> box(hornsculpt::Vec3 size) {
    auto mesh = std::make_shared<hornsculpt::SurfaceMesh>();
    for (const int z : {-1, 1})
        for (const int y : {-1, 1})
            for (const int x : {-1, 1})
                mesh->vertices.push_back({x * size.x * .5, y * size.y * .5, z * size.z * .5});
    for (const auto indices : std::array<std::array<std::uint32_t, 3>, 12>{{{0, 2, 1},
                                                                            {1, 2, 3},
                                                                            {4, 5, 6},
                                                                            {5, 7, 6},
                                                                            {0, 1, 4},
                                                                            {1, 5, 4},
                                                                            {2, 6, 3},
                                                                            {3, 6, 7},
                                                                            {0, 4, 2},
                                                                            {2, 4, 6},
                                                                            {1, 3, 5},
                                                                            {3, 7, 5}}})
        mesh->triangles.push_back({indices, 1, hornsculpt::BoundaryCondition::driven});
    return mesh;
}
QRect silhouette(hornsculpt::MeshViewport &viewport) {
    viewport.repaint();
    const auto image = viewport.grabFramebuffer();
    QRect bounds;
    const int axes = static_cast<int>(100 * viewport.devicePixelRatioF());
    for (int y = 0; y < image.height(); ++y)
        for (int x = 0; x < image.width(); ++x) {
            if (x < axes && y > image.height() - axes)
                continue;
            const auto color = image.pixelColor(x, y);
            if (color.blue() > color.red() * 1.7 && color.blue() > color.green() * 1.25)
                bounds |= QRect(x, y, 1, 1);
        }
    if (bounds.width() < 10 || bounds.height() < 10)
        throw std::runtime_error("Einpassung: Modellsilhouette fehlt im OpenGL-Bild.");
    return bounds;
}
void checkPerspectiveFit(hornsculpt::MeshViewport &viewport) {
    viewport.setColorMode(hornsculpt::MeshColorMode::boundary);
    viewport.setWireframe(false);
    for (const auto size : {hornsculpt::Vec3{.2, 4, .3}, hornsculpt::Vec3{5, .2, 3}}) {
        auto mesh = box(size);
        viewport.setMesh(mesh);
        drag(viewport, {37, -16});
        viewport.fitView();
        // Resize after fitting exercises the real QOpenGLWidget dock-resize path.
        for (const auto resolution : {QSize(420, 900), QSize(1200, 300)}) {
            viewport.resize(resolution);
            QApplication::processEvents();
            const auto bounds = silhouette(viewport);
            const auto framebuffer = viewport.grabFramebuffer().size();
            const QRect safe(static_cast<int>(framebuffer.width() * .06),
                             static_cast<int>(framebuffer.height() * .06),
                             static_cast<int>(framebuffer.width() * .88),
                             static_cast<int>(framebuffer.height() * .88));
            if (!safe.contains(bounds))
                throw std::runtime_error("Perspektivische Einpassung schneidet das Modell ab.");
        }
        const QPointF center(viewport.width() * .5, viewport.height() * .5);
        QWheelEvent wheel(center, viewport.mapToGlobal(center.toPoint()), {}, {0, -360},
                          Qt::NoButton, Qt::NoModifier, Qt::NoScrollPhase, false);
        QApplication::sendEvent(&viewport, &wheel);
        drag(viewport, {31, 18}, Qt::ShiftModifier);
        const auto manual = silhouette(viewport);
        viewport.setMesh(std::make_shared<hornsculpt::SurfaceMesh>(*mesh));
        if (silhouette(viewport) != manual)
            throw std::runtime_error("Live-Geometrie setzt die manuelle Kamera zurück.");
        viewport.resize(700, 600);
        QApplication::processEvents();
        viewport.resize(1200, 300);
        QApplication::processEvents();
        if (silhouette(viewport) != manual)
            throw std::runtime_error("Dock-Resize setzt die manuelle Kamera zurück.");
    }
}
void checkFaceColors(hornsculpt::MeshViewport &viewport, const hornsculpt::SurfaceMesh &mesh) {
    using hornsculpt::MeshColorMode;
    static_assert(static_cast<int>(MeshColorMode::boundary) == 0 &&
                  static_cast<int>(MeshColorMode::elementSize) == 1 &&
                  static_cast<int>(MeshColorMode::minimumAngle) == 2 &&
                  static_cast<int>(MeshColorMode::aspectRatio) == 3 &&
                  static_cast<int>(MeshColorMode::orientation) == 4 &&
                  static_cast<int>(MeshColorMode::faceColors) == 5);
    const auto originalMesh = mesh;
    const hornsculpt::FaceColors blue{{55, {210, 40, 180}}, {56, {40, 80, 220}}};
    const hornsculpt::FaceColors green{{55, {210, 40, 180}}, {56, {40, 200, 80}}};
    const QPointF center(viewport.width() * .5, viewport.height() * .5);
    const auto pixel = [&] {
        viewport.repaint();
        const auto image = viewport.grabFramebuffer();
        return image.pixelColor(image.width() / 2, image.height() / 2);
    };
    viewport.setWireframe(false);
    viewport.setSelectedFaces({});
    viewport.setColorMode(MeshColorMode::faceColors);
    viewport.setFaceColors(blue);
    const auto assigned = pixel();
    if (assigned.blue() < assigned.red() * 2 || assigned.blue() < assigned.green() * 1.5)
        throw std::runtime_error("Zugewiesene RGB-Flächenfarbe fehlt im Framebuffer.");
    viewport.setFaceColors({{55, {210, 40, 180}}, {56, {220, 50, 40}}});
    const auto changed = pixel();
    if (changed.red() < changed.green() * 2 || changed.red() < changed.blue() * 2)
        throw std::runtime_error("Änderung der Flächenfarbe erreicht den GPU-Puffer nicht.");
    viewport.setFaceColors({{55, {210, 40, 180}}});
    const auto neutral = pixel();
    if (std::abs(neutral.red() - neutral.green()) > 2 ||
        std::abs(neutral.red() - neutral.blue()) > 2)
        throw std::runtime_error("Nicht zugewiesene Fläche ist im Farbmodus nicht neutralgrau.");
    viewport.setFaceColors({{56, {0, 0, 0}}});
    const auto black = pixel();
    if (black.red() > 1 || black.green() > 1 || black.blue() > 1)
        throw std::runtime_error("RGB-Schwarz wurde als fehlende Flächenfarbe behandelt.");
    viewport.setFaceColors(green);
    viewport.setMesh(std::make_shared<hornsculpt::SurfaceMesh>(mesh));
    const auto replacement = pixel();
    if (replacement.green() < replacement.red() * 2 ||
        replacement.green() < replacement.blue() * 1.5)
        throw std::runtime_error("Farben stabiler Flächen-IDs gehen beim Netzaustausch verloren.");
    viewport.setSelectedFaces({56});
    const auto selected = pixel();
    if (selected.red() < selected.green() * 1.15 || selected.green() < selected.blue() * 1.5)
        throw std::runtime_error("Auswahlmarkierung überlagert die freie Flächenfarbe nicht.");
    viewport.setSelectedFaces({});
    viewport.setHiddenFaces({56});
    const auto revealed = pixel();
    if (revealed.red() < revealed.green() * 2 || revealed.blue() < revealed.green() * 2 ||
        viewport.pickFace(center) != 55)
        throw std::runtime_error("Versteckte farbige Fläche verdeckt Darstellung oder Picking.");
    viewport.setHiddenFaces({55, 56});
    (void)pixel();
    if (viewport.pickFace(center))
        throw std::runtime_error("Verborgene farbige Flächen bleiben auswählbar.");
    viewport.setHiddenFaces({});
    for (const auto mode :
         {MeshColorMode::boundary, MeshColorMode::elementSize, MeshColorMode::minimumAngle,
          MeshColorMode::aspectRatio, MeshColorMode::orientation}) {
        viewport.setColorMode(mode);
        viewport.setFaceColors(blue);
        const auto before = pixel();
        viewport.setFaceColors(green);
        if (pixel() != before)
            throw std::runtime_error(
                "Freie Farben verändern Rand-, Qualitäts- oder Normalenfarben.");
    }
    if (mesh != originalMesh)
        throw std::runtime_error("Flächenfarben verändern das physische Quellnetz.");
    viewport.setFaceColors({});
    viewport.setColorMode(MeshColorMode::boundary);
    viewport.setWireframe(true);
}
void checkFrontView(hornsculpt::MeshViewport &viewport) {
    auto mesh = box({2, 2, 2});
    for (std::size_t i = 0; i < mesh->triangles.size(); ++i)
        mesh->triangles[i].faceId = static_cast<std::uint32_t>(i / 2 + 1);
    const auto original = *mesh;
    viewport.setMesh(mesh);
    viewport.setHiddenFaces({});
    const QPointF center(viewport.width() * .5, viewport.height() * .5);
    for (const auto &[direction, expected] :
         std::array<std::pair<hornsculpt::Vec3, std::uint32_t>, 9>{{{{0, 0, 1}, 2},
                                                                    {{0, 0, -1}, 1},
                                                                    {{2, 0, 0}, 6},
                                                                    {{-1, 0, 0}, 5},
                                                                    {{0, 1, 0}, 4},
                                                                    {{0, -1, 0}, 3},
                                                                    {{1, 2, 3}, 2},
                                                                    {{3, 2, 1}, 6},
                                                                    {{-1, -3, 2}, 3}}}) {
        viewport.viewFront(direction);
        if (viewport.pickFace(center) != expected)
            throw std::runtime_error(
                "Vorderansicht bildet die physische Vorwärtsachse falsch auf den Bildschirm ab.");
    }
    viewport.viewFront({0, 0, 1});
    for (const auto invalid :
         {hornsculpt::Vec3{}, hornsculpt::Vec3{std::numeric_limits<double>::quiet_NaN(), 0, 1},
          hornsculpt::Vec3{0, std::numeric_limits<double>::infinity(), 1}}) {
        bool rejected = false;
        try {
            viewport.viewFront(invalid);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        if (!rejected || viewport.pickFace(center) != 2)
            throw std::runtime_error(
                "Ungültige Vorderachse wird akzeptiert oder verändert die bestehende Ansicht.");
    }
    if (*mesh != original)
        throw std::runtime_error("Vorderansicht verändert das physische Netz.");
}
void checkInterfaceOverlay(hornsculpt::MeshViewport &viewport, const QString &screenshot = {}) {
    auto source = box({2, 2, 2});
    viewport.setMesh(source);
    viewport.setHiddenFaces({});
    viewport.setWireframe(false);
    viewport.fitView();
    viewport.viewFront({0, 0, 1});
    const QPointF center(viewport.width() * .5, viewport.height() * .5);
    const auto picked = viewport.pickFace(center);
    if (!picked)
        throw std::runtime_error("Interface-Test: physische Fläche nicht auswählbar.");
    auto interfaceMesh = std::make_shared<hornsculpt::SurfaceMesh>();
    interfaceMesh->vertices = {{-.5, -.5, 1.01}, {.5, -.5, 1.01}, {0, .5, 1.01}};
    interfaceMesh->triangles = {{{0, 1, 2}, 999, hornsculpt::BoundaryCondition::wall}};
    viewport.setInterfaceMesh(interfaceMesh);
    viewport.setInterfaceVisible(true);
    viewport.repaint();
    if (viewport.pickFace(center) != picked)
        throw std::runtime_error("Interface-Overlay verändert die Flächenauswahl.");
    const auto visible = viewport.grabFramebuffer();
    if (!screenshot.isEmpty() && !visible.save(screenshot))
        throw std::runtime_error("Interface-Screenshot konnte nicht gespeichert werden.");
    std::size_t cyan = 0;
    const int top = static_cast<int>(60 * viewport.devicePixelRatioF());
    for (int y = top; y < visible.height(); ++y)
        for (int x = 0; x < visible.width(); ++x) {
            const auto c = visible.pixelColor(x, y);
            if (c.red() < 70 && c.green() > 110 && c.blue() > 120 && c.blue() > c.red() + 70)
                ++cyan;
        }
    if (cyan == 0)
        throw std::runtime_error("Interface-Drahtnetz fehlt im OpenGL-Bild.");
    viewport.setInterfaceVisible(false);
    viewport.repaint();
    if (viewport.pickFace(center) != picked)
        throw std::runtime_error("Interface-Sichtbarkeit verändert die Flächenauswahl.");
    const auto hidden = viewport.grabFramebuffer();
    std::size_t hiddenCyan = 0;
    for (int y = top; y < hidden.height(); ++y)
        for (int x = 0; x < hidden.width(); ++x) {
            const auto c = hidden.pixelColor(x, y);
            if (c.red() < 70 && c.green() > 110 && c.blue() > 120 && c.blue() > c.red() + 70)
                ++hiddenCyan;
        }
    if (hiddenCyan >= cyan)
        throw std::runtime_error("Interface-Drahtnetz bleibt sichtbar.");
    viewport.setInterfaceMesh({});
    viewport.setInterfaceVisible(true);
}
// Compare rendered contributions rather than a particular driver-dependent blend color.
std::size_t changedPixels(const QImage &a, const QImage &b) {
    if (a.size() != b.size())
        throw std::runtime_error("Framebuffergröße änderte sich während des Tests.");
    std::size_t changed = 0;
    for (int y = 0; y < a.height(); ++y)
        for (int x = 0; x < a.width(); ++x) {
            const auto left = a.pixelColor(x, y), right = b.pixelColor(x, y);
            if (std::max({std::abs(left.red() - right.red()),
                          std::abs(left.green() - right.green()),
                          std::abs(left.blue() - right.blue())}) > 8)
                ++changed;
        }
    return changed;
}
void checkAcousticGuides(hornsculpt::MeshViewport &viewport, const QString &screenshot) {
    using namespace hornsculpt;
    auto source = box({2, 2, 2});
    const auto original = *source;
    viewport.setMesh(source);
    viewport.setHiddenFaces({});
    viewport.setSelectedFaces({});
    viewport.setWireframe(false);
    viewport.setColorMode(MeshColorMode::boundary);
    viewport.viewFront({1, 1, 2});
    viewport.fitView();
    const QPointF center(viewport.width() * .5, viewport.height() * .5);
    const auto picked = viewport.pickFace(center);
    if (!picked)
        throw std::runtime_error("Hilfsebenen-Test: Modell nicht auswählbar.");
    const auto image = [&] {
        viewport.repaint();
        return viewport.grabFramebuffer();
    };
    const auto baseline = image();
    const OrthogonalSymmetryPlanes planes{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 1e-6};
    const ExteriorBaffle baffle{{0, 0, 0}, {0, 0, 1}, 1e-6};
    viewport.setExteriorBaffle(baffle);
    const auto baffleImage = image();
    if (changedPixels(baseline, baffleImage) < 100)
        throw std::runtime_error("Schallwandebene fehlt im OpenGL-Bild.");
    viewport.setExteriorBaffle({});
    viewport.setSymmetryPlanes(planes);
    const auto symmetryImage = image();
    if (changedPixels(baseline, symmetryImage) < 100)
        throw std::runtime_error("Symmetrieebenen fehlen im OpenGL-Bild.");
    // Moving one plane while keeping the other fixed must affect the projection.
    auto moved = planes;
    moved.secondNormal = {0, 0, 1};
    viewport.setSymmetryPlanes(moved);
    if (changedPixels(symmetryImage, image()) < 100)
        throw std::runtime_error("Die zweite Symmetrieebene ist nicht sichtbar.");
    moved = planes;
    moved.firstNormal = {0, 0, 1};
    viewport.setSymmetryPlanes(moved);
    if (changedPixels(symmetryImage, image()) < 100)
        throw std::runtime_error("Die erste Symmetrieebene ist nicht sichtbar.");
    viewport.setSymmetryPlanes(planes);
    viewport.setExteriorBaffle(baffle);
    const auto combined = image();
    if (!screenshot.isEmpty() && !combined.save(screenshot))
        throw std::runtime_error("Hilfsebenen-Screenshot konnte nicht gespeichert werden.");
    viewport.setAcousticGuidesVisible(false);
    if (changedPixels(baseline, image()) != 0)
        throw std::runtime_error(
            "Ausgeblendete Hilfsebenen bleiben sichtbar oder verändern die Kamera.");
    viewport.setAcousticGuidesVisible(true);
    viewport.setSymmetryPlanes({});
    viewport.setExteriorBaffle({});
    if (changedPixels(baseline, image()) != 0)
        throw std::runtime_error("Gelöschte Hilfsebenen hinterlassen Bildartefakte.");

    auto interfaceMesh = std::make_shared<SurfaceMesh>();
    interfaceMesh->vertices = {{-.6, -.6, 1.01}, {.6, -.6, 1.01}, {0, .6, 1.01}};
    interfaceMesh->triangles = {{{0, 1, 2}, 999, BoundaryCondition::wall}};
    const auto interfaceOriginal = *interfaceMesh;
    const std::vector<InterfaceOrientationLabels> labels{{{0}, "Interface", "Innen", "Außen"}};
    viewport.setInterfaceMesh(interfaceMesh, labels);
    viewport.setInterfaceVisible(true);
    viewport.setInterfaceNormalsVisible(false);
    const auto withoutNormals = image();
    viewport.setInterfaceNormalsVisible(true);
    const auto forward = image();
    if (changedPixels(withoutNormals, forward) < 5)
        throw std::runtime_error("Interface-Normale und Gebietsbezeichnungen fehlen.");
    auto reversed = std::make_shared<SurfaceMesh>(*interfaceMesh);
    std::swap(reversed->triangles[0].vertices[1], reversed->triangles[0].vertices[2]);
    viewport.setInterfaceMesh(reversed, labels);
    const auto backward = image();
    if (changedPixels(forward, backward) < 5)
        throw std::runtime_error("Interface-Normale folgt nicht der Dreiecksorientierung.");
    viewport.setInterfaceNormalsVisible(false);
    if (changedPixels(withoutNormals, image()) != 0)
        throw std::runtime_error("Ausgeblendete Interface-Normalen hinterlassen Artefakte.");
    viewport.setInterfaceMesh({});
    viewport.setInterfaceNormalsVisible(true);
    if (changedPixels(baseline, image()) != 0 || viewport.pickFace(center) != picked)
        throw std::runtime_error(
            "Gelöschtes Interface verändert Kamera, Bild oder Flächenauswahl.");
    if (*source != original || *interfaceMesh != interfaceOriginal)
        throw std::runtime_error("Akustische Hilfen verändern die Quellgeometrie.");
}
void checkPressureField(hornsculpt::MeshViewport &viewport) {
    using namespace hornsculpt;
    auto mesh = std::make_shared<SurfaceMesh>();
    mesh->vertices = {{-.25, -.25, 0}, {.25, -.25, 0}, {.25, .25, 0}, {-.25, .25, 0}};
    mesh->triangles = {{{0, 1, 2}, 1}, {{0, 2, 3}, 1}};
    const auto original = *mesh;
    viewport.setMesh(mesh);
    viewport.setHiddenFaces({});
    viewport.setSelectedFaces({});
    viewport.setWireframe(false);
    viewport.setColorMode(MeshColorMode::faceColors);
    viewport.setFaceColors({{1, {35, 70, 220}}});
    auto plane = std::make_shared<PressureFieldPlane>();
    plane->columns = plane->rows = 2;
    plane->minimumLateral = plane->minimumAxial = -1;
    plane->maximumLateral = plane->maximumAxial = 1;
    plane->frame = {{0, 0, -.1}, {0, 1, 0}, {1, 0, 0}};
    plane->frequency = 1000;
    plane->revision = 7;
    plane->current = true;
    for (const double y : {-1., 1.})
        for (const double x : {-1., 1.}) {
            plane->positions.push_back({x, y, -.1});
            plane->pressure.push_back({1, 0});
            plane->valid.push_back(true);
        }
    PressureFieldAppearance appearance{0, 40, 100, 0};
    viewport.setPressureField(plane, appearance);
    viewport.setPressureFieldVisible(true);
    viewport.viewFront({0, 0, 1});
    const auto pixel = [&] {
        viewport.repaint();
        const auto image = viewport.grabFramebuffer();
        return image.pixelColor(image.width() / 2, image.height() / 2);
    };
    const auto behind = pixel();
    if (behind.blue() < 2 * behind.red())
        throw std::runtime_error(
            "Feld hinter CAD wird ohne korrekte Tiefenprüfung darüber gezeichnet.");
    auto front = std::make_shared<PressureFieldPlane>(*plane);
    front->frame.origin.z = .1;
    for (auto &p : front->positions)
        p.z = .1;
    viewport.setPressureField(front, appearance);
    const auto ahead = pixel();
    if (ahead.red() < ahead.blue() * 2 || ahead.red() < ahead.green() * 1.2)
        throw std::runtime_error(
            "Echte 3D-Feldebene vor CAD fehlt oder verwendet falsche Pegelfarben.");
    auto masked = std::make_shared<PressureFieldPlane>(*front);
    masked->valid[0] = false; // Both triangles touch this excluded sample.
    viewport.setPressureField(masked, appearance);
    if (pixel().blue() < pixel().red() * 2)
        throw std::runtime_error("Die Feldebene übermalt eine ungültige Probenmaske.");
    auto wrapped = std::make_shared<PressureFieldPlane>(*front);
    for (std::size_t i = 0; i < wrapped->pressure.size(); ++i)
        wrapped->pressure[i] = std::polar(1., (i % 2 ? -179. : 179.) * std::acos(-1.) / 180.);
    appearance.quantity = 1;
    viewport.setPressureField(wrapped, appearance);
    const auto phase = pixel();
    if (phase.red() < phase.blue() * 2 || phase.red() < phase.green() * 2)
        throw std::runtime_error(
            "±180°-Phasenübergang wird ohne komplexe Interpolation falsch eingefärbt.");
    auto zero = std::make_shared<PressureFieldPlane>(*front);
    std::fill(zero->pressure.begin(), zero->pressure.end(), std::complex<double>{});
    viewport.setPressureField(zero, appearance);
    const auto zeroPhase = pixel();
    const auto colorDistance = [](QColor a, QColor b) {
        return std::max({std::abs(a.red() - b.red()), std::abs(a.green() - b.green()),
                         std::abs(a.blue() - b.blue())});
    };
    if (colorDistance(zeroPhase, behind) > 3)
        throw std::runtime_error("Ein Nullfeld erhält eine erfundene Phase statt einer Lücke.");
    viewport.setPressureField(zero, {0, 40, 100, 0});
    if (colorDistance(pixel(), pressureFieldColor(40, {0, 40, 100, 0})) > 3)
        throw std::runtime_error("Der gültige Nullpegel erhält nicht die minimale Skalenfarbe.");
    for (const double rotation : {0., 180.}) {
        viewport.setPressureField(zero, {2, 40, 94, rotation});
        if (colorDistance(pixel(), pressureFieldColor(0, {2, 40, 94, rotation})) > 3)
            throw std::runtime_error(
                "Der gültige Null-Realteil erhält nicht die korrekte Skalenfarbe.");
    }
    viewport.setPressureField(front, {2, 40, 94, 0});
    const auto positive = pixel();
    viewport.setPressureField(front, {2, 40, 94, 180});
    const auto negative = pixel();
    if (positive.red() < positive.blue() * 2 || negative.blue() < negative.red() * 2)
        throw std::runtime_error(
            "Phasenregler rotiert den Realteil nicht aus denselben Feldphasoren.");
    viewport.focusPressureField();
    drag(viewport, {31, -18});
    viewport.setPressureFieldVisible(false);
    viewport.repaint();
    const auto before = viewport.grabFramebuffer();
    auto nextFrequency = std::make_shared<PressureFieldPlane>(*front);
    nextFrequency->frequency = 2000;
    viewport.setPressureField(nextFrequency, appearance);
    viewport.repaint();
    if (viewport.grabFramebuffer() != before)
        throw std::runtime_error(
            "Frequenzwechsel der Feldebene setzt die manuelle gemeinsame Kamera zurück.");
    if (*mesh != original)
        throw std::runtime_error("Feldanzeige verändert das physische CAD-Netz.");
    viewport.setPressureField({});
    viewport.setPressureFieldVisible(true);
}
} // namespace

int main(int argc, char **argv) {
    QApplication app(argc, argv);
    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addOption(
        {"acoustic-guides-only", "Nur akustische Hilfselemente mit kleinem Modell prüfen."});
    parser.addOption({"interface-only", "Nur Interface-Overlay mit kleinem Modell prüfen."});
    parser.addOption({"output", "JSON-Ausgabe.", "path"});
    parser.addOption({"screenshot", "Viewport-Bild.", "path"});
    parser.addOption({"frames", "Gemessene Frames nach 30 Aufwärmframes.", "count", "240"});
    parser.process(app);
    hornsculpt::MeshViewport viewport;
    viewport.resize(1000, 700);
    viewport.setWindowTitle("HornSculpt · Renderbenchmark");
    std::optional<hornsculpt::DraftGeometry> geometry;
    if (parser.isSet("interface-only") || parser.isSet("acoustic-guides-only"))
        viewport.setMesh(box({2, 2, 2}));
    else {
        geometry.emplace();
        viewport.setMesh(std::make_shared<hornsculpt::SurfaceMesh>(geometry->previewMesh()));
    }
    viewport.show();
    QTimer::singleShot(500, &viewport, [&] {
        try {
            if (parser.isSet("acoustic-guides-only")) {
                checkAcousticGuides(viewport, parser.value("screenshot"));
                const QJsonObject report{{"acousticGuidesPassed", true},
                                         {"sourceTriangles", 12},
                                         {"interfaceTriangles", 1}};
                const auto bytes = QJsonDocument(report).toJson(QJsonDocument::Indented);
                if (!parser.value("output").isEmpty()) {
                    QFile file(parser.value("output"));
                    if (!file.open(QIODevice::WriteOnly) || file.write(bytes) != bytes.size())
                        throw std::runtime_error(
                            "Hilfsebenen-Testbericht konnte nicht gespeichert werden.");
                }
                std::fwrite(bytes.constData(), 1, static_cast<std::size_t>(bytes.size()), stdout);
                app.exit(0);
                return;
            }
            if (parser.isSet("interface-only")) {
                checkInterfaceOverlay(viewport, parser.value("screenshot"));
                {
                    hornsculpt::MeshViewport recreated;
                    recreated.resize(640, 480);
                    recreated.setMesh(box({2, 2, 2}));
                    recreated.show();
                    QApplication::processEvents();
                    checkInterfaceOverlay(recreated);
                }
                const QJsonObject report{{"interfaceOverlayPassed", true},
                                         {"pickingPassed", true},
                                         {"independentContextCleanupPassed", true},
                                         {"sourceTriangles", 12},
                                         {"interfaceTriangles", 1}};
                const auto bytes = QJsonDocument(report).toJson(QJsonDocument::Indented);
                if (!parser.value("output").isEmpty()) {
                    QFile file(parser.value("output"));
                    if (!file.open(QIODevice::WriteOnly) || file.write(bytes) != bytes.size())
                        throw std::runtime_error(
                            "Interface-Testbericht konnte nicht gespeichert werden.");
                }
                std::fwrite(bytes.constData(), 1, static_cast<std::size_t>(bytes.size()), stdout);
                app.exit(0);
                return;
            }
            bool ok = false;
            const auto frames = parser.value("frames").toULongLong(&ok);
            if (!ok)
                throw std::invalid_argument("Ungültige Framezahl.");
            QJsonArray measurements;
            auto preview = measure(viewport, static_cast<std::size_t>(frames));
            preview["case"] = "default-preview";
            measurements.append(preview);
            viewport.setMesh(
                std::make_shared<hornsculpt::SurfaceMesh>(geometry->previewMesh(.0001)));
            viewport.repaint();
            auto fine = measure(viewport, static_cast<std::size_t>(frames));
            fine["case"] = "fine-surface-tessellation";
            measurements.append(fine);
            auto triangle = std::make_shared<hornsculpt::SurfaceMesh>();
            triangle->vertices = {{-1, -1, 0},  {1, -1, 0},  {0, 1, 0},
                                  {-1, -1, .2}, {1, -1, .2}, {0, 1, .2}};
            triangle->triangles = {{{0, 1, 2}, 55, hornsculpt::BoundaryCondition::wall},
                                   {{3, 5, 4}, 56, hornsculpt::BoundaryCondition::wall}};
            viewport.setMesh(triangle);
            viewport.fitView();
            viewport.repaint();
            if (viewport.pickFace(QPointF(viewport.width() * .5, viewport.height() * .5)) != 56 ||
                viewport.pickFace(QPointF(-5, -5)))
                throw std::runtime_error("Face-ID-Picking liefert falschen Treffer.");
            checkFaceColors(viewport, *triangle);
            const QPointF center(viewport.width() * .5, viewport.height() * .5);
            viewport.setHiddenFaces({56});
            viewport.repaint();
            if (viewport.pickFace(center) != 55 || triangle->triangles.size() != 2)
                throw std::runtime_error(
                    "Ausgeblendete Apertur verdeckt die Front oder verändert das Netz.");
            viewport.setHiddenFaces({55, 56});
            viewport.repaint();
            if (viewport.pickFace(center))
                throw std::runtime_error("Verborgene Flächen bleiben auswählbar.");
            viewport.setHiddenFaces({});
            viewport.repaint();
            click(viewport, center);
            if (viewport.selectedFaces() != std::vector<std::uint32_t>{56})
                throw std::runtime_error("Klick ersetzt die Auswahl nicht.");
            click(viewport, center, Qt::ControlModifier);
            if (!viewport.selectedFaces().empty())
                throw std::runtime_error("Strg-Klick entfernt die Fläche nicht.");
            viewport.setSelectedFaces({55, 56, 55});
            click(viewport, QPointF(1, 1), Qt::ControlModifier);
            if (viewport.selectedFaces() != std::vector<std::uint32_t>{55, 56})
                throw std::runtime_error("Strg-Leerklick verändert die Auswahl.");
            click(viewport, QPointF(1, 1));
            if (!viewport.selectedFaces().empty())
                throw std::runtime_error("Leerklick leert die Auswahl nicht.");
            viewport.setSelectedFaces({55});
            QMouseEvent press(QEvent::MouseButtonPress, center,
                              viewport.mapToGlobal(center.toPoint()), Qt::LeftButton,
                              Qt::LeftButton, Qt::ControlModifier);
            QMouseEvent move(QEvent::MouseMove, center + QPointF(20, 0),
                             viewport.mapToGlobal((center + QPointF(20, 0)).toPoint()),
                             Qt::NoButton, Qt::LeftButton, Qt::ControlModifier);
            QMouseEvent release(QEvent::MouseButtonRelease, center,
                                viewport.mapToGlobal(center.toPoint()), Qt::LeftButton,
                                Qt::NoButton, Qt::ControlModifier);
            QApplication::sendEvent(&viewport, &press);
            QApplication::sendEvent(&viewport, &move);
            QApplication::sendEvent(&viewport, &release);
            if (viewport.selectedFaces() != std::vector<std::uint32_t>{55})
                throw std::runtime_error("Kameraziehen verändert die Auswahl.");
            viewport.setMesh(
                std::make_shared<hornsculpt::SurfaceMesh>(geometry->previewMesh(.0001)));
            viewport.repaint();
            if (!viewport.selectedFaces().empty())
                throw std::runtime_error("Ein Modellwechsel übernimmt fremde Flächen-IDs.");
            viewport.setSelectedFaces({1, 2, 3});
            if (viewport.selectedFaces() != std::vector<std::uint32_t>{1, 2, 3})
                throw std::runtime_error("Mehrfachauswahl wurde nicht übernommen.");
            (void)viewport.benchmarkFrames(2, 0);
            const QImage image = viewport.grabFramebuffer();
            std::size_t highlightedPixels{};
            for (int y = 0; y < image.height(); y += 4)
                for (int x = 0; x < image.width(); x += 4) {
                    const auto color = image.pixelColor(x, y);
                    if (color.red() > color.green() * 1.2 && color.green() > color.blue() * 1.4)
                        ++highlightedPixels;
                }
            if (highlightedPixels < 100)
                throw std::runtime_error("GPU-Flächenhervorhebung fehlt im Framebuffer.");
            auto orientedTriangle = std::make_shared<hornsculpt::SurfaceMesh>();
            orientedTriangle->vertices = {{-1, -1, 0}, {1, -1, 0}, {0, 1, 0}};
            orientedTriangle->triangles = {{{0, 1, 2}, 7, hornsculpt::BoundaryCondition::wall}};
            viewport.setMesh(orientedTriangle);
            viewport.setSelectedFaces({});
            viewport.fitView();
            viewport.setColorMode(hornsculpt::MeshColorMode::orientation);
            viewport.setWireframe(false);
            (void)viewport.benchmarkFrames(2, 0);
            const auto frontImage = viewport.grabFramebuffer();
            const auto front =
                frontImage.pixelColor(frontImage.width() / 2, frontImage.height() / 2);
            std::swap(orientedTriangle->triangles[0].vertices[1],
                      orientedTriangle->triangles[0].vertices[2]);
            viewport.setMesh(std::make_shared<hornsculpt::SurfaceMesh>(*orientedTriangle));
            (void)viewport.benchmarkFrames(2, 0);
            const auto backImage = viewport.grabFramebuffer();
            const auto back = backImage.pixelColor(backImage.width() / 2, backImage.height() / 2);
            if (front.green() < front.red() * 1.5 || back.red() < back.green() * 1.5)
                throw std::runtime_error("Normalenfarben folgen nicht der Dreiecksorientierung.");
            checkPerspectiveFit(viewport);
            checkFrontView(viewport);
            checkInterfaceOverlay(viewport);
            checkPressureField(viewport);
            if (!parser.value("screenshot").isEmpty() && !image.save(parser.value("screenshot")))
                throw std::runtime_error("Screenshot konnte nicht gespeichert werden.");
            const QJsonObject report{
                {"measurement", "Unpaced scene and overlay rendering to GPU completion; excludes "
                                "widget composition, presentation and vertical sync."},
                {"wireframe", true},
                {"pickingPassed", true},
                {"highlightPassed", true},
                {"orientationPassed", true},
                {"faceColorsPassed", true},
                {"perspectiveFitPassed", true},
                {"frontViewPassed", true},
                {"pressureField3dPassed", true},
                {"measurements", measurements}};
            const auto bytes = QJsonDocument(report).toJson(QJsonDocument::Indented);
            if (!parser.value("output").isEmpty()) {
                QFile file(parser.value("output"));
                if (!file.open(QIODevice::WriteOnly) || file.write(bytes) != bytes.size())
                    throw std::runtime_error("Benchmarkdatei konnte nicht gespeichert werden.");
            }
            std::fwrite(bytes.constData(), 1, static_cast<std::size_t>(bytes.size()), stdout);
            app.exit(0);
        } catch (const std::exception &error) {
            std::fprintf(stderr, "Renderbenchmark: %s\n", error.what());
            app.exit(1);
        }
    });
    return app.exec();
}
