#include "cad_fixture.hpp"
#include "fixture_geometry.hpp"
#include "project_file.hpp"
#include "result_export.hpp"
#include "simulation_request.hpp"
#include <hornsculpt/axisymmetric_mesh.hpp>
#include <hornsculpt/mesh_orientation.hpp>
#include <hornsculpt/surface_measures.hpp>

#include <QCoreApplication>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QStringList>
#include <QTemporaryDir>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

using namespace hornsculpt;

namespace {
constexpr std::uint64_t largeRevision = (std::uint64_t{1} << 53) + 17;

void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
QJsonObject preWidebandArchive(QJsonObject root, unsigned version) {
    root.insert("version", int(version));
    if (version < 15)
        root.remove("comparisonSnapshots");
    if (version < 13)
        root.remove("slot");
    auto simulation = root["simulation"].toObject();
    const auto historicalSettings = [](QJsonObject settings) {
        auto fmm = settings["fastMultipole"].toObject();
        if (!fmm.isEmpty()) {
            fmm.remove("preconditioner_block_width");
            fmm.insert("leafPanelCapacity", fmm.take("leaf_support_capacity"));
            fmm.insert("interpolationOrder", 6);
            settings.insert("fastMultipole", fmm);
        }
        return settings;
    };
    simulation.insert("settings", historicalSettings(simulation["settings"].toObject()));
    auto request = simulation["request"].toObject();
    if (!request.isEmpty()) {
        request.insert("numerics", historicalSettings(request["numerics"].toObject()));
        simulation.insert("request", request);
    }
    root.insert("simulation", simulation);
    auto result = simulation["result"].toObject();
    if (result.isEmpty())
        return root;
    result.remove("liveDeviceBytes");
    auto frequencies = result["frequencies"].toArray();
    for (qsizetype i = 0; i < frequencies.size(); ++i) {
        auto f = frequencies[i].toObject();
        auto stats = f["fastMultipole"].toObject();
        if (!stats.isEmpty()) {
            const auto history = stats["historicalChebyshevStatistics"].toObject();
            for (const auto *name :
                 {"pressurePressureNearPairs", "pressurePressureFarPairs", "pressureFluxNearPairs",
                  "pressureFluxFarPairs", "historicalChebyshevStatistics", "graphExecution"})
                stats.remove(name);
            stats.insert("directPanelPairCount",
                         history.isEmpty() ? QJsonValue("0") : history["directPanelPairCount"]);
            stats.insert("representedFarPanelPairs",
                         history.isEmpty() ? QJsonValue("0") : history["representedFarPanelPairs"]);
            f.insert("fastMultipole", stats);
        }
        auto excitations = f["excitations"].toArray();
        for (qsizetype j = 0; j < excitations.size(); ++j) {
            auto e = excitations[j].toObject();
            const auto history = e["historicalChebyshevEvidence"].toObject();
            for (const auto *name :
                 {"fastMultipoleAccuracySeed", "fastMultipoleComparisonAccuracySeed",
                  "fastMultipoleVerificationAccuracySeed", "fastMultipoleAccuracyRefinements",
                  "fastMultipoleComparisonLevels", "fastMultipoleVerificationLevels",
                  "historicalChebyshevEvidence", "quadratureEvidence", "fastMultipoleKrylovTarget"})
                e.remove(name);
            e.insert("fastMultipoleComparisonOrder",
                     history.isEmpty() ? QJsonValue(0) : history["comparisonOrder"]);
            e.insert("fastMultipoleAuditOrder",
                     history.isEmpty() ? QJsonValue(0) : history["auditOrder"]);
            excitations[j] = e;
        }
        f.insert("excitations", excitations);
        frequencies[i] = f;
    }
    result.insert("frequencies", frequencies);
    simulation.insert("result", result);
    root.insert("simulation", simulation);
    return root;
}
void close(double actual, double expected, const char *message) {
    require(std::isfinite(actual) && std::isfinite(expected) &&
                std::abs(actual - expected) <= 1e-14 * std::max(1.0, std::abs(expected)),
            message);
}
void sameComplex(std::complex<double> actual, std::complex<double> expected, const char *message) {
    close(actual.real(), expected.real(), message);
    close(actual.imag(), expected.imag(), message);
}
void samePoint(Vec3 actual, Vec3 expected) {
    close(actual.x, expected.x, "Point X changed");
    close(actual.y, expected.y, "Point Y changed");
    close(actual.z, expected.z, "Point Z changed");
}
void sameMesh(const SurfaceMesh &actual, const SurfaceMesh &expected) {
    require(actual.vertices.size() == expected.vertices.size() &&
                actual.triangles.size() == expected.triangles.size() &&
                actual.fixedVertices == expected.fixedVertices,
            "Mesh dimensions changed");
    for (std::size_t index = 0; index < actual.vertices.size(); ++index)
        samePoint(actual.vertices[index], expected.vertices[index]);
    for (std::size_t index = 0; index < actual.triangles.size(); ++index) {
        const auto &a = actual.triangles[index], &b = expected.triangles[index];
        require(a.vertices == b.vertices && a.faceId == b.faceId && a.boundary == b.boundary,
                "Mesh connectivity, physical tags or boundary conditions changed");
    }
}

Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double length(Vec3 point) {
    return std::hypot(point.x, point.y, point.z);
}

SurfaceMesh testMesh(double scale = 1) {
    auto mesh = hornsculpt::testing::sphere(0);
    mesh.fixedVertices = {0, 2};
    for (auto &vertex : mesh.vertices) {
        vertex.x *= scale;
        vertex.y *= scale;
        vertex.z *= scale;
    }
    for (std::size_t index = 0; index < mesh.triangles.size(); ++index) {
        mesh.triangles[index].faceId = static_cast<std::uint32_t>(index + 1);
        mesh.triangles[index].boundary =
            index < 2 ? BoundaryCondition::driven : BoundaryCondition::wall;
    }
    return mesh;
}

std::shared_ptr<CadSnapshot> testCad(const SurfaceMesh &mesh) {
    auto cad = std::make_shared<CadSnapshot>();
    cad->documentId = QStringLiteral("test-cad-document");
    cad->documentName = QStringLiteral("Prüfkörper");
    cad->sessionId = QStringLiteral("offline-test-session");
    cad->revision = largeRevision + 401;
    cad->mesh = mesh;
    cad->acousticOrigin = {.01, -.02, .03};
    cad->acousticAxis = {0, 1, 0};
    cad->objects.push_back({QStringLiteral("object-a"),
                            QStringLiteral("Solid"),
                            QStringLiteral("Synthetic test solid"),
                            QStringLiteral("octahedron"),
                            {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}});
    for (std::size_t index = 0; index < mesh.triangles.size(); ++index) {
        const auto &triangle = mesh.triangles[index];
        const auto a = mesh.vertices[triangle.vertices[0]], b = mesh.vertices[triangle.vertices[1]],
                   c = mesh.vertices[triangle.vertices[2]];
        const auto normal = cross(subtract(b, a), subtract(c, a));
        const double norm = length(normal);
        CadFace face;
        face.id = triangle.faceId;
        face.objectId = QStringLiteral("object-a");
        face.sourceKey = QStringLiteral("Face%1").arg(index + 1);
        face.label = QStringLiteral("Test surface %1").arg(index + 1);
        face.signature = {QStringLiteral("Plane"),
                          QStringLiteral("geometry-%1").arg(index),
                          QStringLiteral("topology-%1").arg(index),
                          norm / 2,
                          {(a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3, (a.z + b.z + c.z) / 3},
                          {normal.x / norm, normal.y / norm, normal.z / norm},
                          3,
                          1};
        face.boundary = triangle.boundary;
        face.normalVelocity =
            index == 0 ? std::complex<double>{.3, .6} : std::complex<double>{-.2, .7};
        face.state = CadFaceState::confirmed;
        cad->faces.push_back(face);
    }
    return cad;
}

SimulationArchive syntheticTestArchive(std::shared_ptr<const SurfaceMesh> mesh) {
    // Structurally valid persistence fixture only. These values were not computed
    // by an acoustic solver and must never be used as numerical reference data.
    SimulationArchive archive;
    archive.settings = {315, 1000, 2, {1.18, 346.2}, .7, 27};
    archive.mesh = std::move(mesh);
    archive.revision = largeRevision;
    archive.current = true;
    BemRequest request;
    request.medium = archive.settings.medium;
    request.frequencies = {315, 1000};
    request.excitations = {{"synthetic test excitation A", {{1, {.3, .6}}, {2, {-.2, .7}}}},
                           {"synthetic test excitation B", {{1, {-.4, .2}}, {2, {.9, -.1}}}}};
    request.observers = {{.3, 0, .5}, {-.2, .4, .5}};
    request.numerics = {4, 6, 1, 1e-4, 6, 5e-5, 3};
    archive.request = request;
    auto result = std::make_shared<BemResult>();
    result->device = "SYNTHETIC PERSISTENCE TEST FIXTURE - no GPU computation";
    result->triangleCount = archive.mesh->triangles.size();
    result->pressureDofCount = archive.mesh->vertices.size();
    result->peakDeviceBytes = 1234567;
    for (const auto frequency : request.frequencies) {
        BemFrequencyResult frame;
        frame.frequency = frequency;
        frame.assemblyMilliseconds = 12.3;
        frame.solveMilliseconds = 1.7;
        frame.factorizationMilliseconds = 1;
        frame.triangularSolveMilliseconds = .4;
        frame.fieldMilliseconds = .5;
        for (std::size_t excitation = 0; excitation < request.excitations.size(); ++excitation) {
            BemExcitationResult output;
            BemBoundaryData boundary;
            output.name = request.excitations[excitation].name;
            output.refinementSteps = 1;
            output.boundaryResidual = 2e-6;
            for (std::size_t panel = 0; panel < archive.mesh->triangles.size(); ++panel)
                boundary.surfacePressure.emplace_back(
                    static_cast<double>(panel + 1) * .125 + static_cast<double>(excitation), -.7);
            for (std::size_t observer = 0; observer < request.observers.size(); ++observer)
                output.observerPressure.emplace_back(static_cast<double>(observer) + .75,
                                                     -.2 - static_cast<double>(excitation));
            for (const auto &face : request.excitations[excitation].faces) {
                const auto &triangle = archive.mesh->triangles[face.faceId - 1];
                const auto a = archive.mesh->vertices[triangle.vertices[0]],
                           b = archive.mesh->vertices[triangle.vertices[1]],
                           c = archive.mesh->vertices[triangle.vertices[2]];
                const double area = length(cross(subtract(b, a), subtract(c, a))) / 2;
                const auto mean = boundary.surfacePressure[face.faceId - 1];
                const auto volume = face.normalVelocity * area;
                output.impedances.push_back({face.faceId, area, mean, volume, mean / volume});
            }
            output.setBoundary(std::move(boundary));
            frame.excitations.push_back(std::move(output));
        }
        result->frequencies.push_back(std::move(frame));
    }
    archive.result = result;
    return archive;
}

Project fullTestProject(ProjectSource source) {
    Project project;
    project.source = source;
    project.revision = project.geometryRevision = project.meshRevision = largeRevision;
    project.meshing.targetEdgeLength = .0065;
    project.meshing.adaptive = {true, .0002, .00003, .45};
    project.meshing.quality = {false, 27.5, 3.25};
    project.draftFaceTargetEdgeLengths = {{5, .002}};
    project.cadFaceTargetEdgeLengths = {{2, .003}};
    auto &p = project.parameters;
    p.horizontalAngle = .72;
    p.verticalAngle = .46;
    p.entry = {.026, .018, 3.5};
    p.length = .19;
    p.profile.throatAngle = .04;
    p.profile.k = 1.1;
    p.profile.terminationScale = .1;
    p.profile.terminationQ = .96;
    p.profile.terminationExponent = 4.2;
    p.verticalProfile = {.03, .87, .08, .94, 5.1};
    p.crossSectionExponent = 3.1;
    p.wallThickness = .006;
    p.guideEnabled = true;
    p.guide.type = CurveType::superformula;
    p.guide.width = .28;
    p.guide.aspectRatio = .75;
    p.guide.rotation = .17;
    p.guide.position = .9;
    p.guide.positionMode = GuidePositionMode::relative;
    p.guide.superellipseExponent = 3.5;
    p.guide.superformula = {1.1, .9, 4, 4, 2.2, 2.4, 2.1};
    p.transition = {TransitionMethod::lut, .11, .83, {{0, 0}, {.2, .06}, {.7, .9}, {1, 1}}};
    p.mouthMorph.enabled = true;
    p.mouthMorph.startFraction = .92;
    p.mouthMorph.endFraction = 1;
    p.mouthMorph.target.width = .31;
    p.mouthMorph.target.aspectRatio = .72;
    p.mouthMorph.target.rotation = .12;
    p.mouthMorph.target.superellipseExponent = 3.7;
    (void)DraftGeometry(p);
    const auto mesh = std::make_shared<SurfaceMesh>(testMesh());
    project.preview = mesh;
    project.mesh = std::make_shared<MeshingResult>(MeshingResult{*mesh, inspectMesh(*mesh)});
    project.cad = testCad(*mesh);
    project.simulation = syntheticTestArchive(mesh);
    project.simulation.settings.beamwidthTargetEnabled = true;
    project.simulation.settings.targetHorizontalDegrees = 100;
    project.simulation.settings.targetVerticalDegrees = 65;
    project.simulation.settings.targetToleranceDegrees = 12;
    project.simulation.measurements.push_back(
        {QStringLiteral("Measured fixture"), false, 82, {{100, -3}, {1000, 0}}});
    return project;
}

QByteArray bytes(const QString &path) {
    QFile file(path);
    require(file.open(QIODevice::ReadOnly), "Test project file cannot be read");
    return file.readAll();
}
void saveBytes(const QString &path, const QByteArray &value) {
    QFile file(path);
    require(file.open(QIODevice::WriteOnly | QIODevice::Truncate),
            "Test fixture cannot be written");
    require(file.write(value) == value.size(), "Test fixture write is incomplete");
}

void expectRejected(const std::function<void()> &operation, const char *message) {
    try {
        operation();
    } catch (const std::exception &) {
        return;
    }
    throw std::runtime_error(message);
}

void sameSettings(const SimulationSettings &a, const SimulationSettings &b) {
    require(a.geometryMode == b.geometryMode && a.axisymmetricSampling == b.axisymmetricSampling,
            "Geometry mode or meridian sampling changed");
    require(a.ctaExcitationIndex == b.ctaExcitationIndex, "CTA selection changed");
    require(a.fieldDisplay == b.fieldDisplay && a.phaseDisplay == b.phaseDisplay,
            "Field or phase display settings changed");
    require(a.angularResponses == b.angularResponses, "Angular response display settings changed");
    require(a.infiniteBaffle == b.infiniteBaffle, "Infinite baffle settings changed");
    require(a.acceleration == b.acceleration && a.fastMultipole == b.fastMultipole,
            "Fast Multipole settings changed");
    require(a.nearMaximumDepth == b.nearMaximumDepth, "Near quadrature depth changed");
    require(a.singularOrder == b.singularOrder, "Singular quadrature order changed");
    close(a.minimumFrequency, b.minimumFrequency, "Minimum frequency changed");
    close(a.maximumFrequency, b.maximumFrequency, "Maximum frequency changed");
    require(a.frequencyCount == b.frequencyCount, "Frequency count changed");
    close(a.medium.density, b.medium.density, "Medium density changed");
    close(a.medium.soundSpeed, b.medium.soundSpeed, "Sound speed changed");
    close(a.velocityScale, b.velocityScale, "Velocity scale changed");
    close(a.phaseDegrees, b.phaseDegrees, "Drive phase changed");
    require(a.normalizeToOnAxis == b.normalizeToOnAxis, "Result normalization setting changed");
    require(a.beamwidthTargetEnabled == b.beamwidthTargetEnabled, "Beamwidth target mode changed");
    close(a.targetHorizontalDegrees, b.targetHorizontalDegrees, "Horizontal target changed");
    close(a.targetVerticalDegrees, b.targetVerticalDegrees, "Vertical target changed");
    close(a.targetToleranceDegrees, b.targetToleranceDegrees, "Target tolerance changed");
    close(a.normalizationAngleDegrees, b.normalizationAngleDegrees, "Normalization angle changed");
    close(a.directivityMinimumDb, b.directivityMinimumDb, "Directivity minimum changed");
    close(a.directivityMaximumDb, b.directivityMaximumDb, "Directivity maximum changed");
    close(a.directivityStepDb, b.directivityStepDb, "Directivity scale step changed");
    require(a.couplingEnabled == b.couplingEnabled && a.couplingNetwork == b.couplingNetwork,
            "The LEM editor graph or coupling mode changed");
    require(a.farField == b.farField, "The observation field mode changed");
    close(a.observationDistance, b.observationDistance, "The observation distance changed");
    require(a.observationFrame.has_value() == b.observationFrame.has_value(),
            "Acoustic frame override mode changed");
    if (a.observationFrame) {
        samePoint(a.observationFrame->origin, b.observationFrame->origin);
        samePoint(a.observationFrame->forward, b.observationFrame->forward);
        samePoint(a.observationFrame->horizontal, b.observationFrame->horizontal);
    }
}

void sameDraft(const DraftParameters &a, const DraftParameters &b) {
    const auto &x = a.profile, &y = b.profile;
    require(a.asymmetricEnabled == b.asymmetricEnabled, "Asymmetric mode changed");
    close(a.leftAngle, b.leftAngle, "Left angle changed");
    close(a.bottomAngle, b.bottomAngle, "Bottom angle changed");
    for (const auto pair : {std::pair{&a.leftProfile, &b.leftProfile},
                            std::pair{&a.bottomProfile, &b.bottomProfile}}) {
        close(pair.first->throatAngle, pair.second->throatAngle, "Opposite throat angle changed");
        close(pair.first->k, pair.second->k, "Opposite OSSE k changed");
        close(pair.first->terminationScale, pair.second->terminationScale,
              "Opposite termination scale changed");
        close(pair.first->terminationQ, pair.second->terminationQ,
              "Opposite termination Q changed");
        close(pair.first->terminationExponent, pair.second->terminationExponent,
              "Opposite termination exponent changed");
    }
    close(a.entry.width, b.entry.width, "Entry width changed");
    close(a.entry.height, b.entry.height, "Entry height changed");
    close(a.entry.exponent, b.entry.exponent, "Entry exponent changed");
    require(a.slot.enabled == b.slot.enabled, "Diffraction slot enabled state changed");
    require(a.slot.contour == b.slot.contour && a.slot.exit == b.slot.exit,
            "Diffraction slot rounding modes changed");
    close(a.slot.cornerRadius, b.slot.cornerRadius, "Slot corner radius changed");
    close(a.slot.exitRadius, b.slot.exitRadius, "Slot exit radius changed");
    close(a.slot.aperture.width, b.slot.aperture.width, "Slot width changed");
    close(a.slot.aperture.height, b.slot.aperture.height, "Slot height changed");
    close(a.slot.aperture.exponent, b.slot.aperture.exponent, "Slot exponent changed");
    close(a.slot.adapterLength, b.slot.adapterLength, "Slot adapter length changed");
    close(a.slot.straightLength, b.slot.straightLength, "Slot straight length changed");
    close(a.slot.transitionLength, b.slot.transitionLength, "Slot transition length changed");
    close(a.length, b.length, "Horn length changed");
    require(a.mouthMorph.enabled == b.mouthMorph.enabled, "Mouth morph state changed");
    close(a.mouthMorph.startFraction, b.mouthMorph.startFraction, "Mouth morph start changed");
    close(a.mouthMorph.endFraction, b.mouthMorph.endFraction, "Mouth morph end changed");
    const auto &mt = a.mouthMorph.target, &nt = b.mouthMorph.target;
    require(mt.type == nt.type, "Mouth target type changed");
    close(mt.width, nt.width, "Mouth target width changed");
    close(mt.aspectRatio, nt.aspectRatio, "Mouth target aspect changed");
    close(mt.rotation, nt.rotation, "Mouth target rotation changed");
    close(mt.superellipseExponent, nt.superellipseExponent, "Mouth target shape changed");
    close(mt.superformula.a, nt.superformula.a, "Mouth superformula a changed");
    close(mt.superformula.b, nt.superformula.b, "Mouth superformula b changed");
    close(mt.superformula.m1, nt.superformula.m1, "Mouth superformula m1 changed");
    close(mt.superformula.m2, nt.superformula.m2, "Mouth superformula m2 changed");
    close(mt.superformula.n1, nt.superformula.n1, "Mouth superformula n1 changed");
    close(mt.superformula.n2, nt.superformula.n2, "Mouth superformula n2 changed");
    close(mt.superformula.n3, nt.superformula.n3, "Mouth superformula n3 changed");
    close(x.throatAngle, y.throatAngle, "Throat angle changed");
    close(x.k, y.k, "OSSE k changed");
    close(x.terminationScale, y.terminationScale, "Termination scale changed");
    close(x.terminationQ, y.terminationQ, "Termination Q changed");
    close(x.terminationExponent, y.terminationExponent, "Termination exponent changed");
    close(a.verticalProfile.throatAngle, b.verticalProfile.throatAngle,
          "Vertical throat angle changed");
    close(a.verticalProfile.k, b.verticalProfile.k, "Vertical OSSE k changed");
    close(a.verticalProfile.terminationScale, b.verticalProfile.terminationScale,
          "Vertical termination scale changed");
    close(a.verticalProfile.terminationQ, b.verticalProfile.terminationQ,
          "Vertical termination Q changed");
    close(a.verticalProfile.terminationExponent, b.verticalProfile.terminationExponent,
          "Vertical termination exponent changed");
    close(a.horizontalAngle, b.horizontalAngle, "Horizontal angle changed");
    close(a.verticalAngle, b.verticalAngle, "Vertical angle changed");
    close(a.crossSectionExponent, b.crossSectionExponent, "Cross section exponent changed");
    close(a.wallThickness, b.wallThickness, "Wall thickness changed");
    require(a.guideEnabled == b.guideEnabled, "Guide enabled state changed");
    const auto &g = a.guide, &h = b.guide;
    require(g.type == h.type && g.positionMode == h.positionMode, "Guide strategy changed");
    close(g.width, h.width, "Guide width changed");
    close(g.aspectRatio, h.aspectRatio, "Guide aspect ratio changed");
    close(g.rotation, h.rotation, "Guide rotation changed");
    close(g.position, h.position, "Guide position changed");
    close(g.superellipseExponent, h.superellipseExponent, "Superellipse exponent changed");
    const auto &sf = g.superformula, &tf = h.superformula;
    close(sf.a, tf.a, "Superformula a changed");
    close(sf.b, tf.b, "Superformula b changed");
    close(sf.m1, tf.m1, "Superformula m1 changed");
    close(sf.m2, tf.m2, "Superformula m2 changed");
    close(sf.n1, tf.n1, "Superformula n1 changed");
    close(sf.n2, tf.n2, "Superformula n2 changed");
    close(sf.n3, tf.n3, "Superformula n3 changed");
    require(a.transition.method == b.transition.method &&
                a.transition.knots.size() == b.transition.knots.size(),
            "Transition strategy changed");
    close(a.transition.startFraction, b.transition.startFraction, "Transition start changed");
    close(a.transition.endFraction, b.transition.endFraction, "Transition end changed");
    for (std::size_t index = 0; index < a.transition.knots.size(); ++index) {
        close(a.transition.knots[index].fraction, b.transition.knots[index].fraction,
              "LUT position changed");
        close(a.transition.knots[index].blend, b.transition.knots[index].blend,
              "LUT blend changed");
    }
}

void sameCad(const CadSnapshot &a, const CadSnapshot &b) {
    require(a.documentId == b.documentId && a.documentName == b.documentName &&
                a.sessionId == b.sessionId && a.revision == b.revision &&
                a.unresolvedPreviousFaces == b.unresolvedPreviousFaces,
            "CAD identity or revision changed");
    samePoint(a.acousticOrigin, b.acousticOrigin);
    samePoint(a.acousticAxis, b.acousticAxis);
    sameMesh(a.mesh, b.mesh);
    require(a.brepObjects.size() == b.brepObjects.size() && bool(a.geometry) == bool(b.geometry),
            "Exact CAD source or native geometry cache disappeared");
    close(a.brepMetresPerUnit, b.brepMetresPerUnit, "BRep source units changed");
    for (std::size_t i = 0; i < a.brepObjects.size(); ++i)
        require(a.brepObjects[i].objectId == b.brepObjects[i].objectId &&
                    a.brepObjects[i].data == b.brepObjects[i].data &&
                    a.brepObjects[i].faceIds == b.brepObjects[i].faceIds,
                "Exact CAD bytes or physical face binding changed");
    require(a.faces.size() == b.faces.size() && a.objects.size() == b.objects.size(),
            "CAD metadata dimensions changed");
    for (std::size_t index = 0; index < a.objects.size(); ++index) {
        const auto &x = a.objects[index], &y = b.objects[index];
        require(x.id == y.id && x.name == y.name && x.label == y.label &&
                    x.topologyFingerprint == y.topologyFingerprint,
                "CAD object identity changed");
        for (unsigned coordinate = 0; coordinate < 16; ++coordinate)
            close(x.placement[coordinate], y.placement[coordinate], "CAD placement changed");
    }
    for (std::size_t index = 0; index < a.faces.size(); ++index) {
        const auto &x = a.faces[index], &y = b.faces[index];
        require(x.id == y.id && x.objectId == y.objectId && x.sourceKey == y.sourceKey &&
                    x.label == y.label && x.boundary == y.boundary && x.state == y.state,
                "CAD face identity, assignment or review state changed");
        sameComplex(x.normalVelocity, y.normalVelocity, "Complex CAD face velocity changed");
        const auto &s = x.signature, &t = y.signature;
        require(s.surfaceType == t.surfaceType && s.geometryKey == t.geometryKey &&
                    s.topologyKey == t.topologyKey && s.edgeCount == t.edgeCount &&
                    s.wireCount == t.wireCount,
                "CAD face signature changed");
        close(s.area, t.area, "CAD face area changed");
        samePoint(s.center, t.center);
        samePoint(s.normal, t.normal);
    }
}

void sameRequest(const BemRequest &a, const BemRequest &b) {
    require(a.interactionPlanFrequencies == b.interactionPlanFrequencies,
            "The immutable interaction planning raster changed");
    require(a.axisymmetric == b.axisymmetric, "The immutable m=0 reduction changed");
    require(a.infiniteBaffle == b.infiniteBaffle, "The immutable baffle domain changed");
    require(a.coupling == b.coupling, "The immutable LEM request graph changed");
    require(a.farField == b.farField,
            "Finite/far-field mode, origin or reference distance changed");
    close(a.medium.density, b.medium.density, "Request density changed");
    close(a.medium.soundSpeed, b.medium.soundSpeed, "Request sound speed changed");
    require(a.frequencies == b.frequencies && a.excitations.size() == b.excitations.size() &&
                a.observers.size() == b.observers.size(),
            "BEM request dimensions changed");
    for (std::size_t index = 0; index < a.observers.size(); ++index)
        samePoint(a.observers[index], b.observers[index]);
    const auto &x = a.numerics, &y = b.numerics;
    require(x.acceleration == y.acceleration && x.fastMultipole == y.fastMultipole,
            "Fast Multipole request changed");
    require(x.quadratureOrder == y.quadratureOrder && x.singularOrder == y.singularOrder &&
                x.nearMaximumDepth == y.nearMaximumDepth &&
                x.maximumRefinementSteps == y.maximumRefinementSteps,
            "BEM integration strategy changed");
    close(x.nearDistanceRatio, y.nearDistanceRatio, "Near distance ratio changed");
    close(x.nearRelativeTolerance, y.nearRelativeTolerance, "Near tolerance changed");
    close(x.relativeResidualTolerance, y.relativeResidualTolerance, "Residual tolerance changed");
    for (std::size_t index = 0; index < a.excitations.size(); ++index) {
        const auto &u = a.excitations[index], &v = b.excitations[index];
        require(u.name == v.name && u.faces.size() == v.faces.size(),
                "Excitation identity changed");
        for (std::size_t face = 0; face < u.faces.size(); ++face) {
            require(u.faces[face].faceId == v.faces[face].faceId, "Excitation face changed");
            sameComplex(u.faces[face].normalVelocity, v.faces[face].normalVelocity,
                        "Complex excitation velocity changed");
        }
    }
}

void sameLemSolution(const lem::Solution &a, const lem::Solution &b) {
    require(a.frequency == b.frequency && a.relativeResidual == b.relativeResidual &&
                a.nodes.size() == b.nodes.size() && a.elements.size() == b.elements.size() &&
                a.drivers.size() == b.drivers.size() && a.ports.size() == b.ports.size(),
            "The LEM solution metadata or dimensions changed");
    for (std::size_t i = 0; i < a.nodes.size(); ++i)
        require(a.nodes[i].id == b.nodes[i].id && a.nodes[i].effort == b.nodes[i].effort,
                "A LEM node potential changed");
    for (std::size_t i = 0; i < a.elements.size(); ++i) {
        const auto &x = a.elements[i], &y = b.elements[i];
        require(x.id == y.id && x.effort == y.effort && x.flow == y.flow &&
                    x.activePower == y.activePower,
                "A LEM element result changed");
    }
    for (std::size_t i = 0; i < a.drivers.size(); ++i) {
        const auto &x = a.drivers[i], &y = b.drivers[i];
        require(x.id == y.id && x.current == y.current && x.velocity == y.velocity &&
                    x.displacement == y.displacement && x.electricalPower == y.electricalPower &&
                    x.copperPower == y.copperPower && x.mechanicalPower == y.mechanicalPower &&
                    x.acousticPower == y.acousticPower,
                "A computed driver quantity changed");
    }
    for (std::size_t i = 0; i < a.ports.size(); ++i) {
        const auto &x = a.ports[i], &y = b.ports[i];
        require(x.id == y.id && x.pressure == y.pressure && x.volumeVelocity == y.volumeVelocity &&
                    x.activePower == y.activePower,
                "A physical BEM port quantity changed");
    }
}

void sameResult(const BemResult &a, const BemResult &b) {
    require(a.geometryMode == b.geometryMode, "The computed geometry mode changed");
    require(a.acceleration == b.acceleration, "BEM backend identity changed");
    require(a.device == b.device && a.triangleCount == b.triangleCount &&
                a.pressureDofCount == b.pressureDofCount &&
                a.peakDeviceBytes == b.peakDeviceBytes && a.liveDeviceBytes == b.liveDeviceBytes &&
                a.frequencies.size() == b.frequencies.size(),
            "BEM result metadata changed");
    for (std::size_t index = 0; index < a.frequencies.size(); ++index) {
        const auto &x = a.frequencies[index], &y = b.frequencies[index];
        close(x.frequency, y.frequency, "Result frequency changed");
        close(x.assemblyMilliseconds, y.assemblyMilliseconds, "Assembly timing changed");
        close(x.solveMilliseconds, y.solveMilliseconds, "Solve timing changed");
        close(x.factorizationMilliseconds, y.factorizationMilliseconds, "LU timing changed");
        close(x.triangularSolveMilliseconds, y.triangularSolveMilliseconds, "RHS timing changed");
        close(x.fieldMilliseconds, y.fieldMilliseconds, "Field timing changed");
        close(x.residualCheckMilliseconds, y.residualCheckMilliseconds, "Residual timing changed");
        close(x.interactionPlanMilliseconds, y.interactionPlanMilliseconds,
              "Host planning timing changed");
        require(x.interactionPlanReused == y.interactionPlanReused,
                "Host plan cache-hit flag changed");
        require(x.fastMultipole.nodeCount == y.fastMultipole.nodeCount &&
                    x.fastMultipole.leafCount == y.fastMultipole.leafCount &&
                    x.fastMultipole.farInteractionCount == y.fastMultipole.farInteractionCount &&
                    x.fastMultipole.pressurePressureNearPairs ==
                        y.fastMultipole.pressurePressureNearPairs &&
                    x.fastMultipole.pressurePressureFarPairs ==
                        y.fastMultipole.pressurePressureFarPairs &&
                    x.fastMultipole.pressureFluxNearPairs ==
                        y.fastMultipole.pressureFluxNearPairs &&
                    x.fastMultipole.pressureFluxFarPairs == y.fastMultipole.pressureFluxFarPairs &&
                    x.fastMultipole.historicalChebyshevStatistics ==
                        y.fastMultipole.historicalChebyshevStatistics &&
                    x.fastMultipole.graphExecution == y.fastMultipole.graphExecution,
                "Fast Multipole hierarchy diagnostics changed");
        require(x.excitations.size() == y.excitations.size(), "Result RHS count changed");
        for (std::size_t excitation = 0; excitation < x.excitations.size(); ++excitation) {
            const auto &u = x.excitations[excitation], &v = y.excitations[excitation];
            require(u.name == v.name && u.refinementSteps == v.refinementSteps &&
                        u.iterations == v.iterations &&
                        u.boundary()->surfacePressure == v.boundary()->surfacePressure &&
                        u.observerPressure == v.observerPressure &&
                        u.boundary()->pressureTrace == v.boundary()->pressureTrace &&
                        u.boundary()->normalDerivative == v.boundary()->normalDerivative &&
                        u.impedances.size() == v.impedances.size(),
                    "Complex pressure result changed");
            close(u.boundaryResidual, v.boundaryResidual, "Boundary residual changed");
            require(u.residualEvaluation == v.residualEvaluation,
                    "Residual evaluation metadata changed");
            close(u.fastMultipoleKrylovResidual, v.fastMultipoleKrylovResidual,
                  "Krylov residual changed");
            close(u.fastMultipoleKrylovTarget, v.fastMultipoleKrylovTarget,
                  "Configured inner Krylov target changed");
            close(u.fastMultipoleOperatorDiscrepancy, v.fastMultipoleOperatorDiscrepancy,
                  "Operator discrepancy changed");
            require(u.historicalChebyshevEvidence == v.historicalChebyshevEvidence &&
                        u.fastMultipoleComparisonLevels == v.fastMultipoleComparisonLevels &&
                        u.fastMultipoleVerificationLevels == v.fastMultipoleVerificationLevels &&
                        u.fastMultipoleAccuracyRefinements == v.fastMultipoleAccuracyRefinements,
                    "Automatic level evidence or historical provenance changed");
            close(u.fastMultipoleAccuracySeed, v.fastMultipoleAccuracySeed,
                  "Accuracy seed changed");
            close(u.fastMultipoleComparisonAccuracySeed, v.fastMultipoleComparisonAccuracySeed,
                  "Comparison seed changed");
            close(u.fastMultipoleVerificationAccuracySeed, v.fastMultipoleVerificationAccuracySeed,
                  "Verification seed changed");
            require(u.farQuadratureComparisonDegree == v.farQuadratureComparisonDegree &&
                        u.farQuadratureVerificationDegree == v.farQuadratureVerificationDegree &&
                        u.farQuadratureRefinements == v.farQuadratureRefinements &&
                        u.fieldQuadratureComparisonDegree == v.fieldQuadratureComparisonDegree &&
                        u.fieldQuadratureVerificationDegree == v.fieldQuadratureVerificationDegree,
                    "Quadrature-rule evidence changed");
            for (const auto [before, after] :
                 {std::pair{u.farQuadratureDiscrepancy, v.farQuadratureDiscrepancy},
                  std::pair{u.farQuadraturePreparationMilliseconds,
                            v.farQuadraturePreparationMilliseconds},
                  std::pair{u.farQuadraturePressureMilliseconds,
                            v.farQuadraturePressureMilliseconds},
                  std::pair{u.farQuadratureFluxMilliseconds, v.farQuadratureFluxMilliseconds},
                  std::pair{u.fieldQuadratureDiscrepancy, v.fieldQuadratureDiscrepancy}})
                close(before, after, "Quadrature comparison metric changed");
            require(bool(u.lemSolution) == bool(v.lemSolution) &&
                        u.solvedVelocities.size() == v.solvedVelocities.size(),
                    "Coupled solution or physical velocity data disappeared");
            if (u.lemSolution)
                sameLemSolution(*u.lemSolution, *v.lemSolution);
            for (std::size_t i = 0; i < u.solvedVelocities.size(); ++i)
                require(u.solvedVelocities[i].faceId == v.solvedVelocities[i].faceId &&
                            u.solvedVelocities[i].normalVelocity ==
                                v.solvedVelocities[i].normalVelocity,
                        "A coupled face velocity changed");
            for (std::size_t impedance = 0; impedance < u.impedances.size(); ++impedance) {
                const auto &s = u.impedances[impedance], &t = v.impedances[impedance];
                require(s.faceId == t.faceId, "Impedance face identity changed");
                close(s.area, t.area, "Impedance area changed");
                sameComplex(s.meanPressure, t.meanPressure, "Impedance mean pressure changed");
                sameComplex(s.volumeVelocity, t.volumeVelocity, "Impedance flow changed");
                sameComplex(s.volumeImpedance, t.volumeImpedance,
                            "Complex volume impedance changed");
            }
        }
    }
}

void sameProject(const Project &actual, const Project &expected) {
    require(actual.draftFaceColors == expected.draftFaceColors &&
                actual.cadFaceColors == expected.cadFaceColors,
            "Display face colors changed");
    require(actual.source == expected.source && actual.revision == expected.revision &&
                actual.geometryRevision == expected.geometryRevision &&
                actual.meshRevision == expected.meshRevision,
            "Project source or exact uint64 revisions changed");
    sameDraft(actual.parameters, expected.parameters);
    close(actual.meshing.targetEdgeLength, expected.meshing.targetEdgeLength,
          "Target mesh size changed");
    require(actual.meshing.quality.enabled == expected.meshing.quality.enabled,
            "Mesh quality enforcement changed");
    close(actual.meshing.quality.minimumAngleDegrees, expected.meshing.quality.minimumAngleDegrees,
          "Mesh minimum-angle constraint changed");
    close(actual.meshing.quality.maximumAspectRatio, expected.meshing.quality.maximumAspectRatio,
          "Mesh aspect-ratio constraint changed");
    require(actual.draftFaceTargetEdgeLengths == expected.draftFaceTargetEdgeLengths &&
                actual.cadFaceTargetEdgeLengths == expected.cadFaceTargetEdgeLengths,
            "Local face refinement changed");
    require(actual.draftReversedFaceIds == expected.draftReversedFaceIds &&
                actual.cadReversedFaceIds == expected.cadReversedFaceIds,
            "Physical face orientation masks changed");
    require(actual.meshing.adaptive.enabled == expected.meshing.adaptive.enabled,
            "Adaptive meshing mode changed");
    close(actual.meshing.adaptive.minimumEdgeLength, expected.meshing.adaptive.minimumEdgeLength,
          "Adaptive minimum mesh size changed");
    close(actual.meshing.adaptive.geometricTolerance, expected.meshing.adaptive.geometricTolerance,
          "Adaptive geometric tolerance changed");
    close(actual.meshing.adaptive.growthRate, expected.meshing.adaptive.growthRate,
          "Adaptive size growth changed");
    require(bool(actual.preview) == bool(expected.preview) &&
                bool(actual.mesh) == bool(expected.mesh) && bool(actual.cad) == bool(expected.cad),
            "Project snapshot presence changed");
    if (expected.preview)
        sameMesh(*actual.preview, *expected.preview);
    if (expected.mesh) {
        sameMesh(actual.mesh->mesh, expected.mesh->mesh);
        require(actual.mesh->quality.watertight(), "Restored mesh quality was not validated");
        const auto quality = inspectMesh(actual.mesh->mesh);
        close(actual.mesh->quality.minimumAngleDegrees, quality.minimumAngleDegrees,
              "Restored mesh quality was not recomputed");
        close(actual.mesh->quality.maximumAspectRatio, quality.maximumAspectRatio,
              "Restored aspect ratio was not recomputed");
    }
    if (expected.cad)
        sameCad(*actual.cad, *expected.cad);
    const auto &a = actual.simulation, &b = expected.simulation;
    require(a.measurements == b.measurements, "Imported measurements changed");
    require(a.axisymmetricSampling == b.axisymmetricSampling,
            "The immutable meridian sampling changed");
    sameSettings(a.settings, b.settings);
    require(a.revision == b.revision && a.current == b.current && bool(a.mesh) == bool(b.mesh) &&
                a.request.has_value() == b.request.has_value() && bool(a.result) == bool(b.result),
            "Simulation archive state changed");
    if (b.mesh)
        sameMesh(*a.mesh, *b.mesh);
    if (b.request)
        sameRequest(*a.request, *b.request);
    if (b.result)
        sameResult(*a.result, *b.result);
}

QJsonValue replaceAt(QJsonValue value, const QStringList &path, qsizetype position,
                     const QJsonValue &replacement) {
    if (position == path.size())
        return replacement;
    if (value.isObject()) {
        auto object = value.toObject();
        require(object.contains(path[position]),
                "Corruption fixture references a missing schema key");
        object[path[position]] =
            replaceAt(object.value(path[position]), path, position + 1, replacement);
        return object;
    }
    require(value.isArray(), "Corruption fixture path is not an object or array");
    auto array = value.toArray();
    bool valid = false;
    const auto index = path[position].toLongLong(&valid);
    require(valid && index >= 0 && index < array.size(),
            "Corruption fixture array index is invalid");
    array[index] = replaceAt(array[index], path, position + 1, replacement);
    return array;
}

void roundTripTests(const QTemporaryDir &directory) {
    for (const auto source : {ProjectSource::draft, ProjectSource::freecad}) {
        const auto project = fullTestProject(source);
        const auto path =
            directory.filePath(source == ProjectSource::draft ? "draft.json" : "cad.json");
        writeProject(path, project);
        const auto restoredCurrent = readProject(path);
        sameProject(restoredCurrent, project);
        require(restoredCurrent.simulation.current &&
                    restoredCurrent.simulation.resultComplete == std::optional{true},
                "A legacy current archive was not recognized as a completed solve");
        const auto root = QJsonDocument::fromJson(bytes(path)).object();
        require(root.value("revision").isString() && root.value("geometryRevision").isString() &&
                    root.value("meshRevision").isString(),
                "uint64 revisions are not represented exactly");
        require(root.value("revision").toString() == QString::number(largeRevision),
                "Revision above 2^53 lost precision in the file");
    }

    {
        auto project = fullTestProject(ProjectSource::draft);
        project.simulation.current = false;
        project.simulation.settings.acousticRegions = {{7, "Kammer", {1, 2}, {4}},
                                                       {8, "Hals", {3}, {4}}};
        const auto path = directory.filePath("regions.json");
        writeProject(path, project);
        require(readProject(path).simulation.settings.acousticRegions ==
                    project.simulation.settings.acousticRegions,
                "Acoustic region definitions did not survive project round trip");
    }
    {
        auto previous = *testCad(testMesh());
        auto current = previous;
        const auto id = previous.faces.front().id;
        std::vector<AcousticRegionDefinition> regions{{7, "Innen", {id}, {}, false}};
        require(!reconciledCadAcousticRegions(regions, &previous, &current)
                     .front()
                     .needsAssignmentReview,
                "Stable CAD region was flagged for review");
        current.documentId += "-different";
        const auto changed = reconciledCadAcousticRegions(regions, &previous, &current);
        require(changed.front().needsAssignmentReview &&
                    changed.front().faceIds == regions.front().faceIds,
                "Document change silently reused a region face ID");
        current = previous;
        current.faces.clear();
        require(reconciledCadAcousticRegions(regions, &previous, &current)
                    .front()
                    .needsAssignmentReview,
                "Missing CAD region face was accepted");
        auto project = fullTestProject(ProjectSource::draft);
        project.simulation.current = false;
        project.simulation.settings.acousticRegions = changed;
        const auto path = directory.filePath("regions-review.json");
        writeProject(path, project);
        require(readProject(path).simulation.settings.acousticRegions.front().needsAssignmentReview,
                "CAD region review obligation did not survive project round trip");
    }
    auto stale = fullTestProject(ProjectSource::freecad);
    stale.revision += 3;
    stale.geometryRevision = stale.meshRevision = stale.revision;
    const auto changed = std::make_shared<SurfaceMesh>(testMesh(1.07));
    stale.preview = changed;
    stale.mesh = std::make_shared<MeshingResult>(MeshingResult{*changed, inspectMesh(*changed)});
    stale.cad = testCad(*changed);
    stale.simulation.current = false;
    const auto stalePath = directory.filePath("stale.json");
    writeProject(stalePath, stale);
    const auto restored = readProject(stalePath);
    sameProject(restored, stale);
    require(restored.simulation.mesh->vertices.front().x !=
                    restored.mesh->mesh.vertices.front().x &&
                restored.simulation.revision < restored.revision && !restored.simulation.current,
            "Stale result was silently attached to the current mesh or marked current");

    for (const auto complete : {false, true}) {
        auto tagged = stale;
        tagged.simulation.resultComplete = complete;
        const auto taggedPath =
            directory.filePath(complete ? "complete-result.json" : "partial-result.json");
        writeProject(taggedPath, tagged);
        require(readProject(taggedPath).simulation.resultComplete == std::optional{complete},
                "Result completeness was not preserved in the archive");
    }
    require(!restored.simulation.resultComplete,
            "A legacy stale result without completeness metadata was guessed complete or partial");
    auto invalidComplete = fullTestProject(ProjectSource::freecad);
    invalidComplete.simulation.resultComplete = false;
    expectRejected(
        [&] { writeProject(directory.filePath("invalid-partial.json"), invalidComplete); },
        "An incomplete result was marked current");

    auto pending = stale;
    pending.simulation.settings.normalizeToOnAxis = false;
    auto cad = std::make_shared<CadSnapshot>(*pending.cad);
    cad->faces[0].state = CadFaceState::needsReview;
    cad->faces[1].state = CadFaceState::newFace;
    cad->unresolvedPreviousFaces = {81, 82};
    pending.cad = cad;
    const auto pendingPath = directory.filePath("pending-cad.json");
    writeProject(pendingPath, pending);
    const auto restoredPending = readProject(pendingPath);
    sameProject(restoredPending, pending);
    require(restoredPending.cad->hasUnresolvedAssignments(),
            "Pending CAD assignments disappeared offline");

    Project parametersOnly;
    const auto plainPath = directory.filePath("parameters-only.json");
    writeProject(plainPath, parametersOnly);
    sameProject(readProject(plainPath), parametersOnly);
}

lem::Network testLemNetwork() {
    // Deliberately generic in-vacuo parameters; this fixture is not a driver specification.
    const lem::DriverParameters parameters{6.4, .0007, 4.98, .004, .0016, .262, .00544};
    auto network = lem::makeDriverNetwork(parameters, {{1}, {2}}, {12, 3}, .0004);
    std::get<lem::BemPort>(network.elements.back().component).polarity = -1;
    return network;
}

Project coupledTestProject(bool silent) {
    auto project = fullTestProject(ProjectSource::freecad);
    auto &archive = project.simulation;
    archive.settings.minimumFrequency = 300;
    archive.settings.maximumFrequency = 600;
    archive.settings.frequencyCount = 2;
    archive.settings.couplingEnabled = true;
    archive.settings.couplingNetwork = testLemNetwork();
    archive.settings.farField = false;
    archive.settings.observationDistance = 1.5;
    if (silent)
        std::get<lem::EffortSource>(archive.settings.couplingNetwork.elements.front().component)
            .amplitude = {};
    auto &request = *archive.request;
    request.frequencies = {300, 600};
    request.excitations.clear();
    request.coupling = archive.settings.couplingNetwork;
    request.farField.reset();
    auto result = std::make_shared<BemResult>();
    result->device =
        "SYNTHETIC BEM PERSISTENCE FIXTURE with actual CPU LEM solution - no GPU computation";
    result->triangleCount = archive.mesh->triangles.size();
    result->pressureDofCount = archive.mesh->vertices.size();
    const lem::AcousticLoad load{{{18000, 3000}, {3000, -1000}, {3000, -1000}, {24000, 7000}}, {}};
    for (const auto frequency : request.frequencies) {
        BemFrequencyResult frame;
        frame.frequency = frequency;
        BemExcitationResult excitation;
        BemBoundaryData boundary;
        excitation.name = "LEM";
        excitation.lemSolution = lem::solve(*request.coupling, frequency, load,
                                            {request.medium.density, request.medium.soundSpeed});
        boundary.surfacePressure.resize(archive.mesh->triangles.size());
        excitation.observerPressure.resize(request.observers.size());
        const auto ports = lem::bemPorts(*request.coupling);
        for (std::size_t i = 0; i < ports.size(); ++i) {
            const auto &port = std::get<lem::BemPort>(ports[i]->component);
            const auto faceId = port.faceIds.front();
            const auto &triangle = archive.mesh->triangles[faceId - 1];
            const auto a = archive.mesh->vertices[triangle.vertices[0]],
                       b = archive.mesh->vertices[triangle.vertices[1]],
                       c = archive.mesh->vertices[triangle.vertices[2]];
            const double area = length(cross(subtract(b, a), subtract(c, a))) / 2;
            const auto &physical = excitation.lemSolution->ports[i];
            boundary.surfacePressure[faceId - 1] = physical.pressure;
            excitation.solvedVelocities.push_back({faceId, physical.volumeVelocity / area});
            if (physical.volumeVelocity != lem::Complex{})
                excitation.impedances.push_back({faceId, area, physical.pressure,
                                                 physical.volumeVelocity,
                                                 physical.pressure / physical.volumeVelocity});
        }
        excitation.setBoundary(std::move(boundary));
        frame.excitations.push_back(std::move(excitation));
        result->frequencies.push_back(std::move(frame));
    }
    archive.result = result;
    return project;
}

void coupledRoundTrip(const QTemporaryDir &directory) {
    const auto path = directory.filePath("coupled.json");
    for (const bool silent : {false, true}) {
        const auto project = coupledTestProject(silent);
        writeProject(path, project);
        sameProject(readProject(path), project);
    }
    auto stale = coupledTestProject(false);
    ++stale.revision;
    stale.simulation.current = false;
    std::get<lem::EffortSource>(
        stale.simulation.settings.couplingNetwork.elements.front().component)
        .amplitude = {5, 0};
    writeProject(path, stale);
    sameProject(readProject(path), stale);
    require(readProject(path).simulation.settings.couplingNetwork !=
                *readProject(path).simulation.request->coupling,
            "An edited circuit overwrote the immutable archived circuit");

    const auto project = coupledTestProject(false);
    writeProject(path, project);
    const auto original = QJsonDocument::fromJson(bytes(path)).object();
    const auto damaged = directory.filePath("bad-coupled.json");
    const auto reject = [&](const QStringList &keys, const QJsonValue &value) {
        saveBytes(damaged, QJsonDocument(replaceAt(original, keys, 0, value).toObject()).toJson());
        expectRejected([&] { (void)readProject(damaged); },
                       "An inconsistent coupled archive was accepted");
    };
    const QStringList solutionPath{"simulation",  "result", "frequencies", "0",
                                   "excitations", "0",      "lemSolution"};
    const auto rejectSolution = [&](QStringList keys, const QJsonValue &value) {
        reject(solutionPath + keys, value);
    };
    rejectSolution({"frequency"}, 299);
    rejectSolution({"relativeResidual"}, 1);
    rejectSolution({"nodes", "0", "id"}, 2);
    rejectSolution({"nodes", "0", "effort"}, QJsonArray{1, 0});
    rejectSolution({"elements", "0", "activePower"}, 999);
    rejectSolution({"drivers", "0", "id"}, 1);
    rejectSolution({"drivers", "0", "current"}, QJsonArray{999, 0});
    rejectSolution({"drivers", "0", "displacement"}, QJsonArray{1, 0});
    rejectSolution({"drivers", "0", "copperPower"}, -1);
    rejectSolution({"drivers", "0", "acousticPower"}, QJsonValue::Null);
    rejectSolution({"ports", "0", "id"}, 2);
    rejectSolution({"ports", "0", "pressure"}, QJsonArray{1, 0});
    rejectSolution({"ports", "0", "volumeVelocity"}, QJsonArray{1, 0});
    rejectSolution({"ports"}, QJsonArray{});
    reject({"simulation", "result", "frequencies", "0", "excitations", "0", "solvedVelocities", "0",
            "normalVelocity"},
           QJsonArray{1, 0});
    reject({"simulation", "request", "coupling", "elements", "3", "component", "faceIds"},
           QJsonArray{3});
    reject({"simulation", "request", "coupling", "elements", "1", "component", "parameters",
            "movingMassInVacuo"},
           0);
    reject({"simulation", "request", "coupling", "nodes", "0", "reference"}, false);
    reject({"simulation", "settings", "observationDistance"}, 0);
    reject({"simulation", "settings", "farField"}, "false");
}

void unfinishedLemGraphRoundTrip(const QTemporaryDir &directory) {
    Project project;
    auto &settings = project.simulation.settings;
    settings.couplingEnabled = true;
    settings.couplingNetwork = testLemNetwork();
    auto &network = settings.couplingNetwork;
    network.nodes.push_back(
        {{100}, lem::Domain::mechanical, true, "mechanical reference", {-10, 17}});
    network.nodes.push_back({{101}, lem::Domain::mechanical, false, "mechanics", {5, 35}});
    network.nodes.push_back({{102}, lem::Domain::acoustic, false, "unwired terminal", {19, 31}});
    network.elements.push_back({{100}, "R", {7, 9}, lem::Resistance{{{101}, {100}}, 4}});
    network.elements.push_back({{101}, "M", {7, 19}, lem::Inertance{{{101}, {100}}, .012}});
    network.elements.push_back({{102}, "C", {7, 29}, lem::Compliance{{{101}, {100}}, .0005}});
    network.elements.push_back(
        {{103}, "flow", {7, 39}, lem::FlowSource{{{100}, {101}}, {.1, -.2}}});
    std::get<lem::BemPort>(network.elements[3].component).faceIds.clear();
    std::get<lem::Driver>(network.elements[1].component).parameters.pistonArea = 0;
    const auto path = directory.filePath("unfinished-lem.json");
    writeProject(path, project);
    sameProject(readProject(path), project);
    expectRejected([&] { (void)makeSimulationRequest(settings, {}); },
                   "An unfinished stored graph was accepted for simulation");
    const auto original = QJsonDocument::fromJson(bytes(path)).object();
    const auto bad = directory.filePath("bad-lem-graph.json");
    const auto reject = [&](const QStringList &keys, const QJsonValue &value) {
        saveBytes(
            bad,
            QJsonDocument(replaceAt(original,
                                    QStringList{"simulation", "settings", "couplingNetwork"} + keys,
                                    0, value)
                              .toObject())
                .toJson());
        expectRejected([&] { (void)readProject(bad); },
                       "An invalid LEM graph structure was accepted");
    };
    reject({"nodes", "0", "id"}, 2);
    reject({"nodes", "0", "domain"}, "other");
    reject({"nodes", "0", "domain"}, "acoustic");
    reject({"elements", "0", "id"}, 2);
    reject({"elements", "1", "component", "front"}, 99999);
    reject({"elements", "1", "component", "parameters", "forceFactor"}, QJsonValue::Null);
    reject({"elements", "3", "component", "polarity"}, 0);
}

lem::Network interfaceNetwork(double volume = 5e-6) {
    lem::Network network;
    network.nodes = {{{1}, lem::Domain::acoustic, false, "Chamber", {}},
                     {{2}, lem::Domain::acoustic, true, "Reference", {}}};
    const lem::Terminals terminals{{1}, {2}};
    network.elements = {{{1},
                         "Virtual interface",
                         {},
                         lem::BemPort{terminals, {3}, 1, lem::BemPortKind::interfacePort}},
                        {{2}, "Chamber assumption", {}, lem::Cavity{terminals, volume}}};
    return network;
}

Project interfaceTestProject(bool silent) {
    auto project = fullTestProject(ProjectSource::freecad);
    auto &archive = project.simulation;
    auto &settings = archive.settings;
    settings.minimumFrequency = settings.maximumFrequency = 1000;
    settings.frequencyCount = 1;
    settings.velocityScale = silent ? 0 : .7;
    settings.phaseDegrees = 37;
    settings.couplingEnabled = true;
    settings.couplingNetwork = interfaceNetwork();
    settings.infiniteBaffle = InfiniteBaffle{{}, {0, 0, 1}, {4}};
    settings.farField = true;
    const std::vector<FaceVelocity> known{{1, {.3, .6}}, {2, {.15, .3}}};
    archive.request = makeSimulationRequest(settings, known);
    const auto &request = *archive.request;
    auto result = std::make_shared<BemResult>();
    result->device =
        "SYNTHETIC INTERFACE PERSISTENCE FIXTURE with actual CPU LEM solution - no GPU computation";
    result->triangleCount = archive.mesh->triangles.size();
    result->pressureDofCount = archive.mesh->vertices.size();
    BemFrequencyResult frame;
    frame.frequency = 1000;
    BemExcitationResult excitation;
    BemBoundaryData boundary;
    excitation.name = "LEM";
    excitation.lemSolution =
        lem::solve(*request.coupling, frame.frequency, {{{18000, 3000}}, {{2, -.5}}},
                   {request.medium.density, request.medium.soundSpeed});
    boundary.surfacePressure.resize(archive.mesh->triangles.size(), {1, -.2});
    excitation.observerPressure.resize(request.observers.size());
    excitation.solvedVelocities = request.excitations.front().faces;
    const auto area = [&](std::uint32_t id) {
        const auto &t = archive.mesh->triangles[id - 1];
        const auto a = archive.mesh->vertices[t.vertices[0]],
                   b = archive.mesh->vertices[t.vertices[1]],
                   c = archive.mesh->vertices[t.vertices[2]];
        return length(cross(subtract(b, a), subtract(c, a))) / 2;
    };
    const auto &port = excitation.lemSolution->ports.front();
    boundary.surfacePressure[2] = port.pressure;
    excitation.solvedVelocities.push_back({3, port.volumeVelocity / area(3)});
    for (const auto &velocity : excitation.solvedVelocities) {
        if (velocity.normalVelocity == lem::Complex{})
            continue;
        const auto pressure = boundary.surfacePressure[velocity.faceId - 1];
        const auto volumeVelocity = velocity.normalVelocity * area(velocity.faceId);
        excitation.impedances.push_back({velocity.faceId, area(velocity.faceId), pressure,
                                         volumeVelocity, pressure / volumeVelocity});
    }
    excitation.setBoundary(std::move(boundary));
    frame.excitations.push_back(std::move(excitation));
    result->frequencies.push_back(std::move(frame));
    archive.result = result;
    return project;
}

void interfaceRoundTrip(const QTemporaryDir &directory) {
    const auto path = directory.filePath("interface.hsproj");
    for (const bool silent : {false, true}) {
        const auto project = interfaceTestProject(silent);
        writeProject(path, project);
        sameProject(readProject(path), project);
        const auto &request = *project.simulation.request;
        require(request.coupling && request.excitations.size() == 1 &&
                    request.excitations.front().faces.size() == 2,
                "Interface coupling lost the known dome/surround excitation");
        sameComplex(request.excitations[0].faces[1].normalVelocity,
                    .5 * request.excitations[0].faces[0].normalVelocity,
                    "The prescribed surround/dome velocity ratio changed");
    }
    auto settings = interfaceTestProject(false).simulation.settings;
    const std::vector<FaceVelocity> known{{1, {1, 0}}};
    const auto expected =
        std::polar(settings.velocityScale, settings.phaseDegrees * std::acos(-1.) / 180);
    sameComplex(
        makeSimulationRequest(settings, known).excitations.front().faces.front().normalVelocity,
        expected, "Interface coupling ignored known velocity scale/phase");
    settings.couplingEnabled = false;
    expectRejected([&] { (void)makeSimulationRequest(settings, known); },
                   "Disabling coupling silently closed a virtual interface");
    auto &port = std::get<lem::BemPort>(settings.couplingNetwork.elements[0].component);
    port.faceIds.clear();
    port.needsAssignmentReview = true;
    expectRejected([&] { (void)makeSimulationRequest(settings, known); },
                   "A lost interface assignment silently became Wall while disabled");
    port.needsAssignmentReview = false;
    require(makeSimulationRequest(settings, known).excitations.size() == 1,
            "A new unbound editor port blocked ordinary BEM");
    settings.couplingEnabled = true;
    settings.couplingNetwork = interfaceNetwork();
    settings.infiniteBaffle.reset();
    expectRejected([&] { (void)makeSimulationRequest(settings, known); },
                   "An interface was accepted by exterior BEM");
    settings.infiniteBaffle = InfiniteBaffle{{}, {0, 0, 1}, {4}};
    auto other = settings.couplingNetwork.elements.front();
    other.id = {3};
    std::get<lem::BemPort>(other.component) = {{{1}, {2}}, {1}, 1};
    settings.couplingNetwork.elements.push_back(other);
    require(makeSimulationRequest(settings, known).excitations.empty(),
            "A mixed baffled network prescribed an already port-driven face twice");
    const std::vector<FaceVelocity> partlyBound{{1, {1, 0}}, {2, {.5, 0}}};
    const auto mixedRequest = makeSimulationRequest(settings, partlyBound);
    require(mixedRequest.excitations.size() == 1 && mixedRequest.excitations[0].faces.size() == 1 &&
                mixedRequest.excitations[0].faces[0].faceId == 2,
            "A mixed baffled network lost or duplicated its remaining known Driven face");
    settings.infiniteBaffle.reset();
    settings.couplingNetwork = testLemNetwork();
    require(makeSimulationRequest(settings, known).excitations.empty(),
            "Old driven-port coupling acquired a prescribed excitation");

    Project pending;
    pending.simulation.settings.couplingEnabled = true;
    pending.simulation.settings.couplingNetwork = interfaceNetwork(0);
    pending.simulation.settings.infiniteBaffle = InfiniteBaffle{{}, {0, 0, 1}, {4}};
    writeProject(path, pending);
    sameProject(readProject(path), pending);
    expectRejected([&] { (void)makeSimulationRequest(pending.simulation.settings, known); },
                   "An unset chamber volume was assigned a usable default");
    auto &pendingPort =
        std::get<lem::BemPort>(pending.simulation.settings.couplingNetwork.elements[0].component);
    pendingPort.faceIds.clear();
    pendingPort.needsAssignmentReview = true;
    pending.simulation.settings.couplingEnabled = false;
    writeProject(path, pending);
    sameProject(readProject(path), pending);

    const auto before = testCad(testMesh());
    auto after = *before;
    ++after.revision;
    std::reverse(after.faces.begin(), after.faces.end());
    const auto assigned = interfaceNetwork();
    require(reconciledCadBemAssignments(assigned, before.get(), &after) == assigned,
            "Reordering uniquely mapped CAD faces lost a virtual interface");
    const auto requireReview = [&](const CadSnapshot *previous, const CadSnapshot *current) {
        const auto reconciled = reconciledCadBemAssignments(assigned, previous, current);
        const auto &changed = std::get<lem::BemPort>(reconciled.elements[0].component);
        require(changed.faceIds.empty() && changed.needsAssignmentReview,
                "Lost/ambiguous CAD interface was reused or its review obligation disappeared");
    };
    for (const auto state : {CadFaceState::needsReview, CadFaceState::newFace}) {
        auto changed = after;
        for (auto &face : changed.faces)
            if (face.id == 3)
                face.state = state;
        requireReview(before.get(), &changed);
    }
    auto changed = after;
    for (auto &face : changed.faces)
        if (face.id == 3)
            face.boundary = BoundaryCondition::driven;
    requireReview(before.get(), &changed);
    changed = after;
    std::erase_if(changed.faces, [](const auto &face) { return face.id == 3; });
    requireReview(before.get(), &changed);
    changed = after;
    changed.documentId += "-new";
    requireReview(before.get(), &changed);
    requireReview(nullptr, &after);
    requireReview(before.get(), nullptr);
    require(reconciledCadBemAssignments(pending.simulation.settings.couplingNetwork, before.get(),
                                        &after) == pending.simulation.settings.couplingNetwork,
            "A later geometry match silently cleared persistent assignment review");

    const auto fixture = interfaceTestProject(false);
    writeProject(path, fixture);
    const auto json = QJsonDocument::fromJson(bytes(path)).object();
    const auto bad = directory.filePath("bad-interface.hsproj");
    const auto reject = [&](QStringList keys, QJsonValue value) {
        saveBytes(bad, QJsonDocument(replaceAt(json, keys, 0, value).toObject()).toJson());
        expectRejected([&] { (void)readProject(bad); },
                       "Invalid BEM↔LEM interface archive accepted");
    };
    for (const auto section : {"settings", "request"}) {
        const QStringList prefix{
            "simulation",
            section,
            section == QStringLiteral("settings") ? "couplingNetwork" : "coupling",
            "elements",
            "0",
            "component"};
        reject(prefix + QStringList{"kind"}, "unknown");
        reject(prefix + QStringList{"kind"}, 1);
        reject(prefix + QStringList{"needsAssignmentReview"}, "false");
        reject(prefix + QStringList{"needsAssignmentReview"}, 0);
    }
    reject({"simulation", "request", "coupling", "elements", "0", "component", "faceIds"},
           QJsonArray{1});
    reject({"simulation", "request", "coupling", "elements", "0", "component", "faceIds"},
           QJsonArray{4});
    reject({"simulation", "request", "coupling", "elements", "0", "component",
            "needsAssignmentReview"},
           true);
    reject({"simulation", "result", "frequencies", "0", "excitations", "0", "solvedVelocities", "0",
            "normalVelocity"},
           QJsonArray{1, 0});
    reject({"simulation", "result", "frequencies", "0", "excitations", "0", "solvedVelocities", "2",
            "normalVelocity"},
           QJsonArray{1, 0});
    reject({"simulation", "settings", "couplingEnabled"}, false);
    reject({"simulation", "settings", "couplingNetwork", "elements", "0", "component", "kind"},
           "drivenSurface");
    reject({"simulation", "settings", "couplingNetwork", "elements", "0", "component",
            "needsAssignmentReview"},
           true);
    reject({"simulation", "settings", "couplingNetwork", "elements", "1", "component", "volume"},
           6e-6);
    auto stale = fixture;
    stale.simulation.current = false;
    std::get<lem::Cavity>(stale.simulation.settings.couplingNetwork.elements[1].component).volume =
        6e-6;
    writeProject(path, stale);
    sameProject(readProject(path), stale);

    // Missing role fields retain the existing driven-surface semantics for old files.
    auto old = coupledTestProject(false);
    writeProject(path, old);
    auto legacy = QJsonDocument::fromJson(bytes(path)).object();
    for (const auto section : {"settings", "request"}) {
        const QStringList prefix{
            "simulation",
            section,
            section == QStringLiteral("settings") ? "couplingNetwork" : "coupling",
            "elements",
            "3",
            "component"};
        auto component =
            legacy["simulation"]
                .toObject()[section]
                .toObject()[section == QStringLiteral("settings") ? "couplingNetwork" : "coupling"]
                .toObject()["elements"]
                .toArray()[3]
                .toObject()["component"]
                .toObject();
        component.remove("kind");
        component.remove("needsAssignmentReview");
        component.remove("velocityWeights");
        legacy = replaceAt(legacy, prefix, 0, component).toObject();
    }
    saveBytes(path, QJsonDocument(legacy).toJson());
    sameProject(readProject(path), old);
}

Project weightedPortTestProject(bool baffled, bool withSlit = false, bool closedChamber = false) {
    // The BEM samples below only exercise archive invariants. The CPU MNA solution
    // is real, but these synthetic fields are not acoustic validation data.
    auto project = interfaceTestProject(false);
    auto &archive = project.simulation;
    auto &settings = archive.settings;
    lem::BemPort driven{{{1}, {2}}, {1, 2}, 1};
    driven.velocityWeights = {{2, .5}};
    settings.couplingNetwork.elements.push_back({{3}, "Modal drive", {}, driven});
    if (withSlit) {
        settings.couplingNetwork.nodes.push_back(
            {{3}, lem::Domain::acoustic, false, "Behind slit", {}});
        std::get<lem::Cavity>(settings.couplingNetwork.elements[1].component).terminals.positive = {
            3};
        settings.couplingNetwork.elements.push_back(
            {{4},
             "Real channel behind interface",
             {},
             lem::AcousticSlit{{{1}, {3}}, .03, 150e-6, .00175, 18.4e-6}});
    }
    if (closedChamber) {
        const auto terminals =
            std::get<lem::Cavity>(settings.couplingNetwork.elements[1].component).terminals;
        settings.couplingNetwork.elements[1].component =
            lem::ClosedSlitCavity{terminals, .03, 150e-6, .00175, 18.4e-6, 1.4, .71};
    }
    if (!baffled) {
        settings.infiniteBaffle.reset();
        settings.couplingNetwork.elements.erase(settings.couplingNetwork.elements.begin());
    }
    archive.request =
        makeSimulationRequest(settings, std::vector<FaceVelocity>{{1, {1, 0}}, {2, {.5, 0}}});
    auto result = std::make_shared<BemResult>();
    result->device =
        "SYNTHETIC WEIGHTED-PORT ARCHIVE with actual CPU LEM solution - no GPU computation";
    result->triangleCount = archive.mesh->triangles.size();
    result->pressureDofCount = archive.mesh->vertices.size();
    BemFrequencyResult frequency;
    frequency.frequency = 1000;
    BemExcitationResult output;
    BemBoundaryData boundary;
    output.name = "LEM";
    const auto ports = lem::bemPorts(settings.couplingNetwork);
    const lem::AcousticLoad load =
        baffled ? lem::AcousticLoad{{{18000, 3000}, {2000, -500}, {2000, -500}, {25000, 4000}},
                                    {{2, -.5}, {1, .3}}}
                : lem::AcousticLoad{{{25000, 4000}}, {{1, .3}}};
    output.lemSolution = lem::solve(settings.couplingNetwork, 1000, load,
                                    {settings.medium.density, settings.medium.soundSpeed});
    boundary.surfacePressure.resize(archive.mesh->triangles.size());
    output.observerPressure.resize(archive.request->observers.size());
    const auto area = [&](std::uint32_t id) {
        const auto &t = archive.mesh->triangles[id - 1];
        const auto a = archive.mesh->vertices[t.vertices[0]],
                   b = archive.mesh->vertices[t.vertices[1]],
                   c = archive.mesh->vertices[t.vertices[2]];
        return length(cross(subtract(b, a), subtract(c, a))) / 2;
    };
    for (std::size_t i = 0; i < ports.size(); ++i) {
        const auto &port = std::get<lem::BemPort>(ports[i]->component);
        const auto weight = [&](std::uint32_t id) {
            return port.velocityWeights.contains(id) ? port.velocityWeights.at(id) : 1.;
        };
        double weightedArea = 0;
        for (const auto id : port.faceIds)
            weightedArea += weight(id) * area(id);
        for (const auto id : port.faceIds) {
            const auto velocity =
                weight(id) * output.lemSolution->ports[i].volumeVelocity / weightedArea;
            auto pressure = output.lemSolution->ports[i].pressure;
            if (port.kind == lem::BemPortKind::drivenSurface)
                pressure += id == 1 ? lem::Complex{.1, -.07} : lem::Complex{-.2, .14};
            boundary.surfacePressure[id - 1] = pressure;
            output.solvedVelocities.push_back({id, velocity});
            if (velocity != lem::Complex{}) {
                const auto volumeVelocity = area(id) * velocity;
                output.impedances.push_back(
                    {id, area(id), pressure, volumeVelocity, pressure / volumeVelocity});
            }
        }
    }
    output.setBoundary(std::move(boundary));
    frequency.excitations.push_back(std::move(output));
    result->frequencies.push_back(std::move(frequency));
    archive.result = result;
    return project;
}

void weightedPortRoundTrip(const QTemporaryDir &directory) {
    const auto path = directory.filePath("weighted-port.hsproj");
    for (const bool baffled : {false, true}) {
        const auto project = weightedPortTestProject(baffled);
        require(project.simulation.request->excitations.empty(),
                "Port-driven diaphragm faces were also prescribed independently");
        writeProject(path, project);
        sameProject(readProject(path), project);
        const auto &values =
            project.simulation.result->frequencies[0].excitations[0].solvedVelocities;
        const auto velocity = [&](std::uint32_t id) {
            return std::find_if(values.begin(), values.end(),
                                [id](const auto &value) { return value.faceId == id; })
                ->normalVelocity;
        };
        sameComplex(velocity(2), .5 * velocity(1),
                    "Stored modal drive did not preserve the surround weight");
    }
    auto project = weightedPortTestProject(true);
    const auto archived = *project.simulation.request->coupling;
    auto &port = std::get<lem::BemPort>(
        project.simulation.settings.couplingNetwork.elements.back().component);
    port.velocityWeights[2] = .7;
    require(!lem::samePhysics(archived, project.simulation.settings.couplingNetwork),
            "Velocity weights are absent from physical snapshot equality");
    expectRejected([&] { writeProject(path, project); },
                   "A changed drive weight kept the archived result current");
    project.simulation.current = false;
    writeProject(path, project);
    sameProject(readProject(path), project);
    require(readProject(path).simulation.request->coupling == archived,
            "An edited weight overwrote the immutable solved snapshot");
    project = weightedPortTestProject(true);
    auto scaled = project;
    const auto rescaleWeights = [](lem::Network &network) {
        std::get<lem::BemPort>(network.elements.back().component).velocityWeights = {{1, 1e300},
                                                                                     {2, .5e300}};
    };
    rescaleWeights(scaled.simulation.settings.couplingNetwork);
    rescaleWeights(*scaled.simulation.request->coupling);
    auto largeResult = std::make_shared<BemResult>(*scaled.simulation.result);
    // A linear field rescaling exposes intermediate products w*p and w*Q
    // overflowing, while all physical results and normalized modal bases stay finite.
    constexpr double fieldScale = 1e40;
    for (auto &frequency : largeResult->frequencies)
        for (auto &excitation : frequency.excitations) {
            auto boundary = *excitation.boundary();
            for (auto &value : boundary.surfacePressure)
                value *= fieldScale;
            excitation.setBoundary(std::move(boundary));
            for (auto &value : excitation.observerPressure)
                value *= fieldScale;
            for (auto &value : excitation.solvedVelocities)
                value.normalVelocity *= fieldScale;
            for (auto &value : excitation.impedances) {
                value.meanPressure *= fieldScale;
                value.volumeVelocity *= fieldScale;
            }
            for (auto &value : excitation.lemSolution->nodes)
                value.effort *= fieldScale;
            for (auto &value : excitation.lemSolution->elements) {
                value.effort *= fieldScale;
                value.flow *= fieldScale;
                value.activePower = std::real(value.effort * std::conj(value.flow));
            }
            for (auto &value : excitation.lemSolution->ports) {
                value.pressure *= fieldScale;
                value.volumeVelocity *= fieldScale;
                value.activePower = std::real(value.pressure * std::conj(value.volumeVelocity));
            }
        }
    scaled.simulation.result = largeResult;
    writeProject(path, scaled);
    sameProject(readProject(path), scaled);
    writeProject(path, project);
    const auto original = QJsonDocument::fromJson(bytes(path)).object();
    const auto bad = directory.filePath("bad-weighted-port.hsproj");
    const auto reject = [&](QStringList keys, const QJsonValue &value) {
        saveBytes(bad, QJsonDocument(replaceAt(original, keys, 0, value).toObject()).toJson());
        expectRejected([&] { (void)readProject(bad); },
                       "An invalid weighted-port archive was accepted");
    };
    for (const auto section : {"settings", "request"}) {
        const QStringList prefix{
            "simulation",
            section,
            section == QStringLiteral("settings") ? "couplingNetwork" : "coupling",
            "elements",
            "2",
            "component",
            "velocityWeights"};
        for (const auto &value :
             {QJsonValue(QJsonArray{QJsonArray{2, 0}}), QJsonValue(QJsonArray{QJsonArray{2, -1}}),
              QJsonValue(QJsonArray{QJsonArray{2, QJsonValue::Null}}),
              QJsonValue(QJsonArray{QJsonArray{2, "0.5"}}),
              QJsonValue(QJsonArray{QJsonArray{3, .5}}),
              QJsonValue(QJsonArray{QJsonArray{2, .5}, QJsonArray{2, 1}}),
              QJsonValue(QJsonObject{{"2", .5}})})
            reject(prefix, value);
        reject({"simulation", section,
                section == QStringLiteral("settings") ? "couplingNetwork" : "coupling", "elements",
                "0", "component", "velocityWeights"},
               QJsonArray{QJsonArray{3, .5}});
    }
    reject({"simulation", "request", "excitations"},
           QJsonArray{QJsonObject{
               {"name", "Duplicate drive"},
               {"faces",
                QJsonArray{QJsonObject{{"faceId", 1}, {"normalVelocity", QJsonArray{1, 0}}}}}}});
    reject({"simulation", "result", "frequencies", "0", "excitations", "0", "solvedVelocities", "2",
            "normalVelocity"},
           QJsonArray{1, 0});

    const auto before = testCad(testMesh());
    auto after = *before;
    ++after.revision;
    std::erase_if(after.faces, [](const auto &face) { return face.id == 2; });
    const auto mapped = reconciledCadBemAssignments(archived, before.get(), &after);
    const auto &retained = std::get<lem::BemPort>(mapped.elements.back().component);
    require(retained.faceIds == std::vector<std::uint32_t>{1} && retained.velocityWeights.empty() &&
                retained.needsAssignmentReview,
            "Lost CAD faces retained weights or silently cleared their review status");
}

void acousticSlitRoundTrip(const QTemporaryDir &directory) {
    const auto path = directory.filePath("acoustic-slit.hsproj");
    const auto fixture = weightedPortTestProject(true, true);
    writeProject(path, fixture);
    sameProject(readProject(path), fixture);
    const auto &solution = *fixture.simulation.result->frequencies[0].excitations[0].lemSolution;
    const auto dissipative =
        std::find_if(solution.elements.begin(), solution.elements.end(),
                     [](const auto &element) { return element.id.value == 4; });
    require(dissipative != solution.elements.end() && dissipative->activePower > 0,
            "The slit persistence fixture does not exercise an actual viscous CPU LEM result");
    const auto original = QJsonDocument::fromJson(bytes(path)).object();
    const auto bad = directory.filePath("bad-acoustic-slit.hsproj");
    const auto reject = [&](QStringList keys, const QJsonValue &value) {
        saveBytes(bad, QJsonDocument(replaceAt(original, keys, 0, value).toObject()).toJson());
        expectRejected([&] { (void)readProject(bad); },
                       "Malformed typed acoustic slit was accepted");
    };
    for (const auto section : {"settings", "request"}) {
        const QStringList prefix{
            "simulation",
            section,
            section == QStringLiteral("settings") ? "couplingNetwork" : "coupling",
            "elements",
            "3",
            "component"};
        for (const auto *parameter : {"length", "area", "gap", "dynamicViscosity"}) {
            for (const auto &value : {QJsonValue(QJsonValue::Null), QJsonValue("0.03"),
                                      QJsonValue(true), QJsonValue(QJsonArray{})})
                reject(prefix + QStringList{parameter}, value);
            auto component =
                original["simulation"]
                    .toObject()[section]
                    .toObject()[section == QStringLiteral("settings") ? "couplingNetwork"
                                                                      : "coupling"]
                    .toObject()["elements"]
                    .toArray()[3]
                    .toObject()["component"]
                    .toObject();
            component.remove(parameter);
            reject(prefix, component);
            if (section == QStringLiteral("request"))
                reject(prefix + QStringList{parameter}, 0);
        }
        reject(prefix + QStringList{"terminals"}, QJsonArray{1, 999});
    }
    reject({"simulation", "request", "coupling", "nodes", "2", "domain"}, "mechanical");
    auto edited = fixture;
    std::get<lem::AcousticSlit>(edited.simulation.settings.couplingNetwork.elements[3].component)
        .gap = .002;
    expectRejected([&] { writeProject(path, edited); },
                   "Changing a real channel kept the solved archive current");
    edited.simulation.current = false;
    writeProject(path, edited);
    sameProject(readProject(path), edited);

    Project pending;
    pending.simulation.settings.couplingEnabled = true;
    pending.simulation.settings.couplingNetwork.nodes = {
        {{1}, lem::Domain::acoustic, false, "Front", {}},
        {{2}, lem::Domain::acoustic, true, "Reference", {}}};
    pending.simulation.settings.couplingNetwork.elements = {
        {{1}, "Unspecified real channel", {}, lem::AcousticSlit{{{1}, {2}}}}};
    writeProject(path, pending);
    sameProject(readProject(path), pending);
    expectRejected([&] { (void)makeSimulationRequest(pending.simulation.settings, {}); },
                   "An unfinished slit acquired usable implicit geometry");
}

void closedSlitCavityRoundTrip(const QTemporaryDir &directory) {
    const auto path = directory.filePath("closed-slit-cavity.hsproj");
    const auto fixture = weightedPortTestProject(true, true, true);
    writeProject(path, fixture);
    sameProject(readProject(path), fixture);
    const auto &solution = *fixture.simulation.result->frequencies[0].excitations[0].lemSolution;
    const auto chamber = std::find_if(solution.elements.begin(), solution.elements.end(),
                                      [](const auto &element) { return element.id.value == 2; });
    require(chamber != solution.elements.end() && chamber->activePower > 0,
            "Closed 1D slit fixture has no actual dissipative CPU LEM solution");
    const auto original = QJsonDocument::fromJson(bytes(path)).object();
    const auto bad = directory.filePath("bad-closed-slit.hsproj");
    const auto reject = [&](QStringList keys, const QJsonValue &value) {
        saveBytes(bad, QJsonDocument(replaceAt(original, keys, 0, value).toObject()).toJson());
        expectRejected([&] { (void)readProject(bad); },
                       "Malformed closed 1D slit archive was accepted");
    };
    for (const auto section : {"settings", "request"}) {
        const auto networkKey =
            section == QStringLiteral("settings") ? "couplingNetwork" : "coupling";
        const QStringList prefix{"simulation", section, networkKey, "elements", "1", "component"};
        const auto component = original["simulation"]
                                   .toObject()[section]
                                   .toObject()[networkKey]
                                   .toObject()["elements"]
                                   .toArray()[1]
                                   .toObject()["component"]
                                   .toObject();
        for (const auto *parameter :
             {"length", "area", "gap", "dynamicViscosity", "heatCapacityRatio", "prandtlNumber"}) {
            for (const auto &value : {QJsonValue(QJsonValue::Null), QJsonValue("1.4"),
                                      QJsonValue(true), QJsonValue(QJsonArray{})})
                reject(prefix + QStringList{parameter}, value);
            auto missing = component;
            missing.remove(parameter);
            reject(prefix, missing);
            if (section == QStringLiteral("request"))
                reject(prefix + QStringList{parameter}, 0);
        }
        auto unknown = component;
        unknown.insert("volume", 5e-6);
        reject(prefix, unknown);
    }
    reject({"simulation", "request", "coupling", "elements", "1", "component", "heatCapacityRatio"},
           .999);
    reject({"simulation", "request", "coupling", "nodes", "2", "domain"}, "electrical");
    auto edited = fixture;
    std::get<lem::ClosedSlitCavity>(
        edited.simulation.settings.couplingNetwork.elements[1].component)
        .length = .035;
    expectRejected([&] { writeProject(path, edited); },
                   "Changing the wave-carrying chamber kept old results current");
    edited.simulation.current = false;
    writeProject(path, edited);
    sameProject(readProject(path), edited);
    Project pending;
    pending.simulation.settings.couplingEnabled = true;
    pending.simulation.settings.couplingNetwork.nodes = {
        {{1}, lem::Domain::acoustic, false, "Input", {}},
        {{2}, lem::Domain::acoustic, true, "Reference", {}}};
    pending.simulation.settings.couplingNetwork.elements = {
        {{1}, "Unset 1D chamber", {}, lem::ClosedSlitCavity{{{1}, {2}}}}};
    writeProject(path, pending);
    sameProject(readProject(path), pending);
    expectRejected([&] { (void)makeSimulationRequest(pending.simulation.settings, {}); },
                   "Unset 1D chamber dimensions acquired implicit values");
}

void orientationRoundTrip(const QTemporaryDir &directory) {
    auto project = fullTestProject(ProjectSource::freecad);
    ++project.revision;
    project.geometryRevision = project.revision;
    project.simulation.current = false;
    project.draftReversedFaceIds = {1, 5};
    project.cadReversedFaceIds = {2};
    project.preview = std::make_shared<SurfaceMesh>(
        applyFaceOrientation(project.cad->mesh, project.cadReversedFaceIds));
    const auto path = directory.filePath("orientation-conflict.json");
    writeProject(path, project);
    sameProject(readProject(path), project);
    require(!inspectMesh(*readProject(path).preview).watertight(),
            "Orientation conflict must remain visible and reversible after reopening.");
    auto contradictory = project;
    contradictory.meshRevision = contradictory.revision;
    expectRejected([&] { writeProject(path, contradictory); },
                   "A conflicting orientation was accepted as a current BEM mesh.");
    sameProject(readProject(path), project);
    project.cadReversedFaceIds.clear();
    project.preview = std::make_shared<SurfaceMesh>(project.cad->mesh);
    project.meshRevision = project.revision;
    writeProject(path, project);
    sameProject(readProject(path), project);
    require(inspectMesh(*readProject(path).preview).watertight(),
            "Reverting orientation did not restore the valid source.");
}

void exactCadRoundTrip(const QTemporaryDir &directory) {
    Project project;
    project.source = ProjectSource::freecad;
    project.revision = project.geometryRevision = 7;
    auto cad = std::make_shared<CadSnapshot>(testing::exactCadFixture());
    cad->faces.front().boundary = BoundaryCondition::driven;
    cad->faces.front().normalVelocity = {.2, -.3};
    for (auto &triangle : cad->mesh.triangles)
        if (triangle.faceId == cad->faces.front().id)
            triangle.boundary = BoundaryCondition::driven;
    project.cad = cad;
    project.preview = std::make_shared<SurfaceMesh>(cad->mesh);
    project.cadFaceTargetEdgeLengths = {{cad->faces.front().id, .002}};
    const auto path = directory.filePath("exact-offline-cad.json");
    writeProject(path, project);
    const auto restored = readProject(path);
    sameProject(restored, project);
    auto regenerated = *restored.cad;
    prepareCadGeometry(regenerated, true);
    require(regenerated.geometry && inspectMesh(regenerated.mesh).watertight(),
            "Offline exact source cannot regenerate a native preview without FreeCAD");
    const auto original = QJsonDocument::fromJson(bytes(path)).object();
    const auto damaged = directory.filePath("bad-exact-cad.json");
    for (const auto &replacement :
         {replaceAt(original, {"cad", "brepObjects", "0", "data"}, 0, "@@@"),
          replaceAt(original, {"cad", "brepObjects", "0", "sha256"}, 0, QString(64, '0')),
          replaceAt(original, {"cad", "brepObjects", "0", "faceIds"}, 0, QJsonArray{1, 1, 3, 4}),
          replaceAt(original, {"cad", "faces", "0", "signature", "area"}, 0, .1),
          replaceAt(original, {"cad", "faces", "0", "signature", "center"}, 0, QJsonArray{1, 2, 3}),
          replaceAt(original, {"cad", "brepMetresPerUnit"}, 0, 1)}) {
        saveBytes(damaged, QJsonDocument(replacement.toObject()).toJson());
        expectRejected([&] { (void)readProject(damaged); },
                       "Corrupt exact CAD geometry or face binding was accepted");
    }
    auto legacy = preWidebandArchive(original, 6);
    legacy.remove("draftReversedFaces");
    legacy.remove("cadReversedFaces");
    legacy.insert("version", 6);
    legacy.remove("slot");
    auto offlineCad = legacy["cad"].toObject();
    offlineCad.remove("brepObjects");
    offlineCad.remove("brepMetresPerUnit");
    legacy.insert("cad", offlineCad);
    saveBytes(damaged, QJsonDocument(legacy).toJson());
    const auto legacyProject = readProject(damaged);
    require(legacyProject.cad && legacyProject.cad->brepObjects.empty() &&
                !legacyProject.cad->geometry,
            "Legacy mesh-only CAD archive must not claim exact CAD geometry");
    sameMesh(legacyProject.cad->mesh, cad->mesh);
}

void farFieldRoundTrip(const QTemporaryDir &directory) {
    auto project = fullTestProject(ProjectSource::draft);
    auto &settings = project.simulation.settings;
    settings.normalizationAngleDegrees = 23.5;
    settings.directivityMinimumDb = -42;
    settings.directivityMaximumDb = 3;
    settings.directivityStepDb = 1.5;
    project.simulation.request->farField = FarFieldObservation{{.1, .2, .3}, 2.5};
    const auto path = directory.filePath("far-field.json");
    writeProject(path, project);
    sameProject(readProject(path), project);
    const auto json = QJsonDocument::fromJson(bytes(path)).object();
    const auto damaged = directory.filePath("bad-far-field.json");
    for (const auto &value :
         {replaceAt(json, {"simulation", "request", "farField", "referenceDistance"}, 0, 0),
          replaceAt(json, {"simulation", "request", "farField", "origin"}, 0,
                    QJsonArray{.3, 0, .5}),
          replaceAt(json, {"simulation", "request", "farField", "origin"}, 0, QJsonArray{0, 0}),
          replaceAt(json, {"simulation", "request", "farField"}, 0, QJsonValue::Null)}) {
        saveBytes(damaged, QJsonDocument(value.toObject()).toJson());
        expectRejected([&] { (void)readProject(damaged); },
                       "Invalid far-field request was accepted");
    }
    auto legacy = preWidebandArchive(json, 7);
    legacy.remove("draftReversedFaces");
    legacy.remove("cadReversedFaces");
    legacy.insert("version", 7);
    legacy.remove("slot");
    auto simulation = legacy["simulation"].toObject();
    auto oldSettings = simulation["settings"].toObject();
    oldSettings.remove("normAngle");
    oldSettings.remove("DbScale");
    simulation.insert("settings", oldSettings);
    auto request = simulation["request"].toObject();
    request.remove("farField");
    simulation.insert("request", request);
    legacy.insert("simulation", simulation);
    saveBytes(damaged, QJsonDocument(legacy).toJson());
    auto expected = project;
    expected.simulation.settings.normalizationAngleDegrees = 0;
    expected.simulation.settings.directivityMinimumDb = -30;
    expected.simulation.settings.directivityMaximumDb = 3;
    expected.simulation.settings.directivityStepDb = 3;
    expected.simulation.request->farField.reset();
    sameProject(readProject(damaged), expected);
    settings.directivityMinimumDb = 0;
    settings.directivityMaximumDb = 999;
    settings.directivityStepDb = 1;
    writeProject(path, project);
    sameProject(readProject(path), project);
}

void faceColorsRoundTrip(const QTemporaryDir &directory) {
    const FaceColor orange{255, 145, 0};
    const auto path = directory.filePath("face-colors.hsproj");
    for (const auto source : {ProjectSource::draft, ProjectSource::freecad}) {
        auto project = fullTestProject(source);
        project.draftFaceColors = {{1, orange}, {5, {0, 0, 0}}};
        project.cadFaceColors = {{2, {0, 255, 255}}, {3, orange}};
        writeProject(path, project);
        const auto restored = readProject(path);
        sameProject(restored, project);
        require(restored.simulation.current && restored.revision == project.revision,
                "Display colors invalidated current acoustic data");
        const auto json = QJsonDocument::fromJson(bytes(path)).object();
        auto old = json;
        old.remove("draftFaceColors");
        old.remove("cadFaceColors");
        saveBytes(path, QJsonDocument(old).toJson());
        auto expected = project;
        expected.draftFaceColors.clear();
        expected.cadFaceColors.clear();
        sameProject(readProject(path), expected);

        const auto reject = [&](const char *field, QJsonValue value) {
            saveBytes(path, QJsonDocument(replaceAt(json, {field}, 0, value).toObject()).toJson());
            expectRejected([&] { (void)readProject(path); },
                           "Malformed display face colors were accepted");
        };
        for (const auto field : {"draftFaceColors", "cadFaceColors"}) {
            for (const auto &value :
                 {QJsonValue(QJsonValue::Null), QJsonValue(QJsonObject{}), QJsonValue(true),
                  QJsonValue(QJsonArray{QJsonArray{0, QJsonArray{1, 2, 3}}}),
                  QJsonValue(QJsonArray{QJsonArray{1, QJsonArray{1, 2, 3}},
                                        QJsonArray{1, QJsonArray{4, 5, 6}}}),
                  QJsonValue(QJsonArray{QJsonArray{1.5, QJsonArray{1, 2, 3}}}),
                  QJsonValue(QJsonArray{QJsonArray{1, QJsonArray{1, 2}}}),
                  QJsonValue(QJsonArray{QJsonArray{1, QJsonArray{1, 2, 3, 255}}}),
                  QJsonValue(QJsonArray{QJsonArray{1, QJsonArray{256, 2, 3}}}),
                  QJsonValue(QJsonArray{QJsonArray{1, QJsonArray{-1, 2, 3}}}),
                  QJsonValue(QJsonArray{QJsonArray{1, QJsonArray{1.5, 2, 3}}}),
                  QJsonValue(QJsonArray{QJsonArray{1, QJsonArray{true, 2, 3}}}),
                  QJsonValue(QJsonArray{QJsonArray{1, QJsonArray{"1", 2, 3}}}),
                  QJsonValue(QJsonArray{QJsonArray{1, QJsonValue::Null}}),
                  QJsonValue(QJsonArray{QJsonArray{1, QJsonArray{1, 2, 3}, 4}})})
                reject(field, value);
        }
        reject("draftFaceColors", QJsonArray{QJsonArray{9, QJsonArray{255, 145, 0}}});
        reject("cadFaceColors", QJsonArray{QJsonArray{999, QJsonArray{255, 145, 0}}});
        auto invalid = project;
        invalid.draftFaceColors[9] = orange;
        expectRejected([&] { writeProject(path, invalid); }, "Unknown draft color ID was saved");
        invalid = project;
        invalid.cadFaceColors[999] = orange;
        expectRejected([&] { writeProject(path, invalid); }, "Unknown CAD color ID was saved");
    }
    Project plain;
    plain.cadFaceColors = {{1, orange}};
    expectRejected([&] { writeProject(path, plain); },
                   "CAD colors without a CAD identity source were saved");

    const auto before = testCad(testMesh());
    auto after = *before;
    ++after.revision;
    std::reverse(after.faces.begin(), after.faces.end());
    for (auto &face : after.faces)
        face.sourceKey = "Face999";
    const FaceColors assigned{{1, {0, 255, 0}}, {3, orange}, {999, {1, 2, 3}}};
    const FaceColors expected{{1, {0, 255, 0}}, {3, orange}};
    require(retainedCadFaceColors(assigned, before.get(), &after) == expected,
            "Mapped CAD colors followed face order/names or retained an unknown ID");
    auto changed = after;
    for (auto &face : changed.faces)
        if (face.id == 1)
            face.objectId = "different-object";
        else if (face.id == 3)
            face.state = CadFaceState::needsReview;
    require(retainedCadFaceColors(assigned, before.get(), &changed).empty(),
            "Ambiguous or object-reused CAD IDs inherited colors");
    changed = after;
    std::erase_if(changed.faces, [](const auto &face) { return face.id == 1; });
    for (auto &face : changed.faces)
        if (face.id == 3)
            face.state = CadFaceState::newFace;
    require(retainedCadFaceColors(assigned, before.get(), &changed).empty(),
            "Removed or new topology silently reused old colors");
    changed = after;
    changed.documentId = "different-document";
    require(retainedCadFaceColors(assigned, before.get(), &changed).empty() &&
                retainedCadFaceColors(assigned, nullptr, &after).empty() &&
                retainedCadFaceColors(assigned, before.get(), nullptr).empty(),
            "Colors crossed a missing or different document identity");
    require(assigned.at(3) == orange, "CAD reconciliation mutated its input colors");
}

void angularResponsesRoundTrip(const QTemporaryDir &directory) {
    const AngularResponseSettings defaults;
    require(defaults.anglesDegrees == std::vector<double>({0, 15, 30, 45, 60}) &&
                !defaults.normalizeToAngle && defaults.referenceAngleDegrees == 0,
            "Angular response defaults changed");
    validateAngularResponseSettings(defaults);
    AngularResponseSettings maximum;
    maximum.anglesDegrees.clear();
    for (unsigned i = 0; i < 12; ++i)
        maximum.anglesDegrees.push_back(-180 + 360. * i / 11);
    maximum.referenceAngleDegrees = 180;
    validateAngularResponseSettings(maximum);
    require(maximum != defaults, "Angular response equality ignored the selected angles");

    auto project = fullTestProject(ProjectSource::draft);
    const std::vector<FaceVelocity> velocity{{1, {1, 0}}};
    const auto request = makeSimulationRequest(project.simulation.settings, velocity);
    project.simulation.settings.angularResponses = {{-180, -45.5, 0, 17.25, 180}, true, -12.5};
    project.simulation.settings.fieldDisplay = {true, -.04, .04, -.035, .1, 65,
                                                81,   1,    30,  115,   25, .00005};
    project.simulation.settings.fieldDisplay.frequencyHz = 1000;
    project.simulation.settings.fieldDisplay.excitationIndex = 1;
    project.simulation.settings.ctaExcitationIndex = 1;
    project.simulation.settings.phaseDisplay = {30, -15, .0001, true, true};
    sameRequest(makeSimulationRequest(project.simulation.settings, velocity), request);
    const auto path = directory.filePath("angular-responses.json");
    writeProject(path, project);
    const auto saved = bytes(path);
    const auto restored = readProject(path);
    sameProject(restored, project);
    require(restored.simulation.current && restored.simulation.result,
            "Display-only angular responses invalidated a current result");

    const auto json = QJsonDocument::fromJson(saved).object();
    const auto simulation = json["simulation"].toObject();
    const auto settings = simulation["settings"].toObject();
    require(settings["angularResponses"].isObject() &&
                settings["angularResponses"].toObject().size() == 3 &&
                !simulation["request"].toObject().contains("angularResponses"),
            "Angular responses were omitted or added to the acoustic request");

    auto oldFieldSettings = settings;
    auto oldField = oldFieldSettings["fieldDisplay"].toObject();
    oldField.remove("frequencyHz");
    oldField.remove("excitationIndex");
    oldFieldSettings.insert("fieldDisplay", oldField);
    oldFieldSettings.remove("ctaExcitationIndex");
    auto oldFieldSimulation = simulation;
    oldFieldSimulation.insert("settings", oldFieldSettings);
    auto oldFieldProject = json;
    oldFieldProject.insert("simulation", oldFieldSimulation);
    const auto oldFieldPath = directory.filePath("old-field-selection.json");
    saveBytes(oldFieldPath, QJsonDocument(oldFieldProject).toJson());
    const auto oldFieldLoaded = readProject(oldFieldPath);
    require(oldFieldLoaded.simulation.settings.fieldDisplay.frequencyHz == 0 &&
                oldFieldLoaded.simulation.settings.fieldDisplay.excitationIndex == -1 &&
                oldFieldLoaded.simulation.settings.ctaExcitationIndex == -1 &&
                oldFieldLoaded.simulation.current,
            "Old project selection fallback damaged the result");
    auto oldSettings = settings;
    oldSettings.remove("angularResponses");
    auto oldSimulation = simulation;
    oldSimulation.insert("settings", oldSettings);
    auto old = json;
    old.insert("simulation", oldSimulation);
    const auto damaged = directory.filePath("bad-angular-responses.json");
    saveBytes(damaged, QJsonDocument(old).toJson());
    auto expected = project;
    expected.simulation.settings.angularResponses = defaults;
    sameProject(readProject(damaged), expected);

    Project plain;
    const auto defaultsPath = directory.filePath("angular-response-defaults.json");
    writeProject(defaultsPath, plain);
    const auto defaultsJson = QJsonDocument::fromJson(bytes(defaultsPath)).object();
    require(
        defaultsJson["simulation"].toObject()["settings"].toObject()["angularResponses"].isObject(),
        "Default angular response settings were not written");
    sameProject(readProject(defaultsPath), plain);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();
    std::vector<AngularResponseSettings> invalid{
        {{}, false, 0},         {{0, 0}, false, 0},      {{0, -0.}, false, 0},
        {{nan}, false, 0},      {{infinity}, false, 0},  {{-infinity}, false, 0},
        {{-180.01}, false, 0},  {{180.01}, false, 0},    {{0}, false, nan},
        {{0}, false, infinity}, {{0}, false, -infinity}, {{0}, false, -180.01},
        {{0}, false, 180.01}};
    maximum.anglesDegrees.push_back(1);
    invalid.push_back(maximum);
    for (const auto &bad : invalid) {
        expectRejected([&] { validateAngularResponseSettings(bad); },
                       "Invalid angular response settings passed common validation");
        auto broken = project;
        broken.simulation.settings.angularResponses = bad;
        expectRejected([&] { writeProject(path, broken); },
                       "Invalid angular response settings were written");
        require(bytes(path) == saved, "An invalid display setting overwrote a valid project");
    }

    const QStringList prefix{"simulation", "settings", "angularResponses"};
    const auto reject = [&](QStringList keys, QJsonValue value) {
        saveBytes(damaged, QJsonDocument(replaceAt(json, keys, 0, value).toObject()).toJson());
        expectRejected([&] { (void)readProject(damaged); },
                       "Malformed angular response JSON was accepted");
    };
    for (const auto &value : {QJsonValue(QJsonValue::Null), QJsonValue(true), QJsonValue(0),
                              QJsonValue("angles"), QJsonValue(QJsonArray{})})
        reject(prefix, value);
    for (const auto &value :
         {QJsonValue(QJsonArray{}), QJsonValue(QJsonArray{0, 0}), QJsonValue(QJsonArray{0, -0.}),
          QJsonValue(QJsonArray{181}), QJsonValue(QJsonArray{-181}),
          QJsonValue(QJsonArray{QJsonValue::Null}), QJsonValue(QJsonArray{"15"}),
          QJsonValue(QJsonArray{true}), QJsonValue(QJsonArray{QJsonObject{}}), QJsonValue("0,15"),
          QJsonValue(QJsonValue::Null)})
        reject(prefix + QStringList{"anglesDegrees"}, value);
    QJsonArray tooMany;
    for (unsigned i = 0; i < 13; ++i)
        tooMany.append(int(i));
    reject(prefix + QStringList{"anglesDegrees"}, tooMany);
    for (const auto &value : {QJsonValue("true"), QJsonValue(1), QJsonValue(QJsonValue::Null)})
        reject(prefix + QStringList{"normalizeToAngle"}, value);
    for (const auto &value : {QJsonValue(-181), QJsonValue(181), QJsonValue("0"), QJsonValue(false),
                              QJsonValue(QJsonValue::Null)})
        reject(prefix + QStringList{"referenceAngleDegrees"}, value);
    const auto encoded = settings["angularResponses"].toObject();
    for (const auto key : {"anglesDegrees", "normalizeToAngle", "referenceAngleDegrees"}) {
        auto missing = encoded;
        missing.remove(key);
        reject(prefix, missing);
    }
    auto unknown = encoded;
    unknown.insert("automatic", true);
    reject(prefix, unknown);
}
void savedFieldTracesRoundTrip(const QTemporaryDir &directory) {
    auto project = fullTestProject(ProjectSource::draft);
    auto result = std::make_shared<BemResult>(*project.simulation.result);
    for (auto &f : result->frequencies)
        for (auto &e : f.excitations) {
            // Persistence fixture only: a constant P1 field makes its exact
            // triangle means explicit without pretending this is a BEM solve.
            const std::complex<double> pressure{.75, -.125};
            auto boundary = *e.boundary();
            boundary.pressureTrace.assign(project.simulation.mesh->vertices.size(), pressure);
            boundary.normalDerivative.assign(project.simulation.mesh->triangles.size(), {.4, -.2});
            boundary.surfacePressure.assign(project.simulation.mesh->triangles.size(), pressure);
            e.setBoundary(std::move(boundary));
            for (auto &z : e.impedances) {
                z.meanPressure = pressure;
                z.volumeImpedance = pressure / z.volumeVelocity;
            }
        }
    // Copies share immutable boundary data but not observer arrays. Replacing
    // one copy's boundary must not modify an already published snapshot.
    auto copied = *result;
    auto &originalExcitation = result->frequencies[0].excitations[0];
    auto &copiedExcitation = copied.frequencies[0].excitations[0];
    require(copiedExcitation.boundary() == originalExcitation.boundary(),
            "Result snapshots did not share immutable boundary data.");
    require(result->frequencies[0].excitations[0].boundary() !=
                result->frequencies[1].excitations[0].boundary(),
            "Different frequencies shared one boundary block.");
    auto replacement = *copiedExcitation.boundary();
    replacement.pressureTrace[0] += 1.;
    copiedExcitation.setBoundary(std::move(replacement));
    require(copiedExcitation.boundary() != originalExcitation.boundary() &&
                copiedExcitation.boundary()->pressureTrace[0] !=
                    originalExcitation.boundary()->pressureTrace[0],
            "Replacing a copied boundary changed the original snapshot.");
    copiedExcitation.observerPressure[0] += 1.;
    require(copiedExcitation.observerPressure[0] != originalExcitation.observerPressure[0],
            "Snapshot observer data was not independent.");
    project.simulation.result = result;
    const auto path = directory.filePath("saved-field-traces.json");
    writeProject(path, project);
    const auto restored = readProject(path);
    for (std::size_t f = 0; f < result->frequencies.size(); ++f)
        for (std::size_t e = 0; e < result->frequencies[f].excitations.size(); ++e) {
            const auto &a = result->frequencies[f].excitations[e];
            const auto &b = restored.simulation.result->frequencies[f].excitations[e];
            require(a.boundary()->pressureTrace == b.boundary()->pressureTrace &&
                        a.boundary()->normalDerivative == b.boundary()->normalDerivative,
                    "Saved P1/DP0 Cauchy traces changed on project reload.");
        }
    auto broken = std::make_shared<BemResult>(*result);
    auto missingTrace = *broken->frequencies[0].excitations[0].boundary();
    missingTrace.normalDerivative.clear();
    broken->frequencies[0].excitations[0].setBoundary(std::move(missingTrace));
    project.simulation.result = broken;
    expectRejected([&] { writeProject(path, project); }, "A half-present Cauchy trace was saved.");
    broken = std::make_shared<BemResult>(*result);
    auto badTrace = *broken->frequencies[0].excitations[0].boundary();
    badTrace.pressureTrace[0] += 1.;
    broken->frequencies[0].excitations[0].setBoundary(std::move(badTrace));
    project.simulation.result = broken;
    expectRejected([&] { writeProject(path, project); },
                   "P1 trace incompatible with saved surface means was accepted.");
}

void axisymmetricRoundTrip(const QTemporaryDir &directory) {
    // Persistence-only data on an actual periodic integration mesh, not a solver reference.
    AxisymmetricMeridian meridian;
    meridian.points = {{0, -.005}, {.01, -.005}, {.015, 0}, {0, 0}};
    meridian.segments = {{0, 1, 1, BoundaryCondition::driven},
                         {1, 2, 2, BoundaryCondition::wall},
                         {2, 3, 3, BoundaryCondition::wall}};
    auto generated = buildAxisymmetricMesh(meridian, 16);
    Project project;
    project.revision = project.geometryRevision = largeRevision;
    project.preview = std::make_shared<SurfaceMesh>(std::move(generated.mesh));
    auto &a = project.simulation;
    a.revision = project.revision;
    a.current = true;
    a.mesh = project.preview;
    a.settings.minimumFrequency = a.settings.maximumFrequency = 1000;
    a.settings.frequencyCount = 1;
    a.settings.geometryMode = BemGeometryMode::axisymmetric;
    a.settings.axisymmetricSampling = {.003, .00007, 2e-8, 16};
    a.settings.infiniteBaffle = InfiniteBaffle{{}, {0, 0, 1}, {3}};
    a.axisymmetricSampling = a.settings.axisymmetricSampling;
    const std::vector<FaceVelocity> velocities{{1, {1, .25}}};
    a.request = makeSimulationRequest(a.settings, velocities);
    a.request->axisymmetric =
        AxisymmetricBem{meridian.origin,        meridian.axis,      16,
                        generated.pressureDofs, generated.fluxDofs, generated.testPanels};
    auto result = std::make_shared<BemResult>();
    result->device = "SYNTHETIC 2D PERSISTENCE FIXTURE - no GPU computation";
    result->geometryMode = BemGeometryMode::axisymmetric;
    result->triangleCount = a.mesh->triangles.size();
    result->pressureDofCount = meridian.points.size();
    BemExcitationResult excitation;
    BemBoundaryData boundary;
    excitation.name = a.request->excitations.front().name;
    const std::complex<double> pressure{.7, -.2};
    for (const auto dof : generated.pressureDofs)
        boundary.pressureTrace.push_back(pressure + std::complex<double>{.05 * dof, .03 * dof});
    boundary.normalDerivative.assign(result->triangleCount, {.4, .1});
    excitation.observerPressure.assign(a.request->observers.size(), {.01, -.002});
    double area = 0;
    const auto measures = surfaceMeasures(*a.mesh, a.request->axisymmetric);
    std::complex<double> integratedPressure{};
    bool hasNonCartesianMean = false;
    for (std::size_t ti = 0; ti < a.mesh->triangles.size(); ++ti) {
        std::complex<double> mean{}, cartesianMean{};
        const auto &triangle = a.mesh->triangles[ti];
        for (unsigned local = 0; local < 3; ++local) {
            const auto value = boundary.pressureTrace[triangle.vertices[local]];
            mean += (measures.pressureIntegrals[ti][local] / measures.areas[ti]) * value;
            cartesianMean += value / 3.;
        }
        boundary.surfacePressure.push_back(mean);
        hasNonCartesianMean |= std::abs(mean - cartesianMean) > 1e-5;
        if (triangle.faceId == 1) {
            area += measures.areas[ti];
            integratedPressure += measures.areas[ti] * mean;
        }
    }
    require(hasNonCartesianMean,
            "The continuous persistence fixture never exercises radial weights.");
    const auto flow = area * velocities.front().normalVelocity;
    const auto meanPressure = integratedPressure / area;
    excitation.impedances.push_back({1, area, meanPressure, flow, meanPressure / flow});
    excitation.setBoundary(std::move(boundary));
    result->frequencies.push_back({1000, {excitation}});
    a.result = result;
    const auto path = directory.filePath("axisymmetric.hsproj");
    writeProject(path, project);
    const auto restored = readProject(path);
    sameProject(restored, project);
    require(!restored.mesh && restored.meshRevision == 0 && restored.simulation.current &&
                restored.simulation.result->pressureDofCount == 4 &&
                restored.simulation.result->frequencies[0]
                        .excitations[0]
                        .boundary()
                        ->surfacePressure.size() == restored.simulation.mesh->triangles.size(),
            "2D results require a 3D quality mesh or lost expanded surface fields.");
    const auto continuousCsv = resultCsv(restored.simulation);
    require(continuousCsv.contains("surface_contract=continuousRevolution") &&
                continuousCsv.contains("carrier_triangles="),
            "Continuous ring CSV labels its carrier cells as integration triangles.");
    const auto json = QJsonDocument::fromJson(bytes(path)).object();
    const auto damaged = directory.filePath("bad-axisymmetric.hsproj");
    const auto reject = [&](QStringList keys, QJsonValue value) {
        saveBytes(damaged, QJsonDocument(replaceAt(json, keys, 0, value).toObject()).toJson());
        expectRejected([&] { (void)readProject(damaged); },
                       "Invalid axisymmetric snapshot was accepted.");
    };
    for (const auto prefix : {QStringList{"simulation", "settings", "axisymmetricSampling"},
                              QStringList{"simulation", "axisymmetricSampling"}}) {
        for (const auto field : {"maximumLength", "chordTolerance", "rotationalTolerance"}) {
            reject(prefix + QStringList{field}, 0);
            reject(prefix + QStringList{field}, -1);
            reject(prefix + QStringList{field}, QJsonValue::Null);
            reject(prefix + QStringList{field}, "0.002");
        }
        for (const auto sectors : {0, 7, 10, 4097})
            reject(prefix + QStringList{"sectors"}, sectors);
        reject(prefix + QStringList{"sectors"}, "16");
        reject(prefix, QJsonValue::Null);
    }
    const QStringList reduction{"simulation", "request", "axisymmetric"};
    reject(reduction + QStringList{"surfaceContract"}, "unknown");
    reject(reduction + QStringList{"surfaceContract"}, QJsonValue::Null);
    reject(reduction + QStringList{"surfaceContract"}, 1);
    auto historicalFixture = project;
    historicalFixture.simulation.request->axisymmetric->surfaceContract =
        AxisymmetricSurfaceContract::triangulatedOrbit;
    auto historicalResult = std::make_shared<BemResult>(*result);
    const auto historicalMeasures =
        surfaceMeasures(*a.mesh, historicalFixture.simulation.request->axisymmetric);
    double historicalArea = 0;
    std::complex<double> historicalIntegral{};
    auto &historicalExcitation = historicalResult->frequencies[0].excitations[0];
    auto historicalBoundary = *historicalExcitation.boundary();
    for (std::size_t ti = 0; ti < a.mesh->triangles.size(); ++ti) {
        std::complex<double> mean{};
        const auto &triangle = a.mesh->triangles[ti];
        for (unsigned local = 0; local < 3; ++local)
            mean +=
                (historicalMeasures.pressureIntegrals[ti][local] / historicalMeasures.areas[ti]) *
                historicalBoundary.pressureTrace[triangle.vertices[local]];
        historicalBoundary.surfacePressure[ti] = mean;
        if (triangle.faceId == 1) {
            historicalArea += historicalMeasures.areas[ti];
            historicalIntegral += historicalMeasures.areas[ti] * mean;
        }
    }
    historicalExcitation.setBoundary(std::move(historicalBoundary));
    auto &historicalImpedance = historicalResult->frequencies[0].excitations[0].impedances[0];
    historicalImpedance.area = historicalArea;
    historicalImpedance.meanPressure = historicalIntegral / historicalArea;
    historicalImpedance.volumeVelocity = historicalArea * velocities[0].normalVelocity;
    historicalImpedance.volumeImpedance =
        historicalImpedance.meanPressure / historicalImpedance.volumeVelocity;
    historicalFixture.simulation.result = historicalResult;
    const auto historicalPath = directory.filePath("historical-axisymmetric.hsproj");
    writeProject(historicalPath, historicalFixture);
    auto historicalJson = QJsonDocument::fromJson(bytes(historicalPath)).object();
    auto historicalSimulation = historicalJson["simulation"].toObject();
    auto historicalRequest = historicalSimulation["request"].toObject();
    auto historicalReduction = historicalRequest["axisymmetric"].toObject();
    historicalReduction.remove("surfaceContract");
    historicalRequest["axisymmetric"] = historicalReduction;
    historicalSimulation["request"] = historicalRequest;
    historicalJson["simulation"] = historicalSimulation;
    saveBytes(historicalPath, QJsonDocument(historicalJson).toJson());
    auto historical = readProject(historicalPath);
    sameProject(historical, historicalFixture);
    const auto historicalCsv = resultCsv(historical.simulation);
    require(historicalCsv.contains("surface_contract=triangulatedOrbit") &&
                historicalCsv.contains("integration_triangles="),
            "Historical polygonal CSV changed its surface interpretation.");
    require(
        historical.simulation.request->axisymmetric->surfaceContract ==
                AxisymmetricSurfaceContract::triangulatedOrbit &&
            historical.simulation.result->frequencies[0].excitations[0].boundary()->pressureTrace ==
                restored.simulation.result->frequencies[0].excitations[0].boundary()->pressureTrace,
        "Historical triangular traces were reinterpreted as a continuous revolution.");
    writeProject(historicalPath, historical);
    sameProject(readProject(historicalPath), historical);
    reject(reduction + QStringList{"pressureDofs"}, QJsonArray{});
    reject(reduction + QStringList{"fluxDofs"}, QJsonArray{});
    reject(reduction + QStringList{"testPanels"}, QJsonArray{0, 0});
    reject(reduction + QStringList{"pressureDofs", "0"}, -1);
    reject(reduction + QStringList{"pressureDofs", "0"}, 1);
    reject(reduction + QStringList{"fluxDofs", "0"}, 999);
    reject(reduction + QStringList{"axis"}, QJsonArray{0, 0, 2});
    reject(reduction + QStringList{"origin"}, QJsonArray{.001, 0, 0});
    reject(reduction + QStringList{"sectors"}, 8);
    reject({"simulation", "result", "geometryMode"}, "threeDimensional");
    reject({"simulation", "result", "pressureDofCount"}, double(a.mesh->vertices.size()));
    reject({"simulation", "settings", "geometryMode"}, "threeDimensional");
    reject({"simulation", "settings", "geometryMode"}, "2D");
    reject({"simulation", "settings", "axisymmetricSampling", "maximumLength"}, .004);
    reject({"simulation", "request", "numerics", "acceleration"}, "fastMultipole");
    reject({"simulation", "request", "infiniteBaffle"}, QJsonValue::Null);
    reject({"simulation", "solvedMesh", "vertices", "1", "0"}, .091);
    auto invalid = a.settings.axisymmetricSampling;
    invalid.rotationalTolerance = std::numeric_limits<double>::quiet_NaN();
    expectRejected([&] { validateAxisymmetricSamplingSettings(invalid); },
                   "Nonfinite 2D geometry tolerance was accepted.");
    auto unsupported = a.settings;
    unsupported.infiniteBaffle.reset();
    expectRejected([&] { (void)makeSimulationRequest(unsupported, velocities); },
                   "A reduced free-field request without supported air cavity was accepted.");
    unsupported = a.settings;
    unsupported.acceleration = BemAcceleration::fastMultipole;
    expectRejected([&] { (void)makeSimulationRequest(unsupported, velocities); },
                   "2D FMM silently selected a different solver.");

    a.current = false;
    a.settings.geometryMode = BemGeometryMode::threeDimensional;
    a.settings.axisymmetricSampling.maximumLength = .004;
    writeProject(path, project);
    sameProject(readProject(path), project);
    require(readProject(path).simulation.axisymmetricSampling->maximumLength == .003,
            "Editing 2D inputs rewrote the computed integration mesh snapshot.");

    auto legacy = fullTestProject(ProjectSource::draft);
    writeProject(path, legacy);
    auto old = QJsonDocument::fromJson(bytes(path)).object();
    auto simulation = old["simulation"].toObject();
    auto settings = simulation["settings"].toObject();
    settings.remove("geometryMode");
    settings.remove("axisymmetricSampling");
    simulation["settings"] = settings;
    auto oldResult = simulation["result"].toObject();
    oldResult.remove("geometryMode");
    simulation["result"] = oldResult;
    old["simulation"] = simulation;
    saveBytes(path, QJsonDocument(old).toJson());
    sameProject(readProject(path), legacy);
}

void interactionPlanningRoundTrip(const QTemporaryDir &directory) {
    auto project = fullTestProject(ProjectSource::draft);
    project.simulation.request->interactionPlanFrequencies = {315, 500, 1000};
    auto result = std::make_shared<BemResult>(*project.simulation.result);
    result->frequencies[0].interactionPlanMilliseconds = 7.25;
    result->frequencies[1].interactionPlanReused = true;
    project.simulation.result = result;
    const auto path = directory.filePath("interaction-planning.hsproj");
    writeProject(path, project);
    sameProject(readProject(path), project);
    const auto csv = resultCsv(project.simulation);
    require(csv.contains("cpu_build_ms=7.25") && csv.contains("host_cache_hit=true"),
            "Host planning time and cache state are absent from the result CSV.");
    const auto json = QJsonDocument::fromJson(bytes(path)).object();
    const auto damaged = directory.filePath("bad-interaction-planning.hsproj");
    const auto reject = [&](QStringList keys, QJsonValue value) {
        saveBytes(damaged, QJsonDocument(replaceAt(json, keys, 0, value).toObject()).toJson());
        expectRejected([&] { (void)readProject(damaged); },
                       "Invalid interaction planning archive accepted.");
    };
    for (const auto value :
         {QJsonValue(QJsonValue::Null), QJsonValue("315,1000"), QJsonValue(QJsonArray{}),
          QJsonValue(QJsonArray{315, 315, 1000}), QJsonValue(QJsonArray{315, 500}),
          QJsonValue(QJsonArray{315, 0, 1000}), QJsonValue(QJsonArray{315, -1, 1000}),
          QJsonValue(QJsonArray{315, "1000"})})
        reject({"simulation", "request", "interactionPlanFrequencies"}, value);
    const QStringList timing{"simulation", "result", "frequencies", "0"};
    for (const auto value :
         {QJsonValue(-1), QJsonValue(1000), QJsonValue("7.25"), QJsonValue(QJsonValue::Null)})
        reject(timing + QStringList{"interactionPlanMilliseconds"}, value);
    reject(timing + QStringList{"interactionPlanReused"}, true);
    reject(timing + QStringList{"interactionPlanReused"}, 1);

    project.simulation.current = false;
    project.simulation.request->frequencies.resize(1);
    auto partial = std::make_shared<BemResult>(*result);
    partial->frequencies.resize(1);
    project.simulation.result = partial;
    writeProject(path, project);
    sameProject(readProject(path), project);
    require(readProject(path).simulation.request->interactionPlanFrequencies.size() == 3,
            "A partial-frequency archive discarded the complete planning raster.");

    auto old = json;
    auto simulation = old["simulation"].toObject();
    auto request = simulation["request"].toObject();
    request.remove("interactionPlanFrequencies");
    simulation["request"] = request;
    auto output = simulation["result"].toObject();
    auto frequencies = output["frequencies"].toArray();
    for (qsizetype i = 0; i < frequencies.size(); ++i) {
        auto f = frequencies[i].toObject();
        f.remove("interactionPlanMilliseconds");
        f.remove("interactionPlanReused");
        frequencies[i] = f;
    }
    output["frequencies"] = frequencies;
    simulation["result"] = output;
    old["simulation"] = simulation;
    saveBytes(path, QJsonDocument(old).toJson());
    sameProject(readProject(path), fullTestProject(ProjectSource::draft));
}

void infiniteBaffleRoundTrip(const QTemporaryDir &directory) {
    auto project = fullTestProject(ProjectSource::draft);
    const InfiniteBaffle baffle{{.001, -.002, .03}, {0, 0, 1}, {3, 4}};
    project.simulation.settings.infiniteBaffle = baffle;
    project.simulation.request->infiniteBaffle = baffle;
    const auto path = directory.filePath("infinite-baffle.hsproj");
    writeProject(path, project);
    sameProject(readProject(path), project);
    const auto json = QJsonDocument::fromJson(bytes(path)).object();
    const auto damaged = directory.filePath("bad-infinite-baffle.hsproj");
    const auto reject = [&](QStringList keys, QJsonValue value) {
        saveBytes(damaged, QJsonDocument(replaceAt(json, keys, 0, value).toObject()).toJson());
        expectRejected([&] { (void)readProject(damaged); },
                       "Invalid infinite-baffle domain was accepted");
    };
    for (const auto section : {QStringLiteral("settings"), QStringLiteral("request")}) {
        const QStringList prefix{"simulation", section, "infiniteBaffle"};
        reject(prefix + QStringList{"normal"}, QJsonArray{0, 0, 0});
        reject(prefix + QStringList{"normal"}, QJsonArray{0, 0, 2});
        reject(prefix + QStringList{"origin"}, QJsonArray{0, 0});
        reject(prefix + QStringList{"origin", "0"}, QJsonValue::Null);
        reject(prefix + QStringList{"apertureFaceIds"}, QJsonArray{});
        reject(prefix + QStringList{"apertureFaceIds"}, QJsonArray{3, 3});
        reject(prefix + QStringList{"apertureFaceIds"}, QJsonArray{0});
        reject(prefix + QStringList{"apertureFaceIds"}, QJsonArray{3.5});
        reject(prefix, QJsonValue::Null);
        auto unknownKey =
            json["simulation"].toObject()[section].toObject()["infiniteBaffle"].toObject();
        unknownKey.insert("unexpected", true);
        reject(prefix, unknownKey);
    }
    reject({"simulation", "request", "infiniteBaffle", "apertureFaceIds"}, QJsonArray{1});
    reject({"simulation", "request", "infiniteBaffle", "apertureFaceIds"}, QJsonArray{999});
    reject({"simulation", "request", "numerics", "acceleration"}, "fastMultipole");

    auto &settings = project.simulation.settings;
    const std::vector<FaceVelocity> velocities{{1, {1, 0}}};
    require(makeSimulationRequest(settings, velocities).infiniteBaffle == baffle,
            "The common request builder lost the infinite-baffle domain");
    settings.observationFrame = AcousticFrame{{.002, -.001, .01}, {1, 0, 0}, {0, 1, 0}};
    const auto rotated = makeSimulationRequest(settings, velocities);
    require(rotated.observers.size() == 147 && rotated.infiniteBaffle == baffle &&
                rotated.observers.front().x == 1.002 &&
                rotated.farField->origin == settings.observationFrame->origin,
            "The baffle replaced the independently selected world observation directions");
    settings.observationFrame.reset();
    settings.acceleration = BemAcceleration::fastMultipole;
    expectRejected([&] { (void)makeSimulationRequest(settings, velocities); },
                   "Unsupported baffle FMM silently selected another solver");
    settings.acceleration = BemAcceleration::dense;
    settings.couplingEnabled = true;
    expectRejected([&] { (void)makeSimulationRequest(settings, velocities); },
                   "Unsupported baffle LEM coupling was accepted");
    settings.couplingEnabled = false;
    settings.farField = false;
    expectRejected([&] { (void)makeSimulationRequest(settings, velocities); },
                   "A finite H/V circle with points inside the baffle was accepted");
    settings.farField = true;

    project.simulation.current = false;
    settings.infiniteBaffle->origin.z += .01;
    writeProject(path, project);
    const auto stale = readProject(path);
    sameProject(stale, project);
    require(stale.simulation.request->infiniteBaffle == baffle &&
                stale.simulation.settings.infiniteBaffle != baffle,
            "Editing the baffle overwrote the archived result's domain");
    project.simulation.current = true;
    expectRejected([&] { writeProject(path, project); },
                   "An archive claimed a current result for a different baffle domain");
}

void fastMultipoleRoundTrip(const QTemporaryDir &directory) {
    auto project = fullTestProject(ProjectSource::draft);
    auto &settings = project.simulation.settings;
    settings.acceleration = BemAcceleration::fastMultipole;
    settings.fastMultipole = {24, 2.75, 275, 37};
    auto &numerics = project.simulation.request->numerics;
    numerics.acceleration = settings.acceleration;
    numerics.fastMultipole = settings.fastMultipole;
    auto result = std::make_shared<BemResult>(*project.simulation.result);
    result->acceleration = BemAcceleration::fastMultipole;
    for (auto &frequency : result->frequencies) {
        frequency.fastMultipole = {};
        frequency.fastMultipole.nodeCount = 3;
        frequency.fastMultipole.leafCount = 2;
        frequency.fastMultipole.farInteractionCount = 2;
        frequency.fastMultipole.historicalChebyshevStatistics =
            BemHistoricalChebyshevStatistics{32, 32};
        frequency.residualCheckMilliseconds = 2.75;
        for (auto &excitation : frequency.excitations) {
            excitation.iterations = 37;
            excitation.residualEvaluation = BemResidualEvaluation::higherOrderFastMultipole;
            excitation.fastMultipoleKrylovResidual = 1e-7;
            excitation.fastMultipoleOperatorDiscrepancy = 3e-6;
            excitation.historicalChebyshevEvidence = BemHistoricalChebyshevEvidence{7, 8};
        }
    }
    project.simulation.result = result;
    const auto path = directory.filePath("fast-multipole.json");
    writeProject(path, project);
    sameProject(readProject(path), project);
    for (const unsigned width : {64U, 128U, 256U, 512U, 1024U, 2048U}) {
        settings.fastMultipole.preconditionerBlockWidth = width;
        numerics.fastMultipole.preconditionerBlockWidth = width;
        const auto selectedRequest =
            makeSimulationRequest(settings, std::vector<FaceVelocity>{{1, {1, 0}}});
        require(selectedRequest.numerics.fastMultipole.preconditionerBlockWidth == width,
                "Desktop/study request lost the selected Jacobi-block width");
        writeProject(path, project);
        sameProject(readProject(path), project);
    }
    settings.fastMultipole.preconditionerBlockWidth = 32;
    numerics.fastMultipole.preconditionerBlockWidth = 32;
    const auto request = makeSimulationRequest(settings, std::vector<FaceVelocity>{{1, {1, 0}}});
    require(request.numerics.acceleration == settings.acceleration &&
                request.numerics.fastMultipole == settings.fastMultipole,
            "Desktop/study request lost the selected solver settings");
    settings.acceleration = BemAcceleration::dense;
    project.simulation.current = false;
    writeProject(path, project);
    sameProject(readProject(path), project);
    const auto json = QJsonDocument::fromJson(bytes(path)).object();
    const auto damaged = directory.filePath("bad-fast-multipole.json");
    auto missingWidth = json;
    auto oldSimulation = missingWidth["simulation"].toObject();
    const auto removeWidth = [](QJsonObject container) {
        auto fmm = container["fastMultipole"].toObject();
        fmm.remove("preconditioner_block_width");
        container.insert("fastMultipole", fmm);
        return container;
    };
    oldSimulation.insert("settings", removeWidth(oldSimulation["settings"].toObject()));
    auto requestWithoutWidth = oldSimulation["request"].toObject();
    requestWithoutWidth.insert("numerics", removeWidth(requestWithoutWidth["numerics"].toObject()));
    oldSimulation.insert("request", requestWithoutWidth);
    missingWidth.insert("simulation", oldSimulation);
    saveBytes(damaged, QJsonDocument(missingWidth).toJson());
    sameProject(readProject(damaged), project);
    for (const auto width : {0, 31, 33, 129, 255, 257, 513, 1023, 1025, 2047, 2049}) {
        auto mutation = replaceAt(
            json, {"simulation", "settings", "fastMultipole", "preconditioner_block_width"}, 0,
            width);
        saveBytes(damaged, QJsonDocument(mutation.toObject()).toJson());
        expectRejected([&] { (void)readProject(damaged); },
                       "Unsupported spatial preconditioner block width accepted.");
    }
    const auto currentFmm =
        json["simulation"].toObject()["settings"].toObject()["fastMultipole"].toObject();
    require(currentFmm.contains("leaf_support_capacity") &&
                !currentFmm.contains("interpolationOrder") &&
                !currentFmm.contains("leafPanelCapacity"),
            "New settings retained an obsolete manual order or panel-only leaf contract.");
    saveBytes(damaged, QJsonDocument(preWidebandArchive(json, 11)).toJson());
    sameProject(readProject(damaged), project);
    // Both historical p values migrate to the same automatic solver settings.
    auto historical = QJsonValue(preWidebandArchive(json, 11));
    historical = replaceAt(historical,
                           {"simulation", "settings", "fastMultipole", "interpolationOrder"}, 0, 3);
    historical = replaceAt(
        historical, {"simulation", "request", "numerics", "fastMultipole", "interpolationOrder"}, 0,
        3);
    saveBytes(damaged, QJsonDocument(historical.toObject()).toJson());
    sameProject(readProject(damaged), project);
    const auto frequencyJson =
        json["simulation"].toObject()["result"].toObject()["frequencies"].toArray()[0].toObject();
    for (const auto *field : {"fastMultipole", "residualCheckMilliseconds"}) {
        auto missing = frequencyJson;
        missing.remove(field);
        const auto mutation =
            replaceAt(json, {"simulation", "result", "frequencies", "0"}, 0, missing);
        saveBytes(damaged, QJsonDocument(mutation.toObject()).toJson());
        expectRejected([&] { (void)readProject(damaged); },
                       "FMM result without required execution metadata was accepted");
    }
    auto missingIterations = frequencyJson["excitations"].toArray()[0].toObject();
    missingIterations.remove("iterations");
    const auto missingIterationMutation =
        replaceAt(json, {"simulation", "result", "frequencies", "0", "excitations", "0"}, 0,
                  missingIterations);
    saveBytes(damaged, QJsonDocument(missingIterationMutation.toObject()).toJson());
    expectRejected([&] { (void)readProject(damaged); },
                   "FMM result without its iteration counter was accepted");
    for (const auto *field : {"fastMultipoleKrylovResidual", "fastMultipoleOperatorDiscrepancy",
                              "fastMultipoleAccuracySeed", "fastMultipoleComparisonAccuracySeed",
                              "fastMultipoleVerificationAccuracySeed",
                              "fastMultipoleAccuracyRefinements", "fastMultipoleComparisonLevels",
                              "fastMultipoleVerificationLevels", "historicalChebyshevEvidence"}) {
        auto missing = frequencyJson["excitations"].toArray()[0].toObject();
        missing.remove(field);
        saveBytes(damaged, QJsonDocument(replaceAt(json,
                                                   {"simulation", "result", "frequencies", "0",
                                                    "excitations", "0"},
                                                   0, missing)
                                             .toObject())
                               .toJson());
        expectRejected([&] { (void)readProject(damaged); }, "Incomplete audit metadata accepted");
    }
    for (const auto &mutation :
         {replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "excitations", "0",
                     "residualEvaluation"},
                    0, "exactFmm"),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "excitations", "0",
                     "historicalChebyshevEvidence", "auditOrder"},
                    0, 7),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "excitations", "0",
                     "fastMultipoleOperatorDiscrepancy"},
                    0, numerics.relativeResidualTolerance)}) {
        saveBytes(damaged, QJsonDocument(mutation.toObject()).toJson());
        expectRejected([&] { (void)readProject(damaged); }, "Invalid audit metadata accepted");
    }
    for (const auto &value :
         {replaceAt(json, {"simulation", "settings", "acceleration"}, 0, "automatic"),
          replaceAt(json, {"simulation", "settings", "fastMultipole", "leaf_support_capacity"}, 0,
                    2),
          replaceAt(json,
                    {"simulation", "request", "numerics", "fastMultipole", "leaf_support_capacity"},
                    0, 0),
          replaceAt(json, {"simulation", "request", "numerics", "fastMultipole", "separationRatio"},
                    0, 1),
          replaceAt(json,
                    {"simulation", "request", "numerics", "fastMultipole", "maximumIterations"}, 0,
                    0),
          replaceAt(json, {"simulation", "request", "numerics", "fastMultipole", "restartLength"},
                    0, 1),
          replaceAt(json, {"simulation", "result", "acceleration"}, 0, "dense"),
          replaceAt(json, {"simulation", "result", "frequencies", "0", "residualCheckMilliseconds"},
                    0, -1),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "fastMultipole", "leafCount"}, 0,
                    QString("4")),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "fastMultipole", "leafCount"}, 0,
                    QString("0")),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "fastMultipole",
                     "historicalChebyshevStatistics", "directPanelPairCount"},
                    0, QString("31")),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "fastMultipole",
                     "historicalChebyshevStatistics", "representedFarPanelPairs"},
                    0, QString("33")),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "fastMultipole",
                     "historicalChebyshevStatistics", "directPanelPairCount"},
                    0, QString::number(std::numeric_limits<qulonglong>::max())),
          replaceAt(json, {"simulation", "result", "frequencies", "0", "fastMultipole"}, 0,
                    QJsonObject{{"nodeCount", "0"},
                                {"leafCount", "0"},
                                {"farInteractionCount", "0"},
                                {"directPanelPairCount", "0"},
                                {"representedFarPanelPairs", "0"}}),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "excitations", "0", "iterations"},
                    0, 100000)}) {
        saveBytes(damaged, QJsonDocument(value.toObject()).toJson());
        expectRejected([&] { (void)readProject(damaged); },
                       "Invalid Fast Multipole archive was accepted");
    }
    // A genuine zero-RHS solve may finish without Krylov iterations. Timing is
    // mandatory metadata but a zero duration is not proof that no audit ran.
    auto zeroCounters = replaceAt(
        json, {"simulation", "result", "frequencies", "0", "excitations", "0", "iterations"}, 0, 0);
    zeroCounters =
        replaceAt(zeroCounters,
                  {"simulation", "result", "frequencies", "0", "residualCheckMilliseconds"}, 0, 0);
    saveBytes(damaged, QJsonDocument(zeroCounters.toObject()).toJson());
    const auto zeroCounterResult = readProject(damaged).simulation.result;
    require(zeroCounterResult->frequencies.front().excitations.front().iterations == 0 &&
                zeroCounterResult->frequencies.front().residualCheckMilliseconds == 0,
            "Valid zero FMM counters were rejected or changed");
    auto directFmm = project;
    auto directFmmResult = std::make_shared<BemResult>(*project.simulation.result);
    for (auto &frequency : directFmmResult->frequencies)
        for (auto &excitation : frequency.excitations) {
            excitation.residualEvaluation = BemResidualEvaluation::directGalerkin;
            excitation.historicalChebyshevEvidence = BemHistoricalChebyshevEvidence{};
            excitation.fastMultipoleOperatorDiscrepancy = 0;
            excitation.fastMultipoleKrylovResidual = 2e-6;
            excitation.refinementSteps = 4;
            excitation.iterations = 1125;
        }
    directFmm.simulation.result = directFmmResult;
    const auto directPath = directory.filePath("fmm-direct-audit.json");
    writeProject(directPath, directFmm);
    sameProject(readProject(directPath), directFmm);
    auto baseline = fullTestProject(ProjectSource::draft);
    writeProject(path, baseline);
    auto old = QJsonDocument::fromJson(bytes(path)).object();
    auto simulation = old["simulation"].toObject();
    auto oldSettings = simulation["settings"].toObject();
    oldSettings.remove("acceleration");
    oldSettings.remove("fastMultipole");
    simulation.insert("settings", oldSettings);
    auto oldRequest = simulation["request"].toObject();
    auto oldNumerics = oldRequest["numerics"].toObject();
    oldNumerics.remove("acceleration");
    oldNumerics.remove("fastMultipole");
    oldRequest.insert("numerics", oldNumerics);
    simulation.insert("request", oldRequest);
    auto oldResult = simulation["result"].toObject();
    oldResult.remove("acceleration");
    auto frequencies = oldResult["frequencies"].toArray();
    for (qsizetype i = 0; i < frequencies.size(); ++i) {
        auto frequency = frequencies[i].toObject();
        frequency.remove("fastMultipole");
        frequency.remove("residualCheckMilliseconds");
        auto excitations = frequency["excitations"].toArray();
        for (qsizetype j = 0; j < excitations.size(); ++j) {
            auto excitation = excitations[j].toObject();
            excitation.remove("iterations");
            excitations[j] = excitation;
        }
        frequency.insert("excitations", excitations);
        frequencies[i] = frequency;
    }
    oldResult.insert("frequencies", frequencies);
    simulation.insert("result", oldResult);
    old.insert("simulation", simulation);
    saveBytes(damaged, QJsonDocument(old).toJson());
    sameProject(readProject(damaged), baseline);

    auto directArchive = preWidebandArchive(QJsonDocument::fromJson(bytes(path)).object(), 10);
    auto archivedSimulation = directArchive["simulation"].toObject();
    auto archivedResult = archivedSimulation["result"].toObject();
    auto archivedFrequencies = archivedResult["frequencies"].toArray();
    for (qsizetype i = 0; i < archivedFrequencies.size(); ++i) {
        auto frequency = archivedFrequencies[i].toObject();
        auto excitations = frequency["excitations"].toArray();
        for (qsizetype j = 0; j < excitations.size(); ++j) {
            auto excitation = excitations[j].toObject();
            for (const auto *field :
                 {"residualEvaluation", "fastMultipoleKrylovResidual", "fastMultipoleKrylovTarget",
                  "fastMultipoleOperatorDiscrepancy", "fastMultipoleComparisonOrder",
                  "fastMultipoleAuditOrder"})
                excitation.remove(field);
            excitations[j] = excitation;
        }
        frequency.insert("excitations", excitations);
        archivedFrequencies[i] = frequency;
    }
    archivedResult.insert("frequencies", archivedFrequencies);
    archivedSimulation.insert("result", archivedResult);
    directArchive.insert("simulation", archivedSimulation);
    saveBytes(damaged, QJsonDocument(directArchive).toJson());
    sameProject(readProject(damaged), baseline);
    directArchive.insert("version", 11);
    saveBytes(damaged, QJsonDocument(directArchive).toJson());
    sameProject(readProject(damaged), baseline);
}

void widebandEvidenceRoundTrip(const QTemporaryDir &directory) {
    // Metadata-only synthetic fixture; it does not claim a numerical solve.
    auto project = fullTestProject(ProjectSource::draft);
    project.simulation.settings.acceleration = BemAcceleration::fastMultipole;
    project.simulation.request->numerics.acceleration = BemAcceleration::fastMultipole;
    auto result = std::make_shared<BemResult>(*project.simulation.result);
    result->acceleration = BemAcceleration::fastMultipole;
    result->liveDeviceBytes = 456789;
    const auto n = result->pressureDofCount, t = result->triangleCount;
    for (auto &frequency : result->frequencies) {
        auto &s = frequency.fastMultipole;
        s.nodeCount = 3;
        s.leafCount = 2;
        s.farInteractionCount = 2;
        s.pressurePressureNearPairs = n * n - 12;
        s.pressurePressureFarPairs = 12;
        s.pressureFluxNearPairs = n * t - 16;
        s.pressureFluxFarPairs = 16;
        s.graphExecution = BemFastMultipoleGraphEvidence{};
        auto &graph = *s.graphExecution;
        graph.captureAttempts = 7;
        graph.instantiations = 6;
        graph.launches = 13;
        graph.reuses = 7;
        graph.executedBatches = 41;
        graph.canonicalChunks = 9;
        graph.canonicalBatches = 27;
        graph.invalidatedExecutables = 4;
        graph.budgetBypasses = 2;
        graph.admissionClosures = 1;
        graph.captureMilliseconds = 1.25;
        graph.instantiationMilliseconds = 2.5;
        graph.maximumCachedExecutables = 4;
        graph.maximumCachedNodes = 108;
        graph.peakOwnedHostMetadataBytes = 16384;
        graph.memoryObservationAvailable = true;
        graph.minimumObservedFreeDeviceBytes = 3000000000;
        graph.totalDeviceBytes = 8000000000;
        graph.requiredFreeReserveBytes = 1600000000;
        for (auto &e : frequency.excitations) {
            e.iterations = 37;
            e.residualEvaluation = BemResidualEvaluation::higherOrderFastMultipole;
            e.fastMultipoleKrylovResidual = 1e-7;
            e.fastMultipoleKrylovTarget = 2.5e-7;
            e.fastMultipoleOperatorDiscrepancy = 3e-6;
            e.fastMultipoleAccuracySeed = 1e-5;
            e.fastMultipoleComparisonAccuracySeed = 2.5e-6;
            e.fastMultipoleVerificationAccuracySeed = 6.25e-7;
            e.fastMultipoleAccuracyRefinements = 1;
            e.farQuadratureComparisonDegree = 8;
            e.farQuadratureVerificationDegree = 12;
            e.farQuadratureRefinements = 2;
            e.farQuadratureDiscrepancy = 2e-7;
            e.farQuadraturePreparationMilliseconds = 1.25;
            e.farQuadraturePressureMilliseconds = 2.5;
            e.farQuadratureFluxMilliseconds = 3.75;
            e.fieldQuadratureComparisonDegree = 12;
            e.fieldQuadratureVerificationDegree = 14;
            e.fieldQuadratureDiscrepancy = 1e-7;
            e.fastMultipoleComparisonLevels = {
                {2, BemFastMultipoleRepresentation::regularSpherical, 5, 10, 10},
                {3, BemFastMultipoleRepresentation::planeWave, 8, 16, 16}};
            e.fastMultipoleVerificationLevels = {
                {2, BemFastMultipoleRepresentation::regularSpherical, 6, 12, 12},
                {3, BemFastMultipoleRepresentation::planeWave, 9, 18, 18}};
        }
    }
    project.simulation.result = result;
    const auto path = directory.filePath("wideband-evidence.json"),
               bad = directory.filePath("bad-wideband-evidence.json");
    writeProject(path, project);
    sameProject(readProject(path), project);
    const auto json = QJsonDocument::fromJson(bytes(path)).object();
    require(json["version"].toInt() == 16,
            "Automatic evidence was written in a format an old reader can misinterpret.");
    for (const auto *field : {"residualEvaluation", "quadratureEvidence"}) {
        auto missing = json["simulation"]
                           .toObject()["result"]
                           .toObject()["frequencies"]
                           .toArray()[0]
                           .toObject()["excitations"]
                           .toArray()[0]
                           .toObject();
        missing.remove(field);
        saveBytes(bad, QJsonDocument(replaceAt(json,
                                               {"simulation", "result", "frequencies", "0",
                                                "excitations", "0"},
                                               0, missing)
                                         .toObject())
                           .toJson());
        expectRejected([&] { (void)readProject(bad); },
                       "Format 16 accepted incomplete numerical provenance");
    }
    auto earlierArchive = json;
    earlierArchive.insert("version", 14);
    earlierArchive.remove("comparisonSnapshots");
    auto earlierSimulation = earlierArchive["simulation"].toObject();
    auto earlierResult = earlierSimulation["result"].toObject();
    auto earlierFrequencies = earlierResult["frequencies"].toArray();
    for (qsizetype i = 0; i < earlierFrequencies.size(); ++i) {
        auto f = earlierFrequencies[i].toObject();
        auto stats = f["fastMultipole"].toObject();
        stats.remove("graphExecution");
        f.insert("fastMultipole", stats);
        auto excitations = f["excitations"].toArray();
        for (qsizetype j = 0; j < excitations.size(); ++j) {
            auto e = excitations[j].toObject();
            e.remove("quadratureEvidence");
            e.remove("fastMultipoleKrylovTarget");
            excitations[j] = e;
        }
        f.insert("excitations", excitations);
        earlierFrequencies[i] = f;
    }
    earlierResult.insert("frequencies", earlierFrequencies);
    earlierSimulation.insert("result", earlierResult);
    earlierArchive.insert("simulation", earlierSimulation);
    auto earlierExpected = project;
    auto earlierData = std::make_shared<BemResult>(*result);
    for (auto &f : earlierData->frequencies) {
        f.fastMultipole.graphExecution.reset();
        for (auto &e : f.excitations) {
            e.farQuadratureComparisonDegree = e.farQuadratureVerificationDegree =
                e.farQuadratureRefinements = 0;
            e.fastMultipoleKrylovTarget = 0;
            e.farQuadratureDiscrepancy = e.farQuadraturePreparationMilliseconds =
                e.farQuadraturePressureMilliseconds = e.farQuadratureFluxMilliseconds = 0;
            e.fieldQuadratureComparisonDegree = e.fieldQuadratureVerificationDegree = 0;
            e.fieldQuadratureDiscrepancy = 0;
        }
    }
    earlierExpected.simulation.result = earlierData;
    saveBytes(bad, QJsonDocument(earlierArchive).toJson());
    sameProject(readProject(bad), earlierExpected);
    for (const auto mutation :
         {replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "excitations", "0",
                     "fastMultipoleKrylovTarget"},
                    0, -1.0),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "excitations", "0",
                     "quadratureEvidence", "farVerificationDegree"},
                    0, 6),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "excitations", "0",
                     "quadratureEvidence", "farRefinements"},
                    0, 3),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "excitations", "0",
                     "quadratureEvidence", "fieldDiscrepancy"},
                    0, 1.0)}) {
        saveBytes(bad, QJsonDocument(mutation.toObject()).toJson());
        expectRejected([&] { (void)readProject(bad); }, "Invalid quadrature evidence accepted");
    }
    for (const auto mutation :
         {replaceAt(json, {"simulation", "result", "liveDeviceBytes"}, 0, QString("99999999")),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "fastMultipole",
                     "pressurePressureNearPairs"},
                    0, QString("64")),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "excitations", "0",
                     "fastMultipoleVerificationAccuracySeed"},
                    0, 2.5e-6),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "excitations", "0",
                     "fastMultipoleComparisonLevels", "0", "representation"},
                    0, "chebyshev"),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "excitations", "0",
                     "fastMultipoleVerificationLevels", "0", "momentBandlimit"},
                    0, 4),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "excitations", "0",
                     "fastMultipoleVerificationLevels", "1", "level"},
                    0, 2),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "excitations", "0",
                     "historicalChebyshevEvidence"},
                    0, QJsonObject{{"comparisonOrder", 7}, {"auditOrder", 8}})}) {
        saveBytes(bad, QJsonDocument(mutation.toObject()).toJson());
        expectRejected([&] { (void)readProject(bad); },
                       "Contradictory automatic/historical FMM evidence was accepted");
    }
    // An observed zero free-memory value is different from no observation.
    const auto zeroFree = replaceAt(json,
                                    {"simulation", "result", "frequencies", "0", "fastMultipole",
                                     "graphExecution", "minimumObservedFreeDeviceBytes"},
                                    0, "0");
    saveBytes(bad, QJsonDocument(zeroFree.toObject()).toJson());
    const auto observedZero =
        readProject(bad).simulation.result->frequencies.front().fastMultipole.graphExecution;
    require(observedZero && observedZero->memoryObservationAvailable &&
                observedZero->minimumObservedFreeDeviceBytes == 0,
            "Observed zero free VRAM was confused with unavailable memory evidence.");
    // All-zero execution is a real recorded NoFar case, unlike an absent block.
    auto noFar = std::make_shared<BemResult>(*result);
    for (auto &f : noFar->frequencies) {
        auto &stats = f.fastMultipole;
        stats.farInteractionCount = 0;
        stats.pressurePressureNearPairs = n * n;
        stats.pressurePressureFarPairs = 0;
        stats.pressureFluxNearPairs = n * t;
        stats.pressureFluxFarPairs = 0;
        stats.graphExecution = BemFastMultipoleGraphEvidence{};
        for (auto &e : f.excitations) {
            e.residualEvaluation = BemResidualEvaluation::directGalerkin;
            e.fastMultipoleComparisonAccuracySeed = e.fastMultipoleVerificationAccuracySeed = 0;
            e.fastMultipoleComparisonLevels.clear();
            e.fastMultipoleVerificationLevels.clear();
            e.fastMultipoleOperatorDiscrepancy = 0;
            e.farQuadratureComparisonDegree = e.farQuadratureVerificationDegree =
                e.farQuadratureRefinements = 0;
            e.farQuadratureDiscrepancy = e.farQuadraturePreparationMilliseconds =
                e.farQuadraturePressureMilliseconds = e.farQuadratureFluxMilliseconds = 0;
        }
    }
    auto noFarProject = project;
    noFarProject.simulation.result = noFar;
    writeProject(path, noFarProject);
    sameProject(readProject(path), noFarProject);
    for (const auto &mutation :
         {replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "fastMultipole", "graphExecution",
                     "reuses"},
                    0, "14"),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "fastMultipole", "graphExecution",
                     "instantiations"},
                    0, "8"),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "fastMultipole", "graphExecution",
                     "captureMilliseconds"},
                    0, -1),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "fastMultipole", "graphExecution",
                     "memoryObservationAvailable"},
                    0, false),
          replaceAt(json,
                    {"simulation", "result", "frequencies", "0", "fastMultipole", "graphExecution",
                     "minimumObservedFreeDeviceBytes"},
                    0, "9000000000")}) {
        saveBytes(bad, QJsonDocument(mutation.toObject()).toJson());
        expectRejected([&] { (void)readProject(bad); },
                       "Inconsistent graph execution evidence accepted");
    }
    auto denseProject = fullTestProject(ProjectSource::draft);
    writeProject(path, denseProject);
    const auto denseJson = QJsonDocument::fromJson(bytes(path)).object();
    const auto graphBlock = json["simulation"]
                                .toObject()["result"]
                                .toObject()["frequencies"]
                                .toArray()[0]
                                .toObject()["fastMultipole"]
                                .toObject()["graphExecution"];
    saveBytes(bad, QJsonDocument(replaceAt(denseJson,
                                           {"simulation", "result", "frequencies", "0",
                                            "fastMultipole", "graphExecution"},
                                           0, graphBlock)
                                     .toObject())
                       .toJson());
    expectRejected([&] { (void)readProject(bad); },
                   "Dense result accepted FMM graph execution evidence");
    auto direct = std::make_shared<BemResult>(*result);
    for (auto &frequency : direct->frequencies)
        for (auto &e : frequency.excitations) {
            e.residualEvaluation = BemResidualEvaluation::directGalerkin;
            e.fastMultipoleComparisonAccuracySeed = e.fastMultipoleVerificationAccuracySeed = 0;
            e.fastMultipoleComparisonLevels.clear();
            e.fastMultipoleVerificationLevels.clear();
            e.fastMultipoleOperatorDiscrepancy = 0;
        }
    project.simulation.result = direct;
    writeProject(path, project);
    sameProject(readProject(path), project);
}
void previousFormatRoundTrip(const QTemporaryDir &directory) {
    // Format 14 predates wideband evidence even though the isolated FMM branch
    // used format 12. Detect the archived field contract, not that branch number.
    for (bool multipole : {false, true}) {
        auto project = fullTestProject(ProjectSource::draft);
        if (multipole) {
            project.simulation.settings.acceleration = BemAcceleration::fastMultipole;
            project.simulation.request->numerics.acceleration = BemAcceleration::fastMultipole;
            auto result = std::make_shared<BemResult>(*project.simulation.result);
            result->acceleration = BemAcceleration::fastMultipole;
            for (auto &frequency : result->frequencies) {
                frequency.fastMultipole.nodeCount = 3;
                frequency.fastMultipole.leafCount = 2;
                frequency.fastMultipole.farInteractionCount = 2;
                frequency.fastMultipole.historicalChebyshevStatistics =
                    BemHistoricalChebyshevStatistics{32, 32};
                for (auto &excitation : frequency.excitations)
                    excitation.historicalChebyshevEvidence = BemHistoricalChebyshevEvidence{};
            }
            project.simulation.result = result;
        }
        const auto path =
            directory.filePath(multipole ? "previous-fmm.json" : "previous-dense.json");
        writeProject(path, project);
        auto previous = preWidebandArchive(QJsonDocument::fromJson(bytes(path)).object(), 14);
        auto simulation = previous["simulation"].toObject();
        auto settings = simulation["settings"].toObject();
        settings.remove("singularOrder");
        simulation.insert("settings", settings);
        auto result = simulation["result"].toObject();
        auto frequencies = result["frequencies"].toArray();
        for (qsizetype fi = 0; fi < frequencies.size(); ++fi) {
            auto frequency = frequencies[fi].toObject();
            auto excitations = frequency["excitations"].toArray();
            for (qsizetype ei = 0; ei < excitations.size(); ++ei) {
                auto excitation = excitations[ei].toObject();
                for (const auto *field :
                     {"residualEvaluation", "fastMultipoleKrylovResidual",
                      "fastMultipoleOperatorDiscrepancy", "fastMultipoleComparisonOrder",
                      "fastMultipoleAuditOrder"})
                    excitation.remove(field);
                excitations[ei] = excitation;
            }
            frequency.insert("excitations", excitations);
            frequencies[fi] = frequency;
        }
        result.insert("frequencies", frequencies);
        simulation.insert("result", result);
        previous.insert("simulation", simulation);
        saveBytes(path, QJsonDocument(previous).toJson());
        const auto restored = readProject(path);
        sameProject(restored, project);
        // Rewriting an old result preserves its historical provenance in format 16.
        writeProject(path, restored);
        sameProject(readProject(path), project);
    }
}
void silentExcitationRoundTrip(const QTemporaryDir &directory) {
    // Persistence fixture only; numerical zero-RHS behavior is tested on CUDA.
    auto project = fullTestProject(ProjectSource::draft);
    for (auto &excitation : project.simulation.request->excitations)
        for (auto &face : excitation.faces)
            face.normalVelocity = {};
    auto result = std::make_shared<BemResult>(*project.simulation.result);
    for (auto &frequency : result->frequencies)
        for (auto &excitation : frequency.excitations) {
            auto boundary = *excitation.boundary();
            std::fill(boundary.surfacePressure.begin(), boundary.surfacePressure.end(),
                      std::complex<double>{});
            std::fill(boundary.pressureTrace.begin(), boundary.pressureTrace.end(),
                      std::complex<double>{});
            std::fill(boundary.normalDerivative.begin(), boundary.normalDerivative.end(),
                      std::complex<double>{});
            excitation.setBoundary(std::move(boundary));
            std::fill(excitation.observerPressure.begin(), excitation.observerPressure.end(),
                      std::complex<double>{});
            excitation.impedances.clear();
            excitation.boundaryResidual = 0;
            excitation.refinementSteps = 0;
        }
    project.simulation.result = result;
    const auto path = directory.filePath("silent-explicit-face.json");
    writeProject(path, project);
    sameProject(readProject(path), project);
    const auto json = QJsonDocument::fromJson(bytes(path)).object();
    const auto empty =
        replaceAt(json, {"simulation", "request", "excitations", "0", "faces"}, 0, QJsonArray{});
    saveBytes(path, QJsonDocument(empty.toObject()).toJson());
    expectRejected([&] { (void)readProject(path); },
                   "Empty face assignment was confused with a valid zero RHS");
}
void nearQuadratureDepthRoundTrip(const QTemporaryDir &directory) {
    auto project = fullTestProject(ProjectSource::draft);
    const auto path = directory.filePath("near-quadrature-depth.json");
    for (const unsigned depth : {0u, 8u, 10u}) {
        project.simulation.settings.nearMaximumDepth = depth;
        project.simulation.request->numerics.nearMaximumDepth = depth;
        writeProject(path, project);
        sameProject(readProject(path), project);
    }
    for (const unsigned order : {2u, 8u, 20u}) {
        project.simulation.settings.singularOrder = order;
        project.simulation.request->numerics.singularOrder = order;
        writeProject(path, project);
        sameProject(readProject(path), project);
    }
    const auto json = QJsonDocument::fromJson(bytes(path)).object();
    const auto damaged = directory.filePath("bad-near-quadrature-depth.json");
    for (const double depth : {-1.0, 11.0, 1.5}) {
        const auto replacement =
            replaceAt(json, {"simulation", "settings", "nearMaximumDepth"}, 0, depth);
        saveBytes(damaged, QJsonDocument(replacement.toObject()).toJson());
        expectRejected([&] { (void)readProject(damaged); },
                       "Invalid near quadrature depth was accepted");
    }
    for (const double order : {1.0, 21.0, 3.5}) {
        const auto replacement =
            replaceAt(json, {"simulation", "settings", "singularOrder"}, 0, order);
        saveBytes(damaged, QJsonDocument(replacement.toObject()).toJson());
        expectRejected([&] { (void)readProject(damaged); },
                       "Invalid singular quadrature order was accepted");
    }
    const auto baseline = fullTestProject(ProjectSource::draft);
    writeProject(path, baseline);
    auto old = QJsonDocument::fromJson(bytes(path)).object();
    auto simulation = old["simulation"].toObject();
    auto settings = simulation["settings"].toObject();
    settings.remove("nearMaximumDepth");
    settings.remove("singularOrder");
    simulation.insert("settings", settings);
    old.insert("simulation", simulation);
    saveBytes(path, QJsonDocument(old).toJson());
    sameProject(readProject(path), baseline);
}

void observationFrameRoundTrip(const QTemporaryDir &directory) {
    for (const auto source : {ProjectSource::draft, ProjectSource::freecad}) {
        auto project = fullTestProject(source);
        project.simulation.current = false;
        auto &settings = project.simulation.settings;
        settings.observationFrame = normalizedAcousticFrame({{-.2, .3, .4}, {0, 3, 4}, {4, 3, 4}});
        const auto path = directory.filePath("observation-frame.json");
        writeProject(path, project);
        sameProject(readProject(path), project);
        const auto json = QJsonDocument::fromJson(bytes(path)).object();
        const auto damaged = directory.filePath("bad-observation-frame.json");
        const auto frame =
            json["simulation"].toObject()["settings"].toObject()["observationFrame"].toObject();
        auto missing = frame;
        missing.remove("forward");
        auto extra = frame;
        extra.insert("vertical", QJsonArray{0, 1, 0});
        for (const QJsonValue bad :
             {QJsonValue(QJsonValue::Null), QJsonValue(missing), QJsonValue(extra)}) {
            const auto replacement =
                replaceAt(json, {"simulation", "settings", "observationFrame"}, 0, bad);
            saveBytes(damaged, QJsonDocument(replacement.toObject()).toJson());
            expectRejected([&] { (void)readProject(damaged); },
                           "Invalid frame schema was accepted");
        }
        for (const auto &replacement :
             {replaceAt(json, {"simulation", "settings", "observationFrame", "origin"}, 0,
                        QJsonArray{0, 0}),
              replaceAt(json, {"simulation", "settings", "observationFrame", "origin"}, 0,
                        QJsonArray{QJsonValue::Null, 0, 0}),
              replaceAt(json, {"simulation", "settings", "observationFrame", "forward"}, 0,
                        QJsonArray{0, 0, 0}),
              replaceAt(json, {"simulation", "settings", "observationFrame", "horizontal"}, 0,
                        QJsonArray{0, 0, 0}),
              replaceAt(json, {"simulation", "settings", "observationFrame", "horizontal"}, 0,
                        QJsonArray{0, 3, 4}),
              replaceAt(json, {"simulation", "settings", "observationFrame", "horizontal"}, 0,
                        QJsonArray{1e-10, 3, 4})}) {
            saveBytes(damaged, QJsonDocument(replacement.toObject()).toJson());
            expectRejected([&] { (void)readProject(damaged); },
                           "Invalid observation frame was accepted");
        }
        // Direction magnitudes do not carry physical units; valid non-unit axes
        // are normalized, including the forward component of the horizontal hint.
        auto nonUnit = replaceAt(json, {"simulation", "settings", "observationFrame", "forward"}, 0,
                                 QJsonArray{0, 3, 4});
        nonUnit = replaceAt(nonUnit, {"simulation", "settings", "observationFrame", "horizontal"},
                            0, QJsonArray{4, 3, 4});
        saveBytes(damaged, QJsonDocument(nonUnit.toObject()).toJson());
        sameProject(readProject(damaged), project);

        auto legacy = json;
        auto simulation = legacy["simulation"].toObject();
        auto oldSettings = simulation["settings"].toObject();
        oldSettings.remove("observationFrame");
        simulation.insert("settings", oldSettings);
        legacy.insert("simulation", simulation);
        saveBytes(damaged, QJsonDocument(legacy).toJson());
        auto expected = project;
        expected.simulation.settings.observationFrame.reset();
        const auto restored = readProject(damaged);
        sameProject(restored, expected);
        const Vec3 origin =
            source == ProjectSource::freecad ? restored.cad->acousticOrigin : Vec3{};
        const Vec3 axis =
            source == ProjectSource::freecad ? restored.cad->acousticAxis : Vec3{0, 0, 1};
        const auto effective = simulationFrame(restored.simulation.settings, origin, axis);
        const auto sourceFrame = defaultAcousticFrame(origin, axis);
        samePoint(effective.origin, sourceFrame.origin);
        samePoint(effective.forward, sourceFrame.forward);
        samePoint(effective.horizontal, sourceFrame.horizontal);

        settings.observationFrame->horizontal = settings.observationFrame->forward;
        expectRejected([&] { writeProject(path, project); }, "Parallel frame was written");
        settings.observationFrame = AcousticFrame{{std::numeric_limits<double>::infinity(), 0, 0}};
        expectRejected([&] { writeProject(path, project); }, "Non-finite frame was written");
    }
}

void comparisonSnapshotsRoundTrip(const QTemporaryDir &directory) {
    auto project = fullTestProject(ProjectSource::draft);
    auto first = project;
    auto second = project;
    second.simulation.current = false;
    project.comparisonSnapshots = {
        {QStringLiteral("Referenz"), std::make_shared<const Project>(first)},
        {QStringLiteral("Variante"), std::make_shared<const Project>(second)}};
    const auto path = directory.filePath("comparison-snapshots.hsproj");
    writeProject(path, project);
    auto loaded = readProject(path);
    require(loaded.comparisonSnapshots.size() == 2 &&
                loaded.comparisonSnapshots[0].name == QStringLiteral("Referenz") &&
                loaded.comparisonSnapshots[1].name == QStringLiteral("Variante") &&
                loaded.comparisonSnapshots[0].project->comparisonSnapshots.empty() &&
                !loaded.comparisonSnapshots[1].project->simulation.current,
            "Comparison names, stale state or non-recursive snapshot did not round-trip.");
    sameProject(*loaded.comparisonSnapshots[0].project, first);
    loaded.comparisonSnapshots[0].name = QStringLiteral("Umbenannt");
    writeProject(path, loaded);
    require(readProject(path).comparisonSnapshots[0].name == QStringLiteral("Umbenannt"),
            "Renamed comparison snapshot did not persist.");
    loaded.comparisonSnapshots[1].name = QStringLiteral("Umbenannt");
    expectRejected([&] { writeProject(path, loaded); }, "Duplicate comparison names accepted");
    auto nested = first;
    nested.comparisonSnapshots = {
        {QStringLiteral("Innen"), std::make_shared<const Project>(first)}};
    loaded.comparisonSnapshots = {
        {QStringLiteral("Außen"), std::make_shared<const Project>(nested)}};
    expectRejected([&] { writeProject(path, loaded); }, "Nested comparison snapshots accepted");
    loaded.comparisonSnapshots.clear();
    for (int i = 0; i < 65; ++i)
        loaded.comparisonSnapshots.push_back(
            {QString::number(i), std::make_shared<const Project>(first)});
    expectRejected([&] { writeProject(path, loaded); }, "Too many comparison snapshots accepted");
    QFile oversized(directory.filePath("oversized-project.hsproj"));
    require(oversized.open(QIODevice::WriteOnly) && oversized.resize(256LL * 1024 * 1024 + 1),
            "Could not prepare oversized project fixture");
    oversized.close();
    expectRejected([&] { readProject(oversized.fileName()); }, "Oversized project accepted");
}

void largeFrequencyRoundTrip(const QTemporaryDir &directory) {
    auto project = fullTestProject(ProjectSource::draft);
    constexpr std::size_t count = 5000;
    project.simulation.settings.minimumFrequency = 200;
    project.simulation.settings.maximumFrequency = 5199;
    project.simulation.settings.frequencyCount = count;
    project.simulation.request->frequencies.clear();
    auto result = std::make_shared<BemResult>(*project.simulation.result);
    const auto frame = result->frequencies.front();
    result->frequencies.clear();
    for (std::size_t index = 0; index < count; ++index) {
        const double frequency = 200 + static_cast<double>(index);
        project.simulation.request->frequencies.push_back(frequency);
        auto next = frame;
        next.frequency = frequency;
        result->frequencies.push_back(std::move(next));
    }
    project.simulation.result = result;
    const auto path = directory.filePath("5000-frequencies.json");
    writeProject(path, project);
    const auto restored = readProject(path);
    sameProject(restored, project);
    require(restored.simulation.settings.frequencyCount == count &&
                restored.simulation.request->frequencies.size() == count &&
                restored.simulation.result->frequencies.size() == count,
            "Large frequency sweep was capped or truncated during persistence");
}

void legacyFormatTest(const QTemporaryDir &directory) {
    const auto project = fullTestProject(ProjectSource::draft);
    const auto currentPath = directory.filePath("legacy-source.json");
    writeProject(currentPath, project);
    auto current = QJsonDocument::fromJson(bytes(currentPath)).object();
    require(current["version"].toInt() == 16, "Projects must use integrated FMM format 16");
    current.remove("slot");
    current = preWidebandArchive(current, 11);
    auto version8 = current;
    version8.insert("version", 8);
    version8.remove("draftReversedFaces");
    version8.remove("cadReversedFaces");
    const auto version8Path = directory.filePath("format-8.json");
    saveBytes(version8Path, QJsonDocument(version8).toJson());
    sameProject(readProject(version8Path), project);
    current.remove("draftReversedFaces");
    current.remove("cadReversedFaces");
    require(!current["meshing"].toObject().contains("entryAngleToleranceDegrees"),
            "Current projects must not retain obsolete entry meshing settings");
    auto oldCad = current["cad"].toObject();
    oldCad.remove("brepObjects");
    oldCad.remove("brepMetresPerUnit");
    current.insert("cad", oldCad);
    auto version6 = current;
    version6.insert("version", 6);
    const auto version6Path = directory.filePath("format-6.json");
    saveBytes(version6Path, QJsonDocument(version6).toJson());
    sameProject(readProject(version6Path), project);
    auto version5 = current;
    version5.insert("version", 5);
    const auto version5Path = directory.filePath("format-5.json");
    for (const double oldAngle : {.05, .7, 5.0}) {
        auto sizing = current["meshing"].toObject();
        sizing.insert("entryAngleToleranceDegrees", oldAngle);
        version5.insert("meshing", sizing);
        saveBytes(version5Path, QJsonDocument(version5).toJson());
        sameProject(readProject(version5Path), project);
    }
    for (const auto &sizes : {std::map<std::uint32_t, double>{{6, .002}},
                              std::map<std::uint32_t, double>{{1, .003}, {6, .002}},
                              std::map<std::uint32_t, double>{{1, .001}, {6, .002}}}) {
        QJsonArray legacySizes;
        for (const auto &[id, size] : sizes)
            legacySizes.append(QJsonArray{static_cast<double>(id), size});
        auto oldFaces = version5;
        oldFaces.insert("draftFaceSizes", legacySizes);
        saveBytes(version5Path, QJsonDocument(oldFaces).toJson());
        auto expected = project;
        expected.draftFaceTargetEdgeLengths = {
            {1, sizes.contains(1) ? std::min(sizes.at(1), sizes.at(6)) : sizes.at(6)}};
        sameProject(readProject(version5Path), expected);
    }
    for (const QJsonValue oldAngle : {QJsonValue(0), QJsonValue(5.01), QJsonValue("0.7")}) {
        auto sizing = current["meshing"].toObject();
        sizing.insert("entryAngleToleranceDegrees", oldAngle);
        version5.insert("meshing", sizing);
        saveBytes(version5Path, QJsonDocument(version5).toJson());
        expectRejected([&] { (void)readProject(version5Path); },
                       "Malformed obsolete settings in a legacy project were accepted");
    }
    auto version4 = current;
    version4.insert("version", 4);
    auto oldMeshing = version4["meshing"].toObject();
    oldMeshing.remove("quality");
    oldMeshing.insert("entryAngleToleranceDegrees", .7);
    version4.insert("meshing", oldMeshing);
    const auto version4Path = directory.filePath("format-4.json");
    saveBytes(version4Path, QJsonDocument(version4).toJson());
    auto version4Expected = project;
    version4Expected.meshing.quality = {};
    sameProject(readProject(version4Path), version4Expected);
    auto oldProfile = current["profile"].toObject();
    oldProfile.remove("vertical");
    oldProfile.insert("throatRadius", project.parameters.entry.width / 2);
    oldProfile.insert("length", project.parameters.length);
    auto sharedExpected = project;
    sharedExpected.parameters.mouthMorph = {};
    sharedExpected.parameters.entry.height = sharedExpected.parameters.entry.width;
    sharedExpected.parameters.entry.exponent = 2;
    sharedExpected.meshing = {};
    sharedExpected.meshing.targetEdgeLength = project.meshing.targetEdgeLength;
    sharedExpected.draftFaceTargetEdgeLengths.clear();
    sharedExpected.cadFaceTargetEdgeLengths.clear();
    sharedExpected.simulation.measurements.clear();
    sharedExpected.simulation.settings.beamwidthTargetEnabled = false;
    sharedExpected.simulation.settings.targetHorizontalDegrees = 90;
    sharedExpected.simulation.settings.targetVerticalDegrees = 60;
    sharedExpected.simulation.settings.targetToleranceDegrees = 10;
    sharedExpected.parameters.verticalProfile =
        static_cast<const OsseDirectionalParameters &>(project.parameters.profile);
    auto oldCurrent = current;
    for (const auto key :
         {"length", "entry", "meshing", "draftFaceSizes", "cadFaceSizes", "mouthMorph"})
        oldCurrent.remove(key);
    oldCurrent.insert("version", 3);
    oldCurrent.insert("targetEdgeLength", project.meshing.targetEdgeLength);
    oldCurrent.insert("cadTargetEdgeLength", .0075);
    oldCurrent.insert("cadAdaptiveSizing", current["meshing"].toObject()["adaptiveSizing"]);
    auto oldSimulation = oldCurrent["simulation"].toObject();
    oldSimulation.remove("measurements");
    auto oldSettings = oldSimulation["settings"].toObject();
    oldSettings.remove("beamwidthTarget");
    oldSimulation.insert("settings", oldSettings);
    oldCurrent.insert("simulation", oldSimulation);
    oldCurrent.insert("profile", oldProfile);
    const auto sharedPath = directory.filePath("shared-directions-format-3.json");
    saveBytes(sharedPath, QJsonDocument(oldCurrent).toJson());
    sameProject(readProject(sharedPath), sharedExpected);
    auto version2 = oldCurrent;
    version2.insert("version", 2);
    version2.remove("cadAdaptiveSizing");
    const auto version2Path = directory.filePath("format-2.json");
    saveBytes(version2Path, QJsonDocument(version2).toJson());
    auto version2Expected = sharedExpected;
    sameProject(readProject(version2Path), version2Expected);
    auto cadV3 = oldCurrent;
    cadV3.insert("source", "freecad");
    const auto cadV3Path = directory.filePath("cad-format-3.json");
    saveBytes(cadV3Path, QJsonDocument(cadV3).toJson());
    auto cadExpected = sharedExpected;
    cadExpected.source = ProjectSource::freecad;
    cadExpected.meshing.targetEdgeLength = .0075;
    cadExpected.meshing.adaptive = project.meshing.adaptive;
    sameProject(readProject(cadV3Path), cadExpected);
    auto withoutFixed = version2;
    auto unconstrained = [](QJsonValue value) {
        auto mesh = value.toObject();
        mesh.remove("fixedVertices");
        return mesh;
    };
    withoutFixed.insert("preview", unconstrained(withoutFixed["preview"]));
    withoutFixed.insert("mesh", unconstrained(withoutFixed["mesh"]));
    auto cad = withoutFixed["cad"].toObject();
    cad.insert("mesh", unconstrained(cad["mesh"]));
    withoutFixed.insert("cad", cad);
    auto simulation = withoutFixed["simulation"].toObject();
    simulation.insert("solvedMesh", unconstrained(simulation["solvedMesh"]));
    withoutFixed.insert("simulation", simulation);
    const auto withoutFixedPath = directory.filePath("without-fixed-vertices.json");
    saveBytes(withoutFixedPath, QJsonDocument(withoutFixed).toJson());
    const auto unrestricted = readProject(withoutFixedPath);
    require(unrestricted.preview->fixedVertices.empty() &&
                unrestricted.mesh->mesh.fixedVertices.empty() &&
                unrestricted.cad->mesh.fixedVertices.empty() &&
                unrestricted.simulation.mesh->fixedVertices.empty(),
            "Meshes without optional constraints acquired fixed vertices on load");
    // Exactly the public format-1 keys written by the original Draft editor.
    QJsonObject legacy;
    for (const auto key : {"format", "units", "profileType", "profile", "horizontalAngle",
                           "verticalAngle", "crossSectionExponent", "wallThickness",
                           "targetEdgeLength", "guideEnabled", "transition", "guide"})
        legacy.insert(key, oldCurrent.value(key));
    legacy.insert("profile", oldProfile);
    legacy.insert("version", 1);
    const auto path = directory.filePath("legacy-v1.json");
    saveBytes(path, QJsonDocument(legacy).toJson());
    const auto restored = readProject(path);
    auto legacyParameters = sharedExpected.parameters;
    legacyParameters.leftProfile = legacyParameters.profile;
    legacyParameters.bottomProfile = legacyParameters.verticalProfile;
    legacyParameters.leftAngle = legacyParameters.horizontalAngle;
    legacyParameters.bottomAngle = legacyParameters.verticalAngle;
    sameDraft(restored.parameters, legacyParameters);
    close(restored.meshing.targetEdgeLength, project.meshing.targetEdgeLength,
          "Legacy mesh size changed");
    require(restored.source == ProjectSource::draft && !restored.cad && !restored.mesh &&
                !restored.simulation.result && !restored.simulation.current,
            "Legacy Draft load invented CAD, mesh or simulation state");
}

void diffractionSlotRoundTrip(const QTemporaryDir &directory) {
    Project project;
    auto &p = project.parameters;
    p.slot.aperture = {.011234567891234, .036789123456789, 4.25};
    p.slot.adapterLength = .023456789123456;
    p.slot.straightLength = .004567891234567;
    p.slot.cornerRadius = .002345678912345;
    p.slot.exitRadius = .018123456789123;
    p.transition.startFraction = .3;
    project.draftFaceTargetEdgeLengths = {{6, .002}, {7, .0025}, {8, .003}};
    project.draftFaceColors = {{6, {255, 145, 0}}, {7, {140, 145, 150}}, {8, {0, 145, 255}}};
    project.draftReversedFaceIds = {6, 7, 8};
    const auto path = directory.filePath("diffraction-slot.hsproj");
    for (const bool enabled : {false, true}) {
        p.slot.enabled = enabled;
        for (const double transitionLength : {0., .006123456789123}) {
            p.slot.transitionLength = transitionLength;
            writeProject(path, project);
            sameProject(readProject(path), project);
        }
    }
    for (const auto contour : {SlotContour::superellipse, SlotContour::roundedRectangle}) {
        p.slot.contour = contour;
        p.slot.exit = SlotExit::circularFillet;
        p.slot.exitRadius = 0;
        writeProject(path, project);
        sameProject(readProject(path), project);
    }
    auto circular = project;
    circular.parameters.slot.contour = SlotContour::roundedRectangle;
    circular.parameters.slot.cornerRadius = .002;
    circular.parameters.slot.exit = SlotExit::circularFillet;
    circular.parameters.slot.exitRadius = .002;
    circular.parameters.slot.straightLength = .010;
    circular.parameters.profile.k = circular.parameters.verticalProfile.k = 0;
    circular.parameters.profile.terminationScale =
        circular.parameters.verticalProfile.terminationScale = 0;
    circular.parameters.horizontalAngle = .5;
    circular.parameters.verticalAngle = .3;
    for (const double radius : {.002, .018}) {
        circular.parameters.slot.exitRadius = radius;
        writeProject(path, circular);
        sameProject(readProject(path), circular);
    }
    p.slot.contour = SlotContour::superellipse;
    p.slot.exit = SlotExit::axialBlend;
    p.slot.exitRadius = .018123456789123;
    writeProject(path, project);
    const auto original = QJsonDocument::fromJson(bytes(path)).object();
    const auto reject = [&](const QStringList &keys, const QJsonValue &value) {
        const auto damaged = replaceAt(original, keys, 0, value);
        saveBytes(path, QJsonDocument(damaged.toObject()).toJson());
        expectRejected([&] { (void)readProject(path); }, "Invalid slot definition was accepted");
    };
    reject({"slot"}, QJsonValue::Null);
    reject({"slot", "enabled"}, "true");
    reject({"slot", "contour"}, 2);
    reject({"slot", "exit"}, 2);
    reject({"slot", "cornerRadius"}, -.001);
    reject({"slot", "exitRadius"}, -.001);
    reject({"slot", "cornerRadius"}, "2");
    reject({"slot", "exitRadius"}, QJsonValue::Null);
    for (const auto field : {"width", "height", "exponent"}) {
        reject({"slot", "aperture", field}, 0);
        reject({"slot", "aperture", field}, QJsonValue::Null);
    }
    reject({"slot", "aperture", "exponent"}, 1);
    reject({"slot", "adapterLength"}, 0);
    reject({"slot", "straightLength"}, -.001);
    reject({"slot", "transitionLength"}, -.001);
    reject({"slot", "transitionLength"}, .2);
    reject({"transition", "startFraction"}, 0);
    for (const auto field : {"enabled", "aperture", "adapterLength", "straightLength",
                             "transitionLength", "contour", "cornerRadius", "exit", "exitRadius"}) {
        auto missing = original["slot"].toObject();
        missing.remove(field);
        reject({"slot"}, missing);
    }
    auto missingSlot = original;
    missingSlot.remove("slot");
    saveBytes(path, QJsonDocument(missingSlot).toJson());
    expectRejected([&] { (void)readProject(path); }, "Format 14 omitted required slot parameters");
    auto format13 = original;
    format13.insert("version", 13);
    format13.remove("comparisonSnapshots");
    auto oldSlot = format13["slot"].toObject();
    for (const auto field : {"contour", "cornerRadius", "exit", "exitRadius"})
        oldSlot.remove(field);
    format13.insert("slot", oldSlot);
    saveBytes(path, QJsonDocument(format13).toJson());
    auto expected13 = project;
    expected13.parameters.slot.cornerRadius = DiffractionSlotParameters{}.cornerRadius;
    expected13.parameters.slot.exitRadius = DiffractionSlotParameters{}.exitRadius;
    sameProject(readProject(path), expected13);
    auto old = missingSlot;
    old.insert("version", 12);
    old.remove("comparisonSnapshots");
    saveBytes(path, QJsonDocument(old).toJson());
    auto expected = project;
    expected.parameters.slot = {};
    sameProject(readProject(path), expected);
    p.slot.enabled = false;
    p.slot.aperture.width = 0;
    expectRejected([&] { writeProject(path, project); }, "Invalid dormant slot dimensions saved");
}

void asymmetricDraftRoundTrip(const QTemporaryDir &directory) {
    Project project;
    auto &p = project.parameters;
    p.horizontalAngle = .70;
    p.verticalAngle = .50;
    p.profile = {.02, 1.03, .07, .95, 3.9};
    p.verticalProfile = {.04, .90, .09, .92, 4.7};
    p.leftAngle = .35;
    p.bottomAngle = .25;
    p.leftProfile = {.03, .85, .08, .93, 3.7};
    p.bottomProfile = {.01, 1.15, .12, .96, 4.3};
    const auto path = directory.filePath("asymmetric-draft.hsproj");
    for (const bool enabled : {false, true}) {
        p.asymmetricEnabled = enabled;
        writeProject(path, project);
        sameDraft(readProject(path).parameters, p);
    }
    auto original = QJsonDocument::fromJson(bytes(path)).object();
    auto old = original;
    old.insert("version", 10);
    old.remove("comparisonSnapshots");
    old.remove("slot");
    old.remove("asymmetry");
    saveBytes(path, QJsonDocument(old).toJson());
    auto expected = p;
    expected.asymmetricEnabled = false;
    expected.leftProfile = expected.profile;
    expected.bottomProfile = expected.verticalProfile;
    expected.leftAngle = expected.horizontalAngle;
    expected.bottomAngle = expected.verticalAngle;
    sameDraft(readProject(path).parameters, expected);
    for (const auto field :
         {"enabled", "leftAngle", "bottomAngle", "leftProfile", "bottomProfile"}) {
        auto invalid = original;
        auto asymmetry = invalid["asymmetry"].toObject();
        asymmetry.remove(field);
        invalid.insert("asymmetry", asymmetry);
        saveBytes(path, QJsonDocument(invalid).toJson());
        expectRejected([&] { (void)readProject(path); }, "Incomplete asymmetry was accepted");
    }
    auto incomplete = original;
    incomplete.remove("asymmetry");
    saveBytes(path, QJsonDocument(incomplete).toJson());
    expectRejected([&] { (void)readProject(path); }, "Format 11 without asymmetry was accepted");
    p.asymmetricEnabled = false;
    p.leftProfile.terminationQ = 0;
    expectRejected([&] { writeProject(path, project); }, "Invalid dormant profile was saved");
    p.leftProfile = expected.leftProfile;
    p.bottomAngle = -1;
    expectRejected([&] { writeProject(path, project); }, "Invalid dormant angle was saved");
}

void relaxedMeshRoundTrip(const QTemporaryDir &directory) {
    auto project = fullTestProject(ProjectSource::draft);
    project.simulation = {};
    auto flattened = std::make_shared<MeshingResult>(*project.mesh);
    for (auto &vertex : flattened->mesh.vertices) {
        vertex.x *= .0001;
        vertex.y *= .0001;
    }
    flattened->quality = inspectMesh(flattened->mesh);
    require(flattened->quality.watertight() && flattened->quality.minimumAngleDegrees < 20 &&
                flattened->quality.maximumAspectRatio > 4,
            "Relaxed-quality archive fixture must be closed and violate default quality limits");
    project.mesh = flattened;
    project.preview = std::make_shared<SurfaceMesh>(flattened->mesh);
    const auto path = directory.filePath("relaxed-quality.json");
    for (const bool enabled : {true, false}) {
        project.meshing.quality = {enabled, .001, 1000000};
        writeProject(path, project);
        sameProject(readProject(path), project);
    }
    const auto saved = bytes(path);
    auto open = std::make_shared<MeshingResult>(*flattened);
    open->mesh.triangles.pop_back();
    project.mesh = open;
    project.preview = std::make_shared<SurfaceMesh>(open->mesh);
    expectRejected([&] { writeProject(path, project); },
                   "Disabling optional quality limits must not permit an open BEM boundary");
    require(bytes(path) == saved, "Structurally invalid mesh damaged a saved relaxed-quality mesh");
}

void invalidWriteAndAtomicTests(const QTemporaryDir &directory) {
    const auto valid = fullTestProject(ProjectSource::draft);
    const auto path = directory.filePath("atomic.json");
    writeProject(path, valid);
    const auto original = bytes(path);
    auto invalid = valid;
    invalid.meshing.targetEdgeLength = std::numeric_limits<double>::quiet_NaN();
    expectRejected([&] { writeProject(path, invalid); }, "Nonfinite project data was written");
    require(bytes(path) == original, "Rejected write damaged the previously saved project");
    auto missing = valid;
    missing.simulation.mesh.reset();
    expectRejected([&] { writeProject(path, missing); }, "Result without solved mesh was written");
    missing = valid;
    missing.simulation.request.reset();
    expectRejected([&] { writeProject(path, missing); },
                   "Result without original request was written");
    missing = valid;
    missing.simulation.result.reset();
    expectRejected([&] { writeProject(path, missing); },
                   "Incomplete request/mesh archive was written");
    require(bytes(path) == original, "Archive validation failure changed the existing project");
#ifdef _WIN32
    struct FileLock {
        HANDLE value;
        ~FileLock() {
            if (value != INVALID_HANDLE_VALUE)
                CloseHandle(value);
        }
    } lock{CreateFileW(reinterpret_cast<LPCWSTR>(path.utf16()), GENERIC_READ, FILE_SHARE_READ,
                       nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr)};
    require(lock.value != INVALID_HANDLE_VALUE, "Atomic-write test cannot lock its existing file");
    auto changed = valid;
    changed.meshing.targetEdgeLength = .009;
    expectRejected([&] { writeProject(path, changed); },
                   "Locked-file atomic commit unexpectedly succeeded");
    require(bytes(path) == original,
            "Failed atomic commit replaced or truncated the existing file");
#endif
    sameProject(readProject(path), valid);
}

void damagedFileTests(const QTemporaryDir &directory) {
    const auto source = directory.filePath("damage-source.json");
    writeProject(source, fullTestProject(ProjectSource::freecad));
    const auto original = QJsonDocument::fromJson(bytes(source)).object();
    const auto damaged = directory.filePath("damaged.json");
    auto rejectValue = [&](const QStringList &path, const QJsonValue &value) {
        const auto object = replaceAt(original, path, 0, value).toObject();
        saveBytes(damaged, QJsonDocument(object).toJson());
        expectRejected([&] { (void)readProject(damaged); }, "Corrupt project JSON was accepted");
    };
    rejectValue({"profile", "vertical"}, QJsonValue::Null);
    rejectValue({"profile", "vertical"}, QJsonObject{});
    rejectValue({"profile", "vertical", "k"}, "0.87");
    rejectValue({"profile", "vertical", "terminationScale"}, QJsonValue::Null);
    auto extraVertical = original["profile"].toObject()["vertical"].toObject();
    extraVertical.insert("unknown", 1);
    rejectValue({"profile", "vertical"}, extraVertical);
    auto extraDirection = original["profile"].toObject();
    extraDirection.insert("unexpectedDirection", QJsonObject{});
    rejectValue({"profile"}, extraDirection);
    rejectValue({"source"}, "unknown-source");
    rejectValue({"entry", "width"}, 0);
    rejectValue({"entry", "height"}, -.01);
    rejectValue({"entry", "exponent"}, QJsonValue::Null);
    rejectValue({"meshing", "adaptiveSizing", "enabled"}, "true");
    rejectValue({"meshing", "adaptiveSizing", "minimumEdgeLength"}, -.001);
    rejectValue({"meshing", "adaptiveSizing", "minimumEdgeLength"}, .1);
    rejectValue({"meshing", "adaptiveSizing", "geometricTolerance"}, 0);
    rejectValue({"meshing", "adaptiveSizing", "growthRate"}, 1.1);
    auto obsoleteMeshing = original["meshing"].toObject();
    obsoleteMeshing.insert("entryAngleToleranceDegrees", .7);
    rejectValue({"meshing"}, obsoleteMeshing);
    rejectValue({"meshing", "quality", "enabled"}, "false");
    rejectValue({"meshing", "quality", "minimumAngleDegrees"}, 0);
    rejectValue({"meshing", "quality", "minimumAngleDegrees"}, 60);
    rejectValue({"meshing", "quality", "maximumAspectRatio"}, 1);
    rejectValue({"meshing", "quality", "maximumAspectRatio"}, 1000001);
    rejectValue({"meshing", "quality"}, QJsonValue::Null);
    rejectValue({"draftFaceSizes"}, QJsonArray{QJsonValue(QJsonArray{9, .002})});
    rejectValue({"draftFaceSizes"}, QJsonArray{QJsonValue(QJsonArray{10, .002})});
    rejectValue({"cadFaceSizes"}, QJsonArray{QJsonValue(QJsonArray{999999, .002})});
    rejectValue({"draftReversedFaces"}, QJsonArray{0});
    rejectValue({"draftReversedFaces"}, QJsonArray{1, 1});
    rejectValue({"draftReversedFaces"}, QJsonArray{9});
    rejectValue({"cadReversedFaces"}, QJsonArray{1, 1});
    rejectValue({"cadReversedFaces"}, QJsonArray{999999});
    rejectValue({"cadReversedFaces"}, QJsonArray{1});
    rejectValue({"draftFaceSizes"}, QJsonArray{QJsonValue(QJsonArray{5, .02})});
    rejectValue({"draftFaceSizes"}, QJsonArray{QJsonValue(QJsonArray{5, .0000001})});
    rejectValue({"revision"}, static_cast<double>(largeRevision));
    rejectValue({"revision"}, "09007199254741009");
    rejectValue({"revision"}, "18446744073709551616");
    rejectValue({"geometryRevision"}, QString::number(largeRevision + 1));
    rejectValue({"simulation", "revision"}, QString::number(largeRevision + 1));
    rejectValue({"simulation", "revision"}, QString::number(largeRevision - 1));
    rejectValue({"simulation", "solvedMesh"}, QJsonValue::Null);
    rejectValue({"simulation", "request"}, QJsonValue::Null);
    rejectValue({"simulation", "result"}, QJsonValue::Null);
    rejectValue({"simulation", "settings", "normalizeToOnAxis"}, "true");
    rejectValue({"simulation", "settings", "normalizeToOnAxis"}, QJsonValue::Null);
    rejectValue({"simulation", "settings", "normAngle"}, 180.1);
    rejectValue({"simulation", "settings", "normAngle"}, -180.1);
    rejectValue({"simulation", "settings", "normAngle"}, QJsonValue::Null);
    rejectValue({"simulation", "settings", "DbScale", "minimum"}, 3);
    rejectValue({"simulation", "settings", "DbScale", "maximum"}, -30);
    rejectValue({"simulation", "settings", "DbScale", "step"}, 0);
    rejectValue({"simulation", "settings", "DbScale", "step"}, -.1);
    rejectValue({"simulation", "settings", "DbScale", "step"}, .001);
    rejectValue({"simulation", "settings", "DbScale"}, QJsonValue::Null);
    rejectValue({"simulation", "request", "excitations", "0", "faces", "0", "faceId"}, 999999);
    rejectValue({"simulation", "request", "excitations", "0", "faces", "0", "faceId"}, 3);
    rejectValue({"simulation", "request", "excitations", "0", "faces", "0", "normalVelocity", "1"},
                QJsonValue::Null);
    rejectValue({"simulation", "request", "frequencies", "0"}, -10);
    rejectValue({"simulation", "result", "frequencies", "0", "frequency"}, 999);
    rejectValue({"simulation", "result", "frequencies", "0", "excitations", "0", "surfacePressure"},
                QJsonArray{QJsonValue(QJsonArray{1, 2})});
    rejectValue(
        {"simulation", "result", "frequencies", "0", "excitations", "0", "observerPressure"},
        QJsonArray{});
    rejectValue({"simulation", "result", "frequencies", "0", "excitations", "0", "surfacePressure",
                 "0", "1"},
                QJsonValue::Null);
    rejectValue({"simulation", "result", "frequencies"}, QJsonArray{});
    rejectValue({"simulation", "result", "triangleCount"}, 9);
    rejectValue({"simulation", "result", "pressureDofCount"}, 7);
    rejectValue({"simulation", "result", "frequencies", "0", "excitations"}, QJsonArray{});
    rejectValue(
        {"simulation", "result", "frequencies", "0", "excitations", "0", "boundaryResidual"}, 1);
    rejectValue({"simulation", "result", "frequencies", "0", "excitations", "0", "impedances", "0",
                 "faceId"},
                999999);
    rejectValue({"simulation", "result", "frequencies", "0", "excitations", "0", "impedances", "0",
                 "volumeImpedance"},
                QJsonArray{1, 2});
    rejectValue({"mesh", "vertices", "0", "0"}, QJsonValue::Null);
    rejectValue({"mesh", "triangles", "0", "0"}, 999999);
    rejectValue({"mesh", "triangles", "0", "4"}, 99);
    rejectValue({"mesh", "triangles", "0", "3"}, 0);
    rejectValue({"mesh", "fixedVertices"}, QJsonArray{0, 0});
    rejectValue({"mesh", "fixedVertices"}, QJsonArray{-1});
    rejectValue({"mesh", "fixedVertices"}, QJsonArray{999999});
    rejectValue({"mesh", "fixedVertices"}, QJsonArray{.5});
    rejectValue({"mesh", "fixedVertices"}, QJsonValue::Null);
    rejectValue({"cad", "faces", "0", "objectId"}, "unknown-object");
    rejectValue({"cad", "faces", "0", "state"}, 99);
    rejectValue({"cad", "faces", "0", "normalVelocity", "0"}, QJsonValue::Null);
    rejectValue({"cad", "mesh", "triangles", "0", "3"}, 999999);
    rejectValue({"cad", "objects", "0", "placement", "0"}, 2);
    rejectValue({"cad"}, QJsonValue::Null);
    rejectValue({"version"}, 999);
    saveBytes(damaged, "{\"format\":\"HornSculpt\",\"version\":2,");
    expectRejected([&] { (void)readProject(damaged); }, "Truncated JSON was accepted");
}

void ctaExportContract() {
    SimulationArchive a;
    SimulationSettings settings;
    const std::vector<FaceVelocity> velocity{{1, {1, 0}}};
    BemRequest request = makeSimulationRequest(settings, velocity);
    request.observers.resize(75);
    // Use a valid legacy ±90° raster generated from the full request.
    const auto full = makeSimulationRequest(settings, velocity);
    request.observers.clear();
    request.observers.push_back(full.observers[0]);
    for (int i = 18; i <= 54; ++i)
        request.observers.push_back(full.observers[1 + i]);
    for (int i = 18; i <= 54; ++i)
        request.observers.push_back(full.observers[74 + i]);
    a.request = request;
    a.revision = 23;
    a.current = false;
    a.resultComplete = false;
    auto result = std::make_shared<BemResult>();
    BemFrequencyResult f;
    f.frequency = 1000;
    BemExcitationResult first;
    first.name = "A";
    first.observerPressure.assign(75, {.002, 0});
    BemExcitationResult second;
    second.name = "B";
    second.observerPressure.assign(75, {.004, 0});
    f.excitations = {first, second};
    result->frequencies.push_back(f);
    a.result = result;
    const auto one = ctaResults(*a.request, *a.result, 0),
               two = ctaResults(*a.request, *a.result, 1);
    require(one.size() == 1 && two.size() == 1 && one[0].curves.onAxisDb && two[0].curves.onAxisDb,
            "CTA excitation selection missing");
    close(*two[0].curves.onAxisDb - *one[0].curves.onAxisDb, 20 * std::log10(2.),
          "CTA selected wrong excitation");
    require(!two[0].curves.soundPowerDb, "Incomplete rear angles fabricated power");
    auto malformed = a;
    malformed.request->observers[5].x += .1;
    expectRejected([&] { (void)ctaCsv(malformed, 1); }, "Deformed CTA observer raster accepted");
    malformed = a;
    std::swap(malformed.request->observers[4], malformed.request->observers[5]);
    expectRejected([&] { (void)ctaCsv(malformed, 1); }, "Reordered CTA observer raster accepted");
    malformed = a;
    malformed.request->farField->origin.x += .1;
    expectRejected([&] { (void)ctaCsv(malformed, 1); }, "Shifted far-field origin accepted");
    for (const double scale : {1e-4, 1e4}) {
        auto scaled = a;
        for (auto &point : scaled.request->observers) {
            point.x *= scale;
            point.y *= scale;
            point.z *= scale;
        }
        require(!ctaCsv(scaled, 1).isEmpty(), "A scaled valid directional raster was rejected");
    }
    const auto csv = ctaCsv(a, 1);
    require(csv.contains("excitation=B") && csv.contains("revision=23; current=false") &&
                csv.contains("completeness=partial") &&
                csv.contains("70_direction_approximation") &&
                csv.contains("reference_pressure_Pa_RMS=2e-5") && csv.contains(";;;"),
            "CTA export omitted metadata or empty incomplete values");
}

void resultMemoryEstimateContract() {
    ResultMemoryDimensions d{100, 60, 0, 147, 3, 10, 2, 1};
    const auto base = estimateResultBytes(d);
    d.exteriorTriangles = 17;
    const auto withRayleigh = estimateResultBytes(d);
    require(withRayleigh - base == 10.L * 2 * 17 * sizeof(std::complex<double>),
            "Rayleigh boundary data missing from result memory estimate.");
    d.snapshots = 4;
    const auto fourSnapshots = estimateResultBytes(d);
    d.exteriorTriangles = 0;
    require(fourSnapshots - estimateResultBytes(d) == 10.L * 2 * 17 * sizeof(std::complex<double>),
            "Shared boundary data was multiplied by snapshot count.");
    require(fourSnapshots > withRayleigh, "Independent snapshot data was not counted per owner.");
    d.frequencies = std::numeric_limits<std::size_t>::max();
    require(estimateResultBytes(d) >
                static_cast<long double>(std::numeric_limits<std::size_t>::max()),
            "Large result estimate overflowed or wrapped.");
}

} // namespace

int main(int argc, char **argv) {
    QCoreApplication application(argc, argv);
    try {
        if (argc == 3 && QString::fromLocal8Bit(argv[1]) == "--export-cad") {
            auto fixture = fullTestProject(ProjectSource::freecad);
            fixture.parameters = DraftParameters{};
            fixture.simulation = {};
            writeProject(QString::fromLocal8Bit(argv[2]), fixture);
            std::cout << "Exported offline CAD geometry fixture without acoustic results.\n";
            return 0;
        }
        const bool legacyOnly = argc == 2 && QString::fromLocal8Bit(argv[1]) == "--legacy-only";
        const bool regionsOnly =
            argc == 2 && QString::fromLocal8Bit(argv[1]) == "--subdomain-settings-only";
        const bool axisymmetricOnly =
            argc == 2 && QString::fromLocal8Bit(argv[1]) == "--axisymmetric-only";
        require(argc == 1 || axisymmetricOnly || legacyOnly || regionsOnly,
                "Unexpected project test arguments");
        QTemporaryDir directory;
        require(directory.isValid(), "Temporary project directory cannot be created");
        if (legacyOnly) {
            legacyFormatTest(directory);
            previousFormatRoundTrip(directory);
            return 0;
        }
        if (regionsOnly) {
            roundTripTests(directory);
            return 0;
        }
        if (axisymmetricOnly) {
            axisymmetricRoundTrip(directory);
            interactionPlanningRoundTrip(directory);
            std::cout << "Axisymmetric surface contracts, weighted traces, planning and archive "
                         "rejection tests passed.\n";
            return 0;
        }
        const auto runStage = [](const char *name, auto test) {
            try {
                test();
            } catch (const std::exception &e) {
                throw std::runtime_error(std::string(name) + ": " + e.what());
            }
        };
        runStage("resultMemoryEstimateContract", [&] { resultMemoryEstimateContract(); });
        runStage("roundTripTests", [&] { roundTripTests(directory); });
        runStage("coupledRoundTrip", [&] { coupledRoundTrip(directory); });
        runStage("unfinishedLemGraphRoundTrip", [&] { unfinishedLemGraphRoundTrip(directory); });
        runStage("interfaceRoundTrip", [&] { interfaceRoundTrip(directory); });
        runStage("weightedPortRoundTrip", [&] { weightedPortRoundTrip(directory); });
        runStage("acousticSlitRoundTrip", [&] { acousticSlitRoundTrip(directory); });
        runStage("closedSlitCavityRoundTrip", [&] { closedSlitCavityRoundTrip(directory); });
        runStage("orientationRoundTrip", [&] { orientationRoundTrip(directory); });
        runStage("exactCadRoundTrip", [&] { exactCadRoundTrip(directory); });
        runStage("farFieldRoundTrip", [&] { farFieldRoundTrip(directory); });
        runStage("faceColorsRoundTrip", [&] { faceColorsRoundTrip(directory); });
        runStage("angularResponsesRoundTrip", [&] { angularResponsesRoundTrip(directory); });
        runStage("ctaExportContract", [&] { ctaExportContract(); });
        runStage("savedFieldTracesRoundTrip", [&] { savedFieldTracesRoundTrip(directory); });
        runStage("axisymmetricRoundTrip", [&] { axisymmetricRoundTrip(directory); });
        runStage("interactionPlanningRoundTrip", [&] { interactionPlanningRoundTrip(directory); });
        runStage("infiniteBaffleRoundTrip", [&] { infiniteBaffleRoundTrip(directory); });
        runStage("fastMultipoleRoundTrip", [&] { fastMultipoleRoundTrip(directory); });
        runStage("nearQuadratureDepthRoundTrip", [&] { nearQuadratureDepthRoundTrip(directory); });
        runStage("observationFrameRoundTrip", [&] { observationFrameRoundTrip(directory); });
        runStage("comparisonSnapshotsRoundTrip", [&] { comparisonSnapshotsRoundTrip(directory); });
        runStage("largeFrequencyRoundTrip", [&] { largeFrequencyRoundTrip(directory); });
        runStage("legacyFormatTest", [&] { legacyFormatTest(directory); });
        runStage("asymmetricDraftRoundTrip", [&] { asymmetricDraftRoundTrip(directory); });
        runStage("diffractionSlotRoundTrip", [&] { diffractionSlotRoundTrip(directory); });
        runStage("relaxedMeshRoundTrip", [&] { relaxedMeshRoundTrip(directory); });
        runStage("invalidWriteAndAtomicTests", [&] { invalidWriteAndAtomicTests(directory); });
        runStage("damagedFileTests", [&] { damagedFileTests(directory); });
        runStage("previousFormatRoundTrip", [&] { previousFormatRoundTrip(directory); });
        runStage("widebandEvidenceRoundTrip", [&] { widebandEvidenceRoundTrip(directory); });
        runStage("silentExcitationRoundTrip", [&] { silentExcitationRoundTrip(directory); });
        std::cout << "Draft/CAD/mesh/synthetic archive persistence and rejection tests passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
