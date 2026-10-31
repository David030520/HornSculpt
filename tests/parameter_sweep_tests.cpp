#include "fixture_geometry.hpp"
#include "sweep_model.hpp"
#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonArray>
#include <QJsonObject>
#include <QTemporaryDir>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
template <class F> void rejects(F operation) {
    try {
        operation();
    } catch (const std::exception &) {
        return;
    }
    throw std::runtime_error("An invalid study definition was accepted.");
}
QByteArray contents(const QString &path) {
    QFile file(path);
    require(file.open(QIODevice::ReadOnly), "Cannot read study fixture.");
    return file.readAll();
}
void save(const QString &path, const QByteArray &bytes) {
    QFile file(path);
    require(file.open(QIODevice::WriteOnly | QIODevice::Truncate), "Cannot write study fixture.");
    require(file.write(bytes) == bytes.size(), "Incomplete study fixture.");
}
void gridAndBindings() {
    const std::vector<SweepAxis> axes{{"length", .16, .20, 3}, {"horizontalAngle", .5, .7, 2}};
    const auto grid = makeSweepGrid(axes);
    require(grid.size() == 6, "The two-axis Cartesian product is incomplete.");
    constexpr double lengths[]{.16, .18, .20, .16, .18, .20};
    for (std::size_t i = 0; i < grid.size(); ++i) {
        require(std::abs(grid[i].values[0] - lengths[i]) < 1e-14 &&
                    grid[i].values[1] == (i < 3 ? .5 : .7),
                "Raster values or deterministic order changed.");
    }
    const auto descending = makeSweepGrid(std::vector<SweepAxis>{{"length", .2, .1, 2}});
    require(descending[0].values[0] == .2 && descending[1].values[0] == .1,
            "Descending endpoints were reordered.");
    const auto single = makeSweepGrid(std::vector<SweepAxis>{{"length", .18, .18, 1}});
    require(single.size() == 1 && single[0].values[0] == .18,
            "Single-point studies are incorrect.");
    DraftParameters parameters;
    std::set<QString> seen;
    for (const auto &p : sweepParameters()) {
        require(seen.insert(p.key).second && !p.label.isEmpty() && std::isfinite(p.displayScale) &&
                    p.displayScale > 0,
                "A parameter descriptor is ambiguous or invalid.");
        if (p.key == "meshTargetEdgeLength")
            continue;
        const auto value = sweepParameterValue(parameters, p.key);
        setSweepParameterValue(parameters, p.key, value);
        require(sweepParameterValue(parameters, p.key) == value,
                "A parameter binding changed an unchanged value.");
    }
    setSweepParameterValue(parameters, "verticalProfile.k", .73);
    require(parameters.verticalProfile.k == .73 && parameters.profile.k == 1,
            "A vertical parameter modified its horizontal counterpart.");
    setSweepParameterValue(parameters, "leftProfile.k", .65);
    setSweepParameterValue(parameters, "bottomProfile.terminationScale", .21);
    setSweepParameterValue(parameters, "leftAngle", .4);
    setSweepParameterValue(parameters, "bottomAngle", .3);
    require(parameters.leftProfile.k == .65 && parameters.profile.k == 1 &&
                parameters.bottomProfile.terminationScale == .21 &&
                parameters.verticalProfile.terminationScale != .21 && parameters.leftAngle == .4 &&
                parameters.bottomAngle == .3 && parameters.horizontalAngle != .4 &&
                parameters.verticalAngle != .3,
            "Independent side sweep bindings modified another direction.");
    rejects([&] { setSweepParameterValue(parameters, "unknown", 1); });
    rejects([&] {
        setSweepParameterValue(parameters, "length", std::numeric_limits<double>::quiet_NaN());
    });
    for (const auto invalid : std::vector<std::vector<SweepAxis>>{
             {},
             {{"length", .1, .2, 0}},
             {{"length", .1, .2, 1}},
             {{"unknown", 1, 2, 2}},
             {{"length", .1, .2, 2}, {"length", .1, .2, 2}},
             {{"length", .1, .2, 2}, {"horizontalAngle", .3, .4, 2}, {"verticalAngle", .4, .5, 2}},
             {{"length", .1, .2, 1000000000}, {"horizontalAngle", .3, .4, 1000000000}}})
        rejects([&] { (void)makeSweepGrid(invalid); });
}
void cadConvergence(const QString &directory) {
    // A synthetic closed octahedron verifies snapshot/BC persistence, not acoustics.
    auto cad = std::make_shared<CadSnapshot>();
    cad->documentId = "study-test-document";
    cad->sessionId = "study-test-session";
    cad->revision = 17;
    cad->mesh = testing::sphere(0);
    cad->acousticOrigin = {.2, .3, .4};
    cad->acousticAxis = {0, 1, 0};
    cad->objects.push_back({"solid",
                            "Solid",
                            "Test solid",
                            "octahedron",
                            {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}});
    CadFace face;
    face.id = 1;
    face.objectId = "solid";
    face.sourceKey = "synthetic-surface";
    face.signature = {"mesh", "octahedron", "eight-triangles", .1, {}, {}, 12, 1};
    face.boundary = BoundaryCondition::driven;
    face.normalVelocity = {.4, -.7};
    face.state = CadFaceState::confirmed;
    cad->faces.push_back(face);
    Project baseline;
    baseline.source = ProjectSource::freecad;
    baseline.cad = cad;
    baseline.cadFaceTargetEdgeLengths = {{1, .002}};
    auto study = createSweepStudy(directory, baseline, {{"meshTargetEdgeLength", .008, .004, 2}},
                                  {100, 333, 1700});
    const auto loaded = readSweepStudy(sweepManifestPath(study));
    const auto fine = sweepVariantInputs(loaded, 1);
    auto solvedInputs = fine;
    solvedInputs.simulation.settings.minimumFrequency = 100;
    solvedInputs.simulation.settings.maximumFrequency = 1700;
    solvedInputs.simulation.settings.frequencyCount = 3;
    require(sweepVariantMatches(loaded, 1, solvedInputs),
            "CAD convergence did not preserve its offline baseline inputs.");
    require(fine.source == ProjectSource::freecad && fine.cad && fine.cad->revision == 17 &&
                fine.cad->mesh.triangles.size() == 8 && fine.cad->acousticOrigin.x == .2 &&
                fine.cad->acousticAxis.y == 1 &&
                fine.cad->faces[0].normalVelocity == std::complex<double>(.4, -.7) &&
                fine.cad->faces[0].state == CadFaceState::confirmed &&
                fine.cadFaceTargetEdgeLengths == baseline.cadFaceTargetEdgeLengths &&
                fine.meshing.targetEdgeLength == .004,
            "CAD convergence changed the offline geometry, placement or boundary assignment.");
    baseline.cadFaceTargetEdgeLengths[1] = .006;
    rejects([&] {
        (void)createSweepStudy(directory, baseline, {{"meshTargetEdgeLength", .008, .004, 2}},
                               {100});
    });
    baseline.cadFaceTargetEdgeLengths.clear();
    rejects([&] { (void)createSweepStudy(directory, baseline, {{"length", .1, .2, 2}}); });
}
void filesAndSnapshots() {
    QTemporaryDir directory;
    require(directory.isValid(), "Temporary directory failed.");
    const auto untouched = directory.filePath("user-file.txt");
    save(untouched, "keep");
    Project baseline;
    baseline.draftFaceTargetEdgeLengths = {{1, .003}};
    baseline.simulation.settings.minimumFrequency = 100;
    baseline.simulation.settings.maximumFrequency = 2000;
    baseline.simulation.settings.frequencyCount = 7;
    baseline.simulation.settings.acceleration = BemAcceleration::fastMultipole;
    baseline.simulation.settings.fastMultipole.leafSupportCapacity = 64;
    baseline.simulation.settings.beamwidthTargetEnabled = true;
    const auto reference = std::make_shared<const Project>(baseline);
    baseline.comparisonSnapshots = {{QStringLiteral("Referenz"), reference}};
    const std::vector<SweepAxis> axes{{"length", .16, .20, 3}, {"horizontalAngle", .5, .7, 2}};
    auto study = createSweepStudy(directory.path(), baseline, axes);
    require(baseline.comparisonSnapshots.size() == 1 &&
                baseline.comparisonSnapshots.front().project == reference &&
                study.baseline.comparisonSnapshots.empty(),
            "Study copied or modified the source project's comparison snapshots.");
    const auto path = sweepManifestPath(study);
    require(QFileInfo::exists(path) &&
                QFileInfo::exists(QDir(study.directory).filePath("baseline.hsproj")),
            "The manifest or independent input snapshot is missing.");
    const auto variant = sweepVariantInputs(study, 4);
    require(sweepVariantMatches(study, 4, variant), "Valid variant inputs did not match.");
    auto wrongInput = variant;
    wrongInput.parameters.verticalAngle += .01;
    require(!sweepVariantMatches(study, 4, wrongInput), "Unrelated geometry mismatch was accepted.");
    wrongInput = variant;
    wrongInput.simulation.settings.velocityScale *= 2;
    require(!sweepVariantMatches(study, 4, wrongInput), "Changed excitation was accepted.");
    require(variant.comparisonSnapshots.empty() &&
                readProject(QDir(study.directory).filePath("baseline.hsproj")).comparisonSnapshots.empty(),
            "Historical comparison snapshots leaked into the study or variant.");
    require(std::abs(variant.parameters.length - .18) < 1e-14 &&
                variant.parameters.horizontalAngle == .7,
            "A variant was not derived from its two raster coordinates.");
    require(variant.parameters.verticalAngle == baseline.parameters.verticalAngle &&
                variant.draftFaceTargetEdgeLengths == baseline.draftFaceTargetEdgeLengths &&
                variant.simulation.settings.frequencyCount == 7 &&
                variant.simulation.settings.acceleration == BemAcceleration::fastMultipole &&
                variant.simulation.settings.fastMultipole.leafSupportCapacity == 64 &&
                variant.simulation.settings.beamwidthTargetEnabled,
            "Common mesh, solver or unrelated geometry settings were changed.");
    require(!variant.mesh && !variant.preview && !variant.simulation.result &&
                !variant.simulation.current && !variant.geometryRevision && !variant.meshRevision,
            "A parameter variant inherited results for another geometry.");
    study.variants[1].state = SweepState::failed;
    study.variants[1].error = "Test fixture: deliberately invalid variant";
    study.variants[2].state = SweepState::cancelled;
    writeSweepStudy(study);
    const auto loaded = readSweepStudy(path);
    require(loaded.axes.size() == 2 && loaded.variants.size() == 6 &&
                loaded.variants[1].state == SweepState::failed &&
                loaded.variants[1].error == study.variants[1].error &&
                loaded.variants[2].state == SweepState::cancelled,
            "Per-variant failure or cancellation state did not round-trip.");
    const auto original = contents(path);
    auto corrupt = study;
    corrupt.variants[0].projectFile = "../user-file.txt";
    rejects([&] { writeSweepStudy(corrupt); });
    require(contents(path) == original && contents(untouched) == "keep",
            "A rejected study write changed existing files.");
    auto missing = study;
    missing.variants[0].state = SweepState::completed;
    writeSweepStudy(missing);
    rejects([&] { (void)readSweepStudy(path); });
    writeSweepStudy(study);
    auto json = QJsonDocument::fromJson(original).object();
    json.insert("baseline", "../user-file.txt");
    save(path, QJsonDocument(json).toJson());
    rejects([&] { (void)readSweepStudy(path); });
    writeSweepStudy(study);
    const auto another = createSweepStudy(directory.path(), baseline, axes);
    require(another.directory != study.directory && contents(untouched) == "keep",
            "Starting another study reused user files.");
    const auto convergence = createSweepStudy(
        directory.path(), baseline, {{"meshTargetEdgeLength", .008, .004, 3}}, {100, 500, 1000});
    const auto fine = sweepVariantInputs(convergence, 2);
    auto completedInputs = fine;
    completedInputs.simulation.settings.minimumFrequency = 100;
    completedInputs.simulation.settings.maximumFrequency = 1000;
    completedInputs.simulation.settings.frequencyCount = 3;
    require(sweepVariantMatches(convergence, 2, completedInputs),
            "Convergence frequency override did not match its variant.");
    require(fine.meshing.targetEdgeLength == .004 &&
                std::abs(sweepVariantInputs(convergence, 1).meshing.targetEdgeLength - .006) < 1e-14 &&
                fine.draftFaceTargetEdgeLengths == baseline.draftFaceTargetEdgeLengths,
            "Convergence changed local face settings or the fine target.");
    require(readSweepStudy(sweepManifestPath(convergence)).frequencyOverride ==
                std::vector<double>({100, 500, 1000}),
            "Selected convergence frequencies were replaced by a full logarithmic sweep.");
    for (const auto axis : std::vector<SweepAxis>{{"meshTargetEdgeLength", .008, .008, 1},
                                                  {"meshTargetEdgeLength", .004, .008, 2},
                                                  {"length", .1, .2, 2}}) {
        rejects([&] { (void)createSweepStudy(directory.path(), baseline, {axis}, {100}); });
        auto invalid = convergence;
        invalid.axes = {axis};
        invalid.variants = makeSweepGrid(invalid.axes);
        rejects([&] { (void)sweepVariantInputs(invalid, 0); });
    }
    const auto convergencePath = sweepManifestPath(convergence);
    const auto convergenceBytes = contents(convergencePath);
    auto invalidWrite = convergence;
    invalidWrite.axes[0].parameter = "length";
    invalidWrite.variants = makeSweepGrid(invalidWrite.axes);
    rejects([&] { writeSweepStudy(invalidWrite); });
    require(contents(convergencePath) == convergenceBytes,
            "Rejected convergence write changed its saved manifest.");
    auto convergenceJson = QJsonDocument::fromJson(convergenceBytes).object();
    auto convergenceAxes = convergenceJson["axes"].toArray();
    auto alteredAxis = convergenceAxes[0].toObject();
    alteredAxis["parameter"] = "length";
    convergenceAxes[0] = alteredAxis;
    convergenceJson["axes"] = convergenceAxes;
    save(convergencePath, QJsonDocument(convergenceJson).toJson());
    rejects([&] { (void)readSweepStudy(convergencePath); });
    save(convergencePath, convergenceBytes);
    rejects([&] {
        (void)createSweepStudy(directory.path(), baseline,
                               {{"meshTargetEdgeLength", .008, .004, 2}}, {100, 100});
    });
    rejects(
        [&] { (void)createSweepStudy(directory.path(), baseline, {{"guide.width", .2, .3, 2}}); });
    auto inactiveConflict = baseline;
    inactiveConflict.cadFaceTargetEdgeLengths = {{37, .006}};
    rejects([&] {
        (void)createSweepStudy(directory.path(), inactiveConflict,
                               {{"meshTargetEdgeLength", .008, .004, 2}}, {100});
    });
    inactiveConflict.cadFaceTargetEdgeLengths[37] = .00001;
    rejects([&] { validateSweepMeshSettings(inactiveConflict); });
    inactiveConflict.cadFaceTargetEdgeLengths[37] = std::numeric_limits<double>::infinity();
    rejects([&] { validateSweepMeshSettings(inactiveConflict); });
    cadConvergence(directory.path());
}
void asymmetricSnapshots() {
    QTemporaryDir directory;
    require(directory.isValid(), "Temporary directory failed.");
    Project baseline;
    const auto symmetric = sweepParameters(baseline.parameters);
    for (const auto &parameter : symmetric) {
        require(parameter.key != "leftAngle" && parameter.key != "bottomAngle" &&
                    !parameter.key.startsWith("leftProfile.") &&
                    !parameter.key.startsWith("bottomProfile.") &&
                    !parameter.label.contains("rechts") && !parameter.label.contains("oben"),
                "Symmetric catalog exposes individual sides or directional labels.");
    }
    const std::vector<SweepAxis> axes{{"leftAngle", .3, .4, 2}, {"bottomProfile.k", .75, .9, 2}};
    rejects([&] { (void)createSweepStudy(directory.path(), baseline, axes); });
    SweepStudy invalidStudy;
    invalidStudy.baseline = baseline;
    invalidStudy.axes = axes;
    invalidStudy.variants = makeSweepGrid(axes);
    rejects([&] { (void)sweepVariantInputs(invalidStudy, 0); });
    const auto symmetricStudy = createSweepStudy(directory.path(), baseline,
        {{"horizontalAngle", .4, .5, 2}, {"verticalAngle", .3, .35, 2}});
    const auto symmetricRestored = readSweepStudy(sweepManifestPath(symmetricStudy));
    const auto symmetricVariant = sweepVariantInputs(symmetricRestored, 3);
    const DraftGeometry symmetricGeometry(symmetricVariant.parameters);
    const double z = symmetricVariant.parameters.length;
    constexpr double pi = std::numbers::pi;
    require(std::abs(symmetricGeometry.radius(z, 0) - symmetricGeometry.radius(z, pi)) < 1e-10 &&
                std::abs(symmetricGeometry.radius(z, pi / 2) -
                         symmetricGeometry.radius(z, 3 * pi / 2)) < 1e-10 &&
                !sweepVariantLabel(symmetricRestored, 3).contains("rechts") &&
                !sweepVariantLabel(symmetricRestored, 3).contains("oben"),
            "Symmetric study changed opposing geometry or used directional labels.");
    baseline.parameters.asymmetricEnabled = true;
    const auto asymmetric = sweepParameters(baseline.parameters);
    require(std::any_of(asymmetric.begin(), asymmetric.end(), [](const auto &p) {
                return p.key == "leftAngle" && p.label.contains("links");
            }) && std::any_of(asymmetric.begin(), asymmetric.end(), [](const auto &p) {
                return p.key == "horizontalAngle" && p.label.contains("rechts");
            }) && std::any_of(asymmetric.begin(), asymmetric.end(), [](const auto &p) {
                return p.key == "verticalAngle" && p.label.contains("oben");
            }) && std::any_of(asymmetric.begin(), asymmetric.end(), [](const auto &p) {
                return p.key == "bottomAngle" && p.label.contains("unten");
            }), "Asymmetric catalog did not expose all four named directions.");
    baseline.parameters.leftProfile.terminationQ = .94;
    baseline.parameters.bottomAngle = .35;
    const auto study = createSweepStudy(directory.path(), baseline, axes);
    const auto restored = readSweepStudy(sweepManifestPath(study));
    const auto variant = sweepVariantInputs(restored, 3);
    require(variant.parameters.asymmetricEnabled && variant.parameters.leftAngle == .4 &&
                variant.parameters.bottomProfile.k == .9 &&
                variant.parameters.leftProfile.terminationQ == .94 &&
                variant.parameters.bottomAngle == .35 &&
                variant.parameters.horizontalAngle == baseline.parameters.horizontalAngle &&
                variant.parameters.verticalProfile.k == baseline.parameters.verticalProfile.k,
            "Four-direction study did not preserve independent baseline and variant inputs.");
    const DraftGeometry geometry(variant.parameters);
    require(std::abs(geometry.radius(variant.parameters.length, 0) -
                     geometry.radius(variant.parameters.length, pi)) > 1e-8 &&
                std::abs(geometry.radius(variant.parameters.length, pi / 2) -
                         geometry.radius(variant.parameters.length, 3 * pi / 2)) > 1e-8,
            "Asymmetric variant did not produce independently shaped sides.");
}
void slotSnapshots() {
    QTemporaryDir directory;
    require(directory.isValid(), "Temporary directory failed.");
    Project baseline;
    const std::vector<SweepAxis> axes{{"slot.aperture.width", .01, .014, 2},
                                      {"slot.transitionLength", 0, .006, 2}};
    rejects([&] { (void)createSweepStudy(directory.path(), baseline, axes); });
    baseline.parameters.slot.enabled = true;
    baseline.parameters.slot.aperture.height = .037123456789;
    baseline.parameters.transition.startFraction = .3;
    const auto study = createSweepStudy(directory.path(), baseline, axes);
    const auto restored = readSweepStudy(sweepManifestPath(study));
    for (std::size_t index = 0; index < restored.variants.size(); ++index) {
        const auto variant = sweepVariantInputs(restored, index);
        require(variant.parameters.slot.enabled &&
                    variant.parameters.slot.aperture.width == (index % 2 == 0 ? .01 : .014) &&
                    variant.parameters.slot.transitionLength == (index < 2 ? 0 : .006) &&
                    variant.parameters.slot.aperture.height ==
                        baseline.parameters.slot.aperture.height &&
                    variant.parameters.slot.adapterLength ==
                        baseline.parameters.slot.adapterLength &&
                    variant.parameters.entry.width == baseline.parameters.entry.width &&
                    variant.parameters.length == baseline.parameters.length,
                "Slot studies lost independent inputs or changed the driver entry/total length.");
        (void)DraftGeometry(variant.parameters);
    }
}
void physicalSlotRadiusSnapshots() {
    QTemporaryDir directory;
    require(directory.isValid(), "Temporary directory failed.");
    Project baseline;
    auto &p = baseline.parameters;
    p.slot.enabled = true;
    p.transition.startFraction = .25;
    p.transition.endFraction = .65;
    rejects([&] {
        (void)createSweepStudy(directory.path(), baseline, {{"slot.cornerRadius", .001, .003, 2}});
    });
    rejects([&] {
        (void)createSweepStudy(directory.path(), baseline, {{"slot.exitRadius", .001, .003, 2}});
    });
    p.slot.contour = SlotContour::roundedRectangle;
    p.slot.cornerRadius = .002;
    p.slot.exit = SlotExit::circularFillet;
    p.slot.exitRadius = .002;
    p.slot.straightLength = .010;
    p.slot.transitionLength = .007123456789123;
    p.profile.k = p.verticalProfile.k = 0;
    p.profile.terminationScale = p.verticalProfile.terminationScale = 0;
    p.horizontalAngle = .5;
    p.verticalAngle = .3;
    rejects([&] {
        (void)createSweepStudy(directory.path(), baseline,
                               {{"slot.aperture.exponent", 2, 4, 2}});
    });
    rejects([&] {
        (void)createSweepStudy(directory.path(), baseline,
                               {{"slot.transitionLength", .001, .003, 2}});
    });
    const std::vector<SweepAxis> axes{{"slot.cornerRadius", .001, .003, 2},
                                      {"slot.exitRadius", .001, .003, 2}};
    const auto study = createSweepStudy(directory.path(), baseline, axes);
    const auto restored = readSweepStudy(sweepManifestPath(study));
    for (std::size_t index = 0; index < restored.variants.size(); ++index) {
        const auto variant = sweepVariantInputs(restored, index);
        const auto &slot = variant.parameters.slot;
        require(slot.contour == SlotContour::roundedRectangle &&
                    slot.exit == SlotExit::circularFillet &&
                    slot.cornerRadius == (index % 2 == 0 ? .001 : .003) &&
                    slot.exitRadius == (index < 2 ? .001 : .003) &&
                    slot.transitionLength == p.slot.transitionLength &&
                    slot.aperture.exponent == p.slot.aperture.exponent &&
                    variant.parameters.profile.k == 0 &&
                    variant.parameters.horizontalAngle == p.horizontalAngle,
                "Radius studies lost their independent dimensions, modes or dormant settings.");
        (void)DraftGeometry(variant.parameters);
    }
}
} // namespace
int main(int argc, char **argv) {
    QCoreApplication application(argc, argv);
    try {
        gridAndBindings();
        filesAndSnapshots();
        asymmetricSnapshots();
        slotSnapshots();
        physicalSlotRadiusSnapshots();
        std::cout
            << "Parameter raster, immutable inputs and atomic study persistence tests passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
