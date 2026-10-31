#include "project_assembly.hpp"
#include "project_json.hpp"
#include "step_import.hpp"
#include <QCoreApplication>
#include <QTemporaryDir>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
template <class F> void rejects(F action, const char *message) {
    bool rejected = false;
    try {
        action();
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, message);
}
} // namespace
int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    try {
        Project p;
        p.parameters.length = .04;
        p.parameters.asymmetricEnabled = true;
        p.parameters.leftAngle = .3;
        p.meshing.targetEdgeLength = .008;
        p.meshing.adaptive.geometricTolerance = .001;
        p.meshing.quality.enabled = false;
        DraftGeometry g(p.parameters);
        const auto local = g.previewMesh(.001);
        p.assembly.instances[0].position = {-.3, 0, 0};
        auto second = p.assembly.instances[0];
        second.id = 7;
        second.name = "Right";
        second.position = {.3, 0, 0};
        second.mirrorHorizontal = true;
        second.mouthAnchored = true;
        p.assembly.instances.push_back(second);
        auto composed = assemblyPreview(local, p);
        const auto q = inspectMesh(composed), single = inspectMesh(local);
        require(q.watertight() && q.selfIntersectionPairs == 0 && q.connectedComponents == 2,
                "Mirroring must preserve closed outward components");
        require(std::abs(q.signedVolume - 2 * single.signedVolume) < 1e-10,
                "Mirroring changed volume/sign");
        require(composed.triangles[local.triangles.size()].faceId ==
                    instanceFaceId(7, local.triangles[0].faceId),
                "Instance IDs are not stable");
        require(instancePoint(second, {0, 0, .04}, .04) == second.position, "Mouth anchor moved");
        require(instancePoint(second, {0, 0, .08}, .08) == second.position,
                "Length edit moved mouth anchor");
        auto bad = p.assembly;
        bad.instances[1].id = 1;
        rejects([&] { validateDraftAssembly(bad); }, "Duplicate instance ID accepted");
        bad = p.assembly;
        bad.instances[1].axes[0].x = 2;
        rejects([&] { validateDraftAssembly(bad); }, "Scaled instance accepted");
        p.preview = std::make_shared<SurfaceMesh>(composed);
        p.geometryRevision = p.revision;
        QTemporaryDir temp;
        writeProject(temp.filePath("assembly.hsproj"), p);
        auto loaded = readProject(temp.filePath("assembly.hsproj"));
        require(loaded.assembly == p.assembly && *loaded.preview == composed,
                "Assembly archive lost placement");
        require(project_json::decodeAssembly(project_json::encodeAssembly(p.assembly)) ==
                    p.assembly,
                "Assembly JSON roundtrip");
        p.cad = readStepSnapshot(QStringLiteral(HS_TEST_CAD), .001);
        rejects([&] { (void)appendAssemblyCad(*p.cad, *p.cad); },
                "Overlapping duplicate CAD imports accepted");
        p.assembly.includeCad = true;
        p.assembly.cadTargetEdgeLength = .02;
        p.preview = std::make_shared<SurfaceMesh>(assemblyPreview(local, p));
        const auto result = meshDraftAssembly(p, local);
        require(result.quality.watertight() && result.quality.connectedComponents == 3,
                "Combined Draft/CAD mesh is not three closed components");
        bool cadWall = false, driveLeft = false, driveRight = false;
        for (auto t : result.mesh.triangles) {
            cadWall |= t.faceId >= assemblyCadFaceOffset && t.boundary == BoundaryCondition::wall;
            driveLeft |= t.faceId == 5 && t.boundary == BoundaryCondition::driven;
            driveRight |= t.faceId == 53 && t.boundary == BoundaryCondition::driven;
        }
        require(cadWall && driveLeft && driveRight, "Solver source identities lost");
        p.mesh = std::make_shared<MeshingResult>(result);
        p.meshRevision = p.revision;
        writeProject(temp.filePath("mixed.hsproj"), p);
        loaded = readProject(temp.filePath("mixed.hsproj"));
        require(loaded.cad && loaded.assembly == p.assembly && loaded.mesh->mesh == result.mesh,
                "Mixed archive roundtrip");
        p.assembly.instances[1] = p.assembly.instances[0];
        p.assembly.instances[1].id = 7;
        rejects([&] { (void)meshDraftAssembly(p, local); }, "Overlapping horn instances accepted");
        std::cout
            << "Linked instances, reflection, stable IDs, CAD meshing and persistence passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
