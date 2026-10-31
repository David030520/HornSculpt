#include "fixture_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <hornsculpt/acoustic_partition.hpp>
#include <hornsculpt/symmetry_sector.hpp>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
template <class F> void rejected(F f) {
    try {
        f();
    } catch (const std::invalid_argument &) {
        return;
    }
    throw std::runtime_error("Invalid symmetry sector was accepted");
}
SurfaceMesh sector() {
    auto mesh = testing::sphere(0);
    mesh.triangles = {mesh.triangles[0], mesh.triangles[4]};
    // Remove unused nodes; keep the positive X/Y quarter sphere.
    mesh.vertices = {{.1, 0, 0}, {0, .1, 0}, {0, 0, .1}, {0, 0, -.1}};
    mesh.triangles = {{{2, 0, 1}, 1, BoundaryCondition::wall},
                      {{3, 1, 0}, 1, BoundaryCondition::wall}};
    return mesh;
}
OrthogonalSymmetryPlanes planes() {
    return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 1e-9};
}
void valid() {
    const auto mesh = sector();
    const auto old = mesh;
    auto report = validateSymmetrySector(mesh, planes());
    require(report.firstPlaneEdges == 2 && report.secondPlaneEdges == 2 &&
                !report.intersectionAxisEdges && mesh == old,
            "Valid sector edge mapping or immutability failed");
    auto rotated = mesh;
    for (auto &p : rotated.vertices) {
        const auto x = p.x;
        p.x = -p.y + 3;
        p.y = x + 4;
        p.z += 5;
    }
    auto transform = OrthogonalSymmetryPlanes{{3, 4, 5}, {0, 1, 0}, {-1, 0, 0}, 1e-9};
    auto actual = validateSymmetrySector(rotated, transform);
    require(actual.firstPlaneEdges == 2 && actual.secondPlaneEdges == 2,
            "Rotated plane mapping failed");
    auto near = mesh;
    near.vertices[1].x = -5e-10;
    report = validateSymmetrySector(near, planes());
    require(report.firstPlaneEdges == 2 && near.vertices[1].x == -5e-10,
            "Tolerance moved or rejected seam node");
}
void partition() {
    auto mesh = sector();
    const auto inner = mesh.triangles.size();
    auto outer = sector();
    const auto offset = static_cast<std::uint32_t>(mesh.vertices.size());
    for (auto v : outer.vertices)
        mesh.vertices.push_back({2 * v.x, 2 * v.y, 2 * v.z});
    for (auto t : outer.triangles) {
        for (auto &v : t.vertices)
            v += offset;
        std::swap(t.vertices[1], t.vertices[2]);
        t.faceId = 2;
        mesh.triangles.push_back(t);
    }
    MeshingOptions options;
    options.targetEdgeLength = .1;
    options.quality.enabled = false;
    auto symmetry = planes();
    const std::vector<AcousticRegionDefinition> regions{{7, "Shell", {1}, {2}}};
    auto result = partitionAcousticSurface(mesh, regions, options, {}, &symmetry);
    require(result.mesh == mesh && result.model.interfaces.size() == 1 &&
                result.model.interfaces.front().triangleIndices.size() == inner,
            "Symmetric explicit partition altered source or interface");
    validateSymmetryAcousticPartition(result.mesh, result.model, symmetry);
    AcousticModelDefinition explicitModel;
    explicitModel.domains = {{7, "Shell", BemDomainKind::interior, {1}, {}, false},
                             {99, "Outside", BemDomainKind::exterior, {}, {}, false}};
    explicitModel.interfaces = {{42, "Mouth", 7, 99, {2}, false}};
    const auto named = partitionAcousticSurface(mesh, explicitModel, options, {}, &symmetry);
    require(named.model.domains.back().id == 99 && named.model.interfaces[0].id == 42 &&
                named.model.interfaces[0].name == "Mouth",
            "Explicit domain or interface identity was lost");
    std::reverse(explicitModel.domains.begin(), explicitModel.domains.end());
    const auto reordered = partitionAcousticSurface(mesh, explicitModel, options, {}, &symmetry);
    require(reordered.model.interfaces[0].firstDomainId == 7 &&
                reordered.model.interfaces[0].secondDomainId == 99,
            "Interface changed when domains were reordered");
    auto reversedMesh = mesh;
    for (auto &triangle : reversedMesh.triangles)
        if (triangle.faceId == 2)
            std::swap(triangle.vertices[1], triangle.vertices[2]);
    explicitModel.interfaces[0].firstDomainId = 99;
    explicitModel.interfaces[0].secondDomainId = 7;
    const auto outsideFirst =
        partitionAcousticSurface(reversedMesh, explicitModel, options, {}, &symmetry);
    require(outsideFirst.model.interfaces[0].firstDomainId == 99 &&
                outsideFirst.model.interfaces[0].secondDomainId == 7,
            "Exterior first interface side was lost");
    auto flipped = result.model;
    for (auto &panel : flipped.domains[0].panels)
        if (panel.triangleIndex < inner)
            panel.normalSign = -panel.normalSign;
    rejected([&] { validateSymmetryAcousticPartition(result.mesh, flipped, symmetry); });
    auto bad = result.model;
    bad.interfaces[0].triangleIndices.pop_back();
    rejected([&] { validateSymmetryAcousticPartition(result.mesh, bad, symmetry); });
    rejected([&] {
        partitionAcousticSurface(mesh, std::vector<AcousticRegionDefinition>{{7, "Shell", {1}, {}}},
                                 options, {}, &symmetry);
    });
    auto hole = mesh;
    hole.triangles.erase(hole.triangles.begin());
    rejected([&] { partitionAcousticSurface(hole, regions, options, {}, &symmetry); });
}
void sharedWall() {
    auto mesh = sector();
    for (auto &triangle : mesh.triangles)
        std::swap(triangle.vertices[1], triangle.vertices[2]);
    mesh.triangles[1].faceId = 2;
    AcousticRegionDefinition region{7, "Inside", {}, {2}};
    region.twoSidedFaceIds = {1};
    MeshingOptions options;
    options.targetEdgeLength = .1;
    options.quality.enabled = false;
    const auto symmetry = planes();
    auto result = partitionAcousticSurface(mesh, std::vector<AcousticRegionDefinition>{region},
                                           options, {}, &symmetry);
    require(result.mesh == mesh && result.model.interfaces.size() == 1 &&
                result.model.twoSidedSurfaces.size() == 1 &&
                result.model.twoSidedSurfaces.front().triangleIndices ==
                    std::vector<std::uint32_t>{0},
            "Two-sided wall was not retained independently of interface");
    validateSymmetryAcousticPartition(result.mesh, result.model, symmetry);
    auto invalidModel = result.model;
    invalidModel.twoSidedSurfaces[0].triangleIndices = {1};
    rejected([&] { validateSymmetryAcousticPartition(result.mesh, invalidModel, symmetry); });
    invalidModel = result.model;
    for (auto &panel : invalidModel.domains[1].panels)
        if (panel.triangleIndex == 0)
            panel.normalSign = -panel.normalSign;
    rejected([&] { validateSymmetryAcousticPartition(result.mesh, invalidModel, symmetry); });
    auto driven = mesh;
    driven.triangles[0].boundary = BoundaryCondition::driven;
    const auto drivenResult = partitionAcousticSurface(
        driven, std::vector<AcousticRegionDefinition>{region}, options, {}, &symmetry);
    require(drivenResult.model.twoSidedSurfaces.size() == 1 &&
                drivenResult.mesh.triangles[0].boundary == BoundaryCondition::driven,
            "Explicit two-sided driven surface was not retained");
    region.interfaceFaceIds.push_back(1);
    rejected([&] {
        partitionAcousticSurface(mesh, std::vector<AcousticRegionDefinition>{region}, options, {},
                                 &symmetry);
    });
}
void baffle() {
    auto mesh = sector();
    mesh.triangles.resize(1);
    mesh.vertices.resize(3);
    const ExteriorBaffle plane{{0, 0, 0}, {0, 0, 1}, 1e-9};
    rejected([&] { validateSymmetrySector(mesh, planes()); });
    (void)validateSymmetrySector(mesh, planes(), &plane);
    auto bad = plane;
    bad.normal = {1, 0, 0};
    rejected([&] { validateSymmetrySector(mesh, planes(), &bad); });
    bad = plane;
    bad.origin.z = .01;
    rejected([&] { validateSymmetrySector(mesh, planes(), &bad); });
    mesh.vertices[0].z = .01;
    rejected([&] { validateSymmetrySector(mesh, planes(), &plane); });
}
void fullBaffle() {
    SurfaceMesh mesh;
    mesh.vertices = {{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}, {0, 0, 1}};
    mesh.triangles = {{{0, 1, 4}, 1, BoundaryCondition::wall},
                      {{1, 2, 4}, 1, BoundaryCondition::wall},
                      {{2, 3, 4}, 1, BoundaryCondition::wall},
                      {{3, 0, 4}, 1, BoundaryCondition::wall}};
    ExteriorBaffle plane{{0, 0, 0}, {0, 0, 1}, 1e-9};
    validateBaffleBoundary(mesh, plane);
    auto hole = mesh;
    hole.triangles.pop_back();
    rejected([&] { validateBaffleBoundary(hole, plane); });
    mesh.vertices[0].z = -.01;
    rejected([&] { validateBaffleBoundary(mesh, plane); });
}
void invalid() {
    auto mesh = sector();
    auto p = planes();
    p.secondNormal = {1, 0, 0};
    rejected([&] { validateSymmetrySector(mesh, p); });
    p = planes();
    p.firstNormal = {2, 0, 0};
    rejected([&] { validateSymmetrySector(mesh, p); });
    p = planes();
    p.planeDistanceTolerance = 0;
    rejected([&] { validateSymmetrySector(mesh, p); });
    p = planes();
    p.origin.x = std::numeric_limits<double>::quiet_NaN();
    rejected([&] { validateSymmetrySector(mesh, p); });
    auto outside = mesh;
    outside.vertices[0].y = -1e-5;
    rejected([&] { validateSymmetrySector(outside, planes()); });
    auto opening = mesh;
    opening.vertices[1].x = .01;
    rejected([&] { validateSymmetrySector(opening, planes()); });
    auto bad = mesh;
    std::swap(bad.triangles[1].vertices[1], bad.triangles[1].vertices[2]);
    rejected([&] { validateSymmetrySector(bad, planes()); });
    auto coplanar = mesh;
    coplanar.vertices[0] = {0, 0, 0};
    rejected([&] { validateSymmetrySector(coplanar, planes()); });
}
} // namespace
int main() {
    try {
        valid();
        invalid();
        baffle();
        fullBaffle();
        partition();
        sharedWall();
        std::cout << "Symmetry sector tests PASS.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
