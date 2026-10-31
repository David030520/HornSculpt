#include <cmath>
#include <hornsculpt/gmsh_mesh.hpp>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace hornsculpt;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
const std::string msh22 = R"($MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
1
2 7 "Driven face"
$EndPhysicalNames
$Nodes
4
10 0 0 0
20 1000 0 0
30 0 1000 0
40 0 0 1000
$EndNodes
$Elements
4
21 2 2 7 4 10 30 20
22 2 2 7 4 10 20 40
23 2 2 7 5 20 30 40
24 2 2 0 6 30 10 40
$EndElements
)";
const std::string msh41 = R"($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
2
2 7 "first group"
2 8 "other group"
$EndPhysicalNames
$Entities
0 0 1 0
4 0 0 0 1000 1000 0 2 7 8 0
$EndEntities
$Nodes
1 3 10 30
2 4 0 3
10 20 30
0 0 0
1000 0 0
0 1000 0
$EndNodes
$Elements
1 1 41 41
2 4 2 1
41 10 20 30
$EndElements
)";
GmshSurfaceImport load(const std::string &s, double scale = .001) {
    std::istringstream in(s);
    return importGmsh(in, {scale});
}
void replace(std::string &s, const std::string &a, const std::string &b) {
    auto pos = s.find(a);
    require(pos != std::string::npos, "fixture replacement missing");
    s.replace(pos, a.size(), b);
}
void reject(const std::string &s, double scale = .001) {
    bool failed = false;
    try {
        (void)load(s, scale);
    } catch (const std::invalid_argument &) {
        failed = true;
    }
    require(failed, "Invalid input was accepted");
}
void v22() {
    auto x = load(msh22);
    require(x.mesh.vertices.size() == 4 && x.mesh.triangles.size() == 4 &&
                x.triangleMetadata.size() == 4,
            "2.2 dimensions");
    require(x.mesh.vertices[1].x == 1. &&
                x.mesh.triangles[0].vertices == std::array<std::uint32_t, 3>{0, 2, 1},
            "scale or winding");
    require(x.mesh.triangles[2].faceId == 5 && x.triangleMetadata[2].elementTag == 23 &&
                x.triangleMetadata[2].physicalTags == std::vector<std::uint32_t>{7} &&
                x.triangleMetadata[3].physicalTags.empty(),
            "2.2 tags");
    require(x.physicalGroups.size() == 1 && x.physicalGroups[0].name == "Driven face",
            "2.2 physical name");
    require(x.diagnostics.boundaryEdges == 0 && x.diagnostics.nonManifoldEdges == 0 &&
                x.diagnostics.inconsistentEdges == 0 && x.diagnostics.degenerateTriangles == 0,
            "closed tetra diagnostics");
}
void v41() {
    auto x = load(msh41);
    require(x.mesh.triangles.size() == 1 && x.mesh.triangles[0].faceId == 4 &&
                x.triangleMetadata[0].elementTag == 41 &&
                x.triangleMetadata[0].physicalTags == std::vector<std::uint32_t>({7, 8}),
            "4.1 entity/group mapping");
    require(x.physicalGroups.size() == 2 && x.physicalGroups[1].name == "other group" &&
                x.diagnostics.boundaryEdges == 3,
            "4.1 names/open surface");
}
void mixedElements() {
    auto s = msh22;
    replace(s, "$Nodes\n4", "$Nodes\n5");
    replace(s, "$EndNodes", "99 2 2 2\n$EndNodes");
    replace(s, "$Elements\n4", "$Elements\n6");
    replace(s, "$EndElements", "25 15 2 0 8 99\n26 1 2 0 9 10 99\n$EndElements");
    auto x = load(s);
    require(x.mesh.triangles.size() == 4 && x.diagnostics.skippedPoints == 1 &&
                x.diagnostics.skippedLines == 1 && x.mesh.vertices.size() == 4 &&
                x.nodeTags == std::vector<std::uint64_t>({10, 20, 30, 40}),
            "2.2 auxiliary element counts and node compaction");
    s = msh41;
    replace(s, "0 0 1 0\n4", "1 1 1 0\n1 0 0 0 0\n2 0 0 0 1000 0 0 0 0\n4");
    replace(s, "1 1 41 41\n2 4 2 1\n41 10 20 30",
            "3 3 41 43\n0 1 15 1\n42 10\n1 2 1 1\n43 10 20\n2 4 2 1\n41 10 20 30");
    x = load(s);
    require(x.mesh.triangles.size() == 1 && x.diagnostics.skippedPoints == 1 &&
                x.diagnostics.skippedLines == 1,
            "4.1 auxiliary entity blocks");
}
void quads() {
    auto s = msh41;
    replace(s, "1 3 10 30", "1 4 10 40");
    replace(s, "2 4 0 3\n10 20 30\n0 0 0\n1000 0 0\n0 1000 0",
            "2 4 0 4\n10 20 30 40\n0 0 0\n1000 0 0\n1000 1000 0\n0 1000 0");
    replace(s, "2 4 2 1\n41 10 20 30", "2 4 3 1\n41 10 20 30 40");
    auto x = load(s);
    require(x.mesh.triangles.size() == 2 && x.triangleMetadata[0].elementTag == 41 &&
                x.triangleMetadata[1].elementTag == 41 &&
                x.mesh.triangles[0].vertices == std::array<std::uint32_t, 3>{0, 1, 2} &&
                x.mesh.triangles[1].vertices == std::array<std::uint32_t, 3>{0, 2, 3},
            "quad source and orientation");
    auto bad = s;
    replace(bad, "0 1000 0\n$EndNodes", "0 1000 1\n$EndNodes");
    reject(bad);
    bad = s;
    replace(bad, "1000 1000 0\n0 1000 0", "200 200 0\n0 1000 0");
    reject(bad);
}
void cancelledImport() {
    std::stop_source source;
    source.request_stop();
    std::istringstream in(msh22);
    bool cancelled = false;
    try {
        (void)importGmsh(in, {.001, source.get_token()});
    } catch (const std::runtime_error &) {
        cancelled = true;
    }
    require(cancelled, "Cancelled import was not stopped");
}
void invalid() {
    auto s = msh22;
    replace(s, "24 2 2 0 6", "24 2 2 0 5");
    reject(s);
    s = msh22;
    replace(s, "10 0 0 0", "-1 0 0 0");
    reject(s);
    s = msh22;
    replace(s, "21 2 2 7 4", "-1 2 2 7 4");
    reject(s);
    s = msh22;
    replace(s, "$Nodes\n4", "$Nodes\n-1");
    reject(s);
    s = msh22;
    replace(s, "21 2 2 7 4", "21 2 2000001 7 4");
    reject(s);
    s = msh22;
    replace(s, "2.2 0 8", "2.2 1 8");
    reject(s);
    s = msh22;
    replace(s, "2 7 \"Driven face\"", "2 7 Driven");
    reject(s);
    s = msh22;
    replace(s, "21 2 2 7 4 10 30 20", "21 2 3 7 4 1 10 30 20");
    reject(s);
    s = msh41;
    replace(s, "2 7 8 0", "2 7 7 0");
    reject(s);
    s = msh22;
    replace(s, "2.2 0 8", "4.0 0 8");
    reject(s);
    s = msh22;
    replace(s, "21 2 2 7 4 10 30 20", "21 9 2 7 4 10 30 20");
    reject(s);
    s = msh22;
    replace(s, "21 2 2 7 4 10 30 20", "21 2 2 7 4 10 30 999");
    reject(s);
    s = msh22;
    replace(s, "20 1000 0 0", "10 1000 0 0");
    reject(s);
    s = msh22;
    replace(s, "22 2 2 7 4", "21 2 2 7 4");
    reject(s);
    s = msh41;
    replace(s, "2 4 2 1", "2 999 2 1");
    reject(s);
    s = msh41;
    replace(s, "1 3 10 30", "1 3 10 40");
    reject(s);
    s = msh41;
    replace(s, "41 10 20 30", "41 10 20");
    reject(s);
    s = msh41;
    s += "$Periodic\n";
    reject(s);
    reject(msh22, 0);
}
} // namespace
int main() {
    try {
        v22();
        cancelledImport();
        v41();
        mixedElements();
        quads();
        invalid();
        std::cout << "Gmsh import tests PASS.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
