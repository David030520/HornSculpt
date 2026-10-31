#include "acoustic_partition_fixture.hpp"
#include "gmsh_source.hpp"
#include <QCoreApplication>
#include <QJsonArray>
#include <QJsonDocument>
#include <iostream>
#include <sstream>
#include <stdexcept>
using namespace hornsculpt;
void require(bool v, const char *m) {
    if (!v)
        throw std::runtime_error(m);
}
int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    try {
        std::istringstream in("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n1\n2 7 "
                              "\"Mund\"\n$EndPhysicalNames\n"
                              "$Nodes\n3\n10 0 0 0\n20 1000 0 0\n30 0 1000 0\n$EndNodes\n"
                              "$Elements\n1\n40 2 2 7 5 10 20 30\n$EndElements\n");
        MeshImportSelection selection{importGmsh(in, {.001}), .001};
        auto snapshot = makeGmshCadSnapshot(selection, QStringLiteral("mesh.msh"));
        require(snapshot->sessionId == "gmsh-import-offline" && snapshot->faces.size() == 1 &&
                    snapshot->faces[0].label.contains("Mund") &&
                    snapshot->faces[0].signature.area == .5,
                "Polygon snapshot data");
        const auto encoded = encodeGmshMetadata(selection);
        const auto decoded = decodeGmshMetadata(encoded, snapshot->mesh);
        require(decoded.lengthScaleToMetres == .001 &&
                    decoded.imported.nodeTags == selection.imported.nodeTags &&
                    decoded.imported.triangleMetadata[0].elementTag == 40 &&
                    decoded.imported.triangleMetadata[0].physicalTags ==
                        std::vector<std::uint32_t>{7} &&
                    decoded.imported.physicalGroups[0].name == "Mund",
                "Metadata roundtrip");
        auto broken = encoded;
        broken["nodes"] = QJsonArray{};
        bool rejected = false;
        try {
            (void)decodeGmshMetadata(broken, snapshot->mesh);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Invalid metadata accepted");
        broken = encoded;
        auto triangles = broken["triangles"].toArray();
        auto triangleJson = triangles[0].toObject();
        triangleJson["physical"] = QJsonArray{QStringLiteral("7"), QStringLiteral("7")};
        triangles[0] = triangleJson;
        broken["triangles"] = triangles;
        rejected = false;
        try {
            (void)decodeGmshMetadata(broken, snapshot->mesh);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Duplicate physical tags accepted");
        auto cup = testing::acousticCup();
        for (auto &triangle : cup.triangles)
            triangle.boundary = BoundaryCondition::wall;
        const auto center = static_cast<std::uint32_t>(cup.vertices.size());
        cup.vertices.push_back({0, 0, .08});
        for (std::uint32_t i = 0; i < 16; ++i)
            cup.triangles.push_back(
                {{center, 64 + (i + 1) % 16, 64 + i}, 6, BoundaryCondition::wall});
        MeshingOptions options;
        options.targetEdgeLength = .02;
        options.quality.enabled = false;
        const std::vector<AcousticRegionDefinition> regions{{7, "Innen", {2, 3}, {6}}};
        const auto checked = validateGmshMesh(cup, regions, options);
        require(checked.mesh == cup && !checked.quality.watertight(),
                "Explicit Gmsh interface did not preserve the source mesh.");
        rejected = false;
        try {
            (void)validateGmshMesh(cup, {}, options);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "A three-sided source was accepted as a single exterior domain.");
        std::cout << "Gmsh source tests PASS.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
