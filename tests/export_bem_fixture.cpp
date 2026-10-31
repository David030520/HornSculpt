#include "fixture_geometry.hpp"
#include <QCommandLineParser>
#include <QCoreApplication>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <cmath>
#include <cstdio>
#include <hornsculpt/draft_geometry.hpp>
#include <hornsculpt/meshing.hpp>
#include <stdexcept>

using namespace hornsculpt;
namespace {
QJsonArray point(Vec3 p) {
    return {p.x, p.y, p.z};
}
QJsonArray observers(Vec3 origin) {
    QJsonArray positions;
    positions.append(point({origin.x, origin.y, origin.z + 1}));
    for (int plane = 0; plane < 2; ++plane)
        for (int degrees = -90; degrees <= 90; degrees += 5) {
            const double angle = degrees * std::acos(-1.0) / 180;
            positions.append(
                point({origin.x + (plane == 0 ? std::sin(angle) : 0),
                       origin.y + (plane == 1 ? std::sin(angle) : 0), origin.z + std::cos(angle)}));
        }
    if (positions.size() != 75)
        throw std::logic_error("Observer count mismatch.");
    for (const auto position : positions)
        if (!position.isArray() || position.toArray().size() != 3)
            throw std::logic_error("Each observer must have three coordinates.");
    return positions;
}
QJsonObject exportFixture(SurfaceMesh mesh, bool horn) {
    const auto quality = inspectMesh(mesh);
    if (!quality.watertight() || quality.signedVolume <= 0)
        throw std::runtime_error("Fixture geometry is invalid.");
    for (auto &triangle : mesh.triangles)
        if (triangle.boundary == BoundaryCondition::driven) {
            Vec3 center{};
            for (const auto vertex : triangle.vertices) {
                const auto p = mesh.vertices[vertex];
                center.x += p.x;
                center.z += p.z;
            }
            // Split the same physical driver into two tags, without changing any geometry.
            triangle.faceId = (horn ? center.x : center.z) >= 0 ? 50U : 51U;
        }
    QJsonArray vertices, triangles, tags, first, second;
    for (const auto p : mesh.vertices)
        vertices.append(point(p));
    for (const auto &triangle : mesh.triangles) {
        const auto v = triangle.vertices;
        triangles.append(QJsonArray{static_cast<double>(v[0]), static_cast<double>(v[1]),
                                    static_cast<double>(v[2])});
        tags.append(static_cast<double>(triangle.faceId));
        const bool driven = triangle.boundary == BoundaryCondition::driven;
        first.append(QJsonArray{driven ? 1.0 : 0.0, 0.0});
        second.append(QJsonArray{driven ? (triangle.faceId == 50 ? 1.0 : -1.0) : 0.0, 0.0});
    }
    QJsonArray faceVelocitiesUniform{QJsonObject{{"face_id", 50}, {"real", 1}, {"imag", 0}},
                                     QJsonObject{{"face_id", 51}, {"real", 1}, {"imag", 0}}};
    QJsonArray faceVelocitiesAntisymmetric{QJsonObject{{"face_id", 50}, {"real", 1}, {"imag", 0}},
                                           QJsonObject{{"face_id", 51}, {"real", -1}, {"imag", 0}}};
    return {
        {"schema", "hornsculpt-solver-input-1"},
        {"case", horn ? "horn-osse-8mm" : "sphere"},
        {"phasor", "exp(+iwt)"},
        {"amplitude", "RMS"},
        {"normal", "outward into fluid"},
        {"rho", 1.21},
        {"c", 343.0},
        {"vertices_m", vertices},
        {"triangles0", triangles},
        {"tags", tags},
        {"driven_tags", QJsonArray{50, 51}},
        {"driven_velocity", QJsonArray{1, 0}},
        {"normal_velocities_m_s", QJsonArray{first, second}},
        {"excitations",
         QJsonArray{QJsonObject{{"name", "uniform"}, {"face_velocities", faceVelocitiesUniform}},
                    QJsonObject{{"name", "antisymmetric"},
                                {"face_velocities", faceVelocitiesAntisymmetric}}}},
        {"frequencies_hz",
         horn ? QJsonArray{400, 1000, 2000, 4000} : QJsonArray{100, 500, 1000, 1715}},
        {"observations_m", observers(horn ? Vec3{0, 0, .18} : Vec3{})}};
}
} // namespace
int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addOption({"shape", "Geometrie: sphere oder horn.", "shape", "sphere"});
    parser.addOption({"refinement", "Anzahl Kugelverfeinerungen (0..5).", "level", "3"});
    parser.addOption({"output", "JSON-Ausgabedatei.", "path"});
    parser.process(app);
    try {
        if (!parser.isSet("output"))
            throw std::invalid_argument("--output is required.");
        const auto shape = parser.value("shape");
        if (shape != "sphere" && shape != "horn")
            throw std::invalid_argument("Unknown fixture shape.");
        bool parsed{};
        const unsigned level = parser.value("refinement").toUInt(&parsed);
        if (!parsed || level > 5)
            throw std::invalid_argument("Invalid sphere refinement.");
        const bool horn = shape == "horn";
        auto mesh = horn ? remeshIsotropic(DraftGeometry{}.previewMesh(), MeshingOptions{.008}).mesh
                         : testing::sphere(level);
        const auto fixture = exportFixture(std::move(mesh), horn);
        QSaveFile output(parser.value("output"));
        if (!output.open(QIODevice::WriteOnly))
            throw std::runtime_error(output.errorString().toStdString());
        const auto bytes = QJsonDocument(fixture).toJson(QJsonDocument::Compact);
        if (output.write(bytes) != bytes.size() || !output.commit())
            throw std::runtime_error("Fixture write failed.");
        std::fprintf(stdout, "Exported %lld triangles and two spatially distinct excitations.\n",
                     static_cast<long long>(fixture["triangles0"].toArray().size()));
        return 0;
    } catch (const std::exception &error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
