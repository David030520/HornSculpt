#include "scaled_fixture_geometry.hpp"
#include <QCommandLineParser>
#include <QCoreApplication>
#include <QCryptographicHash>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <cmath>
#include <cstdio>
#include <hornsculpt/draft_geometry.hpp>
#include <hornsculpt/meshing.hpp>
#include <limits>
#include <map>
#include <set>
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

QJsonObject refineExistingFixture(const QString &path, std::size_t target) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly))
        throw std::runtime_error("Cannot read the source fixture.");
    const auto bytes = file.readAll();
    const auto document = QJsonDocument::fromJson(bytes);
    if (!document.isObject())
        throw std::invalid_argument("Source fixture must be a JSON object.");
    auto output = document.object();
    auto index = [](const QJsonValue &value) {
        const auto number = value.toDouble(-1);
        if (!value.isDouble() || !std::isfinite(number) || number < 0 ||
            number > std::numeric_limits<std::uint32_t>::max() || std::floor(number) != number)
            throw std::invalid_argument("Invalid fixture index.");
        return static_cast<std::uint32_t>(number);
    };
    SurfaceMesh source;
    for (const auto &value : output["vertices_m"].toArray()) {
        const auto p = value.toArray();
        if (p.size() != 3)
            throw std::invalid_argument("Invalid fixture vertex.");
        for (const auto coordinate : p)
            if (!coordinate.isDouble() || !std::isfinite(coordinate.toDouble()))
                throw std::invalid_argument("Invalid fixture coordinate.");
        source.vertices.push_back({p[0].toDouble(), p[1].toDouble(), p[2].toDouble()});
    }
    std::set<std::uint32_t> driven;
    for (const auto &tag : output["driven_tags"].toArray())
        driven.insert(index(tag));
    const auto panels = output["triangles0"].toArray(), tags = output["tags"].toArray();
    if (panels.empty() || panels.size() != tags.size() || driven.empty())
        throw std::invalid_argument("Source fixture needs panels, tags and explicit driven tags.");
    for (qsizetype i = 0; i < panels.size(); ++i) {
        const auto vertices = panels[i].toArray();
        if (vertices.size() != 3)
            throw std::invalid_argument("Invalid fixture triangle.");
        const auto tag = index(tags[i]);
        source.triangles.push_back(
            {{index(vertices[0]), index(vertices[1]), index(vertices[2])},
             tag,
             driven.contains(tag) ? BoundaryCondition::driven : BoundaryCondition::wall});
    }
    const auto original = inspectMesh(source);
    if (!original.watertight() || original.selfIntersectionPairs != 0 || original.signedVolume <= 0)
        throw std::invalid_argument("Source fixture is not an outward closed boundary.");
    auto areas = [](const SurfaceMesh &mesh) {
        std::map<std::uint32_t, double> result;
        for (const auto &triangle : mesh.triangles) {
            const auto a = mesh.vertices[triangle.vertices[0]],
                       b = mesh.vertices[triangle.vertices[1]],
                       c = mesh.vertices[triangle.vertices[2]];
            const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
            result[triangle.faceId] += .5 * std::hypot(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z,
                                                       u.x * v.y - u.y * v.x);
        }
        return result;
    };
    const auto originalAreas = areas(source);
    // Identity projection keeps every new edge midpoint on its original
    // polyhedral facets; it does not invent a more accurate curved surface.
    const auto refined = testing::refinePolyhedralFixture(std::move(source), target);
    const auto quality = inspectMesh(refined);
    const auto refinedAreas = areas(refined);
    constexpr double relativeIntegralPreservationTolerance = 1e-10;
    if (!quality.watertight() || quality.selfIntersectionPairs != 0 ||
        quality.minimumAngleDegrees < 20 || quality.maximumAspectRatio > 4 ||
        quality.connectedComponents != original.connectedComponents ||
        std::abs(quality.signedVolume / original.signedVolume - 1) >
            relativeIntegralPreservationTolerance ||
        originalAreas.size() != refinedAreas.size())
        throw std::runtime_error(
            "Refined fixture violates geometry, topology or quality preservation.");
    for (const auto &[face, area] : originalAreas)
        if (std::abs(refinedAreas.at(face) / area - 1) > relativeIntegralPreservationTolerance)
            throw std::runtime_error("Refinement changed a physical face or driven area.");
    QJsonArray vertices, triangles, faceTags, columns;
    for (const auto p : refined.vertices)
        vertices.append(point(p));
    for (const auto &triangle : refined.triangles) {
        triangles.append(QJsonArray{int(triangle.vertices[0]), int(triangle.vertices[1]),
                                    int(triangle.vertices[2])});
        faceTags.append(int(triangle.faceId));
    }
    const auto excitations = output["excitations"].toArray();
    const auto originalColumns = output["normal_velocities_m_s"].toArray();
    if (excitations.empty() || originalColumns.size() != excitations.size())
        throw std::invalid_argument("Refinement needs explicit face excitations and DP0 columns.");
    for (qsizetype e = 0; e < excitations.size(); ++e) {
        std::map<std::uint32_t, QJsonArray> velocity;
        for (const auto &value : excitations[e].toObject()["face_velocities"].toArray()) {
            const auto face = value.toObject();
            velocity.emplace(index(face["face_id"]), QJsonArray{face["real"], face["imag"]});
        }
        auto value = [&](std::uint32_t tag) {
            const auto found = velocity.find(tag);
            return found == velocity.end() ? QJsonArray{0, 0} : found->second;
        };
        const auto old = originalColumns[e].toArray();
        if (old.size() != panels.size())
            throw std::invalid_argument("Invalid source DP0 column length.");
        for (qsizetype p = 0; p < tags.size(); ++p)
            if (old[p].toArray() != value(index(tags[p])))
                throw std::invalid_argument(
                    "Source DP0 data do not match preserved face velocities.");
        QJsonArray column;
        for (const auto &triangle : refined.triangles)
            column.append(value(triangle.faceId));
        columns.append(column);
    }
    output.insert("vertices_m", vertices);
    output.insert("triangles0", triangles);
    output.insert("tags", faceTags);
    output.insert("normal_velocities_m_s", columns);
    output.insert("case", output["case"].toString() + "-polyhedral-refined");
    output.insert(
        "refinement",
        QJsonObject{{"source_sha256",
                     QString::fromLatin1(
                         QCryptographicHash::hash(bytes, QCryptographicHash::Sha256).toHex())},
                    {"source_triangles", panels.size()},
                    {"target_triangles", static_cast<qint64>(target)},
                    {"method", "Uniform 1-to-4 subdivision and conforming longest-edge bisection "
                               "on the unchanged source polyhedron"},
                    {"minimum_angle_degrees", quality.minimumAngleDegrees},
                    {"maximum_aspect_ratio", quality.maximumAspectRatio},
                    {"minimum_edge_m", quality.minimumEdgeLength},
                    {"median_edge_m", quality.medianEdgeLength},
                    {"maximum_edge_m", quality.maximumEdgeLength},
                    {"signed_volume_m3", quality.signedVolume}});
    return output;
}
} // namespace
int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addOption({"shape", "Geometrie: sphere oder horn.", "shape", "sphere"});
    parser.addOption({"input", "Bestehendes Polyedernetz geometrieerhaltend verfeinern.", "path"});
    parser.addOption({"refinement", "Anzahl Kugelverfeinerungen (0..5).", "level", "3"});
    parser.addOption({"triangles", "Exakte gerade Dreieckszahl für Kugel oder --input (mindestens 8).", "count"});
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
        std::size_t triangles{};
        if (parser.isSet("triangles")) {
            triangles = parser.value("triangles").toUInt(&parsed);
            if (!parsed || (!parser.isSet("input") && horn) || triangles < 8 ||
                triangles > 2000000 || triangles % 2 != 0)
                throw std::invalid_argument("Invalid exact fixture triangle count.");
        }
        if (parser.isSet("input") &&
            (!triangles || parser.value("input") == parser.value("output")))
            throw std::invalid_argument(
                "Source refinement needs --triangles and a distinct output.");
        QJsonObject fixture;
        if (parser.isSet("input"))
            fixture = refineExistingFixture(parser.value("input"), triangles);
        else {
            auto mesh =
                horn ? remeshIsotropic(DraftGeometry{}.previewMesh(), MeshingOptions{.008}).mesh
                : triangles ? testing::scaledSphere(triangles)
                            : testing::sphere(level);
            fixture = exportFixture(std::move(mesh), horn);
        }
        QSaveFile output(parser.value("output"));
        if (!output.open(QIODevice::WriteOnly))
            throw std::runtime_error(output.errorString().toStdString());
        const auto bytes = QJsonDocument(fixture).toJson(QJsonDocument::Compact);
        if (output.write(bytes) != bytes.size() || !output.commit())
            throw std::runtime_error("Fixture write failed.");
        std::fprintf(stdout, "Exported %lld triangles and %lld excitations.\n",
                     static_cast<long long>(fixture["triangles0"].toArray().size()),
                     static_cast<long long>(fixture["excitations"].toArray().size()));
        return 0;
    } catch (const std::exception &error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
