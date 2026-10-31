#include "fixture_geometry.hpp"
#include <galerkin/interaction_plan.hpp>
#include <galerkin/space.hpp>
#include <QCoreApplication>
#include <QCryptographicHash>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QProcess>
#include <QTemporaryDir>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
void checkExporter(const QString &exporter) {
    auto mesh = testing::sphere(4);
    QJsonArray vertices, triangles;
    for (const auto p : mesh.vertices) vertices.append(QJsonArray{p.x, p.y, p.z});
    for (const auto &t : mesh.triangles)
        triangles.append(QJsonArray{static_cast<qint64>(t.vertices[0]),
                                    static_cast<qint64>(t.vertices[1]),
                                    static_cast<qint64>(t.vertices[2])});
    const QJsonObject fixture{{"schema", "hornsculpt-solver-input-1"},
                              {"vertices_m", vertices}, {"triangles0", triangles},
                              {"c", 343.0}, {"frequencies_hz", QJsonArray{100, 1715}}};
    const auto bytes = QJsonDocument(fixture).toJson(QJsonDocument::Compact);
    QTemporaryDir dir;
    require(dir.isValid(), "Cannot create test directory");
    const auto input = dir.filePath("fixture.json"), output = dir.filePath("plan.json");
    QFile f(input);
    require(f.open(QIODevice::WriteOnly) && f.write(bytes) == bytes.size(), "Cannot write fixture");
    f.close();
    QProcess process;
    process.start(exporter, {"--input", input, "--output", output});
    require(process.waitForFinished(120000) && process.exitStatus() == QProcess::NormalExit &&
                process.exitCode() == 0, "Near exporter failed");
    QFile result(output);
    require(result.open(QIODevice::ReadOnly), "Cannot read exported plan");
    const auto plan = QJsonDocument::fromJson(result.readAll()).object();
    require(plan["fixture_sha256"].toString().toLatin1() ==
                QCryptographicHash::hash(bytes, QCryptographicHash::Sha256).toHex(),
            "World fixture hash changed");
    require(plan["coordinate_convention"].toString() == "fixture_vertices_minus_first_vertex",
            "Coordinate convention missing");
    const auto origin = mesh.vertices.front();
    const auto exportedOrigin = plan["origin_m"].toArray();
    require(exportedOrigin.size() == 3 && exportedOrigin[0].toDouble() == origin.x &&
                exportedOrigin[1].toDouble() == origin.y && exportedOrigin[2].toDouble() == origin.z,
            "World-space origin mismatch");
    const auto world = galerkin::buildSpace(mesh);
    for (auto &v : mesh.vertices) { v.x -= origin.x; v.y -= origin.y; v.z -= origin.z; }
    const auto centered = galerkin::buildSpace(mesh);
    galerkin::PairQuadratureSettings settings;
    settings.nearDistanceRatio = 1;
    settings.relativeTolerance = 1e-4;
    settings.maximumSubdivisionDepth = 6;
    const double k = 2 * std::numbers::pi * 1715.0 / 343.0;
    const auto expected = galerkin::makeInteractionPlan(centered, settings, k);
    const auto uncentered = galerkin::makeInteractionPlan(world, settings, k);
    const auto pairs = plan["triangle_pairs0"].toArray();
    require(static_cast<std::size_t>(pairs.size()) == expected.near.size(), "Pair count mismatch");
    bool different = uncentered.near.size() != expected.near.size();
    for (qsizetype i = 0; i < pairs.size(); ++i) {
        const auto pair = pairs[i].toArray();
        const auto &job = expected.near[static_cast<std::size_t>(i)];
        require(pair.size() == 2 && pair[0].toInteger() == job.test &&
                    pair[1].toInteger() == job.trial, "Exporter does not use centered plan");
        if (!different && static_cast<std::size_t>(i) < uncentered.near.size()) {
            const auto &old = uncentered.near[static_cast<std::size_t>(i)];
            different = old.test != job.test || old.trial != job.trial;
        }
    }
    require(different, "Sphere fixture does not exercise a translation-sensitive near boundary");
}
} // namespace
int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    try {
        require(argc == 2, "Exporter path required");
        checkExporter(QString::fromLocal8Bit(argv[1]));
        std::cout << "Near exporter uses the solver's centered coordinate convention.\n";
        return 0;
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
