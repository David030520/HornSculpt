#include <QCommandLineParser>
#include <QCoreApplication>
#include <QCryptographicHash>
#include <QElapsedTimer>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <galerkin/interaction_plan.hpp>
#include <galerkin/space.hpp>
#include <hornsculpt/meshing.hpp>
#include <limits>
#include <numbers>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::invalid_argument(message);
}
double scalar(const QJsonValue &value) {
    require(value.isDouble() && std::isfinite(value.toDouble()), "Expected a finite scalar.");
    return value.toDouble();
}
std::uint32_t node(const QJsonValue &value) {
    const auto index = scalar(value);
    require(index >= 0 && index <= std::numeric_limits<std::uint32_t>::max() &&
                std::floor(index) == index,
            "Invalid triangle index.");
    return static_cast<std::uint32_t>(index);
}
SurfaceMesh readMesh(const QJsonObject &fixture) {
    SurfaceMesh mesh;
    require(fixture["vertices_m"].isArray() && fixture["triangles0"].isArray(),
            "Missing fixture mesh.");
    for (const auto value : fixture["vertices_m"].toArray()) {
        const auto p = value.toArray();
        require(p.size() == 3, "Expected three vertex coordinates.");
        mesh.vertices.push_back({scalar(p[0]), scalar(p[1]), scalar(p[2])});
    }
    for (const auto value : fixture["triangles0"].toArray()) {
        const auto t = value.toArray();
        require(t.size() == 3, "Expected three triangle nodes.");
        // Tags and velocities do not affect geometric pair membership; the input
        // hash still identifies their unchanged original fixture definitions.
        mesh.triangles.push_back({{node(t[0]), node(t[1]), node(t[2])}, 1});
    }
    require(inspectMesh(mesh).watertight(), "Expected a closed, non-intersecting fixture mesh.");
    return mesh;
}
} // namespace

int main(int argc, char **argv) {
    QCoreApplication application(argc, argv);
    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addOption({"input", "Unchanged exported fixture JSON.", "path"});
    parser.addOption({"output", "New directed near-pair plan JSON.", "path"});
    parser.process(application);
    try {
        require(parser.isSet("input") && parser.isSet("output"), "Input and output are required.");
        require(!QFile::exists(parser.value("output")),
                "Refusing to replace a previous near plan.");
        QFile input(parser.value("input"));
        require(input.open(QIODevice::ReadOnly), "Cannot open fixture.");
        const auto bytes = input.readAll();
        QJsonParseError error;
        const auto document = QJsonDocument::fromJson(bytes, &error);
        require(error.error == QJsonParseError::NoError && document.isObject(),
                "Invalid fixture JSON.");
        const auto fixture = document.object();
        require(fixture["schema"].toString() == "hornsculpt-solver-input-1",
                "Invalid fixture schema.");
        const double soundSpeed = scalar(fixture["c"]);
        require(soundSpeed > 0 && fixture["frequencies_hz"].isArray(),
                "Invalid frequency/medium data.");
        double maximumFrequency = 0;
        for (const auto value : fixture["frequencies_hz"].toArray()) {
            const auto frequency = scalar(value);
            require(frequency > 0, "Frequency must be positive.");
            maximumFrequency = std::max(maximumFrequency, frequency);
        }
        require(maximumFrequency > 0, "No frequencies.");
        const double maximumK = 2 * std::numbers::pi * maximumFrequency / soundSpeed;
        auto mesh = readMesh(fixture);
        // Match CudaBemSolver's coordinate convention before classifying pairs:
        // subtraction can change floating-point decisions at the strict near boundary.
        const auto origin = mesh.vertices.front();
        for (auto &vertex : mesh.vertices) {
            vertex.x -= origin.x;
            vertex.y -= origin.y;
            vertex.z -= origin.z;
        }
        const auto space = galerkin::buildSpace(mesh);
        galerkin::PairQuadratureSettings settings;
        settings.nearDistanceRatio = 1;
        settings.relativeTolerance = 1e-4;
        settings.maximumSubdivisionDepth = 6;
        QElapsedTimer elapsed;
        elapsed.start();
        const auto plan = galerkin::makeInteractionPlan(space, settings, maximumK);
        const auto planningMilliseconds = elapsed.nsecsElapsed() * 1e-6;
        std::vector<std::pair<std::uint32_t, std::uint32_t>> directedPairs;
        for (const auto &job : plan.near)
            directedPairs.emplace_back(job.test, job.trial);
        if (!plan.wholePanelDegrees.empty())
            for (std::uint32_t trial = 1; trial < space.panels.size(); ++trial)
                for (std::uint32_t test = 0; test < trial; ++test)
                    if (plan.wholePanelDegrees[galerkin::unorderedPanelPairIndex(test, trial)]) {
                        directedPairs.emplace_back(test, trial);
                        directedPairs.emplace_back(trial, test);
                    }
        std::sort(directedPairs.begin(), directedPairs.end());
        QJsonArray pairs;
        for (const auto &[test, trial] : directedPairs)
            pairs.append(QJsonArray{static_cast<qint64>(test), static_cast<qint64>(trial)});
        const QJsonObject report{
            {"schema", "hornsculpt-near-interaction-plan-1"},
            {"fixture_sha256",
             QString::fromLatin1(
                 QCryptographicHash::hash(bytes, QCryptographicHash::Sha256).toHex())},
            {"triangle_pairs0", pairs},
            {"coordinate_convention", "fixture_vertices_minus_first_vertex"},
            {"origin_m", QJsonArray{origin.x, origin.y, origin.z}},
            {"distance_ratio", settings.nearDistanceRatio},
            {"max_wavenumber_m_inv", maximumK},
            {"relative_tolerance", settings.relativeTolerance},
            {"maximum_subdivision_depth", static_cast<int>(settings.maximumSubdivisionDepth)},
            {"planning_ms", planningMilliseconds},
            {"touching_jobs", static_cast<qint64>(plan.touching.size())},
            {"directed_near_jobs", static_cast<qint64>(directedPairs.size())},
            {"shared_patches", static_cast<qint64>(plan.patches.size())},
            {"maximum_estimated_relative_error", plan.maximumEstimatedRelativeError}};
        QSaveFile output(parser.value("output"));
        require(output.open(QIODevice::WriteOnly), "Cannot open plan output.");
        const auto serialized = QJsonDocument(report).toJson(QJsonDocument::Compact);
        require(output.write(serialized) == serialized.size() && output.commit(),
                "Cannot write plan.");
        std::printf("Exported %zu directed near pairs, %zu shared patches in %.3f ms.\n",
                    directedPairs.size(), plan.patches.size(), planningMilliseconds);
        return 0;
    } catch (const std::exception &error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
