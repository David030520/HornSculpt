#include <QCommandLineParser>
#include <QCoreApplication>
#include <QCryptographicHash>
#include <QDateTime>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cuda_runtime_api.h>
#include <hornsculpt/bem.hpp>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <set>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::invalid_argument(message);
}
double number(const QJsonObject &object, const QString &name) {
    const auto value = object.value(name);
    require(value.isDouble() && std::isfinite(value.toDouble()),
            "Missing or invalid scalar fixture value.");
    return value.toDouble();
}
QJsonArray array(const QJsonObject &object, const QString &name) {
    require(object.value(name).isArray(), "Missing fixture array.");
    return object.value(name).toArray();
}
Vec3 position(const QJsonValue &value) {
    const auto coordinates = value.toArray();
    require(coordinates.size() == 3, "Position must contain three coordinates.");
    for (const auto &coordinate : coordinates)
        require(coordinate.isDouble() && std::isfinite(coordinate.toDouble()),
                "Invalid position coordinate.");
    return {coordinates[0].toDouble(), coordinates[1].toDouble(), coordinates[2].toDouble()};
}
std::uint32_t index(const QJsonValue &value) {
    const double scalar = value.toDouble(-1);
    require(value.isDouble() && scalar >= 0 &&
                scalar <= std::numeric_limits<std::uint32_t>::max() && std::floor(scalar) == scalar,
            "Invalid unsigned fixture index.");
    return static_cast<std::uint32_t>(scalar);
}
void synchronize() {
    const auto status = cudaDeviceSynchronize();
    if (status != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(status));
}
template <class Start> double elapsed(Start start) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
        .count();
}
QJsonArray complexValue(std::complex<double> value) {
    return {value.real(), value.imag()};
}

struct Fixture {
    SurfaceMesh mesh;
    BemRequest request;
    QByteArray inputHash;
};
Fixture readFixture(const QString &path) {
    QFile file(path);
    require(file.open(QIODevice::ReadOnly), "Fixture cannot be opened.");
    const auto bytes = file.readAll();
    QJsonParseError parseError;
    const auto document = QJsonDocument::fromJson(bytes, &parseError);
    require(parseError.error == QJsonParseError::NoError && document.isObject(),
            "Invalid fixture JSON.");
    const auto root = document.object();
    require(root.value("phasor").toString() == "exp(+iwt)",
            "Expected positive-time phasor fixture.");
    Fixture fixture;
    fixture.inputHash = QCryptographicHash::hash(bytes, QCryptographicHash::Sha256).toHex();
    fixture.request.medium = {number(root, "rho"), number(root, "c")};
    for (const auto &vertex : array(root, "vertices_m"))
        fixture.mesh.vertices.push_back(position(vertex));
    std::set<std::uint32_t> driven;
    if (root.contains("driven_tags"))
        for (const auto &tag : array(root, "driven_tags"))
            driven.insert(index(tag));
    const auto tags = array(root, "tags"), triangles = array(root, "triangles0");
    require(tags.size() == triangles.size(), "Triangle and tag counts differ.");
    if (driven.empty())
        for (const auto &tag : tags)
            driven.insert(index(tag));
    for (qsizetype i = 0; i < triangles.size(); ++i) {
        const auto vertices = triangles[i].toArray();
        require(vertices.size() == 3, "Triangle must contain three indices.");
        const std::array indices{index(vertices[0]), index(vertices[1]), index(vertices[2])};
        for (const auto vertex : indices)
            require(vertex < fixture.mesh.vertices.size(), "Triangle index is out of range.");
        const auto tag = index(tags[i]);
        fixture.mesh.triangles.push_back(
            {indices, tag,
             driven.contains(tag) ? BoundaryCondition::driven : BoundaryCondition::wall});
    }
    for (const auto &point : array(root, "observations_m"))
        fixture.request.observers.push_back(position(point));
    for (const auto &frequency : array(root, "frequencies_hz")) {
        require(frequency.isDouble() && frequency.toDouble() > 0, "Frequency must be positive.");
        fixture.request.frequencies.push_back(frequency.toDouble());
    }
    if (root.contains("excitations")) {
        for (const auto &value : array(root, "excitations")) {
            const auto definition = value.toObject();
            BemExcitation excitation;
            excitation.name = definition.value("name").toString().toStdString();
            for (const auto &faceValue : array(definition, "face_velocities")) {
                const auto face = faceValue.toObject();
                excitation.faces.push_back(
                    {index(face.value("face_id")), {number(face, "real"), number(face, "imag")}});
            }
            fixture.request.excitations.push_back(std::move(excitation));
        }
    } else {
        const auto velocity = array(root, "driven_velocity");
        require(velocity.size() == 2, "Invalid legacy complex velocity.");
        const std::complex<double> unit{velocity[0].toDouble(), velocity[1].toDouble()};
        for (const auto factor : {std::complex<double>{1, 0}, std::polar(.6, .71)}) {
            BemExcitation excitation;
            excitation.name =
                fixture.request.excitations.empty() ? "uniform" : "proportional-phase";
            for (const auto tag : driven)
                excitation.faces.push_back({tag, unit * factor});
            fixture.request.excitations.push_back(std::move(excitation));
        }
    }
    require(!fixture.request.excitations.empty(), "Fixture has no excitation.");
    if (root.contains("normal_velocities_m_s")) {
        const auto columns = array(root, "normal_velocities_m_s");
        require(static_cast<std::size_t>(columns.size()) == fixture.request.excitations.size(),
                "DP0 RHS count mismatch.");
        for (qsizetype column = 0; column < columns.size(); ++column) {
            const auto values = columns[column].toArray();
            require(values.size() == triangles.size(), "DP0 column length mismatch.");
            std::map<std::uint32_t, std::complex<double>> byFace;
            for (const auto &face :
                 fixture.request.excitations[static_cast<std::size_t>(column)].faces)
                byFace.emplace(face.faceId, face.normalVelocity);
            for (qsizetype row = 0; row < values.size(); ++row) {
                const auto complex = values[row].toArray();
                require(complex.size() == 2, "Invalid DP0 velocity.");
                const std::complex<double> expected{complex[0].toDouble(), complex[1].toDouble()};
                const auto &triangle = fixture.mesh.triangles[static_cast<std::size_t>(row)];
                const auto found = byFace.find(triangle.faceId);
                const auto actual = found == byFace.end() ? std::complex<double>{} : found->second;
                require(std::abs(expected - actual) < 1e-12,
                        "Face excitation does not reproduce fixture DP0 data.");
                require(triangle.boundary == BoundaryCondition::driven || std::abs(expected) == 0,
                        "Nonzero Wall velocity.");
            }
        }
    }
    return fixture;
}

QJsonObject frequencyData(const BemFrequencyResult &solved) {
    QJsonArray pressures, impedances, residuals, iterations, refinements;
    for (const auto &excitation : solved.excitations) {
        QJsonArray observerPressure, faceImpedances;
        for (const auto p : excitation.observerPressure)
            observerPressure.append(complexValue(p));
        for (const auto &face : excitation.impedances)
            faceImpedances.append(
                QJsonObject{{"face_id", static_cast<qint64>(face.faceId)},
                            {"area_m2", face.area},
                            {"mean_pressure_pa", complexValue(face.meanPressure)},
                            {"volume_velocity_m3_s", complexValue(face.volumeVelocity)},
                            {"volume_impedance_pa_s_m3", complexValue(face.volumeImpedance)}});
        pressures.append(observerPressure);
        impedances.append(faceImpedances);
        residuals.append(excitation.boundaryResidual);
        iterations.append(static_cast<int>(excitation.iterations));
        refinements.append(static_cast<int>(excitation.refinementSteps));
    }
    return {{"frequency_hz", solved.frequency},
            {"assembly_ms", solved.assemblyMilliseconds},
            {"factorization_ms", solved.factorizationMilliseconds},
            {"solve_ms", solved.solveMilliseconds},
            {"triangular_solve_ms", solved.triangularSolveMilliseconds},
            {"field_ms", solved.fieldMilliseconds},
            {"residual_check_ms", solved.residualCheckMilliseconds},
            {"download_ms", QJsonValue::Null},
            {"observer_pressure_pa", pressures},
            {"impedances", impedances},
            {"boundary_relative_l2_residuals", residuals},
            {"krylov_iterations", iterations},
            {"refinement_steps", refinements},
            {"fast_multipole_statistics",
             QJsonObject{{"node_count", static_cast<qint64>(solved.fastMultipole.nodeCount)},
                         {"leaf_count", static_cast<qint64>(solved.fastMultipole.leafCount)},
                         {"far_interaction_count",
                          static_cast<qint64>(solved.fastMultipole.farInteractionCount)},
                         {"direct_panel_pair_count",
                          static_cast<qint64>(solved.fastMultipole.directPanelPairCount)},
                         {"represented_far_panel_pairs",
                          static_cast<qint64>(solved.fastMultipole.representedFarPanelPairs)}}}};
}

QJsonObject accuracySettings(const BemNumerics &settings) {
    return {
        {"near_relative_tolerance", settings.nearRelativeTolerance},
        {"near_maximum_depth", static_cast<int>(settings.nearMaximumDepth)},
        {"relative_residual_tolerance", settings.relativeResidualTolerance},
        {"maximum_refinement_steps", static_cast<int>(settings.maximumRefinementSteps)},
        {"fast_multipole",
         QJsonObject{
             {"interpolation_order", static_cast<int>(settings.fastMultipole.interpolationOrder)},
             {"leaf_panel_capacity", static_cast<int>(settings.fastMultipole.leafPanelCapacity)},
             {"separation_ratio", settings.fastMultipole.separationRatio},
             {"maximum_iterations", static_cast<int>(settings.fastMultipole.maximumIterations)},
             {"restart_length", static_cast<int>(settings.fastMultipole.restartLength)}}}};
}

void writeReport(const QString &path, QJsonObject report) {
    QFile executable(QCoreApplication::applicationFilePath());
    require(executable.open(QIODevice::ReadOnly), "Cannot identify the measured executable.");
    QCryptographicHash executableHash(QCryptographicHash::Sha256);
    require(executableHash.addData(&executable), "Cannot hash the measured executable.");
    report.insert("executable_sha256", QString::fromLatin1(executableHash.result().toHex()));
    QSaveFile output(path);
    require(output.open(QIODevice::WriteOnly), "Benchmark output cannot be opened.");
    const auto bytes = QJsonDocument(report).toJson(QJsonDocument::Indented);
    require(output.write(bytes) == bytes.size() && output.commit(),
            "Benchmark output cannot be written.");
}
} // namespace

int main(int argc, char **argv) {
    QCoreApplication application(argc, argv);
    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addOption({"input", "Exact exported fixture JSON.", "path"});
    parser.addOption({"output", "New measurement JSON.", "path"});
    parser.addOption({"cold-only", "One first complete frequency sweep with only the first RHS; no "
                                   "warmup or repeated solve."});
    parser.addOption({"fast-multipole", "Use hierarchical Fast Multipole and restarted GMRES; "
                                        "without this option use dense CUDA LU."});
    parser.addOption({"fmm-order",
                      "Chebyshev interpolation order (3..8), requires --fast-multipole.", "order",
                      QString::number(FastMultipoleSettings{}.interpolationOrder)});
    parser.addOption({"fmm-leaf-capacity",
                      "Maximum panels per FMM leaf (4..256), requires --fast-multipole.", "count",
                      QString::number(FastMultipoleSettings{}.leafPanelCapacity)});
    parser.addOption({"fmm-separation-ratio",
                      "FMM separation ratio (1.5..10), requires --fast-multipole.", "ratio",
                      QString::number(FastMultipoleSettings{}.separationRatio, 'g', 17)});
    parser.addOption({"warmup", "Warmup calls per frequency and RHS scenario.", "count", "1"});
    parser.addOption(
        {"repetitions", "Measured calls per frequency and RHS scenario.", "count", "2"});
    parser.addOption({"near-distance-ratio",
                      "Near correction distance ratio; use 0 only for the explicit baseline "
                      "without near correction.",
                      "ratio", "1"});
    parser.process(application);
    try {
        require(parser.isSet("input") && parser.isSet("output"), "Input and output are required.");
        require(!QFile::exists(parser.value("output")),
                "Refusing to replace previous measurements.");
        const bool coldOnly = parser.isSet("cold-only");
        require(!coldOnly || (!parser.isSet("warmup") && !parser.isSet("repetitions")),
                "Cold-only cannot be combined with warmup or repetition options.");
        bool warmupOk = false, repetitionsOk = false, nearOk = false, orderOk = false,
             leafCapacityOk = false, separationOk = false;
        const int warmup = parser.value("warmup").toInt(&warmupOk),
                  repetitions = parser.value("repetitions").toInt(&repetitionsOk);
        const double nearRatio = parser.value("near-distance-ratio").toDouble(&nearOk);
        const bool fastMultipole = parser.isSet("fast-multipole");
        const int multipoleOrder = parser.value("fmm-order").toInt(&orderOk);
        const int leafCapacity = parser.value("fmm-leaf-capacity").toInt(&leafCapacityOk);
        const double separationRatio = parser.value("fmm-separation-ratio").toDouble(&separationOk);
        require(!parser.isSet("fmm-order") || fastMultipole,
                "--fmm-order requires --fast-multipole.");
        require(!parser.isSet("fmm-leaf-capacity") || fastMultipole,
                "--fmm-leaf-capacity requires --fast-multipole.");
        require(!parser.isSet("fmm-separation-ratio") || fastMultipole,
                "--fmm-separation-ratio requires --fast-multipole.");
        require(orderOk && multipoleOrder >= 3 && multipoleOrder <= 8,
                "Fast Multipole interpolation order must be an integer from 3 to 8.");
        require(leafCapacityOk && leafCapacity >= 4 && leafCapacity <= 256,
                "Fast Multipole leaf capacity must be an integer from 4 to 256.");
        require(separationOk && std::isfinite(separationRatio) && separationRatio >= 1.5 &&
                    separationRatio <= 10,
                "Fast Multipole separation ratio must be finite and between 1.5 and 10.");
        require(warmupOk && repetitionsOk && warmup >= 1 && repetitions >= 2,
                "Invalid repetition counts.");
        require(nearOk && std::isfinite(nearRatio) && nearRatio >= 0,
                "Invalid near correction distance ratio.");
        auto fixture = readFixture(parser.value("input"));
        fixture.request.numerics.nearDistanceRatio = nearRatio;
        fixture.request.numerics.acceleration =
            fastMultipole ? BemAcceleration::fastMultipole : BemAcceleration::dense;
        fixture.request.numerics.fastMultipole.interpolationOrder =
            static_cast<unsigned>(multipoleOrder);
        fixture.request.numerics.fastMultipole.leafPanelCapacity =
            static_cast<unsigned>(leafCapacity);
        fixture.request.numerics.fastMultipole.separationRatio = separationRatio;
        const QString acceleration = fastMultipole ? "fastMultipole" : "dense";
        const QString precision = fastMultipole ? "Mixed ComplexF32/ComplexF64" : "ComplexF32";
        const auto constructorStart = std::chrono::steady_clock::now();
        CudaBemSolver solver(fixture.mesh);
        synchronize();
        const double constructorMilliseconds = elapsed(constructorStart);
        if (coldOnly) {
            auto request = fixture.request;
            request.excitations.resize(1);
            const auto solveStart = std::chrono::steady_clock::now();
            const auto result = solver.solve(request, {}, [](const BemProgress &progress) {
                if (progress.completedFrequencies > 0)
                    std::cout << "Own cold BM completed-frequency-hz=" << progress.frequency
                              << " completed=" << progress.completedFrequencies << '/'
                              << progress.totalFrequencies << std::endl;
            });
            require(result.acceleration == request.numerics.acceleration,
                    "Solver returned a different acceleration mode than requested.");
            synchronize();
            const double solveMilliseconds = elapsed(solveStart);
            const double fullSweepMilliseconds = elapsed(constructorStart);
            require(result.frequencies.size() == request.frequencies.size(),
                    "Cold sweep frequency count mismatch.");
            QJsonArray frequencies;
            for (const auto &frequency : result.frequencies) {
                require(frequency.excitations.size() == 1,
                        "Cold sweep must contain exactly one RHS.");
                frequencies.append(frequencyData(frequency));
            }
            writeReport(
                parser.value("output"),
                QJsonObject{
                    {"schema", "hornsculpt-solver-cold-benchmark-1"},
                    {"role", "own-cold-solver"},
                    {"backend", "HornSculpt CUDA P1/DP0 Burton-Miller"},
                    {"precision", precision},
                    {"acceleration", acceleration},
                    {"accuracy_qualified", false},
                    {"performance_qualified", false},
                    {"input_sha256", QString::fromLatin1(fixture.inputHash)},
                    {"input_path", parser.value("input")},
                    {"triangles", static_cast<qint64>(fixture.mesh.triangles.size())},
                    {"vertices", static_cast<qint64>(fixture.mesh.vertices.size())},
                    {"p1_dofs", static_cast<qint64>(result.pressureDofCount)},
                    {"dp0_dofs", static_cast<qint64>(fixture.mesh.triangles.size())},
                    {"quadrature_order", static_cast<int>(request.numerics.quadratureOrder)},
                    {"singular_order", static_cast<int>(request.numerics.singularOrder)},
                    {"accuracy_settings", accuracySettings(request.numerics)},
                    {"near_distance_ratio", nearRatio},
                    {"near_correction_enabled", nearRatio > 0},
                    {"symmetry", "off"},
                    {"rhs_count", 1},
                    {"warmup_sweeps", 0},
                    {"solve_calls", 1},
                    {"geometry_constructor_ms", constructorMilliseconds},
                    {"first_complete_solve_ms", solveMilliseconds},
                    {"constructor_and_first_sweep_ms", fullSweepMilliseconds},
                    {"cuda_device_name", QString::fromStdString(result.device)},
                    {"reported_peak_device_bytes", static_cast<qint64>(result.peakDeviceBytes)},
                    {"timing_note", "One solver construction and one complete solve of all fixture "
                                    "frequencies and first RHS. External process wall timing "
                                    "includes CLI startup, fixture import and JSON output."},
                    {"utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODateWithMs)},
                    {"frequencies", frequencies}});
            std::cout << "Own cold BM: one complete sweep, first RHS, constructor+solve-ms="
                      << fullSweepMilliseconds << std::endl;
            return 0;
        }
        QJsonArray samples;
        std::size_t pressureDofs{}, reportedPeakBytes{};
        QString device;
        std::set<std::size_t> rhsCounts{1, fixture.request.excitations.size()};
        for (const auto rhsCount : rhsCounts)
            for (const auto frequency : fixture.request.frequencies)
                for (int repetition = 0; repetition < warmup + repetitions; ++repetition) {
                    auto request = fixture.request;
                    request.frequencies = {frequency};
                    request.excitations.resize(rhsCount);
                    synchronize();
                    const auto start = std::chrono::steady_clock::now();
                    const auto result = solver.solve(request);
                    require(result.acceleration == request.numerics.acceleration,
                            "Solver returned a different acceleration mode than requested.");
                    synchronize();
                    const double totalMilliseconds = elapsed(start);
                    require(result.frequencies.size() == 1,
                            "Solver returned an unexpected frequency count.");
                    const auto &solved = result.frequencies.front();
                    pressureDofs = result.pressureDofCount;
                    reportedPeakBytes = std::max(reportedPeakBytes, result.peakDeviceBytes);
                    device = QString::fromStdString(result.device);
                    auto sample = frequencyData(solved);
                    sample.insert("rhs_count", static_cast<qint64>(rhsCount));
                    sample.insert("repetition", repetition - warmup + 1);
                    sample.insert("warmup", repetition < warmup);
                    sample.insert("warm_ready_total_ms", totalMilliseconds);
                    samples.append(sample);
                    std::cout << "Own BM rhs=" << rhsCount << " f=" << frequency
                              << " rep=" << repetition << " total-ms=" << totalMilliseconds
                              << std::endl;
                }
        const QJsonObject report{
            {"schema", "hornsculpt-solver-benchmark-1"},
            {"role", "own-solver"},
            {"backend", "HornSculpt CUDA P1/DP0 Burton-Miller"},
            {"precision", precision},
            {"acceleration", acceleration},
            {"accuracy_qualified", false},
            {"performance_qualified", false},
            {"cuda_device_name", device},
            {"input_sha256", QString::fromLatin1(fixture.inputHash)},
            {"input_path", parser.value("input")},
            {"triangles", static_cast<qint64>(fixture.mesh.triangles.size())},
            {"vertices", static_cast<qint64>(fixture.mesh.vertices.size())},
            {"p1_dofs", static_cast<qint64>(pressureDofs)},
            {"dp0_dofs", static_cast<qint64>(fixture.mesh.triangles.size())},
            {"quadrature_order", static_cast<int>(fixture.request.numerics.quadratureOrder)},
            {"singular_order", static_cast<int>(fixture.request.numerics.singularOrder)},
            {"accuracy_settings", accuracySettings(fixture.request.numerics)},
            {"symmetry", "off"},
            {"near_distance_ratio", nearRatio},
            {"geometry_constructor_ms", constructorMilliseconds},
            {"near_correction_enabled", nearRatio > 0},
            {"warmup_repetitions", warmup},
            {"measured_repetitions", repetitions},
            {"iteration_order", "RHS count, frequency, warmup/repetitions"},
            {"configuration_prepare_scope",
             fastMultipole
                 ? "Measured calls repeat the same frequency and RHS count. All operator "
                   "preparation performed by solve is included; no setup time is subtracted. "
                   "Raw warmup wall times are retained."
                 : "First warmup prepares the configuration; measured calls repeat the same "
                   "frequency "
                   "and RHS count without plan changes. Raw warmup wall times are retained."},
            {"reported_peak_device_bytes", static_cast<qint64>(reportedPeakBytes)},
            {"memory_note", "Solver-reported owned allocation accounting; not independently "
                            "profiled process peak."},
            {"timing_note", "Device-synchronized complete single-frequency solve; geometry "
                            "retained. Download is not separately exposed by the solver API."},
            {"utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODateWithMs)},
            {"samples", samples}};
        writeReport(parser.value("output"), report);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Solver benchmark: " << error.what() << '\n';
        return 1;
    }
}
