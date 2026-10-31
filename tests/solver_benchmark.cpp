#include "project_file.hpp"
#include "project_json.hpp"
#include "simulation_request.hpp"
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
#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <condition_variable>
#include <cuda/support_near_operator.hpp>
#include <cuda_runtime_api.h>
#include <hornsculpt/bem.hpp>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numbers>
#include <set>
#include <stdexcept>
#include <thread>

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

Fixture readProjectFixture(const QString &path) {
    QFile input(path);
    require(input.open(QIODevice::ReadOnly), "Project cannot be opened.");
    const auto hash = QCryptographicHash::hash(input.readAll(), QCryptographicHash::Sha256).toHex();
    const auto project = readProject(path);
    require(project.mesh != nullptr, "Project has no stored mesh.");
    Fixture fixture;
    fixture.mesh = project.mesh->mesh;
    fixture.inputHash = hash;
    std::vector<FaceVelocity> velocities;
    std::set<std::uint32_t> faces;
    for (const auto &triangle : fixture.mesh.triangles)
        if (triangle.boundary == BoundaryCondition::driven && faces.insert(triangle.faceId).second)
            velocities.push_back({triangle.faceId, {1, 0}});
    const auto origin = project.source == ProjectSource::freecad && project.cad
                            ? project.cad->acousticOrigin
                            : Vec3{0, 0, DraftGeometry(project.parameters).length()};
    const auto axis = project.source == ProjectSource::freecad && project.cad
                          ? project.cad->acousticAxis
                          : Vec3{0, 0, 1};
    fixture.request = makeSimulationRequest(project.simulation.settings, velocities, origin, axis);
    return fixture;
}

QJsonValue graphExecutionEvidence(const std::optional<BemFastMultipoleGraphEvidence> &value) {
    if (!value)
        return QJsonValue::Null;
    const auto &e = *value;
    return QJsonObject{
        {"capture_attempts", qint64(e.captureAttempts)},
        {"instantiations", qint64(e.instantiations)},
        {"launches", qint64(e.launches)},
        {"reuses", qint64(e.reuses)},
        {"executed_batches", qint64(e.executedBatches)},
        {"canonical_chunks", qint64(e.canonicalChunks)},
        {"canonical_batches", qint64(e.canonicalBatches)},
        {"invalidated_executables", qint64(e.invalidatedExecutables)},
        {"budget_bypasses", qint64(e.budgetBypasses)},
        {"admission_closures", qint64(e.admissionClosures)},
        {"maximum_cached_executables", qint64(e.maximumCachedExecutables)},
        {"maximum_cached_nodes", qint64(e.maximumCachedNodes)},
        {"peak_owned_host_metadata_bytes", qint64(e.peakOwnedHostMetadataBytes)},
        {"minimum_observed_free_device_bytes", qint64(e.minimumObservedFreeDeviceBytes)},
        {"total_device_bytes", qint64(e.totalDeviceBytes)},
        {"required_free_reserve_bytes", qint64(e.requiredFreeReserveBytes)},
        {"capture_milliseconds", e.captureMilliseconds},
        {"instantiation_milliseconds", e.instantiationMilliseconds},
        {"memory_observation_available", e.memoryObservationAvailable},
        {"scope", "frequency; all RHS/seeds/rules; owned metadata excludes opaque CUDA graph "
                  "driver storage"}};
}
QJsonObject frequencyData(const BemFrequencyResult &solved) {
    QJsonArray pressures, impedances, residuals, iterations, refinements, evaluations, lemSolutions,
        solvedVelocities, krylovResiduals, krylovTargets, discrepancies, evidence,
        quadratureEvidence;
    const auto levels = [](const std::vector<BemFastMultipoleLevelEvidence> &values) {
        QJsonArray result;
        for (const auto &v : values)
            result.append(QJsonObject{
                {"level", int(v.level)},
                {"representation", v.representation == BemFastMultipoleRepresentation::planeWave
                                       ? "planeWave"
                                       : "regularSpherical"},
                {"moment_bandlimit", int(v.momentBandlimit)},
                {"translation_order", int(v.translationOrder)},
                {"quadrature_order", int(v.quadratureOrder)}});
        return result;
    };
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
                            {"specific_impedance_pa_s_m",
                             face.volumeVelocity != std::complex<double>{} && face.area > 0
                                 ? QJsonValue(complexValue(face.volumeImpedance * face.area))
                                 : QJsonValue(QJsonValue::Null)},
                            {"volume_impedance_pa_s_m3", complexValue(face.volumeImpedance)}});
        pressures.append(observerPressure);
        impedances.append(faceImpedances);
        residuals.append(excitation.boundaryResidual);
        iterations.append(static_cast<int>(excitation.iterations));
        refinements.append(static_cast<int>(excitation.refinementSteps));
        evaluations.append(excitation.residualEvaluation == BemResidualEvaluation::directGalerkin
                               ? "directGalerkin"
                               : "higherOrderFastMultipole");
        krylovResiduals.append(excitation.fastMultipoleKrylovResidual);
        krylovTargets.append(excitation.fastMultipoleKrylovTarget);
        discrepancies.append(excitation.fastMultipoleOperatorDiscrepancy);
        evidence.append(QJsonObject{
            {"accuracy_seed", excitation.fastMultipoleAccuracySeed},
            {"comparison_accuracy_seed", excitation.fastMultipoleComparisonAccuracySeed},
            {"verification_accuracy_seed", excitation.fastMultipoleVerificationAccuracySeed},
            {"accuracy_refinements", int(excitation.fastMultipoleAccuracyRefinements)},
            {"comparison_levels", levels(excitation.fastMultipoleComparisonLevels)},
            {"verification_levels", levels(excitation.fastMultipoleVerificationLevels)}});
        quadratureEvidence.append(QJsonObject{
            {"far_comparison_degree", int(excitation.farQuadratureComparisonDegree)},
            {"far_verification_degree", int(excitation.farQuadratureVerificationDegree)},
            {"far_refinements", int(excitation.farQuadratureRefinements)},
            {"far_discrepancy", excitation.farQuadratureDiscrepancy},
            {"far_preparation_ms", excitation.farQuadraturePreparationMilliseconds},
            {"far_pressure_ms", excitation.farQuadraturePressureMilliseconds},
            {"far_flux_ms", excitation.farQuadratureFluxMilliseconds},
            {"field_comparison_degree", int(excitation.fieldQuadratureComparisonDegree)},
            {"field_verification_degree", int(excitation.fieldQuadratureVerificationDegree)},
            {"field_discrepancy", excitation.fieldQuadratureDiscrepancy}});
        lemSolutions.append(excitation.lemSolution ? QJsonValue(project_json::encodeLemSolution(
                                                         *excitation.lemSolution))
                                                   : QJsonValue(QJsonValue::Null));
        QJsonArray velocities;
        for (const auto &velocity : excitation.solvedVelocities)
            velocities.append(
                QJsonObject{{"face_id", static_cast<qint64>(velocity.faceId)},
                            {"normal_velocity_m_s", complexValue(velocity.normalVelocity)}});
        solvedVelocities.append(velocities);
    }
    return {
        {"frequency_hz", solved.frequency},
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
        {"residual_evaluation", evaluations},
        {"fast_multipole_krylov_residuals", krylovResiduals},
        {"inner_left_preconditioned_relative_residual_targets", krylovTargets},
        {"fast_multipole_operator_discrepancies", discrepancies},
        {"fast_multipole_evidence", evidence},
        {"quadrature_evidence", quadratureEvidence},
        {"lem_solutions", lemSolutions},
        {"solved_face_velocities_m_s", solvedVelocities},
        {"krylov_iterations", iterations},
        {"refinement_steps", refinements},
        {"fast_multipole_statistics",
         QJsonObject{
             {"graph_execution", graphExecutionEvidence(solved.fastMultipole.graphExecution)},
             {"node_count", static_cast<qint64>(solved.fastMultipole.nodeCount)},
             {"leaf_count", static_cast<qint64>(solved.fastMultipole.leafCount)},
             {"far_interaction_count",
              static_cast<qint64>(solved.fastMultipole.farInteractionCount)},
             {"pressure_pressure_near_pairs",
              qint64(solved.fastMultipole.pressurePressureNearPairs)},
             {"pressure_pressure_far_pairs", qint64(solved.fastMultipole.pressurePressureFarPairs)},
             {"pressure_flux_near_pairs", qint64(solved.fastMultipole.pressureFluxNearPairs)},
             {"pressure_flux_far_pairs", qint64(solved.fastMultipole.pressureFluxFarPairs)}}}};
}

QJsonObject accuracySettings(const BemNumerics &settings) {
    return {
        {"near_relative_tolerance", settings.nearRelativeTolerance},
        {"near_maximum_depth", static_cast<int>(settings.nearMaximumDepth)},
        {"singular_order", static_cast<int>(settings.singularOrder)},
        {"relative_residual_tolerance", settings.relativeResidualTolerance},
        {"maximum_refinement_steps", static_cast<int>(settings.maximumRefinementSteps)},
        {"fast_multipole",
         QJsonObject{
             {"order_policy", "automatic per level; basis seed=0.1*relative residual tolerance"},
             {"leaf_support_capacity",
              static_cast<int>(settings.fastMultipole.leafSupportCapacity)},
             {"separation_ratio", settings.fastMultipole.separationRatio},
             {"preconditioner_block_width", int(settings.fastMultipole.preconditionerBlockWidth)},
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
int measureNearOnly(const Fixture &fixture, const QString &inputPath, const QString &outputPath,
                    std::chrono::steady_clock::time_point start, unsigned timeoutSeconds,
                    bool profile) {
    using namespace wideband;
    require(fixture.request.frequencies.size() == 1, "Near-only requires one explicit frequency.");
    std::stop_source stop;
    std::atomic<bool> timedOut{};
    std::mutex mutex;
    std::condition_variable_any wake;
    std::jthread watchdog([&](std::stop_token watcherStop) {
        std::unique_lock lock(mutex);
        wake.wait_until(lock, watcherStop, start + std::chrono::seconds(timeoutSeconds),
                        [] { return false; });
        if (!watcherStop.stop_requested()) {
            timedOut = true;
            stop.request_stop();
        }
    });
    SupportNearConstructionStats stats;
    double spaceMilliseconds{}, treeMilliseconds{}, planMilliseconds{}, gpuMilliseconds{};
    std::size_t hostPlanBytes{}, hostPeakBytes{};
    QString error;
    std::unique_ptr<SupportNearOperator> op;
    try {
        const auto spaceStart = std::chrono::steady_clock::now();
        const auto space = galerkin::buildSpace(fixture.mesh);
        spaceMilliseconds = elapsed(spaceStart);
        SupportTreeSettings treeSettings;
        treeSettings.leafCapacity = fixture.request.numerics.fastMultipole.leafSupportCapacity;
        treeSettings.separationRatio = fixture.request.numerics.fastMultipole.separationRatio;
        treeSettings.nearDistanceRatio = fixture.request.numerics.nearDistanceRatio;
        const auto treeStart = std::chrono::steady_clock::now();
        const auto tree = makeSupportTree(space, treeSettings, stop.get_token());
        treeMilliseconds = elapsed(treeStart);
        std::set<std::uint32_t> active;
        for (const auto &panel : space.panels)
            if (panel.boundary == BoundaryCondition::driven)
                active.insert(panel.faceId);
        const std::vector<std::uint32_t> faces(active.begin(), active.end());
        const auto planStart = std::chrono::steady_clock::now();
        const auto plan =
            makeSupportNearPlan(space, tree, faces, 1024ull * 1024 * 1024, stop.get_token());
        planMilliseconds = elapsed(planStart);
        hostPlanBytes = plan.allocatedBytes();
        hostPeakBytes = plan.peakHostBytes;
        SupportNearOperatorSettings settings;
        settings.profileConstruction = profile;
        const auto &n = fixture.request.numerics;
        settings.quadrature = {n.quadratureOrder, n.singularOrder, n.nearDistanceRatio,
                               n.nearRelativeTolerance, n.nearMaximumDepth};
        const auto gpuStart = std::chrono::steady_clock::now();
        double lastPrint{};
        try {
            op = std::make_unique<SupportNearOperator>(
                space, tree, plan,
                2 * std::numbers::pi * fixture.request.frequencies[0] /
                    fixture.request.medium.soundSpeed,
                settings, stop.get_token(), nullptr,
                [&](const SupportNearConstructionStats &progress) {
                    stats = progress;
                    if (elapsed(start) - lastPrint > 5000) {
                        lastPrint = elapsed(start);
                        std::cout << "Near-only elapsed_ms=" << lastPrint
                                  << " visited_panel_jobs=" << stats.maskedPanelJobs
                                  << " batches=" << stats.batches
                                  << " max_kernel_ms=" << stats.maximumKernelMilliseconds
                                  << std::endl;
                    }
                });
            stats = op->constructionStats();
        } catch (...) {
            gpuMilliseconds = elapsed(gpuStart);
            throw;
        }
        gpuMilliseconds = elapsed(gpuStart);
    } catch (const std::exception &e) {
        error = QString::fromUtf8(e.what());
    }
    watchdog.request_stop();
    wake.notify_all();
    const bool complete = bool(op) && error.isEmpty() && !timedOut;
    writeReport(
        outputPath,
        {{"scope", "Actual wideband near-only preparation; no solve or acoustic result"},
         {"status", complete   ? "complete"
                    : timedOut ? "timeout"
                               : "failed"},
         {"error", error},
         {"input_path", inputPath},
         {"input_sha256", QString::fromLatin1(fixture.inputHash)},
         {"triangles", qint64(fixture.mesh.triangles.size())},
         {"frequency_hz", fixture.request.frequencies[0]},
         {"requested_near_distance_ratio", fixture.request.numerics.nearDistanceRatio},
         {"support_tree_near_distance_ratio", fixture.request.numerics.nearDistanceRatio},
         {"leaf_support_capacity", int(fixture.request.numerics.fastMultipole.leafSupportCapacity)},
         {"separation_ratio", fixture.request.numerics.fastMultipole.separationRatio},
         {"total_elapsed_ms_including_fixture_read", elapsed(start)},
         {"space_ms", spaceMilliseconds},
         {"tree_ms", treeMilliseconds},
         {"near_plan_ms", planMilliseconds},
         {"gpu_construction_ms", gpuMilliseconds},
         {"near_plan_host_bytes", qint64(hostPlanBytes)},
         {"near_plan_peak_host_bytes", qint64(hostPeakBytes)},
         {"allocated_stationary_device_payload_bytes", qint64(stats.liveDeviceBytes)},
         {"allocated_peak_device_payload_bytes", qint64(stats.peakDeviceBytes)},
         {"iterator_workspace_bytes", qint64(stats.iteratorWorkspaceBytes)},
         {"classification_workspace_bytes", qint64(stats.classificationWorkspaceBytes)},
         {"qualified_regular_degree4_jobs", qint64(stats.regularRuleJobs[0])},
         {"qualified_regular_degree6_jobs", qint64(stats.regularRuleJobs[1])},
         {"qualified_regular_degree8_jobs", qint64(stats.regularRuleJobs[2])},
         {"visited_masked_panel_jobs", qint64(stats.maskedPanelJobs)},
         {"completed_quadrature_jobs", qint64(stats.quadratureJobs)},
         {"completed_batches", qint64(stats.batches)},
         {"maximum_batch_jobs", qint64(stats.maximumBatchJobs)},
         {"maximum_batch_patches", qint64(stats.maximumBatchPatches)},
         {"maximum_batch_green_samples", qint64(stats.maximumBatchGreenSamples)},
         {"maximum_measured_kernel_ms", stats.maximumKernelMilliseconds},
         {"sum_cuda_event_kernel_ms", stats.kernelMilliseconds},
         {"batch_upload_kernel_and_sync_wall_ms", stats.batchUploadAndSyncWallMilliseconds},
         {"profile_cpu_timing_enabled", profile},
         {"iterator_cpu_ms", stats.iteratorMilliseconds},
         {"classification_cpu_ms", stats.classificationMilliseconds},
         {"adaptive_queue_decision_cpu_ms", stats.adaptiveDecisionMilliseconds},
         {"adaptive_estimator_gpu_kernel_ms", stats.estimatorKernelMilliseconds},
         {"adaptive_estimator_transfer_kernel_sync_wall_ms", stats.estimatorWallMilliseconds},
         {"adaptive_estimator_rounds", qint64(stats.adaptiveRounds)},
         {"adaptive_estimator_patch_evaluations", qint64(stats.adaptivePatchEvaluations)},
         {"maximum_estimator_batch_patches", qint64(stats.maximumEstimatorBatchPatches)},
         {"adaptive_quadrature_calls", qint64(stats.adaptiveQuadratureCalls)},
         {"adaptive_queue_and_transfer_host_peak_bytes", qint64(stats.adaptiveWorkspacePeakBytes)},
         {"maximum_measured_single_pair_workspace_bytes",
          qint64(stats.maximumSinglePairWorkspaceBytes)},
         {"iterator_batch_bytes", qint64(stats.iteratorBatchBytes)},
         {"adaptive_timing_scope",
          "CPU queue decisions and CUDA estimator kernel time are separate; estimator wall "
          "includes transfer, kernel and synchronization. Host peak includes shared states and "
          "reserved bounded transfer scratch."},
         {"memory_scope",
          "Owned successfully allocated payloads; excludes CUDA context, allocator and "
          "event internals. Stored peak remains valid after timeout frees buffers."},
         {"timeout_seconds", int(timeoutSeconds)}});
    return complete ? 0 : timedOut ? 124 : 1;
}
} // namespace

int main(int argc, char **argv) {
    QCoreApplication application(argc, argv);
    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addOption({"input", "Exact exported fixture JSON.", "path"});
    parser.addOption(
        {"project", "Read stored mesh and actual excitation/LEM settings from hsproj.", "path"});
    parser.addOption({"output", "New measurement JSON.", "path"});
    parser.addOption(
        {"validate-input-only", "Parse and summarize the fixture on CPU; no CUDA or solve."});
    parser.addOption({"cold-only", "One first complete frequency sweep with only the first RHS; no "
                                   "warmup or repeated solve."});
    parser.addOption({"all-rhs", "Retain every prescribed excitation in a cold reference sweep."});
    parser.addOption({"cold-timeout-seconds",
                      "Cooperative complete cold-run limit (1..600 seconds). Requires --cold-only.",
                      "seconds", "600"});
    parser.addOption({"near-only", "Measure only the actual wideband near construction; "
                                   "cooperative bounded limit, no solve."});
    parser.addOption(
        {"near-only-seconds", "Near-only cooperative limit (1..180 seconds).", "seconds", "90"});
    parser.addOption(
        {"profile-near", "Measure detailed near-construction CPU phases (adds timing overhead)."});
    parser.addOption(
        {"frequency", "Solve only this frequency in Hz (fixture bytes retain their hash).", "hz"});
    parser.addOption(
        {"near-maximum-depth", "Maximum adaptive near integration depth (0..10).", "depth"});
    parser.addOption({"singular-order", "Touching-panel integration order (2..20).", "order"});
    parser.addOption({"fast-multipole", "Use hierarchical Fast Multipole and restarted GMRES; "
                                        "without this option use dense CUDA LU."});
    parser.addOption(
        {"fmm-leaf-capacity",
         "Maximum complete supports per channel in a FMM leaf (4..256), requires --fast-multipole.",
         "count", QString::number(FastMultipoleSettings{}.leafSupportCapacity)});
    parser.addOption({"fmm-separation-ratio",
                      "FMM separation ratio (1.5..10), requires --fast-multipole.", "ratio",
                      QString::number(FastMultipoleSettings{}.separationRatio, 'g', 17)});
    parser.addOption({"fmm-block-width",
                      "Near-CSR block-Jacobi width in pressure DOFs (32..2048, powers of two).",
                      "count", QString::number(FastMultipoleSettings{}.preconditionerBlockWidth)});
    parser.addOption({"warmup", "Warmup calls per frequency and RHS scenario.", "count", "1"});
    parser.addOption(
        {"repetitions", "Measured calls per frequency and RHS scenario.", "count", "2"});
    parser.addOption({"near-distance-ratio",
                      "Near correction distance ratio; use 0 only for the explicit baseline "
                      "without near correction.",
                      "ratio", "1"});
    parser.process(application);
    try {
        require((parser.isSet("input") != parser.isSet("project")) && parser.isSet("output"),
                "Exactly one input fixture/project and an output are required.");
        require(!QFile::exists(parser.value("output")),
                "Refusing to replace previous measurements.");
        const bool coldOnly = parser.isSet("cold-only");
        require(!parser.isSet("all-rhs") || coldOnly, "--all-rhs requires --cold-only.");
        require(!coldOnly || (!parser.isSet("warmup") && !parser.isSet("repetitions")),
                "Cold-only cannot be combined with warmup or repetition options.");
        bool coldTimeoutOk = false;
        const auto coldTimeoutSeconds = parser.value("cold-timeout-seconds").toUInt(&coldTimeoutOk);
        require(coldTimeoutOk && coldTimeoutSeconds >= 1 && coldTimeoutSeconds <= 600,
                "Cold timeout must be 1..600 seconds.");
        require(!parser.isSet("cold-timeout-seconds") || coldOnly,
                "--cold-timeout-seconds requires --cold-only.");
        bool warmupOk = false, repetitionsOk = false, nearOk = false, leafCapacityOk = false,
             blockWidthOk = false, separationOk = false;
        const int warmup = parser.value("warmup").toInt(&warmupOk),
                  repetitions = parser.value("repetitions").toInt(&repetitionsOk);
        const double nearRatio = parser.value("near-distance-ratio").toDouble(&nearOk);
        const bool fastMultipole = parser.isSet("fast-multipole");
        const int leafCapacity = parser.value("fmm-leaf-capacity").toInt(&leafCapacityOk);
        const auto blockWidth = parser.value("fmm-block-width").toUInt(&blockWidthOk);
        require(!parser.isSet("fmm-block-width") || fastMultipole,
                "--fmm-block-width requires --fast-multipole.");
        require(blockWidthOk && (blockWidth == 32 || blockWidth == 64 || blockWidth == 128 ||
                                 blockWidth == 256 || blockWidth == 512 || blockWidth == 1024 ||
                                 blockWidth == 2048),
                "Fast Multipole block width must be a power of two from 32 to 2048.");
        const double separationRatio = parser.value("fmm-separation-ratio").toDouble(&separationOk);
        require(!parser.isSet("fmm-leaf-capacity") || fastMultipole,
                "--fmm-leaf-capacity requires --fast-multipole.");
        require(!parser.isSet("fmm-separation-ratio") || fastMultipole,
                "--fmm-separation-ratio requires --fast-multipole.");
        require(leafCapacityOk && leafCapacity >= 4 && leafCapacity <= 256,
                "Fast Multipole leaf capacity must be an integer from 4 to 256.");
        require(separationOk && std::isfinite(separationRatio) && separationRatio >= 1.5 &&
                    separationRatio <= 10,
                "Fast Multipole separation ratio must be finite and between 1.5 and 10.");
        require(warmupOk && repetitionsOk && warmup >= 1 && repetitions >= 2,
                "Invalid repetition counts.");
        require(nearOk && std::isfinite(nearRatio) && nearRatio >= 0,
                "Invalid near correction distance ratio.");
        const auto preparationStart = std::chrono::steady_clock::now();
        const auto inputPath = parser.value(parser.isSet("project") ? "project" : "input");
        auto fixture =
            parser.isSet("project") ? readProjectFixture(inputPath) : readFixture(inputPath);
        if (parser.isSet("validate-input-only")) {
            require(parser.isSet("input"), "Input validation requires an explicit fixture.");
            const auto finitePosition = [](const Vec3 &value) {
                return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
            };
            const bool finite =
                std::all_of(fixture.mesh.vertices.begin(), fixture.mesh.vertices.end(),
                            finitePosition) &&
                std::all_of(fixture.request.observers.begin(), fixture.request.observers.end(),
                            finitePosition) &&
                std::all_of(fixture.request.frequencies.begin(), fixture.request.frequencies.end(),
                            [](double value) { return std::isfinite(value); }) &&
                std::all_of(fixture.request.excitations.begin(), fixture.request.excitations.end(),
                            [](const BemExcitation &excitation) {
                                return std::all_of(
                                    excitation.faces.begin(), excitation.faces.end(),
                                    [](const auto &face) {
                                        return std::isfinite(face.normalVelocity.real()) &&
                                               std::isfinite(face.normalVelocity.imag());
                                    });
                            });
            require(finite, "Imported fixture contains non-finite values.");
            writeReport(parser.value("output"),
                        {{"schema", "solver-benchmark-input-validation-1"},
                         {"status", "fixture_parser_passed_not_solve"},
                         {"input_sha256", QString::fromLatin1(fixture.inputHash)},
                         {"vertices", static_cast<qint64>(fixture.mesh.vertices.size())},
                         {"triangles", static_cast<qint64>(fixture.mesh.triangles.size())},
                         {"rhs_count", static_cast<qint64>(fixture.request.excitations.size())},
                         {"observers", static_cast<qint64>(fixture.request.observers.size())},
                         {"frequencies", static_cast<qint64>(fixture.request.frequencies.size())},
                         {"parsed_values_finite", finite},
                         {"cuda_or_solver_called", false},
                         {"parser", "existing readFixture without numerical changes"}});
            return 0;
        }
        if (parser.isSet("frequency")) {
            bool ok{};
            const auto frequency = parser.value("frequency").toDouble(&ok);
            require(ok && std::isfinite(frequency) && frequency > 0,
                    "Frequency must be finite and positive.");
            fixture.request.frequencies = {frequency};
        }
        if (parser.isSet("near-maximum-depth")) {
            bool ok{};
            const auto depth = parser.value("near-maximum-depth").toUInt(&ok);
            require(ok && depth <= 10, "Near maximum depth must be from 0 to 10.");
            fixture.request.numerics.nearMaximumDepth = depth;
        }
        fixture.request.numerics.nearDistanceRatio = nearRatio;
        if (parser.isSet("singular-order")) {
            bool ok{};
            const auto order = parser.value("singular-order").toUInt(&ok);
            require(ok && order >= 2 && order <= 20, "Singular order must be from 2 to 20.");
            fixture.request.numerics.singularOrder = order;
        }
        fixture.request.numerics.acceleration =
            fastMultipole ? BemAcceleration::fastMultipole : BemAcceleration::dense;
        fixture.request.numerics.fastMultipole.leafSupportCapacity =
            static_cast<unsigned>(leafCapacity);
        fixture.request.numerics.fastMultipole.separationRatio = separationRatio;
        fixture.request.numerics.fastMultipole.preconditionerBlockWidth = blockWidth;
        if (parser.isSet("near-only")) {
            bool ok{};
            const auto seconds = parser.value("near-only-seconds").toUInt(&ok);
            require(ok && seconds >= 1 && seconds <= 180,
                    "Near-only duration must be 1..180 seconds.");
            return measureNearOnly(fixture, inputPath, parser.value("output"), preparationStart,
                                   seconds, parser.isSet("profile-near"));
        }
        const QString acceleration = fastMultipole ? "fastMultipole" : "dense";
        const QString precision = fastMultipole ? "Mixed ComplexF32/ComplexF64" : "ComplexF32";
        const auto constructorStart = std::chrono::steady_clock::now();
        CudaBemSolver solver(fixture.mesh);
        synchronize();
        const double constructorMilliseconds = elapsed(constructorStart);
        if (coldOnly) {
            auto request = fixture.request;
            if (!request.coupling && !parser.isSet("all-rhs"))
                request.excitations.resize(1);
            const auto solveStart = std::chrono::steady_clock::now();
            const auto phaseName = [](BemSolvePhase phase) -> const char * {
                switch (phase) {
                case BemSolvePhase::preparing:
                    return "preparing";
                case BemSolvePhase::nearConstruction:
                    return "near-construction";
                case BemSolvePhase::farPreparation:
                    return "far-preparation";
                case BemSolvePhase::krylovSolve:
                    return "krylov-solve";
                case BemSolvePhase::residualAudit:
                    return "residual-audit";
                case BemSolvePhase::fieldEvaluation:
                    return "field-evaluation";
                case BemSolvePhase::frequencyComplete:
                    return "frequency-complete";
                }
                throw std::logic_error("Unknown solver progress phase.");
            };
            QJsonArray phaseTimings;
            QJsonArray operatorEvents;
            std::optional<BemProgress> lastProgress;
            auto phaseStart = solveStart;
            unsigned lastPrintedIteration = 0;
            std::stop_source coldStop;
            std::atomic_bool coldTimedOut{};
            std::mutex coldMutex;
            std::condition_variable_any coldWake;
            std::jthread coldWatcher([&](std::stop_token watcherStop) {
                std::unique_lock lock(coldMutex);
                coldWake.wait_until(lock, watcherStop,
                                    preparationStart + std::chrono::seconds(coldTimeoutSeconds),
                                    [] { return false; });
                if (!watcherStop.stop_requested()) {
                    coldTimedOut = true;
                    coldStop.request_stop();
                }
            });
            BemResult result;
            try {
                result =
                    solver.solve(request, coldStop.get_token(), [&](const BemProgress &progress) {
                        const bool changed =
                            !lastProgress || lastProgress->phase != progress.phase ||
                            lastProgress->frequency != progress.frequency ||
                            lastProgress->rightHandSideIndex != progress.rightHandSideIndex;
                        if (changed && lastProgress) {
                            phaseTimings.append(
                                QJsonObject{{"frequency_hz", lastProgress->frequency},
                                            {"phase", phaseName(lastProgress->phase)},
                                            {"wall_ms", elapsed(phaseStart)},
                                            {"rhs_index", qint64(lastProgress->rightHandSideIndex)},
                                            {"rhs_count", qint64(lastProgress->rightHandSideCount)},
                                            {"completed_iterations", int(lastProgress->iterations)},
                                            {"pressure_operator_applications",
                                             qint64(lastProgress->pressureOperatorApplications)},
                                            {"face_flux_operator_applications",
                                             qint64(lastProgress->faceFluxOperatorApplications)}});
                            phaseStart = std::chrono::steady_clock::now();
                        }
                        const bool firstPressure =
                            progress.pressureOperatorApplications > 0 &&
                            (!lastProgress || lastProgress->pressureOperatorApplications == 0);
                        const bool firstFlux =
                            progress.faceFluxOperatorApplications > 0 &&
                            (!lastProgress || lastProgress->faceFluxOperatorApplications == 0);
                        const bool operatorChanged =
                            (lastProgress && (lastProgress->pressureOperatorApplications !=
                                                  progress.pressureOperatorApplications ||
                                              lastProgress->faceFluxOperatorApplications !=
                                                  progress.faceFluxOperatorApplications)) ||
                            firstPressure || firstFlux;
                        if (operatorChanged)
                            operatorEvents.append(
                                QJsonObject{{"elapsed_ms", elapsed(solveStart)},
                                            {"frequency_hz", progress.frequency},
                                            {"phase", phaseName(progress.phase)},
                                            {"rhs_index", qint64(progress.rightHandSideIndex)},
                                            {"completed_iterations", int(progress.iterations)},
                                            {"pressure_operator_applications",
                                             qint64(progress.pressureOperatorApplications)},
                                            {"face_flux_operator_applications",
                                             qint64(progress.faceFluxOperatorApplications)}});
                        if (changed || firstPressure || firstFlux || operatorChanged ||
                            progress.iterations >= lastPrintedIteration + 10) {
                            std::cout << "Own cold BM elapsed_ms=" << elapsed(solveStart)
                                      << " frequency_hz=" << progress.frequency
                                      << " rhs=" << progress.rightHandSideIndex << '/'
                                      << progress.rightHandSideCount
                                      << " phase=" << phaseName(progress.phase)
                                      << " iterations=" << progress.iterations
                                      << " Ap=" << progress.pressureOperatorApplications
                                      << " Bq=" << progress.faceFluxOperatorApplications
                                      << " completed=" << progress.completedFrequencies << '/'
                                      << progress.totalFrequencies << std::endl;
                            lastPrintedIteration = progress.iterations;
                        }
                        lastProgress = progress;
                    });
            } catch (const std::exception &failure) {
                coldWatcher.request_stop();
                if (lastProgress)
                    phaseTimings.append(
                        QJsonObject{{"frequency_hz", lastProgress->frequency},
                                    {"phase", phaseName(lastProgress->phase)},
                                    {"wall_ms", elapsed(phaseStart)},
                                    {"rhs_index", qint64(lastProgress->rightHandSideIndex)},
                                    {"rhs_count", qint64(lastProgress->rightHandSideCount)}});
                writeReport(
                    parser.value("output"),
                    QJsonObject{
                        {"schema", "hornsculpt-solver-cold-benchmark-1"},
                        {"role", "own-cold-solver"},
                        {"status", coldTimedOut ? "timeout" : "failed"},
                        {"error", failure.what()},
                        {"completed_result", false},
                        {"input_sha256", QString::fromLatin1(fixture.inputHash)},
                        {"input_path", inputPath},
                        {"triangles", qint64(fixture.mesh.triangles.size())},
                        {"frequencies_hz", QJsonArray::fromVariantList([&] {
                             QVariantList v;
                             for (auto f : request.frequencies)
                                 v.append(f);
                             return v;
                         }())},
                        {"rhs_count",
                         qint64(request.coupling
                                    ? lastProgress ? lastProgress->rightHandSideCount : 0
                                    : request.excitations.size())},
                        {"timeout_seconds", int(coldTimeoutSeconds)},
                        {"total_elapsed_ms_including_fixture_read", elapsed(preparationStart)},
                        {"geometry_constructor_ms", constructorMilliseconds},
                        {"accuracy_settings", accuracySettings(request.numerics)},
                        {"requested_near_distance_ratio", request.numerics.nearDistanceRatio},
                        {"support_tree_near_distance_ratio",
                         fastMultipole ? QJsonValue(request.numerics.nearDistanceRatio)
                                       : QJsonValue(QJsonValue::Null)},
                        {"leaf_support_capacity",
                         fastMultipole
                             ? QJsonValue(int(request.numerics.fastMultipole.leafSupportCapacity))
                             : QJsonValue(QJsonValue::Null)},
                        {"separation_ratio",
                         fastMultipole ? QJsonValue(request.numerics.fastMultipole.separationRatio)
                                       : QJsonValue(QJsonValue::Null)},
                        {"phase_wall_times", phaseTimings},
                        {"operator_application_events", operatorEvents},
                        {"last_phase", lastProgress ? phaseName(lastProgress->phase) : "preparing"},
                        {"last_rhs_index",
                         lastProgress ? qint64(lastProgress->rightHandSideIndex) : 0},
                        {"completed_iterations", lastProgress ? int(lastProgress->iterations) : 0},
                        {"pressure_operator_applications",
                         lastProgress ? qint64(lastProgress->pressureOperatorApplications) : 0},
                        {"face_flux_operator_applications",
                         lastProgress ? qint64(lastProgress->faceFluxOperatorApplications) : 0}});
                std::cerr << "Own cold BM " << (coldTimedOut ? "timeout" : "failure") << ": "
                          << failure.what() << std::endl;
                return coldTimedOut ? 124 : 1;
            }
            coldWatcher.request_stop();
            require(result.acceleration == request.numerics.acceleration,
                    "Solver returned a different acceleration mode than requested.");
            synchronize();
            const double solveMilliseconds = elapsed(solveStart);
            const double fullSweepMilliseconds = elapsed(constructorStart);
            require(result.frequencies.size() == request.frequencies.size(),
                    "Cold sweep frequency count mismatch.");
            QJsonArray frequencies;
            for (const auto &frequency : result.frequencies) {
                require(frequency.excitations.size() ==
                            (request.coupling ? 1 : request.excitations.size()),
                        "Cold sweep returned an unexpected RHS count.");
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
                    {"input_path", inputPath},
                    {"coupled_network", request.coupling.has_value()},
                    {"triangles", static_cast<qint64>(fixture.mesh.triangles.size())},
                    {"vertices", static_cast<qint64>(fixture.mesh.vertices.size())},
                    {"p1_dofs", static_cast<qint64>(result.pressureDofCount)},
                    {"dp0_dofs", static_cast<qint64>(fixture.mesh.triangles.size())},
                    {"quadrature_order", static_cast<int>(request.numerics.quadratureOrder)},
                    {"singular_order", static_cast<int>(request.numerics.singularOrder)},
                    {"accuracy_settings", accuracySettings(request.numerics)},
                    {"requested_near_distance_ratio", request.numerics.nearDistanceRatio},
                    {"support_tree_near_distance_ratio",
                     fastMultipole ? QJsonValue(request.numerics.nearDistanceRatio)
                                   : QJsonValue(QJsonValue::Null)},
                    {"leaf_support_capacity",
                     fastMultipole
                         ? QJsonValue(int(request.numerics.fastMultipole.leafSupportCapacity))
                         : QJsonValue(QJsonValue::Null)},
                    {"separation_ratio",
                     fastMultipole ? QJsonValue(request.numerics.fastMultipole.separationRatio)
                                   : QJsonValue(QJsonValue::Null)},
                    {"near_distance_ratio", nearRatio},
                    {"near_correction_enabled", nearRatio > 0},
                    {"symmetry", "off"},
                    {"rhs_count", qint64(request.coupling ? 1 : request.excitations.size())},
                    {"warmup_sweeps", 0},
                    {"solve_calls", 1},
                    {"geometry_constructor_ms", constructorMilliseconds},
                    {"first_complete_solve_ms", solveMilliseconds},
                    {"phase_wall_times", phaseTimings},
                    {"operator_application_events", operatorEvents},
                    {"progress_timing_scope",
                     "Coarse solver callbacks; wall intervals may contain nested preparation, "
                     "callback and synchronization costs. No per-panel clocks."},
                    {"operator_counter_scope", "Cumulative within the final completed frequency; "
                                               "phase snapshots identify frequency and RHS."},
                    {"pressure_operator_applications",
                     lastProgress ? qint64(lastProgress->pressureOperatorApplications) : 0},
                    {"face_flux_operator_applications",
                     lastProgress ? qint64(lastProgress->faceFluxOperatorApplications) : 0},
                    {"constructor_and_first_sweep_ms", fullSweepMilliseconds},
                    {"total_elapsed_ms_including_fixture_read", elapsed(preparationStart)},
                    {"timeout_seconds", int(coldTimeoutSeconds)},
                    {"status", "complete"},
                    {"completed_result", true},
                    {"cuda_device_name", QString::fromStdString(result.device)},
                    {"reported_peak_device_bytes", static_cast<qint64>(result.peakDeviceBytes)},
                    {"timing_note",
                     "One solver construction and one complete solve of the requested "
                     "frequencies and excitations. External process wall timing "
                     "includes CLI startup, fixture import and JSON output."},
                    {"utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODateWithMs)},
                    {"frequencies", frequencies}});
            std::cout << "Own cold BM: one complete requested sweep, constructor+solve-ms="
                      << fullSweepMilliseconds << std::endl;
            return 0;
        }
        QJsonArray samples;
        std::size_t pressureDofs{}, reportedPeakBytes{};
        QString device;
        std::set<std::size_t> rhsCounts{1};
        if (!fixture.request.coupling)
            rhsCounts.insert(fixture.request.excitations.size());
        for (const auto rhsCount : rhsCounts)
            for (const auto frequency : fixture.request.frequencies)
                for (int repetition = 0; repetition < warmup + repetitions; ++repetition) {
                    auto request = fixture.request;
                    request.frequencies = {frequency};
                    if (!request.coupling)
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
            {"input_path", inputPath},
            {"coupled_network", fixture.request.coupling.has_value()},
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
