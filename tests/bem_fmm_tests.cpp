#include "fixture_geometry.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <cuda_runtime_api.h>
#include <fstream>
#include <future>
#include <galerkin/space.hpp>
#include <hornsculpt/bem.hpp>
#include <iomanip>
#include <iostream>
#include <locale>
#include <numbers>
#include <span>
#include <sstream>
#include <stdexcept>
#include <stop_token>
#include <string>
#include <string_view>
#include <thread>

using namespace hornsculpt;

namespace {
using Complex = std::complex<double>;
constexpr Complex imaginary{0, 1};
constexpr Complex mixture{.3, .4};
// These compare the same discretization. Physical sphere gates additionally
// include the faceted geometry error and are deliberately separate.
constexpr double parityLimit = 1e-3;
constexpr double sphereErrorLimit = .04;

void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}

double relativeError(std::span<const Complex> actual, std::span<const Complex> expected) {
    require(actual.size() == expected.size() && !actual.empty(),
            "Comparison requires equal nonempty vectors.");
    double difference = 0, reference = 0;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        require(std::isfinite(actual[i].real()) && std::isfinite(actual[i].imag()),
                "FMM produced nonfinite pressure.");
        difference += std::norm(actual[i] - expected[i]);
        reference += std::norm(expected[i]);
    }
    require(reference > 0, "Comparison has zero reference norm.");
    return std::sqrt(difference / reference);
}

double compare(const BemResult &actual, const BemResult &expected, double limit) {
    require(actual.triangleCount == expected.triangleCount &&
                actual.pressureDofCount == expected.pressureDofCount &&
                actual.frequencies.size() == expected.frequencies.size(),
            "FMM changed the mesh, pressure space or frequency count.");
    double worst = 0;
    for (std::size_t f = 0; f < actual.frequencies.size(); ++f) {
        const auto &a = actual.frequencies[f], &b = expected.frequencies[f];
        require(a.frequency == b.frequency && a.excitations.size() == b.excitations.size(),
                "FMM changed frequency or excitation ordering.");
        for (std::size_t e = 0; e < a.excitations.size(); ++e) {
            const auto &x = a.excitations[e], &y = b.excitations[e];
            require(x.name == y.name, "FMM changed the excitation name.");
            const double surface = relativeError(x.boundary()->surfacePressure, y.boundary()->surfacePressure);
            const double field = relativeError(x.observerPressure, y.observerPressure);
            worst = std::max({worst, surface, field});
            require(x.impedances.size() == y.impedances.size(),
                    "FMM changed the face impedance count.");
            double pressureNorm = 0, pressureDifference = 0;
            for (std::size_t face = 0; face < x.impedances.size(); ++face) {
                const auto &xa = x.impedances[face], &ya = y.impedances[face];
                require(xa.faceId == ya.faceId && xa.area == ya.area &&
                            std::abs(xa.volumeVelocity - ya.volumeVelocity) <=
                                limit * std::max(1e-18, std::abs(ya.volumeVelocity)),
                        "FMM changed the driven boundary or its physical normalization.");
                pressureNorm += std::norm(ya.meanPressure) * ya.area;
                pressureDifference += std::norm(xa.meanPressure - ya.meanPressure) * ya.area;
                if (std::abs(xa.volumeVelocity) > 0)
                    require(std::abs(xa.volumeImpedance * xa.volumeVelocity - xa.meanPressure) <
                                1e-9 * std::max(1.0, std::abs(xa.meanPressure)),
                            "FMM impedance does not equal pressure divided by volume flow.");
            }
            const double impedance = std::sqrt(pressureDifference / pressureNorm);
            worst = std::max(worst, impedance);
            std::cout << a.frequency << " Hz, " << x.name << ": surface/field/port L2 " << surface
                      << '/' << field << '/' << impedance << ", residual " << x.boundaryResidual
                      << '\n';
        }
    }
    require(worst < limit, "Fast-multipole versus dense complex parity gate failed.");
    return worst;
}

BemRequest sphereRequest() {
    BemRequest value;
    value.frequencies = {100, 500, 1000, 1715};
    value.excitations = {{"uniform", {{1, 1.0}}}, {"complex phase", {{1, std::polar(.6, .71)}}}};
    value.observers = {{0, 0, .3}, {1, 0, 0}, {0, 0, 3}};
    return value;
}

BemRequest accelerated(BemRequest value, double residualTolerance = 0) {
    value.numerics.acceleration = BemAcceleration::fastMultipole;
    if (residualTolerance > 0)
        value.numerics.relativeResidualTolerance = residualTolerance;
    value.numerics.fastMultipole.leafSupportCapacity = 4;
    return value;
}

void checkResidual(const BemResult &result, const BemRequest &request) {
    for (const auto &frequency : result.frequencies)
        for (const auto &excitation : frequency.excitations) {
            require(std::isfinite(excitation.boundaryResidual) &&
                        excitation.boundaryResidual <= request.numerics.relativeResidualTolerance,
                    "FMM did not meet its requested final audit residual.");
            const bool higherOrder =
                excitation.residualEvaluation == BemResidualEvaluation::higherOrderFastMultipole;
            require(excitation.fastMultipoleAccuracySeed > 0 &&
                        !excitation.historicalChebyshevEvidence,
                    "New wideband solve lacks automatic seed evidence or retained historical "
                    "interpolation metadata.");
            require((higherOrder && excitation.fastMultipoleComparisonAccuracySeed > 0 &&
                     excitation.fastMultipoleVerificationAccuracySeed > 0 &&
                     excitation.fastMultipoleVerificationAccuracySeed <
                         excitation.fastMultipoleComparisonAccuracySeed &&
                     excitation.fastMultipoleComparisonAccuracySeed <
                         excitation.fastMultipoleAccuracySeed &&
                     excitation.fastMultipoleComparisonLevels !=
                         excitation.fastMultipoleVerificationLevels) ||
                        (!higherOrder && excitation.fastMultipoleComparisonAccuracySeed == 0 &&
                         excitation.fastMultipoleVerificationAccuracySeed == 0 &&
                         excitation.fastMultipoleComparisonLevels.empty() &&
                         excitation.fastMultipoleVerificationLevels.empty() &&
                         excitation.fastMultipoleOperatorDiscrepancy == 0),
                    "Wideband audit seeds or actual independent level evidence are inconsistent.");
            require(std::isfinite(excitation.fastMultipoleKrylovResidual) &&
                        std::isfinite(excitation.fastMultipoleOperatorDiscrepancy) &&
                        excitation.fastMultipoleOperatorDiscrepancy <=
                            fastMultipoleOrderComparisonBudgetFraction *
                                request.numerics.relativeResidualTolerance,
                    "FMM accepted an unconverged interpolation-order comparison.");
        }
}

void checkFastExecution(const BemResult &result) {
    require(result.acceleration == BemAcceleration::fastMultipole,
            "Requested FMM was silently replaced by a dense solve.");
    for (const auto &frequency : result.frequencies) {
        const auto &stats = frequency.fastMultipole;
        require(stats.nodeCount > 1 && stats.leafCount > 1 && stats.farInteractionCount > 0 &&
                    stats.pressurePressureFarPairs > 0,
                "FMM parity fixture did not execute multipole translations.");
        require(stats.pressurePressureNearPairs + stats.pressurePressureFarPairs ==
                        result.pressureDofCount * result.pressureDofCount &&
                    stats.pressureFluxNearPairs + stats.pressureFluxFarPairs ==
                        result.pressureDofCount * result.triangleCount,
                "Wideband complete-support interaction partition loses or duplicates DOF pairs.");
        require(!stats.historicalChebyshevStatistics,
                "New wideband execution retained historical panel-pair provenance.");
        for (const auto &excitation : frequency.excitations)
            require(excitation.iterations > 0, "FMM did not report actual Krylov iterations.");
        std::cout << frequency.frequency << " Hz: M2L " << stats.farInteractionCount
                  << ", P1/P1 near/far " << stats.pressurePressureNearPairs << '/'
                  << stats.pressurePressureFarPairs << ", P1/DP0 near/far "
                  << stats.pressureFluxNearPairs << '/' << stats.pressureFluxFarPairs << '\n';
    }
}

Complex analytic(double frequency, double radius) {
    const double k = 2 * std::numbers::pi * frequency / 343.0;
    return 1.21 * 343.0 * (imaginary * k * .1) / (1.0 + imaginary * k * .1) * (.1 / radius) *
           std::exp(-imaginary * k * (radius - .1));
}

void pinnedReference(const BemResult &result) {
    std::ifstream input(HS_PARITY_REFERENCE_PATH);
    require(static_cast<bool>(input), "Pinned BEAT reference is missing.");
    std::string line;
    std::size_t f = 0;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#')
            continue;
        std::istringstream row(line);
        row.imbue(std::locale::classic());
        double frequency{};
        row >> frequency;
        require(f < result.frequencies.size() && result.frequencies[f].frequency == frequency,
                "Pinned BEAT fixture frequency mismatch.");
        std::array<Complex, 3> expected{};
        for (auto &value : expected) {
            char first{}, second{};
            double real{}, imag{};
            row >> first >> real >> second >> imag;
            require(static_cast<bool>(row) && first == ',' && second == ',',
                    "Malformed pinned BEAT fixture.");
            value = {real, imag};
        }
        const auto &actual = result.frequencies[f++].excitations.front().observerPressure;
        const double error = relativeError(actual, expected);
        require(error <= .01, "FMM failed the pinned BEAT complex L2 gate.");
        for (std::size_t p = 0; p < expected.size(); ++p) {
            const auto ratio = actual[p] / expected[p];
            require(std::abs(20 * std::log10(std::abs(ratio))) <= .1 &&
                        std::abs(std::arg(ratio)) * 180 / std::numbers::pi <= 1,
                    "FMM failed the pinned BEAT level or phase gate.");
        }
        std::cout << "Pinned BEAT, " << frequency << " Hz: complex L2 " << error << '\n';
    }
    require(f == result.frequencies.size(), "Pinned BEAT fixture is incomplete.");
}

void checkSphere() {
    CudaBemSolver solver(testing::sphere(3));
    const auto input = sphereRequest(), fastInput = accelerated(input);
    const auto dense = solver.solve(input);
    const auto fast = solver.solve(fastInput);
    checkFastExecution(fast);
    compare(fast, dense, parityLimit);
    checkResidual(fast, fastInput);
    pinnedReference(fast);
    for (const auto &frequency : fast.frequencies) {
        const auto &unit = frequency.excitations.front();
        for (std::size_t p = 0; p < 3; ++p) {
            require(std::abs(unit.observerPressure[p] /
                                 analytic(frequency.frequency, std::array{.3, 1.0, 3.0}[p]) -
                             1.0) < sphereErrorLimit,
                    "FMM failed the independent pulsating-sphere formula.");
            require(
                std::abs(frequency.excitations[1].observerPressure[p] / unit.observerPressure[p] -
                         std::polar(.6, .71)) < 5e-4,
                "FMM lost complex source amplitude or phase.");
        }
        require(unit.impedances.front().volumeImpedance.real() > 0,
                "Pulsating sphere radiation resistance must be positive.");
        require(std::abs(unit.impedances.front().meanPressure / analytic(frequency.frequency, .1) -
                         1.0) < 2 * sphereErrorLimit,
                "FMM failed the independent sphere surface impedance.");
    }
    // A reused solver must change its actual execution path when switched off.
    auto again = input;
    again.frequencies = {500};
    BemResult reference = dense;
    reference.frequencies = {dense.frequencies[1]};
    compare(solver.solve(again), reference, 1e-5);
    std::cout << "512T peak device bytes: dense " << dense.peakDeviceBytes << ", FMM "
              << fast.peakDeviceBytes << " (small meshes include FMM setup overhead).\n";
}

BemRequest independentRequest(SurfaceMesh &mesh) {
    BemRequest request;
    request.frequencies = {500, 1715, 2452.940460658903};
    request.observers = {{0, 0, .3}, {0, 0, -.3}, {.3, 0, 0}, {.3, .4, .5}};
    request.excitations = {{"radial", {}}, {"translation-z", {}}, {"mixed", {}}};
    for (std::size_t panel = 0; panel < mesh.triangles.size(); ++panel) {
        auto &triangle = mesh.triangles[panel];
        triangle.faceId = static_cast<std::uint32_t>(panel + 1);
        const auto a = mesh.vertices[triangle.vertices[0]], b = mesh.vertices[triangle.vertices[1]],
                   c = mesh.vertices[triangle.vertices[2]];
        const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
        const Vec3 normal{u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
        const double nz = normal.z / std::hypot(normal.x, normal.y, normal.z);
        request.excitations[0].faces.push_back({triangle.faceId, 1.0});
        request.excitations[1].faces.push_back({triangle.faceId, nz});
        request.excitations[2].faces.push_back({triangle.faceId, 1.0 + mixture * nz});
    }
    return request;
}

void checkIndependentRhs() {
    auto mesh = testing::sphere(3);
    const auto input = independentRequest(mesh), fastInput = accelerated(input);
    CudaBemSolver solver(std::move(mesh));
    const auto dense = solver.solve(input), fast = solver.solve(fastInput);
    checkFastExecution(fast);
    compare(fast, dense, parityLimit);
    checkResidual(fast, fastInput);
    for (const auto &frequency : fast.frequencies) {
        const auto &radial = frequency.excitations[0].observerPressure,
                   &dipole = frequency.excitations[1].observerPressure,
                   &mixed = frequency.excitations[2].observerPressure;
        std::vector<Complex> expected;
        for (std::size_t p = 0; p < radial.size(); ++p)
            expected.push_back(radial[p] + mixture * dipole[p]);
        require(relativeError(mixed, expected) < 5e-4,
                "Independent FMM right-hand sides violate superposition.");
        require(std::abs(dipole[2]) < 5e-4 * std::abs(dipole[0]) &&
                    std::abs(dipole[0] + dipole[1]) < 5e-4 * std::abs(dipole[0]),
                "FMM breaks the dipole null or antisymmetry.");
    }
}

BemRequest irregularRequest(SurfaceMesh &mesh) {
    auto request = independentRequest(mesh);
    request.frequencies = {913};
    request.excitations.resize(1);
    request.excitations.front().name = "deterministic complex spatial drive";
    for (std::size_t panel = 0; panel < mesh.triangles.size(); ++panel) {
        const double argument = static_cast<double>(panel + 1);
        request.excitations.front().faces[panel].normalVelocity = {std::sin(argument * .417) + .2,
                                                                   std::cos(argument * .731)};
    }
    return request;
}

void checkSparseRhs() {
    auto mesh = testing::sphere(3);
    for (auto &triangle : mesh.triangles) {
        double z = 0;
        for (const auto vertex : triangle.vertices)
            z += mesh.vertices[vertex].z / 3;
        triangle.faceId = z > .07 ? 2u : (z < -.07 ? 3u : 1u);
        triangle.boundary =
            triangle.faceId == 1 ? BoundaryCondition::wall : BoundaryCondition::driven;
    }
    BemRequest input;
    input.frequencies = {500};
    input.observers = {{0, 0, .3}, {0, 0, -.3}, {.3, .4, .5}};
    // Omitted ports and explicit signed zeros must both remain silent. A tiny
    // nonzero source guards against replacing the exact-zero test by a cutoff.
    input.excitations = {{"north", {{2, 1.0}, {3, {0.0, -0.0}}}},
                         {"south", {{3, mixture}}},
                         {"combined", {{2, 1.0}, {3, mixture}}},
                         {"tiny", {{2, {1e-12, -1e-12}}}}};
    CudaBemSolver solver(std::move(mesh));
    const auto fastInput = accelerated(input);
    const auto dense = solver.solve(input), fast = solver.solve(fastInput);
    checkFastExecution(fast);
    compare(fast, dense, parityLimit);
    checkResidual(fast, fastInput);
    const auto &results = fast.frequencies.front().excitations;
    std::vector<Complex> sum, tiny;
    for (std::size_t i = 0; i < input.observers.size(); ++i) {
        sum.push_back(results[0].observerPressure[i] + results[1].observerPressure[i]);
        tiny.push_back(results[0].observerPressure[i] * Complex{1e-12, -1e-12});
    }
    require(relativeError(results[2].observerPressure, sum) < parityLimit,
            "Sparse mixed Wall/Driven excitations violate superposition.");
    require(relativeError(results[3].observerPressure, tiny) < parityLimit,
            "Exact-zero flux skipping changed a tiny nonzero source.");
}

void checkConvergence() {
    auto mesh = testing::sphere(3);
    const auto input = irregularRequest(mesh);
    CudaBemSolver solver(std::move(mesh));
    const auto dense = solver.solve(input);
    const auto lowInput = accelerated(input, 1e-3), highInput = accelerated(input, 1e-6);
    const auto low = solver.solve(lowInput), high = solver.solve(highInput);
    checkFastExecution(low);
    checkFastExecution(high);
    const double lowError = compare(low, dense, .05);
    const double highError = compare(high, dense, parityLimit);
    require(highError < .5 * lowError || highError < 1e-4,
            "Tightening wideband residual/seed budgets did not improve complex pressure parity.");
    checkResidual(high, highInput);
    std::cout << "Wideband convergence residual budgets 1e-3/1e-6: " << lowError << '/' << highError
              << '\n';
}

template <class Exception, class Function> void rejects(Function function, const char *message) {
    bool rejected = false;
    try {
        function();
    } catch (const Exception &) {
        rejected = true;
    }
    require(rejected, message);
}

void checkRunningCancellation();

void checkLifecycle() {
    auto mesh = testing::sphere(3);
    const auto input = accelerated(irregularRequest(mesh));
    CudaBemSolver solver(mesh);
    auto invalid = input;
    invalid.numerics.fastMultipole.separationRatio = 1;
    rejects<std::invalid_argument>([&] { solver.solve(invalid); },
                                   "FMM accepted an inadmissible separation ratio.");
    invalid = input;
    invalid.numerics.fastMultipole.leafSupportCapacity = 0;
    rejects<std::invalid_argument>([&] { solver.solve(invalid); },
                                   "FMM accepted zero leaf capacity.");
    invalid = input;
    invalid.numerics.fastMultipole.maximumIterations = 1;
    invalid.numerics.fastMultipole.restartLength = 2;
    bool unconverged = false;
    try {
        solver.solve(invalid);
    } catch (const std::runtime_error &error) {
        const std::string message = error.what();
        require(message.find("1/1 iterations") != std::string::npos &&
                    message.find("unpreconditioned relative FMM equation residual") !=
                        std::string::npos &&
                    message.find("preconditioned relative residual") != std::string::npos,
                "GMRES failure omitted its actual residual or iteration budget.");
        unconverged = true;
    }
    require(unconverged, "Unconverged FMM solve returned successful results.");
    std::stop_source cancelled;
    cancelled.request_stop();
    rejects<std::runtime_error>([&] { solver.solve(input, cancelled.get_token()); },
                                "Pre-cancelled FMM solve returned successful results.");
    checkRunningCancellation();
    auto sweep = input;
    sweep.frequencies = {500, 1000};
    std::stop_source during;
    unsigned completed = 0;
    rejects<std::runtime_error>(
        [&] {
            solver.solve(sweep, during.get_token(), [&](const BemProgress &progress) {
                if (progress.phase == BemSolvePhase::frequencyComplete) {
                    completed = static_cast<unsigned>(progress.completedFrequencies);
                    during.request_stop();
                }
            });
        },
        "FMM ignored cancellation between sweep frequencies.");
    require(completed == 1, "FMM published a second frequency after cancellation.");
    auto high = input;
    high.frequencies = {1500};
    const auto reuseSpace = galerkin::buildSpace(mesh);
    double maximumKh = 0;
    for (const auto &panel : reuseSpace.panels)
        maximumKh = std::max(maximumKh, 2 * std::numbers::pi * high.frequencies.front() /
                                            high.medium.soundSpeed * panel.maximumEdgeLength);
    require(maximumKh <= 2 * std::numbers::pi / 6,
            "Lifecycle reuse fixture needs at least six panel edges per wavelength.");
    auto direct = high;
    direct.numerics.acceleration = BemAcceleration::dense;
    // Exercise a real frequency/plan and backend change on a resolved mesh.
    // High-frequency algebra is covered by the independent wideband gates.
    const auto denseBefore = solver.solve(direct);
    const auto acceleratedResult = solver.solve(high);
    compare(acceleratedResult, denseBefore, parityLimit);
    const auto denseAfter = solver.solve(direct);
    compare(denseAfter, denseBefore, parityLimit);
    constexpr double offset = 1e9;
    for (auto &point : mesh.vertices) {
        point.x += offset;
        point.y -= offset;
        point.z += offset;
    }
    auto localMesh = mesh;
    for (auto &point : localMesh.vertices) {
        point.x -= offset;
        point.y += offset;
        point.z -= offset;
    }
    auto world = input;
    for (auto &point : world.observers) {
        point.x += offset;
        point.y -= offset;
        point.z += offset;
    }
    auto local = world;
    for (auto &point : local.observers) {
        point.x -= offset;
        point.y += offset;
        point.z -= offset;
    }
    CudaBemSolver worldSolver(std::move(mesh)), localSolver(std::move(localMesh));
    compare(worldSolver.solve(world), localSolver.solve(local), 2e-4);
}

void cudaCheck(cudaError_t status) {
    if (status != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(status));
}

void checkRunningCancellation() {
    using namespace std::chrono;
    auto mesh = testing::sphere(4);
    auto input = accelerated(irregularRequest(mesh));
    input.frequencies = {500};
    CudaBemSolver solver(std::move(mesh));
    int device{};
    cudaCheck(cudaGetDevice(&device));
    std::size_t baseline{}, total{};
    cudaCheck(cudaMemGetInfo(&baseline, &total));
    std::stop_source stop;
    std::promise<void> ready;
    auto entered = ready.get_future();
    std::atomic<unsigned> published{};
    auto worker = std::async(std::launch::async, [&] {
        ready.set_value();
        try {
            cudaCheck(cudaSetDevice(device));
            solver.solve(input, stop.get_token(), [&](const BemProgress &progress) {
                if (progress.phase == BemSolvePhase::frequencyComplete)
                    ++published;
            });
            return std::string{"FMM returned a completed result after cancellation."};
        } catch (const std::runtime_error &error) {
            const std::string message = error.what();
            return message.find("cancelled") != std::string::npos ? std::string{} : message;
        }
    });
    entered.get();
    // The solver is constructed before this baseline. New device allocations
    // therefore witness the first solve entering its CUDA preparation, without
    // a product-only progress hook or a timing assumption about CPU startup.
    bool deviceWorkStarted = false;
    const auto startDeadline = steady_clock::now() + seconds(5);
    while (steady_clock::now() < startDeadline &&
           worker.wait_for(milliseconds(0)) != std::future_status::ready) {
        std::size_t available{};
        cudaCheck(cudaMemGetInfo(&available, &total));
        if (available < baseline && baseline - available >= 1024 * 1024) {
            deviceWorkStarted = true;
            break;
        }
        std::this_thread::sleep_for(milliseconds(2));
    }
    // Allow the first construction kernels to start after their allocations.
    std::this_thread::sleep_for(milliseconds(25));
    const bool stillRunning = worker.wait_for(milliseconds(0)) != std::future_status::ready;
    const auto cancelStart = steady_clock::now();
    stop.request_stop();
    const auto outcome = worker.get();
    const auto cancellationSeconds = duration<double>(steady_clock::now() - cancelStart).count();
    require(deviceWorkStarted && stillRunning,
            "Cancellation fixture did not reach an active first CUDA FMM solve.");
    require(outcome.empty(), "Running FMM cancellation failed: " + outcome);
    require(published.load() == 0, "Cancelled first FMM solve published a frequency result.");
    require(cancellationSeconds < 15,
            "Running FMM solve exceeded the 15-second cancellation return budget.");
    std::cout << "First-solve CUDA cancellation returned in " << cancellationSeconds << " s\n";
}

// Detailed pressure/flux actions, high-order rotations, LF/HF handoffs,
// exact spatial preconditioning and independent direct solve/audit checks live
// in the corresponding wideband/support CUDA tests. This public regression
// retains asymmetric geometry and arbitrary complex face excitation parity.
void checkOperator() {
    auto mesh = testing::sphere(2);
    for (auto &point : mesh.vertices)
        point = {1.27 * point.x + .13 * point.z, .71 * point.y, .95 * point.z};
    const auto request = irregularRequest(mesh);
    CudaBemSolver solver(mesh);
    const auto direct = solver.solve(request);
    const auto fastRequest = accelerated(request);
    const auto fast = solver.solve(fastRequest);
    checkFastExecution(fast);
    checkResidual(fast, fastRequest);
    compare(fast, direct, parityLimit);
    auto silent = fastRequest;
    for (auto &face : silent.excitations.front().faces)
        face.normalVelocity = {};
    const auto zero = solver.solve(silent);
    for (const auto &frequency : zero.frequencies)
        for (const auto &excitation : frequency.excitations) {
            require(excitation.boundaryResidual == 0 && excitation.iterations == 0 &&
                        excitation.fastMultipoleKrylovResidual == 0 &&
                        excitation.fastMultipoleOperatorDiscrepancy == 0,
                    "Zero-flow public wideband solve invented a residual or iteration.");
            for (const auto pressure : excitation.boundary()->surfacePressure)
                require(pressure == Complex{},
                        "Zero-flow public wideband solve retained earlier pressure.");
        }
    compare(solver.solve(fastRequest), fast, parityLimit);
}

void checkLem() {
    auto mesh = testing::sphere(3);
    for (auto &triangle : mesh.triangles) {
        const double z = mesh.vertices[triangle.vertices[0]].z +
                         mesh.vertices[triangle.vertices[1]].z +
                         mesh.vertices[triangle.vertices[2]].z;
        triangle.faceId = z > 0 ? 1U : 2U;
    }
    const auto space = galerkin::buildSpace(mesh);
    std::array<double, 2> areas{};
    for (const auto &panel : space.panels)
        areas[panel.faceId - 1] += panel.area;
    BemRequest basis;
    basis.frequencies = {300, 700};
    basis.observers = {{0, 0, .5}, {0, 0, -.7}, {.6, .2, 0}};
    basis.excitations = {{"North", {{1, 1.0 / areas[0]}}}, {"South", {{2, 1.0 / areas[1]}}}};
    CudaBemSolver solver(mesh);
    const auto basisResult = solver.solve(basis);
    const lem::DriverParameters driver{6, .0003, 1.5, .0003, .02, .05, .03};
    const std::array<Complex, 2> voltage{{{2, .3}, {4, -2}}};
    auto coupled = basis;
    coupled.excitations.clear();
    coupled.coupling = lem::makeDriverNetwork(driver, {{1}, {2}}, voltage[0], 1.0);
    coupled.coupling->nodes.push_back(
        {{100}, lem::Domain::electrical, false, "Second amplifier", {}});
    coupled.coupling->elements.push_back(
        {{100}, "Second voltage", {}, lem::EffortSource{{{100}, {1}}, voltage[1]}});
    std::size_t driverIndex = 0;
    for (auto &element : coupled.coupling->elements)
        if (auto *part = std::get_if<lem::Driver>(&element.component))
            if (driverIndex++ == 1)
                part->electrical.positive = {100};
    require(driverIndex == 2, "Coupled FMM fixture does not contain two drivers.");
    const auto dense = solver.solve(coupled);
    const auto fastInput = accelerated(coupled);
    const auto fast = solver.solve(fastInput);
    compare(fast, dense, parityLimit);
    checkFastExecution(fast);
    checkResidual(fast, fastInput);
    const auto close = [](Complex actual, Complex expected, double tolerance, const char *message) {
        require(std::abs(actual - expected) <= tolerance * std::max(1e-15, std::abs(expected)),
                message);
    };
    for (std::size_t f = 0; f < fast.frequencies.size(); ++f) {
        const auto &output = fast.frequencies[f].excitations.front();
        require(output.lemSolution.has_value(), "FMM coupled output lost its LEM solution.");
        const auto &solution = *output.lemSolution;
        require(solution.drivers.size() == 2 && solution.ports.size() == 2 &&
                    solution.relativeResidual <= lem::maximumRelativeResidual,
                "FMM coupled output lost driver/port data or violates its MNA equations.");
        std::array<Complex, 4> impedance{};
        for (std::size_t column = 0; column < 2; ++column)
            for (std::size_t p = 0; p < space.panels.size(); ++p) {
                const auto row = space.panels[p].faceId - 1;
                impedance[2 * row + column] +=
                    basisResult.frequencies[f].excitations[column].boundary()->surfacePressure[p] *
                    space.panels[p].area / areas[row];
            }
        const double omega = 2 * std::numbers::pi * fast.frequencies[f].frequency;
        const Complex ze = driver.resistance + imaginary * omega * driver.inductance;
        const Complex zm =
            driver.mechanicalResistance + imaginary * omega * driver.movingMass +
            1.0 / (imaginary * omega * driver.compliance) +
            driver.pistonArea * driver.pistonArea * 1.21 * 343 * 343 / (imaginary * omega);
        const Complex unloaded = zm + driver.forceFactor * driver.forceFactor / ze;
        const double areaSquared = driver.pistonArea * driver.pistonArea;
        const Complex a = unloaded + areaSquared * impedance[0], b = areaSquared * impedance[1],
                      c = areaSquared * impedance[2], d = unloaded + areaSquared * impedance[3],
                      force0 = driver.forceFactor * voltage[0] / ze,
                      force1 = driver.forceFactor * voltage[1] / ze;
        const std::array<Complex, 2> expectedVelocity{(d * force0 - b * force1) / (a * d - b * c),
                                                      (a * force1 - c * force0) / (a * d - b * c)};
        if (f == 0)
            require(std::abs(expectedVelocity[0] / (force0 / a) - 1.0) > .05,
                    "FMM LEM fixture does not detect omitted mutual acoustic loading.");
        double driverAcousticPower = 0, portPower = 0;
        for (std::size_t port = 0; port < 2; ++port) {
            const auto &actual = solution.drivers[port];
            close(actual.velocity, expectedVelocity[port], parityLimit,
                  "FMM loaded diaphragm velocity fails independent two-port equations.");
            close(actual.current,
                  (voltage[port] - driver.forceFactor * expectedVelocity[port]) / ze, parityLimit,
                  "FMM driver current loses back-EMF coupling.");
            close(actual.displacement, actual.velocity / (imaginary * omega), 1e-12,
                  "FMM driver displacement uses a different phasor convention.");
            close(solution.ports[port].volumeVelocity, driver.pistonArea * actual.velocity, 1e-11,
                  "FMM diaphragm area transformer changes the volume flow.");
            close(actual.electricalPower,
                  actual.copperPower + actual.mechanicalPower + actual.acousticPower, 1e-9,
                  "FMM coupled driver violates active-power balance.");
            Complex mean{};
            for (std::size_t p = 0; p < space.panels.size(); ++p)
                if (space.panels[p].faceId == port + 1)
                    mean += output.boundary()->surfacePressure[p] * space.panels[p].area / areas[port];
            close(mean, solution.ports[port].pressure, parityLimit,
                  "FMM reconstructed boundary pressure fails the LEM feedback loop.");
            close(solution.ports[port].activePower,
                  std::real(solution.ports[port].pressure *
                            std::conj(solution.ports[port].volumeVelocity)),
                  1e-11, "FMM acoustic port power is inconsistent with RMS pressure and flow.");
            driverAcousticPower += actual.acousticPower;
            portPower += solution.ports[port].activePower;
        }
        close(driverAcousticPower, portPower, 1e-9,
              "FMM coupled network does not conserve acoustic active power.");
    }
    auto prescribed = accelerated(basis);
    prescribed.frequencies = {basis.frequencies.front()};
    prescribed.excitations = {
        {"Actual", fast.frequencies.front().excitations.front().solvedVelocities}};
    const auto actualVelocityResult = solver.solve(prescribed);
    require(
        relativeError(actualVelocityResult.frequencies.front().excitations.front().observerPressure,
                      fast.frequencies.front().excitations.front().observerPressure) < parityLimit,
        "FMM coupled field differs from a direct actual-velocity solve.");
    auto zero = fastInput;
    zero.frequencies = {300};
    for (auto &element : zero.coupling->elements)
        if (auto *source = std::get_if<lem::EffortSource>(&element.component))
            source->amplitude = {};
    const auto zeroResult = solver.solve(zero);
    const auto &silence = zeroResult.frequencies.front().excitations.front();
    require(silence.boundaryResidual == 0 && silence.impedances.empty() &&
                silence.solvedVelocities.size() == 2 && silence.lemSolution.has_value(),
            "Zero-drive FMM coupling lost its finite zero-flow metadata.");
    for (const auto pressure : silence.boundary()->surfacePressure)
        require(pressure == Complex{}, "Zero-drive FMM surface pressure is nonzero.");
    for (const auto pressure : silence.observerPressure)
        require(pressure == Complex{}, "Zero-drive FMM field is nonzero.");
}
} // namespace

int main(int argc, char **argv) {
    try {
        std::cout << std::setprecision(10) << "CUDA device: " << CudaBemSolver::deviceDescription()
                  << '\n';
        const std::string_view mode = argc == 2 ? argv[1] : "--sphere";
        if (mode == "--sphere")
            checkSphere();
        else if (mode == "--independent-rhs")
            checkIndependentRhs();
        else if (mode == "--sparse-rhs")
            checkSparseRhs();
        else if (mode == "--convergence")
            checkConvergence();
        else if (mode == "--lifecycle")
            checkLifecycle();
        else if (mode == "--operator")
            checkOperator();
        else if (mode == "--lem")
            checkLem();
        else
            throw std::invalid_argument("Unknown FMM test mode.");
        std::cout << "Fast-multipole BEM " << mode << " passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
