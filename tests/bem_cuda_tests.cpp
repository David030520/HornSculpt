#include "cuda/workspace_size.hpp"
#include "fixture_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cuda_runtime_api.h>
#include <fstream>
#include <galerkin/interaction_plan.hpp>
#include <galerkin/space.hpp>
#include <hornsculpt/bem.hpp>
#include <hornsculpt/draft_geometry.hpp>
#include <hornsculpt/meshing.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string_view>

using namespace hornsculpt;
using hornsculpt::testing::sphere;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
void checkMemoryPreflight() {
    require(cudaFree(nullptr) == cudaSuccess, "CUDA context initialization failed.");
    std::size_t free{}, total{};
    require(cudaMemGetInfo(&free, &total) == cudaSuccess && total > 0,
            "CUDA device capacity query failed.");
    constexpr std::size_t rhsCount = 1024, maximumObservers = 100000;
    SurfaceMesh mesh;
    for (unsigned refinement = 2; refinement <= 8; ++refinement) {
        mesh = sphere(refinement);
        if (cuda_detail::denseWorkspaceSize<std::complex<float>, std::complex<double>>(
                mesh.vertices.size(), mesh.triangles.size(), rhsCount, maximumObservers)
                .minimumBytes > total)
            break;
    }
    const auto needed = [&](std::size_t observers) {
        return cuda_detail::denseWorkspaceSize<std::complex<float>, std::complex<double>>(
                   mesh.vertices.size(), mesh.triangles.size(), rhsCount, observers)
            .minimumBytes;
    };
    require(needed(maximumObservers) > total,
            "The test fixture cannot exceed this device's physical capacity safely.");
    std::size_t low = 0, high = maximumObservers;
    while (low < high) {
        const auto middle = low + (high - low) / 2;
        if (needed(middle) > total)
            high = middle;
        else
            low = middle + 1;
    }
    const auto minimumBytes = needed(low);
    BemRequest input;
    input.frequencies = {1000};
    input.farField = FarFieldObservation{};
    input.observers.assign(low, {0, 0, 1});
    input.excitations.assign(rhsCount, {"Capacity check", {{1, {1, 0}}}});
    CudaBemSolver solver(std::move(mesh));
    bool rejectedBeforeQuadrature = false;
    try {
        (void)solver.solve(input);
    } catch (const std::runtime_error &error) {
        rejectedBeforeQuadrature =
            std::string_view(error.what()).find("GPU-Gesamtkapazität") != std::string_view::npos;
    }
    require(rejectedBeforeQuadrature,
            "An impossible dense workspace was not rejected before near-pair planning.");
    require(cudaMemGetInfo(&free, &total) == cudaSuccess,
            "CUDA device query failed after capacity rejection.");
    std::cout << "CUDA capacity preflight: minimum " << minimumBytes << " bytes > total " << total
              << " bytes; rejected before near-pair planning without a dense allocation.\n";
}
void checkCoarseDraft() {
    DraftParameters parameters;
    parameters.entry.width = parameters.entry.height = .004;
    parameters.length = .020;
    parameters.wallThickness = .001;
    parameters.horizontalAngle = 20 * std::numbers::pi / 180;
    parameters.verticalAngle = 15 * std::numbers::pi / 180;
    MeshingOptions options;
    options.targetEdgeLength = .008;
    options.quality.enabled = false;
    const auto source = DraftGeometry(parameters).previewMesh(options.adaptive.geometricTolerance);
    auto remeshed = remeshIsotropic(source, options);
    const auto &quality = remeshed.quality;
    require(quality.watertight() && quality.selfIntersectionPairs == 0 &&
                quality.connectedComponents == 1 && quality.eulerCharacteristic == 2 &&
                quality.signedVolume > 0,
            "Coarse Draft regression did not produce a valid closed exterior.");
    const auto triangleCount = remeshed.mesh.triangles.size();
    const auto vertexCount = remeshed.mesh.vertices.size();
    BemRequest input;
    input.frequencies = {100};
    const auto driven = static_cast<std::uint32_t>(DraftFace::driven);
    input.excitations = {{"Uniform entry", {{driven, {1, 0}}}}};
    input.observers = {{0, 0, parameters.length + 1}};
    input.farField = FarFieldObservation{{0, 0, parameters.length}, 1};
    CudaBemSolver solver(std::move(remeshed.mesh));
    const auto result = solver.solve(input);
    require(result.triangleCount == triangleCount && result.pressureDofCount == vertexCount &&
                result.frequencies.size() == 1 && result.frequencies.front().frequency == 100 &&
                result.frequencies.front().excitations.size() == 1,
            "Coarse Draft solve changed its pressure space or request.");
    const auto &excitation = result.frequencies.front().excitations.front();
    require(std::isfinite(excitation.boundaryResidual) &&
                excitation.boundaryResidual <= input.numerics.relativeResidualTolerance &&
                excitation.boundary()->surfacePressure.size() == triangleCount &&
                excitation.observerPressure.size() == 1 && excitation.impedances.size() == 1,
            "Coarse Draft solve failed the unchanged residual or result dimensions.");
    for (const auto pressure : excitation.boundary()->surfacePressure)
        require(std::isfinite(pressure.real()) && std::isfinite(pressure.imag()),
                "Coarse Draft surface pressure is not finite.");
    const auto &impedance = excitation.impedances.front();
    const auto pressure = excitation.observerPressure.front();
    require(impedance.faceId == driven && impedance.area > 0 &&
                std::isfinite(impedance.volumeImpedance.real()) &&
                std::isfinite(impedance.volumeImpedance.imag()) && std::isfinite(pressure.real()) &&
                std::isfinite(pressure.imag()) &&
                std::abs(impedance.volumeVelocity - std::complex<double>(impedance.area, 0)) <=
                    1e-12 * impedance.area,
            "Coarse Draft far-field pressure or explicit driven flow is invalid.");
    std::cout << "Coarse Draft: " << triangleCount << " triangles, minimum angle "
              << quality.minimumAngleDegrees << ", aspect " << quality.maximumAspectRatio
              << ", 100 Hz residual " << excitation.boundaryResidual
              << ", unchanged near tolerance/depth " << input.numerics.nearRelativeTolerance << '/'
              << input.numerics.nearMaximumDepth << " PASS.\n";
}
std::complex<double> analytic(double frequency, double radius) {
    const double k = 2 * std::numbers::pi * frequency / 343.0;
    const std::complex<double> i{0, 1};
    return 1.21 * 343.0 * (i * k * .1) / (1.0 + i * k * .1) * (.1 / radius) *
           std::exp(-i * k * (radius - .1));
}
BemRequest request() {
    BemRequest value;
    value.frequencies = {100, 500, 1000, 1715};
    value.excitations = {{"unit", {{1, {1, 0}}}}, {"phase", {{1, std::polar(.6, .71)}}}};
    value.observers = {{0, 0, .3}, {1, 0, 0}, {0, 0, 3}};
    return value;
}
double checkSphere(const BemResult &result, double maximumError,
                   double residualTolerance = BemNumerics{}.relativeResidualTolerance) {
    require(result.pressureDofCount == result.triangleCount / 2 + 2,
            "Closed spherical P1 pressure space has an incorrect number of DOFs");
    double worst = 0;
    for (const auto &frequency : result.frequencies) {
        const auto &unit = frequency.excitations[0];
        for (const auto &excitation : frequency.excitations) {
            require(excitation.fieldQuadratureComparisonDegree > 0 &&
                        excitation.fieldQuadratureVerificationDegree >
                            excitation.fieldQuadratureComparisonDegree &&
                        excitation.fieldQuadratureVerificationDegree <= 14 &&
                        std::isfinite(excitation.fieldQuadratureDiscrepancy) &&
                        excitation.fieldQuadratureDiscrepancy <= .025 * residualTolerance,
                    "Dense observer field lacks a genuine independent rule comparison.");
            require(excitation.fastMultipoleAccuracySeed == 0 &&
                        excitation.fastMultipoleComparisonLevels.empty() &&
                        excitation.fastMultipoleVerificationLevels.empty() &&
                        excitation.farQuadratureComparisonDegree == 0,
                    "Dense field quadrature must not invent FMM expansion evidence.");
        }
        for (std::size_t j = 0; j < 3; ++j) {
            const auto expected = analytic(frequency.frequency, std::array{.3, 1.0, 3.0}[j]);
            const double error = std::abs(unit.observerPressure[j] / expected - 1.0);
            worst = std::max(worst, error);
            std::cout << result.triangleCount << " triangles, " << frequency.frequency
                      << " Hz, observer " << j << ", complex relative error " << error
                      << ", refinement " << unit.refinementSteps << ", boundary "
                      << unit.boundaryResidual << ", assembly/solve/field ms "
                      << frequency.assemblyMilliseconds << " / " << frequency.solveMilliseconds
                      << " / " << frequency.fieldMilliseconds << '\n';
            require(error < maximumError, "Pulsating sphere pressure fails the analytic reference");
            require(
                std::abs(frequency.excitations[1].observerPressure[j] / unit.observerPressure[j] -
                         std::polar(.6, .71)) < 3e-4,
                "Complex excitation amplitude or phase is inconsistent");
        }
        const auto &impedance = unit.impedances.front();
        require(impedance.volumeImpedance.real() > 0, "Radiation resistance must be positive");
        require(std::abs(impedance.meanPressure / analytic(frequency.frequency, .1) - 1.0) <
                    maximumError * 2,
                "Complex specific sphere impedance fails the analytic reference");
        require(std::abs(impedance.meanPressure / impedance.volumeVelocity -
                         impedance.volumeImpedance) < 1e-8,
                "Volume impedance normalization is inconsistent");
    }
    return worst;
}
void checkMappedQuadrature() {
    CudaBemSolver device(sphere(3), {CudaCacheResidency::device});
    CudaBemSolver mapped(sphere(3), {CudaCacheResidency::mappedHost});
    auto input = request();
    // Reuse the registration, then replace its source plan with a higher maximum k.
    for (const auto frequencies : {std::vector<double>{500, 1000}, std::vector<double>{500, 1000},
                                   std::vector<double>{750, 1250}}) {
        input.frequencies = frequencies;
        const auto expected = device.solve(input), actual = mapped.solve(input);
        require(actual.peakDeviceBytes < expected.peakDeviceBytes,
                "Mapped quadrature did not remove its nonempty device patch allocation.");
        for (std::size_t f = 0; f < expected.frequencies.size(); ++f)
            for (std::size_t e = 0; e < input.excitations.size(); ++e) {
                const auto &a = actual.frequencies[f].excitations[e];
                const auto &b = expected.frequencies[f].excitations[e];
                require(a.boundaryResidual <= input.numerics.relativeResidualTolerance,
                        "Mapped quadrature violated the original residual contract.");
                double error = 0, scale = 0;
                for (std::size_t j = 0; j < a.boundary()->surfacePressure.size(); ++j) {
                    error += std::norm(a.boundary()->surfacePressure[j] -
                                       b.boundary()->surfacePressure[j]);
                    scale += std::norm(b.boundary()->surfacePressure[j]);
                }
                require(std::sqrt(error / scale) < 2e-5,
                        "Mapped patch storage changed the boundary-pressure solution.");
                for (std::size_t j = 0; j < a.observerPressure.size(); ++j)
                    require(std::abs(a.observerPressure[j] / b.observerPressure[j] - 1.) < 2e-5,
                            "Mapped patch storage changed the complex observer field.");
                require(std::abs(a.impedances.front().volumeImpedance /
                                     b.impedances.front().volumeImpedance -
                                 1.) < 2e-5,
                        "Mapped patch storage changed the complex impedance.");
            }
    }
    std::cout << "Mapped/device CUDA quadrature parity, multiple excitations, cache reuse and "
                 "replacement passed.\n";
}
void checkCompactQuadrature() {
    const auto mesh = sphere(2);
    auto input = request();
    input.frequencies = {2500};
    galerkin::PairQuadratureSettings settings;
    settings.relativeTolerance = input.numerics.nearRelativeTolerance;
    settings.maximumSubdivisionDepth = input.numerics.nearMaximumDepth;
    const double k = 2 * std::numbers::pi * input.frequencies.front() / input.medium.soundSpeed;
    const auto plan = galerkin::makeInteractionPlan(galerkin::buildSpace(mesh), settings, k);
    const auto compactCount =
        std::count_if(plan.wholePanelDegrees.begin(), plan.wholePanelDegrees.end(),
                      [](auto d) { return d != 0; });
    require(compactCount > 0, "CUDA regression does not exercise compact full-panel rules.");
    CudaBemSolver solver(mesh);
    const auto actual = solver.solve(input);
    auto strict = input;
    strict.numerics.nearRelativeTolerance = 1e-6;
    const auto reference = solver.solve(strict);
    double worst{};
    for (std::size_t e = 0; e < input.excitations.size(); ++e) {
        const auto &a = actual.frequencies.front().excitations[e];
        const auto &b = reference.frequencies.front().excitations[e];
        require(a.boundaryResidual <= input.numerics.relativeResidualTolerance &&
                    b.boundaryResidual <= strict.numerics.relativeResidualTolerance,
                "Compact CUDA quadrature failed the unchanged solve residual.");
        for (std::size_t j = 0; j < a.observerPressure.size(); ++j)
            worst = std::max(worst, std::abs(a.observerPressure[j] / b.observerPressure[j] - 1.));
        worst = std::max(worst, std::abs(a.impedances.front().volumeImpedance /
                                             b.impedances.front().volumeImpedance -
                                         1.));
    }
    require(worst < 1e-4, "Compact CUDA rules disagree with stricter FP64 adaptive quadrature.");
    std::cout << "Compact CUDA rules: " << compactCount << " pairs, " << mesh.triangles.size()
              << " triangles, complex field/impedance difference " << worst << '\n';
}

void checkTightResidual() {
    // A single existing sphere fixture exercises the supported residual floor
    // without a frequency search or a geometry/discretization tolerance change.
    auto input = request();
    input.frequencies = {1715};
    input.numerics.relativeResidualTolerance = 1e-7;
    CudaBemSolver solver(sphere(4));
    const auto result = solver.solve(input);
    require(result.triangleCount == 2048 && result.pressureDofCount == 1026 &&
                result.frequencies.size() == 1 &&
                result.frequencies.front().excitations.size() == 2,
            "Tight residual regression changed the fixture or excitation dimensions.");
    for (const auto &excitation : result.frequencies.front().excitations) {
        require(std::isfinite(excitation.boundaryResidual) &&
                    excitation.boundaryResidual <= input.numerics.relativeResidualTolerance,
                "The requested tight residual was not met for every complex RHS.");
        require(excitation.refinementSteps <= input.numerics.maximumRefinementSteps,
                "Tight residual solve exceeded the unchanged refinement limit.");
        std::cout << "Tight residual 2048T/1026DOF, " << excitation.name << ": "
                  << std::setprecision(12) << excitation.boundaryResidual << ", refinement "
                  << excitation.refinementSteps << '\n';
    }
    checkSphere(result, .015, input.numerics.relativeResidualTolerance);
    std::cout << "Tight CUDA residual, analytic pressure, impedance and complex phase PASS.\n";
}
void checkThinBody() {
    SurfaceMesh mesh;
    mesh.vertices = {{0, 0, 0},   {.1, 0, 0},   {.1, .1, 0},   {0, .1, 0},
                     {0, 0, .01}, {.1, 0, .01}, {.1, .1, .01}, {0, .1, .01}};
    for (const auto vertices : std::array<std::array<std::uint32_t, 3>, 12>{{{0, 2, 1},
                                                                             {0, 3, 2},
                                                                             {4, 5, 6},
                                                                             {4, 6, 7},
                                                                             {0, 1, 5},
                                                                             {0, 5, 4},
                                                                             {1, 2, 6},
                                                                             {1, 6, 5},
                                                                             {2, 3, 7},
                                                                             {2, 7, 6},
                                                                             {3, 0, 4},
                                                                             {3, 4, 7}}})
        mesh.triangles.push_back({vertices, 1, BoundaryCondition::driven});
    BemRequest input;
    input.frequencies = {20, 1000};
    input.observers = {{.05, .05, 1.005}, {.05, .05, -.995}};
    input.farField = FarFieldObservation{{.05, .05, .005}, 1};
    input.excitations = {{"unit", {{1, {1, 0}}}}, {"phase", {{1, std::polar(.6, .71)}}}};
    CudaBemSolver solver(std::move(mesh));
    const auto result = solver.solve(input);
    auto strict = input;
    strict.numerics.nearRelativeTolerance = 1e-6;
    strict.numerics.nearMaximumDepth = 8;
    const auto reference = solver.solve(strict);
    double worst = 0;
    for (std::size_t f = 0; f < result.frequencies.size(); ++f) {
        const auto &frequency = result.frequencies[f];
        const auto &unit = frequency.excitations[0];
        require(std::abs(unit.observerPressure[0] / unit.observerPressure[1] - 1.) < .003,
                "Thin-body quadrature loses reflection symmetry");
        for (std::size_t e = 0; e < frequency.excitations.size(); ++e) {
            const auto &actual = frequency.excitations[e];
            const auto &expected = reference.frequencies[f].excitations[e];
            require(actual.boundaryResidual <= input.numerics.relativeResidualTolerance &&
                        actual.refinementSteps <= input.numerics.maximumRefinementSteps,
                    "Thin-body solve exceeds the unchanged residual/refinement limits");
            for (std::size_t j = 0; j < input.observers.size(); ++j) {
                const double difference =
                    std::abs(actual.observerPressure[j] / expected.observerPressure[j] - 1.);
                require(difference < 1e-4,
                        "Thin-body GPU field does not converge with stricter near quadrature");
                worst = std::max(worst, difference);
            }
            require(std::abs(actual.impedances.front().volumeImpedance /
                                 expected.impedances.front().volumeImpedance -
                             1.) < 1e-4,
                    "Thin-body GPU impedance does not converge with stricter near quadrature");
        }
        for (std::size_t j = 0; j < input.observers.size(); ++j)
            require(
                std::abs(frequency.excitations[1].observerPressure[j] / unit.observerPressure[j] -
                         std::polar(.6, .71)) < 3e-4,
                "Thin-body solve loses complex excitation linearity");
        std::cout << "Thin-body 100x100x10 mm: " << frequency.frequency
                  << " Hz, unchanged default near depth/tolerance 6/1e-4, residual "
                  << unit.boundaryResidual << '\n';
    }
    std::cout << "Thin-body GPU near-quadrature convergence, phase, impedance and symmetry PASS; "
                 "maximum complex field difference "
              << worst << '\n';
}
void exportReference(const char *path, const SurfaceMesh &mesh, const BemRequest &inputs,
                     const BemResult &results) {
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot write the external parity fixture.");
    out << std::setprecision(17)
        << "{\n\"phasor\":\"exp(+iwt)\",\n\"rho\":1.21,\"c\":343,\"driven_velocity\":[1,0],"
           "\n\"vertices_m\":[";
    for (std::size_t i = 0; i < mesh.vertices.size(); ++i) {
        const auto p = mesh.vertices[i];
        if (i)
            out << ',';
        out << '[' << p.x << ',' << p.y << ',' << p.z << ']';
    }
    out << "],\n\"triangles0\":[";
    for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
        const auto v = mesh.triangles[i].vertices;
        if (i)
            out << ',';
        out << '[' << v[0] << ',' << v[1] << ',' << v[2] << ']';
    }
    out << "],\n\"tags\":[";
    for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
        if (i)
            out << ',';
        out << mesh.triangles[i].faceId;
    }
    out << "],\n\"frequencies_hz\":[";
    for (std::size_t i = 0; i < inputs.frequencies.size(); ++i) {
        if (i)
            out << ',';
        out << inputs.frequencies[i];
    }
    out << "],\n\"observations_m\":[";
    for (std::size_t i = 0; i < inputs.observers.size(); ++i) {
        const auto p = inputs.observers[i];
        if (i)
            out << ',';
        out << '[' << p.x << ',' << p.y << ',' << p.z << ']';
    }
    out << "],\n\"hornsculpt_pressure\":[";
    for (std::size_t f = 0; f < results.frequencies.size(); ++f) {
        if (f)
            out << ',';
        out << '[';
        const auto &values = results.frequencies[f].excitations[0].observerPressure;
        for (std::size_t j = 0; j < values.size(); ++j) {
            if (j)
                out << ',';
            out << '[' << values[j].real() << ',' << values[j].imag() << ']';
        }
        out << ']';
    }
    out << "]\n}\n";
    if (!out)
        throw std::runtime_error("Could not finish the external parity fixture.");
}
void checkExternalReference(const BemResult &results) {
    std::ifstream input(HS_PARITY_REFERENCE_PATH);
    if (!input)
        throw std::runtime_error("Pinned external reference data is missing.");
    std::string line;
    std::size_t frequencyIndex = 0;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#')
            continue;
        std::istringstream row(line);
        row.imbue(std::locale::classic());
        double frequency{};
        row >> frequency;
        require(frequencyIndex < results.frequencies.size(),
                "External reference has unexpected frequencies");
        const auto &actual = results.frequencies[frequencyIndex++];
        require(actual.frequency == frequency, "External reference frequency mismatch");
        double errorSquared = 0, referenceSquared = 0;
        for (std::size_t j = 0; j < 3; ++j) {
            char comma{};
            double real{}, imag{};
            row >> comma;
            require(comma == ',', "Invalid reference CSV");
            row >> real;
            row >> comma;
            require(comma == ',', "Invalid reference CSV");
            row >> imag;
            require(static_cast<bool>(row), "Invalid reference pressure");
            const std::complex<double> expected{real, imag};
            const auto pressure = actual.excitations[0].observerPressure[j],
                       ratio = pressure / expected;
            errorSquared += std::norm(pressure - expected);
            referenceSquared += std::norm(expected);
            require(std::abs(20 * std::log10(std::abs(ratio))) <= .1,
                    "External reference level gate failed");
            require(std::abs(std::arg(ratio)) * 180 / std::numbers::pi <= 1,
                    "External reference phase gate failed");
        }
        const double relative = std::sqrt(errorSquared / referenceSquared);
        require(relative <= .01, "External reference complex L2 gate failed");
        std::cout << "Pinned BEAT CPU reference, " << frequency << " Hz: complex L2 " << relative
                  << '\n';
    }
    require(frequencyIndex == results.frequencies.size(),
            "External reference frequencies are incomplete");
}
void checkDipoleFarField() {
    auto mesh = sphere(4);
    BemRequest input;
    input.frequencies = {1715};
    input.observers = {{0, 0, 1}, {0, 0, -1}, {1, 0, 0}};
    input.farField = FarFieldObservation{};
    input.excitations.push_back({"translation", {}});
    for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
        auto &t = mesh.triangles[i];
        t.faceId = static_cast<std::uint32_t>(i + 1);
        const auto a = mesh.vertices[t.vertices[0]], b = mesh.vertices[t.vertices[1]],
                   c = mesh.vertices[t.vertices[2]];
        const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
        const Vec3 n{u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
        input.excitations.front().faces.push_back({t.faceId, {n.z / std::hypot(n.x, n.y, n.z), 0}});
    }
    CudaBemSolver solver(std::move(mesh));
    const auto result = solver.solve(input);
    const auto &p = result.frequencies.front().excitations.front().observerPressure;
    const double k = 2 * std::numbers::pi * 1715 / input.medium.soundSpeed, z = k * .1;
    const std::complex<double> i{0, 1};
    const auto derivative = std::exp(-i * z) * (i / z + 2. / (z * z) - 2. * i / (z * z * z));
    const auto expected =
        i * input.medium.density * input.medium.soundSpeed / k / derivative * std::exp(-i * k);
    const double error = std::abs(p[0] / expected - 1.);
    require(error < .02 && std::abs(p[1] / p[0] + 1.) < 1e-4 && std::abs(p[2] / p[0]) < 1e-4,
            "Asymptotic dipole pattern fails its analytic phase, reversal or transverse null.");
    input.farField.reset();
    input.observers = {{0, 0, 1}, {0, 0, 10}, {0, 0, 1000}};
    const auto finite = solver.solve(input);
    double previous = 1;
    for (std::size_t j = 0; j < input.observers.size(); ++j) {
        const double r = input.observers[j].z;
        const auto normalized = finite.frequencies.front().excitations.front().observerPressure[j] *
                                r * std::exp(i * k * (r - 1));
        const double difference = std::abs(normalized / p[0] - 1.);
        require(difference < previous * .2,
                "Dipole finite-radius field does not converge to its asymptote.");
        previous = difference;
    }
    require(previous < 5e-5, "Dipole far-field limit differs from the finite-radius GPU solution.");
    std::cout << "Asymptotic dipole analytic error " << error << ", finite-radius limit error "
              << previous << '\n';
}
void checkAsymptoticFarField() {
    auto input = request();
    input.frequencies = {100, 1000, 1715};
    input.observers = {{0, 0, .01}, {.01, 0, 0}, {.006, .008, 0}};
    input.farField = FarFieldObservation{};
    CudaBemSolver coarse(sphere(3)), fine(sphere(4));
    const auto coarseResult = coarse.solve(input), farResult = fine.solve(input);
    double coarseError = 0, fineError = 0;
    for (std::size_t f = 0; f < input.frequencies.size(); ++f)
        for (std::size_t j = 0; j < input.observers.size(); ++j) {
            const auto expected = analytic(input.frequencies[f], 1);
            const auto &unit = farResult.frequencies[f].excitations[0];
            coarseError = std::max(
                coarseError,
                std::abs(coarseResult.frequencies[f].excitations[0].observerPressure[j] / expected -
                         1.));
            fineError = std::max(fineError, std::abs(unit.observerPressure[j] / expected - 1.));
            require(std::abs(farResult.frequencies[f].excitations[1].observerPressure[j] /
                                 unit.observerPressure[j] -
                             std::polar(.6, .71)) < 3e-4,
                    "Asymptotic field lost complex multi-RHS linearity.");
        }
    require(fineError < .015 && fineError < .7 * coarseError,
            "Asymptotic sphere pressure did not converge to the analytic reference.");
    auto distant = input;
    distant.frequencies = {1715};
    distant.observers = {{0, 0, 1}, {0, 0, 10}, {0, 0, 100}, {0, 0, 10000}};
    distant.farField.reset();
    const auto finite = fine.solve(distant);
    const double k = 2 * std::numbers::pi * 1715 / input.medium.soundSpeed;
    const auto reference = farResult.frequencies.back().excitations[0].observerPressure[0];
    double firstError = 0, lastError = 0;
    for (std::size_t j = 0; j < distant.observers.size(); ++j) {
        const double r = distant.observers[j].z;
        const auto normalized = finite.frequencies[0].excitations[0].observerPressure[j] * r *
                                std::exp(std::complex<double>(0, k * (r - 1)));
        const double error = std::abs(normalized / reference - 1.);
        if (j == 0)
            firstError = error;
        lastError = error;
        std::cout << "finite-to-asymptotic radius " << r << " m, relative error " << error << '\n';
    }
    require(lastError < 5e-5 && lastError < .1 * firstError + 1e-5,
            "Finite GPU observations do not approach the asymptotic GPU field.");
    const auto &farSurface = farResult.frequencies.back().excitations[0];
    const auto &finiteSurface = finite.frequencies[0].excitations[0];
    require(std::abs(farSurface.impedances[0].volumeImpedance /
                         finiteSurface.impedances[0].volumeImpedance -
                     1.) < 1e-4,
            "Far-field mode changed the solved radiation impedance.");
    for (std::size_t j = 0; j < farSurface.boundary()->surfacePressure.size(); ++j)
        require(std::abs(farSurface.boundary()->surfacePressure[j] /
                             finiteSurface.boundary()->surfacePressure[j] -
                         1.) < 1e-4,
                "Far-field mode changed surface pressure.");
    auto scaled = input;
    scaled.frequencies = {1715};
    scaled.farField->referenceDistance = 2;
    for (auto &observer : scaled.observers) {
        observer.x *= 1e6;
        observer.y *= 1e6;
        observer.z *= 1e6;
    }
    const auto scaledResult = fine.solve(scaled);
    for (std::size_t j = 0; j < scaled.observers.size(); ++j)
        require(std::abs(scaledResult.frequencies[0].excitations[0].observerPressure[j] /
                             farResult.frequencies.back().excitations[0].observerPressure[j] * 2. *
                             std::exp(std::complex<double>(0, k)) -
                         1.) < 1e-4,
                "Asymptotic reference distance or direction normalization is inconsistent.");
    scaled.farField->referenceDistance = 1e15;
    const auto enormous = fine.solve(scaled);
    for (std::size_t j = 0; j < scaled.observers.size(); ++j)
        require(
            std::abs(enormous.frequencies[0].excitations[0].observerPressure[j] /
                         farResult.frequencies.back().excitations[0].observerPressure[j] *
                         scaled.farField->referenceDistance *
                         std::exp(std::complex<double>(0, k * scaled.farField->referenceDistance)) *
                         std::exp(std::complex<double>(0, -k)) -
                     1.) < 1e-4,
            "An enormous reference distance erased the aperture phase.");
    auto shiftedOrigin = input;
    shiftedOrigin.frequencies = {1715};
    shiftedOrigin.farField->origin = {.02, -.01, .03};
    shiftedOrigin.observers = {{.02, -.01, 1.03}};
    const auto shifted = fine.solve(shiftedOrigin);
    require(std::abs(shifted.frequencies[0].excitations[0].observerPressure[0] / reference *
                         std::exp(std::complex<double>(0, k * .03)) -
                     1.) < 1e-4,
            "Changing the far-field phase origin retained a stale observer cache.");
    const auto rejectsRequest = [&](const BemRequest &bad) {
        bool rejected = false;
        try {
            (void)fine.solve(bad);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Invalid far-field observation was accepted.");
    };
    auto invalid = input;
    invalid.observers = {{0, 0, 0}};
    rejectsRequest(invalid);
    invalid = input;
    invalid.farField->referenceDistance = 0;
    rejectsRequest(invalid);
    invalid = input;
    invalid.farField->origin.x = std::numeric_limits<double>::infinity();
    rejectsRequest(invalid);
    // These same coordinates are legitimate direction markers but are inside
    // the body in finite mode. Mode switches must invalidate observer caching.
    invalid = input;
    invalid.farField.reset();
    rejectsRequest(invalid);
    auto translated = sphere(3);
    const Vec3 offset{1e9, -1e9, 1e9};
    for (auto &p : translated.vertices) {
        p.x += offset.x;
        p.y += offset.y;
        p.z += offset.z;
    }
    auto rounded = translated;
    for (auto &p : rounded.vertices) {
        p.x -= offset.x;
        p.y -= offset.y;
        p.z -= offset.z;
    }
    auto world = input;
    world.frequencies = {1000};
    world.farField->origin = offset;
    for (auto &p : world.observers) {
        p.x += offset.x;
        p.y += offset.y;
        p.z += offset.z;
    }
    auto local = world;
    local.farField->origin = {};
    for (auto &p : local.observers) {
        p.x -= offset.x;
        p.y -= offset.y;
        p.z -= offset.z;
    }
    CudaBemSolver translatedSolver(std::move(translated)), roundedSolver(std::move(rounded));
    const auto worldResult = translatedSolver.solve(world),
               localResult = roundedSolver.solve(local);
    for (std::size_t e = 0; e < world.excitations.size(); ++e)
        for (std::size_t j = 0; j < world.observers.size(); ++j)
            require(
                std::abs(worldResult.frequencies[0].excitations[e].observerPressure[j] /
                             localResult.frequencies[0].excitations[e].observerPressure[j] -
                         1.) < 1e-4,
                "Common translation changed the asymptotic field on identically rounded geometry.");
    std::cout << "Asymptotic GPU field: coarse/fine analytic errors " << coarseError << " / "
              << fineError << ", distance/phase/cache/translation checks passed.\n";
    checkDipoleFarField();
}
} // namespace

int main(int argc, char **argv) {
    try {
        std::cout << "CUDA device: " << CudaBemSolver::deviceDescription() << '\n';
        if (argc == 2 && std::string_view(argv[1]) == "--compact-quadrature") {
            checkCompactQuadrature();
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--mapped-quadrature") {
            checkMappedQuadrature();
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--memory-preflight") {
            checkMemoryPreflight();
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--tight-residual") {
            checkTightResidual();
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--thin-body") {
            checkThinBody();
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--coarse-draft") {
            checkCoarseDraft();
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--asymptotic") {
            checkAsymptoticFarField();
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--low-frequency") {
            auto lowInput = request();
            lowInput.frequencies = {20, 50};
            CudaBemSolver lowSolver(sphere(4));
            checkSphere(lowSolver.solve(lowInput), .015);
            std::cout << "20/50 Hz analytic CUDA pressure and complex impedance tests passed.\n";
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--far-phase") {
            auto farInput = request();
            farInput.frequencies = {1715};
            farInput.observers = {{0, 0, 1000}, {0, 0, 1e6}};
            CudaBemSolver farSolver(sphere(4));
            const auto farResult = farSolver.solve(farInput);
            const double k =
                2 * std::numbers::pi * farInput.frequencies.front() / farInput.medium.soundSpeed;
            double error = 0;
            for (const auto &excitation : farResult.frequencies.front().excitations) {
                const auto ratio = excitation.observerPressure[1] / excitation.observerPressure[0];
                error = std::max(
                    error,
                    std::abs(ratio * 1000.0 * std::exp(std::complex<double>(0, k * (1e6 - 1000))) -
                             1.0));
            }
            std::cout << "Far-observer phase/amplitude ratio relative error " << error << '\n';
            require(error < 1e-3, "Large-distance observer phase fails radial sphere propagation");
            return 0;
        }
        auto inputs = request();
        auto coarseInputs = inputs;
        coarseInputs.frequencies = {100, 500};
        CudaBemSolver coarse(sphere(2));
        const auto coarseResult = coarse.solve(coarseInputs);
        const auto coarseError = checkSphere(coarseResult, .16);
        auto strictNearInput = coarseInputs;
        strictNearInput.numerics.nearRelativeTolerance = 1e-6;
        const auto strictNearResult = coarse.solve(strictNearInput);
        double nearPrecisionError = 0;
        for (std::size_t f = 0; f < coarseInputs.frequencies.size(); ++f)
            for (std::size_t e = 0; e < coarseInputs.excitations.size(); ++e) {
                const auto &mixed = coarseResult.frequencies[f].excitations[e],
                           &strict = strictNearResult.frequencies[f].excitations[e];
                for (std::size_t j = 0; j < coarseInputs.observers.size(); ++j)
                    nearPrecisionError = std::max(
                        nearPrecisionError,
                        std::abs(mixed.observerPressure[j] / strict.observerPressure[j] - 1.0));
                nearPrecisionError = std::max(
                    nearPrecisionError, std::abs(mixed.impedances.front().volumeImpedance /
                                                     strict.impedances.front().volumeImpedance -
                                                 1.0));
            }
        require(
            nearPrecisionError < 1e-4,
            "Mixed and stricter double near quadrature disagree in complex pressure or impedance");
        std::cout << "Mixed versus strict double near pressure/impedance relative error "
                  << nearPrecisionError << '\n';
        auto alternateQuadrature = coarseInputs;
        alternateQuadrature.numerics.singularOrder = 6;
        const auto higherOrder = coarse.solve(alternateQuadrature);
        for (std::size_t f = 0; f < coarseInputs.frequencies.size(); ++f)
            for (std::size_t j = 0; j < coarseInputs.observers.size(); ++j)
                require(
                    std::abs(higherOrder.frequencies[f].excitations[0].observerPressure[j] /
                                 coarseResult.frequencies[f].excitations[0].observerPressure[j] -
                             1.0) < .003,
                    "Singular Gauss orders four and six did not converge");
        alternateQuadrature.numerics.quadratureOrder = 2;
        checkSphere(coarse.solve(alternateQuadrature), .16);
        CudaBemSolver fine(sphere(3));
        const auto fineResult = fine.solve(inputs);
        const auto fineError = checkSphere(fineResult, .035);
        checkExternalReference(fineResult);
        auto nearInput = inputs;
        nearInput.frequencies = {inputs.frequencies.back()};
        nearInput.observers = {{0, 0, .101}};
        const auto nearResult = fine.solve(nearInput);
        const auto nearPressure = nearResult.frequencies.front().excitations[0].observerPressure[0];
        require(std::abs(nearPressure / analytic(nearInput.frequencies.front(), .101) - 1.0) < .035,
                "Adaptive near field fails the analytic sphere reference");
        require(std::abs(nearResult.frequencies.front().excitations[1].observerPressure[0] /
                             nearPressure -
                         std::polar(.6, .71)) < 3e-4,
                "Adaptive near field loses complex excitation linearity");
        auto warmInput = inputs;
        warmInput.frequencies = {inputs.frequencies.back()};
        const auto gain = std::polar(.8, -.24);
        for (auto &excitation : warmInput.excitations)
            for (auto &face : excitation.faces)
                face.normalVelocity *= gain;
        const auto warmResult = fine.solve(warmInput);
        for (std::size_t e = 0; e < warmInput.excitations.size(); ++e)
            for (std::size_t j = 0; j < warmInput.observers.size(); ++j)
                require(
                    std::abs(warmResult.frequencies.front().excitations[e].observerPressure[j] /
                                 fineResult.frequencies.back().excitations[e].observerPressure[j] -
                             gain) < 3e-4,
                    "Persistent workspace retained a previous excitation RHS");
        auto resizedInput = warmInput;
        resizedInput.excitations.resize(1);
        resizedInput.observers.clear();
        const auto resizedResult = fine.solve(resizedInput);
        require(resizedResult.frequencies.front().excitations.front().observerPressure.empty(),
                "Zero-observer workspace retained a previous field buffer");
        require(std::abs(resizedResult.frequencies.front()
                                 .excitations.front()
                                 .impedances.front()
                                 .volumeImpedance /
                             warmResult.frequencies.front()
                                 .excitations.front()
                                 .impedances.front()
                                 .volumeImpedance -
                         1.0) < 3e-4,
                "Workspace resizing changed complex impedance");
        const auto restoredResult = fine.solve(warmInput);
        for (std::size_t e = 0; e < warmInput.excitations.size(); ++e)
            for (std::size_t j = 0; j < warmInput.observers.size(); ++j)
                require(
                    std::abs(restoredResult.frequencies.front().excitations[e].observerPressure[j] /
                                 warmResult.frequencies.front().excitations[e].observerPressure[j] -
                             1.0) < 3e-4,
                    "Restored workspace dimensions changed the field");
        double fineLowFrequencyError = 0;
        for (const auto &f : fineResult.frequencies)
            if (f.frequency <= 500)
                fineLowFrequencyError =
                    std::max(fineLowFrequencyError, std::abs(f.excitations[0].observerPressure[0] /
                                                                 analytic(f.frequency, .3) -
                                                             1.0));
        require(fineLowFrequencyError < coarseError * .85,
                "Sphere mesh refinement did not reduce the acoustic error");
        CudaBemSolver refined(sphere(4));
        const auto refinedResult = refined.solve(inputs);
        const auto refinedError = checkSphere(refinedResult, .015);
        require(refinedError < fineError * .7,
                "Third sphere refinement did not improve the resonant-frequency reference");
        if (argc == 2)
            exportReference(argv[1], sphere(3), inputs, fineResult);
        bool rejected = false;
        auto invalid = inputs;
        invalid.frequencies = {0};
        try {
            fine.solve(invalid);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Nonpositive frequency was accepted");
        const auto boundaryMesh = sphere(3);
        const auto face = boundaryMesh.triangles.front().vertices;
        const auto a = boundaryMesh.vertices[face[0]], b = boundaryMesh.vertices[face[1]],
                   c = boundaryMesh.vertices[face[2]];
        for (const auto point : std::array<Vec3, 4>{
                 Vec3{0, 0, 0}, a, Vec3{(a.x + b.x) / 2, (a.y + b.y) / 2, (a.z + b.z) / 2},
                 Vec3{.2 * a.x + .3 * b.x + .5 * c.x, .2 * a.y + .3 * b.y + .5 * c.y,
                      .2 * a.z + .3 * b.z + .5 * c.z}}) {
            invalid = inputs;
            invalid.observers = {point};
            rejected = false;
            try {
                fine.solve(invalid);
            } catch (const std::invalid_argument &) {
                rejected = true;
            }
            require(rejected, "Interior or boundary observer was accepted");
        }
        std::stop_source stop;
        stop.request_stop();
        rejected = false;
        try {
            fine.solve(inputs, stop.get_token());
        } catch (const std::runtime_error &) {
            rejected = true;
        }
        require(rejected, "Cancelled GPU request was accepted");
        for (const double offset : {1e9, 1e12}) {
            auto translated = sphere(2);
            for (auto &vertex : translated.vertices) {
                vertex.x += offset;
                vertex.y -= offset;
                vertex.z += offset;
            }
            auto roundedLocal = translated;
            for (auto &vertex : roundedLocal.vertices) {
                vertex.x -= offset;
                vertex.y += offset;
                vertex.z -= offset;
            }
            auto worldInput = inputs;
            worldInput.frequencies = {500};
            for (auto &observer : worldInput.observers) {
                observer.x += offset;
                observer.y -= offset;
                observer.z += offset;
            }
            auto localInput = worldInput;
            for (auto &observer : localInput.observers) {
                observer.x -= offset;
                observer.y += offset;
                observer.z -= offset;
            }
            CudaBemSolver translatedSolver(std::move(translated)),
                localSolver(std::move(roundedLocal));
            const auto worldResult = translatedSolver.solve(worldInput);
            const auto localResult = localSolver.solve(localInput);
            for (std::size_t e = 0; e < worldInput.excitations.size(); ++e) {
                const auto &world = worldResult.frequencies.front().excitations[e],
                           &local = localResult.frequencies.front().excitations[e];
                for (std::size_t j = 0; j < worldInput.observers.size(); ++j)
                    require(std::abs(world.observerPressure[j] / local.observerPressure[j] - 1.0) <
                                1e-4,
                            "Common CAD translation changed pressure on the same rounded geometry");
                require(std::abs(world.impedances.front().volumeImpedance /
                                     local.impedances.front().volumeImpedance -
                                 1.0) < 1e-4,
                        "Common CAD translation changed complex impedance");
            }
        }
        checkAsymptoticFarField();
        std::cout << "Analytic CUDA sphere, phase, impedance, refinement and cancellation tests "
                     "passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
