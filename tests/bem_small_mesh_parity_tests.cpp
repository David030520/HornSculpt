#include "bem_support_operator_oracle.hpp"
#include "fixture_geometry.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cuda_runtime_api.h>
#include <fmm/support_level_plan.hpp>
#include <hornsculpt/bem.hpp>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
using namespace hornsculpt;
using Complex = std::complex<double>;
using Clock = std::chrono::steady_clock;
constexpr Complex mixture{.27, -.31};
// Same-discretization comparisons allow FP32 dense assembly/factorization and
// the unchanged 5e-5 original-equation residual. These are not mesh-error gates.
constexpr double pressureParityTolerance = 1e-3;
constexpr double superpositionTolerance = 5e-4;
constexpr double levelToleranceDb = .02;
constexpr double phaseToleranceRadians = .1 * std::numbers::pi / 180;
constexpr double observableAmplitudeFraction = .05;
constexpr double residualTolerance = 5e-5;

void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}
bool finite(Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}
Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
double length(Vec3 v) {
    return std::hypot(v.x, v.y, v.z);
}
double area(const SurfaceMesh &mesh, const Triangle &triangle) {
    const auto a = mesh.vertices[triangle.vertices[0]],
               b = subtract(mesh.vertices[triangle.vertices[1]], a),
               c = subtract(mesh.vertices[triangle.vertices[2]], a);
    return .5 * length({b.y * c.z - b.z * c.y, b.z * c.x - b.x * c.z, b.x * c.y - b.y * c.x});
}
std::vector<double> panelAreas(const SurfaceMesh &mesh) {
    std::vector<double> result;
    for (const auto &triangle : mesh.triangles)
        result.push_back(area(mesh, triangle));
    return result;
}
double relative(std::span<const Complex> actual, std::span<const Complex> expected,
                std::span<const double> weights = {}) {
    require(actual.size() == expected.size() && !actual.empty() &&
                (weights.empty() || weights.size() == actual.size()),
            "Parity comparison requires matching nonempty vectors and weights.");
    double difference{}, norm{};
    for (std::size_t i = 0; i < actual.size(); ++i) {
        require(finite(actual[i]) && finite(expected[i]), "Nonfinite parity pressure.");
        const double weight = weights.empty() ? 1 : weights[i];
        difference += weight * std::norm(actual[i] - expected[i]);
        norm += weight * std::norm(expected[i]);
    }
    require(norm > 0, "Parity reference has zero norm.");
    return std::sqrt(difference / norm);
}

struct Fixture {
    SurfaceMesh mesh;
    Vec3 origin{};
    std::array<std::uint32_t, 2> ports{17, 29};
};

Fixture gradedEllipsoid() {
    Fixture result{testing::sphere(3), {.037, -.024, .019}};
    constexpr double grading = .45, rotation = .31;
    for (auto &point : result.mesh.vertices) {
        const auto u = Vec3{10 * point.x, 10 * point.y, 10 * point.z};
        // A smooth sphere automorphism redistributes existing conforming nodes;
        // the affine ellipsoid/shear then changes shape without opening seams.
        const double denominator = 1 + grading * u.z,
                     factor = std::sqrt(1 - grading * grading) / denominator;
        const Vec3 v{factor * u.x, factor * u.y, (u.z + grading) / denominator};
        const double x = .11 * v.x + .02 * v.z, y = .075 * v.y;
        point = add(result.origin, {std::cos(rotation) * x - std::sin(rotation) * y,
                                    std::sin(rotation) * x + std::cos(rotation) * y, .16 * v.z});
    }
    for (auto &triangle : result.mesh.triangles) {
        Vec3 center{};
        for (const auto vertex : triangle.vertices)
            center = add(center, subtract(result.mesh.vertices[vertex], result.origin));
        if (center.z > .10) {
            triangle.faceId = result.ports[0];
            triangle.boundary = BoundaryCondition::driven;
        } else if (center.x < -.10) {
            triangle.faceId = result.ports[1];
            triangle.boundary = BoundaryCondition::driven;
        } else {
            triangle.faceId = 90;
            triangle.boundary = BoundaryCondition::wall;
        }
    }
    return result;
}

Fixture hornShell() {
    Fixture result;
    result.origin = {.027, -.018, .031};
    constexpr unsigned segments = 16, intervals = 7;
    constexpr double throatRadius = .015, mouthRadius = .12, wallThickness = .005, hornLength = .22,
                     verticalRatio = .78;
    const auto vertex = [](unsigned side, unsigned axial, unsigned angle) {
        return static_cast<std::uint32_t>((side * (intervals + 1) + axial) * segments +
                                          angle % segments);
    };
    for (unsigned side = 0; side < 2; ++side)
        for (unsigned axial = 0; axial <= intervals; ++axial) {
            const double t = static_cast<double>(axial) / intervals,
                         radius = throatRadius + (mouthRadius - throatRadius) * t * t +
                                  side * wallThickness,
                         z = side ? -wallThickness + (hornLength + wallThickness) * t
                                  : hornLength * t;
            for (unsigned angle = 0; angle < segments; ++angle) {
                const double phi = 2 * std::numbers::pi * angle / segments;
                result.mesh.vertices.push_back(
                    add(result.origin,
                        {radius * std::cos(phi), verticalRatio * radius * std::sin(phi), z}));
            }
        }
    const auto wall = [&](std::array<std::uint32_t, 3> indices) {
        result.mesh.triangles.push_back({indices, 90, BoundaryCondition::wall});
    };
    for (unsigned side = 0; side < 2; ++side)
        for (unsigned axial = 0; axial < intervals; ++axial)
            for (unsigned angle = 0; angle < segments; ++angle) {
                const auto a = vertex(side, axial, angle), b = vertex(side, axial, angle + 1),
                           c = vertex(side, axial + 1, angle),
                           d = vertex(side, axial + 1, angle + 1);
                if (side) {
                    wall({a, b, d});
                    wall({a, d, c});
                } else {
                    wall({a, d, b});
                    wall({a, c, d});
                }
            }
    for (unsigned angle = 0; angle < segments; ++angle) {
        const auto a = vertex(0, intervals, angle), b = vertex(0, intervals, angle + 1),
                   c = vertex(1, intervals, angle), d = vertex(1, intervals, angle + 1);
        wall({a, c, d});
        wall({a, d, b});
    }
    // The throat diaphragm faces +Z into the air cavity; the separate rear
    // surface faces -Z. The mouth remains open, with a physical annular rim.
    const auto throatCenter = static_cast<std::uint32_t>(result.mesh.vertices.size());
    result.mesh.vertices.push_back(result.origin);
    const auto rearCenter = static_cast<std::uint32_t>(result.mesh.vertices.size());
    result.mesh.vertices.push_back(add(result.origin, {0, 0, -wallThickness}));
    for (unsigned angle = 0; angle < segments; ++angle) {
        result.mesh.triangles.push_back(
            {{throatCenter, vertex(0, 0, angle), vertex(0, 0, angle + 1)},
             result.ports[angle < segments / 2 ? 0 : 1],
             BoundaryCondition::driven});
        wall({rearCenter, vertex(1, 0, angle + 1), vertex(1, 0, angle)});
    }
    require(result.mesh.triangles.size() == 512, "Horn fixture changed its declared panel count.");
    return result;
}

Fixture separatedComponents() {
    Fixture result;
    result.origin = {.6, 0, 0};
    result.ports = {1, 5};
    constexpr std::array<Vec3, 4> tetra{
        {{0, 0, 0}, {.00525, 0, 0}, {0, .00425, 0}, {0, 0, .00325}}};
    // Electrically large group supports and electrically small closed leaves
    // exercise the HF/LF handoff without acoustically coarse individual panels.
    for (unsigned side = 0; side < 2; ++side)
        for (unsigned x = 0; x < 2; ++x)
            for (unsigned y = 0; y < 2; ++y) {
                const Vec3 center{side * 1.2 + (x ? .05 : -.05), y ? .05 : -.05, 0};
                const auto first = static_cast<std::uint32_t>(result.mesh.vertices.size());
                for (const auto point : tetra)
                    result.mesh.vertices.push_back(add(center, point));
                for (const std::array<std::uint32_t, 3> indices :
                     {std::array<std::uint32_t, 3>{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}})
                    result.mesh.triangles.push_back(
                        {{first + indices[0], first + indices[1], first + indices[2]},
                         result.ports[side],
                         BoundaryCondition::driven});
            }
    return result;
}

BemRequest request(const Fixture &fixture, std::vector<double> frequencies, bool farField) {
    BemRequest result;
    result.medium = {1.17, 346.1};
    result.frequencies = std::move(frequencies);
    const std::array<Complex, 2> first{{{.8, -.3}, {-.2, .1}}}, second{{{-.1, .5}, {.6, -.4}}};
    result.excitations = {{"independent complex A", {}},
                          {"independent complex B", {}},
                          {"complex superposition", {}}};
    for (unsigned port = 0; port < fixture.ports.size(); ++port) {
        result.excitations[0].faces.push_back({fixture.ports[port], first[port]});
        result.excitations[1].faces.push_back({fixture.ports[port], second[port]});
        result.excitations[2].faces.push_back(
            {fixture.ports[port], first[port] + mixture * second[port]});
    }
    // Both cuts share the same explicit origin and cover rear as well as front
    // radiation. Finite samples are outside the complete fixture boundary.
    constexpr double azimuth = .37;
    for (unsigned cut = 0; cut < 2; ++cut)
        for (int degrees = -180; degrees <= 180; degrees += 30) {
            const double angle = degrees * std::numbers::pi / 180,
                         phi = azimuth + cut * std::numbers::pi / 2, radius = farField ? 1. : 2.;
            result.observers.push_back(
                add(fixture.origin,
                    {radius * std::sin(angle) * std::cos(phi),
                     radius * std::sin(angle) * std::sin(phi), radius * std::cos(angle)}));
        }
    if (farField)
        result.farField = FarFieldObservation{fixture.origin, 1.};
    result.numerics.relativeResidualTolerance = residualTolerance;
    result.numerics.fastMultipole.leafSupportCapacity = 4;
    return result;
}

struct TimedResult {
    BemResult result;
    double milliseconds{};
};
TimedResult solve(const SurfaceMesh &mesh, const BemRequest &input) {
    const auto start = Clock::now();
    CudaBemSolver solver(mesh);
    auto result = solver.solve(input);
    return {std::move(result),
            std::chrono::duration<double, std::milli>(Clock::now() - start).count()};
}

void checkExecution(const BemResult &result, const BemRequest &input) {
    require(result.acceleration == input.numerics.acceleration &&
                result.frequencies.size() == input.frequencies.size(),
            "Public solve changed the requested execution path or frequency count.");
    for (const auto &frequency : result.frequencies) {
        if (result.acceleration == BemAcceleration::fastMultipole) {
            const auto &stats = frequency.fastMultipole;
            require(
                stats.nodeCount > 1 && stats.leafCount > 1 && stats.farInteractionCount > 0 &&
                    stats.pressurePressureFarPairs > 0 && stats.pressureFluxFarPairs > 0,
                "Small-mesh parity did not execute genuine pressure and flux far interactions.");
            require(stats.pressurePressureNearPairs + stats.pressurePressureFarPairs ==
                            result.pressureDofCount * result.pressureDofCount &&
                        stats.pressureFluxNearPairs + stats.pressureFluxFarPairs ==
                            result.pressureDofCount * result.triangleCount,
                    "Public FMM interaction partition lost or duplicated original DOF pairs.");
            std::cout << frequency.frequency << " Hz: M2L=" << stats.farInteractionCount
                      << ", P1/P1 near/far=" << stats.pressurePressureNearPairs << '/'
                      << stats.pressurePressureFarPairs
                      << ", P1/DP0 near/far=" << stats.pressureFluxNearPairs << '/'
                      << stats.pressureFluxFarPairs << '\n';
        }
        require(frequency.excitations.size() == input.excitations.size(),
                "Public solve changed the independent RHS count.");
        for (const auto &excitation : frequency.excitations) {
            require(std::isfinite(excitation.boundaryResidual) &&
                        excitation.boundaryResidual <= residualTolerance,
                    "Public solution failed the unchanged original-equation residual target.");
            if (result.acceleration == BemAcceleration::fastMultipole)
                require(excitation.iterations > 0 && excitation.fastMultipoleAccuracySeed > 0,
                        "Public FMM solve did not report actual Krylov work and accuracy seed.");
        }
    }
}

std::array<double, 2> fieldLevelAndPhase(std::span<const Complex> actual,
                                         std::span<const Complex> expected) {
    double referencePower{};
    for (const auto value : expected)
        referencePower += std::norm(value);
    const double amplitudeFloor =
        observableAmplitudeFraction * std::sqrt(referencePower / expected.size());
    std::array<double, 2> worst{};
    for (std::size_t observer = 0; observer < expected.size(); ++observer) {
        if (std::abs(expected[observer]) < amplitudeFloor)
            continue;
        const auto ratio = actual[observer] / expected[observer];
        require(finite(ratio) && std::abs(ratio) > 0, "Invalid observable H/V pressure ratio.");
        worst[0] = std::max(worst[0], std::abs(20 * std::log10(std::abs(ratio))));
        worst[1] = std::max(worst[1], std::abs(std::arg(ratio)));
    }
    return worst;
}

void compareResults(const BemResult &actual, const BemResult &expected,
                    std::span<const double> areas) {
    require(actual.triangleCount == expected.triangleCount &&
                actual.pressureDofCount == expected.pressureDofCount &&
                actual.triangleCount == areas.size(),
            "FMM changed the same-discretization parity mesh or pressure space.");
    for (std::size_t f = 0; f < expected.frequencies.size(); ++f) {
        const auto &af = actual.frequencies[f], &ef = expected.frequencies[f];
        require(af.frequency == ef.frequency, "FMM reordered parity frequencies.");
        for (std::size_t e = 0; e < ef.excitations.size(); ++e) {
            const auto &a = af.excitations[e], &b = ef.excitations[e];
            require(a.name == b.name && a.impedances.size() == b.impedances.size(),
                    "FMM reordered excitations or physical ports.");
            const double surface = relative(a.boundary()->surfacePressure, b.boundary()->surfacePressure, areas),
                         field = relative(a.observerPressure, b.observerPressure);
            const auto levelPhase = fieldLevelAndPhase(a.observerPressure, b.observerPressure);
            double portPower{}, portArea{};
            for (const auto &port : b.impedances) {
                portPower += port.area * std::norm(port.meanPressure);
                portArea += port.area;
            }
            require(portPower > 0 && portArea > 0, "Parity has no active physical port pressure.");
            const double portPressureFloor = 1e-3 * std::sqrt(portPower / portArea);
            double worstPort{};
            for (std::size_t p = 0; p < b.impedances.size(); ++p) {
                const auto &x = a.impedances[p], &y = b.impedances[p];
                require(x.faceId == y.faceId && x.area == y.area && x.area > 0 &&
                            std::abs(x.volumeVelocity - y.volumeVelocity) <=
                                1e-12 * std::max(1e-18, std::abs(y.volumeVelocity)) &&
                            finite(x.meanPressure) && finite(x.volumeImpedance),
                        "FMM changed physical face identity, area, prescribed flow or units.");
                const double denominator = std::max(portPressureFloor, std::abs(y.meanPressure));
                worstPort =
                    std::max(worstPort, std::abs(x.meanPressure - y.meanPressure) / denominator);
                require(std::abs(x.volumeVelocity) > 0 &&
                            std::abs(x.volumeImpedance * x.volumeVelocity - x.meanPressure) <=
                                1e-12 * std::max(1., std::abs(x.meanPressure)),
                        "Physical port impedance does not satisfy Z=p/Q.");
                // Since both paths prescribe exactly the same complex Q, this
                // is also a complex impedance comparison with a shared floor.
                worstPort = std::max(worstPort, std::abs(x.volumeImpedance - y.volumeImpedance) *
                                                    std::abs(y.volumeVelocity) / denominator);
            }
            std::cout << af.frequency << " Hz, " << a.name
                      << ": surface(area-L2)/field/port=" << surface << '/' << field << '/'
                      << worstPort << ", observable level dB/phase degrees=" << levelPhase[0] << '/'
                      << levelPhase[1] * 180 / std::numbers::pi
                      << ", dense/FMM residual=" << b.boundaryResidual << '/' << a.boundaryResidual
                      << ", audit="
                      << (a.residualEvaluation == BemResidualEvaluation::directGalerkin
                              ? "direct original Galerkin"
                              : "higher-order FMM estimate")
                      << '\n';
            require(surface <= pressureParityTolerance && field <= pressureParityTolerance &&
                        worstPort <= pressureParityTolerance,
                    "Complex surface/field/port same-discretization parity failed.");
            require(levelPhase[0] <= levelToleranceDb && levelPhase[1] <= phaseToleranceRadians,
                    "H/V field exceeds the declared level/phase gates above the amplitude mask.");
        }
    }
}

void checkSuperposition(const BemResult &result, std::span<const double> areas) {
    for (const auto &frequency : result.frequencies) {
        const auto &a = frequency.excitations[0], &b = frequency.excitations[1],
                   &mixed = frequency.excitations[2];
        for (unsigned field = 0; field < 2; ++field) {
            const auto &first = field ? a.observerPressure : a.boundary()->surfacePressure,
                       &second = field ? b.observerPressure : b.boundary()->surfacePressure,
                       &actual = field ? mixed.observerPressure : mixed.boundary()->surfacePressure;
            std::vector<Complex> expected;
            for (std::size_t i = 0; i < first.size(); ++i)
                expected.push_back(first[i] + mixture * second[i]);
            const double error =
                relative(actual, expected, field ? std::span<const double>{} : areas);
            std::cout << frequency.frequency << " Hz, "
                      << (result.acceleration == BemAcceleration::dense ? "Dense" : "FMM")
                      << (field ? " field" : " surface") << " complex superposition=" << error
                      << '\n';
            require(error <= superpositionTolerance,
                    "Independent complex port RHS violate surface/field superposition.");
        }
    }
}

void describeMesh(std::string_view name, const Fixture &fixture, const BemRequest &input) {
    double largestEdge{};
    for (const auto &triangle : fixture.mesh.triangles)
        for (unsigned edge = 0; edge < 3; ++edge)
            largestEdge = std::max(
                largestEdge,
                length(subtract(fixture.mesh.vertices[triangle.vertices[edge]],
                                fixture.mesh.vertices[triangle.vertices[(edge + 1) % 3]])));
    std::cout << name << ": " << fixture.mesh.triangles.size() << " triangles, "
              << fixture.mesh.vertices.size() << " P1 DOFs, h_max=" << largestEdge << " m"
              << ", rho=" << input.medium.density << " kg/m^3, c=" << input.medium.soundSpeed
              << " m/s, quadrature/singular=" << input.numerics.quadratureOrder << '/'
              << input.numerics.singularOrder
              << ", near tolerance/depth=" << input.numerics.nearRelativeTolerance << '/'
              << input.numerics.nearMaximumDepth << ", residual gate=" << residualTolerance << '\n';
    for (const auto frequency : input.frequencies)
        std::cout << "  " << frequency
                  << " Hz, h_max/lambda=" << largestEdge * frequency / input.medium.soundSpeed
                  << ", k*h_max="
                  << 2 * std::numbers::pi * largestEdge * frequency / input.medium.soundSpeed
                  << '\n';
    std::cout << "Same mesh/quadrature parity only; no physical mesh-convergence claim.\n";
}

struct PairResults {
    BemResult dense, fast;
};
PairResults checkFixture(std::string_view name, const Fixture &fixture, BemRequest input) {
    describeMesh(name, fixture, input);
    const auto weights = panelAreas(fixture.mesh);
    auto fastInput = input;
    fastInput.numerics.acceleration = BemAcceleration::fastMultipole;
    const auto dense = solve(fixture.mesh, input), fast = solve(fixture.mesh, fastInput);
    checkExecution(dense.result, input);
    checkExecution(fast.result, fastInput);
    compareResults(fast.result, dense.result, weights);
    checkSuperposition(dense.result, weights);
    checkSuperposition(fast.result, weights);
    std::cout << "Constructor+solve wall milliseconds dense/FMM=" << dense.milliseconds << '/'
              << fast.milliseconds
              << ", peak device bytes dense/FMM=" << dense.result.peakDeviceBytes << '/'
              << fast.result.peakDeviceBytes << '\n';
    for (std::size_t f = 0; f < dense.result.frequencies.size(); ++f) {
        const auto &d = dense.result.frequencies[f], &a = fast.result.frequencies[f];
        std::cout << d.frequency
                  << " Hz assembly/solve/audit/field milliseconds dense=" << d.assemblyMilliseconds
                  << '/' << d.solveMilliseconds << '/' << d.residualCheckMilliseconds << '/'
                  << d.fieldMilliseconds << ", FMM=" << a.assemblyMilliseconds << '/'
                  << a.solveMilliseconds << '/' << a.residualCheckMilliseconds << '/'
                  << a.fieldMilliseconds << '\n';
    }
    return {dense.result, fast.result};
}

void independentResidual(const Fixture &fixture, const BemRequest &input,
                         const PairResults &solutions) {
    const auto space = galerkin::buildSpace(fixture.mesh);
    galerkin::PairQuadratureSettings settings;
    settings.regularOrder = input.numerics.quadratureOrder;
    settings.singularOrder = input.numerics.singularOrder;
    settings.nearDistanceRatio = input.numerics.nearDistanceRatio;
    settings.relativeTolerance = input.numerics.nearRelativeTolerance;
    settings.maximumSubdivisionDepth = input.numerics.nearMaximumDepth;
    const auto n = space.vertices.size();
    for (std::size_t f = 0; f < input.frequencies.size(); ++f) {
        const double omega = 2 * std::numbers::pi * input.frequencies[f],
                     k = omega / input.medium.soundSpeed;
        const auto reference = test::support_operator_oracle::directReference(space, k, settings);
        for (const auto *solution : {&solutions.dense, &solutions.fast})
            for (std::size_t e = 0; e < input.excitations.size(); ++e) {
                const auto &surface = solution->frequencies[f].excitations[e].boundary()->surfacePressure;
                std::vector<Complex> pressure(n);
                // A tetrahedron's four face means determine its four P1 nodal
                // values exactly: p_i=sum(face_means)-3*opposite_face_mean.
                // This recovers the public solution without exposing internals.
                for (std::size_t component = 0; component < n / 4; ++component) {
                    Complex sum{};
                    for (unsigned face = 0; face < 4; ++face)
                        sum += surface[4 * component + face];
                    for (unsigned node = 0; node < 4; ++node)
                        pressure[4 * component + node] =
                            sum - 3. * surface[4 * component + 3 - node];
                }
                std::array<Complex, 2> flux{};
                for (const auto &face : input.excitations[e].faces)
                    flux[face.faceId == 5 ? 0 : 1] =
                        Complex{0, -input.medium.density * omega} * face.normalVelocity;
                double residualPower{}, rightPower{};
                for (std::size_t row = 0; row < n; ++row) {
                    Complex right{}, action{};
                    for (unsigned face = 0; face < 2; ++face)
                        right += reference.b[2 * row + face] * flux[face];
                    for (std::size_t column = 0; column < n; ++column)
                        action += reference.a[n * row + column] * pressure[column];
                    residualPower += std::norm(action - right);
                    rightPower += std::norm(right);
                }
                require(rightPower > 0, "Independent direct residual has zero RHS.");
                const double measured = std::sqrt(residualPower / rightPower);
                std::cout << input.frequencies[f] << " Hz, " << input.excitations[e].name << ", "
                          << (solution->acceleration == BemAcceleration::dense ? "Dense" : "FMM")
                          << " independent CPU original weak-form residual=" << measured << '\n';
                require(std::isfinite(measured) && measured <= residualTolerance,
                        "Independent direct CPU Galerkin residual exceeds the unchanged target.");
            }
    }
}

void checkWideband() {
    const auto fixture = separatedComponents();
    auto input = request(fixture, {300, 55 * 346.1 / (2 * std::numbers::pi)}, true);
    input.numerics.nearDistanceRatio = .5;
    input.numerics.fastMultipole.leafSupportCapacity = 4;
    const auto space = galerkin::buildSpace(fixture.mesh);
    wideband::SupportTreeSettings settings;
    settings.leafCapacity = input.numerics.fastMultipole.leafSupportCapacity;
    settings.nearDistanceRatio = input.numerics.nearDistanceRatio;
    settings.separationRatio = input.numerics.fastMultipole.separationRatio;
    const auto tree = wideband::makeSupportTree(space, settings);
    std::vector<std::uint32_t> faces;
    for (const auto &panel : space.panels)
        faces.push_back(panel.faceId == 1 ? 0u : 1u);
    // Inspect the exact public-request expansion plan before solving; this
    // protects the test's LF/HF scope even when the cost policy uses a direct
    // final residual audit and therefore exports no audit-level evidence.
    const auto plan = wideband::makeSupportExpansionPlan(
        tree, faces, 55, .1 * input.numerics.relativeResidualTolerance);
    bool hf = false, lf = false, handoff = false;
    for (std::size_t index = 0; index < plan.levels.size(); ++index) {
        const auto &level = plan.levels[index];
        if (level.expansion.representation == wideband::LevelRepresentation::inactive)
            continue;
        std::cout << "Wideband planned level=" << index
                  << ", k*a=" << level.expansion.electricalRadius
                  << ", radius=" << level.expansion.activeRadius
                  << ", separation=" << level.policySeparation << ", representation="
                  << (level.expansion.representation == wideband::LevelRepresentation::planeWave
                          ? "HF"
                          : "LF")
                  << ", moment order=" << level.expansion.signatureOrder << '\n';
        hf |= level.expansion.representation == wideband::LevelRepresentation::planeWave;
        lf |= level.expansion.representation == wideband::LevelRepresentation::regularSpherical;
    }
    for (std::size_t node = 0; node < tree.nodes.size(); ++node) {
        const auto &child = tree.nodes[node];
        if (child.parent == wideband::noSupportNode || !plan.activeChannels[node])
            continue;
        handoff |= plan.levels[child.level].expansion.representation ==
                       wideband::LevelRepresentation::regularSpherical &&
                   plan.levels[tree.nodes[child.parent].level].expansion.representation ==
                       wideband::LevelRepresentation::planeWave;
    }
    require(hf && lf && handoff,
            "Wideband public parity fixture must retain active HF, LF and their handoff.");
    const auto solutions = checkFixture("Separated closed components, LF/HF", fixture, input);
    independentResidual(fixture, input, solutions);
}
} // namespace

int main(int argc, char **argv) {
    int devices{};
    if (cudaGetDeviceCount(&devices) != cudaSuccess || devices == 0) {
        std::cout << "CUDA unavailable; no small-mesh parity calculation executed.\n";
        return 77;
    }
    try {
        std::cout << std::setprecision(9) << CudaBemSolver::deviceDescription() << '\n';
        const std::string_view mode = argc == 2 ? argv[1] : "--all";
        require(argc <= 2, "Expected one small-mesh parity mode.");
        require(mode == "--all" || mode == "--graded" || mode == "--horn" || mode == "--wideband",
                "Unknown small-mesh parity mode.");
        if (mode == "--all" || mode == "--graded") {
            const auto fixture = gradedEllipsoid();
            // Fixed before evaluating parity errors: 2500 Hz keeps k*h below
            // the public 3.25 quadrature qualification limit on this 512T mesh.
            checkFixture("Graded asymmetric ellipsoid, finite H/V", fixture,
                         request(fixture, {300, 2500}, false));
        }
        if (mode == "--all" || mode == "--horn") {
            const auto fixture = hornShell();
            checkFixture("Elliptical horn shell, asymptotic H/V", fixture,
                         request(fixture, {300, 2500}, true));
        }
        if (mode == "--all" || mode == "--wideband")
            checkWideband();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Small-mesh parity failed: " << error.what() << '\n';
        return 1;
    }
}
