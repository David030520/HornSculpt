#include "../modules/bem/lem_coupling.hpp"
#include "fixture_geometry.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <hornsculpt/bem.hpp>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>

using namespace hornsculpt;
namespace {
using Complex = std::complex<double>;
constexpr Complex imaginary{0, 1};
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
void close(Complex actual, Complex expected, double tolerance, const char *message) {
    if (std::abs(actual - expected) > tolerance * std::max(1e-12, std::abs(expected))) {
        std::cerr << message << ": actual " << actual << ", expected " << expected << '\n';
        throw std::runtime_error(message);
    }
}
double triangleArea(const SurfaceMesh &mesh, const Triangle &triangle) {
    const auto a = mesh.vertices[triangle.vertices[0]], b = mesh.vertices[triangle.vertices[1]],
               c = mesh.vertices[triangle.vertices[2]];
    const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
    return .5 * std::hypot(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
}
lem::Network flowNetwork(const std::vector<std::vector<std::uint32_t>> &groups,
                         const std::vector<Complex> &flows, int lastPolarity = 1) {
    lem::Network network;
    network.nodes.push_back({{1}, lem::Domain::acoustic, true, "Reference", {}});
    for (std::size_t i = 0; i < groups.size(); ++i) {
        const auto id = static_cast<std::uint32_t>(i + 2);
        network.nodes.push_back({{id}, lem::Domain::acoustic, false, "Port", {}});
        const int polarity = i + 1 == groups.size() ? lastPolarity : 1;
        network.elements.push_back({{2 * id}, "Flow", {}, lem::FlowSource{{{1}, {id}}, flows[i]}});
        network.elements.push_back(
            {{2 * id + 1}, "Boundary", {}, lem::BemPort{{{id}, {1}}, groups[i], polarity}});
    }
    return network;
}
void checkClosureContract() {
    bem_detail::CouplingBasis basis;
    basis.areas = {.03};
    basis.panelPorts = {0, 0};
    basis.panelVelocities = {1 / .03, 1 / .03};
    const std::array<double, 2> areas{.01, .02};
    // The port mean vanishes although each panel carries finite pressure.
    std::array<Complex, 2> pressure{{{3, 6}, {-1.5, -3}}};
    lem::Solution solution;
    solution.ports.push_back({{1}, {}, {}, 0});
    bem_detail::validateCouplingClosure(basis, solution, areas, pressure);
    solution.ports[0].pressure = {1e-5, 0};
    bem_detail::validateCouplingClosure(basis, solution, areas, pressure);
    solution.ports[0].pressure = {.1, 0};
    bool rejected = false;
    try {
        bem_detail::validateCouplingClosure(basis, solution, areas, pressure);
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    require(rejected, "Closure check accepted an inconsistent pressure feedback");
    pressure = {};
    solution.ports[0].pressure = {};
    bem_detail::validateCouplingClosure(basis, solution, areas, pressure);
}
void checkModalBasisContract() {
    SurfaceMesh mesh;
    mesh.vertices = {{0, 0, 0}, {.2, 0, 0}, {0, .1, 0}, {.4, 0, 0}};
    mesh.triangles = {{{0, 1, 2}, 1, BoundaryCondition::driven},
                      {{0, 2, 3}, 4, BoundaryCondition::driven}};
    const std::array<double, 2> areas{.01, .02};
    auto network = flowNetwork({{1, 4}}, {{.002, -.001}});
    auto &port = std::get<lem::BemPort>(network.elements[1].component);
    port.velocityWeights = {{4, .5}};
    auto basis = bem_detail::makeCouplingBasis(network, mesh, areas);
    close(basis.areas[0], .02, 1e-14, "Weighted port uses the wrong modal area.");
    close(basis.excitations[0].faces[0].normalVelocity, 50, 1e-14, "Dome unit-Q shape is wrong.");
    close(basis.excitations[0].faces[1].normalVelocity, 25, 1e-14,
          "Surround unit-Q shape is wrong.");
    const Complex q{.002, -.001};
    const std::array<Complex, 2> pressure{Complex{2, 3}, Complex{6, -1}};
    const Complex modal = .5 * (pressure[0] + pressure[1]);
    lem::Solution solution;
    solution.ports.push_back({{1}, modal, q, std::real(modal * std::conj(q))});
    const auto velocities = bem_detail::coupledVelocities(basis, solution);
    close(velocities[1].normalVelocity, .5 * velocities[0].normalVelocity, 1e-14,
          "Weighted physical velocities do not preserve their ratio.");
    close(areas[0] * velocities[0].normalVelocity + areas[1] * velocities[1].normalVelocity, q,
          1e-14, "Modal volume flow is not the integral of physical velocity.");
    const auto boundaryPower = areas[0] * pressure[0] * std::conj(velocities[0].normalVelocity) +
                               areas[1] * pressure[1] * std::conj(velocities[1].normalVelocity);
    close(boundaryPower, modal * std::conj(q), 1e-14,
          "Modal pressure and Q do not preserve complex power.");
    bem_detail::validateCouplingClosure(basis, solution, areas, pressure);
    auto scaledNetwork = network;
    std::get<lem::BemPort>(scaledNetwork.elements[1].component).velocityWeights = {{1, 7},
                                                                                   {4, 3.5}};
    const auto scaled = bem_detail::makeCouplingBasis(scaledNetwork, mesh, areas);
    for (std::size_t i = 0; i < 2; ++i)
        close(scaled.panelVelocities[i], basis.panelVelocities[i], 1e-14,
              "Uniform mode-weight scaling changed its physical unit-Q basis.");
    port.velocityWeights.clear();
    const auto uniform = bem_detail::makeCouplingBasis(network, mesh, areas);
    for (auto velocity : uniform.panelVelocities)
        close(velocity, 1 / .03, 1e-14, "Missing weights changed the previous uniform default.");
}
void checkPulsatingSphere() {
    auto mesh = testing::sphere(3);
    double area = 0;
    for (const auto &triangle : mesh.triangles)
        area += triangleArea(mesh, triangle);
    CudaBemSolver solver(mesh);
    BemRequest request;
    request.frequencies = {300};
    request.observers = {{0, 0, 1}};
    request.coupling = flowNetwork({{1}}, {{area, 0}});
    const auto result = solver.solve(request);
    const auto &output = result.frequencies.front().excitations.front();
    const double k = 2 * std::numbers::pi * 300 / 343;
    const auto analytic = 1.21 * 343 * (imaginary * k * .1) / (1.0 + imaginary * k * .1) * .1 *
                          std::exp(-imaginary * k * .9);
    close(output.observerPressure.front(), analytic, .025,
          "Coupled sphere does not match analytic finite-distance pressure");
    close(output.solvedVelocities.front().normalVelocity, {1, 0}, 1e-12,
          "Coupled unit flow did not yield unit velocity");
    require(output.lemSolution.has_value() && output.boundaryResidual <= 5e-5,
            "Coupled sphere has no valid LEM or BEM residual");
}
void checkTwoPorts() {
    auto mesh = testing::sphere(2);
    std::array<double, 2> areas{};
    for (auto &triangle : mesh.triangles) {
        const auto centerZ = mesh.vertices[triangle.vertices[0]].z +
                             mesh.vertices[triangle.vertices[1]].z +
                             mesh.vertices[triangle.vertices[2]].z;
        triangle.faceId = centerZ > 0 ? 1 : 2;
        areas[triangle.faceId - 1] += triangleArea(mesh, triangle);
    }
    CudaBemSolver solver(mesh);
    BemRequest basis;
    basis.frequencies = {300};
    basis.observers = {{0, 0, .5}, {0, 0, -.7}, {.6, .2, 0}};
    basis.excitations = {{"North", {{1, {1 / areas[0], 0}}}}, {"South", {{2, {1 / areas[1], 0}}}}};
    const auto reference = solver.solve(basis);
    const auto &columns = reference.frequencies.front().excitations;
    std::array<Complex, 4> z{};
    for (std::size_t j = 0; j < 2; ++j)
        for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
            const auto row = mesh.triangles[t].faceId - 1;
            z[2 * row + j] +=
                columns[j].boundary()->surfacePressure[t] * triangleArea(mesh, mesh.triangles[t]) / areas[row];
        }
    close(z[1], z[2], 2e-4, "Radiation matrix violates reciprocity on symmetric mesh");
    require(z[0].real() > 0 && z[3].real() > 0 &&
                z[0].real() * z[3].real() - std::pow(.5 * (z[1].real() + z[2].real()), 2) > 0,
            "Radiation resistance matrix is not passive");

    // Unequal drive makes omission of mutual loading observable even with a
    // geometrically symmetric north/south radiator.
    const std::array<Complex, 2> voltages{{{2, .3}, {4, -2}}};
    const double omega = 2 * std::numbers::pi * 300;
    for (const double movingMass : {.02, .0003}) {
        const lem::DriverParameters driver{6, .0003, 1.5, movingMass, .02, .05, .03};
        BemRequest coupled = basis;
        coupled.excitations.clear();
        coupled.coupling = lem::makeDriverNetwork(driver, {{1}, {2}}, voltages[0], 1.0);
        coupled.coupling->nodes.push_back(
            {{100}, lem::Domain::electrical, false, "Second amplifier", {}});
        coupled.coupling->elements.push_back(
            {{100}, "Second voltage", {}, lem::EffortSource{{{100}, {1}}, voltages[1]}});
        std::size_t driverIndex = 0;
        for (auto &element : coupled.coupling->elements)
            if (auto *part = std::get_if<lem::Driver>(&element.component))
                if (driverIndex++ == 1)
                    part->electrical.positive = {100};
        require(driverIndex == 2, "Driver fixture has an unexpected driver count");
        const auto result = solver.solve(coupled);
        const auto &output = result.frequencies.front().excitations.front();
        const auto &solution = *output.lemSolution;
        const Complex ze = driver.resistance + imaginary * omega * driver.inductance;
        const Complex zm =
            driver.mechanicalResistance + imaginary * omega * driver.movingMass +
            1.0 / (imaginary * omega * driver.compliance) +
            driver.pistonArea * driver.pistonArea * 1.21 * 343 * 343 / (imaginary * omega);
        const auto base = zm + driver.forceFactor * driver.forceFactor / ze;
        const double sd2 = driver.pistonArea * driver.pistonArea;
        const Complex a = base + sd2 * z[0], b = sd2 * z[1], c = sd2 * z[2], d = base + sd2 * z[3];
        const Complex f0 = driver.forceFactor * voltages[0] / ze,
                      f1 = driver.forceFactor * voltages[1] / ze;
        const std::array<Complex, 2> velocity{(d * f0 - b * f1) / (a * d - b * c),
                                              (a * f1 - c * f0) / (a * d - b * c)};
        for (std::size_t i = 0; i < 2; ++i) {
            close(solution.drivers[i].velocity, velocity[i], 3e-5,
                  "Loaded driver differs from independent two-port electromechanical equations");
            close(solution.ports[i].volumeVelocity, driver.pistonArea * velocity[i], 3e-5,
                  "Diaphragm-to-BEM transformer area is wrong");
            close(solution.drivers[i].current,
                  (voltages[i] - driver.forceFactor * velocity[i]) / ze, 3e-5,
                  "Back EMF is not coupled to the loaded driver motion");
            Complex measuredPressure{};
            for (std::size_t t = 0; t < mesh.triangles.size(); ++t)
                if (mesh.triangles[t].faceId == i + 1)
                    measuredPressure += output.boundary()->surfacePressure[t] *
                                        triangleArea(mesh, mesh.triangles[t]) / areas[i];
            close(measuredPressure, solution.ports[i].pressure, 8e-5,
                  "Reconstructed boundary pressure does not close the LEM feedback loop");
        }
        for (std::size_t i = 0; i < basis.observers.size(); ++i)
            close(output.observerPressure[i],
                  columns[0].observerPressure[i] * driver.pistonArea * velocity[0] +
                      columns[1].observerPressure[i] * driver.pistonArea * velocity[1],
                  8e-5, "Coupled observer field violates independent BEM superposition");
        require(output.boundaryResidual <= coupled.numerics.relativeResidualTolerance &&
                    solution.relativeResidual < 1e-10,
                "Coupled solution failed its original BEM or MNA equation residual");
        BemRequest prescribed = basis;
        prescribed.excitations = {{"Actual", output.solvedVelocities}};
        const auto independent = solver.solve(prescribed);
        for (std::size_t i = 0; i < basis.observers.size(); ++i)
            close(output.observerPressure[i],
                  independent.frequencies.front().excitations.front().observerPressure[i], 8e-5,
                  "Coupled field differs from independent actual-velocity BEM solve");
        if (movingMass < .001)
            require(std::abs(velocity[0] / (f0 / a) - 1.0) > .05,
                    "Strong-load fixture would not detect missing mutual coupling");
        std::cout << "Two-port loaded driver M=" << movingMass << " kg, BEM residual "
                  << output.boundaryResidual << ", MNA residual " << solution.relativeResidual
                  << '\n';
    }

    BemRequest flows = basis;
    flows.excitations.clear();
    flows.coupling = flowNetwork({{1}, {2}}, {{.001, .0002}, {.0003, -.0001}}, -1);
    const auto polarityResult = solver.solve(flows);
    const auto &polarity = polarityResult.frequencies.front().excitations.front();
    close(polarity.lemSolution->ports[1].volumeVelocity, {-.0003, .0001}, 1e-12,
          "BEM port polarity is applied twice or lost");
    for (std::size_t i = 0; i < basis.observers.size(); ++i)
        close(polarity.observerPressure[i],
              columns[0].observerPressure[i] * Complex(.001, .0002) -
                  columns[1].observerPressure[i] * Complex(.0003, -.0001),
              8e-5, "Port polarity does not reach the physical pressure field");
    flows.coupling = flowNetwork({{1}, {2}}, {{}, {}});
    const auto zeroResult = solver.solve(flows);
    const auto &zero = zeroResult.frequencies.front().excitations.front();
    require(zero.boundaryResidual == 0 && zero.impedances.empty() &&
                zero.solvedVelocities.size() == 2,
            "Zero drive did not preserve finite zero-flow result metadata");
    for (const auto pressure : zero.boundary()->surfacePressure)
        require(pressure == Complex{}, "Zero-drive surface pressure is nonzero");
    for (const auto pressure : zero.observerPressure)
        require(pressure == Complex{}, "Zero-drive observer pressure is nonzero");

    BemRequest weighted = basis;
    weighted.excitations.clear();
    const Complex modeFlow{.001, .0003};
    weighted.coupling = flowNetwork({{1, 2}}, {modeFlow});
    std::get<lem::BemPort>(weighted.coupling->elements[1].component).velocityWeights = {{2, .5}};
    const double modalArea = areas[0] + .5 * areas[1];
    for (const auto acceleration : {BemAcceleration::dense, BemAcceleration::fastMultipole}) {
        weighted.numerics.acceleration = acceleration;
        const auto modeResult = solver.solve(weighted);
        const auto &mode = modeResult.frequencies[0].excitations[0];
        close(mode.solvedVelocities[0].normalVelocity, modeFlow / modalArea, 1e-12,
              "Exterior weighted port loses its first face velocity.");
        close(mode.solvedVelocities[1].normalVelocity, .5 * modeFlow / modalArea, 1e-12,
              "Exterior weighted port loses its velocity ratio.");
        Complex pressure{}, power{};
        for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
            const auto row = mesh.triangles[t].faceId - 1;
            const double weight = row ? .5 : 1;
            const auto area = triangleArea(mesh, mesh.triangles[t]);
            pressure += weight * area * mode.boundary()->surfacePressure[t] / modalArea;
            power += area * mode.boundary()->surfacePressure[t] * std::conj(weight * modeFlow / modalArea);
        }
        close(mode.lemSolution->ports[0].pressure, pressure, relativePortPressureClosureTolerance,
              "Exterior modal feedback pressure is not the power-dual mean.");
        close(power, pressure * std::conj(modeFlow), 1e-12,
              "Exterior weighted mode does not conserve complex acoustic power.");
        require(power.real() > 0 &&
                    mode.boundaryResidual <= weighted.numerics.relativeResidualTolerance,
                "Exterior weighted mode violates passivity or original residual.");
        for (std::size_t o = 0; o < basis.observers.size(); ++o)
            close(mode.observerPressure[o],
                  columns[0].observerPressure[o] * areas[0] * modeFlow / modalArea +
                      columns[1].observerPressure[o] * .5 * areas[1] * modeFlow / modalArea,
                  8e-5,
                  "Dense/FMM weighted exterior field differs from physical independent solves.");
        std::cout << "Exterior modal acceleration=" << static_cast<int>(acceleration)
                  << " BEM residual=" << mode.boundaryResidual << '\n';
    }

    const auto rejects = [&](const BemRequest &invalid) {
        bool rejected = false;
        try {
            (void)solver.solve(invalid);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Invalid BEM port bindings were accepted");
    };
    auto invalid = flows;
    invalid.excitations = basis.excitations;
    rejects(invalid);
    invalid = flows;
    invalid.coupling = flowNetwork({{1}, {99}}, {{1, 0}, {1, 0}});
    rejects(invalid);
    invalid = flows;
    std::get<lem::BemPort>(invalid.coupling->elements[1].component).kind =
        lem::BemPortKind::interfacePort;
    rejects(invalid);
    invalid.coupling = flowNetwork({{1}, {1}}, {{1, 0}, {1, 0}});
    rejects(invalid);
    invalid.coupling = flowNetwork({{1}}, {{1, 0}});
    rejects(invalid);
    for (auto &triangle : mesh.triangles)
        if (triangle.faceId == 2)
            triangle.boundary = BoundaryCondition::wall;
    CudaBemSolver wallSolver(mesh);
    bool wallRejected = false;
    try {
        (void)wallSolver.solve(flows);
    } catch (const std::invalid_argument &) {
        wallRejected = true;
    }
    require(wallRejected, "A Wall face was accepted as a driven BEM port");
    std::stop_source stop;
    stop.request_stop();
    bool cancelled = false;
    try {
        (void)solver.solve(flows, stop.get_token());
    } catch (const std::runtime_error &) {
        cancelled = true;
    }
    require(cancelled, "Coupled BEM solve ignored cancellation");
}
} // namespace
int main() {
    try {
        std::cout << "CUDA device: " << CudaBemSolver::deviceDescription() << '\n';
        checkClosureContract();
        checkModalBasisContract();
        checkPulsatingSphere();
        checkTwoPorts();
        std::cout << "LEM-BEM coupling: analytic sphere, mutual loads, feedback, polarity, zero "
                     "and validation PASS\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
