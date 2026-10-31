#include <algorithm>
#include <cmath>
#include <complex>
#include <hornsculpt/axisymmetric_mesh.hpp>
#include <hornsculpt/bem.hpp>
#include <hornsculpt/field_evaluation.hpp>
#include <hornsculpt/meshing.hpp>
#include <hornsculpt/surface_measures.hpp>
#include <iostream>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string_view>
#include <vector>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
// A closed air cylinder. Its bottom piston normal points out of air (-Z),
// whereas the public positive piston velocity enters air (+Z).
SurfaceMesh cylinder(unsigned sectors, unsigned radial, double depth) {
    constexpr double radius = .025;
    SurfaceMesh mesh;
    for (unsigned layer = 0; layer < 2; ++layer) {
        const double z = layer == 0 ? -depth : 0;
        mesh.vertices.push_back({0, 0, z});
        for (unsigned r = 1; r <= radial; ++r)
            for (unsigned i = 0; i < sectors; ++i) {
                const double angle = 2 * std::numbers::pi * i / sectors;
                mesh.vertices.push_back({radius * r / radial * std::cos(angle),
                                         radius * r / radial * std::sin(angle), z});
            }
    }
    const unsigned perLayer = 1 + sectors * radial;
    auto node = [&](unsigned layer, unsigned r, unsigned i) {
        return layer * perLayer + (r == 0 ? 0 : 1 + (r - 1) * sectors + i % sectors);
    };
    auto tri = [&](unsigned a, unsigned b, unsigned c, unsigned face, bool driven) {
        mesh.triangles.push_back(
            {{a, b, c}, face, driven ? BoundaryCondition::driven : BoundaryCondition::wall});
    };
    for (unsigned layer = 0; layer < 2; ++layer)
        for (unsigned r = 1; r <= radial; ++r)
            for (unsigned i = 0; i < sectors; ++i) {
                const unsigned a = node(layer, r - 1, i), b = node(layer, r, i),
                               c = node(layer, r, i + 1);
                const unsigned face = layer == 0 ? 1 : 3;
                if (layer == 0)
                    tri(a, c, b, face, true);
                else
                    tri(a, b, c, face, false);
                if (r > 1) {
                    const auto d = node(layer, r - 1, i + 1);
                    if (layer == 0)
                        tri(a, d, c, face, true);
                    else
                        tri(a, c, d, face, false);
                }
            }
    for (unsigned i = 0; i < sectors; ++i) {
        const auto a = node(0, radial, i), b = node(0, radial, i + 1), c = node(1, radial, i + 1),
                   d = node(1, radial, i);
        tri(a, b, c, 2, false);
        tri(a, c, d, 2, false);
    }
    return mesh;
}
BemRequest request() {
    BemRequest value;
    value.infiniteBaffle = InfiniteBaffle{{0, 0, 0}, {0, 0, 1}, {3}};
    value.frequencies = {500, 2000};
    value.excitations = {{"Piston", {{1, {1, 0}}}}, {"Complex piston", {{1, {.3, -.4}}}}};
    value.farField = FarFieldObservation{};
    value.observers = {{0, 0, 1}, {1, 0, 1}, {0, 0, -1}};
    value.numerics.singularOrder = 6;
    value.numerics.nearRelativeTolerance = 1e-4;
    value.numerics.nearMaximumDepth = 8;
    return value;
}
double triangleArea(const SurfaceMesh &mesh, const Triangle &t) {
    const auto a = mesh.vertices[t.vertices[0]], b = mesh.vertices[t.vertices[1]],
               c = mesh.vertices[t.vertices[2]];
    const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
    return .5 * std::hypot(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
}

// Independent mapping for the older closed cylindrical fixture.
AxisymmetricBem cylinderReduction(const SurfaceMesh &mesh, unsigned sectors, unsigned radial) {
    AxisymmetricBem result;
    result.sectors = sectors;
    const unsigned perLayer = 1 + sectors * radial;
    for (unsigned vertex = 0; vertex < mesh.vertices.size(); ++vertex) {
        const unsigned layer = vertex / perLayer, local = vertex % perLayer;
        result.pressureDofs.push_back(layer * (radial + 1) +
                                      (local ? 1 + (local - 1) / sectors : 0));
    }
    std::map<std::pair<unsigned, unsigned>, unsigned> strips;
    for (unsigned t = 0; t < mesh.triangles.size(); ++t) {
        std::set<unsigned> nodes, angles;
        for (auto vertex : mesh.triangles[t].vertices) {
            nodes.insert(result.pressureDofs[vertex]);
            const auto local = vertex % perLayer;
            if (local)
                angles.insert((local - 1) % sectors);
        }
        require(nodes.size() == 2, "Cylinder fixture is not a meridian strip mesh.");
        const auto key = std::pair{*nodes.begin(), *nodes.rbegin()};
        const auto next = static_cast<unsigned>(strips.size());
        result.fluxDofs.push_back(strips.emplace(key, next).first->second);
        if (angles == std::set<unsigned>{0, 1})
            result.testPanels.push_back(t);
    }
    return result;
}

void checkTrace(const SurfaceMesh &mesh, const BemExcitationResult &value,
                const AcousticMedium &medium, double frequency, unsigned apertureFace = 3,
                const std::optional<AxisymmetricBem> &axisymmetric = {}) {
    require(value.boundary()->pressureTrace.size() == mesh.vertices.size() &&
                value.boundary()->normalDerivative.size() == mesh.triangles.size(),
            "Solved Cauchy trace does not cover the full mesh.");
    std::map<unsigned, std::complex<double>> velocities;
    const auto measures = surfaceMeasures(mesh, axisymmetric);
    for (const auto &face : value.solvedVelocities)
        velocities.emplace(face.faceId, face.normalVelocity);
    for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
        const auto &triangle = mesh.triangles[t];
        std::complex<double> mean{};
        for (unsigned local = 0; local < 3; ++local)
            if (axisymmetric &&
                axisymmetric->surfaceContract == AxisymmetricSurfaceContract::continuousRevolution)
                mean += value.boundary()->pressureTrace[triangle.vertices[local]] *
                        (measures.pressureIntegrals[t][local] / measures.areas[t]);
            else
                mean += value.boundary()->pressureTrace[triangle.vertices[local]] / 3.0;
        require(mean == value.boundary()->surfacePressure[t], "Triangle mean loses its solved P1 trace.");
        if (triangle.faceId != apertureFace) {
            const auto found = velocities.find(triangle.faceId);
            const auto q =
                found == velocities.end()
                    ? std::complex<double>{}
                    : std::complex<double>{0, 2 * std::numbers::pi * frequency * medium.density} *
                          found->second;
            require(std::abs(q - value.boundary()->normalDerivative[t]) <= 1e-14 * std::max(1.0, std::abs(q)),
                    "Physical derivative trace uses an incorrect air-cavity normal.");
        }
    }
}
SurfaceMesh segmentedCylinder() {
    auto mesh = cylinder(16, 3, .006);
    for (auto &t : mesh.triangles) {
        if (t.faceId != 1)
            continue;
        double radius{}, x{};
        for (auto node : t.vertices) {
            const auto p = mesh.vertices[node];
            radius = std::max(radius, std::hypot(p.x, p.y));
            x += p.x;
        }
        if (radius > .025 * 2 / 3 + 1e-12) {
            t.faceId = x > 0 ? 5 : 6;
            t.boundary = BoundaryCondition::wall;
        } else if (radius > .025 / 3 + 1e-12)
            t.faceId = 4;
    }
    return mesh;
}
struct Sample {
    double pistonError{};
    std::vector<std::complex<double>> pressure;
};
Sample run(unsigned sectors, unsigned radial, double depth, unsigned singularOrder = 6) {
    auto mesh = cylinder(sectors, radial, depth);
    const auto quality = inspectMesh(mesh);
    require(quality.watertight() && quality.signedVolume > 0,
            "Cylinder fixture is not a closed air volume.");
    auto input = request();
    input.numerics.singularOrder = singularOrder;
    const auto result = CudaBemSolver(mesh).solve(input);
    require(result.frequencies.size() == 2 && result.pressureDofCount == mesh.vertices.size(),
            "Baffled result shape mismatch.");
    Sample sample;
    for (std::size_t f = 0; f < result.frequencies.size(); ++f) {
        const auto &frequency = result.frequencies[f];
        const auto &one = frequency.excitations[0], &complex = frequency.excitations[1];
        checkTrace(mesh, one, input.medium, frequency.frequency);
        require(one.boundaryResidual < 1e-10, "Baffled original block residual is too large.");
        require(one.boundary()->surfacePressure.size() == mesh.triangles.size(),
                "Baffled triangle pressures are missing.");
        require(one.observerPressure[2] == std::complex<double>{},
                "Baffle rear radiation was not screened.");
        require(one.impedances.size() == 1 && one.impedances[0].volumeImpedance.real() > 0,
                "Baffled piston has nonpositive radiation resistance.");
        for (std::size_t o = 0; o < one.observerPressure.size(); ++o)
            require(std::abs(complex.observerPressure[o] -
                             std::complex<double>{.3, -.4} * one.observerPressure[o]) <
                        1e-10 * std::max(1.0, std::abs(one.observerPressure[o])),
                    "Baffled independent complex RHS linearity failed.");
        const double k = 2 * std::numbers::pi * frequency.frequency / input.medium.soundSpeed;
        // Equal-area circular piston removes the polygon's known area error.
        const double area = one.impedances[0].area;
        const std::complex<double> exact =
            std::complex<double>{0, input.medium.density * 2 * std::numbers::pi *
                                        frequency.frequency * area / (2 * std::numbers::pi)} *
            std::exp(std::complex<double>{0, -k});
        const double error = std::abs(one.observerPressure[0] / exact - 1.0);
        sample.pistonError = std::max(sample.pistonError, error);
        sample.pressure.push_back(one.observerPressure[0]);
        std::cout << "cylinder " << sectors << '/' << radial << " h=" << depth
                  << " f=" << frequency.frequency << " piston error=" << error
                  << " residual=" << one.boundaryResidual
                  << " Z=" << one.impedances[0].volumeImpedance << '\n';
    }
    return sample;
}
double difference(const Sample &a, const Sample &b) {
    double largest{};
    for (std::size_t i = 0; i < a.pressure.size(); ++i)
        largest = std::max(largest, std::abs(a.pressure[i] / b.pressure[i] - 1.0));
    return largest;
}
void finiteField() {
    const auto mesh = cylinder(24, 3, .001);
    auto input = request();
    input.frequencies = {2000};
    CudaBemSolver solver(mesh);
    const auto far = solver.solve(input).frequencies[0].excitations[0].observerPressure[0];
    input.farField.reset();
    input.observers = {{0, 0, 10}, {0, 0, .0005}};
    const auto finite = solver.solve(input).frequencies[0].excitations[0];
    const double k = 2 * std::numbers::pi * 2000 / input.medium.soundSpeed;
    const auto normalized =
        10.0 * finite.observerPressure[0] * std::exp(std::complex<double>{0, k * 9});
    require(std::abs(normalized / far - 1.0) < .002,
            "Finite baffled field does not approach its asymptotic field.");
    require(std::isfinite(finite.observerPressure[1].real()) &&
                std::abs(finite.observerPressure[1]) > 0,
            "Adaptive near-aperture Rayleigh field is missing.");
    // Independently integrate the actual saved CUDA Cauchy trace. This guards
    // aperture q/k reconstruction and the -2 Rayleigh factor at both a remote
    // and a genuinely near-aperture observation point.
    const auto reconstructed = evaluatePressureField(mesh, input, input.frequencies[0], finite,
                                                     input.observers, {1e-4, 1e-5, 10});
    require(reconstructed.size() == finite.observerPressure.size(),
            "Saved-trace field lost a finite CUDA observation point.");
    for (std::size_t i = 0; i < reconstructed.size(); ++i)
        require(reconstructed[i].status == FieldPointStatus::air &&
                    std::abs(reconstructed[i].pressure - finite.observerPressure[i]) <=
                        3e-4 * std::max(std::abs(finite.observerPressure[i]), 1e-12),
                "Saved CUDA Cauchy trace does not reproduce its finite Rayleigh field.");
}
void rigidPlacement() {
    auto mesh = cylinder(12, 2, .005);
    // Preserve exactly the same binary geometry across the rigid placement.
    // Otherwise a translation changes the last coordinate bits as well as pose.
    for (auto &p : mesh.vertices) {
        p.x = std::ldexp(std::round(std::ldexp(p.x, 30)), -30);
        p.y = std::ldexp(std::round(std::ldexp(p.y, 30)), -30);
        p.z = std::ldexp(std::round(std::ldexp(p.z, 30)), -30);
    }
    auto input = request();
    // Keep this pose test away from exact near/far cutoff equalities.
    input.numerics.nearDistanceRatio = 1.3;
    const auto original = CudaBemSolver(mesh).solve(input);
    auto lowerInput = input, higherInput = input;
    lowerInput.numerics.singularOrder = 4;
    higherInput.numerics.singularOrder = 8;
    const auto lower = CudaBemSolver(mesh).solve(lowerInput),
               higher = CudaBemSolver(mesh).solve(higherInput);
    double lowerChange{}, higherChange{};
    for (std::size_t f = 0; f < original.frequencies.size(); ++f) {
        const auto reference = original.frequencies[f].excitations[0].observerPressure[0];
        lowerChange = std::max(
            lowerChange,
            std::abs(lower.frequencies[f].excitations[0].observerPressure[0] / reference - 1.0));
        higherChange = std::max(
            higherChange,
            std::abs(higher.frequencies[f].excitations[0].observerPressure[0] / reference - 1.0));
    }
    std::cout << "Independent singular-order changes 4/6=" << lowerChange << " 6/8=" << higherChange
              << '\n';
    require(higherChange < lowerChange && higherChange < .003,
            "Baffled singular quadrature did not converge independently of the mesh.");
    const Vec3 translation{.125, -.25, .5};
    auto placed = [&](Vec3 p) {
        return Vec3{p.z + translation.x, p.y + translation.y, -p.x + translation.z};
    };
    for (auto &p : mesh.vertices)
        p = placed(p);
    for (auto &p : input.observers)
        p = placed(p);
    input.infiniteBaffle->origin = translation;
    input.infiniteBaffle->normal = {1, 0, 0};
    input.farField->origin = translation;
    const auto transformed = CudaBemSolver(mesh).solve(input);
    for (std::size_t f = 0; f < original.frequencies.size(); ++f)
        for (std::size_t r = 0; r < original.frequencies[f].excitations.size(); ++r) {
            const auto &a = original.frequencies[f].excitations[r],
                       &b = transformed.frequencies[f].excitations[r];
            for (std::size_t o = 0; o < a.observerPressure.size(); ++o) {
                if (std::abs(a.observerPressure[o] - b.observerPressure[o]) >=
                    1e-9 * std::max(1.0, std::abs(a.observerPressure[o])))
                    std::cerr << "Placement mismatch f=" << f << " rhs=" << r << " observer=" << o
                              << " original=" << a.observerPressure[o]
                              << " transformed=" << b.observerPressure[o] << " relative="
                              << std::abs(a.observerPressure[o] - b.observerPressure[o]) /
                                     std::max(1.0, std::abs(a.observerPressure[o]))
                              << '\n';
                require(std::abs(a.observerPressure[o] - b.observerPressure[o]) <
                            1e-9 * std::max(1.0, std::abs(a.observerPressure[o])),
                        "Baffled plane/phase origin rotation or translation changed the field.");
            }
            if (std::abs(a.impedances[0].volumeImpedance / b.impedances[0].volumeImpedance - 1.0) >=
                1e-9)
                std::cerr << "Placement impedance relative="
                          << std::abs(a.impedances[0].volumeImpedance /
                                          b.impedances[0].volumeImpedance -
                                      1.0)
                          << '\n';
            require(std::abs(a.impedances[0].volumeImpedance / b.impedances[0].volumeImpedance -
                             1.0) < 1e-9,
                    "Baffled rigid placement changed the physical port impedance.");
        }
}
void interfaceCoupling() {
    using Complex = std::complex<double>;
    constexpr Complex imaginary{0, 1};
    auto mesh = segmentedCylinder();
    auto panelArea = [&](const Triangle &t) { return triangleArea(mesh, t); };
    std::map<std::uint32_t, double> areas;
    for (auto &t : mesh.triangles) {
        areas[t.faceId] += panelArea(t);
    }
    require(areas[1] > 0 && areas[4] > 0 && areas[5] > 0 && areas[6] > 0,
            "Interface fixture misses a complete driven annulus or interface group.");
    auto prescribedMesh = mesh;
    for (auto &t : prescribedMesh.triangles)
        if (t.faceId == 5 || t.faceId == 6)
            t.boundary = BoundaryCondition::driven;
    auto input = request();
    input.frequencies = {400, 1900};
    input.excitations = {{"Dome and whole surround", {{1, {1, 0}}, {4, {.5, 0}}}},
                         {"Unit interface 1", {{5, {1 / areas[5], 0}}}},
                         {"Unit interface 2", {{6, {1 / areas[6], 0}}}}};
    const auto reference = CudaBemSolver(prescribedMesh).solve(input);
    auto mean = [&](const BemExcitationResult &result, std::uint32_t face) {
        Complex pressure{};
        for (std::size_t t = 0; t < mesh.triangles.size(); ++t)
            if (mesh.triangles[t].faceId == face)
                pressure += panelArea(mesh.triangles[t]) * result.boundary()->surfacePressure[t] / areas[face];
        return pressure;
    };
    auto close = [&](Complex actual, Complex expected, double tolerance, const char *message) {
        if (!(std::abs(actual - expected) <= tolerance * std::max(1e-20, std::abs(expected)))) {
            std::cerr << message << " actual=" << actual << " expected=" << expected << '\n';
            throw std::runtime_error(message);
        }
    };
    auto network = [&](double volumeScale, bool losses) {
        lem::Network value;
        value.nodes = {{{1}, lem::Domain::acoustic, true, {}, {}},
                       {{2}, lem::Domain::acoustic, false, {}, {}},
                       {{3}, lem::Domain::acoustic, false, {}, {}}};
        for (std::uint32_t i = 0; i < 2; ++i) {
            const lem::Terminals terminals{{i + 2}, {1}};
            value.elements.push_back({{3 * i + 1},
                                      "Rear volume",
                                      {},
                                      lem::Cavity{terminals, volumeScale * (i ? 7e-6 : 5e-6)}});
            value.elements.push_back(
                {{3 * i + 2},
                 "Interface",
                 {},
                 lem::BemPort{terminals, {5 + i}, 1, lem::BemPortKind::interfacePort}});
            if (losses)
                value.elements.push_back(
                    {{3 * i + 3}, "Rear loss", {}, lem::Resistance{terminals, i ? 3e7 : 2e7}});
        }
        return value;
    };
    input.excitations.resize(1);
    input.coupling = network(1, true);
    CudaBemSolver solver(mesh);
    const auto coupled = solver.solve(input);
    double mutualEffect{};
    for (std::size_t f = 0; f < coupled.frequencies.size(); ++f) {
        const auto &columns = reference.frequencies[f].excitations;
        const auto &result = coupled.frequencies[f].excitations.at(0);
        require(coupled.frequencies[f].excitations.size() == 1 && result.name == "LEM" &&
                    result.lemSolution && result.solvedVelocities.size() == 4,
                "Coupled output omits prescribed or interface velocities.");
        const auto &ports = result.lemSolution->ports;
        const Complex p0[]{mean(columns[0], 5), mean(columns[0], 6)},
            z[]{mean(columns[1], 5), mean(columns[2], 5), mean(columns[1], 6), mean(columns[2], 6)};
        const double omega = 2 * std::numbers::pi * input.frequencies[f];
        const double bulk =
            input.medium.density * input.medium.soundSpeed * input.medium.soundSpeed;
        const Complex y[]{1.0 / 2e7 + imaginary * omega * 5e-6 / bulk,
                          1.0 / 3e7 + imaginary * omega * 7e-6 / bulk};
        const Complex a = 1.0 + y[0] * z[0], b = y[0] * z[1], c = y[1] * z[2],
                      d = 1.0 + y[1] * z[3], r0 = -y[0] * p0[0], r1 = -y[1] * p0[1];
        const Complex q[]{(d * r0 - b * r1) / (a * d - b * c), (a * r1 - c * r0) / (a * d - b * c)};
        mutualEffect = std::max(mutualEffect, std::abs(q[0] / (r0 / a) - 1.0));
        close(z[1], z[2], .004, "Baffled interface mutual impedance violates reciprocity.");
        double rearPower{};
        for (std::size_t i = 0; i < 2; ++i) {
            close(ports[i].volumeVelocity, q[i], 1e-9, "Analytic rear-load feedback Q mismatch.");
            const auto pressure = p0[i] + z[2 * i] * q[0] + z[2 * i + 1] * q[1];
            close(ports[i].pressure, pressure, 1e-9, "Per-frequency blocked pressure is lost.");
            close(mean(result, 5 + static_cast<std::uint32_t>(i)), pressure, 1e-9,
                  "Reconstructed BEM and LEM interface pressures differ.");
            close(q[i], -y[i] * pressure, 1e-9, "Interface flux continuity has the wrong sign.");
            rearPower += ports[i].activePower + std::norm(pressure) / (i ? 3e7 : 2e7);
        }
        require(std::abs(rearPower) < 1e-10 && result.boundaryResidual < 1e-10 &&
                    result.lemSolution->relativeResidual < 1e-10,
                "Coupled original block residual or passive rear power does not close.");
        double frontPower{};
        for (const auto &velocity : result.solvedVelocities) {
            if (velocity.faceId == 1)
                close(velocity.normalVelocity, {1, 0}, 1e-12, "Prescribed dome changed.");
            if (velocity.faceId == 4)
                close(velocity.normalVelocity, {.5, 0}, 1e-12, "Whole surround changed.");
            frontPower += std::real(mean(result, velocity.faceId) *
                                    std::conj(velocity.normalVelocity * areas[velocity.faceId]));
        }
        require(frontPower > 0, "Passive front BEM domain creates acoustic power.");
        for (std::size_t o = 0; o < input.observers.size(); ++o)
            close(result.observerPressure[o],
                  columns[0].observerPressure[o] + columns[1].observerPressure[o] * q[0] +
                      columns[2].observerPressure[o] * q[1],
                  1e-9,
                  "Coupled Rayleigh field does not include known drive and both interface flows.");
        std::cout << "interface f=" << input.frequencies[f] << " Q=" << q[0] << ',' << q[1]
                  << " front power=" << frontPower << " residual=" << result.boundaryResidual
                  << '\n';
    }
    require(mutualEffect > .005, "Interface fixture does not exercise material mutual coupling.");
    auto phased = input;
    constexpr Complex phase{.3, -.4};
    for (auto &face : phased.excitations[0].faces)
        face.normalVelocity *= phase;
    const auto phaseResult = solver.solve(phased);
    for (std::size_t f = 0; f < phaseResult.frequencies.size(); ++f) {
        const auto &actual = phaseResult.frequencies[f].excitations[0];
        const auto &original = coupled.frequencies[f].excitations[0];
        for (std::size_t i = 0; i < 2; ++i) {
            close(actual.lemSolution->ports[i].volumeVelocity,
                  phase * original.lemSolution->ports[i].volumeVelocity, 1e-9,
                  "Complex prescribed phase is lost in interface flow feedback.");
            close(actual.lemSolution->ports[i].pressure,
                  phase * original.lemSolution->ports[i].pressure, 1e-9,
                  "Complex prescribed phase is lost in interface pressure feedback.");
        }
        close(actual.observerPressure[0], phase * original.observerPressure[0], 1e-9,
              "Complex prescribed phase is lost in coupled Rayleigh radiation.");
    }
    input.coupling = network(1e-7, false);
    const auto tiny = solver.solve(input);
    for (std::size_t f = 0; f < tiny.frequencies.size(); ++f) {
        const auto &output = tiny.frequencies[f].excitations[0];
        for (std::size_t i = 0; i < 2; ++i)
            close(output.lemSolution->ports[i].pressure,
                  mean(reference.frequencies[f].excitations[0], 5 + static_cast<std::uint32_t>(i)),
                  1e-6, "Vanishing sealed volume does not approach blocked pressure.");
        close(output.observerPressure[0],
              reference.frequencies[f].excitations[0].observerPressure[0], 1e-6,
              "Vanishing sealed volume does not approach rigid-interface field.");
    }
    // Interface-only forcing also works without a prescribed physical surface.
    auto flowNetwork = network(1, false);
    flowNetwork.elements.erase(
        std::remove_if(
            flowNetwork.elements.begin(), flowNetwork.elements.end(),
            [](const auto &e) { return std::holds_alternative<lem::Cavity>(e.component); }),
        flowNetwork.elements.end());
    const Complex forced[]{Complex{2e-5, -1e-5}, Complex{-1e-5, 3e-5}};
    for (std::uint32_t i = 0; i < 2; ++i)
        flowNetwork.elements.push_back(
            {{20 + i}, "Flow source", {}, lem::FlowSource{{{1}, {i + 2}}, forced[i]}});
    input.excitations.clear();
    input.coupling = flowNetwork;
    const auto forcedResult = solver.solve(input);
    for (std::size_t f = 0; f < forcedResult.frequencies.size(); ++f) {
        const auto &output = forcedResult.frequencies[f].excitations[0];
        for (std::size_t i = 0; i < 2; ++i)
            close(output.lemSolution->ports[i].volumeVelocity, forced[i], 1e-12,
                  "Interface-only source loses physical volume flow.");
        close(output.observerPressure[0],
              reference.frequencies[f].excitations[1].observerPressure[0] * forced[0] +
                  reference.frequencies[f].excitations[2].observerPressure[0] * forced[1],
              1e-9, "Interface-only complex source does not reach Rayleigh field.");
    }
    auto zeroPrescribed = input;
    zeroPrescribed.excitations = {{"Silent physical surfaces", {{1, {}}, {4, {}}}}};
    const auto withZeroSurfaces = solver.solve(zeroPrescribed);
    for (std::size_t f = 0; f < withZeroSurfaces.frequencies.size(); ++f) {
        const auto &actual = withZeroSurfaces.frequencies[f].excitations[0];
        require(actual.solvedVelocities.size() == 4,
                "Zero prescribed velocities are absent from the coupled metadata.");
        close(actual.observerPressure[0],
              forcedResult.frequencies[f].excitations[0].observerPressure[0], 1e-9,
              "Zero prescribed velocity changes an interface-driven field.");
    }
    for (auto &e : flowNetwork.elements)
        if (auto *source = std::get_if<lem::FlowSource>(&e.component))
            source->amplitude = {};
    input.coupling = flowNetwork;
    const auto zero = solver.solve(input).frequencies[0].excitations[0];
    require(zero.boundaryResidual == 0 && zero.impedances.empty() &&
                zero.solvedVelocities.size() == 2,
            "Zero interface drive loses finite zero-flow metadata.");
    for (auto p : zero.boundary()->surfacePressure)
        require(p == Complex{}, "Zero interface drive creates surface pressure.");
    for (auto p : zero.observerPressure)
        require(p == Complex{}, "Zero interface drive creates radiation.");

    auto rejects = [&](BemRequest bad) {
        bool rejected{};
        try {
            (void)solver.solve(bad);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Invalid subdomain interface binding was accepted.");
    };
    auto bad = input;
    bad.excitations = {{"one", {{1, {1, 0}}}}, {"two", {{4, {.5, 0}}}}};
    rejects(bad);
    bad = input;
    std::get<lem::BemPort>(bad.coupling->elements[0].component).faceIds = {3};
    rejects(bad);
    bad = input;
    std::get<lem::BemPort>(bad.coupling->elements[0].component).faceIds = {1};
    rejects(bad);
    bad = input;
    std::get<lem::BemPort>(bad.coupling->elements[0].component).kind =
        lem::BemPortKind::drivenSurface;
    rejects(bad);
}
void modalDriverCoupling() {
    using Complex = std::complex<double>;
    constexpr Complex imaginary{0, 1};
    auto mesh = segmentedCylinder();
    std::map<std::uint32_t, double> areas;
    for (auto &t : mesh.triangles) {
        // This test uses one complete axisymmetric interface. The separate
        // interfaceCoupling case retains its two independent azimuthal faces.
        if (t.faceId == 6)
            t.faceId = 5;
        if (t.faceId == 2)
            t.boundary = BoundaryCondition::driven;
        areas[t.faceId] += triangleArea(mesh, t);
    }
    double sw = areas[1] + .5 * areas[4], interfaceArea = areas[5];
    auto integrationAreas = surfaceMeasures(mesh).areas;
    auto independentMesh = mesh;
    for (auto &t : independentMesh.triangles)
        if (t.faceId == 5)
            t.boundary = BoundaryCondition::driven;
    auto basisRequest = request();
    basisRequest.frequencies = {400, 1900, 6000};
    // Rotated DOF numbering can permute the finite singular quadrature. Use
    // order eight for the tighter 1e-8 H/V covariance check, independently of
    // the production default and the existing order 4/6/8 convergence fixture.
    basisRequest.numerics.singularOrder = 8;
    for (double angle : {30.0, 60.0}) {
        const double radians = angle * std::numbers::pi / 180;
        const double sine = std::sin(radians), cosine = std::cos(radians);
        for (auto direction : {Vec3{sine, 0, cosine}, Vec3{-sine, 0, cosine}, Vec3{0, sine, cosine},
                               Vec3{0, -sine, cosine}})
            basisRequest.observers.push_back(direction);
    }
    basisRequest.excitations = {{"Modal Q", {{1, {1 / sw, 0}}, {4, {.5 / sw, 0}}}},
                                {"Interface Q", {{5, {1 / interfaceArea, 0}}}},
                                {"Independent known source", {{2, {.02, -.01}}}}};
    auto modalMean = [&](const BemExcitationResult &value, bool interfacePort) {
        Complex sum{};
        for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
            const auto face = mesh.triangles[t].faceId;
            double density{};
            if (interfacePort && face == 5)
                density = 1 / interfaceArea;
            if (!interfacePort && (face == 1 || face == 4))
                density = (face == 4 ? .5 : 1) / sw;
            sum += density * integrationAreas[t] * value.boundary()->surfacePressure[t];
        }
        return sum;
    };
    auto close = [&](Complex actual, Complex expected, double tolerance, const char *message) {
        if (!(std::abs(actual - expected) <= tolerance * std::max(1e-20, std::abs(expected)))) {
            std::cerr << message << " actual=" << actual << " expected=" << expected << '\n';
            throw std::runtime_error(message);
        }
    };
    const lem::DriverParameters parameters{5.6, 50e-6, 1.8, .00008, .0002, .04, sw};
    const Complex voltage{2.83, .2};
    lem::Network network;
    network.nodes = {
        {{1}, lem::Domain::electrical, true, {}, {}}, {{2}, lem::Domain::electrical, false, {}, {}},
        {{3}, lem::Domain::acoustic, true, {}, {}},   {{4}, lem::Domain::acoustic, false, {}, {}},
        {{5}, lem::Domain::acoustic, false, {}, {}},  {{6}, lem::Domain::acoustic, false, {}, {}}};
    lem::BemPort driven{{{4}, {3}}, {1, 4}, 1};
    driven.velocityWeights = {{4, .5}};
    const lem::AcousticSlit slit{{{5}, {6}}, .002, 6e-5, .0003, 1.84e-5};
    const lem::ClosedSlitCavity closed{{{6}, {3}}, .03, 5e-6 / .03, .00175, 1.84e-5, 1.4, .71};
    network.elements = {
        {{1}, "Voltage", {}, lem::EffortSource{{{2}, {1}}, voltage}},
        {{2}, "Motor", {}, lem::Driver{{{2}, {1}}, {4}, {3}, parameters}},
        {{3}, "Dome and whole surround", {}, driven},
        {{4}, "Interface", {}, lem::BemPort{{{5}, {3}}, {5}, 1, lem::BemPortKind::interfacePort}},
        {{5}, "Series slit", {}, slit},
        {{6}, "Closed rear slit", {}, closed}};
    for (bool axisymmetric : {false, true}) {
        auto modeBasisRequest = basisRequest;
        if (axisymmetric)
            modeBasisRequest.axisymmetric = cylinderReduction(mesh, 16, 3);
        integrationAreas = surfaceMeasures(mesh, modeBasisRequest.axisymmetric).areas;
        areas.clear();
        for (std::size_t t = 0; t < mesh.triangles.size(); ++t)
            areas[mesh.triangles[t].faceId] += integrationAreas[t];
        sw = areas[1] + .5 * areas[4];
        interfaceArea = areas[5];
        modeBasisRequest.excitations[0].faces = {{1, {1 / sw, 0}}, {4, {.5 / sw, 0}}};
        modeBasisRequest.excitations[1].faces = {{5, {1 / interfaceArea, 0}}};
        std::get<lem::Driver>(network.elements[1].component).parameters.pistonArea = sw;
        const auto basis = CudaBemSolver(independentMesh).solve(modeBasisRequest);
        auto input = modeBasisRequest;
        input.excitations.clear();
        input.coupling = network;
        CudaBemSolver solver(mesh);
        for (bool known : {false, true}) {
            input.excitations = known ? std::vector<BemExcitation>{basisRequest.excitations[2]}
                                      : std::vector<BemExcitation>{};
            const auto solved = solver.solve(input);
            for (std::size_t f = 0; f < solved.frequencies.size(); ++f) {
                const auto &columns = basis.frequencies[f].excitations;
                const auto &actual = solved.frequencies[f].excitations[0];
                checkTrace(mesh, actual, input.medium, solved.frequencies[f].frequency, 3,
                           input.axisymmetric);
                require(solved.geometryMode == (axisymmetric ? BemGeometryMode::axisymmetric
                                                             : BemGeometryMode::threeDimensional),
                        "Mixed motor result loses its selected Galerkin space.");
                const auto &solution = *actual.lemSolution;
                const Complex zdd = modalMean(columns[0], false),
                              zdi = modalMean(columns[1], false), zid = modalMean(columns[0], true),
                              zii = modalMean(columns[1], true),
                              p0d = known ? modalMean(columns[2], false) : Complex{},
                              p0i = known ? modalMean(columns[2], true) : Complex{};
                const double omega = 2 * std::numbers::pi * input.frequencies[f],
                             bulk = input.medium.density * input.medium.soundSpeed *
                                    input.medium.soundSpeed;
                const auto chi =
                    (slit.gap / 2) *
                    std::sqrt(imaginary * omega * input.medium.density / slit.dynamicViscosity);
                const auto lineChi =
                    (closed.gap / 2) *
                    std::sqrt(imaginary * omega * input.medium.density / closed.dynamicViscosity);
                const auto thermalChi = lineChi * std::sqrt(closed.prandtlNumber);
                const Complex lineDensity =
                                  input.medium.density / (1.0 - std::tanh(lineChi) / lineChi),
                              lineBeta = (1.0 + (closed.heatCapacityRatio - 1) *
                                                    std::tanh(thermalChi) / thermalChi) /
                                         bulk,
                              lineK = omega * std::sqrt(lineDensity * lineBeta),
                              lineZc = std::sqrt(lineDensity / lineBeta) / closed.area,
                              lineZ = -imaginary * lineZc / std::tan(lineK * closed.length);
                const Complex slitZ = imaginary * omega * input.medium.density * slit.length /
                                      slit.area / (1.0 - std::tanh(chi) / chi),
                              volumeY = 1.0 / lineZ, y = 1.0 / (slitZ + 1.0 / volumeY),
                              beta = -y * zid * sw / (1.0 + y * zii),
                              gamma = -y * p0i / (1.0 + y * zii),
                              ze =
                                  parameters.resistance + imaginary * omega * parameters.inductance,
                              zm = parameters.mechanicalResistance +
                                   imaginary * omega * parameters.movingMass +
                                   1.0 / (imaginary * omega * parameters.compliance);
                const Complex velocity =
                    (parameters.forceFactor * voltage / ze - sw * (p0d + zdi * gamma)) /
                    (zm + parameters.forceFactor * parameters.forceFactor / ze +
                     sw * (zdd * sw + zdi * beta));
                const Complex qd = sw * velocity, qi = beta * velocity + gamma,
                              pd = p0d + zdd * qd + zdi * qi, pi = p0i + zid * qd + zii * qi;
                close(solution.drivers[0].velocity, velocity, 1e-8,
                      "Electrical driver misses weighted BEM feedback.");
                close(solution.ports[0].volumeVelocity, qd, 1e-8,
                      "Modal port loses integrated physical Q.");
                close(solution.ports[1].volumeVelocity, qi, 1e-8,
                      "Passive interface Q violates the analytic load.");
                close(solution.ports[0].pressure, pd, 1e-8,
                      "Weighted motor pressure is not power-dual.");
                close(solution.ports[1].pressure, pi, 1e-8,
                      "Mixed interface pressure feedback is wrong.");
                close(qi, -y * pi, 1e-8, "Passive interface flow continuity has the wrong sign.");
                const auto cavityPressure = -qi / volumeY;
                const auto slitResult =
                    std::find_if(solution.elements.begin(), solution.elements.end(),
                                 [](const auto &e) { return e.id.value == 5; });
                const auto cavityNode = std::find_if(solution.nodes.begin(), solution.nodes.end(),
                                                     [](const auto &n) { return n.id.value == 6; });
                require(slitResult != solution.elements.end() && cavityNode != solution.nodes.end(),
                        "Series slit or separated cavity pressure result is absent.");
                close(slitResult->flow, -qi, 1e-8, "BEM interface and series slit flows differ.");
                close(cavityNode->effort, cavityPressure, 1e-8,
                      "Slit pressure drop was omitted before cavity pressure.");
                close(pi - cavityPressure, slitZ * (-qi), 1e-8,
                      "Coupled slit does not obey its Stokes impedance.");
                close(modalMean(actual, false), pd, 1e-8,
                      "Actual modal pressure differs from its MNA feedback.");
                close(modalMean(actual, true), pi, 1e-8,
                      "Actual interface pressure differs from its MNA feedback.");
                std::map<std::uint32_t, Complex> velocities;
                Complex flowDriven{}, flowInterface{}, boundaryPower{};
                for (const auto &v : actual.solvedVelocities)
                    velocities.emplace(v.faceId, v.normalVelocity);
                close(velocities.at(1), velocity, 1e-8,
                      "Dome velocity is not the mechanical modal velocity.");
                close(velocities.at(4), .5 * velocity, 1e-8,
                      "Whole surround is not half-speed in phase.");
                for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
                    const auto face = mesh.triangles[t].faceId;
                    const auto found = velocities.find(face);
                    if (found == velocities.end())
                        continue;
                    const auto area = integrationAreas[t];
                    if (face == 1 || face == 4)
                        flowDriven += area * found->second;
                    if (face == 5)
                        flowInterface += area * found->second;
                    if (face != 2)
                        boundaryPower +=
                            area * actual.boundary()->surfacePressure[t] * std::conj(found->second);
                }
                close(flowDriven, qd, 1e-8, "Integrated driven modal flux is inconsistent.");
                close(flowInterface, qi, 1e-8, "Integrated interface flux is inconsistent.");
                close(boundaryPower, pd * std::conj(qd) + pi * std::conj(qi), 1e-8,
                      "Modal and actual boundary complex powers differ.");
                const auto &driver = solution.drivers[0];
                const double losses = driver.copperPower + driver.mechanicalPower +
                                      (lineZ.real() + slitZ.real()) * std::norm(qi);
                require(std::abs(driver.electricalPower - losses - boundaryPower.real()) < 1e-9 &&
                            actual.boundaryResidual < 1e-10 && solution.relativeResidual < 1e-10,
                        "Electrical-to-acoustic power or original block equations do not close.");
                if (!known)
                    require(boundaryPower.real() > 0, "Passive mixed-port BEM creates energy.");
                for (std::size_t o = 0; o < input.observers.size(); ++o)
                    close(actual.observerPressure[o],
                          columns[0].observerPressure[o] * qd +
                              columns[1].observerPressure[o] * qi +
                              (known ? columns[2].observerPressure[o] : Complex{}),
                          1e-8,
                          "Voltage-driven modal field differs from independent physical basis "
                          "solves.");
                // Dome, whole surround, complete interface and known wall source
                // have no azimuthal variation.
                // Full 3D quadrature has finite local integration errors and its
                // canonical touching rules depend on DOF numbering. m=0 binds
                // the rotational space directly and retains the strict check.
                const double symmetryTolerance =
                    axisymmetric ? 1e-8
                                 : std::min(1e-6, .01 * input.numerics.nearRelativeTolerance);
                double symmetryDifference{};
                for (std::size_t group = 3; group < input.observers.size(); group += 4)
                    for (std::size_t direction = 1; direction < 4; ++direction) {
                        const double scale = std::max(std::abs(actual.observerPressure[0]),
                                                      std::abs(actual.observerPressure[group]));
                        const double difference =
                            std::abs(actual.observerPressure[group] -
                                     actual.observerPressure[group + direction]);
                        symmetryDifference = std::max(symmetryDifference, difference / scale);
                        if (difference > symmetryTolerance * scale)
                            std::cerr << "modal symmetry m0=" << axisymmetric << " known=" << known
                                      << " f=" << input.frequencies[f] << " group=" << group
                                      << " direction=" << direction
                                      << " qOrder=" << input.numerics.quadratureOrder
                                      << " sOrder=" << input.numerics.singularOrder
                                      << " nearTolerance=" << input.numerics.nearRelativeTolerance
                                      << " reference=" << actual.observerPressure[group]
                                      << " other=" << actual.observerPressure[group + direction]
                                      << " normalized difference=" << difference / scale << '\n';
                        require(difference <= symmetryTolerance * scale,
                                "Axisymmetric modal motor has unequal H/V or signed-angle fields.");
                    }
                std::cout << "modal motor m0=" << axisymmetric << " known=" << known
                          << " f=" << input.frequencies[f] << " velocity=" << velocity
                          << " max H/V normalized difference=" << symmetryDifference
                          << " residual=" << actual.boundaryResidual << '\n';
            }
        }
        auto invalid = input;
        invalid.excitations = {{"Overlap", {{1, {.01, 0}}}}};
        bool rejected{};
        try {
            (void)solver.solve(invalid);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Known excitation overlaps a weighted Driven network port.");
    }
}

AxisymmetricMesh rotationalCylinder(unsigned sectors, unsigned radial, double depth) {
    AxisymmetricMeridian meridian;
    meridian.axis = {0, 0, 1};
    meridian.radial = {1, 0, 0};
    for (unsigned r = 0; r <= radial; ++r) {
        meridian.points.push_back({.025 * r / radial, -depth});
        if (r)
            meridian.segments.push_back({r - 1, r, 1, BoundaryCondition::driven});
    }
    meridian.points.push_back({.025, 0});
    meridian.segments.push_back({radial, radial + 1, 2, BoundaryCondition::wall});
    for (unsigned r = radial; r-- > 0;) {
        const auto next = static_cast<unsigned>(meridian.points.size());
        meridian.points.push_back({.025 * r / radial, 0});
        meridian.segments.push_back({next - 1, next, 3, BoundaryCondition::wall});
    }
    return buildAxisymmetricMesh(meridian, sectors);
}

AxisymmetricBem reduction(const AxisymmetricMesh &mesh, unsigned sectors) {
    AxisymmetricBem result;
    result.sectors = sectors;
    result.pressureDofs = mesh.pressureDofs;
    result.fluxDofs = mesh.fluxDofs;
    result.testPanels = mesh.testPanels;
    return result;
}

void planCacheLifecycle() {
    const auto geometry = rotationalCylinder(16, 2, .006);
    auto input = request();
    input.axisymmetric = reduction(geometry, 16);
    input.interactionPlanFrequencies = input.frequencies;
    CudaBemSolver coldSolver(geometry.mesh);
    const auto cold = coldSolver.solve(input);
    require(!cold.frequencies[0].interactionPlanReused &&
                cold.frequencies[0].interactionPlanMilliseconds > 0,
            "Cold baffled solve has no separately measured host planning time.");
    const auto compare = [](const BemFrequencyResult &a, const BemFrequencyResult &b) {
        require(a.frequency == b.frequency && a.excitations.size() == b.excitations.size(),
                "Cached solve changed frequency or excitation dimensions.");
        for (std::size_t e = 0; e < a.excitations.size(); ++e) {
            const auto &x = a.excitations[e], &y = b.excitations[e];
            const auto equal = [](const auto &one, const auto &two) {
                require(one.size() == two.size(), "Cached result field dimension changed.");
                for (std::size_t i = 0; i < one.size(); ++i)
                    require(std::abs(one[i] - two[i]) <= 1e-12 * std::max(1., std::abs(one[i])),
                            "Warm host-plan reuse changed the actual numerical result.");
            };
            equal(x.boundary()->surfacePressure, y.boundary()->surfacePressure);
            equal(x.boundary()->pressureTrace, y.boundary()->pressureTrace);
            equal(x.boundary()->normalDerivative, y.boundary()->normalDerivative);
            equal(x.observerPressure, y.observerPressure);
            require(std::abs(x.boundaryResidual - y.boundaryResidual) <= 1e-12,
                    "Warm plan reuse changed the original-system residual.");
        }
    };
    const auto warm = coldSolver.solve(input);
    require(warm.frequencies[0].interactionPlanReused &&
                warm.frequencies[0].interactionPlanMilliseconds == 0,
            "An unchanged baffled sweep rebuilt its host plan.");
    for (std::size_t f = 0; f < cold.frequencies.size(); ++f)
        compare(cold.frequencies[f], warm.frequencies[f]);

    CudaBemSolver scheduledSolver(geometry.mesh);
    for (std::size_t f = 0; f < input.frequencies.size(); ++f) {
        auto scheduled = input;
        scheduled.frequencies = {input.frequencies[f]};
        const auto one = scheduledSolver.solve(scheduled);
        require(one.frequencies[0].interactionPlanReused == (f != 0),
                "A scheduled frequency did not retain the explicit full-sweep plan.");
        compare(cold.frequencies[f], one.frequencies[0]);
    }
    auto coupled = input;
    coupled.frequencies = {500};
    coupled.excitations.clear();
    coupled.coupling =
        lem::makeDriverNetwork({6.4, .0007, 4.98, .004, .0016, .262, .00544}, {{1}}, {1, 0}, .0004);
    const auto motor = scheduledSolver.solve(coupled);
    require(motor.frequencies[0].interactionPlanReused,
            "Changing the LEM network rebuilt a geometry/quadrature-only plan.");
    std::get<lem::EffortSource>(coupled.coupling->elements.front().component).amplitude *= 2;
    const auto doubled = scheduledSolver.solve(coupled);
    require(doubled.frequencies[0].interactionPlanReused,
            "Changing a LEM source rebuilt the unchanged host plan.");
    for (std::size_t o = 0; o < input.observers.size(); ++o)
        require(
            std::abs(doubled.frequencies[0].excitations[0].observerPressure[o] -
                     2. * motor.frequencies[0].excitations[0].observerPressure[o]) <=
                1e-10 *
                    std::max(1., std::abs(motor.frequencies[0].excitations[0].observerPressure[o])),
            "Cached LEM source variation reused an old field rather than only its plan.");
    auto bad = input;
    bad.numerics.singularOrder = 21;
    bool rejected{};
    try {
        (void)scheduledSolver.solve(bad);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected && scheduledSolver.solve(input).frequencies[0].interactionPlanReused,
            "Failed plan replacement destroyed the preceding valid entry.");
    std::stop_source stop;
    stop.request_stop();
    rejected = false;
    try {
        (void)scheduledSolver.solve(input, stop.get_token());
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    require(rejected && scheduledSolver.solve(input).frequencies[0].interactionPlanReused,
            "Cancellation destroyed the preceding valid entry.");
    auto changed = input;
    changed.frequencies = {1000};
    changed.interactionPlanFrequencies.clear();
    const auto newFrequency = scheduledSolver.solve(changed);
    require(newFrequency.frequencies[0].interactionPlanReused,
            "A continuous geometric plan was unnecessarily tied to its old frequency.");
    const auto independentFrequency = CudaBemSolver(geometry.mesh).solve(changed);
    compare(independentFrequency.frequencies[0], newFrequency.frequencies[0]);
    changed = input;
    changed.numerics.singularOrder += 2;
    require(!scheduledSolver.solve(changed).frequencies[0].interactionPlanReused,
            "A different quadrature order reused the previous plan.");
    std::cout << "Baffled plan cache cold_ms=" << cold.frequencies[0].interactionPlanMilliseconds
              << " warm_ms=" << warm.frequencies[0].interactionPlanMilliseconds
              << "; actual fields/residuals, scheduled raster and LEM source variation passed.\n";
}

void axisymmetricSpace() {
    double firstDifference{}, lastDifference{};
    for (unsigned level = 0; level < 3; ++level) {
        const unsigned sectors = 16u << level, radial = 2u << level;
        const auto geometry = rotationalCylinder(sectors, radial, .006);
        auto input = request();
        input.frequencies = {500, 2000, 6000};
        input.numerics.nearDistanceRatio = 1.3;
        const double sine = std::sin(std::numbers::pi / 3), cosine = .5;
        for (auto direction : {Vec3{sine, 0, cosine}, Vec3{-sine, 0, cosine}, Vec3{0, sine, cosine},
                               Vec3{0, -sine, cosine}})
            input.observers.push_back(direction);
        CudaBemSolver solver(geometry.mesh);
        const auto full = solver.solve(input);
        input.axisymmetric = reduction(geometry, sectors);
        const auto reduced = solver.solve(input);
        require(reduced.geometryMode == BemGeometryMode::axisymmetric &&
                    reduced.pressureDofCount == 2 * (radial + 1) &&
                    reduced.triangleCount == geometry.mesh.triangles.size() &&
                    reduced.peakDeviceBytes < full.peakDeviceBytes,
                "Axisymmetric mode did not reduce the actual GPU operator space.");
        double difference{};
        for (std::size_t f = 0; f < reduced.frequencies.size(); ++f) {
            const auto &one = reduced.frequencies[f].excitations[0];
            const auto &complex = reduced.frequencies[f].excitations[1];
            checkTrace(geometry.mesh, one, input.medium, input.frequencies[f], 3,
                       input.axisymmetric);
            require(one.boundaryResidual < 1e-10 &&
                        one.observerPressure[2] == std::complex<double>{},
                    "Reduced original BM/Rayleigh residual or screening is wrong.");
            for (std::size_t observer = 0; observer < input.observers.size(); ++observer) {
                require(std::abs(complex.observerPressure[observer] -
                                 std::complex<double>{.3, -.4} * one.observerPressure[observer]) <=
                            1e-10 * std::max(1.0, std::abs(one.observerPressure[observer])),
                        "Reduced complex excitation loses its phase.");
                difference = std::max(
                    difference,
                    std::abs(one.observerPressure[observer] -
                             full.frequencies[f].excitations[0].observerPressure[observer]) /
                        std::abs(full.frequencies[f].excitations[0].observerPressure[0]));
            }
            for (unsigned observer = 4; observer < 7; ++observer)
                require(std::abs(one.observerPressure[observer] - one.observerPressure[3]) <=
                            1e-8 * std::abs(one.observerPressure[0]),
                        "Reduced solved field differs between H/V signed angles.");
        }
        if (level == 0)
            firstDifference = difference;
        lastDifference = difference;
        std::cout << "m0 vs 3D N=" << sectors << " radial=" << radial
                  << " pressure DOFs=" << reduced.pressureDofCount << '/' << full.pressureDofCount
                  << " field difference/axis=" << difference
                  << " GPU bytes=" << reduced.peakDeviceBytes << '/' << full.peakDeviceBytes
                  << " assembly ms=" << reduced.frequencies[0].assemblyMilliseconds << '/'
                  << full.frequencies[0].assemblyMilliseconds << '\n';
    }
    require(lastDifference < firstDifference && lastDifference < .02,
            "Meridian-strip space did not converge against independent full 3D BEM.");

    const auto shallow = rotationalCylinder(32, 4, .001);
    auto input = request();
    input.numerics.nearDistanceRatio = 1.3;
    input.axisymmetric = reduction(shallow, 32);
    const auto original = CudaBemSolver(shallow.mesh).solve(input);
    for (std::size_t f = 0; f < original.frequencies.size(); ++f) {
        const auto &one = original.frequencies[f].excitations[0];
        const double omega = 2 * std::numbers::pi * input.frequencies[f],
                     k = omega / input.medium.soundSpeed;
        const auto exact =
            std::complex<double>{0, input.medium.density * omega * one.impedances[0].area /
                                        (2 * std::numbers::pi)} *
            std::exp(std::complex<double>{0, -k});
        require(std::abs(one.observerPressure[0] / exact - 1.0) < .06,
                "Reduced shallow cavity misses the analytical baffled piston.");
    }
    auto placedMesh = shallow.mesh;
    const Vec3 translation{.125, -.25, .5};
    const auto placed = [&](Vec3 p) {
        return Vec3{p.z + translation.x, p.y + translation.y, -p.x + translation.z};
    };
    for (auto &p : placedMesh.vertices)
        p = placed(p);
    for (auto &p : input.observers)
        p = placed(p);
    input.axisymmetric->origin = translation;
    input.axisymmetric->axis = {1, 0, 0};
    input.infiniteBaffle->origin = translation;
    input.infiniteBaffle->normal = {1, 0, 0};
    input.farField->origin = translation;
    const auto moved = CudaBemSolver(placedMesh).solve(input);
    for (std::size_t f = 0; f < moved.frequencies.size(); ++f)
        require(std::abs(moved.frequencies[f].excitations[0].observerPressure[0] /
                             original.frequencies[f].excitations[0].observerPressure[0] -
                         1.0) < 1e-7,
                "Reduced frame changed its field after rigid rotation and translation.");

    auto expectRejected = [&](SurfaceMesh mesh, BemRequest bad) {
        bool rejected{};
        try {
            (void)CudaBemSolver(std::move(mesh)).solve(bad);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Invalid periodic geometry or non-axisymmetric mode was accepted.");
    };
    input = request();
    input.axisymmetric = reduction(shallow, 32);
    auto bad = input;
    bad.axisymmetric->pressureDofs[1] = bad.axisymmetric->pressureDofs[0];
    expectRejected(shallow.mesh, bad);
    bad = input;
    bad.axisymmetric->testPanels.push_back(bad.axisymmetric->testPanels.front());
    expectRejected(shallow.mesh, bad);
    auto asymmetric = shallow.mesh;
    asymmetric.vertices[1].x += 1e-7;
    expectRejected(asymmetric, input);
    asymmetric = shallow.mesh;
    asymmetric.triangles.front().faceId = 4;
    bad = input;
    bad.excitations = {{"Azimuthal velocity", {{1, {1, 0}}, {4, {.5, 0}}}}};
    expectRejected(asymmetric, bad);
}
void rejectInvalid() {
    const auto mesh = cylinder(12, 2, .005);
    auto expect = [&](BemRequest input) {
        bool rejected = false;
        try {
            (void)CudaBemSolver(mesh).solve(input);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Baffled unsupported request was accepted.");
    };
    auto input = request();
    input.infiniteBaffle->origin.z = -.001;
    expect(input);
    input = request();
    input.infiniteBaffle->normal = {0, 0, -1};
    expect(input);
    input = request();
    input.infiniteBaffle->apertureFaceIds = {1};
    expect(input);
    input = request();
    input.numerics.acceleration = BemAcceleration::fastMultipole;
    expect(input);
    input = request();
    input.farField.reset();
    expect(input);
    input = request();
    input.coupling = lem::Network{};
    expect(input);
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 2 && std::string_view(argv[1]) == "--plan-cache-only") {
            planCacheLifecycle();
            return 0;
        }
        require(argc == 1, "Unknown baffled test argument.");
        rejectInvalid();
        const auto coarse = run(16, 2, .006), medium = run(24, 3, .006), fine = run(32, 4, .006);
        const double first = difference(coarse, medium), second = difference(medium, fine);
        std::cout << "Independent fixed-depth mesh differences " << first << " -> " << second
                  << '\n';
        require(second < first && second < .03,
                "Baffled fixed-depth mesh sequence did not converge.");
        const auto shallow = run(24, 3, .001);
        require(shallow.pistonError < .06 && shallow.pistonError < fine.pistonError,
                "Shallow cavity did not converge toward the analytical baffled piston.");
        // At 50 micrometres the opposite, geometrically disjoint disk panels
        // invoke the projected near rule. This exercises its actual CUDA map,
        // not just its CPU planner, against the equal-area piston limit.
        const auto veryShallow = run(16, 2, .00005, 8);
        require(veryShallow.pistonError < .01 && veryShallow.pistonError < shallow.pistonError,
                "Projected GPU near integration misses the very shallow piston limit.");
        finiteField();
        interfaceCoupling();
        modalDriverCoupling();
        axisymmetricSpace();
        rigidPlacement();
        std::cout << "Baffled GPU interior BM / Rayleigh tests passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
