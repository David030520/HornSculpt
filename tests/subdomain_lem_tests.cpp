#include "fixture_geometry.hpp"
#include <cmath>
#include <complex>
#include <hornsculpt/bem.hpp>
#include <hornsculpt/surface_measures.hpp>
#include <iostream>
#include <map>
#include <stdexcept>
using namespace hornsculpt;
namespace {
using Complex = std::complex<double>;
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
void close(Complex a, Complex b, double tolerance, const char *message) {
    if (std::abs(a - b) > tolerance * std::max({1e-12, std::abs(a), std::abs(b)}))
        throw std::runtime_error(message);
}
struct Fixture {
    SurfaceMesh mesh;
    BemSubdomainModel model;
};
Fixture makeFixture() {
    Fixture f;
    f.mesh = testing::sphere(1);
    auto interface = testing::sphere(1);
    auto second = testing::sphere(1);
    const auto innerCount = static_cast<std::uint32_t>(f.mesh.triangles.size());
    const auto add = [&](SurfaceMesh part, double scale, Vec3 shift, std::uint32_t face,
                         BoundaryCondition bc) {
        const auto offset = static_cast<std::uint32_t>(f.mesh.vertices.size());
        for (auto p : part.vertices)
            f.mesh.vertices.push_back(
                {scale * p.x + shift.x, scale * p.y + shift.y, scale * p.z + shift.z});
        for (auto t : part.triangles) {
            for (auto &v : t.vertices)
                v += offset;
            t.faceId = face;
            t.boundary = bc;
            f.mesh.triangles.push_back(t);
        }
    };
    add(interface, 2, {0, 0, 0}, 2, BoundaryCondition::wall);
    const auto interfaceEnd = static_cast<std::uint32_t>(f.mesh.triangles.size());
    add(second, 1, {.6, 0, 0}, 3, BoundaryCondition::driven);
    BemSubdomain shell{1, "shell", BemDomainKind::interior, {}},
        exterior{2, "outside", BemDomainKind::exterior, {}};
    BemInterface opening{3, "opening", 1, 2, {}};
    for (std::uint32_t t = 0; t < innerCount; ++t)
        shell.panels.push_back({t, -1});
    for (auto t = innerCount; t < interfaceEnd; ++t) {
        shell.panels.push_back({t, 1});
        exterior.panels.push_back({t, -1});
        opening.triangleIndices.push_back(t);
    }
    for (auto t = interfaceEnd; t < f.mesh.triangles.size(); ++t)
        exterior.panels.push_back({t, -1});
    f.model.domains = {shell, exterior};
    f.model.interfaces = {opening};
    return f;
}
lem::Network network() {
    lem::Network n;
    n.nodes = {{{1}, lem::Domain::acoustic, true, "reference", {}},
               {{2}, lem::Domain::acoustic, false, "front", {}},
               {{3}, lem::Domain::acoustic, false, "rear", {}}};
    n.elements.push_back({{1}, "front flow", {}, lem::FlowSource{{{1}, {2}}, {.001, 0}}});
    n.elements.push_back({{2}, "front port", {}, lem::BemPort{{{2}, {1}}, {1}, 1}});
    n.elements.push_back({{3}, "rear flow", {}, lem::FlowSource{{{1}, {3}}, {.0004, .0002}}});
    n.elements.push_back({{4}, "rear port", {}, lem::BemPort{{{3}, {1}}, {3}, 1}});
    return n;
}
std::vector<Complex> panelPressure(const Fixture &f, const BemExcitationResult &e) {
    std::vector<Complex> p(f.mesh.triangles.size());
    for (const auto &trace : e.subdomainTraces()) {
        std::map<std::uint32_t, Complex> nodes;
        for (std::size_t i = 0; i < trace.vertexIndices.size(); ++i)
            nodes.emplace(trace.vertexIndices[i], trace.pressureTrace[i]);
        for (const auto t : trace.triangleIndices) {
            if (f.mesh.triangles[t].faceId == 2)
                continue;
            for (auto v : f.mesh.triangles[t].vertices)
                p[t] += nodes.at(v) / 3.;
        }
    }
    return p;
}
void run() {
    const auto f = makeFixture();
    const auto measures = surfaceMeasures(f.mesh);
    double areas[2]{};
    for (std::size_t t = 0; t < f.mesh.triangles.size(); ++t) {
        if (f.mesh.triangles[t].faceId == 1)
            areas[0] += measures.areas[t];
        if (f.mesh.triangles[t].faceId == 3)
            areas[1] += measures.areas[t];
    }
    CudaBemSolver solver(f.mesh, f.model);
    BemRequest reference;
    reference.subdomains = f.model;
    reference.frequencies = {300};
    reference.excitations = {{"basis one", {{1, {1 / areas[0], 0}}}},
                             {"basis two", {{3, {1 / areas[1], 0}}}}};
    const auto basis = solver.solve(reference);
    lem::AcousticLoad load;
    load.impedance.resize(4);
    for (std::size_t j = 0; j < 2; ++j) {
        const auto p = panelPressure(f, basis.frequencies[0].excitations[j]);
        for (std::size_t t = 0; t < p.size(); ++t) {
            const auto face = f.mesh.triangles[t].faceId;
            if (face == 1 || face == 3) {
                const auto row = face == 1 ? 0u : 1u;
                load.impedance[row * 2 + j] += p[t] * measures.areas[t] / areas[row];
            }
        }
    }
    require(std::abs(load.impedance[1]) > 1e-9 && std::abs(load.impedance[2]) > 1e-9,
            "Missing mutual loading");
    const auto n = network();
    const auto expected = lem::solve(n, 300, load, {1.21, 343});
    BemRequest coupled;
    coupled.subdomains = f.model;
    coupled.frequencies = {300};
    coupled.coupling = n;
    const auto result = solver.solve(coupled);
    const auto &actual = result.frequencies[0].excitations[0];
    require(actual.name == "LEM" && actual.lemSolution && actual.lemSolution->ports.size() == 2,
            "Missing coupled result");
    const auto copied = actual;
    require(!actual.subdomainTraces().empty() &&
                &actual.subdomainTraces() == &copied.subdomainTraces(),
            "Copies do not share immutable subdomain traces");
    auto detached = copied;
    detached.setSubdomainTraces({});
    require(detached.subdomainTraces().empty() && !actual.subdomainTraces().empty(),
            "Replacing traces modified an existing result");
    for (unsigned i = 0; i < 2; ++i) {
        close(actual.lemSolution->ports[i].pressure, expected.ports[i].pressure, 1e-6,
              "Wrong loaded port pressure");
        close(actual.lemSolution->ports[i].volumeVelocity, expected.ports[i].volumeVelocity, 1e-6,
              "Wrong loaded flow");
    }
    require(actual.boundaryResidual <= coupled.numerics.relativeResidualTolerance,
            "Final residual not accepted");
    // One electromechanical diaphragm injects and withdraws flow into two
    // different BEM domains; pressure feedback must change its motion.
    lem::Network driver;
    driver.nodes = {{{1}, lem::Domain::electrical, true, "electric reference", {}},
                    {{2}, lem::Domain::electrical, false, "drive", {}},
                    {{3}, lem::Domain::acoustic, true, "acoustic reference", {}},
                    {{4}, lem::Domain::acoustic, false, "front", {}},
                    {{5}, lem::Domain::acoustic, false, "rear", {}}};
    driver.elements.push_back({{1}, "voltage", {}, lem::EffortSource{{{2}, {1}}, {2, 0}}});
    driver.elements.push_back(
        {{2},
         "driver",
         {},
         lem::Driver{{{2}, {1}}, {4}, {5}, {6, .0002, 4, .003, .0005, 1, .002}}});
    driver.elements.push_back({{3}, "front", {}, lem::BemPort{{{4}, {3}}, {1}, 1}});
    driver.elements.push_back({{4}, "rear", {}, lem::BemPort{{{5}, {3}}, {3}, 1}});
    const auto expectedDriver = lem::solve(driver, 300, load, {1.21, 343});
    lem::AcousticLoad noLoad;
    noLoad.impedance.resize(4);
    const auto unloaded = lem::solve(driver, 300, noLoad, {1.21, 343});
    require(std::abs(expectedDriver.drivers[0].velocity - unloaded.drivers[0].velocity) > 1e-6,
            "Two-sided acoustic loading did not change diaphragm motion");
    coupled.coupling = driver;
    const auto driverResult = solver.solve(coupled).frequencies[0].excitations[0];
    require(driverResult.lemSolution && driverResult.lemSolution->drivers.size() == 1,
            "Missing diaphragm result");
    close(driverResult.lemSolution->drivers[0].velocity, expectedDriver.drivers[0].velocity, 1e-6,
          "Incorrect pressure feedback on diaphragm");
    close(driverResult.lemSolution->ports[0].volumeVelocity,
          -driverResult.lemSolution->ports[1].volumeVelocity, 1e-6,
          "Diaphragm flows are not opposite");
    auto invalid = n;
    std::get<lem::BemPort>(invalid.elements[1].component).faceIds = {1, 3};
    coupled.coupling = invalid;
    bool rejected = false;
    try {
        (void)solver.solve(coupled);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "Cross-domain port was accepted");
}
} // namespace
int main() {
    try {
        run();
        std::cout << "Subdomain LEM tests PASS\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
