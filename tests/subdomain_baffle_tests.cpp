#include "fixture_geometry.hpp"
#include "subdomain_baffle_fixture.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <hornsculpt/bem.hpp>
#include <iostream>
#include <map>
#include <stdexcept>
using namespace hornsculpt;
namespace {
struct Fixture {
    SurfaceMesh mesh;
    BemSubdomainModel model;
};
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
Fixture fixture(bool sector, bool mirrored = false) {
    Fixture f;
    const auto inner = testing::sphere(1);
    const auto outer = testing::sphere(1);
    const auto append = [&](const SurfaceMesh &sphere, double scale, std::uint32_t face,
                            BoundaryCondition boundary, double center) {
        std::map<std::uint32_t, std::uint32_t> indices;
        for (const auto &source : sphere.triangles) {
            if (sector && std::any_of(source.vertices.begin(), source.vertices.end(), [&](auto v) {
                    return sphere.vertices[v].x < -1e-12 || sphere.vertices[v].y < -1e-12;
                }))
                continue;
            auto triangle = source;
            for (auto &v : triangle.vertices) {
                const auto [it, added] =
                    indices.emplace(v, static_cast<std::uint32_t>(f.mesh.vertices.size()));
                if (added) {
                    const auto q = sphere.vertices[v];
                    f.mesh.vertices.push_back({scale * q.x, scale * q.y, scale * q.z + center});
                }
                v = it->second;
            }
            triangle.faceId = face;
            triangle.boundary = boundary;
            f.mesh.triangles.push_back(triangle);
        }
    };
    BemSubdomain exterior{2, "outside", BemDomainKind::exterior, {}};
    for (int copy = 0; copy < (mirrored ? 2 : 1); ++copy) {
        const double center = copy ? -3.0 : 3.0;
        const auto begin = static_cast<std::uint32_t>(f.mesh.triangles.size());
        append(inner, 1, 1, BoundaryCondition::driven, center);
        const auto middle = static_cast<std::uint32_t>(f.mesh.triangles.size());
        append(outer, 2, 2, BoundaryCondition::wall, center);
        BemSubdomain shell{
            static_cast<std::uint32_t>(1 + copy * 3), "shell", BemDomainKind::interior, {}};
        BemInterface aperture{
            static_cast<std::uint32_t>(3 + copy * 3), "interface", shell.id, 2, {}};
        for (auto t = begin; t < middle; ++t)
            shell.panels.push_back({t, -1});
        for (auto t = middle; t < f.mesh.triangles.size(); ++t) {
            shell.panels.push_back({t, 1});
            exterior.panels.push_back({t, -1});
            aperture.triangleIndices.push_back(t);
        }
        f.model.domains.push_back(shell);
        f.model.interfaces.push_back(aperture);
    }
    f.model.domains.push_back(exterior);
    return f;
}
BemRequest request(const Fixture &f) {
    BemRequest r;
    r.subdomains = f.model;
    r.frequencies = {300};
    r.excitations = {{"source", {{1, {1, 0}}}}};
    r.farField = FarFieldObservation{{0, 0, 0}, 1};
    r.observers = {{0, 0, 1}, {1, 0, 0}, {-1, 0, 0}, {0, -1, 0}};
    return r;
}
void compare(std::complex<double> actual, std::complex<double> expected, double limit,
             const char *message) {
    if (std::abs(actual - expected) > limit * std::max(1e-15, std::abs(expected)))
        throw std::runtime_error(message);
}
Fixture openRimFixture(bool sector) {
    Fixture f;
    BemSubdomain inside{1, "inside", BemDomainKind::interior, {}};
    BemSubdomain outside{2, "outside", BemDomainKind::exterior, {}};
    BemInterface interface{3, "spatial interface", 1, 2, {}};
    f.mesh = testing::openRimSourceMesh(sector);
    for (std::size_t index = 0; index < f.mesh.triangles.size(); ++index) {
        const auto face = f.mesh.triangles[index].faceId;
        const auto triangle = static_cast<std::uint32_t>(index);
        inside.panels.push_back({triangle, -1});
        if (face == 2) {
            outside.panels.push_back({triangle, 1});
            interface.triangleIndices.push_back(triangle);
        }
    }
    f.model.domains = {inside, outside};
    f.model.interfaces = {interface};
    return f;
}
void runOpenRim() {
    const auto full = openRimFixture(false), quarter = openRimFixture(true);
    const OrthogonalSymmetryPlanes planes{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 1e-11};
    const ExteriorBaffle baffle{{0, 0, 0}, {0, 0, 1}, 1e-11};
    auto fullRequest = request(full);
    fullRequest.frequencies = {180, 550, 1200};
    fullRequest.exteriorBaffle = baffle;
    fullRequest.observers = {{0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {1, 0, 0}, {0, 0, -1}};
    fullRequest.observerValidity = {1, 1, 1, 1, 0};
    CudaBemSolver fullSolver(full.mesh, full.model, baffle);
    auto invalidCoupling = fullRequest;
    invalidCoupling.coupling = lem::Network{};
    bool rejected = false;
    try {
        (void)fullSolver.solve(invalidCoupling);
    } catch (const std::invalid_argument &error) {
        rejected =
            std::string(error.what()).find("LEM coupling with symmetry or an exterior baffle") !=
            std::string::npos;
    }
    require(rejected, "Exterior baffle with LEM must be explicitly rejected");
    auto invalidMask = fullRequest;
    invalidMask.observerValidity.back() = 1;
    rejected = false;
    try {
        (void)fullSolver.solve(invalidMask);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "An unmasked rear-hemisphere direction must be rejected");
    invalidMask = fullRequest;
    invalidMask.observerValidity.pop_back();
    rejected = false;
    try {
        (void)fullSolver.solve(invalidMask);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "A wrong-sized observer mask must be rejected");
    const auto reference = fullSolver.solve(fullRequest);
    auto reducedRequest = request(quarter);
    reducedRequest.frequencies = fullRequest.frequencies;
    reducedRequest.observers = fullRequest.observers;
    reducedRequest.observerValidity = fullRequest.observerValidity;
    reducedRequest.symmetry = planes;
    reducedRequest.exteriorBaffle = baffle;
    CudaBemSolver reducedSolver(quarter.mesh, quarter.model, planes, baffle);
    const auto reduced = reducedSolver.solve(reducedRequest);
    for (std::size_t f = 0; f < reference.frequencies.size(); ++f) {
        const auto &a = reference.frequencies[f].excitations[0];
        const auto &b = reduced.frequencies[f].excitations[0];
        require(a.observerValidity == fullRequest.observerValidity &&
                    b.observerValidity == fullRequest.observerValidity,
                "Observer validity was not preserved");
        for (std::size_t i = 0; i < a.observerPressure.size(); ++i) {
            if (!a.observerValidity[i])
                continue;
            compare(b.observerPressure[i], a.observerPressure[i], 0.001,
                    "Open-rim baffle far field differs");
        }
        compare(b.impedances[0].meanPressure, a.impedances[0].meanPressure, 0.001,
                "Open-rim baffle pressure differs");
        compare(b.impedances[0].volumeVelocity, a.impedances[0].volumeVelocity, 1e-10,
                "Open-rim source volume differs");
        compare(b.impedances[0].volumeImpedance, a.impedances[0].volumeImpedance, 0.001,
                "Open-rim impedance differs");
    }
}
void run() {
    const auto full = fixture(false, true), quarter = fixture(true);
    OrthogonalSymmetryPlanes planes{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 1e-11};
    auto fullRequest = request(full);
    CudaBemSolver fullSolver(full.mesh, full.model);
    const auto fullResult = fullSolver.solve(fullRequest);
    auto reducedRequest = request(quarter);
    reducedRequest.symmetry = planes;
    ExteriorBaffle baffle{{0, 0, 0}, {0, 0, 1}, 1e-11};
    reducedRequest.exteriorBaffle = baffle;
    CudaBemSolver reducedSolver(quarter.mesh, quarter.model, planes, baffle);
    const auto reduced = reducedSolver.solve(reducedRequest);
    require(reduced.pressureDofCount < fullResult.pressureDofCount, "No pressure DOF reduction");
    const auto &a = fullResult.frequencies[0].excitations[0],
               &b = reduced.frequencies[0].excitations[0];
    require(a.observerPressure.size() == b.observerPressure.size() && b.impedances.size() == 1,
            "Incomplete symmetry result");
    for (std::size_t i = 0; i < a.observerPressure.size(); ++i)
        compare(b.observerPressure[i], a.observerPressure[i], .015,
                "Symmetric image field differs from full mesh");
    compare(b.impedances[0].meanPressure, a.impedances[0].meanPressure, .015,
            "Symmetry pressure differs");
    compare(2.0 * b.impedances[0].volumeVelocity, a.impedances[0].volumeVelocity, 1e-10,
            "Baffle volume-flow scaling differs");
    compare(b.impedances[0].volumeImpedance, 2.0 * a.impedances[0].volumeImpedance, .015,
            "Baffle impedance differs");
}
} // namespace
int main() {
    try {
        run();
        runOpenRim();
        std::cout << "Subdomain baffle PASS\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
