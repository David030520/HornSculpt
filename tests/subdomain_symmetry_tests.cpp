#include "fixture_geometry.hpp"
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
Fixture fixture(bool sector) {
    Fixture f;
    const auto inner = testing::sphere(1);
    const auto outer = testing::sphere(1);
    const auto append = [&](const SurfaceMesh &sphere, double scale, std::uint32_t face,
                            BoundaryCondition boundary) {
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
                    f.mesh.vertices.push_back({scale * q.x, scale * q.y, scale * q.z});
                }
                v = it->second;
            }
            triangle.faceId = face;
            triangle.boundary = boundary;
            f.mesh.triangles.push_back(triangle);
        }
    };
    append(inner, 1, 1, BoundaryCondition::driven);
    const auto middle = static_cast<std::uint32_t>(f.mesh.triangles.size());
    append(outer, 2, 2, BoundaryCondition::wall);
    BemSubdomain shell{1, "shell", BemDomainKind::interior, {}},
        exterior{2, "outside", BemDomainKind::exterior, {}};
    BemInterface aperture{3, "interface", 1, 2, {}};
    for (std::uint32_t t = 0; t < middle; ++t)
        shell.panels.push_back({t, -1});
    for (auto t = middle; t < f.mesh.triangles.size(); ++t) {
        shell.panels.push_back({t, 1});
        exterior.panels.push_back({t, -1});
        aperture.triangleIndices.push_back(t);
    }
    f.model.domains = {shell, exterior};
    f.model.interfaces = {aperture};
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
void run() {
    const auto full = fixture(false), quarter = fixture(true);
    OrthogonalSymmetryPlanes planes{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 1e-11};
    auto fullRequest = request(full);
    CudaBemSolver fullSolver(full.mesh, full.model);
    const auto fullResult = fullSolver.solve(fullRequest);
    auto reducedRequest = request(quarter);
    reducedRequest.symmetry = planes;
    CudaBemSolver reducedSolver(quarter.mesh, quarter.model, planes);
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
    compare(b.impedances[0].volumeVelocity, a.impedances[0].volumeVelocity, 1e-10,
            "Symmetry volume-flow scaling differs");
    compare(b.impedances[0].volumeImpedance, a.impedances[0].volumeImpedance, .015,
            "Symmetry impedance differs");
}
} // namespace
int main() {
    try {
        run();
        std::cout << "Subdomain symmetry PASS\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
