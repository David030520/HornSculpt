#include "acoustic_partition_fixture.hpp"
#include "fixture_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <hornsculpt/acoustic_partition.hpp>
#include <hornsculpt/bem.hpp>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

using namespace hornsculpt;
namespace {
using Complex = std::complex<double>;
constexpr double quadratureTolerance = .01;
constexpr double energyTolerance = .05;
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
double relative(Complex a, Complex b) {
    return std::abs(a - b) / std::max(std::abs(b), 1e-30);
}
const BemExcitationResult &excitation(const BemResult &r, std::size_t f = 0) {
    require(f < r.frequencies.size() && r.frequencies[f].excitations.size() == 1,
            "Missing frequency result");
    return r.frequencies[f].excitations[0];
}
Vec3 sub(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double area(const SurfaceMesh &m, const Triangle &t) {
    auto v = cross(sub(m.vertices[t.vertices[1]], m.vertices[t.vertices[0]]),
                   sub(m.vertices[t.vertices[2]], m.vertices[t.vertices[0]]));
    return std::hypot(v.x, v.y, v.z) / 2;
}
const BemSubdomainTrace &domain(const BemExcitationResult &r, std::uint32_t id) {
    auto it = std::find_if(r.subdomainTraces().begin(), r.subdomainTraces().end(),
                           [id](const auto &x) { return x.domainId == id; });
    require(it != r.subdomainTraces().end(), "Missing domain trace");
    return *it;
}
Complex interfaceFlow(const SurfaceMesh &mesh, const BemSubdomainTrace &trace,
                      const BemInterface &iface) {
    require(trace.triangleIndices.size() == trace.normalDerivative.size(),
            "Missing triangle flux trace");
    Complex value{};
    for (auto id : iface.triangleIndices) {
        auto it = std::find(trace.triangleIndices.begin(), trace.triangleIndices.end(), id);
        require(it != trace.triangleIndices.end(), "Interface panel missing from domain");
        value +=
            area(mesh, mesh.triangles[id]) *
            trace.normalDerivative[static_cast<std::size_t>(it - trace.triangleIndices.begin())];
    }
    return value;
}
void cupQuadrature() {
    MeshingOptions opts;
    opts.targetEdgeLength = .02;
    auto partition = partitionAcousticSurface(
        testing::acousticCup(), std::vector<AcousticRegionDefinition>{{7, "Horn air", {2, 3}, {}}},
        opts);
    require(partition.mesh.triangles.size() < 1000 && partition.model.interfaces.size() == 1,
            "Unexpected cup fixture");
    BemRequest request;
    request.subdomains = partition.model;
    request.frequencies = {500};
    request.excitations = {{"drive", {{3, {1, 0}}}}};
    request.observers = {{0, 0, .2}, {.2, 0, .2}, {0, 0, .04}};
    request.numerics.quadratureOrder = 4;
    request.numerics.singularOrder = 4;
    CudaBemSolver solver(partition.mesh, partition.model);
    const auto coarse = solver.solve(request);
    request.numerics.singularOrder = 6;
    const auto fine = solver.solve(request);
    request.numerics.nearRelativeTolerance = 1e-6;
    const auto regularTight = solver.solve(request);
    request.numerics.quadratureOrder = 4;
    request.numerics.nearRelativeTolerance = 1e-4;
    const auto &a = excitation(coarse), &b = excitation(fine);
    for (std::size_t j = 0; j < request.observers.size(); ++j) {
        const auto error =
            relative(excitation(regularTight).observerPressure.at(j), b.observerPressure.at(j));
        std::cout << "Cup regular/near tolerance 1e-4/1e-6, observer " << j
                  << ", relative complex difference " << error << '\n';
        require(error < quadratureTolerance,
                "Qualified/adaptive quadrature tolerance comparison exceeds 1%");
    }
    require(a.observerPressure.size() == 3 && b.observerPressure.size() == 3, "Missing cup fields");
    for (std::size_t j = 0; j < 3; ++j) {
        auto error = relative(a.observerPressure[j], b.observerPressure[j]);
        std::cout << "Cup singular quadrature 4/6, observer " << j
                  << ", relative complex difference " << error << '\n';
        require(error < quadratureTolerance, "Cup quadrature 4/6 pressure difference exceeds 1%");
    }
    request.observers = {{0, 0, 1}, {1, 0, 0}};
    request.farField = FarFieldObservation{{0, 0, 0}, 1};
    request.numerics.singularOrder = 4;
    const auto coarseFar = solver.solve(request);
    request.numerics.singularOrder = 6;
    const auto fineFar = solver.solve(request);
    for (std::size_t j = 0; j < 2; ++j) {
        const auto error = relative(excitation(coarseFar).observerPressure[j],
                                    excitation(fineFar).observerPressure[j]);
        std::cout << "Cup far field quadrature 4/6, direction " << j
                  << ", relative complex difference " << error << '\n';
        require(error < quadratureTolerance, "Cup quadrature 4/6 far field difference exceeds 1%");
    }
    const auto &iface = partition.model.interfaces.front();
    auto q1 = interfaceFlow(partition.mesh, domain(b, iface.firstDomainId), iface);
    auto q2 = interfaceFlow(partition.mesh, domain(b, iface.secondDomainId), iface);
    const auto mismatch = std::abs(q1 + q2) / std::max({std::abs(q1), std::abs(q2), 1e-30});
    std::cout << "Cup interface integrated derivative mismatch " << mismatch << '\n';
    require(std::abs(q1) > 1e-12 && mismatch < 1e-8,
            "Integrated interface flux is not equal and opposite");
}
struct SphereFixture {
    SurfaceMesh mesh;
    BemSubdomainModel model;
};
SphereFixture sphereFixture() {
    SphereFixture f;
    f.mesh = testing::sphere(3);
    auto outer = testing::sphere(3);
    const auto nv = static_cast<std::uint32_t>(f.mesh.vertices.size());
    const auto nt = static_cast<std::uint32_t>(f.mesh.triangles.size());
    for (auto p : outer.vertices)
        f.mesh.vertices.push_back({1.5 * p.x, 1.5 * p.y, 1.5 * p.z});
    for (auto t : outer.triangles) {
        for (auto &v : t.vertices)
            v += nv;
        t.faceId = 2;
        t.boundary = BoundaryCondition::wall;
        f.mesh.triangles.push_back(t);
    }
    BemSubdomain shell{1, "air shell", BemDomainKind::interior, {}};
    BemSubdomain outside{2, "outside", BemDomainKind::exterior, {}};
    BemInterface interface{3, "interface", 1, 2, {}};
    for (std::uint32_t t = 0; t < nt; ++t)
        shell.panels.push_back({t, -1});
    for (std::uint32_t t = nt; t < f.mesh.triangles.size(); ++t) {
        shell.panels.push_back({t, 1});
        outside.panels.push_back({t, -1});
        interface.triangleIndices.push_back(t);
    }
    f.model.domains = {shell, outside};
    f.model.interfaces = {interface};
    return f;
}
void sphereEnergy() {
    const auto f = sphereFixture();
    require(f.mesh.triangles.size() <= 1024, "Sphere fixture too large");
    BemRequest req;
    req.subdomains = f.model;
    // Frequencies correspond to ka approximately 1, pi and 2pi at a=0.1 m.
    req.frequencies = {343. / (2 * std::numbers::pi * .1), 1715, 3430};
    req.excitations = {{"drive", {{1, {1, 0}}}}};
    req.farField = FarFieldObservation{{0, 0, 0}, 1};
    constexpr unsigned directions = 128;
    constexpr double golden = 2.39996322972865332;
    for (unsigned j = 0; j < directions; ++j) {
        const double z = 1 - 2 * (j + .5) / directions, radial = std::sqrt(1 - z * z),
                     phi = golden * j;
        req.observers.push_back({radial * std::cos(phi), radial * std::sin(phi), z});
    }
    CudaBemSolver split(f.mesh, f.model);
    const auto splitResult = split.solve(req);
    auto singleReq = req;
    singleReq.subdomains.reset();
    CudaBemSolver single(testing::sphere(3));
    const auto singleResult = single.solve(singleReq);
    for (std::size_t i = 0; i < req.frequencies.size(); ++i) {
        const auto &a = excitation(splitResult, i), &b = excitation(singleResult, i);
        require(a.impedances.size() == 1 && a.observerPressure.size() == directions &&
                    b.observerPressure.size() == directions,
                "Missing sphere output");
        double power = 0, worst = 0;
        for (unsigned j = 0; j < directions; ++j) {
            power += std::norm(a.observerPressure[j]) * 4 * std::numbers::pi / directions /
                     (req.medium.density * req.medium.soundSpeed);
            worst = std::max(worst, relative(a.observerPressure[j], b.observerPressure[j]));
        }
        const auto &port = a.impedances.front();
        const double input = std::real(port.meanPressure * std::conj(port.volumeVelocity));
        require(input > 0 && std::isfinite(input), "Active source power is not positive");
        const double imbalance = std::abs(power - input) / input;
        std::cout << "Sphere f=" << req.frequencies[i] << " Hz, far/input power " << power << '/'
                  << input << ", relative energy imbalance " << imbalance << ", split difference "
                  << worst << '\n';
        require(imbalance < energyTolerance, "Active power and far-field power differ by over 5%");
        require(worst < .05, "Split/unsplit sphere difference exceeds 5%");
    }
}
} // namespace
int main() {
    try {
        std::cout << "CUDA device: " << CudaBemSolver::deviceDescription() << '\n';
        cupQuadrature();
        sphereEnergy();
        std::cout << "Subdomain accuracy tests PASS.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Subdomain accuracy test failed: " << e.what() << '\n';
        return 1;
    }
}
