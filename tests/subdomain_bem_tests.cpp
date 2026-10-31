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
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
struct Fixture {
    SurfaceMesh mesh;
    BemSubdomainModel model;
};
Fixture fixture(unsigned level) {
    Fixture f;
    f.mesh = testing::sphere(level);
    auto outer = testing::sphere(level);
    const auto nv = static_cast<std::uint32_t>(f.mesh.vertices.size());
    const auto nt = static_cast<std::uint32_t>(f.mesh.triangles.size());
    for (const auto v : outer.vertices)
        f.mesh.vertices.push_back({2 * v.x, 2 * v.y, 2 * v.z});
    for (auto t : outer.triangles) {
        for (auto &v : t.vertices)
            v += nv;
        t.faceId = 2;
        t.boundary = BoundaryCondition::wall;
        f.mesh.triangles.push_back(t);
    }
    BemSubdomain shell{1, "shell", BemDomainKind::interior, {}};
    BemSubdomain exterior{2, "exterior", BemDomainKind::exterior, {}};
    BemInterface iface{3, "interface", 1, 2, {}};
    for (std::uint32_t i = 0; i < nt; ++i)
        shell.panels.push_back({i, -1});
    for (std::uint32_t i = nt; i < f.mesh.triangles.size(); ++i) {
        shell.panels.push_back({i, 1});
        exterior.panels.push_back({i, -1});
        iface.triangleIndices.push_back(i);
    }
    f.model.domains = {shell, exterior};
    f.model.interfaces = {iface};
    return f;
}
BemRequest request(const BemSubdomainModel &model) {
    BemRequest r;
    r.subdomains = model;
    r.frequencies = {500};
    r.excitations = {{"unit", {{1, {1, 0}}}}};
    r.observers = {{0, 0, .3}, {1, 0, 0}};
    r.numerics.relativeResidualTolerance = 5e-5;
    return r;
}
Complex analytic(double radius) {
    const Complex i{0, 1};
    const double k = 2 * std::numbers::pi * 500 / 343.;
    return 1.21 * 343. * (i * k * .1) / (1. + i * k * .1) * (.1 / radius) *
           std::exp(-i * k * (radius - .1));
}
const BemSubdomainTrace &trace(const BemExcitationResult &r, std::uint32_t id) {
    auto it = std::find_if(r.subdomainTraces().begin(), r.subdomainTraces().end(),
                           [id](const auto &x) { return x.domainId == id; });
    require(it != r.subdomainTraces().end(), "Missing domain trace");
    return *it;
}
void checkInterface(const BemExcitationResult &r) {
    const auto &a = trace(r, 1), &b = trace(r, 2);
    require(a.interfaceFluxes.size() == 1 && b.interfaceFluxes.size() == 1,
            "Missing P1 interface flux");
    const auto &qa = a.interfaceFluxes.front(), &qb = b.interfaceFluxes.front();
    require(qa.interfaceId == 3 && qb.interfaceId == 3 && !qa.vertexIndices.empty() &&
                qa.vertexIndices.size() == qa.normalDerivative.size() &&
                qb.vertexIndices.size() == qb.normalDerivative.size() &&
                a.vertexIndices.size() == a.pressureTrace.size() &&
                b.vertexIndices.size() == b.pressureTrace.size(),
            "Invalid interface trace dimensions");
    double signal = 0;
    for (std::size_t j = 0; j < qa.vertexIndices.size(); ++j) {
        auto vertex = qa.vertexIndices[j];
        auto q = std::find(qb.vertexIndices.begin(), qb.vertexIndices.end(), vertex);
        auto pa = std::find(a.vertexIndices.begin(), a.vertexIndices.end(), vertex);
        auto pb = std::find(b.vertexIndices.begin(), b.vertexIndices.end(), vertex);
        require(q != qb.vertexIndices.end() && pa != a.vertexIndices.end() &&
                    pb != b.vertexIndices.end(),
                "Missing interface vertex");
        auto x = qa.normalDerivative[j],
             y = qb.normalDerivative[static_cast<std::size_t>(q - qb.vertexIndices.begin())];
        auto p = a.pressureTrace[static_cast<std::size_t>(pa - a.vertexIndices.begin())];
        auto other = b.pressureTrace[static_cast<std::size_t>(pb - b.vertexIndices.begin())];
        require(std::abs(x + y) <= 1e-8 * (1 + std::abs(x)), "Interface flux sign mismatch");
        require(std::abs(p - other) <= 1e-8 * (1 + std::abs(p)),
                "Interface pressure discontinuity");
        signal += std::norm(x);
    }
    require(signal > 1e-16, "Interface flux silently zero");
}
BemExcitationResult solve(const Fixture &f) {
    CudaBemSolver solver(f.mesh, f.model);
    auto out = solver.solve(request(f.model));
    require(out.frequencies.size() == 1 && out.frequencies.front().excitations.size() == 1,
            "Invalid result dimensions");
    auto r = out.frequencies.front().excitations.front();
    require(std::isfinite(r.boundaryResidual) && r.boundaryResidual <= 5e-5 &&
                r.observerPressure.size() == 2,
            "Residual or observer failure");
    checkInterface(r);
    return r;
}
void rejected(const Fixture &f, const char *message) {
    bool ok = false;
    try {
        CudaBemSolver solver(f.mesh, f.model);
        (void)solver.solve(request(f.model));
    } catch (const std::invalid_argument &) {
        ok = true;
    }
    require(ok, message);
}
void invalidModels() {
    const auto original = fixture(0);
    auto bad = original;
    bad.model.domains[0].panels[0].normalSign = 0;
    rejected(bad, "Accepted zero normal sign");
    bad = original;
    bad.model.domains[1].id = 1;
    rejected(bad, "Accepted duplicate domain ID");
    bad = original;
    bad.model.interfaces[0].secondDomainId = 999;
    rejected(bad, "Accepted unknown domain");
    bad = original;
    bad.model.domains[0].panels.pop_back();
    rejected(bad, "Accepted missing ownership");
    bad = original;
    bad.model.interfaces[0].triangleIndices[0] = 999999;
    rejected(bad, "Accepted unknown triangle");
    bad = original;
    bad.model.domains[1].panels[0].normalSign = 1;
    rejected(bad, "Accepted same interface normals");
    bad = original;
    auto third = bad.model.domains[1];
    third.id = 4;
    bad.model.domains.push_back(third);
    bad.model.interfaces.push_back({4, "overlap", 1, 4, bad.model.interfaces[0].triangleIndices});
    rejected(bad, "Accepted multi-domain interface junction");
}
void convergence() {
    std::vector<double> errors;
    for (unsigned level = 0; level <= 2; ++level) {
        auto f = fixture(level);
        auto s = solve(f);
        auto req = request(f.model);
        req.subdomains.reset();
        CudaBemSolver single(testing::sphere(level));
        auto out = single.solve(req);
        const auto &base = out.frequencies.front().excitations.front();
        double worst = 0;
        for (std::size_t i = 0; i < 2; ++i) {
            worst =
                std::max(worst, std::abs(s.observerPressure[i] / analytic(i == 0 ? .3 : 1.) - 1.));
            if (level == 2)
                require(std::abs(s.observerPressure[i] / base.observerPressure[i] - 1.) < .15,
                        "Split/unsplit mismatch");
        }
        errors.push_back(worst);
        std::cout << "Sphere level " << level << ", triangles " << f.mesh.triangles.size()
                  << ", error " << worst << ", residual " << s.boundaryResidual << '\n';
    }
    require(errors[1] < errors[0] && errors[2] < errors[1] && errors[2] < .15,
            "Sphere did not converge");
}
void orientationAndOrder() {
    auto f = fixture(1);
    const auto reference = solve(f);
    for (auto index : f.model.interfaces[0].triangleIndices) {
        std::swap(f.mesh.triangles[index].vertices[1], f.mesh.triangles[index].vertices[2]);
        for (auto &d : f.model.domains)
            for (auto &p : d.panels)
                if (p.triangleIndex == index)
                    p.normalSign = -p.normalSign;
    }
    std::swap(f.model.interfaces[0].firstDomainId, f.model.interfaces[0].secondDomainId);
    std::reverse(f.model.interfaces[0].triangleIndices.begin(),
                 f.model.interfaces[0].triangleIndices.end());
    for (auto &d : f.model.domains)
        std::reverse(d.panels.begin(), d.panels.end());
    std::reverse(f.model.domains.begin(), f.model.domains.end());
    auto changed = solve(f);
    for (std::size_t i = 0; i < reference.observerPressure.size(); ++i)
        require(std::abs(changed.observerPressure[i] / reference.observerPressure[i] - 1.) < 1e-5,
                "Orientation/order changed solution");
}
void cupSplitInvariance() {
    const auto physical = testing::acousticCup();
    MeshingOptions meshing;
    meshing.targetEdgeLength = .02;
    const auto divided = partitionAcousticSurface(
        physical, std::vector<AcousticRegionDefinition>{{7, "Hornluft", {2, 3}}}, meshing);
    BemRequest req;
    req.frequencies = {500};
    req.excitations = {{"Kolben", {{3, {1, 0}}}}};
    req.observers = {{0, 0, .2}, {.2, 0, .2}, {0, 0, .04}};
    CudaBemSolver original(physical);
    const auto baseline = original.solve(req);
    req.subdomains = divided.model;
    CudaBemSolver coupled(divided.mesh, divided.model);
    const auto output = coupled.solve(req);
    const auto &a = baseline.frequencies.front().excitations.front();
    const auto &b = output.frequencies.front().excitations.front();
    double worst = 0;
    for (std::size_t i = 0; i < req.observers.size(); ++i)
        worst = std::max(worst, std::abs(b.observerPressure[i] / a.observerPressure[i] - 1.));
    std::cout << "Cup interface-rim split comparison " << worst << ", residual "
              << b.boundaryResidual << '\n';
    require(worst < .05, "Cup with interface rim changed physical pressure by more than 5%");
    require(b.impedances.size() == 1 && b.impedances.front().faceId == 3,
            "Subdomain physical source impedance missing");
    auto segmented = physical;
    for (auto &t : segmented.triangles)
        if (t.faceId == 2) {
            double z = 0;
            for (auto vertex : t.vertices)
                z += segmented.vertices[vertex].z / 3;
            if (z > .04)
                t.faceId = 6;
        }
    const auto chain = partitionAcousticSurface(
        segmented,
        std::vector<AcousticRegionDefinition>{{7, "Hals", {2, 3}, {}}, {8, "Glocke", {6}, {}}},
        meshing);
    req.subdomains = chain.model;
    CudaBemSolver threeDomains(chain.mesh, chain.model);
    const auto threeResult = threeDomains.solve(req);
    const auto &threePressure =
        threeResult.frequencies.front().excitations.front().observerPressure;
    double chainDifference = 0;
    for (std::size_t i = 0; i < req.observers.size(); ++i)
        chainDifference =
            std::max(chainDifference, std::abs(threePressure[i] / a.observerPressure[i] - 1.));
    std::cout << "Three-domain horn split comparison " << chainDifference << '\n';
    require(chainDifference < .05,
            "Throat/bell/exterior decomposition changed physical pressure by more than 5%");
    const auto volume = std::abs(a.impedances.front().volumeVelocity);
    require(std::abs(b.impedances.front().volumeVelocity - a.impedances.front().volumeVelocity) <
                1e-10 * volume,
            "Subdomain decomposition changed physical driven volume flow");
}
} // namespace
int main() {
    try {
        std::cout << "CUDA device: " << CudaBemSolver::deviceDescription() << '\n';
        invalidModels();
        convergence();
        orientationAndOrder();
        cupSplitInvariance();
        std::cout << "Subdomain tests PASS.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Subdomain BEM test failed: " << e.what() << '\n';
        return 1;
    }
}
