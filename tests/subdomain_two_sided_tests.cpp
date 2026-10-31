#include "fixture_geometry.hpp"
#include <algorithm>
#include <complex>
#include <hornsculpt/acoustic_dof_map.hpp>
#include <hornsculpt/bem.hpp>
#include <iostream>
#include <map>
#include <numbers>
#include <stdexcept>
using namespace hornsculpt;
namespace {
using Complex = std::complex<double>;
void require(bool x, const char *m) {
    if (!x)
        throw std::runtime_error(m);
}
struct Fixture {
    SurfaceMesh mesh;
    BemSubdomainModel model;
    std::uint32_t south{};
};
Fixture fixture(bool driven) {
    Fixture f;
    f.mesh = testing::sphere(1);
    auto outer = testing::sphere(1);
    const auto offset = static_cast<std::uint32_t>(f.mesh.vertices.size());
    for (std::uint32_t i = 0; i < outer.vertices.size(); ++i) {
        const auto v = outer.vertices[i];
        f.mesh.vertices.push_back({2 * v.x, 2 * v.y, 2 * v.z});
        if (v.z < -.099999)
            f.south = offset + i;
    }
    const auto inner = static_cast<std::uint32_t>(f.mesh.triangles.size());
    BemSubdomain shell{1, "inside", BemDomainKind::interior, {}},
        ext{2, "outside", BemDomainKind::exterior, {}};
    BemInterface interface{3, "opening", 1, 2, {}};
    BemTwoSidedSurface shared{1, 2, {}};
    for (std::uint32_t t = 0; t < inner; ++t)
        shell.panels.push_back({t, -1});
    for (auto triangle : outer.triangles) {
        double z = 0;
        for (auto &v : triangle.vertices) {
            z += outer.vertices[v].z;
            v += offset;
        }
        const auto index = static_cast<std::uint32_t>(f.mesh.triangles.size());
        const bool opening = z >= 0;
        triangle.faceId = opening ? 2 : 3;
        triangle.boundary =
            !opening && driven ? BoundaryCondition::driven : BoundaryCondition::wall;
        f.mesh.triangles.push_back(triangle);
        shell.panels.push_back({index, 1});
        ext.panels.push_back({index, -1});
        if (opening)
            interface.triangleIndices.push_back(index);
        else
            shared.triangleIndices.push_back(index);
    }
    f.model.domains = {shell, ext};
    f.model.interfaces = {interface};
    f.model.twoSidedSurfaces = {shared};
    return f;
}
const BemSubdomainTrace &trace(const BemExcitationResult &result, std::uint32_t id) {
    const auto it = std::find_if(result.subdomainTraces().begin(), result.subdomainTraces().end(),
                                 [&](const auto &v) { return v.domainId == id; });
    if (it == result.subdomainTraces().end())
        throw std::runtime_error("Missing domain trace");
    return *it;
}
Complex pressure(const BemSubdomainTrace &t, std::uint32_t v) {
    const auto it = std::find(t.vertexIndices.begin(), t.vertexIndices.end(), v);
    if (it == t.vertexIndices.end())
        throw std::runtime_error("Missing pressure vertex");
    return t.pressureTrace[static_cast<std::size_t>(it - t.vertexIndices.begin())];
}
Complex derivative(const BemSubdomainTrace &t, std::uint32_t triangle) {
    const auto it = std::find(t.triangleIndices.begin(), t.triangleIndices.end(), triangle);
    if (it == t.triangleIndices.end())
        throw std::runtime_error("Missing flux triangle");
    return t.normalDerivative[static_cast<std::size_t>(it - t.triangleIndices.begin())];
}
void run(bool driven) {
    const auto f = fixture(driven);
    const auto dofs = makeAcousticPressureDofMap(f.mesh, f.model);
    require(dofs.domains[0].at(f.south) != dofs.domains[1].at(f.south),
            "Wall-side pressure DOFs were welded");
    const auto vertex = f.mesh.triangles[f.model.interfaces[0].triangleIndices.front()].vertices[0];
    require(dofs.domains[0].at(vertex) == dofs.domains[1].at(vertex),
            "Interface pressure DOFs were not welded");
    BemRequest request;
    request.subdomains = f.model;
    request.frequencies = {300};
    request.excitations = {{"drive", {{driven ? 3u : 1u, {1, 0}}}}};
    CudaBemSolver solver(f.mesh, f.model);
    const auto result = solver.solve(request);
    const auto &excitation = result.frequencies[0].excitations[0];
    const auto &a = trace(excitation, 1), &b = trace(excitation, 2);
    require(std::abs(pressure(a, f.south) - pressure(b, f.south)) > 1e-4,
            "Two physical sides have an artificial common pressure");
    for (auto t : f.model.interfaces[0].triangleIndices)
        for (auto v : f.mesh.triangles[t].vertices)
            require(std::abs(pressure(a, v) - pressure(b, v)) <
                        1e-9 * (1 + std::abs(pressure(a, v))),
                    "Transparent interface pressure is discontinuous");
    if (driven) {
        const auto t = f.model.twoSidedSurfaces[0].triangleIndices.front();
        const Complex expected{0, 1.21 * 2 * std::numbers::pi * 300};
        require(std::abs(derivative(a, t) - expected) < 1e-8 * std::abs(expected),
                "First membrane side has wrong flux");
        require(std::abs(derivative(b, t) + expected) < 1e-8 * std::abs(expected),
                "Second membrane side has wrong flux");
        require(excitation.impedances.size() == 1 && excitation.impedances[0].faceId == 3,
                "Two-sided membrane impedance missing");
    }
}
} // namespace
int main() {
    try {
        run(false);
        run(true);
        std::cout << "Two-sided subdomain tests PASS\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
