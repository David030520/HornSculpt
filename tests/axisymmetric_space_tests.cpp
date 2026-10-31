#include "axisymmetric_reduction.hpp"
#include "galerkin/axisymmetric_space.hpp"
#include <hornsculpt/axisymmetric_mesh.hpp>
#include <hornsculpt/surface_measures.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
void close(double value, double expected, const char *message) {
    require(std::isfinite(value) && std::isfinite(expected) &&
                std::abs(value - expected) <= 2e-12 * std::max(std::abs(expected), 1e-8),
            message);
}
void closeComplex(std::complex<double> value, std::complex<double> expected, const char *message) {
    close(value.real(), expected.real(), message);
    close(value.imag(), expected.imag(), message);
}
template <class F> void rejects(F function) {
    bool failed{};
    try {
        function();
    } catch (const std::invalid_argument &) {
        failed = true;
    }
    require(failed, "Invalid rotational surface contract was accepted.");
}
AxisymmetricMeridian frustum(double r0, double r1, double height) {
    AxisymmetricMeridian m;
    m.points = {{0, 0}, {r0, 0}, {r1, height}, {0, height}};
    m.segments = {{0, 1, 10, BoundaryCondition::driven},
                  {1, 2, 20, BoundaryCondition::wall},
                  {2, 3, 30, BoundaryCondition::wall}};
    return m;
}
AxisymmetricBem reduction(const AxisymmetricMesh &mesh, const AxisymmetricMeridian &m,
                          unsigned sectors) {
    AxisymmetricBem result;
    result.origin = m.origin;
    result.axis = m.axis;
    result.sectors = sectors;
    result.pressureDofs = mesh.pressureDofs;
    result.fluxDofs = mesh.fluxDofs;
    result.testPanels = mesh.testPanels;
    result.surfaceContract = AxisymmetricSurfaceContract::continuousRevolution;
    return result;
}
void checkMass(const galerkin::AxisymmetricElement &e, double r0, double r1, double length) {
    const double factor = 2 * std::numbers::pi * length;
    const std::array<double, 4> expected{factor * (3 * r0 + r1) / 12, factor * (r0 + r1) / 12,
                                         factor * (r0 + r1) / 12, factor * (r0 + 3 * r1) / 12};
    const std::array<double, 2> mixed{factor * (2 * r0 + r1) / 6, factor * (r0 + 2 * r1) / 6};
    close(e.length, length, "Meridian arclength changed.");
    close(e.area, std::numbers::pi * (r0 + r1) * length, "Exact revolved area is incorrect.");
    for (unsigned i = 0; i < 4; ++i)
        close(e.massP1P1[i], expected[i], "Exact P1/P1 moment is incorrect.");
    for (unsigned i = 0; i < 2; ++i) {
        close(e.massP1DP0[i], mixed[i], "Exact P1/DP0 moment is incorrect.");
        close(e.massP1P1[2 * i] + e.massP1P1[2 * i + 1], mixed[i],
              "Consistent mass row sum differs from mixed mass.");
    }
    close(e.massP1DP0[0] + e.massP1DP0[1], e.area, "Mixed masses do not conserve surface area.");
}
void checkCellMoments(const AxisymmetricMesh &mesh, const AxisymmetricBem &r,
                      const galerkin::AxisymmetricSpace &space, const SurfaceMeasures &measures) {
    require(measures.areas.size() == mesh.mesh.triangles.size() &&
                measures.pressureIntegrals.size() == measures.areas.size(),
            "Surface measure dimensions differ from the full mesh.");
    std::vector<double> area(space.elements.size());
    std::vector<std::array<double, 2>> mixed(space.elements.size());
    std::vector<std::complex<double>> pressure(space.elements.size());
    for (std::size_t t = 0; t < mesh.mesh.triangles.size(); ++t) {
        const auto strip = r.fluxDofs[t];
        const auto &e = space.elements[strip];
        area[strip] += measures.areas[t];
        close(std::accumulate(measures.pressureIntegrals[t].begin(),
                              measures.pressureIntegrals[t].end(), 0.),
              measures.areas[t], "A continuous half-cell does not conserve its area.");
        for (unsigned local = 0; local < 3; ++local) {
            const auto dof = r.pressureDofs[mesh.mesh.triangles[t].vertices[local]];
            const double weight = measures.pressureIntegrals[t][local];
            require(weight > 0, "Continuous nodal integral is not positive.");
            mixed[strip][dof == e.first ? 0 : 1] += weight;
            pressure[strip] += weight * std::complex<double>{1. + dof, .3 * dof - .2};
        }
    }
    for (std::size_t strip = 0; strip < space.elements.size(); ++strip) {
        const auto &e = space.elements[strip];
        close(area[strip], e.area, "Full ring-cell sum differs from physical revolved area.");
        for (unsigned endpoint = 0; endpoint < 2; ++endpoint)
            close(mixed[strip][endpoint], e.massP1DP0[endpoint],
                  "Full-cell nodal integrals differ from radial P1 masses.");
        const std::complex<double> p0{1. + e.first, .3 * e.first - .2},
            p1{1. + e.second, .3 * e.second - .2};
        closeComplex(pressure[strip], e.massP1DP0[0] * p0 + e.massP1DP0[1] * p1,
                     "Pressure integration lost power-conjugate radial weighting.");
    }
}
std::array<double, 3> quadratureHalfCell(double r0, double r1, double length, unsigned sectors,
                                         bool second) {
    // Independent tensor Gauss integration of the parameter-triangle basis.
    const double offset = std::sqrt(3. / 5) / 2;
    const std::array<double, 3> nodes{.5 - offset, .5, .5 + offset},
        weights{5. / 18, 4. / 9, 5. / 18};
    std::array<double, 3> result{};
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            const double s = nodes[i], t = second ? 1 - s + s * nodes[j] : (1 - s) * nodes[j];
            const double jacobian = (second ? s : 1 - s) * length * (r0 + (r1 - r0) * s) * 2 *
                                    std::numbers::pi / sectors;
            const std::array<double, 3> shape = second
                                                    ? std::array<double, 3>{1 - s, s + t - 1, 1 - t}
                                                    : std::array<double, 3>{1 - s - t, t, s};
            for (unsigned k = 0; k < 3; ++k)
                result[k] += weights[i] * weights[j] * jacobian * shape[k];
        }
    return result;
}
void exactCylinderAndFrustum() {
    for (const auto dimensions :
         {std::array<double, 3>{.03, .03, .02}, std::array<double, 3>{.02, .04, .03}}) {
        const auto m = frustum(dimensions[0], dimensions[1], dimensions[2]);
        const double sideLength = std::hypot(dimensions[1] - dimensions[0], dimensions[2]);
        for (unsigned sectors : {8u, 32u}) {
            const auto mesh = buildAxisymmetricMesh(m, sectors);
            const auto r = reduction(mesh, m, sectors);
            const auto space = galerkin::buildAxisymmetricSpace(mesh.mesh, {}, r);
            const auto measures = surfaceMeasures(mesh.mesh, r);
            require(space.elements.size() == 3 && space.meridian.points.size() == 4,
                    "Meridian extraction changed the DOF layout.");
            const auto &bottom = space.elements[0], &side = space.elements[1],
                       &top = space.elements[2];
            require(bottom.first == 0 && bottom.second == 1 && bottom.faceId == 10 &&
                        bottom.boundary == BoundaryCondition::driven && side.first == 1 &&
                        side.second == 2 && side.faceId == 20 && top.faceId == 30,
                    "Meridian orientation, face identities or boundary conditions changed.");
            checkMass(bottom, 0, dimensions[0], dimensions[0]);
            checkMass(side, dimensions[0], dimensions[1], sideLength);
            checkMass(top, dimensions[1], 0, dimensions[1]);
            close(bottom.radialNormal, 0, "Bottom normal has a radial component.");
            close(bottom.axialNormal, -1, "Bottom air normal does not point outwards.");
            close(top.axialNormal, 1, "Top air normal does not point outwards.");
            close(side.radialNormal, dimensions[2] / sideLength,
                  "Revolved side radial normal is wrong.");
            close(side.axialNormal, (dimensions[0] - dimensions[1]) / sideLength,
                  "Revolved side axial normal is wrong.");
            checkCellMoments(mesh, r, space, measures);
            close(std::accumulate(space.lumpedMass.begin(), space.lumpedMass.end(), 0.),
                  bottom.area + side.area + top.area,
                  "Lumped pressure masses do not conserve total area.");
            unsigned half = 0;
            for (auto t : mesh.testPanels)
                if (r.fluxDofs[t] == 1) {
                    const auto exact = quadratureHalfCell(dimensions[0], dimensions[1], sideLength,
                                                          sectors, half++ != 0);
                    for (unsigned i = 0; i < 3; ++i)
                        close(measures.pressureIntegrals[t][i], exact[i],
                              "Half-cell weights differ from independent parameter Gauss "
                              "integration.");
                }
            require(half == 2, "Frustum strip lacks a complete parameter-cell diagonal.");
            for (auto t : mesh.testPanels)
                if (r.fluxDofs[t] != 1)
                    for (double weight : measures.pressureIntegrals[t])
                        close(weight, measures.areas[t] / 3,
                              "Collapsed polar fan did not use its exact sector limit.");
            if (dimensions[0] != dimensions[1]) {
                const double mean = (2 * side.massP1DP0[0] + 7 * side.massP1DP0[1]) / side.area;
                require(std::abs(mean - 4.5) > .1,
                        "A tapered revolved pressure incorrectly uses endpoint arithmetic means.");
            }
            auto old = r;
            old.surfaceContract = AxisymmetricSurfaceContract::triangulatedOrbit;
            const auto legacy = surfaceMeasures(mesh.mesh, old),
                       cartesian = surfaceMeasures(mesh.mesh);
            require(legacy.areas == cartesian.areas &&
                        legacy.pressureIntegrals == cartesian.pressureIntegrals,
                    "Old polygonal traces were silently reinterpreted as continuous surfaces.");
            require(std::accumulate(legacy.areas.begin(), legacy.areas.end(), 0.) <
                        bottom.area + side.area + top.area,
                    "Finite polygonal facets unexpectedly have exact continuous revolution area.");
            rejects([&] { (void)galerkin::buildAxisymmetricSpace(mesh.mesh, {}, old); });
        }
    }
}
void conePoleAndPlacement() {
    AxisymmetricMeridian cone;
    cone.points = {{0, 0}, {.03, 0}, {0, .04}};
    cone.segments = {{0, 1, 11, BoundaryCondition::wall}, {1, 2, 12, BoundaryCondition::driven}};
    const auto coneMesh = buildAxisymmetricMesh(cone, 16);
    const auto coneReduction = reduction(coneMesh, cone, 16);
    const auto coneSpace = galerkin::buildAxisymmetricSpace(coneMesh.mesh, {}, coneReduction);
    checkMass(coneSpace.elements[1], .03, 0, .05);
    close(coneSpace.elements[1].radialNormal, .8, "Cone radial normal is incorrect.");
    close(coneSpace.elements[1].axialNormal, .6, "Cone axial normal is incorrect.");
    checkCellMoments(coneMesh, coneReduction, coneSpace,
                     surfaceMeasures(coneMesh.mesh, coneReduction));

    auto m = frustum(.02, .04, .03);
    const auto original = buildAxisymmetricMesh(m, 16);
    const auto originalReduction = reduction(original, m, 16);
    const auto originalMeasures = surfaceMeasures(original.mesh, originalReduction);
    m.origin = {.2, -.3, .17};
    m.axis = {1 / std::sqrt(14.), 2 / std::sqrt(14.), 3 / std::sqrt(14.)};
    m.radial = {2 / std::sqrt(5.), -1 / std::sqrt(5.), 0};
    auto moved = buildAxisymmetricMesh(m, 16);
    const auto movedReduction = reduction(moved, m, 16);
    const Vec3 shift{.13, -.29, .09};
    for (auto &p : moved.mesh.vertices) {
        p.x -= shift.x;
        p.y -= shift.y;
        p.z -= shift.z;
    }
    const auto source = galerkin::buildSpace(moved.mesh);
    const auto space = galerkin::buildAxisymmetricSpace(source, shift, movedReduction);
    const auto measures = surfaceMeasures(moved.mesh, movedReduction, shift);
    for (std::size_t i = 0; i < m.points.size(); ++i) {
        close(space.meridian.points[i].radius, m.points[i].radius,
              "Rigid placement changed a meridian radius.");
        require(std::abs(space.meridian.points[i].axial - m.points[i].axial) < 2e-14,
                "Rigid placement changed axial meridian coordinates.");
    }
    for (std::size_t t = 0; t < measures.areas.size(); ++t) {
        close(measures.areas[t], originalMeasures.areas[t],
              "Rigid placement changed continuous cell area.");
        for (unsigned i = 0; i < 3; ++i)
            close(measures.pressureIntegrals[t][i], originalMeasures.pressureIntegrals[t][i],
                  "Rigid placement changed nodal surface moments.");
    }
    checkCellMoments(moved, movedReduction, space, measures);
}
void invalidContracts() {
    const auto m = frustum(.02, .04, .03);
    auto mesh = buildAxisymmetricMesh(m, 16);
    const auto r = reduction(mesh, m, 16);
    auto bad = r;
    bad.pressureDofs[1] = bad.pressureDofs[0];
    rejects([&] { (void)surfaceMeasures(mesh.mesh, bad); });
    bad = r;
    bad.fluxDofs[0] = 1;
    rejects([&] { (void)galerkin::buildAxisymmetricSpace(mesh.mesh, {}, bad); });
    bad = r;
    bad.testPanels.push_back(bad.testPanels.front());
    rejects([&] { (void)surfaceMeasures(mesh.mesh, bad); });
    bad = r;
    bad.surfaceContract = static_cast<AxisymmetricSurfaceContract>(99);
    rejects([&] { (void)surfaceMeasures(mesh.mesh, bad); });
    auto changed = mesh.mesh;
    changed.triangles[0].faceId = 999;
    rejects([&] { (void)surfaceMeasures(changed, r); });
    changed = mesh.mesh;
    for (auto &t : changed.triangles)
        std::swap(t.vertices[0], t.vertices[1]);
    rejects([&] { (void)galerkin::buildAxisymmetricSpace(changed, {}, r); });
    rejects([&] {
        (void)surfaceMeasures(mesh.mesh, r, {std::numeric_limits<double>::quiet_NaN(), 0, 0});
    });
}
void expandedTraceAliases() {
    const auto meridian = frustum(.015, .025, .03);
    const auto mesh = buildAxisymmetricMesh(meridian, 16);
    const auto r = reduction(mesh, meridian, 16);
    std::vector<std::complex<double>> pressure, derivative;
    for (const auto id : r.pressureDofs)
        pressure.emplace_back(.1 * id, -.2 * id);
    for (const auto id : r.fluxDofs)
        derivative.emplace_back(.03 * id, .07 * id);
    bem_detail::validateAxisymmetricTrace(r, pressure, derivative);
    auto changed = pressure;
    changed[1] += 1.;
    rejects([&] { bem_detail::validateAxisymmetricTrace(r, changed, derivative); });
    auto changedDerivative = derivative;
    changedDerivative[0] += std::complex<double>{0, .001};
    rejects([&] { bem_detail::validateAxisymmetricTrace(r, pressure, changedDerivative); });
    changed.pop_back();
    rejects([&] { bem_detail::validateAxisymmetricTrace(r, changed, derivative); });
}
} // namespace
int main() {
    try {
        exactCylinderAndFrustum();
        conePoleAndPlacement();
        invalidContracts();
        expandedTraceAliases();
        std::cout << "Continuous rotational area, P1 masses, half-cell moments, pole limits, "
                     "placement and archive contracts passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
