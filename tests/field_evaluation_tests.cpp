#include "fixture_geometry.hpp"
#include "galerkin/space.hpp"
#include <cmath>
#include <hornsculpt/field_evaluation.hpp>
#include <iostream>
#include <numbers>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
using Complex = std::complex<double>;
Complex monopole(Vec3 p, double k) {
    const double r = std::hypot(p.x, p.y, p.z);
    return std::polar(1 / r, -k * r);
}
BemExcitationResult monopoleTrace(const SurfaceMesh &mesh, double k) {
    const auto space = galerkin::buildSpace(mesh);
    BemExcitationResult trace;
    BemBoundaryData boundary;
    for (const auto p : mesh.vertices)
        boundary.pressureTrace.push_back(monopole(p, k));
    for (const auto &t : space.panels) {
        const auto p = t.centroid;
        const double r = std::hypot(p.x, p.y, p.z);
        boundary.normalDerivative.push_back(
            monopole(p, k) * Complex{-1 / r, -k} *
            ((p.x * t.normal.x + p.y * t.normal.y + p.z * t.normal.z) / r));
    }
    trace.setBoundary(std::move(boundary));
    return trace;
}
void exteriorMonopole() {
    BemRequest request;
    const double frequency = 350, k = 2 * std::numbers::pi * frequency / request.medium.soundSpeed;
    const auto mesh = testing::sphere(4);
    const auto trace = monopoleTrace(mesh, k);
    const auto coarseMesh = testing::sphere(3);
    const auto coarse =
        evaluatePressureField(coarseMesh, request, frequency, monopoleTrace(coarseMesh, k),
                              {{0, 0, .3}, {.2, .15, .4}}, {1e-5, 1e-5, 8});
    const auto field = evaluatePressureField(
        mesh, request, frequency, trace, {{0, 0, .3}, {.2, .15, .4}, {0, 0, 0}, mesh.vertices[0]},
        {1e-5, 1e-5, 8});
    for (unsigned i = 0; i < 2; ++i) {
        require(field[i].status == FieldPointStatus::air, "Exterior air was masked.");
        const auto exact = monopole(field[i].position, k);
        const double relativeError = std::abs(field[i].pressure - exact) / std::abs(exact);
        require(relativeError < .003 && std::abs(field[i].pressure - exact) <
                                            .35 * std::abs(coarse[i].pressure - exact),
                "Exterior layer potential disagrees with the analytic monopole.");
    }
    require(field[2].status == FieldPointStatus::excluded,
            "Solid interior was shown as sound field.");
    require(field[3].status == FieldPointStatus::boundary, "Surface point was not masked.");
    auto incomplete = trace;
    auto incompleteData = *trace.boundary();
    incompleteData.normalDerivative.clear();
    incomplete.setBoundary(std::move(incompleteData));
    bool rejected = false;
    try {
        evaluatePressureField(mesh, request, frequency, incomplete, {{0, 0, .3}});
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "Missing saved derivative trace was silently reconstructed.");
    std::stop_source cancelled;
    cancelled.request_stop();
    bool stopped = false;
    try {
        evaluatePressureField(mesh, request, frequency, trace, {{0, 0, .3}}, {},
                              cancelled.get_token());
    } catch (const std::runtime_error &) {
        stopped = true;
    }
    require(stopped, "Field evaluation ignored cancellation.");
}
SurfaceMesh cylinder(unsigned sectors) {
    SurfaceMesh mesh;
    mesh.vertices = {{0, 0, 0}, {0, 0, -.05}};
    for (unsigned layer = 0; layer < 2; ++layer)
        for (unsigned s = 0; s < sectors; ++s) {
            const double a = 2 * std::numbers::pi * s / sectors;
            mesh.vertices.push_back({.1 * std::cos(a), .1 * std::sin(a), layer == 0 ? 0. : -.05});
        }
    for (unsigned s = 0; s < sectors; ++s) {
        const auto a = 2 + s, b = 2 + (s + 1) % sectors, c = a + sectors, d = b + sectors;
        mesh.triangles.push_back({{0, a, b}, 2});
        mesh.triangles.push_back({{1, d, c}, 1});
        mesh.triangles.push_back({{a, c, d}, 1});
        mesh.triangles.push_back({{a, d, b}, 1});
    }
    return mesh;
}
void apertureAndInterior() {
    const auto mesh = cylinder(128);
    const auto space = galerkin::buildSpace(mesh);
    BemRequest request;
    request.infiniteBaffle = InfiniteBaffle{{}, {0, 0, 1}, {2}};
    const double frequency = 700, k = 2 * std::numbers::pi * frequency / request.medium.soundSpeed;
    const Complex q{1.3, -.4};
    BemExcitationResult trace;
    BemBoundaryData boundary;
    boundary.pressureTrace.resize(mesh.vertices.size());
    for (const auto &t : mesh.triangles)
        boundary.normalDerivative.push_back(t.faceId == 2 ? q : Complex{});
    trace.setBoundary(boundary);
    const auto field = evaluatePressureField(
        mesh, request, frequency, trace, {{0, 0, .15}, {0, 0, 1}, {.2, 0, -.03}}, {1e-5, 1e-6, 8});
    for (unsigned i = 0; i < 2; ++i) {
        const double z = field[i].position.z;
        const auto exact = -q * (std::polar(1., -k * z) - std::polar(1., -k * std::hypot(z, .1))) *
                           Complex{0, -1 / k};
        require(field[i].status == FieldPointStatus::air &&
                    std::abs(field[i].pressure - exact) / std::abs(exact) < .0007,
                "Baffle aperture field disagrees with exact disk Rayleigh integral.");
    }
    require(field[2].status == FieldPointStatus::excluded,
            "Back half-space outside cavity was not masked.");
    // Independent plane-wave Cauchy data on the complete closed cavity. The
    // representation must reproduce the interior wave, including its sign.
    const double lowFrequency = 40,
                 lowK = 2 * std::numbers::pi * lowFrequency / request.medium.soundSpeed;
    boundary.pressureTrace.clear();
    boundary.normalDerivative.clear();
    for (const auto p : mesh.vertices)
        boundary.pressureTrace.push_back(std::polar(1., -lowK * p.z));
    for (const auto &t : space.panels)
        boundary.normalDerivative.push_back(Complex{0, -lowK * t.normal.z} *
                                         std::polar(1., -lowK * t.centroid.z));
    trace.setBoundary(boundary);
    const auto inside =
        evaluatePressureField(mesh, request, lowFrequency, trace, {{0, 0, -.025}}, {1e-5, 1e-6, 8});
    require(inside[0].status == FieldPointStatus::air &&
                std::abs(inside[0].pressure - std::polar(1., lowK * .025)) < .0004,
            "Interior cavity layer representation has incorrect pressure or orientation.");
    const auto unresolved =
        evaluatePressureField(mesh, request, frequency, trace, {{0, 0, .00002}}, {1e-6, 1e-12, 0});
    require(unresolved[0].status == FieldPointStatus::quadratureNotConverged,
            "Unresolved near field was accepted as calculated pressure.");
}
void continuousRevolution() {
    constexpr unsigned sectors = 32;
    const auto mesh = cylinder(sectors);
    AxisymmetricBem reduction;
    reduction.sectors = sectors;
    reduction.pressureDofs = {3, 0};
    for (unsigned layer = 0; layer < 2; ++layer)
        for (unsigned i = 0; i < sectors; ++i)
            reduction.pressureDofs.push_back(layer ? 1 : 2);
    for (unsigned i = 0; i < sectors; ++i) {
        reduction.fluxDofs.insert(reduction.fluxDofs.end(), {2, 0, 1, 1});
        if (i == 0)
            reduction.testPanels = {0, 1, 2, 3};
    }
    BemRequest request;
    request.axisymmetric = reduction;
    request.infiniteBaffle = InfiniteBaffle{{}, {0, 0, 1}, {2}};
    const double frequency = 700, k = 2 * std::numbers::pi * frequency / request.medium.soundSpeed;
    const Complex q{1.3, -.4};
    BemExcitationResult trace;
    BemBoundaryData boundary;
    boundary.pressureTrace.resize(mesh.vertices.size());
    for (const auto &t : mesh.triangles)
        boundary.normalDerivative.push_back(t.faceId == 2 ? q : Complex{});
    trace.setBoundary(boundary);
    const auto field = evaluatePressureField(
        mesh, request, frequency, trace,
        {{0, 0, .15}, {0, 0, 1}, {.2, 0, -.03}, {0, 0, -.05}, {.2, 0, 0}}, {1e-5, 1e-6, 8});
    for (unsigned i = 0; i < 2; ++i) {
        const double z = field[i].position.z;
        const auto exact = -q * (std::polar(1., -k * z) - std::polar(1., -k * std::hypot(z, .1))) *
                           Complex{0, -1 / k};
        require(field[i].status == FieldPointStatus::air &&
                    std::abs(field[i].pressure / exact - 1.) < 2e-6,
                "Continuous ring field misses the analytic disk Rayleigh integral.");
    }
    require(field[2].status == FieldPointStatus::excluded &&
                field[3].status == FieldPointStatus::boundary &&
                field[4].status == FieldPointStatus::boundary,
            "Continuous ring field misclassifies solid, physical pole or infinite baffle.");
    auto historical = request;
    historical.axisymmetric->surfaceContract = AxisymmetricSurfaceContract::triangulatedOrbit;
    const auto archived =
        evaluatePressureField(mesh, historical, frequency, trace, {{0, 0, 1}}, {1e-5, 1e-6, 8});
    require(archived[0].status == FieldPointStatus::air &&
                std::abs(archived[0].pressure / field[1].pressure - 1.) > 1e-4,
            "Historical polygonal traces were silently reinterpreted as continuous rings.");
    const double lowFrequency = 40,
                 lowK = 2 * std::numbers::pi * lowFrequency / request.medium.soundSpeed;
    boundary.pressureTrace.clear();
    boundary.normalDerivative.clear();
    for (const auto p : mesh.vertices)
        boundary.pressureTrace.push_back(std::polar(1., -lowK * p.z));
    for (const auto &t : mesh.triangles) {
        if (t.faceId == 2)
            boundary.normalDerivative.push_back(Complex{0, -lowK});
        else if (t.vertices[0] == 1)
            boundary.normalDerivative.push_back(Complex{0, lowK} * std::polar(1., lowK * .05));
        else
            boundary.normalDerivative.push_back({});
    }
    trace.setBoundary(boundary);
    const auto inside =
        evaluatePressureField(mesh, request, lowFrequency, trace,
                              {{0, 0, -.025}, {.03, 0, -.025}, {0, 0, -.03}}, {1e-5, 1e-6, 8});
    for (const auto &value : inside)
        require(value.status == FieldPointStatus::air &&
                    std::abs(value.pressure - std::polar(1., -lowK * value.position.z)) < .0004,
                "Continuous interior representation masks the axis or changes its physical sign.");
    auto corrupted = trace;
    auto corruptedData = *trace.boundary();
    corruptedData.pressureTrace[3] *= 1.01;
    corrupted.setBoundary(std::move(corruptedData));
    bool rejected = false;
    try {
        evaluatePressureField(mesh, request, lowFrequency, corrupted, {{0, 0, -.025}});
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "Continuous field accepted different pressures within one ring DOF.");
}
} // namespace
int main() {
    try {
        exteriorMonopole();
        apertureAndInterior();
        continuousRevolution();
        std::cout
            << "Analytic exterior, Rayleigh aperture, interior, masks and cancellation passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
