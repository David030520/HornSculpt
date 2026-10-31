#include "baffle_exterior.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <hornsculpt/bem.hpp>
#include <hornsculpt/field_evaluation.hpp>
#include <hornsculpt/surface_measures.hpp>
#include <iostream>
#include <numbers>
#include <stdexcept>
using namespace hornsculpt;
using Complex = std::complex<double>;
constexpr double pi = std::numbers::pi;
void require(bool b, const char *m) {
    if (!b)
        throw std::runtime_error(m);
}
void close(Complex a, Complex b, double tol, const char *m) {
    if (std::abs(a - b) > tol * std::max(std::abs(b), 1e-20)) {
        std::cerr << m << " actual=" << a << " expected=" << b << '\n';
        throw std::runtime_error(m);
    }
}
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

SurfaceMesh disk(unsigned sectors, unsigned radial) {
    const auto closed = cylinder(sectors, radial, .02);
    SurfaceMesh m;
    const unsigned layer = 1 + sectors * radial;
    m.vertices.assign(closed.vertices.begin() + layer, closed.vertices.end());
    for (auto &v : m.vertices)
        v.x += .07;
    for (auto t : closed.triangles)
        if (t.faceId == 3) {
            for (auto &n : t.vertices)
                n -= layer;
            t.faceId = 4;
            t.boundary = BoundaryCondition::driven;
            m.triangles.push_back(t);
        }
    return m;
}
Complex pistonImpedance(double frequency, double area, AcousticMedium medium) {
    const double a = std::sqrt(area / pi), ka = 2 * pi * frequency * a / medium.soundSpeed,
                 z = 2 * ka;
    // H1(z), independently integrated after t=sin(theta).
    double sum = 0;
    constexpr int n = 1024;
    const double step = pi / (2 * n);
    for (int j = 0; j <= n; ++j) {
        double t = j * step, v = std::cos(t) * std::cos(t) * std::sin(z * std::sin(t));
        sum += (j == 0 || j == n ? 1 : j % 2 ? 4 : 2) * v;
    }
    const double h1 = 2 * z / pi * sum * step / 3;
    return medium.density * medium.soundSpeed / area *
           Complex{1 - std::cyl_bessel_j(1, z) / ka, h1 / ka};
}
Complex mean(const BemExcitationResult &e, unsigned face) {
    for (const auto &z : e.impedances)
        if (z.faceId == face)
            return z.meanPressure;
    throw std::runtime_error("Expected source pressure missing.");
}
int main() {
    try {
        auto mesh = cylinder(16, 2, .02);
        auto exterior = disk(16, 2);
        const auto areas = surfaceMeasures(exterior).areas;
        double area = 0;
        for (auto a : areas)
            area += a;
        BemRequest input;
        input.medium = {1.205, 343.32};
        input.frequencies = {500, 2000};
        input.infiniteBaffle = InfiniteBaffle{{}, {0, 0, 1}, {3}, exterior};
        input.numerics.singularOrder = 6;
        input.numerics.nearRelativeTolerance = 1e-4;
        input.farField = FarFieldObservation{};
        input.observers = {{0, 0, 1}, {.7, 0, 1}};
        input.excitations = {{"Rear Q", {{1, {1 / area, 0}}, {4, {0, 0}}}},
                             {"Front Q", {{1, {0, 0}}, {4, {1 / area, 0}}}}};
        bem_detail::validateBaffleExteriorMesh(mesh, exterior, {}, {0, 0, 1});
        const auto extSpace = galerkin::buildSpace(exterior);
        for (auto f : input.frequencies) {
            const auto op = bem_detail::baffleExteriorOperator(
                galerkin::buildSpace(mesh), {}, extSpace, input.numerics,
                2 * pi * f / input.medium.soundSpeed);
            Complex integral{};
            for (std::size_t i = 0; i < areas.size(); ++i)
                for (std::size_t j = 0; j < areas.size(); ++j) {
                    integral += areas[i] * op.self[i * areas.size() + j];
                    close(areas[i] * op.self[i * areas.size() + j],
                          areas[j] * op.self[j * areas.size() + i], 1e-12,
                          "Exterior self operator violates reciprocity.");
                }
            auto z = Complex{0, 4 * pi * f * input.medium.density} * integral / (area * area);
            close(z, pistonImpedance(f, area, input.medium), .025,
                  "Exterior disk misses analytical piston load.");
            std::cout << "piston " << f << " Hz Z=" << z
                      << " reference=" << pistonImpedance(f, area, input.medium) << '\n';
        }
        CudaBemSolver solver(mesh);
        const auto basis = solver.solve(input);
        std::vector<lem::AcousticLoad> loads;
        for (const auto &f : basis.frequencies) {
            const auto &a = f.excitations[0], &b = f.excitations[1];
            // a unit-Q piston permits direct pressure comparison.
            auto integralMean = [&](const BemExcitationResult &e, unsigned face) {
                Complex p{};
                double s = 0;
                if (face == 4) {
                    for (std::size_t i = 0; i < areas.size(); ++i) {
                        p += areas[i] * e.boundary()->exteriorSurfacePressure[i];
                        s += areas[i];
                    }
                } else {
                    auto m = surfaceMeasures(mesh);
                    for (std::size_t i = 0; i < mesh.triangles.size(); ++i)
                        if (mesh.triangles[i].faceId == face) {
                            p += m.areas[i] * e.boundary()->surfacePressure[i];
                            s += m.areas[i];
                        }
                }
                return p / s;
            };
            auto z11 = integralMean(a, 1), z12 = integralMean(b, 1), z21 = integralMean(a, 4),
                 z22 = integralMean(b, 4);
            close(z12, z21, .015, "Full cavity/exterior mutual load violates reciprocity.");
            require(z11.real() > 0 && z22.real() > 0 &&
                        z11.real() * z22.real() > std::norm((z12 + std::conj(z21)) * .5),
                    "Exterior/cavity load is not passive.");
            require(std::abs(z12) > .001 * std::abs(z22), "Mutual load was omitted.");
            loads.push_back({{z11, z12, z21, z22}, {}});
            std::cout << "mutual " << f.frequency << " Hz Z12=" << z12 << " Z21=" << z21 << '\n';
        }
        lem::Network net;
        net.nodes = {{{1}, lem::Domain::electrical, true, {}, {}},
                     {{2}, lem::Domain::electrical, false, {}, {}},
                     {{3}, lem::Domain::acoustic, true, {}, {}},
                     {{4}, lem::Domain::acoustic, false, {}, {}},
                     {{5}, lem::Domain::acoustic, false, {}, {}}};
        net.elements = {
            {{1}, "Voltage", {}, lem::EffortSource{{{2}, {1}}, {1, .2}}},
            {{2},
             "Driver",
             {},
             lem::Driver{{{2}, {1}}, {5}, {4}, {5.6, 50e-6, 1.8, .00008, .0002, .04, area}}},
            {{3}, "Rear", {}, lem::BemPort{{{4}, {3}}, {1}, 1}},
            {{4}, "Front", {}, lem::BemPort{{{5}, {3}}, {4}, 1}}};
        input.excitations.clear();
        input.coupling = net;
        const auto coupled = solver.solve(input);
        for (std::size_t i = 0; i < coupled.frequencies.size(); ++i) {
            const auto &e = coupled.frequencies[i].excitations[0];
            auto expected = lem::solve(net, input.frequencies[i], loads[i],
                                       {input.medium.density, input.medium.soundSpeed});
            close(e.lemSolution->drivers[0].current, expected.drivers[0].current, 1e-9,
                  "Full coupled electrical current differs from independent modal closure.");
            for (unsigned o = 0; o < input.observers.size(); ++o)
                close(e.observerPressure[o],
                      basis.frequencies[i].excitations[0].observerPressure[o] *
                              expected.ports[0].volumeVelocity +
                          basis.frequencies[i].excitations[1].observerPressure[o] *
                              expected.ports[1].volumeVelocity,
                      1e-9, "Coupled far field omits an exterior term.");
            require(e.boundaryResidual < input.numerics.relativeResidualTolerance,
                    "Coupled block residual failed.");
        }
        input.frequencies = {500};
        input.farField.reset();
        input.observers = {{.01, 0, .2}, {.04, 0, .15}};
        const auto finite = solver.solve(input).frequencies[0].excitations[0];
        const auto evaluated = evaluatePressureField(mesh, input, 500, finite, input.observers);
        for (unsigned o = 0; o < input.observers.size(); ++o)
            close(finite.observerPressure[o], evaluated[o].pressure, 1e-12,
                  "Saved exterior source field changed.");
        auto invalid = exterior;
        for (auto &t : invalid.triangles)
            std::swap(t.vertices[1], t.vertices[2]);
        bool rejected = false;
        try {
            bem_detail::validateBaffleExteriorMesh(mesh, invalid, {}, {0, 0, 1});
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Reversed external source orientation was accepted.");
        std::cout << "Full exterior baffle coupling tests passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
