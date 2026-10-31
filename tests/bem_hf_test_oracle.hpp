#pragma once

#include <galerkin/pair_quadrature.hpp>
#include <galerkin/space.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <numbers>
#include <vector>

namespace hornsculpt::test {
// Independent CPU action oracle for a complete, geometrically separated BEM
// block. Neither production plane-wave formulas nor production triangle rules
// are used. Arrays for each RHS use the full Space's global DOF numbering.
struct WidebandFarBlockFixture {
    galerkin::Space space;
    Vec3 sourceCenter{}, targetCenter{.37, .11, .07};
    std::vector<std::uint32_t> sourcePressure{0, 1, 2, 3}, targetPressure{4, 5, 6, 7},
        sourceFlux{0, 1, 2, 3};
    std::vector<galerkin::TriangleQuadraturePoint> triangleRule;
    std::vector<std::vector<std::complex<double>>> pressure, flux, expectedPressure, expectedFlux;
    double wavenumber{};
};
namespace wideband_oracle {
using C = std::complex<double>;
inline Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline std::vector<std::array<double, 2>> gauss(unsigned count) {
    std::vector<std::array<double, 2>> result;
    for (unsigned index = 0; index < count; ++index) {
        double x = std::cos(std::numbers::pi * (index + .75) / (count + .5)), derivative{};
        for (unsigned iteration = 0; iteration < 40; ++iteration) {
            double p0 = 1, p1 = x;
            for (unsigned degree = 2; degree <= count; ++degree) {
                const double p = ((2 * degree - 1) * x * p1 - (degree - 1) * p0) / degree;
                p0 = p1;
                p1 = p;
            }
            derivative = count * (x * p1 - p0) / (x * x - 1);
            const double step = p1 / derivative;
            x -= step;
            if (std::abs(step) < 3e-16)
                break;
        }
        result.push_back({(1 + x) / 2, 1 / ((1 - x * x) * derivative * derivative)});
    }
    return result;
}
inline Vec3 difference(const galerkin::Panel &x, const double *bx, const galerkin::Panel &y,
                       const double *by) {
    Vec3 result{};
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            const auto delta = subtract(x.vertices[i], y.vertices[j]);
            const double weight = bx[i] * by[j];
            result.x += weight * delta.x;
            result.y += weight * delta.y;
            result.z += weight * delta.z;
        }
    return result;
}
} // namespace wideband_oracle

inline WidebandFarBlockFixture makeWidebandFarBlockFixture(double k, unsigned gaussOrder = 12,
                                                           Vec3 targetCenter = {.37, .11, .07}) {
    using namespace wideband_oracle;
    WidebandFarBlockFixture result;
    result.wavenumber = k;
    result.targetCenter = targetCenter;
    SurfaceMesh mesh;
    const std::array<Vec3, 4> tetra{{{0, 0, 0}, {.021, 0, 0}, {0, .017, 0}, {0, 0, .013}}};
    for (const Vec3 offset : {result.sourceCenter, result.targetCenter}) {
        const auto first = static_cast<std::uint32_t>(mesh.vertices.size());
        for (const auto p : tetra)
            mesh.vertices.push_back({p.x + offset.x, p.y + offset.y, p.z + offset.z});
        for (const std::array<std::uint32_t, 3> ids :
             {std::array<std::uint32_t, 3>{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}})
            mesh.triangles.push_back({{first + ids[0], first + ids[1], first + ids[2]}, first + 1});
    }
    result.space = galerkin::buildSpace(mesh);
    const auto rule = gauss(gaussOrder);
    for (const auto &[u, wu] : rule)
        for (const auto &[v, wv] : rule)
            result.triangleRule.push_back({{(1 - u) * (1 - v), u, (1 - u) * v}, wu * wv * (1 - u)});
    result.pressure = {{{1.1, .2}, {-.3, .7}, {.2, -.8}, {-.5, -.1}, {}, {}, {}, {}},
                       {{-.7, .3}, {.8, -.2}, {.1, .5}, {-.2, -.9}, {}, {}, {}, {}}};
    result.flux = {{{.6, -.4}, {-.3, .2}, {.1, .7}, {-.8, -.2}, {}, {}, {}, {}},
                   {{-.2, .8}, {.4, -.6}, {-.5, .1}, {.3, .2}, {}, {}, {}, {}}};
    result.expectedPressure.assign(2, std::vector<C>(8));
    result.expectedFlux.assign(2, std::vector<C>(8));
    constexpr C imaginary{0, 1};
    for (unsigned targetPanel = 4; targetPanel < 8; ++targetPanel)
        for (unsigned sourcePanel = 0; sourcePanel < 4; ++sourcePanel) {
            const auto &test = result.space.panels[targetPanel],
                       &trial = result.space.panels[sourcePanel];
            for (const auto &x : result.triangleRule)
                for (const auto &y : result.triangleRule) {
                    const Vec3 delta = difference(test, x.barycentric, trial, y.barycentric);
                    const double radius = std::sqrt(dot(delta, delta)),
                                 weight = x.weight * y.weight * 4 * test.area * trial.area;
                    const C green = std::exp(-imaginary * (k * radius)) /
                                    (4 * std::numbers::pi * radius),
                            radial = green * (1.0 + imaginary * (k * radius)) / (radius * radius),
                            sourceDerivative = radial * dot(delta, trial.normal),
                            targetDerivative = -radial * dot(delta, test.normal);
                    for (unsigned rhs = 0; rhs < 2; ++rhs) {
                        C p{}, curlX{}, curlY{}, curlZ{};
                        for (unsigned j = 0; j < 3; ++j) {
                            const C value = result.pressure[rhs][trial.nodes[j]];
                            p += value * y.barycentric[j];
                            curlX += value * trial.curls[j].x;
                            curlY += value * trial.curls[j].y;
                            curlZ += value * trial.curls[j].z;
                        }
                        for (unsigned i = 0; i < 3; ++i) {
                            const C weakH =
                                green *
                                (test.curls[i].x * curlX + test.curls[i].y * curlY +
                                 test.curls[i].z * curlZ -
                                 k * k * x.barycentric[i] * dot(test.normal, trial.normal) * p);
                            result.expectedPressure[rhs][test.nodes[i]] +=
                                weight *
                                (-x.barycentric[i] * sourceDerivative * p - imaginary * weakH / k);
                            result.expectedFlux[rhs][test.nodes[i]] +=
                                weight * x.barycentric[i] * result.flux[rhs][sourcePanel] *
                                (-green + imaginary * targetDerivative / k);
                        }
                    }
                }
        }
    return result;
}
} // namespace hornsculpt::test
