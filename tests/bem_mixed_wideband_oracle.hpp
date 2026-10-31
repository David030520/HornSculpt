#pragma once

#include "bem_hf_test_oracle.hpp"
#include <numeric>
#include <span>

namespace hornsculpt::test {
// Two separated groups of four closed complete P1 components. The full
// support tree has electrically large coarse boxes and small leaf supports;
// direct actions below use independent Duffy triangle integration of the
// original curl-regularized weak form, with no scalar IBP/MLFMA translation.
inline WidebandFarBlockFixture
makeMixedWidebandFixture(double k, unsigned gaussOrder,
                         std::span<const unsigned char> farPressureMask = {},
                         std::span<const unsigned char> farFluxMask = {}) {
    using namespace wideband_oracle;
    const bool masked = !farPressureMask.empty();
    if (masked != !farFluxMask.empty() ||
        (!farPressureMask.empty() && farPressureMask.size() != 32 * 32) ||
        (!farFluxMask.empty() && farFluxMask.size() != 32 * 32))
        throw std::invalid_argument(
            "Mixed wideband oracle needs both global P1/P1 and P1/DP0 masks or neither.");
    WidebandFarBlockFixture result;
    result.wavenumber = k;
    SurfaceMesh mesh;
    const std::array<Vec3, 4> tetra{{{0, 0, 0}, {.00525, 0, 0}, {0, .00425, 0}, {0, 0, .00325}}};
    for (unsigned side = 0; side < 2; ++side)
        for (unsigned x = 0; x < 2; ++x)
            for (unsigned y = 0; y < 2; ++y) {
                const Vec3 center{side * 1.2 + (x ? .05 : -.05), y ? .05 : -.05, 0};
                const auto first = static_cast<std::uint32_t>(mesh.vertices.size());
                for (const auto point : tetra)
                    mesh.vertices.push_back(
                        {point.x + center.x, point.y + center.y, point.z + center.z});
                for (const std::array<std::uint32_t, 3> indices :
                     {std::array<std::uint32_t, 3>{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}})
                    mesh.triangles.push_back(
                        {{first + indices[0], first + indices[1], first + indices[2]},
                         side ? 5u : 1u});
            }
    result.space = galerkin::buildSpace(mesh);
    for (auto &panel : result.space.panels)
        panel.boundary = BoundaryCondition::driven;
    const unsigned sourceCount = masked ? 32u : 16u, firstTarget = masked ? 0u : 16u;
    result.sourcePressure.resize(sourceCount);
    result.sourceFlux.resize(sourceCount);
    result.targetPressure.resize(32 - firstTarget);
    std::iota(result.sourcePressure.begin(), result.sourcePressure.end(), 0);
    std::iota(result.sourceFlux.begin(), result.sourceFlux.end(), 0);
    std::iota(result.targetPressure.begin(), result.targetPressure.end(), firstTarget);
    const auto rule = gauss(gaussOrder);
    for (const auto &[u, wu] : rule)
        for (const auto &[v, wv] : rule)
            result.triangleRule.push_back({{(1 - u) * (1 - v), u, (1 - u) * v}, wu * wv * (1 - u)});
    result.pressure.assign(2, std::vector<C>(32));
    result.flux.assign(2, std::vector<C>(32));
    result.expectedPressure.assign(2, std::vector<C>(32));
    result.expectedFlux.assign(2, std::vector<C>(32));
    for (unsigned rhs = 0; rhs < 2; ++rhs)
        for (unsigned j = 0; j < sourceCount; ++j) {
            result.pressure[rhs][j] = {.2 + std::cos(.61 * j + .37 * rhs),
                                       -.3 + std::sin(.89 * j - .41 * rhs)};
            // The unmasked 16->16 fixture preserves its original excitation.
            // Full masked actions drive both physical face groups with
            // different face-constant complex Neumann fluxes.
            result.flux[rhs][j] =
                j < 16 ? (rhs ? C{-.4, .8} : C{.7, -.2}) : (rhs ? C{.55, .25} : C{-.35, .45});
        }
    constexpr C imaginary{0, 1};
    for (unsigned ti = firstTarget; ti < 32; ++ti)
        for (unsigned si = 0; si < sourceCount; ++si) {
            const auto &t = result.space.panels[ti], &s = result.space.panels[si];
            if (masked) {
                bool pressurePresent = false, fluxPresent = false;
                for (unsigned i = 0; i < 3; ++i) {
                    fluxPresent |= farFluxMask[t.nodes[i] * 32 + si] != 0;
                    for (unsigned j = 0; j < 3; ++j)
                        pressurePresent |= farPressureMask[t.nodes[i] * 32 + s.nodes[j]] != 0;
                }
                // In particular, coincident/touching panels have no actual
                // far products. Skip before evaluating any singular Green
                // sample: multiplying a NaN by a zero mask is insufficient.
                if (!pressurePresent && !fluxPresent)
                    continue;
                if (galerkin::classifyPair(t, s).topology != galerkin::PairTopology::disjoint)
                    throw std::invalid_argument(
                        "Mixed far oracle mask includes touching panel supports.");
            }
            for (const auto &x : result.triangleRule)
                for (const auto &y : result.triangleRule) {
                    const auto delta = difference(t, x.barycentric, s, y.barycentric);
                    const double radius = std::sqrt(dot(delta, delta)),
                                 weight = x.weight * y.weight * 4 * t.area * s.area;
                    const C green = std::exp(-imaginary * (k * radius)) /
                                    (4 * std::numbers::pi * radius),
                            radial = green * (1. + imaginary * (k * radius)) / (radius * radius),
                            sourceDerivative = radial * dot(delta, s.normal),
                            targetDerivative = -radial * dot(delta, t.normal);
                    for (unsigned rhs = 0; rhs < 2; ++rhs) {
                        for (unsigned i = 0; i < 3; ++i) {
                            C p{}, curlX{}, curlY{}, curlZ{};
                            for (unsigned j = 0; j < 3; ++j) {
                                if (masked && !farPressureMask[t.nodes[i] * 32 + s.nodes[j]])
                                    continue;
                                const auto value = result.pressure[rhs][s.nodes[j]];
                                p += value * y.barycentric[j];
                                curlX += value * s.curls[j].x;
                                curlY += value * s.curls[j].y;
                                curlZ += value * s.curls[j].z;
                            }
                            const C weakH =
                                green * (t.curls[i].x * curlX + t.curls[i].y * curlY +
                                         t.curls[i].z * curlZ -
                                         k * k * x.barycentric[i] * dot(t.normal, s.normal) * p);
                            result.expectedPressure[rhs][t.nodes[i]] +=
                                weight *
                                (-x.barycentric[i] * sourceDerivative * p - imaginary * weakH / k);
                            if (!masked || farFluxMask[t.nodes[i] * 32 + si])
                                result.expectedFlux[rhs][t.nodes[i]] +=
                                    weight * x.barycentric[i] * result.flux[rhs][si] *
                                    (-green + imaginary * targetDerivative / k);
                        }
                    }
                }
        }
    return result;
}
} // namespace hornsculpt::test
