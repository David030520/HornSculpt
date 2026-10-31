#pragma once
#include "bem_hf_test_oracle.hpp"
#include <array>
#include <complex>
#include <galerkin/pair_quadrature.hpp>
#include <galerkin/regular_qualification.hpp>
#include <galerkin/space.hpp>
#include <numbers>
#include <span>
#include <stdexcept>
namespace hornsculpt::test::support_operator_oracle {
using Complex = std::complex<double>;
inline double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline Vec3 delta(const galerkin::Panel &t, const galerkin::Panel &s,
                  const galerkin::PairQuadraturePoint &p) {
    Vec3 result{};
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            const auto weight = p.test[i] * p.trial[j];
            result.x += (t.vertices[i].x - s.vertices[j].x) * weight;
            result.y += (t.vertices[i].y - s.vertices[j].y) * weight;
            result.z += (t.vertices[i].z - s.vertices[j].z) * weight;
        }
    return result;
}
struct Matrices {
    std::vector<Complex> a, b;
};
inline std::vector<galerkin::TriangleQuadraturePoint> independentRegularRule(unsigned gaussOrder) {
    const auto gauss = wideband_oracle::gauss(gaussOrder);
    std::vector<galerkin::TriangleQuadraturePoint> result;
    result.reserve(gauss.size() * gauss.size());
    for (const auto &[u, wu] : gauss)
        for (const auto &[v, wv] : gauss)
            result.push_back({{1 - u, u * (1 - v), u * v}, wu * wv * u});
    return result;
}
template <class Integrate>
void forIndependentRegularPair(std::span<const galerkin::TriangleQuadraturePoint> rule,
                               Integrate &&integrate) {
    for (const auto &x : rule)
        for (const auto &y : rule) {
            galerkin::PairQuadraturePoint point{};
            point.weight = x.weight * y.weight;
            for (unsigned i = 0; i < 3; ++i) {
                point.test[i] = x.barycentric[i];
                point.trial[i] = y.barycentric[i];
            }
            integrate(point);
        }
}
inline Matrices directReference(const galerkin::Space &space, double k,
                                galerkin::PairQuadratureSettings settings) {
    const auto n = space.vertices.size(), p = space.panels.size();
    Matrices result;
    result.a.resize(n * n);
    result.b.resize(n * 2);
    constexpr Complex imaginary{0, 1};
    const auto regular12 = independentRegularRule(12), regular16 = independentRegularRule(16);
    std::vector<galerkin::RegularPanelBounds> bounds;
    for (const auto &panel : space.panels)
        bounds.push_back(galerkin::makeRegularPanelBounds(panel));
    galerkin::AdaptiveQuadratureWorkspace workspace(256 * 1024 * 1024);
    // Independent complex weak-form algebra. The established CPU quadrature
    // supplies identical integration domains; no CUDA moment helper is used.
    for (std::size_t ti = 0; ti < p; ++ti)
        for (std::size_t si = 0; si < p; ++si) {
            const auto &t = space.panels[ti], &s = space.panels[si];
            std::array<Complex, 9> a{};
            std::array<Complex, 3> b{};
            const auto integrate = [&](const galerkin::PairQuadraturePoint &point) {
                const auto d = delta(t, s, point);
                const auto r = std::sqrt(dot(d, d));
                const auto g = std::exp(-imaginary * k * r) / (4 * std::numbers::pi * r);
                const auto radial = g * (1. + imaginary * k * r) / (r * r);
                const auto sourceDerivative = ti == si ? Complex{} : radial * dot(d, s.normal);
                const auto targetDerivative = ti == si ? Complex{} : -radial * dot(d, t.normal);
                const auto weight = point.weight * 4 * t.area * s.area;
                for (unsigned i = 0; i < 3; ++i) {
                    b[i] += weight * point.test[i] * (-g + imaginary * targetDerivative / k);
                    for (unsigned j = 0; j < 3; ++j) {
                        const auto basis = point.test[i] * point.trial[j];
                        const auto weakH = (dot(t.curls[i], s.curls[j]) -
                                            k * k * dot(t.normal, s.normal) * basis) *
                                           g;
                        a[3 * i + j] +=
                            weight * (-basis * sourceDerivative - imaginary * weakH / k);
                    }
                }
            };
            const bool disjoint =
                galerkin::classifyPair(t, s).topology == galerkin::PairTopology::disjoint;
            if (disjoint && !galerkin::isNearPair(t, s, settings.nearDistanceRatio) &&
                galerkin::qualifiedRegularDegree(bounds[ti], bounds[si], k,
                                                 settings.relativeTolerance,
                                                 settings.regularOrder) != 0) {
                forIndependentRegularPair(regular16, integrate);
                const auto fineA = a;
                const auto fineB = b;
                a = {};
                b = {};
                forIndependentRegularPair(regular12, integrate);
                const auto refined = [](const auto &low, const auto &high) {
                    double difference{}, norm{};
                    for (std::size_t i = 0; i < low.size(); ++i) {
                        difference += std::norm(low[i] - high[i]);
                        norm += std::norm(high[i]);
                    }
                    return norm > 0 ? difference <= 1e-18 * norm : difference == 0;
                };
                if (!refined(a, fineA) || !refined(b, fineB))
                    throw std::runtime_error(
                        "Independent regular Duffy12/16 oracle did not converge.");
                a = fineA;
                b = fineB;
            } else {
                galerkin::PairQuadratureRule rule;
                if (disjoint) {
                    galerkin::AdaptivePairQuadratureState state(t, s, settings, k, workspace);
                    while (!state.done()) {
                        const auto tasks = state.requests();
                        std::array<galerkin::PairPatchEstimate, 4> estimates;
                        for (std::size_t i = 0; i < tasks.size(); ++i)
                            estimates[i] = galerkin::estimatePairPatchCpu(t, s, tasks[i], k);
                        state.commit(std::span(estimates).first(tasks.size()));
                    }
                    rule = state.takeResult();
                } else
                    rule = galerkin::makePairQuadrature(t, s, settings, k);
                for (const auto &point : rule.points)
                    integrate(point);
                for (const auto &patch : rule.patches) {
                    const auto quadrature = galerkin::triangleRule(patch.degree);
                    for (const auto &x : quadrature)
                        for (const auto &y : quadrature)
                            integrate(galerkin::mapPairPatch(x, y, patch));
                }
            }
            if (ti == si)
                for (unsigned i = 0; i < 3; ++i) {
                    b[i] += imaginary * t.area / (6 * k);
                    for (unsigned j = 0; j < 3; ++j)
                        a[3 * i + j] += .5 * t.area * (i == j ? 2. : 1.) / 12.;
                }
            const auto face = s.boundary == BoundaryCondition::wall ? UINT32_MAX
                              : s.faceId == 5                       ? 0u
                                                                    : 1u;
            for (unsigned i = 0; i < 3; ++i) {
                const auto row = t.nodes[i];
                for (unsigned j = 0; j < 3; ++j)
                    result.a[std::size_t(row) * n + s.nodes[j]] += a[3 * i + j];
                if (face != UINT32_MAX)
                    result.b[std::size_t(row) * 2 + face] += b[i];
            }
        }
    return result;
}
} // namespace hornsculpt::test::support_operator_oracle
