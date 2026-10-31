#include "bem_hf_test_oracle.hpp"
#include <algorithm>
#include <galerkin/regular_qualification.hpp>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
using namespace hornsculpt;
using C = std::complex<double>;
using Rule = std::vector<galerkin::TriangleQuadraturePoint>;
constexpr C imaginary{0, 1};
constexpr double budget = 1.25e-6;
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
// Independent differentiated free-space Green function. The scalar form is
// evaluated by direct double quadrature, without an angular/FMM expansion or
// surface integration by parts. The reference retains the curl weak form.
std::array<std::vector<C>, 2> scalarAction(const test::WidebandFarBlockFixture &fixture,
                                           const Rule &rule, unsigned rhs) {
    using namespace test::wideband_oracle;
    std::array<std::vector<C>, 2> action{std::vector<C>(8), std::vector<C>(8)};
    const double k = fixture.wavenumber;
    for (unsigned target = 4; target < 8; ++target)
        for (unsigned source = 0; source < 4; ++source) {
            const auto &xPanel = fixture.space.panels[target],
                       &yPanel = fixture.space.panels[source];
            for (const auto &x : rule)
                for (const auto &y : rule) {
                    const auto delta = difference(xPanel, x.barycentric, yPanel, y.barycentric);
                    const double r2 = dot(delta, delta), r = std::sqrt(r2),
                                 weight = 4 * xPanel.area * yPanel.area * x.weight * y.weight,
                                 dnX = dot(delta, xPanel.normal) / r,
                                 dnY = dot(delta, yPanel.normal) / r;
                    const C green = std::exp(-imaginary * (k * r)) / (4 * std::numbers::pi * r),
                            radial = (1. + imaginary * (k * r)) / r,
                            sourceDerivative = green * radial * dnY,
                            targetDerivative = -green * radial * dnX,
                            hessian =
                                green * ((3 / r2 + imaginary * (3 * k / r) - k * k) * dnX * dnY -
                                         (1 / r2 + imaginary * (k / r)) *
                                             dot(xPanel.normal, yPanel.normal));
                    C pressure{};
                    for (unsigned j = 0; j < 3; ++j)
                        pressure += fixture.pressure[rhs][yPanel.nodes[j]] * y.barycentric[j];
                    for (unsigned i = 0; i < 3; ++i) {
                        action[0][xPanel.nodes[i]] += weight * x.barycentric[i] * pressure *
                                                      (-sourceDerivative - imaginary * hessian / k);
                        action[1][xPanel.nodes[i]] += weight * x.barycentric[i] *
                                                      fixture.flux[rhs][source] *
                                                      (-green + imaginary * targetDerivative / k);
                    }
                }
        }
    return action;
}
double relative(std::span<const C> actual, std::span<const C> expected) {
    double error{}, norm{};
    for (unsigned i = 4; i < 8; ++i) {
        error += std::norm(actual[i] - expected[i]);
        norm += std::norm(expected[i]);
    }
    require(norm > 0, "Far quadrature reference action is degenerate.");
    return std::sqrt(error / norm);
}
using Moments = std::array<C, 36>;
Rule independentRule(unsigned degree) {
    Rule result;
    const auto gauss = test::wideband_oracle::gauss(degree);
    for (const auto &[u, wu] : gauss)
        for (const auto &[v, wv] : gauss)
            result.push_back({{1 - u, u * (1 - v), u * v}, wu * wv * u});
    return result;
}
Moments independentMoments(const galerkin::Panel &testPanel, const galerkin::Panel &sourcePanel,
                           double k, const Rule &rule) {
    using namespace test::wideband_oracle;
    Moments result{};
    const double h = std::max(testPanel.maximumEdgeLength, sourcePanel.maximumEdgeLength);
    for (const auto &x : rule)
        for (const auto &y : rule) {
            const auto delta = difference(testPanel, x.barycentric, sourcePanel, y.barycentric);
            const double r2 = dot(delta, delta), r = std::sqrt(r2);
            const C green = std::exp(-imaginary * k * r) / r;
            const C derivative = -green * (1. + imaginary * k * r) * h / r2;
            const std::array<C, 4> components{green, derivative * delta.x, derivative * delta.y,
                                              derivative * delta.z};
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j)
                    for (unsigned component = 0; component < 4; ++component)
                        result[4 * (3 * i + j) + component] += x.weight * y.weight *
                                                               x.barycentric[i] * y.barycentric[j] *
                                                               components[component];
        }
    return result;
}
double momentError(const Moments &actual, const Moments &expected) {
    double error{}, norm{};
    for (unsigned i = 0; i < 36; ++i) {
        error = std::max(error, std::abs(actual[i] - expected[i]));
        norm = std::max(norm, std::abs(expected[i]));
    }
    require(norm > 0, "Regular qualification reference is zero.");
    return error / norm;
}
void diagnoseRegularQualification() {
    using namespace test::wideband_oracle;
    const auto coarseRule = independentRule(12), finalRule = independentRule(16);
    const std::array<unsigned, 4> degrees{4, 6, 8, 14};
    const std::array<Rule, 4> rules{galerkin::triangleRule(4), galerkin::triangleRule(6),
                                    galerkin::triangleRule(8), galerkin::triangleRule(14)};
    double maximumReferenceDifference{};
    std::array<std::size_t, 4> selected{};
    for (const double kh : {0., .25, .375, .5, .75, 1., 1.5, 2., 2.5, 3.25})
        for (const double gapRatio : {0., .5, 1., 2., 3., 4., 6., 8., 16.}) {
            std::array<double, 4> errors{};
            for (unsigned shape = 0; shape < 5; ++shape) {
                galerkin::Panel x{}, y{};
                x.vertices = {Vec3{-.5, -.25, 0}, Vec3{.5, -.25, 0}, Vec3{0, .5, 0}};
                // Holdouts include near-minimum-angle triangles in two
                // different planes, not only the original broad triangles.
                if (shape >= 3)
                    x.vertices = {Vec3{-.5, -.06333333333333333, 0},
                                  Vec3{.5, -.06333333333333333, 0}, Vec3{0, .12666666666666667, 0}};
                if (shape >= 3)
                    y.vertices = {Vec3{-.5, 0, -.06333333333333333},
                                  Vec3{.5, 0, -.06333333333333333}, Vec3{0, 0, .12666666666666667}};
                else if (shape == 0)
                    y.vertices = {Vec3{-.5, -.25, 0}, Vec3{.5, -.25, 0}, Vec3{0, .5, 0}};
                else if (shape == 1)
                    y.vertices = {Vec3{-.4, 0, -.2}, Vec3{.4, .1, -.2}, Vec3{0, -.1, .4}};
                else
                    y.vertices = {Vec3{-.1, -.05, 0}, Vec3{.1, -.05, 0}, Vec3{0, .1, 0}};
                const auto radiusAndEdge = [](const auto &vertices) {
                    double radius{}, edge{};
                    for (unsigned i = 0; i < 3; ++i) {
                        radius = std::max(radius, std::sqrt(dot(vertices[i], vertices[i])));
                        const auto d = subtract(vertices[i], vertices[(i + 1) % 3]);
                        edge = std::max(edge, std::sqrt(dot(d, d)));
                    }
                    return std::array<double, 2>{radius, edge};
                };
                const auto xt = radiusAndEdge(x.vertices), ys = radiusAndEdge(y.vertices);
                x.maximumEdgeLength = xt[1];
                y.maximumEdgeLength = ys[1];
                const double h = std::max(xt[1], ys[1]), distance = xt[0] + ys[0] + gapRatio * h;
                const Vec3 direction = shape == 0   ? Vec3{0, 0, 1}
                                       : shape == 4 ? Vec3{0, .6, .8}
                                                    : Vec3{.6, 0, .8};
                for (auto &vertex : y.vertices) {
                    vertex.x += distance * direction.x;
                    vertex.y += distance * direction.y;
                    vertex.z += distance * direction.z;
                }
                x.centroid = {};
                y.centroid = {distance * direction.x, distance * direction.y,
                              distance * direction.z};
                const auto expected = independentMoments(x, y, kh / h, finalRule);
                maximumReferenceDifference =
                    std::max(maximumReferenceDifference,
                             momentError(independentMoments(x, y, kh / h, coarseRule), expected));
                for (unsigned i = 0; i < degrees.size(); ++i)
                    errors[i] =
                        std::max(errors[i],
                                 momentError(independentMoments(x, y, kh / h, rules[i]), expected));
                const auto xb = galerkin::makeRegularPanelBounds(x),
                           yb = galerkin::makeRegularPanelBounds(y);
                const auto degree = galerkin::qualifiedRegularDegree(xb, yb, kh / h, 1e-4);
                ++selected[degree == 4 ? 0 : degree == 6 ? 1 : degree == 8 ? 2 : 3];
                if (degree) {
                    const auto actual =
                        independentMoments(x, y, kh / h, galerkin::triangleRule(degree));
                    require(momentError(actual, expected) < .2e-4,
                            "Qualified regular rule exceeded the unchanged near moment budget.");
                }
                require(galerkin::qualifiedRegularDegree(xb, yb, kh / h, 1e-8) == 0,
                        "Strict tolerances must use adaptive integration outside the envelope.");
                const auto strictDegree = galerkin::qualifiedRegularDegree(xb, yb, kh / h, 1e-5);
                if (strictDegree)
                    require(momentError(independentMoments(x, y, kh / h,
                                                           galerkin::triangleRule(strictDegree)),
                                        expected) < .2e-5,
                            "Stricter qualified envelope exceeded its unchanged moment budget.");
            }
            std::cout << "Regular qualification kh=" << kh << ", sphereGap/h=" << gapRatio;
            for (unsigned i = 0; i < degrees.size(); ++i)
                std::cout << ", degree" << degrees[i] << "=" << errors[i];
            std::cout << '\n';
        }
    std::cout << "Regular qualification independent Duffy12/16=" << maximumReferenceDifference
              << '\n';
    require(maximumReferenceDifference < 1e-9,
            "Regular qualification oracle failed independent quadrature refinement.");
    require(std::all_of(selected.begin(), selected.end(), [](auto count) { return count > 0; }),
            "Regular qualification must exercise all rules and adaptive fallback.");
    require(
        galerkin::minimumEndpointDegree(.01) == 4 && galerkin::minimumEndpointDegree(.75) == 6 &&
            galerkin::minimumEndpointDegree(1.75) == 8 && galerkin::minimumEndpointDegree(3) == 12,
        "Electrical endpoint policy must preserve cheap well-resolved rules.");
    bool coarseRejected = false;
    try {
        static_cast<void>(galerkin::minimumEndpointDegree(3.251));
    } catch (const std::runtime_error &) {
        coarseRejected = true;
    }
    require(coarseRejected, "Unqualified electrical panel size must fail explicitly.");
}
void checkElectricalPanelQuadrature() {
    // Budgets are fixed independently of the measured results. This is an
    // electrically coarse panel diagnostic, not a claim of mesh convergence.
    constexpr double k = 2 * std::numbers::pi * 20000 / 343;
    const auto reference = test::makeWidebandFarBlockFixture(k, 16),
               independent = test::makeWidebandFarBlockFixture(k, 12);
    double weakRefinement{}, scalarRefinement{}, maximumLowError{};
    double h{};
    for (const auto &panel : reference.space.panels)
        h = std::max(h, panel.maximumEdgeLength);
    std::cout << std::setprecision(12) << "20kHz far-panel diagnostic: k*hmax=" << k * h
              << ", independent quadrature budget=" << budget << '\n';
    for (unsigned rhs = 0; rhs < 2; ++rhs) {
        weakRefinement =
            std::max({weakRefinement,
                      relative(independent.expectedPressure[rhs], reference.expectedPressure[rhs]),
                      relative(independent.expectedFlux[rhs], reference.expectedFlux[rhs])});
        const auto refinedScalar = scalarAction(reference, independent.triangleRule, rhs),
                   finalScalar = scalarAction(reference, reference.triangleRule, rhs);
        scalarRefinement =
            std::max({scalarRefinement, relative(refinedScalar[0], reference.expectedPressure[rhs]),
                      relative(refinedScalar[1], reference.expectedFlux[rhs]),
                      relative(finalScalar[0], reference.expectedPressure[rhs]),
                      relative(finalScalar[1], reference.expectedFlux[rhs])});
        for (const unsigned degree : {4u, 8u, 14u}) {
            const auto rule = galerkin::triangleRule(degree);
            const auto actual = scalarAction(reference, rule, rhs);
            const double a = relative(actual[0], reference.expectedPressure[rhs]),
                         b = relative(actual[1], reference.expectedFlux[rhs]);
            std::cout << "20kHz panel degree=" << degree << " rhs=" << rhs
                      << " points=" << rule.size() << " scalar/original-weak A/B=" << a << '/' << b
                      << " budgetPassed=" << (a <= budget && b <= budget) << '\n';
            if (degree == 4)
                maximumLowError = std::max({maximumLowError, a, b});
        }
    }
    std::cout << "20kHz original-weak Duffy12/16=" << weakRefinement
              << ", independent scalar/original-weak=" << scalarRefinement << '\n';
    require(weakRefinement < 1e-9 && scalarRefinement < 1e-8,
            "20kHz independent original-weak/scalar quadrature failed to converge.");
    require(maximumLowError > budget, "Electrical-panel fixture failed to expose the low-rule "
                                      "error missed by a seed-only audit.");
}
void checkQualificationHoldouts() {
    using namespace test::wideband_oracle;
    const auto referenceRule = independentRule(16), checkRule = independentRule(12);
    std::size_t qualified{}, adaptive{};
    for (const double ratio : {1e-4, .02, 1., 50., 1e4})
        for (const double kh : {.499999, .500001, 1.999999, 2.000001})
            for (const double gap : {.499999, .500001, .999999, 1.000001, 2.999999, 3.000001}) {
                const double xt = 1 / std::max(1., ratio), ys = ratio * xt;
                galerkin::Panel a{}, b{};
                a.vertices = {Vec3{-.5 * xt, -.25 * xt, 0}, Vec3{.5 * xt, -.25 * xt, 0},
                              Vec3{0, .5 * xt, 0}};
                b.vertices = {Vec3{-.5 * ys, 0, -.25 * ys}, Vec3{.5 * ys, 0, -.25 * ys},
                              Vec3{0, 0, .5 * ys}};
                a.maximumEdgeLength = xt;
                b.maximumEdgeLength = ys;
                const auto ba = galerkin::makeRegularPanelBounds(a),
                           bb = galerkin::makeRegularPanelBounds(b);
                const double distance = ba.radius + bb.radius + gap;
                // A non-coordinate-aligned rigid rotation plus a large common
                // translation. Geometric quality limits are not assumed on.
                const auto transform = [](Vec3 v) {
                    return Vec3{1e6 + .6 * v.x - .8 * v.z, -1e6 + .8 * v.x + .6 * v.z, 1e6 + v.y};
                };
                for (auto &v : a.vertices)
                    v = transform(v);
                for (auto &v : b.vertices)
                    v = transform({v.x, v.y + .6 * distance, v.z + .8 * distance});
                a.centroid = transform({});
                b.centroid = transform({0, .6 * distance, .8 * distance});
                const auto degree =
                    galerkin::qualifiedRegularDegree(galerkin::makeRegularPanelBounds(a),
                                                     galerkin::makeRegularPanelBounds(b), kh, 1e-4);
                if (!degree) {
                    ++adaptive;
                    continue;
                }
                ++qualified;
                const auto reference = independentMoments(a, b, kh, referenceRule);
                require(momentError(independentMoments(a, b, kh, checkRule), reference) < 1e-9,
                        "Adversarial regular oracle did not converge independently.");
                require(
                    momentError(independentMoments(a, b, kh, galerkin::triangleRule(degree)),
                                reference) < .2e-4,
                    "Qualified rule failed an independent threshold/size/rigid-motion holdout.");
            }
    galerkin::Panel slender{};
    slender.vertices = {Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{.5, 1e-6, 0}};
    slender.centroid = {.5, 1e-6 / 3, 0};
    slender.maximumEdgeLength = 1;
    const auto bounds = galerkin::makeRegularPanelBounds(slender);
    require(!bounds.qualifiedShape &&
                galerkin::qualifiedRegularDegree(bounds, bounds, .1, 1e-4) == 0,
            "Valid slender panels must remain available through adaptive integration.");
    require(qualified > 0 && adaptive > 0,
            "Holdouts must exercise qualification and genuine adaptive fallbacks.");
    std::cout << "Regular adversarial holdouts qualified/adaptive=" << qualified << '/' << adaptive
              << '\n';
}
} // namespace
int main() {
    try {
        checkElectricalPanelQuadrature();
        checkQualificationHoldouts();
        diagnoseRegularQualification();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
