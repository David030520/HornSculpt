#include "bem_hf_test_oracle.hpp"
#include "reference/scaled_projection.hpp"
#include <fmm/scaled_translation.hpp>
#include <fmm/support_level_plan.hpp>
#include <fmm/wideband_transfer.hpp>
#include <galerkin/space.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

using namespace hornsculpt;
using namespace hornsculpt::galerkin;
using namespace hornsculpt::test::reference;
namespace {
using C = std::complex<double>;
using CV = std::array<C, 3>;
constexpr C imaginary{0, 1};
constexpr double pi = std::numbers::pi;
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3 sub(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 scale(Vec3 a, double b) {
    return {a.x * b, a.y * b, a.z * b};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
double length(Vec3 a) {
    return std::sqrt(dot(a, a));
}
double axis(Vec3 a, unsigned i) {
    return i == 0 ? a.x : i == 1 ? a.y : a.z;
}
C dot(const CV &a, const CV &b) {
    C value{};
    for (unsigned i = 0; i < 3; ++i)
        value += a[i] * b[i];
    return value;
}
C dot(const CV &a, Vec3 b) {
    return a[0] * b.x + a[1] * b.y + a[2] * b.z;
}
double norm(const CV &value) {
    return std::sqrt(std::norm(value[0]) + std::norm(value[1]) + std::norm(value[2]));
}
void close(C actual, C expected, double absoluteTolerance, const char *message) {
    if (std::abs(actual - expected) > absoluteTolerance) {
        std::cerr << std::setprecision(17) << message << ": " << actual << " versus " << expected
                  << ", tolerance " << absoluteTolerance << '\n';
        throw std::runtime_error(message);
    }
}
// Newton/Legendre quadrature is independent of the production triangle rules
// and of the angular-grid implementation. Weights integrate [0,1].
std::vector<std::array<double, 2>> gauss(unsigned count) {
    std::vector<std::array<double, 2>> result;
    for (unsigned index = 0; index < count; ++index) {
        double x = std::cos(pi * (index + .75) / (count + .5)), derivative{};
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
SurfaceMesh tetrahedron(Vec3 offset = {}, double size = 1) {
    SurfaceMesh mesh;
    mesh.vertices = {offset, add(offset, {.021 * size, 0, 0}), add(offset, {0, .017 * size, 0}),
                     add(offset, {0, 0, .013 * size})};
    mesh.triangles = {{{0, 2, 1}, 1}, {{0, 1, 3}, 1}, {{0, 3, 2}, 1}, {{1, 2, 3}, 1}};
    return mesh;
}
struct Sample {
    Vec3 localPoint, normal, curl;
    double phi{}, weight{};
};
// A complete pressure basis support is its entire vertex star. In particular,
// its contribution is not truncated at the boundary of an open panel cluster.
std::vector<Sample> supportSamples(const Space &space, unsigned node, Vec3 center, unsigned degree,
                                   bool openPanel = false) {
    const auto rule = gauss(degree);
    std::vector<Sample> result;
    for (auto incidence = space.incidenceOffsets[node];
         incidence < space.incidenceOffsets[node + 1]; ++incidence) {
        const auto entry = space.incidences[incidence];
        const auto &panel = space.panels[entry.panel];
        const Vec3 a = sub(panel.vertices[0], center),
                   b = sub(panel.vertices[1], panel.vertices[0]),
                   c = sub(panel.vertices[2], panel.vertices[0]);
        for (const auto &[u, wu] : rule)
            for (const auto &[v, wv] : rule) {
                const std::array bary{(1 - u) * (1 - v), u, (1 - u) * v};
                result.push_back({add(a, add(scale(b, bary[1]), scale(c, bary[2]))), panel.normal,
                                  panel.curls[entry.local], bary[entry.local],
                                  wu * wv * (1 - u) * 2 * panel.area});
            }
        if (openPanel)
            break;
    }
    return result;
}
struct Moments {
    C scalar{};
    CV normal{}, curl{};
    double absoluteCurl{}, absoluteScalar{};
};
Moments moments(const std::vector<Sample> &samples, Vec3 direction, double k, double sign) {
    Moments result;
    for (const auto &sample : samples) {
        const C phase = std::exp(imaginary * (sign * k * dot(direction, sample.localPoint)));
        result.scalar += phase * (sample.weight * sample.phi);
        for (unsigned i = 0; i < 3; ++i) {
            result.normal[i] += phase * (sample.weight * sample.phi * axis(sample.normal, i));
            result.curl[i] += phase * (sample.weight * axis(sample.curl, i));
        }
        result.absoluteCurl += sample.weight * length(sample.curl);
        result.absoluteScalar += sample.weight * sample.phi;
    }
    return result;
}
double ibpError(const Moments &value, Vec3 d, double k, double sign) {
    const CV crossMoment{value.normal[1] * d.z - value.normal[2] * d.y,
                         value.normal[2] * d.x - value.normal[0] * d.z,
                         value.normal[0] * d.y - value.normal[1] * d.x};
    CV defect{};
    for (unsigned i = 0; i < 3; ++i)
        defect[i] = value.curl[i] + imaginary * (sign * k) * crossMoment[i];
    return norm(defect) / (value.absoluteCurl + k * value.absoluteScalar);
}
void completeSupportIdentity() {
    const auto space = buildSpace(tetrahedron({}, 10));
    const Vec3 d = scale({.371, .529, .763}, 1 / length({.371, .529, .763}));
    for (const double k : {.1, 1., 10., 100., 2 * pi * 20000 / 343}) {
        const auto source = moments(supportSamples(space, 0, {}, 64), d, k, 1);
        const auto test = moments(supportSamples(space, 1, {}, 64), d, k, -1);
        const auto open = moments(supportSamples(space, 0, {}, 64, true), d, k, 1);
        require(ibpError(source, d, k, 1) < 2e-13, "Complete source support IBP failed.");
        require(ibpError(test, d, k, -1) < 2e-13, "Complete test support IBP failed.");
        require(ibpError(open, d, k, 1) > 1e-5,
                "Open panel counterexample did not expose the missing boundary term.");
        const C h = dot(test.curl, source.curl) - k * k * dot(test.normal, source.normal);
        const C scalarH = -k * k * dot(test.normal, d) * dot(source.normal, d);
        const double absoluteScale = (test.absoluteCurl + k * test.absoluteScalar) *
                                     (source.absoluteCurl + k * source.absoluteScalar);
        close(h, scalarH, 2e-13 * absoluteScale, "Scalar hypersingular factorization failed.");
        const C a = -imaginary * k * test.scalar * dot(source.normal, d) - imaginary * h / k;
        const C scalarA =
            imaginary * k * (dot(test.normal, d) - test.scalar) * dot(source.normal, d);
        close(a, scalarA, 2e-13 * absoluteScale / k, "Scalar Burton-Miller factorization failed.");
    }
}
struct AngularSample {
    Vec3 direction;
    C weightedTransfer;
};
// Validate the production radial recurrence against independent C++ special
// functions before using its angular rule in direct Green/Galerkin oracles.
std::vector<AngularSample> diagonalRule(Vec3 separation, double k, unsigned cutoff) {
    const double distance = length(separation), z = k * distance;
    const Vec3 unit = scale(separation, 1 / distance);
    const auto sphere = wideband::makeSphereRule(cutoff + 1);
    const auto production = wideband::diagonalTranslation(sphere, cutoff, k, separation);
    std::vector<AngularSample> result;
    for (std::size_t q = 0; q < sphere.samples.size(); ++q) {
        const Vec3 d = sphere.samples[q].direction;
        C transfer{}, phase{1, 0};
        for (unsigned degree = 0; degree <= cutoff; ++degree) {
            const C hankel{std::sph_bessel(degree, z), -std::sph_neumann(degree, z)};
            transfer += static_cast<double>(2 * degree + 1) * phase * hankel *
                        std::legendre(degree, dot(d, unit));
            phase *= -imaginary;
        }
        const C independent = -imaginary * k * transfer / (16 * pi * pi);
        close(production[q], independent, 3e-12 * (1 + std::abs(independent)),
              "Production diagonal radial recurrence differs from independent special functions.");
        result.push_back({d, sphere.samples[q].weight * production[q]});
    }
    return result;
}
struct KernelValues {
    C single{}, sourceNormal{}, testNormal{}, hypersingular{};
};
KernelValues closedKernel(Vec3 displacement, Vec3 nx, Vec3 ny, double k) {
    const double radius = length(displacement);
    const Vec3 unit = scale(displacement, 1 / radius);
    const C green = std::exp(-imaginary * (k * radius)) / (4 * pi * radius);
    const C radial = (1.0 + imaginary * (k * radius)) / radius;
    // H=-d_nx d_ny G = nx.Hessian(G).ny, independently differentiated
    // from exp(-ikr)/r; no surface-curl integration by parts is used here.
    return {green, green * radial * dot(unit, ny), -green * radial * dot(unit, nx),
            green * ((3. / (radius * radius) + imaginary * (3 * k / radius) - k * k) *
                         dot(unit, nx) * dot(unit, ny) -
                     (1. / (radius * radius) + imaginary * (k / radius)) * dot(nx, ny))};
}
void diagonalGreen() {
    const Vec3 separation{.37, .11, .07}, nx = scale({.2, -.5, .7}, 1 / length({.2, -.5, .7})),
                                          ny = scale({-.3, .8, .1}, 1 / length({-.3, .8, .1}));
    for (const double k : {200., 2 * pi * 20000 / 343}) {
        const auto rule = diagonalRule(separation, k, 30);
        for (const Vec3 relative : {Vec3{}, Vec3{.009, -.007, .004}, Vec3{-.012, .011, -.003}}) {
            KernelValues actual;
            for (const auto &entry : rule) {
                const C value = entry.weightedTransfer *
                                std::exp(-imaginary * (k * dot(entry.direction, relative)));
                const double dx = dot(entry.direction, nx), dy = dot(entry.direction, ny);
                actual.single += value;
                actual.sourceNormal += imaginary * k * dy * value;
                actual.testNormal -= imaginary * k * dx * value;
                actual.hypersingular -= k * k * dx * dy * value;
            }
            const auto expected = closedKernel(add(separation, relative), nx, ny, k);
            const double reference = std::abs(expected.single);
            close(actual.single, expected.single, 2e-11 * reference, "Diagonal Green failed.");
            close(actual.sourceNormal, expected.sourceNormal, 2e-11 * k * reference,
                  "Diagonal source derivative failed.");
            close(actual.testNormal, expected.testNormal, 2e-11 * k * reference,
                  "Diagonal test derivative failed.");
            close(actual.hypersingular, expected.hypersingular, 2e-11 * k * k * reference,
                  "Diagonal hypersingular derivative failed.");
        }
    }
}
void automaticNearMacKernel(bool actualGeometry = false) {
    using namespace hornsculpt::wideband;
    const double k = 2 * pi * 20000 / 343, activeRadius = .21, separationRatio = 1.50001,
                 budget = 5e-6;
    const auto plan = makeWidebandLevelPlan(activeRadius, k,
                                            actualGeometry ? 2 * separationRatio - 1 : 1.5, budget);
    require(plan.representation == LevelRepresentation::planeWave,
            "Automatic HF test did not exercise the diagonal representation.");
    require(makeWidebandLevelPlan(.001, k, 1.5, budget).representation ==
                LevelRepresentation::regularSpherical,
            "Automatic policy applied an unstable diagonal form in the LF regime.");
    const Vec3 unit = scale({.8, .5, .3}, 1 / length({.8, .5, .3})),
               separation = scale(unit, 2 * activeRadius * separationRatio),
               nx = scale({.2, -.5, .7}, 1 / length({.2, -.5, .7})),
               ny = scale({-.3, .8, .1}, 1 / length({-.3, .8, .1}));
    double maximumError{}, maximumDiscrepancy{};
    std::array<std::vector<KernelValues>, 2> results;
    const std::array<Vec3, 4> offsets{
        Vec3{}, scale(unit, -1.96 * activeRadius), scale(unit, 1.96 * activeRadius),
        scale({-.5, .8, 0}, 1.8 * activeRadius / length({-.5, .8, 0}))};
    for (unsigned refinement = 0; refinement < 2; ++refinement) {
        const auto rule = makeSphereRule(plan.quadratureOrder + 4 * refinement);
        const auto transfer =
            diagonalTranslation(rule, plan.translationOrder + 4 * refinement, k, separation);
        for (const auto offset : offsets) {
            KernelValues actual;
            for (std::size_t q = 0; q < rule.samples.size(); ++q) {
                const auto d = rule.samples[q].direction;
                const C value = rule.samples[q].weight * transfer[q] *
                                std::exp(-imaginary * (k * dot(d, offset)));
                actual.single += value;
                actual.sourceNormal += imaginary * k * dot(d, ny) * value;
                actual.testNormal -= imaginary * k * dot(d, nx) * value;
                actual.hypersingular -= k * k * dot(d, nx) * dot(d, ny) * value;
            }
            const auto exact = closedKernel(add(separation, offset), nx, ny, k);
            const double magnitude = std::abs(exact.single);
            maximumError = std::max(
                {maximumError, std::abs(actual.single - exact.single) / magnitude,
                 std::abs(actual.sourceNormal - exact.sourceNormal) / (k * magnitude),
                 std::abs(actual.testNormal - exact.testNormal) / (k * magnitude),
                 std::abs(actual.hypersingular - exact.hypersingular) / (k * k * magnitude)});
            results[refinement].push_back(actual);
        }
    }
    for (std::size_t j = 0; j < offsets.size(); ++j) {
        const auto &a = results[0][j], &b = results[1][j];
        const double magnitude =
            std::abs(closedKernel(add(separation, offsets[j]), nx, ny, k).single);
        maximumDiscrepancy =
            std::max({maximumDiscrepancy, std::abs(a.single - b.single) / magnitude,
                      std::abs(a.sourceNormal - b.sourceNormal) / (k * magnitude),
                      std::abs(a.testNormal - b.testNormal) / (k * magnitude),
                      std::abs(a.hypersingular - b.hypersingular) / (k * k * magnitude)});
    }
    std::cout << "Automatic near-MAC HF kernel: ka=" << k * activeRadius
              << " ratio=" << separationRatio
              << " moment/translation/quadrature=" << plan.signatureOrder << '/'
              << plan.translationOrder << '/' << plan.quadratureOrder
              << " q=" << plan.angularSamples << " exact error=" << maximumError
              << " refined discrepancy=" << maximumDiscrepancy << '\n';
    require(maximumError < budget && maximumDiscrepancy < .25 * budget,
            "Automatic HF bandwidth fails the near-MAC direct kernel or refinement budget.");
}
C independentHarmonic(unsigned l, int m, Vec3 d) {
    const auto absoluteM = static_cast<unsigned>(std::abs(m));
    const double normalization =
        std::sqrt((2. * l + 1) / (4 * pi) *
                  std::exp(std::lgamma(l - absoluteM + 1.) - std::lgamma(l + absoluteM + 1.)));
    const C positive = normalization * (absoluteM % 2 ? -1. : 1.) *
                       std::assoc_legendre(l, absoluteM, d.z) *
                       std::exp(imaginary * (absoluteM * std::atan2(d.y, d.x)));
    return m < 0 ? (absoluteM % 2 ? -1. : 1.) * std::conj(positive) : positive;
}
std::size_t harmonicIndex(unsigned l, int m) {
    return static_cast<std::size_t>(l) * l + l + m;
}
void angularTransforms() {
    const auto childRule = wideband::makeSphereRule(8), parentRule = wideband::makeSphereRule(16);
    double weight = 0;
    for (const auto &sample : childRule.samples) {
        weight += sample.weight;
        const auto actual = wideband::sphericalHarmonics(8, sample.direction);
        for (unsigned l = 0; l <= 8; ++l)
            for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m)
                close(
                    actual[harmonicIndex(l, m)], independentHarmonic(l, m, sample.direction), 3e-13,
                    "Normalized spherical harmonic disagrees with independent Legendre functions.");
    }
    close(weight, 4 * pi, 4e-13, "Sphere rule does not integrate its surface area.");
    const Vec3 source{.012, -.008, .005}, shift{.011, -.009, .007};
    const double k = 12;
    std::vector<C> outgoing, parentIncoming;
    for (const auto &sample : childRule.samples)
        outgoing.push_back(std::exp(imaginary * (k * dot(sample.direction, source))));
    for (std::size_t q = 0; q < parentRule.samples.size(); ++q) {
        const auto d = parentRule.samples[q].direction;
        parentIncoming.push_back(C{.2 + std::sin(1.3 * static_cast<double>(q)),
                                   -.1 + std::cos(.7 * static_cast<double>(q))} +
                                 independentHarmonic(3, -2, d));
    }
    const auto plan = wideband::makeAngularTransferPlan(childRule, parentRule, 8);
    const auto resampled = wideband::resampleSpectrum(outgoing, plan);
    const auto up = wideband::outgoingToParent(outgoing, childRule, parentRule, 8, k, shift);
    const auto down = wideband::incomingToChild(parentIncoming, parentRule, childRule, 8, k, shift);
    C left{}, right{};
    for (std::size_t q = 0; q < up.size(); ++q) {
        close(resampled[q],
              std::exp(imaginary * (k * dot(parentRule.samples[q].direction, source))), 3e-13,
              "Separable resampling changed a resolved point-source spectrum.");
        close(up[q],
              std::exp(imaginary * (k * dot(parentRule.samples[q].direction, add(source, shift)))),
              3e-13, "Outgoing parent shift has an incorrect sign or normalization.");
        left += std::conj(up[q]) * parentIncoming[q] * parentRule.samples[q].weight;
    }
    for (std::size_t q = 0; q < down.size(); ++q)
        right += std::conj(outgoing[q]) * down[q] * childRule.samples[q].weight;
    close(left, right, 2e-12, "Angular down-transfer is not the weighted Hermitian adjoint.");
    const std::vector<C> constant(parentRule.samples.size(), C{1, 0});
    const auto shiftedIncoming =
        wideband::incomingToChild(constant, parentRule, childRule, 8, k, shift);
    for (std::size_t q = 0; q < shiftedIncoming.size(); ++q)
        close(shiftedIncoming[q],
              std::exp(-imaginary * (k * dot(childRule.samples[q].direction, shift))), 3e-13,
              "Incoming child shift has an incorrect physical phase.");
    // LF/HF handoff is tested against actual regular waves and a closed-form
    // point-source phase, including a very small argument where a diagonal
    // Hankel translation must not be used.
    const Vec3 direction = scale(source, 1 / length(source));
    for (const double handoffK : {1e-7, 1., 100.}) {
        const auto radial = wideband::regularRadialTable(8, handoffK * length(source));
        std::vector<C> multipoles(81);
        for (unsigned l = 0; l <= 8; ++l)
            for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m)
                multipoles[harmonicIndex(l, m)] =
                    radial.values[l] * std::conj(independentHarmonic(l, m, direction));
        const auto spectrum = wideband::multipoleToSpectrum(multipoles, 8, childRule);
        const double truncationTolerance = handoffK == 100 ? 2e-6 : 3e-13;
        for (std::size_t q = 0; q < spectrum.size(); ++q)
            close(spectrum[q],
                  std::exp(imaginary * (handoffK * dot(childRule.samples[q].direction, source))),
                  truncationTolerance, "Outgoing regular-wave to spectrum handoff failed.");
    }
    std::vector<C> mode;
    for (const auto &sample : childRule.samples)
        mode.push_back(independentHarmonic(3, -2, sample.direction));
    const auto local = wideband::spectrumToLocal(mode, childRule, 8);
    for (unsigned l = 0; l <= 8; ++l)
        for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m)
            close(local[harmonicIndex(l, m)], l == 3 && m == -2 ? C{0, 4 * pi} : C{}, 3e-12,
                  "Incoming spectrum-to-local handoff failed.");
}
void regularRadialReference() {
    const auto origin = wideband::regularRadialTable(32, 0);
    for (unsigned l = 0; l <= 32; ++l) {
        close(origin.values[l], l == 0 ? 1. : 0., 0, "Regular radial origin value failed.");
        close(origin.derivatives[l], l == 1 ? 1. / 3 : 0., 0,
              "Regular radial origin derivative failed.");
    }
    const auto tiny = wideband::regularRadialTable(8, 1e-300);
    close(tiny.derivatives[1], 1. / 3, 1e-15, "Tiny-argument regular derivative was lost.");
    close(tiny.derivatives[2], 2e-300 / 15, 1e-313,
          "Tiny-argument higher regular derivative was lost.");
    for (const double x : {.001, .1, .999, 1.001, pi, 7.725251836937707, 10., 50., 100.}) {
        const auto actual = wideband::regularRadialTable(64, x);
        for (unsigned l = 0; l <= 64; ++l) {
            const double expected = std::sph_bessel(l, x), next = std::sph_bessel(l + 1, x),
                         previous = l ? std::sph_bessel(l - 1, x) : expected,
                         derivative = l == 0 ? -next : l * expected / x - next;
            // Neighboring regular modes supply a meaningful absolute scale at
            // roots of j_l or its derivative without hiding tiny evanescent modes.
            const double modeScale =
                std::max({std::abs(expected), std::abs(previous), std::abs(next), 1e-250});
            close(actual.values[l], expected, 5e-12 * modeScale,
                  "Regular radial recursion disagrees with independent special functions.");
            close(actual.derivatives[l], derivative,
                  5e-12 * std::max(std::abs(derivative), modeScale),
                  "Regular radial derivative disagrees with independent special functions.");
        }
    }
}
void fullSupportGalerkin() {
    const Vec3 separation{.37, .11, .07};
    const auto sourceSpace = buildSpace(tetrahedron()),
               testSpace = buildSpace(tetrahedron(separation));
    const auto source = supportSamples(sourceSpace, 0, {}, 12),
               test = supportSamples(testSpace, 1, separation, 12);
    for (const double k : {200., 2 * pi * 20000 / 343}) {
        C directA{}, directB{}, normalH{}, weakH{};
        for (const auto &x : test)
            for (const auto &y : source) {
                const Vec3 displacement = add(separation, sub(x.localPoint, y.localPoint));
                const auto kernel = closedKernel(displacement, x.normal, y.normal, k);
                const double weight = x.weight * y.weight;
                const C weak = kernel.single * (dot(x.curl, y.curl) -
                                                k * k * x.phi * y.phi * dot(x.normal, y.normal));
                weakH += weight * weak;
                normalH += weight * x.phi * y.phi * kernel.hypersingular;
                directA += weight * (-x.phi * y.phi * kernel.sourceNormal - imaginary * weak / k);
                // A DP0 source panel is not integrated by parts. Retaining one
                // complete panel suffices to exercise its distinct source space.
                if (y.normal.z < -.5)
                    directB +=
                        weight * x.phi * (-kernel.single + imaginary * kernel.testNormal / k);
            }
        C diagonalA{}, diagonalB{};
        for (const auto &entry : diagonalRule(separation, k, 30)) {
            const auto x = moments(test, entry.direction, k, -1),
                       y = moments(source, entry.direction, k, 1);
            C flux{};
            for (const auto &sample : source)
                if (sample.normal.z < -.5)
                    flux += sample.weight *
                            std::exp(imaginary * (k * dot(entry.direction, sample.localPoint)));
            const C tested = dot(x.normal, entry.direction) - x.scalar;
            diagonalA +=
                entry.weightedTransfer * imaginary * k * tested * dot(y.normal, entry.direction);
            diagonalB += entry.weightedTransfer * tested * flux;
        }
        const double aScale =
            k * sourceSpace.surfaceArea * testSpace.surfaceArea / (4 * pi * length(separation));
        close(weakH, normalH, 3e-11 * k * aScale,
              "Weak and independently differentiated full-support H disagree.");
        close(diagonalA, directA, 3e-11 * aScale,
              "Scalar diagonal pressure operator disagrees with direct Galerkin integration.");
        close(diagonalB, directB, 3e-11 * aScale / k,
              "Scalar diagonal flux operator disagrees with direct Galerkin integration.");
        require(std::abs(directA) > aScale * 1e-4 && std::abs(directB) > aScale / k * 1e-4,
                "Galerkin reference is too small to exercise the operator signs.");
    }
}
void farBlockOracleConvergence() {
    for (const double k : {200., 2 * pi * 20000 / 343}) {
        const auto coarse = hornsculpt::test::makeWidebandFarBlockFixture(k, 12),
                   fine = hornsculpt::test::makeWidebandFarBlockFixture(k, 16);
        for (unsigned rhs = 0; rhs < 2; ++rhs) {
            double pressureNorm{}, fluxNorm{}, pressureError{}, fluxError{};
            for (const auto target : fine.targetPressure) {
                pressureNorm += std::norm(fine.expectedPressure[rhs][target]);
                fluxNorm += std::norm(fine.expectedFlux[rhs][target]);
                pressureError += std::norm(coarse.expectedPressure[rhs][target] -
                                           fine.expectedPressure[rhs][target]);
                fluxError +=
                    std::norm(coarse.expectedFlux[rhs][target] - fine.expectedFlux[rhs][target]);
            }
            require(pressureNorm > 0 && fluxNorm > 0, "Full far-block oracle is vacuous.");
            require(std::sqrt(pressureError / pressureNorm) < 1e-10,
                    "Full far-block pressure oracle lacks quadrature convergence.");
            require(std::sqrt(fluxError / fluxNorm) < 1e-10,
                    "Full far-block flux oracle lacks quadrature convergence.");
        }
    }
}
void scaledLowFrequencyTranslations() {
    using namespace hornsculpt::wideband;
    const Vec3 sourcePoint{.012, -.006, .004}, upShift{.015, .008, -.006},
        separation{.37, .11, .07}, childShift{.008, .003, -.002}, targetPoint{.003, -.002, .004};
    const C strength{.7, -.4};
    for (const double k : {0., 1e-7, 10.}) {
        const auto basis = scaledRegularWaves(10, k, .025, sourcePoint);
        ScaledMultipole multipole{10, .025, std::vector<C>(basis.waves.size())};
        for (std::size_t j = 0; j < basis.waves.size(); ++j)
            multipole.coefficients[j] = strength * std::conj(basis.waves[j].value);
        const auto sphere = makeSphereRule(24);
        const auto outgoingSpectrum = scaledMultipoleToSpectrum(multipole, k, sphere);
        for (std::size_t q = 0; q < outgoingSpectrum.size(); ++q)
            close(outgoingSpectrum[q],
                  strength *
                      std::exp(imaginary * (k * dot(sphere.samples[q].direction, sourcePoint))),
                  5e-13 * std::abs(strength), "Scaled outgoing LF/HF handoff failed.");
        const std::vector<C> incomingConstant(sphere.samples.size(), strength);
        const auto incomingLocal = spectrumToScaledLocal(incomingConstant, sphere, 10, k, .025);
        const auto incomingValue = evaluateScaledLocal(incomingLocal, k, targetPoint);
        const double incomingArgument = k * length(targetPoint);
        const double incomingJ0 =
            incomingArgument == 0 ? 1 : std::sin(incomingArgument) / incomingArgument;
        close(incomingValue.value, 4 * pi * strength * incomingJ0, 5e-13 * std::abs(strength),
              "Scaled incoming HF/LF handoff failed.");
        const auto parent = scaledMultipoleToParent(multipole, k, upShift, .06, 10, 24);
        const auto local = scaledMultipoleToLocal(parent, k, separation, .04, 10, 24),
                   refinedLocal = scaledMultipoleToLocal(parent, k, separation, .04, 10, 32),
                   child = scaledLocalToChild(local, k, childShift, .018, 10, 24);
        const Vec3 worldPoint = add(separation, add(childShift, targetPoint)),
                   displacement = sub(worldPoint, add(upShift, sourcePoint));
        const double radius = length(displacement);
        const C expected = strength * std::exp(-imaginary * (k * radius)) / (4 * pi * radius),
                radial = -expected * (1.0 + imaginary * (k * radius)) / (radius * radius);
        const auto actual = evaluateScaledLocal(child, k, targetPoint),
                   refined = evaluateScaledLocal(refinedLocal, k, add(childShift, targetPoint)),
                   outgoing = evaluateScaledMultipole(parent, k, worldPoint);
        for (const auto value : {actual.value, refined.value, outgoing.value})
            close(value, expected, 2e-11 * std::abs(expected),
                  "Stable LF M2M/M2L/L2L disagrees with direct outgoing Green.");
        for (unsigned axisIndex = 0; axisIndex < 3; ++axisIndex) {
            const C derivative = radial * axis(displacement, axisIndex);
            close(actual.gradient[axisIndex], derivative, 2e-10 * std::abs(expected) / radius,
                  "Stable LF hierarchy has an incorrect gradient.");
            close(actual.gradient[axisIndex], refined.gradient[axisIndex],
                  2e-10 * std::abs(expected) / radius,
                  "Stable LF projection lacks independent angular refinement convergence.");
        }
    }
}
void scaledAxialProjectionReference() {
    using namespace hornsculpt::wideband;
    const double k = 10;
    ScaledMultipole source{8, .02, std::vector<C>(81)};
    ScaledLocal parent{8, .06, std::vector<C>(81)};
    for (std::size_t j = 0; j < source.coefficients.size(); ++j) {
        const double degree = std::floor(std::sqrt(static_cast<double>(j)));
        source.coefficients[j] = std::pow(.2, degree) * C{std::sin(.7 * j), std::cos(1.1 * j)};
        parent.coefficients[j] = std::pow(.2, degree) * C{std::cos(.6 * j), std::sin(1.3 * j)};
    }
    const auto up = makeScaledAxialPlan(ScaledAxialKind::multipoleToParent, 8, 8, k, .02, .06, .015,
                                        25),
               across = makeScaledAxialPlan(ScaledAxialKind::multipoleToLocal, 8, 8, k, .02, .04,
                                            .2, 25),
               down =
                   makeScaledAxialPlan(ScaledAxialKind::localToChild, 8, 8, k, .06, .02, .015, 25);
    const auto expectedUp = scaledMultipoleToParent(source, k, {0, 0, .015}, .06, 8, 24);
    const auto expectedAcross = scaledMultipoleToLocal(source, k, {0, 0, .2}, .04, 8, 24);
    const auto expectedDown = scaledLocalToChild(parent, k, {0, 0, .015}, .02, 8, 24);
    const std::array actual{applyScaledAxialPlan(up, source.coefficients),
                            applyScaledAxialPlan(across, source.coefficients),
                            applyScaledAxialPlan(down, parent.coefficients)};
    const std::array expected{expectedUp.coefficients, expectedAcross.coefficients,
                              expectedDown.coefficients};
    for (std::size_t operation = 0; operation < actual.size(); ++operation) {
        double error{}, magnitude{};
        for (std::size_t j = 0; j < actual[operation].size(); ++j) {
            error += std::norm(actual[operation][j] - expected[operation][j]);
            magnitude += std::norm(expected[operation][j]);
        }
        require(std::sqrt(error / magnitude) < 3e-11,
                "Axial m-block plan disagrees with independent full-sphere projection.");
    }
}
void scaledTransitionAndZeros() {
    using namespace hornsculpt::wideband;
    constexpr double radius = .04;
    const C strength{.7, -.4};
    for (const double electrical :
         {1., pi * (1 - 1e-9), pi * (1 + 1e-9), 4.493409457909064 * (1 - 1e-9),
          4.493409457909064 * (1 + 1e-9), 10., 20., 30., 32., 38., 40., 100.}) {
        const double k = electrical / radius;
        const unsigned order = electrical > 40
                                   ? makeWidebandLevelPlan(radius, k, 1.5, 1e-10).signatureOrder
                                   : static_cast<unsigned>(std::ceil(electrical)) + 20,
                       polar = 2 * order + 9;
        ScaledMultipole monopole{0, radius, {strength / std::sqrt(4 * pi)}};
        const auto plan = makeScaledAxialPlan(ScaledAxialKind::multipoleToLocal, 0, order, k,
                                              radius, radius, .2, polar),
                   refinedPlan = makeScaledAxialPlan(ScaledAxialKind::multipoleToLocal, 0, order, k,
                                                     radius, radius, .2, polar + 8);
        const ScaledLocal local{order, radius, applyScaledAxialPlan(plan, monopole.coefficients)},
            refined{order, radius, applyScaledAxialPlan(refinedPlan, monopole.coefficients)};
        const Vec3 point{.019, -.015, .021}, displacement = add({0, 0, .2}, point);
        const double distance = length(displacement);
        const C exact = strength * std::exp(-imaginary * (k * distance)) / (4 * pi * distance),
                derivative = -exact * (1.0 + imaginary * (k * distance)) / (distance * distance);
        const auto actual = evaluateScaledLocal(local, k, point),
                   other = evaluateScaledLocal(refined, k, point);
        double error = std::abs(actual.value - exact) / std::abs(exact),
               discrepancy = std::abs(actual.value - other.value) / std::abs(exact);
        for (unsigned component = 0; component < 3; ++component) {
            error = std::max(error, std::abs(actual.gradient[component] -
                                             derivative * axis(displacement, component)) /
                                        ((k + 1 / distance) * std::abs(exact)));
            discrepancy = std::max(
                discrepancy, std::abs(actual.gradient[component] - other.gradient[component]) /
                                 ((k + 1 / distance) * std::abs(exact)));
        }
        const auto upward = makeScaledAxialPlan(ScaledAxialKind::multipoleToParent, 0, order, k,
                                                .01, radius, .015, polar);
        const ScaledMultipole translated{order, radius,
                                         applyScaledAxialPlan(upward, monopole.coefficients)};
        const auto upField = evaluateScaledMultipole(translated, k, displacement);
        const Vec3 shifted = sub(displacement, {0, 0, .015});
        const double shiftedDistance = length(shifted);
        const C shiftedExact =
            strength * std::exp(-imaginary * (k * shiftedDistance)) / (4 * pi * shiftedDistance);
        error = std::max(error, std::abs(upField.value - shiftedExact) / std::abs(shiftedExact));
        std::cout << "Scaled axial transition ka=" << electrical << " p=" << order
                  << " polar=" << polar << " Green/gradient error=" << error
                  << " mu-refinement=" << discrepancy << '\n';
        require(error < 3e-10 && discrepancy < 3e-11,
                "Scaled axial transition or Bessel-zero Cauchy projection failed.");
    }
}
std::pair<std::array<C, 8>, std::array<C, 8>>
scaledFarBlockAction(const hornsculpt::test::WidebandFarBlockFixture &fixture, unsigned rhs,
                     unsigned order, double radius, bool localTranslation,
                     const hornsculpt::wideband::ScaledAxialPlan *axial = nullptr) {
    using namespace hornsculpt::wideband;
    const double k = fixture.wavenumber;
    ScaledMultipole pressure{order, radius, std::vector<C>((order + 1) * (order + 1))},
        flux = pressure;
    for (unsigned panelId = 0; panelId < 4; ++panelId) {
        const auto &panel = fixture.space.panels[panelId];
        for (const auto &sample : fixture.triangleRule) {
            Vec3 point{};
            C p{};
            for (unsigned j = 0; j < 3; ++j) {
                point = add(point, scale(sub(panel.vertices[j], fixture.sourceCenter),
                                         sample.barycentric[j]));
                p += fixture.pressure[rhs][panel.nodes[j]] * sample.barycentric[j];
            }
            const auto regular = scaledRegularWaves(order, k, radius, point);
            const double weight = sample.weight * 2 * panel.area;
            for (std::size_t j = 0; j < regular.waves.size(); ++j) {
                C derivative{};
                for (unsigned a = 0; a < 3; ++a)
                    derivative += axis(panel.normal, a) * regular.waves[j].gradient[a];
                pressure.coefficients[j] += weight * p * std::conj(derivative);
                flux.coefficients[j] +=
                    weight * fixture.flux[rhs][panelId] * std::conj(regular.waves[j].value);
            }
        }
    }
    ScaledLocal aLocal, bLocal;
    if (localTranslation) {
        if (axial) {
            aLocal = {axial->targetOrder, radius,
                      applyScaledAxialPlan(*axial, pressure.coefficients)};
            bLocal = {axial->targetOrder, radius, applyScaledAxialPlan(*axial, flux.coefficients)};
        } else {
            aLocal = scaledMultipoleToLocal(pressure, k, fixture.targetCenter, radius, order,
                                            2 * order + 4);
            bLocal =
                scaledMultipoleToLocal(flux, k, fixture.targetCenter, radius, order, 2 * order + 4);
        }
    }
    std::array<C, 8> a{}, b{};
    for (unsigned panelId = 4; panelId < 8; ++panelId) {
        const auto &panel = fixture.space.panels[panelId];
        for (const auto &sample : fixture.triangleRule) {
            Vec3 point{};
            for (unsigned j = 0; j < 3; ++j)
                point = add(point, scale(sub(panel.vertices[j], fixture.targetCenter),
                                         sample.barycentric[j]));
            const auto pa = localTranslation ? evaluateScaledLocal(aLocal, k, point)
                                             : evaluateScaledMultipole(
                                                   pressure, k, add(fixture.targetCenter, point)),
                       pb =
                           localTranslation
                               ? evaluateScaledLocal(bLocal, k, point)
                               : evaluateScaledMultipole(flux, k, add(fixture.targetCenter, point));
            C na{}, nb{};
            for (unsigned j = 0; j < 3; ++j) {
                na += pa.gradient[j] * axis(panel.normal, j);
                nb += pb.gradient[j] * axis(panel.normal, j);
            }
            for (unsigned j = 0; j < 3; ++j) {
                const double weight = sample.weight * 2 * panel.area * sample.barycentric[j];
                a[panel.nodes[j]] += weight * (-pa.value + imaginary * na / k);
                b[panel.nodes[j]] += weight * (-pb.value + imaginary * nb / k);
            }
        }
    }
    return {a, b};
}
void separateProbeNearMacGalerkin(bool actualGeometry = false) {
    using namespace hornsculpt::wideband;
    constexpr double k = 200, radius = .1, physicalRadius = .021, mac = 1.50001,
                     distance = 2 * physicalRadius * mac, probe = distance / 1.5;
    const auto fixture = hornsculpt::test::makeWidebandFarBlockFixture(k, 14, {0, 0, distance});
    auto policy = makeWidebandLevelPlan(radius, k,
                                        actualGeometry ? distance / physicalRadius - 1 : 1.5, 5e-6);
    enforceSupportDiagonalStability(policy, k, distance);
    require(policy.representation == LevelRepresentation::regularSpherical,
            "Separate-probe regression must exercise the actual LF order policy.");
    const unsigned order = policy.signatureOrder;
    const auto plan = makeScaledAxialPlan(ScaledAxialKind::multipoleToLocal, order, order, k,
                                          radius, radius, distance, 2 * order + 33, probe),
               refined =
                   makeScaledAxialPlan(ScaledAxialKind::multipoleToLocal, order + 2, order + 2, k,
                                       radius, radius, distance, 2 * order + 53, probe);
    require(2 * radius > distance && 2 * physicalRadius * 1.5 < distance &&
                plan.projectionRadius < distance,
            "Separate-probe regression must overlap normalization spheres and obey physical MAC.");
    for (unsigned rhs = 0; rhs < 2; ++rhs) {
        const auto actual = scaledFarBlockAction(fixture, rhs, order, radius, true, &plan),
                   finer = scaledFarBlockAction(fixture, rhs, order + 2, radius, true, &refined);
        double normA{}, normB{}, errorA{}, errorB{}, deltaA{}, deltaB{};
        for (const auto target : fixture.targetPressure) {
            normA += std::norm(fixture.expectedPressure[rhs][target]);
            normB += std::norm(fixture.expectedFlux[rhs][target]);
            errorA += std::norm(actual.first[target] - fixture.expectedPressure[rhs][target]);
            errorB += std::norm(actual.second[target] - fixture.expectedFlux[rhs][target]);
            deltaA += std::norm(actual.first[target] - finer.first[target]);
            deltaB += std::norm(actual.second[target] - finer.second[target]);
        }
        std::cout << "Separate-probe near-MAC LF block p=" << order << " rhs=" << rhs
                  << " relative A/B=" << std::sqrt(errorA / normA) << '/'
                  << std::sqrt(errorB / normB) << " refinement=" << std::sqrt(deltaA / normA) << '/'
                  << std::sqrt(deltaB / normB) << '\n';
        require(std::sqrt(errorA / normA) < 1e-8 && std::sqrt(errorB / normB) < 1e-8 &&
                    std::sqrt(deltaA / normA) < 1e-8 && std::sqrt(deltaB / normB) < 1e-8,
                "Independent Cauchy probe failed complete near-MAC Galerkin/refinement gate.");
    }
}
void asymmetricNearMacProbe(bool actualGeometry = false) {
    using namespace hornsculpt::wideband;
    constexpr double radius = .1, k = 200, distance = .063;
    const double sourceRadius = length({.003, -.002, .035}),
                 targetRadius = length({.001, .001, -.004}),
                 convergence = std::min((distance - targetRadius) / sourceRadius,
                                        (distance - sourceRadius) / targetRadius);
    auto policy = makeWidebandLevelPlan(radius, k, actualGeometry ? convergence : 1.5, 5e-6);
    enforceSupportDiagonalStability(policy, k, distance);
    const unsigned order = policy.signatureOrder;
    const auto plan = makeScaledAxialPlan(ScaledAxialKind::multipoleToLocal, order, order, k,
                                          radius, radius, distance, 2 * order + 33),
               refined =
                   makeScaledAxialPlan(ScaledAxialKind::multipoleToLocal, order + 2, order + 2, k,
                                       radius, radius, distance, 2 * order + 53);
    const std::array configurations{std::array{Vec3{.003, -.002, .035}, Vec3{.001, .001, -.004}},
                                    std::array{Vec3{.001, .001, .004}, Vec3{.003, -.002, -.035}}};
    constexpr C amplitude{.7, -.4};
    for (const auto &points : configurations) {
        const auto sourcePoint = points[0], targetPoint = points[1];
        require(1.5 * (length(sourcePoint) + length(targetPoint)) < distance,
                "Asymmetric-probe fixture must satisfy the actual physical MAC.");
        const auto translate = [&](unsigned p, const ScaledAxialPlan &translation) {
            const auto sourceBasis = scaledRegularWaves(p, k, radius, sourcePoint);
            std::vector<C> moments;
            for (const auto &wave : sourceBasis.waves)
                moments.push_back(amplitude * std::conj(wave.value));
            const ScaledLocal local{p, radius, applyScaledAxialPlan(translation, moments)};
            return evaluateScaledLocal(local, k, targetPoint);
        };
        const auto actual = translate(order, plan), other = translate(order + 2, refined);
        const Vec3 delta = sub(add({0, 0, distance}, targetPoint), sourcePoint);
        const double separation = length(delta);
        const C exact = amplitude * std::exp(-imaginary * (k * separation)) / (4 * pi * separation),
                derivative =
                    -exact * (1.0 + imaginary * (k * separation)) / (separation * separation);
        double error = std::abs(actual.value - exact) / std::abs(exact),
               discrepancy = std::abs(actual.value - other.value) / std::abs(exact);
        for (unsigned component = 0; component < 3; ++component) {
            const double norm = (k + 1 / separation) * std::abs(exact);
            error = std::max(
                error,
                std::abs(actual.gradient[component] - derivative * axis(delta, component)) / norm);
            discrepancy =
                std::max(discrepancy,
                         std::abs(actual.gradient[component] - other.gradient[component]) / norm);
        }
        std::cout << "Asymmetric near-MAC LF probe source/target radius=" << length(sourcePoint)
                  << '/' << length(targetPoint) << " Green/gradient error=" << error
                  << " p/mu refinement=" << discrepancy << '\n';
        require(error < 1.25e-6 && discrepancy < 1.25e-6,
                "Asymmetric source/target Cauchy probe failed physical field/refinement gate.");
    }
}
void scaledLowFrequencyGalerkin() {
    for (const double k : {.01, 1., 10.}) {
        const auto fixture = hornsculpt::test::makeWidebandFarBlockFixture(k, 12);
        for (unsigned rhs = 0; rhs < 2; ++rhs) {
            const auto [a, b] = scaledFarBlockAction(fixture, rhs, 10, .04, true);
            double normA{}, normB{}, errorA{}, errorB{};
            for (const auto target : fixture.targetPressure) {
                normA += std::norm(fixture.expectedPressure[rhs][target]);
                normB += std::norm(fixture.expectedFlux[rhs][target]);
                errorA += std::norm(a[target] - fixture.expectedPressure[rhs][target]);
                errorB += std::norm(b[target] - fixture.expectedFlux[rhs][target]);
            }
            std::cout << "LF full P1/DP0 block k=" << k << " rhs=" << rhs
                      << " relative A/B=" << std::sqrt(errorA / normA) << '/'
                      << std::sqrt(errorB / normB) << '\n';
            require(std::sqrt(errorA / normA) < 1e-8 && std::sqrt(errorB / normB) < 1e-8,
                    "Stable LF complete Galerkin block disagrees with independent weak-form "
                    "quadrature.");
        }
    }
}
void actualConvergenceGeometry() {
    using namespace hornsculpt::wideband;
    auto illConditioned = makeWidebandLevelPlan(.1, 200, 2.00002, 5e-6);
    const auto retainedOrder = illConditioned.signatureOrder;
    require(illConditioned.representation == LevelRepresentation::planeWave,
            "Conditioning regression did not exercise the geometric HF decision.");
    const double indicator = enforceSupportDiagonalStability(illConditioned, 200, .06300042);
    require(illConditioned.representation == LevelRepresentation::regularSpherical &&
                illConditioned.signatureOrder == retainedOrder &&
                indicator * std::numeric_limits<double>::epsilon() > .01 * 5e-6,
            "Ill-conditioned common-level diagonal was accepted or silently truncated.");
    auto stable = makeWidebandLevelPlan(.21, 2 * pi * 20000 / 343, 2.00002, 5e-6);
    enforceSupportDiagonalStability(stable, 2 * pi * 20000 / 343, .6300042);
    require(stable.representation == LevelRepresentation::planeWave,
            "Conditioning guard unnecessarily removed the stable physical HF regime.");
    bool rejected = false;
    try {
        auto outsideDomain = makeWidebandLevelPlan(1, 200, 2.00002, 5e-6);
        enforceSupportDiagonalStability(outsideDomain, 200, .06300042);
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    require(rejected, "Conditioning guard silently forced an unsupported LF electrical domain.");
    SupportTree tree;
    tree.separationRatio = 1.5;
    tree.levels.resize(3);
    tree.nodes.resize(6);
    for (unsigned level = 0; level < 3; ++level)
        for (unsigned body = 0; body < 2; ++body) {
            const unsigned id = 2 * level + body;
            auto &node = tree.nodes[id];
            node.level = level;
            node.center = {body * .300001, 0, level * .01};
            node.parent = level ? id - 2 : noSupportNode;
            node.pressureRadius = .1 - level * .03;
            node.fluxRadius = .12 - level * .03;
            node.pressureCount = node.fluxCount = 1;
            tree.levels[level].nodes.push_back(id);
        }
    tree.fluxLeaf = {4};
    tree.pressurePressure.farInteractions = {{1, 0}};
    tree.pressureFlux.farInteractions = {{1, 0}};
    const std::array<std::uint32_t, 1> absent{UINT32_MAX}, requested{0};
    for (double k : {.1, 10., 50., 200., 1000.})
        for (double seed : {5e-6, 1e-10}) {
            const auto pressure = makeSupportExpansionPlan(tree, absent, k, seed),
                       both = makeSupportExpansionPlan(tree, requested, k, seed);
            const double balanced = (.300001 - .1) / .1, flux = (.300001 - .1) / .12;
            require(
                std::abs(pressure.levels[0].sourceSeparation - balanced) < 1e-14 &&
                    std::abs(pressure.levels[0].targetSeparation - balanced) < 1e-14,
                "Balanced actual convergence ratio must be near2, not MAC1.5 or center ratio3.");
            require(std::abs(both.levels[0].policySeparation - flux) < 1e-14 &&
                        pressure.levels[0].maximumActiveRadius == .1 &&
                        both.levels[0].maximumActiveRadius == .12,
                    "Inactive wall flux altered a pressure plan, or a requested zero-flow face was "
                    "pruned.");
            for (unsigned l = 1; l < 3; ++l) {
                require(pressure.levels[l].policySeparation <=
                                pressure.levels[l - 1].policySeparation &&
                            pressure.levels[l].expansion.activeRadius <=
                                pressure.levels[l - 1].expansion.activeRadius,
                        "Inherited convergence geometry/radius envelope is not monotone.");
                require(
                    !(pressure.levels[l - 1].expansion.representation ==
                          LevelRepresentation::regularSpherical &&
                      pressure.levels[l].expansion.representation ==
                          LevelRepresentation::planeWave),
                    "Actual convergence policy introduced an unsupported LF-parent/HF-child edge.");
            }
        }
    // An independently rooted fine-level far subtree can have a larger raw
    // radius than the active coarse branch. Envelope only active levels.
    auto disjoint = tree;
    disjoint.nodes.resize(8);
    for (unsigned body = 0; body < 2; ++body) {
        auto &node = disjoint.nodes[6 + body];
        node.level = 2;
        node.center = {1. + .6 * body, 0, 0};
        node.pressureRadius = .15;
        node.pressureCount = 1;
        disjoint.levels[2].nodes.push_back(6 + body);
    }
    disjoint.pressurePressure.farInteractions.push_back({7, 6});
    const auto envelope = makeSupportExpansionPlan(disjoint, absent, 100, 5e-6);
    require(envelope.levels[0].maximumActiveRadius == .1 &&
                envelope.levels[2].maximumActiveRadius == .15 &&
                envelope.levels[0].expansion.activeRadius == .15 &&
                envelope.levels[1].expansion.activeRadius == .15,
            "Disjoint active fine roots escaped the conservative common-level radius envelope.");
    for (unsigned l = 1; l < 3; ++l)
        require(!(envelope.levels[l - 1].expansion.representation ==
                      LevelRepresentation::regularSpherical &&
                  envelope.levels[l].expansion.representation == LevelRepresentation::planeWave),
                "Disjoint active fine roots caused an LF-parent/HF-child representation.");
    // A strongly unequal pair approaches the public MAC, rather than the
    // balanced ratio2. No universal doubling of the MAC is admissible.
    tree.levels.resize(1);
    tree.nodes.resize(2);
    tree.fluxLeaf = {0};
    tree.nodes[0].pressureRadius = .02;
    tree.nodes[1].pressureRadius = .08;
    tree.nodes[1].center.x = .150001;
    tree.pressureFlux.farInteractions.clear();
    const auto asymmetric = makeSupportExpansionPlan(tree, absent, 10, 5e-6);
    require(std::abs(asymmetric.levels[0].policySeparation - (.150001 - .02) / .08) < 1e-14,
            "Unequal support radii lost the limiting target convergence ratio.");
    std::stop_source cancellation;
    cancellation.request_stop();
    bool stopped = false;
    try {
        (void)makeSupportExpansionPlan(tree, absent, 10, 5e-6, cancellation.get_token());
    } catch (const std::runtime_error &) {
        stopped = true;
    }
    require(stopped, "Support level planning ignored cancellation.");
    std::cout << "Actual source/target convergence geometry, active-face pruning and monotone "
                 "LF/HF policy passed.\n";
}
} // namespace
int main() {
    try {
        completeSupportIdentity();
        angularTransforms();
        regularRadialReference();
        diagonalGreen();
        automaticNearMacKernel();
        fullSupportGalerkin();
        farBlockOracleConvergence();
        scaledLowFrequencyTranslations();
        scaledAxialProjectionReference();
        scaledTransitionAndZeros();
        scaledLowFrequencyGalerkin();
        separateProbeNearMacGalerkin();
        asymmetricNearMacProbe();
        actualConvergenceGeometry();
        automaticNearMacKernel(true);
        separateProbeNearMacGalerkin(true);
        asymmetricNearMacProbe(true);
        std::cout << "Full-support IBP, open-panel counterexample, outgoing diagonal Green "
                     "derivatives and scalar Burton-Miller/direct Galerkin references passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
