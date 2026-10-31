#include "ring_kernel_reference.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <galerkin/axisymmetric_quadrature.hpp>
#include <galerkin/pair_quadrature.hpp>
#include <iostream>
#include <numbers>
#include <stdexcept>

using namespace hornsculpt;
using namespace hornsculpt::galerkin;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
std::complex<double> complex(RingComplex x) {
    return {x.real, x.imag};
}
RingGeometry geometry(double r, double t, double z, double nr, double nz, double mr, double mz) {
    return {r, t, r - t, z, nr, nz, mr, mz, nr * (r - t) + nz * z, mr * (t - r) - mz * z};
}
void ellipticTests() {
    for (double q : {1.0, .9, .5, .1, .001, 1e-6}) {
        const auto actual = ring_detail::completeElliptic(q);
        const double modulus = std::sqrt(1 - q);
        require(!actual.status &&
                    std::abs(actual.first / std::comp_ellint_1(modulus) - 1) < 8e-12 &&
                    std::abs(actual.second / std::comp_ellint_2(modulus) - 1) < 8e-12,
                "Complementary elliptic evaluation differs from independent standard functions.");
    }
    const auto tiny = ring_detail::completeElliptic(1e-40);
    require(!tiny.status && std::abs(tiny.first - (std::log(4.) - .5 * std::log(1e-40))) < 2e-14 &&
                tiny.second == 1,
            "Complementary tiny-q limit lost precision.");
}
void logRuleTests() {
    for (unsigned n : {4u, 6u, 8u, 12u, 16u, 20u}) {
        const auto rule = ringRule(n);
        double harmonic = 0;
        for (unsigned power = 0; power < n; ++power) {
            harmonic += 1.0 / (power + 1);
            double first = 0, last = 0, regular = 0;
            for (const auto &p : rule) {
                const double value = std::pow(p.position, power);
                first += p.logFirst * value;
                last += p.logLast * value;
                regular += p.weight * value;
            }
            require(std::abs(first + 1.0 / ((power + 1.) * (power + 1.))) < 2e-13 &&
                        std::abs(last + harmonic / (power + 1.)) < 2e-13 &&
                        std::abs(regular - 1.0 / (power + 1.)) < 2e-13,
                    "Meridian logarithmic product weights fail analytic polynomial moments.");
        }
    }
    for (unsigned n : {8u, 16u}) {
        const auto rule = gaussLegendre01(n);
        for (unsigned i = 0; i < n; ++i) {
            const auto stored = ring_detail::gauss(n, i);
            require(std::abs(stored.position - rule[i].position) < 2e-14 &&
                        std::abs(stored.weight - rule[i].weight) < 2e-14,
                    "Shared angular Gauss literals differ from generated nodes.");
        }
    }
}
void kernelTests() {
    const std::vector<RingGeometry> cases{
        geometry(.04, .041, .002, .6, .8, -.8, .6),
        geometry(.085, .085001, -.00000677, .011, -std::sqrt(1 - .011 * .011), 0, 1),
        geometry(.025, .0250000000025, .000000000025, 1, 0, 1, 0),
        geometry(0, .025, .01, 0, 1, 0, 1),
        geometry(1e-12, .025, .01, 0, 1, 0, 1),
        geometry(.001, .08, -.02, .8, .6, .6, -.8)};
    double worst = 0;
    for (const auto &g : cases)
        for (double k : {0.0, 73.2733, 293.093}) {
            const auto actual = ringKernels(g, k, 1e-9);
            require(!actual.status, "Ring kernel did not attain its angular error budget.");
            const auto reference = testing::referenceRingKernel(g, k);
            const auto stricter = testing::referenceRingKernel(g, k, 32, 2e-12L);
            const double size = ring_detail::scale(g);
            for (unsigned c = 0; c < 4; ++c) {
                const double floor = c < 2 ? 1 / size : 1 / (size * size);
                const long double scale = std::max(std::abs(stricter[c]), (long double)floor);
                require(std::abs(reference[c] - stricter[c]) / scale < 2e-9L,
                        "Independent angular reference has not converged.");
                const std::complex<long double> value{actual.value[c].real, actual.value[c].imag};
                const double error = static_cast<double>(std::abs(value - stricter[c]) / scale);
                require(error < 3e-8,
                        "Analytic/dynamic ring kernel fails independent 3D integration.");
                worst = std::max(worst, error);
            }
            auto reversed = g;
            std::swap(reversed.testRadius, reversed.sourceRadius);
            std::swap(reversed.testRadialNormal, reversed.sourceRadialNormal);
            std::swap(reversed.testAxialNormal, reversed.sourceAxialNormal);
            std::swap(reversed.testNormalDistance, reversed.sourceNormalDistance);
            reversed.deltaRadius = -g.deltaRadius;
            reversed.deltaZ = -g.deltaZ;
            const auto other = ringKernels(reversed, k, 1e-9);
            require(!other.status, "Reversed ring kernel failed.");
            for (unsigned c = 0; c < 4; ++c) {
                const unsigned match = c < 2 ? c : 5 - c;
                require(std::abs(complex(actual.value[c]) - complex(other.value[match])) <=
                            2e-13 * std::max(1.0, ringMagnitude(actual.value[c])),
                        "Ring kernel reciprocity or derivative orientation changed.");
            }
            if (g.testRadius * g.sourceRadius > 0) {
                const auto split = ringStaticSplit(g);
                const auto statics = ringKernels(g, 0);
                require(!split.status && !statics.status, "Static logarithmic split failed.");
                for (unsigned c = 0; c < 4; ++c) {
                    const double reconstructed =
                        split.regular[c] + split.logarithmic[c] * split.logFactor;
                    require(std::abs(reconstructed - statics.value[c].real) <=
                                3e-14 * std::max(1.0, std::abs(statics.value[c].real)),
                            "Static logarithmic subtraction changed the kernel.");
                }
            }
        }
    std::cout << "Ring kernels versus independent clustered3D-psi reference: " << worst << '\n';
}

AxisymmetricSpace fixture(std::vector<MeridianPoint> points,
                          std::vector<std::array<unsigned, 2>> edges) {
    AxisymmetricSpace result;
    result.meridian.points = std::move(points);
    result.lumpedMass.resize(result.meridian.points.size());
    for (const auto &edge : edges) {
        const auto a = result.meridian.points[edge[0]], b = result.meridian.points[edge[1]];
        AxisymmetricElement e;
        e.first = edge[0];
        e.second = edge[1];
        e.faceId = 1;
        const double dr = b.radius - a.radius, dz = b.axial - a.axial;
        e.length = std::hypot(dr, dz);
        e.radialNormal = dz / e.length;
        e.axialNormal = -dr / e.length;
        e.area = std::numbers::pi * (a.radius + b.radius) * e.length;
        e.massP1DP0 = {std::numbers::pi * e.length * (2 * a.radius + b.radius) / 3,
                       std::numbers::pi * e.length * (a.radius + 2 * b.radius) / 3};
        result.lumpedMass[e.first] += e.massP1DP0[0];
        result.lumpedMass[e.second] += e.massP1DP0[1];
        result.elements.push_back(e);
    }
    return result;
}
std::vector<RingPairPatch> roots(const RingInteractionPlan &plan, unsigned a, unsigned b) {
    const auto range = plan.pairs[a * plan.segments.size() + b];
    return {plan.roots.begin() + range.first, plan.roots.begin() + range.first + range.count};
}
void pairTests() {
    const std::vector<AxisymmetricSpace> cases{
        fixture({{.025, -.006}, {.025, 0}}, {{0, 1}}), fixture({{.025, 0}, {0, 0}}, {{0, 1}}),
        fixture({{0, -.03}, {.006, -.028}}, {{0, 1}}),
        fixture({{.08518012928817169, -.000012165773387238084},
                 {.0862725888227, 0},
                 {.08492457962234534, 0}},
                {{0, 1}, {1, 2}})};
    BemNumerics settings;
    settings.singularOrder = 8;
    settings.nearRelativeTolerance = 1e-4;
    settings.nearMaximumDepth = 8;
    for (unsigned c = 0; c < cases.size(); ++c) {
        const auto plan = makeRingInteractionPlan(cases[c], settings);
        const unsigned source = c == 3 ? 1 : 0;
        const auto domains = roots(plan, 0, source);
        const auto actual = integrateRingPair(plan.segments[0], plan.segments[source], domains,
                                              plan.coarse, plan.fine, 293.093, 1e-4, 8);
        std::cout << "Ring pair case " << c << " status=" << actual.status
                  << " patches=" << actual.patches << " depth=" << actual.maximumDepth
                  << " estimate=" << actual.estimatedRelativeError << std::endl;
        require(!actual.status && actual.estimatedRelativeError < 1e-4,
                "Controlled self/axis/thin-wedge ring-pair integration failed.");
        const auto high = integrateRingPair(plan.segments[0], plan.segments[source], domains,
                                            ringRule(12), ringRule(16), 293.093, 2e-6, 10);
        std::cout << "Stricter ring pair case " << c << " status=" << high.status
                  << " patches=" << high.patches << " depth=" << high.maximumDepth
                  << " estimate=" << high.estimatedRelativeError << std::endl;
        require(!high.status, "Stricter ring-pair reference failed.");
        const auto difference = ringSubtract(actual.moments, high.moments);
        const double relative =
            ringOperatorNorm(difference, plan.segments[0], plan.segments[source], 293.093) /
            ringOperatorNorm(high.moments, plan.segments[0], plan.segments[source], 293.093);
        require(relative < 1e-4, "Ring-pair operator refinement exceeds its accuracy budget.");
    }
    // Closed-cylinder Laplace representation: (1/2 M + D)1=0. This verifies
    // the complete surface factor, normal derivatives and all corner/pole cells.
    const auto cylinder =
        fixture({{0, -.006}, {.025, -.006}, {.025, 0}, {0, 0}}, {{0, 1}, {1, 2}, {2, 3}});
    settings.nearDistanceRatio = 100;
    const auto plan = makeRingInteractionPlan(cylinder, settings);
    std::vector<std::complex<double>> rows(cylinder.meridian.points.size());
    for (unsigned a = 0; a < 3; ++a)
        for (unsigned b = 0; b < 3; ++b) {
            const auto integrated =
                integrateRingPair(plan.segments[a], plan.segments[b], roots(plan, a, b),
                                  ringRule(8), ringRule(12), 0, 2e-6, 10);
            require(!integrated.status, "Closed-cylinder Laplace pair failed.");
            const unsigned nodes[]{plan.segments[a].first, plan.segments[a].second};
            for (unsigned i = 0; i < 2; ++i)
                for (unsigned j = 0; j < 2; ++j)
                    rows[nodes[i]] += complex(integrated.moments.value[12 + 2 * i + j]);
        }
    double worst = 0;
    for (unsigned i = 0; i < rows.size(); ++i)
        worst = std::max(worst, std::abs(rows[i] / cylinder.lumpedMass[i] + .5));
    std::cout << "Closed-cylinder Laplace constant-trace residual " << worst << '\n';
    require(worst < 2e-4, "Meridian Galerkin jump/normal/surface factor is inconsistent.");
}
} // namespace
int main() {
    try {
        ellipticTests();
        logRuleTests();
        kernelTests();
        pairTests();
        std::cout << "Continuous revolution kernel and pair tests passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
