#pragma once

#include "galerkin/ring_kernels.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace hornsculpt::testing {
namespace ring_reference_detail {
using Real = long double;
using Complex = std::complex<Real>;
using Values = std::array<Complex, 4>;
struct GaussPoint {
    Real position{}, weight{};
};

// Generated independently of the production quadrature tables. On MSVC,
// long double is FP64: this is an independent method, not extra precision.
inline std::vector<GaussPoint> gauss(unsigned order) {
    std::vector<GaussPoint> rule(order);
    const Real pi = std::numbers::pi_v<Real>, eps = std::numeric_limits<Real>::epsilon();
    for (unsigned i = 0; i < (order + 1) / 2; ++i) {
        Real x = std::cos(pi * (i + .75L) / (order + .5L));
        Real derivative{};
        bool converged = false;
        for (unsigned iteration = 0; iteration < 64; ++iteration) {
            Real p0 = 1, p1 = x;
            for (unsigned n = 2; n <= order; ++n) {
                const Real pn = ((2 * n - 1) * x * p1 - (n - 1) * p0) / n;
                p0 = p1;
                p1 = pn;
            }
            derivative = order * (x * p1 - p0) / (x * x - 1);
            const Real step = p1 / derivative;
            x -= step;
            if (std::abs(step) <= 8 * eps * std::max(Real(1), std::abs(x))) {
                converged = true;
                break;
            }
        }
        if (!converged)
            throw std::runtime_error("Independent Gauss root did not converge.");
        // Re-evaluate the derivative at the accepted root.
        Real p0 = 1, p1 = x;
        for (unsigned n = 2; n <= order; ++n) {
            const Real pn = ((2 * n - 1) * x * p1 - (n - 1) * p0) / n;
            p0 = p1;
            p1 = pn;
        }
        derivative = order * (x * p1 - p0) / (x * x - 1);
        const Real weight = 1 / ((1 - x * x) * derivative * derivative);
        rule[i] = {(1 - x) / 2, weight};
        rule[order - 1 - i] = {(1 + x) / 2, weight};
    }
    return rule;
}
struct Phase {
    Complex potential, gradient;
};
inline Phase phase(Real x) {
    if (std::abs(x) >= .5L) {
        const Complex exponential{std::cos(x), -std::sin(x)};
        return {exponential, Complex{-1, -x} * exponential};
    }
    // d(exp(-ikR)/R)/dR has numerator (-1-ikR)exp(-ikR).
    // Summing (n-1)(-ix)^n/n! directly retains the cubic imaginary term;
    // forming sin(x)-x*cos(x) would lose it near the static limit.
    Complex term{1, 0}, exponential{1, 0}, gradient{-1, 0};
    const Real eps = std::numeric_limits<Real>::epsilon();
    for (unsigned n = 1; n <= 128; ++n) {
        term *= Complex{0, -x} / Real(n);
        exponential += term;
        gradient += Real(n - 1) * term;
        if (n >= 8 && std::abs(term) <= eps * eps)
            return {exponential, gradient};
    }
    throw std::runtime_error("Independent direct Green series did not converge.");
}
struct Integrand {
    const galerkin::RingGeometry &geometry;
    Real k{}, width{}, extent{};
    bool sinhCluster{};
    Values operator()(Real coordinate) const {
        const Real pi = std::numbers::pi_v<Real>;
        const Real psi = sinhCluster ? width * std::sinh(extent * coordinate) : pi * coordinate;
        const Real jacobian = sinhCluster ? width * extent * std::cosh(extent * coordinate) : pi;
        const Real sh = std::sin(psi / 2), oneMinusCos = 2 * sh * sh;
        const Real r = geometry.testRadius, t = geometry.sourceRadius;
        const Real delta = std::hypot(Real(geometry.deltaRadius), Real(geometry.deltaZ));
        const Real distance = std::hypot(delta, 2 * std::sqrt(r) * std::sqrt(t) * sh);
        const auto wave = phase(k * distance);
        // These are explicit stable geometric dot products, not derivatives
        // of an elliptic-integral expression. They also respect an exact
        // same-segment normal distance supplied by the geometry contract.
        const Real nx =
            Real(geometry.testNormalDistance) + Real(geometry.testRadialNormal) * t * oneMinusCos;
        const Real ny = Real(geometry.sourceNormalDistance) +
                        Real(geometry.sourceRadialNormal) * r * oneMinusCos;
        const Real weight = jacobian / (2 * pi), inverse = 1 / distance;
        const Complex green = wave.potential * inverse;
        const Complex gradient = wave.gradient * inverse * inverse * inverse;
        return {green * weight, green * std::cos(psi) * weight, gradient * nx * weight,
                gradient * ny * weight};
    }
};
inline Values integrate(const Integrand &integrand, const std::vector<GaussPoint> &rule, Real first,
                        Real last) {
    Values sum{}, compensation{};
    for (const auto p : rule) {
        const auto value = integrand(first + (last - first) * p.position);
        for (unsigned c = 0; c < 4; ++c) {
            const auto term = value[c] * ((last - first) * p.weight) - compensation[c];
            const auto next = sum[c] + term;
            compensation[c] = (next - sum[c]) - term;
            sum[c] = next;
        }
    }
    return sum;
}
} // namespace ring_reference_detail

// Noncoincident point-ring reference by direct three-dimensional Green and
// gradient integration. No K/E, static split, production Gauss table or angular
// adaptation is reused. A failed convergence budget throws, never relaxes.
inline std::array<std::complex<long double>, 4>
referenceRingKernel(const galerkin::RingGeometry &g, double k, unsigned order = 24,
                    long double relativeTolerance = 2e-11L) {
    using namespace ring_reference_detail;
    const Real delta = std::hypot(Real(g.deltaRadius), Real(g.deltaZ));
    if (!(g.testRadius >= 0 && g.sourceRadius >= 0 && delta > 0 && k >= 0 &&
          relativeTolerance > 0) ||
        !std::isfinite(delta) || !std::isfinite(k) || !std::isfinite(relativeTolerance) ||
        order < 4 || order > 128)
        throw std::invalid_argument("Invalid independent ring reference request.");
    Integrand integrand{g, Real(k)};
    const Real radialScale = std::sqrt(Real(g.testRadius)) * std::sqrt(Real(g.sourceRadius));
    if (radialScale > 0) {
        integrand.width = delta / radialScale;
        integrand.sinhCluster = integrand.width < std::numbers::pi_v<Real>;
        if (integrand.sinhCluster)
            integrand.extent = std::asinh(std::numbers::pi_v<Real> / integrand.width);
    }
    const auto lowRule = gauss(order), highRule = gauss(2 * order);
    const auto pilot = integrate(integrand, highRule, 0, 1);
    const Real size = std::hypot(delta, 2 * radialScale);
    std::array<Real, 4> scales{};
    for (unsigned c = 0; c < 4; ++c)
        scales[c] = std::max(std::abs(pilot[c]), c < 2 ? 1 / size : 1 / (size * size));
    struct Interval {
        Real first{}, last{};
        unsigned depth{};
    };
    std::vector<Interval> pending{{0, 1, 0}};
    Values result{}, compensation{};
    std::size_t evaluations = highRule.size();
    while (!pending.empty()) {
        const auto interval = pending.back();
        pending.pop_back();
        const auto low = integrate(integrand, lowRule, interval.first, interval.last),
                   high = integrate(integrand, highRule, interval.first, interval.last);
        evaluations += lowRule.size() + highRule.size();
        bool accepted = true;
        for (unsigned c = 0; c < 4; ++c) {
            if (!std::isfinite(high[c].real()) || !std::isfinite(high[c].imag()))
                throw std::runtime_error("Independent ring reference is non-finite.");
            accepted &= std::abs(high[c] - low[c]) <=
                        relativeTolerance * scales[c] * (interval.last - interval.first);
        }
        if (accepted) {
            for (unsigned c = 0; c < 4; ++c) {
                const auto term = high[c] - compensation[c], next = result[c] + term;
                compensation[c] = (next - result[c]) - term;
                result[c] = next;
            }
        } else {
            if (interval.depth >= 28 || evaluations > 4000000)
                throw std::runtime_error(
                    "Independent direct-angle reference did not meet its error budget.");
            const Real mid = (interval.first + interval.last) / 2;
            pending.push_back({mid, interval.last, interval.depth + 1});
            pending.push_back({interval.first, mid, interval.depth + 1});
        }
    }
    return result;
}
} // namespace hornsculpt::testing
