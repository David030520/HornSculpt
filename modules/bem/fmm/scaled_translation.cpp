#include "scaled_translation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace hornsculpt::wideband {
namespace {
constexpr Complex imaginary{0, 1};
std::size_t index(unsigned l, int m) {
    return static_cast<std::size_t>(l) * l + l + m;
}
std::size_t count(unsigned order) {
    if (order > maximumScaledSphericalOrder)
        throw std::invalid_argument(
            "Scaled spherical reference order exceeds its bounded CPU contract.");
    return static_cast<std::size_t>(order + 1) * (order + 1);
}
double norm(Vec3 v) {
    return std::hypot(v.x, v.y, v.z);
}
Vec3 scale(Vec3 v, double s) {
    return {v.x * s, v.y * s, v.z * s};
}
Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
double axis(Vec3 v, unsigned i) {
    return i == 0 ? v.x : i == 1 ? v.y : v.z;
}
void domain(double k, double radius) {
    if (!std::isfinite(k) || k < 0 || !std::isfinite(radius) || radius <= 0 ||
        k * radius > maximumScaledElectricalRadius)
        throw std::invalid_argument("Scaled spherical waves require k>=0, a>0 and supported ka.");
}
void expansion(unsigned order, double radius, std::span<const Complex> values, double k) {
    domain(k, radius);
    if (values.size() != count(order))
        throw std::invalid_argument(
            "Scaled spherical expansion has an inconsistent coefficient count.");
    for (const auto value : values)
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::invalid_argument(
                "Scaled spherical expansion contains a non-finite coefficient.");
}
ScalarWave multiply(const ScalarWave &a, const ScalarWave &b) {
    ScalarWave result;
    result.value = a.value * b.value;
    for (unsigned i = 0; i < 3; ++i)
        result.gradient[i] = a.gradient[i] * b.value + a.value * b.gradient[i];
    return result;
}
ScalarWave times(ScalarWave a, double factor) {
    a.value *= factor;
    for (auto &g : a.gradient)
        g *= factor;
    return a;
}
ScalarWave minus(ScalarWave a, const ScalarWave &b) {
    a.value -= b.value;
    for (unsigned i = 0; i < 3; ++i)
        a.gradient[i] -= b.gradient[i];
    return a;
}
// Normalized regular solid harmonics, generated as polynomials with exact
// first-derivative propagation; no polar-angle division occurs at either pole.
std::vector<ScalarWave> solidHarmonics(unsigned order, Vec3 u) {
    std::vector<ScalarWave> result(count(order));
    const ScalarWave xy{{u.x, u.y}, {Complex{1, 0}, Complex{0, 1}, Complex{}}},
        z{{u.z, 0}, {Complex{}, Complex{}, Complex{1, 0}}},
        r2{{u.x * u.x + u.y * u.y + u.z * u.z, 0},
           {Complex{2 * u.x, 0}, Complex{2 * u.y, 0}, Complex{2 * u.z, 0}}};
    ScalarWave diagonal;
    diagonal.value = 1 / std::sqrt(4 * std::numbers::pi);
    for (unsigned m = 0; m <= order; ++m) {
        if (m)
            diagonal = times(multiply(xy, diagonal), -std::sqrt((2. * m + 1) / (2. * m)));
        result[index(m, static_cast<int>(m))] = diagonal;
        if (m < order)
            result[index(m + 1, static_cast<int>(m))] =
                times(multiply(z, diagonal), std::sqrt(2. * m + 3));
        for (unsigned l = m + 2; l <= order; ++l) {
            const double denominator = static_cast<double>(l) * l - static_cast<double>(m) * m,
                         a = std::sqrt((4. * l * l - 1) / denominator),
                         b = std::sqrt((2. * l + 1) *
                                       ((l - 1.) * (l - 1.) - static_cast<double>(m) * m) /
                                       ((2. * l - 3) * denominator));
            result[index(l, static_cast<int>(m))] =
                minus(times(multiply(z, result[index(l - 1, static_cast<int>(m))]), a),
                      times(multiply(r2, result[index(l - 2, static_cast<int>(m))]), b));
        }
        if (m)
            for (unsigned l = m; l <= order; ++l) {
                auto negative = result[index(l, static_cast<int>(m))];
                const double sign = m % 2 ? -1. : 1.;
                negative.value = sign * std::conj(negative.value);
                for (auto &g : negative.gradient)
                    g = sign * std::conj(g);
                result[index(l, -static_cast<int>(m))] = negative;
            }
    }
    return result;
}
// DLMF 10.52.1 and 10.53.1: https://dlmf.nist.gov/10.53.E1
// Jhat_l(z)=(2l+1)!!*j_l(z)/z^l has a convergent small-argument
// series with leading term one. Miller recursion avoids power-series
// cancellation away from zero, including the LF/HF transition region.
std::vector<double> normalizedRegular(unsigned order, double z) {
    if (!std::isfinite(z) || z < 0 || z > 2 * maximumScaledElectricalRadius)
        throw std::invalid_argument(
            "Scaled regular reference projection exceeds its stable LF domain.");
    std::vector<double> result(order + 1);
    const auto radial = z > 1 ? regularRadialTable(order, z) : RegularRadialTable{};
    double normalization = 1;
    for (unsigned l = 0; l <= order; ++l) {
        if (l)
            normalization *= z / (2. * l + 1);
        // Once the first series correction is <=1, its alternating terms
        // decrease and cancellation is harmless. This branch also avoids a
        // tiny j_l divided by an underflowed leading power at large l/small z.
        if (z > 1 && z * z > 2 * (2. * l + 3)) {
            result[l] = radial.values[l] / normalization;
            continue;
        }
        double term = 1, sum = 1;
        for (unsigned n = 1; n < 128; ++n) {
            term *= -z * z / (2. * n * (2. * l + 2 * n + 1));
            sum += term;
            if (std::abs(term) < std::numeric_limits<double>::epsilon() * std::abs(sum))
                break;
        }
        result[l] = sum;
    }
    return result;
}
void addWeighted(ScalarWave &result, const ScalarWave &basis, Complex coefficient) {
    result.value += coefficient * basis.value;
    for (unsigned i = 0; i < 3; ++i)
        result.gradient[i] += coefficient * basis.gradient[i];
}
std::vector<Complex> outgoingRadial(unsigned order, double k, double radius, double distance) {
    const double t = radius / distance, ka = k * radius;
    std::vector<Complex> result(order + 2);
    result[0] = t * std::exp(-imaginary * (k * distance));
    result[1] = t * (1.0 + imaginary * (k * distance)) * result[0];
    for (unsigned l = 1; l <= order; ++l)
        result[l + 1] = t * result[l] - ka * ka / ((2. * l - 1) * (2. * l + 1)) * result[l - 1];
    return result;
}
} // namespace

ScaledWaveBasis scaledRegularWaves(unsigned order, double k, double radius, Vec3 point) {
    domain(k, radius);
    const auto solid = solidHarmonics(order, scale(point, 1 / radius));
    const auto radial = normalizedRegular(order + 1, k * norm(point));
    ScaledWaveBasis result{std::vector<ScalarWave>(solid.size())};
    for (unsigned l = 0; l <= order; ++l)
        for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m) {
            const auto j = index(l, m);
            result.waves[j].value = radial[l] * solid[j].value;
            for (unsigned i = 0; i < 3; ++i)
                result.waves[j].gradient[i] =
                    radial[l] * solid[j].gradient[i] / radius -
                    (k * k * axis(point, i) / (2. * l + 3)) * radial[l + 1] * solid[j].value;
        }
    return result;
}

ScaledWaveBasis scaledOutgoingWaves(unsigned order, double k, double radius, Vec3 point) {
    domain(k, radius);
    const double distance = norm(point);
    if (!std::isfinite(distance) || distance <= 0)
        throw std::invalid_argument(
            "Scaled outgoing basis requires a finite point outside its singular center.");
    const Vec3 unit = scale(point, 1 / distance);
    const auto solid = solidHarmonics(order, unit);
    // Rescale DLMF 10.51.1-2 with the limit in 10.52.2:
    // https://dlmf.nist.gov/10.51.E1 ; https://dlmf.nist.gov/10.52.E2
    // O_l=(a/r)^(l+1)*Hhat_l(kr). This recurrence never forms
    // h_l(kr)~(kr)^(-l-1), nor a product of huge and tiny k-powers.
    const auto outgoing = outgoingRadial(order, k, radius, distance);
    ScaledWaveBasis result{std::vector<ScalarWave>(count(order))};
    for (unsigned l = 0; l <= order; ++l) {
        const Complex derivative =
            (l / distance) * outgoing[l] - (2. * l + 1) / radius * outgoing[l + 1];
        const double normalization = 1 / (radius * (2. * l + 1));
        for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m) {
            const auto j = index(l, m);
            ScalarWave wave;
            wave.value = normalization * outgoing[l] * solid[j].value;
            for (unsigned i = 0; i < 3; ++i)
                wave.gradient[i] =
                    normalization *
                    (derivative * axis(unit, i) * solid[j].value +
                     outgoing[l] / distance *
                         (solid[j].gradient[i] - l * axis(unit, i) * solid[j].value));
            result.waves[j] = wave;
        }
    }
    return result;
}

ScalarWave evaluateScaledMultipole(const ScaledMultipole &source, double k, Vec3 point) {
    expansion(source.order, source.radius, source.coefficients, k);
    const auto basis = scaledOutgoingWaves(source.order, k, source.radius, point);
    ScalarWave result;
    for (std::size_t j = 0; j < basis.waves.size(); ++j)
        addWeighted(result, basis.waves[j], source.coefficients[j]);
    return result;
}

ScalarWave evaluateScaledLocal(const ScaledLocal &local, double k, Vec3 point) {
    expansion(local.order, local.radius, local.coefficients, k);
    const auto basis = scaledRegularWaves(local.order, k, local.radius, point);
    ScalarWave result;
    for (std::size_t j = 0; j < basis.waves.size(); ++j)
        addWeighted(result, basis.waves[j], local.coefficients[j]);
    return result;
}

std::vector<Complex> scaledMultipoleToSpectrum(const ScaledMultipole &source, double k,
                                               const SphereRule &rule) {
    expansion(source.order, source.radius, source.coefficients, k);
    auto coefficients = source.coefficients;
    double radialScale = 1;
    for (unsigned l = 0; l <= source.order; ++l) {
        if (l)
            radialScale *= k * source.radius / (2. * l + 1);
        for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m)
            coefficients[index(l, m)] *= radialScale;
    }
    return multipoleToSpectrum(coefficients, source.order, rule);
}

ScaledLocal spectrumToScaledLocal(std::span<const Complex> incoming, const SphereRule &rule,
                                  unsigned order, double k, double radius) {
    domain(k, radius);
    count(order);
    ScaledLocal result{order, radius, spectrumToLocal(incoming, rule, order)};
    double radialScale = 1;
    for (unsigned l = 0; l <= order; ++l) {
        if (l)
            radialScale *= k * radius / (2. * l + 1);
        for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m)
            result.coefficients[index(l, m)] *= radialScale;
    }
    return result;
}

ScaledAxialPlan makeScaledAxialPlan(ScaledAxialKind kind, unsigned sourceOrder,
                                    unsigned targetOrder, double k, double sourceRadius,
                                    double targetRadius, double distance, unsigned polarCount,
                                    double projectionRadius, std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Scaled axial plan construction cancelled.");
    domain(k, sourceRadius);
    domain(k, targetRadius);
    count(sourceOrder);
    count(targetOrder);
    if (kind != ScaledAxialKind::multipoleToParent && kind != ScaledAxialKind::multipoleToLocal &&
        kind != ScaledAxialKind::localToChild)
        throw std::invalid_argument("Scaled axial translation kind is invalid.");
    if (!std::isfinite(distance) || distance < 0 ||
        polarCount < 2 * std::max(sourceOrder, targetOrder) + 1 || polarCount > 1024)
        throw std::invalid_argument("Scaled axial plan needs finite separation and independently "
                                    "resolved polar quadrature.");
    if (!std::isfinite(projectionRadius) || projectionRadius < 0)
        throw std::invalid_argument(
            "Scaled axial projection radius must be finite and nonnegative.");
    const double probe = projectionRadius > 0                         ? projectionRadius
                         : kind == ScaledAxialKind::multipoleToParent ? sourceRadius * 1.01
                         : kind == ScaledAxialKind::multipoleToLocal
                             ? std::min(targetRadius, distance * .5)
                             : targetRadius;
    if (!(probe > 0) || (kind == ScaledAxialKind::multipoleToLocal && probe >= distance))
        throw std::invalid_argument(
            "Scaled axial Cauchy probe intersects the outgoing singular center.");
    ScaledAxialPlan result{kind, sourceOrder, targetOrder, polarCount, {0}, {}};
    const unsigned maximumM = std::min(sourceOrder, targetOrder);
    for (unsigned m = 0; m <= maximumM; ++m)
        result.blockOffsets.push_back(result.blockOffsets.back() +
                                      static_cast<std::size_t>(targetOrder - m + 1) *
                                          (sourceOrder - m + 1));
    result.coefficients.resize(result.blockOffsets.back());
    result.projectionRadius = probe;
    const auto sphere = makeSphereRule(polarCount - 1);
    const auto targetRadial = normalizedRegular(targetOrder + 1, k * probe);
    for (unsigned polar = 0; polar < polarCount; ++polar) {
        if (stop.stop_requested())
            throw std::runtime_error("Scaled axial plan construction cancelled.");
        const auto &sample = sphere.samples[static_cast<std::size_t>(polar) * 2 * polarCount];
        const Vec3 direction = sample.direction;
        // The phi integral is analytic. makeSphereRule's phi=0 weight times
        // its azimuth count already equals 2*pi times the polar Gauss weight.
        const double weight = sample.weight * (2 * polarCount);
        const double traceRadius = probe;
        const Vec3 point = add({0, 0, distance}, scale(direction, traceRadius));
        const auto angular = sphericalHarmonics(std::max(sourceOrder, targetOrder), direction);
        ScaledWaveBasis sourceBasis, targetBasis;
        if (kind == ScaledAxialKind::multipoleToParent) {
            sourceBasis =
                scaledOutgoingWaves(sourceOrder, k, sourceRadius, scale(direction, traceRadius));
            targetBasis = scaledRegularWaves(targetOrder, k, targetRadius, point);
        } else if (kind == ScaledAxialKind::multipoleToLocal) {
            sourceBasis = scaledOutgoingWaves(sourceOrder, k, sourceRadius, point);
        } else {
            sourceBasis = scaledRegularWaves(sourceOrder, k, sourceRadius, point);
        }
        for (unsigned m = 0; m <= maximumM; ++m) {
            if (stop.stop_requested())
                throw std::runtime_error("Scaled axial plan construction cancelled.");
            for (unsigned l = m; l <= targetOrder; ++l) {
                const auto row = index(l, static_cast<int>(m));
                const double radialDerivative = l * targetRadial[l] - k * k * probe * probe /
                                                                          (2. * l + 3) *
                                                                          targetRadial[l + 1],
                             denominator = targetRadial[l] * targetRadial[l] +
                                           radialDerivative * radialDerivative;
                for (unsigned n = m; n <= sourceOrder; ++n) {
                    const auto column = index(n, static_cast<int>(m));
                    Complex sourceNormal{};
                    for (unsigned component = 0; component < 3; ++component)
                        sourceNormal += sourceBasis.waves[column].gradient[component] *
                                        axis(direction, component);
                    Complex value{};
                    if (kind == ScaledAxialKind::multipoleToParent) {
                        Complex targetNormal{};
                        for (unsigned component = 0; component < 3; ++component)
                            targetNormal += targetBasis.waves[row].gradient[component] *
                                            axis(direction, component);
                        value = traceRadius * traceRadius *
                                (sourceBasis.waves[column].value * std::conj(targetNormal) -
                                 sourceNormal * std::conj(targetBasis.waves[row].value));
                    } else {
                        if (!(denominator > 0) || !std::isfinite(denominator))
                            throw std::runtime_error("Scaled axial Cauchy projection has an "
                                                     "unresolved radial normalization.");
                        value = std::pow(targetRadius / probe, l) * std::conj(angular[row]) *
                                (targetRadial[l] * sourceBasis.waves[column].value +
                                 radialDerivative * probe * sourceNormal) /
                                denominator;
                    }
                    const auto entry = result.blockOffsets[m] +
                                       static_cast<std::size_t>(l - m) * (sourceOrder - m + 1) + n -
                                       m;
                    result.coefficients[entry] += weight * value;
                }
            }
        }
    }
    for (const auto value : result.coefficients)
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::runtime_error("Scaled axial plan has a non-finite translation coefficient.");
    return result;
}

std::vector<Complex> applyScaledAxialPlan(const ScaledAxialPlan &plan,
                                          std::span<const Complex> source) {
    if (source.size() != count(plan.sourceOrder) ||
        plan.blockOffsets.size() != std::min(plan.sourceOrder, plan.targetOrder) + 2 ||
        plan.blockOffsets.front() != 0 || plan.blockOffsets.back() != plan.coefficients.size())
        throw std::invalid_argument("Scaled axial application has inconsistent dimensions.");
    std::vector<Complex> result(count(plan.targetOrder));
    for (const auto value : source)
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::invalid_argument("Scaled axial source contains a non-finite coefficient.");
    for (unsigned m = 0; m <= std::min(plan.sourceOrder, plan.targetOrder); ++m) {
        const auto expected =
            static_cast<std::size_t>(plan.targetOrder - m + 1) * (plan.sourceOrder - m + 1);
        if (plan.blockOffsets[m] > plan.blockOffsets[m + 1] ||
            plan.blockOffsets[m + 1] > plan.coefficients.size() ||
            plan.blockOffsets[m + 1] - plan.blockOffsets[m] != expected)
            throw std::invalid_argument("Scaled axial application has an invalid m-block span.");
        for (const int signedM : {static_cast<int>(m), -static_cast<int>(m)}) {
            if (m == 0 && signedM != 0)
                continue;
            for (unsigned l = m; l <= plan.targetOrder; ++l) {
                Complex sum{};
                for (unsigned n = m; n <= plan.sourceOrder; ++n) {
                    const auto entry =
                        plan.blockOffsets[m] +
                        static_cast<std::size_t>(l - m) * (plan.sourceOrder - m + 1) + n - m;
                    sum += plan.coefficients[entry] * source[index(n, signedM)];
                }
                result[index(l, signedM)] = sum;
            }
            if (m == 0)
                break;
        }
    }
    for (const auto value : result)
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::runtime_error("Scaled axial application produced a non-finite coefficient.");
    return result;
}

} // namespace hornsculpt::wideband
