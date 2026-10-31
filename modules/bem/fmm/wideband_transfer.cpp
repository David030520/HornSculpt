#include "wideband_transfer.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace hornsculpt::wideband {
namespace {
constexpr double pi = std::numbers::pi;
constexpr unsigned maximumOrder = 1024;
void validateOrder(unsigned order) {
    if (order > maximumOrder)
        throw std::invalid_argument("Wideband angular order exceeds the supported CPU contract.");
}
std::size_t coefficientCount(unsigned order) {
    validateOrder(order);
    return static_cast<std::size_t>(order + 1) * (order + 1);
}
std::size_t index(unsigned l, int m) {
    return static_cast<std::size_t>(l) * l + l + m;
}
bool finite(Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
void validateRule(const SphereRule &rule) {
    validateOrder(rule.order);
    const auto n = static_cast<std::size_t>(rule.order + 1);
    if (rule.samples.size() != 2 * n * n)
        throw std::invalid_argument("Wideband sphere rule has an inconsistent sample count.");
    for (const auto &sample : rule.samples)
        if (!std::isfinite(sample.weight) || sample.weight <= 0 ||
            !std::isfinite(dot(sample.direction, sample.direction)) ||
            std::abs(dot(sample.direction, sample.direction) - 1) > 1e-12)
            throw std::invalid_argument("Wideband sphere rule has an invalid direction or weight.");
}
void validateValues(std::span<const Complex> values, const SphereRule &rule) {
    validateRule(rule);
    if (values.size() != rule.samples.size())
        throw std::invalid_argument("Wideband spectrum does not match its angular rule.");
    for (const auto value : values)
        if (!finite(value))
            throw std::invalid_argument("Wideband spectrum contains a non-finite value.");
}
Complex phase(double argument) {
    if (!std::isfinite(argument))
        throw std::invalid_argument("Wideband translation phase is non-finite.");
    return {std::cos(argument), std::sin(argument)};
}
std::vector<Complex> project(std::span<const Complex> values, const SphereRule &rule,
                             unsigned order) {
    validateValues(values, rule);
    if (order > rule.order)
        throw std::invalid_argument("Wideband projection exceeds its angular sampling order.");
    std::vector<Complex> coefficients(coefficientCount(order));
    for (std::size_t q = 0; q < values.size(); ++q) {
        const auto harmonics = sphericalHarmonics(order, rule.samples[q].direction);
        const auto weighted = values[q] * rule.samples[q].weight;
        for (std::size_t j = 0; j < coefficients.size(); ++j)
            coefficients[j] += weighted * std::conj(harmonics[j]);
    }
    return coefficients;
}
std::vector<Complex> evaluate(std::span<const Complex> coefficients, unsigned order,
                              const SphereRule &rule) {
    validateRule(rule);
    if (coefficients.size() != coefficientCount(order))
        throw std::invalid_argument("Wideband harmonic coefficients have an invalid size.");
    std::vector<Complex> result(rule.samples.size());
    for (std::size_t q = 0; q < result.size(); ++q) {
        const auto harmonics = sphericalHarmonics(order, rule.samples[q].direction);
        for (std::size_t j = 0; j < coefficients.size(); ++j)
            result[q] += coefficients[j] * harmonics[j];
    }
    return result;
}
std::size_t triangularIndex(unsigned l, unsigned m) {
    return static_cast<std::size_t>(l) * (l + 1) / 2 + m;
}
std::vector<double> associatedLegendre(const AngularTransferPlan &plan, double mu) {
    const unsigned order = plan.bandlimit;
    std::vector<double> result(plan.recurrenceFirst.size());
    const double radial = std::sqrt(std::max(0., 1 - mu * mu));
    double diagonal = 1 / std::sqrt(4 * pi);
    for (unsigned m = 0; m <= order; ++m) {
        if (m)
            diagonal *= plan.diagonalScale[m] * radial;
        result[triangularIndex(m, m)] = diagonal;
        if (m < order)
            result[triangularIndex(m + 1, m)] = std::sqrt(2. * m + 3) * mu * diagonal;
        for (unsigned l = m + 2; l <= order; ++l) {
            const auto entry = triangularIndex(l, m);
            result[entry] = plan.recurrenceFirst[entry] * mu * result[triangularIndex(l - 1, m)] -
                            plan.recurrenceSecond[entry] * result[triangularIndex(l - 2, m)];
        }
    }
    return result;
}
std::vector<PolarSample> polarSamples(const SphereRule &rule) {
    validateRule(rule);
    const unsigned azimuthCount = 2 * (rule.order + 1);
    std::vector<PolarSample> result;
    for (unsigned polar = 0; polar <= rule.order; ++polar) {
        const auto &first = rule.samples[static_cast<std::size_t>(polar) * azimuthCount];
        const double mu = first.direction.z, radius = std::sqrt(std::max(0., 1 - mu * mu));
        for (unsigned azimuth = 0; azimuth < azimuthCount; ++azimuth) {
            const auto &sample =
                rule.samples[static_cast<std::size_t>(polar) * azimuthCount + azimuth];
            const double angle = 2 * pi * azimuth / azimuthCount;
            if (std::abs(sample.direction.z - mu) > 1e-13 ||
                std::abs(sample.direction.x - radius * std::cos(angle)) > 1e-13 ||
                std::abs(sample.direction.y - radius * std::sin(angle)) > 1e-13 ||
                std::abs(sample.weight - first.weight) > 1e-13 * first.weight)
                throw std::invalid_argument(
                    "Wideband transfer requires a polar-major canonical sphere rule.");
        }
        result.push_back({mu, first.weight * azimuthCount / (2 * pi)});
    }
    return result;
}
} // namespace

SphereRule makeSphereRule(unsigned order) {
    validateOrder(order);
    SphereRule result{order, {}};
    const unsigned polarCount = order + 1, azimuthCount = 2 * polarCount;
    result.samples.reserve(static_cast<std::size_t>(polarCount) * azimuthCount);
    for (unsigned j = 0; j < polarCount; ++j) {
        double mu = std::cos(pi * (j + .75) / (polarCount + .5)), derivative = 0;
        bool converged = false;
        for (unsigned iteration = 0; iteration < 100; ++iteration) {
            double previous = 0, current = 1;
            for (unsigned l = 1; l <= polarCount; ++l) {
                const double next = ((2 * l - 1) * mu * current - (l - 1) * previous) / l;
                previous = current;
                current = next;
            }
            derivative = polarCount * (mu * current - previous) / (mu * mu - 1);
            const double correction = current / derivative;
            mu -= correction;
            if (std::abs(correction) <= 4 * std::numeric_limits<double>::epsilon()) {
                converged = true;
                break;
            }
        }
        if (!converged)
            throw std::runtime_error("Wideband Gauss-Legendre rule did not converge.");
        const double weight = 2 / ((1 - mu * mu) * derivative * derivative),
                     radial = std::sqrt(std::max(0., 1 - mu * mu));
        for (unsigned m = 0; m < azimuthCount; ++m) {
            const double phi = 2 * pi * m / azimuthCount;
            result.samples.push_back({{radial * std::cos(phi), radial * std::sin(phi), mu},
                                      weight * 2 * pi / azimuthCount});
        }
    }
    return result;
}

std::vector<Complex> sphericalHarmonics(unsigned order, Vec3 direction) {
    const double norm2 = dot(direction, direction);
    if (!std::isfinite(norm2) || std::abs(norm2 - 1) > 1e-12)
        throw std::invalid_argument("Spherical harmonic direction must be a finite unit vector.");
    std::vector<Complex> result(coefficientCount(order));
    const double mu = std::clamp(direction.z, -1., 1.),
                 radial = std::hypot(direction.x, direction.y);
    const Complex azimuth =
        radial == 0 ? Complex{1, 0} : Complex{direction.x / radial, direction.y / radial};
    Complex diagonal{1 / std::sqrt(4 * pi), 0};
    for (unsigned m = 0; m <= order; ++m) {
        if (m)
            diagonal *= -std::sqrt((2. * m + 1) / (2. * m)) * radial * azimuth;
        result[index(m, static_cast<int>(m))] = diagonal;
        if (m < order)
            result[index(m + 1, static_cast<int>(m))] = std::sqrt(2. * m + 3) * mu * diagonal;
        for (unsigned l = m + 2; l <= order; ++l) {
            const double denominator = static_cast<double>(l) * l - static_cast<double>(m) * m;
            const double a = std::sqrt((4. * l * l - 1) / denominator),
                         b = std::sqrt((2. * l + 1) *
                                       ((l - 1.) * (l - 1.) - static_cast<double>(m) * m) /
                                       ((2. * l - 3) * denominator));
            result[index(l, static_cast<int>(m))] =
                a * mu * result[index(l - 1, static_cast<int>(m))] -
                b * result[index(l - 2, static_cast<int>(m))];
        }
        if (m)
            for (unsigned l = m; l <= order; ++l)
                result[index(l, -static_cast<int>(m))] =
                    (m % 2 ? -1. : 1.) * std::conj(result[index(l, static_cast<int>(m))]);
    }
    return result;
}

RegularRadialTable regularRadialTable(unsigned order, double x) {
    validateOrder(order);
    if (!std::isfinite(x) || x < 0)
        throw std::invalid_argument(
            "Regular spherical waves require a finite nonnegative argument.");
    RegularRadialTable result;
    result.values.resize(order + 1);
    result.derivatives.resize(order + 1);
    if (x == 0) {
        result.values[0] = 1;
        if (order)
            result.derivatives[1] = 1. / 3;
        return result;
    }
    std::vector<double> radial(order + 2);
    if (x <= 1) {
        double leading = 1;
        for (unsigned l = 0; l <= order + 1; ++l) {
            const double leadingOverX = l ? leading / (2. * l + 1) : 0;
            if (l)
                leading *= x / (2. * l + 1);
            double term = 1, sum = 1, weighted = 0;
            for (unsigned n = 1; n < 128; ++n) {
                term *= -x * x / (2. * n * (2. * l + 2 * n + 1));
                sum += term;
                weighted += 2 * n * term;
                if (std::abs(term) <= std::numeric_limits<double>::epsilon() * std::abs(sum))
                    break;
            }
            radial[l] = leading * sum;
            if (l && l <= order)
                result.derivatives[l] = leadingOverX * (l * sum + weighted);
        }
    } else {
        const double j0 = std::sin(x) / x, j1 = (j0 - std::cos(x)) / x;
        if (x > order + 32.) {
            radial[0] = j0;
            radial[1] = j1;
            for (unsigned l = 1; l <= order; ++l)
                radial[l + 1] = ((2. * l + 1) / x) * radial[l] - radial[l - 1];
        } else {
            // Begin sufficiently past the turning region, then normalize by the
            // larger of j0,j1 so roots of either normalization mode are harmless.
            const auto top =
                std::max(order + 48u, static_cast<unsigned>(std::ceil(x + 48 + 8 * std::cbrt(x))));
            std::vector<double> miller(top + 2);
            miller[top] = 1;
            for (unsigned l = top; l > 0; --l) {
                miller[l - 1] = ((2. * l + 1) / x) * miller[l] - miller[l + 1];
                if (std::abs(miller[l - 1]) > 1e200)
                    for (std::size_t j = l - 1; j < miller.size(); ++j)
                        miller[j] *= 1e-200;
            }
            const double normalization =
                std::abs(j0) > std::abs(j1) ? j0 / miller[0] : j1 / miller[1];
            for (unsigned l = 0; l <= order + 1; ++l)
                radial[l] = normalization * miller[l];
        }
        for (unsigned l = 1; l <= order; ++l)
            result.derivatives[l] = (l / x) * radial[l] - radial[l + 1];
    }
    result.derivatives[0] = -radial[1];
    std::copy_n(radial.begin(), order + 1, result.values.begin());
    for (unsigned l = 0; l <= order; ++l)
        if (!std::isfinite(result.values[l]) || !std::isfinite(result.derivatives[l]))
            throw std::runtime_error(
                "Regular spherical-wave recursion produced a non-finite value.");
    return result;
}

WidebandLevelPlan makeWidebandLevelPlan(double activeRadius, double k, double separation,
                                        double accuracy) {
    if (!std::isfinite(activeRadius) || activeRadius < 0 || !std::isfinite(k) || k <= 0 ||
        !std::isfinite(separation) || separation <= 1 || !std::isfinite(accuracy) ||
        accuracy <= 0 || accuracy >= 1)
        throw std::invalid_argument(
            "Wideband level planning requires finite positive physical scales and accuracy.");
    WidebandLevelPlan result;
    result.activeRadius = activeRadius;
    result.electricalRadius = k * activeRadius;
    result.accuracySeed = accuracy;
    if (!std::isfinite(result.electricalRadius))
        throw std::invalid_argument("Wideband active electrical radius is non-finite.");
    if (activeRadius == 0)
        return result;
    // Eq. 18 blends geometric and oscillatory asymptotics. The representation
    // switch p-p_hi<=2 avoids applying the divergent real-direction diagonal
    // form to the low-frequency regime; see the discussion following Eq. 18.
    // https://doi.org/10.1121/1.3021297
    const double low = 1 - (std::log(accuracy) + 1.5 * std::log1p(-1 / separation)) /
                               std::log(separation),
                 high = result.electricalRadius + .5 * std::pow(3 * -std::log(accuracy), 2. / 3) *
                                                      std::cbrt(result.electricalRadius),
                 largest = std::max(low, high),
                 order = largest *
                         std::pow(std::pow(low / largest, 4) + std::pow(high / largest, 4), .25);
    // The angular grid must also fit the supported recurrence range. Failing
    // this preflight never lowers the order to make a request appear feasible.
    if (!std::isfinite(order) || std::ceil(order) > maximumOrder / 2)
        throw std::runtime_error(
            "Wideband full angular order exceeds the supported resource contract.");
    result.signatureOrder = static_cast<unsigned>(std::ceil(order));
    result.translationOrder = result.quadratureOrder = 2 * result.signatureOrder;
    result.sphericalCoefficients = coefficientCount(result.signatureOrder);
    const auto samples = static_cast<std::size_t>(result.quadratureOrder + 1);
    result.angularSamples = 2 * samples * samples;
    result.representation =
        order - high <= 2 ? LevelRepresentation::planeWave : LevelRepresentation::regularSpherical;
    if (result.representation == LevelRepresentation::regularSpherical &&
        (result.electricalRadius > maximumScaledElectricalRadius ||
         result.signatureOrder > maximumScaledSphericalOrder))
        throw std::runtime_error("Wideband LF seed requires a spherical domain beyond the "
                                 "supported numerical/resource contract.");
    return result;
}

std::vector<Complex> multipoleToSpectrum(std::span<const Complex> multipoles, unsigned order,
                                         const SphereRule &rule) {
    if (multipoles.size() != coefficientCount(order))
        throw std::invalid_argument("Wideband multipoles have an invalid size.");
    auto coefficients = std::vector<Complex>(multipoles.begin(), multipoles.end());
    Complex factor{4 * pi, 0};
    for (unsigned l = 0; l <= order; ++l) {
        for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m) {
            if (!finite(coefficients[index(l, m)]))
                throw std::invalid_argument("Wideband multipole is non-finite.");
            coefficients[index(l, m)] *= factor;
        }
        factor *= Complex{0, 1};
    }
    return evaluate(coefficients, order, rule);
}

std::vector<Complex> spectrumToLocal(std::span<const Complex> incoming, const SphereRule &rule,
                                     unsigned order) {
    auto result = project(incoming, rule, order);
    Complex factor{4 * pi, 0};
    for (unsigned l = 0; l <= order; ++l) {
        for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m)
            result[index(l, m)] *= factor;
        factor *= Complex{0, -1};
    }
    return result;
}

std::vector<Complex> resampleSpectrum(std::span<const Complex> values, const SphereRule &source,
                                      const SphereRule &target, unsigned bandlimit) {
    validateValues(values, source);
    return resampleSpectrum(values, makeAngularTransferPlan(source, target, bandlimit));
}

AngularTransferPlan makeAngularTransferPlan(const SphereRule &source, const SphereRule &target,
                                            unsigned bandlimit) {
    validateOrder(bandlimit);
    if (bandlimit > source.order || bandlimit > target.order)
        throw std::invalid_argument(
            "Wideband transfer bandlimit exceeds its source or target sampling.");
    AngularTransferPlan result;
    result.sourceAzimuthCount = 2 * (source.order + 1);
    result.targetAzimuthCount = 2 * (target.order + 1);
    result.bandlimit = bandlimit;
    result.sourcePolar = polarSamples(source);
    result.targetPolar = polarSamples(target);
    const std::size_t count = static_cast<std::size_t>(bandlimit + 1) * (bandlimit + 2) / 2;
    result.recurrenceFirst.resize(count);
    result.recurrenceSecond.resize(count);
    result.diagonalScale.resize(bandlimit + 1);
    for (unsigned m = 0; m <= bandlimit; ++m) {
        if (m)
            result.diagonalScale[m] = -std::sqrt((2. * m + 1) / (2. * m));
        for (unsigned l = m + 2; l <= bandlimit; ++l) {
            const double denominator = static_cast<double>(l) * l - static_cast<double>(m) * m;
            const auto entry = triangularIndex(l, m);
            result.recurrenceFirst[entry] = std::sqrt((4. * l * l - 1) / denominator);
            result.recurrenceSecond[entry] =
                std::sqrt((2. * l + 1) * ((l - 1.) * (l - 1.) - static_cast<double>(m) * m) /
                          ((2. * l - 3) * denominator));
        }
    }
    return result;
}

std::vector<Complex> resampleSpectrum(std::span<const Complex> values,
                                      const AngularTransferPlan &plan) {
    validateOrder(plan.bandlimit);
    const unsigned order = plan.bandlimit;
    const std::size_t count = static_cast<std::size_t>(order + 1) * (order + 2) / 2;
    if (plan.sourcePolar.size() <= order || plan.targetPolar.size() <= order ||
        plan.sourceAzimuthCount != 2 * plan.sourcePolar.size() ||
        plan.targetAzimuthCount != 2 * plan.targetPolar.size() ||
        plan.recurrenceFirst.size() != count || plan.recurrenceSecond.size() != count ||
        plan.diagonalScale.size() != order + 1 ||
        values.size() != plan.sourcePolar.size() * plan.sourceAzimuthCount)
        throw std::invalid_argument(
            "Wideband separable transfer has an inconsistent plan or spectrum.");
    for (const auto value : values)
        if (!finite(value))
            throw std::invalid_argument("Wideband spectrum contains a non-finite value.");
    for (const auto *polar : {&plan.sourcePolar, &plan.targetPolar})
        for (const auto sample : *polar)
            if (!std::isfinite(sample.cosine) || std::abs(sample.cosine) > 1 ||
                !std::isfinite(sample.weight) || sample.weight <= 0)
                throw std::invalid_argument(
                    "Wideband separable transfer has an invalid polar rule.");
    for (const auto *array : {&plan.recurrenceFirst, &plan.recurrenceSecond, &plan.diagonalScale})
        for (const auto value : *array)
            if (!std::isfinite(value))
                throw std::invalid_argument(
                    "Wideband separable transfer has a non-finite recurrence.");
    std::vector<Complex> coefficients(coefficientCount(order));
    for (std::size_t polar = 0; polar < plan.sourcePolar.size(); ++polar) {
        const auto associated = associatedLegendre(plan, plan.sourcePolar[polar].cosine);
        for (int m = -static_cast<int>(order); m <= static_cast<int>(order); ++m) {
            Complex fourier{};
            for (unsigned azimuth = 0; azimuth < plan.sourceAzimuthCount; ++azimuth)
                fourier += values[polar * plan.sourceAzimuthCount + azimuth] *
                           phase(-2 * pi * m * azimuth / plan.sourceAzimuthCount);
            fourier *= 2 * pi * plan.sourcePolar[polar].weight / plan.sourceAzimuthCount;
            const unsigned absoluteM = static_cast<unsigned>(std::abs(m));
            const double sign = m < 0 && absoluteM % 2 ? -1. : 1.;
            for (unsigned l = absoluteM; l <= order; ++l)
                coefficients[index(l, m)] +=
                    fourier * (sign * associated[triangularIndex(l, absoluteM)]);
        }
    }
    std::vector<Complex> result(plan.targetPolar.size() * plan.targetAzimuthCount);
    for (std::size_t polar = 0; polar < plan.targetPolar.size(); ++polar) {
        const auto associated = associatedLegendre(plan, plan.targetPolar[polar].cosine);
        for (int m = -static_cast<int>(order); m <= static_cast<int>(order); ++m) {
            const unsigned absoluteM = static_cast<unsigned>(std::abs(m));
            const double sign = m < 0 && absoluteM % 2 ? -1. : 1.;
            Complex fourier{};
            for (unsigned l = absoluteM; l <= order; ++l)
                fourier +=
                    coefficients[index(l, m)] * (sign * associated[triangularIndex(l, absoluteM)]);
            for (unsigned azimuth = 0; azimuth < plan.targetAzimuthCount; ++azimuth)
                result[polar * plan.targetAzimuthCount + azimuth] +=
                    fourier * phase(2 * pi * m * azimuth / plan.targetAzimuthCount);
        }
    }
    return result;
}

std::vector<Complex> outgoingToParent(std::span<const Complex> child, const SphereRule &childRule,
                                      const SphereRule &parentRule, unsigned childBandlimit,
                                      double wavenumber, Vec3 childMinusParent) {
    auto result = resampleSpectrum(child, childRule, parentRule, childBandlimit);
    for (std::size_t q = 0; q < result.size(); ++q)
        result[q] *= phase(wavenumber * dot(parentRule.samples[q].direction, childMinusParent));
    return result;
}

std::vector<Complex> incomingToChild(std::span<const Complex> parent, const SphereRule &parentRule,
                                     const SphereRule &childRule, unsigned childBandlimit,
                                     double wavenumber, Vec3 childMinusParent) {
    validateValues(parent, parentRule);
    std::vector<Complex> shifted(parent.begin(), parent.end());
    for (std::size_t q = 0; q < shifted.size(); ++q)
        shifted[q] *= phase(-wavenumber * dot(parentRule.samples[q].direction, childMinusParent));
    return resampleSpectrum(shifted, parentRule, childRule, childBandlimit);
}

std::vector<Complex> diagonalTranslation(const SphereRule &rule, unsigned expansionOrder,
                                         double wavenumber, Vec3 targetMinusSource) {
    validateRule(rule);
    validateOrder(expansionOrder);
    const double distance = std::sqrt(dot(targetMinusSource, targetMinusSource)),
                 z = wavenumber * distance;
    if (!std::isfinite(distance) || !std::isfinite(z) || wavenumber <= 0 || distance <= 0 || z == 0)
        throw std::invalid_argument("Wideband diagonal translation requires finite positive kR.");
    const Vec3 direction{targetMinusSource.x / distance, targetMinusSource.y / distance,
                         targetMinusSource.z / distance};
    std::vector<Complex> radial(expansionOrder + 1);
    radial[0] = Complex{std::sin(z), std::cos(z)} / z;
    if (expansionOrder)
        radial[1] = (Complex{1 / z, 1}) * radial[0];
    for (unsigned l = 1; l < expansionOrder; ++l)
        radial[l + 1] = ((2. * l + 1) / z) * radial[l] - radial[l - 1];
    for (const auto value : radial)
        if (!finite(value))
            throw std::runtime_error(
                "Unstable diagonal translation: use the low-frequency representation.");
    std::vector<Complex> result(rule.samples.size());
    const Complex prefactor{0, -wavenumber / (16 * pi * pi)};
    for (std::size_t q = 0; q < result.size(); ++q) {
        const double mu = std::clamp(dot(rule.samples[q].direction, direction), -1., 1.);
        double previous = 0, current = 1;
        Complex factor{1, 0}, sum{};
        for (unsigned l = 0; l <= expansionOrder; ++l) {
            sum += (2. * l + 1) * factor * radial[l] * current;
            const double next = ((2. * l + 1) * mu * current - l * previous) / (l + 1.);
            previous = current;
            current = next;
            factor *= Complex{0, -1};
        }
        result[q] = prefactor * sum;
        if (!finite(result[q]))
            throw std::runtime_error("Wideband diagonal translation is non-finite.");
    }
    return result;
}

} // namespace hornsculpt::wideband
