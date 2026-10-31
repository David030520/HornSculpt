#pragma once
#include <algorithm>
#include <cmath>
#include <fmm/scaled_translation.hpp>
#include <limits>
#include <stdexcept>

namespace hornsculpt::test::reference {
using namespace wideband;
namespace projection_detail {
inline std::size_t count(unsigned order) {
    return static_cast<std::size_t>(order + 1) * (order + 1);
}
inline std::size_t index(unsigned l, int m) {
    return static_cast<std::size_t>(l) * l + l + m;
}
inline double norm(Vec3 v) {
    return std::hypot(v.x, v.y, v.z);
}
inline Vec3 scale(Vec3 v, double s) {
    return {v.x * s, v.y * s, v.z * s};
}
inline Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline double axis(Vec3 v, unsigned i) {
    return i == 0 ? v.x : i == 1 ? v.y : v.z;
}
inline void domain(double k, double radius) {
    if (!std::isfinite(k) || k < 0 || !std::isfinite(radius) || radius <= 0 || k * radius > 1)
        throw std::invalid_argument("Full-sphere projection oracle requires 0<=ka<=1.");
}
inline void expansion(unsigned order, double radius, std::span<const Complex> values, double k) {
    domain(k, radius);
    if (values.size() != count(order))
        throw std::invalid_argument("Projection oracle coefficient count is invalid.");
    for (const auto value : values)
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::invalid_argument("Projection oracle coefficients are non-finite.");
}
// Independent small-argument series for the test-only full-sphere projection.
// The production axial builder instead uses its stable turning-region radial
// primitive and analytically integrates azimuth before constructing m blocks.
inline std::vector<double> normalizedRegular(unsigned order, double z) {
    if (!std::isfinite(z) || z < 0 || z > 1)
        throw std::invalid_argument("Projection oracle radial series requires 0<=z<=1.");
    std::vector<double> result(order + 1);
    for (unsigned l = 0; l <= order; ++l) {
        double term = 1, sum = 1;
        for (unsigned n = 1; n < 128; ++n) {
            term *= -z * z / (2. * n * (2. * l + 2 * n + 1));
            sum += term;
            if (std::abs(term) <= std::numeric_limits<double>::epsilon() * std::abs(sum))
                break;
        }
        result[l] = sum;
    }
    return result;
}
template <class Field>
ScaledLocal projectLocal(unsigned order, double k, double radius, unsigned angularOrder,
                         Field field) {
    domain(k, radius);
    if (angularOrder < 2 * order)
        throw std::invalid_argument(
            "Scaled local projection requires independently oversampled angular quadrature.");
    ScaledLocal result{order, radius, std::vector<Complex>(count(order))};
    const auto radial = normalizedRegular(order + 1, k * radius);
    for (const auto &sample : makeSphereRule(angularOrder).samples) {
        const auto harmonic = sphericalHarmonics(order, sample.direction);
        const auto value = field(scale(sample.direction, radius));
        Complex normalDerivative{};
        for (unsigned component = 0; component < 3; ++component)
            normalDerivative +=
                value.gradient[component] * (radius * axis(sample.direction, component));
        for (unsigned l = 0; l <= order; ++l) {
            const double derivative =
                             l * radial[l] - k * k * radius * radius / (2. * l + 3) * radial[l + 1],
                         denominator = radial[l] * radial[l] + derivative * derivative;
            if (!(denominator > 0) || !std::isfinite(denominator))
                throw std::runtime_error(
                    "Scaled Cauchy projection has an unresolved radial normalization.");
            for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m)
                result.coefficients[index(l, m)] +=
                    sample.weight * std::conj(harmonic[index(l, m)]) *
                    (radial[l] * value.value + derivative * normalDerivative) / denominator;
        }
    }
    return result;
}

} // namespace projection_detail
inline ScaledMultipole scaledMultipoleToParent(const ScaledMultipole &child, double k, Vec3 shift,
                                               double parentRadius, unsigned parentOrder,
                                               unsigned angularOrder) {
    projection_detail::expansion(child.order, child.radius, child.coefficients, k);
    projection_detail::domain(k, parentRadius);
    if (projection_detail::norm(shift) + child.radius > parentRadius ||
        angularOrder < 2 * std::max(child.order, parentOrder))
        throw std::invalid_argument(
            "Scaled M2M reference needs nested enclosing spheres and oversampling.");
    ScaledMultipole result{parentOrder, parentRadius,
                           std::vector<Complex>(projection_detail::count(parentOrder))};
    const double traceRadius = child.radius * 1.01;
    for (const auto &sample : makeSphereRule(angularOrder).samples) {
        const Vec3 tracePoint = projection_detail::scale(sample.direction, traceRadius);
        const auto outgoing = evaluateScaledMultipole(child, k, tracePoint);
        Complex outgoingNormal{};
        for (unsigned component = 0; component < 3; ++component)
            outgoingNormal +=
                outgoing.gradient[component] * projection_detail::axis(sample.direction, component);
        const auto basis = scaledRegularWaves(parentOrder, k, parentRadius,
                                              projection_detail::add(shift, tracePoint));
        for (std::size_t j = 0; j < result.coefficients.size(); ++j) {
            Complex regularNormal{};
            for (unsigned component = 0; component < 3; ++component)
                regularNormal += basis.waves[j].gradient[component] *
                                 projection_detail::axis(sample.direction, component);
            result.coefficients[j] += traceRadius * traceRadius * sample.weight *
                                      (outgoing.value * std::conj(regularNormal) -
                                       outgoingNormal * std::conj(basis.waves[j].value));
        }
    }
    return result;
}

inline ScaledLocal scaledMultipoleToLocal(const ScaledMultipole &source, double k, Vec3 separation,
                                          double targetRadius, unsigned targetOrder,
                                          unsigned angularOrder) {
    projection_detail::expansion(source.order, source.radius, source.coefficients, k);
    if (projection_detail::norm(separation) <= source.radius + targetRadius)
        throw std::invalid_argument("Scaled M2L reference requires disjoint enclosing spheres.");
    return projection_detail::projectLocal(
        targetOrder, k, targetRadius, angularOrder, [&](Vec3 point) {
            return evaluateScaledMultipole(source, k, projection_detail::add(separation, point));
        });
}

inline ScaledLocal scaledLocalToChild(const ScaledLocal &parent, double k, Vec3 shift,
                                      double childRadius, unsigned childOrder,
                                      unsigned angularOrder) {
    projection_detail::expansion(parent.order, parent.radius, parent.coefficients, k);
    if (projection_detail::norm(shift) + childRadius > parent.radius)
        throw std::invalid_argument("Scaled L2L reference requires nested enclosing spheres.");
    return projection_detail::projectLocal(
        childOrder, k, childRadius, angularOrder, [&](Vec3 point) {
            return evaluateScaledLocal(parent, k, projection_detail::add(shift, point));
        });
}

} // namespace hornsculpt::test::reference
