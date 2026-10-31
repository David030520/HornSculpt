#include <hornsculpt/osse_profile.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace hornsculpt {
namespace {
bool finite(double value) {
    return std::isfinite(value);
}
} // namespace

OsseProfile::OsseProfile(OsseParameters parameters) : parameters_(parameters) {
    const auto &p = parameters_;
    if (!finite(p.throatRadius) || p.throatRadius <= 0.0 || !finite(p.length) || p.length <= 0.0 ||
        !finite(p.throatAngle) || p.throatAngle < 0.0 || p.throatAngle >= std::numbers::pi / 2.0 ||
        !finite(p.k) || p.k < 0.0 || !finite(p.terminationScale) || p.terminationScale < 0.0 ||
        !finite(p.terminationQ) || p.terminationQ <= 0.0 || p.terminationQ > 1.0 ||
        !finite(p.terminationExponent) || p.terminationExponent < 2.0) {
        throw std::invalid_argument("OS-SE: ungueltige Profilparameter.");
    }
}

double OsseProfile::termination(double z) const {
    const auto &p = parameters_;
    if (p.terminationScale == 0.0 || z == 0.0)
        return 0.0;
    const double power = std::pow(p.terminationQ * z / p.length, p.terminationExponent);
    const double rise =
        power == 1.0 ? 1.0 : -std::expm1(std::log1p(-power) / p.terminationExponent);
    return p.terminationScale * p.length / p.terminationQ * rise;
}

double OsseProfile::radius(double z, double angle) const {
    if (!finite(z) || z < 0.0 || z > parameters_.length || !finite(angle) || angle < 0.0 ||
        angle >= std::numbers::pi / 2.0) {
        throw std::invalid_argument("OS-SE: Position oder Halbwinkel ungueltig.");
    }
    const auto &p = parameters_;
    const double throat = p.k * p.throatRadius;
    const double axial = z * std::tan(angle);
    const double square =
        throat * throat + 2.0 * throat * z * std::tan(p.throatAngle) + axial * axial;
    return std::sqrt(square) + p.throatRadius * (1.0 - p.k) + termination(z);
}

double OsseProfile::angleThrough(double z, double targetRadius) const {
    if (!finite(z) || z <= 0.0 || z > parameters_.length || !finite(targetRadius) ||
        targetRadius <= 0.0) {
        throw std::invalid_argument("OS-SE: Leitkurvenposition oder Radius ungueltig.");
    }
    const auto &p = parameters_;
    const double base = targetRadius - p.throatRadius * (1.0 - p.k) - termination(z);
    const double throat = p.k * p.throatRadius;
    const double minimumSquare = throat * throat + 2.0 * throat * z * std::tan(p.throatAngle);
    const double angleSquare = base * base - minimumSquare;
    const double roundoff =
        64.0 * std::numeric_limits<double>::epsilon() * std::max(base * base, minimumSquare);
    if (base < 0.0 || angleSquare < -roundoff) {
        throw std::invalid_argument("Leitkurve mit diesem OS-SE-Profil nicht erreichbar.");
    }
    return std::atan(std::sqrt(std::max(0.0, angleSquare)) / z);
}

double OsseProfile::slope(double z, double angle) const {
    // Reuse domain validation, including the q=1 mouth endpoint.
    (void)radius(z, angle);
    const auto &p = parameters_;
    const double tangent = std::tan(angle);
    const double throat = p.k * p.throatRadius;
    const double root = std::sqrt(throat * throat + 2.0 * throat * z * std::tan(p.throatAngle) +
                                  z * z * tangent * tangent);
    double result =
        root == 0.0 ? tangent : (throat * std::tan(p.throatAngle) + z * tangent * tangent) / root;
    if (p.terminationScale > 0.0 && z > 0.0) {
        const double t = p.terminationQ * z / p.length;
        if (t == 1.0)
            return std::numeric_limits<double>::infinity();
        result +=
            p.terminationScale * std::pow(t, p.terminationExponent - 1.0) *
            std::pow(1.0 - std::pow(t, p.terminationExponent), 1.0 / p.terminationExponent - 1.0);
    }
    return result;
}

double OsseProfile::length() const noexcept {
    return parameters_.length;
}

} // namespace hornsculpt
