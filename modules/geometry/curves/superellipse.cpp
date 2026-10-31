#include <hornsculpt/radial_curve.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace hornsculpt {

Superellipse::Superellipse(double a, double b, double exponent)
    : a_(a), b_(b), exponent_(exponent) {
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(exponent) || a <= 0.0 ||
        b <= 0.0 || exponent < 2.0) {
        throw std::invalid_argument(
            "Superellipse: Halbachsen positiv, Exponent >= 2 erforderlich.");
    }
}

double Superellipse::radius(double phi) const {
    if (!std::isfinite(phi))
        throw std::invalid_argument("Superellipse: Winkel muss endlich sein.");
    const double x = std::abs(std::cos(phi)) / a_;
    const double y = std::abs(std::sin(phi)) / b_;
    // Scaling keeps high exponents finite without changing the superellipse.
    const double largest = std::max(x, y);
    const double normalizedX = x / largest, normalizedY = y / largest;
    double norm;
    if (exponent_ == 2.0 || exponent_ == 4.0 || exponent_ == 6.0) {
        const double x2 = normalizedX * normalizedX, y2 = normalizedY * normalizedY;
        if (exponent_ == 2.0)
            norm = std::sqrt(x2 + y2);
        else if (exponent_ == 4.0)
            norm = std::sqrt(std::sqrt(x2 * x2 + y2 * y2));
        else
            norm = std::cbrt(std::sqrt(x2 * x2 * x2 + y2 * y2 * y2));
    } else {
        norm = std::pow(std::pow(normalizedX, exponent_) + std::pow(normalizedY, exponent_),
                        1.0 / exponent_);
    }
    const double result = 1.0 / (largest * norm);
    if (!std::isfinite(result) || result <= 0.0)
        throw std::invalid_argument("Superellipse erzeugt keinen endlichen positiven Radius.");
    return result;
}

double Superellipse::area() const {
    // Four times integral_0^a b*(1-(x/a)^n)^(1/n) dx, evaluated
    // by the beta integral. n>=2 keeps both gamma arguments in [1,2].
    const double gamma = std::tgamma(1 + 1 / exponent_);
    const double result = 4 * a_ * b_ * gamma * gamma / std::tgamma(1 + 2 / exponent_);
    if (!std::isfinite(result) || result <= 0)
        throw std::invalid_argument("Superellipse erzeugt keine endliche positive Flaeche.");
    return result;
}

} // namespace hornsculpt
