#include <hornsculpt/radial_curve.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace hornsculpt {
RoundedRectangle::RoundedRectangle(double a, double b, double cornerRadius)
    : a_(a), b_(b), radius_(cornerRadius) {
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(cornerRadius) || a <= 0 ||
        b <= 0 || cornerRadius < 0 || cornerRadius > std::min(a, b))
        throw std::invalid_argument(
            "Rechteck: positive Halbbreite/-hoehe und Eckradius zwischen 0 und "
            "halber kleinster Seitenlaenge erforderlich.");
}

double RoundedRectangle::radius(double phi) const {
    if (!std::isfinite(phi))
        throw std::invalid_argument("Rechteck: endlicher Winkel erforderlich.");
    const double c = std::abs(std::cos(phi)), s = std::abs(std::sin(phi));
    if (a_ * s <= (b_ - radius_) * c)
        return a_ / c;
    if (b_ * c <= (a_ - radius_) * s)
        return b_ / s;
    const double x = a_ - radius_, y = b_ - radius_;
    const double along = x * c + y * s;
    const double across = x * s - y * c;
    return along + std::sqrt(std::max(0., (radius_ - across) * (radius_ + across)));
}

double RoundedRectangle::angularDerivative(double phi) const {
    const double r = radius(phi);
    const double c = std::cos(phi), s = std::sin(phi);
    double nx{}, ny{};
    if (a_ * std::abs(s) <= (b_ - radius_) * std::abs(c))
        nx = std::copysign(1., c);
    else if (b_ * std::abs(c) <= (a_ - radius_) * std::abs(s))
        ny = std::copysign(1., s);
    else {
        nx = r * c - std::copysign(a_ - radius_, c);
        ny = r * s - std::copysign(b_ - radius_, s);
    }
    return r * (nx * s - ny * c) / (nx * c + ny * s);
}

double RoundedRectangle::area() const {
    return 4 * a_ * b_ - (4 - std::numbers::pi) * radius_ * radius_;
}
} // namespace hornsculpt
