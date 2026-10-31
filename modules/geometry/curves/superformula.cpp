#include <hornsculpt/radial_curve.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace hornsculpt {

Superformula::Superformula(SuperformulaParameters parameters) : parameters_(parameters) {
    const auto &p = parameters_;
    for (const double value : {p.a, p.b, p.n2, p.n3}) {
        if (!std::isfinite(value) || value <= 0.0)
            throw std::invalid_argument("Superformula: a, b, n2 und n3 muessen positiv sein.");
    }
    if (!std::isfinite(p.n1) || p.n1 == 0.0)
        throw std::invalid_argument("Superformula: n1 muss endlich und ungleich null sein.");
    if (!std::isfinite(p.m1) || !std::isfinite(p.m2) || p.m1 < 0.0 || p.m2 < 0.0)
        throw std::invalid_argument(
            "Superformula: Frequenzen muessen endlich und nichtnegativ sein.");
    // Equal endpoint radii alone do not establish a periodic curve.
    constexpr std::array probes{-0.371, 0.0, 0.173, 0.791, 1.411, 2.337, 4.781};
    for (const double phi : probes) {
        const double r = radius(phi);
        const double periodic = radius(phi + 2.0 * std::numbers::pi);
        if (std::abs(r - periodic) > 1e-8 * std::max(r, periodic))
            throw std::invalid_argument(
                "Superformula: m1/m2 erzeugen keine bei 360 Grad periodische Kontur.");
        constexpr double derivativeStep = 1e-5;
        const double derivative =
            (radius(phi + derivativeStep) - radius(phi - derivativeStep)) / (2.0 * derivativeStep);
        const double periodicDerivative = (radius(phi + 2.0 * std::numbers::pi + derivativeStep) -
                                           radius(phi + 2.0 * std::numbers::pi - derivativeStep)) /
                                          (2.0 * derivativeStep);
        if (std::abs(derivative - periodicDerivative) >
            1e-6 * std::max({1.0, std::abs(derivative), std::abs(periodicDerivative)}))
            throw std::invalid_argument("Superformula besitzt keine periodische Tangente.");
    }
}

double Superformula::radius(double phi) const {
    const auto &p = parameters_;
    const double x = std::pow(std::abs(std::cos(p.m1 * phi / 4.0)) / p.a, p.n2);
    const double y = std::pow(std::abs(std::sin(p.m2 * phi / 4.0)) / p.b, p.n3);
    const double result = std::pow(x + y, -1.0 / p.n1);
    if (!std::isfinite(result) || result <= 0.0)
        throw std::invalid_argument("Superformula erzeugt keinen endlichen positiven Radius.");
    return result;
}

} // namespace hornsculpt
