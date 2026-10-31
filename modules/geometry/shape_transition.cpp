#include <hornsculpt/shape_transition.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace hornsculpt {

ShapeTransition::ShapeTransition(ShapeTransitionParameters parameters)
    : parameters_(std::move(parameters)) {
    const auto &p = parameters_;
    if (p.method == TransitionMethod::smoothRange) {
        if (!std::isfinite(p.startFraction) || !std::isfinite(p.endFraction) ||
            p.startFraction < 0.0 || p.endFraction > 1.0 || p.startFraction >= p.endFraction)
            throw std::invalid_argument("Formuebergang: 0 <= Beginn < Ende <= 1 erforderlich.");
        return;
    }
    const auto &knots = p.knots;
    if (knots.size() < 2 || knots.front().fraction != 0.0 || knots.back().fraction != 1.0 ||
        knots.front().blend != 0.0 || knots.back().blend != 1.0)
        throw std::invalid_argument("LUT muss mit (0,0) beginnen und mit (1,1) enden.");
    for (std::size_t i = 0; i < knots.size(); ++i) {
        if (!std::isfinite(knots[i].fraction) || !std::isfinite(knots[i].blend) ||
            knots[i].fraction < 0.0 || knots[i].fraction > 1.0 || knots[i].blend < 0.0 ||
            knots[i].blend > 1.0 ||
            (i > 0 &&
             (knots[i].fraction <= knots[i - 1].fraction || knots[i].blend < knots[i - 1].blend)))
            throw std::invalid_argument(
                "LUT: Positionen streng steigend und Formanteile monoton 0..1 erforderlich.");
    }
    std::vector<double> width(knots.size() - 1), secant(knots.size() - 1);
    for (std::size_t i = 0; i + 1 < knots.size(); ++i) {
        width[i] = knots[i + 1].fraction - knots[i].fraction;
        secant[i] = (knots[i + 1].blend - knots[i].blend) / width[i];
    }
    slopes_.resize(knots.size());
    // Zero endpoint slopes preserve the entry tangent and a settled mouth shape.
    slopes_.front() = slopes_.back() = 0.0;
    for (std::size_t i = 1; i + 1 < knots.size(); ++i) {
        if (secant[i - 1] == 0.0 || secant[i] == 0.0) {
            slopes_[i] = 0.0;
        } else {
            const double w1 = 2.0 * width[i] + width[i - 1];
            const double w2 = width[i] + 2.0 * width[i - 1];
            slopes_[i] = (w1 + w2) / (w1 / secant[i - 1] + w2 / secant[i]);
        }
    }
}

double ShapeTransition::blend(double fraction) const {
    if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0)
        throw std::invalid_argument("Formuebergangsposition ausserhalb 0..1.");
    if (parameters_.method == TransitionMethod::smoothRange) {
        if (fraction <= parameters_.startFraction)
            return 0.0;
        if (fraction >= parameters_.endFraction)
            return 1.0;
        const double t = (fraction - parameters_.startFraction) /
                         (parameters_.endFraction - parameters_.startFraction);
        return t * t * t * (t * (6.0 * t - 15.0) + 10.0);
    }
    const auto &knots = parameters_.knots;
    if (fraction == 1.0)
        return 1.0;
    const auto upper = std::upper_bound(
        knots.begin(), knots.end(), fraction,
        [](double value, const TransitionKnot &knot) { return value < knot.fraction; });
    const auto i = static_cast<std::size_t>(upper - knots.begin() - 1);
    const double width = knots[i + 1].fraction - knots[i].fraction;
    const double t = (fraction - knots[i].fraction) / width;
    const double t2 = t * t, t3 = t2 * t;
    return (2.0 * t3 - 3.0 * t2 + 1.0) * knots[i].blend + (t3 - 2.0 * t2 + t) * width * slopes_[i] +
           (-2.0 * t3 + 3.0 * t2) * knots[i + 1].blend + (t3 - t2) * width * slopes_[i + 1];
}

double ShapeTransition::derivative(double fraction) const {
    (void)blend(fraction);
    if (parameters_.method == TransitionMethod::smoothRange) {
        if (fraction <= parameters_.startFraction || fraction >= parameters_.endFraction)
            return 0.0;
        const double width = parameters_.endFraction - parameters_.startFraction;
        const double t = (fraction - parameters_.startFraction) / width;
        return 30.0 * t * t * (t - 1.0) * (t - 1.0) / width;
    }
    if (fraction == 1.0)
        return slopes_.back();
    const auto &knots = parameters_.knots;
    const auto upper = std::upper_bound(
        knots.begin(), knots.end(), fraction,
        [](double value, const TransitionKnot &knot) { return value < knot.fraction; });
    const auto i = static_cast<std::size_t>(upper - knots.begin() - 1);
    const double width = knots[i + 1].fraction - knots[i].fraction;
    const double t = (fraction - knots[i].fraction) / width;
    return (6.0 * t * t - 6.0 * t) * knots[i].blend / width +
           (3.0 * t * t - 4.0 * t + 1.0) * slopes_[i] +
           (-6.0 * t * t + 6.0 * t) * knots[i + 1].blend / width +
           (3.0 * t * t - 2.0 * t) * slopes_[i + 1];
}

} // namespace hornsculpt
