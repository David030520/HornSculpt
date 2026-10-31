#include <hornsculpt/mouth_morph.hpp>

#include <cmath>
#include <stdexcept>

namespace hornsculpt {

MouthMorph::MouthMorph(const MouthMorphParameters &parameters)
    : target_(parameters.target),
      transition_(ShapeTransitionParameters{TransitionMethod::smoothRange, parameters.startFraction,
                                            parameters.endFraction}) {}

RadialSample MouthMorph::apply(RadialSample base, double baseMouthRadius, double axialFraction,
                               double length, double phi) const {
    if (!std::isfinite(base.radius) || base.radius <= 0.0 || std::isnan(base.slope) ||
        !std::isfinite(baseMouthRadius) || baseMouthRadius <= 0.0 || !std::isfinite(length) ||
        length <= 0.0 || !std::isfinite(phi))
        throw std::invalid_argument("Mundmorphing: ungueltiges Grundprofil oder Achslaenge.");
    const double blend = transition_.blend(axialFraction);
    const double derivative = transition_.derivative(axialFraction) / length;
    if (blend == 0.0 && derivative == 0.0)
        return base;
    const double target = target_.radius(phi), displacement = target - baseMouthRadius;
    const double radius = axialFraction == 1.0 ? target : base.radius + blend * displacement;
    const double slope = derivative == 0.0 ? base.slope : base.slope + derivative * displacement;
    if (!std::isfinite(radius) || radius <= 0.0 || std::isnan(slope))
        throw std::invalid_argument("Mundmorphing erzeugt keinen positiven endlichen Querschnitt.");
    return {radius, slope};
}

} // namespace hornsculpt
