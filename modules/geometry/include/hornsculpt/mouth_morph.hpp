#pragma once

#include <hornsculpt/longitudinal_profile.hpp>
#include <hornsculpt/physical_contour.hpp>
#include <hornsculpt/shape_transition.hpp>

namespace hornsculpt {

struct MouthMorphParameters {
    bool enabled{false};
    double startFraction{.75};
    double endFraction{1.0};
    ContourParameters target;
};

// C2 radial displacement to a physical mouth contour. Its zero endpoint slopes
// preserve the input profile's tangents at the entry and mouth.
class MouthMorph {
  public:
    explicit MouthMorph(const MouthMorphParameters &parameters);
    RadialSample apply(RadialSample base, double baseMouthRadius, double axialFraction,
                       double length, double phi) const;

  private:
    PhysicalContour target_;
    ShapeTransition transition_;
};

} // namespace hornsculpt
