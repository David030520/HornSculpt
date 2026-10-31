#pragma once

#include <hornsculpt/radial_curve.hpp>

#include <memory>

namespace hornsculpt {

// Physical width/height and rotation are independent of the curve's shape parameters.
class PhysicalContour final : public RadialCurve {
  public:
    explicit PhysicalContour(const ContourParameters &parameters);
    double radius(double phi) const override;

  private:
    std::shared_ptr<const RadialCurve> curve_;
    double scaleX_{}, scaleY_{}, rotation_{};
};

} // namespace hornsculpt
