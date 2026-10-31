#pragma once

namespace hornsculpt {

struct RadialSample {
    double radius{};
    double slope{};
};

// A radial longitudinal profile in metres and radians. angleThrough is the
// provider's inverse; the composition layer does not assume an OS-SE equation.
class LongitudinalProfile {
  public:
    virtual ~LongitudinalProfile() = default;
    virtual double radius(double z, double angle) const = 0;
    virtual double slope(double z, double angle) const;
    virtual double angleThrough(double z, double radius) const = 0;
    virtual double length() const noexcept = 0;
};

} // namespace hornsculpt
