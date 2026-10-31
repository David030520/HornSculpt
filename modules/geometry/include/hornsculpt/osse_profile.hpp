#pragma once

#include <hornsculpt/longitudinal_profile.hpp>

namespace hornsculpt {

struct OsseDirectionalParameters {
    double throatAngle{0.0};
    double k{1.0};
    double terminationScale{0.12};
    double terminationQ{0.98};
    double terminationExponent{4.0};
};

struct OsseParameters : OsseDirectionalParameters {
    double throatRadius{0.0127};
    double length{0.18};
};

// OS-SE radius and inverse from Batik, OS-SE Waveguide (2020), equations 1–2.
class OsseProfile final : public LongitudinalProfile {
  public:
    explicit OsseProfile(OsseParameters parameters);
    double radius(double z, double angle) const override;
    double slope(double z, double angle) const override;
    double angleThrough(double z, double radius) const override;
    double length() const noexcept override;

  private:
    double termination(double z) const;
    OsseParameters parameters_;
};

} // namespace hornsculpt
