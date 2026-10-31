#pragma once

namespace hornsculpt {

class RadialCurve {
  public:
    virtual ~RadialCurve() = default;
    virtual double radius(double phi) const = 0;
};

class Superellipse final : public RadialCurve {
  public:
    Superellipse(double a, double b, double exponent);
    double radius(double phi) const override;
    // Enclosed area in the squared units of the supplied half-axes.
    double area() const;

  private:
    double a_, b_, exponent_;
};

// Exact centred rectangle with circular corner arcs. Dimensions are half-axes
// and a physical radius in their common length unit.
class RoundedRectangle final : public RadialCurve {
  public:
    RoundedRectangle(double a, double b, double cornerRadius);
    double radius(double phi) const override;
    double angularDerivative(double phi) const;
    double area() const;

  private:
    double a_, b_, radius_;
};

struct SuperformulaParameters {
    double a{1.0};
    double b{1.0};
    double m1{4.0};
    double m2{4.0};
    double n1{2.0};
    double n2{2.0};
    double n3{2.0};
};

class Superformula final : public RadialCurve {
  public:
    explicit Superformula(SuperformulaParameters parameters);
    double radius(double phi) const override;

  private:
    SuperformulaParameters parameters_;
};

enum class CurveType { superellipse, superformula };
enum class GuidePositionMode { relative, absolute };

struct ContourParameters {
    CurveType type{CurveType::superellipse};
    double width{0.3};
    double aspectRatio{0.7};
    double rotation{0.0};
    double superellipseExponent{4.0};
    // ATH 4.8.2, section 6.4.2: a nondegenerate example with independent
    // shape controls. The generic SuperformulaParameters default stays circular.
    SuperformulaParameters superformula{.95, .82, 4, 4, .8, 7, 1.9};
};

struct GuideParameters : ContourParameters {
    GuidePositionMode positionMode{GuidePositionMode::relative};
    double position{1.0};
};

} // namespace hornsculpt
