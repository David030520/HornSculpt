#pragma once

#include <vector>

namespace hornsculpt {

enum class TransitionMethod { smoothRange, lut };

struct TransitionKnot {
    double fraction{};
    double blend{};
};

struct ShapeTransitionParameters {
    TransitionMethod method{TransitionMethod::smoothRange};
    double startFraction{0.0};
    double endFraction{0.65};
    std::vector<TransitionKnot> knots{{0.0, 0.0}, {1.0, 1.0}};
};

// Monotone axial blending with zero endpoint slopes; independent of profile type.
class ShapeTransition {
  public:
    explicit ShapeTransition(ShapeTransitionParameters parameters);
    double blend(double fraction) const;
    double derivative(double fraction) const;

  private:
    ShapeTransitionParameters parameters_;
    std::vector<double> slopes_;
};

} // namespace hornsculpt
