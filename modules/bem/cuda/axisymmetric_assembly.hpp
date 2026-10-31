#pragma once

#include "../galerkin/axisymmetric_quadrature.hpp"
#include <cuComplex.h>
#include <memory>

namespace hornsculpt::cuda_detail {
struct RingObservation {
    double radius{}, axial{}, directionR{}, directionZ{}, originPhase{};
    unsigned blocked{};
};

// Continuous-revolution assembly only. The shared baffled solver retains LU,
// original-system auditing and the power-conjugate LEM block elimination.
class AxisymmetricAssembly {
  public:
    AxisymmetricAssembly(const galerkin::RingInteractionPlan &plan,
                         const std::vector<int> &elementApertures,
                         const std::vector<unsigned> &representativePanels,
                         const std::vector<unsigned> &faceIds);
    ~AxisymmetricAssembly();
    AxisymmetricAssembly(const AxisymmetricAssembly &) = delete;
    AxisymmetricAssembly &operator=(const AxisymmetricAssembly &) = delete;
    std::size_t bytes() const;
    void assemble(unsigned pressureCount, unsigned dimension, unsigned panelCount,
                  unsigned rhsCount, double k, double tolerance, unsigned maximumDepth,
                  const double *rowScales, const cuDoubleComplex *knownFlux,
                  cuDoubleComplex *matrix, cuDoubleComplex *rhs, std::stop_token stop);
    void field(const std::vector<RingObservation> &observations, bool farField,
               double referenceDistance, unsigned pressureCount, unsigned dimension,
               unsigned rhsCount, double k, double tolerance, unsigned maximumDepth,
               const cuDoubleComplex *solution, cuDoubleComplex *output, std::stop_token stop);

  private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};
} // namespace hornsculpt::cuda_detail
