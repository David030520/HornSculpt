#pragma once

#include "../galerkin/interaction_plan.hpp"
#include "../galerkin/space.hpp"
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <hornsculpt/bem.hpp>
#include <memory>
#include <stop_token>

namespace hornsculpt::cuda_detail {

// Matrix-free weak P1/DP0 Burton-Miller. Applications return the original
// unscaled Galerkin equation. Inputs and outputs reside on the supplied CUDA
// stream; distinct input/output buffers are required. A consolidated sparse
// near operator retains singular/near integration; separated blocks use FMM.
class FastMultipoleOperator {
  public:
    FastMultipoleOperator(const galerkin::Space &space, const galerkin::InteractionPlan &plan,
                          const BemNumerics &numerics, double wavenumber, cudaStream_t stream,
                          std::stop_token stop = {});
    ~FastMultipoleOperator();
    FastMultipoleOperator(const FastMultipoleOperator &) = delete;
    FastMultipoleOperator &operator=(const FastMultipoleOperator &) = delete;
    void applyPressure(const cuDoubleComplex *pressure, cuDoubleComplex *output);
    void applyFlux(const cuFloatComplex *flux, cuDoubleComplex *output);
    // Independent audit: exact regular quadrature on every far cluster product,
    // without allocating a global dense matrix. This costs quadratic arithmetic.
    void applyPressureDirect(const cuDoubleComplex *pressure, cuDoubleComplex *output);
    void applyFluxDirect(const cuFloatComplex *flux, cuDoubleComplex *output);
    // Apply spatial block-Jacobi in place. Each block contains at most 32
    // pressure DOFs and uses pivoted FP64 LU of the exact sparse near operator.
    void precondition(cuDoubleComplex *values);
    // Peak of this operator's explicit device allocations, including temporary
    // construction buffers released before the first operator application.
    std::size_t deviceBytes() const noexcept;
    const FastMultipoleStatistics &statistics() const noexcept;

  private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};
} // namespace hornsculpt::cuda_detail
