#pragma once

#include <cuComplex.h>
#include <cuda_runtime_api.h>
#include <hornsculpt/bem.hpp>
#include <memory>
#include <stop_token>

namespace hornsculpt::cuda_detail {
class FastMultipoleOperator;

struct IterativeSolution {
    double boundaryResidual{};
    double residualCheckMilliseconds{};
    unsigned iterations{}, refinementSteps{};
};

// Restarted, twice-orthogonalized complex GMRES with a spatial block-Jacobi left
// preconditioner. Acceptance uses the direct original Galerkin operator.
class FastMultipoleLinearSolver {
  public:
    FastMultipoleLinearSolver(FastMultipoleOperator &operation, std::size_t pressureCount,
                              const BemNumerics &options, cudaStream_t stream,
                              std::stop_token stop);
    ~FastMultipoleLinearSolver();
    FastMultipoleLinearSolver(const FastMultipoleLinearSolver &) = delete;
    FastMultipoleLinearSolver &operator=(const FastMultipoleLinearSolver &) = delete;
    IterativeSolution solve(const cuFloatComplex *flux, cuDoubleComplex *pressure,
                            bool useInitialGuess = false);
    std::size_t deviceBytes() const;

  private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};
} // namespace hornsculpt::cuda_detail
