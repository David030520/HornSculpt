#pragma once

#include <cuComplex.h>
#include <cuda_runtime_api.h>
#include <functional>
#include <hornsculpt/bem.hpp>
#include <memory>
#include <stop_token>

namespace hornsculpt::cuda_detail {
inline constexpr double deviceGmresResidualFraction = .1;
// Restarted complex GMRES with twice classical Gram-Schmidt. Both callbacks
// operate on the supplied stream; they must complete before reusing external
// scratch on another stream. The preconditioner applies in place. Geometry,
// operator representation and outer residual verification belong to the caller.
class DeviceGmresSolver {
  public:
    using PressureApplication = std::function<void(const cuDoubleComplex *, cuDoubleComplex *)>;
    using Preconditioner = std::function<void(cuDoubleComplex *)>;
    using IterationProgress = std::function<void()>;
    DeviceGmresSolver(PressureApplication pressure, Preconditioner preconditioner,
                      std::size_t pressureCount, const BemNumerics &options, cudaStream_t stream,
                      std::stop_token stop = {}, IterationProgress progress = {});
    ~DeviceGmresSolver();
    DeviceGmresSolver(const DeviceGmresSolver &) = delete;
    DeviceGmresSolver &operator=(const DeviceGmresSolver &) = delete;
    unsigned solve(const cuDoubleComplex *right, cuDoubleComplex *solution,
                   double &relativeResidual);
    std::size_t deviceBytes() const noexcept;

  private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};
} // namespace hornsculpt::cuda_detail
