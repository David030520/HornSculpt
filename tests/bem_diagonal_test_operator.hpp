#pragma once
#include "cuda/wideband_operator.hpp"

namespace hornsculpt::wideband {
// Test-only diagonal executor for composing isolated angular/handoff gates.
// The production hierarchy owns its resident spectra and diagonal launches.
struct SpectralInteraction {
    std::uint32_t source{}, translationClass{};
};
struct DiagonalLevelPlan {
    std::uint32_t sourceNodes{}, targetNodes{}, angularSamples{}, translationClasses{};
    // CSR by target node. Each translation class denotes an actual dyadic
    // center offset at this level. The scalar coefficient includes the
    // Helmholtz prefactor, but not angular quadrature weights.
    std::vector<std::size_t> targetOffsets;
    std::vector<SpectralInteraction> interactions;
    // class-major: translationCoefficients[class*angularSamples+direction].
    std::vector<Complex> translationCoefficients;
};

// One real, reusable GPU M2L level: O(q*numberOfFarPairs), with translation
// reuse by class. Source and local spectra are node-major. This does not
// claim a complete broadband solve; near integration, hierarchy transfers and
// discretization/audit gates are separate parts of the operator contract.
class DiagonalM2LOperator {
  public:
    explicit DiagonalM2LOperator(const DiagonalLevelPlan &plan, std::stop_token stop = {},
                                 cudaStream_t stream = nullptr);
    ~DiagonalM2LOperator();
    DiagonalM2LOperator(const DiagonalM2LOperator &) = delete;
    DiagonalM2LOperator &operator=(const DiagonalM2LOperator &) = delete;
    std::vector<Complex> apply(std::span<const Complex> source, std::stop_token stop = {});
    // Device arrays keep the constructor's fixed dimensions and stay on CUDA.
    // Only scalar failure status is read back; calls complete before returning.
    void applyDevice(const cuDoubleComplex *source, cuDoubleComplex *output,
                     std::stop_token stop = {});
    cudaStream_t stream() const noexcept;
    std::size_t deviceBytes() const noexcept;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace hornsculpt::wideband
