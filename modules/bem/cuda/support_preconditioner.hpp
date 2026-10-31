#pragma once

#include "support_near_operator.hpp"

namespace hornsculpt::wideband {
// Spatial principal blocks extracted from the exact near matrix. Entries
// outside its CSR pattern are zero; this is a preconditioner, not a replacement
// for the complete BEM operator. Pressure order is a complete DOF permutation.
class SupportBlockJacobi {
  public:
    SupportBlockJacobi(SupportDeviceCsr nearPressure, std::span<const std::uint32_t> pressureOrder,
                       unsigned blockWidth = 32, std::stop_token stop = {},
                       cudaStream_t stream = nullptr);
    ~SupportBlockJacobi();
    SupportBlockJacobi(const SupportBlockJacobi &) = delete;
    SupportBlockJacobi &operator=(const SupportBlockJacobi &) = delete;
    void applyDevice(cuDoubleComplex *values, std::stop_token stop = {});
    std::size_t deviceBytes() const noexcept;
    std::size_t peakDeviceBytes() const noexcept;

  private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};
} // namespace hornsculpt::wideband
