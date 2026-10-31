#pragma once

#include "support_near_operator.hpp"

namespace hornsculpt::wideband {
inline constexpr std::size_t maximumDirectAuditPanels = 20000;
// Matrix-free original weak Galerkin action for an explicitly bounded small
// mesh. Reuses the exact near operator; every remaining panel contribution is
// integrated directly, with DOF masks from the same geometric partition.
// It owns no second near matrix, dense matrix, or multipole hierarchy.
class SupportDirectAudit {
  public:
    SupportDirectAudit(const galerkin::Space &space, const SupportTree &tree,
                       const SupportNearPlan &plan, std::shared_ptr<SupportNearOperator> near,
                       double wavenumber, unsigned regularOrder, std::stop_token stop = {});
    ~SupportDirectAudit();
    SupportDirectAudit(const SupportDirectAudit &) = delete;
    SupportDirectAudit &operator=(const SupportDirectAudit &) = delete;
    void applyPressureDevice(const cuDoubleComplex *input, cuDoubleComplex *output,
                             std::stop_token stop = {});
    void applyFaceFluxDevice(const cuDoubleComplex *input, cuDoubleComplex *output,
                             std::stop_token stop = {});
    std::size_t ownDeviceBytes() const noexcept;

  private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};
} // namespace hornsculpt::wideband
