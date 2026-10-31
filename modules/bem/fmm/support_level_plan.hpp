#pragma once

#include "support_tree.hpp"
#include "wideband_transfer.hpp"
#include <span>

namespace hornsculpt::wideband {
struct SupportLevelExpansionPlan {
    WidebandLevelPlan expansion;
    // Actual admissible source/target convergence ratios, including inherited
    // bounds below a far endpoint. They are not the public sum-radius MAC.
    double sourceSeparation{}, targetSeparation{}, policySeparation{};
    double maximumActiveRadius{};
    double minimumFarDistance{}, diagonalConditionIndicator{};
};
struct SupportExpansionPlan {
    std::vector<SupportLevelExpansionPlan> levels;
    // Channels needed by PP and genuinely requested face-flux far actions.
    // A requested face remains active even when its current flux is zero.
    std::vector<unsigned> activeChannels;
};
// For a far pair: delta_s=(distance-r_target)/r_source and conversely delta_t.
// See Gumerov/Duraiswami, JCP 225 (2007), section5.3, Eq.123 (b/a), and
// JASA125 (2009), Eq.18. Nested support balls preserve these lower bounds
// below each endpoint. A cumulative minimum of delta and a coarseward maximum
// of radius preserve the common-level HF-parent/LF-child transition. Neither
// operation changes the geometric partition or certifies BM operator error.
SupportExpansionPlan
makeSupportExpansionPlan(const SupportTree &tree,
                         std::span<const std::uint32_t> panelToActiveFaceColumn, double wavenumber,
                         double accuracySeed, std::stop_token stop = {});
// Radial-conditioning check, independent of the choice of geometric delta.
// By |P_l|<=1, kD*sum_l(2l+1)|h_l^(2)(kD)| bounds the continuous absolute
// translation / |G(D)|. It is not a complete rounding-error bound (angular
// quadrature and summation remain separately tested). If the coefficient-
// roundoff indicator exceeds1% of the seed, retain the full order but use
// stable spherical translation instead of an ill-conditioned diagonal.
double enforceSupportDiagonalStability(WidebandLevelPlan &plan, double wavenumber,
                                       double minimumFarDistance);
} // namespace hornsculpt::wideband
