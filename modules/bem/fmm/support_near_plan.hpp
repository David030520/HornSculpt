#pragma once

#include "support_tree.hpp"
#include <span>

namespace hornsculpt::wideband {

struct SupportCsr {
    // Sorted global P1 columns for A, activeFaceIds column indices for B.
    std::vector<std::uint32_t> rowOffsets, columns;
};

struct SupportNearPlan {
    SupportCsr pressure, flux;
    std::vector<std::uint32_t> activeFaceIds;
    // UINT32_MAX for inactive or wall panels. Source velocities are constant
    // within each active face, with arbitrary complex values in each RHS.
    std::vector<std::uint32_t> panelToActiveFaceColumn;
    // Complete, deduplicated incident panels of P1 supports, indexed by node.
    // Interior nodes have empty spans; each leaf span is sorted by panel ID.
    std::vector<std::uint32_t> nodePanelOffsets, supportPanels;
    // Target-leaf adjacency into tree.pressureFlux.nearInteractions. Retaining
    // this compact block graph avoids expanding inactive wall panel columns.
    std::vector<std::uint32_t> fluxBlockOffsets, fluxBlockIndices;
    std::size_t peakHostBytes{};
    std::size_t allocatedBytes() const noexcept;
};

// The host budget includes owned vector capacities and construction scratch,
// excluding the caller-owned Space/SupportTree and allocator bookkeeping.
// No global panel-pair table is constructed. CSR indices must fit uint32_t.
SupportNearPlan makeSupportNearPlan(const galerkin::Space &space, const SupportTree &tree,
                                    std::span<const std::uint32_t> activeFaceIds,
                                    std::size_t maximumHostBytes, std::stop_token stop = {});

struct SupportNearJob {
    std::uint32_t testPanel{}, sourcePanel{};
    // Local P1/P1 bit (3*testLocal+sourceLocal), and P1/DP0 bit testLocal.
    // All near contributions of one panel pair are combined into these masks.
    std::uint16_t pressureMask{};
    std::uint8_t fluxMask{};
};

// Non-owning iterator; Space, Tree and NearPlan must outlive it. The immutable
// plan can be reused for any pressure vector and multiple complex face RHSs.
class SupportNearJobIterator {
  public:
    SupportNearJobIterator(const galerkin::Space &space, const SupportTree &tree,
                           const SupportNearPlan &plan,
                           std::size_t maximumWorkspaceBytes = 64 * 1024 * 1024);
    // Caller storage is the only job buffer. The explicit byte budget is
    // checked before any write; no hidden batch-sized allocation is made.
    std::size_t next(std::span<SupportNearJob> output, std::size_t maximumBatchBytes,
                     std::stop_token stop = {});
    bool done() const noexcept;
    // Actual owned capacities: one stamp, combined mask and candidate ID per
    // source panel. This is O(T), independent of the near panel-pair count.
    std::size_t workspaceBytes() const noexcept;

  private:
    const galerkin::Space &space_;
    const SupportTree &tree_;
    const SupportNearPlan &plan_;
    void prepare(std::stop_token stop);
    std::uint16_t &mask(std::uint32_t source);
    std::vector<std::uint32_t> stamps_, candidates_;
    std::vector<std::uint16_t> masks_;
    std::size_t testPanel_{}, count_{}, cursor_{};
    bool prepared_{}, done_{};
};

} // namespace hornsculpt::wideband
