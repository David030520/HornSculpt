#include "support_near_plan.hpp"
#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace hornsculpt::wideband {
namespace {
template <class T> std::size_t bytes(const std::vector<T> &v) {
    return v.capacity() * sizeof(T);
}
std::size_t addBytes(std::size_t a, std::size_t b) {
    if (b > std::numeric_limits<std::size_t>::max() - a)
        throw std::overflow_error("Support near buffer size overflow.");
    return a + b;
}
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Support near construction cancelled.");
}
void prefix(std::vector<std::uint32_t> &offsets) {
    std::uint64_t total = 0;
    for (auto &offset : offsets) {
        total += offset;
        if (total > UINT32_MAX)
            throw std::overflow_error("Support near CSR exceeds 32-bit indexing.");
        offset = static_cast<std::uint32_t>(total);
    }
}
class Builder {
  public:
    Builder(const galerkin::Space &geometry, const SupportTree &hierarchy, std::size_t budget,
            std::stop_token token)
        : space(geometry), tree(hierarchy), maximumBytes(budget), stop(token) {}
    SupportNearPlan build(std::span<const std::uint32_t> faces) {
        cancelled(stop);
        if (faces.size() > UINT32_MAX)
            throw std::overflow_error("Active near face column index overflow.");
        if (space.vertices.size() != tree.pressureOrder.size() ||
            space.panels.size() != tree.fluxOrder.size() || tree.nodes.empty() ||
            tree.pressureLeaf.size() != space.vertices.size() ||
            tree.fluxLeaf.size() != space.panels.size())
            throw std::invalid_argument("Support near space and hierarchy disagree.");
        allocate(plan.activeFaceIds, faces.size());
        std::copy(faces.begin(), faces.end(), plan.activeFaceIds.begin());
        allocate(plan.panelToActiveFaceColumn, space.panels.size());
        std::fill(plan.panelToActiveFaceColumn.begin(), plan.panelToActiveFaceColumn.end(),
                  UINT32_MAX);
        // Active face counts are small, explicit RHS dimensions. This avoids a
        // global panel-pair hash and retains the caller's exact column order.
        allocate(marks, faces.size());
        for (std::size_t i = 0; i < faces.size(); ++i) {
            if (std::find(faces.begin(), faces.begin() + i, faces[i]) != faces.begin() + i)
                throw std::invalid_argument("Active near face IDs are duplicated.");
        }
        for (std::size_t p = 0; p < space.panels.size(); ++p) {
            cancelled(stop);
            const auto found = std::find(faces.begin(), faces.end(), space.panels[p].faceId);
            if (found == faces.end())
                continue;
            if (space.panels[p].boundary != BoundaryCondition::driven)
                throw std::invalid_argument("Active near face includes a wall panel.");
            const auto column = static_cast<std::uint32_t>(found - faces.begin());
            plan.panelToActiveFaceColumn[p] = column;
            marks[column] = 1;
        }
        if (std::find(marks.begin(), marks.end(), 0u) != marks.end())
            throw std::invalid_argument("Active near face has no driven panels.");
        releaseScratch();
        allocate(plan.pressure.rowOffsets, space.vertices.size() + 1);
        allocate(plan.flux.rowOffsets, space.vertices.size() + 1);
        fluxCsr();
        releaseScratch();
        pressureCsr();
        releaseScratch();
        supportPanels();
        releaseScratch();
        record();
        return std::move(plan);
    }

  private:
    std::size_t used() const {
        auto result = plan.allocatedBytes();
        for (auto n : {bytes(offsets), bytes(indices), bytes(cursor), bytes(columns), bytes(marks)})
            result = addBytes(result, n);
        return result;
    }
    void record() {
        const auto current = used();
        if (current > maximumBytes)
            throw std::runtime_error("Support near host buffer budget exceeded.");
        plan.peakHostBytes = std::max(plan.peakHostBytes, current);
    }
    template <class T>
    void allocate(std::vector<T> &v, std::size_t count, bool reserveOnly = false) {
        if (!v.empty() || v.capacity())
            throw std::logic_error("Support near buffer allocated twice.");
        if (count > std::numeric_limits<std::size_t>::max() / sizeof(T) ||
            addBytes(used(), count * sizeof(T)) > maximumBytes)
            throw std::runtime_error("Support near host buffer budget exceeded before allocation.");
        if (reserveOnly)
            v.reserve(count);
        else
            v.resize(count);
        record();
    }
    void releaseScratch() {
        std::vector<std::uint32_t>().swap(offsets);
        std::vector<std::uint32_t>().swap(indices);
        std::vector<std::uint32_t>().swap(cursor);
        std::vector<std::uint32_t>().swap(columns);
        std::vector<std::uint32_t>().swap(marks);
    }
    void adjacency(const SupportPartition &partition) {
        if (partition.nearInteractions.size() > UINT32_MAX)
            throw std::overflow_error("Support near block index overflow.");
        allocate(offsets, tree.nodes.size() + 1);
        allocate(indices, partition.nearInteractions.size());
        allocate(cursor, tree.nodes.size());
        for (const auto block : partition.nearInteractions) {
            if (block.target >= tree.nodes.size() || block.source >= tree.nodes.size() ||
                !tree.nodes[block.target].leaf() || !tree.nodes[block.source].leaf())
                throw std::invalid_argument("Support near block is not a leaf Cartesian product.");
            ++offsets[std::size_t(block.target) + 1];
        }
        prefix(offsets);
        std::copy_n(offsets.begin(), cursor.size(), cursor.begin());
        for (std::uint32_t i = 0; i < partition.nearInteractions.size(); ++i)
            indices[cursor[partition.nearInteractions[i].target]++] = i;
    }
    std::size_t faceColumns(std::size_t node) {
        const auto epoch = static_cast<std::uint32_t>(node + 1);
        std::size_t count = 0;
        for (auto b = offsets[node]; b < offsets[node + 1]; ++b) {
            const auto &source = tree.nodes[tree.pressureFlux.nearInteractions[indices[b]].source];
            for (std::size_t j = source.firstFlux;
                 j < std::size_t(source.firstFlux) + source.fluxCount; ++j) {
                const auto column = plan.panelToActiveFaceColumn[tree.fluxOrder[j]];
                if (column != UINT32_MAX && marks[column] != epoch) {
                    marks[column] = epoch;
                    columns[count++] = column;
                }
            }
        }
        std::sort(columns.begin(), columns.begin() + count);
        return count;
    }
    void fluxCsr() {
        adjacency(tree.pressureFlux);
        allocate(columns, plan.activeFaceIds.size());
        allocate(marks, plan.activeFaceIds.size());
        for (std::size_t n = 0; n < tree.nodes.size(); ++n) {
            cancelled(stop);
            const auto &node = tree.nodes[n];
            if (!node.leaf() || !node.pressureCount)
                continue;
            const auto count = faceColumns(n);
            for (std::size_t j = node.firstPressure;
                 j < std::size_t(node.firstPressure) + node.pressureCount; ++j)
                plan.flux.rowOffsets[std::size_t(tree.pressureOrder[j]) + 1] =
                    static_cast<std::uint32_t>(count);
        }
        prefix(plan.flux.rowOffsets);
        allocate(plan.flux.columns, plan.flux.rowOffsets.back());
        std::fill(marks.begin(), marks.end(), 0);
        for (std::size_t n = 0; n < tree.nodes.size(); ++n) {
            cancelled(stop);
            const auto &node = tree.nodes[n];
            if (!node.leaf() || !node.pressureCount)
                continue;
            const auto count = faceColumns(n);
            for (std::size_t j = node.firstPressure;
                 j < std::size_t(node.firstPressure) + node.pressureCount; ++j)
                std::copy_n(columns.begin(), count,
                            plan.flux.columns.begin() +
                                plan.flux.rowOffsets[tree.pressureOrder[j]]);
        }
        plan.fluxBlockOffsets = std::move(offsets);
        plan.fluxBlockIndices = std::move(indices);
    }
    void pressureCsr() {
        adjacency(tree.pressurePressure);
        allocate(columns, space.vertices.size());
        for (std::size_t n = 0; n < tree.nodes.size(); ++n) {
            cancelled(stop);
            const auto &node = tree.nodes[n];
            if (!node.leaf() || !node.pressureCount)
                continue;
            std::uint64_t count = 0;
            for (auto b = offsets[n]; b < offsets[n + 1]; ++b)
                count += tree.nodes[tree.pressurePressure.nearInteractions[indices[b]].source]
                             .pressureCount;
            if (count > space.vertices.size())
                throw std::invalid_argument("Support near source DOFs are duplicated.");
            for (std::size_t j = node.firstPressure;
                 j < std::size_t(node.firstPressure) + node.pressureCount; ++j)
                plan.pressure.rowOffsets[std::size_t(tree.pressureOrder[j]) + 1] =
                    static_cast<std::uint32_t>(count);
        }
        prefix(plan.pressure.rowOffsets);
        if (plan.pressure.rowOffsets.back() != tree.pressurePressure.representedNearPairs)
            throw std::invalid_argument("Support near CSR and represented pair count disagree.");
        allocate(plan.pressure.columns, plan.pressure.rowOffsets.back());
        for (std::size_t n = 0; n < tree.nodes.size(); ++n) {
            cancelled(stop);
            const auto &node = tree.nodes[n];
            if (!node.leaf() || !node.pressureCount)
                continue;
            std::size_t count = 0;
            for (auto b = offsets[n]; b < offsets[n + 1]; ++b) {
                const auto &source =
                    tree.nodes[tree.pressurePressure.nearInteractions[indices[b]].source];
                std::copy_n(tree.pressureOrder.begin() + source.firstPressure, source.pressureCount,
                            columns.begin() + count);
                count += source.pressureCount;
            }
            std::sort(columns.begin(), columns.begin() + count);
            if (std::adjacent_find(columns.begin(), columns.begin() + count) !=
                columns.begin() + count)
                throw std::invalid_argument("Support near CSR duplicates a source DOF.");
            for (std::size_t j = node.firstPressure;
                 j < std::size_t(node.firstPressure) + node.pressureCount; ++j)
                std::copy_n(columns.begin(), count,
                            plan.pressure.columns.begin() +
                                plan.pressure.rowOffsets[tree.pressureOrder[j]]);
        }
    }
    void supportPanels() {
        allocate(plan.nodePanelOffsets, tree.nodes.size() + 1);
        allocate(plan.supportPanels, space.incidences.size(), true);
        allocate(marks, space.panels.size());
        for (std::size_t n = 0; n < tree.nodes.size(); ++n) {
            cancelled(stop);
            const auto &node = tree.nodes[n];
            const auto first = plan.supportPanels.size();
            if (node.leaf())
                for (std::size_t j = node.firstPressure;
                     j < std::size_t(node.firstPressure) + node.pressureCount; ++j) {
                    const auto vertex = tree.pressureOrder[j];
                    for (auto k = space.incidenceOffsets[vertex];
                         k < space.incidenceOffsets[vertex + 1]; ++k) {
                        const auto panel = space.incidences[k].panel;
                        if (marks[panel] != n + 1) {
                            marks[panel] = static_cast<std::uint32_t>(n + 1);
                            plan.supportPanels.push_back(panel);
                        }
                    }
                }
            if (plan.supportPanels.size() > UINT32_MAX)
                throw std::overflow_error("Support panel reference index overflow.");
            std::sort(plan.supportPanels.begin() + first, plan.supportPanels.end());
            plan.nodePanelOffsets[n + 1] = static_cast<std::uint32_t>(plan.supportPanels.size());
        }
    }
    const galerkin::Space &space;
    const SupportTree &tree;
    std::size_t maximumBytes;
    std::stop_token stop;
    SupportNearPlan plan;
    std::vector<std::uint32_t> offsets, indices, cursor, columns, marks;
};
} // namespace

std::size_t SupportNearPlan::allocatedBytes() const noexcept {
    return bytes(pressure.rowOffsets) + bytes(pressure.columns) + bytes(flux.rowOffsets) +
           bytes(flux.columns) + bytes(activeFaceIds) + bytes(panelToActiveFaceColumn) +
           bytes(nodePanelOffsets) + bytes(supportPanels) + bytes(fluxBlockOffsets) +
           bytes(fluxBlockIndices);
}
SupportNearPlan makeSupportNearPlan(const galerkin::Space &space, const SupportTree &tree,
                                    std::span<const std::uint32_t> faces,
                                    std::size_t maximumHostBytes, std::stop_token stop) {
    return Builder(space, tree, maximumHostBytes, stop).build(faces);
}
SupportNearJobIterator::SupportNearJobIterator(const galerkin::Space &space,
                                               const SupportTree &tree, const SupportNearPlan &plan,
                                               std::size_t maximumWorkspaceBytes)
    : space_(space), tree_(tree), plan_(plan) {
    if (plan.panelToActiveFaceColumn.size() != space.panels.size() ||
        plan.nodePanelOffsets.size() != tree.nodes.size() + 1 ||
        plan.fluxBlockOffsets.size() != tree.nodes.size() + 1 ||
        plan.fluxBlockIndices.size() != tree.pressureFlux.nearInteractions.size() ||
        space.panels.size() > UINT32_MAX)
        throw std::invalid_argument("Support near iterator plan and geometry disagree.");
    constexpr auto bytesPerPanel = 2 * sizeof(std::uint32_t) + sizeof(std::uint16_t);
    if (space.panels.size() > maximumWorkspaceBytes / bytesPerPanel)
        throw std::runtime_error(
            "Support near iterator workspace budget exceeded before allocation.");
    stamps_.resize(space.panels.size());
    candidates_.resize(space.panels.size());
    masks_.resize(space.panels.size());
    if (workspaceBytes() > maximumWorkspaceBytes)
        throw std::runtime_error("Support near iterator actual workspace exceeds budget.");
}
bool SupportNearJobIterator::done() const noexcept {
    return done_;
}
std::size_t SupportNearJobIterator::workspaceBytes() const noexcept {
    return bytes(stamps_) + bytes(candidates_) + bytes(masks_);
}
std::uint16_t &SupportNearJobIterator::mask(std::uint32_t source) {
    const auto epoch = static_cast<std::uint32_t>(testPanel_ + 1);
    if (stamps_[source] != epoch) {
        stamps_[source] = epoch;
        masks_[source] = 0;
        candidates_[count_++] = source;
    }
    return masks_[source];
}
void SupportNearJobIterator::prepare(std::stop_token stop) {
    count_ = cursor_ = 0;
    const auto &target = space_.panels[testPanel_];
    // A CSR already encodes the exact P1 near ownership. Expanding only each
    // source vertex's incident panels unions the nine local basis slots without
    // searching or storing global panel pairs. Legitimate local addends survive.
    for (unsigned i = 0; i < 3; ++i) {
        cancelled(stop);
        const auto vertex = target.nodes[i];
        for (auto j = plan_.pressure.rowOffsets[vertex];
             j < plan_.pressure.rowOffsets[std::size_t(vertex) + 1]; ++j) {
            if ((j & 4095u) == 0)
                cancelled(stop);
            const auto sourceVertex = plan_.pressure.columns[j];
            for (auto k = space_.incidenceOffsets[sourceVertex];
                 k < space_.incidenceOffsets[std::size_t(sourceVertex) + 1]; ++k) {
                const auto incidence = space_.incidences[k];
                mask(incidence.panel) |=
                    static_cast<std::uint16_t>(1u << (3 * i + incidence.local));
            }
        }
        // B's face CSR alone cannot identify which panels of a face are near.
        // Use the exact leaf block ownership and omit inactive flux columns.
        const auto leaf = tree_.pressureLeaf[vertex];
        for (auto b = plan_.fluxBlockOffsets[leaf];
             b < plan_.fluxBlockOffsets[std::size_t(leaf) + 1]; ++b) {
            const auto &source =
                tree_.nodes[tree_.pressureFlux.nearInteractions[plan_.fluxBlockIndices[b]].source];
            for (auto k = source.firstFlux; k < std::size_t(source.firstFlux) + source.fluxCount;
                 ++k) {
                if ((k & 4095u) == 0)
                    cancelled(stop);
                const auto panel = tree_.fluxOrder[k];
                if (plan_.panelToActiveFaceColumn[panel] != UINT32_MAX)
                    mask(panel) |= static_cast<std::uint16_t>(1u << (9 + i));
            }
        }
    }
    prepared_ = true;
}
std::size_t SupportNearJobIterator::next(std::span<SupportNearJob> output, std::size_t budget,
                                         std::stop_token stop) {
    cancelled(stop);
    if (output.empty())
        throw std::invalid_argument("Support near job batch is empty.");
    if (output.size() > budget / sizeof(SupportNearJob))
        throw std::runtime_error("Support near job buffer budget exceeded.");
    std::size_t written = 0;
    while (written < output.size() && !done_) {
        cancelled(stop);
        if (testPanel_ == space_.panels.size()) {
            done_ = true;
            break;
        }
        if (!prepared_) {
            try {
                prepare(stop);
            } catch (...) {
                // A stopped partial panel must be restartable: its epoch marks
                // cannot hide candidates after the incomplete list is reset.
                std::fill(stamps_.begin(), stamps_.end(), 0u);
                throw;
            }
        }
        if (cursor_ == count_) {
            ++testPanel_;
            prepared_ = false;
            continue;
        }
        const auto source = candidates_[cursor_++];
        output[written++] = {static_cast<std::uint32_t>(testPanel_), source,
                             static_cast<std::uint16_t>(masks_[source] & 0x1ffu),
                             static_cast<std::uint8_t>(masks_[source] >> 9)};
    }
    return written;
}
} // namespace hornsculpt::wideband
