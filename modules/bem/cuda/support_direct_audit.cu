#include "galerkin_pair_moments.cuh"
#include "support_direct_audit.hpp"
#include "wideband_runtime.cuh"
#include <algorithm>

namespace hornsculpt::wideband {
namespace {
using namespace detail;
using cuda_detail::Panel;
struct LeafAdjacency {
    const std::uint32_t *offsets{}, *columns{};
};
__device__ bool contains(LeafAdjacency graph, unsigned row, unsigned column) {
    auto begin = graph.offsets[row], end = graph.offsets[row + 1];
    const auto last = end;
    while (begin < end) {
        const auto middle = begin + (end - begin) / 2;
        if (graph.columns[middle] < column)
            begin = middle + 1;
        else
            end = middle;
    }
    return begin < last && graph.columns[begin] == column;
}
template <bool fluxMode>
__global__ void directFar(const Panel *panels, unsigned panelCount, const unsigned *pressureLeaf,
                          const unsigned *fluxLeaf, const unsigned *faceColumns,
                          LeafAdjacency adjacency, const galerkin::TriangleQuadraturePoint *rule,
                          unsigned ruleCount, double k, const Z *input, Z *output,
                          std::size_t first, int *error) {
    constexpr unsigned sourceBatch = 256;
    const auto chunks = (panelCount + sourceBatch - 1) / sourceBatch;
    const auto job = first + blockIdx.x;
    const unsigned testId = static_cast<unsigned>(job / chunks),
                   chunk = static_cast<unsigned>(job % chunks);
    const auto &test = panels[testId];
    Z sum[3]{};
    for (unsigned trialId = chunk * sourceBatch + threadIdx.x;
         trialId < panelCount && trialId < (chunk + 1) * sourceBatch; trialId += blockDim.x) {
        const auto &trial = panels[trialId];
        unsigned mask = 0;
        if constexpr (fluxMode) {
            if (faceColumns[trialId] == UINT32_MAX)
                continue;
            for (unsigned i = 0; i < 3; ++i)
                if (!contains(adjacency, pressureLeaf[test.nodes[i]], fluxLeaf[trialId]))
                    mask |= 1u << i;
        } else {
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j)
                    if (!contains(adjacency, pressureLeaf[test.nodes[i]],
                                  pressureLeaf[trial.nodes[j]]))
                        mask |= 1u << (3 * i + j);
        }
        if (!mask)
            continue;
        // Higher rules contain up to4096 samples: promote each canonical
        // mixed sample before accumulation instead of a long FP32 sum.
        const auto moments = cuda_detail::integrateRegularPairThread<true>(
            test, trial, rule, ruleCount, k, error, 0, 1);
        for (unsigned i = 0; i < 3; ++i) {
            if constexpr (fluxMode) {
                if (mask & (1u << i)) {
                    const auto a = cuda_detail::fluxCoefficient(moments, test, i, k, false),
                               b = input[faceColumns[trialId]];
                    sum[i].x += a.x * b.x - a.y * b.y;
                    sum[i].y += a.x * b.y + a.y * b.x;
                }
            } else {
                for (unsigned j = 0; j < 3; ++j)
                    if (mask & (1u << (3 * i + j))) {
                        const auto a = cuda_detail::pressureCoefficient(moments, test, trial, i, j,
                                                                        k),
                                   b = input[trial.nodes[j]];
                        sum[i].x += a.x * b.x - a.y * b.y;
                        sum[i].y += a.x * b.y + a.y * b.x;
                    }
            }
        }
    }
    for (unsigned i = 0; i < 3; ++i) {
        for (unsigned step = 16; step; step /= 2) {
            sum[i].x += __shfl_down_sync(0xffffffffu, sum[i].x, step);
            sum[i].y += __shfl_down_sync(0xffffffffu, sum[i].y, step);
        }
        if (!(threadIdx.x & 31)) {
            atomicAdd(&output[test.nodes[i]].x, sum[i].x);
            atomicAdd(&output[test.nodes[i]].y, sum[i].y);
        }
    }
}
Panel packedPanel(const galerkin::Panel &source) {
    Panel value{};
    for (unsigned j = 0; j < 3; ++j) {
        value.vertices[j] = source.vertices[j];
        value.nodes[j] = source.nodes[j];
        value.curls[j] = source.curls[j];
    }
    value.normal = source.normal;
    value.area = source.area;
    value.centroid = source.centroid;
    value.edge = source.maximumEdgeLength;
    return value;
}
struct Adjacency {
    DeviceBuffer<std::uint32_t> offsets, columns;
    std::size_t allocatedBytes{};
    void build(const SupportPartition &partition, std::size_t nodes, cudaStream_t stream) {
        std::vector<std::vector<std::uint32_t>> rows(nodes);
        for (const auto pair : partition.nearInteractions)
            rows.at(pair.target).push_back(pair.source);
        std::vector<std::uint32_t> rowOffsets{0}, entries;
        for (auto &row : rows) {
            std::sort(row.begin(), row.end());
            row.erase(std::unique(row.begin(), row.end()), row.end());
            entries.insert(entries.end(), row.begin(), row.end());
            if (entries.size() > UINT32_MAX)
                throw std::overflow_error("Direct audit leaf adjacency exceeds uint32 indexing.");
            rowOffsets.push_back(static_cast<std::uint32_t>(entries.size()));
        }
        offsets.allocate(rowOffsets.size());
        columns.allocate(entries.size());
        allocatedBytes = product(rowOffsets.size() + entries.size(), sizeof(std::uint32_t));
        offsets.upload(rowOffsets, stream);
        columns.upload(entries, stream);
        checked(cudaStreamSynchronize(stream), "Direct audit leaf adjacency upload");
    }
    LeafAdjacency view() const {
        return {offsets.data(), columns.data()};
    }
    std::size_t bytes() const {
        return allocatedBytes;
    }
};
} // namespace
class SupportDirectAudit::Implementation {
  public:
    Implementation(const galerkin::Space &space, const SupportTree &tree,
                   const SupportNearPlan &plan, std::shared_ptr<SupportNearOperator> operation,
                   double wavenumber, unsigned regularOrder, std::stop_token stop)
        : near(std::move(operation)), k(wavenumber) {
        cancelled(stop);
        if (!near || space.panels.empty() || space.panels.size() > maximumDirectAuditPanels ||
            space.vertices.size() != tree.pressureLeaf.size() ||
            space.panels.size() != tree.fluxLeaf.size() ||
            space.panels.size() != plan.panelToActiveFaceColumn.size() || !std::isfinite(k) ||
            k <= 0 ||
            (regularOrder != 2 && regularOrder != 4 && regularOrder != 6 && regularOrder != 8 &&
             regularOrder != 12 && regularOrder != 14))
            throw std::invalid_argument(
                "Direct audit requires a valid shared partition and at most 20000 panels.");
        stream = near->stream();
        panelCount = static_cast<unsigned>(space.panels.size());
        pressureCount = space.vertices.size();
        for (const auto column : plan.panelToActiveFaceColumn)
            if (column != UINT32_MAX && column >= near->faceFluxMatrix().columnCount)
                throw std::invalid_argument(
                    "Direct audit face column exceeds the shared RHS basis.");
        for (const auto leaf : tree.pressureLeaf)
            if (leaf >= tree.nodes.size())
                throw std::invalid_argument(
                    "Direct audit pressure leaf is outside the shared tree.");
        for (const auto leaf : tree.fluxLeaf)
            if (leaf >= tree.nodes.size())
                throw std::invalid_argument("Direct audit flux leaf is outside the shared tree.");
        const auto quadrature = galerkin::triangleRule(regularOrder);
        ruleCount = static_cast<unsigned>(quadrature.size());
        allocatedBytes =
            sum(product(space.panels.size(), sizeof(Panel)),
                sum(product(tree.pressureLeaf.size() + tree.fluxLeaf.size() +
                                plan.panelToActiveFaceColumn.size(),
                            sizeof(unsigned)),
                    sum(product(quadrature.size(), sizeof(galerkin::TriangleQuadraturePoint)),
                        sizeof(int))));
        const auto graphEstimate =
            product(2 * (tree.nodes.size() + 1) + tree.pressurePressure.nearInteractions.size() +
                        tree.pressureFlux.nearInteractions.size(),
                    sizeof(std::uint32_t));
        std::size_t free{}, total{};
        checked(cudaMemGetInfo(&free, &total), "Direct audit memory query");
        if (sum(allocatedBytes, graphEstimate) > free - free / 5)
            throw std::runtime_error("Direct audit buffers exceed available GPU memory budget.");
        std::vector<Panel> host;
        host.reserve(panelCount);
        for (const auto &panel : space.panels)
            host.push_back(packedPanel(panel));
        panels.allocate(host.size());
        panels.upload(host, stream);
        pressureLeaf.allocate(tree.pressureLeaf.size());
        pressureLeaf.upload(tree.pressureLeaf, stream);
        fluxLeaf.allocate(tree.fluxLeaf.size());
        fluxLeaf.upload(tree.fluxLeaf, stream);
        faceColumns.allocate(plan.panelToActiveFaceColumn.size());
        faceColumns.upload(plan.panelToActiveFaceColumn, stream);
        rule.allocate(quadrature.size());
        rule.upload(quadrature, stream);
        error.allocate(1);
        pressureNear.build(tree.pressurePressure, tree.nodes.size(), stream);
        fluxNear.build(tree.pressureFlux, tree.nodes.size(), stream);
        allocatedBytes = sum(allocatedBytes, sum(pressureNear.bytes(), fluxNear.bytes()));
        checked(cudaStreamSynchronize(stream), "Direct audit preparation");
        cancelled(stop);
    }
    template <bool fluxMode> void apply(const Z *input, Z *output, std::stop_token stop) {
        cancelled(stop);
        if (!input || !output || input == output)
            throw std::invalid_argument("Direct audit needs distinct nonnull device vectors.");
        if constexpr (fluxMode)
            near->applyFaceFluxDevice(input, output, stop);
        else
            near->applyPressureDevice(input, output, stop);
        checked(cudaMemsetAsync(error.data(), 0, sizeof(int), stream), "Direct audit status reset");
        const auto jobs = product(panelCount, (panelCount + 255) / 256);
        for (std::size_t first = 0; first < jobs; first += 2048) {
            cancelled(stop);
            const auto count = static_cast<unsigned>(std::min(std::size_t{2048}, jobs - first));
            directFar<fluxMode><<<count, threads, 0, stream>>>(
                panels.data(), panelCount, pressureLeaf.data(), fluxLeaf.data(), faceColumns.data(),
                fluxMode ? fluxNear.view() : pressureNear.view(), rule.data(), ruleCount, k, input,
                output, first, error.data());
            checked(cudaGetLastError(), "Direct audit far launch");
            int status{};
            checked(cudaMemcpyAsync(&status, error.data(), sizeof(status), cudaMemcpyDeviceToHost,
                                    stream),
                    "Direct audit status download");
            checked(cudaStreamSynchronize(stream), "Direct audit bounded batch");
            if (status)
                throw std::runtime_error("Direct Galerkin audit integration failed with status " +
                                         std::to_string(status) + ".");
        }
        verifySpectrum(output, pressureCount, error.data(), stream, stop);
    }
    std::size_t bytes() const {
        return allocatedBytes;
    }
    std::shared_ptr<SupportNearOperator> near;
    cudaStream_t stream{};
    double k{};
    unsigned panelCount{}, ruleCount{};
    std::size_t pressureCount{}, allocatedBytes{};
    DeviceBuffer<Panel> panels;
    DeviceBuffer<unsigned> pressureLeaf, fluxLeaf, faceColumns;
    DeviceBuffer<galerkin::TriangleQuadraturePoint> rule;
    DeviceBuffer<int> error;
    Adjacency pressureNear, fluxNear;
};
SupportDirectAudit::SupportDirectAudit(const galerkin::Space &space, const SupportTree &tree,
                                       const SupportNearPlan &plan,
                                       std::shared_ptr<SupportNearOperator> near, double k,
                                       unsigned order, std::stop_token stop)
    : implementation_(
          std::make_unique<Implementation>(space, tree, plan, std::move(near), k, order, stop)) {}
SupportDirectAudit::~SupportDirectAudit() = default;
void SupportDirectAudit::applyPressureDevice(const Z *in, Z *out, std::stop_token stop) {
    implementation_->apply<false>(in, out, stop);
}
void SupportDirectAudit::applyFaceFluxDevice(const Z *in, Z *out, std::stop_token stop) {
    implementation_->apply<true>(in, out, stop);
}
std::size_t SupportDirectAudit::ownDeviceBytes() const noexcept {
    return implementation_->bytes();
}
} // namespace hornsculpt::wideband
