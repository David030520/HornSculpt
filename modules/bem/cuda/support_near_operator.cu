#include "../galerkin/regular_qualification.hpp"
#include "adaptive_patch_estimator.hpp"
#include "galerkin_pair_moments.cuh"
#include "regular_pair_contraction.cuh"
#include "support_near_operator.hpp"
#include "wideband_runtime.cuh"
#include <algorithm>
#include <array>
#include <chrono>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>

namespace hornsculpt::wideband {
namespace {
using namespace detail;
using cuda_detail::PairMoments;
using cuda_detail::Panel;
using Clock = std::chrono::steady_clock;
double milliseconds(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
struct MaskedQuadratureJob {
    cuda_detail::PairIntegrationJob integration;
    std::uint32_t faceColumn{UINT32_MAX};
    std::uint16_t pressureMask{};
    std::uint8_t fluxMask{};
    std::uint8_t regularRuleIndex{};
};
struct RegularRulesView {
    const galerkin::TriangleQuadraturePoint *points{};
    unsigned offsets[3]{}, counts[3]{};
};
struct MutableCsr {
    const std::uint32_t *offsets{}, *columns{};
    Z *values{};
};
class KernelTimer {
  public:
    KernelTimer() {
        checked(cudaEventCreate(&start), "Near kernel start event");
        try {
            checked(cudaEventCreate(&end), "Near kernel end event");
        } catch (...) {
            cudaEventDestroy(start);
            throw;
        }
    }
    ~KernelTimer() {
        cudaEventDestroy(start);
        cudaEventDestroy(end);
    }
    cudaEvent_t start{}, end{};
};
template <class T> class PinnedArray {
  public:
    static_assert(std::is_nothrow_default_constructible_v<T> &&
                  std::is_trivially_destructible_v<T>);
    ~PinnedArray() {
        if (data_)
            cudaFreeHost(data_);
    }
    PinnedArray() = default;
    PinnedArray(const PinnedArray &) = delete;
    PinnedArray &operator=(const PinnedArray &) = delete;
    void allocate(std::size_t count) {
        if (data_)
            throw std::logic_error("Near pinned staging is already allocated.");
        if (!count)
            return;
        void *storage{};
        checked(cudaHostAlloc(&storage, product(count, sizeof(T)), cudaHostAllocDefault),
                "Near pinned staging allocation");
        data_ = static_cast<T *>(storage);
        count_ = count;
        std::uninitialized_value_construct_n(data_, count_);
    }
    T *data() const noexcept {
        return data_;
    }
    std::size_t size() const noexcept {
        return count_;
    }
    std::size_t bytes() const {
        return product(count_, sizeof(T));
    }

  private:
    T *data_{};
    std::size_t count_{};
};
class CompletionEvent {
  public:
    CompletionEvent() {
        checked(cudaEventCreateWithFlags(&value_, cudaEventDisableTiming | cudaEventBlockingSync),
                "Near batch completion event");
    }
    ~CompletionEvent() {
        if (value_)
            cudaEventDestroy(value_);
    }
    CompletionEvent(const CompletionEvent &) = delete;
    CompletionEvent &operator=(const CompletionEvent &) = delete;
    operator cudaEvent_t() const noexcept {
        return value_;
    }

  private:
    cudaEvent_t value_{};
};
struct ConstructionBatch {
    PinnedArray<MaskedQuadratureJob> jobs;
    PinnedArray<galerkin::PairQuadraturePatch> patches;
    PinnedArray<int> status;
    KernelTimer regularTimer, heavyTimer;
    CompletionEvent completion;
    std::size_t jobCount{}, patchCount{}, samples{}, ordinal{};
    unsigned regularJobs{}, heavyJobs{};
    bool preparationStarted{};
    void clear() noexcept {
        jobCount = patchCount = samples = 0;
        regularJobs = heavyJobs = 0;
        preparationStarted = false;
    }
};
// One active device batch and one prepared pinned host batch. The borrowed
// stream preserves all existing kernel/scatter order; no concurrent CSR writes
// or estimator stream are introduced. Declare this owner after every borrowed
// panel/rule/estimator buffer so its destructor drains before they can die.
class ConstructionPipeline {
  public:
    ConstructionPipeline(std::size_t maximumJobs, std::size_t maximumPatches,
                         std::size_t maximumPinnedBytes, cudaStream_t stream,
                         SupportNearConstructionStats &statistics,
                         const std::function<void(const SupportNearConstructionStats &)> &progress)
        : stream_(stream), stats_(statistics), progress_(progress) {
        const auto slotBytes =
            sum(sum(product(maximumJobs, sizeof(MaskedQuadratureJob)),
                    product(maximumPatches, sizeof(galerkin::PairQuadraturePatch))),
                sizeof(int));
        if (product(slots_.size(), slotBytes) > maximumPinnedBytes)
            throw std::runtime_error("Near pinned construction staging exceeds its host budget.");
        // Construction only establishes ownership; it never queues a transfer.
        // A throwing constructor therefore has no in-flight pinned references.
        for (auto &slot : slots_) {
            slot.jobs.allocate(maximumJobs);
            slot.patches.allocate(maximumPatches);
            slot.status.allocate(1);
            stats_.pinnedConstructionHostBytes =
                sum(stats_.pinnedConstructionHostBytes,
                    sum(sum(slot.jobs.bytes(), slot.patches.bytes()), slot.status.bytes()));
        }
        stats_.maximumPinnedConstructionHostBytes = maximumPinnedBytes;
        if (stats_.pinnedConstructionHostBytes > maximumPinnedBytes)
            throw std::runtime_error("Actual near pinned capacities exceed their host budget.");
        jobs_.allocate(maximumJobs);
        patches_.allocate(maximumPatches);
        stats_.deviceConstructionBatchBytes =
            sum(product(maximumJobs, sizeof(MaskedQuadratureJob)),
                product(maximumPatches, sizeof(galerkin::PairQuadraturePatch)));
    }
    ~ConstructionPipeline() {
        // Normal finish checks completion and status explicitly. During any
        // exception/cancellation, drain before slots_ or device buffers die.
        cudaStreamSynchronize(stream_);
    }
    ConstructionPipeline(const ConstructionPipeline &) = delete;
    ConstructionPipeline &operator=(const ConstructionPipeline &) = delete;
    ConstructionBatch &prepared() noexcept {
        return slots_[preparing_];
    }
    void beginPreparation() {
        auto &batch = prepared();
        if (batch.preparationStarted)
            return;
        batch.preparationStarted = true;
        if (active_) {
            const auto status = cudaEventQuery(slots_[*active_].completion);
            if (status == cudaErrorNotReady)
                ++stats_.preparationsStartedWhileAssemblyPending;
            else
                checked(status, "Near pending assembly observation");
        }
    }
    void observeResidentBatches() {
        stats_.maximumResidentBatches =
            std::max(stats_.maximumResidentBatches,
                     std::size_t(active_.has_value()) + std::size_t(prepared().jobCount != 0));
    }
    template <class Enqueue> void submit(Enqueue enqueue, int *error, std::stop_token stop) {
        auto &batch = prepared();
        if (!batch.jobCount)
            return;
        cancelled(stop);
        const auto heavyBegin = std::partition(
            batch.jobs.data(), batch.jobs.data() + batch.jobCount, [](const auto &job) {
                return job.integration.permutation.topology == galerkin::PairTopology::disjoint &&
                       job.integration.patchCount == 0;
            });
        batch.regularJobs = static_cast<unsigned>(heavyBegin - batch.jobs.data());
        batch.heavyJobs = static_cast<unsigned>(batch.jobCount - batch.regularJobs);
        // No shared device buffer is overwritten until the previous completion
        // and its status snapshot have both been checked in submission order.
        retire(stop);
        const auto start = Clock::now();
        checked(cudaMemcpyAsync(jobs_.data(), batch.jobs.data(),
                                product(batch.jobCount, sizeof(MaskedQuadratureJob)),
                                cudaMemcpyHostToDevice, stream_),
                "Near pinned job upload");
        if (batch.patchCount)
            checked(
                cudaMemcpyAsync(patches_.data(), batch.patches.data(),
                                product(batch.patchCount, sizeof(galerkin::PairQuadraturePatch)),
                                cudaMemcpyHostToDevice, stream_),
                "Near pinned patch upload");
        enqueue(batch, jobs_.data(), patches_.data());
        // The assembly status remains sticky throughout construction. The
        // estimator owns a separate status buffer on this same stream.
        checked(cudaMemcpyAsync(batch.status.data(), error, sizeof(int), cudaMemcpyDeviceToHost,
                                stream_),
                "Near pinned batch status snapshot");
        checked(cudaEventRecord(batch.completion, stream_), "Near batch completion record");
        batch.ordinal = ++stats_.submittedBatches;
        active_ = preparing_;
        preparing_ = 1 - preparing_;
        prepared().clear();
        stats_.maximumResidentBatches = std::max(stats_.maximumResidentBatches, std::size_t{1});
        const auto elapsed = milliseconds(start);
        stats_.batchSubmissionMilliseconds += elapsed;
        stats_.batchUploadAndSyncWallMilliseconds += elapsed;
        if (progress_)
            progress_(stats_);
        cancelled(stop);
    }
    void finish(std::stop_token stop) {
        retire(stop);
        checked(cudaStreamSynchronize(stream_), "Near construction final drain");
        cancelled(stop);
    }

  private:
    void retire(std::stop_token stop) {
        if (!active_)
            return;
        cancelled(stop);
        auto &batch = slots_[*active_];
        const auto start = Clock::now();
        checked(cudaEventSynchronize(batch.completion), "Near batch completion wait");
        const auto elapsed = milliseconds(start);
        stats_.batchRetirementMilliseconds += elapsed;
        stats_.batchUploadAndSyncWallMilliseconds += elapsed;
        if (*batch.status.data())
            throw std::runtime_error("Support near CUDA integration/scatter failed in batch " +
                                     std::to_string(batch.ordinal) + ", status " +
                                     std::to_string(*batch.status.data()) + ".");
        const auto record = [&](const KernelTimer &timer) {
            float kernelMilliseconds{};
            checked(cudaEventElapsedTime(&kernelMilliseconds, timer.start, timer.end),
                    "Near completed batch kernel time");
            stats_.maximumKernelMilliseconds =
                std::max(stats_.maximumKernelMilliseconds, double(kernelMilliseconds));
            stats_.kernelMilliseconds += kernelMilliseconds;
        };
        if (batch.regularJobs)
            record(batch.regularTimer);
        if (batch.heavyJobs)
            record(batch.heavyTimer);
        ++stats_.batches;
        stats_.quadratureJobs += batch.jobCount;
        stats_.patches += batch.patchCount;
        stats_.maximumBatchJobs = std::max(stats_.maximumBatchJobs, batch.jobCount);
        stats_.maximumBatchPatches = std::max(stats_.maximumBatchPatches, batch.patchCount);
        stats_.maximumBatchGreenSamples = std::max(stats_.maximumBatchGreenSamples, batch.samples);
        active_.reset();
        if (progress_)
            progress_(stats_);
        cancelled(stop);
    }
    cudaStream_t stream_{};
    SupportNearConstructionStats &stats_;
    const std::function<void(const SupportNearConstructionStats &)> &progress_;
    std::array<ConstructionBatch, 2> slots_;
    DeviceBuffer<MaskedQuadratureJob> jobs_;
    DeviceBuffer<galerkin::PairQuadraturePatch> patches_;
    std::optional<std::size_t> active_;
    std::size_t preparing_{};
};
__device__ void addSparse(MutableCsr matrix, unsigned row, unsigned column, Z value, int *error) {
    auto first = matrix.offsets[row], last = matrix.offsets[row + 1];
    const auto end = last;
    while (first < last) {
        const auto middle = first + (last - first) / 2;
        if (matrix.columns[middle] < column)
            first = middle + 1;
        else
            last = middle;
    }
    if (first == end || matrix.columns[first] != column) {
        atomicOr(error, 4);
        return;
    }
    if (!isfinite(value.x) || !isfinite(value.y)) {
        atomicOr(error, 2);
        return;
    }
    atomicAdd(&matrix.values[first].x, value.x);
    atomicAdd(&matrix.values[first].y, value.y);
}
__global__ void assemble(const Panel *panels, const MaskedQuadratureJob *jobs,
                         cuda_detail::PairQuadratureView quadrature, double k, MutableCsr a,
                         MutableCsr b, int *error) {
    const auto job = jobs[blockIdx.x];
    const auto &test = panels[job.integration.test], &trial = panels[job.integration.trial];
    auto moments =
        cuda_detail::integratePairThread(test, trial, job.integration, quadrature, k, error);
    cuda_detail::reduceWarp(moments);
    __shared__ PairMoments partial[4];
    if ((threadIdx.x & 31) == 0)
        partial[threadIdx.x / 32] = moments;
    __syncthreads();
    if (threadIdx.x < 32) {
        PairMoments total{};
        if (threadIdx.x < 4)
            total = partial[threadIdx.x];
        cuda_detail::reduceWarp(total);
        if (threadIdx.x == 0)
            for (unsigned i = 0; i < 3; ++i) {
                const bool coincident = job.integration.test == job.integration.trial;
                if (job.fluxMask & (1u << i))
                    addSparse(b, test.nodes[i], job.faceColumn,
                              cuda_detail::fluxCoefficient(total, test, i, k, coincident), error);
                for (unsigned j = 0; j < 3; ++j)
                    if (job.pressureMask & (1u << (3 * i + j))) {
                        auto value = cuda_detail::pressureCoefficient(total, test, trial, i, j, k);
                        if (coincident)
                            value.x += test.area * (i == j ? 2. : 1.) / 24.;
                        addSparse(a, test.nodes[i], trial.nodes[j], value, error);
                    }
            }
    }
}
// Four independent pairs share a block, with private tensor-contraction scratch
// per warp. Inactive trailing warps never participate in another pair.
__global__ void assembleRegularWarp(const Panel *panels, const MaskedQuadratureJob *jobs,
                                    unsigned jobCount, RegularRulesView rules, double k,
                                    MutableCsr a, MutableCsr b, int *error) {
    const unsigned index = blockIdx.x * 4 + threadIdx.x / 32;
    if (index >= jobCount)
        return;
    const auto job = jobs[index];
    const auto &test = panels[job.integration.test], &trial = panels[job.integration.trial];
    const auto *regular = rules.points + rules.offsets[job.regularRuleIndex];
    const auto regularCount = rules.counts[job.regularRuleIndex];
    __shared__ cuda_detail::RegularPairScratch scratch[4];
    const auto value = cuda_detail::contractRegularPairWarp(test, trial, regular, regularCount, k,
                                                            scratch[threadIdx.x / 32], error);
    const unsigned lane = threadIdx.x & 31u;
    if (lane < 9 && (job.pressureMask & (1u << lane)))
        addSparse(a, test.nodes[lane / 3], trial.nodes[lane % 3], value, error);
    else if (lane >= 9 && lane < 12 && (job.fluxMask & (1u << (lane - 9))))
        addSparse(b, test.nodes[lane - 9], job.faceColumn, value, error);
}

__global__ void applyCsr(SupportDeviceCsr matrix, const Z *input, Z *output, int *error) {
    const auto row = blockIdx.x * 4 + threadIdx.x / 32;
    if (row >= matrix.rows)
        return;
    Z value{};
    for (auto j = matrix.rowOffsets[row] + (threadIdx.x & 31); j < matrix.rowOffsets[row + 1];
         j += 32) {
        const auto x = input[matrix.columns[j]], a = matrix.values[j];
        value.x += a.x * x.x - a.y * x.y;
        value.y += a.x * x.y + a.y * x.x;
    }
    for (unsigned shift = 16; shift; shift /= 2) {
        value.x += __shfl_down_sync(0xffffffffu, value.x, shift);
        value.y += __shfl_down_sync(0xffffffffu, value.y, shift);
    }
    if ((threadIdx.x & 31) == 0) {
        if (!isfinite(value.x) || !isfinite(value.y))
            atomicOr(error, 2);
        output[row] = value;
    }
}
struct Matrix {
    DeviceBuffer<std::uint32_t> offsets, columns;
    DeviceBuffer<Z> values;
    unsigned rows{}, columnCount{}, nonzeros{};
    void create(const SupportCsr &host, unsigned columnSize, cudaStream_t stream) {
        rows = static_cast<unsigned>(host.rowOffsets.size() - 1);
        columnCount = columnSize;
        nonzeros = static_cast<unsigned>(host.columns.size());
        offsets.allocate(host.rowOffsets.size());
        offsets.upload(host.rowOffsets, stream);
        columns.allocate(host.columns.size());
        columns.upload(host.columns, stream);
        values.allocate(host.columns.size());
        if (nonzeros)
            checked(cudaMemsetAsync(values.data(), 0, product(nonzeros, sizeof(Z)), stream),
                    "Near matrix clear");
    }
    SupportDeviceCsr view() const {
        return {offsets.data(), columns.data(), values.data(), rows, columnCount, nonzeros};
    }
    MutableCsr writable() const {
        return {offsets.data(), columns.data(), values.data()};
    }
};
Panel panel(const galerkin::Panel &source) {
    Panel value{};
    value.centroid = source.centroid;
    value.normal = source.normal;
    value.area = source.area;
    value.edge = source.maximumEdgeLength;
    for (unsigned i = 0; i < 3; ++i) {
        value.vertices[i] = source.vertices[i];
        value.curls[i] = source.curls[i];
        value.nodes[i] = source.nodes[i];
    }
    return value;
}
std::size_t matrixBytes(const SupportCsr &csr) {
    return sum(product(csr.rowOffsets.size(), 4), product(csr.columns.size(), 4 + sizeof(Z)));
}
void validate(const SupportCsr &csr, std::size_t rows, std::size_t columns) {
    if (csr.rowOffsets.size() != rows + 1 || csr.rowOffsets.front() != 0 ||
        csr.rowOffsets.back() != csr.columns.size() || csr.columns.size() > UINT32_MAX)
        throw std::invalid_argument("Near CSR dimensions are invalid.");
    for (std::size_t i = 0; i < rows; ++i) {
        if (csr.rowOffsets[i] > csr.rowOffsets[i + 1] || csr.rowOffsets[i + 1] > csr.columns.size())
            throw std::invalid_argument("Near CSR row offsets are not monotone.");
        const auto first = csr.columns.begin() + csr.rowOffsets[i],
                   last = csr.columns.begin() + csr.rowOffsets[i + 1];
        if (!std::is_sorted(first, last) || std::adjacent_find(first, last) != last ||
            (first != last && *(last - 1) >= columns))
            throw std::invalid_argument("Near CSR columns are invalid.");
    }
}
} // namespace

class SupportNearOperator::Implementation {
  public:
    Implementation(const galerkin::Space &space, const SupportTree &tree,
                   const SupportNearPlan &plan, double k, SupportNearOperatorSettings settings,
                   std::stop_token stop, cudaStream_t borrowed,
                   const std::function<void(const SupportNearConstructionStats &)> &progress) {
        cancelled(stop);
        const auto &q = settings.quadrature;
        if (space.vertices.empty() || !(k > 0) || !std::isfinite(k) ||
            (q.regularOrder != 2 && q.regularOrder != 4) || q.singularOrder < 2 ||
            q.singularOrder > 20 || !std::isfinite(q.nearDistanceRatio) ||
            q.nearDistanceRatio < 0 || !std::isfinite(q.relativeTolerance) ||
            q.relativeTolerance <= 0 || q.maximumSubdivisionDepth > 16 ||
            q.nearDistanceRatio > tree.nearDistanceRatio || !settings.maximumBatchJobs ||
            !settings.maximumBatchGreenSamples)
            throw std::invalid_argument(
                "Near operator physical, quadrature or support-distance contract is invalid.");
        validate(plan.pressure, space.vertices.size(), space.vertices.size());
        validate(plan.flux, space.vertices.size(), plan.activeFaceIds.size());
        const auto minimumBatch =
            sum(sum(sizeof(MaskedQuadratureJob), sizeof(galerkin::PairQuadraturePatch)),
                sizeof(SupportNearJob));
        if (settings.maximumConstructionBatchBytes < minimumBatch)
            throw std::invalid_argument("Near construction batch cannot hold one job and patch.");
        const auto maskBatchCount = std::min<std::size_t>(
            4096, std::max<std::size_t>(1, (settings.maximumConstructionBatchBytes - minimumBatch) /
                                               4 / sizeof(SupportNearJob)));
        const auto maskBatchBytes = product(maskBatchCount, sizeof(SupportNearJob));
        const auto payloadBudget = settings.maximumConstructionBatchBytes - maskBatchBytes;
        const auto maximumJobs = std::min<std::size_t>(
            settings.maximumBatchJobs,
            (payloadBudget - sizeof(galerkin::PairQuadraturePatch)) / sizeof(MaskedQuadratureJob));
        const auto maximumPatches =
            (payloadBudget - product(maximumJobs, sizeof(MaskedQuadratureJob))) /
            sizeof(galerkin::PairQuadraturePatch);
        std::vector<galerkin::TriangleQuadraturePoint> regular;
        RegularRulesView ruleView;
        for (unsigned i = 0; i < 3; ++i) {
            const auto rule = galerkin::triangleRule(i == 0 ? 4 : i == 1 ? 6 : 8);
            ruleView.offsets[i] = static_cast<unsigned>(regular.size());
            ruleView.counts[i] = static_cast<unsigned>(rule.size());
            regular.insert(regular.end(), rule.begin(), rule.end());
        }
        auto nearRules = galerkin::triangleRule(8);
        stats.classificationWorkspaceBytes =
            product(space.panels.size(), sizeof(galerkin::RegularPanelBounds));
        if (stats.classificationWorkspaceBytes > settings.maximumClassificationWorkspaceBytes)
            throw std::runtime_error("Near panel-qualification workspace exceeds its host budget.");
        auto bounds = std::make_unique<galerkin::RegularPanelBounds[]>(space.panels.size());
        for (std::size_t i = 0; i < space.panels.size(); ++i) {
            cancelled(stop);
            bounds[i] = galerkin::makeRegularPanelBounds(space.panels[i]);
        }
        const auto baseCount = static_cast<unsigned>(nearRules.size());
        const auto high = galerkin::triangleRule(14);
        nearRules.insert(nearRules.end(), high.begin(), high.end());
        const auto gauss = galerkin::gaussLegendre01(q.singularOrder);
        liveBytes =
            sum(sum(matrixBytes(plan.pressure), matrixBytes(plan.flux)),
                sum(product(2 * space.vertices.size() + plan.activeFaceIds.size(), sizeof(Z)),
                    sizeof(int)));
        const auto temporary =
            sum(sum(product(space.panels.size(), sizeof(Panel)),
                    product(maximumJobs, sizeof(MaskedQuadratureJob))),
                sum(product(maximumPatches, sizeof(galerkin::PairQuadraturePatch)),
                    sum(product(regular.size() + nearRules.size(),
                                sizeof(galerkin::TriangleQuadraturePoint)),
                        product(gauss.size(), sizeof(galerkin::GaussPoint)))));
        peakBytes = sum(liveBytes, temporary);
        std::size_t available{}, total{};
        checked(cudaMemGetInfo(&available, &total), "Near memory preflight");
        if (peakBytes > settings.maximumDeviceBytes || peakBytes > available - available / 5)
            throw std::runtime_error(
                "Support near device buffer budget exceeded before allocation.");
        if (borrowed)
            streamValue = borrowed;
        else {
            ownedStream = std::make_unique<Stream>();
            streamValue = *ownedStream;
        }
        AdaptivePairBatch adaptiveBatch(space, k, q, settings.maximumAdaptiveWorkspaceBytes,
                                        settings.maximumDeviceBytes - peakBytes,
                                        settings.maximumActiveAdaptivePairs, streamValue, stop);
        peakBytes = sum(peakBytes, adaptiveBatch.deviceBytes());
        if (peakBytes > settings.maximumDeviceBytes || peakBytes > available - available / 5)
            throw std::runtime_error(
                "Combined near and adaptive estimator device budget exceeded before allocation.");
        a.create(plan.pressure, static_cast<unsigned>(space.vertices.size()), streamValue);
        b.create(plan.flux, static_cast<unsigned>(plan.activeFaceIds.size()), streamValue);
        inputPressure.allocate(space.vertices.size());
        inputFace.allocate(plan.activeFaceIds.size());
        output.allocate(space.vertices.size());
        error.allocate(1);
        checked(cudaMemsetAsync(error.data(), 0, sizeof(int), streamValue), "Near error clear");
        DeviceBuffer<Panel> panels;
        panels.allocate(space.panels.size());
        std::array<Panel, 256> hostPanels{};
        for (std::size_t first = 0; first < space.panels.size(); first += hostPanels.size()) {
            cancelled(stop);
            const auto count = std::min(hostPanels.size(), space.panels.size() - first);
            for (std::size_t j = 0; j < count; ++j)
                hostPanels[j] = panel(space.panels[first + j]);
            checked(cudaMemcpyAsync(panels.data() + first, hostPanels.data(),
                                    product(count, sizeof(Panel)), cudaMemcpyHostToDevice,
                                    streamValue),
                    "Near panel upload");
            checked(cudaStreamSynchronize(streamValue), "Near panel upload completion");
        }
        DeviceBuffer<galerkin::TriangleQuadraturePoint> regularDevice, nearDevice;
        DeviceBuffer<galerkin::GaussPoint> gaussDevice;
        regularDevice.allocate(regular.size());
        regularDevice.upload(regular, streamValue);
        ruleView.points = regularDevice.data();
        nearDevice.allocate(nearRules.size());
        nearDevice.upload(nearRules, streamValue);
        gaussDevice.allocate(gauss.size());
        gaussDevice.upload(gauss, streamValue);
        checked(cudaStreamSynchronize(streamValue), "Near static quadrature upload completion");
        // Keep the existing batch partition unchanged. The explicit total host
        // staging bound scales with that same user-selected batch contract.
        ConstructionPipeline pipeline(maximumJobs, maximumPatches,
                                      product(2, settings.maximumConstructionBatchBytes),
                                      streamValue, stats, progress);
        // These are now allocated payloads, not a preflight prediction. CUDA
        // context/allocator bookkeeping and event internals are excluded.
        stats.liveDeviceBytes = liveBytes;
        stats.peakDeviceBytes = peakBytes;
        if (progress)
            progress(stats);
        const auto flush = [&] {
            if (!pipeline.prepared().jobCount)
                return;
            pipeline.submit(
                [&](ConstructionBatch &batch, const MaskedQuadratureJob *jobs,
                    const galerkin::PairQuadraturePatch *patches) {
                    if (batch.regularJobs) {
                        checked(cudaEventRecord(batch.regularTimer.start, streamValue),
                                "Near regular start");
                        assembleRegularWarp<<<(batch.regularJobs + 3) / 4, 128, 0, streamValue>>>(
                            panels.data(), jobs, batch.regularJobs, ruleView, k, a.writable(),
                            b.writable(), error.data());
                        checked(cudaGetLastError(), "Near regular warp launch");
                        checked(cudaEventRecord(batch.regularTimer.end, streamValue),
                                "Near regular end");
                    }
                    if (batch.heavyJobs) {
                        checked(cudaEventRecord(batch.heavyTimer.start, streamValue),
                                "Near kernel start");
                        assemble<<<batch.heavyJobs, 128, 0, streamValue>>>(
                            panels.data(), jobs + batch.regularJobs,
                            {regularDevice.data(), ruleView.counts[0], gaussDevice.data(),
                             q.singularOrder, patches, nearDevice.data(), baseCount,
                             static_cast<unsigned>(high.size())},
                            k, a.writable(), b.writable(), error.data());
                        checked(cudaGetLastError(), "Near singular/adaptive launch");
                        checked(cudaEventRecord(batch.heavyTimer.end, streamValue),
                                "Near kernel end");
                    }
                },
                error.data(), stop);
        };
        const auto append = [&](MaskedQuadratureJob job, const galerkin::PairQuadratureRule &rule) {
            if (rule.patches.empty()) {
                const auto sampleCount =
                    job.integration.permutation.topology == galerkin::PairTopology::disjoint
                        ? product(ruleView.counts[job.regularRuleIndex],
                                  ruleView.counts[job.regularRuleIndex])
                        : product(product(product(q.singularOrder, q.singularOrder),
                                          product(q.singularOrder, q.singularOrder)),
                                  galerkin::pairSectorCount(job.integration.permutation.topology));
                if (sampleCount > settings.maximumBatchGreenSamples)
                    throw std::runtime_error(
                        "One near quadrature job exceeds the Green sample budget.");
                if (pipeline.prepared().jobCount == maximumJobs ||
                    sampleCount > settings.maximumBatchGreenSamples - pipeline.prepared().samples)
                    flush();
                pipeline.beginPreparation();
                auto &batch = pipeline.prepared();
                batch.jobs.data()[batch.jobCount++] = job;
                batch.samples += sampleCount;
                pipeline.observeResidentBatches();
            } else
                for (std::size_t first = 0; first < rule.patches.size();) {
                    const auto limit =
                        std::min({std::size_t(32), maximumPatches, rule.patches.size() - first});
                    std::size_t count = 0, sampleCount = 0;
                    while (count < limit) {
                        const auto degreeCount =
                            rule.patches[first + count].degree == 14 ? high.size() : baseCount;
                        const auto nextSamples = galerkin::nearPatchPointCount(
                            rule.patches[first + count], static_cast<unsigned>(degreeCount));
                        if (nextSamples > settings.maximumBatchGreenSamples)
                            throw std::runtime_error(
                                "One near patch exceeds the Green sample budget.");
                        if (nextSamples > settings.maximumBatchGreenSamples - sampleCount)
                            break;
                        sampleCount += nextSamples;
                        ++count;
                    }
                    if (pipeline.prepared().jobCount == maximumJobs ||
                        pipeline.prepared().patchCount + count > maximumPatches ||
                        sampleCount >
                            settings.maximumBatchGreenSamples - pipeline.prepared().samples)
                        flush();
                    pipeline.beginPreparation();
                    auto &batch = pipeline.prepared();
                    job.integration.firstPatch = static_cast<unsigned>(batch.patchCount);
                    job.integration.patchCount = static_cast<unsigned>(count);
                    std::copy_n(rule.patches.begin() + first, count,
                                batch.patches.data() + batch.patchCount);
                    batch.patchCount += count;
                    batch.jobs.data()[batch.jobCount++] = job;
                    batch.samples += sampleCount;
                    pipeline.observeResidentBatches();
                    first += count;
                }
        };
        SupportNearJobIterator iterator(space, tree, plan, settings.maximumIteratorWorkspaceBytes);
        stats.iteratorWorkspaceBytes = iterator.workspaceBytes();
        std::vector<SupportNearJob> masks(maskBatchCount);
        stats.iteratorBatchBytes = masks.capacity() * sizeof(SupportNearJob);
        if (stats.iteratorBatchBytes + maximumJobs * sizeof(MaskedQuadratureJob) +
                maximumPatches * sizeof(galerkin::PairQuadraturePatch) >
            settings.maximumConstructionBatchBytes)
            throw std::runtime_error(
                "Actual near job/patch/mask capacities exceed the construction budget.");
        const auto preparedJob = [&](const SupportNearJob &mask,
                                     galerkin::PairPermutation permutation) {
            MaskedQuadratureJob job;
            job.integration.test = mask.testPanel;
            job.integration.trial = mask.sourcePanel;
            job.integration.permutation = permutation;
            job.pressureMask = mask.pressureMask;
            job.fluxMask = mask.fluxMask;
            job.faceColumn = plan.panelToActiveFaceColumn[mask.sourcePanel];
            return job;
        };
        const auto recordAdaptiveStats = [&] {
            const auto measured = adaptiveBatch.statistics();
            stats.adaptiveDecisionMilliseconds = measured.decisionMilliseconds;
            stats.estimatorKernelMilliseconds = measured.estimatorKernelMilliseconds;
            stats.estimatorWallMilliseconds = measured.estimatorWallMilliseconds;
            stats.adaptiveWorkspacePeakBytes = measured.peakWorkspaceBytes;
            stats.maximumSinglePairWorkspaceBytes = measured.maximumSinglePairWorkspaceBytes;
            stats.adaptiveRounds = measured.rounds;
            stats.adaptivePatchEvaluations = measured.evaluatedPatches;
            stats.maximumEstimatorBatchPatches = measured.maximumRoundPatches;
        };
        const auto consume = [&](const SupportNearJob &mask,
                                 const galerkin::PairQuadratureRule &rule) {
            ++stats.adaptiveQuadratureCalls;
            recordAdaptiveStats();
            append(preparedJob(mask, galerkin::classifyPair(space.panels[mask.testPanel],
                                                            space.panels[mask.sourcePanel])),
                   rule);
        };
        const galerkin::PairQuadratureRule emptyRule;
        while (!iterator.done()) {
            pipeline.beginPreparation();
            const auto iteratorStart =
                settings.profileConstruction ? Clock::now() : Clock::time_point{};
            const auto count = iterator.next(masks, stats.iteratorBatchBytes, stop);
            if (settings.profileConstruction)
                stats.iteratorMilliseconds += milliseconds(iteratorStart);
            for (std::size_t i = 0; i < count; ++i) {
                cancelled(stop);
                const auto mask = masks[i];
                ++stats.maskedPanelJobs;
                const auto &test = space.panels[mask.testPanel],
                           &trial = space.panels[mask.sourcePanel];
                const auto classificationStart =
                    settings.profileConstruction ? Clock::now() : Clock::time_point{};
                const auto permutation = galerkin::classifyPair(test, trial);
                unsigned regularDegree = 0;
                bool adaptive = false;
                if (permutation.topology == galerkin::PairTopology::disjoint) {
                    const auto &testBounds = bounds[mask.testPanel],
                               &trialBounds = bounds[mask.sourcePanel];
                    regularDegree = galerkin::qualifiedRegularDegree(
                        testBounds, trialBounds, k, q.relativeTolerance, q.regularOrder);
                    adaptive = regularDegree == 0;
                    if (!adaptive && q.nearDistanceRatio > 0) {
                        const Vec3 delta{testBounds.center.x - trialBounds.center.x,
                                         testBounds.center.y - trialBounds.center.y,
                                         testBounds.center.z - trialBounds.center.z};
                        const double distanceSquared =
                            delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                        const double bound =
                            testBounds.radius + trialBounds.radius +
                            q.nearDistanceRatio *
                                std::max(testBounds.maximumEdge, trialBounds.maximumEdge);
                        // A conservative sphere gap avoids exact triangle distances
                        // only when non-nearness is proved. User near settings remain binding.
                        if (!(distanceSquared > bound * bound))
                            adaptive = galerkin::isNearPair(test, trial, q.nearDistanceRatio);
                    }
                }
                if (settings.profileConstruction)
                    stats.classificationMilliseconds += milliseconds(classificationStart);
                if (adaptive) {
                    while (!adaptiveBatch.submit(mask))
                        adaptiveBatch.advance(consume);
                } else {
                    auto job = preparedJob(mask, permutation);
                    if (permutation.topology == galerkin::PairTopology::disjoint) {
                        job.regularRuleIndex = static_cast<std::uint8_t>(regularDegree == 4   ? 0
                                                                         : regularDegree == 6 ? 1
                                                                                              : 2);
                        ++stats.regularRuleJobs[job.regularRuleIndex];
                    }
                    append(job, emptyRule);
                }
            }
        }
        adaptiveBatch.drain(consume);
        recordAdaptiveStats();
        flush();
        pipeline.finish(stop);
        checkError();
        cancelled(stop);
    }
    void checkError() const {
        checked(cudaGetLastError(), "Near CUDA launch");
        int value{};
        checked(cudaMemcpyAsync(&value, error.data(), sizeof(value), cudaMemcpyDeviceToHost,
                                streamValue),
                "Near error download");
        checked(cudaStreamSynchronize(streamValue), "Near CUDA completion");
        if (value)
            throw std::runtime_error("Support near CUDA integration/scatter failed, status " +
                                     std::to_string(value) + ".");
    }
    void apply(const Matrix &matrix, const Z *values, Z *destination, std::stop_token stop) {
        cancelled(stop);
        if ((matrix.nonzeros && !values) || !destination || values == destination)
            throw std::invalid_argument("Near device vectors are null or alias.");
        checked(cudaMemsetAsync(error.data(), 0, sizeof(int), streamValue),
                "Near apply error clear");
        applyCsr<<<(matrix.rows + 3) / 4, 128, 0, streamValue>>>(matrix.view(), values, destination,
                                                                 error.data());
        checkError();
        cancelled(stop);
    }
    std::vector<Complex> hostApply(bool flux, std::span<const Complex> values,
                                   std::stop_token stop) {
        cancelled(stop);
        const auto &matrix = flux ? b : a;
        if (values.size() != matrix.columnCount)
            throw std::invalid_argument("Near input vector dimension is invalid.");
        const auto host = packed(values);
        auto &input = flux ? inputFace : inputPressure;
        input.upload(host, streamValue);
        apply(matrix, input.data(), output.data(), stop);
        const auto result = output.download(streamValue);
        std::vector<Complex> converted;
        converted.reserve(result.size());
        for (auto z : result)
            converted.emplace_back(z.x, z.y);
        return converted;
    }
    std::unique_ptr<Stream> ownedStream;
    cudaStream_t streamValue{};
    Matrix a, b;
    DeviceBuffer<Z> inputPressure, inputFace, output;
    DeviceBuffer<int> error;
    std::size_t liveBytes{}, peakBytes{};
    SupportNearConstructionStats stats;
};
SupportNearOperator::SupportNearOperator(
    const galerkin::Space &space, const SupportTree &tree, const SupportNearPlan &plan, double k,
    SupportNearOperatorSettings settings, std::stop_token stop, cudaStream_t stream,
    std::function<void(const SupportNearConstructionStats &)> progress)
    : impl_(std::make_unique<Implementation>(space, tree, plan, k, settings, stop, stream,
                                             progress)) {}
SupportNearOperator::~SupportNearOperator() = default;
std::vector<Complex> SupportNearOperator::applyPressure(std::span<const Complex> values,
                                                        std::stop_token stop) {
    return impl_->hostApply(false, values, stop);
}
std::vector<Complex> SupportNearOperator::applyFaceFlux(std::span<const Complex> values,
                                                        std::stop_token stop) {
    return impl_->hostApply(true, values, stop);
}
void SupportNearOperator::applyPressureDevice(const cuDoubleComplex *values,
                                              cuDoubleComplex *output, std::stop_token stop) {
    impl_->apply(impl_->a, values, output, stop);
}
void SupportNearOperator::applyFaceFluxDevice(const cuDoubleComplex *values,
                                              cuDoubleComplex *output, std::stop_token stop) {
    impl_->apply(impl_->b, values, output, stop);
}
SupportDeviceCsr SupportNearOperator::pressureMatrix() const noexcept {
    return impl_->a.view();
}
SupportDeviceCsr SupportNearOperator::faceFluxMatrix() const noexcept {
    return impl_->b.view();
}
cudaStream_t SupportNearOperator::stream() const noexcept {
    return impl_->streamValue;
}
std::size_t SupportNearOperator::deviceBytes() const noexcept {
    return impl_->liveBytes;
}
std::size_t SupportNearOperator::peakDeviceBytes() const noexcept {
    return impl_->peakBytes;
}
SupportNearConstructionStats SupportNearOperator::constructionStats() const noexcept {
    return impl_->stats;
}
} // namespace hornsculpt::wideband
