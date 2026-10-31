#include "adaptive_patch_estimator.hpp"
#include "wideband_runtime.cuh"
#include "../galerkin/near_projection.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <type_traits>

namespace hornsculpt::wideband {
namespace {
using namespace detail;
struct Geometry {
    Vec3 vertices[3], normal;
    double edge, area;
};
struct Evaluation {
    Z integral[36];
    double error, innerError, outerError;
};
template <class T> class PinnedMirror {
  public:
    static_assert(std::is_nothrow_default_constructible_v<T> &&
                  std::is_trivially_destructible_v<T>);
    PinnedMirror() = default;
    ~PinnedMirror() {
        if (values_)
            cudaFreeHost(values_);
    }
    PinnedMirror(const PinnedMirror &) = delete;
    PinnedMirror &operator=(const PinnedMirror &) = delete;
    void allocate(std::size_t count) {
        if (values_ || !count)
            throw std::logic_error("Estimator pinned mirror allocation is invalid.");
        void *storage{};
        checked(cudaHostAlloc(&storage, product(count, sizeof(T)), cudaHostAllocDefault),
                "Estimator pinned mirror allocation");
        values_ = static_cast<T *>(storage);
        count_ = count;
        std::uninitialized_value_construct_n(values_, count_);
    }
    T *data() const noexcept {
        return values_;
    }
    std::size_t bytes() const noexcept {
        return count_ * sizeof(T);
    }

  private:
    T *values_{};
    std::size_t count_{};
};
std::size_t estimatorHostBytes(std::size_t count) {
    return sum(product(count, sizeof(AdaptivePatchTask) + sizeof(Evaluation) +
                                  sizeof(AdaptivePatchMoments)),
               sizeof(int));
}
std::size_t greenSamples(const galerkin::PairQuadraturePatch &patch) {
    const unsigned low = patch.degree == 8 ? 12 : 49, high = patch.degree == 8 ? 16 : 64;
    if (!patch.projection)
        return low * low + high * high;
    return sum(sum(galerkin::nearPatchPointCount(patch, low),
                   galerkin::nearPatchPointCount(patch, high)),
               galerkin::nearPatchPointCount(patch, high, 14));
}
__device__ Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
__device__ Vec3 local(const Geometry &g, const double *b) {
    const auto x = subtract(g.vertices[1], g.vertices[0]),
               y = subtract(g.vertices[2], g.vertices[0]);
    return {x.x * b[1] + y.x * b[2], x.y * b[1] + y.y * b[2], x.z * b[1] + y.z * b[2]};
}
__device__ Z warpSum(Z v) {
    for (unsigned d = 16; d; d /= 2) {
        v.x += __shfl_down_sync(0xffffffffu, v.x, d);
        v.y += __shfl_down_sync(0xffffffffu, v.y, d);
    }
    return v;
}
__global__ void estimate(const Geometry *panels, const AdaptivePatchTask *tasks,
                         const galerkin::TriangleQuadraturePoint *rules, double k, Z *moments,
                         int *error) {
    const auto task = tasks[blockIdx.x / 3];
    const unsigned role = blockIdx.x % 3;
    if (role == 2 && !task.patch.projection)
        return;
    const bool fine = role != 0;
    const unsigned innerOrder = role == 2 ? 14 : 20;
    const unsigned offset = task.patch.degree == 8 ? (fine ? 12 : 0) : (fine ? 77 : 28);
    const unsigned count = task.patch.degree == 8 ? (fine ? 16 : 12) : (fine ? 64 : 49);
    const auto &test = panels[task.testPanel], &trial = panels[task.sourcePanel];
    const auto origin = subtract(test.vertices[0], trial.vertices[0]);
    const double scale = fmax(test.edge, trial.edge);
    Z values[36]{};
    const unsigned points = galerkin::nearPatchPointCount(task.patch, count, innerOrder);
    for (unsigned index = threadIdx.x; index < points; index += blockDim.x) {
        const auto point = galerkin::mapNearPatchPoint(rules + offset, count, task.patch,
                                                       index, test, trial, innerOrder);
        const auto tx = local(test, point.test), ty = local(trial, point.trial);
        const Vec3 d{origin.x + tx.x - ty.x, origin.y + tx.y - ty.y, origin.z + tx.z - ty.z};
        const double r2 = d.x * d.x + d.y * d.y + d.z * d.z, r = sqrt(r2);
        if (!(r > 0) || !isfinite(r)) {
            atomicOr(error, 1);
            continue;
        }
        double sine, cosine;
        sincos(k * r, &sine, &cosine);
        const Z g = make_cuDoubleComplex(cosine / r, -sine / r);
        const double factor = -scale / r2;
        const Z derivative =
            make_cuDoubleComplex((g.x - g.y * k * r) * factor, (g.y + g.x * k * r) * factor);
        const Z kernels[4]{g, make_cuDoubleComplex(derivative.x * d.x, derivative.y * d.x),
                           make_cuDoubleComplex(derivative.x * d.y, derivative.y * d.y),
                           make_cuDoubleComplex(derivative.x * d.z, derivative.y * d.z)};
#pragma unroll
        for (unsigned i = 0; i < 3; ++i)
#pragma unroll
            for (unsigned j = 0; j < 3; ++j) {
                const double weight = point.weight * point.test[i] * point.trial[j];
#pragma unroll
                for (unsigned c = 0; c < 4; ++c) {
                    auto &v = values[c * 9 + 3 * i + j];
                    v.x += kernels[c].x * weight;
                    v.y += kernels[c].y * weight;
                }
            }
    }
    __shared__ Z partial[4];
#pragma unroll
    for (unsigned c = 0; c < 36; ++c) {
        const auto sum = warpSum(values[c]);
        if (!(threadIdx.x & 31))
            partial[threadIdx.x / 32] = sum;
        __syncthreads();
        if (threadIdx.x < 32) {
            const auto total =
                warpSum(threadIdx.x < 4 ? partial[threadIdx.x] : make_cuDoubleComplex(0, 0));
            if (!threadIdx.x)
                moments[std::size_t(blockIdx.x) * 36 + c] = total;
        }
        __syncthreads();
    }
}
__global__ void finish(const Z *moments, const AdaptivePatchTask *tasks, Evaluation *output, unsigned count, int *error) {
    const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= count)
        return;
    double maximum = 0, innerMaximum = 0;
    for (unsigned j = 0; j < 36; ++j) {
        const auto coarse = moments[(std::size_t(i) * 3) * 36 + j],
                   fine = moments[(std::size_t(i) * 3 + 1) * 36 + j];
        if (!isfinite(fine.x) || !isfinite(fine.y) || !isfinite(coarse.x) || !isfinite(coarse.y))
            atomicOr(error, 2);
        output[i].integral[j] = fine;
        maximum = fmax(maximum, hypot(fine.x - coarse.x, fine.y - coarse.y));
        if (tasks[i].patch.projection) {
            const auto inner = moments[(std::size_t(i) * 3 + 2) * 36 + j];
            if (!isfinite(inner.x) || !isfinite(inner.y))
                atomicOr(error, 2);
            innerMaximum = fmax(innerMaximum, hypot(fine.x - inner.x, fine.y - inner.y));
        }
    }
    if (!isfinite(maximum))
        atomicOr(error, 2);
    output[i].outerError = maximum;
    output[i].innerError = innerMaximum;
    output[i].error = maximum + innerMaximum;
}
class Timer {
  public:
    Timer() {
        checked(cudaEventCreate(&start), "Estimator start event");
        try {
            checked(cudaEventCreate(&end), "Estimator end event");
        } catch (...) {
            cudaEventDestroy(start);
            throw;
        }
    }
    ~Timer() {
        cudaEventDestroy(start);
        cudaEventDestroy(end);
    }
    cudaEvent_t start{}, end{};
};
} // namespace
struct AdaptivePatchEstimator::Implementation {
    Stream stream;
    const galerkin::Space &space;
    std::size_t panelCount{};
    DeviceBuffer<Geometry> geometry;
    DeviceBuffer<AdaptivePatchTask> tasks;
    DeviceBuffer<galerkin::TriangleQuadraturePoint> rules;
    DeviceBuffer<Z> moments;
    DeviceBuffer<Evaluation> output;
    DeviceBuffer<int> error;
    PinnedMirror<AdaptivePatchTask> hostTasks;
    PinnedMirror<Evaluation> hostEvaluated;
    PinnedMirror<int> hostStatus;
    std::vector<AdaptivePatchMoments> hostResult;
    std::size_t capacity, bytes{}, chunks{}, largestChunkSamples{};
    double kernelMs{};
    Timer timer;
    Implementation(const galerkin::Space &geometrySpace, double k_, std::size_t count,
                   std::size_t budget, cudaStream_t borrowed)
        : stream(borrowed), space(geometrySpace), capacity(count), k(k_) {
        if (space.panels.empty() || space.panels.size() > UINT32_MAX || !count ||
            count > UINT32_MAX / 3 || !(k >= 0) || !std::isfinite(k))
            throw std::invalid_argument(
                "Adaptive estimator geometry, frequency or batch capacity is invalid.");
        std::vector<galerkin::TriangleQuadraturePoint> hostRules;
        for (unsigned degree : {6u, 8u, 12u, 14u}) {
            const auto rule = galerkin::triangleRule(degree);
            hostRules.insert(hostRules.end(), rule.begin(), rule.end());
        }
        bytes =
            sum(product(space.panels.size(), sizeof(Geometry)),
                sum(product(count, sizeof(AdaptivePatchTask) + 108 * sizeof(Z) + sizeof(Evaluation)),
                    sum(product(hostRules.size(), sizeof(galerkin::TriangleQuadraturePoint)),
                        sizeof(int))));
        std::size_t available{}, total{};
        checked(cudaMemGetInfo(&available, &total), "Estimator memory preflight");
        if (bytes > budget || bytes > available - available / 5)
            throw std::runtime_error(
                "Adaptive estimator device budget exceeded before allocation.");
        panelCount = space.panels.size();
        geometry.allocate(panelCount);
        std::array<Geometry, 256> chunk;
        for (std::size_t first = 0; first < panelCount; first += chunk.size()) {
            const auto n = std::min(chunk.size(), panelCount - first);
            for (std::size_t i = 0; i < n; ++i) {
                for (unsigned j = 0; j < 3; ++j)
                    chunk[i].vertices[j] = space.panels[first + i].vertices[j];
                chunk[i].edge = space.panels[first + i].maximumEdgeLength;
                chunk[i].normal = space.panels[first + i].normal;
                chunk[i].area = space.panels[first + i].area;
            }
            checked(cudaMemcpyAsync(geometry.data() + first, chunk.data(), n * sizeof(Geometry),
                                    cudaMemcpyHostToDevice, stream),
                    "Estimator geometry upload");
            checked(cudaStreamSynchronize(stream), "Estimator geometry chunk lifetime");
        }
        rules.allocate(hostRules.size());
        rules.upload(hostRules, stream);
        tasks.allocate(capacity);
        moments.allocate(product(capacity, 108));
        output.allocate(capacity);
        error.allocate(1);
        hostTasks.allocate(capacity);
        hostEvaluated.allocate(capacity);
        hostStatus.allocate(1);
        hostResult.resize(capacity);
        if (hostResult.capacity() != capacity)
            throw std::runtime_error(
                "Estimator host capacity exceeds its fixed payload reservation.");
        checked(cudaStreamSynchronize(stream), "Estimator geometry upload");
    }
    double k;
    ~Implementation() {
        // The same stream can be borrowed by near assembly. Drain before any
        // pinned or device mirror is released during error/cancellation unwind.
        cudaStreamSynchronize(stream);
    }
    std::span<const AdaptivePatchMoments> evaluateView(std::span<const AdaptivePatchTask> values,
                                                       std::stop_token stop) {
        try {
            return evaluateImpl(values, stop);
        } catch (...) {
            checked(cudaStreamSynchronize(stream), "Estimator failure drain");
            throw;
        }
    }
    std::span<const AdaptivePatchMoments> evaluateImpl(std::span<const AdaptivePatchTask> values,
                                                       std::stop_token stop) {
        cancelled(stop);
        kernelMs = 0;
        chunks = largestChunkSamples = 0;
        if (values.size() > capacity)
            throw std::invalid_argument("Adaptive estimator batch exceeds its fixed capacity.");
        if (values.empty())
            return {};
        for (const auto &v : values) {
            if (v.testPanel >= panelCount || v.sourcePanel >= panelCount ||
                (v.patch.degree != 8 && v.patch.degree != 14) || !(v.patch.areaFraction > 0) ||
                !std::isfinite(v.patch.areaFraction) || v.patch.projection > 1 ||
                (v.patch.projection && (!v.patch.projectionIntervals ||
                                        v.patch.projectionIntervals > 128)))
                throw std::invalid_argument("Adaptive estimator pair/patch contract is invalid.");
            if (galerkin::classifyPair(space.panels[v.testPanel], space.panels[v.sourcePanel])
                    .topology != galerkin::PairTopology::disjoint)
                throw std::invalid_argument("Adaptive estimator requires a disjoint panel pair.");
            for (const auto *bary : {v.patch.test, v.patch.trial})
                for (unsigned i = 0; i < 3; ++i) {
                    double sum = 0;
                    for (unsigned j = 0; j < 3; ++j) {
                        const auto x = bary[3 * i + j];
                        if (!std::isfinite(x) || x < 0 || x > 1)
                            throw std::invalid_argument(
                                "Invalid adaptive patch barycentric coordinate.");
                        sum += x;
                    }
                    if (std::abs(sum - 1) > 1e-12)
                        throw std::invalid_argument(
                            "Adaptive patch barycentric row is not normalized.");
                }
        }
        std::copy(values.begin(), values.end(), hostTasks.data());
        for (std::size_t first = 0; first < values.size();) {
            cancelled(stop);
            std::size_t count = 0, samples = 0;
            while (first + count < values.size()) {
                const auto next = greenSamples(values[first + count].patch);
                if (next > maximumChunkGreenSamples)
                    throw std::runtime_error("One adaptive patch exceeds the estimator Green-sample chunk budget.");
                if (next > maximumChunkGreenSamples - samples)
                    break;
                samples += next;
                ++count;
            }
            checked(cudaMemcpyAsync(tasks.data(), hostTasks.data() + first,
                                    product(count, sizeof(AdaptivePatchTask)),
                                    cudaMemcpyHostToDevice, stream),
                    "Estimator pinned task upload");
            checked(cudaMemsetAsync(error.data(), 0, sizeof(int), stream), "Estimator error clear");
            checked(cudaEventRecord(timer.start, stream), "Estimator timer start");
            estimate<<<static_cast<unsigned>(3 * count), 128, 0, stream>>>(
                geometry.data(), tasks.data(), rules.data(), k, moments.data(), error.data());
            finish<<<static_cast<unsigned>((count + 127) / 128), 128, 0, stream>>>(
                moments.data(), tasks.data(), output.data(), static_cast<unsigned>(count), error.data());
            checked(cudaGetLastError(), "Adaptive estimator kernels");
            checked(cudaEventRecord(timer.end, stream), "Estimator timer end");
            checked(cudaMemcpyAsync(hostEvaluated.data() + first, output.data(),
                                    product(count, sizeof(Evaluation)), cudaMemcpyDeviceToHost,
                                    stream),
                    "Estimator pinned moments download");
            checked(cudaMemcpyAsync(hostStatus.data(), error.data(), sizeof(int),
                                    cudaMemcpyDeviceToHost, stream),
                    "Estimator pinned status download");
            checked(cudaStreamSynchronize(stream), "Estimator chunk completion");
            cancelled(stop);
            if (*hostStatus.data())
                throw std::runtime_error("Adaptive estimator encountered an invalid Green sample.");
            float elapsed{};
            checked(cudaEventElapsedTime(&elapsed, timer.start, timer.end),
                    "Estimator kernel time");
            kernelMs += elapsed;
            ++chunks;
            largestChunkSamples = std::max(largestChunkSamples, samples);
            first += count;
        }
        auto &result = hostResult;
        for (std::size_t i = 0; i < values.size(); ++i) {
            for (unsigned j = 0; j < 36; ++j)
                result[i].integral[j] = {hostEvaluated.data()[i].integral[j].x,
                                         hostEvaluated.data()[i].integral[j].y};
            result[i].absoluteError = hostEvaluated.data()[i].error;
            result[i].innerError = hostEvaluated.data()[i].innerError;
            result[i].outerError = hostEvaluated.data()[i].outerError;
        }
        return std::span(result).first(values.size());
    }
};
AdaptivePatchEstimator::AdaptivePatchEstimator(const galerkin::Space &space, double k,
                                               std::size_t count, std::size_t budget,
                                               cudaStream_t stream)
    : impl_(std::make_unique<Implementation>(space, k, count, budget, stream)) {}
AdaptivePatchEstimator::~AdaptivePatchEstimator() = default;
std::vector<AdaptivePatchMoments>
AdaptivePatchEstimator::evaluate(std::span<const AdaptivePatchTask> values, std::stop_token stop) {
    const auto view = impl_->evaluateView(values, stop);
    return {view.begin(), view.end()};
}
std::span<const AdaptivePatchMoments>
AdaptivePatchEstimator::evaluateView(std::span<const AdaptivePatchTask> tasks,
                                     std::stop_token stop) {
    return impl_->evaluateView(tasks, stop);
}
std::size_t AdaptivePatchEstimator::deviceBytes() const noexcept {
    return impl_->bytes;
}
std::size_t AdaptivePatchEstimator::hostBytes() const noexcept {
    return impl_->hostTasks.bytes() + impl_->hostEvaluated.bytes() + impl_->hostStatus.bytes() +
           impl_->hostResult.capacity() * sizeof(AdaptivePatchMoments);
}
std::size_t AdaptivePatchEstimator::lastChunkCount() const noexcept {
    return impl_->chunks;
}
std::size_t AdaptivePatchEstimator::lastMaximumChunkGreenSamples() const noexcept {
    return impl_->largestChunkSamples;
}
double AdaptivePatchEstimator::lastKernelMilliseconds() const noexcept {
    return impl_->kernelMs;
}

struct AdaptivePairBatch::Implementation {
    using Clock = std::chrono::steady_clock;
    struct Slot {
        SupportNearJob job;
        std::unique_ptr<galerkin::AdaptivePairQuadratureState> state;
        std::size_t first{}, count{};
    };
    const galerkin::Space &space;
    galerkin::PairQuadratureSettings settings;
    double k;
    std::stop_token stop;
    galerkin::AdaptiveQuadratureWorkspace workspace;
    std::vector<Slot> slots;
    std::vector<AdaptivePatchTask> tasks;
    std::unique_ptr<AdaptivePatchEstimator> estimator;
    AdaptivePairBatchStats stats;
    std::size_t active{}, next{};
    static double ms(Clock::time_point start) {
        return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    }
    Implementation(const galerkin::Space &geometry, double frequency,
                   galerkin::PairQuadratureSettings q, std::size_t hostBytes,
                   std::size_t deviceBudget, std::size_t capacity, cudaStream_t stream,
                   std::stop_token token)
        : space(geometry), settings(q), k(frequency), stop(token), workspace(hostBytes) {
        if (!capacity || capacity > 4096)
            throw std::invalid_argument("Adaptive pair batch capacity must be 1..4096.");
        const auto taskCount = product(capacity, 4);
        // Reserve the CPU task list and all pinned task/evaluation/status
        // mirrors plus the converted result vector before any allocation.
        // Allocator metadata and immutable rule tables remain outside payload.
        workspace.acquire(
            sum(sizeof(Implementation), sum(product(capacity, sizeof(Slot)),
                                            sum(product(taskCount, sizeof(AdaptivePatchTask)),
                                                estimatorHostBytes(taskCount)))));
        slots.resize(capacity);
        tasks.reserve(taskCount);
        if (slots.capacity() != capacity || tasks.capacity() != taskCount)
            throw std::runtime_error(
                "Unexpected adaptive batch vector capacity exceeds its reserved payload.");
        estimator =
            std::make_unique<AdaptivePatchEstimator>(space, k, taskCount, deviceBudget, stream);
        if (estimator->hostBytes() != estimatorHostBytes(taskCount))
            throw std::runtime_error("Estimator host mirrors exceed their payload reservation.");
    }
    ~Implementation() {
        for (auto &slot : slots)
            if (slot.state) {
                slot.state.reset();
                workspace.release(sizeof(galerkin::AdaptivePairQuadratureState));
            }
    }
    bool submit(const SupportNearJob &job) {
        cancelled(stop);
        if (active == slots.size())
            return false;
        if (job.testPanel >= space.panels.size() || job.sourcePanel >= space.panels.size())
            throw std::invalid_argument("Adaptive pair batch panel index is invalid.");
        for (std::size_t j = 0; j < slots.size(); ++j) {
            auto &slot = slots[(next + j) % slots.size()];
            if (slot.state)
                continue;
            const auto start = Clock::now();
            workspace.acquire(sizeof(galerkin::AdaptivePairQuadratureState));
            try {
                slot.state = std::make_unique<galerkin::AdaptivePairQuadratureState>(
                    space.panels[job.testPanel], space.panels[job.sourcePanel], settings, k,
                    workspace, stop);
            } catch (...) {
                workspace.release(sizeof(galerkin::AdaptivePairQuadratureState));
                throw;
            }
            stats.decisionMilliseconds += ms(start);
            slot.job = job;
            ++active;
            next = (next + j + 1) % slots.size();
            return true;
        }
        throw std::logic_error("Adaptive batch free-slot count is inconsistent.");
    }
    void advance(const Consumer &consume) {
        cancelled(stop);
        if (!active)
            return;
        auto start = Clock::now();
        tasks.clear();
        for (auto &slot : slots)
            if (slot.state) {
                slot.first = tasks.size();
                const auto pending = slot.state->requests();
                slot.count = pending.size();
                for (const auto &patch : pending)
                    tasks.push_back({slot.job.testPanel, slot.job.sourcePanel, patch});
            }
        stats.decisionMilliseconds += ms(start);
        start = Clock::now();
        const auto evaluated = estimator->evaluateView(tasks, stop);
        stats.estimatorWallMilliseconds += ms(start);
        stats.estimatorKernelMilliseconds += estimator->lastKernelMilliseconds();
        ++stats.rounds;
        stats.evaluatedPatches += tasks.size();
        stats.maximumRoundPatches = std::max(stats.maximumRoundPatches, tasks.size());
        for (auto &slot : slots)
            if (slot.state) {
                start = Clock::now();
                slot.state->commit(std::span(evaluated).subspan(slot.first, slot.count));
                stats.decisionMilliseconds += ms(start);
                if (slot.state->done()) {
                    stats.maximumSinglePairWorkspaceBytes =
                        std::max(stats.maximumSinglePairWorkspaceBytes,
                                 slot.state->result().peakWorkspaceBytes);
                    consume(slot.job, slot.state->result());
                    slot.state.reset();
                    workspace.release(sizeof(galerkin::AdaptivePairQuadratureState));
                    --active;
                    ++stats.completedPairs;
                }
            }
        stats.peakWorkspaceBytes = workspace.peakBytes();
    }
};
AdaptivePairBatch::AdaptivePairBatch(const galerkin::Space &space, double k,
                                     galerkin::PairQuadratureSettings settings,
                                     std::size_t hostBytes, std::size_t deviceBytes,
                                     std::size_t active, cudaStream_t stream, std::stop_token stop)
    : impl_(std::make_unique<Implementation>(space, k, settings, hostBytes, deviceBytes, active,
                                             stream, stop)) {}
AdaptivePairBatch::~AdaptivePairBatch() = default;
bool AdaptivePairBatch::submit(const SupportNearJob &job) {
    return impl_->submit(job);
}
void AdaptivePairBatch::advance(const Consumer &consume) {
    impl_->advance(consume);
}
void AdaptivePairBatch::drain(const Consumer &consume) {
    while (impl_->active)
        impl_->advance(consume);
}
std::size_t AdaptivePairBatch::deviceBytes() const noexcept {
    return impl_->estimator->deviceBytes();
}
AdaptivePairBatchStats AdaptivePairBatch::statistics() const noexcept {
    auto result = impl_->stats;
    result.liveWorkspaceBytes = impl_->workspace.liveBytes();
    result.peakWorkspaceBytes = impl_->workspace.peakBytes();
    return result;
}
} // namespace hornsculpt::wideband
