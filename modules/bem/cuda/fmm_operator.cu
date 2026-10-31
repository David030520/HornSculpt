#include "../fmm/tree.hpp"
#include "fmm_operator.hpp"
#include "panel_integrals.cuh"
#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cublas_v2.h>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace hornsculpt::cuda_detail {
namespace {
using Z = cuDoubleComplex;
constexpr unsigned channels = 7;
constexpr unsigned threads = 128;
constexpr unsigned preconditionerWidth = 32;

void checked(cudaError_t status, const char *operation) {
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
}
void checked(cublasStatus_t status, const char *operation) {
    if (status != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error(std::string(operation) + ": cuBLAS " + std::to_string(status));
}
class BlasHandle {
  public:
    BlasHandle() = default;
    ~BlasHandle() {
        if (handle_)
            cublasDestroy(handle_);
    }
    BlasHandle(const BlasHandle &) = delete;
    BlasHandle &operator=(const BlasHandle &) = delete;
    void create(cudaStream_t stream) {
        checked(cublasCreate(&handle_), "FMM block preconditioner cuBLAS creation");
        checked(cublasSetStream(handle_, stream), "FMM block preconditioner cuBLAS stream");
    }
    cublasHandle_t get() const noexcept {
        return handle_;
    }

  private:
    cublasHandle_t handle_{};
};
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Fast Multipole calculation cancelled.");
}
std::size_t productSize(std::size_t a, std::size_t b) {
    if (b && a > std::numeric_limits<std::size_t>::max() / b)
        throw std::overflow_error("Fast Multipole storage size overflow.");
    return a * b;
}
template <class T> class Buffer {
  public:
    ~Buffer() {
        if (pointer_)
            cudaFree(pointer_);
    }
    Buffer() = default;
    Buffer(const Buffer &) = delete;
    Buffer &operator=(const Buffer &) = delete;
    void allocate(std::size_t count) {
        if (pointer_)
            throw std::logic_error("Fast Multipole buffer already allocated.");
        bytes_ = productSize(count, sizeof(T));
        if (bytes_)
            checked(cudaMalloc(reinterpret_cast<void **>(&pointer_), bytes_), "FMM allocation");
    }
    void assign(const std::vector<T> &values, cudaStream_t stream) {
        allocate(values.size());
        if (bytes_)
            checked(
                cudaMemcpyAsync(pointer_, values.data(), bytes_, cudaMemcpyHostToDevice, stream),
                "FMM upload");
    }
    void clear(cudaStream_t stream) {
        if (bytes_)
            checked(cudaMemsetAsync(pointer_, 0, bytes_, stream), "FMM clear");
    }
    void release() {
        if (pointer_)
            checked(cudaFree(pointer_), "FMM buffer release");
        pointer_ = nullptr;
        bytes_ = 0;
    }
    T *data() const noexcept {
        return pointer_;
    }
    std::size_t bytes() const noexcept {
        return bytes_;
    }
    std::vector<T> download(cudaStream_t stream) const {
        std::vector<T> result(bytes_ / sizeof(T));
        if (bytes_)
            checked(
                cudaMemcpyAsync(result.data(), pointer_, bytes_, cudaMemcpyDeviceToHost, stream),
                "FMM download");
        checked(cudaStreamSynchronize(stream), "FMM download synchronization");
        return result;
    }

  private:
    T *pointer_{};
    std::size_t bytes_{};
};

struct Node {
    Vec3 center;
    double halfWidth;
    unsigned parent, child[2], firstPanel, panelCount;
};
struct DirectJob {
    unsigned test{}, trial{}, firstPatch{}, patchCount{}, swapBarycentric{};
    galerkin::PairPermutation permutation;
};
struct AuditJob {
    unsigned test{}, trialNode{};
};
struct SparseMatrix {
    const std::size_t *offsets;
    const unsigned *columns;
    Z *values;
};

__global__ void extractPreconditionerBlocks(SparseMatrix matrix, const unsigned *permutation,
                                            unsigned pressureCount, std::size_t entries, Z *blocks,
                                            int *error) {
    const auto index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= entries)
        return;
    const auto block = index / (preconditionerWidth * preconditionerWidth);
    const unsigned row = index % preconditionerWidth,
                   column = (index / preconditionerWidth) % preconditionerWidth;
    const auto permutedRow = block * preconditionerWidth + row,
               permutedColumn = block * preconditionerWidth + column;
    Z value{};
    // Identity padding decouples the incomplete final block from real DOFs.
    if (permutedRow >= pressureCount || permutedColumn >= pressureCount)
        value.x = row == column ? 1.0 : 0.0;
    else {
        const unsigned globalRow = permutation[permutedRow],
                       globalColumn = permutation[permutedColumn];
        auto first = matrix.offsets[globalRow], last = matrix.offsets[globalRow + 1];
        const auto end = last;
        while (first < last) {
            const auto middle = first + (last - first) / 2;
            if (matrix.columns[middle] < globalColumn)
                first = middle + 1;
            else
                last = middle;
        }
        if (first < end && matrix.columns[first] == globalColumn)
            value = matrix.values[first];
    }
    if (!isfinite(value.x) || !isfinite(value.y))
        atomicExch(error, 2);
    blocks[index] = value;
}

__global__ void gatherPreconditionerRhs(const Z *values, const unsigned *permutation,
                                        unsigned pressureCount, std::size_t paddedCount, Z *rhs) {
    const auto index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index < paddedCount)
        rhs[index] = index < pressureCount ? values[permutation[index]] : Z{};
}

__global__ void scatterPreconditionerSolution(const Z *rhs, const unsigned *permutation,
                                              unsigned pressureCount, Z *values, int *error) {
    const auto index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index < pressureCount) {
        const auto value = rhs[index];
        if (!isfinite(value.x) || !isfinite(value.y))
            atomicExch(error, 2);
        values[permutation[index]] = value;
    }
}

__global__ void validatePreconditionerFactors(const Z *factors, std::size_t count, int *error) {
    const auto index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index < count && (!isfinite(factors[index].x) || !isfinite(factors[index].y)))
        atomicExch(error, 2);
}

__device__ Z zadd(Z a, Z b) {
    return make_cuDoubleComplex(a.x + b.x, a.y + b.y);
}
__device__ Z zscale(Z a, double b) {
    return make_cuDoubleComplex(a.x * b, a.y * b);
}
__device__ Z zproduct(Z a, Z b) {
    return make_cuDoubleComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}
__device__ Z zimaginary(Z a, double factor) {
    return make_cuDoubleComplex(-factor * a.y, factor * a.x);
}
__device__ void atomicZ(Z *output, Z value) {
    atomicAdd(&output->x, value.x);
    atomicAdd(&output->y, value.y);
}
__device__ Z converted(DoubleComplex z) {
    return make_cuDoubleComplex(z.real, z.imag);
}
template <bool mixed> struct MultipoleSum;
template <> struct MultipoleSum<false> {
    Z sum{};
    __device__ void addProduct(Z kernel, Z moment) {
        sum = zadd(sum, zproduct(kernel, moment));
    }
    __device__ Z value() const {
        return sum;
    }
};
template <> struct MultipoleSum<true> {
    float real{}, imag{}, realCorrection{}, imagCorrection{};
    __device__ void addProduct(Z kernel, Z moment) {
        const float kr = static_cast<float>(kernel.x), ki = static_cast<float>(kernel.y),
                    mr = static_cast<float>(moment.x), mi = static_cast<float>(moment.y);
        const float productReal = kr * mr - ki * mi, productImag = kr * mi + ki * mr;
        // Compensated FP32 amplitudes preserve small contributions across large
        // interaction lists. Geometry, propagation phase and translations remain
        // FP64; strict residual requests select the FP64 accumulator above.
        const float nextReal = productReal - realCorrection,
                    nextImag = productImag - imagCorrection;
        const float totalReal = real + nextReal, totalImag = imag + nextImag;
        realCorrection = (totalReal - real) - nextReal;
        imagCorrection = (totalImag - imag) - nextImag;
        real = totalReal;
        imag = totalImag;
    }
    __device__ Z value() const {
        return make_cuDoubleComplex(static_cast<double>(real) - realCorrection,
                                    static_cast<double>(imag) - imagCorrection);
    }
};
__device__ void sparseAdd(SparseMatrix matrix, unsigned row, unsigned column, Z value, int *error) {
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
        atomicExch(error, 4);
        return;
    }
    atomicZ(matrix.values + first, value);
}
__device__ double coordinate(Vec3 v, unsigned axis) {
    return axis == 0 ? v.x : axis == 1 ? v.y : v.z;
}
__device__ Vec3 interpolationPoint(const Node &node, const double *abscissae, unsigned index,
                                   unsigned order) {
    const unsigned x = index % order, y = (index / order) % order, z = index / (order * order);
    return {node.center.x + node.halfWidth * abscissae[x],
            node.center.y + node.halfWidth * abscissae[y],
            node.center.z + node.halfWidth * abscissae[z]};
}
__device__ PairMoments promotedMoments(const MixedMoments &mixed) {
    PairMoments result{};
    result.green = promoted(mixed.green);
    for (unsigned i = 0; i < 9; ++i)
        result.matrix[i] = promoted(mixed.matrix[i]);
    for (unsigned i = 0; i < 3; ++i)
        result.rhs[i] = promoted(mixed.rhs[i]);
    return result;
}
__device__ void regularMoments(PairMoments &moments, const Panel &test, const Panel &trial,
                               const galerkin::TriangleQuadraturePoint *rule, unsigned count,
                               double k, int *error) {
    // Match the dense solver's regular FP32 integral kernels. Origin subtraction
    // remains FP64; contractions with the current pressure and global sums stay
    // FP64. The direct audit still evaluates every far pair independently.
    MixedMoments mixed{};
    for (unsigned x = 0; x < count; ++x)
        for (unsigned y = 0; y < count; ++y) {
            galerkin::PairQuadraturePoint sample;
            for (unsigned i = 0; i < 3; ++i) {
                sample.test[i] = rule[x].barycentric[i];
                sample.trial[i] = rule[y].barycentric[i];
            }
            sample.weight = rule[x].weight * rule[y].weight;
            mixedSample(mixed, pairDelta(test, trial, sample.test, sample.trial), test, trial,
                        sample.test, sample.trial,
                        static_cast<float>(sample.weight * 4 * test.area * trial.area),
                        static_cast<float>(k), 1.f / static_cast<float>(k), false, error);
        }
    moments = promotedMoments(mixed);
}
__device__ Z pressureCoefficient(const PairMoments &moments, const Panel &test, const Panel &trial,
                                 unsigned i, unsigned j, double k) {
    const double curl = dot(test.curls[i], trial.curls[j]) / k;
    return make_cuDoubleComplex(moments.matrix[3 * i + j].real + moments.green.imag * curl,
                                moments.matrix[3 * i + j].imag - moments.green.real * curl);
}

__global__ void buildDirect(const Panel *panels, const DirectJob *jobs, std::size_t first,
                            const galerkin::TriangleQuadraturePoint *regular, unsigned regularCount,
                            const galerkin::GaussPoint *gauss, unsigned singularOrder,
                            const galerkin::PairQuadraturePatch *patches,
                            const galerkin::TriangleQuadraturePoint *nearRule,
                            unsigned baseNearCount, unsigned highNearCount, double k,
                            SparseMatrix pressureMatrix, SparseMatrix fluxMatrix, int *error) {
    const auto index = first + blockIdx.x;
    const auto job = jobs[index];
    const auto &test = panels[job.test];
    const auto &trial = panels[job.trial];
    PairMoments moments{};
    if (job.permutation.topology != galerkin::PairTopology::disjoint) {
        const unsigned perSector = singularOrder * singularOrder * singularOrder * singularOrder;
        const unsigned count = perSector * galerkin::pairSectorCount(job.permutation.topology);
        for (unsigned sampleIndex = threadIdx.x; sampleIndex < count; sampleIndex += blockDim.x) {
            unsigned tuple = sampleIndex % perSector;
            const unsigned w = tuple % singularOrder;
            tuple /= singularOrder;
            const unsigned v = tuple % singularOrder;
            tuple /= singularOrder;
            const unsigned u = tuple % singularOrder;
            tuple /= singularOrder;
            const unsigned r = tuple;
            auto sample = galerkin::canonicalPairTransform(
                job.permutation.topology, sampleIndex / perSector, gauss[r].position,
                gauss[u].position, gauss[v].position, gauss[w].position);
            sample.weight *= gauss[r].weight * gauss[u].weight * gauss[v].weight * gauss[w].weight;
            sample = galerkin::permutePairPoint(sample, job.permutation);
            pairSample(moments, test, trial, sample,
                       pairDelta(test, trial, sample.test, sample.trial), k,
                       job.permutation.topology == galerkin::PairTopology::coincident, error);
        }
    } else if (job.patchCount) {
        for (unsigned patchIndex = 0; patchIndex < job.patchCount; ++patchIndex) {
            const auto &patch = patches[job.firstPatch + patchIndex];
            const bool promoted = patch.degree == 14;
            const auto *rule = nearRule + (promoted ? baseNearCount : 0);
            const unsigned count = promoted ? highNearCount : baseNearCount;
            const unsigned points = galerkin::nearPatchPointCount(patch, count);
            for (unsigned sampleIndex = threadIdx.x; sampleIndex < points;
                 sampleIndex += blockDim.x) {
                auto sample = galerkin::mapNearPatchPoint(rule, count, patch, sampleIndex,
                                                          job.swapBarycentric ? trial : test,
                                                          job.swapBarycentric ? test : trial);
                if (job.swapBarycentric)
                    for (unsigned i = 0; i < 3; ++i) {
                        const double value = sample.test[i];
                        sample.test[i] = sample.trial[i];
                        sample.trial[i] = value;
                    }
                pairSample(moments, test, trial, sample,
                           pairDelta(test, trial, sample.test, sample.trial), k, false, error);
            }
        }
    } else {
        MixedMoments mixed{};
        for (unsigned sampleIndex = threadIdx.x; sampleIndex < regularCount * regularCount;
             sampleIndex += blockDim.x) {
            galerkin::PairQuadraturePoint sample;
            const auto &x = regular[sampleIndex / regularCount];
            const auto &y = regular[sampleIndex % regularCount];
            for (unsigned i = 0; i < 3; ++i) {
                sample.test[i] = x.barycentric[i];
                sample.trial[i] = y.barycentric[i];
            }
            sample.weight = x.weight * y.weight;
            mixedSample(mixed, pairDelta(test, trial, sample.test, sample.trial), test, trial,
                        sample.test, sample.trial,
                        static_cast<float>(sample.weight * 4 * test.area * trial.area),
                        static_cast<float>(k), 1.f / static_cast<float>(k), false, error);
        }
        moments = promotedMoments(mixed);
    }
    reduceWarp(moments);
    __shared__ PairMoments partial[4];
    if ((threadIdx.x & 31) == 0)
        partial[threadIdx.x / 32] = moments;
    __syncthreads();
    if (threadIdx.x < 32) {
        PairMoments sum{};
        if (threadIdx.x < 4)
            sum = partial[threadIdx.x];
        reduceWarp(sum);
        if (threadIdx.x == 0) {
            for (unsigned i = 0; i < 3; ++i) {
                auto b = converted(sum.rhs[i]);
                if (job.test == job.trial)
                    b.y += test.area / (6 * k);
                sparseAdd(fluxMatrix, test.nodes[i], job.trial, b, error);
                for (unsigned j = 0; j < 3; ++j) {
                    auto a = pressureCoefficient(sum, test, trial, i, j, k);
                    if (job.test == job.trial)
                        a.x += test.area * (i == j ? 2.0 : 1.0) / 24.0;
                    sparseAdd(pressureMatrix, test.nodes[i], trial.nodes[j], a, error);
                    if (!isfinite(a.x) || !isfinite(a.y))
                        atomicExch(error, 2);
                }
                if (!isfinite(b.x) || !isfinite(b.y))
                    atomicExch(error, 2);
            }
        }
    }
}

template <bool fluxMode>
__global__ void applySparse(SparseMatrix matrix, unsigned pressureCount, const Z *pressure,
                            const cuFloatComplex *flux, Z *output) {
    const unsigned row = blockIdx.x * 4 + threadIdx.x / 32;
    if (row >= pressureCount)
        return;
    Z sum{};
    for (auto entry = matrix.offsets[row] + (threadIdx.x & 31u); entry < matrix.offsets[row + 1];
         entry += 32) {
        const auto column = matrix.columns[entry];
        if constexpr (fluxMode) {
            const auto q = flux[column];
            sum = zadd(sum, zproduct(matrix.values[entry], make_cuDoubleComplex(q.x, q.y)));
        } else {
            sum = zadd(sum, zproduct(matrix.values[entry], pressure[column]));
        }
    }
    for (unsigned offset = 16; offset; offset /= 2) {
        sum.x += __shfl_down_sync(0xffffffffu, sum.x, offset);
        sum.y += __shfl_down_sync(0xffffffffu, sum.y, offset);
    }
    if ((threadIdx.x & 31) == 0)
        output[row] = sum;
}

template <bool fluxMode>
__global__ void directFarAudit(const Panel *panels, const Node *nodes, const unsigned *panelOrder,
                               const AuditJob *jobs, std::size_t first,
                               const galerkin::TriangleQuadraturePoint *rule, unsigned count,
                               double k, const Z *pressure, const cuFloatComplex *flux, Z *output,
                               int *error) {
    const auto job = jobs[first + blockIdx.x];
    const auto &test = panels[job.test];
    const auto &source = nodes[job.trialNode];
    Z sum[3]{};
    for (unsigned index = threadIdx.x; index < source.panelCount; index += blockDim.x) {
        const unsigned trialIndex = panelOrder[source.firstPanel + index];
        if constexpr (fluxMode) {
            // A prescribed zero flux contributes exactly zero to this linear
            // audit. Skip only exact zeros, including inactive driven ports;
            // nonzero sources retain the full original quadrature contract.
            const auto q = flux[trialIndex];
            if (q.x == 0 && q.y == 0)
                continue;
        }
        const auto &trial = panels[trialIndex];
        PairMoments moments{};
        regularMoments(moments, test, trial, rule, count, k, error);
        for (unsigned i = 0; i < 3; ++i) {
            if constexpr (fluxMode) {
                const auto q = flux[trialIndex];
                sum[i] = zadd(sum[i],
                              zproduct(converted(moments.rhs[i]), make_cuDoubleComplex(q.x, q.y)));
            } else {
                for (unsigned j = 0; j < 3; ++j)
                    sum[i] =
                        zadd(sum[i], zproduct(pressureCoefficient(moments, test, trial, i, j, k),
                                              pressure[trial.nodes[j]]));
            }
        }
    }
    for (unsigned i = 0; i < 3; ++i) {
        for (unsigned offset = 16; offset; offset /= 2) {
            sum[i].x += __shfl_down_sync(0xffffffffu, sum[i].x, offset);
            sum[i].y += __shfl_down_sync(0xffffffffu, sum[i].y, offset);
        }
        if ((threadIdx.x & 31) == 0)
            atomicZ(output + test.nodes[i], sum[i]);
    }
}

template <bool fluxMode>
__global__ void particleToMultipole(const Panel *panels, const Node *nodes,
                                    const unsigned *leafNodes, const unsigned *panelOrder,
                                    unsigned leafCount, const double *weights, unsigned order,
                                    const galerkin::TriangleQuadraturePoint *rule, unsigned count,
                                    const Z *pressure, const cuFloatComplex *flux, Z *moments) {
    const unsigned coefficients = order * order * order;
    const auto flat = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (flat >= static_cast<std::size_t>(leafCount) * coefficients)
        return;
    const unsigned leaf = leafNodes[flat / coefficients], point = flat % coefficients;
    const unsigned ix = point % order, iy = point / order % order, iz = point / (order * order);
    const auto &node = nodes[leaf];
    Z result[channels]{};
    for (unsigned member = 0; member < node.panelCount; ++member) {
        const unsigned panelId = panelOrder[node.firstPanel + member];
        const auto &panel = panels[panelId];
        Z curl[3]{};
        if constexpr (!fluxMode)
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned axis = 0; axis < 3; ++axis)
                    curl[axis] = zadd(curl[axis], zscale(pressure[panel.nodes[i]],
                                                         coordinate(panel.curls[i], axis)));
        for (unsigned sample = 0; sample < count; ++sample) {
            const auto *w =
                weights + (static_cast<std::size_t>(panelId) * count + sample) * 3 * order;
            const double factor =
                w[ix] * w[order + iy] * w[2 * order + iz] * (2 * panel.area * rule[sample].weight);
            if constexpr (fluxMode) {
                const auto q = flux[panelId];
                result[0] = zadd(result[0], make_cuDoubleComplex(factor * q.x, factor * q.y));
            } else {
                Z p{};
                for (unsigned i = 0; i < 3; ++i)
                    p = zadd(p, zscale(pressure[panel.nodes[i]], rule[sample].barycentric[i]));
                for (unsigned axis = 0; axis < 3; ++axis) {
                    result[axis] =
                        zadd(result[axis], zscale(p, factor * coordinate(panel.normal, axis)));
                    result[3 + axis] = zadd(result[3 + axis], zscale(curl[axis], factor));
                }
            }
        }
    }
    auto *out = moments + (static_cast<std::size_t>(leaf) * coefficients + point) * channels;
    for (unsigned channel = 0; channel < channels; ++channel)
        out[channel] = result[channel];
}

// Separability reduces each tensor translation from O(p^6) to three O(p^4)
// contractions. The child-indexed scratch arrays prevent inter-level races.
template <unsigned axis, bool upward, bool fluxMode>
__global__ void translateAxis(const Node *nodes, const unsigned *level, unsigned levelCount,
                              unsigned order, const double *transfers, const Z *input, Z *output) {
    const unsigned coefficients = order * order * order;
    const auto flat = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (flat >= static_cast<std::size_t>(levelCount) * coefficients)
        return;
    const unsigned child = level[flat / coefficients], point = flat % coefficients;
    const unsigned stride = axis == 0 ? 1 : axis == 1 ? order : order * order;
    const unsigned coordinate = point / stride % order;
    const unsigned sourceNode = !upward && axis == 0 ? nodes[child].parent : child;
    const Z *source = input + static_cast<std::size_t>(sourceNode) * coefficients * channels;
    const double *transfer =
        transfers + (static_cast<std::size_t>(child) * 3 + axis) * order * order;
    Z result[channels]{};
    constexpr unsigned count = upward ? (fluxMode ? 1u : 6u) : (fluxMode ? 4u : 7u);
    for (unsigned j = 0; j < order; ++j) {
        const unsigned sourcePoint = point - coordinate * stride + j * stride;
        const double weight =
            upward ? transfer[coordinate * order + j] : transfer[j * order + coordinate];
        for (unsigned channel = 0; channel < count; ++channel)
            result[channel] =
                zadd(result[channel], zscale(source[sourcePoint * channels + channel], weight));
    }
    auto *out = output + (static_cast<std::size_t>(child) * coefficients + point) * channels;
    for (unsigned channel = 0; channel < channels; ++channel)
        if constexpr (!upward && axis == 2)
            out[channel] = zadd(out[channel], result[channel]);
        else
            out[channel] = result[channel];
}

__global__ void gatherMultipoles(const Node *nodes, const unsigned *level, unsigned levelCount,
                                 unsigned coefficients, const Z *translated, Z *moments) {
    const auto flat = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (flat >= static_cast<std::size_t>(levelCount) * coefficients)
        return;
    const unsigned parent = level[flat / coefficients], point = flat % coefficients;
    const auto &node = nodes[parent];
    if (node.child[0] == UINT32_MAX)
        return;
    const auto *first =
        translated + (static_cast<std::size_t>(node.child[0]) * coefficients + point) * channels;
    const auto *second =
        translated + (static_cast<std::size_t>(node.child[1]) * coefficients + point) * channels;
    auto *out = moments + (static_cast<std::size_t>(parent) * coefficients + point) * channels;
    for (unsigned channel = 0; channel < channels; ++channel)
        out[channel] = zadd(first[channel], second[channel]);
}

template <bool fluxMode, bool mixed>
__global__ void multipoleToLocal(const Node *nodes, unsigned nodeCount, const std::size_t *offsets,
                                 const unsigned *sources, unsigned order, const double *abscissae,
                                 double k, const Z *moments, Z *locals, std::size_t firstPoint,
                                 std::size_t pointCount, int *error) {
    const unsigned coefficients = order * order * order;
    const auto batchPoint = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto flat = firstPoint + batchPoint;
    if (batchPoint >= pointCount || flat >= static_cast<std::size_t>(nodeCount) * coefficients)
        return;
    const unsigned target = flat / coefficients, point = flat % coefficients;
    const Vec3 x = interpolationPoint(nodes[target], abscissae, point, order);
    MultipoleSum<mixed> result[channels]{};
    for (auto interaction = offsets[target]; interaction < offsets[target + 1]; ++interaction) {
        const unsigned source = sources[interaction];
        const Z *multipole = moments + static_cast<std::size_t>(source) * coefficients * channels;
        for (unsigned q = 0; q < coefficients; ++q) {
            const auto delta = sub(x, interpolationPoint(nodes[source], abscissae, q, order));
            const double r2 = dot(delta, delta), inverseR = rsqrt(r2), r = r2 * inverseR;
            double sine, cosine;
            sincos(k * r, &sine, &cosine);
            const Z green = make_cuDoubleComplex(cosine * inverseR * .07957747154594766788,
                                                 -sine * inverseR * .07957747154594766788);
            // grad_x G = -(1+ikr) G (x-y)/r^2 for exp(+i omega t).
            const Z radial = make_cuDoubleComplex(-(green.x - k * r * green.y) / r2,
                                                  -(green.y + k * r * green.x) / r2);
            const Z *m = multipole + q * channels;
            if constexpr (fluxMode) {
                result[0].addProduct(green, m[0]);
                for (unsigned axis = 0; axis < 3; ++axis)
                    result[1 + axis].addProduct(zscale(radial, coordinate(delta, axis)), m[0]);
            } else {
                for (unsigned axis = 0; axis < 3; ++axis) {
                    result[0].addProduct(zscale(radial, coordinate(delta, axis)), m[axis]);
                    result[1 + axis].addProduct(green, m[axis]);
                    result[4 + axis].addProduct(green, m[3 + axis]);
                }
            }
        }
    }
    auto *out = locals + flat * channels;
    for (unsigned channel = 0; channel < channels; ++channel) {
        out[channel] = result[channel].value();
        if (!isfinite(out[channel].x) || !isfinite(out[channel].y))
            atomicExch(error, 2);
    }
}

template <bool fluxMode>
__global__ void localToParticle(const Panel *panels, unsigned panelCount, const unsigned *panelLeaf,
                                unsigned order, const double *weights,
                                const galerkin::TriangleQuadraturePoint *rule, unsigned count,
                                double k, const Z *locals, Z *output, int *error) {
    const unsigned panelId = blockIdx.x * blockDim.x + threadIdx.x;
    if (panelId >= panelCount)
        return;
    const auto &panel = panels[panelId];
    const unsigned coefficients = order * order * order;
    const Z *source =
        locals + static_cast<std::size_t>(panelLeaf[panelId]) * coefficients * channels;
    Z result[3]{};
    for (unsigned sample = 0; sample < count; ++sample) {
        Z value[channels]{};
        const auto *w = weights + (static_cast<std::size_t>(panelId) * count + sample) * 3 * order;
        for (unsigned q = 0; q < coefficients; ++q) {
            const unsigned ix = q % order, iy = q / order % order, iz = q / (order * order);
            const double weight = w[ix] * w[order + iy] * w[2 * order + iz];
            for (unsigned channel = 0; channel < (fluxMode ? 4u : 7u); ++channel)
                value[channel] =
                    zadd(value[channel], zscale(source[q * channels + channel], weight));
        }
        Z normal{};
        for (unsigned axis = 0; axis < 3; ++axis)
            normal = zadd(normal, zscale(value[1 + axis], coordinate(panel.normal, axis)));
        for (unsigned i = 0; i < 3; ++i) {
            const double weight = 2 * panel.area * rule[sample].weight;
            const double phi = rule[sample].barycentric[i];
            if constexpr (fluxMode) {
                // -S q + (i/k) K' q; jump mass is in the direct self block.
                result[i] =
                    zadd(result[i], zscale(zadd(zscale(value[0], -1), zimaginary(normal, 1 / k)),
                                           weight * phi));
            } else {
                Z curl{};
                for (unsigned axis = 0; axis < 3; ++axis)
                    curl = zadd(curl, zscale(value[4 + axis], coordinate(panel.curls[i], axis)));
                const Z integrand = zadd(zscale(zadd(value[0], zimaginary(normal, k)), phi),
                                         zimaginary(curl, -1 / k));
                result[i] = zadd(result[i], zscale(integrand, weight));
            }
        }
    }
    for (unsigned i = 0; i < 3; ++i) {
        if (!isfinite(result[i].x) || !isfinite(result[i].y))
            atomicExch(error, 2);
        atomicZ(output + panel.nodes[i], result[i]);
    }
}

unsigned gridSize(std::size_t count) {
    const auto blocks = (count + threads - 1) / threads;
    if (blocks > static_cast<std::size_t>(INT_MAX))
        throw std::overflow_error("Fast Multipole CUDA grid is too large.");
    return static_cast<unsigned>(blocks);
}
} // namespace

class FastMultipoleOperator::Implementation {
  public:
    Implementation(const galerkin::Space &space, const galerkin::InteractionPlan &plan,
                   const BemNumerics &numerics, double kValue, cudaStream_t cudaStream,
                   std::stop_token token);
    template <bool fluxMode>
    void apply(const Z *pressure, const cuFloatComplex *flux, Z *output, bool direct);
    void precondition(Z *values);
    void checkError() {
        checked(cudaGetLastError(), "FMM kernel launch");
        const auto code = error.download(stream).front();
        if (code == 4)
            throw std::logic_error("FMM sparse near graph omitted a panel contribution.");
        if (code)
            throw std::runtime_error(code == 1 ? "FMM integration encountered an invalid distance."
                                               : "FMM operation produced a non-finite value.");
        cancelled(stop);
    }
    std::size_t deviceBytes() const noexcept;
    std::size_t allocationBytes() const noexcept;
    unsigned order{}, coefficientCount{}, panelCount{}, pressureCount{}, nodeCount{}, leafCount{};
    unsigned regularCount{}, baseNearCount{}, highNearCount{};
    unsigned blockCount{};
    std::size_t directCount{}, auditCount{};
    std::size_t peakBytes{};
    double k{};
    bool mixedMultipoles{};
    cudaStream_t stream;
    std::stop_token stop;
    FastMultipoleStatistics stats;
    std::vector<std::pair<std::size_t, unsigned>> levelRanges;
    BlasHandle blas;
    Buffer<Panel> panels;
    Buffer<Node> nodes;
    Buffer<unsigned> leaves, panelOrder, panelLeaf, levelNodes, farSources, pressureColumns,
        fluxColumns, pressurePermutation;
    Buffer<std::size_t> pressureOffsets, fluxOffsets, farOffsets;
    Buffer<DirectJob> directJobs;
    Buffer<AuditJob> auditJobs;
    Buffer<galerkin::TriangleQuadraturePoint> regularRule, nearRules;
    Buffer<galerkin::GaussPoint> gauss;
    Buffer<galerkin::PairQuadraturePatch> patches;
    Buffer<double> abscissae, weights, transfers;
    Buffer<Z> moments, locals, scratchFirst, scratchSecond, pressureValues, fluxValues,
        blockFactors, blockRhs;
    Buffer<Z *> blockPointers, rhsPointers;
    Buffer<int> error, blockPivots, blockInfo;
};

FastMultipoleOperator::Implementation::Implementation(const galerkin::Space &space,
                                                      const galerkin::InteractionPlan &plan,
                                                      const BemNumerics &numerics, double kValue,
                                                      cudaStream_t cudaStream,
                                                      std::stop_token token)
    : order(numerics.fastMultipole.interpolationOrder), coefficientCount(order * order * order),
      // Below 1e-5 the inner Krylov target approaches FP32 amplitude precision.
      // That regime keeps all M2L products and sums in FP64 as well.
      k(kValue), mixedMultipoles(numerics.relativeResidualTolerance >= 1e-5), stream(cudaStream),
      stop(token) {
    cancelled(stop);
    if (!(k > 0) || !std::isfinite(k) || space.vertices.empty() ||
        space.vertices.size() > UINT32_MAX)
        throw std::invalid_argument("Invalid Fast Multipole frequency or space.");
    std::size_t freeBeforeTree{}, deviceTotal{};
    checked(cudaMemGetInfo(&freeBeforeTree, &deviceTotal), "FMM hierarchy memory query");
    // Count mandatory expansion, direct and audit storage before any pair list
    // is materialized. The final exact preflight also includes the sparse graphs.
    fmm::TreeMemoryBudget treeBudget;
    treeBudget.maximumBytes = freeBeforeTree - freeBeforeTree / 5;
    treeBudget.nodeBytes = productSize(productSize(channels, coefficientCount), 4 * sizeof(Z));
    treeBudget.directPairBytes = sizeof(DirectJob);
    treeBudget.farPairBytes = sizeof(unsigned);
    treeBudget.auditJobBytes = sizeof(AuditJob);
    treeBudget.fixedBytes = productSize(plan.patches.size(), sizeof(galerkin::PairQuadraturePatch));
    const auto tree = fmm::makeTree(space, numerics.fastMultipole.leafPanelCapacity, order,
                                    numerics.fastMultipole.separationRatio,
                                    numerics.nearDistanceRatio, k, stop, treeBudget);
    panelCount = static_cast<unsigned>(space.panels.size());
    pressureCount = static_cast<unsigned>(space.vertices.size());
    blockCount = pressureCount / preconditionerWidth + (pressureCount % preconditionerWidth != 0);
    const auto paddedPressureCount = productSize(blockCount, preconditionerWidth);
    const auto blockEntries = productSize(paddedPressureCount, preconditionerWidth);
    // First occurrence along the spatial tree traversal assigns each shared P1
    // node to exactly one block, independent of arbitrary mesh node numbering.
    std::vector<unsigned> permutation;
    permutation.reserve(pressureCount);
    std::vector<unsigned char> assigned(pressureCount);
    for (const auto panel : tree.panelOrder) {
        if ((permutation.size() & 255u) == 0)
            cancelled(stop);
        for (const auto node : space.panels[panel].nodes) {
            if (node >= pressureCount)
                throw std::logic_error("FMM block preconditioner encountered an invalid DOF.");
            if (!assigned[node]) {
                assigned[node] = 1;
                permutation.push_back(node);
            }
        }
    }
    if (permutation.size() != pressureCount)
        throw std::logic_error("FMM block preconditioner did not cover all pressure DOFs.");
    nodeCount = static_cast<unsigned>(tree.nodes.size());
    leafCount = static_cast<unsigned>(tree.leafNodes.size());
    directCount = tree.directInteractions.size();
    stats = {nodeCount, leafCount, tree.farInteractions.size(), directCount,
             tree.representedFarPanelPairs};
    const auto nearLess = [](const auto &a, const auto &b) {
        return a.test == b.test ? a.trial < b.trial : a.test < b.test;
    };
    // Verify the conservative MAC against the actual adaptive integration plan.
    for (const auto &job : plan.near)
        if (!std::binary_search(tree.directInteractions.begin(), tree.directInteractions.end(),
                                fmm::Interaction{job.test, job.trial}, nearLess))
            throw std::logic_error("Fast Multipole hierarchy omitted an adaptive near pair.");
    for (const auto &job : plan.touching)
        if (!std::binary_search(tree.directInteractions.begin(), tree.directInteractions.end(),
                                fmm::Interaction{job.test, job.trial}, nearLess))
            throw std::logic_error("Fast Multipole hierarchy omitted a singular pair.");
    auto regular = galerkin::triangleRule(numerics.quadratureOrder);
    auto high = galerkin::triangleRule(8);
    const auto promoted = galerkin::triangleRule(14);
    regularCount = static_cast<unsigned>(regular.size());
    baseNearCount = static_cast<unsigned>(high.size());
    highNearCount = static_cast<unsigned>(promoted.size());
    high.insert(high.end(), promoted.begin(), promoted.end());
    const auto gaussHost = galerkin::gaussLegendre01(numerics.singularOrder);
    const auto chebyshev = fmm::chebyshevNodes(order);
    std::vector<Panel> packedPanels;
    packedPanels.reserve(panelCount);
    for (const auto &source : space.panels) {
        Panel panel{};
        panel.centroid = source.centroid;
        panel.normal = source.normal;
        panel.area = source.area;
        panel.edge = source.maximumEdgeLength;
        for (unsigned i = 0; i < 3; ++i) {
            panel.vertices[i] = source.vertices[i];
            panel.curls[i] = source.curls[i];
            panel.nodes[i] = source.nodes[i];
        }
        packedPanels.push_back(panel);
    }
    std::vector<Node> packedNodes;
    packedNodes.reserve(nodeCount);
    for (const auto &node : tree.nodes)
        packedNodes.push_back({node.center,
                               node.halfWidth,
                               node.parent,
                               {node.children[0], node.children[1]},
                               node.firstPanel,
                               node.panelCount});
    std::vector<DirectJob> jobs;
    jobs.reserve(directCount);
    std::vector<std::size_t> dOffsets(static_cast<std::size_t>(panelCount) + 1);
    for (const auto &interaction : tree.directInteractions) {
        DirectJob job;
        job.test = interaction.test;
        job.trial = interaction.trial;
        job.permutation = galerkin::classifyPair(space.panels[job.test], space.panels[job.trial]);
        const auto near =
            std::lower_bound(plan.near.begin(), plan.near.end(), interaction, nearLess);
        if (near != plan.near.end() && near->test == job.test && near->trial == job.trial) {
            job.firstPatch = near->firstPatch;
            job.patchCount = near->patchCount;
            job.swapBarycentric = near->swapBarycentric;
        }
        jobs.push_back(job);
        ++dOffsets[static_cast<std::size_t>(job.test) + 1];
    }
    for (std::size_t i = 1; i < dOffsets.size(); ++i)
        dOffsets[i] += dOffsets[i - 1];
    // Consolidate all panel contributions into the actual P1 and DP0 sparsity
    // graphs. Keeping twelve complex entries per ordered panel pair would
    // duplicate shared pressure DOFs many times and waste device memory.
    std::vector<std::size_t> pOffsets(static_cast<std::size_t>(pressureCount) + 1),
        qOffsets(static_cast<std::size_t>(pressureCount) + 1);
    std::vector<unsigned> pColumns, qColumns, rowPanels, rowNodes;
    for (unsigned row = 0; row < pressureCount; ++row) {
        if ((row & 127u) == 0)
            cancelled(stop);
        rowPanels.clear();
        rowNodes.clear();
        for (auto incidence = space.incidenceOffsets[row];
             incidence < space.incidenceOffsets[row + 1]; ++incidence) {
            const auto test = space.incidences[incidence].panel;
            for (auto pair = dOffsets[test]; pair < dOffsets[test + 1]; ++pair)
                rowPanels.push_back(jobs[pair].trial);
        }
        std::sort(rowPanels.begin(), rowPanels.end());
        rowPanels.erase(std::unique(rowPanels.begin(), rowPanels.end()), rowPanels.end());
        for (const auto trial : rowPanels)
            for (const auto node : space.panels[trial].nodes)
                rowNodes.push_back(node);
        std::sort(rowNodes.begin(), rowNodes.end());
        rowNodes.erase(std::unique(rowNodes.begin(), rowNodes.end()), rowNodes.end());
        const auto columns = pColumns.size() + qColumns.size() + rowNodes.size() + rowPanels.size();
        const auto graphBytes = productSize(columns, sizeof(unsigned) + sizeof(Z));
        if (graphBytes > treeBudget.maximumBytes - tree.constructionMinimumBytes)
            throw std::runtime_error(
                "Fast Multipole sparse-graph counting exceeds the GPU memory budget before column "
                "materialization: nodes=" +
                std::to_string(nodeCount) + ", direct=" + std::to_string(directCount) +
                ", far=" + std::to_string(tree.farInteractions.size()) +
                ", audit=" + std::to_string(tree.auditJobCount) +
                ", graphBytes=" + std::to_string(graphBytes) + ". Mesh was not changed.");
        qColumns.insert(qColumns.end(), rowPanels.begin(), rowPanels.end());
        qOffsets[static_cast<std::size_t>(row) + 1] = qColumns.size();
        pColumns.insert(pColumns.end(), rowNodes.begin(), rowNodes.end());
        pOffsets[static_cast<std::size_t>(row) + 1] = pColumns.size();
    }
    std::vector<std::size_t> fOffsets(static_cast<std::size_t>(nodeCount) + 1);
    std::vector<unsigned> sources;
    std::vector<AuditJob> audits;
    sources.reserve(tree.farInteractions.size());
    audits.reserve(tree.auditJobCount);
    for (const auto &interaction : tree.farInteractions) {
        ++fOffsets[static_cast<std::size_t>(interaction.test) + 1];
        sources.push_back(interaction.trial);
        const auto &target = tree.nodes[interaction.test];
        for (unsigned i = 0; i < target.panelCount; ++i)
            audits.push_back({tree.panelOrder[target.firstPanel + i], interaction.trial});
    }
    auditCount = audits.size();
    for (std::size_t i = 1; i < fOffsets.size(); ++i)
        fOffsets[i] += fOffsets[i - 1];
    std::vector<unsigned> levels;
    for (const auto &level : tree.levels) {
        levelRanges.emplace_back(levels.size(), static_cast<unsigned>(level.size()));
        levels.insert(levels.end(), level.begin(), level.end());
    }
    const auto interpolationSize = productSize(productSize(panelCount, regularCount), 3 * order);
    std::vector<double> interpolation(interpolationSize);
    const auto hostCoordinate = [](Vec3 p, unsigned axis) {
        return axis == 0 ? p.x : axis == 1 ? p.y : p.z;
    };
    for (unsigned panelId = 0; panelId < panelCount; ++panelId) {
        if ((panelId & 255u) == 0)
            cancelled(stop);
        const auto &panel = space.panels[panelId];
        const auto &leaf = tree.nodes[tree.panelLeaf[panelId]];
        for (unsigned sample = 0; sample < regularCount; ++sample) {
            Vec3 point{};
            for (unsigned i = 0; i < 3; ++i) {
                point.x += panel.vertices[i].x * regular[sample].barycentric[i];
                point.y += panel.vertices[i].y * regular[sample].barycentric[i];
                point.z += panel.vertices[i].z * regular[sample].barycentric[i];
            }
            for (unsigned axis = 0; axis < 3; ++axis) {
                const double coordinate =
                    (hostCoordinate(point, axis) - hostCoordinate(leaf.center, axis)) /
                    leaf.halfWidth;
                const auto w = fmm::interpolationWeights(chebyshev, coordinate);
                const auto first =
                    ((static_cast<std::size_t>(panelId) * regularCount + sample) * 3 + axis) *
                    order;
                std::copy(w.begin(), w.end(), interpolation.begin() + first);
            }
        }
    }
    std::vector<double> translation(productSize(nodeCount, 3 * order * order));
    for (unsigned child = 1; child < nodeCount; ++child) {
        const auto &node = tree.nodes[child];
        const auto &parent = tree.nodes[node.parent];
        for (unsigned axis = 0; axis < 3; ++axis)
            for (unsigned j = 0; j < order; ++j) {
                const double x =
                    (hostCoordinate(node.center, axis) - hostCoordinate(parent.center, axis) +
                     node.halfWidth * chebyshev[j]) /
                    parent.halfWidth;
                const auto w = fmm::interpolationWeights(chebyshev, x);
                for (unsigned i = 0; i < order; ++i)
                    translation[(static_cast<std::size_t>(child) * 3 + axis) * order * order +
                                i * order + j] = w[i];
            }
    }
    std::size_t bytes = 0;
    const auto add = [&bytes](std::size_t count, std::size_t width) {
        const auto value = productSize(count, width);
        if (value > std::numeric_limits<std::size_t>::max() - bytes)
            throw std::overflow_error("Fast Multipole workspace size overflow.");
        bytes += value;
    };
    const auto expansionCount = productSize(productSize(nodeCount, coefficientCount), channels);
    add(expansionCount, 4 * sizeof(Z));
    add(directCount, sizeof(DirectJob));
    add(panelCount, sizeof(Panel) + 2 * sizeof(unsigned));
    add(nodeCount, sizeof(Node) + 2 * sizeof(unsigned));
    add(pOffsets.size() + qOffsets.size() + fOffsets.size(), sizeof(std::size_t));
    add(pColumns.size() + qColumns.size(), sizeof(unsigned) + sizeof(Z));
    add(sources.size(), sizeof(unsigned));
    add(audits.size(), sizeof(AuditJob));
    add(interpolation.size() + translation.size() + chebyshev.size(), sizeof(double));
    add(regular.size() + high.size(), sizeof(galerkin::TriangleQuadraturePoint));
    add(gaussHost.size(), sizeof(galerkin::GaussPoint));
    add(plan.patches.size(), sizeof(galerkin::PairQuadraturePatch));
    add(pressureCount, sizeof(unsigned));
    add(blockEntries, sizeof(Z));
    add(paddedPressureCount, sizeof(Z) + sizeof(int));
    add(blockCount, 2 * sizeof(Z *) + sizeof(int));
    add(1, sizeof(int));
    std::size_t free{}, total{};
    checked(cudaMemGetInfo(&free, &total), "FMM memory query");
    if (bytes > free - free / 5)
        throw std::runtime_error("Fast Multipole near blocks and expansions exceed available GPU "
                                 "memory; mesh was not changed.");
    panels.assign(packedPanels, stream);
    nodes.assign(packedNodes, stream);
    leaves.assign(tree.leafNodes, stream);
    panelOrder.assign(tree.panelOrder, stream);
    panelLeaf.assign(tree.panelLeaf, stream);
    levelNodes.assign(levels, stream);
    directJobs.assign(jobs, stream);
    pressureOffsets.assign(pOffsets, stream);
    fluxOffsets.assign(qOffsets, stream);
    pressureColumns.assign(pColumns, stream);
    fluxColumns.assign(qColumns, stream);
    pressureValues.allocate(pColumns.size());
    fluxValues.allocate(qColumns.size());
    pressureValues.clear(stream);
    fluxValues.clear(stream);
    farOffsets.assign(fOffsets, stream);
    farSources.assign(sources, stream);
    auditJobs.assign(audits, stream);
    regularRule.assign(regular, stream);
    nearRules.assign(high, stream);
    gauss.assign(gaussHost, stream);
    patches.assign(plan.patches, stream);
    abscissae.assign(chebyshev, stream);
    weights.assign(interpolation, stream);
    transfers.assign(translation, stream);
    moments.allocate(expansionCount);
    locals.allocate(expansionCount);
    scratchFirst.allocate(expansionCount);
    scratchSecond.allocate(expansionCount);
    pressurePermutation.assign(permutation, stream);
    blockFactors.allocate(blockEntries);
    blockRhs.allocate(paddedPressureCount);
    blockPivots.allocate(paddedPressureCount);
    blockInfo.allocate(blockCount);
    std::vector<Z *> factorAddresses(blockCount), rhsAddresses(blockCount);
    for (unsigned block = 0; block < blockCount; ++block) {
        factorAddresses[block] = blockFactors.data() + static_cast<std::size_t>(block) *
                                                           preconditionerWidth *
                                                           preconditionerWidth;
        rhsAddresses[block] =
            blockRhs.data() + static_cast<std::size_t>(block) * preconditionerWidth;
    }
    blockPointers.assign(factorAddresses, stream);
    rhsPointers.assign(rhsAddresses, stream);
    error.allocate(1);
    error.clear(stream);
    for (std::size_t first = 0; first < directCount; first += 32768) {
        cancelled(stop);
        const auto count = static_cast<unsigned>(std::min<std::size_t>(32768, directCount - first));
        buildDirect<<<count, threads, 0, stream>>>(
            panels.data(), directJobs.data(), first, regularRule.data(), regularCount, gauss.data(),
            numerics.singularOrder, patches.data(), nearRules.data(), baseNearCount, highNearCount,
            k, {pressureOffsets.data(), pressureColumns.data(), pressureValues.data()},
            {fluxOffsets.data(), fluxColumns.data(), fluxValues.data()}, error.data());
        checked(cudaStreamSynchronize(stream), "FMM near block construction");
    }
    checkError();
    extractPreconditionerBlocks<<<gridSize(blockEntries), threads, 0, stream>>>(
        {pressureOffsets.data(), pressureColumns.data(), pressureValues.data()},
        pressurePermutation.data(), pressureCount, blockEntries, blockFactors.data(), error.data());
    checkError();
    blas.create(stream);
    checked(cublasZgetrfBatched(blas.get(), preconditionerWidth, blockPointers.data(),
                                preconditionerWidth, blockPivots.data(), blockInfo.data(),
                                static_cast<int>(blockCount)),
            "FMM block preconditioner LU factorization");
    const auto factorInfo = blockInfo.download(stream);
    for (unsigned block = 0; block < blockCount; ++block) {
        if (factorInfo[block] < 0)
            throw std::runtime_error("FMM block preconditioner LU argument is invalid.");
        if (factorInfo[block] > 0)
            throw std::runtime_error("FMM block preconditioner is singular in spatial block " +
                                     std::to_string(block) + ".");
    }
    validatePreconditionerFactors<<<gridSize(blockEntries), threads, 0, stream>>>(
        blockFactors.data(), blockEntries, error.data());
    checkError();
    // Temporary construction buffers are never needed during Krylov iterations.
    // Retain the observed allocation peak separately from the live byte count.
    peakBytes = allocationBytes();
    directJobs.release();
    gauss.release();
    patches.release();
    nearRules.release();
    blockInfo.release();
}

void FastMultipoleOperator::Implementation::precondition(Z *values) {
    cancelled(stop);
    if (!values)
        throw std::invalid_argument("FMM block preconditioner requires a CUDA vector.");
    error.clear(stream);
    const auto paddedCount = productSize(blockCount, preconditionerWidth);
    gatherPreconditionerRhs<<<gridSize(paddedCount), threads, 0, stream>>>(
        values, pressurePermutation.data(), pressureCount, paddedCount, blockRhs.data());
    checked(cudaGetLastError(), "FMM block preconditioner input permutation");
    int info{};
    checked(cublasZgetrsBatched(blas.get(), CUBLAS_OP_N, preconditionerWidth, 1,
                                blockPointers.data(), preconditionerWidth, blockPivots.data(),
                                rhsPointers.data(), preconditionerWidth, &info,
                                static_cast<int>(blockCount)),
            "FMM block preconditioner triangular solve");
    if (info)
        throw std::runtime_error("FMM block preconditioner triangular solve argument is invalid.");
    scatterPreconditionerSolution<<<gridSize(pressureCount), threads, 0, stream>>>(
        blockRhs.data(), pressurePermutation.data(), pressureCount, values, error.data());
    checkError();
}

template <bool fluxMode>
void FastMultipoleOperator::Implementation::apply(const Z *pressure, const cuFloatComplex *flux,
                                                  Z *output, bool direct) {
    cancelled(stop);
    if (!output || (fluxMode ? !flux : !pressure) || (!fluxMode && output == pressure))
        throw std::invalid_argument(
            "Fast Multipole application requires distinct CUDA input/output buffers.");
    checked(cudaMemsetAsync(output, 0, productSize(pressureCount, sizeof(Z)), stream),
            "FMM output clear");
    error.clear(stream);
    const SparseMatrix sparse =
        fluxMode
            ? SparseMatrix{fluxOffsets.data(), fluxColumns.data(), fluxValues.data()}
            : SparseMatrix{pressureOffsets.data(), pressureColumns.data(), pressureValues.data()};
    applySparse<fluxMode><<<(pressureCount + 3u) / 4u, threads, 0, stream>>>(
        sparse, pressureCount, pressure, flux, output);
    if (direct) {
        for (std::size_t first = 0; first < auditCount; first += 32768) {
            cancelled(stop);
            const auto count =
                static_cast<unsigned>(std::min<std::size_t>(32768, auditCount - first));
            directFarAudit<fluxMode><<<count, threads, 0, stream>>>(
                panels.data(), nodes.data(), panelOrder.data(), auditJobs.data(), first,
                regularRule.data(), regularCount, k, pressure, flux, output, error.data());
            checked(cudaStreamSynchronize(stream), "FMM direct audit batch");
        }
    } else if (stats.farInteractionCount) {
        particleToMultipole<fluxMode>
            <<<gridSize(static_cast<std::size_t>(leafCount) * coefficientCount), threads, 0,
               stream>>>(panels.data(), nodes.data(), leaves.data(), panelOrder.data(), leafCount,
                         weights.data(), order, regularRule.data(), regularCount, pressure, flux,
                         moments.data());
        for (std::size_t depth = levelRanges.size() - 1; depth > 0; --depth) {
            const auto [first, count] = levelRanges[depth];
            const auto blocks = gridSize(static_cast<std::size_t>(count) * coefficientCount);
            translateAxis<0, true, fluxMode><<<blocks, threads, 0, stream>>>(
                nodes.data(), levelNodes.data() + first, count, order, transfers.data(),
                moments.data(), scratchFirst.data());
            translateAxis<1, true, fluxMode><<<blocks, threads, 0, stream>>>(
                nodes.data(), levelNodes.data() + first, count, order, transfers.data(),
                scratchFirst.data(), scratchSecond.data());
            translateAxis<2, true, fluxMode><<<blocks, threads, 0, stream>>>(
                nodes.data(), levelNodes.data() + first, count, order, transfers.data(),
                scratchSecond.data(), scratchFirst.data());
            const auto [parentFirst, parentCount] = levelRanges[depth - 1];
            gatherMultipoles<<<gridSize(static_cast<std::size_t>(parentCount) * coefficientCount),
                               threads, 0, stream>>>(nodes.data(), levelNodes.data() + parentFirst,
                                                     parentCount, coefficientCount,
                                                     scratchFirst.data(), moments.data());
        }
        const auto targetPoints = static_cast<std::size_t>(nodeCount) * coefficientCount;
        // A bounded group of target interpolation points limits uninterruptible
        // GPU work. No partial expansion is exposed to the solver on cancellation.
        constexpr std::size_t maximumBatchPoints = 16384;
        for (std::size_t first = 0; first < targetPoints; first += maximumBatchPoints) {
            cancelled(stop);
            const auto count = std::min(maximumBatchPoints, targetPoints - first);
            if (mixedMultipoles)
                multipoleToLocal<fluxMode, true><<<gridSize(count), threads, 0, stream>>>(
                    nodes.data(), nodeCount, farOffsets.data(), farSources.data(), order,
                    abscissae.data(), k, moments.data(), locals.data(), first, count, error.data());
            else
                multipoleToLocal<fluxMode, false><<<gridSize(count), threads, 0, stream>>>(
                    nodes.data(), nodeCount, farOffsets.data(), farSources.data(), order,
                    abscissae.data(), k, moments.data(), locals.data(), first, count, error.data());
            checked(cudaStreamSynchronize(stream), "FMM translation batch");
        }
        for (std::size_t depth = 1; depth < levelRanges.size(); ++depth) {
            const auto [first, count] = levelRanges[depth];
            const auto blocks = gridSize(static_cast<std::size_t>(count) * coefficientCount);
            translateAxis<0, false, fluxMode><<<blocks, threads, 0, stream>>>(
                nodes.data(), levelNodes.data() + first, count, order, transfers.data(),
                locals.data(), scratchFirst.data());
            translateAxis<1, false, fluxMode><<<blocks, threads, 0, stream>>>(
                nodes.data(), levelNodes.data() + first, count, order, transfers.data(),
                scratchFirst.data(), scratchSecond.data());
            translateAxis<2, false, fluxMode><<<blocks, threads, 0, stream>>>(
                nodes.data(), levelNodes.data() + first, count, order, transfers.data(),
                scratchSecond.data(), locals.data());
        }
        localToParticle<fluxMode><<<gridSize(panelCount), threads, 0, stream>>>(
            panels.data(), panelCount, panelLeaf.data(), order, weights.data(), regularRule.data(),
            regularCount, k, locals.data(), output, error.data());
    }
    checkError();
}
std::size_t FastMultipoleOperator::Implementation::allocationBytes() const noexcept {
    return panels.bytes() + nodes.bytes() + leaves.bytes() + panelOrder.bytes() +
           panelLeaf.bytes() + levelNodes.bytes() + farSources.bytes() + pressureOffsets.bytes() +
           fluxOffsets.bytes() + pressureColumns.bytes() + fluxColumns.bytes() +
           pressureValues.bytes() + fluxValues.bytes() + farOffsets.bytes() + directJobs.bytes() +
           auditJobs.bytes() + regularRule.bytes() + nearRules.bytes() + gauss.bytes() +
           patches.bytes() + abscissae.bytes() + weights.bytes() + transfers.bytes() +
           moments.bytes() + locals.bytes() + scratchFirst.bytes() + scratchSecond.bytes() +
           pressurePermutation.bytes() + blockFactors.bytes() + blockRhs.bytes() +
           blockPointers.bytes() + rhsPointers.bytes() + blockPivots.bytes() + blockInfo.bytes() +
           error.bytes();
}
std::size_t FastMultipoleOperator::Implementation::deviceBytes() const noexcept {
    return std::max(peakBytes, allocationBytes());
}
FastMultipoleOperator::FastMultipoleOperator(const galerkin::Space &space,
                                             const galerkin::InteractionPlan &plan,
                                             const BemNumerics &numerics, double wavenumber,
                                             cudaStream_t stream, std::stop_token stop)
    : implementation_(
          std::make_unique<Implementation>(space, plan, numerics, wavenumber, stream, stop)) {}
FastMultipoleOperator::~FastMultipoleOperator() = default;
void FastMultipoleOperator::applyPressure(const Z *pressure, Z *output) {
    implementation_->apply<false>(pressure, nullptr, output, false);
}
void FastMultipoleOperator::applyFlux(const cuFloatComplex *flux, Z *output) {
    implementation_->apply<true>(nullptr, flux, output, false);
}
void FastMultipoleOperator::applyPressureDirect(const Z *pressure, Z *output) {
    implementation_->apply<false>(pressure, nullptr, output, true);
}
void FastMultipoleOperator::applyFluxDirect(const cuFloatComplex *flux, Z *output) {
    implementation_->apply<true>(nullptr, flux, output, true);
}
void FastMultipoleOperator::precondition(Z *values) {
    implementation_->precondition(values);
}
std::size_t FastMultipoleOperator::deviceBytes() const noexcept {
    return implementation_->deviceBytes();
}
const FastMultipoleStatistics &FastMultipoleOperator::statistics() const noexcept {
    return implementation_->stats;
}
} // namespace hornsculpt::cuda_detail
