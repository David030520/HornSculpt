#include "../galerkin/exterior_domain.hpp"
#include "../galerkin/interaction_plan.hpp"
#include "../galerkin/regular_qualification.hpp"
#include "../galerkin/space.hpp"
#include "../lem_coupling.hpp"
#include "../planning_frequencies.hpp"
#include "baffled_solver.hpp"
#include "galerkin_pair_moments.cuh"
#include "subdomain_solver.hpp"
#include "support_fmm_solve.hpp"
#include "workspace_size.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <hornsculpt/bem.hpp>
#include <hornsculpt/meshing.hpp>
#include <iomanip>
#include <limits>
#include <map>
#include <mutex>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace hornsculpt {
namespace {
using namespace cuda_detail;
using Complex = cuFloatComplex;
using AccurateComplex = cuDoubleComplex;
using Clock = std::chrono::steady_clock;
void cudaCheck(cudaError_t status, const char *operation) {
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
}
void blasCheck(cublasStatus_t status, const char *operation) {
    if (status != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error(std::string(operation) + ": cuBLAS status " +
                                 std::to_string(status));
}
void solverCheck(cusolverStatus_t status, const char *operation) {
    if (status != CUSOLVER_STATUS_SUCCESS)
        throw std::runtime_error(std::string(operation) + ": cuSOLVER status " +
                                 std::to_string(status));
}
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("BEM calculation cancelled.");
}
double milliseconds(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

template <class T> class DeviceArray {
  public:
    explicit DeviceArray(std::size_t count) : count_(count) {
        if (count_ > std::numeric_limits<std::size_t>::max() / sizeof(T))
            throw std::invalid_argument("CUDA allocation size overflow.");
        if (count_)
            cudaCheck(cudaMalloc(reinterpret_cast<void **>(&data_), count_ * sizeof(T)),
                      "CUDA allocation");
    }
    ~DeviceArray() {
        if (data_)
            cudaFree(data_);
    }
    DeviceArray(const DeviceArray &) = delete;
    DeviceArray &operator=(const DeviceArray &) = delete;
    T *data() noexcept {
        return data_;
    }
    const T *data() const noexcept {
        return data_;
    }
    std::size_t size() const noexcept {
        return count_;
    }
    std::size_t bytes() const noexcept {
        return count_ * sizeof(T);
    }
    void upload(const std::vector<T> &values, cudaStream_t stream) {
        if (values.size() != count_)
            throw std::logic_error("CUDA upload size mismatch.");
        if (count_)
            cudaCheck(
                cudaMemcpyAsync(data_, values.data(), bytes(), cudaMemcpyHostToDevice, stream),
                "CUDA upload");
    }
    std::vector<T> download(cudaStream_t stream) const {
        std::vector<T> result(count_);
        if (count_)
            cudaCheck(
                cudaMemcpyAsync(result.data(), data_, bytes(), cudaMemcpyDeviceToHost, stream),
                "CUDA download");
        cudaCheck(cudaStreamSynchronize(stream), "CUDA synchronization");
        return result;
    }
    void clear(cudaStream_t stream) {
        if (count_)
            cudaCheck(cudaMemsetAsync(data_, 0, bytes(), stream), "CUDA buffer clear");
    }

  private:
    T *data_{};
    std::size_t count_{};
};
// The source vector outlives this registration and is immutable until release.
// Only the address changes: all assembly kernels read the same patch coefficients.
class QuadraturePatchArray {
  public:
    using Patch = galerkin::PairQuadraturePatch;
    QuadraturePatchArray(const std::vector<Patch> &values, cudaStream_t stream, bool mapped) {
        if (!mapped) {
            device_ = std::make_unique<DeviceArray<Patch>>(values.size());
            device_->upload(values, stream);
            address_ = device_->data();
        } else if (!values.empty()) {
            int device{}, supported{}, canMap{};
            cudaCheck(cudaGetDevice(&device), "CUDA mapped-cache device query");
            cudaCheck(cudaDeviceGetAttribute(&supported, cudaDevAttrHostRegisterSupported, device),
                      "CUDA host registration query");
            cudaCheck(cudaDeviceGetAttribute(&canMap, cudaDevAttrCanMapHostMemory, device),
                      "CUDA host mapping query");
            if (!supported || !canMap)
                throw std::runtime_error("CUDA device does not support mapped quadrature storage.");
            auto *host = const_cast<Patch *>(values.data());
            cudaCheck(cudaHostRegister(host, checkedSizeProduct(values.size(), sizeof(Patch)),
                                       cudaHostRegisterMapped),
                      "CUDA quadrature host registration");
            void *address{};
            const auto status = cudaHostGetDevicePointer(&address, host, 0);
            if (status != cudaSuccess) {
                cudaHostUnregister(host);
                cudaCheck(status, "CUDA quadrature mapped pointer");
            }
            registered_ = host;
            address_ = static_cast<const Patch *>(address);
        }
    }
    ~QuadraturePatchArray() {
        if (registered_)
            cudaHostUnregister(registered_);
    }
    QuadraturePatchArray(const QuadraturePatchArray &) = delete;
    QuadraturePatchArray &operator=(const QuadraturePatchArray &) = delete;
    const Patch *data() const noexcept {
        return address_;
    }
    std::size_t bytes() const noexcept {
        return device_ ? device_->bytes() : 0;
    }

  private:
    std::unique_ptr<DeviceArray<Patch>> device_;
    Patch *registered_{};
    const Patch *address_{};
};
class Stream {
  public:
    Stream() {
        cudaCheck(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking), "CUDA stream");
    }
    ~Stream() {
        cudaStreamDestroy(stream_);
    }
    operator cudaStream_t() const noexcept {
        return stream_;
    }
    void synchronize() const {
        cudaCheck(cudaStreamSynchronize(stream_), "CUDA synchronization");
    }

  private:
    cudaStream_t stream_{};
};
class Blas {
  public:
    explicit Blas(cudaStream_t stream) {
        blasCheck(cublasCreate(&handle_), "cuBLAS initialization");
        const auto status = cublasSetStream(handle_, stream);
        if (status != CUBLAS_STATUS_SUCCESS) {
            cublasDestroy(handle_);
            blasCheck(status, "cuBLAS stream");
        }
    }
    ~Blas() {
        cublasDestroy(handle_);
    }
    operator cublasHandle_t() const noexcept {
        return handle_;
    }

  private:
    cublasHandle_t handle_{};
};
class DenseSolver {
  public:
    explicit DenseSolver(cudaStream_t stream) {
        solverCheck(cusolverDnCreate(&handle_), "cuSOLVER initialization");
        const auto status = cusolverDnSetStream(handle_, stream);
        if (status != CUSOLVER_STATUS_SUCCESS) {
            cusolverDnDestroy(handle_);
            solverCheck(status, "cuSOLVER stream");
        }
    }
    ~DenseSolver() {
        cusolverDnDestroy(handle_);
    }
    operator cusolverDnHandle_t() const noexcept {
        return handle_;
    }

  private:
    cusolverDnHandle_t handle_{};
};

class SolveWorkspace {
  public:
    SolveWorkspace(int pressureCount, int rhsCount, int observations, const WorkspaceSize &size,
                   cusolverDnHandle_t solver, bool accelerated)
        : excitations(rhsCount), observerCount(observations), fastMultipole(accelerated),
          tileCount(static_cast<unsigned>(size.tileCount)), matrix(size.matrixElements),
          factor(matrix.size()), rhs(size.rhsElements), correction(rhs.size()),
          flux(size.fluxElements), fields(size.fieldElements), pressure(rhs.size()),
          originalResidual(rhs.size()), fieldPartials(size.fieldPartialElements),
          observers(observations), pivots(accelerated ? 0 : pressureCount), info(1),
          integrationError(1) {
        if (accelerated)
            return;
        int required{};
        solverCheck(cusolverDnCgetrf_bufferSize(solver, pressureCount, pressureCount, factor.data(),
                                                pressureCount, &required),
                    "CUDA LU workspace query");
        std::size_t free{}, total{};
        cudaCheck(cudaMemGetInfo(&free, &total), "CUDA LU memory query");
        if (required < 0 || checkedSizeProduct(static_cast<std::size_t>(required),
                                               sizeof(Complex)) > deviceMemoryBudget(free))
            throw std::runtime_error("CUDA LU workspace exceeds available memory.");
        lu = std::make_unique<DeviceArray<Complex>>(static_cast<std::size_t>(required));
    }
    bool matches(std::size_t rhsCount, std::size_t observations, bool accelerated) const {
        return rhsCount == static_cast<std::size_t>(excitations) &&
               observations == static_cast<std::size_t>(observerCount) &&
               fastMultipole == accelerated;
    }
    std::size_t bytes() const {
        return matrix.bytes() + factor.bytes() + rhs.bytes() + pressure.bytes() +
               correction.bytes() + originalResidual.bytes() + flux.bytes() + fields.bytes() +
               fieldPartials.bytes() + observers.bytes() + pivots.bytes() + info.bytes() +
               integrationError.bytes() + (lu ? lu->bytes() : 0);
    }
    int excitations, observerCount;
    bool fastMultipole;
    unsigned tileCount;
    DeviceArray<Complex> matrix, factor, rhs, correction, flux, fields;
    DeviceArray<AccurateComplex> pressure, originalResidual;
    DeviceArray<DoubleComplex> fieldPartials;
    DeviceArray<Vec3> observers;
    DeviceArray<int> pivots, info, integrationError;
    std::unique_ptr<DeviceArray<Complex>> lu;
};

__global__ void assembleRegular(const Panel *__restrict__ panels, int fluxCount, int pressureCount,
                                const galerkin::TriangleQuadraturePoint *__restrict__ rules,
                                const galerkin::RegularPanelBounds *__restrict__ bounds,
                                const double *__restrict__ rowScale,
                                const Complex *__restrict__ flux, int excitationCount, double k,
                                const unsigned *__restrict__ nearOffsets,
                                const unsigned *__restrict__ nearTrials,
                                const std::uint8_t *__restrict__ wholePanelDegrees,
                                const galerkin::TriangleQuadraturePoint *__restrict__ nearRules,
                                Complex *__restrict__ matrix, Complex *__restrict__ rhs,
                                int *__restrict__ error, unsigned firstTrial, double tolerance,
                                unsigned minimumDegree) {
    const unsigned testIndex = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned trialIndex = firstTrial + blockIdx.y;
    if (testIndex >= static_cast<unsigned>(fluxCount))
        return;
    const auto &test = panels[testIndex];
    const auto &trial = panels[trialIndex];
    bool touching = false;
#pragma unroll
    for (unsigned i = 0; i < 3; ++i)
#pragma unroll
        for (unsigned j = 0; j < 3; ++j)
            touching |= test.nodes[i] == trial.nodes[j];
    if (touching)
        return;
    unsigned begin = nearOffsets[testIndex], end = nearOffsets[testIndex + 1];
    const auto stop = end;
    while (begin < end) {
        const unsigned middle = begin + (end - begin) / 2;
        if (nearTrials[middle] < trialIndex)
            begin = middle + 1;
        else
            end = middle;
    }
    if (begin < stop && nearTrials[begin] == trialIndex)
        return;
    const auto acceptedDegree =
        wholePanelDegrees
            ? wholePanelDegrees[galerkin::unorderedPanelPairIndex(testIndex, trialIndex)]
            : 0;
    if (acceptedDegree) {
        const auto *rule = nearRules + (acceptedDegree == 14 ? 16 : 0);
        const unsigned count = acceptedDegree == 14 ? 64 : 16;
        PairMoments moments{};
        for (unsigned index = 0; index < count * count; ++index) {
            const auto &x = rule[index / count], &y = rule[index % count];
            galerkin::PairQuadraturePoint sample{};
            for (unsigned i = 0; i < 3; ++i) {
                sample.test[i] = x.barycentric[i];
                sample.trial[i] = y.barycentric[i];
            }
            sample.weight = x.weight * y.weight;
            const auto delta = pairDelta(test, trial, sample.test, sample.trial);
            if (tolerance < 1e-5)
                pairSample(moments, test, trial, sample, delta, k, false, error);
            else {
                MixedMoments one{};
                mixedSample(one, delta, test, trial, sample.test, sample.trial,
                            static_cast<float>(sample.weight) *
                                static_cast<float>(4 * test.area * trial.area),
                            k, 1.f / static_cast<float>(k), false, error);
                const auto increment = promotedMoments(one);
                moments.green = plus(moments.green, increment.green);
                for (unsigned i = 0; i < 9; ++i)
                    moments.matrix[i] = plus(moments.matrix[i], increment.matrix[i]);
                for (unsigned i = 0; i < 3; ++i)
                    moments.rhs[i] = plus(moments.rhs[i], increment.rhs[i]);
            }
        }
        scatterPair(moments, test, trial, pressureCount, fluxCount, trialIndex, excitationCount, k,
                    rowScale, flux, matrix, rhs, error);
        return;
    }
    const auto degree = galerkin::qualifiedRegularDegreeUnchecked(
        bounds[testIndex], bounds[trialIndex], k, tolerance, minimumDegree);
    if (!degree) {
        atomicExch(error, 4);
        return;
    }
    const unsigned offset = degree == 4   ? 0
                            : degree == 6 ? 6
                                          : 18,
                   count = degree == 4   ? 6
                           : degree == 6 ? 12
                                         : 16;
    const auto moments =
        integrateRegularPairThread<true>(test, trial, rules + offset, count, k, error, 0, 1);
    scatterPair(moments, test, trial, pressureCount, fluxCount, trialIndex, excitationCount, k,
                rowScale, flux, matrix, rhs, error);
}
template <class Moments>
__device__ void completeBlock(Moments value, Moments *partial, const Panel &test,
                              const Panel &trial, int pressureCount, int fluxCount, int trialIndex,
                              int excitationCount, double k, const double *rowScale,
                              const Complex *flux, Complex *matrix, Complex *rhs, int *error) {
    reduceWarp(value);
    if ((threadIdx.x & 31) == 0)
        partial[threadIdx.x / 32] = value;
    __syncthreads();
    if (threadIdx.x < 32) {
        Moments total{};
        if (threadIdx.x < blockDim.x / 32)
            total = partial[threadIdx.x];
        reduceWarp(total);
        if (threadIdx.x == 0) {
            if constexpr (std::is_same_v<Moments, MixedMoments>)
                scatterMixed(total, test, trial, pressureCount, fluxCount, trialIndex,
                             excitationCount, static_cast<float>(k), rowScale, flux, matrix, rhs,
                             error);
            else
                scatterPair(total, test, trial, pressureCount, fluxCount, trialIndex,
                            excitationCount, k, rowScale, flux, matrix, rhs, error);
        }
    }
}
template <unsigned fixedOrder>
__global__ void assembleSingular(const Panel *panels, const PairJob *jobs,
                                 const galerkin::GaussPoint *gauss, unsigned runtimeOrder,
                                 int fluxCount, int pressureCount, const double *rowScale,
                                 const Complex *flux, int excitationCount, double k,
                                 Complex *matrix, Complex *rhs, int *error, std::size_t first) {
    const auto &job = jobs[first + blockIdx.x];
    const auto &test = panels[job.test];
    const auto &trial = panels[job.trial];
    const unsigned order = fixedOrder ? fixedOrder : runtimeOrder;
    const unsigned sectorPoints = order * order * order;
    const unsigned pointCount = galerkin::pairSectorCount(job.permutation.topology) * sectorPoints;
    MixedMoments moments{};
    for (unsigned sampleIndex = threadIdx.x; sampleIndex < pointCount; sampleIndex += blockDim.x) {
        unsigned tuple = sampleIndex % sectorPoints;
        const unsigned w = tuple % order;
        tuple /= order;
        const unsigned v = tuple % order;
        tuple /= order;
        const unsigned u = tuple;
        auto sample = galerkin::canonicalPairTransform(
            job.permutation.topology, sampleIndex / sectorPoints, 1.0, gauss[u].position,
            gauss[v].position, gauss[w].position);
        sample.weight *= gauss[u].weight * gauss[v].weight * gauss[w].weight;
        sample = galerkin::permutePairPoint(sample, job.permutation);
        const auto delta = pairDelta(test, trial, sample.test, sample.trial);
#pragma unroll
        for (unsigned radial = 0; radial < order; ++radial) {
            const double r = gauss[radial].position;
            double testBasis[3], trialBasis[3];
#pragma unroll
            for (unsigned i = 0; i < 3; ++i) {
                testBasis[i] = r * sample.test[i];
                trialBasis[i] = r * sample.trial[i];
            }
            testBasis[job.permutation.test[0]] += 1 - r;
            trialBasis[job.permutation.trial[0]] += 1 - r;
            mixedSample(moments, mul(delta, r), test, trial, testBasis, trialBasis,
                        static_cast<float>(sample.weight * r * r * r * gauss[radial].weight) *
                            static_cast<float>(4 * test.area * trial.area),
                        k, 1.f / static_cast<float>(k),
                        job.permutation.topology == galerkin::PairTopology::coincident, error);
        }
    }
    __shared__ MixedMoments partial[4];
    completeBlock(moments, partial, test, trial, pressureCount, fluxCount, job.trial,
                  excitationCount, k, rowScale, flux, matrix, rhs, error);
}
template <bool mixed>
__global__ void
assembleNear(const Panel *panels, const NearJob *jobs, const galerkin::PairQuadraturePatch *patches,
             const galerkin::TriangleQuadraturePoint *rules, int baseRuleCount, int highRuleCount,
             int fluxCount, int pressureCount, const double *rowScale, const Complex *flux,
             int excitationCount, double k, Complex *matrix, Complex *rhs, int *error,
             std::size_t first) {
    const auto &job = jobs[first + blockIdx.x];
    const auto &test = panels[job.test];
    const auto &trial = panels[job.trial];
    using Moments = std::conditional_t<mixed, MixedMoments, PairMoments>;
    Moments moments{};
    extern __shared__ NearSample samples[];
    auto *testSamples = samples;
    auto *trialSamples = samples + highRuleCount;
    for (unsigned patchIndex = 0; patchIndex < job.patchCount; ++patchIndex) {
        const auto &patch = patches[job.firstPatch + patchIndex];
        const bool promoted = patch.degree == 14;
        const auto *rule = rules + (promoted ? baseRuleCount : 0);
        const int ruleCount = promoted ? highRuleCount : baseRuleCount;
        if (patch.projection) {
            const unsigned points =
                galerkin::nearPatchPointCount(patch, static_cast<unsigned>(ruleCount));
            for (unsigned index = threadIdx.x; index < points; index += blockDim.x) {
                auto sample = galerkin::mapNearPatchPoint(
                    rule, static_cast<unsigned>(ruleCount), patch, index,
                    job.swapBarycentric ? trial : test, job.swapBarycentric ? test : trial);
                if (job.swapBarycentric)
                    for (unsigned i = 0; i < 3; ++i) {
                        const double swap = sample.test[i];
                        sample.test[i] = sample.trial[i];
                        sample.trial[i] = swap;
                    }
                const auto delta = pairDelta(test, trial, sample.test, sample.trial);
                if constexpr (mixed)
                    mixedSample(moments, delta, test, trial, sample.test, sample.trial,
                                static_cast<float>(sample.weight * 4 * test.area * trial.area),
                                static_cast<float>(k), 1.f / static_cast<float>(k), false, error);
                else
                    pairSample(moments, test, trial, sample, delta, k, false, error);
            }
            __syncthreads();
            continue;
        }
        const std::size_t perPatch = static_cast<std::size_t>(ruleCount) * ruleCount;
        if (threadIdx.x < static_cast<unsigned>(ruleCount)) {
            auto sample = galerkin::mapPairPatch(rule[threadIdx.x], rule[threadIdx.x], patch);
            if (job.swapBarycentric)
                for (unsigned i = 0; i < 3; ++i) {
                    const double value = sample.test[i];
                    sample.test[i] = sample.trial[i];
                    sample.trial[i] = value;
                }
            for (unsigned i = 0; i < 3; ++i) {
                testSamples[threadIdx.x].basis[i] = sample.test[i];
                trialSamples[threadIdx.x].basis[i] = sample.trial[i];
            }
            testSamples[threadIdx.x].local = localPoint(test, sample.test);
            trialSamples[threadIdx.x].local =
                add(sub(trial.vertices[0], test.vertices[0]), localPoint(trial, sample.trial));
        }
        __syncthreads();
        for (std::size_t point = threadIdx.x; point < perPatch; point += blockDim.x) {
            const auto x = point / ruleCount, y = point % ruleCount;
            const auto &firstSample = testSamples[x];
            const auto &secondSample = trialSamples[y];
            const auto delta = sub(firstSample.local, secondSample.local);
            const double weight = rule[x].weight * rule[y].weight * patch.areaFraction;
            if constexpr (mixed)
                mixedSample(moments, delta, test, trial, firstSample.basis, secondSample.basis,
                            static_cast<float>(weight) *
                                static_cast<float>(4 * test.area * trial.area),
                            k, 1.f / static_cast<float>(k), false, error);
            else {
                galerkin::PairQuadraturePoint sample{};
                for (unsigned i = 0; i < 3; ++i) {
                    sample.test[i] = firstSample.basis[i];
                    sample.trial[i] = secondSample.basis[i];
                }
                sample.weight = weight;
                pairSample(moments, test, trial, sample, delta, k, false, error);
            }
        }
        __syncthreads();
    }
    __shared__ Moments partial[4];
    completeBlock(moments, partial, test, trial, pressureCount, fluxCount, job.trial,
                  excitationCount, k, rowScale, flux, matrix, rhs, error);
}
__global__ void addMass(const Panel *panels, int panelCount, int pressureCount,
                        const double *rowScale, const Complex *flux, int excitations, double k,
                        Complex *matrix, Complex *rhs, int *error) {
    const unsigned index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= static_cast<unsigned>(panelCount))
        return;
    const auto &panel = panels[index];
    for (unsigned i = 0; i < 3; ++i) {
        const auto row = panel.nodes[i];
        for (unsigned j = 0; j < 3; ++j)
            atomicComplex(matrix + row + static_cast<std::size_t>(panel.nodes[j]) * pressureCount,
                          {rowScale[row] * panel.area * (i == j ? 2.0 : 1.0) / 24.0, 0}, error);
        const double coefficient = panel.area / (6.0 * k) * rowScale[row];
        for (int excitation = 0; excitation < excitations; ++excitation) {
            const auto q = flux[static_cast<std::size_t>(excitation) * panelCount + index];
            if (q.x == 0 && q.y == 0)
                continue;
            atomicComplex(rhs + row + static_cast<std::size_t>(excitation) * pressureCount,
                          {-coefficient * q.y, coefficient * q.x}, error);
        }
    }
}
__global__ void initializePressure(const Complex *initial, AccurateComplex *pressure,
                                   std::size_t count) {
    const auto index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index < count)
        pressure[index] = make_cuDoubleComplex(initial[index].x, initial[index].y);
}
__global__ void updatePressure(const Complex *correction, AccurateComplex *pressure,
                               std::size_t count) {
    const auto index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index < count) {
        pressure[index].x += static_cast<double>(correction[index].x);
        pressure[index].y += static_cast<double>(correction[index].y);
    }
}
// Mixed-precision refinement: FP32 LU is the correction solver, while b-Ap and
// the accumulated solution remain FP64. Recomputing b-Ap with FP32 GEMM can
// stall at its summation error even when the factors still give useful updates.
__global__ void trueResidual(const Complex *matrix, const Complex *rhs,
                             const AccurateComplex *pressure, const double *mass,
                             Complex *correction, AccurateComplex *original, int count) {
    const unsigned row = blockIdx.x * 32 + threadIdx.x;
    const auto offset = static_cast<std::size_t>(blockIdx.y) * count;
    double real = 0, imag = 0;
    if (row < static_cast<unsigned>(count))
        for (unsigned column = threadIdx.y; column < static_cast<unsigned>(count); column += 8) {
            const auto a = matrix[row + static_cast<std::size_t>(column) * count];
            const auto p = pressure[offset + column];
            real += static_cast<double>(a.x) * p.x - static_cast<double>(a.y) * p.y;
            imag += static_cast<double>(a.x) * p.y + static_cast<double>(a.y) * p.x;
        }
    __shared__ AccurateComplex partial[8][32];
    partial[threadIdx.y][threadIdx.x] = make_cuDoubleComplex(real, imag);
    __syncthreads();
    if (threadIdx.y == 0 && row < static_cast<unsigned>(count)) {
        for (unsigned lane = 1; lane < 8; ++lane) {
            real += partial[lane][threadIdx.x].x;
            imag += partial[lane][threadIdx.x].y;
        }
        const auto index = offset + row;
        real = static_cast<double>(rhs[index].x) - real;
        imag = static_cast<double>(rhs[index].y) - imag;
        correction[index] = make_cuFloatComplex(static_cast<float>(real), static_cast<float>(imag));
        original[index] = make_cuDoubleComplex(real * mass[row], imag * mass[row]);
    }
}
__device__ void storeFieldPartial(DoubleComplex value, DoubleComplex *partials) {
    for (unsigned offset = 16; offset > 0; offset /= 2) {
        value.real += __shfl_down_sync(0xffffffffU, value.real, offset);
        value.imag += __shfl_down_sync(0xffffffffU, value.imag, offset);
    }
    __shared__ DoubleComplex partial[4];
    if ((threadIdx.x & 31) == 0)
        partial[threadIdx.x / 32] = value;
    __syncthreads();
    if (threadIdx.x < 32) {
        DoubleComplex sum{};
        if (threadIdx.x < 4)
            sum = partial[threadIdx.x];
        for (unsigned offset = 16; offset > 0; offset /= 2) {
            sum.real += __shfl_down_sync(0xffffffffU, sum.real, offset);
            sum.imag += __shfl_down_sync(0xffffffffU, sum.imag, offset);
        }
        if (threadIdx.x == 0) {
            partials[blockIdx.x] = sum;
        }
    }
}

__global__ void evaluateField(const Panel *panels, int panelCount, const Vec3 *observers,
                              int observerCount, const galerkin::TriangleQuadraturePoint *rule,
                              int ruleCount, double k, const AccurateComplex *pressure,
                              int pressureCount, const Complex *flux, DoubleComplex *partials,
                              int excitations, unsigned tileCount, double tolerance,
                              unsigned maximumDepth, int *error) {
    const unsigned index = blockIdx.x / tileCount, tile = blockIdx.x % tileCount;
    if (index >= static_cast<unsigned>(observerCount * excitations))
        return;
    const unsigned excitation = index / observerCount, observer = index % observerCount;
    DoubleComplex value{};
    const unsigned end = min((tile + 1) * 256u, static_cast<unsigned>(panelCount));
    for (unsigned panel = tile * 256 + threadIdx.x; panel < end; panel += blockDim.x) {
        const auto offset = sub(observers[observer], panels[panel].centroid);
        const auto *nodal = pressure + static_cast<std::size_t>(excitation) * pressureCount;
        const auto q = flux[static_cast<std::size_t>(excitation) * panelCount + panel];
        if (dot(offset, offset) > 4 * panels[panel].edge * panels[panel].edge)
            value = plus(value, regularField(observers[observer], panels[panel], rule,
                                             static_cast<unsigned>(ruleCount), nodal, q, k, error));
        else
            value =
                plus(value, integrateNearField(observers[observer], panels[panel], rule, ruleCount,
                                               nodal, q, k, tolerance, maximumDepth, error));
    }
    storeFieldPartial(value, partials);
}

__global__ void evaluateFarField(const Panel *panels, int panelCount, const Vec3 *directions,
                                 int observerCount, Vec3 observationOrigin,
                                 double referenceDistance,
                                 const galerkin::TriangleQuadraturePoint *rule, int ruleCount,
                                 double k, const AccurateComplex *pressure, int pressureCount,
                                 const Complex *flux, DoubleComplex *partials, int excitations,
                                 unsigned tileCount, int *error) {
    const unsigned index = blockIdx.x / tileCount, tile = blockIdx.x % tileCount;
    if (index >= static_cast<unsigned>(observerCount * excitations))
        return;
    const unsigned excitation = index / observerCount, observer = index % observerCount;
    DoubleComplex value{};
    const unsigned end = min((tile + 1) * 256u, static_cast<unsigned>(panelCount));
    for (unsigned panel = tile * 256 + threadIdx.x; panel < end; panel += blockDim.x) {
        const auto *nodal = pressure + static_cast<std::size_t>(excitation) * pressureCount;
        const auto q = flux[static_cast<std::size_t>(excitation) * panelCount + panel];
        value = plus(value, asymptoticField(directions[observer], observationOrigin, panels[panel],
                                            rule, ruleCount, nodal, q, k, error));
    }
    double sine, cosine;
    sincos(k * referenceDistance, &sine, &cosine);
    value = product(value, {cosine, -sine});
    value.real /= referenceDistance;
    value.imag /= referenceDistance;
    storeFieldPartial(value, partials);
}
__global__ void finishField(const DoubleComplex *partials, Complex *field, std::size_t outputCount,
                            unsigned tileCount, int *error) {
    const auto index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= outputCount)
        return;
    DoubleComplex sum{};
    for (unsigned tile = 0; tile < tileCount; ++tile)
        sum = plus(sum, partials[index * tileCount + tile]);
    field[index] = make_cuFloatComplex(static_cast<float>(sum.real), static_cast<float>(sum.imag));
    if (!isfinite(field[index].x) || !isfinite(field[index].y))
        atomicExch(error, 2);
}
double vectorNorm(cublasHandle_t handle, int count, const AccurateComplex *vector) {
    double result{};
    blasCheck(cublasDznrm2(handle, count, vector, 1, &result), "CUDA residual norm");
    return result;
}
bool inside(const std::vector<Panel> &panels, Vec3 point) {
    double solidAngle = 0;
    for (const auto &panel : panels) {
        const auto edge0 = sub(panel.vertices[1], panel.vertices[0]);
        const auto edge1 = sub(panel.vertices[2], panel.vertices[0]);
        const auto relative = sub(point, panel.vertices[0]);
        const double coordinateScale = std::max(length(point), length(panel.vertices[0]));
        const double tolerance = std::max(
            panel.edge * 1e-10, coordinateScale * 32 * std::numeric_limits<double>::epsilon());
        if (std::abs(dot(relative, panel.normal)) <= tolerance) {
            const double aa = dot(edge0, edge0), ab = dot(edge0, edge1), bb = dot(edge1, edge1),
                         pa = dot(relative, edge0), pb = dot(relative, edge1),
                         determinant = aa * bb - ab * ab;
            const double u = (bb * pa - ab * pb) / determinant,
                         v = (aa * pb - ab * pa) / determinant;
            if (u >= -1e-10 && v >= -1e-10 && u + v <= 1 + 1e-10)
                return true;
        }
        const auto a = sub(panel.vertices[0], point), b = sub(panel.vertices[1], point),
                   c = sub(panel.vertices[2], point);
        const double denominator = length(a) * length(b) * length(c) + dot(a, b) * length(c) +
                                   dot(b, c) * length(a) + dot(c, a) * length(b);
        solidAngle += 2 * std::atan2(dot(a, cross(b, c)), denominator);
    }
    return std::abs(solidAngle) > 2 * std::numbers::pi;
}
void checkIntegrationError(const DeviceArray<int> &error, cudaStream_t stream) {
    const auto code = error.download(stream).front();
    if (code == 1)
        throw std::runtime_error("BEM quadrature encountered zero or invalid source distance.");
    if (code == 2)
        throw std::runtime_error("BEM calculation overflows complex FP32 storage.");
    if (code == 3)
        throw std::runtime_error("BEM field integration did not reach its near-field tolerance.");
    if (code == 4)
        throw std::runtime_error("BEM dense regular pair escaped its qualified adaptive plan.");
}
std::vector<Panel> pack(const galerkin::Space &space) {
    std::vector<Panel> panels;
    panels.reserve(space.panels.size());
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
        panels.push_back(panel);
    }
    return panels;
}
std::complex<double> hostComplex(Complex value) {
    return {value.x, value.y};
}
std::complex<double> hostComplex(AccurateComplex value) {
    return {value.x, value.y};
}
} // namespace

class CudaBemSolver::Implementation {
  public:
    explicit Implementation(SurfaceMesh value, CudaStorageOptions storageOptions,
                            std::optional<BemSubdomainModel> model = {},
                            std::optional<OrthogonalSymmetryPlanes> symmetry = {},
                            std::optional<ExteriorBaffle> baffle = {})
        : storage(storageOptions), domainModel(std::move(model)),
          symmetryModel(std::move(symmetry)), baffleModel(std::move(baffle)), origin{},
          mesh(centered(std::move(value), origin)),
          space(domainModel ? galerkin::buildSpace(mesh) : validatedSpace(mesh)),
          baffledPlanCache(mesh), panels(pack(space)), stream(), blas(stream), solver(stream),
          devicePanels(panels.size()), deviceMass(space.lumpedMass.size()),
          deviceRowScale(space.lumpedMass.size()) {
        if (storage.quadraturePatches != CudaCacheResidency::automatic &&
            storage.quadraturePatches != CudaCacheResidency::device &&
            storage.quadraturePatches != CudaCacheResidency::mappedHost)
            throw std::invalid_argument("Invalid CUDA quadrature storage policy.");
        devicePanels.upload(panels, stream);
        deviceMass.upload(space.lumpedMass, stream);
        auto scale = space.lumpedMass;
        for (auto &entry : scale)
            entry = 1.0 / entry;
        deviceRowScale.upload(scale, stream);
        stream.synchronize();
    }
    ~Implementation() {
        // Complete all reads before any mapped host range or source vector is released.
        cudaStreamSynchronize(stream);
    }
    static SurfaceMesh centered(SurfaceMesh value, Vec3 &origin) {
        if (value.vertices.empty())
            throw std::invalid_argument("BEM mesh has no vertices.");
        origin = value.vertices.front();
        for (auto &vertex : value.vertices) {
            vertex = sub(vertex, origin);
            if (!std::isfinite(vertex.x) || !std::isfinite(vertex.y) || !std::isfinite(vertex.z))
                throw std::invalid_argument("BEM mesh coordinates are not finite.");
        }
        return value;
    }
    static galerkin::Space validatedSpace(const SurfaceMesh &mesh) {
        const auto quality = inspectMesh(mesh);
        if (!quality.watertight())
            throw std::invalid_argument(
                "BEM requires a closed, outward-oriented, non-intersecting mesh.");
        galerkin::validateExteriorComponents(mesh);
        return galerkin::buildSpace(mesh);
    }
    void releaseDenseQuadrature() {
        deviceSingular.reset();
        deviceNear.reset();
        deviceNearOffsets.reset();
        deviceNearTrials.reset();
        deviceWholePanelDegrees.reset();
        devicePatches.reset();
        deviceNearRules.reset();
        deviceGauss.reset();
    }
    void prepareQuadrature(const BemNumerics &options, double maximumK, std::stop_token stop,
                           std::size_t workspaceBytes) {
        cancelled(stop);
        const bool dense = options.acceleration == BemAcceleration::dense;
        const bool samePlan = planValid && cachedAcceleration == options.acceleration &&
                              cached.quadratureOrder == options.quadratureOrder &&
                              cached.singularOrder == options.singularOrder &&
                              cached.nearDistanceRatio == options.nearDistanceRatio &&
                              cached.nearRelativeTolerance == options.nearRelativeTolerance &&
                              cached.nearMaximumDepth == options.nearMaximumDepth &&
                              cachedMaximumK == maximumK;
        if (samePlan && cacheValid && cachedAcceleration == options.acceleration) {
            std::size_t free{}, total{};
            cudaCheck(cudaMemGetInfo(&free, &total), "CUDA cached workspace memory query");
            if (workspaceBytes <= deviceMemoryBudget(free))
                return;
        }
        stream.synchronize();
        cacheValid = false;
        // Release incompatible caches before the capacity query and allocation.
        // FMM owns its construction jobs; retaining the dense copies wastes GPU
        // memory even after the FMM operator has released its temporary patches.
        releaseDenseQuadrature();
        const bool rebuildRegular = dense && (!deviceRegular || !deviceRegularBounds ||
                                              cached.quadratureOrder != options.quadratureOrder);
        if (rebuildRegular || !dense) {
            deviceRegular.reset();
            deviceRegularBounds.reset();
        }
        if (!samePlan) {
            planValid = false;
            galerkin::PairQuadratureSettings settings;
            settings.regularOrder = options.quadratureOrder;
            settings.singularOrder = options.singularOrder;
            settings.nearDistanceRatio = options.nearDistanceRatio;
            settings.relativeTolerance = options.nearRelativeTolerance;
            settings.maximumSubdivisionDepth = options.nearMaximumDepth;
            if (dense)
                interactionPlan = galerkin::makeInteractionPlan(space, settings, maximumK, stop);
            else
                interactionPlan = {};
            cached = options;
            cachedMaximumK = maximumK;
            planValid = true;
        }
        const auto &plan = interactionPlan;
        std::vector<PairJob> singular;
        std::vector<NearJob> near;
        std::vector<galerkin::TriangleQuadraturePoint> high, nearRules;
        std::vector<galerkin::GaussPoint> gauss;
        if (dense)
            high = galerkin::triangleRule(8);
        if (dense) {
            singular.reserve(plan.touching.size());
            for (const auto &job : plan.touching)
                singular.push_back({job.test, job.trial, job.permutation});
            near.reserve(plan.near.size());
            for (const auto &job : plan.near)
                near.push_back(
                    {job.test, job.trial, job.firstPatch, job.patchCount, job.swapBarycentric});
            nearRules = high;
            const auto promotedNearRule = galerkin::triangleRule(14);
            nearRules.insert(nearRules.end(), promotedNearRule.begin(), promotedNearRule.end());
            gauss = galerkin::gaussLegendre01(options.singularOrder);
        }
        std::vector<galerkin::TriangleQuadraturePoint> regularRule;
        std::vector<galerkin::RegularPanelBounds> regularBounds;
        if (rebuildRegular) {
            for (const unsigned degree : {4u, 6u, 8u}) {
                const auto rule = galerkin::triangleRule(degree);
                regularRule.insert(regularRule.end(), rule.begin(), rule.end());
            }
            regularBounds.reserve(space.panels.size());
            for (const auto &panel : space.panels)
                regularBounds.push_back(galerkin::makeRegularPanelBounds(panel));
        }
        std::size_t bytes{};
        const auto addBytes = [&bytes](std::size_t count, std::size_t width) {
            bytes = checkedSizeSum(bytes, checkedSizeProduct(count, width));
        };
        if (dense) {
            addBytes(singular.size(), sizeof(PairJob));
            addBytes(near.size(), sizeof(NearJob));
            addBytes(plan.patches.size(), sizeof(galerkin::PairQuadraturePatch));
            addBytes(plan.nearOffsets.size(), sizeof(unsigned));
            addBytes(plan.nearTrials.size(), sizeof(unsigned));
            addBytes(plan.wholePanelDegrees.size(), sizeof(std::uint8_t));
            addBytes(nearRules.size(), sizeof(galerkin::TriangleQuadraturePoint));
            addBytes(gauss.size(), sizeof(galerkin::GaussPoint));
        }
        if (rebuildRegular) {
            addBytes(regularRule.size(), sizeof(galerkin::TriangleQuadraturePoint));
            addBytes(regularBounds.size(), sizeof(galerkin::RegularPanelBounds));
        }
        std::size_t free{}, total{};
        cudaCheck(cudaMemGetInfo(&free, &total), "CUDA quadrature memory query");
        const auto patchBytes =
            dense ? checkedSizeProduct(plan.patches.size(), sizeof(galerkin::PairQuadraturePatch))
                  : 0;
        const bool mappedPatches =
            dense && (storage.quadraturePatches == CudaCacheResidency::mappedHost ||
                      (storage.quadraturePatches == CudaCacheResidency::automatic &&
                       !deviceCacheAndWorkspaceFit(free, bytes, workspaceBytes)));
        if (mappedPatches)
            bytes -= patchBytes;
        if (bytes > deviceMemoryBudget(free))
            throw std::runtime_error(
                "CUDA near-pair cache exceeds available device memory; mesh was not changed.");
        if (dense) {
            deviceSingular = std::make_unique<DeviceArray<PairJob>>(singular.size());
            deviceSingular->upload(singular, stream);
            deviceNear = std::make_unique<DeviceArray<NearJob>>(near.size());
            deviceNear->upload(near, stream);
            deviceNearOffsets = std::make_unique<DeviceArray<unsigned>>(plan.nearOffsets.size());
            deviceNearOffsets->upload(plan.nearOffsets, stream);
            deviceNearTrials = std::make_unique<DeviceArray<unsigned>>(plan.nearTrials.size());
            deviceNearTrials->upload(plan.nearTrials, stream);
            deviceWholePanelDegrees =
                std::make_unique<DeviceArray<std::uint8_t>>(plan.wholePanelDegrees.size());
            deviceWholePanelDegrees->upload(plan.wholePanelDegrees, stream);
            devicePatches =
                std::make_unique<QuadraturePatchArray>(plan.patches, stream, mappedPatches);
            deviceNearRules =
                std::make_unique<DeviceArray<galerkin::TriangleQuadraturePoint>>(nearRules.size());
            deviceNearRules->upload(nearRules, stream);
            deviceGauss = std::make_unique<DeviceArray<galerkin::GaussPoint>>(gauss.size());
            deviceGauss->upload(gauss, stream);
        }
        if (rebuildRegular) {
            deviceRegular = std::make_unique<DeviceArray<galerkin::TriangleQuadraturePoint>>(
                regularRule.size());
            deviceRegular->upload(regularRule, stream);
            deviceRegularBounds =
                std::make_unique<DeviceArray<galerkin::RegularPanelBounds>>(regularBounds.size());
            deviceRegularBounds->upload(regularBounds, stream);
        }
        stream.synchronize();
        cancelled(stop);
        cachedAcceleration = options.acceleration;
        cacheValid = true;
    }
    std::size_t cacheBytes() const {
        const auto bytes = [](const auto &array) { return array ? array->bytes() : 0; };
        return devicePanels.bytes() + deviceMass.bytes() + deviceRowScale.bytes() +
               bytes(deviceSingular) + bytes(deviceNear) + bytes(deviceNearOffsets) +
               bytes(deviceNearTrials) + bytes(deviceWholePanelDegrees) + bytes(devicePatches) +
               bytes(deviceRegular) + bytes(deviceRegularBounds) + bytes(deviceNearRules) +
               bytes(deviceGauss);
    }
    CudaStorageOptions storage;
    std::optional<BemSubdomainModel> domainModel;
    std::optional<OrthogonalSymmetryPlanes> symmetryModel;
    std::optional<ExteriorBaffle> baffleModel;
    // Declared before the registration so reverse destruction unregisters first.
    galerkin::InteractionPlan interactionPlan;
    Vec3 origin;
    SurfaceMesh mesh;
    galerkin::Space space;
    cuda_detail::BaffledPlanCache baffledPlanCache;
    std::vector<Panel> panels;
    Stream stream;
    Blas blas;
    DenseSolver solver;
    DeviceArray<Panel> devicePanels;
    DeviceArray<double> deviceMass, deviceRowScale;
    std::unique_ptr<DeviceArray<PairJob>> deviceSingular;
    std::unique_ptr<DeviceArray<NearJob>> deviceNear;
    std::unique_ptr<DeviceArray<unsigned>> deviceNearOffsets, deviceNearTrials;
    std::unique_ptr<DeviceArray<std::uint8_t>> deviceWholePanelDegrees;
    std::unique_ptr<QuadraturePatchArray> devicePatches;
    std::unique_ptr<DeviceArray<galerkin::TriangleQuadraturePoint>> deviceRegular;
    std::unique_ptr<DeviceArray<galerkin::RegularPanelBounds>> deviceRegularBounds;
    std::unique_ptr<DeviceArray<galerkin::TriangleQuadraturePoint>> deviceNearRules;
    std::unique_ptr<DeviceArray<galerkin::GaussPoint>> deviceGauss;
    BemNumerics cached;
    std::shared_ptr<const wideband::SupportTree> supportTree;
    std::shared_ptr<const wideband::SupportNearPlan> supportNearPlan;
    wideband::SupportTreeSettings supportTreeSettings;
    double cachedMaximumK{};
    BemAcceleration cachedAcceleration{BemAcceleration::dense};
    bool planValid{}, cacheValid{};
    std::unique_ptr<SolveWorkspace> workspace;
    bool observersValid{};
    std::vector<Vec3> cachedWorldObservers, cachedLocalObservers;
    std::optional<FarFieldObservation> cachedFarField;
    std::mutex mutex;
};

CudaBemSolver::CudaBemSolver(SurfaceMesh mesh)
    : CudaBemSolver(std::move(mesh), CudaStorageOptions{}) {}
CudaBemSolver::CudaBemSolver(SurfaceMesh mesh, CudaStorageOptions storage)
    : implementation_(std::make_unique<Implementation>(std::move(mesh), storage)) {}
CudaBemSolver::CudaBemSolver(SurfaceMesh mesh, BemSubdomainModel model, CudaStorageOptions storage)
    : implementation_(
          std::make_unique<Implementation>(std::move(mesh), storage, std::move(model))) {}
CudaBemSolver::CudaBemSolver(SurfaceMesh mesh, BemSubdomainModel model,
                             OrthogonalSymmetryPlanes symmetry, CudaStorageOptions storage)
    : implementation_(std::make_unique<Implementation>(std::move(mesh), storage, std::move(model),
                                                       std::move(symmetry))) {}
CudaBemSolver::CudaBemSolver(SurfaceMesh mesh, BemSubdomainModel model,
                             OrthogonalSymmetryPlanes symmetry, ExteriorBaffle baffle,
                             CudaStorageOptions storage)
    : implementation_(std::make_unique<Implementation>(std::move(mesh), storage, std::move(model),
                                                       std::move(symmetry), std::move(baffle))) {}
CudaBemSolver::CudaBemSolver(SurfaceMesh mesh, BemSubdomainModel model, ExteriorBaffle baffle,
                             CudaStorageOptions storage)
    : implementation_(std::make_unique<Implementation>(std::move(mesh), storage, std::move(model),
                                                       std::nullopt, std::move(baffle))) {}
CudaBemSolver::~CudaBemSolver() = default;
std::string CudaBemSolver::deviceDescription() {
    int count{};
    cudaCheck(cudaGetDeviceCount(&count), "CUDA device detection");
    if (!count)
        throw std::runtime_error("No CUDA device is available.");
    int selected{};
    cudaCheck(cudaGetDevice(&selected), "CUDA device selection");
    cudaDeviceProp properties{};
    cudaCheck(cudaGetDeviceProperties(&properties, selected), "CUDA device properties");
    return properties.name;
}
BemResult CudaBemSolver::solve(const BemRequest &request, std::stop_token stop,
                               std::function<void(const BemProgress &)> progress) {
    auto &state = *implementation_;
    std::unique_lock lock(state.mutex);
    cancelled(stop);
    if (request.symmetry != state.symmetryModel)
        throw std::invalid_argument("BEM request and solver symmetry snapshots differ.");
    if (request.exteriorBaffle != state.baffleModel)
        throw std::invalid_argument("BEM request and solver exterior baffle snapshots differ.");
    if (!request.observerValidity.empty() &&
        (!request.subdomains || !request.exteriorBaffle || !request.farField ||
         request.observerValidity.size() != request.observers.size()))
        throw std::invalid_argument(
            "Observer validity requires one entry per exterior-baffle far-field direction.");
    if (request.subdomains != state.domainModel)
        throw std::invalid_argument("BEM request and solver subdomain snapshots differ.");
    if (request.subdomains) {
        auto result = cuda_detail::solveSubdomains(state.mesh, state.origin, request, stop,
                                                   std::move(progress));
        result.peakDeviceBytes +=
            state.cacheBytes() + (state.workspace ? state.workspace->bytes() : 0);
        return result;
    }
    if (request.infiniteBaffle) {
        auto result = cuda_detail::solveBaffledCavity(state.mesh, state.origin, request,
                                                      state.baffledPlanCache, stop, progress);
        result.peakDeviceBytes +=
            state.cacheBytes() + (state.workspace ? state.workspace->bytes() : 0);
        return result;
    }
    if (request.axisymmetric)
        throw std::invalid_argument(
            "Axisymmetric BEM currently requires the recessed-air infinite-baffle domain.");
    const auto &options = request.numerics;
    const bool fastMultipole = options.acceleration == BemAcceleration::fastMultipole;
    const auto &fmm = options.fastMultipole;
    if ((options.acceleration != BemAcceleration::dense && !fastMultipole) ||
        (fastMultipole &&
         (fmm.leafSupportCapacity < 4 || fmm.leafSupportCapacity > 256 ||
          !std::isfinite(fmm.separationRatio) || fmm.separationRatio < 1.5 ||
          fmm.separationRatio > 10 || fmm.maximumIterations < 1 || fmm.maximumIterations > 10000 ||
          fmm.restartLength < 2 || fmm.restartLength > 300 ||
          (fmm.preconditionerBlockWidth != 32 && fmm.preconditionerBlockWidth != 64 &&
           fmm.preconditionerBlockWidth != 128 && fmm.preconditionerBlockWidth != 256 &&
           fmm.preconditionerBlockWidth != 512 && fmm.preconditionerBlockWidth != 1024 &&
           fmm.preconditionerBlockWidth != 2048))))
        throw std::invalid_argument("Invalid Fast Multipole or GMRES settings.");
    if (!std::isfinite(request.medium.density) || request.medium.density <= 0 ||
        !std::isfinite(request.medium.soundSpeed) || request.medium.soundSpeed <= 0 ||
        request.frequencies.empty() || (!request.coupling && request.excitations.empty()) ||
        (options.quadratureOrder != 2 && options.quadratureOrder != 4) ||
        options.singularOrder < 2 || options.singularOrder > 20 ||
        !std::isfinite(options.nearDistanceRatio) || options.nearDistanceRatio < 0 ||
        options.nearDistanceRatio > 10 || !std::isfinite(options.nearRelativeTolerance) ||
        options.nearRelativeTolerance < 1e-10 || options.nearRelativeTolerance > .01 ||
        options.nearMaximumDepth > 10 || !std::isfinite(options.relativeResidualTolerance) ||
        options.relativeResidualTolerance < 1e-7 || options.relativeResidualTolerance > .01 ||
        options.maximumRefinementSteps > 10)
        throw std::invalid_argument(
            "Invalid BEM medium, frequencies, excitation or Galerkin settings.");
    std::optional<bem_detail::CouplingBasis> couplingBasis;
    std::vector<double> couplingPanelAreas;
    if (request.coupling) {
        if (!request.excitations.empty())
            throw std::invalid_argument("LEM coupling replaces prescribed BEM excitations.");
        couplingPanelAreas.reserve(state.panels.size());
        for (const auto &panel : state.panels)
            couplingPanelAreas.push_back(panel.area);
        couplingBasis =
            bem_detail::makeCouplingBasis(*request.coupling, state.mesh, couplingPanelAreas);
    }
    const auto &excitations = couplingBasis ? couplingBasis->excitations : request.excitations;
    const double maximumK = bem_detail::planningWavenumber(
        request.frequencies, request.interactionPlanFrequencies, request.medium.soundSpeed);
    Vec3 farOrigin{};
    if (request.farField) {
        const auto &far = *request.farField;
        farOrigin = sub(far.origin, state.origin);
        if (!std::isfinite(far.origin.x) || !std::isfinite(far.origin.y) ||
            !std::isfinite(far.origin.z) || !std::isfinite(farOrigin.x) ||
            !std::isfinite(farOrigin.y) || !std::isfinite(farOrigin.z) ||
            !std::isfinite(far.referenceDistance) || far.referenceDistance <= 0 ||
            !std::isfinite(maximumK * far.referenceDistance))
            throw std::invalid_argument(
                "Far-field origin and positive reference distance must be finite.");
    }
    const auto samePoint = [](Vec3 a, Vec3 b) { return a.x == b.x && a.y == b.y && a.z == b.z; };
    if (!state.observersValid || state.cachedFarField != request.farField ||
        !std::equal(request.observers.begin(), request.observers.end(),
                    state.cachedWorldObservers.begin(), state.cachedWorldObservers.end(),
                    samePoint)) {
        std::vector<Vec3> local;
        local.reserve(request.observers.size());
        for (const auto worldPoint : request.observers) {
            auto point =
                sub(worldPoint, request.farField ? request.farField->origin : state.origin);
            if (!std::isfinite(worldPoint.x) || !std::isfinite(worldPoint.y) ||
                !std::isfinite(worldPoint.z) || !std::isfinite(point.x) ||
                !std::isfinite(point.y) || !std::isfinite(point.z))
                throw std::invalid_argument("BEM observation coordinates must be finite.");
            if (request.farField) {
                const double distance = std::hypot(point.x, point.y, point.z);
                if (!(distance > 0) || !std::isfinite(distance))
                    throw std::invalid_argument(
                        "A far-field direction must differ from its origin.");
                point = {point.x / distance, point.y / distance, point.z / distance};
            } else if (inside(state.panels, point))
                throw std::invalid_argument(
                    "BEM observation points must be finite and outside the solid.");
            local.push_back(point);
        }
        auto world = request.observers;
        state.cachedWorldObservers = std::move(world);
        state.cachedLocalObservers = std::move(local);
        state.cachedFarField = request.farField;
        state.observersValid = true;
    }
    const auto &localObservers = state.cachedLocalObservers;
    const auto panelCount = state.panels.size(), pressureCount = state.space.vertices.size(),
               excitationCount = excitations.size();
    if (pressureCount > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        panelCount > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        excitationCount > 1024 || request.observers.size() > 100000)
        throw std::invalid_argument("BEM request exceeds indexing limits.");
    const int np = static_cast<int>(pressureCount), nt = static_cast<int>(panelCount),
              nrhs = static_cast<int>(excitationCount);
    const auto workspaceSize = solveWorkspaceSize<Complex, DoubleComplex>(
        pressureCount, panelCount, excitationCount, request.observers.size(), fastMultipole);
    const bool needsWorkspace =
        !state.workspace ||
        !state.workspace->matches(excitationCount, request.observers.size(), fastMultipole);
    const auto mib = [](std::size_t bytes, bool roundUp = false) {
        constexpr std::size_t unit = 1024 * 1024;
        return std::to_string(bytes / unit + (roundUp && bytes % unit != 0));
    };
    if (needsWorkspace) {
        std::size_t free{}, total{};
        cudaCheck(cudaMemGetInfo(&free, &total), "CUDA capacity query");
        // Total capacity is independent of caches that can be released below.
        // Reject only a provably impossible request before any near-pair planning.
        if (workspaceSize.minimumBytes > total)
            throw std::runtime_error("BEM: mindestens " + mib(workspaceSize.minimumBytes, true) +
                                     " MiB + LU; GPU-Gesamtkapazität " + mib(total) + " MiB.");
    }
    std::map<std::uint32_t, BoundaryCondition> conditions;
    for (const auto &triangle : state.mesh.triangles) {
        const auto [it, inserted] = conditions.emplace(triangle.faceId, triangle.boundary);
        if (!inserted && it->second != triangle.boundary)
            throw std::invalid_argument("Conflicting boundary conditions on one face.");
    }
    std::vector<std::vector<std::complex<double>>> velocities;
    std::map<std::uint32_t, std::size_t> activeFaceColumns;
    for (const auto &excitation : excitations) {
        std::map<std::uint32_t, std::complex<double>> values;
        for (const auto &face : excitation.faces) {
            if (!conditions.contains(face.faceId) ||
                conditions.at(face.faceId) != BoundaryCondition::driven ||
                !std::isfinite(face.normalVelocity.real()) ||
                !std::isfinite(face.normalVelocity.imag()) ||
                !values.emplace(face.faceId, face.normalVelocity).second)
                throw std::invalid_argument(
                    "Excitation must assign each known Driven face at most once.");
            // Preserve all requested ports, including a currently zero flow.
            activeFaceColumns.emplace(face.faceId, 0);
        }
        if (values.empty())
            throw std::invalid_argument("Each excitation needs at least one assigned Driven face.");
        auto &perPanel = velocities.emplace_back();
        perPanel.reserve(panelCount);
        for (const auto &triangle : state.mesh.triangles) {
            const auto found = values.find(triangle.faceId);
            perPanel.push_back(found == values.end() ? std::complex<double>{} : found->second);
        }
    }
    if (state.workspace && needsWorkspace) {
        state.stream.synchronize();
        state.workspace.reset();
    }
    std::vector<std::uint32_t> activeFaceIds;
    for (auto &[faceId, column] : activeFaceColumns) {
        column = activeFaceIds.size();
        activeFaceIds.push_back(faceId);
    }
    std::vector<std::vector<std::complex<double>>> faceVelocities(
        excitationCount, std::vector<std::complex<double>>(activeFaceIds.size()));
    for (std::size_t e = 0; e < excitationCount; ++e)
        for (const auto &face : excitations[e].faces)
            faceVelocities[e][activeFaceColumns.at(face.faceId)] = face.normalVelocity;
    const bool hasNonzeroExcitation =
        std::any_of(faceVelocities.begin(), faceVelocities.end(), [](const auto &values) {
            return std::any_of(values.begin(), values.end(), [](auto v) { return v != 0.; });
        });
    bool hasActiveFar = false;
    if (fastMultipole) {
        if (!state.supportTree ||
            state.supportTreeSettings.leafCapacity != fmm.leafSupportCapacity ||
            state.supportTreeSettings.separationRatio != fmm.separationRatio ||
            state.supportTreeSettings.nearDistanceRatio != options.nearDistanceRatio) {
            wideband::SupportTreeSettings settings;
            settings.leafCapacity = fmm.leafSupportCapacity;
            settings.separationRatio = fmm.separationRatio;
            settings.nearDistanceRatio = options.nearDistanceRatio;
            state.supportTree = std::make_shared<wideband::SupportTree>(
                wideband::makeSupportTree(state.space, settings, stop));
            state.supportTreeSettings = settings;
            state.supportNearPlan.reset();
        }
        if (!state.supportNearPlan || state.supportNearPlan->activeFaceIds != activeFaceIds) {
            // Explicitly bound host CSR indices and construction scratch; the
            // complete support tree is independent of frequency and RHS values.
            constexpr std::size_t maximumNearHostBytes = std::size_t{2} * 1024 * 1024 * 1024;
            state.supportNearPlan =
                std::make_shared<wideband::SupportNearPlan>(wideband::makeSupportNearPlan(
                    state.space, *state.supportTree, activeFaceIds, maximumNearHostBytes, stop));
        }
        hasActiveFar = !state.supportTree->pressurePressure.farInteractions.empty();
        if (!hasActiveFar) {
            std::vector<std::size_t> activePrefix(panelCount + 1);
            for (std::size_t i = 0; i < panelCount; ++i)
                activePrefix[i + 1] =
                    activePrefix[i] +
                    (state.supportNearPlan
                         ->panelToActiveFaceColumn[state.supportTree->fluxOrder[i]] != UINT32_MAX);
            for (const auto pair : state.supportTree->pressureFlux.farInteractions) {
                const auto &node = state.supportTree->nodes[pair.source];
                hasActiveFar |=
                    activePrefix[node.firstFlux + node.fluxCount] > activePrefix[node.firstFlux];
            }
        }
    } else {
        state.supportNearPlan.reset();
        state.supportTree.reset();
    }
    state.prepareQuadrature(options, maximumK, stop,
                            needsWorkspace ? workspaceSize.minimumBytes : 0);
    cancelled(stop);
    const std::size_t rhsElements = workspaceSize.rhsElements;
    if (!state.workspace) {
        std::size_t free{}, total{};
        cudaCheck(cudaMemGetInfo(&free, &total), "CUDA memory query");
        if (workspaceSize.minimumBytes > deviceMemoryBudget(free))
            throw std::runtime_error("BEM: mindestens " + mib(workspaceSize.minimumBytes, true) +
                                     " MiB + LU; verfügbar " + mib(free) + " MiB, Budget " +
                                     mib(deviceMemoryBudget(free)) + " MiB.");
        state.workspace =
            std::make_unique<SolveWorkspace>(np, nrhs, static_cast<int>(request.observers.size()),
                                             workspaceSize, state.solver, fastMultipole);
    }
    auto &buffers = *state.workspace;
    auto &matrix = buffers.matrix;
    auto &factor = buffers.factor;
    auto &rhs = buffers.rhs;
    auto &pressure = buffers.pressure;
    auto &correction = buffers.correction;
    auto &originalResidual = buffers.originalResidual;
    auto &flux = buffers.flux;
    auto &fields = buffers.fields;
    auto &observers = buffers.observers;
    auto &pivots = buffers.pivots;
    auto &info = buffers.info;
    auto &integrationError = buffers.integrationError;
    DeviceArray<AccurateComplex> faceFlux(fastMultipole ? activeFaceIds.size() * excitationCount
                                                        : 0);
    std::vector<AccurateComplex> faceSource(faceFlux.size());
    observers.upload(localObservers, state.stream);
    BemResult result;
    result.acceleration = options.acceleration;
    result.device = deviceDescription();
    result.triangleCount = panelCount;
    result.pressureDofCount = pressureCount;
    result.peakDeviceBytes = state.cacheBytes() + buffers.bytes() + faceFlux.bytes();
    for (const auto frequency : request.frequencies) {
        cancelled(stop);
        BemProgress frequencyProgress{result.frequencies.size(), request.frequencies.size(),
                                      frequency};
        frequencyProgress.rightHandSideCount = excitationCount + (couplingBasis ? 1 : 0);
        const auto report = [&](BemSolvePhase phase) {
            frequencyProgress.phase = phase;
            if (progress)
                progress(frequencyProgress);
        };
        report(BemSolvePhase::preparing);
        BemFrequencyResult frame;
        frame.frequency = frequency;
        const double omega = 2 * std::numbers::pi * frequency,
                     k = omega / request.medium.soundSpeed;
        std::vector<Complex> source(panelCount * excitationCount);
        for (std::size_t e = 0; e < excitationCount; ++e)
            for (std::size_t panel = 0; panel < panelCount; ++panel) {
                const auto q =
                    std::complex<double>(0, -request.medium.density * omega) * velocities[e][panel];
                if (!std::isfinite(q.real()) || !std::isfinite(q.imag()) ||
                    std::abs(q) > std::numeric_limits<float>::max())
                    throw std::invalid_argument("BEM flux exceeds FP32 storage.");
                source[e * panelCount + panel] =
                    make_cuFloatComplex(static_cast<float>(q.real()), static_cast<float>(q.imag()));
            }
        const auto assemblyStart = Clock::now();
        double maximumKh{};
        for (const auto &panel : state.space.panels)
            maximumKh = std::max(maximumKh, k * panel.maximumEdgeLength);
        const unsigned endpointDegree =
            fastMultipole
                ? (hasActiveFar && hasNonzeroExcitation
                       ? galerkin::minimumEndpointDegree(maximumKh, options.quadratureOrder)
                       : 4)
                : options.quadratureOrder;
        std::unique_ptr<wideband::SupportFmmOperator> fastOperator;
        std::unique_ptr<wideband::SupportFmmLinearSolver> fastSolver;
        flux.upload(source, state.stream);
        matrix.clear(state.stream);
        rhs.clear(state.stream);
        integrationError.clear(state.stream);
        if (fastMultipole) {
            for (std::size_t e = 0; e < excitationCount; ++e)
                for (std::size_t face = 0; face < activeFaceIds.size(); ++face) {
                    const auto q = std::complex<double>(0, -request.medium.density * omega) *
                                   faceVelocities[e][face];
                    faceSource[e * activeFaceIds.size() + face] =
                        make_cuDoubleComplex(q.real(), q.imag());
                }
            faceFlux.upload(faceSource, state.stream);
            wideband::SupportNearOperatorSettings nearSettings;
            nearSettings.quadrature.regularOrder = options.quadratureOrder;
            nearSettings.quadrature.singularOrder = options.singularOrder;
            nearSettings.quadrature.nearDistanceRatio = options.nearDistanceRatio;
            nearSettings.quadrature.relativeTolerance = options.nearRelativeTolerance;
            nearSettings.quadrature.maximumSubdivisionDepth = options.nearMaximumDepth;
            const auto persistent = state.cacheBytes() + buffers.bytes() + faceFlux.bytes();
            report(BemSolvePhase::nearConstruction);
            auto near = std::make_shared<wideband::SupportNearOperator>(
                state.space, *state.supportTree, *state.supportNearPlan, k, nearSettings, stop,
                state.stream);
            result.peakDeviceBytes =
                std::max(result.peakDeviceBytes, persistent + near->peakDeviceBytes());
            auto pc = std::make_shared<wideband::SupportBlockJacobi>(
                near->pressureMatrix(), state.supportTree->pressureOrder,
                fmm.preconditionerBlockWidth, stop, state.stream);
            result.peakDeviceBytes = std::max(
                result.peakDeviceBytes, persistent + near->deviceBytes() + pc->peakDeviceBytes());
            wideband::SupportFmmOperatorSettings settings;
            settings.role = wideband::SupportFmmOperatorRole::distancePrimary;
            settings.accuracySeed = .1 * options.relativeResidualTolerance;
            settings.panelRule = galerkin::triangleRule(endpointDegree);
            report(BemSolvePhase::farPreparation);
            fastOperator = std::make_unique<wideband::SupportFmmOperator>(
                state.space, state.supportTree, state.supportNearPlan, near, pc, k,
                std::move(settings), stop);
            result.peakDeviceBytes =
                std::max(result.peakDeviceBytes, persistent + fastOperator->peakDeviceBytes());
            wideband::SupportDirectAuditPolicy directPolicy;
            if (panelCount <= wideband::maximumDirectAuditPanels) {
                int deviceId{};
                cudaCheck(cudaGetDevice(&deviceId), "Direct audit cost device");
                int throughputRatio{};
                cudaCheck(cudaDeviceGetAttribute(&throughputRatio,
                                                 cudaDevAttrSingleToDoublePrecisionPerfRatio,
                                                 deviceId),
                          "Direct audit cost throughput");
                const auto ratio = static_cast<unsigned>(std::max(1, throughputRatio));
                directPolicy.preferDirect = [tree = state.supportTree, ratio](auto levels,
                                                                              unsigned degree) {
                    const auto points =
                        static_cast<unsigned>(galerkin::triangleRule(degree).size());
                    return wideband::preferSupportDirectAudit(*tree, levels, points, ratio);
                };
                directPolicy.makeDirect = [&state, near, k, stop](unsigned order) {
                    return std::make_unique<wideband::SupportDirectAudit>(
                        state.space, *state.supportTree, *state.supportNearPlan, near, k, order,
                        stop);
                };
            }
            fastSolver = std::make_unique<wideband::SupportFmmLinearSolver>(
                *fastOperator, pressureCount, options, stop, std::move(directPolicy),
                [&](BemSolvePhase phase, unsigned count, std::size_t pressureCount,
                    std::size_t fluxCount) {
                    frequencyProgress.iterations = count;
                    frequencyProgress.pressureOperatorApplications = pressureCount;
                    frequencyProgress.faceFluxOperatorApplications = fluxCount;
                    report(phase);
                },
                endpointDegree);
            frame.fastMultipole = fastOperator->statistics();
            result.peakDeviceBytes =
                std::max(result.peakDeviceBytes,
                         state.cacheBytes() + buffers.bytes() + faceFlux.bytes() +
                             fastOperator->deviceBytes() + fastSolver->deviceBytes());
        } else {
            constexpr std::size_t batchPairs = 4194304, batchJobs = 4096;
            const auto batchTrials =
                std::min<std::size_t>(65535, std::max<std::size_t>(1, batchPairs / panelCount));
            for (std::size_t first = 0; first < panelCount; first += batchTrials) {
                const auto count = std::min(batchTrials, panelCount - first);
                const dim3 grid(static_cast<unsigned>((panelCount + 127) / 128),
                                static_cast<unsigned>(count));
                assembleRegular<<<grid, 128, 0, state.stream>>>(
                    state.devicePanels.data(), nt, np, state.deviceRegular->data(),
                    state.deviceRegularBounds->data(), state.deviceRowScale.data(), flux.data(),
                    nrhs, k, state.deviceNearOffsets->data(), state.deviceNearTrials->data(),
                    state.deviceWholePanelDegrees->data(), state.deviceNearRules->data(),
                    matrix.data(), rhs.data(), integrationError.data(),
                    static_cast<unsigned>(first), options.nearRelativeTolerance,
                    options.quadratureOrder);
                cudaCheck(cudaGetLastError(), "CUDA regular Galerkin assembly");
                state.stream.synchronize();
                cancelled(stop);
            }
            for (std::size_t first = 0; first < state.deviceSingular->size(); first += batchJobs) {
                const auto count = std::min(batchJobs, state.deviceSingular->size() - first);
                const auto launch = [&]<unsigned order>() {
                    assembleSingular<order><<<static_cast<unsigned>(count), 128, 0, state.stream>>>(
                        state.devicePanels.data(), state.deviceSingular->data(),
                        state.deviceGauss->data(), options.singularOrder, nt, np,
                        state.deviceRowScale.data(), flux.data(), nrhs, k, matrix.data(),
                        rhs.data(), integrationError.data(), first);
                };
                if (options.singularOrder == 4)
                    launch.template operator()<4>();
                else
                    launch.template operator()<0>();
                cudaCheck(cudaGetLastError(), "CUDA singular Galerkin assembly");
                state.stream.synchronize();
                cancelled(stop);
            }
            for (std::size_t first = 0; first < state.deviceNear->size(); first += batchJobs) {
                const auto count = std::min(batchJobs, state.deviceNear->size() - first);
                const auto baseCount = galerkin::triangleRule(8).size();
                const auto promotedCount = state.deviceNearRules->size() - baseCount;
                const auto launch = [&]<bool mixed>() {
                    assembleNear<mixed><<<static_cast<unsigned>(count), 128,
                                          2 * promotedCount * sizeof(NearSample), state.stream>>>(
                        state.devicePanels.data(), state.deviceNear->data(),
                        state.devicePatches->data(), state.deviceNearRules->data(),
                        static_cast<int>(baseCount), static_cast<int>(promotedCount), nt, np,
                        state.deviceRowScale.data(), flux.data(), nrhs, k, matrix.data(),
                        rhs.data(), integrationError.data(), first);
                };
                if (options.nearRelativeTolerance >= 1e-5)
                    launch.template operator()<true>();
                else
                    launch.template operator()<false>();
                cudaCheck(cudaGetLastError(), "CUDA near Galerkin assembly");
                state.stream.synchronize();
                cancelled(stop);
            }
            addMass<<<static_cast<unsigned>((panelCount + 127) / 128), 128, 0, state.stream>>>(
                state.devicePanels.data(), nt, np, state.deviceRowScale.data(), flux.data(), nrhs,
                k, matrix.data(), rhs.data(), integrationError.data());
            cudaCheck(cudaGetLastError(), "CUDA Galerkin identity assembly");
            checkIntegrationError(integrationError, state.stream);
        }
        auto rhsHost = rhs.download(state.stream);
        std::vector<double> rhsNorm(excitationCount, 0);
        for (std::size_t e = 0; e < excitationCount; ++e)
            for (std::size_t row = 0; row < pressureCount; ++row)
                rhsNorm[e] =
                    std::hypot(rhsNorm[e], std::abs(hostComplex(rhsHost[e * pressureCount + row])) *
                                               state.space.lumpedMass[row]);
        for (const auto value : rhsNorm)
            if (!std::isfinite(value))
                throw std::runtime_error("Galerkin right-hand side is not finite.");
        frame.assemblyMilliseconds = milliseconds(assemblyStart);
        const auto solveStart = Clock::now();
        auto triangularStart = Clock::now();
        if (!fastMultipole) {
            cudaCheck(cudaMemcpyAsync(factor.data(), matrix.data(), matrix.bytes(),
                                      cudaMemcpyDeviceToDevice, state.stream),
                      "CUDA factorization copy");
            state.stream.synchronize();
            const auto factorStart = Clock::now();
            solverCheck(cusolverDnCgetrf(state.solver, np, np, factor.data(), np,
                                         buffers.lu->data(), pivots.data(), info.data()),
                        "CUDA pivoted LU");
            const int factorInfo = info.download(state.stream).front();
            if (factorInfo != 0)
                throw std::runtime_error("CUDA Burton-Miller factorization failed at pivot " +
                                         std::to_string(factorInfo));
            frame.factorizationMilliseconds = milliseconds(factorStart);
            cudaCheck(cudaMemcpyAsync(correction.data(), rhs.data(), rhs.bytes(),
                                      cudaMemcpyDeviceToDevice, state.stream),
                      "CUDA excitation RHS copy");
            state.stream.synchronize();
            triangularStart = Clock::now();
            solverCheck(cusolverDnCgetrs(state.solver, CUBLAS_OP_N, np, nrhs, factor.data(), np,
                                         pivots.data(), correction.data(), np, info.data()),
                        "CUDA multi-RHS solve");
            if (info.download(state.stream).front() != 0)
                throw std::runtime_error("CUDA triangular solve failed.");
            initializePressure<<<static_cast<unsigned>((rhsElements + 255) / 256), 256, 0,
                                 state.stream>>>(correction.data(), pressure.data(), rhsElements);
            cudaCheck(cudaGetLastError(), "CUDA accurate solution initialization");
            frame.triangularSolveMilliseconds = milliseconds(triangularStart);
        }
        std::vector<double> relative(excitationCount);
        std::vector<wideband::SupportFmmSolution> fastEvidence(excitationCount);
        std::vector<unsigned> iterations(excitationCount), refinementCounts(excitationCount);
        unsigned refinement = 0;
        bool fastInitialGuess = false;
        const auto refine = [&](int activeCount) {
            if (fastMultipole) {
                for (int e = 0; e < activeCount; ++e) {
                    frequencyProgress.rightHandSideIndex =
                        fastInitialGuess ? excitationCount : static_cast<std::size_t>(e);
                    auto solved = fastSolver->solve(
                        faceFlux.data() + static_cast<std::size_t>(e) * activeFaceIds.size(),
                        pressure.data() + static_cast<std::size_t>(e) * pressureCount,
                        fastInitialGuess);
                    relative[e] = solved.boundaryResidual;
                    iterations[e] += solved.iterations;
                    refinementCounts[e] += solved.refinementSteps;
                    frame.residualCheckMilliseconds += solved.residualCheckMilliseconds;
                    fastEvidence[e] = std::move(solved);
                }
                return;
            }
            const auto activeElements = pressureCount * static_cast<std::size_t>(activeCount);
            for (;;) {
                cancelled(stop);
                trueResidual<<<dim3(static_cast<unsigned>((pressureCount + 31) / 32), activeCount),
                               dim3(32, 8), 0, state.stream>>>(
                    matrix.data(), rhs.data(), pressure.data(), state.deviceMass.data(),
                    correction.data(), originalResidual.data(), np);
                cudaCheck(cudaGetLastError(), "CUDA FP64 physical residual");
                state.stream.synchronize();
                double largest = 0;
                for (int e = 0; e < activeCount; ++e) {
                    const auto residual =
                        vectorNorm(state.blas, np, originalResidual.data() + e * pressureCount);
                    relative[e] =
                        rhsNorm[e] > 0
                            ? residual / rhsNorm[e]
                            : (residual == 0 ? 0 : std::numeric_limits<double>::infinity());
                    if (!std::isfinite(relative[e]))
                        throw std::runtime_error("Burton-Miller residual is not finite.");
                    largest = std::max(largest, relative[e]);
                }
                if (largest <= options.relativeResidualTolerance)
                    break;
                if (refinement >= options.maximumRefinementSteps) {
                    std::ostringstream message;
                    message << std::setprecision(9) << "Burton-Miller residual " << largest
                            << " exceeds tolerance " << options.relativeResidualTolerance << " at "
                            << frequency << " Hz after " << refinement << " refinement steps.";
                    throw std::runtime_error(message.str());
                }
                state.stream.synchronize();
                triangularStart = Clock::now();
                solverCheck(cusolverDnCgetrs(state.solver, CUBLAS_OP_N, np, activeCount,
                                             factor.data(), np, pivots.data(), correction.data(),
                                             np, info.data()),
                            "CUDA refinement solve");
                if (info.download(state.stream).front() != 0)
                    throw std::runtime_error("CUDA refinement triangular solve failed.");
                frame.triangularSolveMilliseconds += milliseconds(triangularStart);
                updatePressure<<<static_cast<unsigned>((activeElements + 255) / 256), 256, 0,
                                 state.stream>>>(correction.data(), pressure.data(),
                                                 activeElements);
                cudaCheck(cudaGetLastError(), "CUDA accurate refinement update");
                ++refinement;
            }
        };
        refine(nrhs);

        std::optional<lem::Solution> lemSolution;
        auto actualExcitations = excitations;
        if (couplingBasis) {
            cancelled(stop);
            auto nodalBasis = pressure.download(state.stream);
            lem::AcousticLoad load;
            load.impedance.resize(excitationCount * excitationCount);
            // The unit-Q velocity shape is dual to modal pressure: p_modal * conj(Q)
            // is the boundary's complex acoustic power. Include every mutual term.
            for (std::size_t j = 0; j < excitationCount; ++j)
                for (std::size_t panel = 0; panel < panelCount; ++panel) {
                    const int port = couplingBasis->panelPorts[panel];
                    if (port < 0)
                        continue;
                    std::complex<double> mean{};
                    for (const auto node : state.panels[panel].nodes)
                        mean += hostComplex(nodalBasis[j * pressureCount + node]) / 3.0;
                    load.impedance[static_cast<std::size_t>(port) * excitationCount + j] +=
                        mean * (state.panels[panel].area * couplingBasis->panelVelocities[panel]);
                }
            lemSolution = lem::solve(*request.coupling, frequency, load,
                                     {request.medium.density, request.medium.soundSpeed}, stop);
            actualExcitations = {
                {"LEM", bem_detail::coupledVelocities(*couplingBasis, *lemSolution)}};

            // Reconstruct the coupled field and Bq in the same discrete basis.
            // The retained LU factors also refine this combined equation; its
            // residual is measured directly, never inferred from basis residuals.
            std::vector<AccurateComplex> actualPressure(rhsElements);
            std::vector<Complex> actualRhs(rhsElements);
            rhsNorm[0] = 0;
            for (std::size_t row = 0; row < pressureCount; ++row) {
                std::complex<double> p{}, b{};
                for (std::size_t j = 0; j < excitationCount; ++j) {
                    const auto flow = lemSolution->ports[j].volumeVelocity;
                    p += hostComplex(nodalBasis[j * pressureCount + row]) * flow;
                    b += hostComplex(rhsHost[j * pressureCount + row]) * flow;
                }
                if (!std::isfinite(p.real()) || !std::isfinite(p.imag()) ||
                    !std::isfinite(b.real()) || !std::isfinite(b.imag()) ||
                    std::abs(b) > std::numeric_limits<float>::max())
                    throw std::runtime_error("Coupled BEM pressure or RHS exceeds storage.");
                actualPressure[row] = make_cuDoubleComplex(p.real(), p.imag());
                actualRhs[row] =
                    make_cuFloatComplex(static_cast<float>(b.real()), static_cast<float>(b.imag()));
                rhsNorm[0] = std::hypot(rhsNorm[0], std::abs(hostComplex(actualRhs[row])) *
                                                        state.space.lumpedMass[row]);
            }
            std::fill(source.begin(), source.end(), make_cuFloatComplex(0, 0));
            for (std::size_t panel = 0; panel < panelCount; ++panel) {
                const int port = couplingBasis->panelPorts[panel];
                if (port < 0)
                    continue;
                const auto q = std::complex<double>(0, -request.medium.density * omega) *
                               lemSolution->ports[port].volumeVelocity *
                               couplingBasis->panelVelocities[panel];
                if (!std::isfinite(q.real()) || !std::isfinite(q.imag()) ||
                    std::abs(q) > std::numeric_limits<float>::max())
                    throw std::runtime_error("Coupled BEM flux exceeds FP32 storage.");
                source[panel] =
                    make_cuFloatComplex(static_cast<float>(q.real()), static_cast<float>(q.imag()));
            }
            pressure.upload(actualPressure, state.stream);
            rhs.upload(actualRhs, state.stream);
            flux.upload(source, state.stream);
            if (fastMultipole) {
                std::fill(faceSource.begin(), faceSource.end(), make_cuDoubleComplex(0, 0));
                for (const auto &face : actualExcitations.front().faces) {
                    const auto q = std::complex<double>(0, -request.medium.density * omega) *
                                   face.normalVelocity;
                    faceSource[activeFaceColumns.at(face.faceId)] =
                        make_cuDoubleComplex(q.real(), q.imag());
                }
                faceFlux.upload(faceSource, state.stream);
                // The single coupled output represents every unit-port solve.
                // Count all basis work and subsequent actual-field corrections.
                for (std::size_t e = 1; e < excitationCount; ++e) {
                    iterations[0] += iterations[e];
                    refinementCounts[0] += refinementCounts[e];
                }
            }
            fastInitialGuess = true;
            refine(1);
        }
        state.stream.synchronize();
        frame.solveMilliseconds = milliseconds(solveStart);
        if (fastMultipole) {
            frame.fastMultipole.graphExecution = fastSolver->graphExecutionEvidence();
            const auto persistent =
                state.cacheBytes() + buffers.bytes() + faceFlux.bytes() + fastSolver->deviceBytes();
            result.liveDeviceBytes = std::max(
                result.liveDeviceBytes, persistent + fastSolver->stationaryOperatorDeviceBytes());
            result.peakDeviceBytes = std::max(result.peakDeviceBytes,
                                              persistent + fastSolver->peakOperatorDeviceBytes());
        } else {
            result.liveDeviceBytes =
                std::max(result.liveDeviceBytes, state.cacheBytes() + buffers.bytes());
        }
        const auto outputCount = actualExcitations.size();
        const int fieldRhs = static_cast<int>(outputCount);
        for (std::size_t e = 0; e < outputCount; ++e) {
            BemExcitationResult output;
            output.name = actualExcitations[e].name;
            output.refinementSteps = fastMultipole ? refinementCounts[e] : refinement;
            output.iterations = iterations[e];
            output.boundaryResidual = relative[e];
            if (fastMultipole) {
                const auto &evidence = fastEvidence[e];
                output.residualEvaluation = evidence.residualEvaluation;
                output.fastMultipoleKrylovResidual = evidence.krylovResidual;
                output.fastMultipoleKrylovTarget = evidence.krylovTarget;
                output.fastMultipoleOperatorDiscrepancy = evidence.operatorDiscrepancy;
                output.fastMultipoleAccuracySeed = evidence.accuracySeed;
                output.fastMultipoleComparisonAccuracySeed = evidence.comparisonAccuracySeed;
                output.fastMultipoleVerificationAccuracySeed = evidence.verificationAccuracySeed;
                output.fastMultipoleAccuracyRefinements = evidence.accuracyRefinements;
                output.fastMultipoleComparisonLevels = evidence.comparisonLevels;
                output.fastMultipoleVerificationLevels = evidence.verificationLevels;
                output.farQuadratureComparisonDegree = evidence.farQuadratureComparisonDegree;
                output.farQuadratureVerificationDegree = evidence.farQuadratureVerificationDegree;
                output.farQuadratureRefinements = evidence.farQuadratureRefinements;
                output.farQuadratureDiscrepancy = evidence.farQuadratureDiscrepancy;
                output.farQuadraturePreparationMilliseconds =
                    evidence.farQuadraturePreparationMilliseconds;
                output.farQuadraturePressureMilliseconds =
                    evidence.farQuadraturePressureMilliseconds;
                output.farQuadratureFluxMilliseconds = evidence.farQuadratureFluxMilliseconds;
            }
            output.solvedVelocities = actualExcitations[e].faces;
            if (lemSolution)
                output.lemSolution = std::move(lemSolution);
            frame.excitations.push_back(std::move(output));
        }
        report(BemSolvePhase::fieldEvaluation);
        const auto fieldStart = Clock::now();
        std::vector<Complex> observed;
        if (!request.observers.empty() && !hasNonzeroExcitation) {
            observed.resize(outputCount * request.observers.size());
        } else if (!request.observers.empty()) {
            const auto evaluate = [&](unsigned degree) {
                cancelled(stop);
                const auto rule = galerkin::triangleRule(degree);
                DeviceArray<galerkin::TriangleQuadraturePoint> deviceRule(rule.size());
                deviceRule.upload(rule, state.stream);
                result.peakDeviceBytes = std::max(
                    result.peakDeviceBytes, state.cacheBytes() + buffers.bytes() +
                                                faceFlux.bytes() + deviceRule.bytes() +
                                                (fastOperator ? fastOperator->deviceBytes() : 0) +
                                                (fastSolver ? fastSolver->deviceBytes() : 0));
                if (request.farField)
                    evaluateFarField<<<static_cast<unsigned>(outputCount *
                                                             request.observers.size() *
                                                             buffers.tileCount),
                                       128, 0, state.stream>>>(
                        state.devicePanels.data(), nt, observers.data(),
                        static_cast<int>(request.observers.size()), farOrigin,
                        request.farField->referenceDistance, deviceRule.data(),
                        static_cast<int>(rule.size()), k, pressure.data(), np, flux.data(),
                        buffers.fieldPartials.data(), fieldRhs, buffers.tileCount,
                        integrationError.data());
                else
                    evaluateField<<<static_cast<unsigned>(outputCount * request.observers.size() *
                                                          buffers.tileCount),
                                    128, 0, state.stream>>>(
                        state.devicePanels.data(), nt, observers.data(),
                        static_cast<int>(request.observers.size()), deviceRule.data(),
                        static_cast<int>(rule.size()), k, pressure.data(), np, flux.data(),
                        buffers.fieldPartials.data(), fieldRhs, buffers.tileCount,
                        std::min(options.nearRelativeTolerance,
                                 .025 * options.relativeResidualTolerance),
                        options.nearMaximumDepth, integrationError.data());
                cudaCheck(cudaGetLastError(), "CUDA P1 field evaluation");
                finishField<<<static_cast<unsigned>((outputCount * request.observers.size() + 127) /
                                                    128),
                              128, 0, state.stream>>>(buffers.fieldPartials.data(), fields.data(),
                                                      outputCount * request.observers.size(),
                                                      buffers.tileCount, integrationError.data());
                cudaCheck(cudaGetLastError(), "CUDA P1 field tile reduction");
                checkIntegrationError(integrationError, state.stream);
                return fields.download(state.stream);
            };
            unsigned degree = galerkin::minimumEndpointDegree(maximumKh, options.quadratureOrder);
            observed = evaluate(degree);
            for (;;) {
                const auto next = galerkin::refinedEndpointDegree(degree);
                auto refined = evaluate(next);
                bool passed = true;
                for (std::size_t e = 0; e < outputCount; ++e) {
                    double differenceSquared{}, normSquared{};
                    for (std::size_t point = 0; point < request.observers.size(); ++point) {
                        const auto index = e * request.observers.size() + point;
                        const auto higher = hostComplex(refined[index]);
                        differenceSquared += std::norm(higher - hostComplex(observed[index]));
                        normSquared += std::norm(higher);
                    }
                    const double indicator =
                        normSquared > 0          ? std::sqrt(differenceSquared / normSquared)
                        : differenceSquared == 0 ? 0
                                                 : std::numeric_limits<double>::infinity();
                    auto &output = frame.excitations[e];
                    output.fieldQuadratureComparisonDegree = degree;
                    output.fieldQuadratureVerificationDegree = next;
                    output.fieldQuadratureDiscrepancy = indicator;
                    passed &= indicator <= .025 * options.relativeResidualTolerance;
                }
                observed = std::move(refined);
                if (passed)
                    break;
                if (degree == 12)
                    throw std::runtime_error(
                        "BEM complex observer-field quadrature 12/14 did not converge within "
                        "0.025*residualTolerance; refine the mesh or observer separation. "
                        "Equation residuals do not bound field quadrature.");
                degree = next;
            }
        }
        const auto nodal = pressure.download(state.stream);
        frame.fieldMilliseconds = milliseconds(fieldStart);
        for (std::size_t e = 0; e < outputCount; ++e) {
            auto &output = frame.excitations[e];
            BemBoundaryData boundary;
            boundary.surfacePressure.reserve(panelCount);
            boundary.normalDerivative.reserve(panelCount);
            boundary.pressureTrace.reserve(pressureCount);
            std::map<std::uint32_t, std::complex<double>> actualVelocities;
            for (const auto &face : actualExcitations[e].faces)
                actualVelocities.emplace(face.faceId, face.normalVelocity);
            for (std::size_t node = 0; node < pressureCount; ++node)
                boundary.pressureTrace.push_back(hostComplex(nodal[e * pressureCount + node]));
            for (const auto &panel : state.space.panels) {
                std::complex<double> average{};
                for (const auto node : panel.nodes)
                    average += hostComplex(nodal[e * pressureCount + node]) / 3.0;
                if (!std::isfinite(average.real()) || !std::isfinite(average.imag()))
                    throw std::runtime_error("Burton-Miller pressure is not finite.");
                boundary.surfacePressure.push_back(average);
                const auto found = actualVelocities.find(panel.faceId);
                boundary.normalDerivative.push_back(
                    found == actualVelocities.end()
                        ? std::complex<double>{}
                        : std::complex<double>{0, -request.medium.density * omega} * found->second);
            }
            if (couplingBasis)
                bem_detail::validateCouplingClosure(*couplingBasis, *output.lemSolution,
                                                    couplingPanelAreas, boundary.surfacePressure);
            for (std::size_t point = 0; point < request.observers.size(); ++point)
                output.observerPressure.push_back(
                    hostComplex(observed[e * request.observers.size() + point]));
            for (const auto face : actualExcitations[e].faces) {
                if (std::abs(face.normalVelocity) == 0)
                    continue;
                double area = 0;
                std::complex<double> integrated{};
                for (std::size_t panel = 0; panel < panelCount; ++panel)
                    if (state.mesh.triangles[panel].faceId == face.faceId) {
                        area += state.panels[panel].area;
                        integrated += boundary.surfacePressure[panel] * state.panels[panel].area;
                    }
                const auto mean = integrated / area, volume = face.normalVelocity * area;
                output.impedances.push_back({face.faceId, area, mean, volume, mean / volume});
            }
            output.setBoundary(std::move(boundary));
        }
        cancelled(stop);
        result.frequencies.push_back(std::move(frame));
        frequencyProgress.completedFrequencies = result.frequencies.size();
        report(BemSolvePhase::frequencyComplete);
    }
    return result;
}
} // namespace hornsculpt
