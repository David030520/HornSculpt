#include "wideband_operator.hpp"
#include "wideband_runtime.cuh"
#include <algorithm>
#include <climits>
#include <cufft.h>
#include <numbers>

namespace hornsculpt::wideband {
namespace {
using namespace detail;
constexpr unsigned modeTile = 8;
constexpr double pi = std::numbers::pi;
void checkedFft(cufftResult status, const char *operation) {
    if (status != CUFFT_SUCCESS)
        throw std::runtime_error(std::string(operation) + ": cuFFT " + std::to_string(status));
}
class FftPlan {
  public:
    FftPlan() {
        checkedFft(cufftCreate(&handle_), "Wideband FFT handle");
    }
    ~FftPlan() {
        cufftDestroy(handle_);
    }
    FftPlan(const FftPlan &) = delete;
    FftPlan &operator=(const FftPlan &) = delete;
    std::size_t prepare(unsigned length, std::size_t batches, cudaStream_t stream) {
        if (!length || length > INT_MAX || !batches || batches > INT_MAX)
            throw std::overflow_error(
                "Wideband FFT dimensions exceed the supported integer range.");
        checkedFft(cufftSetAutoAllocation(handle_, 0), "Wideband explicit FFT workspace");
        int n = static_cast<int>(length);
        checkedFft(cufftMakePlanMany(handle_, 1, &n, nullptr, 1, n, nullptr, 1, n, CUFFT_Z2Z,
                                     static_cast<int>(batches), &workBytes_),
                   "Wideband FFT plan");
        checkedFft(cufftSetStream(handle_, stream), "Wideband FFT stream");
        return workBytes_;
    }
    void workspace(void *data) {
        if (workBytes_)
            checkedFft(cufftSetWorkArea(handle_, data), "Wideband FFT workspace binding");
    }
    void execute(Z *values, int direction) {
        checkedFft(cufftExecZ2Z(handle_, values, values, direction), "Wideband azimuth FFT");
    }

  private:
    cufftHandle handle_{};
    std::size_t workBytes_{};
};
__device__ std::size_t triangular(unsigned l, unsigned m) {
    return static_cast<std::size_t>(l) * (l + 1) / 2 + m;
}
__device__ std::size_t harmonic(unsigned l, int m) {
    return static_cast<std::size_t>(l) * l + l + m;
}
__device__ double diagonal(unsigned m, double mu, const double *scales) {
    const double radius = sqrt(fmax(0., 1 - mu * mu));
    double result = .28209479177387814347;
    for (unsigned k = 1; k <= m; ++k)
        result *= scales[k] * radius;
    return result;
}
__global__ void associatedTile(const PolarSample *polar, unsigned polarCount, unsigned order,
                               const double *first, const double *second, const double *scales,
                               unsigned firstMode, double *table, std::size_t firstPoint,
                               std::size_t count) {
    const auto within = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (within >= count)
        return;
    const auto flat = firstPoint + within;
    const unsigned theta = static_cast<unsigned>(flat % polarCount),
                   m = firstMode + static_cast<unsigned>(flat / polarCount);
    const double mu = polar[theta].cosine;
    double previous = 0, current = diagonal(m, mu, scales);
    auto *row = table + flat * (order + 1);
    for (unsigned l = m; l <= order; ++l) {
        row[l] = current;
        if (l == order)
            break;
        const double next = l == m ? sqrt(2. * m + 3) * mu * current
                                   : first[triangular(l + 1, m)] * mu * current -
                                         second[triangular(l + 1, m)] * previous;
        previous = current;
        current = next;
    }
}
__global__ void projectTile(const Z *fourier, unsigned azimuthCount, const PolarSample *polar,
                            unsigned polarCount, unsigned order, const double *table,
                            unsigned firstMode, unsigned modes, Z *coefficients,
                            std::size_t firstPoint, std::size_t count) {
    const auto within = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (within >= count)
        return;
    auto flat = firstPoint + within;
    const unsigned l = static_cast<unsigned>(flat % (order + 1));
    flat /= order + 1;
    const unsigned localMode = static_cast<unsigned>(flat % modes), m = firstMode + localMode;
    if (l < m)
        return;
    const auto batch = flat / modes;
    Z positive{}, negative{};
    const double sign = m % 2 ? -1. : 1.;
    for (unsigned theta = 0; theta < polarCount; ++theta) {
        const double p =
            table[(static_cast<std::size_t>(localMode) * polarCount + theta) * (order + 1) + l] *
            (2 * pi * polar[theta].weight / azimuthCount);
        const auto base = (batch * polarCount + theta) * azimuthCount;
        const auto plus = fourier[base + m];
        positive.x += p * plus.x;
        positive.y += p * plus.y;
        if (m) {
            const auto minus = fourier[base + azimuthCount - m];
            negative.x += sign * p * minus.x;
            negative.y += sign * p * minus.y;
        }
    }
    const auto base = batch * static_cast<std::size_t>(order + 1) * (order + 1);
    coefficients[base + harmonic(l, static_cast<int>(m))] = positive;
    if (m)
        coefficients[base + harmonic(l, -static_cast<int>(m))] = negative;
}
__global__ void evaluateTile(const Z *coefficients, unsigned azimuthCount, const PolarSample *polar,
                             unsigned polarCount, unsigned order, const double *first,
                             const double *second, const double *scales, unsigned firstMode,
                             unsigned modes, Z *fourier, std::size_t firstPoint,
                             std::size_t count) {
    const auto within = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (within >= count)
        return;
    auto flat = firstPoint + within;
    const unsigned theta = static_cast<unsigned>(flat % polarCount);
    flat /= polarCount;
    const unsigned m = firstMode + static_cast<unsigned>(flat % modes);
    const auto batch = flat / modes;
    const double mu = polar[theta].cosine, sign = m % 2 ? -1. : 1.;
    double previous = 0, current = diagonal(m, mu, scales);
    Z positive{}, negative{};
    const auto coefficientBase = batch * static_cast<std::size_t>(order + 1) * (order + 1);
    for (unsigned l = m; l <= order; ++l) {
        const auto plus = coefficients[coefficientBase + harmonic(l, static_cast<int>(m))];
        positive.x += current * plus.x;
        positive.y += current * plus.y;
        if (m) {
            const auto minus = coefficients[coefficientBase + harmonic(l, -static_cast<int>(m))];
            negative.x += sign * current * minus.x;
            negative.y += sign * current * minus.y;
        }
        if (l == order)
            break;
        const double next = l == m ? sqrt(2. * m + 3) * mu * current
                                   : first[triangular(l + 1, m)] * mu * current -
                                         second[triangular(l + 1, m)] * previous;
        previous = current;
        current = next;
    }
    const auto outputBase = (batch * polarCount + theta) * azimuthCount;
    fourier[outputBase + m] = positive;
    if (m)
        fourier[outputBase + azimuthCount - m] = negative;
}
__global__ void shiftSpectrum(Z *values, const PolarSample *polar, unsigned polarCount,
                              unsigned azimuthCount, Vec3 phaseShift, std::size_t firstPoint,
                              std::size_t count) {
    const auto within = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (within >= count)
        return;
    const auto flat = firstPoint + within;
    const auto point = flat % (static_cast<std::size_t>(polarCount) * azimuthCount);
    const unsigned azimuth = static_cast<unsigned>(point % azimuthCount),
                   theta = static_cast<unsigned>(point / azimuthCount);
    const double mu = polar[theta].cosine, radius = sqrt(fmax(0., 1 - mu * mu));
    double sine, cosine;
    sincos(2 * pi * azimuth / azimuthCount, &sine, &cosine);
    const double angle = radius * (cosine * phaseShift.x + sine * phaseShift.y) + mu * phaseShift.z;
    sincos(angle, &sine, &cosine);
    const auto value = values[flat];
    values[flat] =
        make_cuDoubleComplex(value.x * cosine - value.y * sine, value.x * sine + value.y * cosine);
}

template <bool outgoing>
__global__ void scaleSpherical(const Z *input, Z *output, unsigned order, double ka,
                               std::size_t firstPoint, std::size_t count) {
    const auto within = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (within >= count)
        return;
    const auto flat = firstPoint + within,
               coefficient = flat % (static_cast<std::size_t>(order + 1) * (order + 1));
    const unsigned degree = static_cast<unsigned>(sqrt(static_cast<double>(coefficient)));
    // A scaled source moment is normalized by (2l+1)!!/(ka)^l.
    // Rayleigh's identity supplies 4*pi*i^l for the outgoing spectrum and
    // 4*pi*(-i)^l for the incoming local coefficients.
    double factor = 4 * pi;
    for (unsigned l = 1; l <= degree; ++l)
        factor *= ka / (2. * l + 1);
    const auto value = input[flat];
    const unsigned quadrant = outgoing ? degree % 4 : (4 - degree % 4) % 4;
    if (quadrant == 0)
        output[flat] = make_cuDoubleComplex(factor * value.x, factor * value.y);
    else if (quadrant == 1)
        output[flat] = make_cuDoubleComplex(-factor * value.y, factor * value.x);
    else if (quadrant == 2)
        output[flat] = make_cuDoubleComplex(-factor * value.x, -factor * value.y);
    else
        output[flat] = make_cuDoubleComplex(factor * value.y, -factor * value.x);
}

} // namespace

struct AngularTransferOperator::Impl {
    Stream stream;
    DeviceBuffer<std::byte> workspace;
    FftPlan forward, inverse;
    DeviceBuffer<PolarSample> sourcePolar, targetPolar;
    DeviceBuffer<double> recurrenceFirst, recurrenceSecond, diagonalScale, associated;
    DeviceBuffer<Z> source, target, harmonics;
    DeviceBuffer<int> error;
    unsigned order{}, batches{}, sourceTheta{}, targetTheta{}, sourceAzimuth{}, targetAzimuth{};
    std::size_t sourceCount{}, targetCount{}, bytes{};
    Impl(const AngularTransferPlan &plan, unsigned batchCount, std::stop_token stop,
         cudaStream_t borrowed)
        : stream(borrowed) {
        cancelled(stop);
        const auto orderCount = static_cast<std::size_t>(plan.bandlimit) + 1;
        const auto recurrenceCount = product(orderCount, orderCount + 1) / 2;
        if (!batchCount || plan.bandlimit > 1024 || plan.sourcePolar.size() <= plan.bandlimit ||
            plan.targetPolar.size() <= plan.bandlimit ||
            plan.sourceAzimuthCount != 2 * plan.sourcePolar.size() ||
            plan.targetAzimuthCount != 2 * plan.targetPolar.size() ||
            plan.recurrenceFirst.size() != recurrenceCount ||
            plan.recurrenceSecond.size() != recurrenceCount ||
            plan.diagonalScale.size() != orderCount)
            throw std::invalid_argument("Wideband GPU angular plan has inconsistent dimensions.");
        for (const auto *polar : {&plan.sourcePolar, &plan.targetPolar})
            for (const auto point : *polar)
                if (!std::isfinite(point.cosine) || std::abs(point.cosine) > 1 ||
                    !std::isfinite(point.weight) || point.weight <= 0)
                    throw std::invalid_argument("Wideband GPU polar rule is invalid.");
        for (const auto *array :
             {&plan.recurrenceFirst, &plan.recurrenceSecond, &plan.diagonalScale})
            for (const auto value : *array)
                if (!std::isfinite(value))
                    throw std::invalid_argument("Wideband GPU Legendre coefficient is non-finite.");
        order = plan.bandlimit;
        batches = batchCount;
        sourceTheta = static_cast<unsigned>(plan.sourcePolar.size());
        targetTheta = static_cast<unsigned>(plan.targetPolar.size());
        sourceAzimuth = plan.sourceAzimuthCount;
        targetAzimuth = plan.targetAzimuthCount;
        sourceCount = product(product(batches, sourceTheta), sourceAzimuth);
        targetCount = product(product(batches, targetTheta), targetAzimuth);
        const auto harmonicCount = product(batches, product(orderCount, orderCount));
        const auto associatedCount =
            product(product(std::min(modeTile, order + 1), sourceTheta), orderCount);
        const auto forwardWork =
                       forward.prepare(sourceAzimuth, product(batches, sourceTheta), stream),
                   inverseWork =
                       inverse.prepare(targetAzimuth, product(batches, targetTheta), stream);
        const auto workBytes = std::max(forwardWork, inverseWork);
        bytes = sum(product(sum(sourceCount, sum(targetCount, harmonicCount)), sizeof(Z)),
                    product(associatedCount, sizeof(double)));
        bytes = sum(bytes, product(sum(plan.sourcePolar.size(), plan.targetPolar.size()),
                                   sizeof(PolarSample)));
        bytes = sum(bytes, product(sum(product(2, recurrenceCount), orderCount), sizeof(double)));
        bytes = sum(sum(bytes, workBytes), sizeof(int));
        std::size_t free = 0, total = 0;
        checked(cudaMemGetInfo(&free, &total), "Wideband angular GPU memory query");
        if (bytes > free - free / 5)
            throw std::runtime_error(
                "Wideband angular transfer exceeds the available GPU memory budget.");
        sourcePolar.allocate(plan.sourcePolar.size());
        targetPolar.allocate(plan.targetPolar.size());
        recurrenceFirst.allocate(recurrenceCount);
        recurrenceSecond.allocate(recurrenceCount);
        diagonalScale.allocate(orderCount);
        associated.allocate(associatedCount);
        source.allocate(sourceCount);
        target.allocate(targetCount);
        harmonics.allocate(harmonicCount);
        workspace.allocate(workBytes);
        error.allocate(1);
        forward.workspace(workspace.data());
        inverse.workspace(workspace.data());
        sourcePolar.upload(plan.sourcePolar, stream);
        targetPolar.upload(plan.targetPolar, stream);
        recurrenceFirst.upload(plan.recurrenceFirst, stream);
        recurrenceSecond.upload(plan.recurrenceSecond, stream);
        diagonalScale.upload(plan.diagonalScale, stream);
        checked(cudaStreamSynchronize(stream), "Wideband angular plan preparation");
        cancelled(stop);
    }
    template <class Function>
    void tiles(std::size_t points, unsigned workPerPoint, std::stop_token stop, Function launch) {
        const auto capacity = std::max(
            std::size_t{512}, std::min(maximumBatchPoints, std::size_t{8000000} / workPerPoint));
        for (std::size_t first = 0; first < points; first += capacity) {
            cancelled(stop);
            const auto count = std::min(capacity, points - first);
            launch(static_cast<unsigned>((count + threads - 1) / threads), first, count);
            checked(cudaGetLastError(), "Wideband angular kernel launch");
            checked(cudaStreamSynchronize(stream), "Wideband angular batch");
        }
    }
    void applyDevice(const Z *input, Z *output, int phaseStage, double k, Vec3 displacement,
                     std::stop_token stop) {
        cancelled(stop);
        const double sign = phaseStage < 0 ? -1. : 1.;
        const Vec3 phaseShift{sign * k * displacement.x, sign * k * displacement.y,
                              sign * k * displacement.z};
        if (!std::isfinite(k) || k < 0 || !std::isfinite(phaseShift.x) ||
            !std::isfinite(phaseShift.y) || !std::isfinite(phaseShift.z))
            throw std::invalid_argument("Wideband angular center shift is non-finite or invalid.");
        if (!input || !output)
            throw std::invalid_argument("Wideband angular device pointer is null.");
        if (input != source.data())
            checked(cudaMemcpyAsync(source.data(), input, product(sourceCount, sizeof(Z)),
                                    cudaMemcpyDeviceToDevice, stream),
                    "Wideband device angular input");
        if (phaseStage < 0)
            tiles(sourceCount, 1, stop, [&](unsigned blocks, std::size_t first, std::size_t count) {
                shiftSpectrum<<<blocks, threads, 0, stream>>>(source.data(), sourcePolar.data(),
                                                              sourceTheta, sourceAzimuth,
                                                              phaseShift, first, count);
            });
        forward.execute(source.data(), CUFFT_FORWARD);
        checked(cudaMemsetAsync(target.data(), 0, product(targetCount, sizeof(Z)), stream),
                "Wideband target Fourier initialization");
        for (unsigned firstMode = 0; firstMode <= order; firstMode += modeTile) {
            const unsigned modes = std::min(modeTile, order - firstMode + 1);
            tiles(product(modes, sourceTheta), order + 1, stop,
                  [&](unsigned blocks, std::size_t first, std::size_t count) {
                      associatedTile<<<blocks, threads, 0, stream>>>(
                          sourcePolar.data(), sourceTheta, order, recurrenceFirst.data(),
                          recurrenceSecond.data(), diagonalScale.data(), firstMode,
                          associated.data(), first, count);
                  });
            tiles(product(product(batches, modes), order + 1), sourceTheta, stop,
                  [&](unsigned blocks, std::size_t first, std::size_t count) {
                      projectTile<<<blocks, threads, 0, stream>>>(
                          source.data(), sourceAzimuth, sourcePolar.data(), sourceTheta, order,
                          associated.data(), firstMode, modes, harmonics.data(), first, count);
                  });
            tiles(product(product(batches, modes), targetTheta), order + 1, stop,
                  [&](unsigned blocks, std::size_t first, std::size_t count) {
                      evaluateTile<<<blocks, threads, 0, stream>>>(
                          harmonics.data(), targetAzimuth, targetPolar.data(), targetTheta, order,
                          recurrenceFirst.data(), recurrenceSecond.data(), diagonalScale.data(),
                          firstMode, modes, target.data(), first, count);
                  });
        }
        inverse.execute(target.data(), CUFFT_INVERSE);
        if (phaseStage > 0)
            tiles(targetCount, 1, stop, [&](unsigned blocks, std::size_t first, std::size_t count) {
                shiftSpectrum<<<blocks, threads, 0, stream>>>(target.data(), targetPolar.data(),
                                                              targetTheta, targetAzimuth,
                                                              phaseShift, first, count);
            });
        if (output != target.data())
            checked(cudaMemcpyAsync(output, target.data(), product(targetCount, sizeof(Z)),
                                    cudaMemcpyDeviceToDevice, stream),
                    "Wideband device angular output");
        verifySpectrum(output, targetCount, error.data(), stream, stop);
    }
    void handoff(const Z *input, Z *output, bool outgoing, double k, double radius,
                 Vec3 displacement, std::stop_token stop) {
        cancelled(stop);
        const double ka = k * radius, sign = outgoing ? 1. : -1.;
        const Vec3 shift{sign * k * displacement.x, sign * k * displacement.y,
                         sign * k * displacement.z};
        if (!input || !output || input == output || order > maximumDeviceScaledOrder ||
            !std::isfinite(k) || k <= 0 || !std::isfinite(radius) || radius <= 0 ||
            !std::isfinite(ka) || ka > maximumScaledElectricalRadius || !std::isfinite(shift.x) ||
            !std::isfinite(shift.y) || !std::isfinite(shift.z))
            throw std::invalid_argument("Wideband LF/HF handoff dimensions or scales are invalid.");
        const auto harmonicCount = product(batches, product(order + 1, order + 1));
        if (outgoing) {
            tiles(harmonicCount, order + 1, stop,
                  [&](unsigned blocks, std::size_t first, std::size_t count) {
                      scaleSpherical<true><<<blocks, threads, 0, stream>>>(input, harmonics.data(),
                                                                           order, ka, first, count);
                  });
            checked(cudaMemsetAsync(target.data(), 0, product(targetCount, sizeof(Z)), stream),
                    "Wideband handoff Fourier initialization");
            for (unsigned firstMode = 0; firstMode <= order; firstMode += modeTile) {
                const unsigned modes = std::min(modeTile, order - firstMode + 1);
                tiles(product(product(batches, modes), targetTheta), order + 1, stop,
                      [&](unsigned blocks, std::size_t first, std::size_t count) {
                          evaluateTile<<<blocks, threads, 0, stream>>>(
                              harmonics.data(), targetAzimuth, targetPolar.data(), targetTheta,
                              order, recurrenceFirst.data(), recurrenceSecond.data(),
                              diagonalScale.data(), firstMode, modes, target.data(), first, count);
                      });
            }
            inverse.execute(target.data(), CUFFT_INVERSE);
            tiles(targetCount, 1, stop, [&](unsigned blocks, std::size_t first, std::size_t count) {
                shiftSpectrum<<<blocks, threads, 0, stream>>>(target.data(), targetPolar.data(),
                                                              targetTheta, targetAzimuth, shift,
                                                              first, count);
            });
            checked(cudaMemcpyAsync(output, target.data(), product(targetCount, sizeof(Z)),
                                    cudaMemcpyDeviceToDevice, stream),
                    "Wideband handoff spectrum output");
            verifySpectrum(output, targetCount, error.data(), stream, stop);
        } else {
            checked(cudaMemcpyAsync(source.data(), input, product(sourceCount, sizeof(Z)),
                                    cudaMemcpyDeviceToDevice, stream),
                    "Wideband handoff spectrum input");
            tiles(sourceCount, 1, stop, [&](unsigned blocks, std::size_t first, std::size_t count) {
                shiftSpectrum<<<blocks, threads, 0, stream>>>(source.data(), sourcePolar.data(),
                                                              sourceTheta, sourceAzimuth, shift,
                                                              first, count);
            });
            forward.execute(source.data(), CUFFT_FORWARD);
            for (unsigned firstMode = 0; firstMode <= order; firstMode += modeTile) {
                const unsigned modes = std::min(modeTile, order - firstMode + 1);
                tiles(product(modes, sourceTheta), order + 1, stop,
                      [&](unsigned blocks, std::size_t first, std::size_t count) {
                          associatedTile<<<blocks, threads, 0, stream>>>(
                              sourcePolar.data(), sourceTheta, order, recurrenceFirst.data(),
                              recurrenceSecond.data(), diagonalScale.data(), firstMode,
                              associated.data(), first, count);
                      });
                tiles(product(product(batches, modes), order + 1), sourceTheta, stop,
                      [&](unsigned blocks, std::size_t first, std::size_t count) {
                          projectTile<<<blocks, threads, 0, stream>>>(
                              source.data(), sourceAzimuth, sourcePolar.data(), sourceTheta, order,
                              associated.data(), firstMode, modes, harmonics.data(), first, count);
                      });
            }
            tiles(harmonicCount, order + 1, stop,
                  [&](unsigned blocks, std::size_t first, std::size_t count) {
                      scaleSpherical<false><<<blocks, threads, 0, stream>>>(
                          harmonics.data(), output, order, ka, first, count);
                  });
            verifySpectrum(output, harmonicCount, error.data(), stream, stop);
        }
    }
    std::vector<Complex> apply(std::span<const Complex> input, int phaseStage, double k,
                               Vec3 displacement, std::stop_token stop) {
        if (input.size() != sourceCount)
            throw std::invalid_argument("Wideband batched angular spectrum has an invalid size.");
        const auto packedInput = packed(input);
        source.upload(packedInput, stream);
        applyDevice(source.data(), target.data(), phaseStage, k, displacement, stop);
        const auto output = target.download(stream);
        std::vector<Complex> result;
        result.reserve(output.size());
        for (const auto value : output)
            result.emplace_back(value.x, value.y);
        return result;
    }
};
AngularTransferOperator::AngularTransferOperator(const AngularTransferPlan &plan,
                                                 unsigned batchCount, std::stop_token stop,
                                                 cudaStream_t stream)
    : impl_(std::make_unique<Impl>(plan, batchCount, stop, stream)) {}
AngularTransferOperator::~AngularTransferOperator() = default;
std::vector<Complex> AngularTransferOperator::apply(std::span<const Complex> source,
                                                    std::stop_token stop) {
    return impl_->apply(source, 0, 0, {}, stop);
}
std::vector<Complex> AngularTransferOperator::outgoingToParent(std::span<const Complex> source,
                                                               double wavenumber,
                                                               Vec3 childMinusParent,
                                                               std::stop_token stop) {
    return impl_->apply(source, 1, wavenumber, childMinusParent, stop);
}
std::vector<Complex> AngularTransferOperator::incomingToChild(std::span<const Complex> source,
                                                              double wavenumber,
                                                              Vec3 childMinusParent,
                                                              std::stop_token stop) {
    return impl_->apply(source, -1, wavenumber, childMinusParent, stop);
}
void AngularTransferOperator::applyDevice(const cuDoubleComplex *source, cuDoubleComplex *output,
                                          std::stop_token stop) {
    impl_->applyDevice(source, output, 0, 0, {}, stop);
}
void AngularTransferOperator::outgoingToParentDevice(const cuDoubleComplex *source,
                                                     cuDoubleComplex *output, double k, Vec3 shift,
                                                     std::stop_token stop) {
    impl_->applyDevice(source, output, 1, k, shift, stop);
}
void AngularTransferOperator::incomingToChildDevice(const cuDoubleComplex *source,
                                                    cuDoubleComplex *output, double k, Vec3 shift,
                                                    std::stop_token stop) {
    impl_->applyDevice(source, output, -1, k, shift, stop);
}
cudaStream_t AngularTransferOperator::stream() const noexcept {
    return impl_->stream;
}

std::size_t AngularTransferOperator::deviceBytes() const noexcept {
    return impl_->bytes;
}
void AngularTransferOperator::scaledOutgoingToParentDevice(const Z *multipole, Z *spectrum,
                                                           double k, double radius, Vec3 shift,
                                                           std::stop_token stop) {
    impl_->handoff(multipole, spectrum, true, k, radius, shift, stop);
}
void AngularTransferOperator::incomingToScaledChildDevice(const Z *spectrum, Z *local, double k,
                                                          double radius, Vec3 shift,
                                                          std::stop_token stop) {
    impl_->handoff(spectrum, local, false, k, radius, shift, stop);
}
} // namespace hornsculpt::wideband
