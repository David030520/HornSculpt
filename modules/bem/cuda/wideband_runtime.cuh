#pragma once

#include "../fmm/wideband_transfer.hpp"
#include <algorithm>
#include <cmath>
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <limits>
#include <stdexcept>
#include <stop_token>
#include <string>
#include <utility>

namespace hornsculpt::wideband::detail {
using Z = cuDoubleComplex;
constexpr unsigned threads = 128;
constexpr std::size_t maximumBatchPoints = 262144;
constexpr unsigned maximumBatchInteractions = 32;
inline void checked(cudaError_t status, const char *operation) {
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
}
inline void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Wideband calculation cancelled.");
}
inline std::size_t product(std::size_t a, std::size_t b) {
    if (b && a > std::numeric_limits<std::size_t>::max() / b)
        throw std::overflow_error("Wideband level allocation size overflow.");
    return a * b;
}
inline std::size_t sum(std::size_t a, std::size_t b) {
    if (a > std::numeric_limits<std::size_t>::max() - b)
        throw std::overflow_error("Wideband level allocation size overflow.");
    return a + b;
}
class Stream {
  public:
    explicit Stream(cudaStream_t borrowed = nullptr) : value_(borrowed), owned_(!borrowed) {
        if (owned_)
            checked(cudaStreamCreateWithFlags(&value_, cudaStreamNonBlocking),
                    "Wideband CUDA stream");
    }
    ~Stream() {
        if (owned_ && value_)
            cudaStreamDestroy(value_);
    }
    Stream(const Stream &) = delete;
    Stream &operator=(const Stream &) = delete;
    operator cudaStream_t() const noexcept {
        return value_;
    }

  private:
    cudaStream_t value_{};
    bool owned_{};
};
template <class T> class DeviceBuffer {
  public:
    DeviceBuffer() = default;
    ~DeviceBuffer() {
        if (data_)
            cudaFree(data_);
    }
    DeviceBuffer(const DeviceBuffer &) = delete;
    DeviceBuffer &operator=(const DeviceBuffer &) = delete;
    DeviceBuffer(DeviceBuffer &&other) noexcept
        : data_(std::exchange(other.data_, nullptr)), count_(std::exchange(other.count_, 0)) {}
    DeviceBuffer &operator=(DeviceBuffer &&other) noexcept {
        if (this != &other) {
            if (data_)
                cudaFree(data_);
            data_ = std::exchange(other.data_, nullptr);
            count_ = std::exchange(other.count_, 0);
        }
        return *this;
    }
    void allocate(std::size_t count) {
        count_ = count;
        if (count)
            checked(cudaMalloc(reinterpret_cast<void **>(&data_), product(count, sizeof(T))),
                    "Wideband level allocation");
    }
    void upload(std::span<const T> values, cudaStream_t stream) {
        if (values.size() != count_)
            throw std::invalid_argument("Wideband upload has an inconsistent size.");
        if (count_)
            checked(cudaMemcpyAsync(data_, values.data(), product(count_, sizeof(T)),
                                    cudaMemcpyHostToDevice, stream),
                    "Wideband level upload");
    }
    std::vector<T> download(cudaStream_t stream) const {
        std::vector<T> result(count_);
        if (count_)
            checked(cudaMemcpyAsync(result.data(), data_, product(count_, sizeof(T)),
                                    cudaMemcpyDeviceToHost, stream),
                    "Wideband level download");
        checked(cudaStreamSynchronize(stream), "Wideband level download completion");
        return result;
    }
    T *data() const noexcept {
        return data_;
    }

  private:
    T *data_{};
    std::size_t count_{};
};
inline std::vector<Z> packed(std::span<const Complex> values) {
    std::vector<Z> result;
    result.reserve(values.size());
    for (const auto value : values) {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::invalid_argument("Wideband input contains a non-finite coefficient.");
        result.push_back(make_cuDoubleComplex(value.real(), value.imag()));
    }
    return result;
}

static __global__ void finiteSpectrum(const Z *values, std::size_t first, std::size_t count,
                                      int *error) {
    const auto within = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (within < count) {
        const auto value = values[first + within];
        if (!isfinite(value.x) || !isfinite(value.y))
            atomicExch(error, 1);
    }
}
// Append without resetting a caller-owned status: an invalid intermediate
// must remain visible even if a later stage overwrites its output.
inline void appendSpectrumCheck(const Z *values, std::size_t count, int *error, cudaStream_t stream,
                                std::stop_token stop) {
    for (std::size_t first = 0; first < count; first += maximumBatchPoints) {
        cancelled(stop);
        const auto points = std::min(maximumBatchPoints, count - first);
        finiteSpectrum<<<static_cast<unsigned>((points + threads - 1) / threads), threads, 0,
                         stream>>>(values, first, points, error);
        checked(cudaGetLastError(), "Wideband spectrum finite check");
    }
}
inline void completeSpectrumCheck(int *error, cudaStream_t stream, std::stop_token stop) {
    cancelled(stop);
    int status = 0;
    checked(cudaMemcpyAsync(&status, error, sizeof(status), cudaMemcpyDeviceToHost, stream),
            "Wideband scalar status download");
    checked(cudaStreamSynchronize(stream), "Wideband spectrum completion");
    if (status)
        throw std::runtime_error("Wideband device spectrum contains a non-finite coefficient.");
}
inline void verifySpectrum(const Z *values, std::size_t count, int *error, cudaStream_t stream,
                           std::stop_token stop) {
    checked(cudaMemsetAsync(error, 0, sizeof(int), stream), "Wideband spectrum status reset");
    appendSpectrumCheck(values, count, error, stream, stop);
    completeSpectrumCheck(error, stream, stop);
}

} // namespace hornsculpt::wideband::detail
