#include "support_preconditioner.hpp"
#include "wideband_runtime.cuh"
#include <cublas_v2.h>
#include <limits>

namespace hornsculpt::wideband {
namespace {
using namespace detail;
void checked(cublasStatus_t status, const char *operation) {
    if (status != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error(std::string(operation) + ": cuBLAS " + std::to_string(status));
}
class Blas {
  public:
    explicit Blas(cudaStream_t stream) {
        checked(cublasCreate(&value_), "Support preconditioner cuBLAS creation");
        const auto status = cublasSetStream(value_, stream);
        if (status != CUBLAS_STATUS_SUCCESS) {
            cublasDestroy(value_);
            checked(status, "Support preconditioner cuBLAS stream");
        }
    }
    ~Blas() {
        cublasDestroy(value_);
    }
    operator cublasHandle_t() const noexcept {
        return value_;
    }

  private:
    cublasHandle_t value_{};
};
__global__ void extractBlocks(SupportDeviceCsr matrix, const std::uint32_t *order, unsigned width,
                              Z *factors, Z **factorPointers, Z *right, Z **rightPointers,
                              std::size_t first, std::size_t count) {
    const auto local = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (local >= count)
        return;
    const auto flat = first + local, entries = static_cast<std::size_t>(width) * width,
               block = flat / entries;
    const unsigned row = static_cast<unsigned>(flat % width),
                   column = static_cast<unsigned>((flat / width) % width);
    const auto sortedRow = block * width + row, sortedColumn = block * width + column;
    Z value{};
    if (sortedRow < matrix.rows && sortedColumn < matrix.rows) {
        const auto globalRow = order[sortedRow], globalColumn = order[sortedColumn];
        auto begin = matrix.rowOffsets[globalRow], end = matrix.rowOffsets[globalRow + 1];
        const auto last = end;
        while (begin < end) {
            const auto middle = begin + (end - begin) / 2;
            if (matrix.columns[middle] < globalColumn)
                begin = middle + 1;
            else
                end = middle;
        }
        if (begin < last && matrix.columns[begin] == globalColumn)
            value = matrix.values[begin];
    } else if (row == column)
        value.x = 1;
    factors[flat] = value;
    if (!row && !column) {
        factorPointers[block] = factors + block * entries;
        rightPointers[block] = right + block * width;
    }
}
__global__ void permute(const Z *input, Z *output, const std::uint32_t *order, std::size_t count,
                        std::size_t padded, bool inverse) {
    const auto index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= padded)
        return;
    if (inverse) {
        if (index < count)
            output[order[index]] = input[index];
    } else
        output[index] = index < count ? input[order[index]] : Z{};
}
} // namespace
class SupportBlockJacobi::Implementation {
  public:
    Implementation(SupportDeviceCsr matrix, std::span<const std::uint32_t> permutation,
                   unsigned selectedWidth, std::stop_token stop, cudaStream_t borrowed)
        : stream(borrowed), blas(stream), width(selectedWidth), count(matrix.rows) {
        cancelled(stop);
        if (!count || matrix.columnCount != count || permutation.size() != count ||
            !matrix.rowOffsets || !matrix.columns || !matrix.values || !width || width > 2048 ||
            (width & (width - 1)))
            throw std::invalid_argument(
                "Support preconditioner dimensions or block width are invalid.");
        std::vector<bool> visited(count);
        for (const auto index : permutation) {
            if (index >= count || visited[index])
                throw std::invalid_argument(
                    "Support preconditioner pressure order is not a permutation.");
            visited[index] = true;
        }
        blocks = count / width + (count % width != 0);
        if (blocks > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::invalid_argument(
                "Support preconditioner block count exceeds cuBLAS limits.");
        padded = product(blocks, width);
        const auto entries = product(padded, width);
        bytes = sum(
            sum(product(sum(entries, padded), sizeof(Z)), product(count, sizeof(std::uint32_t))),
            sum(product(padded, sizeof(int)), product(product(2, blocks), sizeof(Z *))));
        peak = sum(bytes, sum(product(blocks, sizeof(int)), sizeof(int)));
        std::size_t free{}, total{};
        detail::checked(cudaMemGetInfo(&free, &total), "Support preconditioner memory query");
        if (peak > free - free / 5)
            throw std::runtime_error("Support preconditioner exceeds available GPU memory budget.");
        order.allocate(count);
        order.upload(permutation, stream);
        factors.allocate(entries);
        right.allocate(padded);
        pivots.allocate(padded);
        factorPointers.allocate(blocks);
        rightPointers.allocate(blocks);
        error.allocate(1);
        DeviceBuffer<int> factorInfo;
        factorInfo.allocate(blocks);
        for (std::size_t first = 0; first < entries; first += maximumBatchPoints) {
            cancelled(stop);
            const auto batch = std::min(maximumBatchPoints, entries - first);
            extractBlocks<<<static_cast<unsigned>((batch + threads - 1) / threads), threads, 0,
                            stream>>>(matrix, order.data(), width, factors.data(),
                                      factorPointers.data(), right.data(), rightPointers.data(),
                                      first, batch);
            detail::checked(cudaGetLastError(), "Support preconditioner extraction launch");
            detail::checked(cudaStreamSynchronize(stream), "Support preconditioner extraction");
        }
        checked(cublasZgetrfBatched(blas, static_cast<int>(width), factorPointers.data(),
                                    static_cast<int>(width), pivots.data(), factorInfo.data(),
                                    static_cast<int>(blocks)),
                "Support preconditioner LU factorization");
        const auto info = factorInfo.download(stream);
        for (std::size_t block = 0; block < blocks; ++block)
            if (info[block])
                throw std::runtime_error("Support preconditioner LU failed in spatial block " +
                                         std::to_string(block) + " with status " +
                                         std::to_string(info[block]) + ".");
        verifySpectrum(factors.data(), entries, error.data(), stream, stop);
        cancelled(stop);
    }
    void apply(Z *values, std::stop_token stop) {
        cancelled(stop);
        if (!values)
            throw std::invalid_argument("Support preconditioner requires a device vector.");
        const auto grid = static_cast<unsigned>((padded + threads - 1) / threads);
        permute<<<grid, threads, 0, stream>>>(values, right.data(), order.data(), count, padded,
                                              false);
        detail::checked(cudaGetLastError(), "Support preconditioner input permutation");
        int info{};
        checked(cublasZgetrsBatched(blas, CUBLAS_OP_N, static_cast<int>(width), 1,
                                    factorPointers.data(), static_cast<int>(width), pivots.data(),
                                    rightPointers.data(), static_cast<int>(width), &info,
                                    static_cast<int>(blocks)),
                "Support preconditioner triangular solve");
        if (info)
            throw std::runtime_error(
                "Support preconditioner triangular solve has invalid arguments.");
        verifySpectrum(right.data(), padded, error.data(), stream, stop);
        permute<<<grid, threads, 0, stream>>>(right.data(), values, order.data(), count, padded,
                                              true);
        detail::checked(cudaGetLastError(), "Support preconditioner output permutation");
        detail::checked(cudaStreamSynchronize(stream), "Support preconditioner completion");
        cancelled(stop);
    }
    Stream stream;
    Blas blas;
    unsigned width{};
    std::size_t count{}, blocks{}, padded{}, bytes{}, peak{};
    DeviceBuffer<std::uint32_t> order;
    DeviceBuffer<Z> factors, right;
    DeviceBuffer<int> pivots, error;
    DeviceBuffer<Z *> factorPointers, rightPointers;
};
SupportBlockJacobi::SupportBlockJacobi(SupportDeviceCsr matrix,
                                       std::span<const std::uint32_t> pressureOrder, unsigned width,
                                       std::stop_token stop, cudaStream_t stream)
    : implementation_(
          std::make_unique<Implementation>(matrix, pressureOrder, width, stop, stream)) {}
SupportBlockJacobi::~SupportBlockJacobi() = default;
void SupportBlockJacobi::applyDevice(cuDoubleComplex *values, std::stop_token stop) {
    implementation_->apply(values, stop);
}
std::size_t SupportBlockJacobi::deviceBytes() const noexcept {
    return implementation_->bytes + sizeof(int);
}
std::size_t SupportBlockJacobi::peakDeviceBytes() const noexcept {
    return implementation_->peak;
}
} // namespace hornsculpt::wideband
