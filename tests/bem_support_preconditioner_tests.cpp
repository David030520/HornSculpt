#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cuda/support_preconditioner.hpp>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string_view>

namespace {
using namespace hornsculpt::wideband;
using C = std::complex<double>;
void check(cudaError_t status) {
    if (status != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(status));
}
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class T> class Device {
  public:
    explicit Device(std::span<const T> values) {
        check(cudaMalloc(reinterpret_cast<void **>(&data), values.size_bytes()));
        try {
            check(cudaMemcpy(data, values.data(), values.size_bytes(), cudaMemcpyHostToDevice));
            check(cudaStreamSynchronize(nullptr));
        } catch (...) {
            cudaFree(data);
            throw;
        }
    }
    ~Device() {
        cudaFree(data);
    }
    Device(const Device &) = delete;
    Device &operator=(const Device &) = delete;
    T *data{};
};
void checkBlocks(unsigned width, bool sparsePivot = false) {
    // Exercise complete physical blocks and a real eleven-DOF tail for every
    // width. A smaller matrix would only test identity padding at large widths.
    const unsigned n = width <= 128 ? 139 : width + 11;
    require(n > width && n % width == 11 && std::gcd(n, 41u) == 1,
            "Preconditioner fixture needs a full block, tail and coprime permutation.");
    std::vector<std::uint32_t> permutation(n), inverse(n), offsets(n + 1), columns;
    std::iota(permutation.begin(), permutation.end(), 0);
    // Coprime stride interleaves distant global indices, so applying blocks in
    // unpermuted numbering cannot accidentally pass the test.
    for (unsigned i = 0; i < n; ++i) {
        permutation[i] = (i * 41) % n;
        inverse[permutation[i]] = i;
    }
    std::vector<cuDoubleComplex> entries, right(n);
    std::vector<C> expected(n);
    for (unsigned i = 0; i < n; ++i)
        expected[i] = {.1 + std::sin(i * .71), -.2 + std::cos(i * .37)};
    for (unsigned row = 0; row < n; ++row) {
        C rhs{};
        offsets[row] = static_cast<std::uint32_t>(columns.size());
        for (unsigned column = 0; column < n; ++column) {
            C value = row == column ? C{4 + row * .001, .3}
                                    : C{.001 * std::sin(row + 3. * column),
                                        .001 * std::cos(2. * row - column)};
            if (sparsePivot) {
                const auto sortedRow = inverse[row], sortedColumn = inverse[column];
                // Both the full block and the real tail require a row pivot.
                // Most other in-block entries are absent from CSR, not small
                // dense substitutes; extraction must supply exact zeros.
                if (row == column && sortedRow % width == 0)
                    value = {};
                else if (sortedRow % width == 0 && sortedColumn == sortedRow + 1)
                    value = {2, .25};
                else if (sortedColumn % width == 0 && sortedRow == sortedColumn + 1)
                    value = {1, -.125};
                else if (row != column && (row + 3 * column) % 29 != 0)
                    continue;
            }
            columns.push_back(column);
            entries.push_back(make_cuDoubleComplex(value.real(), value.imag()));
            if (inverse[row] / width == inverse[column] / width)
                rhs += value * expected[column];
        }
        right[row] = make_cuDoubleComplex(rhs.real(), rhs.imag());
    }
    offsets[n] = static_cast<std::uint32_t>(columns.size());
    Device<std::uint32_t> deviceOffsets(offsets), deviceColumns(columns);
    Device<cuDoubleComplex> deviceEntries(entries), vector(right);
    SupportBlockJacobi operation({deviceOffsets.data, deviceColumns.data, deviceEntries.data, n, n,
                                  static_cast<std::uint32_t>(entries.size())},
                                 permutation, width);
    operation.applyDevice(vector.data);
    check(
        cudaMemcpy(right.data(), vector.data, n * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
    double error{}, magnitude{};
    for (unsigned i = 0; i < n; ++i) {
        error += std::norm(C{right[i].x, right[i].y} - expected[i]);
        magnitude += std::norm(expected[i]);
    }
    const double relative = std::sqrt(error / magnitude);
    std::cout << "Spatial block-Jacobi width=" << width << " rows=" << n
              << " sparse-pivot=" << sparsePivot << " relative error=" << relative
              << " live/peak bytes=" << operation.deviceBytes() << '/'
              << operation.peakDeviceBytes() << '\n';
    require(relative < 2e-14, "Spatial block-Jacobi permutation/padding/LU solve failed.");
    require(operation.deviceBytes() <= operation.peakDeviceBytes(),
            "Preconditioner residency exceeds its peak.");
    const std::size_t blocks = (n + width - 1) / width, padded = blocks * width,
                      live = (padded * width + padded) * sizeof(cuDoubleComplex) +
                             n * sizeof(std::uint32_t) + padded * sizeof(int) +
                             2 * blocks * sizeof(cuDoubleComplex *) + sizeof(int);
    require(operation.deviceBytes() == live &&
                operation.peakDeviceBytes() == live + blocks * sizeof(int),
            "Preconditioner factor/pivot/padding capacity accounting changed.");
    if (sparsePivot) {
        std::stop_source cancellation;
        cancellation.request_stop();
        bool applyCancelled = false, constructionCancelled = false;
        try {
            operation.applyDevice(vector.data, cancellation.get_token());
        } catch (const std::runtime_error &error) {
            applyCancelled = std::string(error.what()).find("cancelled") != std::string::npos;
        }
        try {
            SupportBlockJacobi cancelledOperation({deviceOffsets.data, deviceColumns.data,
                                                   deviceEntries.data, n, n,
                                                   static_cast<std::uint32_t>(entries.size())},
                                                  permutation, width, cancellation.get_token());
        } catch (const std::runtime_error &error) {
            constructionCancelled =
                std::string(error.what()).find("cancelled") != std::string::npos;
        }
        require(applyCancelled && constructionCancelled,
                "Cancelled spatial preconditioner construction/apply did not stop.");
        std::vector<cuDoubleComplex> afterCancellation(n);
        check(cudaMemcpy(afterCancellation.data(), vector.data, n * sizeof(cuDoubleComplex),
                         cudaMemcpyDeviceToHost));
        for (unsigned i = 0; i < n; ++i)
            require(afterCancellation[i].x == right[i].x && afterCancellation[i].y == right[i].y,
                    "Pre-cancelled preconditioner application changed its input.");
    }
    check(cudaMemset(vector.data, 0, n * sizeof(cuDoubleComplex)));
    check(cudaStreamSynchronize(nullptr));
    operation.applyDevice(vector.data);
    check(
        cudaMemcpy(right.data(), vector.data, n * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
    for (const auto value : right)
        require(value.x == 0 && value.y == 0, "Preconditioner retained a previous RHS.");
}
void checkPivotAndSingularity() {
    const std::vector<std::uint32_t> offsets{0, 2, 4, 5}, columns{0, 1, 0, 1, 2},
        permutation{0, 1, 2};
    const std::vector<cuDoubleComplex> values{{0, 0}, {2, 0}, {1, 0}, {3, 0}, {4, 0}};
    Device<std::uint32_t> deviceOffsets(offsets), deviceColumns(columns);
    Device<cuDoubleComplex> deviceValues(values);
    const SupportDeviceCsr matrix{
        deviceOffsets.data, deviceColumns.data, deviceValues.data, 3, 3, 5};
    SupportBlockJacobi operation(matrix, permutation, 2);
    const std::vector<cuDoubleComplex> right{{4, -2}, {7, -1}, {12, 8}};
    Device<cuDoubleComplex> vector(right);
    operation.applyDevice(vector.data);
    std::vector<cuDoubleComplex> result(3);
    check(cudaMemcpy(result.data(), vector.data, result.size() * sizeof(cuDoubleComplex),
                     cudaMemcpyDeviceToHost));
    const std::array<C, 3> expected{{{1, 2}, {2, -1}, {3, 2}}};
    for (unsigned i = 0; i < 3; ++i)
        require(std::abs(C{result[i].x, result[i].y} - expected[i]) < 1e-13,
                "Preconditioner failed pivoted complex triangular solve.");
    check(cudaMemset(deviceValues.data, 0, values.size() * sizeof(cuDoubleComplex)));
    check(cudaStreamSynchronize(nullptr));
    bool rejected = false;
    try {
        SupportBlockJacobi singular(matrix, permutation, 2);
    } catch (const std::runtime_error &error) {
        rejected = std::string(error.what()).find("LU failed") != std::string::npos;
    }
    require(rejected, "Singular spatial block was not reported clearly.");
}
void checkInvalidDescriptors() {
    const std::vector<std::uint32_t> offsets{0, 1, 2, 3}, columns{0, 1, 2}, permutation{0, 1, 2};
    const std::vector<cuDoubleComplex> values{{1, 0}, {2, 0}, {3, 0}};
    Device<std::uint32_t> deviceOffsets(offsets), deviceColumns(columns);
    Device<cuDoubleComplex> deviceValues(values);
    const SupportDeviceCsr matrix{
        deviceOffsets.data, deviceColumns.data, deviceValues.data, 3, 3, 3};
    const auto invalid = [&](SupportDeviceCsr candidate, std::span<const std::uint32_t> order,
                             unsigned width) {
        bool rejected = false;
        try {
            SupportBlockJacobi operation(candidate, order, width);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Invalid preconditioner width/CSR/permutation was accepted.");
    };
    for (const unsigned width : {0u, 3u, 4096u})
        invalid(matrix, permutation, width);
    auto missing = matrix;
    missing.rowOffsets = nullptr;
    invalid(missing, permutation, 2048);
    missing = matrix;
    missing.columns = nullptr;
    invalid(missing, permutation, 2048);
    missing = matrix;
    missing.values = nullptr;
    invalid(missing, permutation, 2048);
    auto rectangular = matrix;
    rectangular.columnCount = 2;
    invalid(rectangular, permutation, 2048);
    const std::vector<std::uint32_t> duplicate{0, 1, 1}, outside{0, 1, 3}, shortOrder{0, 1};
    invalid(matrix, duplicate, 2048);
    invalid(matrix, outside, 2048);
    invalid(matrix, shortOrder, 2048);
}
} // namespace
int main(int argc, char **argv) {
    try {
        unsigned selectedWidth{};
        if (argc != 1) {
            require(
                argc == 3 && std::string_view(argv[1]) == "--width",
                "Usage: bem_support_preconditioner_tests [--width 32|64|128|256|512|1024|2048]");
            const std::string_view argument(argv[2]);
            const auto parsed =
                std::from_chars(argument.data(), argument.data() + argument.size(), selectedWidth);
            require(parsed.ec == std::errc{} && parsed.ptr == argument.data() + argument.size() &&
                        (selectedWidth == 32 || selectedWidth == 64 || selectedWidth == 128 ||
                         selectedWidth == 256 || selectedWidth == 512 || selectedWidth == 1024 ||
                         selectedWidth == 2048),
                    "Invalid isolated preconditioner test width.");
        }
        int count{};
        if (cudaGetDeviceCount(&count) != cudaSuccess || !count)
            return 77;
        // A selected width permits bounded staged qualification without running
        // earlier widths again. CTest's unchanged no-argument run covers all.
        if (selectedWidth)
            checkBlocks(selectedWidth, true);
        else {
            for (const unsigned width : {32, 64, 128, 256})
                checkBlocks(width);
            for (const unsigned width : {256, 512, 1024, 2048})
                checkBlocks(width, true);
        }
        checkPivotAndSingularity();
        checkInvalidDescriptors();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
