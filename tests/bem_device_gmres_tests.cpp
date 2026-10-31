#include <algorithm>
#include <cmath>
#include <complex>
#include <cublas_v2.h>
#include <cuda/device_gmres.hpp>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
using C = std::complex<double>;
using Z = cuDoubleComplex;
void check(cudaError_t status) {
    if (status != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(status));
}
void check(cublasStatus_t status) {
    if (status != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error("GMRES test cuBLAS failure.");
}
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
class Device {
  public:
    explicit Device(std::size_t count) {
        check(cudaMalloc(reinterpret_cast<void **>(&data), count * sizeof(Z)));
    }
    ~Device() {
        cudaFree(data);
    }
    Device(const Device &) = delete;
    Device &operator=(const Device &) = delete;
    Z *data{};
};
class Runtime {
  public:
    Runtime() {
        check(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
        try {
            check(cublasCreate(&blas));
            check(cublasSetStream(blas, stream));
        } catch (...) {
            if (blas)
                cublasDestroy(blas);
            cudaStreamDestroy(stream);
            throw;
        }
    }
    ~Runtime() {
        cublasDestroy(blas);
        cudaStreamDestroy(stream);
    }
    cudaStream_t stream{};
    cublasHandle_t blas{};
};
std::vector<C> denseReference(const std::vector<Z> &matrix, const std::vector<Z> &right) {
    const auto n = right.size();
    require(matrix.size() == n * n, "Dense GMRES reference dimensions disagree.");
    std::vector<C> a(n * n), result(n);
    for (std::size_t row = 0; row < n; ++row) {
        result[row] = {right[row].x, right[row].y};
        for (std::size_t column = 0; column < n; ++column) {
            const auto value = matrix[column * n + row];
            a[row * n + column] = {value.x, value.y};
        }
    }
    // Independent CPU Gaussian elimination with partial row pivoting. Neither
    // the CUDA matvec nor any Arnoldi/GMRES implementation is used here.
    for (std::size_t column = 0; column < n; ++column) {
        auto pivot = column;
        for (std::size_t row = column + 1; row < n; ++row)
            if (std::abs(a[row * n + column]) > std::abs(a[pivot * n + column]))
                pivot = row;
        const double magnitude = std::abs(a[pivot * n + column]);
        require(std::isfinite(magnitude) && magnitude > 0,
                "Dense GMRES reference is singular or nonfinite.");
        if (pivot != column) {
            for (std::size_t j = 0; j < n; ++j)
                std::swap(a[column * n + j], a[pivot * n + j]);
            std::swap(result[column], result[pivot]);
        }
        for (std::size_t row = column + 1; row < n; ++row) {
            const auto factor = a[row * n + column] / a[column * n + column];
            a[row * n + column] = 0;
            for (std::size_t j = column + 1; j < n; ++j)
                a[row * n + j] -= factor * a[column * n + j];
            result[row] -= factor * result[column];
        }
    }
    for (std::size_t remaining = n; remaining > 0; --remaining) {
        const auto row = remaining - 1;
        for (std::size_t column = row + 1; column < n; ++column)
            result[row] -= a[row * n + column] * result[column];
        result[row] /= a[row * n + row];
    }
    return result;
}
void checkRestartedSolve() {
    constexpr unsigned n = 17;
    Runtime runtime;
    std::vector<Z> matrix(n * n), right(n), actual(n);
    std::vector<C> expected(n);
    for (unsigned i = 0; i < n; ++i)
        expected[i] = {std::sin(.7 * i) + .1, std::cos(.3 * i) - .2};
    for (unsigned row = 0; row < n; ++row) {
        C rhs{};
        for (unsigned column = 0; column < n; ++column) {
            C value{.006 * std::sin(row + 2. * column), .004 * std::cos(3. * row - column)};
            if (row == column)
                value += C{3.5 + .03 * row, .2};
            if (column == row + 1)
                value += C{1.3, .7};
            if (row == column + 1)
                value += C{-.2, .6};
            matrix[column * n + row] = make_cuDoubleComplex(value.real(), value.imag());
            rhs += value * expected[column];
        }
        right[row] = make_cuDoubleComplex(rhs.real(), rhs.imag());
    }
    const auto reference = denseReference(matrix, right);
    double referenceError{}, referenceNorm{};
    for (unsigned i = 0; i < n; ++i) {
        referenceError += std::norm(reference[i] - expected[i]);
        referenceNorm += std::norm(expected[i]);
    }
    require(std::sqrt(referenceError / referenceNorm) < 2e-14,
            "Independent dense solve differs from the manufactured complex solution.");
    Device deviceMatrix(matrix.size()), deviceRight(n), solution(n);
    check(cudaMemcpyAsync(deviceMatrix.data, matrix.data(), matrix.size() * sizeof(Z),
                          cudaMemcpyHostToDevice, runtime.stream));
    check(cudaMemcpyAsync(deviceRight.data, right.data(), n * sizeof(Z), cudaMemcpyHostToDevice,
                          runtime.stream));
    const Z one{1, 0}, zero{};
    unsigned applications = 0, conditionings = 0, completedIterations = 0, residualApplications = 0,
             lastResidualIteration = 0;
    bool checkRestartActions = true;
    const auto apply = [&](const Z *input, Z *output) {
        ++applications;
        if (checkRestartActions && input == solution.data) {
            require(completedIterations > lastResidualIteration,
                    "Restart recomputed the residual of an unchanged GMRES solution.");
            lastResidualIteration = completedIterations;
            ++residualApplications;
        }
        check(cublasZgemv(runtime.blas, CUBLAS_OP_N, n, n, &one, deviceMatrix.data, n, input, 1,
                          &zero, output, 1));
    };
    const auto condition = [&](Z *values) {
        ++conditionings;
        const Z factor{.7, .1};
        check(cublasZscal(runtime.blas, n, &factor, values, 1));
    };
    hornsculpt::BemNumerics options;
    options.relativeResidualTolerance = 1e-7;
    options.fastMultipole.restartLength = 4;
    options.fastMultipole.maximumIterations = 120;
    hornsculpt::cuda_detail::DeviceGmresSolver solver(apply, condition, n, options, runtime.stream,
                                                      {}, [&] { ++completedIterations; });
    double residual{};
    const auto iterations = solver.solve(deviceRight.data, solution.data, residual);
    check(cudaMemcpyAsync(actual.data(), solution.data, n * sizeof(Z), cudaMemcpyDeviceToHost,
                          runtime.stream));
    check(cudaStreamSynchronize(runtime.stream));
    double error{}, norm{};
    for (unsigned i = 0; i < n; ++i) {
        error += std::norm(C{actual[i].x, actual[i].y} - reference[i]);
        norm += std::norm(reference[i]);
    }
    const double relative = std::sqrt(error / norm);
    std::cout << "Device GMRES restarted complex solve iterations=" << iterations
              << ", residual=" << residual << ", solution error=" << relative
              << ", applications/preconditionings=" << applications << '/' << conditionings
              << ", completed restart cycles=" << residualApplications << '\n';
    require(iterations > 4 && iterations < 120 && relative < 2e-8 && residual < 1e-8 &&
                residualApplications > 1 && completedIterations == iterations &&
                applications == iterations + residualApplications &&
                conditionings == 1 + iterations + residualApplications,
            "Device GMRES failed restarted nonnormal complex equation.");
    checkRestartActions = false;
    check(cudaMemsetAsync(deviceRight.data, 0, n * sizeof(Z), runtime.stream));
    require(solver.solve(deviceRight.data, solution.data, residual) == 0 && residual == 0,
            "Device GMRES zero RHS retained previous solve.");
    options.fastMultipole.maximumIterations = 1;
    hornsculpt::cuda_detail::DeviceGmresSolver bounded(apply, condition, n, options,
                                                       runtime.stream);
    check(cudaMemcpyAsync(deviceRight.data, right.data(), n * sizeof(Z), cudaMemcpyHostToDevice,
                          runtime.stream));
    bool rejected = false;
    try {
        bounded.solve(deviceRight.data, solution.data, residual);
    } catch (const std::runtime_error &error) {
        const std::string message = error.what();
        rejected = message.find("1/1 iterations") != std::string::npos &&
                   message.find("relative FMM equation residual") != std::string::npos;
    }
    require(rejected, "Device GMRES must report actual exhausted iteration budget/residual.");
}
} // namespace
int main() {
    int count{};
    if (cudaGetDeviceCount(&count) != cudaSuccess || !count)
        return 77;
    try {
        checkRestartedSolve();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
