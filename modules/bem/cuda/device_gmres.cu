#include "device_gmres.hpp"
#include "workspace_size.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cublas_v2.h>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
namespace hornsculpt::cuda_detail {
namespace {
using C = std::complex<double>;
using Z = cuDoubleComplex;
void check(cudaError_t code, const char *context) {
    if (code != cudaSuccess)
        throw std::runtime_error(std::string(context) + ": " + cudaGetErrorString(code));
}
void check(cublasStatus_t code, const char *context) {
    if (code != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error(std::string(context) + ": cuBLAS " + std::to_string(code));
}
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("BEM calculation cancelled.");
}
Z device(C value) {
    return make_cuDoubleComplex(value.real(), value.imag());
}
C host(Z value) {
    return {value.x, value.y};
}
class Buffer {
  public:
    explicit Buffer(std::size_t count) : count_(count) {
        check(cudaMalloc(reinterpret_cast<void **>(&data_), checkedSizeProduct(count, sizeof(Z))),
              "GMRES allocation");
    }
    ~Buffer() {
        cudaFree(data_);
    }
    Buffer(const Buffer &) = delete;
    Buffer &operator=(const Buffer &) = delete;
    Z *data() {
        return data_;
    }
    std::size_t bytes() const {
        return count_ * sizeof(Z);
    }

  private:
    Z *data_{};
    std::size_t count_{};
};
__global__ void difference(const Z *right, const Z *left, Z *output, int count) {
    const int i = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (i < count)
        output[i] = cuCsub(right[i], left[i]);
}
std::string scientific(double value) {
    std::ostringstream output;
    output << std::scientific << std::setprecision(6) << value;
    return output.str();
}
std::size_t checkedBasisSize(std::size_t count, unsigned restart) {
    if (count == 0 || count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("GMRES pressure dimension is invalid.");
    auto elements = checkedSizeProduct(count, static_cast<std::size_t>(restart) + 4);
    const auto projections = checkedSizeProduct(2, restart);
    if (projections > std::numeric_limits<std::size_t>::max() - elements)
        throw std::overflow_error("GMRES workspace size overflow.");
    elements += projections;
    std::size_t free{}, total{};
    check(cudaMemGetInfo(&free, &total), "GMRES capacity query");
    if (checkedSizeProduct(elements, sizeof(Z)) > deviceMemoryBudget(free))
        throw std::runtime_error("GMRES workspace exceeds available device memory.");
    return checkedSizeProduct(count, static_cast<std::size_t>(restart) + 1);
}
} // namespace

class DeviceGmresSolver::Implementation {
  public:
    Implementation(PressureApplication pressure, Preconditioner conditioning, std::size_t size,
                   const BemNumerics &settings, cudaStream_t selectedStream, std::stop_token token,
                   IterationProgress callback)
        : applyPressure(std::move(pressure)), precondition(std::move(conditioning)),
          n(static_cast<int>(size)), options(settings), stream(selectedStream), stop(token),
          progress(std::move(callback)),
          basis(checkedBasisSize(size, settings.fastMultipole.restartLength)), product(size),
          temporary(size), gmresResidual(size),
          projections(checkedSizeProduct(2, settings.fastMultipole.restartLength)) {
        if (!applyPressure || !precondition || !settings.fastMultipole.maximumIterations ||
            !settings.fastMultipole.restartLength ||
            !std::isfinite(settings.relativeResidualTolerance) ||
            settings.relativeResidualTolerance <= 0)
            throw std::invalid_argument(
                "Device GMRES settings or operation callbacks are invalid.");
        check(cublasCreate(&blas), "GMRES cuBLAS creation");
        try {
            check(cublasSetStream(blas, stream), "GMRES cuBLAS stream");
        } catch (...) {
            cublasDestroy(blas);
            throw;
        }
    }
    ~Implementation() {
        cublasDestroy(blas);
    }
    void synchronize() {
        check(cudaStreamSynchronize(stream), "GMRES synchronization");
    }
    void copy(const Z *source, Z *destination) {
        check(cudaMemcpyAsync(destination, source, static_cast<std::size_t>(n) * sizeof(Z),
                              cudaMemcpyDeviceToDevice, stream),
              "GMRES vector copy");
    }
    void clear(Z *value) {
        check(cudaMemsetAsync(value, 0, static_cast<std::size_t>(n) * sizeof(Z), stream),
              "GMRES vector clear");
    }
    void axpy(C coefficient, const Z *input, Z *output) {
        const Z alpha = device(coefficient);
        check(cublasZaxpy(blas, n, &alpha, input, 1, output, 1), "GMRES vector update");
    }
    void scale(C coefficient, Z *value) {
        const Z alpha = device(coefficient);
        check(cublasZscal(blas, n, &alpha, value, 1), "GMRES normalization");
    }
    double norm(const Z *value) {
        double result{};
        check(cublasDznrm2(blas, n, value, 1, &result), "GMRES vector norm");
        synchronize();
        if (!std::isfinite(result))
            throw std::runtime_error("GMRES encountered a non-finite vector.");
        return result;
    }
    void applyPreconditioner(Z *value) {
        precondition(value);
    }
    void orthogonalize(unsigned columns, C *h) {
        const Z one = device(1), zero = device(0), minusOne = device(-1);
        // Twice classical Gram-Schmidt keeps both projection passes on the GPU.
        // Reorthogonalization controls cancellation without one host transfer
        // for every existing basis vector in the Arnoldi sequence.
        for (unsigned pass = 0; pass < 2; ++pass) {
            auto *coefficients = projections.data() + pass * columns;
            check(cublasZgemv(blas, CUBLAS_OP_C, n, static_cast<int>(columns), &one, basis.data(),
                              n, product.data(), 1, &zero, coefficients, 1),
                  "GMRES basis projection");
            check(cublasZgemv(blas, CUBLAS_OP_N, n, static_cast<int>(columns), &minusOne,
                              basis.data(), n, coefficients, 1, &one, product.data(), 1),
                  "GMRES reorthogonalization");
        }
        std::vector<Z> values(2 * columns);
        check(cudaMemcpyAsync(values.data(), projections.data(), values.size() * sizeof(Z),
                              cudaMemcpyDeviceToHost, stream),
              "GMRES projection download");
        synchronize();
        for (unsigned i = 0; i < columns; ++i)
            h[i] = host(values[i]) + host(values[columns + i]);
    }
    void subtract(const Z *right, const Z *left, Z *output) {
        difference<<<(n + 255) / 256, 256, 0, stream>>>(right, left, output, n);
        check(cudaGetLastError(), "GMRES residual subtraction");
    }

    unsigned gmres(const Z *right, Z *solution, double &relativeResidual) {
        // A higher-order equation is audited separately. The tighter Krylov
        // threshold leaves room for expansion and mixed-sample arithmetic
        // errors; original panel geometry and phases remain FP64.
        const auto &settings = options.fastMultipole;
        const unsigned restart = std::min(settings.restartLength, settings.maximumIterations);
        const unsigned stride = restart + 1;
        std::vector<C> h(static_cast<std::size_t>(stride) * restart), g(stride), s(restart);
        std::vector<double> c(restart);
        clear(solution);
        copy(right, temporary.data());
        const double unpreconditionedRightNorm = norm(right);
        applyPreconditioner(temporary.data());
        const double rightNorm = norm(temporary.data());
        if (rightNorm == 0) {
            if (unpreconditionedRightNorm != 0)
                throw std::runtime_error("GMRES preconditioner annihilated a nonzero RHS.");
            relativeResidual = 0;
            return 0;
        }
        const double tolerance =
            options.relativeResidualTolerance * deviceGmresResidualFraction * rightNorm;
        copy(temporary.data(), gmresResidual.data());
        double beta = rightNorm;
        double achievedPreconditionedResidual = 1;
        unsigned iterations = 0;
        while (iterations < settings.maximumIterations) {
            cancelled(stop);
            if (beta <= tolerance) {
                applyPressure(solution, product.data());
                subtract(right, product.data(), gmresResidual.data());
                relativeResidual = norm(gmresResidual.data()) / unpreconditionedRightNorm;
                return iterations;
            }
            copy(gmresResidual.data(), basis.data());
            scale(1.0 / beta, basis.data());
            std::fill(h.begin(), h.end(), C{});
            std::fill(g.begin(), g.end(), C{});
            g[0] = beta;
            unsigned used = 0;
            bool breakdown = false;
            for (unsigned j = 0; j < restart && iterations < settings.maximumIterations; ++j) {
                cancelled(stop);
                applyPressure(basis.data() + static_cast<std::size_t>(j) * n, product.data());
                applyPreconditioner(product.data());
                const double before = norm(product.data());
                orthogonalize(j + 1, h.data() + static_cast<std::size_t>(j) * stride);
                const double next = norm(product.data());
                h[static_cast<std::size_t>(j) * stride + j + 1] = next;
                breakdown = next <= 64 * std::numeric_limits<double>::epsilon() * before;
                if (!breakdown) {
                    Z *v = basis.data() + static_cast<std::size_t>(j + 1) * n;
                    copy(product.data(), v);
                    scale(1.0 / next, v);
                }
                for (unsigned i = 0; i < j; ++i) {
                    auto &a = h[static_cast<std::size_t>(j) * stride + i];
                    auto &b = h[static_cast<std::size_t>(j) * stride + i + 1];
                    const C rotated = c[i] * a + s[i] * b;
                    b = -std::conj(s[i]) * a + c[i] * b;
                    a = rotated;
                }
                auto &a = h[static_cast<std::size_t>(j) * stride + j];
                auto &b = h[static_cast<std::size_t>(j) * stride + j + 1];
                const double length = std::hypot(std::abs(a), std::abs(b));
                if (!(length > 0) || !std::isfinite(length))
                    throw std::runtime_error("FMM GMRES Arnoldi factorization broke down.");
                const C phase = std::abs(a) > 0 ? a / std::abs(a) : C{1, 0};
                c[j] = std::abs(a) / length;
                s[j] = phase * std::conj(b) / length;
                a = phase * length;
                b = 0;
                const C rotated = c[j] * g[j] + s[j] * g[j + 1];
                g[j + 1] = -std::conj(s[j]) * g[j] + c[j] * g[j + 1];
                g[j] = rotated;
                ++iterations;
                if (progress)
                    progress();
                used = j + 1;
                if (std::abs(g[j + 1]) <= tolerance || breakdown)
                    break;
            }
            std::vector<C> y(used);
            for (unsigned remaining = used; remaining > 0; --remaining) {
                const unsigned row = remaining - 1;
                C value = g[row];
                for (unsigned column = row + 1; column < used; ++column)
                    value -= h[static_cast<std::size_t>(column) * stride + row] * y[column];
                y[row] = value / h[static_cast<std::size_t>(row) * stride + row];
                if (!std::isfinite(y[row].real()) || !std::isfinite(y[row].imag()))
                    throw std::runtime_error("FMM GMRES triangular solve is not finite.");
            }
            for (unsigned j = 0; j < used; ++j)
                axpy(y[j], basis.data() + static_cast<std::size_t>(j) * n, solution);
            applyPressure(solution, product.data());
            subtract(right, product.data(), gmresResidual.data());
            relativeResidual = norm(gmresResidual.data()) / unpreconditionedRightNorm;
            applyPreconditioner(gmresResidual.data());
            const double preconditionedResidualNorm = norm(gmresResidual.data());
            achievedPreconditionedResidual = preconditionedResidualNorm / rightNorm;
            if (preconditionedResidualNorm <= tolerance)
                return iterations;
            if (breakdown)
                throw std::runtime_error("FMM GMRES stagnated after Arnoldi breakdown.");
            // The verified restart residual is already M^-1(b-A*x). Keep
            // its vector and norm until the next Arnoldi cycle; recomputing
            // them would apply the unchanged operator and PC a second time.
            beta = preconditionedResidualNorm;
        }
        throw std::runtime_error(
            "FMM GMRES did not converge after " + std::to_string(iterations) + "/" +
            std::to_string(settings.maximumIterations) +
            " iterations: unpreconditioned relative FMM equation residual " +
            scientific(relativeResidual) + ", preconditioned relative residual " +
            scientific(achievedPreconditionedResidual) + " (inner target " +
            scientific(options.relativeResidualTolerance * deviceGmresResidualFraction) + ").");
    }

    std::size_t bytes() const {
        return basis.bytes() + product.bytes() + temporary.bytes() + gmresResidual.bytes() +
               projections.bytes();
    }
    PressureApplication applyPressure;
    Preconditioner precondition;
    int n;
    BemNumerics options;
    cudaStream_t stream;
    std::stop_token stop;
    IterationProgress progress;
    Buffer basis, product, temporary, gmresResidual, projections;
    cublasHandle_t blas{};
};
DeviceGmresSolver::DeviceGmresSolver(PressureApplication apply, Preconditioner precondition,
                                     std::size_t count, const BemNumerics &options,
                                     cudaStream_t stream, std::stop_token stop,
                                     IterationProgress progress)
    : implementation_(std::make_unique<Implementation>(std::move(apply), std::move(precondition),
                                                       count, options, stream, stop,
                                                       std::move(progress))) {}
DeviceGmresSolver::~DeviceGmresSolver() = default;
unsigned DeviceGmresSolver::solve(const cuDoubleComplex *right, cuDoubleComplex *solution,
                                  double &relativeResidual) {
    if (!right || !solution || right == solution)
        throw std::invalid_argument(
            "Device GMRES requires distinct nonnull RHS and solution arrays.");
    return implementation_->gmres(right, solution, relativeResidual);
}
std::size_t DeviceGmresSolver::deviceBytes() const noexcept {
    return implementation_->bytes();
}
} // namespace hornsculpt::cuda_detail
