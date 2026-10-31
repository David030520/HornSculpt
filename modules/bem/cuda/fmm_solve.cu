#include "fmm_operator.hpp"
#include "fmm_solve.hpp"
#include "workspace_size.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cublas_v2.h>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace hornsculpt::cuda_detail {
namespace {
using C = std::complex<double>;
using Z = cuDoubleComplex;
using Clock = std::chrono::steady_clock;
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
double elapsed(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
std::size_t checkedBasisSize(std::size_t count, unsigned restart) {
    if (count == 0 || count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("GMRES pressure dimension is invalid.");
    auto elements = checkedSizeProduct(count, restart + 7);
    const auto projections = checkedSizeProduct(2, restart);
    if (projections > std::numeric_limits<std::size_t>::max() - elements)
        throw std::overflow_error("GMRES workspace size overflow.");
    elements += projections;
    std::size_t free{}, total{};
    check(cudaMemGetInfo(&free, &total), "GMRES capacity query");
    if (checkedSizeProduct(elements, sizeof(Z)) > deviceMemoryBudget(free))
        throw std::runtime_error("GMRES workspace exceeds available device memory.");
    return checkedSizeProduct(count, restart + 1);
}
} // namespace

class FastMultipoleLinearSolver::Implementation {
  public:
    Implementation(FastMultipoleOperator &op, std::size_t size, const BemNumerics &settings,
                   cudaStream_t selectedStream, std::stop_token token)
        : operation(op), n(static_cast<int>(size)), options(settings), stream(selectedStream),
          stop(token), basis(checkedBasisSize(size, settings.fastMultipole.restartLength)),
          exactRhs(size), residual(size), product(size), correction(size), temporary(size),
          gmresResidual(size), projections(2 * settings.fastMultipole.restartLength) {
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
        operation.precondition(value);
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

    unsigned gmres(const Z *right, Z *solution) {
        // The original equation is audited separately. The tighter Krylov
        // threshold leaves room for interpolation error and FP32 panel storage.
        const auto &settings = options.fastMultipole;
        const unsigned restart = std::min(settings.restartLength, settings.maximumIterations);
        const unsigned stride = restart + 1;
        std::vector<C> h(static_cast<std::size_t>(stride) * restart), g(stride), s(restart);
        std::vector<double> c(restart);
        clear(solution);
        copy(right, temporary.data());
        applyPreconditioner(temporary.data());
        const double rightNorm = norm(temporary.data());
        if (rightNorm == 0)
            return 0;
        const double tolerance = options.relativeResidualTolerance * .1 * rightNorm;
        unsigned iterations = 0;
        while (iterations < settings.maximumIterations) {
            cancelled(stop);
            if (iterations == 0)
                copy(temporary.data(), gmresResidual.data());
            else {
                operation.applyPressure(solution, product.data());
                subtract(right, product.data(), gmresResidual.data());
                applyPreconditioner(gmresResidual.data());
            }
            const double beta = norm(gmresResidual.data());
            if (beta <= tolerance)
                return iterations;
            copy(gmresResidual.data(), basis.data());
            scale(1.0 / beta, basis.data());
            std::fill(h.begin(), h.end(), C{});
            std::fill(g.begin(), g.end(), C{});
            g[0] = beta;
            unsigned used = 0;
            bool breakdown = false;
            for (unsigned j = 0; j < restart && iterations < settings.maximumIterations; ++j) {
                cancelled(stop);
                operation.applyPressure(basis.data() + static_cast<std::size_t>(j) * n,
                                        product.data());
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
            operation.applyPressure(solution, product.data());
            subtract(right, product.data(), gmresResidual.data());
            applyPreconditioner(gmresResidual.data());
            if (norm(gmresResidual.data()) <= tolerance)
                return iterations;
            if (breakdown)
                throw std::runtime_error("FMM GMRES stagnated after Arnoldi breakdown.");
        }
        throw std::runtime_error("FMM GMRES did not converge within " +
                                 std::to_string(settings.maximumIterations) + " iterations.");
    }

    IterativeSolution solve(const cuFloatComplex *flux, Z *pressure, bool useInitialGuess) {
        cancelled(stop);
        IterativeSolution result;
        auto audit = Clock::now();
        operation.applyFluxDirect(flux, exactRhs.data());
        const double rhsNorm = norm(exactRhs.data());
        result.residualCheckMilliseconds += elapsed(audit);
        if (!useInitialGuess)
            clear(pressure);
        const bool hadInitialGuess = useInitialGuess;
        unsigned solves = 0;
        for (;;) {
            cancelled(stop);
            audit = Clock::now();
            if (!useInitialGuess)
                copy(exactRhs.data(), residual.data());
            else {
                operation.applyPressureDirect(pressure, product.data());
                subtract(exactRhs.data(), product.data(), residual.data());
            }
            const double error = norm(residual.data());
            result.boundaryResidual = rhsNorm > 0 ? error / rhsNorm : error;
            result.residualCheckMilliseconds += elapsed(audit);
            if (error == 0 ||
                (rhsNorm > 0 && result.boundaryResidual <= options.relativeResidualTolerance))
                return result;
            if (solves >= options.maximumRefinementSteps + (hadInitialGuess ? 0u : 1u))
                throw std::runtime_error(
                    "FMM original Galerkin residual " + std::to_string(result.boundaryResidual) +
                    " exceeds tolerance " + std::to_string(options.relativeResidualTolerance) +
                    "; increase interpolation order or tighten separation.");
            result.iterations += gmres(residual.data(), correction.data());
            axpy(1.0, correction.data(), pressure);
            useInitialGuess = true;
            ++solves;
            result.refinementSteps = solves - (hadInitialGuess ? 0u : 1u);
        }
    }
    std::size_t bytes() const {
        return basis.bytes() + exactRhs.bytes() + residual.bytes() + product.bytes() +
               correction.bytes() + temporary.bytes() + gmresResidual.bytes() + projections.bytes();
    }
    FastMultipoleOperator &operation;
    int n;
    BemNumerics options;
    cudaStream_t stream;
    std::stop_token stop;
    Buffer basis, exactRhs, residual, product, correction, temporary, gmresResidual, projections;
    cublasHandle_t blas{};
};

FastMultipoleLinearSolver::FastMultipoleLinearSolver(FastMultipoleOperator &operation,
                                                     std::size_t count, const BemNumerics &options,
                                                     cudaStream_t stream, std::stop_token stop)
    : implementation_(std::make_unique<Implementation>(operation, count, options, stream, stop)) {}
FastMultipoleLinearSolver::~FastMultipoleLinearSolver() = default;
IterativeSolution FastMultipoleLinearSolver::solve(const cuFloatComplex *flux,
                                                   cuDoubleComplex *pressure, bool initial) {
    return implementation_->solve(flux, pressure, initial);
}
std::size_t FastMultipoleLinearSolver::deviceBytes() const {
    return implementation_->bytes();
}
} // namespace hornsculpt::cuda_detail
