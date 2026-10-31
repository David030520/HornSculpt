#include "cuda/scaled_regular.cuh"
#include "cuda/wideband_runtime.cuh"
#include "fmm/scaled_translation.hpp"
#include "fmm/spherical_rotation.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
using namespace hornsculpt;
using namespace hornsculpt::wideband;
using namespace hornsculpt::wideband::detail;
namespace regular = hornsculpt::wideband::detail::scaled;

// The same m-column/shared-radial layout is used by leaf emission and target
// testing. No thread owns an O(p^2) wave array or O(p) radial scratch.
__global__ void evaluateRegular(unsigned order, double k, double radius, double3 point,
                                regular::ScalarWave *output, int *error) {
    extern __shared__ double radial[];
    if (!threadIdx.x &&
        !regular::fillScaledRegularRadial(
            order, k * sqrt(point.x * point.x + point.y * point.y + point.z * point.z), radial))
        atomicExch(error, 1);
    __syncthreads();
    for (unsigned m = threadIdx.x; m <= order; m += blockDim.x) {
        regular::ScalarWave previous{}, current = regular::solidDiagonal(m, point, radius);
        for (unsigned l = m; l <= order; ++l) {
            if (l > m) {
                const auto next = regular::solidNext(l, m, point, radius, previous, current);
                previous = current;
                current = next;
            }
            const auto value =
                regular::combineSolidRadial(current, l, point, k, radial[l], radial[l + 1]);
            output[static_cast<std::size_t>(l) * l + l + m] = value;
            if (m)
                output[static_cast<std::size_t>(l) * l + l - m] = regular::negativeOrder(value, m);
        }
    }
}
Complex complex(cuDoubleComplex value) {
    return {value.x, value.y};
}
void compareRegular(unsigned order, double argument, Vec3 point, double radius = .1) {
    const double length = std::hypot(point.x, point.y, point.z),
                 k = length ? argument / length : 100;
    const auto count = static_cast<std::size_t>(order + 1) * (order + 1);
    DeviceBuffer<regular::ScalarWave> output;
    DeviceBuffer<int> status;
    Stream stream;
    output.allocate(count);
    status.allocate(1);
    checked(cudaMemsetAsync(status.data(), 0, sizeof(int), stream), "Regular test status");
    evaluateRegular<<<1, 128, (order + 2) * sizeof(double), stream>>>(
        order, k, radius, make_double3(point.x, point.y, point.z), output.data(), status.data());
    checked(cudaGetLastError(), "Scaled regular test launch");
    const auto actual = output.download(stream);
    if (status.download(stream).front())
        throw std::runtime_error("CUDA regular radial helper rejected a valid point.");
    const auto expected = scaledRegularWaves(order, k, radius, point);
    double error = 0, norm = 0, gradientError = 0, gradientNorm = 0;
    for (std::size_t i = 0; i < count; ++i) {
        error += std::norm(complex(actual[i].value) - expected.waves[i].value);
        norm += std::norm(expected.waves[i].value);
        for (unsigned component = 0; component < 3; ++component) {
            gradientError += std::norm(complex(actual[i].gradient[component]) -
                                       expected.waves[i].gradient[component]);
            gradientNorm += std::norm(expected.waves[i].gradient[component]);
        }
    }
    // Independent derivative check: values at displaced physical coordinates,
    // rather than another analytical derivative of the same recurrence.
    double finiteError = 0, finiteNorm = 0;
    const double step = radius * 1e-6;
    for (unsigned component = 0; component < 3; ++component) {
        Vec3 plus = point, minus = point;
        if (component == 0) {
            plus.x += step;
            minus.x -= step;
        }
        if (component == 1) {
            plus.y += step;
            minus.y -= step;
        }
        if (component == 2) {
            plus.z += step;
            minus.z -= step;
        }
        const auto upper = scaledRegularWaves(order, k, radius, plus),
                   lower = scaledRegularWaves(order, k, radius, minus);
        for (std::size_t i = 0; i < count; ++i) {
            const auto derivative = (upper.waves[i].value - lower.waves[i].value) / (2 * step);
            finiteError += std::norm(complex(actual[i].gradient[component]) - derivative);
            finiteNorm += std::norm(derivative);
        }
    }
    const double relative = std::sqrt(error / norm),
                 relativeGradient = std::sqrt(gradientError / gradientNorm),
                 relativeFinite = std::sqrt(finiteError / finiteNorm);
    std::cout << "CUDA scaled regular p=" << order << ", kr=" << argument << ", value=" << relative
              << ", gradient=" << relativeGradient << ", finite difference=" << relativeFinite
              << '\n';
    if (relative > 3e-12 || relativeGradient > 3e-12 || relativeFinite > 3e-6)
        throw std::runtime_error("CUDA scaled regular values/physical gradients failed.");
}
} // namespace

void checkScaledRegularCuda() {
    compareRegular(52, 0, {0, 0, 0});
    compareRegular(36, .01, {.017, -.031, .024});
    compareRegular(52, 1, {0, 0, .071});
    compareRegular(52, 3.141592653589793, {0, 0, -.071});
    compareRegular(52, 32, {.037, -.051, .023});
    compareRegular(64, 64, {-.027, .041, .031});
    compareRegular(maximumDeviceScaledOrder, 128, {.125, 0, 0}, .125);
}
