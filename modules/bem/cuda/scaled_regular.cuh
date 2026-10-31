#pragma once

#include "../fmm/wideband_transfer.hpp"
#include <cfloat>
#include <cuComplex.h>
#include <cuda_runtime.h>

namespace hornsculpt::wideband::detail::scaled {
struct ScalarWave {
    cuDoubleComplex value{};
    cuDoubleComplex gradient[3]{};
};

// DLMF 10.51.1 and 10.53.1. The caller supplies order+2 doubles, shared by
// all m columns at this quadrature point. No per-thread radial array is used.
// Jhat_l(z)=(2l+1)!!*j_l(z)/z^l has the finite limit one at z=0.
__device__ inline bool fillScaledRegularRadial(unsigned order, double z, double *values) {
    if (!values || order > maximumScaledSphericalOrder || !isfinite(z) || z < 0 ||
        z > 2 * maximumScaledElectricalRadius)
        return false;
    const unsigned degree = order + 1;
    if (z > 1) {
        const double j0 = sin(z) / z, j1 = (j0 - cos(z)) / z;
        if (z > degree + 32.) {
            values[0] = j0;
            values[1] = j1;
            for (unsigned l = 1; l < degree; ++l)
                values[l + 1] = ((2. * l + 1) / z) * values[l] - values[l - 1];
        } else {
            const unsigned top =
                max(degree + 48u, static_cast<unsigned>(ceil(z + 48 + 8 * cbrt(z))));
            double current = 1, next = 0;
            for (unsigned l = top; l; --l) {
                double previous = ((2. * l + 1) / z) * current - next;
                if (fabs(previous) > 1e200) {
                    previous *= 1e-200;
                    current *= 1e-200;
                    next *= 1e-200;
                    for (unsigned j = l; j <= degree; ++j)
                        values[j] *= 1e-200;
                }
                if (l - 1 <= degree)
                    values[l - 1] = previous;
                next = current;
                current = previous;
            }
            // j0 and j1 never vanish together; choose the better-conditioned
            // normalization rather than dividing by a prescribed Bessel mode.
            const double normalization = fabs(j0) > fabs(j1) ? j0 / values[0] : j1 / values[1];
            for (unsigned l = 0; l <= degree; ++l)
                values[l] *= normalization;
        }
    }
    double leading = 1;
    for (unsigned l = 0; l <= degree; ++l) {
        if (l)
            leading *= z / (2. * l + 1);
        if (z > 1 && z * z > 2 * (2. * l + 3)) {
            values[l] /= leading;
        } else {
            // Here the alternating series terms decrease from the start, so
            // neither cancellation nor underflow of z^l affects normalization.
            double term = 1, sum = 1;
            for (unsigned n = 1; n < 128; ++n) {
                term *= -z * z / (2. * n * (2. * l + 2 * n + 1));
                sum += term;
                if (fabs(term) <= DBL_EPSILON * fabs(sum))
                    break;
            }
            values[l] = sum;
        }
        if (!isfinite(values[l]))
            return false;
    }
    return true;
}

__device__ inline cuDoubleComplex add(cuDoubleComplex a, cuDoubleComplex b) {
    return make_cuDoubleComplex(a.x + b.x, a.y + b.y);
}
__device__ inline cuDoubleComplex multiply(cuDoubleComplex a, cuDoubleComplex b) {
    return make_cuDoubleComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}
__device__ inline cuDoubleComplex times(cuDoubleComplex a, double factor) {
    return make_cuDoubleComplex(a.x * factor, a.y * factor);
}
__device__ inline ScalarWave multiply(const ScalarWave &a, const ScalarWave &b) {
    ScalarWave result;
    result.value = multiply(a.value, b.value);
    for (unsigned component = 0; component < 3; ++component)
        result.gradient[component] =
            add(multiply(a.gradient[component], b.value), multiply(a.value, b.gradient[component]));
    return result;
}
__device__ inline ScalarWave times(ScalarWave wave, double factor) {
    wave.value = times(wave.value, factor);
    for (auto &gradient : wave.gradient)
        gradient = times(gradient, factor);
    return wave;
}

// S_lm=(r/a)^l Y_lm. Gradients are with respect to the physical point, so
// the final radial combination does not introduce a second factor 1/a.
// Polynomial recurrences have finite derivatives at both poles and r=0.
__device__ inline ScalarWave solidDiagonal(unsigned m, double3 point, double radius) {
    const double inverse = 1 / radius;
    ScalarWave xy;
    xy.value = make_cuDoubleComplex(point.x * inverse, point.y * inverse);
    xy.gradient[0] = make_cuDoubleComplex(inverse, 0);
    xy.gradient[1] = make_cuDoubleComplex(0, inverse);
    ScalarWave result;
    result.value = make_cuDoubleComplex(.282094791773878143474039725780386293, 0);
    for (unsigned degree = 1; degree <= m; ++degree)
        result = times(multiply(xy, result), -sqrt((2. * degree + 1) / (2. * degree)));
    return result;
}

// Compute degree l>m from S_(l-2,m), S_(l-1,m). For l=m+1 only current
// (the diagonal) is read. A caller can stream an entire m column in O(1) state.
__device__ inline ScalarWave solidNext(unsigned l, unsigned m, double3 point, double radius,
                                       const ScalarWave &previous, const ScalarWave &current) {
    const double inverse = 1 / radius;
    ScalarWave z;
    z.value = make_cuDoubleComplex(point.z * inverse, 0);
    z.gradient[2] = make_cuDoubleComplex(inverse, 0);
    if (l == m + 1)
        return times(multiply(z, current), sqrt(2. * m + 3));
    ScalarWave r2;
    const double inverseSquared = inverse * inverse;
    r2.value = make_cuDoubleComplex(
        (point.x * point.x + point.y * point.y + point.z * point.z) * inverseSquared, 0);
    r2.gradient[0] = make_cuDoubleComplex(2 * point.x * inverseSquared, 0);
    r2.gradient[1] = make_cuDoubleComplex(2 * point.y * inverseSquared, 0);
    r2.gradient[2] = make_cuDoubleComplex(2 * point.z * inverseSquared, 0);
    const double denominator = static_cast<double>(l) * l - static_cast<double>(m) * m,
                 first = sqrt((4. * l * l - 1) / denominator),
                 second = sqrt((2. * l + 1) * ((l - 1.) * (l - 1.) - static_cast<double>(m) * m) /
                               ((2. * l - 3) * denominator));
    auto result = times(multiply(z, current), first);
    const auto correction = times(multiply(r2, previous), -second);
    result.value = add(result.value, correction.value);
    for (unsigned component = 0; component < 3; ++component)
        result.gradient[component] =
            add(result.gradient[component], correction.gradient[component]);
    return result;
}

// dJhat_l(kr)/dx_i=-k^2*x_i*Jhat_(l+1)/(2l+3).
__device__ inline ScalarWave combineSolidRadial(const ScalarWave &solid, unsigned l, double3 point,
                                                double k, double radial, double nextRadial) {
    ScalarWave result;
    result.value = times(solid.value, radial);
    const double factor = -k * k * nextRadial / (2. * l + 3);
    result.gradient[0] =
        add(times(solid.gradient[0], radial), times(solid.value, factor * point.x));
    result.gradient[1] =
        add(times(solid.gradient[1], radial), times(solid.value, factor * point.y));
    result.gradient[2] =
        add(times(solid.gradient[2], radial), times(solid.value, factor * point.z));
    return result;
}
__device__ inline ScalarWave negativeOrder(ScalarWave wave, unsigned m) {
    const double sign = m % 2 ? -1 : 1;
    wave.value = make_cuDoubleComplex(sign * wave.value.x, -sign * wave.value.y);
    for (auto &gradient : wave.gradient)
        gradient = make_cuDoubleComplex(sign * gradient.x, -sign * gradient.y);
    return wave;
}
} // namespace hornsculpt::wideband::detail::scaled
