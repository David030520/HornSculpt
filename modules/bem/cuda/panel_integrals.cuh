#pragma once

#include "../galerkin/near_projection.hpp"
#include "../galerkin/pair_quadrature.hpp"
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <hornsculpt/mesh.hpp>

namespace hornsculpt::cuda_detail {
struct Panel {
    Vec3 vertices[3], centroid, normal, curls[3];
    unsigned nodes[3];
    double area, edge;
};
struct PairJob {
    unsigned test{}, trial{};
    galerkin::PairPermutation permutation;
};
struct NearJob {
    unsigned test{}, trial{}, firstPatch{}, patchCount{}, swapBarycentric{};
};
struct NearSample {
    double basis[3];
    Vec3 local;
};
struct DoubleComplex {
    double real{}, imag{};
};
__host__ __device__ inline Vec3 sub(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
__host__ __device__ inline Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
__host__ __device__ inline Vec3 mul(Vec3 a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}
__host__ __device__ inline double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
__host__ __device__ inline Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
__host__ __device__ inline double length(Vec3 a) {
    return sqrt(dot(a, a));
}
__device__ inline DoubleComplex plus(DoubleComplex a, DoubleComplex b) {
    return {a.real + b.real, a.imag + b.imag};
}
__device__ inline DoubleComplex times(DoubleComplex a, double s) {
    return {a.real * s, a.imag * s};
}
__device__ inline DoubleComplex product(DoubleComplex a, DoubleComplex b) {
    return {a.real * b.real - a.imag * b.imag, a.real * b.imag + a.imag * b.real};
}
__device__ inline double magnitude(DoubleComplex a) {
    return hypot(a.real, a.imag);
}
__device__ inline DoubleComplex promoted(cuFloatComplex a) {
    return {a.x, a.y};
}
__device__ inline DoubleComplex promoted(cuDoubleComplex a) {
    return {a.x, a.y};
}
__device__ inline Vec3 localPoint(const Panel &panel, const double *bary) {
    return add(mul(sub(panel.vertices[1], panel.vertices[0]), bary[1]),
               mul(sub(panel.vertices[2], panel.vertices[0]), bary[2]));
}
__device__ inline Vec3 pairDelta(const Panel &test, const Panel &trial, const double *testBasis,
                                 const double *trialBasis) {
    return add(sub(test.vertices[0], trial.vertices[0]),
               sub(localPoint(test, testBasis), localPoint(trial, trialBasis)));
}
struct PairMoments {
    DoubleComplex matrix[9], rhs[3], green;
};
struct MixedMoments {
    cuFloatComplex matrix[9]{}, rhs[3]{}, green{};
};
template <class T> __device__ inline float floatDot(T a, T b) {
    return static_cast<float>(a.x) * static_cast<float>(b.x) +
           static_cast<float>(a.y) * static_cast<float>(b.y) +
           static_cast<float>(a.z) * static_cast<float>(b.z);
}
struct MixedKernelSample {
    cuFloatComplex green{}, pressure{}, flux{};
};

// Evaluate scalar kernels before contracting bases. Dense assembly and the
// tensor-contracted near path thereby share one distance and phase contract.
template <class PanelType>
__device__ inline MixedKernelSample
mixedKernelSample(Vec3 delta, const PanelType &test, const PanelType &trial, float physicalWeight,
                  double waveNumber, float inverseK, bool coincident, int *error) {
    // Preserve the original FP64 displacement and wavenumber for the phase.
    // Above the measured small-phase domain, rounding k or r to FP32 causes
    // an O(epsilon_float*k*r) error that no quadrature refinement can remove.
    const float k = static_cast<float>(waveNumber);
    const float dx = static_cast<float>(delta.x), dy = static_cast<float>(delta.y),
                dz = static_cast<float>(delta.z), radiusSquared = dx * dx + dy * dy + dz * dz;
    if (!(radiusSquared > 0) || !isfinite(radiusSquared)) {
        atomicExch(error, 1);
        return {};
    }
    const float inverseRadius = rsqrtf(radiusSquared), radius = radiusSquared * inverseRadius,
                inverseRadiusSquared = inverseRadius * inverseRadius;
    float phase = k * radius, sine, cosine;
    if (phase <= 1.f) {
        sincosf(phase, &sine, &cosine);
    } else {
        const double radius64 = sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z),
                     phase64 = waveNumber * radius64;
        if (!(radius64 > 0) || !isfinite(phase64)) {
            atomicExch(error, 1);
            return {};
        }
        double sine64, cosine64;
        sincos(phase64, &sine64, &cosine64);
        sine = static_cast<float>(sine64);
        cosine = static_cast<float>(cosine64);
        phase = static_cast<float>(phase64);
    }
    const float scale = physicalWeight * (0.07957747154594766788f * inverseRadius);
    const float gr = cosine * scale, gi = -sine * scale;
    const float radialReal = gr - gi * phase, radialImag = gi + gr * phase;
    const float sourceFactor = coincident ? 0.f
                                          : (dx * static_cast<float>(trial.normal.x) +
                                             dy * static_cast<float>(trial.normal.y) +
                                             dz * static_cast<float>(trial.normal.z)) *
                                                inverseRadiusSquared,
                targetFactor = coincident ? 0.f
                                          : -(dx * static_cast<float>(test.normal.x) +
                                              dy * static_cast<float>(test.normal.y) +
                                              dz * static_cast<float>(test.normal.z)) *
                                                inverseRadiusSquared,
                normalProduct = floatDot(test.normal, trial.normal);
    const float ar = -radialReal * sourceFactor - k * normalProduct * gi,
                ai = -radialImag * sourceFactor + k * normalProduct * gr,
                br = -gr - radialImag * targetFactor * inverseK,
                bi = -gi + radialReal * targetFactor * inverseK;
    return {{gr, gi}, {ar, ai}, {br, bi}};
}

template <class PanelType, class Basis>
__device__ inline void mixedSample(MixedMoments &result, Vec3 delta, const PanelType &test,
                                   const PanelType &trial, const Basis *testBasis,
                                   const Basis *trialBasis, float physicalWeight, double waveNumber,
                                   float inverseK, bool coincident, int *error) {
    const auto sample = mixedKernelSample(delta, test, trial, physicalWeight, waveNumber, inverseK,
                                          coincident, error);
    result.green.x += sample.green.x;
    result.green.y += sample.green.y;
#pragma unroll
    for (unsigned i = 0; i < 3; ++i) {
        const float basis = static_cast<float>(testBasis[i]);
        result.rhs[i].x += sample.flux.x * basis;
        result.rhs[i].y += sample.flux.y * basis;
#pragma unroll
        for (unsigned j = 0; j < 3; ++j) {
            const float product = basis * static_cast<float>(trialBasis[j]);
            result.matrix[3 * i + j].x += sample.pressure.x * product;
            result.matrix[3 * i + j].y += sample.pressure.y * product;
        }
    }
}
__device__ inline void pairSample(PairMoments &result, const Panel &test, const Panel &trial,
                                  const galerkin::PairQuadraturePoint &sample, Vec3 displacement,
                                  double k, bool coincident, int *error) {
    const double radius = length(displacement);
    if (!(radius > 0.0) || !isfinite(radius)) {
        atomicExch(error, 1);
        return;
    }
    double sine, cosine;
    sincos(k * radius, &sine, &cosine);
    const double scale =
        sample.weight * (4.0 * test.area * trial.area) / (12.56637061435917295385 * radius);
    const DoubleComplex green{cosine * scale, -sine * scale};
    result.green = plus(result.green, green);
    const DoubleComplex radial{green.real - green.imag * k * radius,
                               green.imag + green.real * k * radius};
    const double inverseRadiusSquared = 1.0 / (radius * radius);
    const auto d =
        times(radial, coincident ? 0.0 : dot(displacement, trial.normal) * inverseRadiusSquared);
    const auto adjoint =
        times(radial, coincident ? 0.0 : -dot(displacement, test.normal) * inverseRadiusSquared);
    const double normalProduct = dot(test.normal, trial.normal);
    // alpha=-i/k; combine D and the k^2 part of H before gathering the curl term.
    const DoubleComplex a{-d.real - k * normalProduct * green.imag,
                          -d.imag + k * normalProduct * green.real};
    const DoubleComplex b{-green.real - adjoint.imag / k, -green.imag + adjoint.real / k};
    for (unsigned i = 0; i < 3; ++i) {
        result.rhs[i] = plus(result.rhs[i], times(b, sample.test[i]));
        for (unsigned j = 0; j < 3; ++j)
            result.matrix[3 * i + j] =
                plus(result.matrix[3 * i + j], times(a, sample.test[i] * sample.trial[j]));
    }
}
__device__ inline void atomicComplex(cuFloatComplex *output, DoubleComplex value, int *error) {
    const float real = static_cast<float>(value.real), imag = static_cast<float>(value.imag);
    if (!isfinite(real) || !isfinite(imag)) {
        atomicExch(error, 2);
        return;
    }
    atomicAdd(&output->x, real);
    atomicAdd(&output->y, imag);
}
__device__ inline void atomicMixed(cuFloatComplex *output, cuFloatComplex value, int *error) {
    if (!isfinite(value.x) || !isfinite(value.y)) {
        atomicExch(error, 2);
        return;
    }
    atomicAdd(&output->x, value.x);
    atomicAdd(&output->y, value.y);
}
template <class PanelType>
__device__ inline void scatterMixed(const MixedMoments &moments, const PanelType &test,
                                    const PanelType &trial, int pressureCount, int fluxCount,
                                    int trialIndex, int excitationCount, float k,
                                    const double *rowScale, const cuFloatComplex *flux,
                                    cuFloatComplex *matrix, cuFloatComplex *rhs, int *error) {
    const float inverseK = 1.f / k;
    for (unsigned i = 0; i < 3; ++i) {
        const unsigned row = test.nodes[i];
        const float scale = static_cast<float>(rowScale[row]);
        for (unsigned j = 0; j < 3; ++j) {
            const float curl = floatDot(test.curls[i], trial.curls[j]) * inverseK;
            atomicMixed(
                matrix + row + static_cast<std::size_t>(trial.nodes[j]) * pressureCount,
                make_cuFloatComplex((moments.matrix[3 * i + j].x + moments.green.y * curl) * scale,
                                    (moments.matrix[3 * i + j].y - moments.green.x * curl) * scale),
                error);
        }
        for (int excitation = 0; excitation < excitationCount; ++excitation) {
            const auto q = flux[static_cast<std::size_t>(excitation) * fluxCount + trialIndex];
            if (q.x == 0 && q.y == 0)
                continue;
            const auto b = moments.rhs[i];
            atomicMixed(rhs + row + static_cast<std::size_t>(excitation) * pressureCount,
                        make_cuFloatComplex((b.x * q.x - b.y * q.y) * scale,
                                            (b.x * q.y + b.y * q.x) * scale),
                        error);
        }
    }
}
__device__ inline void reduceWarp(MixedMoments &value) {
    for (int offset = 16; offset > 0; offset /= 2) {
        for (unsigned i = 0; i < 9; ++i) {
            value.matrix[i].x += __shfl_down_sync(0xffffffffU, value.matrix[i].x, offset);
            value.matrix[i].y += __shfl_down_sync(0xffffffffU, value.matrix[i].y, offset);
        }
        for (unsigned i = 0; i < 3; ++i) {
            value.rhs[i].x += __shfl_down_sync(0xffffffffU, value.rhs[i].x, offset);
            value.rhs[i].y += __shfl_down_sync(0xffffffffU, value.rhs[i].y, offset);
        }
        value.green.x += __shfl_down_sync(0xffffffffU, value.green.x, offset);
        value.green.y += __shfl_down_sync(0xffffffffU, value.green.y, offset);
    }
}
__device__ inline void scatterPair(const PairMoments &moments, const Panel &test,
                                   const Panel &trial, int pressureCount, int fluxCount,
                                   int trialIndex, int excitationCount, double k,
                                   const double *rowScale, const cuFloatComplex *flux,
                                   cuFloatComplex *matrix, cuFloatComplex *rhs, int *error) {
    for (unsigned i = 0; i < 3; ++i) {
        const unsigned row = test.nodes[i];
        for (unsigned j = 0; j < 3; ++j) {
            const double curlProduct = dot(test.curls[i], trial.curls[j]);
            const auto value = plus(moments.matrix[3 * i + j],
                                    DoubleComplex{moments.green.imag * curlProduct / k,
                                                  -moments.green.real * curlProduct / k});
            atomicComplex(matrix + row + static_cast<std::size_t>(trial.nodes[j]) * pressureCount,
                          times(value, rowScale[row]), error);
        }
        for (int excitation = 0; excitation < excitationCount; ++excitation) {
            const auto q = flux[static_cast<std::size_t>(excitation) * fluxCount + trialIndex];
            if (q.x == 0 && q.y == 0)
                continue;
            atomicComplex(rhs + row + static_cast<std::size_t>(excitation) * pressureCount,
                          times(product(moments.rhs[i], promoted(q)), rowScale[row]), error);
        }
    }
}
__device__ inline void reduceWarp(PairMoments &value) {
    for (int offset = 16; offset > 0; offset /= 2) {
        for (unsigned i = 0; i < 9; ++i) {
            value.matrix[i].real += __shfl_down_sync(0xffffffffU, value.matrix[i].real, offset);
            value.matrix[i].imag += __shfl_down_sync(0xffffffffU, value.matrix[i].imag, offset);
        }
        for (unsigned i = 0; i < 3; ++i) {
            value.rhs[i].real += __shfl_down_sync(0xffffffffU, value.rhs[i].real, offset);
            value.rhs[i].imag += __shfl_down_sync(0xffffffffU, value.rhs[i].imag, offset);
        }
        value.green.real += __shfl_down_sync(0xffffffffU, value.green.real, offset);
        value.green.imag += __shfl_down_sync(0xffffffffU, value.green.imag, offset);
    }
}
struct FieldPatch {
    Vec3 bary[3];
    double fraction;
    unsigned depth;
};
struct FieldIntegral {
    DoubleComplex value;
    double absoluteScale{};
};
__device__ inline FieldPatch childPatch(const FieldPatch &patch, unsigned child) {
    const auto ab = mul(add(patch.bary[0], patch.bary[1]), .5);
    const auto ac = mul(add(patch.bary[0], patch.bary[2]), .5);
    const auto bc = mul(add(patch.bary[1], patch.bary[2]), .5);
    FieldPatch result{};
    result.fraction = patch.fraction * .25;
    result.depth = patch.depth + 1;
    if (child == 0) {
        result.bary[0] = patch.bary[0];
        result.bary[1] = ab;
        result.bary[2] = ac;
    } else if (child == 1) {
        result.bary[0] = ab;
        result.bary[1] = patch.bary[1];
        result.bary[2] = bc;
    } else if (child == 2) {
        result.bary[0] = ac;
        result.bary[1] = bc;
        result.bary[2] = patch.bary[2];
    } else {
        result.bary[0] = ab;
        result.bary[1] = bc;
        result.bary[2] = ac;
    }
    return result;
}
__device__ inline FieldIntegral fieldPatch(Vec3 observer, const Panel &panel,
                                           const FieldPatch &patch,
                                           const galerkin::TriangleQuadraturePoint *rule, int count,
                                           const cuDoubleComplex *pressure, cuFloatComplex flux,
                                           double k, int *error) {
    FieldIntegral result{};
    for (int sample = 0; sample < count; ++sample) {
        const auto lambda = add(add(mul(patch.bary[0], rule[sample].barycentric[0]),
                                    mul(patch.bary[1], rule[sample].barycentric[1])),
                                mul(patch.bary[2], rule[sample].barycentric[2]));
        const double bary[3]{lambda.x, lambda.y, lambda.z};
        const auto displacement = sub(sub(observer, panel.vertices[0]), localPoint(panel, bary));
        const double radius = length(displacement);
        if (!(radius > 0.0)) {
            atomicExch(error, 1);
            return {};
        }
        double sine, cosine;
        sincos(k * radius, &sine, &cosine);
        const double scale = 2.0 * panel.area * patch.fraction * rule[sample].weight /
                             (12.56637061435917295385 * radius);
        const DoubleComplex green{cosine * scale, -sine * scale};
        const double normalFactor = dot(displacement, panel.normal) / (radius * radius);
        const DoubleComplex dual{(green.real - green.imag * k * radius) * normalFactor,
                                 (green.imag + green.real * k * radius) * normalFactor};
        DoubleComplex p{};
        for (unsigned local = 0; local < 3; ++local)
            p = plus(p, times(promoted(pressure[panel.nodes[local]]), bary[local]));
        const auto dp = product(dual, p), sq = product(green, promoted(flux));
        result.value = plus(result.value, {dp.real - sq.real, dp.imag - sq.imag});
        result.absoluteScale += magnitude(dp) + magnitude(sq);
    }
    return result;
}
__device__ inline DoubleComplex regularField(Vec3 observer, const Panel &panel,
                                             const galerkin::TriangleQuadraturePoint *rule,
                                             unsigned count, const cuDoubleComplex *pressure,
                                             cuFloatComplex flux, double k, int *error) {
    const FieldPatch whole{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, 1.0, 0};
    // The same canonical FP64 geometry/integrand is used by both fixed-rule
    // evaluation and adaptive refinement. No rounded point cache survives.
    return fieldPatch(observer, panel, whole, rule, static_cast<int>(count), pressure, flux, k,
                      error)
        .value;
}

__device__ inline DoubleComplex asymptoticField(Vec3 direction, Vec3 observationOrigin,
                                                const Panel &panel,
                                                const galerkin::TriangleQuadraturePoint *rule,
                                                int count, const cuDoubleComplex *pressure,
                                                cuFloatComplex flux, double k, int *error) {
    DoubleComplex result{};
    const auto originDelta = sub(panel.vertices[0], observationOrigin);
    const double normalFactor = k * dot(direction, panel.normal);
    for (int sample = 0; sample < count; ++sample) {
        const auto &point = rule[sample];
        const auto position = add(originDelta, localPoint(panel, point.barycentric));
        // Integrate F(d) first. Applying exp(-ik R)/R separately preserves the
        // small aperture phase even when the reference distance is enormous.
        const double phase = k * dot(direction, position);
        if (!isfinite(phase)) {
            atomicExch(error, 1);
            return {};
        }
        double sine, cosine;
        sincos(phase, &sine, &cosine);
        const double scale = 2 * panel.area * point.weight / 12.56637061435917295385;
        DoubleComplex p{};
        for (unsigned i = 0; i < 3; ++i)
            p = plus(p, times(promoted(pressure[panel.nodes[i]]), point.barycentric[i]));
        const DoubleComplex integrand{-normalFactor * p.imag - flux.x,
                                      normalFactor * p.real - flux.y};
        result = plus(result, product({cosine * scale, sine * scale}, integrand));
    }
    return result;
}
static __device__ __noinline__ DoubleComplex
integrateNearField(Vec3 observer, const Panel &panel, const galerkin::TriangleQuadraturePoint *rule,
                   int count, const cuDoubleComplex *pressure, cuFloatComplex flux, double k,
                   double tolerance, unsigned maximumDepth, int *error) {
    const FieldPatch whole{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, 1.0, 0};
    FieldPatch stack[34];
    unsigned pending = 1;
    stack[0] = whole;
    DoubleComplex result{};
    while (pending) {
        const auto patch = stack[--pending];
        const auto coarse =
            fieldPatch(observer, panel, patch, rule, count, pressure, flux, k, error);
        FieldIntegral fine{};
        for (unsigned child = 0; child < 4; ++child) {
            const auto value = fieldPatch(observer, panel, childPatch(patch, child), rule, count,
                                          pressure, flux, k, error);
            fine.value = plus(fine.value, value.value);
            fine.absoluteScale += value.absoluteScale;
        }
        const double estimate =
            magnitude({fine.value.real - coarse.value.real, fine.value.imag - coarse.value.imag});
        if (estimate <= tolerance * fmax(fine.absoluteScale, 1e-300))
            result = plus(result, fine.value);
        else {
            if (patch.depth >= maximumDepth || pending + 4 > 34) {
                atomicExch(error, 3);
                return {};
            }
            for (unsigned child = 0; child < 4; ++child)
                stack[pending++] = childPatch(patch, child);
        }
    }
    return result;
}
} // namespace hornsculpt::cuda_detail
