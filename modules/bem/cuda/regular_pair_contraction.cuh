#pragma once

#include "panel_integrals.cuh"

namespace hornsculpt::cuda_detail {

// Degrees 4/6/8 have at most sixteen triangle points. One warp owns this
// scratch, so unrelated pairs in the same block never synchronize together.
struct RegularPairScratch {
    cuDoubleComplex source[5][16];
};

// Contract the tensor rule in two stages. At a fixed test point, only three
// pressure source moments, one flux and one Green moment are needed. This
// avoids thirteen complex per-lane accumulators and their warp reductions.
// Scalar kernels retain mixedKernelSample's FP64 phase and FP32 products;
// source moments and the final test contraction accumulate in FP64.
// All 32 lanes participate, including inactive test-point subgroups.
// Lanes 0..8 return A_ij including curl; lanes 9..11 return B_i. Mass jumps
// are excluded: this contract only integrates disjoint regular panel pairs.
__device__ inline cuDoubleComplex
contractRegularPairWarp(const Panel &test, const Panel &trial,
                        const galerkin::TriangleQuadraturePoint *rule, unsigned count, double k,
                        RegularPairScratch &scratch, int *error) {
    const unsigned lane = threadIdx.x & 31u, group = lane / 4u, sourceLane = lane & 3u;
    if (!count || count > 16) {
        if (!lane)
            atomicOr(error, 8);
        return {};
    }
    for (unsigned first = 0; first < count; first += 8) {
        const unsigned x = first + group;
        cuDoubleComplex moments[5]{};
        if (x < count) {
            const auto &testPoint = rule[x];
            for (unsigned y = sourceLane; y < count; y += 4) {
                const auto &trialPoint = rule[y];
                const auto kernel = mixedKernelSample(
                    pairDelta(test, trial, testPoint.barycentric, trialPoint.barycentric), test,
                    trial,
                    static_cast<float>(testPoint.weight * trialPoint.weight * 4 * test.area *
                                       trial.area),
                    k, 1.f / static_cast<float>(k), false, error);
#pragma unroll
                for (unsigned j = 0; j < 3; ++j) {
                    const double basis = static_cast<float>(trialPoint.barycentric[j]);
                    moments[j].x += static_cast<double>(kernel.pressure.x) * basis;
                    moments[j].y += static_cast<double>(kernel.pressure.y) * basis;
                }
                moments[3].x += kernel.flux.x;
                moments[3].y += kernel.flux.y;
                moments[4].x += kernel.green.x;
                moments[4].y += kernel.green.y;
            }
        }
#pragma unroll
        for (unsigned shift = 2; shift; shift /= 2)
#pragma unroll
            for (unsigned j = 0; j < 5; ++j) {
                moments[j].x += __shfl_down_sync(0xffffffffu, moments[j].x, shift, 4);
                moments[j].y += __shfl_down_sync(0xffffffffu, moments[j].y, shift, 4);
            }
        if (!sourceLane && x < count)
#pragma unroll
            for (unsigned j = 0; j < 5; ++j)
                scratch.source[j][x] = moments[j];
    }
    __syncwarp();
    cuDoubleComplex result{};
    if (lane < 12) {
        const unsigned row = lane < 9 ? lane / 3 : lane - 9;
        const unsigned sourceMoment = lane < 9 ? lane % 3 : 3;
        cuDoubleComplex green{};
        for (unsigned x = 0; x < count; ++x) {
            const double basis = static_cast<float>(rule[x].barycentric[row]);
            const auto moment = scratch.source[sourceMoment][x];
            result.x += basis * moment.x;
            result.y += basis * moment.y;
            if (lane < 9) {
                green.x += scratch.source[4][x].x;
                green.y += scratch.source[4][x].y;
            }
        }
        if (lane < 9) {
            const double curl = dot(test.curls[row], trial.curls[lane % 3]) / k;
            result.x += green.y * curl;
            result.y -= green.x * curl;
        }
    }
    return result;
}

} // namespace hornsculpt::cuda_detail
