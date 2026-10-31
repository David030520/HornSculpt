#pragma once

#include "panel_integrals.cuh"

namespace hornsculpt::cuda_detail {
struct PairIntegrationJob {
    unsigned test{}, trial{}, firstPatch{}, patchCount{}, swapBarycentric{};
    galerkin::PairPermutation permutation;
};
struct PairQuadratureView {
    const galerkin::TriangleQuadraturePoint *regular{};
    unsigned regularCount{};
    const galerkin::GaussPoint *gauss{};
    unsigned singularOrder{};
    const galerkin::PairQuadraturePatch *patches{};
    const galerkin::TriangleQuadraturePoint *nearRule{};
    unsigned baseNearCount{}, highNearCount{};
};
__device__ inline PairMoments promotedMoments(const MixedMoments &mixed) {
    PairMoments result{};
    result.green = promoted(mixed.green);
    for (unsigned i = 0; i < 9; ++i)
        result.matrix[i] = promoted(mixed.matrix[i]);
    for (unsigned i = 0; i < 3; ++i)
        result.rhs[i] = promoted(mixed.rhs[i]);
    return result;
}
__device__ inline cuDoubleComplex pressureCoefficient(const PairMoments &moments, const Panel &test,
                                                      const Panel &trial, unsigned i, unsigned j,
                                                      double k) {
    const double curl = dot(test.curls[i], trial.curls[j]) / k;
    return make_cuDoubleComplex(moments.matrix[3 * i + j].real + moments.green.imag * curl,
                                moments.matrix[3 * i + j].imag - moments.green.real * curl);
}

// The lane/stride contract also permits one regular pair per warp. Products
// retain the dense solver's FP32 integral convention; the warp specialization
// promotes each sample before accumulation so reduced parallelism introduces
// no additional FP32 summation of several samples in one lane.
template <bool promoteSamples = false>
__device__ inline PairMoments
integrateRegularPairThread(const Panel &test, const Panel &trial,
                           const galerkin::TriangleQuadraturePoint *regular, unsigned regularCount,
                           double k, int *error, unsigned first, unsigned stride) {
    MixedMoments mixed{};
    PairMoments moments{};
    for (unsigned sampleIndex = first; sampleIndex < regularCount * regularCount;
         sampleIndex += stride) {
        galerkin::PairQuadraturePoint sample;
        const auto &x = regular[sampleIndex / regularCount];
        const auto &y = regular[sampleIndex % regularCount];
        for (unsigned i = 0; i < 3; ++i) {
            sample.test[i] = x.barycentric[i];
            sample.trial[i] = y.barycentric[i];
        }
        sample.weight = x.weight * y.weight;
        MixedMoments one{};
        auto &accumulator = promoteSamples ? one : mixed;
        mixedSample(accumulator, pairDelta(test, trial, sample.test, sample.trial), test, trial,
                    sample.test, sample.trial,
                    static_cast<float>(sample.weight * 4 * test.area * trial.area), k,
                    1.f / static_cast<float>(k), false, error);
        if constexpr (promoteSamples) {
            const auto increment = promotedMoments(one);
            moments.green = plus(moments.green, increment.green);
            for (unsigned i = 0; i < 9; ++i)
                moments.matrix[i] = plus(moments.matrix[i], increment.matrix[i]);
            for (unsigned i = 0; i < 3; ++i)
                moments.rhs[i] = plus(moments.rhs[i], increment.rhs[i]);
        }
    }
    if constexpr (promoteSamples)
        return moments;
    else
        return promotedMoments(mixed);
}

// Each CUDA thread integrates its strided subset of the same regular,
// singular or adaptive pair rule. The caller reduces the moments across
// its 128-thread block before applying local basis masks and mass jumps.
__device__ inline PairMoments integratePairThread(const Panel &test, const Panel &trial,
                                                  const PairIntegrationJob &job,
                                                  PairQuadratureView view, double k, int *error) {
    const auto *regular = view.regular;
    const auto regularCount = view.regularCount;
    const auto *gauss = view.gauss;
    const auto singularOrder = view.singularOrder;
    const auto *patches = view.patches;
    const auto *nearRule = view.nearRule;
    const auto baseNearCount = view.baseNearCount, highNearCount = view.highNearCount;
    PairMoments moments{};
    if (job.permutation.topology != galerkin::PairTopology::disjoint) {
        const unsigned perSector = singularOrder * singularOrder * singularOrder * singularOrder;
        const unsigned count = perSector * galerkin::pairSectorCount(job.permutation.topology);
        for (unsigned sampleIndex = threadIdx.x; sampleIndex < count; sampleIndex += blockDim.x) {
            unsigned tuple = sampleIndex % perSector;
            const unsigned w = tuple % singularOrder;
            tuple /= singularOrder;
            const unsigned v = tuple % singularOrder;
            tuple /= singularOrder;
            const unsigned u = tuple % singularOrder;
            tuple /= singularOrder;
            const unsigned r = tuple;
            auto sample = galerkin::canonicalPairTransform(
                job.permutation.topology, sampleIndex / perSector, gauss[r].position,
                gauss[u].position, gauss[v].position, gauss[w].position);
            sample.weight *= gauss[r].weight * gauss[u].weight * gauss[v].weight * gauss[w].weight;
            sample = galerkin::permutePairPoint(sample, job.permutation);
            pairSample(moments, test, trial, sample,
                       pairDelta(test, trial, sample.test, sample.trial), k,
                       job.permutation.topology == galerkin::PairTopology::coincident, error);
        }
    } else if (job.patchCount) {
        for (unsigned patchIndex = 0; patchIndex < job.patchCount; ++patchIndex) {
            const auto &patch = patches[job.firstPatch + patchIndex];
            const bool promoted = patch.degree == 14;
            const auto *rule = nearRule + (promoted ? baseNearCount : 0);
            const unsigned count = promoted ? highNearCount : baseNearCount;
            const unsigned points = galerkin::nearPatchPointCount(patch, count);
            for (unsigned sampleIndex = threadIdx.x; sampleIndex < points;
                 sampleIndex += blockDim.x) {
                auto sample = galerkin::mapNearPatchPoint(rule, count, patch, sampleIndex,
                                                          job.swapBarycentric ? trial : test,
                                                          job.swapBarycentric ? test : trial);
                if (job.swapBarycentric)
                    for (unsigned i = 0; i < 3; ++i) {
                        const double value = sample.test[i];
                        sample.test[i] = sample.trial[i];
                        sample.trial[i] = value;
                    }
                pairSample(moments, test, trial, sample,
                           pairDelta(test, trial, sample.test, sample.trial), k, false, error);
            }
        }
    } else {
        moments = integrateRegularPairThread(test, trial, regular, regularCount, k, error,
                                             threadIdx.x, blockDim.x);
    }
    return moments;
}
__device__ inline cuDoubleComplex fluxCoefficient(const PairMoments &moments, const Panel &test,
                                                  unsigned i, double k, bool coincident) {
    auto value = make_cuDoubleComplex(moments.rhs[i].real, moments.rhs[i].imag);
    if (coincident)
        value.y += test.area / (6 * k);
    return value;
}
} // namespace hornsculpt::cuda_detail
