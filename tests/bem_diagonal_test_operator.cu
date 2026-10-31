#include "bem_diagonal_test_operator.hpp"
#include "cuda/wideband_runtime.cuh"
#include <algorithm>
#include <cmath>
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <limits>
#include <stdexcept>
#include <string>

namespace hornsculpt::wideband {
namespace {
using namespace detail;
__global__ void translateLevel(const std::size_t *offsets, const SpectralInteraction *interactions,
                               const Z *translation, const Z *sources, Z *locals,
                               unsigned angularCount, std::size_t first, std::size_t count,
                               std::size_t interactionFirst) {
    const auto within = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (within >= count)
        return;
    const auto flat = first + within;
    const auto target = flat / angularCount;
    const auto point = flat % angularCount;
    const auto firstPair = offsets[target], pairCount = offsets[target + 1] - firstPair;
    Z result = interactionFirst ? locals[flat] : make_cuDoubleComplex(0, 0);
    if (interactionFirst < pairCount) {
        const auto batchEnd = interactionFirst + maximumBatchInteractions;
        const auto end = pairCount < batchEnd ? pairCount : batchEnd;
        for (auto pair = interactionFirst; pair < end; ++pair) {
            const auto job = interactions[firstPair + pair];
            const auto coefficient =
                translation[static_cast<std::size_t>(job.translationClass) * angularCount + point];
            const auto source =
                sources[static_cast<std::size_t>(job.source) * angularCount + point];
            result.x += coefficient.x * source.x - coefficient.y * source.y;
            result.y += coefficient.x * source.y + coefficient.y * source.x;
        }
    }
    locals[flat] = result;
}
} // namespace

struct DiagonalM2LOperator::Impl {
    Stream stream;
    DeviceBuffer<std::size_t> offsets;
    DeviceBuffer<SpectralInteraction> interactions;
    DeviceBuffer<Z> translations;
    DeviceBuffer<int> error;
    struct HostBuffers {
        DeviceBuffer<Z> sources, locals;
        HostBuffers(std::size_t sourceCount, std::size_t targetCount) {
            sources.allocate(sourceCount);
            locals.allocate(targetCount);
        }
    };
    std::unique_ptr<HostBuffers> host;
    std::size_t hostBytes{};
    std::size_t sourceCount{}, targetCount{}, maximumIncoming{}, bytes{};
    unsigned angularCount{};
    explicit Impl(const DiagonalLevelPlan &plan, std::stop_token stop, cudaStream_t borrowed)
        : stream(borrowed) {
        cancelled(stop);
        if (!plan.sourceNodes || !plan.targetNodes || !plan.angularSamples ||
            plan.targetOffsets.size() != static_cast<std::size_t>(plan.targetNodes) + 1 ||
            plan.targetOffsets.front() != 0 ||
            plan.targetOffsets.back() != plan.interactions.size())
            throw std::invalid_argument("Wideband level dimensions or CSR offsets are invalid.");
        for (unsigned node = 0; node < plan.targetNodes; ++node) {
            if (plan.targetOffsets[node] > plan.targetOffsets[node + 1])
                throw std::invalid_argument("Wideband level CSR offsets are not monotone.");
            maximumIncoming =
                std::max(maximumIncoming, plan.targetOffsets[node + 1] - plan.targetOffsets[node]);
        }
        for (const auto job : plan.interactions)
            if (job.source >= plan.sourceNodes || job.translationClass >= plan.translationClasses)
                throw std::invalid_argument("Wideband level interaction index is invalid.");
        angularCount = plan.angularSamples;
        sourceCount = product(plan.sourceNodes, angularCount);
        targetCount = product(plan.targetNodes, angularCount);
        const auto translationCount = product(plan.translationClasses, angularCount);
        if (translationCount != plan.translationCoefficients.size())
            throw std::invalid_argument("Wideband translation class table has an invalid size.");
        bytes = sum(product(plan.targetOffsets.size(), sizeof(std::size_t)),
                    product(plan.interactions.size(), sizeof(SpectralInteraction)));
        bytes = sum(sum(bytes, product(translationCount, sizeof(Z))), sizeof(int));
        hostBytes = product(sum(sourceCount, targetCount), sizeof(Z));
        std::size_t free = 0, total = 0;
        checked(cudaMemGetInfo(&free, &total), "Wideband GPU memory query");
        if (bytes > free - free / 5)
            throw std::runtime_error("Wideband level exceeds the available GPU memory budget.");
        const auto coefficients = packed(plan.translationCoefficients);
        offsets.allocate(plan.targetOffsets.size());
        interactions.allocate(plan.interactions.size());
        translations.allocate(translationCount);
        error.allocate(1);
        offsets.upload(plan.targetOffsets, stream);
        interactions.upload(plan.interactions, stream);
        translations.upload(coefficients, stream);
        checked(cudaStreamSynchronize(stream), "Wideband level preparation");
        cancelled(stop);
    }
    void applyDevice(const Z *input, Z *output, std::stop_token stop) {
        cancelled(stop);
        if (!input || !output || input == output)
            throw std::invalid_argument("Wideband M2L requires distinct non-null device arrays.");
        // The first tile initializes all local coefficients, including nodes
        // with no incoming interactions. Tile bounds retain cancellation points
        // without imposing a platform-dependent kernel wall-time assumption.
        for (std::size_t interactionFirst = 0;
             interactionFirst < std::max(std::size_t{1}, maximumIncoming);
             interactionFirst += maximumBatchInteractions)
            for (std::size_t first = 0; first < targetCount; first += maximumBatchPoints) {
                cancelled(stop);
                const auto count = std::min(maximumBatchPoints, targetCount - first);
                translateLevel<<<static_cast<unsigned>((count + threads - 1) / threads), threads, 0,
                                 stream>>>(offsets.data(), interactions.data(), translations.data(),
                                           input, output, angularCount, first, count,
                                           interactionFirst);
                checked(cudaGetLastError(), "Wideband M2L launch");
                checked(cudaStreamSynchronize(stream), "Wideband M2L batch");
            }
        verifySpectrum(output, targetCount, error.data(), stream, stop);
    }
    std::vector<Complex> apply(std::span<const Complex> input, std::stop_token stop) {
        cancelled(stop);
        if (input.size() != sourceCount)
            throw std::invalid_argument("Wideband source spectrum has an invalid size.");
        if (!host) {
            std::size_t free = 0, total = 0;
            checked(cudaMemGetInfo(&free, &total), "Wideband host-wrapper GPU memory query");
            if (hostBytes > free - free / 5)
                throw std::runtime_error(
                    "Wideband host-wrapper staging exceeds the GPU memory budget.");
            host = std::make_unique<HostBuffers>(sourceCount, targetCount);
        }
        const auto values = packed(input);
        host->sources.upload(values, stream);
        applyDevice(host->sources.data(), host->locals.data(), stop);
        const auto output = host->locals.download(stream);
        std::vector<Complex> result;
        result.reserve(output.size());
        for (const auto value : output)
            result.emplace_back(value.x, value.y);
        return result;
    }
};

DiagonalM2LOperator::DiagonalM2LOperator(const DiagonalLevelPlan &plan, std::stop_token stop,
                                         cudaStream_t stream)
    : impl_(std::make_unique<Impl>(plan, stop, stream)) {}
DiagonalM2LOperator::~DiagonalM2LOperator() = default;
std::vector<Complex> DiagonalM2LOperator::apply(std::span<const Complex> source,
                                                std::stop_token stop) {
    return impl_->apply(source, stop);
}
void DiagonalM2LOperator::applyDevice(const cuDoubleComplex *source, cuDoubleComplex *output,
                                      std::stop_token stop) {
    impl_->applyDevice(source, output, stop);
}
cudaStream_t DiagonalM2LOperator::stream() const noexcept {
    return impl_->stream;
}

std::size_t DiagonalM2LOperator::deviceBytes() const noexcept {
    return impl_->bytes + (impl_->host ? impl_->hostBytes : 0);
}

} // namespace hornsculpt::wideband
