#include "wideband_contraction.cuh"
#include "wideband_operator.hpp"
#include "wideband_runtime.cuh"
#include <algorithm>

namespace hornsculpt::wideband {
namespace {
using namespace detail;
__global__ void updateAzimuthPhases(Z *phases, unsigned order, double phi) {
    const unsigned index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index > 2 * order)
        return;
    const int m = static_cast<int>(index) - static_cast<int>(order);
    double sine, cosine;
    sincos(m * phi, &sine, &cosine);
    phases[index] = make_cuDoubleComplex(cosine, sine);
}
template <bool toAxis, unsigned batchTile>
__global__ void rotateCoefficients(const ContractionRowTile *rowTiles, const std::size_t *offsets,
                                   const double *smallD, const Z *phases, const Z *source,
                                   Z *target, unsigned planOrder, unsigned activeOrder,
                                   unsigned firstBatch, unsigned batchCount, int *status) {
    __shared__ double matrix[contractionRows][contractionColumns + 1];
    __shared__ Z values[batchTile][contractionColumns];
    const auto tile = rowTiles[blockIdx.x];
    const unsigned degree = tile.degree, dimension = 2 * degree + 1,
                   localRow = threadIdx.x / contractionRowLanes, row = tile.firstRow + localRow,
                   lane = threadIdx.x % contractionRowLanes, batchBase = blockIdx.y * batchTile;
    const auto state = static_cast<std::size_t>(activeOrder + 1) * (activeOrder + 1),
               degreeBase = static_cast<std::size_t>(degree) * degree;
    Z result[batchTile]{};
    for (unsigned firstColumn = 0; firstColumn < dimension; firstColumn += contractionColumns) {
        // Read contiguous spans of the native row-major d block in both
        // directions; the forward transpose occurs only in shared memory.
        for (unsigned entry = threadIdx.x; entry < contractionRows * contractionColumns;
             entry += contractionThreads) {
            const unsigned tileRow = toAxis ? entry % contractionRows : entry / contractionColumns,
                           tileColumn =
                               toAxis ? entry / contractionRows : entry % contractionColumns,
                           globalRow = tile.firstRow + tileRow, column = firstColumn + tileColumn;
            double value = 0;
            if (globalRow < dimension && column < dimension) {
                const auto matrixRow = toAxis ? column : globalRow,
                           matrixColumn = toAxis ? globalRow : column;
                value = smallD[offsets[degree] + static_cast<std::size_t>(matrixRow) * dimension +
                               matrixColumn];
            }
            matrix[tileRow][tileColumn] = value;
        }
        for (unsigned entry = threadIdx.x; entry < batchTile * contractionColumns;
             entry += contractionThreads) {
            const unsigned batch = entry / contractionColumns,
                           column = firstColumn + entry % contractionColumns;
            Z value{};
            if (batchBase + batch < batchCount && column < dimension) {
                value = source[(static_cast<std::size_t>(firstBatch) + batchBase + batch) * state +
                               degreeBase + column];
                if constexpr (toAxis) {
                    const auto phase = phases[planOrder + column - degree];
                    value = make_cuDoubleComplex(value.x * phase.x - value.y * phase.y,
                                                 value.x * phase.y + value.y * phase.x);
                }
            }
            values[batch][entry % contractionColumns] = value;
        }
        __syncthreads();
        const unsigned columnsInTile = dimension - firstColumn < contractionColumns
                                           ? dimension - firstColumn
                                           : contractionColumns;
        if (row < dimension)
            for (unsigned column = lane; column < columnsInTile; column += contractionRowLanes) {
                const double coefficient = matrix[localRow][column];
#pragma unroll
                for (unsigned batch = 0; batch < batchTile; ++batch) {
                    const auto value = values[batch][column];
                    result[batch].x += coefficient * value.x;
                    result[batch].y += coefficient * value.y;
                }
            }
        __syncthreads();
    }
#pragma unroll
    for (unsigned batch = 0; batch < batchTile; ++batch) {
        auto value = contractRowSum(result[batch]);
        if (lane == 0 && row < dimension && batchBase + batch < batchCount) {
            if constexpr (!toAxis) {
                const auto phase = phases[planOrder + row - degree];
                value = make_cuDoubleComplex(value.x * phase.x + value.y * phase.y,
                                             value.y * phase.x - value.x * phase.y);
            }
            const auto flat = (static_cast<std::size_t>(firstBatch) + batchBase + batch) * state +
                              degreeBase + row;
            storeContractedValue(target, flat, value, status);
        }
    }
}

template <bool toAxis, unsigned batchTile>
void launchRotation(const ContractionRowTile *rowTiles, unsigned tileCount,
                    const std::size_t *offsets, const double *smallD, const Z *phases,
                    const Z *source, Z *target, unsigned planOrder, unsigned activeOrder,
                    unsigned firstBatch, unsigned batchCount, int *status, cudaStream_t stream) {
    const dim3 blocks(tileCount, (batchCount + batchTile - 1) / batchTile);
    rotateCoefficients<toAxis, batchTile><<<blocks, contractionThreads, 0, stream>>>(
        rowTiles, offsets, smallD, phases, source, target, planOrder, activeOrder, firstBatch,
        batchCount, status);
}
} // namespace
struct SphericalRotationOperator::Impl {
    Stream stream;
    DeviceBuffer<std::size_t> offsets;
    DeviceBuffer<double> smallD;
    DeviceBuffer<Z> phases;
    DeviceBuffer<int> error;
    // Allocated once with the native matrix; graph cache class generations
    // protect their coupled owner from destruction or pointer reuse.
    DeviceBuffer<ContractionRowTile> rowTiles;
    std::vector<unsigned> degreeTileCounts;
    struct HostBuffers {
        DeviceBuffer<Z> source, target;
        explicit HostBuffers(std::size_t count) {
            source.allocate(count);
            target.allocate(count);
        }
    };
    std::unique_ptr<HostBuffers> host;
    unsigned order{}, batchCapacity{};
    std::size_t count{}, bytes{}, hostBytes{};
    Impl(const SphericalRotationPlan &plan, unsigned batches, std::stop_token stop,
         cudaStream_t borrowed)
        : stream(borrowed) {
        cancelled(stop);
        if (!batches || plan.order > maximumDeviceScaledOrder ||
            plan.degreeOffsets.size() != plan.order + 2 || plan.degreeOffsets.front() != 0 ||
            plan.degreeOffsets.back() != plan.smallD.size() ||
            plan.azimuthPhases.size() != 2 * plan.order + 1)
            throw std::invalid_argument("Wideband GPU rotation plan dimensions are invalid.");
        for (unsigned l = 0; l <= plan.order; ++l)
            if (plan.degreeOffsets[l] > plan.degreeOffsets[l + 1] ||
                plan.degreeOffsets[l + 1] - plan.degreeOffsets[l] !=
                    static_cast<std::size_t>(2 * l + 1) * (2 * l + 1))
                throw std::invalid_argument("Wideband GPU rotation degree span is invalid.");
        for (const auto value : plan.smallD)
            if (!std::isfinite(value))
                throw std::invalid_argument("Wideband GPU rotation coefficient is non-finite.");
        order = plan.order;
        std::vector<ContractionRowTile> tiles;
        degreeTileCounts.reserve(order + 1);
        for (unsigned degree = 0; degree <= order; ++degree) {
            for (unsigned row = 0; row < 2 * degree + 1; row += contractionRows)
                tiles.push_back({degree, row});
            degreeTileCounts.push_back(static_cast<unsigned>(tiles.size()));
        }
        batchCapacity = batches;
        count = product(batches, product(order + 1, order + 1));
        bytes = sum(sum(product(plan.degreeOffsets.size(), sizeof(std::size_t)),
                        product(plan.smallD.size(), sizeof(double))),
                    sum(sum(product(plan.azimuthPhases.size(), sizeof(Z)), sizeof(int)),
                        product(tiles.size(), sizeof(ContractionRowTile))));
        hostBytes = product(product(2, count), sizeof(Z));
        std::size_t free = 0, total = 0;
        checked(cudaMemGetInfo(&free, &total), "Wideband rotation GPU memory query");
        if (bytes > free - free / 5)
            throw std::runtime_error(
                "Wideband rotation class exceeds the available GPU memory budget.");
        const auto phaseValues = packed(plan.azimuthPhases);
        offsets.allocate(plan.degreeOffsets.size());
        smallD.allocate(plan.smallD.size());
        phases.allocate(phaseValues.size());
        error.allocate(1);
        rowTiles.allocate(tiles.size());
        offsets.upload(plan.degreeOffsets, stream);
        smallD.upload(plan.smallD, stream);
        phases.upload(phaseValues, stream);
        rowTiles.upload(tiles, stream);
        checked(cudaStreamSynchronize(stream), "Wideband rotation class preparation");
        cancelled(stop);
    }
    void applyDevice(const Z *input, Z *output, bool toAxis, unsigned activeOrder,
                     unsigned activeBatches, std::stop_token stop, int *queuedStatus = nullptr) {
        cancelled(stop);
        if (!input || !output || input == output)
            throw std::invalid_argument("Wideband rotation requires distinct device arrays.");
        if (activeOrder > order)
            throw std::invalid_argument("Wideband rotation active order exceeds the stored plan.");
        if (!activeBatches || activeBatches > batchCapacity)
            throw std::invalid_argument("Wideband rotation active batch count exceeds capacity.");
        auto *status = queuedStatus ? queuedStatus : error.data();
        if (!queuedStatus)
            checked(cudaMemsetAsync(status, 0, sizeof(int), stream),
                    "Wideband rotation status reset");
        const auto state = product(activeOrder + 1, activeOrder + 1);
        const auto matrixWork =
            product(activeOrder + 1, product(2 * activeOrder + 1, 2 * activeOrder + 3)) / 3;
        const auto capacity = std::max(
            std::size_t{1}, std::min({maximumBatchPoints / state, std::size_t{8000000} / matrixWork,
                                      maximumContractionBatchGroups * 4}));
        for (std::size_t first = 0; first < activeBatches; first += capacity) {
            cancelled(stop);
            const auto batches = static_cast<unsigned>(std::min(capacity, activeBatches - first)),
                       firstBatch = static_cast<unsigned>(first);
            const auto launch = [&]<bool forward>() {
                if (batches > 2)
                    launchRotation<forward, 4>(rowTiles.data(), degreeTileCounts[activeOrder],
                                               offsets.data(), smallD.data(), phases.data(), input,
                                               output, order, activeOrder, firstBatch, batches,
                                               status, stream);
                else if (batches == 2)
                    launchRotation<forward, 2>(rowTiles.data(), degreeTileCounts[activeOrder],
                                               offsets.data(), smallD.data(), phases.data(), input,
                                               output, order, activeOrder, firstBatch, batches,
                                               status, stream);
                else
                    launchRotation<forward, 1>(rowTiles.data(), degreeTileCounts[activeOrder],
                                               offsets.data(), smallD.data(), phases.data(), input,
                                               output, order, activeOrder, firstBatch, batches,
                                               status, stream);
            };
            if (toAxis)
                launch.template operator()<true>();
            else
                launch.template operator()<false>();
            checked(cudaGetLastError(), "Wideband rotation launch");
            if (!queuedStatus)
                checked(cudaStreamSynchronize(stream), "Wideband rotation batch");
        }
        if (!queuedStatus)
            completeSpectrumCheck(status, stream, stop);
    }
    std::vector<Complex> apply(std::span<const Complex> input, bool toAxis, std::stop_token stop) {
        cancelled(stop);
        if (input.size() != count)
            throw std::invalid_argument("Wideband rotation source has an invalid size.");
        if (!host) {
            std::size_t free = 0, total = 0;
            checked(cudaMemGetInfo(&free, &total), "Wideband rotation staging memory query");
            if (hostBytes > free - free / 5)
                throw std::runtime_error(
                    "Wideband rotation staging exceeds the GPU memory budget.");
            host = std::make_unique<HostBuffers>(count);
        }
        const auto values = packed(input);
        host->source.upload(values, stream);
        applyDevice(host->source.data(), host->target.data(), toAxis, order, batchCapacity, stop);
        const auto output = host->target.download(stream);
        std::vector<Complex> result;
        result.reserve(output.size());
        for (const auto value : output)
            result.emplace_back(value.x, value.y);
        return result;
    }
};
SphericalRotationOperator::SphericalRotationOperator(const SphericalRotationPlan &plan,
                                                     unsigned batches, std::stop_token stop,
                                                     cudaStream_t stream)
    : impl_(std::make_unique<Impl>(plan, batches, stop, stream)) {}
SphericalRotationOperator::~SphericalRotationOperator() = default;
std::vector<Complex> SphericalRotationOperator::toAxis(std::span<const Complex> input,
                                                       std::stop_token stop) {
    return impl_->apply(input, true, stop);
}
std::vector<Complex> SphericalRotationOperator::fromAxis(std::span<const Complex> input,
                                                         std::stop_token stop) {
    return impl_->apply(input, false, stop);
}
void SphericalRotationOperator::toAxisDevice(const cuDoubleComplex *input, cuDoubleComplex *output,
                                             std::stop_token stop) {
    toAxisPrefixDeviceBatch(input, output, impl_->order, impl_->batchCapacity, stop);
}
void SphericalRotationOperator::fromAxisDevice(const cuDoubleComplex *input,
                                               cuDoubleComplex *output, std::stop_token stop) {
    fromAxisPrefixDeviceBatch(input, output, impl_->order, impl_->batchCapacity, stop);
}
void SphericalRotationOperator::toAxisDeviceBatch(const Z *input, Z *output, unsigned activeBatches,
                                                  std::stop_token stop) {
    toAxisPrefixDeviceBatch(input, output, impl_->order, activeBatches, stop);
}
void SphericalRotationOperator::fromAxisDeviceBatch(const Z *input, Z *output,
                                                    unsigned activeBatches, std::stop_token stop) {
    fromAxisPrefixDeviceBatch(input, output, impl_->order, activeBatches, stop);
}
void SphericalRotationOperator::toAxisPrefixDeviceBatch(const Z *input, Z *output,
                                                        unsigned activeOrder,
                                                        unsigned activeBatches,
                                                        std::stop_token stop) {
    impl_->applyDevice(input, output, true, activeOrder, activeBatches, stop);
}
void SphericalRotationOperator::fromAxisPrefixDeviceBatch(const Z *input, Z *output,
                                                          unsigned activeOrder,
                                                          unsigned activeBatches,
                                                          std::stop_token stop) {
    impl_->applyDevice(input, output, false, activeOrder, activeBatches, stop);
}
void SphericalRotationOperator::enqueueToAxisDeviceBatch(const Z *input, Z *output,
                                                         unsigned activeBatches, int *status,
                                                         std::stop_token stop) {
    enqueueToAxisPrefixDeviceBatch(input, output, impl_->order, activeBatches, status, stop);
}
void SphericalRotationOperator::enqueueFromAxisDeviceBatch(const Z *input, Z *output,
                                                           unsigned activeBatches, int *status,
                                                           std::stop_token stop) {
    enqueueFromAxisPrefixDeviceBatch(input, output, impl_->order, activeBatches, status, stop);
}
void SphericalRotationOperator::enqueueToAxisPrefixDeviceBatch(const Z *input, Z *output,
                                                               unsigned activeOrder,
                                                               unsigned activeBatches, int *status,
                                                               std::stop_token stop) {
    if (!status)
        throw std::invalid_argument("Queued rotation requires a device failure status.");
    impl_->applyDevice(input, output, true, activeOrder, activeBatches, stop, status);
}
void SphericalRotationOperator::enqueueFromAxisPrefixDeviceBatch(const Z *input, Z *output,
                                                                 unsigned activeOrder,
                                                                 unsigned activeBatches,
                                                                 int *status,
                                                                 std::stop_token stop) {
    if (!status)
        throw std::invalid_argument("Queued rotation requires a device failure status.");
    impl_->applyDevice(input, output, false, activeOrder, activeBatches, stop, status);
}
SphericalRotationDeviceIdentity SphericalRotationOperator::deviceIdentity() const noexcept {
    return {impl_->offsets.data(), impl_->smallD.data(), impl_->phases.data(), impl_->order,
            impl_->count};
}
cudaStream_t SphericalRotationOperator::stream() const noexcept {
    return impl_->stream;
}
std::size_t SphericalRotationOperator::deviceBytes() const noexcept {
    return impl_->bytes + (impl_->host ? impl_->hostBytes : 0);
}
void SphericalRotationOperator::setAzimuth(double phi, std::stop_token stop) {
    detail::cancelled(stop);
    if (!std::isfinite(phi) || !std::isfinite(impl_->order * phi))
        throw std::invalid_argument("Wideband rotation azimuth or scaled phase is non-finite.");
    updateAzimuthPhases<<<(2 * impl_->order + 1 + detail::threads - 1) / detail::threads,
                          detail::threads, 0, impl_->stream>>>(impl_->phases.data(), impl_->order,
                                                               phi);
    detail::checked(cudaGetLastError(), "Wideband rotation azimuth update launch");
    detail::cancelled(stop);
}
} // namespace hornsculpt::wideband
