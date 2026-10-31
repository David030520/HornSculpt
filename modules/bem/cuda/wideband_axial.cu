#include "wideband_contraction.cuh"
#include "wideband_operator.hpp"
#include "wideband_runtime.cuh"
#include <algorithm>

namespace hornsculpt::wideband {
namespace {
using namespace detail;
template <unsigned batchTile>
__global__ void applyAxial(const ContractionRowTile *rowTiles, const std::size_t *offsets,
                           const Z *coefficients, const Z *source, Z *target, unsigned sourceOrder,
                           unsigned targetOrder, unsigned firstBatch, unsigned batchCount,
                           int *status) {
    __shared__ Z matrix[contractionRows][contractionColumns + 1];
    __shared__ Z values[batchTile][2][contractionColumns];
    const auto tile = rowTiles[blockIdx.x];
    const unsigned m = tile.degree, localRow = threadIdx.x / contractionRowLanes,
                   row = tile.firstRow + localRow, lane = threadIdx.x % contractionRowLanes,
                   batchBase = blockIdx.y * batchTile, rows = targetOrder - m + 1,
                   columns = m <= sourceOrder ? sourceOrder - m + 1 : 0, signs = m == 0 ? 1 : 2;
    const auto sourceState = static_cast<std::size_t>(sourceOrder + 1) * (sourceOrder + 1),
               targetState = static_cast<std::size_t>(targetOrder + 1) * (targetOrder + 1);
    Z result[batchTile][2]{};
    for (unsigned firstColumn = 0; firstColumn < columns; firstColumn += contractionColumns) {
        for (unsigned entry = threadIdx.x; entry < contractionRows * contractionColumns;
             entry += contractionThreads) {
            const unsigned tileRow = entry / contractionColumns,
                           tileColumn = entry % contractionColumns,
                           globalRow = tile.firstRow + tileRow, column = firstColumn + tileColumn;
            Z value{};
            if (globalRow < rows && column < columns)
                value = coefficients[offsets[m] + static_cast<std::size_t>(globalRow) * columns +
                                     column];
            matrix[tileRow][tileColumn] = value;
        }
        for (unsigned entry = threadIdx.x; entry < batchTile * signs * contractionColumns;
             entry += contractionThreads) {
            const unsigned batch = entry / (signs * contractionColumns),
                           sign = entry / contractionColumns % signs,
                           tileColumn = entry % contractionColumns,
                           column = firstColumn + tileColumn, degree = m + column;
            Z value{};
            if (batchBase + batch < batchCount && column < columns) {
                const auto coefficient = static_cast<std::size_t>(degree) * degree + degree;
                const auto index = sign == 0 ? coefficient + m : coefficient - m;
                value = source[(static_cast<std::size_t>(firstBatch) + batchBase + batch) *
                                   sourceState +
                               index];
            }
            values[batch][sign][tileColumn] = value;
        }
        __syncthreads();
        // The axial block is identical for +m/-m. Reuse a matrix load for
        // both azimuthal signs as well as independent source batches.
        const unsigned columnsInTile =
            columns - firstColumn < contractionColumns ? columns - firstColumn : contractionColumns;
        if (row < rows)
            for (unsigned column = lane; column < columnsInTile; column += contractionRowLanes) {
                const auto a = matrix[localRow][column];
#pragma unroll
                for (unsigned batch = 0; batch < batchTile; ++batch)
                    for (unsigned sign = 0; sign < signs; ++sign) {
                        const auto b = values[batch][sign][column];
                        result[batch][sign].x += a.x * b.x - a.y * b.y;
                        result[batch][sign].y += a.x * b.y + a.y * b.x;
                    }
            }
        __syncthreads();
    }
#pragma unroll
    for (unsigned batch = 0; batch < batchTile; ++batch)
        for (unsigned sign = 0; sign < signs; ++sign) {
            const auto value = contractRowSum(result[batch][sign]);
            if (lane == 0 && row < rows && batchBase + batch < batchCount) {
                const unsigned degree = m + row;
                const auto coefficient = static_cast<std::size_t>(degree) * degree + degree,
                           index = sign == 0 ? coefficient + m : coefficient - m,
                           flat = (static_cast<std::size_t>(firstBatch) + batchBase + batch) *
                                      targetState +
                                  index;
                storeContractedValue(target, flat, value, status);
            }
        }
}

template <unsigned batchTile>
void launchAxial(const ContractionRowTile *rowTiles, unsigned tileCount, const std::size_t *offsets,
                 const Z *coefficients, const Z *source, Z *target, unsigned sourceOrder,
                 unsigned targetOrder, unsigned firstBatch, unsigned batchCount, int *status,
                 cudaStream_t stream) {
    const dim3 blocks(tileCount, (batchCount + batchTile - 1) / batchTile);
    applyAxial<batchTile><<<blocks, contractionThreads, 0, stream>>>(
        rowTiles, offsets, coefficients, source, target, sourceOrder, targetOrder, firstBatch,
        batchCount, status);
}
} // namespace

struct ScaledAxialOperator::Impl {
    Stream stream;
    DeviceBuffer<std::size_t> offsets;
    DeviceBuffer<Z> coefficients;
    DeviceBuffer<int> error;
    // Immutable within the same graph-protected class generation as its
    // matrix and offsets. Application never allocates or replaces this table.
    DeviceBuffer<ContractionRowTile> rowTiles;
    struct HostBuffers {
        DeviceBuffer<Z> source, target;
        HostBuffers(std::size_t sourceCount, std::size_t targetCount) {
            source.allocate(sourceCount);
            target.allocate(targetCount);
        }
    };
    std::unique_ptr<HostBuffers> host;
    unsigned sourceOrder{}, targetOrder{};
    unsigned rowTileCount{};
    std::size_t sourceCount{}, targetCount{}, bytes{}, hostBytes{}, matrixWork{};
    Impl(const ScaledAxialPlan &plan, unsigned batches, std::stop_token stop, cudaStream_t borrowed)
        : stream(borrowed) {
        cancelled(stop);
        if (!batches || plan.sourceOrder > maximumDeviceScaledOrder ||
            plan.targetOrder > maximumDeviceScaledOrder ||
            plan.blockOffsets.size() != std::min(plan.sourceOrder, plan.targetOrder) + 2 ||
            plan.blockOffsets.front() != 0 || plan.blockOffsets.back() != plan.coefficients.size())
            throw std::invalid_argument("Wideband GPU axial plan dimensions are invalid.");
        for (unsigned m = 0; m <= std::min(plan.sourceOrder, plan.targetOrder); ++m) {
            const auto count =
                static_cast<std::size_t>(plan.sourceOrder - m + 1) * (plan.targetOrder - m + 1);
            if (plan.blockOffsets[m] > plan.blockOffsets[m + 1] ||
                plan.blockOffsets[m + 1] - plan.blockOffsets[m] != count)
                throw std::invalid_argument("Wideband GPU axial m-block span is inconsistent.");
        }
        sourceOrder = plan.sourceOrder;
        targetOrder = plan.targetOrder;
        std::vector<ContractionRowTile> tiles;
        for (unsigned m = 0; m <= targetOrder; ++m)
            for (unsigned row = 0; row <= targetOrder - m; row += contractionRows)
                tiles.push_back({m, row});
        for (unsigned m = 0; m <= std::min(sourceOrder, targetOrder); ++m)
            matrixWork += product(product(sourceOrder - m + 1, targetOrder - m + 1), m ? 2 : 1);
        rowTileCount = static_cast<unsigned>(tiles.size());
        sourceCount = product(batches, product(sourceOrder + 1, sourceOrder + 1));
        targetCount = product(batches, product(targetOrder + 1, targetOrder + 1));
        bytes = sum(sum(product(plan.blockOffsets.size(), sizeof(std::size_t)),
                        product(plan.coefficients.size(), sizeof(Z))),
                    sum(sizeof(int), product(tiles.size(), sizeof(ContractionRowTile))));
        hostBytes = product(sum(sourceCount, targetCount), sizeof(Z));
        std::size_t free = 0, total = 0;
        checked(cudaMemGetInfo(&free, &total), "Wideband axial GPU memory query");
        if (bytes > free - free / 5)
            throw std::runtime_error(
                "Wideband axial class exceeds the available GPU memory budget.");
        const auto values = packed(plan.coefficients);
        offsets.allocate(plan.blockOffsets.size());
        coefficients.allocate(plan.coefficients.size());
        error.allocate(1);
        rowTiles.allocate(tiles.size());
        offsets.upload(plan.blockOffsets, stream);
        coefficients.upload(values, stream);
        rowTiles.upload(tiles, stream);
        checked(cudaStreamSynchronize(stream), "Wideband axial class preparation");
        cancelled(stop);
    }
    void applyDevice(const Z *input, Z *output, std::size_t activeCount, std::stop_token stop,
                     int *queuedStatus = nullptr) {
        cancelled(stop);
        if (!input || !output || input == output)
            throw std::invalid_argument(
                "Wideband axial application requires distinct device arrays.");
        const auto targetState = product(targetOrder + 1, targetOrder + 1);
        if (!activeCount || activeCount > targetCount || activeCount % targetState)
            throw std::invalid_argument("Wideband axial active batch count exceeds capacity.");
        const auto activeBatches = static_cast<unsigned>(activeCount / targetState);
        auto *status = queuedStatus ? queuedStatus : error.data();
        if (!queuedStatus)
            checked(cudaMemsetAsync(status, 0, sizeof(int), stream), "Wideband axial status reset");
        const auto capacity =
            std::max(std::size_t{1},
                     std::min({maximumBatchPoints / targetState, std::size_t{8000000} / matrixWork,
                               maximumContractionBatchGroups * 2}));
        for (std::size_t first = 0; first < activeBatches; first += capacity) {
            cancelled(stop);
            const auto batches = static_cast<unsigned>(std::min(capacity, activeBatches - first)),
                       firstBatch = static_cast<unsigned>(first);
            if (batches > 1)
                launchAxial<2>(rowTiles.data(), rowTileCount, offsets.data(), coefficients.data(),
                               input, output, sourceOrder, targetOrder, firstBatch, batches, status,
                               stream);
            else
                launchAxial<1>(rowTiles.data(), rowTileCount, offsets.data(), coefficients.data(),
                               input, output, sourceOrder, targetOrder, firstBatch, batches, status,
                               stream);
            checked(cudaGetLastError(), "Wideband scaled axial launch");
            if (!queuedStatus)
                checked(cudaStreamSynchronize(stream), "Wideband scaled axial batch");
        }
        if (!queuedStatus)
            completeSpectrumCheck(status, stream, stop);
    }
    std::vector<Complex> apply(std::span<const Complex> input, std::stop_token stop) {
        cancelled(stop);
        if (input.size() != sourceCount)
            throw std::invalid_argument("Wideband axial source has an invalid size.");
        if (!host) {
            std::size_t free = 0, total = 0;
            checked(cudaMemGetInfo(&free, &total), "Wideband axial staging memory query");
            if (hostBytes > free - free / 5)
                throw std::runtime_error(
                    "Wideband axial staging exceeds the available GPU memory budget.");
            host = std::make_unique<HostBuffers>(sourceCount, targetCount);
        }
        const auto values = packed(input);
        host->source.upload(values, stream);
        applyDevice(host->source.data(), host->target.data(), targetCount, stop);
        const auto output = host->target.download(stream);
        std::vector<Complex> result;
        result.reserve(output.size());
        for (const auto value : output)
            result.emplace_back(value.x, value.y);
        return result;
    }
};
ScaledAxialOperator::ScaledAxialOperator(const ScaledAxialPlan &plan, unsigned batchCount,
                                         std::stop_token stop, cudaStream_t stream)
    : impl_(std::make_unique<Impl>(plan, batchCount, stop, stream)) {}
ScaledAxialOperator::~ScaledAxialOperator() = default;
std::vector<Complex> ScaledAxialOperator::apply(std::span<const Complex> source,
                                                std::stop_token stop) {
    return impl_->apply(source, stop);
}
void ScaledAxialOperator::applyDevice(const cuDoubleComplex *source, cuDoubleComplex *output,
                                      std::stop_token stop) {
    impl_->applyDevice(source, output, impl_->targetCount, stop);
}
void ScaledAxialOperator::applyDeviceBatch(const Z *source, Z *output, unsigned activeBatches,
                                           std::stop_token stop) {
    impl_->applyDevice(source, output,
                       detail::product(activeBatches, detail::product(impl_->targetOrder + 1,
                                                                      impl_->targetOrder + 1)),
                       stop);
}
void ScaledAxialOperator::enqueueDeviceBatch(const Z *source, Z *output, unsigned activeBatches,
                                             int *status, std::stop_token stop) {
    if (!status)
        throw std::invalid_argument("Queued axial translation requires a device failure status.");
    impl_->applyDevice(source, output,
                       detail::product(activeBatches, detail::product(impl_->targetOrder + 1,
                                                                      impl_->targetOrder + 1)),
                       stop, status);
}
ScaledAxialDeviceIdentity ScaledAxialOperator::deviceIdentity() const noexcept {
    return {impl_->offsets.data(), impl_->coefficients.data(), impl_->sourceOrder,
            impl_->targetOrder,    impl_->sourceCount,         impl_->targetCount};
}
cudaStream_t ScaledAxialOperator::stream() const noexcept {
    return impl_->stream;
}
std::size_t ScaledAxialOperator::deviceBytes() const noexcept {
    return impl_->bytes + (impl_->host ? impl_->hostBytes : 0);
}
} // namespace hornsculpt::wideband
