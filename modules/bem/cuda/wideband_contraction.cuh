#pragma once

#include "wideband_runtime.cuh"

namespace hornsculpt::wideband::detail {
// A subgroup contracts one row; a block shares each matrix/source tile among
// sixteen rows and several independent right-hand sides. Keep the tile's
// column extent divisible by the subgroup width so every shuffle is converged.
inline constexpr unsigned contractionRows = 16;
inline constexpr unsigned contractionColumns = 32;
inline constexpr unsigned contractionRowLanes = 8;
inline constexpr unsigned contractionThreads = contractionRows * contractionRowLanes;
// CUDA guarantees at least this many blocks in grid.y. The host launch
// batches are bounded independently from the coefficient-work budget.
inline constexpr std::size_t maximumContractionBatchGroups = 65535;
struct ContractionRowTile {
    unsigned degree{}, firstRow{};
};

inline std::size_t rotationRowTileCount(unsigned order) {
    std::size_t count = 0;
    for (unsigned degree = 0; degree <= order; ++degree)
        count += (2 * degree + contractionRows) / contractionRows;
    return count;
}

inline std::size_t axialRowTileCount(unsigned targetOrder) {
    std::size_t count = 0;
    for (unsigned m = 0; m <= targetOrder; ++m)
        count += (targetOrder - m + contractionRows) / contractionRows;
    return count;
}

__device__ inline Z contractRowSum(Z value) {
    for (unsigned offset = contractionRowLanes / 2; offset; offset /= 2) {
        value.x += __shfl_down_sync(0xffffffffu, value.x, offset, contractionRowLanes);
        value.y += __shfl_down_sync(0xffffffffu, value.y, offset, contractionRowLanes);
    }
    return value;
}

__device__ inline void storeContractedValue(Z *target, std::size_t index, Z value, int *status) {
    target[index] = value;
    // The same sticky scalar is used throughout captured RCR applications.
    // Checking the producing register removes a separate full-spectrum read.
    if (!isfinite(value.x) || !isfinite(value.y))
        atomicExch(status, 1);
}
} // namespace hornsculpt::wideband::detail
