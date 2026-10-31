#pragma once

#include <cstddef>
#include <hornsculpt/mesh.hpp>
#include <limits>
#include <stdexcept>

namespace hornsculpt::cuda_detail {

inline std::size_t checkedSizeProduct(std::size_t a, std::size_t b) {
    if (b && a > std::numeric_limits<std::size_t>::max() / b)
        throw std::overflow_error("BEM workspace size overflow.");
    return a * b;
}
inline std::size_t checkedSizeSum(std::size_t a, std::size_t b) {
    if (a > std::numeric_limits<std::size_t>::max() - b)
        throw std::overflow_error("BEM workspace size overflow.");
    return a + b;
}

// The existing allocation policy reserves 20% of currently free device memory.
inline std::size_t deviceMemoryBudget(std::size_t freeBytes) noexcept {
    return (freeBytes / 5) * 4 + (freeBytes % 5) * 4 / 5;
}

// Reserve the workspace against the free memory remaining after the immutable cache.
inline bool deviceCacheAndWorkspaceFit(std::size_t freeBytes, std::size_t cacheBytes,
                                       std::size_t workspaceBytes) noexcept {
    return cacheBytes <= deviceMemoryBudget(freeBytes) &&
           workspaceBytes <= deviceMemoryBudget(freeBytes - cacheBytes);
}

struct WorkspaceSize {
    std::size_t matrixElements{};
    std::size_t rhsElements{};
    std::size_t fluxElements{};
    std::size_t fieldElements{};
    std::size_t fieldPartialElements{};
    std::size_t tileCount{};
    std::size_t minimumBytes{};
};

// All solve buffers except the vendor-specific LU scratch space. Storage types
// come from the actual backend so the same checked formula also runs on the CPU.
template <class Complex, class PartialComplex>
WorkspaceSize solveWorkspaceSize(std::size_t pressureCount, std::size_t panelCount,
                                 std::size_t rhsCount, std::size_t observerCount,
                                 bool fastMultipole) {
    WorkspaceSize size;
    size.matrixElements = fastMultipole ? 0 : checkedSizeProduct(pressureCount, pressureCount);
    size.rhsElements = checkedSizeProduct(pressureCount, rhsCount);
    size.fluxElements = checkedSizeProduct(panelCount, rhsCount);
    size.fieldElements = checkedSizeProduct(observerCount, rhsCount);
    size.tileCount = panelCount / 256 + (panelCount % 256 != 0);
    size.fieldPartialElements = checkedSizeProduct(size.fieldElements, size.tileCount);
    const auto append = [&](std::size_t elements, std::size_t elementBytes) {
        size.minimumBytes =
            checkedSizeSum(size.minimumBytes, checkedSizeProduct(elements, elementBytes));
    };
    append(size.matrixElements, 2 * sizeof(Complex));
    append(size.rhsElements, 2 * sizeof(Complex));
    append(size.rhsElements, 2 * sizeof(PartialComplex));
    append(size.fluxElements, sizeof(Complex));
    append(size.fieldElements, sizeof(Complex));
    append(size.fieldPartialElements, sizeof(PartialComplex));
    append(observerCount, sizeof(Vec3));
    append(checkedSizeSum(fastMultipole ? 0 : pressureCount, 2), sizeof(int));
    return size;
}

template <class Complex, class PartialComplex>
WorkspaceSize denseWorkspaceSize(std::size_t pressureCount, std::size_t panelCount,
                                 std::size_t rhsCount, std::size_t observerCount) {
    return solveWorkspaceSize<Complex, PartialComplex>(pressureCount, panelCount, rhsCount,
                                                       observerCount, false);
}

} // namespace hornsculpt::cuda_detail
