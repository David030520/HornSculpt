#include "cuda/workspace_size.hpp"
#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace hornsculpt::cuda_detail;
namespace {
using Complex32 = std::array<float, 2>;
using Complex64 = std::array<double, 2>;
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class Function> void rejectsOverflow(Function function) {
    bool rejected = false;
    try {
        function();
    } catch (const std::overflow_error &) {
        rejected = true;
    }
    require(rejected, "An overflowing workspace size was accepted.");
}
} // namespace

int main() {
    try {
        const auto size = denseWorkspaceSize<Complex32, Complex64>(3, 257, 2, 5);
        require(size.matrixElements == 9 && size.rhsElements == 6 && size.fluxElements == 514 &&
                    size.fieldElements == 10 && size.tileCount == 2 &&
                    size.fieldPartialElements == 20 && size.minimumBytes == 5084,
                "Known solve-buffer layout or byte count changed.");
        require(denseWorkspaceSize<Complex32, Complex64>(3, 256, 2, 5).minimumBytes == 4908,
                "The exact field-tile boundary uses an extra tile.");
        require(denseWorkspaceSize<Complex32, Complex64>(3, 0, 2, 0).minimumBytes == 452,
                "An empty field buffer has an incorrect minimum size.");
        const auto fast = solveWorkspaceSize<Complex32, Complex64>(3, 257, 2, 5, true);
        require(fast.matrixElements == 0 && fast.minimumBytes == size.minimumBytes - 156 &&
                    fast.rhsElements == size.rhsElements &&
                    fast.fieldElements == size.fieldElements,
                "Fast Multipole allocated dense matrix/LU buffers or changed field storage.");
        const auto largeFast =
            solveWorkspaceSize<Complex32, Complex64>(90566, 181128, 1, 147, true);
        require(largeFast.matrixElements == 0 && largeFast.minimumBytes < 10000000,
                "Fast Multipole common workspace retains quadratic pressure storage.");
        constexpr auto maximum = std::numeric_limits<std::size_t>::max();
        require(checkedSizeProduct(maximum, 1) == maximum && checkedSizeProduct(maximum, 0) == 0 &&
                    checkedSizeSum(maximum - 1, 1) == maximum,
                "Representable boundary sizes were rejected.");
        rejectsOverflow([&] { (void)checkedSizeProduct(maximum, 2); });
        rejectsOverflow([&] { (void)checkedSizeSum(maximum, 1); });
        rejectsOverflow([&] { (void)denseWorkspaceSize<Complex32, Complex64>(maximum, 1, 1, 1); });
        rejectsOverflow(
            [&] { (void)denseWorkspaceSize<Complex32, Complex64>(0, 257, 1, maximum / 2 + 1); });
        const auto pressure =
            (std::size_t{1} << ((std::numeric_limits<std::size_t>::digits - 4) / 2)) - 1;
        rejectsOverflow([&] { (void)denseWorkspaceSize<Complex32, Complex64>(pressure, 1, 1, 0); });
        require(denseWorkspaceSize<Complex32, Complex64>(0, maximum, 0, 0).tileCount ==
                    maximum / 256 + 1,
                "Tile rounding overflowed at the size limit.");
        for (std::size_t free = 0; free <= 100; ++free)
            require(deviceMemoryBudget(free) == free * 8 / 10,
                    "The existing free-memory budget changed.");
        require(deviceCacheAndWorkspaceFit(1000, 200, 640) &&
                    !deviceCacheAndWorkspaceFit(1000, 200, 641) &&
                    !deviceCacheAndWorkspaceFit(1000, 801, 0) &&
                    !deviceCacheAndWorkspaceFit(1000, maximum, 0) &&
                    deviceCacheAndWorkspaceFit(maximum, 0, deviceMemoryBudget(maximum)),
                "Cache/workspace residency preflight changed its reserve or overflow behavior.");
        require(deviceMemoryBudget(maximum) == maximum - maximum / 5 - (maximum % 5 != 0),
                "The free-memory budget overflowed.");
        if constexpr (sizeof(std::size_t) >= 8)
            require(denseWorkspaceSize<Complex32, Complex64>(90566, 181128, 1, 147).minimumBytes ==
                        131243034080ULL,
                    "The measured Grand pressure space has an incorrect dense lower bound.");
        std::cout << "CPU dense workspace counts, byte bounds, tile rounding and overflow checks "
                     "passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
