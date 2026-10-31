#pragma once
#include <algorithm>
#include <cstddef>
#include <stdexcept>

namespace hornsculpt::wideband {
struct RcrQueuePolicy {
    unsigned batchCapacity{}, maximumChunkBatches{};
    std::size_t workPerLink{};
    static constexpr unsigned maximumQueuedBatches = 32;
    static constexpr std::size_t maximumChunkWork = 64 * 1024 * 1024;
};
// Conservative complex coefficient operations, not a hardware-time promise.
// A rotation has at most (p+1)^2*(2p+1) terms per link. Include both rotations,
// the axial upper bound targetStates*(sourceOrder+1), and both endpoint phases.
// Small electrical orders amortize status readbacks, while high orders reduce
// even the leading active batch so one chunk cannot exceed this work budget.
inline RcrQueuePolicy makeRcrQueuePolicy(unsigned sourceOrder, unsigned targetOrder,
                                         unsigned reservedBatches) {
    if (sourceOrder > 128 || targetOrder > 128 || !reservedBatches || reservedBatches > 256)
        throw std::invalid_argument("LF RCR queue requires supported orders and batch capacity.");
    const auto source = std::size_t(sourceOrder + 1) * (sourceOrder + 1),
               target = std::size_t(targetOrder + 1) * (targetOrder + 1),
               work = source * (2 * sourceOrder + 2) +
                      target * (2 * targetOrder + 2 + sourceOrder + 1);
    const auto capacity = static_cast<unsigned>(
        std::min(std::size_t(reservedBatches), RcrQueuePolicy::maximumChunkWork / work));
    if (!capacity)
        throw std::runtime_error("One LF RCR link exceeds the bounded queue work domain.");
    const auto chunks =
        static_cast<unsigned>(std::min(std::size_t{RcrQueuePolicy::maximumQueuedBatches},
                                       RcrQueuePolicy::maximumChunkWork / (capacity * work)));
    return {capacity, chunks, work};
}
} // namespace hornsculpt::wideband
