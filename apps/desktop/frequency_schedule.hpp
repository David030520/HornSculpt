#pragma once

#include <cstddef>
#include <queue>
#include <utility>
#include <vector>

namespace hornsculpt {
// Cover both ends first, then bisect the widest logarithmic-grid gaps.
// Equal-width gaps at higher frequencies receive priority.
inline std::vector<std::size_t> frequencySchedule(std::size_t count) {
    if (count == 0)
        return {};
    std::vector<std::size_t> order{count - 1};
    if (count == 1)
        return order;
    order.push_back(0);
    using Gap = std::pair<std::size_t, std::pair<std::size_t, std::size_t>>;
    std::priority_queue<Gap> gaps;
    gaps.push({count - 1, {0, count - 1}});
    while (!gaps.empty()) {
        const auto [width, bounds] = gaps.top();
        gaps.pop();
        if (width <= 1)
            continue;
        const auto [low, high] = bounds;
        const auto middle = low + (high - low) / 2;
        order.push_back(middle);
        gaps.push({middle - low, {low, middle}});
        gaps.push({high - middle, {middle, high}});
    }
    return order;
}
} // namespace hornsculpt
