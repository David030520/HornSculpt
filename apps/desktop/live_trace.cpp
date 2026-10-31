#include "live_trace.hpp"

#include <cmath>
#include <iterator>
#include <optional>
#include <stdexcept>

namespace hornsculpt::live {
namespace {
void validateFrequencies(const Trace &trace) {
    for (const auto &[frequency, value] : trace) {
        (void)value;
        if (!std::isfinite(frequency) || frequency <= 0)
            throw std::invalid_argument("Live trace frequencies must be positive and finite.");
    }
}
double fraction(double frequency, double left, double right) {
    return (std::log(frequency) - std::log(left)) / (std::log(right) - std::log(left));
}
std::optional<double> interpolate(const Trace &trace, double frequency) {
    const auto right = trace.lower_bound(frequency);
    if (right != trace.end() && right->first == frequency)
        return std::isfinite(right->second) ? std::optional<double>{right->second} : std::nullopt;
    if (right == trace.begin() || right == trace.end())
        return {};
    const auto left = std::prev(right);
    if (!std::isfinite(left->second) || !std::isfinite(right->second))
        return {};
    return std::lerp(left->second, right->second, fraction(frequency, left->first, right->first));
}
} // namespace

Trace estimate(const Trace &baseline, const Trace &computed) {
    validateFrequencies(baseline);
    validateFrequencies(computed);
    Trace result = baseline;
    for (const auto &[frequency, value] : baseline) {
        if (!std::isfinite(value) || computed.contains(frequency))
            continue;
        const auto right = computed.upper_bound(frequency);
        if (right == computed.begin() || right == computed.end())
            continue;
        const auto left = std::prev(right);
        if (!std::isfinite(left->second) || !std::isfinite(right->second))
            continue;
        const double t = fraction(frequency, left->first, right->first);
        const auto oldLeft = interpolate(baseline, left->first);
        const auto oldRight = interpolate(baseline, right->first);
        result[frequency] = oldLeft && oldRight ? value + std::lerp(left->second - *oldLeft,
                                                                    right->second - *oldRight, t)
                                                : std::lerp(left->second, right->second, t);
    }
    // Exact new values always win, including explicit gaps in a computed channel.
    for (const auto &[frequency, value] : computed)
        result[frequency] = value;
    return result;
}

} // namespace hornsculpt::live
