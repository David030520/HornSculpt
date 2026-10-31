#include <algorithm>
#include <cmath>
#include <hornsculpt/plot_axes.hpp>
#include <stdexcept>

namespace hornsculpt {
QString frequencyLabel(double frequency) {
    if (!std::isfinite(frequency) || frequency < 0)
        throw std::invalid_argument("Frequency labels require finite nonnegative Hz.");
    const int decimals = frequency == 0 ? 0
                         : frequency < 1
                             ? std::min(12, static_cast<int>(std::ceil(-std::log10(frequency))) + 2)
                         : frequency < 1000 ? 2
                                            : 1;
    auto text = QString::number(frequency, 'f', decimals);
    if (text.contains('.')) {
        while (text.endsWith('0'))
            text.chop(1);
        if (text.endsWith('.'))
            text.chop(1);
    }
    return text;
}

std::vector<double> frequencyTicks(double minimum, double maximum, double widthPixels) {
    if (!std::isfinite(minimum) || !std::isfinite(maximum) || minimum <= 0 || maximum < minimum ||
        !std::isfinite(widthPixels) || widthPixels <= 0)
        throw std::invalid_argument("Invalid frequency axis range.");
    if (minimum == maximum)
        return {minimum};
    std::vector<double> candidates{minimum};
    const int firstDecade = static_cast<int>(std::floor(std::log10(minimum)));
    const int lastDecade = static_cast<int>(std::ceil(std::log10(maximum)));
    for (int decade = firstDecade; decade <= lastDecade; ++decade) {
        const double base = std::pow(10.0, decade);
        for (const double multiplier : {1.0, 2.0, 5.0}) {
            const double tick = base * multiplier;
            if (std::isfinite(tick) && tick > minimum && tick < maximum)
                candidates.push_back(tick);
        }
    }
    std::vector<double> ticks{minimum};
    const double span = std::log(maximum) - std::log(minimum);
    for (const double tick : candidates) {
        const double fromPrevious = (std::log(tick) - std::log(ticks.back())) / span * widthPixels;
        const double toEnd = (std::log(maximum) - std::log(tick)) / span * widthPixels;
        if (fromPrevious >= 68 && toEnd >= 58)
            ticks.push_back(tick);
    }
    if (widthPixels >= 58)
        ticks.push_back(maximum);
    return ticks;
}
} // namespace hornsculpt
