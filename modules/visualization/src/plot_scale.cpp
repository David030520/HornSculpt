#include <algorithm>
#include <cmath>
#include <hornsculpt/plot_scale.hpp>
#include <limits>
#include <stdexcept>

namespace hornsculpt {
void validatePlotAxisScale(const PlotAxisScale &s, bool logarithmic) {
    if (s.minimum.has_value() != s.maximum.has_value() ||
        (s.minimum &&
         (!std::isfinite(*s.minimum) || !std::isfinite(*s.maximum) || !(*s.maximum > *s.minimum) ||
          !std::isfinite(*s.maximum - *s.minimum) || (logarithmic && *s.minimum <= 0))) ||
        (s.tickStep && (!std::isfinite(*s.tickStep) || *s.tickStep <= 0)))
        throw std::invalid_argument(
            "Skalierung verlangt endliche, aufsteigende Grenzen und einen positiven Schritt.");
}
void validatePlotDbScale(const PlotDbScale &s) {
    validatePlotAxisScale({s.minimum, s.maximum, s.tickStep});
    if (std::ceil((s.maximum - s.minimum) / s.tickStep) + 1 > 1000)
        throw std::invalid_argument(
            "Die dB-Skala darf einschließlich Enden höchstens 1000 Teilstriche enthalten.");
}
std::vector<double> plotTicks(double lo, double hi, std::optional<double> step, unsigned pixels) {
    validatePlotAxisScale({lo, hi, step});
    double h = step.value_or((hi - lo) / 4);
    if (!step) {
        const double magnitude = std::pow(10., std::floor(std::log10(h)));
        if (magnitude > 0 && std::isfinite(magnitude)) {
            const double fraction = h / magnitude;
            h = magnitude * (fraction <= 1 ? 1 : fraction <= 2 ? 2 : fraction <= 5 ? 5 : 10);
        }
    }
    const double capacity = step ? 999. : std::max(2., static_cast<double>(pixels) / 12);
    const double requested = (hi - lo) / h;
    if (!std::isfinite(requested))
        h = (hi - lo) / capacity;
    else if (requested > capacity)
        h *= std::ceil(requested / capacity);
    if (!(h > 0) || !std::isfinite(h))
        return {lo, hi};
    const double first = std::ceil(lo / h) * h;
    std::vector<double> result;
    for (unsigned i = 0; i <= static_cast<unsigned>(capacity) + 1; ++i) {
        const double tick = first + i * h;
        if (!std::isfinite(tick) || tick > hi + (hi - lo) * 1e-12)
            break;
        if (tick >= lo)
            result.push_back(tick == 0 ? 0 : tick);
    }
    if (result.empty() || result.front() > lo + (hi - lo) * 1e-12)
        result.insert(result.begin(), lo);
    if (result.back() < hi - (hi - lo) * 1e-12)
        result.push_back(hi);
    return result;
}
} // namespace hornsculpt
