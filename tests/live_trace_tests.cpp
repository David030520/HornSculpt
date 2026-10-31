#include "live_trace.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
void near(double actual, double expected, const char *message) {
    require(std::isfinite(actual) && std::abs(actual - expected) < 1e-10, message);
}
} // namespace

int main() {
    using namespace hornsculpt::live;
    try {
        const auto corrected =
            estimate({{1, 3}, {10, 0}, {100, 8}, {1000, 0}, {10000, 4}}, {{10, 2}, {1000, 6}});
        near(corrected.at(100), 12, "Log-midpoint correction must preserve the baseline bump.");
        near(corrected.at(10), 2, "First computed anchor was changed.");
        near(corrected.at(1000), 6, "Last computed anchor was changed.");
        near(corrected.at(1), 3, "Baseline below the new anchor hull must remain unchanged.");
        near(corrected.at(10000), 4, "Baseline above the new anchor hull must remain unchanged.");
        const Trace computed{{10, 2}, {1000, 6}};
        require(estimate({}, computed) == computed,
                "Without baseline, return only actual samples.");
        near(estimate({{100, 80}}, computed).at(100), 4,
             "Unavailable baseline anchors require logarithmic interpolation of new values.");
        const auto shifted = estimate({{10, 0}, {100, 10}, {1000, 20}},
                                      {{std::sqrt(1000.0), 9}, {std::sqrt(100000.0), 23}});
        require(shifted.size() == 5, "Changed frequency grids must retain their union.");
        near(shifted.at(100), 16, "Baseline anchor interpolation must use log-frequency.");
        near(shifted.at(10), 0, "New grids must not cause extrapolation below their anchor hull.");
        near(shifted.at(1000), 20,
             "New grids must not cause extrapolation above their anchor hull.");
        const double gap = std::numeric_limits<double>::quiet_NaN();
        const auto gaps = estimate({{10, 0}, {30, 3}, {100, gap}, {300, 3}, {1000, 0}},
                                   {{10, 2}, {100, gap}, {1000, 6}});
        require(std::isnan(gaps.at(100)), "Computed NaN must remain a gap.");
        near(gaps.at(30), 3, "Do not skip a missing computed anchor on the right.");
        near(gaps.at(300), 3, "Do not skip a missing computed anchor on the left.");
        require(std::isnan(estimate({{10, 0}, {100, gap}, {1000, 0}}, computed).at(100)),
                "Do not invent a missing baseline channel between valid anchors.");
        const auto anchorGap =
            estimate({{10, 0}, {100, gap}, {1000, 8}, {10000, 10}}, {{30, 2}, {10000, 6}});
        const double expected = std::lerp(
            2.0, 6.0, (std::log(1000.0) - std::log(30.0)) / (std::log(10000.0) - std::log(30.0)));
        near(anchorGap.at(1000), expected,
             "Baseline interpolation must not bridge a NaN to manufacture an anchor delta.");
        require(std::isnan(anchorGap.at(100)), "Fallback interpolation must retain baseline gaps.");
        near(estimate({{100, gap}}, {{100, 7}}).at(100), 7,
             "An actual sample must replace a former baseline gap exactly.");
        const Trace baseline{{10, 0}, {100, 5}, {1000, 0}};
        require(estimate(baseline, {{100, 9}}) == Trace({{10, 0}, {100, 9}, {1000, 0}}),
                "A single new anchor must not extrapolate corrections.");
        bool rejected = false;
        try {
            (void)estimate({{0, 1}}, computed);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Nonpositive frequency cannot define a logarithmic coordinate.");
        std::puts("Live trace checks passed.");
    } catch (const std::exception &error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
