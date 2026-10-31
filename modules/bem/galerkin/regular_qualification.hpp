#pragma once

#include "space.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace hornsculpt::galerkin {
#if defined(__CUDACC__)
#define HS_REGULAR_HD __host__ __device__
#else
#define HS_REGULAR_HD
#endif

struct RegularPanelBounds {
    Vec3 center;
    double radius{}, maximumEdge{};
    bool qualifiedShape{};
};

// One O(T) preparation, rather than a triangle-distance calculation for every
// regular pair. The enclosing centroid ball remains conservative under any
// rigid motion. Shape qualification is intentionally narrower than validity of
// a BEM panel; an unqualified panel still has the adaptive integration path.
inline RegularPanelBounds makeRegularPanelBounds(const Panel &panel) {
    RegularPanelBounds result{panel.centroid, 0, panel.maximumEdgeLength, true};
    double edges[3]{};
    for (unsigned i = 0; i < 3; ++i) {
        const auto &a = panel.vertices[i], &b = panel.vertices[(i + 1) % 3];
        edges[i] = std::hypot(a.x - b.x, a.y - b.y, a.z - b.z);
        result.radius =
            std::max(result.radius, std::hypot(a.x - panel.centroid.x, a.y - panel.centroid.y,
                                               a.z - panel.centroid.z));
    }
    // Inflate the computed enclosing radius outwards instead of allowing a
    // rounded-down norm to establish the cheap non-near certificate.
    result.radius = std::nextafter(result.radius * (1 + 8 * std::numeric_limits<double>::epsilon()),
                                   std::numeric_limits<double>::infinity());
    constexpr double minimumAngleCosine = .939692620785908384; // cos(20 degrees)
    for (unsigned i = 0; i < 3; ++i) {
        const double a = edges[i], b = edges[(i + 1) % 3], c = edges[(i + 2) % 3];
        result.qualifiedShape &=
            a > 0 && b > 0 && (a * a + b * b - c * c) <= 2 * a * b * minimumAngleCosine;
    }
    result.qualifiedShape &= std::isfinite(result.radius) && result.maximumEdge > 0;
    return result;
}

// Numerically qualified envelopes, NOT universal analytic error bounds. The
// independent holdouts compare all 36 G/h*grad(G) P1 moments with refined
// FP64 Duffy integration over electrical sizes, orientations, aspect ratios
// and size ratios. A factor four covers the observed rule errors; 1e-6 is
// reserved for the separately tested mixed-sample arithmetic. The qualifying
// error must fit the SAME 0.2*relativeTolerance budget as AdaptiveRule uses.
// Zero means the existing adaptive path, never unconditional regular q4.
HS_REGULAR_HD inline unsigned qualifiedRegularDegreeUnchecked(const RegularPanelBounds &test,
                                                              const RegularPanelBounds &source,
                                                              double wavenumber,
                                                              double relativeTolerance,
                                                              unsigned minimumDegree = 4) {
    if (!test.qualifiedShape || !source.qualifiedShape)
        return 0;
    const double h = test.maximumEdge > source.maximumEdge ? test.maximumEdge : source.maximumEdge,
                 kh = wavenumber * h, dx = test.center.x - source.center.x,
                 dy = test.center.y - source.center.y, dz = test.center.z - source.center.z,
                 distanceSquared = dx * dx + dy * dy + dz * dz,
                 radiusSum = test.radius + source.radius, budget = .2 * relativeTolerance;
    const auto separated = [&](double gap) {
        const double bound = radiusSum + gap * h;
        return distanceSquared >= bound * bound;
    };
    constexpr double arithmeticReserve = 1e-6;
    if (minimumDegree <= 4 &&
        ((kh <= .5 && separated(3) && arithmeticReserve + 4 * 2.3e-6 <= budget) ||
         (kh <= .25 && separated(8) && arithmeticReserve + 4 * 1e-7 <= budget)))
        return 4;
    if (minimumDegree <= 6 &&
        ((kh <= 2 && separated(1) && arithmeticReserve + 4 * 6e-7 <= budget) ||
         (kh <= 1 && separated(2) && arithmeticReserve + 4 * 2.4e-8 <= budget)))
        return 6;
    if (minimumDegree <= 8 && kh <= 2 &&
        ((separated(.5) && arithmeticReserve + 4 * 7.5e-7 <= budget) ||
         (separated(1) && arithmeticReserve + 4 * 3.5e-8 <= budget)))
        return 8;
    return 0;
}
inline unsigned qualifiedRegularDegree(const RegularPanelBounds &test,
                                       const RegularPanelBounds &source, double wavenumber,
                                       double relativeTolerance, unsigned minimumDegree = 4) {
    if (!std::isfinite(wavenumber) || wavenumber < 0 || !std::isfinite(relativeTolerance) ||
        relativeTolerance <= 0)
        throw std::invalid_argument("Regular quadrature qualification needs finite k/tolerance.");
    return qualifiedRegularDegreeUnchecked(test, source, wavenumber, relativeTolerance,
                                           minimumDegree);
}

// Electrical starting rule only. The final fixed-seed A/B action comparison
// separately tests convergence; this does not certify a mesh or an endpoint
// expansion. Electrically coarser panels require refinement instead of quietly
// extrapolating the measured domain of the available rules.
inline unsigned minimumEndpointDegree(double maximumKh, unsigned requestedDegree = 4) {
    if (!std::isfinite(maximumKh) || maximumKh < 0)
        throw std::invalid_argument("Endpoint quadrature requires a finite electrical panel size.");
    if (maximumKh > 3.25)
        throw std::runtime_error("BEM panel quadrature: maximum k*h=" + std::to_string(maximumKh) +
                                 " exceeds the qualified electrical range 3.25; refine the mesh.");
    const unsigned electrical = maximumKh <= .25 ? 4 : maximumKh <= 1 ? 6 : maximumKh <= 2 ? 8 : 12;
    return std::max(electrical, requestedDegree);
}

inline unsigned refinedEndpointDegree(unsigned degree) {
    return degree < 4    ? 4
           : degree < 6  ? 6
           : degree < 8  ? 8
           : degree < 12 ? 12
           : degree < 14 ? 14
                         : 0;
}
} // namespace hornsculpt::galerkin
#undef HS_REGULAR_HD
