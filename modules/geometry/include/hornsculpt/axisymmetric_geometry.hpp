#pragma once

#include <hornsculpt/axisymmetric.hpp>
#include <map>
#include <stop_token>

namespace hornsculpt {
class CadGeometry;
class DraftGeometry;

struct MeridianSampling {
    double maximumLength{.002};       // metres, meridian discretization only
    double chordTolerance{.00005};    // metres, geometric chord deviation
    double rotationalTolerance{1e-8}; // metres, exact-source symmetry check
    std::map<std::uint32_t, double> faceTargetLengths;
    std::map<std::uint32_t, BoundaryCondition> faceConditions;
};

// CAD must consist of complete coaxial analytic/revolved faces with rotational
// trims. The requested frame is checked, never silently estimated or recentered.
AxisymmetricMeridian makeAxisymmetricMeridian(const CadGeometry &geometry, Vec3 origin, Vec3 axis,
                                              const MeridianSampling &sampling = {},
                                              std::stop_token stop = {});
AxisymmetricMeridian makeAxisymmetricMeridian(const DraftGeometry &geometry,
                                              const MeridianSampling &sampling = {},
                                              std::stop_token stop = {});
} // namespace hornsculpt
