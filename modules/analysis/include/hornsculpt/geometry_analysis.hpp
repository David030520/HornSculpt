#pragma once

#include <cstddef>
#include <hornsculpt/draft_geometry.hpp>
#include <vector>
#include <stop_token>

namespace hornsculpt {

struct GeometrySample {
    double z{};
    double volume{};
    double area{};
    double areaDerivative{};
    double horizontalWallAngle{};
    double verticalWallAngle{};
    double leftWallAngle{};
    double bottomWallAngle{};
};

struct GeometryAnalysis {
    std::vector<GeometrySample> samples;
};

// SI units and radians. Volume is the horn air column from throat to z.
// Inlet boundaries add stations. At a crease, equal-z samples contain the
// before limit, undefined (NaN) derivative, and after limit, in that order.
GeometryAnalysis analyzeGeometry(const DraftGeometry &geometry, std::size_t sampleCount = 101,
                                 std::stop_token stop = {});

} // namespace hornsculpt
