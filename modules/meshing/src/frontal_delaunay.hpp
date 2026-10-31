#pragma once

#include "surface_charts.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <stop_token>
#include <string>
#include <vector>

namespace hornsculpt::detail {
using AtlasSegment = std::array<std::uint32_t, 2>;
// Requests coordinated refinement of current shared atlas constraints. It is
// distinct from numerical failures; a chart never splits a seam on its own.
class BoundaryRefinementRequired final : public std::runtime_error {
  public:
    BoundaryRefinementRequired(std::string message, std::vector<AtlasSegment> segments,
                               std::optional<double> violationEnergy = {});
    std::span<const AtlasSegment> segments() const noexcept {
        return segments_;
    }
    std::optional<double> violationEnergy() const noexcept {
        return violationEnergy_;
    }

  private:
    std::vector<AtlasSegment> segments_;
    std::optional<double> violationEnergy_;
};
struct FrontalVertex {
    ChartPoint uv;
    Vec3 position;
    // UINT32_MAX denotes an interior vertex owned only by this chart.
    std::uint32_t atlasVertex{std::numeric_limits<std::uint32_t>::max()};
};
struct FrontalChartMesh {
    std::vector<FrontalVertex> vertices;
    std::vector<std::array<std::uint32_t, 3>> triangles;
    std::size_t insertedInteriorVertices{}, candidateEvaluations{}, qualityFlips{};
};

// Boundary samples are shared atlas constraints, never moved or split locally.
// The front proposes off-centres in the physical surface metric; CGAL supplies
// robust constrained Delaunay insertion. A final bounded physical-quality
// edge optimization can leave the global UV Euclidean-Delaunay configuration.
FrontalChartMesh meshFrontalChart(const SurfaceChart &chart, std::span<const Vec3> atlasVertices,
                                  const MeshingOptions &options, std::stop_token stop = {});
} // namespace hornsculpt::detail
