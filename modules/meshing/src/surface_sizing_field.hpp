#pragma once

#include <hornsculpt/mesh.hpp>
#include <map>
#include <memory>
#include <span>
#include <stop_token>
#include <vector>

namespace hornsculpt::detail {
struct SurfaceEdgeSeed {
    std::array<std::uint32_t, 2> vertices;
    double size;
};

// Grading follows source-surface edges; sampling stays on the requested physical
// patches. In particular, nearby opposing walls do not share a distance field.
class SurfaceSizingField {
  public:
    SurfaceSizingField(const SurfaceMesh &source, std::vector<double> vertexSizes,
                       double growthRate, double targetEdgeLength,
                       std::span<const SurfaceEdgeSeed> edgeSeeds, std::stop_token stop,
                       const std::map<std::uint32_t, double> &faceTargetEdgeLengths = {});
    ~SurfaceSizingField();
    double at(Vec3 position, std::span<const std::uint32_t> faceIds) const;

  private:
    struct Data;
    std::unique_ptr<Data> data_;
};
} // namespace hornsculpt::detail
