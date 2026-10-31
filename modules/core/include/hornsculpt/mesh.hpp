#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace hornsculpt {

// All positions use metres. Face IDs survive tessellation and remeshing.
struct Vec3 {
    double x{};
    double y{};
    double z{};
    bool operator==(const Vec3 &) const = default;
};

enum class BoundaryCondition : std::uint8_t { wall, driven };

struct Triangle {
    std::array<std::uint32_t, 3> vertices{};
    std::uint32_t faceId{};
    BoundaryCondition boundary{BoundaryCondition::wall};
    bool operator==(const Triangle &) const = default;
};

struct SurfaceMesh {
    std::vector<Vec3> vertices;
    std::vector<Triangle> triangles;
    // Preserve these original positions during remeshing; indices refer to vertices.
    // Physical patch boundaries connecting fixed points may still be subdivided.
    std::vector<std::uint32_t> fixedVertices;
    bool operator==(const SurfaceMesh &) const = default;
};

} // namespace hornsculpt
