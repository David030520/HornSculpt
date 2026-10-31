#pragma once
#include "fixture_geometry.hpp"
#include <algorithm>
#include <map>
namespace hornsculpt::testing {
// A quarter of two concentric shells, with an opening on the upper outer shell.
inline SurfaceMesh openRimSourceMesh(bool sector) {
    SurfaceMesh result;
    for (int shell = 0; shell < 2; ++shell) {
        const auto source = sphere(2);
        std::map<std::uint32_t, std::uint32_t> indices;
        for (const auto &element : source.triangles) {
            if (sector &&
                std::any_of(element.vertices.begin(), element.vertices.end(), [&](auto v) {
                    return source.vertices[v].x < -1e-12 || source.vertices[v].y < -1e-12;
                }))
                continue;
            auto triangle = element;
            double z = 0;
            for (auto &v : triangle.vertices) {
                const auto p = source.vertices[v];
                z += p.z;
                const auto [it, added] =
                    indices.emplace(v, static_cast<std::uint32_t>(result.vertices.size()));
                if (added) {
                    const double scale = shell ? 2.0 : .5;
                    result.vertices.push_back({scale * p.x, scale * p.y, scale * p.z});
                }
                v = it->second;
            }
            triangle.faceId = shell ? (z > 0 ? 2u : 3u) : 1u;
            triangle.boundary = shell ? BoundaryCondition::wall : BoundaryCondition::driven;
            if (shell)
                std::swap(triangle.vertices[1], triangle.vertices[2]);
            result.triangles.push_back(triangle);
        }
    }
    return result;
}
} // namespace hornsculpt::testing
