#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <hornsculpt/mesh.hpp>
#include <map>

namespace hornsculpt::testing {

inline SurfaceMesh sphere(unsigned refinements) {
    SurfaceMesh mesh;
    mesh.vertices = {{.1, 0, 0}, {0, .1, 0}, {-.1, 0, 0}, {0, -.1, 0}, {0, 0, .1}, {0, 0, -.1}};
    for (const auto triangle : std::array<std::array<std::uint32_t, 3>, 8>{{{4, 0, 1},
                                                                            {4, 1, 2},
                                                                            {4, 2, 3},
                                                                            {4, 3, 0},
                                                                            {5, 1, 0},
                                                                            {5, 2, 1},
                                                                            {5, 3, 2},
                                                                            {5, 0, 3}}})
        mesh.triangles.push_back({triangle, 1, BoundaryCondition::driven});
    for (unsigned level = 0; level < refinements; ++level) {
        auto next = mesh;
        next.triangles.clear();
        std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> edges;
        auto midpoint = [&](std::uint32_t i, std::uint32_t j) {
            const auto key = std::minmax(i, j);
            if (const auto found = edges.find(key); found != edges.end())
                return found->second;
            const auto a = mesh.vertices[i], b = mesh.vertices[j];
            Vec3 point{a.x + b.x, a.y + b.y, a.z + b.z};
            const double factor =
                .1 / std::sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
            point = {point.x * factor, point.y * factor, point.z * factor};
            const auto index = static_cast<std::uint32_t>(next.vertices.size());
            next.vertices.push_back(point);
            edges.emplace(key, index);
            return index;
        };
        for (const auto &triangle : mesh.triangles) {
            const auto a = triangle.vertices[0], b = triangle.vertices[1], c = triangle.vertices[2];
            const auto ab = midpoint(a, b), bc = midpoint(b, c), ac = midpoint(a, c);
            for (const auto vertices : std::array<std::array<std::uint32_t, 3>, 4>{
                     {{a, ab, ac}, {ab, b, bc}, {ac, bc, c}, {ab, bc, ac}}})
                next.triangles.push_back({vertices, 1, BoundaryCondition::driven});
        }
        mesh = std::move(next);
    }
    return mesh;
}
} // namespace hornsculpt::testing
