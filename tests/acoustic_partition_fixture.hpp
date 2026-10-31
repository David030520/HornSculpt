#pragma once
#include <cmath>
#include <hornsculpt/mesh.hpp>
#include <numbers>
#include <vector>
namespace hornsculpt::testing {
inline SurfaceMesh acousticCup() {
    constexpr unsigned sectors = 16;
    SurfaceMesh mesh;
    const auto ring = [&](double radius, double z) {
        const auto start = static_cast<std::uint32_t>(mesh.vertices.size());
        for (unsigned i = 0; i < sectors; ++i) {
            const auto angle = 2 * std::numbers::pi * i / sectors;
            mesh.vertices.push_back({radius * std::cos(angle), radius * std::sin(angle), z});
        }
        return start;
    };
    std::vector<unsigned> inner, outer;
    for (unsigned layer = 0; layer <= 4; ++layer)
        inner.push_back(ring(.04, .02 * layer));
    for (unsigned layer = 0; layer <= 5; ++layer)
        outer.push_back(ring(.05, -.02 + .02 * layer));
    const auto triangle = [&](unsigned a, unsigned b, unsigned c, unsigned face,
                              BoundaryCondition bc = BoundaryCondition::wall) {
        mesh.triangles.push_back({{a, b, c}, face, bc});
    };
    for (unsigned i = 0; i < sectors; ++i) {
        const auto j = (i + 1) % sectors;
        for (unsigned layer = 0; layer < 4; ++layer) {
            const auto a = inner[layer] + i, b = inner[layer] + j, c = inner[layer + 1] + i,
                       d = inner[layer + 1] + j;
            triangle(a, c, b, 2);
            triangle(b, c, d, 2);
        }
        for (unsigned layer = 0; layer < 5; ++layer) {
            const auto a = outer[layer] + i, b = outer[layer] + j, c = outer[layer + 1] + i,
                       d = outer[layer + 1] + j;
            triangle(a, b, c, 1);
            triangle(b, d, c, 1);
        }
        triangle(inner.back() + i, outer.back() + i, inner.back() + j, 4);
        triangle(inner.back() + j, outer.back() + i, outer.back() + j, 4);
    }
    const auto floor = static_cast<unsigned>(mesh.vertices.size());
    mesh.vertices.push_back({0, 0, 0});
    const auto base = static_cast<unsigned>(mesh.vertices.size());
    mesh.vertices.push_back({0, 0, -.02});
    for (unsigned i = 0; i < sectors; ++i) {
        const auto j = (i + 1) % sectors;
        triangle(floor, inner.front() + i, inner.front() + j, 3, BoundaryCondition::driven);
        triangle(base, outer.front() + j, outer.front() + i, 5);
    }
    return mesh;
}
} // namespace hornsculpt::testing
