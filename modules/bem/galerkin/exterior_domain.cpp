#include "exterior_domain.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace hornsculpt::galerkin {
namespace {
Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 direction(Vec3 value) {
    const double length = std::hypot(value.x, value.y, value.z);
    if (!(length > 0) || !std::isfinite(length))
        throw std::invalid_argument("Exterior components contain an invalid relative position.");
    return {value.x / length, value.y / length, value.z / length};
}

class CompensatedSum {
  public:
    void add(double value) {
        const double total = sum_ + value;
        correction_ +=
            std::abs(sum_) >= std::abs(value) ? (sum_ - total) + value : (value - total) + sum_;
        sum_ = total;
    }
    double value() const {
        return sum_ + correction_;
    }

  private:
    double sum_{};
    double correction_{};
};

struct Component {
    std::vector<std::size_t> panels;
    Vec3 lower, upper;
};

bool contains(const SurfaceMesh &mesh, const Component &component, Vec3 query) {
    if (query.x < component.lower.x || query.x > component.upper.x || query.y < component.lower.y ||
        query.y > component.upper.y || query.z < component.lower.z || query.z > component.upper.z)
        return false;
    CompensatedSum solidAngle;
    for (const auto index : component.panels) {
        const auto &triangle = mesh.triangles[index];
        const auto a = direction(subtract(mesh.vertices[triangle.vertices[0]], query));
        const auto b = direction(subtract(mesh.vertices[triangle.vertices[1]], query));
        const auto c = direction(subtract(mesh.vertices[triangle.vertices[2]], query));
        // Normalizing the rays removes the product of their lengths from both
        // atan2 arguments, avoiding overflow and preserving the solid angle.
        solidAngle.add(2 * std::atan2(dot(a, cross(b, c)), 1 + dot(a, b) + dot(b, c) + dot(c, a)));
    }
    return std::abs(solidAngle.value()) > 2 * std::numbers::pi;
}
} // namespace

void validateExteriorComponents(const SurfaceMesh &mesh) {
    if (mesh.vertices.empty() || mesh.triangles.empty())
        throw std::invalid_argument("Exterior-domain surface mesh is empty.");
    std::vector<std::size_t> parent(mesh.vertices.size());
    std::iota(parent.begin(), parent.end(), 0);
    auto root = [&parent](std::size_t vertex) {
        while (parent[vertex] != vertex) {
            parent[vertex] = parent[parent[vertex]];
            vertex = parent[vertex];
        }
        return vertex;
    };
    for (const auto &triangle : mesh.triangles) {
        for (const auto index : triangle.vertices) {
            if (index >= mesh.vertices.size())
                throw std::invalid_argument("Exterior component has an invalid vertex index.");
            const auto point = mesh.vertices[index];
            if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z))
                throw std::invalid_argument("Exterior component has nonfinite coordinates.");
        }
        const auto first = root(triangle.vertices[0]);
        for (unsigned local = 1; local < 3; ++local)
            parent[root(triangle.vertices[local])] = first;
    }
    std::vector<Component> components;
    std::map<std::size_t, std::size_t> componentIndex;
    for (std::size_t panel = 0; panel < mesh.triangles.size(); ++panel) {
        const auto &triangle = mesh.triangles[panel];
        const auto [entry, inserted] =
            componentIndex.emplace(root(triangle.vertices[0]), components.size());
        if (inserted) {
            const auto point = mesh.vertices[triangle.vertices[0]];
            components.push_back({{}, point, point});
        }
        auto &component = components[entry->second];
        component.panels.push_back(panel);
        for (const auto index : triangle.vertices) {
            const auto point = mesh.vertices[index];
            component.lower = {std::min(component.lower.x, point.x),
                               std::min(component.lower.y, point.y),
                               std::min(component.lower.z, point.z)};
            component.upper = {std::max(component.upper.x, point.x),
                               std::max(component.upper.y, point.y),
                               std::max(component.upper.z, point.z)};
        }
    }
    for (const auto &component : components) {
        CompensatedSum volume;
        for (const auto index : component.panels) {
            const auto &triangle = mesh.triangles[index];
            // A local origin is essential when separate CAD bodies have large
            // placements. A single global origin can hide an inverted body.
            const auto a = subtract(mesh.vertices[triangle.vertices[0]], component.lower);
            const auto b = subtract(mesh.vertices[triangle.vertices[1]], component.lower);
            const auto c = subtract(mesh.vertices[triangle.vertices[2]], component.lower);
            volume.add(dot(a, cross(b, c)) / 6);
        }
        if (!(volume.value() > 0) || !std::isfinite(volume.value()))
            throw std::invalid_argument(
                "Every exterior surface component must have outward winding.");
    }
    for (std::size_t inner = 0; inner < components.size(); ++inner) {
        const auto &triangle = mesh.triangles[components[inner].panels.front()];
        const auto query = mesh.vertices[triangle.vertices[0]];
        for (std::size_t outer = 0; outer < components.size(); ++outer)
            if (inner != outer && contains(mesh, components[outer], query))
                throw std::invalid_argument("Nested surface components and sealed air cavities are "
                                            "not supported by exterior BEM.");
    }
}
} // namespace hornsculpt::galerkin
