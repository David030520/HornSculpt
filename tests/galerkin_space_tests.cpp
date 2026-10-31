#include <galerkin/space.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

using namespace hornsculpt;
using namespace hornsculpt::galerkin;

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
void near(double actual, double expected, double tolerance, const char *message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
        throw std::runtime_error(message);
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 scaled(Vec3 a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}
double norm(Vec3 a) {
    return std::hypot(a.x, a.y, a.z);
}
Vec3 rotate(Vec3 a) {
    // Proper 3D rotation around a non-coordinate axis, independent of cache code.
    const Vec3 u{1.0 / std::sqrt(14.0), 2.0 / std::sqrt(14.0), 3.0 / std::sqrt(14.0)};
    const double angle = 0.731, c = std::cos(angle), s = std::sin(angle);
    const Vec3 cross{u.y * a.z - u.z * a.y, u.z * a.x - u.x * a.z, u.x * a.y - u.y * a.x};
    return add(add(scaled(a, c), scaled(cross, s)), scaled(u, dot(u, a) * (1 - c)));
}
template <class Function> void rejects(Function function, const char *message) {
    bool rejected = false;
    try {
        function();
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, message);
}

SurfaceMesh tetrahedron() {
    SurfaceMesh mesh;
    mesh.vertices = {{0.1, 0.2, 0.3}, {0.4, 0.2, 0.3}, {0.1, 0.6, 0.3}, {0.1, 0.2, 0.8}};
    mesh.triangles = {{{0, 2, 1}, 7, BoundaryCondition::wall},
                      {{0, 1, 3}, 8, BoundaryCondition::driven},
                      {{0, 3, 2}, 9, BoundaryCondition::wall},
                      {{1, 2, 3}, 10, BoundaryCondition::wall}};
    return mesh;
}

void affineReferenceAndMass() {
    SurfaceMesh mesh;
    mesh.vertices = {{0, 0, 0}, {2, 0, 0}, {0, 3, 0}};
    mesh.triangles = {{{0, 1, 2}, 42, BoundaryCondition::driven}};
    const auto space = buildSpace(mesh);
    const auto &p = space.panels.front();
    near(p.area, 3.0, 1e-14, "right triangle area");
    near(norm(subtract(p.normal, {0, 0, 1})), 0.0, 1e-14, "right triangle normal");
    const std::array<Vec3, 3> expected{{{-0.5, -1.0 / 3.0, 0}, {0.5, 0, 0}, {0, 1.0 / 3.0, 0}}};
    for (std::size_t i = 0; i < 3; ++i) {
        near(norm(subtract(p.gradients[i], expected[i])), 0.0, 1e-14,
             "independent barycentric gradients");
        near(dot(p.gradients[i], p.normal), 0.0, 1e-14, "tangent gradient");
        near(dot(p.curls[i], p.normal), 0.0, 1e-14, "tangent surface curl");
        near(p.massP1DP0[i], 1.0, 1e-14, "mixed P1 DP0 integral");
        double row = 0;
        for (std::size_t j = 0; j < 3; ++j) {
            near(p.massP1P1[3 * i + j], i == j ? 0.5 : 0.25, 1e-14, "consistent local mass");
            row += p.massP1P1[3 * i + j];
        }
        near(row, p.massP1DP0[i], 1e-14, "partition of unity mass row");
    }
    near(norm(add(add(p.gradients[0], p.gradients[1]), p.gradients[2])), 0, 1e-14,
         "gradient partition");
    near(norm(add(add(p.curls[0], p.curls[1]), p.curls[2])), 0, 1e-14, "curl partition");
    near(norm(subtract(p.centroid, {2.0 / 3.0, 1, 0})), 0, 1e-14, "triangle centroid");
    require(p.faceId == 42 && p.boundary == BoundaryCondition::driven,
            "physical face metadata preserved");
    const std::array<double, 3> nodalValues{2.0, 8.0, -10.0};
    Vec3 interpolatedGradient{};
    for (std::size_t i = 0; i < 3; ++i)
        interpolatedGradient = add(interpolatedGradient, scaled(p.gradients[i], nodalValues[i]));
    near(norm(subtract(interpolatedGradient, {3, -4, 0})), 0, 1e-14, "exact affine field gradient");
}

void incidenceAndGlobalMass() {
    const auto mesh = tetrahedron();
    const auto space = buildSpace(mesh);
    require(space.vertices.size() == 4 && space.panels.size() == 4, "P1 and DP0 dimensions");
    require(space.incidenceOffsets.size() == 5 && space.incidences.size() == 12, "CSR dimensions");
    for (std::size_t vertex = 0; vertex < 4; ++vertex) {
        require(space.incidenceOffsets[vertex + 1] - space.incidenceOffsets[vertex] == 3,
                "tetra nodal stencil");
        double rowLump = 0.0, diagonal = 0.0;
        for (auto index = space.incidenceOffsets[vertex];
             index < space.incidenceOffsets[vertex + 1]; ++index) {
            const auto &incident = space.incidences[index];
            const auto &panel = space.panels[incident.panel];
            require(panel.nodes[incident.local] == vertex,
                    "incidence points back to pressure vertex");
            rowLump += panel.massP1DP0[incident.local];
            diagonal += panel.massP1P1[4 * incident.local];
        }
        near(rowLump, space.lumpedMass[vertex], 1e-15, "global lumped mass");
        near(diagonal, space.massDiagonal[vertex], 1e-15, "global mass diagonal");
    }
    near(std::accumulate(space.lumpedMass.begin(), space.lumpedMass.end(), 0.0), space.surfaceArea,
         1e-14, "partition integrates whole surface");
    near(2 * std::accumulate(space.massDiagonal.begin(), space.massDiagonal.end(), 0.0),
         space.surfaceArea, 1e-14, "consistent diagonal global integral");
}

void geometricCovariance() {
    const auto mesh = tetrahedron();
    const auto base = buildSpace(mesh);
    for (double scale : {1e-5, 0.25, 7.5, 1e4}) {
        auto transformed = mesh;
        const Vec3 translation{.7, -.9, 1.1};
        for (auto &point : transformed.vertices)
            point = add(rotate(scaled(point, scale)), translation);
        const auto actual = buildSpace(transformed);
        for (std::size_t panel = 0; panel < base.panels.size(); ++panel) {
            const auto &a = actual.panels[panel];
            const auto &b = base.panels[panel];
            const double relativeTolerance = scale < 1e-4 ? 2e-9 : 1e-12;
            near(a.area / (b.area * scale * scale), 1, relativeTolerance, "area scaling");
            near(norm(subtract(a.normal, rotate(b.normal))), 0, relativeTolerance,
                 "normal rotation");
            near(norm(subtract(a.centroid, add(rotate(scaled(b.centroid, scale)), translation))), 0,
                 std::max(1.0, scale) * 1e-12, "centroid rigid covariance");
            for (std::size_t local = 0; local < 3; ++local) {
                near(norm(subtract(scaled(a.gradients[local], scale), rotate(b.gradients[local]))),
                     0, relativeTolerance * 10, "gradient inverse length scaling");
                near(norm(subtract(scaled(a.curls[local], scale), rotate(b.curls[local]))), 0,
                     relativeTolerance * 10, "curl inverse length scaling");
            }
        }
    }
    auto reversed = mesh;
    for (auto &panel : reversed.triangles)
        std::swap(panel.vertices[1], panel.vertices[2]);
    const auto actual = buildSpace(reversed);
    for (std::size_t panel = 0; panel < base.panels.size(); ++panel) {
        near(norm(add(actual.panels[panel].normal, base.panels[panel].normal)), 0, 1e-14,
             "orientation reverses normal");
        for (std::size_t local = 0; local < 3; ++local) {
            const std::size_t original = local == 0 ? 0 : 3 - local;
            near(norm(subtract(actual.panels[panel].gradients[local],
                               base.panels[panel].gradients[original])),
                 0, 1e-13, "gradient follows basis under winding reversal");
            near(norm(add(actual.panels[panel].curls[local], base.panels[panel].curls[original])),
                 0, 1e-13, "surface curl reverses with normal");
        }
    }
}

void invalidInput() {
    rejects([] { (void)buildSpace({}); }, "empty mesh rejected");
    auto mesh = tetrahedron();
    mesh.vertices.push_back({1, 1, 1});
    rejects([&] { (void)buildSpace(mesh); }, "orphan pressure DOF rejected");
    mesh = tetrahedron();
    mesh.triangles[0].vertices[0] = 99;
    rejects([&] { (void)buildSpace(mesh); }, "invalid index rejected");
    mesh = tetrahedron();
    mesh.triangles[0].vertices[1] = mesh.triangles[0].vertices[0];
    rejects([&] { (void)buildSpace(mesh); }, "degenerate panel rejected");
    mesh = tetrahedron();
    mesh.vertices[0].x = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { (void)buildSpace(mesh); }, "nonfinite coordinates rejected");
}
} // namespace

int main() {
    try {
        affineReferenceAndMass();
        incidenceAndGlobalMass();
        geometricCovariance();
        invalidInput();
        std::cout << "Galerkin P1/DP0 spaces, mass integrals and geometric covariance passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
