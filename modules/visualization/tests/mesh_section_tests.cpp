#include <hornsculpt/draft_geometry.hpp>
#include <hornsculpt/mesh_section_view.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
bool near(QPointF a, QPointF b) {
    return std::hypot(a.x() - b.x(), a.y() - b.y()) < 1e-11;
}
bool contains(const std::vector<hornsculpt::SectionSegment> &segments, QPointF point) {
    return std::any_of(segments.begin(), segments.end(), [&](const auto &s) {
        return near(s.first, point) || near(s.second, point);
    });
}
} // namespace

int main() {
    using namespace hornsculpt;
    try {
        SurfaceMesh plane;
        plane.vertices = {{-1, 0, 0}, {2, 0, 0}, {2, 0, 3}, {-1, 0, 3}};
        plane.triangles = {{{0, 1, 2}, 1, BoundaryCondition::wall},
                           {{0, 2, 3}, 1, BoundaryCondition::wall}};
        const auto outline = meshSection(plane, SectionPlane::xz);
        require(outline.size() == 4, "Coplanar patch must omit its tessellation diagonal.");
        require(contains(outline, {0, -1}) && contains(outline, {0, 2}) &&
                    contains(outline, {3, -1}) && contains(outline, {3, 2}),
                "XZ projection changed coordinates or mirrored an asymmetric section.");
        SurfaceMesh crossing;
        crossing.vertices = {{-1, -2, 0}, {1, 2, 0}, {1, -2, 2}};
        crossing.triangles = {{{0, 1, 2}, 5, BoundaryCondition::driven}};
        const auto xz = meshSection(crossing, SectionPlane::xz);
        require(xz.size() == 1 && contains(xz, {0, 0}) && contains(xz, {1, 1}),
                "Crossing-edge interpolation in the horizontal plane is incorrect.");
        require(xz.front().boundary == BoundaryCondition::driven,
                "Section lost the actual driven condition.");
        const auto yz = meshSection(crossing, SectionPlane::yz);
        require(yz.size() == 1 && contains(yz, {0, 0}) && contains(yz, {1, -2}),
                "Vertical plane must project transverse Y, not X.");
        const DraftGeometry geometry;
        const auto preview = geometry.previewMesh();
        for (const auto planeId : {SectionPlane::xz, SectionPlane::yz}) {
            const auto sections = meshSection(preview, planeId);
            const double phi = planeId == SectionPlane::xz ? 0 : std::acos(-1.0) / 2;
            require(
                contains(sections, {geometry.length(), geometry.radius(geometry.length(), phi)}),
                "Section is not taken from the actual final mouth geometry.");
            double throatMinimum = 1, throatMaximum = -1;
            bool driven = false;
            for (const auto &segment : sections)
                if (segment.boundary == BoundaryCondition::driven) {
                    driven = true;
                    for (const auto p : {segment.first, segment.second}) {
                        require(std::abs(p.x()) < 1e-12, "Driven throat section moved axially.");
                        throatMinimum = std::min(throatMinimum, p.y());
                        throatMaximum = std::max(throatMaximum, p.y());
                    }
                }
            require(driven && std::abs(throatMinimum + geometry.radius(0, phi)) < 1e-11 &&
                        std::abs(throatMaximum - geometry.radius(0, phi)) < 1e-11,
                    "The two actual sides of the driven throat must remain visible.");
        }
        std::puts("Mesh section checks passed.");
    } catch (const std::exception &error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
