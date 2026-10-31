#include <algorithm>
#include <cmath>
#include <hornsculpt/draft_instances.hpp>
#include <limits>
#include <set>
#include <stdexcept>

namespace hornsculpt {
namespace {
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
bool finite(Vec3 p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
} // namespace
void validateDraftAssembly(const DraftAssembly &a) {
    if (a.instances.empty() || a.instances.size() > 64 || !std::isfinite(a.cadTargetEdgeLength) ||
        a.cadTargetEdgeLength < .00005 || a.cadTargetEdgeLength > 2)
        throw std::invalid_argument(
            "Baugruppe: 1–64 Hörner und CAD-Zielkante 0,05–2000 mm erforderlich.");
    std::set<std::uint32_t> ids;
    for (const auto &i : a.instances) {
        if (!i.id || i.id > 1000000 || !ids.insert(i.id).second || i.name.empty() ||
            !finite(i.position))
            throw std::invalid_argument("Ungültige oder doppelte Horninstanz.");
        for (const auto v : i.axes)
            if (!finite(v) || std::abs(dot(v, v) - 1) > 1e-10)
                throw std::invalid_argument("Instanzachsen müssen normiert sein.");
        if (std::abs(dot(i.axes[0], i.axes[1])) > 1e-10 ||
            std::abs(dot(i.axes[0], i.axes[2])) > 1e-10 ||
            std::abs(dot(i.axes[1], i.axes[2])) > 1e-10 ||
            std::abs(dot(cross(i.axes[0], i.axes[1]), i.axes[2]) - 1) > 1e-10)
            throw std::invalid_argument("Instanzrahmen muss rechtshändig und orthogonal sein.");
    }
}
bool isPlacedAssembly(const DraftAssembly &a) {
    return a.includeCad || a.instances != DraftAssembly{}.instances;
}
std::uint32_t instanceFaceId(std::uint32_t instance, std::uint32_t face) {
    if (!instance || instance > 1000000 || !face || face > 8)
        throw std::invalid_argument("Ungültige Draft-Flächenidentität.");
    return (instance - 1) * 8 + face;
}
Vec3 instancePoint(const DraftInstance &i, Vec3 p, double length) {
    if (i.mirrorHorizontal)
        p.x = -p.x;
    if (i.mouthAnchored)
        p.z -= length;
    return {i.position.x + i.axes[0].x * p.x + i.axes[1].x * p.y + i.axes[2].x * p.z,
            i.position.y + i.axes[0].y * p.x + i.axes[1].y * p.y + i.axes[2].y * p.z,
            i.position.z + i.axes[0].z * p.x + i.axes[1].z * p.y + i.axes[2].z * p.z};
}
SurfaceMesh composeDraftInstances(const SurfaceMesh &draft, const DraftAssembly &a, double length,
                                  const SurfaceMesh *cad) {
    validateDraftAssembly(a);
    if (a.includeCad && !cad)
        throw std::invalid_argument("CAD-Körper der Baugruppe fehlen.");
    SurfaceMesh result;
    const auto append = [&](const SurfaceMesh &s, const DraftInstance *i) {
        if (result.vertices.size() + s.vertices.size() > std::numeric_limits<std::uint32_t>::max())
            throw std::invalid_argument("Baugruppennetz ist zu groß.");
        const auto offset = static_cast<std::uint32_t>(result.vertices.size());
        for (const auto p : s.vertices)
            result.vertices.push_back(i ? instancePoint(*i, p, length) : p);
        for (auto t : s.triangles) {
            for (auto &v : t.vertices)
                v += offset;
            if (i) {
                t.faceId = instanceFaceId(i->id, t.faceId);
                if (i->mirrorHorizontal)
                    std::swap(t.vertices[1], t.vertices[2]);
            } else {
                if (!t.faceId || t.faceId >= assemblyCadFaceOffset)
                    throw std::invalid_argument("CAD-Flächen-ID außerhalb des Baugruppenbereichs.");
                t.faceId += assemblyCadFaceOffset;
            }
            result.triangles.push_back(t);
        }
        for (auto v : s.fixedVertices)
            result.fixedVertices.push_back(v + offset);
    };
    for (const auto &i : a.instances)
        append(draft, &i);
    if (a.includeCad)
        append(*cad, nullptr);
    return result;
}
} // namespace hornsculpt
