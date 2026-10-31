#include "cad_snapshot.hpp"
#include <algorithm>
#include <cmath>
#include <hornsculpt/cad_geometry.hpp>
#include <hornsculpt/meshing.hpp>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace hornsculpt {
namespace {
void require(bool condition) {
    if (!condition)
        throw std::invalid_argument("Ungültiger archivierter FreeCAD-Snapshot.");
}
bool finite(Vec3 v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}
double squared(Vec3 v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}
void string(const QString &v, bool empty = false) {
    require(v.size() <= 4096 && (empty || !v.isEmpty()));
}
void id(std::uint32_t value) {
    require(value > 0 && value < std::numeric_limits<std::uint32_t>::max());
}
Vec3 transformed(const std::array<double, 16> &m, Vec3 p, bool point) {
    return {m[0] * p.x + m[1] * p.y + m[2] * p.z + (point ? m[3] : 0),
            m[4] * p.x + m[5] * p.y + m[6] * p.z + (point ? m[7] : 0),
            m[8] * p.x + m[9] * p.y + m[10] * p.z + (point ? m[11] : 0)};
}
} // namespace

void prepareCadGeometry(CadSnapshot &s, bool buildPreview) {
    if (s.brepObjects.empty()) {
        require(!buildPreview);
        s.geometry.reset();
        return;
    }
    require(std::isfinite(s.brepMetresPerUnit) && s.brepMetresPerUnit > 0 &&
            s.brepMetresPerUnit <= 1);
    require(s.brepObjects.size() == s.objects.size());
    std::map<QString, const CadObject *> objects;
    for (const auto &object : s.objects) {
        cad_detail::validatePlacement(object.placement);
        require(objects.emplace(object.id, &object).second);
    }
    std::map<std::uint32_t, const CadFace *> faces;
    for (const auto &face : s.faces)
        require(faces.emplace(face.id, &face).second);
    std::set<QString> usedObjects;
    std::set<std::uint32_t> usedFaces;
    std::vector<CadBrepPart> parts;
    for (const auto &part : s.brepObjects) {
        require(objects.contains(part.objectId) && usedObjects.insert(part.objectId).second &&
                !part.data.isEmpty() && part.data.size() <= 48 * 1024 * 1024 &&
                !part.faceIds.empty());
        for (const auto face : part.faceIds)
            require(faces.contains(face) && faces.at(face)->objectId == part.objectId &&
                    usedFaces.insert(face).second);
        parts.push_back({part.data.toStdString(), part.faceIds});
    }
    require(usedFaces.size() == faces.size());
    auto geometry = CadGeometry::fromBrepParts(parts, s.brepMetresPerUnit);
    require(geometry && geometry->faces().size() == faces.size());
    for (const auto &native : geometry->faces()) {
        require(faces.contains(native.faceId));
        const auto &face = *faces.at(native.faceId);
        const auto &signature = face.signature;
        const auto &placement = objects.at(face.objectId)->placement;
        const auto center = transformed(placement, signature.center, true);
        const auto normal = transformed(placement, signature.normal, false);
        const double centerScale =
            std::max({1.0, std::abs(center.x), std::abs(center.y), std::abs(center.z)});
        // Export and native kernel versions may differ; tolerate sub-micrometre roundoff,
        // never a different face. Matching is verified in the global SI frame.
        const auto mismatch = [&](const char *quantity, double supplied, double measured,
                                  double tolerance, const char *unit) {
            std::ostringstream message;
            message << std::setprecision(17) << "FreeCAD-Fläche " << native.faceId << " ("
                    << face.label.toStdString() << "): " << quantity << " passt nicht zum BRep; "
                    << "Metadaten=" << supplied << ' ' << unit << ", BRep=" << measured << ' '
                    << unit << ", Grenze=" << tolerance << ' ' << unit << '.';
            throw std::invalid_argument(message.str());
        };
        const double areaTolerance = 1e-12 + 1e-8 * signature.area;
        if (!std::isfinite(native.area) || std::abs(native.area - signature.area) > areaTolerance)
            mismatch("Flächeninhalt", signature.area, native.area, areaTolerance, "m²");
        const double centerDistance = std::hypot(
            native.center.x - center.x, native.center.y - center.y, native.center.z - center.z);
        if (!std::isfinite(centerDistance) || centerDistance > 1e-8 * centerScale)
            mismatch("Schwerpunktabstand", 0, centerDistance, 1e-8 * centerScale, "m");
        const auto canonicalType = [](QString type) {
            type.remove("Surface");
            if (type == "Cylindrical")
                return QString("Cylinder");
            if (type == "Conical")
                return QString("Cone");
            if (type == "Spherical")
                return QString("Sphere");
            if (type == "Toroidal" || type == "Toroid")
                return QString("Torus");
            if (type == "OfLinearExtrusion" || type == "OfExtrusion")
                return QString("Extrusion");
            return type;
        };
        if (canonicalType(QString::fromStdString(native.surfaceType)) !=
            canonicalType(signature.surfaceType))
            throw std::invalid_argument("FreeCAD-Flaeche " + std::to_string(native.faceId) +
                                        ": Oberflaechentyp " + signature.surfaceType.toStdString() +
                                        " passt nicht zum BRep-Typ " + native.surfaceType + ".");
        if (native.edgeCount != signature.edgeCount)
            mismatch("Kantenanzahl", signature.edgeCount, native.edgeCount, 0, "");
        if (native.wireCount != signature.wireCount)
            mismatch("Randanzahl", signature.wireCount, native.wireCount, 0, "");
        const double normalDistance = std::hypot(
            native.normal.x - normal.x, native.normal.y - normal.y, native.normal.z - normal.z);
        if (!std::isfinite(normalDistance) || normalDistance > 1e-5)
            mismatch("Normalenabstand", 0, normalDistance, 1e-5, "");
    }
    if (buildPreview) {
        // A valid exact BRep may have intersecting preview chords near curved
        // face joins. Refine only this approximation defect, never the source
        // shape or the acceptance tests. Four bounded attempts reach 31.25 um.
        constexpr double initialPreviewChordTolerance = .00025; // metres
        constexpr unsigned maximumPreviewRefinements = 3;
        for (unsigned refinement = 0;; ++refinement) {
            const double chordTolerance =
                std::ldexp(initialPreviewChordTolerance, -static_cast<int>(refinement));
            auto mesh = geometry->previewMesh(chordTolerance);
            for (auto &triangle : mesh.triangles) {
                require(faces.contains(triangle.faceId));
                triangle.boundary = faces.at(triangle.faceId)->boundary;
            }
            const auto quality = inspectMesh(mesh);
            if (quality.watertight() && quality.signedVolume > 0) {
                s.mesh = std::move(mesh);
                break;
            }
            const bool onlyChordIntersections =
                quality.selfIntersectionPairs > 0 && quality.boundaryEdges == 0 &&
                quality.nonManifoldEdges == 0 && quality.inconsistentEdges == 0 &&
                quality.degenerateTriangles == 0 && quality.nonManifoldVertices == 0 &&
                quality.connectedComponents > 0 && quality.signedVolume > 0;
            if (onlyChordIntersections && refinement < maximumPreviewRefinements)
                continue;
            std::ostringstream message;
            message << "FreeCAD-Vorschaunetz ist nicht geschlossen und orientiert: "
                    << "Randkanten=" << quality.boundaryEdges
                    << ", nichtmanifold=" << quality.nonManifoldEdges
                    << ", entartet=" << quality.degenerateTriangles
                    << ", Schnitte=" << quality.selfIntersectionPairs
                    << ", Volumen=" << quality.signedVolume
                    << " m³, Sehnenabweichung=" << chordTolerance << " m.";
            throw std::invalid_argument(message.str());
        }
    }
    s.geometry = std::move(geometry);
}

std::shared_ptr<const CadGeometry> validateCadSnapshot(const CadSnapshot &s) {
    string(s.documentId);
    string(s.documentName, true);
    string(s.sessionId);
    require(s.revision > 0 && s.revision < std::numeric_limits<std::uint64_t>::max());
    require(!s.objects.empty() && s.objects.size() <= 10000 && !s.faces.empty() &&
            s.faces.size() <= 100000);
    require(finite(s.acousticOrigin) && finite(s.acousticAxis) &&
            std::abs(squared(s.acousticAxis) - 1.0) <= 2e-6);
    std::set<QString> objects;
    for (const auto &o : s.objects) {
        string(o.id);
        string(o.name);
        string(o.label, true);
        string(o.topologyFingerprint);
        require(objects.insert(o.id).second);
        for (double v : o.placement)
            require(std::isfinite(v));
        cad_detail::validatePlacement(o.placement);
    }
    std::map<std::uint32_t, BoundaryCondition> faces;
    for (const auto &f : s.faces) {
        id(f.id);
        string(f.objectId);
        string(f.sourceKey, true);
        string(f.label, true);
        require(objects.contains(f.objectId));
        require(f.boundary == BoundaryCondition::wall || f.boundary == BoundaryCondition::driven);
        require(f.state == CadFaceState::confirmed || f.state == CadFaceState::newFace ||
                f.state == CadFaceState::needsReview);
        require(std::isfinite(f.normalVelocity.real()) && std::isfinite(f.normalVelocity.imag()));
        require(faces.emplace(f.id, f.boundary).second);
        const auto &g = f.signature;
        string(g.surfaceType);
        string(g.geometryKey);
        string(g.topologyKey);
        require(std::isfinite(g.area) && g.area > 0 && finite(g.center) && finite(g.normal) &&
                squared(g.normal) <= (1 + 1e-6) * (1 + 1e-6) && g.edgeCount <= 100000 &&
                g.wireCount <= 100000);
    }
    require(s.mesh.vertices.size() >= 3 && s.mesh.vertices.size() <= 1000000 &&
            !s.mesh.triangles.empty() && s.mesh.triangles.size() <= 2000000);
    const auto quality = inspectMesh(s.mesh);
    require(quality.degenerateTriangles == 0);
    std::set<std::uint32_t> used;
    for (const auto &t : s.mesh.triangles) {
        const auto found = faces.find(t.faceId);
        require(found != faces.end() && found->second == t.boundary);
        used.insert(t.faceId);
    }
    require(used.size() == faces.size());
    std::set<std::uint32_t> missing;
    for (auto previous : s.unresolvedPreviousFaces) {
        id(previous);
        require(!faces.contains(previous) && missing.insert(previous).second);
    }
    if (!s.brepObjects.empty()) {
        auto checked = s;
        prepareCadGeometry(checked, false);
        return checked.geometry;
    }
    require(!s.geometry);
    return {};
}

MeshingSource makeCadMeshingSource(const CadSnapshot &snapshot, double geometricTolerance,
                                   std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("CAD-Vernetzung abgebrochen.");
    if (!std::isfinite(geometricTolerance) || geometricTolerance <= 0)
        throw std::invalid_argument("CAD-Detailabweichung muss endlich und positiv sein.");
    if (!snapshot.geometry && !snapshot.brepObjects.empty())
        throw std::invalid_argument("Exakte CAD-Netzquelle wurde nicht vorbereitet.");
    // Native source data and UV bindings are produced together, independently
    // of the archived/display preview's triangulation or detail tolerance.
    auto source = snapshot.geometry ? snapshot.geometry->meshingSource(geometricTolerance, stop)
                                    : MeshingSource{snapshot.mesh, {}};
    std::map<std::uint32_t, BoundaryCondition> conditions;
    for (const auto &face : snapshot.faces)
        if (!conditions.emplace(face.id, face.boundary).second)
            throw std::invalid_argument("CAD-Snapshot enthält doppelte Flächen-IDs.");
    std::set<std::uint32_t> used;
    for (auto &triangle : source.mesh.triangles) {
        if (stop.stop_requested())
            throw std::runtime_error("CAD-Vernetzung abgebrochen.");
        const auto found = conditions.find(triangle.faceId);
        if (found == conditions.end())
            throw std::invalid_argument("CAD-Netzquelle enthält eine unbekannte Fläche.");
        triangle.boundary = found->second;
        used.insert(triangle.faceId);
    }
    if (used.size() != conditions.size())
        throw std::invalid_argument("CAD-Netzquelle enthält nicht alle zugeordneten Flächen.");
    return source;
}
} // namespace hornsculpt
