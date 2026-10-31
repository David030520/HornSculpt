#include "project_assembly.hpp"
#include "galerkin/exterior_domain.hpp"
#include "project_json.hpp"
#include <QUuid>
#include <algorithm>
#include <hornsculpt/cad_geometry.hpp>
#include <hornsculpt/mesh_orientation.hpp>
#include <stdexcept>

namespace hornsculpt {
std::shared_ptr<const CadSnapshot> appendAssemblyCad(const CadSnapshot &current,
                                                     const CadSnapshot &added) {
    auto result = std::make_shared<CadSnapshot>(current);
    result->documentId = QUuid::createUuid().toString(QUuid::WithoutBraces);
    result->documentName = current.documentName + QStringLiteral(" + ") + added.documentName;
    result->sessionId = QStringLiteral("assembly-import");
    result->revision = 1;
    result->brepMetresPerUnit = 1;
    result->brepObjects.clear();
    const auto currentGeometry = validateCadSnapshot(current);
    const auto addedGeometry = validateCadSnapshot(added);
    if (!currentGeometry || !addedGeometry)
        throw std::invalid_argument("Baugruppenimport benötigt exakte CAD-Körper.");
    std::uint32_t offset = 0;
    for (const auto &f : current.faces)
        offset = std::max(offset, f.id);
    for (const auto &f : current.unresolvedPreviousFaces)
        offset = std::max(offset, f);
    const auto prefix = result->documentId + QStringLiteral(":");
    for (auto object : added.objects) {
        object.id = prefix + object.id;
        result->objects.push_back(std::move(object));
    }
    for (auto face : added.faces) {
        if (face.id >= assemblyCadFaceOffset - offset)
            throw std::invalid_argument("Zu viele CAD-Flächen in der Baugruppe.");
        face.id += offset;
        face.objectId = prefix + face.objectId;
        result->faces.push_back(std::move(face));
    }
    for (auto id : added.unresolvedPreviousFaces)
        result->unresolvedPreviousFaces.push_back(id + offset);
    for (const auto &part : currentGeometry->brepParts()) {
        const auto face = std::find_if(current.faces.begin(), current.faces.end(),
                                       [&](const auto &f) { return f.id == part.faceIds.front(); });
        result->brepObjects.push_back(
            {face->objectId, QByteArray::fromStdString(part.data), part.faceIds});
    }
    for (const auto &part : addedGeometry->brepParts()) {
        const auto face = std::find_if(added.faces.begin(), added.faces.end(),
                                       [&](const auto &f) { return f.id == part.faceIds.front(); });
        auto ids = part.faceIds;
        for (auto &id : ids)
            id += offset;
        result->brepObjects.push_back(
            {prefix + face->objectId, QByteArray::fromStdString(part.data), std::move(ids)});
    }
    prepareCadGeometry(*result, true);
    (void)validateCadSnapshot(*result);
    return result;
}
AcousticFrame assemblyAcousticFrame(const DraftAssembly &assembly, double length) {
    validateDraftAssembly(assembly);
    Vec3 origin{}, forward{};
    for (const auto &i : assembly.instances) {
        const auto p = instancePoint(i, {0, 0, length}, length);
        origin.x += p.x;
        origin.y += p.y;
        origin.z += p.z;
        forward.x += i.axes[2].x;
        forward.y += i.axes[2].y;
        forward.z += i.axes[2].z;
    }
    const auto n = static_cast<double>(assembly.instances.size());
    return defaultAcousticFrame({origin.x / n, origin.y / n, origin.z / n}, forward);
}
namespace project_json {
QJsonObject encodeAssembly(const DraftAssembly &a) {
    validateDraftAssembly(a);
    QJsonArray instances;
    for (const auto &i : a.instances) {
        QJsonArray axes;
        for (const auto v : i.axes)
            axes.append(encodeVector(v));
        instances.append(QJsonObject{{"id", static_cast<double>(i.id)},
                                     {"name", QString::fromStdString(i.name)},
                                     {"position", encodeVector(i.position)},
                                     {"axes", axes},
                                     {"mirrorHorizontal", i.mirrorHorizontal},
                                     {"mouthAnchored", i.mouthAnchored}});
    }
    return {{"instances", instances},
            {"includeCad", a.includeCad},
            {"cadTargetEdgeLength", a.cadTargetEdgeLength}};
}
DraftAssembly decodeAssembly(const QJsonValue &value) {
    const auto o = object(value);
    keys(o, {"instances", "includeCad", "cadTargetEdgeLength"});
    DraftAssembly a;
    a.instances.clear();
    a.includeCad = boolean(o["includeCad"]);
    a.cadTargetEdgeLength = positive(o["cadTargetEdgeLength"]);
    for (const auto v : array(o["instances"], 64, 1)) {
        const auto i = object(v);
        keys(i, {"id", "name", "position", "axes", "mirrorHorizontal", "mouthAnchored"});
        DraftInstance instance;
        instance.id = integer(i["id"], 1000000);
        instance.name = text(i["name"]).toStdString();
        instance.position = vector(i["position"]);
        const auto axes = array(i["axes"], 3, 3);
        for (int k = 0; k < 3; k++)
            instance.axes[k] = vector(axes[k]);
        instance.mirrorHorizontal = boolean(i["mirrorHorizontal"]);
        instance.mouthAnchored = boolean(i["mouthAnchored"]);
        a.instances.push_back(std::move(instance));
    }
    validateDraftAssembly(a);
    return a;
}
} // namespace project_json

SurfaceMesh assemblyPreview(const SurfaceMesh &local, const Project &p) {
    auto draft = applyFaceOrientation(local, p.draftReversedFaceIds);
    std::optional<SurfaceMesh> cad;
    if (p.assembly.includeCad) {
        if (!p.cad)
            throw std::invalid_argument("Baugruppe: CAD-Geometrie fehlt.");
        cad = applyFaceOrientation(p.cad->mesh, p.cadReversedFaceIds);
    }
    return composeDraftInstances(draft, p.assembly, p.parameters.length, cad ? &*cad : nullptr);
}

MeshingResult meshDraftAssembly(const Project &p, const SurfaceMesh &local, std::stop_token stop) {
    MeshingOptions options;
    options.targetEdgeLength = p.meshing.targetEdgeLength;
    options.adaptive = p.meshing.adaptive;
    options.quality = p.meshing.quality;
    options.faceTargetEdgeLengths = p.draftFaceTargetEdgeLengths;
    auto draft = remeshIsotropic(local, options, stop).mesh;
    draft = applyFaceOrientation(draft, p.draftReversedFaceIds);
    std::optional<SurfaceMesh> cad;
    if (p.assembly.includeCad) {
        if (!p.cad || p.cad->hasUnresolvedAssignments())
            throw std::invalid_argument("Baugruppe: CAD-Flächen fehlen oder sind ungeklärt.");
        auto source = makeCadMeshingSource(*p.cad, options.adaptive.geometricTolerance, stop);
        options.targetEdgeLength = p.assembly.cadTargetEdgeLength;
        options.faceTargetEdgeLengths = p.cadFaceTargetEdgeLengths;
        options.geometry = p.cad->geometry;
        options.parameterization = source.parameters;
        cad = remeshIsotropic(source.mesh, options, stop).mesh;
        *cad = applyFaceOrientation(*cad, p.cadReversedFaceIds);
    }
    MeshingResult result;
    result.mesh =
        composeDraftInstances(draft, p.assembly, p.parameters.length, cad ? &*cad : nullptr);
    result.quality = inspectMesh(result.mesh);
    if (result.quality.inconsistentEdges &&
        (!p.draftReversedFaceIds.empty() || !p.cadReversedFaceIds.empty()))
        throw std::invalid_argument(
            "Flächennormalen: Flächenauswahl ergibt keine konsistent orientierte geschlossene Hülle.");
    if (!result.quality.watertight() || result.quality.selfIntersectionPairs)
        throw std::invalid_argument(
            "Baugruppe: Körper berühren/schneiden sich oder sind nicht geschlossen.");
    galerkin::validateExteriorComponents(result.mesh);
    return result;
}
} // namespace hornsculpt
