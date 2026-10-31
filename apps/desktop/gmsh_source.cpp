#include "gmsh_source.hpp"
#include "galerkin/exterior_domain.hpp"
#include <QCryptographicHash>
#include <QDataStream>
#include <QFileInfo>
#include <QJsonArray>
#include <QStringList>
#include <QUuid>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace hornsculpt {
MeshingResult validateGmshMesh(const SurfaceMesh &mesh,
                               std::span<const AcousticRegionDefinition> regions,
                               const MeshingOptions &options, std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Netzprüfung abgebrochen.");
    const auto quality = inspectMesh(mesh);
    if (quality.degenerateTriangles)
        throw std::invalid_argument("Importnetz enthält entartete Dreiecke.");
    if (options.quality.enabled &&
        (quality.minimumAngleDegrees < options.quality.minimumAngleDegrees ||
         quality.maximumAspectRatio > options.quality.maximumAspectRatio))
        throw std::invalid_argument(
            "Importnetz verletzt die eingestellten Winkel- oder Seitenverhältnisgrenzen.");
    if (regions.empty()) {
        if (!quality.watertight() || quality.selfIntersectionPairs)
            throw std::invalid_argument(
                "Importnetz ist nicht geschlossen, konsistent orientiert und überschneidungsfrei.");
        galerkin::validateExteriorComponents(mesh);
    } else {
        (void)partitionAcousticSurface(mesh, regions, options, stop);
    }
    return {mesh, quality};
}
namespace {
[[noreturn]] void invalid() {
    throw std::invalid_argument("Ungültige Gmsh-Importmetadaten.");
}
QString idString(std::uint64_t n) {
    return QString::number(static_cast<qulonglong>(n));
}
std::uint64_t parseId(const QJsonValue &value, bool allowZero = false) {
    if (!value.isString())
        invalid();
    bool ok = false;
    const auto text = value.toString();
    auto n = text.toULongLong(&ok);
    if (!ok || (!allowZero && !n) || idString(n) != text)
        invalid();
    return n;
}
std::uint32_t id32(const QJsonValue &value, bool allowZero = false) {
    auto n = parseId(value, allowZero);
    if (n > std::numeric_limits<std::uint32_t>::max())
        invalid();
    return static_cast<std::uint32_t>(n);
}
Vec3 minus(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
void validatePhysicalGroups(const std::vector<GmshTriangleMetadata> &triangles) {
    std::map<std::uint32_t, std::vector<std::uint32_t>> groups;
    for (const auto &triangle : triangles) {
        auto tags = triangle.physicalTags;
        std::sort(tags.begin(), tags.end());
        const auto [it, inserted] = groups.emplace(triangle.entityTag, std::move(tags));
        if (!inserted) {
            auto expected = triangle.physicalTags;
            std::sort(expected.begin(), expected.end());
            if (it->second != expected)
                throw std::invalid_argument("Gmsh: Eine Entity hat unterschiedliche "
                                            "Physical-Gruppen; getrennte Flächen exportieren.");
        }
    }
}
struct FaceData {
    double area{};
    Vec3 center{}, normal{};
    QByteArray content;
    std::set<std::uint32_t> physical;
};
QString meshDigest(const SurfaceMesh &mesh) {
    QByteArray bytes;
    QDataStream data(&bytes, QIODevice::WriteOnly);
    data.setVersion(QDataStream::Qt_6_0);
    data.setByteOrder(QDataStream::BigEndian);
    data << static_cast<quint64>(mesh.vertices.size())
         << static_cast<quint64>(mesh.triangles.size());
    for (const auto &p : mesh.vertices)
        data << p.x << p.y << p.z;
    for (const auto &t : mesh.triangles)
        data << t.vertices[0] << t.vertices[1] << t.vertices[2] << t.faceId;
    return QString::fromLatin1(QCryptographicHash::hash(bytes, QCryptographicHash::Sha256).toHex());
}
} // namespace
std::shared_ptr<const CadSnapshot> makeGmshCadSnapshot(const MeshImportSelection &selection,
                                                       const QString &filename) {
    const auto &input = selection.imported;
    validatePhysicalGroups(input.triangleMetadata);
    if (input.mesh.vertices.size() != input.nodeTags.size() ||
        input.mesh.triangles.size() != input.triangleMetadata.size())
        invalid();
    auto result = std::make_shared<CadSnapshot>();
    result->documentId = QUuid::createUuid().toString(QUuid::WithoutBraces);
    result->documentName = QFileInfo(filename).fileName();
    result->sessionId = QStringLiteral("gmsh-import-offline");
    result->revision = 1;
    result->brepMetresPerUnit = 1;
    result->mesh = input.mesh;
    std::map<std::uint32_t, FaceData> faces;
    for (std::size_t i = 0; i < input.mesh.triangles.size(); ++i) {
        const auto &t = input.mesh.triangles[i];
        const auto &meta = input.triangleMetadata[i];
        if (t.faceId != meta.entityTag || t.boundary != BoundaryCondition::wall)
            invalid();
        for (auto v : t.vertices)
            if (v >= input.mesh.vertices.size())
                invalid();
        auto a = input.mesh.vertices[t.vertices[0]], b = input.mesh.vertices[t.vertices[1]],
             c = input.mesh.vertices[t.vertices[2]];
        auto n = cross(minus(b, a), minus(c, a));
        const double twice = std::hypot(n.x, n.y, n.z);
        if (!(twice > 0) || !std::isfinite(twice))
            invalid();
        auto &f = faces[t.faceId];
        const double area = twice / 2;
        f.area += area;
        f.center.x += area * (a.x + b.x + c.x) / 3;
        f.center.y += area * (a.y + b.y + c.y) / 3;
        f.center.z += area * (a.z + b.z + c.z) / 3;
        f.normal.x += n.x / 2;
        f.normal.y += n.y / 2;
        f.normal.z += n.z / 2;
        f.physical.insert(meta.physicalTags.begin(), meta.physicalTags.end());
        QDataStream data(&f.content, QIODevice::Append);
        data.setVersion(QDataStream::Qt_6_0);
        data.setByteOrder(QDataStream::BigEndian);
        data << static_cast<quint64>(meta.elementTag);
        for (auto v : t.vertices) {
            const auto p = input.mesh.vertices[v];
            data << static_cast<quint64>(input.nodeTags[v]) << p.x << p.y << p.z;
        }
    }
    if (faces.empty())
        invalid();
    QByteArray combined;
    for (const auto &[id, f] : faces) {
        (void)id;
        combined += f.content;
    }
    auto hash = [](const QByteArray &bytes) {
        return QString::fromLatin1(
            QCryptographicHash::hash(bytes, QCryptographicHash::Sha256).toHex());
    };
    CadObject object;
    object.id = QStringLiteral("gmsh-object");
    object.name = object.id;
    object.label = result->documentName;
    object.topologyFingerprint = hash(combined);
    object.placement = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    result->objects.push_back(object);
    for (const auto &[id, f] : faces) {
        CadFace face;
        face.id = id;
        face.objectId = object.id;
        face.sourceKey = object.topologyFingerprint + QStringLiteral(":%1").arg(id);
        QStringList names;
        for (const auto &group : input.physicalGroups)
            if (group.dimension == 2 && f.physical.contains(group.tag))
                names << QString::fromUtf8(group.name.data(), static_cast<int>(group.name.size()));
        face.label = names.isEmpty()
                         ? QStringLiteral("Gmsh Fläche %1").arg(id)
                         : QStringLiteral("Gmsh Fläche %1 · %2").arg(id).arg(names.join(", "));
        face.state = CadFaceState::confirmed;
        face.signature.surfaceType = QStringLiteral("polygonal");
        face.signature.geometryKey = hash(f.content);
        face.signature.topologyKey = face.signature.geometryKey;
        face.signature.area = f.area;
        face.signature.center = {f.center.x / f.area, f.center.y / f.area, f.center.z / f.area};
        const double n = std::hypot(f.normal.x, f.normal.y, f.normal.z);
        face.signature.normal =
            n > 0 ? Vec3{f.normal.x / n, f.normal.y / n, f.normal.z / n} : Vec3{};
        result->faces.push_back(std::move(face));
    }
    (void)validateCadSnapshot(*result);
    return result;
}
QJsonObject encodeGmshMetadata(const MeshImportSelection &selection) {
    const auto &m = selection.imported;
    if (!(selection.lengthScaleToMetres > 0) || !std::isfinite(selection.lengthScaleToMetres) ||
        m.nodeTags.size() != m.mesh.vertices.size() ||
        m.triangleMetadata.size() != m.mesh.triangles.size())
        invalid();
    QJsonArray nodes, triangles, groups;
    for (auto tag : m.nodeTags)
        nodes.append(idString(tag));
    for (const auto &t : m.triangleMetadata) {
        QJsonArray tags;
        for (auto tag : t.physicalTags)
            tags.append(idString(tag));
        triangles.append(QJsonObject{{"element", idString(t.elementTag)},
                                     {"entity", idString(t.entityTag)},
                                     {"physical", tags}});
    }
    for (const auto &g : m.physicalGroups)
        groups.append(QJsonObject{
            {"dimension", g.dimension},
            {"tag", idString(g.tag)},
            {"name", QString::fromUtf8(g.name.data(), static_cast<int>(g.name.size()))}});
    const auto &d = m.diagnostics;
    QJsonObject diagnostics{{"skippedPoints", idString(d.skippedPoints)},
                            {"skippedLines", idString(d.skippedLines)},
                            {"boundaryEdges", idString(d.boundaryEdges)},
                            {"nonManifoldEdges", idString(d.nonManifoldEdges)},
                            {"inconsistentEdges", idString(d.inconsistentEdges)},
                            {"degenerateTriangles", idString(d.degenerateTriangles)}};
    return {{"version", 1},
            {"meshDigest", meshDigest(m.mesh)},
            {"scale", selection.lengthScaleToMetres},
            {"nodes", nodes},
            {"triangles", triangles},
            {"groups", groups},
            {"diagnostics", diagnostics}};
}
MeshImportSelection decodeGmshMetadata(const QJsonObject &json, const SurfaceMesh &mesh) {
    if (json.value("version").toInt() != 1 || !json.value("scale").isDouble() ||
        !json.value("nodes").isArray() || !json.value("triangles").isArray() ||
        !json.value("groups").isArray() || !json.value("diagnostics").isObject())
        invalid();
    if (!json.value("meshDigest").isString() ||
        json.value("meshDigest").toString() != meshDigest(mesh))
        invalid();
    MeshImportSelection result;
    result.lengthScaleToMetres = json.value("scale").toDouble();
    if (!(result.lengthScaleToMetres > 0) || !std::isfinite(result.lengthScaleToMetres))
        invalid();
    auto &m = result.imported;
    m.mesh = mesh;
    std::set<std::uint64_t> seen;
    for (const auto &v : json.value("nodes").toArray()) {
        auto id = parseId(v);
        if (!seen.insert(id).second)
            invalid();
        m.nodeTags.push_back(id);
    }
    if (m.nodeTags.size() != mesh.vertices.size())
        invalid();
    auto triangles = json.value("triangles").toArray();
    if (static_cast<std::size_t>(triangles.size()) != mesh.triangles.size())
        invalid();
    for (qsizetype i = 0; i < triangles.size(); ++i) {
        if (!triangles[i].isObject())
            invalid();
        auto o = triangles[i].toObject();
        GmshTriangleMetadata t;
        t.elementTag = parseId(o.value("element"));
        t.entityTag = id32(o.value("entity"));
        if (t.entityTag != mesh.triangles[static_cast<std::size_t>(i)].faceId ||
            !o.value("physical").isArray())
            invalid();
        std::set<std::uint32_t> uniquePhysical;
        for (const auto &v : o.value("physical").toArray()) {
            auto tag = id32(v);
            if (!uniquePhysical.insert(tag).second)
                invalid();
            t.physicalTags.push_back(tag);
        }
        m.triangleMetadata.push_back(std::move(t));
    }
    std::set<std::pair<int, std::uint32_t>> groups;
    for (const auto &v : json.value("groups").toArray()) {
        if (!v.isObject())
            invalid();
        auto o = v.toObject();
        if (!o.value("dimension").isDouble() || !o.value("name").isString())
            invalid();
        int dim = o.value("dimension").toInt(-1);
        auto tag = id32(o.value("tag"));
        if (dim < 0 || dim > 3 || !groups.emplace(dim, tag).second)
            invalid();
        m.physicalGroups.push_back({dim, tag, o.value("name").toString().toUtf8().toStdString()});
    }
    validatePhysicalGroups(m.triangleMetadata);
    const auto d = json.value("diagnostics").toObject();
    auto count = [&](const char *key) {
        auto n = parseId(d.value(key), true);
        if (n > std::numeric_limits<std::size_t>::max())
            invalid();
        return static_cast<std::size_t>(n);
    };
    m.diagnostics.skippedPoints = count("skippedPoints");
    m.diagnostics.skippedLines = count("skippedLines");
    m.diagnostics.boundaryEdges = count("boundaryEdges");
    m.diagnostics.nonManifoldEdges = count("nonManifoldEdges");
    m.diagnostics.inconsistentEdges = count("inconsistentEdges");
    m.diagnostics.degenerateTriangles = count("degenerateTriangles");
    return result;
}
} // namespace hornsculpt
