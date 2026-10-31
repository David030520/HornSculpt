#pragma once

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <hornsculpt/freecad_link.hpp>
#include <stdexcept>

namespace hornsculpt::testing {
// Exact tetrahedra exported through the real FreeCAD 1.1.4 BRep API, in mm.
inline QJsonObject cadFixture(bool enlarged = false) {
    const QDir directory(QFileInfo(QString::fromUtf8(__FILE__)).absolutePath());
    QFile file(directory.filePath(enlarged ? "fixtures/freecad-tetra-20mm.json"
                                           : "fixtures/freecad-tetra-10mm.json"));
    if (!file.open(QIODevice::ReadOnly))
        throw std::runtime_error("Exact FreeCAD fixture cannot be read.");
    QJsonParseError error;
    const auto document = QJsonDocument::fromJson(file.readAll(), &error);
    if (error.error != QJsonParseError::NoError || !document.isObject())
        throw std::runtime_error("Exact FreeCAD fixture is not valid JSON.");
    return document.object();
}
inline CadSnapshot exactCadFixture() {
    const auto json = cadFixture();
    const auto point = [](QJsonValue value) {
        const auto a = value.toArray();
        return Vec3{a[0].toDouble(), a[1].toDouble(), a[2].toDouble()};
    };
    CadSnapshot result;
    result.documentId = json["documentId"].toString();
    result.documentName = json["documentName"].toString();
    result.sessionId = json["sessionId"].toString();
    result.revision = 1;
    for (const auto &value : json["objects"].toArray()) {
        const auto o = value.toObject();
        CadObject item;
        item.id = o["id"].toString();
        item.name = o["name"].toString();
        item.label = o["label"].toString();
        item.topologyFingerprint = o["topologyFingerprint"].toString();
        const auto placement = o["placement"].toArray();
        for (int i = 0; i < 16; ++i)
            item.placement[i] = placement[i].toDouble();
        result.objects.push_back(std::move(item));
    }
    for (const auto &value : json["faces"].toArray()) {
        const auto f = value.toObject();
        CadFace item;
        item.id = static_cast<std::uint32_t>(result.faces.size() + 1);
        item.objectId = f["objectId"].toString();
        item.sourceKey = f["mappedName"].toString();
        item.label = f["label"].toString();
        item.state = CadFaceState::confirmed;
        item.signature = {f["surfaceType"].toString(),
                          f["geometryKey"].toString(),
                          f["topologyKey"].toString(),
                          f["area"].toDouble(),
                          point(f["center"]),
                          point(f["normal"]),
                          static_cast<std::uint32_t>(f["edgeCount"].toInt()),
                          static_cast<std::uint32_t>(f["wireCount"].toInt())};
        result.faces.push_back(std::move(item));
    }
    for (const auto &value : json["brepParts"].toArray()) {
        const auto p = value.toObject();
        CadBrepObject item;
        item.objectId = p["objectId"].toString();
        item.data = QByteArray::fromBase64(p["data"].toString().toLatin1());
        for (const auto &index : p["faceIndices"].toArray())
            item.faceIds.push_back(static_cast<std::uint32_t>(index.toInt() + 1));
        result.brepObjects.push_back(std::move(item));
    }
    prepareCadGeometry(result, true);
    return result;
}
} // namespace hornsculpt::testing
