#include "project_json.hpp"
#include <QCryptographicHash>

namespace hornsculpt::project_json {
QJsonObject encodeCad(const CadSnapshot &s) {
    QJsonArray faces, objects, missing, breps;
    for (const auto &part : s.brepObjects) {
        QJsonArray ids;
        for (const auto id : part.faceIds)
            ids.append(static_cast<double>(id));
        breps.append(QJsonObject{
            {"objectId", part.objectId},
            {"data", QString::fromLatin1(part.data.toBase64())},
            {"sha256",
             QString::fromLatin1(
                 QCryptographicHash::hash(part.data, QCryptographicHash::Sha256).toHex())},
            {"faceIds", ids}});
    }
    for (const auto &o : s.objects) {
        QJsonArray placement;
        for (auto value : o.placement)
            placement.append(value);
        objects.append(QJsonObject{{"id", o.id},
                                   {"name", o.name},
                                   {"label", o.label},
                                   {"topologyFingerprint", o.topologyFingerprint},
                                   {"placement", placement}});
    }
    for (const auto &f : s.faces) {
        const auto &g = f.signature;
        faces.append(QJsonObject{{"id", double(f.id)},
                                 {"objectId", f.objectId},
                                 {"sourceKey", f.sourceKey},
                                 {"label", f.label},
                                 {"boundary", int(f.boundary)},
                                 {"normalVelocity", encodeComplex(f.normalVelocity)},
                                 {"state", int(f.state)},
                                 {"signature", QJsonObject{{"surfaceType", g.surfaceType},
                                                           {"geometryKey", g.geometryKey},
                                                           {"topologyKey", g.topologyKey},
                                                           {"area", g.area},
                                                           {"center", encodeVector(g.center)},
                                                           {"normal", encodeVector(g.normal)},
                                                           {"edgeCount", double(g.edgeCount)},
                                                           {"wireCount", double(g.wireCount)}}}});
    }
    for (auto id : s.unresolvedPreviousFaces)
        missing.append(double(id));
    return {{"documentId", s.documentId},
            {"documentName", s.documentName},
            {"sessionId", s.sessionId},
            {"revision", encodeWide(s.revision)},
            {"mesh", encodeMesh(s.mesh)},
            {"brepObjects", breps},
            {"brepMetresPerUnit", s.brepMetresPerUnit},
            {"faces", faces},
            {"objects", objects},
            {"acousticOrigin", encodeVector(s.acousticOrigin)},
            {"acousticAxis", encodeVector(s.acousticAxis)},
            {"unresolvedPreviousFaces", missing}};
}

CadSnapshot decodeCad(const QJsonValue &snapshotValue, bool exactGeometryFields) {
    const auto o = object(snapshotValue);
    auto common = o;
    if (exactGeometryFields) {
        if (!common.contains("brepObjects") || !common.contains("brepMetresPerUnit"))
            invalid(QStringLiteral("Exakte CAD-Quellfelder fehlen."));
        common.remove("brepObjects");
        common.remove("brepMetresPerUnit");
    }
    keys(common, {"documentId", "documentName", "sessionId", "revision", "mesh", "faces", "objects",
                  "acousticOrigin", "acousticAxis", "unresolvedPreviousFaces"});
    CadSnapshot s;
    s.documentId = text(o["documentId"]);
    s.documentName = text(o["documentName"], true);
    s.sessionId = text(o["sessionId"]);
    s.revision = wideInteger(o["revision"]);
    s.mesh = decodeMesh(o["mesh"], false);
    s.acousticOrigin = vector(o["acousticOrigin"]);
    s.acousticAxis = vector(o["acousticAxis"]);
    for (const auto &value : array(o["objects"], 10000, 1)) {
        const auto entry = object(value);
        keys(entry, {"id", "name", "label", "topologyFingerprint", "placement"});
        CadObject item;
        item.id = text(entry["id"]);
        item.name = text(entry["name"]);
        item.label = text(entry["label"], true);
        item.topologyFingerprint = text(entry["topologyFingerprint"]);
        const auto placement = array(entry["placement"], 16, 16);
        for (int i = 0; i < 16; ++i)
            item.placement[i] = number(placement[i]);
        s.objects.push_back(std::move(item));
    }
    for (const auto &value : array(o["faces"], 100000, 1)) {
        const auto entry = object(value);
        keys(entry, {"id", "objectId", "sourceKey", "label", "boundary", "normalVelocity", "state",
                     "signature"});
        CadFace item;
        item.id = integer(entry["id"]);
        item.objectId = text(entry["objectId"]);
        item.sourceKey = text(entry["sourceKey"], true);
        item.label = text(entry["label"], true);
        item.boundary = static_cast<BoundaryCondition>(integer(entry["boundary"], 1));
        item.normalVelocity = complex(entry["normalVelocity"]);
        item.state = static_cast<CadFaceState>(integer(entry["state"], 2));
        const auto g = object(entry["signature"]);
        keys(g, {"surfaceType", "geometryKey", "topologyKey", "area", "center", "normal",
                 "edgeCount", "wireCount"});
        item.signature = {text(g["surfaceType"]),
                          text(g["geometryKey"]),
                          text(g["topologyKey"]),
                          positive(g["area"]),
                          vector(g["center"]),
                          vector(g["normal"]),
                          integer(g["edgeCount"], 100000),
                          integer(g["wireCount"], 100000)};
        s.faces.push_back(std::move(item));
    }
    for (const auto &id : array(o["unresolvedPreviousFaces"], 100000))
        s.unresolvedPreviousFaces.push_back(integer(id));
    if (exactGeometryFields) {
        s.brepMetresPerUnit = positive(o["brepMetresPerUnit"]);
        for (const auto &value : array(o["brepObjects"], 10000)) {
            const auto part = object(value);
            keys(part, {"objectId", "data", "sha256", "faceIds"});
            CadBrepObject item;
            item.objectId = text(part["objectId"]);
            if (!part["data"].isString() || part["data"].toString().size() > 64 * 1024 * 1024)
                invalid(QStringLiteral("Ungültige BRep-Datengröße."));
            const auto encoded = part["data"].toString().toLatin1();
            const auto decoded =
                QByteArray::fromBase64Encoding(encoded, QByteArray::AbortOnBase64DecodingErrors);
            if (!decoded || decoded.decoded.isEmpty() || decoded.decoded.toBase64() != encoded)
                invalid(QStringLiteral("Ungültige Base64-BRep-Daten."));
            item.data = decoded.decoded;
            if (text(part["sha256"]).toLatin1() !=
                QCryptographicHash::hash(item.data, QCryptographicHash::Sha256).toHex())
                invalid(QStringLiteral("BRep-Prüfsumme stimmt nicht."));
            for (const auto &id : array(part["faceIds"], 100000, 1))
                item.faceIds.push_back(integer(id));
            s.brepObjects.push_back(std::move(item));
        }
    }
    s.geometry = validateCadSnapshot(s);
    return s;
}
} // namespace hornsculpt::project_json
