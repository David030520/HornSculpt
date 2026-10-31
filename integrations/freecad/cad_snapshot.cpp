#include "cad_snapshot.hpp"

#include <QCryptographicHash>
#include <QHash>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSet>
#include <QtEndian>
#include <algorithm>
#include <cmath>
#include <hornsculpt/meshing.hpp>
#include <limits>
#include <stdexcept>

namespace hornsculpt {

bool CadSnapshot::hasUnresolvedAssignments() const noexcept {
    return !unresolvedPreviousFaces.empty() ||
           std::any_of(faces.begin(), faces.end(),
                       [](const CadFace &face) { return face.state == CadFaceState::needsReview; });
}

namespace cad_detail {
namespace {
[[noreturn]] void fail(const char *message) {
    throw std::invalid_argument(message);
}
QString text(const QJsonObject &object, const char *name) {
    const auto value = object.value(QLatin1String(name));
    if (!value.isString() || value.toString().isEmpty() || value.toString().size() > 4096)
        fail("FreeCAD: Textfeld fehlt oder ist ungueltig.");
    return value.toString();
}
double number(const QJsonValue &value) {
    if (!value.isDouble() || !std::isfinite(value.toDouble()))
        fail("FreeCAD: endliche Zahl erwartet.");
    return value.toDouble();
}
std::uint32_t integer(const QJsonValue &value, std::uint32_t maximum) {
    const double result = number(value);
    if (result < 0.0 || result > maximum || result != std::floor(result))
        fail("FreeCAD: Index ausserhalb des erlaubten Bereichs.");
    return static_cast<std::uint32_t>(result);
}
Vec3 vector(const QJsonValue &value) {
    const auto array = value.toArray();
    if (!value.isArray() || array.size() != 3)
        fail("FreeCAD: Vektor mit drei Komponenten erwartet.");
    return {number(array[0]), number(array[1]), number(array[2])};
}
double norm(Vec3 value) {
    return std::hypot(value.x, value.y, value.z);
}
QString groupingKey(const CadFace &face, const QString &value) {
    return face.objectId + QChar(u'\x1f') + value;
}
bool mappedName(const QString &key) {
    if (key.isEmpty())
        return false;
    if (!key.startsWith(QStringLiteral("Face")))
        return true;
    bool ordinal = false;
    (void)key.mid(4).toUInt(&ordinal);
    return !ordinal;
}
} // namespace

void validatePlacement(const std::array<double, 16> &placement) {
    constexpr double rigidTolerance = 1e-6;
    for (int row = 0; row < 3; ++row) {
        for (int other = 0; other < 3; ++other) {
            double dot = 0.0;
            for (int column = 0; column < 3; ++column)
                dot += placement[4 * row + column] * placement[4 * other + column];
            if (std::abs(dot - (row == other ? 1.0 : 0.0)) > rigidTolerance)
                fail("FreeCAD: Platzierung muss eine starre Transformation sein.");
        }
    }
    if (std::abs(placement[12]) > rigidTolerance || std::abs(placement[13]) > rigidTolerance ||
        std::abs(placement[14]) > rigidTolerance || std::abs(placement[15] - 1.0) > rigidTolerance)
        fail("FreeCAD: ungueltige homogene Platzierung.");
    const double determinant =
        placement[0] * (placement[5] * placement[10] - placement[6] * placement[9]) -
        placement[1] * (placement[4] * placement[10] - placement[6] * placement[8]) +
        placement[2] * (placement[4] * placement[9] - placement[5] * placement[8]);
    if (std::abs(determinant - 1.0) > rigidTolerance)
        fail("FreeCAD: Platzierung muss eine eigentliche Rotation enthalten.");
}
QByteArray frame(const QByteArray &payload) {
    QByteArray result(4, Qt::Uninitialized);
    qToBigEndian(static_cast<quint32>(payload.size()), result.data());
    result.append(payload);
    return result;
}

CadSnapshot decodeSnapshot(const QByteArray &payload, const QString &sessionId) {
    QJsonParseError parseError;
    const auto document = QJsonDocument::fromJson(payload, &parseError);
    if (parseError.error != QJsonParseError::NoError || !document.isObject())
        fail("FreeCAD: ungueltiger JSON-Snapshot.");
    const auto root = document.object();
    if (root.value(QStringLiteral("protocol")).toInt(-1) != 2 ||
        text(root, "type") != QStringLiteral("snapshot") || text(root, "sessionId") != sessionId ||
        text(root, "units") != QStringLiteral("m"))
        fail("FreeCAD: inkompatible Version, Sitzung oder Einheit.");
    CadSnapshot result;
    result.sessionId = sessionId;
    result.documentId = text(root, "documentId");
    result.documentName = text(root, "documentName");
    bool revisionValid = false;
    result.revision = text(root, "revision").toULongLong(&revisionValid);
    if (!revisionValid || result.revision == 0)
        fail("FreeCAD: ungueltige Revision.");
    result.acousticOrigin = vector(root.value(QStringLiteral("acousticOrigin")));
    result.acousticAxis = vector(root.value(QStringLiteral("acousticAxis")));
    if (std::abs(norm(result.acousticAxis) - 1.0) > 1e-6)
        fail("FreeCAD: akustische Achse muss normiert sein.");

    const auto objects = root.value(QStringLiteral("objects")).toArray();
    if (objects.isEmpty() || objects.size() > 10000)
        fail("FreeCAD: Objektanzahl ungueltig.");
    QSet<QString> objectIds;
    for (const auto &value : objects) {
        if (!value.isObject())
            fail("FreeCAD: Objekt erwartet.");
        const auto object = value.toObject();
        CadObject item;
        item.id = text(object, "id");
        item.name = text(object, "name");
        item.label = text(object, "label");
        item.topologyFingerprint = text(object, "topologyFingerprint");
        if (objectIds.contains(item.id))
            fail("FreeCAD: doppelte Objektidentitaet.");
        objectIds.insert(item.id);
        const auto placement = object.value(QStringLiteral("placement")).toArray();
        if (placement.size() != 16)
            fail("FreeCAD: Platzierung mit 16 Werten erwartet.");
        for (qsizetype i = 0; i < 16; ++i)
            item.placement[static_cast<std::size_t>(i)] = number(placement[i]);
        validatePlacement(item.placement);
        result.objects.push_back(std::move(item));
    }
    const auto faces = root.value(QStringLiteral("faces")).toArray();
    if (faces.isEmpty() || faces.size() > 100000)
        fail("FreeCAD: Flaechenanzahl ungueltig.");
    for (const auto &value : faces) {
        if (!value.isObject())
            fail("FreeCAD: Flaechenobjekt erwartet.");
        const auto object = value.toObject();
        CadFace face;
        face.objectId = text(object, "objectId");
        if (!objectIds.contains(face.objectId))
            fail("FreeCAD: Flaeche ohne bekanntes Objekt.");
        const auto sourceKey = object.value(QStringLiteral("mappedName"));
        if (!sourceKey.isString() || sourceKey.toString().size() > 4096)
            fail("FreeCAD: ungueltiger Topologiename.");
        face.sourceKey = sourceKey.toString();
        face.label = text(object, "label");
        face.signature.surfaceType = text(object, "surfaceType");
        face.signature.geometryKey = text(object, "geometryKey");
        face.signature.topologyKey = text(object, "topologyKey");
        face.signature.area = number(object.value(QStringLiteral("area")));
        if (face.signature.area <= 0.0)
            fail("FreeCAD: Flaecheninhalt muss positiv sein.");
        face.signature.center = vector(object.value(QStringLiteral("center")));
        face.signature.normal = vector(object.value(QStringLiteral("normal")));
        if (norm(face.signature.normal) > 1.0 + 1e-6)
            fail("FreeCAD: Flaechennormale ungueltig.");
        face.signature.edgeCount = integer(object.value(QStringLiteral("edgeCount")), 100000);
        face.signature.wireCount = integer(object.value(QStringLiteral("wireCount")), 100000);
        result.faces.push_back(std::move(face));
    }
    if (root.contains("vertices") || root.contains("triangles"))
        fail("FreeCAD: Protokoll 2 erwartet exakte BRep-Geometrie statt eines Quellnetzes.");
    result.brepMetresPerUnit = number(root.value("brepMetresPerUnit"));
    if (result.brepMetresPerUnit != .001)
        fail("FreeCAD: BRep-Quelle muss explizit Millimeter verwenden.");
    const auto parts = root.value("brepParts").toArray();
    if (!root.value("brepParts").isArray() || parts.size() != objects.size())
        fail("FreeCAD: Jedes Objekt benoetigt genau ein BRep-Artefakt.");
    std::vector<bool> usedFaces(result.faces.size(), false);
    QSet<QString> partObjects;
    for (const auto &value : parts) {
        if (!value.isObject())
            fail("FreeCAD: BRep-Objekt erwartet.");
        const auto part = value.toObject();
        CadBrepObject item;
        item.objectId = text(part, "objectId");
        if (!objectIds.contains(item.objectId) || partObjects.contains(item.objectId))
            fail("FreeCAD: BRep-Objektbindung ist nicht eindeutig.");
        partObjects.insert(item.objectId);
        if (!part.value("data").isString())
            fail("FreeCAD: Base64-BRep-Daten erwartet.");
        const auto encoded = part.value("data").toString().toLatin1();
        const auto decoded =
            QByteArray::fromBase64Encoding(encoded, QByteArray::AbortOnBase64DecodingErrors);
        if (!decoded || decoded.decoded.isEmpty() || decoded.decoded.toBase64() != encoded)
            fail("FreeCAD: Ungueltige Base64-BRep-Daten.");
        item.data = decoded.decoded;
        if (text(part, "sha256").toLatin1() !=
            QCryptographicHash::hash(item.data, QCryptographicHash::Sha256).toHex())
            fail("FreeCAD: BRep-Pruefsumme stimmt nicht.");
        const auto indices = part.value("faceIndices").toArray();
        if (!part.value("faceIndices").isArray() || indices.isEmpty())
            fail("FreeCAD: BRep-Flaechenbindung fehlt.");
        for (const auto &index : indices) {
            const auto face = integer(index, static_cast<std::uint32_t>(faces.size() - 1));
            if (usedFaces[face] || result.faces[face].objectId != item.objectId)
                fail("FreeCAD: BRep-Flaechenbindung ist doppelt oder gehoert zum falschen Objekt.");
            usedFaces[face] = true;
            item.faceIds.push_back(face + 1);
        }
        result.brepObjects.push_back(std::move(item));
    }
    if (std::find(usedFaces.begin(), usedFaces.end(), false) != usedFaces.end())
        fail("FreeCAD: Nicht alle Flaechen sind an das exakte BRep gebunden.");
    return result;
}

void reconcileSnapshot(CadSnapshot &current, const CadSnapshot *previous,
                       std::uint32_t &nextFaceId) {
    std::vector<bool> claimed(previous ? previous->faces.size() : 0, false);
    std::vector<int> matches(current.faces.size(), -1);
    QHash<QString, QString> oldTopology, newTopology;
    if (previous) {
        for (const auto &object : previous->objects)
            oldTopology.insert(object.id, object.topologyFingerprint);
        for (const auto &object : current.objects)
            newTopology.insert(object.id, object.topologyFingerprint);
        const auto matchUnique = [&](auto keyFunction) {
            QHash<QString, std::vector<std::size_t>> oldGroups, newGroups;
            for (std::size_t i = 0; i < previous->faces.size(); ++i) {
                if (claimed[i])
                    continue;
                const auto key = keyFunction(previous->faces[i]);
                if (!key.isEmpty())
                    oldGroups[key].push_back(i);
            }
            for (std::size_t i = 0; i < current.faces.size(); ++i) {
                if (matches[i] >= 0)
                    continue;
                const auto key = keyFunction(current.faces[i]);
                if (!key.isEmpty())
                    newGroups[key].push_back(i);
            }
            for (auto group = newGroups.constBegin(); group != newGroups.constEnd(); ++group) {
                const auto candidates = oldGroups.value(group.key());
                if (group.value().size() != 1 || candidates.size() != 1)
                    continue;
                const auto oldIndex = candidates.front(), newIndex = group.value().front();
                claimed[oldIndex] = true;
                matches[newIndex] = static_cast<int>(oldIndex);
            }
        };
        matchUnique([](const CadFace &face) {
            return mappedName(face.sourceKey)
                       ? groupingKey(face, QStringLiteral("name:") + face.sourceKey)
                       : QString();
        });
        matchUnique([](const CadFace &face) {
            return groupingKey(face, QStringLiteral("geometry:") + face.signature.geometryKey);
        });
        matchUnique([&](const CadFace &face) {
            const auto old = oldTopology.constFind(face.objectId),
                       now = newTopology.constFind(face.objectId);
            if (old == oldTopology.constEnd() || now == newTopology.constEnd() ||
                old.value() != now.value())
                return QString();
            return groupingKey(face, QStringLiteral("topology:") + face.signature.topologyKey);
        });
    }
    std::vector<std::uint32_t> faceIds;
    faceIds.reserve(current.faces.size());
    for (std::size_t i = 0; i < current.faces.size(); ++i) {
        auto &face = current.faces[i];
        if (matches[i] >= 0) {
            const auto &old = previous->faces[static_cast<std::size_t>(matches[i])];
            face.id = old.id;
            face.boundary = old.boundary;
            face.normalVelocity = old.normalVelocity;
            face.state = old.state;
        } else {
            if (nextFaceId == std::numeric_limits<std::uint32_t>::max())
                fail("FreeCAD: Flaechenidentitaeten aufgebraucht.");
            face.id = nextFaceId++;
            face.state = previous && oldTopology.contains(face.objectId) ? CadFaceState::needsReview
                                                                         : CadFaceState::newFace;
        }
        faceIds.push_back(face.id);
    }
    if (previous) {
        current.unresolvedPreviousFaces = previous->unresolvedPreviousFaces;
        for (std::size_t i = 0; i < claimed.size(); ++i)
            if (!claimed[i])
                current.unresolvedPreviousFaces.push_back(previous->faces[i].id);
    }
    for (auto &part : current.brepObjects)
        for (auto &id : part.faceIds)
            id = faceIds.at(id - 1);
}

} // namespace cad_detail
} // namespace hornsculpt
