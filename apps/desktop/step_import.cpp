#include "step_import.hpp"
#include <QCryptographicHash>
#include <QFileInfo>
#include <QUuid>
#include <algorithm>
#include <hornsculpt/cad_geometry.hpp>
#include <stdexcept>

namespace hornsculpt {
std::shared_ptr<const CadSnapshot> readStepSnapshot(const QString &path, double tolerance,
                                                    std::stop_token stop) {
    const auto cancelled = [&] {
        if (stop.stop_requested())
            throw std::runtime_error("STEP-Import abgebrochen.");
    };
    cancelled();
    const auto geometry = CadGeometry::fromStep(path.toStdString());
    cancelled();
    auto result = std::make_shared<CadSnapshot>();
    result->documentId = QUuid::createUuid().toString(QUuid::WithoutBraces);
    result->documentName = QFileInfo(path).completeBaseName();
    result->sessionId = "step-import";
    result->revision = 1;
    result->brepMetresPerUnit = 1;
    result->geometry = geometry;
    const auto parts = geometry->brepParts();
    for (std::size_t i = 0; i < parts.size(); ++i) {
        const auto objectId = QStringLiteral("step-%1").arg(i + 1);
        const QByteArray bytes = QByteArray::fromStdString(parts[i].data);
        const auto fingerprint = QString::fromLatin1(
            QCryptographicHash::hash(bytes, QCryptographicHash::Sha256).toHex());
        CadObject object;
        object.id = objectId;
        object.name = objectId;
        object.label = result->documentName;
        object.topologyFingerprint = fingerprint;
        object.placement = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
        result->objects.push_back(std::move(object));
        result->brepObjects.push_back({objectId, bytes, parts[i].faceIds});
        for (const auto id : parts[i].faceIds) {
            const auto measures = geometry->faces();
            const auto found = std::find_if(measures.begin(), measures.end(),
                                            [&](const auto &face) { return face.faceId == id; });
            if (found == measures.end())
                throw std::runtime_error("STEP-Flaechenbindung fehlt.");
            CadFace face;
            face.id = id;
            face.objectId = objectId;
            face.sourceKey = fingerprint + QStringLiteral(":%1").arg(id);
            face.label = QStringLiteral("Fläche %1").arg(id);
            face.state = CadFaceState::confirmed;
            face.signature.surfaceType = QString::fromStdString(found->surfaceType);
            face.signature.geometryKey = face.sourceKey;
            face.signature.topologyKey = face.sourceKey;
            face.signature.area = found->area;
            face.signature.center = found->center;
            face.signature.normal = found->normal;
            face.signature.edgeCount = found->edgeCount;
            face.signature.wireCount = found->wireCount;
            result->faces.push_back(std::move(face));
        }
    }
    result->mesh = geometry->previewMesh(tolerance, stop);
    cancelled();
    (void)validateCadSnapshot(*result);
    return result;
}
} // namespace hornsculpt
