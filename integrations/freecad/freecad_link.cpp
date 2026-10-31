#include <hornsculpt/freecad_link.hpp>

#include "cad_snapshot.hpp"
#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QFutureWatcher>
#include <QHash>
#include <QHostAddress>
#include <QJsonDocument>
#include <QJsonObject>
#include <QPointer>
#include <QSaveFile>
#include <QStandardPaths>
#include <QTcpServer>
#include <QTcpSocket>
#include <QTimer>
#include <QUuid>
#include <QtConcurrentRun>
#include <QtEndian>
#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <utility>

namespace hornsculpt {
namespace {
constexpr quint32 maximumFrameBytes = 64U * 1024U * 1024U;
struct DecodeResult {
    std::optional<CadSnapshot> snapshot;
    QString error;
    std::uint32_t nextFaceId{};
};
} // namespace

struct FreeCadLink::Impl {
    explicit Impl(FreeCadLink *owner) : q(owner), server(owner), deadline(owner), watcher(owner) {
        deadline.setSingleShot(true);
        QObject::connect(&server, &QTcpServer::newConnection, q, [this] { accept(); });
        QObject::connect(&deadline, &QTimer::timeout, q, [this] {
            reject(QStringLiteral(
                "FreeCAD: Zeitlimit fuer Anmeldung oder unvollstaendigen Snapshot."));
        });
        QObject::connect(&watcher, &QFutureWatcher<DecodeResult>::finished, q,
                         [this] { decoded(); });
    }

    FreeCadLink *q;
    QTcpServer server;
    QTimer deadline;
    QFutureWatcher<DecodeResult> watcher;
    QPointer<QTcpSocket> client;
    QByteArray incoming;
    std::deque<QByteArray> pending;
    QString nonce;
    QString discoveryFile;
    const QString instanceId = QUuid::createUuid().toString(QUuid::WithoutBraces);
    QString session;
    QString error;
    bool authenticated{};
    bool sourceDirty{};
    std::uint64_t pendingRevision{};
    QString pendingDocument;
    bool decoding{};
    quint64 generation{};
    quint64 decodingGeneration{};
    quint64 mappingGeneration{};
    quint64 decodingMappingGeneration{};
    QByteArray decodingPayload;
    std::uint32_t nextFaceId{1};
    std::shared_ptr<const CadSnapshot> current;
    QHash<QString, std::shared_ptr<const CadSnapshot>> documents;

    void report(const QString &message) {
        error = message;
        emit q->errorOccurred(message);
    }
    void reject(const QString &message) {
        report(message);
        if (client)
            client->abort();
    }
    void disconnected() {
        const bool wasConnected = authenticated;
        authenticated = false;
        incoming.clear();
        pending.clear();
        deadline.stop();
        ++generation;
        if (client)
            client->deleteLater();
        client = nullptr;
        if (wasConnected)
            emit q->connectionChanged(false);
    }
    void accept() {
        while (server.hasPendingConnections()) {
            auto *socket = server.nextPendingConnection();
            if (client || !socket->peerAddress().isLoopback()) {
                socket->abort();
                socket->deleteLater();
                continue;
            }
            client = socket;
            ++generation;
            authenticated = false;
            session.clear();
            incoming.clear();
            pending.clear();
            socket->setReadBufferSize(maximumFrameBytes + 4ULL);
            QObject::connect(socket, &QTcpSocket::readyRead, q, [this] { receive(); });
            QObject::connect(socket, &QTcpSocket::disconnected, q, [this] { disconnected(); });
            deadline.start(5000);
            if (socket->bytesAvailable())
                receive();
        }
    }
    void receive() {
        if (!client)
            return;
        incoming.append(client->readAll());
        if (incoming.size() > static_cast<qsizetype>(maximumFrameBytes) + 4) {
            reject(QStringLiteral("FreeCAD: Empfangspuffer ueberschreitet das Limit."));
            return;
        }
        while (incoming.size() >= 4) {
            const auto size = qFromBigEndian<quint32>(incoming.constData());
            if (size == 0 || size > maximumFrameBytes || (!authenticated && size > 4096)) {
                reject(QStringLiteral("FreeCAD: ungueltige Nachrichtenlaenge."));
                return;
            }
            if (incoming.size() < static_cast<qsizetype>(size) + 4)
                break;
            QByteArray payload = incoming.mid(4, size);
            incoming.remove(0, static_cast<qsizetype>(size) + 4);
            if (!authenticated) {
                QJsonParseError parseError;
                const auto message = QJsonDocument::fromJson(payload, &parseError);
                const auto root = message.object();
                if (parseError.error != QJsonParseError::NoError || !message.isObject() ||
                    root.value(QStringLiteral("type")) != QStringLiteral("hello") ||
                    root.value(QStringLiteral("protocol")).toInt(-1) != 2 ||
                    root.value(QStringLiteral("token")).toString() != nonce ||
                    root.value(QStringLiteral("sessionId")).toString().isEmpty() ||
                    root.value(QStringLiteral("sessionId")).toString().size() > 128) {
                    reject(QStringLiteral("FreeCAD: Anmeldung oder Protokollversion ungueltig."));
                    return;
                }
                session = root.value(QStringLiteral("sessionId")).toString();
                authenticated = true;
                pendingRevision = 0;
                pendingDocument.clear();
                if (sourceDirty) {
                    sourceDirty = false;
                    emit q->sourceDirtyChanged(false);
                }
                error.clear();
                const auto ready =
                    QJsonDocument(QJsonObject{{QStringLiteral("type"), QStringLiteral("ready")},
                                              {QStringLiteral("protocol"), 2}})
                        .toJson(QJsonDocument::Compact);
                client->write(cad_detail::frame(ready));
                emit q->connectionChanged(true);
            } else {
                if (payload.size() < 4096) {
                    const auto control = QJsonDocument::fromJson(payload).object();
                    if (control.value(QStringLiteral("type")) ==
                        QStringLiteral("geometryPending")) {
                        bool revisionValid = false;
                        const auto revision = control.value(QStringLiteral("revision"))
                                                  .toString()
                                                  .toULongLong(&revisionValid);
                        const auto documentId =
                            control.value(QStringLiteral("documentId")).toString();
                        if (control.value(QStringLiteral("protocol")).toInt(-1) != 2 ||
                            control.value(QStringLiteral("sessionId")).toString() != session ||
                            !revisionValid || revision == 0 || documentId.isEmpty()) {
                            reject(QStringLiteral("FreeCAD: ungueltiger Aenderungsstatus."));
                            return;
                        }
                        const auto prior = documents.value(documentId);
                        if ((prior && prior->sessionId == session && revision <= prior->revision) ||
                            (sourceDirty && pendingDocument == documentId &&
                             revision < pendingRevision)) {
                            report(
                                QStringLiteral("FreeCAD: veralteter Aenderungsstatus verworfen."));
                            continue;
                        }
                        pendingRevision = revision;
                        pendingDocument = documentId;
                        if (!sourceDirty) {
                            sourceDirty = true;
                            emit q->sourceDirtyChanged(true);
                        }
                        continue;
                    }
                }
                if (pending.size() >= 2) {
                    reject(QStringLiteral("FreeCAD: zu viele unverarbeitete Geometriesnapshots."));
                    return;
                }
                pending.push_back(std::move(payload));
                processNext();
            }
        }
        if (incoming.isEmpty())
            deadline.stop();
        else if (!deadline.isActive())
            deadline.start(30000);
    }
    void processNext() {
        if (decoding || pending.empty() || !authenticated)
            return;
        auto payload = std::move(pending.front());
        pending.pop_front();
        const auto expectedSession = session;
        const auto previousDocuments = documents;
        const auto firstFaceId = nextFaceId;
        decoding = true;
        decodingGeneration = generation;
        decodingMappingGeneration = mappingGeneration;
        decodingPayload = payload;
        watcher.setFuture(QtConcurrent::run(
            [payload = std::move(payload), expectedSession, previousDocuments, firstFaceId] {
                DecodeResult result;
                try {
                    result.snapshot = cad_detail::decodeSnapshot(payload, expectedSession);
                    result.nextFaceId = firstFaceId;
                    auto &snapshot = *result.snapshot;
                    const auto previous = previousDocuments.value(snapshot.documentId);
                    cad_detail::reconcileSnapshot(snapshot, previous.get(), result.nextFaceId);
                    prepareCadGeometry(snapshot, true);
                } catch (const std::exception &exception) {
                    result.snapshot.reset();
                    result.error = QString::fromUtf8(exception.what());
                }
                return result;
            }));
    }
    void decoded() {
        decoding = false;
        auto result = watcher.result();
        if (decodingGeneration != generation || !authenticated) {
            decodingPayload.clear();
            processNext();
            return;
        }
        if (decodingMappingGeneration != mappingGeneration) {
            // A face edit during native decoding must not be overwritten by its old copy.
            pending.push_front(std::move(decodingPayload));
            processNext();
            return;
        }
        decodingPayload.clear();
        if (!result.snapshot) {
            reject(result.error);
            return;
        }
        auto snapshot = std::make_shared<CadSnapshot>(std::move(*result.snapshot));
        const auto previous = documents.value(snapshot->documentId);
        if (previous && previous->sessionId == snapshot->sessionId &&
            snapshot->revision <= previous->revision) {
            report(QStringLiteral("FreeCAD: veraltete oder doppelte Revision verworfen."));
            processNext();
            return;
        }
        nextFaceId = result.nextFaceId;
        ++mappingGeneration;
        current = std::move(snapshot);
        documents.insert(current->documentId, current);
        if (sourceDirty && current->documentId == pendingDocument &&
            current->revision >= pendingRevision) {
            sourceDirty = false;
            emit q->sourceDirtyChanged(false);
        }
        emit q->snapshotChanged();
        processNext();
    }
    void replaceCurrent(std::shared_ptr<CadSnapshot> replacement) {
        ++mappingGeneration;
        current = std::move(replacement);
        documents.insert(current->documentId, current);
        emit q->snapshotChanged();
    }
};

FreeCadLink::FreeCadLink(QObject *parent) : QObject(parent), impl_(std::make_unique<Impl>(this)) {}
FreeCadLink::~FreeCadLink() {
    stop();
    // Decoding owns immutable data only; wait before destroying its watcher.
    impl_->watcher.waitForFinished();
}

bool FreeCadLink::listen(quint16 requestedPort) {
    stop();
    impl_->nonce = QUuid::createUuid().toString(QUuid::WithoutBraces);
    impl_->error.clear();
    if (!impl_->server.listen(QHostAddress::LocalHost, requestedPort)) {
        impl_->report(impl_->server.errorString());
        return false;
    }
    const auto base = QStandardPaths::writableLocation(QStandardPaths::GenericDataLocation);
    const QDir directory(base + QStringLiteral("/HornSculpt/freecad-link"));
    if (base.isEmpty() || !QDir().mkpath(directory.absolutePath())) {
        impl_->server.close();
        impl_->report(
            QStringLiteral("FreeCAD: Discovery-Verzeichnis konnte nicht erstellt werden."));
        return false;
    }
    const auto path = directory.filePath(QString::number(QCoreApplication::applicationPid()) + '-' +
                                         impl_->instanceId + QStringLiteral(".json"));
    auto name = QCoreApplication::applicationName();
    if (name.isEmpty())
        name = QStringLiteral("HornSculpt");
    const auto data =
        QJsonDocument(QJsonObject{{QStringLiteral("protocol"), 2},
                                  {QStringLiteral("host"), QStringLiteral("127.0.0.1")},
                                  {QStringLiteral("port"), static_cast<int>(port())},
                                  {QStringLiteral("token"), token()},
                                  {QStringLiteral("pid"), QCoreApplication::applicationPid()},
                                  {QStringLiteral("name"),
                                   name + QStringLiteral(" · ") +
                                       QString::number(QCoreApplication::applicationPid())}})
            .toJson(QJsonDocument::Compact);
    QSaveFile file(path);
    if (!file.open(QIODevice::WriteOnly) ||
        !file.setPermissions(QFileDevice::ReadOwner | QFileDevice::WriteOwner) ||
        file.write(data) != data.size() || !file.commit()) {
        impl_->server.close();
        impl_->report(
            QStringLiteral("FreeCAD: Instanz konnte nicht zur Verbindung veröffentlicht werden."));
        return false;
    }
    impl_->discoveryFile = path;
    return true;
}
void FreeCadLink::stop() {
    if (!impl_->discoveryFile.isEmpty()) {
        QFile::remove(impl_->discoveryFile);
        impl_->discoveryFile.clear();
    }
    impl_->server.close();
    if (impl_->client) {
        auto *socket = impl_->client.data();
        QObject::disconnect(socket, nullptr, this, nullptr);
        socket->abort();
        impl_->disconnected();
    } else {
        impl_->pending.clear();
        impl_->incoming.clear();
        impl_->deadline.stop();
        ++impl_->generation;
    }
}
bool FreeCadLink::isListening() const {
    return impl_->server.isListening();
}
bool FreeCadLink::isConnected() const {
    return impl_->authenticated;
}
bool FreeCadLink::isGeometryPending() const {
    return impl_->sourceDirty;
}
quint16 FreeCadLink::port() const {
    return impl_->server.serverPort();
}
QString FreeCadLink::token() const {
    return impl_->nonce;
}
QString FreeCadLink::discoveryFilePath() const {
    return impl_->discoveryFile;
}
QString FreeCadLink::lastError() const {
    return impl_->error;
}
std::shared_ptr<const CadSnapshot> FreeCadLink::snapshot() const {
    return impl_->current;
}

PreparedCadSnapshot::PreparedCadSnapshot(std::shared_ptr<const CadSnapshot> snapshot)
    : snapshot_(std::move(snapshot)) {}

PreparedCadSnapshot PreparedCadSnapshot::prepare(std::shared_ptr<const CadSnapshot> snapshot) {
    if (!snapshot)
        return PreparedCadSnapshot({});
    auto replacement = std::make_shared<CadSnapshot>(*snapshot);
    replacement->geometry = validateCadSnapshot(*replacement);
    return PreparedCadSnapshot(std::move(replacement));
}

std::shared_ptr<const CadSnapshot> PreparedCadSnapshot::snapshot() const noexcept {
    return snapshot_;
}

void FreeCadLink::restoreSnapshot(std::shared_ptr<const CadSnapshot> snapshot) {
    restoreSnapshot(PreparedCadSnapshot::prepare(std::move(snapshot)));
}

void FreeCadLink::restoreSnapshot(const PreparedCadSnapshot &prepared) {
    const auto replacement = prepared.snapshot();
    QHash<QString, std::shared_ptr<const CadSnapshot>> documents;
    std::uint32_t nextFaceId = 1;
    if (replacement) {
        for (const auto &face : replacement->faces)
            nextFaceId = std::max(nextFaceId, face.id + 1);
        for (const auto missing : replacement->unresolvedPreviousFaces)
            nextFaceId = std::max(nextFaceId, missing + 1);
        documents.insert(replacement->documentId, replacement);
    }
    // Prepare replacements first; failure must not discard the active document.
    const bool wasDirty = impl_->sourceDirty;
    stop();
    impl_->current = std::move(replacement);
    ++impl_->mappingGeneration;
    impl_->documents = std::move(documents);
    impl_->nextFaceId = nextFaceId;
    impl_->sourceDirty = false;
    impl_->pendingRevision = 0;
    impl_->pendingDocument.clear();
    impl_->session.clear();
    impl_->nonce.clear();
    impl_->error.clear();
    if (wasDirty)
        emit sourceDirtyChanged(false);
    emit snapshotChanged();
}

bool FreeCadLink::setBoundaries(const std::vector<std::uint32_t> &faceIds,
                                BoundaryCondition boundary, std::complex<double> normalVelocity) {
    if (!impl_->current ||
        (boundary != BoundaryCondition::wall && boundary != BoundaryCondition::driven) ||
        !std::isfinite(normalVelocity.real()) || !std::isfinite(normalVelocity.imag()) ||
        faceIds.empty())
        return false;
    const std::set<std::uint32_t> selected(faceIds.begin(), faceIds.end());
    if (selected.size() != faceIds.size())
        return false;
    auto replacement = std::make_shared<CadSnapshot>(*impl_->current);
    std::size_t found = 0;
    for (auto &face : replacement->faces) {
        if (!selected.contains(face.id))
            continue;
        face.boundary = boundary;
        face.normalVelocity = normalVelocity;
        face.state = CadFaceState::confirmed;
        ++found;
    }
    if (found != selected.size())
        return false;
    for (auto &triangle : replacement->mesh.triangles)
        if (selected.contains(triangle.faceId))
            triangle.boundary = boundary;
    impl_->replaceCurrent(std::move(replacement));
    return true;
}
bool FreeCadLink::confirmFaces(const std::vector<std::uint32_t> &faceIds) {
    if (!impl_->current || faceIds.empty())
        return false;
    const std::set<std::uint32_t> selected(faceIds.begin(), faceIds.end());
    if (selected.size() != faceIds.size())
        return false;
    auto replacement = std::make_shared<CadSnapshot>(*impl_->current);
    std::size_t found = 0;
    for (auto &face : replacement->faces) {
        if (selected.contains(face.id)) {
            face.state = CadFaceState::confirmed;
            ++found;
        }
    }
    if (found != selected.size())
        return false;
    impl_->replaceCurrent(std::move(replacement));
    return true;
}
bool FreeCadLink::setAssignments(const QString &documentId, const QString &sessionId,
                                 std::uint64_t revision,
                                 const std::vector<CadFaceAssignment> &faces,
                                 const std::vector<std::uint32_t> &missing) {
    if (!impl_->current || impl_->sourceDirty || impl_->current->documentId != documentId ||
        impl_->current->sessionId != sessionId || impl_->current->revision != revision ||
        faces.size() != impl_->current->faces.size())
        return false;
    std::map<std::uint32_t, CadFaceAssignment> edits;
    for (const auto &face : faces)
        if (!face.id ||
            (face.boundary != BoundaryCondition::wall &&
             face.boundary != BoundaryCondition::driven) ||
            !std::isfinite(face.normalVelocity.real()) ||
            !std::isfinite(face.normalVelocity.imag()) ||
            (face.state != CadFaceState::confirmed && face.state != CadFaceState::newFace &&
             face.state != CadFaceState::needsReview) ||
            !edits.emplace(face.id, face).second)
            return false;
    for (const auto &face : impl_->current->faces)
        if (!edits.contains(face.id))
            return false;
    std::set<std::uint32_t> uniqueMissing;
    for (const auto id : missing)
        if (!id || id == std::numeric_limits<std::uint32_t>::max() || edits.contains(id) ||
            !uniqueMissing.insert(id).second)
            return false;
    auto replacement = std::make_shared<CadSnapshot>(*impl_->current);
    for (auto &face : replacement->faces) {
        const auto &edit = edits.at(face.id);
        face.boundary = edit.boundary;
        face.normalVelocity = edit.normalVelocity;
        face.state = edit.state;
    }
    for (auto &triangle : replacement->mesh.triangles)
        triangle.boundary = edits.at(triangle.faceId).boundary;
    replacement->unresolvedPreviousFaces = missing;
    for (const auto id : missing)
        impl_->nextFaceId = std::max(impl_->nextFaceId, id + 1);
    impl_->replaceCurrent(std::move(replacement));
    return true;
}

void FreeCadLink::acknowledgeMissingFaces() {
    if (!impl_->current || impl_->current->unresolvedPreviousFaces.empty())
        return;
    auto replacement = std::make_shared<CadSnapshot>(*impl_->current);
    replacement->unresolvedPreviousFaces.clear();
    impl_->replaceCurrent(std::move(replacement));
}

} // namespace hornsculpt
