#include "cad_fixture.hpp"
#include "project_file.hpp"
#include <hornsculpt/cad_geometry.hpp>
#include <hornsculpt/freecad_link.hpp>
#include <hornsculpt/meshing.hpp>

#include <QCoreApplication>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QFile>
#include <QHostAddress>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QProcess>
#include <QProcessEnvironment>
#include <QTcpSocket>
#include <QTemporaryDir>
#include <QThread>
#include <QtEndian>
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

using namespace hornsculpt;

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
void waitUntil(const std::function<bool()> &condition, int milliseconds = 5000) {
    QElapsedTimer timer;
    timer.start();
    while (!condition() && timer.elapsed() < milliseconds) {
        QCoreApplication::processEvents(QEventLoop::AllEvents, 10);
        QThread::msleep(1);
    }
    require(condition(), "FreeCAD condition timed out");
}
void settle(int milliseconds = 50) {
    QElapsedTimer timer;
    timer.start();
    while (timer.elapsed() < milliseconds) {
        QCoreApplication::processEvents(QEventLoop::AllEvents, 10);
        QThread::msleep(1);
    }
}
QByteArray frame(QByteArray payload) {
    QByteArray bytes(4, Qt::Uninitialized);
    qToBigEndian(static_cast<quint32>(payload.size()), bytes.data());
    bytes += payload;
    return bytes;
}
QByteArray encode(const QJsonObject &object) {
    return frame(QJsonDocument(object).toJson(QJsonDocument::Compact));
}
void connectClient(FreeCadLink &link, QTcpSocket &client, const QString &session,
                   const QString &token = {}) {
    client.connectToHost(QHostAddress::LocalHost, link.port());
    waitUntil([&] { return client.state() == QAbstractSocket::ConnectedState; });
    const auto hello = encode({{QStringLiteral("type"), QStringLiteral("hello")},
                               {QStringLiteral("protocol"), 2},
                               {QStringLiteral("sessionId"), session},
                               {QStringLiteral("token"), token.isEmpty() ? link.token() : token}});
    client.write(hello.left(2));
    settle(10);
    require(!link.isConnected(), "partial authentication is not committed");
    client.write(hello.mid(2));
    if (token.isEmpty())
        waitUntil([&] { return link.isConnected(); });
}

QJsonObject snapshot(const QString &session, unsigned revision, bool changedGeometry = false,
                     bool changedTopology = false, bool reordered = false,
                     const QString &document = QStringLiteral("doc-a")) {
    auto result = testing::cadFixture(changedGeometry);
    auto object = result["objects"].toArray()[0].toObject();
    object.insert("topologyFingerprint", changedTopology ? "changed" : "tetra");
    QJsonArray objects;
    objects.append(object);
    const auto originalFaces = result["faces"].toArray();
    const std::array<unsigned, 4> order =
        reordered ? std::array<unsigned, 4>{2, 0, 3, 1} : std::array<unsigned, 4>{0, 1, 2, 3};
    std::array<unsigned, 4> inverse{};
    QJsonArray faces;
    for (unsigned i = 0; i < 4; ++i) {
        inverse[order[i]] = i;
        auto face = originalFaces[order[i]].toObject();
        face.insert("mappedName", QStringLiteral("Face%1").arg(i + 1));
        face.insert("geometryKey",
                    QStringLiteral("g%1-%2").arg(order[i]).arg(changedGeometry ? revision : 0));
        face.insert("topologyKey", QStringLiteral("t%1").arg(order[i]));
        faces.append(face);
    }
    auto part = result["brepParts"].toArray()[0].toObject();
    QJsonArray indices;
    for (const auto index : inverse)
        indices.append(static_cast<int>(index));
    part.insert("faceIndices", indices);
    QJsonArray parts;
    parts.append(part);
    result.insert("sessionId", session);
    result.insert("documentId", document);
    result.insert("revision", QString::number(revision));
    result.insert("objects", objects);
    result.insert("faces", faces);
    result.insert("brepParts", parts);
    return result;
}

const CadFace &faceByTopology(const CadSnapshot &snapshot, const QString &key) {
    const auto found =
        std::find_if(snapshot.faces.begin(), snapshot.faces.end(),
                     [&](const CadFace &face) { return face.signature.topologyKey == key; });
    require(found != snapshot.faces.end(), "face topology key exists");
    return *found;
}

void protocolAndIdentity() {
    FreeCadLink link;
    require(link.listen(), "loopback server starts");
    require(!link.token().isEmpty() && link.port() != 0, "ephemeral authenticated endpoint");
    const auto discoveryPath = link.discoveryFilePath();
    QFile discovery(discoveryPath);
    require(discovery.open(QIODevice::ReadOnly), "live link publishes discovery descriptor");
    const auto descriptor = QJsonDocument::fromJson(discovery.readAll()).object();
    discovery.close();
    require(descriptor.size() == 6 && descriptor["protocol"].toInt() == 2 &&
                descriptor["host"].toString() == "127.0.0.1" &&
                descriptor["port"].toInt() == link.port() &&
                descriptor["token"].toString() == link.token() &&
                descriptor["pid"].toInteger() == QCoreApplication::applicationPid() &&
                !descriptor["name"].toString().isEmpty(),
            "discovery descriptor identifies the authenticated endpoint");
    QTcpSocket client;
    connectClient(link, client, QStringLiteral("session-a"));
    const auto first = encode(snapshot(QStringLiteral("session-a"), 1));
    client.write(first.left(first.size() / 2));
    settle();
    require(!link.snapshot(), "partial snapshot does not commit");
    client.write(first.mid(first.size() / 2));
    waitUntil([&] { return bool(link.snapshot()); });
    require(link.snapshot()->geometry && link.snapshot()->brepObjects.size() == 1,
            "Exact BRep source is retained independently of its derived preview");
    double extent = 0;
    for (const auto &point : link.snapshot()->mesh.vertices)
        extent = std::max(extent, point.x);
    require(std::abs(extent - .01) < 1e-12, "Native BRep preview uses metres");
    const auto drivenId = faceByTopology(*link.snapshot(), QStringLiteral("t1")).id;
    require(link.setBoundaries({drivenId}, BoundaryCondition::driven, {0.2, -0.3}),
            "assign driven face");
    std::vector<CadFaceAssignment> liveAssignments;
    for (const auto &face : link.snapshot()->faces)
        liveAssignments.push_back({face.id, face.boundary, face.normalVelocity, face.state});
    const auto liveToken = link.token();
    require(link.setAssignments(link.snapshot()->documentId, link.snapshot()->sessionId,
                                link.snapshot()->revision, liveAssignments,
                                link.snapshot()->unresolvedPreviousFaces) &&
                link.isConnected() && link.isListening() && link.token() == liveToken,
            "Live assignment replay disconnected FreeCAD or changed authentication");
    const auto retainedSnapshot = link.snapshot();

    client.write(encode({{QStringLiteral("type"), QStringLiteral("geometryPending")},
                         {QStringLiteral("protocol"), 2},
                         {QStringLiteral("sessionId"), QStringLiteral("session-a")},
                         {QStringLiteral("documentId"), QStringLiteral("doc-a")},
                         {QStringLiteral("revision"), QStringLiteral("2")}}));
    waitUntil([&] { return link.isGeometryPending(); });
    client.write(encode(snapshot(QStringLiteral("session-a"), 2, true, false, true)));
    waitUntil([&] { return link.snapshot()->revision == 2; });
    require(!link.isGeometryPending(), "pending geometry resolves only on committed revision");
    const auto &preserved = faceByTopology(*link.snapshot(), QStringLiteral("t1"));
    require(preserved.id == drivenId && preserved.boundary == BoundaryCondition::driven &&
                preserved.normalVelocity == std::complex<double>(0.2, -0.3),
            "unique semantic face survives geometry and raw ordinal changes");
    require(retainedSnapshot->revision == 1, "published snapshots remain immutable");
    require(!link.snapshot()->hasUnresolvedAssignments(), "unique edit needs no remapping");

    client.write(encode(snapshot(QStringLiteral("session-a"), 1)));
    settle();
    require(link.isConnected() && link.snapshot()->revision == 2,
            "out-of-order revision discarded");
    client.write(encode(snapshot(QStringLiteral("session-a"), 3, true, true, true)));
    waitUntil([&] { return link.snapshot()->revision == 3; });
    require(link.snapshot()->hasUnresolvedAssignments(), "changed topology requires review");
    require(std::find(link.snapshot()->unresolvedPreviousFaces.begin(),
                      link.snapshot()->unresolvedPreviousFaces.end(),
                      drivenId) != link.snapshot()->unresolvedPreviousFaces.end(),
            "lost driven assignment is explicit");
    for (const auto &face : link.snapshot()->faces) {
        require(face.id != drivenId && face.boundary == BoundaryCondition::wall,
                "FaceN does not silently transfer driven assignment");
    }
    const auto review = link.snapshot();
    for (const auto &face : review->faces)
        require(link.confirmFaces({face.id}), "review new face");
    link.acknowledgeMissingFaces();
    require(!link.snapshot()->hasUnresolvedAssignments(), "explicit review resolves mapping");

    const auto beforeDisconnect = link.snapshot();
    const auto incomplete = encode(snapshot(QStringLiteral("session-a"), 4));
    client.write(incomplete.left(incomplete.size() / 3));
    settle(10);
    client.abort();
    waitUntil([&] { return !link.isConnected(); });
    require(link.snapshot() == beforeDisconnect, "disconnect retains complete snapshot");

    QTcpSocket reconnect;
    connectClient(link, reconnect, QStringLiteral("session-b"));
    reconnect.write(encode(snapshot(QStringLiteral("session-b"), 1, true, true, true)));
    waitUntil([&] { return link.snapshot()->sessionId == QStringLiteral("session-b"); });
    require(link.snapshot()->revision == 1, "new authenticated session restarts revision");
    const auto immutable = link.snapshot();
    auto invalid = snapshot(QStringLiteral("session-b"), 2);
    auto parts = invalid["brepParts"].toArray();
    auto part = parts[0].toObject();
    part.insert("faceIndices", QJsonArray{0, 1, 2});
    parts[0] = part;
    invalid.insert("brepParts", parts);
    reconnect.write(encode(invalid));
    waitUntil([&] { return !link.isConnected(); });
    require(link.snapshot() == immutable,
            "Incomplete exact face binding cannot replace a valid snapshot");

    QTcpSocket unauthorized;
    connectClient(link, unauthorized, QStringLiteral("bad"), QStringLiteral("wrong-token"));
    waitUntil([&] { return unauthorized.state() == QAbstractSocket::UnconnectedState; });
    require(link.snapshot() == immutable, "authentication failure cannot change geometry");
    QTcpSocket oversized;
    connectClient(link, oversized, QStringLiteral("session-c"));
    QByteArray size(4, Qt::Uninitialized);
    qToBigEndian(quint32{0xffffffffU}, size.data());
    oversized.write(size);
    waitUntil([&] { return !link.isConnected(); });
    require(link.snapshot() == immutable, "oversized frame rejected without publication");
    link.stop();
    require(!QFile::exists(discoveryPath), "stopping removes the discovery descriptor");
}

void realBridge(const QString &python, const QString &freecadLib, const QString &bridgeScript) {
    FreeCadLink link;
    require(link.listen(), "real bridge loopback server");
    QTemporaryDir directory;
    require(directory.isValid(), "bridge temporary directory");
    const QString scriptPath = directory.filePath(QStringLiteral("cad_live_smoke.py"));
    QFile script(scriptPath);
    require(script.open(QIODevice::WriteOnly), "write FreeCAD integration runner");
    script.write(R"PY(import sys,os,importlib.util
sys.path.insert(0,sys.argv[1])
os.add_dll_directory(os.path.join(os.path.dirname(sys.argv[1]),'bin'))
import FreeCAD as App
import Part
from PySide import QtCore
app=QtCore.QCoreApplication.instance() or QtCore.QCoreApplication([])
spec=importlib.util.spec_from_file_location('hornsculpt_bridge',sys.argv[2])
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
document=App.newDocument('HornSculptLiveTest')
box=document.addObject('Part::Box','Box')
box.Length=10
box.Width=12
box.Height=14
document.recompute()
bridge=module.start(int(sys.argv[3]),sys.argv[4],object_names=['Box'])
def no_op_recompute():
    for unused in range(10):
        document.recompute()
def edit():
    assert bridge.revision == 1, 'Unchanged recomputes published another snapshot'
    box.Length=20
    document.recompute()
QtCore.QTimer.singleShot(700,no_op_recompute)
QtCore.QTimer.singleShot(1500,edit)
QtCore.QTimer.singleShot(3500,app.quit)
app.exec()
bridge.stop()
module.stop()
module.stop()
assert not hasattr(App,'_hornsculpt_live_bridge')
)PY");
    script.close();
    QProcess process;
    auto environment = QProcessEnvironment::systemEnvironment();
    environment.insert("PYTHONDONTWRITEBYTECODE", "1");
    process.setProcessEnvironment(environment);
    process.start(
        python, {scriptPath, freecadLib, bridgeScript, QString::number(link.port()), link.token()});
    waitUntil([&] { return process.state() == QProcess::Running; });
    waitUntil([&] { return bool(link.snapshot()); }, 15000);
    const auto first = link.snapshot();
    require(first->faces.size() == 6 && first->geometry && first->brepObjects.size() == 1 &&
                first->mesh.triangles.size() >= 12,
            "Real FreeCAD box supplies exact BRep and a native derived preview");
    const auto id = first->faces.front().id;
    require(link.setBoundaries({id}, BoundaryCondition::driven, {0.1, 0.2}), "real CAD assignment");
    waitUntil([&] { return link.snapshot()->revision > first->revision; }, 10000);
    const auto second = link.snapshot();
    double maxX = -1.0;
    for (const auto &point : second->mesh.vertices)
        maxX = std::max(maxX, point.x);
    require(std::abs(maxX - 0.02) < 1e-12, "live dimensional edit arrives in metres");
    const auto face = std::find_if(second->faces.begin(), second->faces.end(),
                                   [id](const CadFace &value) { return value.id == id; });
    require(face != second->faces.end() && face->boundary == BoundaryCondition::driven &&
                face->normalVelocity == std::complex<double>(0.1, 0.2),
            "real edit retains unique face assignment");
    require(!second->hasUnresolvedAssignments(), "real edit has unambiguous identity");
    waitUntil([&] { return process.state() == QProcess::NotRunning; }, 10000);
    if (process.exitCode() != 0)
        std::cerr << process.readAllStandardError().constData();
    require(process.exitCode() == 0, "real FreeCAD bridge exits cleanly");
    std::cout << "Real FreeCAD live-edit integration passed.\n";
}
} // namespace

void cadMeshingSnapshot() {
    auto cad = testing::exactCadFixture();
    const auto originalPreview = cad.mesh;
    const auto first = makeCadMeshingSource(cad, .00025);
    require(first.parameters && first.parameters->triangles().size() == first.mesh.triangles.size(),
            "Native CAD meshing lacks its source-corner parameter snapshot");
    const auto firstSource = first.parameters->sourceMesh();
    const auto driven = cad.faces.front().id;
    cad.faces.front().boundary = BoundaryCondition::driven;
    // A display preview is deliberately not a source-coordinate cache.
    cad.mesh.vertices.front().x += 1;
    const auto second = makeCadMeshingSource(cad, .0001);
    require(second.parameters && second.mesh.vertices == originalPreview.vertices &&
                first.parameters->sourceMesh() == firstSource,
            "A later preview or source build changed an immutable parameter snapshot");
    bool sawDriven = false;
    for (const auto &triangle : second.mesh.triangles) {
        require(triangle.boundary == (triangle.faceId == driven ? BoundaryCondition::driven
                                                                : BoundaryCondition::wall),
                "Fresh CAD source lost its assigned physical face conditions");
        sawDriven = sawDriven || triangle.boundary == BoundaryCondition::driven;
    }
    require(sawDriven, "Fresh CAD source lost the assigned driven patch");

    const auto rejects = [](const auto &operation) {
        try {
            operation();
        } catch (const std::invalid_argument &) {
            return true;
        }
        return false;
    };
    auto missing = cad;
    missing.faces.pop_back();
    require(rejects([&] { (void)makeCadMeshingSource(missing, .00025); }),
            "An incomplete CAD face assignment was silently accepted");
    auto unprepared = cad;
    unprepared.geometry.reset();
    require(rejects([&] { (void)makeCadMeshingSource(unprepared, .00025); }),
            "An exact source without its prepared provider fell back to the display mesh");
    std::stop_source stop;
    stop.request_stop();
    bool cancelled = false;
    try {
        (void)makeCadMeshingSource(cad, .00025, stop.get_token());
    } catch (const std::runtime_error &) {
        cancelled = true;
    }
    require(cancelled, "Cancelled CAD source generation continued");
    auto polygon = testing::exactCadFixture();
    polygon.geometry.reset();
    polygon.brepObjects.clear();
    const auto linear = makeCadMeshingSource(polygon, .00025);
    require(!linear.parameters && linear.mesh == polygon.mesh,
            "A polygon-only archive was given invented native parameter data");
}

void archivedSnapshotRestore() {
    FreeCadLink donor;
    require(donor.listen(), "Archive fixture listener failed");
    QTcpSocket source;
    connectClient(donor, source, "archive-source");
    source.write(encode(snapshot("archive-source", 7)));
    waitUntil([&] { return static_cast<bool>(donor.snapshot()); });
    const auto driven = donor.snapshot()->faces.front().id;
    require(donor.setBoundaries({driven}, BoundaryCondition::driven, {.125, -.375}),
            "Archive fixture assignment failed");
    auto archived = std::make_shared<CadSnapshot>(*donor.snapshot());
    archived->unresolvedPreviousFaces.push_back(1000);

    auto cacheless = *archived;
    cacheless.geometry.reset();
    const auto verified = validateCadSnapshot(cacheless);
    require(verified && !cacheless.geometry,
            "Snapshot validation failed to return its provider without mutating the input");
    const auto &reference = cacheless.faces.front();
    const auto lifted = verified->projectSurface(reference.id, reference.signature.center);
    require(std::hypot(lifted.position.x - reference.signature.center.x,
                       lifted.position.y - reference.signature.center.y,
                       lifted.position.z - reference.signature.center.z) < 1e-12,
            "The validated provider does not preserve the source face binding");
    auto changedUnits = *archived;
    changedUnits.brepMetresPerUnit *= 2;
    bool unitsRejected = false;
    try {
        (void)validateCadSnapshot(changedUnits);
    } catch (const std::invalid_argument &) {
        unitsRejected = true;
    }
    require(unitsRejected, "An existing cache hid changed exact-source units");

    FreeCadLink restored;
    auto mutableOriginal = std::make_shared<CadSnapshot>(*archived);
    const auto preparedArchive = PreparedCadSnapshot::prepare(mutableOriginal);
    mutableOriginal->faces.front().boundary = BoundaryCondition::wall;
    mutableOriginal->faces.front().normalVelocity = {};
    mutableOriginal->brepMetresPerUnit *= .5;
    mutableOriginal->brepObjects.front().data.clear();
    restored.restoreSnapshot(preparedArchive);
    FreeCadLink secondRestore;
    secondRestore.restoreSnapshot(preparedArchive);
    require(restored.snapshot()->brepMetresPerUnit == archived->brepMetresPerUnit &&
                restored.snapshot()->brepObjects.front().data == archived->brepObjects.front().data,
            "Mutable input aliases changed a prepared archive's exact geometry");
    require(!restored.isConnected() && !restored.isListening() && !restored.isGeometryPending(),
            "An archive must not claim a live connection or pending live revision");
    require(restored.snapshot()->revision == 7 &&
                restored.snapshot()->unresolvedPreviousFaces == std::vector<std::uint32_t>{1000},
            "Restored CAD revision or unresolved identities changed");
    require(restored.snapshot()->faces.front().id == driven &&
                restored.snapshot()->faces.front().normalVelocity ==
                    std::complex<double>(.125, -.375),
            "Restored complex Driven assignment changed");
    const auto beforeFailure = restored.snapshot();
    require(restored.setBoundaries({driven}, BoundaryCondition::wall) &&
                secondRestore.snapshot()->faces.front().boundary == BoundaryCondition::driven &&
                preparedArchive.snapshot()->faces.front().boundary == BoundaryCondition::driven,
            "Editing one restored link mutated another link or the prepared archive");
    restored.restoreSnapshot(preparedArchive);
    require(!restored.setBoundaries({driven}, static_cast<BoundaryCondition>(2)) &&
                restored.snapshot() == beforeFailure,
            "Unsupported boundary condition changed the CAD archive");
    auto invalid = std::make_shared<CadSnapshot>(*archived);
    invalid->faces.front().id = 0;
    bool rejected = false;
    try {
        restored.restoreSnapshot(invalid);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected && restored.snapshot() == beforeFailure,
            "Invalid archive discarded the active CAD document");

    const std::vector<std::uint32_t> group{archived->faces[0].id, archived->faces[1].id};
    int groupCommits = 0;
    QObject::connect(&restored, &FreeCadLink::snapshotChanged, [&] { ++groupCommits; });
    require(!restored.setBoundaries({group[0], 999999}, BoundaryCondition::driven, {.3, .4}) &&
                restored.snapshot() == beforeFailure && groupCommits == 0,
            "An invalid group edit partially changed a valid face");
    require(!restored.setBoundaries({group[0], group[0]}, BoundaryCondition::driven) &&
                !restored.confirmFaces({group[0], 999999}) && groupCommits == 0,
            "Duplicate or unknown group IDs were accepted");
    require(restored.setBoundaries(group, BoundaryCondition::driven, {.3, -.4}) &&
                groupCommits == 1,
            "A group edit did not commit exactly one snapshot");
    for (const auto &face : restored.snapshot()->faces) {
        const bool selected = std::find(group.begin(), group.end(), face.id) != group.end();
        require(selected ? face.boundary == BoundaryCondition::driven &&
                               face.normalVelocity == std::complex<double>(.3, -.4)
                         : face.boundary == BoundaryCondition::wall,
                "Group assignment changed the wrong CAD faces");
    }
    auto review = std::make_shared<CadSnapshot>(*restored.snapshot());
    for (auto &face : review->faces)
        if (std::find(group.begin(), group.end(), face.id) != group.end())
            face.state = CadFaceState::needsReview;
    restored.restoreSnapshot(review);
    groupCommits = 0;
    require(restored.confirmFaces(group) && groupCommits == 1,
            "Group confirmation did not commit exactly once");
    for (const auto &face : restored.snapshot()->faces)
        if (std::find(group.begin(), group.end(), face.id) != group.end())
            require(face.state == CadFaceState::confirmed &&
                        face.normalVelocity == std::complex<double>(.3, -.4),
                    "Group confirmation replaced an existing boundary definition");

    std::vector<CadFaceAssignment> assignments;
    for (const auto &face : review->faces)
        assignments.push_back({face.id, face.boundary, face.normalVelocity, face.state});
    const auto beforeReplay = restored.snapshot();
    require(!restored.setAssignments(review->documentId, review->sessionId, review->revision + 1,
                                     assignments, {1000}) &&
                restored.snapshot() == beforeReplay,
            "History replay crossed a CAD revision");
    auto duplicateAssignments = assignments;
    duplicateAssignments.back() = duplicateAssignments.front();
    require(!restored.setAssignments(review->documentId, review->sessionId, review->revision,
                                     duplicateAssignments, {1000}) &&
                restored.snapshot() == beforeReplay,
            "Invalid history replay partially committed");
    require(!restored.setAssignments(review->documentId, review->sessionId, review->revision,
                                     assignments, {assignments[0].id}) &&
                restored.snapshot() == beforeReplay,
            "Existing CAD face accepted as missing");
    groupCommits = 0;
    require(restored.setAssignments(review->documentId, review->sessionId, review->revision,
                                    assignments, {1000}) &&
                groupCommits == 1 && restored.snapshot()->hasUnresolvedAssignments(),
            "CAD history replay lost review metadata or committed more than once");
    require(restored.snapshot()->mesh.vertices.size() == beforeReplay->mesh.vertices.size() &&
                restored.snapshot()->mesh.triangles.size() == beforeReplay->mesh.triangles.size(),
            "Assignment history replay replaced geometry");

    require(restored.listen(), "Restored listener failed");
    QTcpSocket later;
    connectClient(restored, later, "archive-later");
    later.write(encode(snapshot("archive-later", 1, false, false, false, "another-document")));
    waitUntil([&] { return restored.snapshot()->documentId == "another-document"; });
    for (const auto &face : restored.snapshot()->faces)
        require(face.id > 1000, "Reconnection reused an archived or missing face identity");

    restored.restoreSnapshot({});
    require(!restored.snapshot() && !restored.isConnected() && !restored.isListening(),
            "New project retained archived CAD state");
    require(restored.listen(), "Reset listener failed");
    QTcpSocket fresh;
    connectClient(restored, fresh, "archive-fresh");
    fresh.write(encode(snapshot("archive-fresh", 1)));
    waitUntil([&] { return static_cast<bool>(restored.snapshot()); });
    const auto first =
        std::min_element(restored.snapshot()->faces.begin(), restored.snapshot()->faces.end(),
                         [](const CadFace &a, const CadFace &b) { return a.id < b.id; });
    require(first->id == 1, "New project retained identities from a previous project");
}

std::shared_ptr<const CadSnapshot> snapshotFile(const QString &path,
                                                const QString &outputProject = {}) {
    QFile file(path);
    require(file.open(QIODevice::ReadOnly), "Snapshot fixture cannot be opened");
    QJsonParseError error;
    const auto document = QJsonDocument::fromJson(file.readAll(), &error);
    require(error.error == QJsonParseError::NoError && document.isObject(),
            "Snapshot fixture is not a JSON object");
    auto payload = document.object();
    require(payload["protocol"].toInt() == 2 && !payload.contains("vertices") &&
                !payload.contains("triangles"),
            "Snapshot fixture must contain an exact protocol-2 source");
    const auto expectedFaces = payload["faces"].toArray().size();
    const auto expectedParts = payload["brepParts"].toArray().size();
    require(expectedFaces > 0 && expectedParts > 0, "Snapshot fixture has no CAD faces");
    const auto session = QStringLiteral("external-snapshot-test");
    payload.insert("sessionId", session);
    payload.insert("revision", "1");
    FreeCadLink link;
    require(link.listen(), "Snapshot fixture listener cannot start");
    QTcpSocket client;
    connectClient(link, client, session);
    client.write(encode(payload));
    waitUntil([&] { return bool(link.snapshot()) || !link.isConnected(); }, 30000);
    if (!link.snapshot())
        throw std::runtime_error("Exact snapshot rejected: " + link.lastError().toStdString());
    const auto original = link.snapshot();
    require(original->geometry &&
                original->faces.size() == static_cast<std::size_t>(expectedFaces) &&
                original->brepObjects.size() == static_cast<std::size_t>(expectedParts),
            "Native snapshot lost an exact CAD part or face");
    const auto quality = inspectMesh(original->mesh);
    require(quality.watertight() && quality.signedVolume > 0,
            "Native snapshot preview is not a closed oriented boundary");
    std::vector<std::uint32_t> selected{original->faces.front().id};
    if (original->faces.size() > 1)
        selected.push_back(original->faces.back().id);
    int commits = 0;
    QObject::connect(&link, &FreeCadLink::snapshotChanged, &link, [&] { ++commits; });
    require(link.setBoundaries(selected, BoundaryCondition::driven, {.2, -.3}) && commits == 1,
            "Group boundary edit did not commit atomically");
    const auto updated = link.snapshot();
    require(updated->geometry == original->geometry && updated->revision == original->revision &&
                updated->mesh.vertices == original->mesh.vertices &&
                updated->mesh.triangles.size() == original->mesh.triangles.size(),
            "Group boundary edit changed the native geometry or its preview");
    for (std::size_t i = 0; i < original->brepObjects.size(); ++i)
        require(updated->brepObjects[i].data == original->brepObjects[i].data &&
                    updated->brepObjects[i].faceIds == original->brepObjects[i].faceIds,
                "Group boundary edit changed the exact source or face binding");
    for (std::size_t i = 0; i < updated->mesh.triangles.size(); ++i) {
        const auto &a = updated->mesh.triangles[i], &b = original->mesh.triangles[i];
        require(a.vertices == b.vertices && a.faceId == b.faceId,
                "Group boundary edit changed mesh connectivity");
    }
    for (const auto &face : updated->faces)
        if (std::find(selected.begin(), selected.end(), face.id) != selected.end())
            require(face.boundary == BoundaryCondition::driven &&
                        face.normalVelocity == std::complex<double>(.2, -.3),
                    "Group boundary edit lost its complex normal velocity");
    if (!outputProject.isEmpty()) {
        Project project;
        project.source = ProjectSource::freecad;
        project.revision = project.geometryRevision = updated->revision;
        project.cad = updated;
        project.preview = std::shared_ptr<const SurfaceMesh>(updated, &updated->mesh);
        writeProject(outputProject, project);
        const auto restored = readProject(outputProject);
        require(restored.cad && restored.cad->geometry && restored.preview &&
                    *restored.preview == updated->mesh && !restored.mesh &&
                    !restored.simulation.result,
                "Exported exact CAD project did not preserve its source-only state");
        std::cout << "Exact CAD project exported without a generated BEM mesh or result.\n";
    }
    std::cout << "Exact protocol-2 snapshot passed: " << expectedFaces << " faces, "
              << expectedParts << " BRep parts, " << original->mesh.triangles.size()
              << " closed preview triangles; group BC preserves exact geometry.\n";
    return original;
}

void extrusionSurfaceType() {
    const QDir directory(QFileInfo(QString::fromUtf8(__FILE__)).absolutePath());
    const auto snapshot =
        snapshotFile(directory.filePath("fixtures/freecad-spline-extrusion.json"));
    auto changed = *snapshot;
    const auto found =
        std::find_if(changed.faces.begin(), changed.faces.end(), [](const auto &face) {
            return face.signature.surfaceType == "SurfaceOfExtrusion";
        });
    require(found != changed.faces.end(), "FreeCAD fixture has no native extrusion face");
    found->signature.surfaceType = "SurfaceOfRevolution";
    bool rejected = false;
    try {
        (void)validateCadSnapshot(changed);
    } catch (const std::invalid_argument &error) {
        const std::string message = error.what();
        rejected = message.find("SurfaceOfRevolution") != std::string::npos &&
                   message.find("SurfaceOfLinearExtrusion") != std::string::npos &&
                   message.find(std::to_string(found->id)) != std::string::npos;
    }
    require(rejected,
            "A different CAD surface type was accepted or its mismatch was not identified");
}


int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    try {
        const auto arguments = app.arguments();
        if ((argc == 3 || argc == 4) && arguments[1] == "--snapshot-file") {
            snapshotFile(arguments[2], argc == 4 ? arguments[3] : QString{});
            return 0;
        }
        require(
            argc == 1 || argc == 4,
            "Use no arguments, FreeCAD's three paths, or --snapshot-file <path> [output.hsproj]");
        protocolAndIdentity();
        cadMeshingSnapshot();
        archivedSnapshotRestore();
        extrusionSurfaceType();

        if (argc == 4)
            realBridge(QString::fromLocal8Bit(argv[1]), QString::fromLocal8Bit(argv[2]),
                       QString::fromLocal8Bit(argv[3]));
        std::cout << "FreeCAD framing, revision, face mapping and reconnect tests passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
