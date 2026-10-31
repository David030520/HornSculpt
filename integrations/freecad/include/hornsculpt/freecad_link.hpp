#pragma once

#include <hornsculpt/mesh.hpp>
#include <hornsculpt/surface_parameters.hpp>

#include <QByteArray>
#include <QObject>
#include <QString>
#include <array>
#include <complex>
#include <cstdint>
#include <memory>
#include <stop_token>
#include <vector>

namespace hornsculpt {
class CadGeometry;

enum class CadFaceState { confirmed, newFace, needsReview };

struct CadFaceAssignment {
    std::uint32_t id{};
    BoundaryCondition boundary{BoundaryCondition::wall};
    std::complex<double> normalVelocity{1, 0};
    CadFaceState state{CadFaceState::newFace};
    bool operator==(const CadFaceAssignment &) const = default;
};

struct CadFaceSignature {
    QString surfaceType;
    QString geometryKey;
    QString topologyKey;
    double area{};
    Vec3 center;
    Vec3 normal;
    std::uint32_t edgeCount{};
    std::uint32_t wireCount{};
};

struct CadFace {
    std::uint32_t id{};
    QString objectId;
    QString sourceKey;
    QString label;
    CadFaceSignature signature;
    BoundaryCondition boundary{BoundaryCondition::wall};
    std::complex<double> normalVelocity{1.0, 0.0};
    CadFaceState state{CadFaceState::newFace};
};

struct CadObject {
    QString id;
    QString name;
    QString label;
    QString topologyFingerprint;
    // Row-major rigid object placement, translation in metres. Mesh is already global.
    std::array<double, 16> placement{};
};

struct CadBrepObject {
    QString objectId;
    QByteArray data;
    // Traversal order is local to this exact artifact; values are reconciled physical IDs.
    std::vector<std::uint32_t> faceIds;
};

struct CadSnapshot {
    QString documentId;
    QString documentName;
    QString sessionId;
    std::uint64_t revision{};
    SurfaceMesh mesh;
    std::vector<CadBrepObject> brepObjects;
    double brepMetresPerUnit{.001};
    std::shared_ptr<const CadGeometry> geometry;
    std::vector<CadFace> faces;
    std::vector<CadObject> objects;
    Vec3 acousticOrigin;
    Vec3 acousticAxis{0.0, 0.0, 1.0};
    std::vector<std::uint32_t> unresolvedPreviousFaces;
    bool hasUnresolvedAssignments() const noexcept;
};

// Validates a mapped snapshot independent of its live or archived origin.
// Returns the newly verified immutable provider for reuse without a second import;
// polygon-only archives return null. The input snapshot is never modified.
[[nodiscard]] std::shared_ptr<const CadGeometry> validateCadSnapshot(const CadSnapshot &snapshot);
// Rebuild the native cache and validate its complete face binding. Archived preview
// data remains unchanged unless explicitly requested for a new live snapshot.
void prepareCadGeometry(CadSnapshot &snapshot, bool buildPreview);
// Fresh source and parameter bindings at the requested detail tolerance.
// Assigned face conditions are copied from this immutable CAD snapshot.
MeshingSource makeCadMeshingSource(const CadSnapshot &snapshot, double geometricTolerance,
                                   std::stop_token stop = {});

// A privately owned snapshot validated before crossing into the GUI thread.
// Public CadSnapshot fields and an existing geometry pointer are not validation proof.
class PreparedCadSnapshot final {
  public:
    static PreparedCadSnapshot prepare(std::shared_ptr<const CadSnapshot> snapshot);
    std::shared_ptr<const CadSnapshot> snapshot() const noexcept;

  private:
    explicit PreparedCadSnapshot(std::shared_ptr<const CadSnapshot> snapshot);
    std::shared_ptr<const CadSnapshot> snapshot_;
};

// Local-only authenticated revision stream. All methods and signals use its Qt thread.
// A snapshot is committed only after the entire frame, geometry and mapping validate.
class FreeCadLink final : public QObject {
    Q_OBJECT
  public:
    explicit FreeCadLink(QObject *parent = nullptr);
    ~FreeCadLink() override;
    bool listen(quint16 requestedPort = 0);
    void stop();
    bool isListening() const;
    bool isConnected() const;
    bool isGeometryPending() const;
    quint16 port() const;
    QString token() const;
    QString discoveryFilePath() const;
    QString lastError() const;
    std::shared_ptr<const CadSnapshot> snapshot() const;
    // Validate a project snapshot before replacing state and disconnecting live input.
    // Passing null clears project-local face identities and pending revisions.
    void restoreSnapshot(std::shared_ptr<const CadSnapshot> snapshot);
    void restoreSnapshot(const PreparedCadSnapshot &snapshot);
    // Group edits commit once; an unknown or repeated ID leaves every face unchanged.
    bool setBoundaries(const std::vector<std::uint32_t> &faceIds, BoundaryCondition boundary,
                       std::complex<double> normalVelocity = {1.0, 0.0});
    bool confirmFaces(const std::vector<std::uint32_t> &faceIds);
    // Atomic edit-history replay on precisely the same live geometry snapshot.
    bool setAssignments(const QString &documentId, const QString &sessionId, std::uint64_t revision,
                        const std::vector<CadFaceAssignment> &faces,
                        const std::vector<std::uint32_t> &missing);
    void acknowledgeMissingFaces();
  signals:
    void connectionChanged(bool connected);
    void snapshotChanged();
    void sourceDirtyChanged(bool dirty);
    void errorOccurred(const QString &message);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace hornsculpt
