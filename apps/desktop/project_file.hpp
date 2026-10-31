#pragma once

#include "angular_response_settings.hpp"
#include "field_display_settings.hpp"
#include "gmsh_import_data.hpp"
#include "measured_response.hpp"
#include "mesh_settings.hpp"
#include <QString>
#include <hornsculpt/acoustic_frame.hpp>
#include <hornsculpt/acoustic_partition.hpp>
#include <hornsculpt/bem.hpp>
#include <hornsculpt/draft_geometry.hpp>
#include <hornsculpt/draft_instances.hpp>
#include <hornsculpt/face_color.hpp>
#include <hornsculpt/freecad_link.hpp>
#include <hornsculpt/meshing.hpp>
#include <map>
#include <memory>
#include <optional>
#include <set>

namespace hornsculpt {

enum class ProjectSource { draft, freecad };

struct AxisymmetricSamplingSettings {
    double maximumLength{.002};
    double chordTolerance{.00005};
    double rotationalTolerance{1e-8};
    unsigned sectors{64};
    bool operator==(const AxisymmetricSamplingSettings &) const = default;
};

struct SimulationSettings {
    double minimumFrequency{20};
    double maximumFrequency{20000};
    std::size_t frequencyCount{17};
    AcousticMedium medium;
    double velocityScale{1};
    double phaseDegrees{};
    bool normalizeToOnAxis{true};
    bool beamwidthTargetEnabled{};
    double targetHorizontalDegrees{90};
    double targetVerticalDegrees{60};
    double targetToleranceDegrees{10};
    double normalizationAngleDegrees{0};
    double directivityMinimumDb{-30};
    double directivityMaximumDb{3};
    double directivityStepDb{3};
    // Empty follows the current geometry's source frame without storing a copy.
    std::optional<AcousticFrame> observationFrame;
    bool farField{true};
    double observationDistance{1}; // m, from the explicit observation frame origin.
    bool couplingEnabled{};
    lem::Network couplingNetwork;
    BemAcceleration acceleration{BemAcceleration::dense};
    FastMultipoleSettings fastMultipole;
    unsigned nearMaximumDepth{6};
    std::optional<InfiniteBaffle> infiniteBaffle;
    AngularResponseSettings angularResponses;
    FieldDisplaySettings fieldDisplay;
    PhaseDisplaySettings phaseDisplay;
    BemGeometryMode geometryMode{BemGeometryMode::threeDimensional};
    AxisymmetricSamplingSettings axisymmetricSampling;
    int ctaExcitationIndex{-1}; // display selection only
    unsigned singularOrder{4};
    AcousticModelDefinition acousticModel;
    bool symmetryXYEnabled{};
    double symmetryPlaneTolerance{1e-6}; // m
};

struct SimulationArchive {
    SimulationSettings settings;
    // These snapshots describe the computed result, even after current inputs change.
    std::shared_ptr<const SurfaceMesh> mesh;
    std::optional<BemRequest> request;
    std::shared_ptr<const BemResult> result;
    std::uint64_t revision{};
    bool current{};
    std::vector<MeasuredResponse> measurements;
    // Sampling belongs to the computed integration mesh, not subsequently edited inputs.
    std::optional<AxisymmetricSamplingSettings> axisymmetricSampling;
    // Missing in older archives: completeness of a stale result cannot be inferred.
    std::optional<bool> resultComplete;
    std::shared_ptr<const SurfaceMesh> sourceMesh;
    AcousticModelDefinition resultModel;
};

struct Project;
// Comparison entries own an immutable project without further comparison entries.
struct NamedProjectSnapshot {
    QString name;
    std::shared_ptr<const Project> project;
};

struct Project {
    DraftParameters parameters;
    DraftAssembly assembly;
    MeshSettings meshing;
    std::map<std::uint32_t, double> draftFaceTargetEdgeLengths;
    std::map<std::uint32_t, double> cadFaceTargetEdgeLengths;
    std::set<std::uint32_t> draftReversedFaceIds;
    std::set<std::uint32_t> cadReversedFaceIds;
    ProjectSource source{ProjectSource::draft};
    std::uint64_t revision{1};
    std::uint64_t geometryRevision{};
    std::uint64_t meshRevision{};
    std::shared_ptr<const SurfaceMesh> preview;
    std::shared_ptr<const MeshingResult> mesh;
    std::shared_ptr<const CadSnapshot> cad;
    std::optional<MeshImportSelection> gmsh;
    SimulationArchive simulation;
    FaceColors draftFaceColors;
    FaceColors cadFaceColors;
    std::vector<NamedProjectSnapshot> comparisonSnapshots;
};

// Reconciled physical IDs come from FreeCadLink, never ordinal FaceN names.
// A new document, missing/object-changed or ambiguous face loses its color.
FaceColors retainedCadFaceColors(const FaceColors &colors, const CadSnapshot *previous,
                                 const CadSnapshot *current);

// Losing an assigned physical surface keeps a persistent review obligation.
// A numeric ID in a new document is never evidence for a BEM-port match.
AcousticModelDefinition reconciledCadAcousticModel(AcousticModelDefinition model,
                                                   const CadSnapshot *previous,
                                                   const CadSnapshot *current);
lem::Network reconciledCadBemAssignments(lem::Network network, const CadSnapshot *previous,
                                         const CadSnapshot *current);

void writeProject(const QString &path, const Project &project);
Project readProject(const QString &path);

} // namespace hornsculpt
