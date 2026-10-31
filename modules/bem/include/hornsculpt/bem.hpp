#pragma once

#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <hornsculpt/acoustic_domains.hpp>
#include <hornsculpt/exterior_baffle.hpp>
#include <hornsculpt/lem.hpp>
#include <hornsculpt/mesh.hpp>
#include <hornsculpt/symmetry_planes.hpp>
#include <memory>
#include <optional>
#include <stop_token>
#include <string>
#include <utility>
#include <vector>

namespace hornsculpt {

// Independent pressure-feedback consistency budget, not a residual-derived
// error estimate. Shared by the solver and persisted-result validation.
inline constexpr double relativePortPressureClosureTolerance = 1e-3;
// Separate interpolation-convergence budget; not a certified remainder bound.
inline constexpr double fastMultipoleOrderComparisonBudgetFraction = .25;

struct AcousticMedium {
    double density{1.21};
    double soundSpeed{343.0};
};

struct FaceVelocity {
    std::uint32_t faceId{};
    // RMS m/s, entering BEM air: out of solid, or opposite the closed air-cavity normal.
    std::complex<double> normalVelocity{1.0, 0.0};
};

struct BemExcitation {
    std::string name;
    std::vector<FaceVelocity> faces;
};

enum class BemAcceleration { dense, fastMultipole };

enum class BemResidualEvaluation { directGalerkin, higherOrderFastMultipole };

enum class BemFastMultipoleRepresentation { regularSpherical, planeWave };
struct BemFastMultipoleLevelEvidence {
    unsigned level{};
    BemFastMultipoleRepresentation representation{BemFastMultipoleRepresentation::regularSpherical};
    unsigned momentBandlimit{}, translationOrder{}, quadratureOrder{};
    bool operator==(const BemFastMultipoleLevelEvidence &) const = default;
};
// Historical result provenance only. These tensor-Chebyshev orders must never
// be interpreted as spherical/plane-wave settings of an automatic wideband solve.
struct BemHistoricalChebyshevEvidence {
    unsigned comparisonOrder{}, auditOrder{};
    bool operator==(const BemHistoricalChebyshevEvidence &) const = default;
};

struct FastMultipoleSettings {
    // Automatic spherical/plane-wave orders derive from the requested residual
    // accuracy. A leaf stops when each channel separately fits this capacity:
    // complete P1 pressure supports and DP0 source panels, not their sum.
    unsigned leafSupportCapacity{16};
    double separationRatio{1.5};
    unsigned maximumIterations{300};
    unsigned restartLength{40};
    // Spatial P1 Near-CSR block-Jacobi width; powers of two from 32 to 2048 DOFs.
    unsigned preconditionerBlockWidth{32};
    bool operator==(const FastMultipoleSettings &) const = default;
};

struct BemHistoricalChebyshevStatistics {
    std::size_t directPanelPairCount{}, representedFarPanelPairs{};
    bool operator==(const BemHistoricalChebyshevStatistics &) const = default;
};
// Recorded once per frequency across every RHS, seed and quadrature attempt.
// Absence means not recorded, including older archives and Dense results.
// Metadata bytes describe owned host payload only; CUDA graph-driver storage
// is opaque. Free-device memory is an observation, not graph-owned memory.
struct BemFastMultipoleGraphEvidence {
    std::uint64_t captureAttempts{};
    std::uint64_t instantiations{};
    std::uint64_t launches{};
    std::uint64_t reuses{};
    std::uint64_t executedBatches{};
    std::uint64_t canonicalChunks{};
    std::uint64_t canonicalBatches{};
    std::uint64_t invalidatedExecutables{};
    std::uint64_t budgetBypasses{};
    std::uint64_t admissionClosures{};
    std::uint64_t maximumCachedExecutables{};
    std::uint64_t maximumCachedNodes{};
    std::uint64_t peakOwnedHostMetadataBytes{};
    double captureMilliseconds{}, instantiationMilliseconds{};
    bool memoryObservationAvailable{};
    std::uint64_t minimumObservedFreeDeviceBytes{};
    std::uint64_t totalDeviceBytes{};
    std::uint64_t requiredFreeReserveBytes{};
    bool operator==(const BemFastMultipoleGraphEvidence &) const = default;
};
struct FastMultipoleStatistics {
    std::size_t nodeCount{}, leafCount{}, farInteractionCount{};
    // Complete, disjoint ordered DOF pairs: near+far=Np*Np and Np*Nt.
    // These are not triangle-pair counts and must not be validated against Nt^2.
    std::uint64_t pressurePressureNearPairs{}, pressurePressureFarPairs{};
    std::uint64_t pressureFluxNearPairs{}, pressureFluxFarPairs{};
    std::optional<BemHistoricalChebyshevStatistics> historicalChebyshevStatistics;
    std::optional<BemFastMultipoleGraphEvidence> graphExecution;
};

struct BemNumerics {
    unsigned quadratureOrder{4};
    unsigned singularOrder{4};
    double nearDistanceRatio{1.0};
    double nearRelativeTolerance{1e-4};
    unsigned nearMaximumDepth{6};
    double relativeResidualTolerance{5e-5};
    unsigned maximumRefinementSteps{3};
    BemAcceleration acceleration{BemAcceleration::dense};
    FastMultipoleSettings fastMultipole;
};

struct FarFieldObservation {
    Vec3 origin{};
    double referenceDistance{1.0}; // metres; includes outgoing phase exp(-ik referenceDistance).
    bool operator==(const FarFieldObservation &) const = default;
};

// A recessed air cavity is closed by these planar interface faces. The plane
// normal points into the radiating half-space; the aperture is not a rigid wall.
struct InfiniteBaffle {
    Vec3 origin{};
    Vec3 normal{0.0, 0.0, 1.0};
    std::vector<std::uint32_t> apertureFaceIds;
    // Separate, planar open patches radiating directly into the same half-space.
    // World coordinates; Driven tags and face IDs disjoint from the cavity mesh.
    SurfaceMesh exteriorDrivenMesh;
    bool operator==(const InfiniteBaffle &) const = default;
};

// Historical triangular traces retain their original surface measure when
// opened. New m=0 calculations use the continuous revolution of the meridian.
enum class AxisymmetricSurfaceContract { triangulatedOrbit, continuousRevolution };

// Explicit m=0 Galerkin space. The full periodic mesh carries the meridian
// topology and reconstructed traces; surfaceContract identifies the physical
// integration surface. Coordinates retain the public world frame.
struct AxisymmetricBem {
    Vec3 origin{};
    Vec3 axis{0.0, 0.0, 1.0};
    unsigned sectors{64};
    std::vector<std::uint32_t> pressureDofs; // One entry per full-mesh vertex.
    std::vector<std::uint32_t> fluxDofs;     // One entry per full-mesh triangle.
    std::vector<std::uint32_t> testPanels;   // All triangles of one angular sector.
    AxisymmetricSurfaceContract surfaceContract{AxisymmetricSurfaceContract::continuousRevolution};
    bool operator==(const AxisymmetricBem &) const = default;
};

enum class BemGeometryMode { threeDimensional, axisymmetric };

struct BemRequest {
    AcousticMedium medium;
    std::optional<BemSubdomainModel> subdomains;
    std::optional<OrthogonalSymmetryPlanes> symmetry;
    std::optional<ExteriorBaffle> exteriorBaffle;
    std::vector<double> frequencies;
    std::vector<BemExcitation> excitations;
    std::vector<Vec3> observers;
    // Empty means all observers are valid. For the exterior-baffle far field,
    // a nonempty 0/1 mask must identify exactly the exterior half-space.
    std::vector<std::uint8_t> observerValidity;
    BemNumerics numerics;
    // World observer points specify directions from origin. Results are the
    // asymptotic outgoing pressure normalized to referenceDistance, not finite
    // field samples at those observer points. Null retains finite observation.
    std::optional<FarFieldObservation> farField;
    // Exterior-solid coupling replaces prescribed excitations. Baffled-air
    // coupling also permits one independent physical Driven excitation; its
    // blocked port pressure is included before solving the LEM. The coupled
    // result is named "LEM". Driven ports use their normalized velocity shape;
    // interface ports use uniform velocity over their face groups.
    // Independently prescribed faces must be disjoint from all network ports.
    std::optional<lem::Network> coupling;
    std::optional<InfiniteBaffle> infiniteBaffle;
    std::optional<AxisymmetricBem> axisymmetric;
    // Empty uses frequencies. Otherwise every solved frequency must belong to
    // this explicit full planning raster; scheduled subset calls keep it intact.
    std::vector<double> interactionPlanFrequencies;
};

struct FaceImpedance {
    std::uint32_t faceId{};
    double area{};
    std::complex<double> meanPressure;
    std::complex<double> volumeVelocity;
    std::complex<double> volumeImpedance; // Pa s/m^3, mean pressure / volume velocity.
};

// Published results own immutable boundary traces. Copies of an excitation
// share these potentially large arrays while observer data remains independent.
struct BemBoundaryData {
    // DP0 mean pressures on mesh triangles, including apertures, in mesh order.
    std::vector<std::complex<double>> surfacePressure;
    // P1 pressure in full mesh-vertex order, reconstructed also for m=0.
    std::vector<std::complex<double>> pressureTrace;
    // DP0 q=dp/dnu in triangle order, with closed-mesh outward normals:
    // -i*rho*omega*v for exterior-solid; +i*rho*omega*v for baffled
    // physical surfaces, and the solved derivative on the aperture.
    std::vector<std::complex<double>> normalDerivative;
    // Rayleigh test-cell means on exteriorDrivenMesh, not cavity P1 unknowns.
    std::vector<std::complex<double>> exteriorSurfacePressure;
};

// Domain-local traces use explicit mesh indices. Flux follows the domain's
// outward air normal. Interface pressure and flux remain spatially resolved.
struct BemSubdomainTrace {
    std::uint32_t domainId{};
    std::vector<std::uint32_t> vertexIndices;
    std::vector<std::complex<double>> pressureTrace;
    std::vector<std::uint32_t> triangleIndices;
    std::vector<std::complex<double>> normalDerivative;
    struct InterfaceFlux {
        std::uint32_t interfaceId{};
        std::vector<std::uint32_t> vertexIndices;
        std::vector<std::complex<double>> normalDerivative;
    };
    std::vector<InterfaceFlux> interfaceFluxes;
};
struct BemExcitationResult {
    std::string name;
    std::vector<std::complex<double>> observerPressure;
    // Empty means all entries are valid. Masked pressure slots are placeholders.
    std::vector<std::uint8_t> observerValidity;
    std::vector<FaceImpedance> impedances;
    unsigned refinementSteps{};
    // Relative L2 residual of the unscaled equation identified below. A
    // higher-order FMM residual is not the exact original Galerkin residual.
    double boundaryResidual{};
    BemResidualEvaluation residualEvaluation{BemResidualEvaluation::directGalerkin};
    // Last inner correction solve's unpreconditioned relative residual; zero
    // when no inner solve was necessary (e.g. a zero coupled excitation).
    double fastMultipoleKrylovResidual{};
    // Configured relative target of the left-preconditioned inner Krylov
    // equation. Not a measured original residual; zero means not recorded.
    double fastMultipoleKrylovTarget{};
    // (||(A_audit-A_comparison)p|| + ||(B_audit-B_comparison)q||) / ||B_audit q||.
    // This is an order-comparison indicator, not a bound on the original
    // operator, mesh discretization, quadrature or pressure error.
    double fastMultipoleOperatorDiscrepancy{};
    // Actual automatic wideband seed and independent final comparison/verification
    // seeds. A seed selects initial level orders; it is not an operator-error bound.
    // Dense solves have zero seeds. Direct audits retain the solve seed, while
    // comparison/verification seeds and their level lists are zero/empty.
    double fastMultipoleAccuracySeed{};
    double fastMultipoleComparisonAccuracySeed{};
    double fastMultipoleVerificationAccuracySeed{};
    unsigned fastMultipoleAccuracyRefinements{};
    std::vector<BemFastMultipoleLevelEvidence> fastMultipoleComparisonLevels;
    std::vector<BemFastMultipoleLevelEvidence> fastMultipoleVerificationLevels;
    std::optional<BemHistoricalChebyshevEvidence> historicalChebyshevEvidence;
    // Independent integration-rule checks at a fixed final accuracy seed.
    // Zero degree pairs mean not evaluated, not demonstrated convergence.
    unsigned farQuadratureComparisonDegree{}, farQuadratureVerificationDegree{};
    unsigned farQuadratureRefinements{};
    // (||Delta Ap|| + ||Delta Bq||) / ||Bq||; separate from FMM seed error.
    double farQuadratureDiscrepancy{};
    double farQuadraturePreparationMilliseconds{}, farQuadraturePressureMilliseconds{},
        farQuadratureFluxMilliseconds{};
    unsigned fieldQuadratureComparisonDegree{}, fieldQuadratureVerificationDegree{};
    double fieldQuadratureDiscrepancy{};
    // Actual physical RMS velocities, including zero-flow coupled faces.
    std::vector<FaceVelocity> solvedVelocities;
    std::optional<lem::Solution> lemSolution;
    unsigned iterations{}; // Actual Krylov steps; zero for dense LU.
    // Null denotes absent boundary data. Legacy archives can contain surface
    // pressure but empty Cauchy traces; field evaluation must reject those.
    // P1 pressures follow full vertex order; DP0 derivatives follow triangle
    // order. Their signs follow the closed-mesh outward normal convention.
    // The value is moved or copied into a newly allocated const object so no
    // mutable shared owner survives publication.
    void setBoundary(BemBoundaryData value) {
        boundary_ = std::make_shared<const BemBoundaryData>(std::move(value));
    }
    const std::shared_ptr<const BemBoundaryData> &boundary() const noexcept {
        return boundary_;
    }
    // Copies share published traces without retaining a mutable owner.
    void setSubdomainTraces(std::vector<BemSubdomainTrace> value) {
        subdomainTraces_ = std::make_shared<const std::vector<BemSubdomainTrace>>(std::move(value));
    }
    const std::vector<BemSubdomainTrace> &subdomainTraces() const noexcept {
        static const std::vector<BemSubdomainTrace> empty;
        return subdomainTraces_ ? *subdomainTraces_ : empty;
    }

  private:
    std::shared_ptr<const BemBoundaryData> boundary_;
    std::shared_ptr<const std::vector<BemSubdomainTrace>> subdomainTraces_;
};

struct BemFrequencyResult {
    double frequency{};
    std::vector<BemExcitationResult> excitations;
    double assemblyMilliseconds{};
    double solveMilliseconds{};
    double factorizationMilliseconds{};
    double triangularSolveMilliseconds{};
    double fieldMilliseconds{};
    double residualCheckMilliseconds{};
    FastMultipoleStatistics fastMultipole;
    // Only the first frequency of one solve call carries its host-plan timing
    // and cache-hit flag. Assembly time already includes this build duration.
    double interactionPlanMilliseconds{};
    bool interactionPlanReused{};
};

struct BemResult {
    std::string device;
    std::size_t triangleCount{};
    std::size_t pressureDofCount{};
    std::size_t peakDeviceBytes{};
    std::vector<BemFrequencyResult> frequencies;
    BemAcceleration acceleration{BemAcceleration::dense};
    BemGeometryMode geometryMode{BemGeometryMode::threeDimensional};
    // Archived measurement: maximum owned stationary buffers held during
    // operator applications (including solver/FFT/RCR/preconditioner workspaces),
    // after temporary construction buffers are released. Not currently live
    // memory after the job ends. Zero means unavailable for historical data.
    std::size_t liveDeviceBytes{};
};

enum class BemSolvePhase {
    preparing,
    nearConstruction,
    farPreparation,
    krylovSolve,
    residualAudit,
    fieldEvaluation,
    frequencyComplete
};
struct BemProgress {
    std::size_t completedFrequencies{};
    std::size_t totalFrequencies{};
    double frequency{};
    BemSolvePhase phase{BemSolvePhase::preparing};
    // Completed Krylov steps and actual A*p / B*q applications, cumulative
    // within the current frequency across excitations, corrections and audits.
    unsigned iterations{};
    std::size_t pressureOperatorApplications{}, faceFluxOperatorApplications{};
    // Zero-based linear RHS index, including coupled-port basis solves when
    // present; count is explicit rather than inferred from the fixture file.
    std::size_t rightHandSideIndex{}, rightHandSideCount{};
};

// Positive-time RMS convention: exp(+i omega t), outgoing exp(-ikr), q=-i rho omega v.
// Continuous P1 pressure / DP0 flux Burton-Miller Galerkin: pivoted CUDA LU or
// hierarchical CUDA Fast Multipole with restarted GMRES and cost-selected audit.
// The closed mesh and GPU geometry data are retained across frequency sweeps.
enum class CudaCacheResidency { automatic, device, mappedHost };
struct CudaStorageOptions {
    // Storage only: the immutable quadrature patches and numerical request are unchanged.
    CudaCacheResidency quadraturePatches{CudaCacheResidency::automatic};
};

class CudaBemSolver {
  public:
    explicit CudaBemSolver(SurfaceMesh mesh);
    CudaBemSolver(SurfaceMesh mesh, CudaStorageOptions storage);
    CudaBemSolver(SurfaceMesh mesh, BemSubdomainModel model, CudaStorageOptions storage = {});
    CudaBemSolver(SurfaceMesh mesh, BemSubdomainModel model, OrthogonalSymmetryPlanes symmetry,
                  CudaStorageOptions storage = {});
    CudaBemSolver(SurfaceMesh mesh, BemSubdomainModel model, OrthogonalSymmetryPlanes symmetry,
                  ExteriorBaffle baffle, CudaStorageOptions storage = {});
    CudaBemSolver(SurfaceMesh mesh, BemSubdomainModel model, ExteriorBaffle baffle,
                  CudaStorageOptions storage = {});
    ~CudaBemSolver();
    CudaBemSolver(const CudaBemSolver &) = delete;
    CudaBemSolver &operator=(const CudaBemSolver &) = delete;
    BemResult solve(const BemRequest &request, std::stop_token stop = {},
                    std::function<void(const BemProgress &)> progress = {});
    static std::string deviceDescription();

  private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};
} // namespace hornsculpt
