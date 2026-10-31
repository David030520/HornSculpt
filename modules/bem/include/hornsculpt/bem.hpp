#pragma once

#include <complex>
#include <cstddef>
#include <functional>
#include <hornsculpt/lem.hpp>
#include <hornsculpt/mesh.hpp>
#include <memory>
#include <optional>
#include <stop_token>
#include <string>
#include <vector>
#include <utility>

namespace hornsculpt {

// Independent pressure-feedback consistency budget, not a residual-derived
// error estimate. Shared by the solver and persisted-result validation.
inline constexpr double relativePortPressureClosureTolerance = 1e-3;

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

struct FastMultipoleSettings {
    // Tensor Chebyshev order, not an error tolerance. Frequency-aware separation
    // and direct original-operator residual checks guard the approximation.
    unsigned interpolationOrder{3};
    unsigned leafPanelCapacity{16};
    double separationRatio{1.5};
    unsigned maximumIterations{300};
    unsigned restartLength{40};
    bool operator==(const FastMultipoleSettings &) const = default;
};

struct FastMultipoleStatistics {
    std::size_t nodeCount{}, leafCount{}, farInteractionCount{};
    std::size_t directPanelPairCount{}, representedFarPanelPairs{};
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
    std::vector<double> frequencies;
    std::vector<BemExcitation> excitations;
    std::vector<Vec3> observers;
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

struct BemExcitationResult {
    std::string name;
    std::vector<std::complex<double>> observerPressure;
    std::vector<FaceImpedance> impedances;
    unsigned refinementSteps{};
    // True relative L2 residual of the original, unscaled Galerkin equation in
    // the selected space (full 3D or the explicit m=0 reduced space).
    double boundaryResidual{};
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
    const std::shared_ptr<const BemBoundaryData> &boundary() const noexcept { return boundary_; }

  private:
    std::shared_ptr<const BemBoundaryData> boundary_;

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
};

struct BemProgress {
    std::size_t completedFrequencies{};
    std::size_t totalFrequencies{};
    double frequency{};
};

// Positive-time RMS convention: exp(+i omega t), outgoing exp(-ikr), q=-i rho omega v.
// Continuous P1 pressure / DP0 flux Burton-Miller Galerkin: pivoted CUDA LU or
// hierarchical CUDA Fast Multipole with restarted GMRES and direct residual audit.
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
