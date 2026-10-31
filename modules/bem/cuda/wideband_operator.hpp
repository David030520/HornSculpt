#pragma once

#include "../fmm/scaled_translation.hpp"
#include "../fmm/spherical_rotation.hpp"
#include "../fmm/support_tree.hpp"
#include "../fmm/wideband_transfer.hpp"
#include "../galerkin/pair_quadrature.hpp"
#include <cstddef>
#include <cstdint>
#include <cuComplex.h>
#include <cuda_runtime_api.h>
#include <memory>
#include <span>
#include <stop_token>
#include <vector>

namespace hornsculpt::wideband {

// Batched separable angular transform: cuFFT in azimuth, normalized Legendre
// projection/evaluation in polar angle. Scratch stores only one bounded m-tile;
// no q-by-q interpolation matrix or full O(p^3) Legendre table is allocated.
class AngularTransferOperator {
  public:
    AngularTransferOperator(const AngularTransferPlan &plan, unsigned batchCount,
                            std::stop_token stop = {}, cudaStream_t stream = nullptr);
    ~AngularTransferOperator();
    AngularTransferOperator(const AngularTransferOperator &) = delete;
    AngularTransferOperator &operator=(const AngularTransferOperator &) = delete;
    std::vector<Complex> apply(std::span<const Complex> source, std::stop_token stop = {});
    // Device arrays keep the constructor's fixed dimensions and stay on CUDA.
    // Only scalar failure status is read back; calls complete before returning.
    void applyDevice(const cuDoubleComplex *source, cuDoubleComplex *output,
                     std::stop_token stop = {});
    cudaStream_t stream() const noexcept;
    std::vector<Complex> outgoingToParent(std::span<const Complex> source, double wavenumber,
                                          Vec3 childMinusParent, std::stop_token stop = {});
    std::vector<Complex> incomingToChild(std::span<const Complex> source, double wavenumber,
                                         Vec3 childMinusParent, std::stop_token stop = {});
    void outgoingToParentDevice(const cuDoubleComplex *source, cuDoubleComplex *output,
                                double wavenumber, Vec3 childMinusParent,
                                std::stop_token stop = {});
    void incomingToChildDevice(const cuDoubleComplex *source, cuDoubleComplex *output,
                               double wavenumber, Vec3 childMinusParent, std::stop_token stop = {});
    // Stable LF/HF boundary: scaled spherical coefficients have the plan's
    // bandlimit, while angular samples retain the independent source/target
    // grids. The source normalization radius belongs to the LF child.
    void scaledOutgoingToParentDevice(const cuDoubleComplex *multipole, cuDoubleComplex *spectrum,
                                      double wavenumber, double childRadius, Vec3 childMinusParent,
                                      std::stop_token stop = {});
    void incomingToScaledChildDevice(const cuDoubleComplex *spectrum, cuDoubleComplex *local,
                                     double wavenumber, double childRadius, Vec3 childMinusParent,
                                     std::stop_token stop = {});
    std::size_t deviceBytes() const noexcept;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Exact non-owning allocation/shape identities for captured graph keys.
struct ScaledAxialDeviceIdentity {
    const void *offsets{}, *coefficients{};
    unsigned sourceOrder{}, targetOrder{};
    std::size_t sourceCapacity{}, targetCapacity{};
};
struct SphericalRotationDeviceIdentity {
    const void *offsets{}, *matrix{}, *phases{};
    unsigned order{};
    std::size_t capacity{};
};
// One distance/radius/order class of the scaled LF translation. Only its
// azimuth-preserving m-blocks are stored: O(p^3), never a dense p^4 matrix.
class ScaledAxialOperator {
  public:
    ScaledAxialOperator(const ScaledAxialPlan &plan, unsigned batchCount, std::stop_token stop = {},
                        cudaStream_t stream = nullptr);
    ~ScaledAxialOperator();
    ScaledAxialOperator(const ScaledAxialOperator &) = delete;
    ScaledAxialOperator &operator=(const ScaledAxialOperator &) = delete;
    std::vector<Complex> apply(std::span<const Complex> source, std::stop_token stop = {});
    void applyDevice(const cuDoubleComplex *source, cuDoubleComplex *output,
                     std::stop_token stop = {});
    // Apply only the leading active batches within the reserved constructor capacity.
    void applyDeviceBatch(const cuDoubleComplex *source, cuDoubleComplex *output,
                          unsigned activeBatches, std::stop_token stop = {});
    // Enqueue on stream() without a host wait. Caller initializes the device
    // status to zero and reads it only after synchronizing the same stream.
    // Each stage ORs finite failures into status; no stage clears it.
    void enqueueDeviceBatch(const cuDoubleComplex *source, cuDoubleComplex *output,
                            unsigned activeBatches, int *status, std::stop_token stop = {});
    ScaledAxialDeviceIdentity deviceIdentity() const noexcept;
    cudaStream_t stream() const noexcept;
    std::size_t deviceBytes() const noexcept;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class SphericalRotationOperator {
  public:
    SphericalRotationOperator(const SphericalRotationPlan &plan, unsigned batchCount,
                              std::stop_token stop = {}, cudaStream_t stream = nullptr);
    ~SphericalRotationOperator();
    SphericalRotationOperator(const SphericalRotationOperator &) = delete;
    SphericalRotationOperator &operator=(const SphericalRotationOperator &) = delete;
    std::vector<Complex> toAxis(std::span<const Complex> source, std::stop_token stop = {});
    std::vector<Complex> fromAxis(std::span<const Complex> source, std::stop_token stop = {});
    // Source and output are distinct. Azimuth phases are fused into the
    // cooperative source tile; device applications retain no state scratch.
    void toAxisDevice(const cuDoubleComplex *source, cuDoubleComplex *output,
                      std::stop_token stop = {});
    void fromAxisDevice(const cuDoubleComplex *source, cuDoubleComplex *output,
                        std::stop_token stop = {});
    // Leading active batches only; trailing storage is neither read nor written.
    void toAxisDeviceBatch(const cuDoubleComplex *source, cuDoubleComplex *output,
                           unsigned activeBatches, std::stop_token stop = {});
    void fromAxisDeviceBatch(const cuDoubleComplex *source, cuDoubleComplex *output,
                             unsigned activeBatches, std::stop_token stop = {});
    // Complete degrees through activeOrder <= the stored plan order. Batches
    // are packed with stride (activeOrder+1)^2; the same native d/phase plan is
    // reused. Trailing reserved storage is neither read nor written.
    void toAxisPrefixDeviceBatch(const cuDoubleComplex *source, cuDoubleComplex *output,
                                 unsigned activeOrder, unsigned activeBatches,
                                 std::stop_token stop = {});
    void fromAxisPrefixDeviceBatch(const cuDoubleComplex *source, cuDoubleComplex *output,
                                   unsigned activeOrder, unsigned activeBatches,
                                   std::stop_token stop = {});
    // Same stream/status ownership contract as ScaledAxialOperator::enqueueDeviceBatch.
    void enqueueToAxisDeviceBatch(const cuDoubleComplex *source, cuDoubleComplex *output,
                                  unsigned activeBatches, int *status, std::stop_token stop = {});
    void enqueueFromAxisDeviceBatch(const cuDoubleComplex *source, cuDoubleComplex *output,
                                    unsigned activeBatches, int *status, std::stop_token stop = {});
    void enqueueToAxisPrefixDeviceBatch(const cuDoubleComplex *source, cuDoubleComplex *output,
                                        unsigned activeOrder, unsigned activeBatches, int *status,
                                        std::stop_token stop = {});
    void enqueueFromAxisPrefixDeviceBatch(const cuDoubleComplex *source, cuDoubleComplex *output,
                                          unsigned activeOrder, unsigned activeBatches, int *status,
                                          std::stop_token stop = {});
    // Enqueue the 2*order+1 azimuth phases on stream(), without a host payload
    // or wait. Rotations/graph replays using these phases must follow on that
    // same stream; a later update follows all earlier uses in stream order.
    // The matrix/phase allocation identities remain unchanged. Reject angles
    // whose largest m*phi would be non-finite before enqueuing device work.
    void setAzimuth(double radians, std::stop_token stop = {});
    SphericalRotationDeviceIdentity deviceIdentity() const noexcept;
    cudaStream_t stream() const noexcept;
    std::size_t deviceBytes() const noexcept;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// One bounded rule bank of immutable integrated LF moments, with shared
// geometry-only panel/slot topology across rule and frequency changes. Snapshot
// objects must remain immutable and outlive the bank. A bank/snapshot change
// requires released leaf leases on its shared stream; a rule change genuinely
// reintegrates the moments and retains only the unchanged topology.
class SupportLeafMomentFamily {
  public:
    explicit SupportLeafMomentFamily(std::size_t maximumBytes);
    ~SupportLeafMomentFamily();
    SupportLeafMomentFamily(const SupportLeafMomentFamily &) = delete;
    SupportLeafMomentFamily &operator=(const SupportLeafMomentFamily &) = delete;
    void prepare(const galerkin::Space &space, const SupportTree &tree,
                 std::span<const std::uint32_t> panelToActiveFaceColumn,
                 std::uint32_t activeFaceCount, double wavenumber,
                 std::span<const galerkin::TriangleQuadraturePoint> rule, cudaStream_t stream);
    std::size_t deviceBytes() const noexcept;
    std::size_t maximumBytes() const noexcept;

  private:
    friend class SupportLeafOperator;
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

struct SupportLeafMomentStatistics {
    // bytes includes topologyBytes once; the latter is a diagnostic subset.
    std::size_t bytes{}, cachedPressureDofs{}, pressureDofs{}, cachedFaceSlots{}, faceSlots{};
    std::size_t topologyBytes{}, topologyHostBytes{}, panelReferences{}, pressureIncidences{},
        fluxPanelReferences{}, constructionPanelReferences{};
    bool topologyReused{};
    double constructionMilliseconds{};
};
struct SupportLeafSettings {
    // These are leaves at one level. Pressure and face inputs retain their
    // global numbering; output states are node-major in precisely this order.
    std::vector<std::uint32_t> nodes;
    std::vector<std::uint32_t> panelToActiveFaceColumn;
    std::uint32_t activeFaceCount{};
    LevelRepresentation representation{LevelRepresentation::regularSpherical};
    unsigned order{};
    double wavenumber{}, normalizationRadius{};
    std::size_t maximumDeviceBytes{std::numeric_limits<std::size_t>::max()};
    // Optional exact integrated LF signatures, including their device indices.
    // A zero budget retains the identical canonical integrand on every apply.
    std::size_t maximumMomentCacheBytes{};
    std::size_t reservedDeviceBytes{};
    // Storage degree may exceed the active hierarchy degree. Every stored
    // degree is genuinely integrated; contractions use order and table stride
    // independently. Miller's start degree may change roundoff in a prefix;
    // canonical-integrand gates bound it rather than asserting bit identity.
    // A family retains tables across serial seed hierarchies.
    unsigned maximumMomentOrder{};
    std::shared_ptr<SupportLeafMomentFamily> momentFamily;
    SphereRule angularRule;
    std::vector<galerkin::TriangleQuadraturePoint> panelRule;
};

// Leaf emission/testing for complete P1 supports. Its LF kernels stream one
// m-column with shared radial factors. A bounded optional cache stores only
// integrated positive-m D/M signatures, never quadrature-point basis tables.
// Flux is face-constant Neumann q, while pressure uses global P1 numbering.
// Device outputs must hold nodes.size()*stateSize() complex values. Testing
// adds to the global pressure result so the exact near contribution survives.
class SupportLeafOperator {
  public:
    SupportLeafOperator(const galerkin::Space &space, const SupportTree &tree,
                        const SupportLeafSettings &settings, std::stop_token stop = {},
                        cudaStream_t stream = nullptr);
    ~SupportLeafOperator();
    SupportLeafOperator(const SupportLeafOperator &) = delete;
    SupportLeafOperator &operator=(const SupportLeafOperator &) = delete;
    void emitPressureDevice(const cuDoubleComplex *pressure, cuDoubleComplex *moments,
                            std::stop_token stop = {});
    void emitFaceFluxDevice(const cuDoubleComplex *flux, cuDoubleComplex *moments,
                            std::stop_token stop = {});
    void testDevice(const cuDoubleComplex *locals, cuDoubleComplex *pressure,
                    std::stop_token stop = {});
    std::size_t stateSize() const noexcept;
    std::size_t deviceBytes() const noexcept;
    std::size_t ownDeviceBytes() const noexcept;
    SupportLeafMomentStatistics momentStatistics() const noexcept;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace hornsculpt::wideband
