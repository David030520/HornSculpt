#include "bem_mixed_wideband_oracle.hpp"
#include "bem_support_operator_oracle.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cuda/support_fmm_operator.hpp>
#include <cuda_runtime_api.h>
#include <fmm/rcr_work.hpp>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace hornsculpt;
using namespace hornsculpt::wideband;
using Complex = std::complex<double>;
using DeviceComplex = cuDoubleComplex;
// The existing complete/mixed weak-form operator gates use this tolerance.
// Isolated far rows below prevent the strong exact near block hiding FMM error.
// A spatial half-block alone is not a near/far partition: complete P1 supports
// and the adaptive tree decide ownership. Subtract the actual shared GPU near
// action, then compare with an independently integrated CPU weak form masked
// by exactly those far DOF products. The full operator keeps its separate gate.
// All cache/graph controls use the same stricter native hierarchy and
// degree-eight endpoint integration; only the primary shortens M2L endpoints.
// No threshold is changed with cache state or numerical role.
constexpr double weakFormTolerance = 2e-6;
constexpr double nativeAccuracySeed = 5e-10;
constexpr unsigned endpointDegree = 8;
constexpr double wavenumber = 20;
constexpr unsigned repetitions = 3;

void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}
void check(cudaError_t status) {
    if (status != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(status));
}
class DeviceArray {
  public:
    explicit DeviceArray(std::size_t size) : size_(size) {
        check(cudaMalloc(reinterpret_cast<void **>(&values), size * sizeof(DeviceComplex)));
    }
    ~DeviceArray() {
        cudaFree(values);
    }
    DeviceArray(const DeviceArray &) = delete;
    DeviceArray &operator=(const DeviceArray &) = delete;
    void upload(std::span<const Complex> input, cudaStream_t stream) {
        require(input.size() == size_, "Queue gate device input dimension mismatch.");
        std::vector<DeviceComplex> packed;
        for (const auto value : input)
            packed.push_back(make_cuDoubleComplex(value.real(), value.imag()));
        check(cudaMemcpyAsync(values, packed.data(), packed.size() * sizeof(DeviceComplex),
                              cudaMemcpyHostToDevice, stream));
        check(cudaStreamSynchronize(stream));
    }
    std::vector<Complex> download(cudaStream_t stream) const {
        std::vector<DeviceComplex> packed(size_);
        check(cudaMemcpyAsync(packed.data(), values, packed.size() * sizeof(DeviceComplex),
                              cudaMemcpyDeviceToHost, stream));
        check(cudaStreamSynchronize(stream));
        std::vector<Complex> result;
        for (const auto value : packed)
            result.emplace_back(value.x, value.y);
        return result;
    }
    DeviceComplex *values{};

  private:
    std::size_t size_{};
};

std::size_t totalMisses(const SupportFmmWorkStatistics &statistics) {
    return std::accumulate(statistics.classCacheMisses.begin(), statistics.classCacheMisses.end(),
                           std::size_t{});
}
double relative(std::span<const Complex> actual, std::span<const Complex> expected,
                std::size_t first, std::size_t count) {
    require(actual.size() == expected.size() && first + count <= actual.size() && count > 0,
            "Queue weak-form comparison dimensions are invalid.");
    double difference{}, reference{};
    for (std::size_t row = first; row < first + count; ++row) {
        require(std::isfinite(actual[row].real()) && std::isfinite(actual[row].imag()),
                "Queue application returned a nonfinite coefficient.");
        difference += std::norm(actual[row] - expected[row]);
        reference += std::norm(expected[row]);
    }
    require(reference > 0, "Queue weak-form oracle has no nonzero comparison field.");
    return std::sqrt(difference / reference);
}
double referenceNorm(std::span<const Complex> values, std::size_t first, std::size_t count) {
    double power{};
    for (std::size_t row = first; row < first + count; ++row)
        power += std::norm(values[row]);
    return std::sqrt(power);
}
std::vector<Complex> action(std::span<const Complex> matrix, std::span<const Complex> input,
                            std::size_t rows) {
    require(matrix.size() == rows * input.size(), "CPU weak-form matrix shape mismatch.");
    std::vector<Complex> result(rows);
    for (std::size_t row = 0; row < rows; ++row)
        for (std::size_t column = 0; column < input.size(); ++column)
            result[row] += matrix[row * input.size() + column] * input[column];
    return result;
}

struct Rhs {
    std::vector<Complex> pressure;
    std::array<Complex, 2> flux;
    std::vector<Complex> farPressure, farFlux;
};
Rhs rightHandSide(unsigned repetition, const test::WidebandFarBlockFixture &reference,
                  std::span<const std::uint32_t> activeFaces) {
    Rhs result;
    result.pressure.resize(32);
    result.farPressure.resize(32);
    result.farFlux.resize(32);
    const Complex first = repetition == 1 ? Complex{} : Complex{1, 0},
                  second = repetition == 0   ? Complex{}
                           : repetition == 1 ? Complex{1, 0}
                                             : Complex{.27, -.31};
    for (unsigned node = 0; node < 32; ++node) {
        result.pressure[node] =
            first * reference.pressure[0][node] + second * reference.pressure[1][node];
        result.farPressure[node] = first * reference.expectedPressure[0][node] +
                                   second * reference.expectedPressure[1][node];
        result.farFlux[node] =
            first * reference.expectedFlux[0][node] + second * reference.expectedFlux[1][node];
    }
    require(activeFaces.size() == 2 && activeFaces[0] == 5 && activeFaces[1] == 1,
            "Queue RHS requires its actual near-plan face order 5,1.");
    // The full masked oracle integrates both independently driven face groups.
    // Obtain their actual panel fluxes in the same order as the near operator;
    // never infer an inactive face from its position in the geometry.
    for (std::size_t face = 0; face < activeFaces.size(); ++face) {
        bool found = false;
        for (std::size_t panel = 0; panel < reference.space.panels.size(); ++panel) {
            if (reference.space.panels[panel].faceId != activeFaces[face])
                continue;
            const auto value = first * reference.flux[0][panel] + second * reference.flux[1][panel];
            require(!found || value == result.flux[face],
                    "Queue far oracle does not have a constant flux within its physical face.");
            result.flux[face] = value;
            found = true;
        }
        require(found, "Queue RHS active face has no physical source panels.");
    }
    return result;
}

struct FarMasks {
    std::vector<unsigned char> pressure, flux;
};
FarMasks exactFarMasks(const SupportTree &tree) {
    constexpr unsigned dofs = 32;
    FarMasks result{std::vector<unsigned char>(dofs * dofs),
                    std::vector<unsigned char>(dofs * dofs)};
    for (const auto pair : tree.pressurePressure.farInteractions) {
        const auto &target = tree.nodes[pair.target], &source = tree.nodes[pair.source];
        for (unsigned row = 0; row < target.pressureCount; ++row)
            for (unsigned column = 0; column < source.pressureCount; ++column) {
                const auto index = dofs * tree.pressureOrder[target.firstPressure + row] +
                                   tree.pressureOrder[source.firstPressure + column];
                require(!result.pressure[index], "Pressure far products overlap in queue fixture.");
                result.pressure[index] = 1;
            }
    }
    for (const auto pair : tree.pressureFlux.farInteractions) {
        const auto &target = tree.nodes[pair.target], &source = tree.nodes[pair.source];
        for (unsigned row = 0; row < target.pressureCount; ++row)
            for (unsigned column = 0; column < source.fluxCount; ++column) {
                const auto index = dofs * tree.pressureOrder[target.firstPressure + row] +
                                   tree.fluxOrder[source.firstFlux + column];
                require(!result.flux[index], "Flux far products overlap in queue fixture.");
                result.flux[index] = 1;
            }
    }
    require(
        static_cast<std::uint64_t>(std::count(result.pressure.begin(), result.pressure.end(), 1)) ==
                tree.pressurePressure.representedFarPairs &&
            static_cast<std::uint64_t>(std::count(result.flux.begin(), result.flux.end(), 1)) ==
                tree.pressureFlux.representedFarPairs,
        "Exact far-mask cardinalities do not match the actual support partition.");
    for (unsigned channel = 0; channel < 2; ++channel) {
        const auto &mask = channel ? result.flux : result.pressure;
        unsigned far{};
        for (unsigned row = 16; row < 32; ++row)
            for (unsigned source = 0; source < 16; ++source)
                far += mask[dofs * row + source];
        std::cout << "Actual opposite-half " << (channel ? "P1/DP0" : "P1/P1")
                  << " near/far DOF products=" << 256 - far << '/' << far << '\n';
        require(far > 0, "Isolated queue oracle has no actual far contribution.");
    }
    return result;
}

void checkGeometryScope(const SupportTree &tree) {
    require(tree.pressurePressure.representedFarPairs > 0 &&
                tree.pressureFlux.representedFarPairs > 0,
            "Queue gate must execute real pressure and flux far interactions.");
    std::set<unsigned> activeLevels;
    std::set<std::array<std::int64_t, 3>> directions;
    for (const auto pair : tree.pressurePressure.farInteractions) {
        const auto &source = tree.nodes[pair.source], &target = tree.nodes[pair.target];
        activeLevels.insert(source.level);
        directions.insert({std::int64_t(target.coordinates[0]) - source.coordinates[0],
                           std::int64_t(target.coordinates[1]) - source.coordinates[1],
                           std::int64_t(target.coordinates[2]) - source.coordinates[2]});
    }
    require(activeLevels.size() > 1 && directions.size() > 4,
            "Queue fixture lost its multiple active levels or direction groups.");
    std::cout << "Queue fixture actual far levels/directions=" << activeLevels.size() << '/'
              << directions.size()
              << ", P1/P1 far pairs=" << tree.pressurePressure.representedFarPairs
              << ", P1/DP0 far pairs=" << tree.pressureFlux.representedFarPairs << '\n';
}

struct CaseResult {
    std::size_t classPayloadBytes{};
    SupportFmmWorkStatistics work;
};
struct DirectedPressureClass {
    std::size_t links{}, drivenCheckedLinks{};
    unsigned sourceOrder{}, targetOrder{}, nativeOrder{};
};
CaseResult checkCase(const galerkin::Space &space, std::shared_ptr<const SupportTree> tree,
                     std::shared_ptr<const SupportNearPlan> plan,
                     std::shared_ptr<SupportNearOperator> near,
                     std::shared_ptr<SupportBlockJacobi> preconditioner,
                     const test::support_operator_oracle::Matrices &reference,
                     const test::WidebandFarBlockFixture &farReference, std::size_t classBudget,
                     bool graphs, bool evict, SupportFmmOperatorRole role) {
    require(
        farReference.sourcePressure.size() == 32 && farReference.sourceFlux.size() == 32 &&
            farReference.targetPressure.size() == 32,
        "Queue multi-link witness requires an independently integrated full 32-row far oracle.");
    SupportFmmOperatorSettings settings;
    settings.role = role;
    settings.accuracySeed = nativeAccuracySeed;
    settings.panelRule = galerkin::triangleRule(endpointDegree);
    settings.maximumClassCacheBytes = classBudget;
    if (!graphs)
        settings.maximumM2lGraphEntries = 0;
    SupportFmmOperator operation(space, tree, plan, near, preconditioner, wavenumber, settings);
    require(operation.activeFaceIds().size() == 2 && operation.activeFaceIds()[0] == 5 &&
                operation.activeFaceIds()[1] == 1,
            "Queue gate changed the independent oracle's explicit face-column order.");
    require(std::all_of(operation.levelEvidence().begin(), operation.levelEvidence().end(),
                        [](const auto &level) {
                            return level.representation ==
                                   BemFastMultipoleRepresentation::regularSpherical;
                        }),
            "The bounded RCR gate must use LF translations on every active level.");
    require(operation.role() == role && operation.accuracySeed() == nativeAccuracySeed,
            "Queue control did not retain its fixed native accuracy contract.");
    for (const auto &level : operation.levelEvidence())
        std::cout << "Queue role=" << static_cast<unsigned>(role) << ", level=" << level.level
                  << ", moment/translation order=" << level.momentBandlimit << '/'
                  << level.translationOrder << ", accuracy seed=" << nativeAccuracySeed
                  << ", endpoint degree=" << endpointDegree << '\n';
    bool shortenedEndpoints = false;
    const auto pressureWitness = rightHandSide(0, farReference, plan->activeFaceIds).pressure;
    require(std::all_of(pressureWitness.begin(), pressureWitness.end(),
                        [](Complex value) { return std::norm(value) > 0; }),
            "Queue multi-link witness must actually drive all 32 source pressure DOFs.");
    // Level and signed integer offset identify the actual directed M2L group.
    // Two driven links in a shortened group exercise compact batch strides.
    // Every target DOF is now independently checked by the full masked far
    // oracle, in addition to the unchanged separate rows-16..31 far gate.
    std::map<std::array<std::int64_t, 4>, DirectedPressureClass> pressureClasses;
    for (const auto *partition : {&tree->pressurePressure, &tree->pressureFlux})
        for (const auto pair : partition->farInteractions) {
            const auto &source = tree->nodes[pair.source], &target = tree->nodes[pair.target];
            const auto level =
                std::find_if(operation.levelEvidence().begin(), operation.levelEvidence().end(),
                             [&](const auto &entry) { return entry.level == source.level; });
            const double distance =
                std::hypot(target.center.x - source.center.x, target.center.y - source.center.y,
                           target.center.z - source.center.z);
            require(level != operation.levelEvidence().end(),
                    "Queue far pair has no native level evidence.");
            const auto sourceOrder = operation.m2lSourceOrder(source.level, distance),
                       targetOrder = operation.m2lTargetOrder(source.level, distance);
            require(sourceOrder == targetOrder && sourceOrder <= level->momentBandlimit,
                    "Queue endpoint orders lost their common conservative envelope.");
            if (role == SupportFmmOperatorRole::nativeAudit)
                require(sourceOrder == level->momentBandlimit,
                        "Queue native control shortened a far translation endpoint order.");
            shortenedEndpoints |= sourceOrder < level->momentBandlimit;
            if (partition == &tree->pressurePressure) {
                const std::array<std::int64_t, 4> key{
                    static_cast<std::int64_t>(source.level),
                    std::int64_t(target.coordinates[0]) - source.coordinates[0],
                    std::int64_t(target.coordinates[1]) - source.coordinates[1],
                    std::int64_t(target.coordinates[2]) - source.coordinates[2]};
                auto [entry, inserted] = pressureClasses.try_emplace(key);
                auto &group = entry->second;
                if (inserted) {
                    group.sourceOrder = sourceOrder;
                    group.targetOrder = targetOrder;
                    group.nativeOrder = level->momentBandlimit;
                } else
                    require(group.sourceOrder == sourceOrder && group.targetOrder == targetOrder &&
                                group.nativeOrder == level->momentBandlimit,
                            "A directed queue pressure class has inconsistent endpoint orders.");
                ++group.links;
                bool drivenSource = false, checkedTarget = false;
                for (unsigned local = 0; local < source.pressureCount; ++local) {
                    const auto dof = tree->pressureOrder[source.firstPressure + local];
                    drivenSource |= std::norm(pressureWitness[dof]) > 0;
                }
                for (unsigned local = 0; local < target.pressureCount; ++local) {
                    const auto dof = tree->pressureOrder[target.firstPressure + local];
                    checkedTarget |= dof < 32;
                }
                group.drivenCheckedLinks += drivenSource && checkedTarget;
            }
        }
    require(shortenedEndpoints == (role == SupportFmmOperatorRole::distancePrimary),
            "Queue primary/native control did not execute its actual endpoint policy.");
    std::size_t shortenedPressureClasses{}, witnessedClasses{}, maximumDrivenCheckedLinks{};
    for (const auto &[key, group] : pressureClasses) {
        if (group.sourceOrder == group.nativeOrder)
            continue;
        ++shortenedPressureClasses;
        maximumDrivenCheckedLinks = std::max(maximumDrivenCheckedLinks, group.drivenCheckedLinks);
        if (group.drivenCheckedLinks < 2)
            continue;
        ++witnessedClasses;
        std::cout << "Queue compact-source witness level=" << key[0] << ", offset=" << key[1] << ','
                  << key[2] << ',' << key[3] << ", links/driven-to-checked-links=" << group.links
                  << '/' << group.drivenCheckedLinks
                  << ", source/target/native order=" << group.sourceOrder << '/'
                  << group.targetOrder << '/' << group.nativeOrder << '\n';
    }
    std::cout << "Queue directed pressure classes/shortened/witnessed=" << pressureClasses.size()
              << '/' << shortenedPressureClasses << '/' << witnessedClasses
              << ", maximum shortened driven-to-checked links=" << maximumDrivenCheckedLinks
              << '\n';
    require((role == SupportFmmOperatorRole::distancePrimary && witnessedClasses > 0) ||
                (role == SupportFmmOperatorRole::nativeAudit && shortenedPressureClasses == 0),
            "Queue primary lacks a shortened pressure class with two driven source links into "
            "independently checked far rows 0..31, or its native control is shortened.");
    const auto initialOwnBytes = operation.ownDeviceBytes();
    DeviceArray pressure(space.vertices.size()), flux(2), output(space.vertices.size()),
        nearOutput(space.vertices.size());
    std::size_t warmMisses{}, secondMisses{};
    SupportFmmWorkStatistics beforeFinal;
    for (unsigned repetition = 0; repetition < repetitions; ++repetition) {
        const auto rhs = rightHandSide(repetition, farReference, plan->activeFaceIds);
        pressure.upload(rhs.pressure, operation.stream());
        flux.upload(rhs.flux, operation.stream());
        if (repetition == repetitions - 1)
            beforeFinal = operation.workStatistics();
        for (unsigned channel = 0; channel < 2; ++channel) {
            if (channel)
                operation.applyFaceFluxDevice(flux.values, output.values);
            else
                operation.applyPressureDevice(pressure.values, output.values);
            require(cudaStreamQuery(operation.stream()) == cudaSuccess,
                    "Completed queue operator returned while its stream still had pending work.");
            const auto actual = output.download(operation.stream());
            const auto expected = channel
                                      ? action(reference.b, rhs.flux, space.vertices.size())
                                      : action(reference.a, rhs.pressure, space.vertices.size());
            if (channel)
                near->applyFaceFluxDevice(flux.values, nearOutput.values);
            else
                near->applyPressureDevice(pressure.values, nearOutput.values);
            const auto exactNear = nearOutput.download(operation.stream());
            const auto &expectedFar = channel ? rhs.farFlux : rhs.farPressure;
            std::vector<Complex> actualFar(actual.size()), nearDifference(actual.size());
            for (std::size_t row = 0; row < actual.size(); ++row) {
                actualFar[row] = actual[row] - exactNear[row];
                nearDifference[row] = expected[row] - expectedFar[row] - exactNear[row];
            }
            const double fullError = relative(actual, expected, 0, actual.size()),
                         farError = relative(actualFar, expectedFar, 16, 16),
                         allFarError = relative(actualFar, expectedFar, 0, actual.size()),
                         nearDifferenceNorm = referenceNorm(nearDifference, 16, 16);
            std::cout << "Queue cache=" << classBudget << ", graphs=" << graphs
                      << ", eviction=" << evict << ", repetition=" << repetition
                      << ", channel=" << (channel ? "Bq" : "Ap")
                      << ", full/exact-far16..31/exact-far0..31 weak-form relative=" << fullError
                      << '/' << farError << '/' << allFarError << ", far/full reference norm="
                      << referenceNorm(expectedFar, 16, 16) /
                             referenceNorm(expected, 0, expected.size())
                      << ", opposite-half CPU-minus-far/GPU-near delta over CPU-full norm="
                      << nearDifferenceNorm / referenceNorm(expected, 16, 16) << '\n';
            require(fullError < weakFormTolerance && farError < weakFormTolerance &&
                        allFarError < weakFormTolerance,
                    "Queued complete/exact-far operator fails independent CPU weak-form action.");
        }
        if (!repetition)
            warmMisses = totalMisses(operation.workStatistics());
        else if (repetition == 1)
            secondMisses = totalMisses(operation.workStatistics());
    }
    const auto observed = operation.workStatistics();
    require(observed.applications == repetitions * 2 && observed.rcrBatchesQueued > 0 &&
                observed.rcrCompletionChecks > 0 &&
                observed.rcrCompletionChecks <= observed.rcrBatchesQueued &&
                observed.maximumRcrChunkBatches <= 32 &&
                observed.maximumRcrChunkWork <= RcrQueuePolicy::maximumChunkWork,
            "Queue work evidence exceeds its bounded application/batch/work contract.");
    const auto warmBatches = observed.rcrBatchesQueued - beforeFinal.rcrBatchesQueued,
               warmChecks = observed.rcrCompletionChecks - beforeFinal.rcrCompletionChecks;
    std::cout << "Queue cache=" << classBudget << ", graphs=" << graphs
              << ", misses warm/second/final=" << warmMisses << '/' << secondMisses << '/'
              << totalMisses(observed) << ", final-pair batches/checks=" << warmBatches << '/'
              << warmChecks << ", graph capture/instantiate/launch/reuse/invalidations="
              << observed.m2lGraphCaptureAttempts << '/' << observed.m2lGraphInstantiations << '/'
              << observed.m2lGraphLaunches << '/' << observed.m2lGraphReuses << '/'
              << observed.m2lGraphInvalidatedExecutables << '\n';
    if (evict)
        require(secondMisses > warmMisses && totalMisses(observed) > secondMisses,
                "Restricted cache did not reconstruct genuinely evicted classes on both reuses.");
    else {
        require(secondMisses == warmMisses && totalMisses(observed) == warmMisses,
                "Ample-cache control did not retain its already prepared translation classes.");
        require(warmBatches > 0 && 2 * warmChecks < warmBatches,
                "Warm multi-direction RCR queue did not at least halve status completion checks.");
        if (graphs)
            require(observed.m2lGraphLaunches > 0 && observed.m2lGraphReuses > 0 &&
                        observed.m2lGraphInvalidatedExecutables == 0,
                    "Ample-cache graph control did not instantiate and reuse actual executables.");
    }
    if (!graphs)
        require(observed.m2lGraphCaptureAttempts == 0 && observed.m2lGraphLaunches == 0 &&
                    observed.m2lGraphDirectBatches > 0,
                "Graph-budget-zero control did not execute the canonical RCR path.");
    require(operation.ownDeviceBytes() >= initialOwnBytes,
            "Queue class payload accounting decreased without releasing state.");
    const auto classPayload = operation.ownDeviceBytes() - initialOwnBytes;
    require(classPayload <= classBudget, "Queue class cache exceeded its explicit payload budget.");
    // A subsequent zero application must overwrite all prior pressure/flux
    // accumulation, including the last queued group and its shared scratch.
    pressure.upload(std::vector<Complex>(space.vertices.size()), operation.stream());
    operation.applyPressureDevice(pressure.values, output.values);
    for (const auto value : output.download(operation.stream()))
        require(value == Complex{}, "Queue reuse retained an earlier complex RHS contribution.");
    operation.releaseDeviceState();
    require(cudaStreamQuery(operation.stream()) == cudaSuccess && operation.ownDeviceBytes() == 0,
            "Queue hierarchy release retained work, captured references or device buffers.");
    return {classPayload, observed};
}

void checkQueueOperator() {
    auto fixture = test::makeMixedWidebandFixture(wavenumber, 4);
    SupportTreeSettings treeSettings;
    treeSettings.leafCapacity = 2;
    treeSettings.nearDistanceRatio = .5;
    auto tree = std::make_shared<SupportTree>(makeSupportTree(fixture.space, treeSettings));
    checkGeometryScope(*tree);
    const auto masks = exactFarMasks(*tree);
    const std::array<std::uint32_t, 2> faces{5, 1};
    auto plan = std::make_shared<SupportNearPlan>(
        makeSupportNearPlan(fixture.space, *tree, faces, 16 * 1024 * 1024));
    const auto coarseFar =
        test::makeMixedWidebandFixture(wavenumber, 12, masks.pressure, masks.flux);
    const auto fineFar = test::makeMixedWidebandFixture(wavenumber, 16, masks.pressure, masks.flux);
    for (unsigned repetition = 0; repetition < repetitions; ++repetition) {
        const auto coarse = rightHandSide(repetition, coarseFar, plan->activeFaceIds),
                   fine = rightHandSide(repetition, fineFar, plan->activeFaceIds);
        const double pressureError = relative(coarse.farPressure, fine.farPressure, 16, 16),
                     fluxError = relative(coarse.farFlux, fine.farFlux, 16, 16),
                     allPressureError = relative(coarse.farPressure, fine.farPressure, 0, 32),
                     allFluxError = relative(coarse.farFlux, fine.farFlux, 0, 32);
        std::cout << "Exact-mask independent Duffy12/16 far quadrature, RHS=" << repetition
                  << ", rows16..31 Ap/Bq relative=" << pressureError << '/' << fluxError
                  << ", rows0..31 Ap/Bq relative=" << allPressureError << '/' << allFluxError
                  << '\n';
        require(pressureError < 2e-9 && fluxError < 2e-9 && allPressureError < 2e-9 &&
                    allFluxError < 2e-9,
                "Independent exact-mask far oracle did not converge in Duffy quadrature.");
    }
    SupportNearOperatorSettings nearSettings;
    nearSettings.quadrature.nearDistanceRatio = .5;
    const auto reference = test::support_operator_oracle::directReference(fixture.space, wavenumber,
                                                                          nearSettings.quadrature);
    auto near = std::make_shared<SupportNearOperator>(fixture.space, *tree, *plan, wavenumber,
                                                      nearSettings);
    auto preconditioner = std::make_shared<SupportBlockJacobi>(
        near->pressureMatrix(), tree->pressureOrder, 32, std::stop_token{}, near->stream());
    constexpr std::size_t ampleBudget = 64 * 1024 * 1024;
    // Size only the cache-pressure experiment from observed owned class bytes.
    // Accuracy, fixture, quadrature and all comparison gates stay fixed. A
    // second/third application must prove actual reconstruction, not infer it
    // from this requested budget. Graph admission can legitimately close under
    // eviction pressure; canonical numerical work must remain correct.
    for (const auto role :
         {SupportFmmOperatorRole::nativeAudit, SupportFmmOperatorRole::distancePrimary}) {
        const auto ample = checkCase(fixture.space, tree, plan, near, preconditioner, reference,
                                     fineFar, ampleBudget, false, false, role);
        checkCase(fixture.space, tree, plan, near, preconditioner, reference, fineFar, ampleBudget,
                  true, false, role);
        require(ample.classPayloadBytes > 4, "Ample queue control recorded no class payload.");
        const std::size_t evictionBudget = ample.classPayloadBytes - ample.classPayloadBytes / 4;
        std::cout << "Measured translation-class payload=" << ample.classPayloadBytes
                  << ", forced eviction budget=" << evictionBudget << " bytes\n";
        checkCase(fixture.space, tree, plan, near, preconditioner, reference, fineFar,
                  evictionBudget, false, true, role);
        checkCase(fixture.space, tree, plan, near, preconditioner, reference, fineFar,
                  evictionBudget, true, true, role);
    }
}
} // namespace

int main() {
    int devices{};
    if (cudaGetDeviceCount(&devices) != cudaSuccess || devices == 0) {
        std::cout << "CUDA unavailable; no RCR queue operator calculation executed.\n";
        return 77;
    }
    try {
        std::cout << std::setprecision(9) << CudaBemSolver::deviceDescription() << '\n';
        checkQueueOperator();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "RCR queue operator gate failed: " << error.what() << '\n';
        return 1;
    }
}
