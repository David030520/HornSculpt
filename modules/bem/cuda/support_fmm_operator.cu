#include "../fmm/rcr_work.hpp"
#include "../fmm/support_level_plan.hpp"
#include "../fmm/translation_plan_cache.hpp"
#include "support_fmm_operator.hpp"
#include "support_m2l_graph.hpp"
#include "wideband_contraction.cuh"
#include "wideband_runtime.cuh"
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <compare>
#include <map>
#include <numeric>
#include <tuple>

namespace hornsculpt::wideband {
namespace {
using namespace detail;
using Direction = std::array<std::int64_t, 3>;
using Key = std::array<std::int64_t, 8>;
using Clock = std::chrono::steady_clock;
double milliseconds(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
struct IndexPair {
    unsigned source{}, target{};
};
struct LinkGroup {
    Direction direction{};
    Vec3 displacement{};
    DeviceBuffer<IndexPair> pairs;
    std::size_t count{};
    unsigned sourceOrder{}, targetOrder{};
};
__global__ void gather(const Z *input, Z *output, const IndexPair *pairs, std::size_t first,
                       unsigned count, unsigned capacity, std::size_t state, bool reverse,
                       std::size_t packed) {
    const auto flat = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (flat >= static_cast<std::size_t>(capacity) * packed)
        return;
    const auto batch = flat / packed, coefficient = flat % packed;
    Z value{};
    if (batch < count) {
        const auto link = pairs[first + batch];
        const auto node = reverse ? link.target : link.source;
        value = input[static_cast<std::size_t>(node) * state + coefficient];
    }
    output[flat] = value;
}
__global__ void scatter(const Z *input, Z *output, const IndexPair *pairs, std::size_t first,
                        unsigned count, std::size_t state, bool reverse, bool add,
                        std::size_t packed) {
    const auto flat = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (flat >= static_cast<std::size_t>(count) * packed)
        return;
    const auto batch = flat / packed, coefficient = flat % packed;
    const auto link = pairs[first + batch];
    const auto node = reverse ? link.source : link.target;
    auto *destination = output + static_cast<std::size_t>(node) * state + coefficient;
    const auto value = input[flat];
    // Within one dyadic offset class each target has a unique source cell.
    // Parent/child groups likewise fix the octant, so no atomics are needed.
    if (add) {
        *destination = make_cuDoubleComplex(destination->x + value.x, destination->y + value.y);
    } else
        *destination = value;
}
__global__ void diagonalProduct(const Z *source, const Z *diagonal, Z *target, std::size_t state,
                                std::size_t total) {
    const auto flat = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (flat >= total)
        return;
    const auto a = source[flat], b = diagonal[flat % state];
    target[flat] = make_cuDoubleComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}
Direction difference(const SupportNode &target, const SupportNode &source) {
    Direction result{};
    for (unsigned axis = 0; axis < 3; ++axis)
        result[axis] =
            static_cast<std::int64_t>(target.coordinates[axis]) - source.coordinates[axis];
    return result;
}
std::uint64_t squared(Direction direction) {
    std::uint64_t result = 0;
    for (const auto value : direction)
        result += static_cast<std::uint64_t>(value * value);
    return result;
}
Key rotationKey(unsigned order, Direction direction) {
    const auto norm = squared(direction),
               z = static_cast<std::uint64_t>(direction[2] * direction[2]),
               divisor = std::gcd(norm, z);
    if (!divisor)
        throw std::invalid_argument("Support rotation direction is zero.");
    return {0,
            order,
            direction[2] < 0 ? -1 : 1,
            static_cast<std::int64_t>(z / divisor),
            static_cast<std::int64_t>(norm / divisor),
            0,
            0,
            0};
}
std::size_t axialEntries(unsigned source, unsigned target) {
    std::size_t result = 0;
    for (unsigned m = 0; m <= std::min(source, target); ++m)
        result = sum(result, product(source - m + 1, target - m + 1));
    return result;
}
} // namespace
struct SupportFmmOperator::Impl {
    struct Level {
        WidebandLevelPlan plan;
        double policySeparation{};
        SphereRule rule;
        std::vector<unsigned> nodes;
        std::vector<LinkGroup> pressureFar, fluxFar, up;
        DeviceBuffer<Z> outgoing, local;
        DeviceBuffer<IndexPair> leafLinks;
        std::unique_ptr<SupportLeafOperator> leaf;
        std::size_t state{}, leafCount{};
    };
    // Equality covers every captured address, source/target shape and scalar.
    // Matrix allocations are immutable within a class generation. Azimuth data
    // is updated before each group in this same stream's execution order.
    struct M2lGraphKey {
        bool flux{};
        int device{};
        unsigned sourceLevel{}, targetLevel{}, sourceOrder{}, targetOrder{}, activeTargetOrder{},
            batchCapacity{}, axialSourceOrder{}, axialTargetOrder{}, sourceRepresentation{},
            targetRepresentation{};
        std::uint64_t hierarchyGeneration{}, streamGeneration{}, sourceMatrixGeneration{},
            targetMatrixGeneration{}, axialGeneration{}, wavenumberBits{}, sourceRadiusBits{},
            targetRadiusBits{}, displacementXBits{}, displacementYBits{}, displacementZBits{};
        std::uintptr_t stream{}, input{}, output{}, first{}, second{}, third{}, error{}, pairs{},
            sourceOffsets{}, sourceMatrix{}, sourcePhases{}, targetOffsets{}, targetMatrix{},
            targetPhases{}, axialOffsets{}, axialCoefficients{};
        std::size_t firstLink{}, linkCount{}, sourceState{}, targetState{},
            sourceRotationCapacity{}, targetRotationCapacity{}, axialSourceCapacity{},
            axialTargetCapacity{};
        Direction direction{};
        auto operator<=>(const M2lGraphKey &) const = default;
    };
    struct M2lGraphEntry {
        M2lGraphKey key;
        M2lChunkGraph executable;
        std::size_t nodes{};
    };
    struct CacheEntry {
        std::shared_ptr<SphericalRotationOperator> rotation;
        std::shared_ptr<ScaledAxialOperator> axial;
        std::shared_ptr<AngularTransferOperator> angular;
        std::shared_ptr<DeviceBuffer<Z>> diagonal;
        std::size_t bytes{}, touch{};
        std::uint64_t generation{};
        bool pinned() const {
            return (rotation && rotation.use_count() > 1) || (axial && axial.use_count() > 1) ||
                   (angular && angular.use_count() > 1) || (diagonal && diagonal.use_count() > 1);
        }
    };
    const galerkin::Space &space;
    std::shared_ptr<const SupportTree> tree;
    std::shared_ptr<const SupportNearPlan> nearPlan;
    std::shared_ptr<SupportNearOperator> near;
    std::shared_ptr<SupportBlockJacobi> preconditioner;
    std::shared_ptr<TranslationPlanCache> cpuPlans;
    std::shared_ptr<SupportLeafMomentFamily> momentFamily;
    std::vector<unsigned> familyOrders;
    std::size_t familyStateReserve{};
    SupportFmmOperatorSettings settings;
    cudaStream_t stream{};
    double k{};
    std::vector<Level> levels;
    std::vector<unsigned> localIndex;
    std::vector<bool> activeFlux;
    std::map<Key, CacheEntry> cache;
    DeviceBuffer<Z> first, second, third, leafState;
    DeviceBuffer<int> error;
    FastMultipoleStatistics statistics;
    SupportFmmWorkStatistics work;
    std::vector<BemFastMultipoleLevelEvidence> evidence;
    std::size_t fixedBytes{}, cacheBytes{}, peak{}, stationary{}, clock{};
    std::size_t combinedPeak{}, combinedStationary{};
    bool released{};
    std::vector<M2lGraphEntry> m2lGraphs;
    std::size_t graphExecutables{}, graphNodes{};
    std::uint64_t hierarchyGeneration{1}, matrixGeneration{};
    unsigned queuedRcrBatches{};
    std::size_t queuedRcrWork{};
    unsigned long long streamGeneration{};
    int device{};
    ~Impl() {
        // Explicit release/eviction paths report drain errors. During stack
        // unwinding, finish the borrowed stream before destroying raw references.
        if (stream)
            cudaStreamSynchronize(stream);
        m2lGraphs.clear();
    }
    Impl(const galerkin::Space &snapshot, std::shared_ptr<const SupportTree> supportTree,
         std::shared_ptr<const SupportNearPlan> plan,
         std::shared_ptr<SupportNearOperator> nearOperator, std::shared_ptr<SupportBlockJacobi> pc,
         double wavenumber, SupportFmmOperatorSettings configuration, std::stop_token stop,
         std::shared_ptr<TranslationPlanCache> sharedPlans = {},
         std::shared_ptr<SupportLeafMomentFamily> sharedMoments = {},
         std::vector<unsigned> sharedOrders = {}, std::size_t sharedStateReserve = 0)
        : space(snapshot), tree(std::move(supportTree)), nearPlan(std::move(plan)),
          near(std::move(nearOperator)), preconditioner(std::move(pc)),
          settings(std::move(configuration)), k(wavenumber) {
        cancelled(stop);
        if (!tree || !nearPlan || !near || !preconditioner || !std::isfinite(k) || k <= 0 ||
            (settings.role != SupportFmmOperatorRole::nativeAudit &&
             settings.role != SupportFmmOperatorRole::distancePrimary) ||
            !std::isfinite(settings.accuracySeed) || settings.accuracySeed <= 0 ||
            settings.accuracySeed >= 1 || !settings.classBatchNodes ||
            settings.classBatchNodes > 256 || settings.panelRule.empty())
            throw std::invalid_argument(
                "Complete support FMM settings or shared state are invalid.");
        if (space.vertices.size() != tree->pressureLeaf.size() ||
            space.panels.size() != tree->fluxLeaf.size() ||
            space.panels.size() != nearPlan->panelToActiveFaceColumn.size() ||
            space.vertices.size() != near->pressureMatrix().rows)
            throw std::invalid_argument(
                "Complete support FMM snapshots have inconsistent dimensions.");
        stream = near->stream();
        checked(cudaGetDevice(&device), "Support FMM device identity");
        checked(cudaStreamGetId(stream, &streamGeneration), "Support FMM stream identity");
        if (settings.maximumM2lGraphEntries > 16384 || settings.maximumM2lGraphExecutables > 8192 ||
            settings.maximumM2lGraphNodes > 524288 ||
            settings.maximumM2lGraphHostBytes > 8 * 1024 * 1024)
            throw std::invalid_argument("Support M2L graph caps exceed their bounded contract.");
        cpuPlans = sharedPlans ? std::move(sharedPlans)
                               : std::make_shared<TranslationPlanCache>(
                                     settings.maximumHostPlanBytes, settings.maximumPlanWorkers);
        statistics.nodeCount = tree->nodes.size();
        for (const auto &node : tree->nodes)
            statistics.leafCount += node.leaf();
        statistics.farInteractionCount = sum(tree->pressurePressure.farInteractions.size(),
                                             tree->pressureFlux.farInteractions.size());
        statistics.pressurePressureNearPairs = tree->pressurePressure.representedNearPairs;
        statistics.pressurePressureFarPairs = tree->pressurePressure.representedFarPairs;
        statistics.pressureFluxNearPairs = tree->pressureFlux.representedNearPairs;
        statistics.pressureFluxFarPairs = tree->pressureFlux.representedFarPairs;
        const auto expansionPlan = makeSupportExpansionPlan(
            *tree, nearPlan->panelToActiveFaceColumn, k, settings.accuracySeed, stop);
        momentFamily =
            sharedMoments
                ? std::move(sharedMoments)
                : std::make_shared<SupportLeafMomentFamily>(settings.maximumLeafMomentCacheBytes);
        // A single bank survives released hierarchies, but a changed rule is a
        // genuinely new integral and may replace it only with no active leases.
        momentFamily->prepare(space, *tree, nearPlan->panelToActiveFaceColumn,
                              static_cast<unsigned>(nearPlan->activeFaceIds.size()), k,
                              settings.panelRule, stream);
        familyOrders = std::move(sharedOrders);
        familyStateReserve = sharedStateReserve;
        if (familyOrders.empty()) {
            familyOrders.resize(tree->levels.size());
            auto includePlan = [&](const SupportExpansionPlan &plan) {
                std::size_t stateBytes = 0, scratch = 0, leaf = 0;
                for (unsigned l = 0; l < plan.levels.size(); ++l) {
                    const auto &p = plan.levels[l].expansion;
                    if (p.representation == LevelRepresentation::inactive)
                        continue;
                    const bool lf = p.representation == LevelRepresentation::regularSpherical;
                    if (lf) {
                        if (p.signatureOrder > maximumDeviceScaledOrder ||
                            p.electricalRadius > maximumScaledElectricalRadius)
                            throw std::runtime_error(
                                "Planned LF audit degree or radius exceeds limit128.");
                        familyOrders[l] = std::max(familyOrders[l], p.signatureOrder);
                    }
                    const auto states = lf ? p.sphericalCoefficients : p.angularSamples;
                    std::size_t nodes = 0, leaves = 0;
                    for (const auto id : tree->levels[l].nodes)
                        if (plan.activeChannels[id]) {
                            ++nodes;
                            leaves += tree->nodes[id].leaf();
                        }
                    stateBytes =
                        sum(stateBytes, product(product(2, product(nodes, states)), sizeof(Z)));
                    scratch = std::max(
                        scratch, product(states, lf ? settings.classBatchNodes
                                                    : std::min(4u, settings.classBatchNodes)));
                    leaf = std::max(leaf, product(states, leaves));
                }
                stateBytes = sum(stateBytes, sum(product(sum(product(3, scratch), leaf), sizeof(Z)),
                                                 sizeof(int)));
                familyStateReserve = std::max(familyStateReserve, stateBytes);
            };
            includePlan(expansionPlan);
            for (const double factor : {.1, .01})
                includePlan(makeSupportExpansionPlan(*tree, nearPlan->panelToActiveFaceColumn, k,
                                                     settings.accuracySeed * factor, stop));
        }
        activeFlux.resize(tree->nodes.size());
        for (std::size_t node = 0; node < tree->nodes.size(); ++node)
            activeFlux[node] = (expansionPlan.activeChannels[node] & fluxOutgoing) != 0;
        levels.resize(tree->levels.size());
        localIndex.assign(tree->nodes.size(), UINT32_MAX);
        std::size_t maximumScratch = 0, maximumLeaf = 0;
        for (unsigned l = 0; l < levels.size(); ++l) {
            auto &level = levels[l];
            level.plan = expansionPlan.levels[l].expansion;
            level.policySeparation = expansionPlan.levels[l].policySeparation;
            if (level.plan.representation == LevelRepresentation::inactive)
                continue;
            if (level.plan.representation == LevelRepresentation::regularSpherical &&
                (level.plan.signatureOrder > maximumDeviceScaledOrder ||
                 level.plan.electricalRadius > maximumScaledElectricalRadius))
                throw std::runtime_error("Support FMM LF degree or electrical radius exceeds the "
                                         "validated RCR limit128.");
            if (level.plan.representation == LevelRepresentation::planeWave)
                level.rule = makeSphereRule(level.plan.quadratureOrder);
            level.state = level.plan.representation == LevelRepresentation::regularSpherical
                              ? level.plan.sphericalCoefficients
                              : level.plan.angularSamples;
            for (const auto id : tree->levels[l].nodes) {
                const auto flags = expansionPlan.activeChannels[id];
                if ((flags & (pressureOutgoing | pressureIncoming)) || activeFlux[id]) {
                    localIndex[id] = static_cast<unsigned>(level.nodes.size());
                    level.nodes.push_back(id);
                    if (tree->nodes[id].leaf())
                        ++level.leafCount;
                }
            }
            if (level.nodes.empty())
                continue;
            evidence.push_back({l,
                                level.plan.representation == LevelRepresentation::regularSpherical
                                    ? BemFastMultipoleRepresentation::regularSpherical
                                    : BemFastMultipoleRepresentation::planeWave,
                                level.plan.signatureOrder, level.plan.translationOrder,
                                level.plan.quadratureOrder});
            const auto states = product(level.nodes.size(), level.state);
            fixedBytes = sum(fixedBytes, product(product(2, states), sizeof(Z)));
            const unsigned batches = level.plan.representation == LevelRepresentation::planeWave
                                         ? std::min(4u, settings.classBatchNodes)
                                         : settings.classBatchNodes;
            maximumScratch = std::max(maximumScratch, product(level.state, batches));
            maximumLeaf = std::max(maximumLeaf, product(level.state, level.leafCount));
        }
        fixedBytes =
            sum(fixedBytes,
                sum(product(sum(product(3, maximumScratch), maximumLeaf), sizeof(Z)), sizeof(int)));
        const auto additionalSeedStates =
            familyStateReserve > fixedBytes ? familyStateReserve - fixedBytes : 0;
        const auto sharedBytes = sharedDeviceBytes();
        if (sharedBytes > settings.maximumDeviceBytes ||
            fixedBytes > settings.maximumDeviceBytes - sharedBytes)
            throw std::runtime_error(
                "Complete support FMM states exceed the configured GPU budget.");
        std::size_t free = 0, total = 0;
        checked(cudaMemGetInfo(&free, &total), "Support FMM state memory query");
        if (fixedBytes > free - free / 5)
            throw std::runtime_error("Complete support FMM states exceed available GPU memory.");
        for (auto &level : levels)
            if (!level.nodes.empty()) {
                level.outgoing.allocate(product(level.nodes.size(), level.state));
                level.local.allocate(product(level.nodes.size(), level.state));
            }
        first.allocate(maximumScratch);
        second.allocate(maximumScratch);
        third.allocate(maximumScratch);
        leafState.allocate(maximumLeaf);
        error.allocate(1);
        for (unsigned l = 0; l < levels.size(); ++l) {
            auto &level = levels[l];
            if (level.nodes.empty())
                continue;
            prepareGroups(level.pressureFar, tree->pressurePressure.farInteractions, l, false,
                          stop);
            prepareGroups(level.fluxFar, tree->pressureFlux.farInteractions, l, true, stop);
            if (l)
                prepareUp(l, stop);
            std::vector<IndexPair> leafLinks;
            SupportLeafSettings leafSettings;
            leafSettings.activeFaceCount = static_cast<unsigned>(nearPlan->activeFaceIds.size());
            leafSettings.panelToActiveFaceColumn = nearPlan->panelToActiveFaceColumn;
            leafSettings.representation = level.plan.representation;
            leafSettings.order = level.plan.signatureOrder;
            leafSettings.wavenumber = k;
            leafSettings.normalizationRadius = level.plan.activeRadius;
            leafSettings.angularRule = level.rule;
            leafSettings.panelRule = settings.panelRule;
            leafSettings.maximumMomentCacheBytes = settings.maximumLeafMomentCacheBytes;
            leafSettings.maximumMomentOrder = std::max(level.plan.signatureOrder, familyOrders[l]);
            leafSettings.momentFamily = momentFamily;
            // Leave room for the largest of the three serial seed states, not
            // three resident hierarchies. No quadrature-point/host mirror exists.
            leafSettings.reservedDeviceBytes =
                sum(settings.maximumClassCacheBytes, additionalSeedStates);
            for (unsigned index = 0; index < level.nodes.size(); ++index)
                if (tree->nodes[level.nodes[index]].leaf()) {
                    leafLinks.push_back({static_cast<unsigned>(leafLinks.size()), index});
                    leafSettings.nodes.push_back(level.nodes[index]);
                }
            if (!leafLinks.empty()) {
                const auto shared = sharedDeviceBytes();
                const auto indexBytes = product(leafLinks.size(), sizeof(IndexPair));
                if (shared > settings.maximumDeviceBytes ||
                    fixedBytes > settings.maximumDeviceBytes - shared ||
                    indexBytes > settings.maximumDeviceBytes - shared - fixedBytes)
                    throw std::runtime_error(
                        "Support leaf indices exceed the configured GPU budget.");
                leafSettings.maximumDeviceBytes = settings.maximumDeviceBytes - shared -
                                                  fixedBytes - indexBytes +
                                                  momentFamily->deviceBytes();
                level.leafLinks.allocate(leafLinks.size());
                level.leafLinks.upload(leafLinks, stream);
                level.leaf =
                    std::make_unique<SupportLeafOperator>(space, *tree, leafSettings, stop, stream);
                fixedBytes = sum(fixedBytes, sum(product(leafLinks.size(), sizeof(IndexPair)),
                                                 level.leaf->ownDeviceBytes()));
                const auto moments = level.leaf->momentStatistics();
                work.leafMomentCacheBytes += moments.bytes;
                work.leafCachedPressureDofs += moments.cachedPressureDofs;
                work.leafCachedFaceSlots += moments.cachedFaceSlots;
                work.leafMomentConstructionMilliseconds += moments.constructionMilliseconds;
            }
        }
        work.leafMomentCacheBytes = momentFamily->deviceBytes();
        if (deviceBytes() > settings.maximumDeviceBytes)
            throw std::runtime_error(
                "Complete support FMM geometry exceeds the configured GPU budget.");
        std::vector<AxialPlanRequest> axialRequests;
        for (unsigned l = 0; l < levels.size(); ++l) {
            if (levels[l].plan.representation != LevelRepresentation::regularSpherical)
                continue;
            for (const auto *groups : {&levels[l].pressureFar, &levels[l].fluxFar})
                for (const auto &group : *groups)
                    axialRequests.push_back(
                        axialRequest(ScaledAxialKind::multipoleToLocal, l, l, group));
            if (l && levels[l - 1].plan.representation == LevelRepresentation::regularSpherical)
                for (const auto &group : levels[l].up) {
                    axialRequests.push_back(
                        axialRequest(ScaledAxialKind::multipoleToParent, l, l - 1, group));
                    axialRequests.push_back(
                        axialRequest(ScaledAxialKind::localToChild, l - 1, l, group));
                }
        }
        const auto planStart = Clock::now();
        cpuPlans->prepare(axialRequests, stop);
        work.cpuPlanPreparationMilliseconds = milliseconds(planStart);
        // Shared near/PC construction peaks belong to their earlier serial
        // preparation phases, not to this seed's resident allocation peak.
        peak = fixedBytes;
        stationary = fixedBytes;
        combinedPeak = combinedStationary = deviceBytes();
        checked(cudaStreamSynchronize(stream), "Support FMM preparation");
        cancelled(stop);
    }
    unsigned m2lEndpointOrder(unsigned l, double distance) const {
        if (l >= levels.size() || !std::isfinite(distance) || distance <= 0)
            throw std::invalid_argument("Support M2L class needs a resident level and distance.");
        const auto &level = levels[l];
        const auto &native = level.plan;
        if (settings.role == SupportFmmOperatorRole::nativeAudit ||
            native.representation != LevelRepresentation::regularSpherical)
            return native.signatureOrder;
        const double radius = native.activeRadius;
        if (!(radius > 0) || !(level.policySeparation > 1))
            throw std::logic_error("Support LF M2L class has no native convergence envelope.");
        // Both complete endpoint balls lie inside this common level envelope:
        // (d-r_target)/r_source and its converse are at least d/R-1.
        if (!(distance > 2 * radius))
            return native.signatureOrder;
        const double delta = distance / radius - 1;
        if (!(delta > level.policySeparation))
            return native.signatureOrder;
        // Eq. 18 is a starting estimate, not an operator error certificate.
        // Budget its seed between the two independently truncated endpoints;
        // uncut seed/rule audits still decide acceptance of the actual solve.
        const auto candidate = makeWidebandLevelPlan(radius, k, delta, .5 * settings.accuracySeed);
        // Retain the native LF representation and full hierarchy states. The
        // compact prefix applies only within this separated M2L class.
        return std::min(native.signatureOrder, candidate.signatureOrder);
    }
    void uploadGroups(std::vector<LinkGroup> &output,
                      std::map<Direction, std::vector<IndexPair>> &groups, unsigned level, bool up,
                      std::stop_token stop) {
        for (auto &[direction, links] : groups) {
            cancelled(stop);
            LinkGroup group;
            group.direction = direction;
            group.count = links.size();
            // Cell centers are dyadic; integer offsets avoid merging distinct
            // translations through rounded floating-point dictionary keys.
            const double step = 2 * tree->nodes[levels[level].nodes.front()].halfWidth;
            if (up)
                group.displacement = {direction[0] * step / 2, direction[1] * step / 2,
                                      direction[2] * step / 2};
            else
                group.displacement = {direction[0] * step, direction[1] * step,
                                      direction[2] * step};
            if (!up)
                group.sourceOrder = group.targetOrder =
                    m2lEndpointOrder(level, std::hypot(group.displacement.x, group.displacement.y,
                                                       group.displacement.z));
            const auto shared = sharedDeviceBytes();
            const auto indexBytes = product(links.size(), sizeof(IndexPair));
            if (shared > settings.maximumDeviceBytes ||
                fixedBytes > settings.maximumDeviceBytes - shared ||
                indexBytes > settings.maximumDeviceBytes - shared - fixedBytes)
                throw std::runtime_error(
                    "Support interaction indices exceed the configured GPU budget.");
            group.pairs.allocate(links.size());
            group.pairs.upload(links, stream);
            fixedBytes = sum(fixedBytes, product(links.size(), sizeof(IndexPair)));
            output.push_back(std::move(group));
        }
    }
    void prepareGroups(std::vector<LinkGroup> &output,
                       const std::vector<SupportInteraction> &interactions, unsigned level,
                       bool flux, std::stop_token stop) {
        std::map<Direction, std::vector<IndexPair>> groups;
        for (const auto pair : interactions)
            if (tree->nodes[pair.target].level == level && (!flux || activeFlux[pair.source])) {
                if (tree->nodes[pair.source].level != level)
                    throw std::logic_error(
                        "Diagonal support M2L requires equal-level interaction endpoints.");
                groups[difference(tree->nodes[pair.target], tree->nodes[pair.source])].push_back(
                    {localIndex[pair.source], localIndex[pair.target]});
            }
        uploadGroups(output, groups, level, false, stop);
    }
    void prepareUp(unsigned l, std::stop_token stop) {
        std::map<Direction, std::vector<IndexPair>> groups;
        for (unsigned i = 0; i < levels[l].nodes.size(); ++i) {
            const auto id = levels[l].nodes[i], parent = tree->nodes[id].parent;
            if (parent == noSupportNode || localIndex[parent] == UINT32_MAX)
                continue;
            Direction direction;
            for (unsigned axis = 0; axis < 3; ++axis)
                direction[axis] = tree->nodes[id].coordinates[axis] % 2 ? 1 : -1;
            groups[direction].push_back({i, localIndex[parent]});
        }
        uploadGroups(levels[l].up, groups, l, true, stop);
    }
    std::size_t sharedDeviceBytes() const {
        return sum(sum(near->deviceBytes(), preconditioner->deviceBytes()),
                   momentFamily->deviceBytes());
    }
    std::size_t deviceBytes() const {
        return sum(sharedDeviceBytes(), sum(fixedBytes, cacheBytes));
    }
    void completeRcrQueue(std::stop_token stop) {
        if (!queuedRcrBatches)
            return;
        completeSpectrumCheck(error.data(), stream, stop);
        ++work.rcrCompletionChecks;
        queuedRcrBatches = 0;
        queuedRcrWork = 0;
    }
    void reserveRcrQueue(unsigned batches, std::size_t operations, std::stop_token stop) {
        if (!batches || batches > RcrQueuePolicy::maximumQueuedBatches || !operations ||
            operations > RcrQueuePolicy::maximumChunkWork)
            throw std::logic_error("One support RCR chunk exceeds the bounded queue contract.");
        if (queuedRcrBatches > RcrQueuePolicy::maximumQueuedBatches - batches ||
            queuedRcrWork > RcrQueuePolicy::maximumChunkWork - operations)
            completeRcrQueue(stop);
        if (!queuedRcrBatches) {
            cancelled(stop);
            checked(cudaMemsetAsync(error.data(), 0, sizeof(int), stream),
                    "Support RCR queue failure status reset");
        }
    }
    void recordRcrQueue(unsigned batches, std::size_t operations) {
        queuedRcrBatches += batches;
        queuedRcrWork += operations;
        work.rcrBatchesQueued += batches;
        work.maximumRcrChunkBatches = std::max(work.maximumRcrChunkBatches, queuedRcrBatches);
        work.maximumRcrChunkWork = std::max(work.maximumRcrChunkWork, queuedRcrWork);
    }
    void room(std::size_t needed) {
        if (queuedRcrBatches)
            throw std::logic_error("Support class allocation requires a completed RCR queue.");
        const auto shared = sum(sharedDeviceBytes(), fixedBytes);
        const auto limit = std::min(
            settings.maximumClassCacheBytes,
            settings.maximumDeviceBytes > shared ? settings.maximumDeviceBytes - shared : 0);
        if (needed > limit)
            throw std::runtime_error(
                "One support translation class exceeds the bounded class cache.");
        while (cacheBytes > limit - needed) {
            auto oldest = cache.end();
            for (auto it = cache.begin(); it != cache.end(); ++it)
                if (!it->second.pinned() &&
                    (oldest == cache.end() || it->second.touch < oldest->second.touch))
                    oldest = it;
            if (oldest == cache.end())
                throw std::runtime_error(
                    "Simultaneously required support classes exceed the class cache.");
            invalidateM2lGraphs(oldest->second.generation);
            cacheBytes -= oldest->second.bytes;
            cache.erase(oldest);
        }
    }
    void installed(CacheEntry &entry, std::size_t bytes) {
        entry.generation = ++matrixGeneration;
        entry.bytes = bytes;
        entry.touch = ++clock;
        cacheBytes = sum(cacheBytes, bytes);
        stationary = std::max(stationary, fixedBytes + cacheBytes);
        peak = std::max(peak, fixedBytes + cacheBytes);
        combinedStationary = std::max(combinedStationary, deviceBytes());
        combinedPeak = std::max(combinedPeak, deviceBytes());
    }
    std::shared_ptr<SphericalRotationOperator> rotation(unsigned order, const LinkGroup &group,
                                                        std::stop_token stop) {
        const auto key = rotationKey(order, group.direction);
        auto it = cache.find(key);
        const bool missed = it == cache.end();
        if (missed)
            ++work.classCacheMisses[0];
        else
            ++work.classCacheHits[0];
        const auto prepared = Clock::now();
        if (it == cache.end()) {
            // In-flight canonical work also borrows matrix/phase buffers. A
            // class miss can evict them, even before any graph exists.
            completeRcrQueue(stop);
            const auto plan = makeSphericalRotationPlan(order, group.displacement);
            const auto estimate =
                sum(product(plan.smallD.size(), sizeof(double)),
                    sum(product(plan.degreeOffsets.size(), sizeof(std::size_t)),
                        sum(sum(product(plan.azimuthPhases.size(), sizeof(Z)), sizeof(int)),
                            product(rotationRowTileCount(order), sizeof(ContractionRowTile)))));
            room(estimate);
            auto object = std::make_shared<SphericalRotationOperator>(
                plan, settings.classBatchNodes, stop, stream);
            auto &entry = cache[key];
            entry.rotation = object;
            installed(entry, object->deviceBytes());
            it = cache.find(key);
        }
        if (missed)
            work.classPreparationMilliseconds += milliseconds(prepared);
        it->second.touch = ++clock;
        auto object = it->second.rotation;
        const auto azimuthStart = Clock::now();
        ++work.azimuthUploads;
        object->setAzimuth(std::atan2(group.displacement.y, group.displacement.x), stop);
        work.azimuthUploadMilliseconds += milliseconds(azimuthStart);
        return object;
    }
    AxialPlanRequest axialRequest(ScaledAxialKind kind, unsigned sourceLevel, unsigned targetLevel,
                                  const LinkGroup &group) const {
        const auto &source = levels[sourceLevel].plan, &target = levels[targetLevel].plan;
        return {
            kind,
            kind == ScaledAxialKind::multipoleToLocal ? group.sourceOrder : source.signatureOrder,
            kind == ScaledAxialKind::multipoleToLocal ? group.targetOrder : target.signatureOrder,
            2 * (source.signatureOrder + target.signatureOrder) + 24,
            k,
            source.activeRadius,
            target.activeRadius,
            std::hypot(group.displacement.x, group.displacement.y, group.displacement.z),
            0};
    }
    Key axialKey(ScaledAxialKind kind, unsigned sourceLevel, unsigned targetLevel,
                 const LinkGroup &group) const {
        const auto &source = levels[sourceLevel].plan, &target = levels[targetLevel].plan;
        const auto sourceOrder =
            kind == ScaledAxialKind::multipoleToLocal ? group.sourceOrder : source.signatureOrder;
        const auto targetOrder =
            kind == ScaledAxialKind::multipoleToLocal ? group.targetOrder : target.signatureOrder;
        return {1,
                static_cast<std::int64_t>(kind),
                sourceLevel,
                targetLevel,
                static_cast<std::int64_t>(squared(group.direction)),
                sourceOrder,
                targetOrder,
                0};
    }
    std::shared_ptr<ScaledAxialOperator> axial(ScaledAxialKind kind, unsigned sourceLevel,
                                               unsigned targetLevel, const LinkGroup &group,
                                               std::stop_token stop) {
        const auto &source = levels[sourceLevel].plan, &target = levels[targetLevel].plan;
        const auto sourceOrder =
            kind == ScaledAxialKind::multipoleToLocal ? group.sourceOrder : source.signatureOrder;
        const auto targetOrder =
            kind == ScaledAxialKind::multipoleToLocal ? group.targetOrder : target.signatureOrder;
        const auto key = axialKey(kind, sourceLevel, targetLevel, group);
        auto it = cache.find(key);
        const bool missed = it == cache.end();
        if (missed)
            ++work.classCacheMisses[1];
        else
            ++work.classCacheHits[1];
        const auto prepared = Clock::now();
        if (it == cache.end()) {
            completeRcrQueue(stop);
            const auto entries = axialEntries(sourceOrder, targetOrder);
            room(sum(sum(product(entries, sizeof(Z)),
                         product(std::min(sourceOrder, targetOrder) + 2, sizeof(std::size_t))),
                     sum(sizeof(int),
                         product(axialRowTileCount(targetOrder), sizeof(ContractionRowTile)))));
            const auto plan =
                cpuPlans->axial(axialRequest(kind, sourceLevel, targetLevel, group), stop);
            auto object = std::make_shared<ScaledAxialOperator>(*plan, settings.classBatchNodes,
                                                                stop, stream);
            auto &entry = cache[key];
            entry.axial = object;
            installed(entry, object->deviceBytes());
            it = cache.find(key);
        }
        if (missed)
            work.classPreparationMilliseconds += milliseconds(prepared);
        it->second.touch = ++clock;
        return it->second.axial;
    }
    std::shared_ptr<AngularTransferOperator> angular(unsigned sourceLevel, unsigned targetLevel,
                                                     unsigned bandlimit, bool handoff,
                                                     std::stop_token stop) {
        const Key key{2, sourceLevel, targetLevel, bandlimit, handoff, 0, 0, 0};
        auto it = cache.find(key);
        const bool missed = it == cache.end();
        if (missed)
            ++work.classCacheMisses[2];
        else
            ++work.classCacheHits[2];
        const auto prepared = Clock::now();
        if (it == cache.end()) {
            const auto &sourceRule = handoff ? levels[std::min(sourceLevel, targetLevel)].rule
                                             : levels[sourceLevel].rule;
            const auto &targetRule = handoff ? sourceRule : levels[targetLevel].rule;
            // cuFFT reports its actual workspace only when the plan is built.
            // Reserve a conservative explicit working estimate, then account
            // the actual allocation and evict other unpinned classes if needed.
            const unsigned batch = std::min(4u, settings.classBatchNodes);
            const auto estimated =
                product(product(batch, sourceRule.samples.size() + targetRule.samples.size() +
                                           product(bandlimit + 1, bandlimit + 1)),
                        sizeof(Z));
            room(sum(product(2, estimated), 4 * 1024 * 1024));
            const auto plan = makeAngularTransferPlan(sourceRule, targetRule, bandlimit);
            auto object = std::make_shared<AngularTransferOperator>(plan, batch, stop, stream);
            peak = std::max(peak, fixedBytes + cacheBytes + object->deviceBytes());
            combinedPeak = std::max(combinedPeak, sum(deviceBytes(), object->deviceBytes()));
            room(object->deviceBytes());
            auto &entry = cache[key];
            entry.angular = object;
            installed(entry, object->deviceBytes());
            it = cache.find(key);
        }
        if (missed)
            work.classPreparationMilliseconds += milliseconds(prepared);
        it->second.touch = ++clock;
        return it->second.angular;
    }
    std::shared_ptr<DeviceBuffer<Z>> diagonal(unsigned l, const LinkGroup &group,
                                              std::stop_token stop) {
        const Key key{3, l, group.direction[0], group.direction[1], group.direction[2], 0, 0, 0};
        auto it = cache.find(key);
        const bool missed = it == cache.end();
        if (missed)
            ++work.classCacheMisses[3];
        else
            ++work.classCacheHits[3];
        const auto prepared = Clock::now();
        if (it == cache.end()) {
            room(product(levels[l].state, sizeof(Z)));
            const auto coefficients = packed(diagonalTranslation(
                levels[l].rule, levels[l].plan.translationOrder, k, group.displacement));
            cancelled(stop);
            auto object = std::make_shared<DeviceBuffer<Z>>();
            object->allocate(coefficients.size());
            object->upload(coefficients, stream);
            checked(cudaStreamSynchronize(stream), "Support diagonal class preparation");
            auto &entry = cache[key];
            entry.diagonal = object;
            installed(entry, product(coefficients.size(), sizeof(Z)));
            it = cache.find(key);
        }
        if (missed)
            work.classPreparationMilliseconds += milliseconds(prepared);
        it->second.touch = ++clock;
        return it->second.diagonal;
    }
    void gathered(const Z *input, const LinkGroup &group, std::size_t start, unsigned count,
                  unsigned batches, std::size_t state, bool reverse, std::stop_token stop,
                  std::size_t prefix = 0) {
        cancelled(stop);
        const auto packed = prefix ? prefix : state;
        if (packed > state)
            throw std::logic_error("Support gather prefix exceeds the native node stride.");
        const auto total = product(batches, packed);
        gather<<<static_cast<unsigned>((total + threads - 1) / threads), threads, 0, stream>>>(
            input, first.data(), group.pairs.data(), start, count, batches, state, reverse, packed);
        checked(cudaGetLastError(), "Support state gather");
    }
    void scattered(Z *output, const LinkGroup &group, std::size_t start, unsigned count,
                   std::size_t state, bool reverse, bool add, std::stop_token stop,
                   std::size_t prefix = 0) {
        cancelled(stop);
        const auto packed = prefix ? prefix : state;
        if (packed > state)
            throw std::logic_error("Support scatter prefix exceeds the native node stride.");
        const auto total = product(count, packed);
        scatter<<<static_cast<unsigned>((total + threads - 1) / threads), threads, 0, stream>>>(
            third.data(), output, group.pairs.data(), start, count, state, reverse, add, packed);
        checked(cudaGetLastError(), "Support state scatter");
    }
    void graphEvidence() {
        work.m2lGraphCachedEntries = m2lGraphs.size();
        work.m2lGraphCachedExecutables = graphExecutables;
        work.m2lGraphCachedNodes = graphNodes;
        work.m2lGraphMaximumCachedExecutables =
            std::max(work.m2lGraphMaximumCachedExecutables, graphExecutables);
        work.m2lGraphMaximumCachedNodes = std::max(work.m2lGraphMaximumCachedNodes, graphNodes);
        work.m2lGraphOwnedHostPayloadBytes = product(m2lGraphs.capacity(), sizeof(M2lGraphEntry));
        work.m2lGraphPeakOwnedHostPayloadBytes =
            std::max(work.m2lGraphPeakOwnedHostPayloadBytes, work.m2lGraphOwnedHostPayloadBytes);
    }
    void closeGraphAdmission() {
        if (!work.m2lGraphAdmissionClosed) {
            work.m2lGraphAdmissionClosed = true;
            ++work.m2lGraphAdmissionClosures;
        }
    }
    void invalidateM2lGraphs(std::uint64_t generation) {
        const auto affected = [&](const M2lGraphEntry &entry) {
            return entry.key.sourceMatrixGeneration == generation ||
                   entry.key.targetMatrixGeneration == generation ||
                   entry.key.axialGeneration == generation;
        };
        const bool executableAffected =
            std::any_of(m2lGraphs.begin(), m2lGraphs.end(), [&](const M2lGraphEntry &entry) {
                return affected(entry) && bool(entry.executable);
            });
        if (executableAffected) {
            checked(cudaStreamSynchronize(stream), "Support M2L graph class eviction drain");
            // A bounded class cache must not induce graph reinstantiation each
            // iteration. Surviving executable hits are still usable.
            closeGraphAdmission();
        }
        for (const auto &entry : m2lGraphs)
            if (affected(entry) && entry.executable) {
                --graphExecutables;
                graphNodes -= entry.nodes;
                ++work.m2lGraphInvalidatedExecutables;
            }
        m2lGraphs.erase(std::remove_if(m2lGraphs.begin(), m2lGraphs.end(), affected),
                        m2lGraphs.end());
        graphEvidence();
    }
    void clearM2lGraphsAfterDrain() {
        work.m2lGraphInvalidatedExecutables += graphExecutables;
        std::vector<M2lGraphEntry>{}.swap(m2lGraphs);
        graphExecutables = graphNodes = 0;
        graphEvidence();
    }
    std::uint64_t classGeneration(const Key &key, const void *object) const {
        // The symbolic class key is already known. Scanning the entire cache
        // for each of three matrices would turn K graph contacts into O(K*C)
        // host pointer walks even when every device plan is already resident.
        const auto found = cache.find(key);
        if (found != cache.end() &&
            (found->second.rotation.get() == object || found->second.axial.get() == object))
            return found->second.generation;
        throw std::logic_error("M2L graph matrix has no live cache generation.");
    }
    M2lGraphKey graphKey(unsigned sourceLevel, unsigned targetLevel, const LinkGroup &group,
                         bool flux, const Z *input, Z *output, std::size_t start, std::size_t count,
                         unsigned batchCapacity, const SphericalRotationOperator &sourceRotation,
                         const SphericalRotationOperator &targetRotation,
                         const ScaledAxialOperator &translation) const {
        const auto source = sourceRotation.deviceIdentity(),
                   target = targetRotation.deviceIdentity();
        const auto axial = translation.deviceIdentity();
        if (source.order != levels[sourceLevel].plan.signatureOrder ||
            target.order != levels[targetLevel].plan.signatureOrder ||
            axial.sourceOrder != group.sourceOrder || axial.sourceOrder > source.order ||
            axial.targetOrder != group.targetOrder || axial.targetOrder > target.order)
            throw std::logic_error("Support M2L graph native/packed orders disagree.");
        const auto address = [](const void *pointer) {
            return reinterpret_cast<std::uintptr_t>(pointer);
        };
        M2lGraphKey key;
        key.flux = flux;
        key.device = device;
        key.sourceLevel = sourceLevel;
        key.targetLevel = targetLevel;
        key.sourceOrder = source.order;
        key.targetOrder = target.order;
        key.activeTargetOrder = axial.targetOrder;
        key.batchCapacity = batchCapacity;
        key.axialSourceOrder = axial.sourceOrder;
        key.axialTargetOrder = axial.targetOrder;
        key.sourceRepresentation = static_cast<unsigned>(levels[sourceLevel].plan.representation);
        key.targetRepresentation = static_cast<unsigned>(levels[targetLevel].plan.representation);
        key.hierarchyGeneration = hierarchyGeneration;
        key.streamGeneration = streamGeneration;
        key.sourceMatrixGeneration =
            classGeneration(rotationKey(source.order, group.direction), &sourceRotation);
        key.targetMatrixGeneration =
            classGeneration(rotationKey(target.order, group.direction), &targetRotation);
        key.axialGeneration = classGeneration(
            axialKey(ScaledAxialKind::multipoleToLocal, sourceLevel, targetLevel, group),
            &translation);
        key.wavenumberBits = std::bit_cast<std::uint64_t>(k);
        key.sourceRadiusBits = std::bit_cast<std::uint64_t>(levels[sourceLevel].plan.activeRadius);
        key.targetRadiusBits = std::bit_cast<std::uint64_t>(levels[targetLevel].plan.activeRadius);
        key.displacementXBits = std::bit_cast<std::uint64_t>(group.displacement.x);
        key.displacementYBits = std::bit_cast<std::uint64_t>(group.displacement.y);
        key.displacementZBits = std::bit_cast<std::uint64_t>(group.displacement.z);
        key.stream = address(stream);
        key.input = address(input);
        key.output = address(output);
        key.first = address(first.data());
        key.second = address(second.data());
        key.third = address(third.data());
        key.error = address(error.data());
        key.pairs = address(group.pairs.data());
        key.sourceOffsets = address(source.offsets);
        key.sourceMatrix = address(source.matrix);
        key.sourcePhases = address(source.phases);
        key.targetOffsets = address(target.offsets);
        key.targetMatrix = address(target.matrix);
        key.targetPhases = address(target.phases);
        key.axialOffsets = address(axial.offsets);
        key.axialCoefficients = address(axial.coefficients);
        key.firstLink = start;
        key.linkCount = count;
        key.sourceState = levels[sourceLevel].state;
        key.targetState = levels[targetLevel].state;
        key.sourceRotationCapacity = source.capacity;
        key.targetRotationCapacity = target.capacity;
        key.axialSourceCapacity = axial.sourceCapacity;
        key.axialTargetCapacity = axial.targetCapacity;
        key.direction = group.direction;
        return key;
    }
    std::size_t graphMemory(bool after) {
        std::size_t free = 0, total = 0;
        checked(cudaMemGetInfo(&free, &total), "Support M2L graph memory observation");
        ++work.m2lGraphMemoryQueries;
        (after ? work.m2lGraphLastFreeAfterBytes : work.m2lGraphLastFreeBeforeBytes) = free;
        work.m2lGraphTotalDeviceBytes = total;
        // This independent admission reserve concerns opaque driver graph
        // storage; it does not replace the class or hierarchy reservations.
        work.m2lGraphRequiredFreeReserveBytes =
            std::max(std::size_t{128} * 1024 * 1024, total / 5 + std::size_t(total % 5 != 0));
        if (work.m2lGraphMemoryQueries == 1 || free < work.m2lGraphMinimumObservedFreeBytes)
            work.m2lGraphMinimumObservedFreeBytes = free;
        return free;
    }
    template <class Enqueue> bool replayM2lGraph(const M2lGraphKey &key, Enqueue enqueue) {
        auto found = std::lower_bound(m2lGraphs.begin(), m2lGraphs.end(), key,
                                      [](const M2lGraphEntry &entry, const M2lGraphKey &wanted) {
                                          return entry.key < wanted;
                                      });
        if (found != m2lGraphs.end() && found->key == key && found->executable) {
            found->executable.launch(stream);
            ++work.m2lGraphLaunches;
            ++work.m2lGraphReuses;
            return true;
        }
        if (work.m2lGraphAdmissionClosed) {
            ++work.m2lGraphBudgetBypasses;
            return false;
        }
        if (m2lGraphs.empty() && m2lGraphs.capacity() == 0) {
            const auto capacity =
                std::min(settings.maximumM2lGraphEntries,
                         settings.maximumM2lGraphHostBytes / sizeof(M2lGraphEntry));
            if (!capacity) {
                closeGraphAdmission();
                ++work.m2lGraphBudgetBypasses;
                return false;
            }
            m2lGraphs.reserve(capacity);
            graphEvidence();
            found = m2lGraphs.end();
        }
        if (found == m2lGraphs.end() || found->key != key) {
            if (m2lGraphs.size() == m2lGraphs.capacity()) {
                closeGraphAdmission();
                ++work.m2lGraphBudgetBypasses;
                return false;
            }
            m2lGraphs.insert(found, M2lGraphEntry{key, {}, 0});
            graphEvidence();
            // Single-contact Bq and independent audit actions do not pay
            // graph capture/instantiation costs that they cannot amortize.
            return false;
        }
        if (graphExecutables == settings.maximumM2lGraphExecutables ||
            graphNodes == settings.maximumM2lGraphNodes) {
            closeGraphAdmission();
            ++work.m2lGraphBudgetBypasses;
            return false;
        }
        const auto freeBefore = graphMemory(false);
        if (freeBefore < work.m2lGraphRequiredFreeReserveBytes) {
            closeGraphAdmission();
            ++work.m2lGraphBudgetBypasses;
            return false;
        }
        M2lChunkGraph::Preparation timing;
        M2lChunkGraph prepared;
        ++work.m2lGraphCaptureAttempts;
        const bool instantiated =
            prepared.prepare(stream, settings.maximumM2lGraphNodes - graphNodes, enqueue, timing);
        work.m2lGraphCaptureMilliseconds += timing.captureMilliseconds;
        work.m2lGraphInstantiationMilliseconds += timing.instantiationMilliseconds;
        if (instantiated)
            ++work.m2lGraphInstantiations;
        const auto freeAfter = instantiated ? graphMemory(true) : 0;
        if (!instantiated || freeAfter < work.m2lGraphRequiredFreeReserveBytes) {
            // The rejected executable has never been launched, hence owns no
            // in-flight references. CUDA errors themselves never enter this path.
            closeGraphAdmission();
            ++work.m2lGraphBudgetBypasses;
            return false;
        }
        found->executable = std::move(prepared);
        found->nodes = timing.nodes;
        ++graphExecutables;
        graphNodes += timing.nodes;
        graphEvidence();
        found->executable.launch(stream);
        ++work.m2lGraphLaunches;
        return true;
    }
    void enqueueRcrBatch(const Z *input, Z *output, unsigned sourceLevel, unsigned targetLevel,
                         const LinkGroup &group, std::size_t start, unsigned count, bool down,
                         SphericalRotationOperator &sourceRotation,
                         ScaledAxialOperator &translation,
                         SphericalRotationOperator &targetRotation, std::stop_token stop) {
        const auto shape = translation.deviceIdentity();
        if (shape.sourceOrder > sourceRotation.deviceIdentity().order ||
            sourceRotation.deviceIdentity().order != levels[sourceLevel].plan.signatureOrder ||
            targetRotation.deviceIdentity().order != levels[targetLevel].plan.signatureOrder ||
            shape.targetOrder > targetRotation.deviceIdentity().order ||
            (sourceLevel != targetLevel &&
             (shape.sourceOrder != sourceRotation.deviceIdentity().order ||
              shape.targetOrder != targetRotation.deviceIdentity().order)))
            throw std::logic_error("Support RCR native source/prefix target shapes disagree.");
        gathered(input, group, start, count, count, levels[sourceLevel].state, down, stop,
                 product(shape.sourceOrder + 1, shape.sourceOrder + 1));
        sourceRotation.enqueueToAxisPrefixDeviceBatch(first.data(), second.data(),
                                                      shape.sourceOrder, count, error.data(), stop);
        translation.enqueueDeviceBatch(second.data(), first.data(), count, error.data(), stop);
        const auto targetOrder = translation.deviceIdentity().targetOrder;
        if (targetOrder < targetRotation.deviceIdentity().order)
            targetRotation.enqueueFromAxisPrefixDeviceBatch(first.data(), third.data(), targetOrder,
                                                            count, error.data(), stop);
        else
            targetRotation.enqueueFromAxisDeviceBatch(first.data(), third.data(), count,
                                                      error.data(), stop);
        scattered(output, group, start, count, levels[targetLevel].state, down, true, stop,
                  product(targetOrder + 1, targetOrder + 1));
    }
    void m2lChunks(unsigned sourceLevel, unsigned targetLevel, const LinkGroup &group, bool flux,
                   const Z *input, Z *output, const RcrQueuePolicy &policy,
                   SphericalRotationOperator &sourceRotation, ScaledAxialOperator &translation,
                   SphericalRotationOperator &targetRotation, std::stop_token stop) {
        for (std::size_t start = 0; start < group.count;) {
            cancelled(stop);
            std::size_t end = start, chunkWork = 0;
            unsigned chunkBatches = 0;
            while (end < group.count && chunkBatches < policy.maximumChunkBatches) {
                const auto count = static_cast<unsigned>(
                    std::min(std::size_t(policy.batchCapacity), group.count - end));
                const auto nextWork = product(count, policy.workPerLink);
                if (chunkBatches && chunkWork > RcrQueuePolicy::maximumChunkWork - nextWork)
                    break;
                chunkWork += nextWork;
                end += count;
                ++chunkBatches;
            }
            reserveRcrQueue(chunkBatches, chunkWork, stop);
            const auto enqueue = [&](bool recording) {
                if (!recording)
                    ++work.m2lGraphDirectChunks;
                for (auto firstLink = start; firstLink < end; firstLink += policy.batchCapacity) {
                    const auto count = static_cast<unsigned>(
                        std::min(std::size_t(policy.batchCapacity), end - firstLink));
                    enqueueRcrBatch(input, output, sourceLevel, targetLevel, group, firstLink,
                                    count, false, sourceRotation, translation, targetRotation,
                                    stop);
                    if (!recording) {
                        ++work.m2lGraphDirectBatches;
                    }
                }
            };
            const auto key =
                graphKey(sourceLevel, targetLevel, group, flux, input, output, start, end - start,
                         policy.batchCapacity, sourceRotation, targetRotation, translation);
            if (replayM2lGraph(key, [&] { enqueue(true); })) {
                work.m2lGraphExecutedBatches += chunkBatches;
            } else
                enqueue(false);
            // Graphs never reset this shared status. Finite later groups must
            // not erase an earlier failure before the bounded queue completes.
            recordRcrQueue(chunkBatches, chunkWork);
            start = end;
        }
    }
    void transfer(unsigned sourceLevel, unsigned targetLevel, const LinkGroup &group, bool down,
                  bool m2l, std::stop_token stop, bool flux = false) {
        auto &source = levels[sourceLevel], &target = levels[targetLevel];
        const bool sourceLf = source.plan.representation == LevelRepresentation::regularSpherical,
                   targetLf = target.plan.representation == LevelRepresentation::regularSpherical;
        const unsigned batches = sourceLf && targetLf ? settings.classBatchNodes
                                                      : std::min(4u, settings.classBatchNodes);
        const auto input = down ? source.local.data() : source.outgoing.data();
        auto output = down || m2l ? target.local.data() : target.outgoing.data();
        if (m2l && sourceLevel != targetLevel)
            throw std::logic_error("Support far interactions must be equal-level.");
        // A direction class is constant over its batches. Pin its plans once;
        // the azimuth update is ordered before their use on this same stream.
        std::shared_ptr<SphericalRotationOperator> sourceRotation, targetRotation;
        std::shared_ptr<ScaledAxialOperator> translation;
        std::shared_ptr<DeviceBuffer<Z>> coefficients;
        std::shared_ptr<AngularTransferOperator> conversion;
        const bool queuedRcr = sourceLf && targetLf;
        const auto policy =
            queuedRcr
                ? makeRcrQueuePolicy(m2l ? group.sourceOrder : source.plan.signatureOrder,
                                     m2l ? group.targetOrder : target.plan.signatureOrder, batches)
                : RcrQueuePolicy{batches, 1, 0};
        try {
            // Angular transfers own their completion/status protocol. Never
            // let one of them reset or consume a pending LF queue's status.
            if (!queuedRcr)
                completeRcrQueue(stop);
            if (sourceLf && targetLf) {
                sourceRotation = rotation(source.plan.signatureOrder, group, stop);
                targetRotation = source.plan.signatureOrder == target.plan.signatureOrder
                                     ? sourceRotation
                                     : rotation(target.plan.signatureOrder, group, stop);
                translation = axial(m2l    ? ScaledAxialKind::multipoleToLocal
                                    : down ? ScaledAxialKind::localToChild
                                           : ScaledAxialKind::multipoleToParent,
                                    sourceLevel, targetLevel, group, stop);
            } else if (m2l)
                coefficients = diagonal(sourceLevel, group, stop);
            else
                conversion = angular(
                    sourceLevel, targetLevel,
                    sourceLf   ? source.plan.signatureOrder
                    : targetLf ? target.plan.signatureOrder
                               : std::min(source.plan.signatureOrder, target.plan.signatureOrder),
                    sourceLf != targetLf, stop);
            if (m2l && queuedRcr) {
                m2lChunks(sourceLevel, targetLevel, group, flux, input, output, policy,
                          *sourceRotation, *translation, *targetRotation, stop);
                return;
            }
            for (std::size_t start = 0; start < group.count; start += policy.batchCapacity) {
                const unsigned count = static_cast<unsigned>(
                    std::min(static_cast<std::size_t>(policy.batchCapacity), group.count - start));
                if (queuedRcr) {
                    const auto nextWork = product(count, policy.workPerLink);
                    reserveRcrQueue(1, nextWork, stop);
                }
                if (sourceLf && targetLf) {
                    enqueueRcrBatch(input, output, sourceLevel, targetLevel, group, start, count,
                                    down, *sourceRotation, *translation, *targetRotation, stop);
                } else if (m2l) {
                    gathered(input, group, start, count, batches, source.state, down, stop);
                    const auto total = product(batches, source.state);
                    diagonalProduct<<<static_cast<unsigned>((total + threads - 1) / threads),
                                      threads, 0, stream>>>(first.data(), coefficients->data(),
                                                            third.data(), source.state, total);
                    checked(cudaGetLastError(), "Support diagonal M2L");
                } else if (sourceLf && !targetLf) {
                    gathered(input, group, start, count, batches, source.state, down, stop);
                    if (down)
                        throw std::logic_error(
                            "A scaled parent cannot have a plane-wave child in the "
                            "monotone wideband plan.");
                    conversion->scaledOutgoingToParentDevice(first.data(), third.data(), k,
                                                             source.plan.activeRadius,
                                                             group.displacement, stop);
                } else if (!sourceLf && targetLf) {
                    gathered(input, group, start, count, batches, source.state, down, stop);
                    if (!down)
                        throw std::logic_error(
                            "A plane-wave child cannot have a scaled parent in the "
                            "monotone wideband plan.");
                    conversion->incomingToScaledChildDevice(first.data(), third.data(), k,
                                                            target.plan.activeRadius,
                                                            group.displacement, stop);
                } else {
                    gathered(input, group, start, count, batches, source.state, down, stop);
                    if (down)
                        conversion->incomingToChildDevice(first.data(), third.data(), k,
                                                          group.displacement, stop);
                    else
                        conversion->outgoingToParentDevice(first.data(), third.data(), k,
                                                           group.displacement, stop);
                }
                if (!queuedRcr)
                    scattered(output, group, start, count, target.state, down, true, stop);
                if (queuedRcr) {
                    recordRcrQueue(1, product(count, policy.workPerLink));
                } else
                    checked(cudaStreamSynchronize(stream), "Support state batch completion");
            }
        } catch (...) {
            checked(cudaStreamSynchronize(stream), "Support transfer failure drain");
            throw;
        }
    }
    void apply(const Z *input, Z *output, bool flux, std::stop_token stop) {
        try {
            applyImpl(input, output, flux, stop);
        } catch (...) {
            checked(cudaStreamSynchronize(stream), "Support application failure drain");
            queuedRcrBatches = 0;
            queuedRcrWork = 0;
            throw;
        }
    }
    void applyImpl(const Z *input, Z *output, bool flux, std::stop_token stop) {
        cancelled(stop);
        int currentDevice = 0;
        unsigned long long currentStream = 0;
        checked(cudaGetDevice(&currentDevice), "Support FMM current device");
        checked(cudaStreamGetId(stream, &currentStream), "Support FMM current stream");
        if (currentDevice != device || currentStream != streamGeneration)
            throw std::runtime_error("Support FMM context/stream generation changed.");
        if (!input || !output || input == output)
            throw std::invalid_argument("Support FMM application needs distinct device arrays.");
        auto phase = Clock::now();
        if (flux)
            near->applyFaceFluxDevice(input, output, stop);
        else
            near->applyPressureDevice(input, output, stop);
        work.nearMilliseconds = milliseconds(phase);
        phase = Clock::now();
        for (auto &level : levels)
            if (!level.nodes.empty()) {
                const auto bytes = product(product(level.nodes.size(), level.state), sizeof(Z));
                checked(cudaMemsetAsync(level.outgoing.data(), 0, bytes, stream),
                        "Support outgoing initialization");
                checked(cudaMemsetAsync(level.local.data(), 0, bytes, stream),
                        "Support local initialization");
                if (level.leaf) {
                    if (flux)
                        level.leaf->emitFaceFluxDevice(input, leafState.data(), stop);
                    else
                        level.leaf->emitPressureDevice(input, leafState.data(), stop);
                    const auto total = product(level.leafCount, level.state);
                    scatter<<<static_cast<unsigned>((total + threads - 1) / threads), threads, 0,
                              stream>>>(leafState.data(), level.outgoing.data(),
                                        level.leafLinks.data(), 0,
                                        static_cast<unsigned>(level.leafCount), level.state, false,
                                        false, level.state);
                    checked(cudaGetLastError(), "Support leaf state injection");
                }
            }
        work.emissionMilliseconds = milliseconds(phase);
        phase = Clock::now();
        for (std::size_t l = levels.size(); l-- > 1;)
            for (const auto &group : levels[l].up)
                transfer(static_cast<unsigned>(l), static_cast<unsigned>(l - 1), group, false,
                         false, stop);
        completeRcrQueue(stop);
        work.upwardMilliseconds = milliseconds(phase);
        phase = Clock::now();
        for (unsigned l = 0; l < levels.size(); ++l)
            for (const auto &group : flux ? levels[l].fluxFar : levels[l].pressureFar)
                transfer(l, l, group, false, true, stop, flux);
        completeRcrQueue(stop);
        work.translationMilliseconds = milliseconds(phase);
        phase = Clock::now();
        for (unsigned l = 1; l < levels.size(); ++l)
            for (const auto &group : levels[l].up)
                transfer(l - 1, l, group, true, false, stop);
        completeRcrQueue(stop);
        work.downwardMilliseconds = milliseconds(phase);
        phase = Clock::now();
        for (auto &level : levels)
            if (level.leaf) {
                const auto total = product(level.leafCount, level.state);
                gather<<<static_cast<unsigned>((total + threads - 1) / threads), threads, 0,
                         stream>>>(level.local.data(), leafState.data(), level.leafLinks.data(), 0,
                                   static_cast<unsigned>(level.leafCount),
                                   static_cast<unsigned>(level.leafCount), level.state, true,
                                   level.state);
                checked(cudaGetLastError(), "Support leaf local extraction");
                level.leaf->testDevice(leafState.data(), output, stop);
            }
        verifySpectrum(output, space.vertices.size(), error.data(), stream, stop);
        work.testingMilliseconds = milliseconds(phase);
        ++work.applications;
        peak = std::max(peak, fixedBytes + cacheBytes);
        combinedStationary = std::max(combinedStationary, deviceBytes());
        combinedPeak = std::max(combinedPeak, deviceBytes());
    }
};
SupportFmmOperator::SupportFmmOperator(const galerkin::Space &space,
                                       std::shared_ptr<const SupportTree> tree,
                                       std::shared_ptr<const SupportNearPlan> plan,
                                       std::shared_ptr<SupportNearOperator> near,
                                       std::shared_ptr<SupportBlockJacobi> pc, double k,
                                       SupportFmmOperatorSettings settings, std::stop_token stop)
    : impl_(std::make_unique<Impl>(space, std::move(tree), std::move(plan), std::move(near),
                                   std::move(pc), k, std::move(settings), stop)) {}
SupportFmmOperator::SupportFmmOperator(SupportFmmOperator &source, double seed,
                                       std::stop_token stop)
    : SupportFmmOperator(source, seed, source.impl_->settings.panelRule,
                         SupportFmmOperatorRole::nativeAudit, stop) {}
SupportFmmOperator::SupportFmmOperator(SupportFmmOperator &source, double seed,
                                       std::vector<galerkin::TriangleQuadraturePoint> panelRule,
                                       SupportFmmOperatorRole role, std::stop_token stop) {
    auto settings = source.impl_->settings;
    settings.accuracySeed = seed;
    settings.role = role;
    settings.panelRule = std::move(panelRule);
    impl_ = std::make_unique<Impl>(source.impl_->space, source.impl_->tree, source.impl_->nearPlan,
                                   source.impl_->near, source.impl_->preconditioner,
                                   source.impl_->k, std::move(settings), stop,
                                   source.impl_->cpuPlans, source.impl_->momentFamily,
                                   source.impl_->familyOrders, source.impl_->familyStateReserve);
    impl_->hierarchyGeneration = source.impl_->hierarchyGeneration + 1;
}
SupportFmmOperator::~SupportFmmOperator() = default;
void SupportFmmOperator::ensureState(std::stop_token stop) {
    if (!impl_->released)
        return;
    auto fresh = std::make_unique<Impl>(impl_->space, impl_->tree, impl_->nearPlan, impl_->near,
                                        impl_->preconditioner, impl_->k, impl_->settings, stop,
                                        impl_->cpuPlans, impl_->momentFamily, impl_->familyOrders,
                                        impl_->familyStateReserve);
    for (unsigned kind = 0; kind < 4; ++kind) {
        fresh->work.classCacheHits[kind] += impl_->work.classCacheHits[kind];
        fresh->work.classCacheMisses[kind] += impl_->work.classCacheMisses[kind];
    }
    fresh->hierarchyGeneration = impl_->hierarchyGeneration + 1;
    fresh->work.m2lGraphCaptureAttempts += impl_->work.m2lGraphCaptureAttempts;
    fresh->work.m2lGraphInstantiations += impl_->work.m2lGraphInstantiations;
    fresh->work.m2lGraphLaunches += impl_->work.m2lGraphLaunches;
    fresh->work.m2lGraphReuses += impl_->work.m2lGraphReuses;
    fresh->work.m2lGraphExecutedBatches += impl_->work.m2lGraphExecutedBatches;
    fresh->work.m2lGraphDirectChunks += impl_->work.m2lGraphDirectChunks;
    fresh->work.m2lGraphDirectBatches += impl_->work.m2lGraphDirectBatches;
    fresh->work.m2lGraphInvalidatedExecutables += impl_->work.m2lGraphInvalidatedExecutables;
    fresh->work.m2lGraphBudgetBypasses += impl_->work.m2lGraphBudgetBypasses;
    fresh->work.m2lGraphAdmissionClosures += impl_->work.m2lGraphAdmissionClosures;
    fresh->work.m2lGraphCaptureMilliseconds += impl_->work.m2lGraphCaptureMilliseconds;
    fresh->work.m2lGraphInstantiationMilliseconds += impl_->work.m2lGraphInstantiationMilliseconds;
    fresh->work.m2lGraphMemoryQueries += impl_->work.m2lGraphMemoryQueries;
    fresh->work.m2lGraphMaximumCachedNodes =
        std::max(fresh->work.m2lGraphMaximumCachedNodes, impl_->work.m2lGraphMaximumCachedNodes);
    fresh->work.m2lGraphMaximumCachedExecutables = std::max(
        fresh->work.m2lGraphMaximumCachedExecutables, impl_->work.m2lGraphMaximumCachedExecutables);
    fresh->work.m2lGraphPeakOwnedHostPayloadBytes =
        std::max(fresh->work.m2lGraphPeakOwnedHostPayloadBytes,
                 impl_->work.m2lGraphPeakOwnedHostPayloadBytes);
    if (impl_->work.m2lGraphMemoryQueries) {
        fresh->work.m2lGraphLastFreeBeforeBytes = impl_->work.m2lGraphLastFreeBeforeBytes;
        fresh->work.m2lGraphLastFreeAfterBytes = impl_->work.m2lGraphLastFreeAfterBytes;
        fresh->work.m2lGraphMinimumObservedFreeBytes = impl_->work.m2lGraphMinimumObservedFreeBytes;
        fresh->work.m2lGraphTotalDeviceBytes = impl_->work.m2lGraphTotalDeviceBytes;
        fresh->work.m2lGraphRequiredFreeReserveBytes = impl_->work.m2lGraphRequiredFreeReserveBytes;
    }
    fresh->work.applications += impl_->work.applications;
    fresh->work.rcrBatchesQueued += impl_->work.rcrBatchesQueued;
    fresh->work.rcrCompletionChecks += impl_->work.rcrCompletionChecks;
    fresh->work.maximumRcrChunkBatches =
        std::max(fresh->work.maximumRcrChunkBatches, impl_->work.maximumRcrChunkBatches);
    fresh->work.maximumRcrChunkWork =
        std::max(fresh->work.maximumRcrChunkWork, impl_->work.maximumRcrChunkWork);
    fresh->work.azimuthUploads += impl_->work.azimuthUploads;
    fresh->work.cpuPlanPreparationMilliseconds += impl_->work.cpuPlanPreparationMilliseconds;
    fresh->work.classPreparationMilliseconds += impl_->work.classPreparationMilliseconds;
    fresh->work.leafMomentConstructionMilliseconds +=
        impl_->work.leafMomentConstructionMilliseconds;
    fresh->work.azimuthUploadMilliseconds += impl_->work.azimuthUploadMilliseconds;
    fresh->peak = std::max(fresh->peak, impl_->peak);
    fresh->stationary = std::max(fresh->stationary, impl_->stationary);
    fresh->combinedStationary = std::max(fresh->combinedStationary, impl_->combinedStationary);
    fresh->combinedPeak = std::max(fresh->combinedPeak, impl_->combinedPeak);
    impl_ = std::move(fresh);
}
void SupportFmmOperator::releaseDeviceState() {
    checked(cudaStreamSynchronize(impl_->stream), "Support hierarchy release");
    impl_->clearM2lGraphsAfterDrain();
    impl_->cache.clear();
    impl_->levels.clear();
    impl_->first = {};
    impl_->second = {};
    impl_->third = {};
    impl_->leafState = {};
    impl_->error = {};
    impl_->fixedBytes = impl_->cacheBytes = 0;
    impl_->released = true;
}
void SupportFmmOperator::applyPressureDevice(const Z *input, Z *output, std::stop_token stop) {
    ensureState(stop);
    impl_->apply(input, output, false, stop);
}
void SupportFmmOperator::applyFaceFluxDevice(const Z *input, Z *output, std::stop_token stop) {
    ensureState(stop);
    impl_->apply(input, output, true, stop);
}
void SupportFmmOperator::preconditionDevice(Z *values, std::stop_token stop) {
    impl_->preconditioner->applyDevice(values, stop);
}
cudaStream_t SupportFmmOperator::stream() const noexcept {
    return impl_->stream;
}
std::span<const std::uint32_t> SupportFmmOperator::activeFaceIds() const noexcept {
    return impl_->nearPlan->activeFaceIds;
}
const FastMultipoleStatistics &SupportFmmOperator::statistics() const noexcept {
    return impl_->statistics;
}
std::span<const BemFastMultipoleLevelEvidence> SupportFmmOperator::levelEvidence() const noexcept {
    return impl_->evidence;
}
double SupportFmmOperator::accuracySeed() const noexcept {
    return impl_->settings.accuracySeed;
}
SupportFmmOperatorRole SupportFmmOperator::role() const noexcept {
    return impl_->settings.role;
}
unsigned SupportFmmOperator::m2lTargetOrder(unsigned level, double centerDistance) const {
    return impl_->m2lEndpointOrder(level, centerDistance);
}
unsigned SupportFmmOperator::m2lSourceOrder(unsigned level, double centerDistance) const {
    return impl_->m2lEndpointOrder(level, centerDistance);
}
bool SupportFmmOperator::hasDistanceTargetGraphs() const noexcept {
    return std::any_of(impl_->m2lGraphs.begin(), impl_->m2lGraphs.end(), [](const auto &entry) {
        return entry.executable && entry.key.activeTargetOrder < entry.key.targetOrder;
    });
}
std::span<const galerkin::TriangleQuadraturePoint> SupportFmmOperator::panelRule() const noexcept {
    return impl_->settings.panelRule;
}
std::size_t SupportFmmOperator::deviceBytes() const noexcept {
    return impl_->deviceBytes();
}
std::size_t SupportFmmOperator::ownDeviceBytes() const noexcept {
    return impl_->fixedBytes + impl_->cacheBytes;
}
std::size_t SupportFmmOperator::sharedDeviceBytes() const noexcept {
    return impl_->sharedDeviceBytes();
}
std::size_t SupportFmmOperator::ownPeakDeviceBytes() const noexcept {
    return impl_->peak;
}
const SupportFmmWorkStatistics &SupportFmmOperator::workStatistics() const noexcept {
    return impl_->work;
}
TranslationPlanCacheStatistics SupportFmmOperator::hostPlanStatistics() const {
    return impl_->cpuPlans->statistics();
}
std::size_t SupportFmmOperator::ownStationaryDeviceBytes() const noexcept {
    return impl_->stationary;
}
std::size_t SupportFmmOperator::stationaryDeviceBytes() const noexcept {
    return std::max(impl_->combinedStationary, impl_->deviceBytes());
}
std::size_t SupportFmmOperator::peakDeviceBytes() const noexcept {
    return std::max(impl_->combinedPeak, impl_->deviceBytes());
}
} // namespace hornsculpt::wideband
