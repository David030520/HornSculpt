#include "scaled_regular.cuh"
#include "wideband_operator.hpp"
#include "wideband_runtime.cuh"
#include <chrono>
#include <exception>
#include <map>
#include <set>

namespace hornsculpt::wideband {
namespace {
using namespace detail;
namespace regular = detail::scaled;
constexpr unsigned absent = UINT32_MAX;
struct LeafPanel {
    double3 vertices[3], normal;
    double area{};
    unsigned pressure[3]{absent, absent, absent}, face{absent};
};
struct LeafTest {
    unsigned node{}, dof{}, local{};
    LeafPanel panel;
};
__device__ double3 position(const LeafPanel &panel, const galerkin::TriangleQuadraturePoint &q) {
    return make_double3(
        panel.vertices[0].x * q.barycentric[0] + panel.vertices[1].x * q.barycentric[1] +
            panel.vertices[2].x * q.barycentric[2],
        panel.vertices[0].y * q.barycentric[0] + panel.vertices[1].y * q.barycentric[1] +
            panel.vertices[2].y * q.barycentric[2],
        panel.vertices[0].z * q.barycentric[0] + panel.vertices[1].z * q.barycentric[1] +
            panel.vertices[2].z * q.barycentric[2]);
}
__device__ double dot(double3 a, double3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
__device__ Z add(Z a, Z b) {
    return regular::add(a, b);
}
__device__ Z times(Z a, double b) {
    return regular::times(a, b);
}
__device__ Z multiply(Z a, Z b) {
    return regular::multiply(a, b);
}
__device__ Z conjugate(Z a) {
    return make_cuDoubleComplex(a.x, -a.y);
}
__device__ Z derivative(const regular::ScalarWave &wave, double3 normal) {
    return add(add(times(wave.gradient[0], normal.x), times(wave.gradient[1], normal.y)),
               times(wave.gradient[2], normal.z));
}
template <bool pressure>
__device__ Z density(const LeafPanel &panel, const galerkin::TriangleQuadraturePoint &q,
                     const Z *values) {
    Z result{};
    if constexpr (pressure) {
        for (unsigned j = 0; j < 3; ++j)
            if (panel.pressure[j] != absent)
                result = add(result, times(values[panel.pressure[j]], q.barycentric[j]));
    } else if (panel.face != absent)
        result = values[panel.face];
    return result;
}
template <bool pressure>
__global__ void emitSpherical(const LeafPanel *panels, const std::size_t *offsets,
                              const galerkin::TriangleQuadraturePoint *rule, unsigned ruleCount,
                              unsigned order, double k, double radius, const Z *input, Z *output,
                              unsigned first, int *error) {
    extern __shared__ double radial[];
    const unsigned node = first + blockIdx.x;
    const auto count = static_cast<std::size_t>(order + 1) * (order + 1), base = node * count;
    for (std::size_t j = threadIdx.x; j < count; j += blockDim.x)
        output[base + j] = {};
    __syncthreads();
    for (auto p = offsets[node]; p < offsets[node + 1]; ++p) {
        const auto panel = panels[p];
        if constexpr (pressure) {
            if (panel.pressure[0] == absent && panel.pressure[1] == absent &&
                panel.pressure[2] == absent)
                continue;
        } else if (panel.face == absent)
            continue;
        for (unsigned q = 0; q < ruleCount; ++q) {
            const auto quadrature = rule[q];
            const auto point = position(panel, quadrature);
            if (!threadIdx.x &&
                !regular::fillScaledRegularRadial(order, k * sqrt(dot(point, point)), radial))
                atomicExch(error, 1);
            __syncthreads();
            const auto weighted = times(density<pressure>(panel, quadrature, input),
                                        2 * panel.area * quadrature.weight);
            for (unsigned m = threadIdx.x; m <= order; m += blockDim.x) {
                regular::ScalarWave previous{}, current = regular::solidDiagonal(m, point, radius);
                for (unsigned l = m; l <= order; ++l) {
                    if (l > m) {
                        const auto next =
                            regular::solidNext(l, m, point, radius, previous, current);
                        previous = current;
                        current = next;
                    }
                    const auto wave =
                        regular::combineSolidRadial(current, l, point, k, radial[l], radial[l + 1]);
                    Z value;
                    if constexpr (pressure)
                        value = derivative(wave, panel.normal);
                    else
                        value = wave.value;
                    const auto positive = base + static_cast<std::size_t>(l) * l + l + m;
                    output[positive] = add(output[positive], multiply(weighted, conjugate(value)));
                    if (m) {
                        const double sign = m % 2 ? -1 : 1;
                        const auto negative = base + static_cast<std::size_t>(l) * l + l - m;
                        output[negative] =
                            add(output[negative], times(multiply(weighted, value), sign));
                    }
                }
            }
            // The shared radial table cannot be overwritten while another
            // m-column is still reading the previous quadrature point.
            __syncthreads();
        }
    }
}
template <bool pressure>
__global__ void emitAngular(const LeafPanel *panels, const std::size_t *offsets,
                            const galerkin::TriangleQuadraturePoint *rule, unsigned ruleCount,
                            const SphereSample *angles, unsigned angleCount, double k,
                            const Z *input, Z *output, std::size_t first, std::size_t count) {
    const auto within = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (within >= count)
        return;
    const auto flat = first + within, node = flat / angleCount, angle = flat % angleCount;
    const auto direction = angles[angle].direction;
    const auto d = make_double3(direction.x, direction.y, direction.z);
    Z result{};
    for (auto p = offsets[node]; p < offsets[node + 1]; ++p) {
        const auto panel = panels[p];
        const double normal = pressure ? dot(panel.normal, d) : 1;
        for (unsigned q = 0; q < ruleCount; ++q) {
            const auto quadrature = rule[q];
            double sine, cosine;
            sincos(k * dot(d, position(panel, quadrature)), &sine, &cosine);
            result = add(result, times(multiply(density<pressure>(panel, quadrature, input),
                                                make_cuDoubleComplex(cosine, sine)),
                                       2 * panel.area * quadrature.weight * normal));
        }
    }
    if constexpr (pressure)
        result = make_cuDoubleComplex(-k * result.y, k * result.x);
    output[flat] = result;
}
__global__ void testSpherical(const LeafTest *incidences, const std::size_t *offsets,
                              const galerkin::TriangleQuadraturePoint *rule, unsigned ruleCount,
                              unsigned order, double k, double radius, const Z *input, Z *output,
                              unsigned first, int *error) {
    extern __shared__ double radial[];
    __shared__ Z partial[threads];
    const unsigned target = first + blockIdx.x;
    const double inverseK = 1 / k;
    Z total{};
    const auto count = static_cast<std::size_t>(order + 1) * (order + 1);
    for (auto incidence = offsets[target]; incidence < offsets[target + 1]; ++incidence) {
        const auto entry = incidences[incidence];
        for (unsigned q = 0; q < ruleCount; ++q) {
            const auto quadrature = rule[q];
            const auto point = position(entry.panel, quadrature);
            if (!threadIdx.x &&
                !regular::fillScaledRegularRadial(order, k * sqrt(dot(point, point)), radial))
                atomicExch(error, 1);
            __syncthreads();
            for (unsigned m = threadIdx.x; m <= order; m += blockDim.x) {
                regular::ScalarWave previous{}, current = regular::solidDiagonal(m, point, radius);
                for (unsigned l = m; l <= order; ++l) {
                    if (l > m) {
                        const auto next =
                            regular::solidNext(l, m, point, radius, previous, current);
                        previous = current;
                        current = next;
                    }
                    auto wave =
                        regular::combineSolidRadial(current, l, point, k, radial[l], radial[l + 1]);
                    auto normal = derivative(wave, entry.panel.normal);
                    auto tested =
                        add(times(wave.value, -1),
                            make_cuDoubleComplex(-normal.y * inverseK, normal.x * inverseK));
                    const auto index = entry.node * count + static_cast<std::size_t>(l) * l + l;
                    Z value = multiply(tested, input[index + m]);
                    if (m) {
                        wave = regular::negativeOrder(wave, m);
                        normal = derivative(wave, entry.panel.normal);
                        tested =
                            add(times(wave.value, -1),
                                make_cuDoubleComplex(-normal.y * inverseK, normal.x * inverseK));
                        value = add(value, multiply(tested, input[index - m]));
                    }
                    total = add(total, times(value, 2 * entry.panel.area * quadrature.weight *
                                                        quadrature.barycentric[entry.local]));
                }
            }
            __syncthreads();
        }
    }
    partial[threadIdx.x] = total;
    __syncthreads();
    for (unsigned stride = blockDim.x / 2; stride; stride /= 2) {
        if (threadIdx.x < stride)
            partial[threadIdx.x] = add(partial[threadIdx.x], partial[threadIdx.x + stride]);
        __syncthreads();
    }
    if (!threadIdx.x) {
        const auto dof = incidences[offsets[target]].dof;
        output[dof] = add(output[dof], partial[0]);
    }
}
__global__ void testAngular(const LeafTest *incidences, const std::size_t *offsets,
                            const galerkin::TriangleQuadraturePoint *rule, unsigned ruleCount,
                            const SphereSample *angles, unsigned angleCount, double k,
                            const Z *input, Z *output, unsigned first) {
    __shared__ Z partial[threads];
    const unsigned target = first + blockIdx.x;
    Z total{};
    for (unsigned angle = threadIdx.x; angle < angleCount; angle += blockDim.x) {
        const auto angular = angles[angle];
        const auto d = make_double3(angular.direction.x, angular.direction.y, angular.direction.z);
        Z tested{};
        for (auto incidence = offsets[target]; incidence < offsets[target + 1]; ++incidence) {
            const auto entry = incidences[incidence];
            for (unsigned q = 0; q < ruleCount; ++q) {
                const auto quadrature = rule[q];
                double sine, cosine;
                sincos(-k * dot(d, position(entry.panel, quadrature)), &sine, &cosine);
                const double weight = 2 * entry.panel.area * quadrature.weight *
                                      quadrature.barycentric[entry.local] *
                                      (dot(entry.panel.normal, d) - 1) * angular.weight;
                tested = add(tested, make_cuDoubleComplex(weight * cosine, weight * sine));
            }
        }
        const auto node = incidences[offsets[target]].node;
        total = add(total,
                    multiply(tested, input[static_cast<std::size_t>(node) * angleCount + angle]));
    }
    partial[threadIdx.x] = total;
    __syncthreads();
    for (unsigned stride = blockDim.x / 2; stride; stride /= 2) {
        if (threadIdx.x < stride)
            partial[threadIdx.x] = add(partial[threadIdx.x], partial[threadIdx.x + stride]);
        __syncthreads();
    }
    if (!threadIdx.x) {
        const auto dof = incidences[offsets[target]].dof;
        output[dof] = add(output[dof], partial[0]);
    }
}
#include "wideband_leaf_moments.cuh"

LeafPanel localPanel(const galerkin::Panel &panel, Vec3 center) {
    LeafPanel result;
    for (unsigned j = 0; j < 3; ++j)
        result.vertices[j] =
            make_double3(panel.vertices[j].x - center.x, panel.vertices[j].y - center.y,
                         panel.vertices[j].z - center.z);
    result.normal = make_double3(panel.normal.x, panel.normal.y, panel.normal.z);
    result.area = panel.area;
    return result;
}
struct LeafMomentTopology {
    std::vector<LeafPanelMomentSlots> hostSlots;
    DeviceBuffer<LeafPanelMomentSlots> slots;
    std::size_t bytes{}, hostBytes{}, pressureIncidences{}, fluxPanels{};
};
struct LeafMomentData {
    DeviceBuffer<Z> normalMoments, massMoments, fluxMoments;
    DeviceBuffer<std::size_t> nodePressureOffsets, nodeFluxOffsets;
    DeviceBuffer<LeafFluxMomentSlot> fluxSlots;
    SupportLeafMomentStatistics statistics;
    std::shared_ptr<LeafMomentTopology> topology;
    unsigned order{}, cachedPressure{}, cachedFlux{};
    double radius{};
};
struct LeafPreparationDrain {
    cudaStream_t stream;
    int exceptions{std::uncaught_exceptions()};
    ~LeafPreparationDrain() {
        if (std::uncaught_exceptions() > exceptions)
            cudaStreamSynchronize(stream);
    }
};
} // namespace
struct SupportLeafMomentFamily::Impl {
    std::size_t maximum{}, bytes{};
    const galerkin::Space *space{};
    const SupportTree *tree{};
    std::vector<std::uint32_t> faceColumns;
    std::vector<galerkin::TriangleQuadraturePoint> rule;
    double k{};
    std::uint32_t activeFaceCount{};
    cudaStream_t stream{};
    std::map<std::vector<std::uint32_t>, std::shared_ptr<LeafMomentData>> entries;
    std::map<std::vector<std::uint32_t>, std::shared_ptr<LeafMomentTopology>> topologies;
    void idle() const {
        for (const auto &[key, entry] : entries)
            if (entry.use_count() != 1)
                throw std::runtime_error(
                    "Cannot replace a leaf moment rule bank with active leases.");
    }
};
SupportLeafMomentFamily::SupportLeafMomentFamily(std::size_t maximumBytes)
    : impl_(std::make_unique<Impl>()) {
    impl_->maximum = maximumBytes;
}
SupportLeafMomentFamily::~SupportLeafMomentFamily() = default;
void SupportLeafMomentFamily::prepare(const galerkin::Space &space, const SupportTree &tree,
                                      std::span<const std::uint32_t> faces,
                                      std::uint32_t activeFaceCount, double k,
                                      std::span<const galerkin::TriangleQuadraturePoint> rule,
                                      cudaStream_t stream) {
    auto &bank = *impl_;
    const bool sameRule = bank.rule.size() == rule.size() &&
                          std::equal(bank.rule.begin(), bank.rule.end(), rule.begin(),
                                     [](const auto &a, const auto &b) {
                                         return a.barycentric[0] == b.barycentric[0] &&
                                                a.barycentric[1] == b.barycentric[1] &&
                                                a.barycentric[2] == b.barycentric[2] &&
                                                a.weight == b.weight;
                                     });
    const bool sameGeometry =
        bank.space == &space && bank.tree == &tree && bank.stream == stream &&
        bank.activeFaceCount == activeFaceCount && bank.faceColumns.size() == faces.size() &&
        std::equal(bank.faceColumns.begin(), bank.faceColumns.end(), faces.begin());
    if (sameGeometry && bank.k == k && sameRule)
        return;
    bank.idle();
    if (!bank.entries.empty() || !bank.topologies.empty())
        checked(cudaStreamSynchronize(bank.stream), "Leaf moment bank replacement");
    bank.entries.clear();
    if (!sameGeometry)
        bank.topologies.clear();
    bank.bytes = 0;
    for (const auto &[key, topology] : bank.topologies)
        bank.bytes = sum(bank.bytes, topology->bytes);
    bank.space = &space;
    bank.tree = &tree;
    bank.k = k;
    bank.activeFaceCount = activeFaceCount;
    bank.stream = stream;
    bank.faceColumns.assign(faces.begin(), faces.end());
    bank.rule.assign(rule.begin(), rule.end());
}
std::size_t SupportLeafMomentFamily::deviceBytes() const noexcept {
    return impl_->bytes;
}
std::size_t SupportLeafMomentFamily::maximumBytes() const noexcept {
    return impl_->maximum;
}
struct SupportLeafOperator::Impl {
    Stream stream;
    DeviceBuffer<LeafPanel> panels;
    DeviceBuffer<LeafTest> incidences;
    DeviceBuffer<std::size_t> panelOffsets, testOffsets;
    DeviceBuffer<galerkin::TriangleQuadraturePoint> rule;
    DeviceBuffer<SphereSample> angles;
    DeviceBuffer<int> error;
    std::shared_ptr<SupportLeafMomentFamily> family;
    std::shared_ptr<LeafMomentData> moments;
    SupportLeafMomentStatistics momentStats;
    unsigned storageOrder{};
    unsigned cachedPressure{}, cachedFlux{};
    unsigned nodes{}, targets{}, order{}, ruleCount{}, angleCount{};
    std::size_t states{}, bytes{};
    double k{}, radius{};
    bool spherical{};
    Impl(const galerkin::Space &space, const SupportTree &tree, const SupportLeafSettings &settings,
         std::stop_token stop, cudaStream_t borrowed)
        : stream(borrowed) {
        cancelled(stop);
        spherical = settings.representation == LevelRepresentation::regularSpherical;
        if (settings.nodes.empty() || settings.nodes.size() >= absent ||
            settings.panelRule.empty() || settings.panelRule.size() >= absent ||
            settings.panelToActiveFaceColumn.size() != space.panels.size() ||
            !std::isfinite(settings.wavenumber) || settings.wavenumber <= 0 ||
            !std::isfinite(settings.normalizationRadius) || settings.normalizationRadius <= 0 ||
            (spherical && (settings.order > maximumDeviceScaledOrder ||
                           settings.wavenumber * settings.normalizationRadius >
                               maximumScaledElectricalRadius)) ||
            (!spherical && settings.representation != LevelRepresentation::planeWave))
            throw std::invalid_argument("Support leaf expansion settings are invalid.");
        for (const auto face : settings.panelToActiveFaceColumn)
            if (face != absent && face >= settings.activeFaceCount)
                throw std::invalid_argument("Support leaf face column exceeds its input basis.");
        if (!spherical &&
            (settings.angularRule.samples.empty() || settings.angularRule.samples.size() >= absent))
            throw std::invalid_argument("Support leaf angular rule is empty or too large.");
        double measure = 0;
        for (const auto &q : settings.panelRule) {
            double sum = 0;
            for (const auto b : q.barycentric) {
                if (!std::isfinite(b) || b < 0 || b > 1)
                    throw std::invalid_argument("Support leaf quadrature point is invalid.");
                sum += b;
            }
            if (std::abs(sum - 1) > 1e-12 || !std::isfinite(q.weight))
                throw std::invalid_argument("Support leaf quadrature weight is invalid.");
            measure += q.weight;
        }
        if (std::abs(measure - .5) > 1e-12)
            throw std::invalid_argument("Support leaf triangle rule has invalid total measure.");
        for (const auto &a : settings.angularRule.samples) {
            if (!std::isfinite(a.weight) || !std::isfinite(a.direction.x) ||
                !std::isfinite(a.direction.y) || !std::isfinite(a.direction.z) ||
                std::abs(std::hypot(a.direction.x, a.direction.y, a.direction.z) - 1) > 1e-12)
                throw std::invalid_argument("Support leaf angular sample is invalid.");
        }
        family = settings.momentFamily;
        if (family)
            family->prepare(space, tree, settings.panelToActiveFaceColumn, settings.activeFaceCount,
                            settings.wavenumber, settings.panelRule, stream);
        storageOrder = std::max(settings.order, settings.maximumMomentOrder);
        if (spherical && storageOrder > maximumDeviceScaledOrder)
            throw std::invalid_argument(
                "Integrated leaf storage degree exceeds the validated limit128.");
        order = settings.order;
        k = settings.wavenumber;
        radius = settings.normalizationRadius;
        nodes = static_cast<unsigned>(settings.nodes.size());
        ruleCount = static_cast<unsigned>(settings.panelRule.size());
        angleCount = static_cast<unsigned>(settings.angularRule.samples.size());
        states = spherical ? product(order + 1, order + 1) : angleCount;
        std::set<unsigned> uniqueNodes;
        std::vector<LeafPanel> source;
        std::vector<LeafTest> tests;
        std::vector<std::size_t> sourceOffsets{0}, targetOffsets{0}, nodeTargets{0};
        for (unsigned local = 0; local < nodes; ++local) {
            cancelled(stop);
            const auto id = settings.nodes[local];
            if (id >= tree.nodes.size() || !tree.nodes[id].leaf() || !uniqueNodes.insert(id).second)
                throw std::invalid_argument(
                    "Support leaf node is invalid, internal or duplicated.");
            const auto &node = tree.nodes[id];
            std::set<unsigned> panelIds;
            for (unsigned j = 0; j < node.pressureCount; ++j) {
                const auto dof = tree.pressureOrder[node.firstPressure + j];
                for (auto t = space.incidenceOffsets[dof]; t < space.incidenceOffsets[dof + 1];
                     ++t) {
                    const auto incidence = space.incidences[t];
                    panelIds.insert(incidence.panel);
                    auto panel = localPanel(space.panels[incidence.panel], node.center);
                    tests.push_back({local, dof, incidence.local, panel});
                }
                targetOffsets.push_back(tests.size());
            }
            for (unsigned j = 0; j < node.fluxCount; ++j) {
                const auto panel = tree.fluxOrder[node.firstFlux + j];
                if (settings.panelToActiveFaceColumn[panel] != absent)
                    panelIds.insert(panel);
            }
            for (const auto panelId : panelIds) {
                const auto &physical = space.panels[panelId];
                auto panel = localPanel(physical, node.center);
                for (unsigned j = 0; j < 3; ++j)
                    if (tree.pressureLeaf[physical.nodes[j]] == id)
                        panel.pressure[j] = physical.nodes[j];
                if (tree.fluxLeaf[panelId] == id)
                    panel.face = settings.panelToActiveFaceColumn[panelId];
                source.push_back(panel);
            }
            sourceOffsets.push_back(source.size());
            nodeTargets.push_back(targetOffsets.size() - 1);
        }
        if (source.size() >= absent || tests.size() >= absent || targetOffsets.size() >= absent)
            throw std::overflow_error("Support leaf target index overflow.");
        targets = static_cast<unsigned>(targetOffsets.size() - 1);
        bytes = sum(
            sum(product(source.size(), sizeof(LeafPanel)), product(tests.size(), sizeof(LeafTest))),
            sum(product(sourceOffsets.size() + targetOffsets.size(), sizeof(std::size_t)),
                sum(product(settings.panelRule.size(), sizeof(galerkin::TriangleQuadraturePoint)),
                    sum(product(settings.angularRule.samples.size(), sizeof(SphereSample)),
                        sizeof(int)))));
        std::size_t free = 0, total = 0;
        checked(cudaMemGetInfo(&free, &total), "Support leaf memory query");
        const auto residentFamily = family ? family->deviceBytes() : 0;
        if (residentFamily > settings.maximumDeviceBytes ||
            bytes > settings.maximumDeviceBytes - residentFamily || bytes > free - free / 5)
            throw std::runtime_error("Support leaf geometry exceeds available GPU memory budget.");
        // Plan exact integrated signatures before any cache allocation. Only
        // small geometry/index vectors live on the host, never moment tables.
        std::vector<LeafFluxMomentSlot> hostFluxSlots;
        std::vector<std::size_t> hostFluxNodeOffsets{0}, hostFluxPanelOffsets{0};
        std::vector<unsigned> hostFluxPanelIds;
        std::shared_ptr<LeafMomentTopology> topology;
        if (spherical && family) {
            const auto found = family->impl_->topologies.find(settings.nodes);
            if (found != family->impl_->topologies.end())
                topology = found->second;
        }
        const auto residentTopology = topology ? topology->bytes : 0;
        bool reused = false, uncachedLease = false;
        if (spherical && family && settings.maximumMomentCacheBytes) {
            auto &bank = *family->impl_;
            const auto found = bank.entries.find(settings.nodes);
            if (found != bank.entries.end()) {
                const auto &entry = found->second;
                if (entry->radius == radius && entry->order >= storageOrder &&
                    entry->statistics.bytes > settings.maximumMomentCacheBytes) {
                    // A retained family is shared state, not permission to
                    // exceed this leaf's smaller optional cache cap. Preserve
                    // other leases and use the same canonical integrand here.
                    uncachedLease = true;
                } else if (entry->radius == radius && entry->order >= storageOrder) {
                    moments = entry;
                    reused = true;
                    storageOrder = entry->order;
                    cachedPressure = entry->cachedPressure;
                    cachedFlux = entry->cachedFlux;
                    momentStats = entry->statistics;
                    momentStats.topologyReused = true;
                    momentStats.constructionPanelReferences = 0;
                    momentStats.constructionMilliseconds = 0;
                } else {
                    if (entry.use_count() != 1)
                        throw std::runtime_error(
                            "Cannot grow leased leaf moment degrees or change their radius.");
                    checked(cudaStreamSynchronize(stream), "Leaf moment degree replacement");
                    bank.bytes -= entry->statistics.bytes - entry->topology->bytes;
                    bank.entries.erase(found);
                    checked(cudaMemGetInfo(&free, &total), "Leaf degree replacement memory query");
                }
            }
        }
        if (!moments)
            moments = std::make_shared<LeafMomentData>();
        if (spherical && settings.maximumMomentCacheBytes && !uncachedLease) {
            for (unsigned node = 0; node < nodes; ++node) {
                std::map<unsigned, std::vector<unsigned>> groups;
                for (auto j = sourceOffsets[node]; j < sourceOffsets[node + 1]; ++j)
                    if (source[j].face != absent)
                        groups[source[j].face].push_back(static_cast<unsigned>(j));
                for (const auto &[face, ids] : groups) {
                    hostFluxSlots.push_back({node, face});
                    hostFluxPanelIds.insert(hostFluxPanelIds.end(), ids.begin(), ids.end());
                    hostFluxPanelOffsets.push_back(hostFluxPanelIds.size());
                }
                hostFluxNodeOffsets.push_back(hostFluxSlots.size());
            }
            if (hostFluxSlots.size() >= absent)
                throw std::overflow_error("Integrated leaf face-slot index overflow.");
            const auto positive = product(storageOrder + 1, storageOrder + 2) / 2;
            const auto pressureBytes = product(product(2, positive), sizeof(Z)),
                       fluxBytes = product(positive, sizeof(Z));
            const auto metadataBound =
                sum(product(nodeTargets.size() + hostFluxNodeOffsets.size(), sizeof(std::size_t)),
                    sum(product(hostFluxSlots.size(), sizeof(LeafFluxMomentSlot)),
                        product(source.size(), sizeof(LeafPanelMomentSlots))));
            const auto resident = family ? family->deviceBytes() : 0;
            const auto configured =
                settings.maximumDeviceBytes > resident ? settings.maximumDeviceBytes - resident : 0;
            const auto allowed = std::min(configured, free - free / 5);
            const auto reserve = sum(bytes, settings.reservedDeviceBytes);
            const auto remaining = family ? family->maximumBytes() - family->deviceBytes()
                                          : settings.maximumMomentCacheBytes;
            // Only this leaf's already resident topology is a credit. Every
            // other leaf/bank byte and the full upload/allocation peak remain
            // charged before admitting a new optional cache payload.
            const auto budget =
                std::min(settings.maximumMomentCacheBytes,
                         sum(std::min(remaining, allowed > reserve ? allowed - reserve : 0),
                             residentTopology));
            if (!reused && budget > metadataBound) {
                const auto available = budget - metadataBound;
                // Leave a small share for repeated face-constant DP0 sources;
                // pressure D/M removes both P2M and L2P quadrature work.
                cachedFlux = static_cast<unsigned>(
                    std::min(hostFluxSlots.size(), available / 8 / fluxBytes));
                cachedPressure = static_cast<unsigned>(
                    std::min(static_cast<std::size_t>(targets),
                             (available - product(cachedFlux, fluxBytes)) / pressureBytes));
                cachedFlux = static_cast<unsigned>(
                    std::min(hostFluxSlots.size(),
                             (available - product(cachedPressure, pressureBytes)) / fluxBytes));
            }
            momentStats.pressureDofs = targets;
            momentStats.faceSlots = hostFluxSlots.size();
            if (cachedPressure || cachedFlux) {
                if (!topology) {
                    topology = std::make_shared<LeafMomentTopology>();
                    topology->hostSlots.resize(source.size());
                    std::map<unsigned, unsigned> pressureTargets;
                    for (unsigned j = 0; j < targets; ++j)
                        pressureTargets.emplace(tests[targetOffsets[j]].dof, j);
                    for (std::size_t p = 0; p < source.size(); ++p)
                        for (unsigned local = 0; local < 3; ++local)
                            if (source[p].pressure[local] != absent) {
                                topology->hostSlots[p].pressure[local] =
                                    pressureTargets.at(source[p].pressure[local]);
                                ++topology->pressureIncidences;
                            }
                    for (unsigned slot = 0; slot < hostFluxSlots.size(); ++slot)
                        for (auto j = hostFluxPanelOffsets[slot];
                             j < hostFluxPanelOffsets[slot + 1]; ++j) {
                            topology->hostSlots[hostFluxPanelIds[j]].flux = slot;
                            ++topology->fluxPanels;
                        }
                    topology->bytes = product(source.size(), sizeof(LeafPanelMomentSlots));
                    topology->hostBytes =
                        sum(sizeof(LeafMomentTopology),
                            product(topology->hostSlots.capacity(), sizeof(LeafPanelMomentSlots)));
                }
                if (topology->hostSlots.size() != source.size())
                    throw std::logic_error("Leaf moment topology and source order disagree.");
                momentStats.topologyBytes = topology->bytes;
                momentStats.topologyHostBytes = topology->hostBytes;
                momentStats.topologyReused = residentTopology != 0;
                momentStats.panelReferences = topology->hostSlots.size();
                momentStats.pressureIncidences = topology->pressureIncidences;
                momentStats.fluxPanelReferences = topology->fluxPanels;
                if (!reused)
                    for (const auto &entry : topology->hostSlots)
                        momentStats.constructionPanelReferences +=
                            entry.pressure[0] < cachedPressure ||
                            entry.pressure[1] < cachedPressure ||
                            entry.pressure[2] < cachedPressure || entry.flux < cachedFlux;
                hostFluxPanelIds.resize(hostFluxPanelOffsets[cachedFlux]);
                hostFluxPanelOffsets.resize(cachedFlux + 1);
                hostFluxSlots.resize(cachedFlux);
                momentStats.bytes =
                    sum(product(cachedPressure, pressureBytes), product(cachedFlux, fluxBytes));
                momentStats.bytes =
                    sum(momentStats.bytes,
                        sum(product(nodeTargets.size() + hostFluxNodeOffsets.size(),
                                    sizeof(std::size_t)),
                            sum(product(hostFluxSlots.size(), sizeof(LeafFluxMomentSlot)),
                                topology->bytes)));
                if (!reused && momentStats.bytes > budget)
                    throw std::runtime_error(
                        "Integrated leaf moment cache exceeds its byte preflight.");
                if (!family)
                    bytes = sum(bytes, momentStats.bytes);
                std::set<unsigned> cachedDofs;
                for (unsigned j = 0; j < cachedPressure; ++j)
                    cachedDofs.insert(tests[targetOffsets[j]].dof);
                for (auto &panel : source)
                    for (auto &dof : panel.pressure)
                        if (cachedDofs.contains(dof))
                            dof = absent;
                for (unsigned slot = 0; slot < cachedFlux; ++slot)
                    for (auto j = hostFluxPanelOffsets[slot]; j < hostFluxPanelOffsets[slot + 1];
                         ++j)
                        source[hostFluxPanelIds[j]].face = absent;
            }
        }
        // Uploads and kernels share this stream; on cancellation/failure drain
        // before temporary host index vectors or unpublished device tables die.
        LeafPreparationDrain drain{stream};
        panels.allocate(source.size());
        incidences.allocate(tests.size());
        panelOffsets.allocate(sourceOffsets.size());
        testOffsets.allocate(targetOffsets.size());
        rule.allocate(settings.panelRule.size());
        angles.allocate(settings.angularRule.samples.size());
        error.allocate(1);
        panels.upload(source, stream);
        incidences.upload(tests, stream);
        panelOffsets.upload(sourceOffsets, stream);
        testOffsets.upload(targetOffsets, stream);
        rule.upload(settings.panelRule, stream);
        angles.upload(settings.angularRule.samples, stream);
        checked(cudaStreamSynchronize(stream), "Support leaf preparation");
        if (!reused && (cachedPressure || cachedFlux)) {
            const auto start = std::chrono::steady_clock::now();
            const unsigned momentThreads = std::min(threads, ((storageOrder + 32u) / 32u) * 32u);
            const auto sharedBytes = product(storageOrder + 2, sizeof(double));
            int device = 0, multiprocessors = 0, blocksPerMultiprocessor = 0;
            checked(cudaGetDevice(&device), "Leaf moment scheduling device");
            checked(
                cudaDeviceGetAttribute(&multiprocessors, cudaDevAttrMultiProcessorCount, device),
                "Leaf moment scheduling multiprocessor count");
            checked(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                        &blocksPerMultiprocessor, buildLeafMoments, static_cast<int>(momentThreads),
                        sharedBytes),
                    "Leaf moment scheduling occupancy");
            if (multiprocessors <= 0 || blocksPerMultiprocessor <= 0)
                throw std::runtime_error("Leaf moment scheduling requires positive occupancy.");
            const auto waveCapacity = static_cast<unsigned>(std::min(
                std::size_t{512}, product(static_cast<std::size_t>(multiprocessors),
                                          static_cast<std::size_t>(blocksPerMultiprocessor))));
            // Occupancy bounds independent leaves, not kernel duration. Short
            // panel slices retain explicit cancellation boundaries at every order.
            constexpr std::size_t panelsPerMomentSlice = 4;
            const auto positive = product(storageOrder + 1, storageOrder + 2) / 2;
            moments->topology = topology;
            if (!residentTopology) {
                topology->slots.allocate(topology->hostSlots.size());
                topology->slots.upload(topology->hostSlots, stream);
            }
            moments->normalMoments.allocate(product(cachedPressure, positive));
            moments->massMoments.allocate(product(cachedPressure, positive));
            moments->fluxMoments.allocate(product(cachedFlux, positive));
            moments->nodePressureOffsets.allocate(nodeTargets.size());
            moments->nodeFluxOffsets.allocate(hostFluxNodeOffsets.size());
            moments->fluxSlots.allocate(hostFluxSlots.size());
            moments->nodePressureOffsets.upload(nodeTargets, stream);
            moments->nodeFluxOffsets.upload(hostFluxNodeOffsets, stream);
            moments->fluxSlots.upload(hostFluxSlots, stream);
            checked(cudaMemsetAsync(error.data(), 0, sizeof(int), stream),
                    "Leaf moment status reset");
            if (cachedPressure) {
                const auto size = product(product(cachedPressure, positive), sizeof(Z));
                checked(cudaMemsetAsync(moments->normalMoments.data(), 0, size, stream),
                        "Leaf normal moment initialization");
                checked(cudaMemsetAsync(moments->massMoments.data(), 0, size, stream),
                        "Leaf mass moment initialization");
            }
            if (cachedFlux)
                checked(cudaMemsetAsync(moments->fluxMoments.data(), 0,
                                        product(product(cachedFlux, positive), sizeof(Z)), stream),
                        "Leaf flux moment initialization");
            for (unsigned first = 0; first < nodes;) {
                cancelled(stop);
                const unsigned count = std::min(waveCapacity, nodes - first);
                std::size_t maximumPanels = 0;
                for (unsigned local = 0; local < count; ++local) {
                    const auto node = first + local;
                    maximumPanels =
                        std::max(maximumPanels, sourceOffsets[node + 1] - sourceOffsets[node]);
                }
                for (std::size_t firstPanel = 0; firstPanel < maximumPanels;) {
                    cancelled(stop);
                    const auto panelCount =
                        std::min(panelsPerMomentSlice, maximumPanels - firstPanel);
                    buildLeafMoments<<<count, momentThreads, sharedBytes, stream>>>(
                        panels.data(), panelOffsets.data(), topology->slots.data(), rule.data(),
                        ruleCount, storageOrder, k, radius, cachedPressure, cachedFlux,
                        moments->normalMoments.data(), moments->massMoments.data(),
                        moments->fluxMoments.data(), first, firstPanel, panelCount, error.data());
                    checked(cudaGetLastError(), "Integrated shared-panel leaf moment construction");
                    checked(cudaStreamSynchronize(stream),
                            "Integrated shared-panel leaf moment panel slice");
                    cancelled(stop);
                    firstPanel += panelCount;
                }
                first += count;
            }
            finish(stop);
            if (cachedPressure) {
                verifySpectrum(moments->normalMoments.data(), product(cachedPressure, positive),
                               error.data(), stream, stop);
                verifySpectrum(moments->massMoments.data(), product(cachedPressure, positive),
                               error.data(), stream, stop);
            }
            if (cachedFlux)
                verifySpectrum(moments->fluxMoments.data(), product(cachedFlux, positive),
                               error.data(), stream, stop);
            momentStats.cachedPressureDofs = cachedPressure;
            momentStats.cachedFaceSlots = cachedFlux;
            momentStats.constructionMilliseconds =
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                    .count();
            moments->statistics = momentStats;
            moments->order = storageOrder;
            moments->radius = radius;
            moments->cachedPressure = cachedPressure;
            moments->cachedFlux = cachedFlux;
            // Publish only complete, finite tables. Failed construction leaves
            // no reusable entry and no host copy of integrated coefficients.
            if (family) {
                auto &bank = *family->impl_;
                // Stage both owning keys before publishing either map node;
                // allocation failure must not leave an unaccounted topology.
                decltype(bank.topologies) pendingTopologies;
                decltype(bank.entries) pendingEntries;
                if (!residentTopology) {
                    const auto entry = pendingTopologies.emplace(settings.nodes, topology).first;
                    topology->hostBytes = sum(topology->hostBytes, product(entry->first.capacity(),
                                                                           sizeof(std::uint32_t)));
                    momentStats.topologyHostBytes = topology->hostBytes;
                    moments->statistics = momentStats;
                }
                pendingEntries.emplace(settings.nodes, moments);
                const auto residentBytes = sum(bank.bytes, momentStats.bytes - residentTopology);
                bank.topologies.merge(pendingTopologies);
                bank.entries.merge(pendingEntries);
                bank.bytes = residentBytes;
            }
        }
    }
    void finish(std::stop_token stop) {
        int status = 0;
        checked(cudaMemcpyAsync(&status, error.data(), sizeof(int), cudaMemcpyDeviceToHost, stream),
                "Support leaf scalar status");
        checked(cudaStreamSynchronize(stream), "Support leaf completion");
        cancelled(stop);
        if (status)
            throw std::runtime_error(
                "Support leaf radial basis exceeded its finite numerical domain.");
    }
    template <bool pressure> void emit(const Z *input, Z *output, std::stop_token stop) {
        cancelled(stop);
        if (!input || !output || input == output)
            throw std::invalid_argument("Support leaf emission requires distinct device arrays.");
        checked(cudaMemsetAsync(error.data(), 0, sizeof(int), stream), "Support leaf status reset");
        if (spherical) {
            // Bounded node batches retain a cancellation boundary. One block
            // streams the complete panel union of one leaf.
            for (unsigned first = 0; first < nodes; first += 64) {
                cancelled(stop);
                const unsigned count = std::min(64u, nodes - first);
                emitSpherical<pressure><<<count, threads, (order + 2) * sizeof(double), stream>>>(
                    panels.data(), panelOffsets.data(), rule.data(), ruleCount, order, k, radius,
                    input, output, first, error.data());
                checked(cudaGetLastError(), "Support leaf spherical emission");
                checked(cudaStreamSynchronize(stream), "Support leaf spherical batch");
            }
        } else {
            const auto total = product(nodes, states);
            for (std::size_t first = 0; first < total; first += 16384) {
                cancelled(stop);
                const auto count = std::min(std::size_t{16384}, total - first);
                emitAngular<pressure>
                    <<<static_cast<unsigned>((count + threads - 1) / threads), threads, 0,
                       stream>>>(panels.data(), panelOffsets.data(), rule.data(), ruleCount,
                                 angles.data(), angleCount, k, input, output, first, count);
                checked(cudaGetLastError(), "Support leaf angular emission");
                checked(cudaStreamSynchronize(stream), "Support leaf angular batch");
            }
        }
        if (spherical && (pressure ? cachedPressure : cachedFlux)) {
            const auto total = product(nodes, states);
            for (std::size_t first = 0; first < total; first += maximumBatchPoints) {
                cancelled(stop);
                const auto count = std::min(maximumBatchPoints, total - first);
                emitIntegratedMoments<pressure>
                    <<<static_cast<unsigned>((count + threads - 1) / threads), threads, 0,
                       stream>>>(pressure ? moments->normalMoments.data()
                                          : moments->fluxMoments.data(),
                                 incidences.data(), testOffsets.data(), moments->fluxSlots.data(),
                                 pressure ? moments->nodePressureOffsets.data()
                                          : moments->nodeFluxOffsets.data(),
                                 pressure ? cachedPressure : cachedFlux, order, storageOrder, input,
                                 output, first, count);
                checked(cudaGetLastError(), "Integrated leaf emission contraction");
                checked(cudaStreamSynchronize(stream), "Integrated leaf emission batch");
            }
        }
        finish(stop);
        verifySpectrum(output, product(nodes, states), error.data(), stream, stop);
    }
    void test(const Z *input, Z *output, std::stop_token stop) {
        cancelled(stop);
        if (!input || !output || input == output)
            throw std::invalid_argument("Support leaf testing requires distinct device arrays.");
        checked(cudaMemsetAsync(error.data(), 0, sizeof(int), stream),
                "Support leaf test status reset");
        for (unsigned first = 0; first < cachedPressure; first += 256) {
            cancelled(stop);
            const unsigned count = std::min(256u, cachedPressure - first);
            testIntegratedMoments<<<count, threads, 0, stream>>>(
                moments->normalMoments.data(), moments->massMoments.data(), incidences.data(),
                testOffsets.data(), order, storageOrder, 1 / k, input, output, first, error.data());
            checked(cudaGetLastError(), "Integrated leaf testing contraction");
            checked(cudaStreamSynchronize(stream), "Integrated leaf testing batch");
        }
        for (unsigned first = cachedPressure; first < targets; first += 64) {
            cancelled(stop);
            const unsigned count = std::min(64u, targets - first);
            if (spherical)
                testSpherical<<<count, threads, (order + 2) * sizeof(double), stream>>>(
                    incidences.data(), testOffsets.data(), rule.data(), ruleCount, order, k, radius,
                    input, output, first, error.data());
            else
                testAngular<<<count, threads, 0, stream>>>(incidences.data(), testOffsets.data(),
                                                           rule.data(), ruleCount, angles.data(),
                                                           angleCount, k, input, output, first);
            checked(cudaGetLastError(), "Support leaf testing launch");
            checked(cudaStreamSynchronize(stream), "Support leaf testing batch");
        }
        finish(stop);
    }
};
SupportLeafOperator::SupportLeafOperator(const galerkin::Space &space, const SupportTree &tree,
                                         const SupportLeafSettings &settings, std::stop_token stop,
                                         cudaStream_t stream)
    : impl_(std::make_unique<Impl>(space, tree, settings, stop, stream)) {}
SupportLeafOperator::~SupportLeafOperator() = default;
void SupportLeafOperator::emitPressureDevice(const Z *pressure, Z *moments, std::stop_token stop) {
    impl_->emit<true>(pressure, moments, stop);
}
void SupportLeafOperator::emitFaceFluxDevice(const Z *flux, Z *moments, std::stop_token stop) {
    impl_->emit<false>(flux, moments, stop);
}
void SupportLeafOperator::testDevice(const Z *locals, Z *pressure, std::stop_token stop) {
    impl_->test(locals, pressure, stop);
}
std::size_t SupportLeafOperator::stateSize() const noexcept {
    return impl_->states;
}
std::size_t SupportLeafOperator::deviceBytes() const noexcept {
    return impl_->bytes + (impl_->family ? impl_->momentStats.bytes : 0);
}
std::size_t SupportLeafOperator::ownDeviceBytes() const noexcept {
    return impl_->bytes;
}
SupportLeafMomentStatistics SupportLeafOperator::momentStatistics() const noexcept {
    return impl_->momentStats;
}
} // namespace hornsculpt::wideband
