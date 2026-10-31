#include "bem_support_operator_oracle.hpp"
#include "fixture_geometry.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cuda/support_near_operator.hpp>
#include <galerkin/pair_quadrature.hpp>
#include <galerkin/regular_qualification.hpp>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <stop_token>

using namespace hornsculpt;
using namespace hornsculpt::wideband;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class F> void rejects(F f) {
    bool failed = false;
    try {
        f();
    } catch (const std::exception &) {
        failed = true;
    }
    require(failed, "Invalid near CUDA operation was accepted.");
}
class TestStream {
  public:
    TestStream() {
        require(cudaStreamCreateWithFlags(&value, cudaStreamNonBlocking) == cudaSuccess,
                "Cannot create borrowed near test stream.");
    }
    ~TestStream() {
        cudaStreamDestroy(value);
    }
    cudaStream_t value{};
};
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 delta(const galerkin::Panel &t, const galerkin::Panel &s,
           const galerkin::PairQuadraturePoint &p) {
    Vec3 result{};
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            const auto w = p.test[i] * p.trial[j];
            result.x += (t.vertices[i].x - s.vertices[j].x) * w;
            result.y += (t.vertices[i].y - s.vertices[j].y) * w;
            result.z += (t.vertices[i].z - s.vertices[j].z) * w;
        }
    return result;
}
struct Matrices {
    std::vector<Complex> a, b;
};
Matrices directReference(const galerkin::Space &space, const SupportTree &tree, double k,
                         galerkin::PairQuadratureSettings settings) {
    const auto n = space.vertices.size(), p = space.panels.size();
    Matrices result;
    result.a.resize(n * n);
    result.b.resize(n * 2);
    std::vector<unsigned char> nearA(n * n), nearB(n * p);
    for (auto pair : tree.pressurePressure.nearInteractions) {
        const auto &a = tree.nodes[pair.target], &b = tree.nodes[pair.source];
        for (std::size_t i = a.firstPressure; i < std::size_t(a.firstPressure) + a.pressureCount;
             ++i)
            for (std::size_t j = b.firstPressure;
                 j < std::size_t(b.firstPressure) + b.pressureCount; ++j)
                nearA[std::size_t(tree.pressureOrder[i]) * n + tree.pressureOrder[j]] = 1;
    }
    for (auto pair : tree.pressureFlux.nearInteractions) {
        const auto &a = tree.nodes[pair.target], &b = tree.nodes[pair.source];
        for (std::size_t i = a.firstPressure; i < std::size_t(a.firstPressure) + a.pressureCount;
             ++i)
            for (std::size_t j = b.firstFlux; j < std::size_t(b.firstFlux) + b.fluxCount; ++j)
                nearB[std::size_t(tree.pressureOrder[i]) * p + tree.fluxOrder[j]] = 1;
    }
    constexpr Complex imaginary{0, 1};
    namespace oracle = hornsculpt::test::support_operator_oracle;
    const auto regular12 = oracle::independentRegularRule(12),
               regular16 = oracle::independentRegularRule(16);
    std::vector<galerkin::RegularPanelBounds> bounds;
    for (const auto &panel : space.panels)
        bounds.push_back(galerkin::makeRegularPanelBounds(panel));
    galerkin::AdaptiveQuadratureWorkspace workspace(64 * 1024 * 1024);
    const auto relativeBlock = [](const auto &a, const auto &b) {
        double difference = 0, norm = 0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            difference += std::norm(a[i] - b[i]);
            norm += std::norm(b[i]);
        }
        return std::sqrt(difference / std::max(norm, 1e-300));
    };
    // Independent complex weak-form algebra. The established CPU quadrature
    // supplies identical integration domains; no CUDA moment helper is used.
    for (std::size_t ti = 0; ti < p; ++ti)
        for (std::size_t si = 0; si < p; ++si) {
            const auto &t = space.panels[ti], &s = space.panels[si];
            std::array<Complex, 9> a{};
            std::array<Complex, 3> b{};
            const auto integrate = [&](const galerkin::PairQuadraturePoint &point) {
                const auto d = delta(t, s, point);
                const auto r = std::sqrt(dot(d, d));
                const auto g = std::exp(-imaginary * k * r) / (4 * std::numbers::pi * r);
                const auto radial = g * (1. + imaginary * k * r) / (r * r);
                const auto sourceDerivative = ti == si ? Complex{} : radial * dot(d, s.normal);
                const auto targetDerivative = ti == si ? Complex{} : -radial * dot(d, t.normal);
                const auto weight = point.weight * 4 * t.area * s.area;
                for (unsigned i = 0; i < 3; ++i) {
                    b[i] += weight * point.test[i] * (-g + imaginary * targetDerivative / k);
                    for (unsigned j = 0; j < 3; ++j) {
                        const auto basis = point.test[i] * point.trial[j];
                        const auto weakH = (dot(t.curls[i], s.curls[j]) -
                                            k * k * dot(t.normal, s.normal) * basis) *
                                           g;
                        a[3 * i + j] +=
                            weight * (-basis * sourceDerivative - imaginary * weakH / k);
                    }
                }
            };
            const bool disjoint =
                galerkin::classifyPair(t, s).topology == galerkin::PairTopology::disjoint;
            const bool regular = disjoint &&
                                 !galerkin::isNearPair(t, s, settings.nearDistanceRatio) &&
                                 galerkin::qualifiedRegularDegree(bounds[ti], bounds[si], k,
                                                                  settings.relativeTolerance,
                                                                  settings.regularOrder) != 0;
            if (regular) {
                oracle::forIndependentRegularPair(regular16, integrate);
                const auto fineA = a;
                const auto fineB = b;
                a = {};
                b = {};
                oracle::forIndependentRegularPair(regular12, integrate);
                require(relativeBlock(a, fineA) < 1e-9 && relativeBlock(b, fineB) < 1e-9,
                        "Independent regular Duffy reference did not converge.");
                a = fineA;
                b = fineB;
            } else {
                galerkin::PairQuadratureRule rule;
                if (disjoint) {
                    galerkin::AdaptivePairQuadratureState state(t, s, settings, k, workspace);
                    while (!state.done()) {
                        const auto tasks = state.requests();
                        std::array<galerkin::PairPatchEstimate, 4> estimates{};
                        require(tasks.size() <= estimates.size(),
                                "Adaptive CPU request capacity changed.");
                        for (std::size_t j = 0; j < tasks.size(); ++j)
                            estimates[j] = galerkin::estimatePairPatchCpu(t, s, tasks[j], k);
                        state.commit(std::span(estimates).first(tasks.size()));
                    }
                    rule = state.takeResult();
                } else
                    rule = galerkin::makePairQuadrature(t, s, settings, k);
                for (const auto &point : rule.points)
                    integrate(point);
                for (const auto &patch : rule.patches) {
                    const auto quadrature = galerkin::triangleRule(patch.degree);
                    for (const auto &x : quadrature)
                        for (const auto &y : quadrature)
                            integrate(galerkin::mapPairPatch(x, y, patch));
                }
            }
            if (ti == si)
                for (unsigned i = 0; i < 3; ++i) {
                    b[i] += imaginary * t.area / (6 * k);
                    for (unsigned j = 0; j < 3; ++j)
                        a[3 * i + j] += .5 * t.area * (i == j ? 2. : 1.) / 12.;
                }
            const auto face = s.boundary == BoundaryCondition::wall ? UINT32_MAX
                              : s.faceId == 2                       ? 0u
                                                                    : 1u;
            for (unsigned i = 0; i < 3; ++i) {
                const auto row = t.nodes[i];
                for (unsigned j = 0; j < 3; ++j)
                    if (nearA[std::size_t(row) * n + s.nodes[j]])
                        result.a[std::size_t(row) * n + s.nodes[j]] += a[3 * i + j];
                if (face != UINT32_MAX && nearB[std::size_t(row) * p + si])
                    result.b[std::size_t(row) * 2 + face] += b[i];
            }
        }
    return result;
}
std::vector<Complex> dense(SupportDeviceCsr matrix) {
    std::vector<std::uint32_t> offsets(matrix.rows + 1), columns(matrix.nonzeros);
    std::vector<cuDoubleComplex> values(matrix.nonzeros);
    require(cudaMemcpy(offsets.data(), matrix.rowOffsets, offsets.size() * 4,
                       cudaMemcpyDeviceToHost) == cudaSuccess,
            "Download near row offsets failed.");
    if (matrix.nonzeros) {
        require(cudaMemcpy(columns.data(), matrix.columns, columns.size() * 4,
                           cudaMemcpyDeviceToHost) == cudaSuccess,
                "Download near columns failed.");
        require(cudaMemcpy(values.data(), matrix.values, values.size() * sizeof(cuDoubleComplex),
                           cudaMemcpyDeviceToHost) == cudaSuccess,
                "Download near values failed.");
    }
    std::vector<Complex> result(std::size_t(matrix.rows) * matrix.columnCount);
    for (unsigned row = 0; row < matrix.rows; ++row)
        for (auto j = offsets[row]; j < offsets[row + 1]; ++j)
            result[std::size_t(row) * matrix.columnCount + columns[j]] = {values[j].x, values[j].y};
    return result;
}
double relative(std::span<const Complex> a, std::span<const Complex> b) {
    double e = 0, n = 0;
    require(a.size() == b.size(), "Reference dimension mismatch.");
    for (std::size_t i = 0; i < a.size(); ++i) {
        e += std::norm(a[i] - b[i]);
        n += std::norm(b[i]);
    }
    return std::sqrt(e / std::max(n, 1e-300));
}
std::vector<Complex> apply(std::span<const Complex> matrix, std::span<const Complex> values,
                           std::size_t rows) {
    std::vector<Complex> result(rows);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < values.size(); ++j)
            result[i] += matrix[i * values.size() + j] * values[j];
    return result;
}
void run() {
    auto mesh = testing::sphere(0);
    const auto other = mesh;
    const auto base = static_cast<std::uint32_t>(mesh.vertices.size());
    for (auto p : other.vertices)
        mesh.vertices.push_back({p.x + .24, p.y, p.z});
    for (auto t : other.triangles) {
        for (auto &v : t.vertices)
            v += base;
        mesh.triangles.push_back(t);
    }
    for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
        mesh.triangles[i].faceId = static_cast<std::uint32_t>(i % 3 + 1);
        mesh.triangles[i].boundary =
            i % 3 == 2 ? BoundaryCondition::wall : BoundaryCondition::driven;
    }
    const auto space = galerkin::buildSpace(mesh);
    SupportTreeSettings treeSettings;
    treeSettings.leafCapacity = 1;
    treeSettings.nearDistanceRatio = .5;
    const auto tree = makeSupportTree(space, treeSettings);
    const std::array<std::uint32_t, 2> faces{2, 1};
    const auto plan = makeSupportNearPlan(space, tree, faces, 16 * 1024 * 1024);
    SupportNearOperatorSettings settings;
    settings.quadrature.nearDistanceRatio = .5;
    settings.quadrature.relativeTolerance = 1e-4;
    settings.quadrature.maximumSubdivisionDepth = 6;
    // This holds one adaptive patch in the bounded serialized device buffer,
    // forcing a single panel/mask integral to cross batch boundaries.
    settings.maximumConstructionBatchBytes = 256;
    settings.maximumBatchJobs = 7;
    constexpr double k = 20;
    const auto reference = directReference(space, tree, k, settings.quadrature);
    SupportNearOperator gpu(space, tree, plan, k, settings);
    const auto a = dense(gpu.pressureMatrix()), b = dense(gpu.faceFluxMatrix());
    const auto errorA = relative(a, reference.a), errorB = relative(b, reference.b);
    std::cout << "Real CUDA support near matrices: A relative=" << errorA
              << ", Bface relative=" << errorB << '\n';
    require(errorA < 3e-6 && errorB < 3e-6,
            "Near CUDA matrix differs from direct CPU weak-form reference.");
    const auto stats = gpu.constructionStats();
    require(stats.submittedBatches == stats.batches && stats.maximumResidentBatches == 2 &&
                stats.pinnedConstructionHostBytes ==
                    2 * (stats.deviceConstructionBatchBytes + sizeof(int)) &&
                stats.pinnedConstructionHostBytes <= stats.maximumPinnedConstructionHostBytes &&
                stats.maximumPinnedConstructionHostBytes ==
                    2 * settings.maximumConstructionBatchBytes,
            "Near pinned slots, device batch or completed FIFO accounting is inconsistent.");
    {
        auto serialSettings = settings;
        serialSettings.maximumActiveAdaptivePairs = 1;
        SupportNearOperator serial(space, tree, plan, k, serialSettings);
        require(relative(dense(serial.pressureMatrix()), a) < 1e-13 &&
                    relative(dense(serial.faceFluxMatrix()), b) < 1e-13,
                "CUDA adaptive batching changed the assembled near matrices.");
        require(stats.adaptiveRounds > 0 && stats.estimatorKernelMilliseconds > 0 &&
                    stats.adaptiveWorkspacePeakBytes <= settings.maximumAdaptiveWorkspaceBytes &&
                    stats.maximumSinglePairWorkspaceBytes > 0,
                "Near CUDA adaptive budget or measured pair workspace is inconsistent.");
    }
    require(stats.batches > 1 && stats.patches > 0 && stats.maximumBatchJobs <= 7 &&
                stats.maximumBatchPatches == 1 && stats.quadratureJobs > stats.maskedPanelJobs,
            "Near test did not exercise bounded adaptive-patch construction.");
    for (unsigned rhs = 0; rhs < 2; ++rhs) {
        std::vector<Complex> p(space.vertices.size());
        for (std::size_t i = 0; i < p.size(); ++i)
            p[i] = {std::sin(.7 * (i + rhs)), std::cos(.3 * (i + 2 * rhs))};
        const std::array<Complex, 2> q{{{.4 + rhs, .3 - rhs}, {-.2 * rhs, .9}}};
        const auto pressure = gpu.applyPressure(p), flux = gpu.applyFaceFlux(q);
        require(relative(pressure, apply(reference.a, p, p.size())) < 3e-6 &&
                    relative(flux, apply(reference.b, q, p.size())) < 3e-6,
                "Near multiple complex RHS action is inaccurate.");
        require(relative(gpu.applyPressure(p), pressure) < 1e-13,
                "Pressure/face-flux/pressure reuse changed the operator.");
    }
    const auto zero = gpu.applyPressure(std::vector<Complex>(space.vertices.size()));
    for (auto z : zero)
        require(z == Complex{}, "Near zero pressure reused stale state.");
    const auto zeroFlux = gpu.applyFaceFlux(std::array<Complex, 2>{});
    for (auto z : zeroFlux)
        require(z == Complex{}, "Near zero face RHS reused stale state.");
    auto strictBudget = settings;
    strictBudget.maximumDeviceBytes = gpu.peakDeviceBytes();
    {
        SupportNearOperator exactBudget(space, tree, plan, k, strictBudget);
        require(exactBudget.peakDeviceBytes() == gpu.peakDeviceBytes(),
                "Exact near device budget changed actual buffer accounting.");
    }
    strictBudget.maximumDeviceBytes = gpu.peakDeviceBytes() - 1;
    rejects([&] { SupportNearOperator rejected(space, tree, plan, k, strictBudget); });
    strictBudget.maximumConstructionBatchBytes = 1;
    strictBudget.maximumDeviceBytes = std::numeric_limits<std::size_t>::max();
    rejects([&] { SupportNearOperator rejected(space, tree, plan, k, strictBudget); });
    auto classificationBudget = settings;
    const auto classificationBytes = gpu.constructionStats().classificationWorkspaceBytes;
    require(classificationBytes > 0, "Panel qualification workspace was not accounted.");
    classificationBudget.maximumClassificationWorkspaceBytes = classificationBytes;
    {
        SupportNearOperator exact(space, tree, plan, k, classificationBudget);
        require(exact.constructionStats().classificationWorkspaceBytes == classificationBytes,
                "Exact panel qualification workspace budget changed payload accounting.");
    }
    classificationBudget.maximumClassificationWorkspaceBytes = classificationBytes - 1;
    rejects([&] { SupportNearOperator rejected(space, tree, plan, k, classificationBudget); });
    std::stop_source stop;
    stop.request_stop();
    rejects(
        [&] { SupportNearOperator rejected(space, tree, plan, k, settings, stop.get_token()); });
    rejects([&] { gpu.applyFaceFlux({}); });
    rejects(
        [&] { gpu.applyPressure(std::vector<Complex>(space.vertices.size()), stop.get_token()); });
    std::array<Complex, 2> invalid{{{std::numeric_limits<double>::quiet_NaN(), 0}, {}}};
    rejects([&] { gpu.applyFaceFlux(invalid); });
    require(gpu.deviceBytes() < gpu.peakDeviceBytes(),
            "Temporary near construction buffers were not released/accounted.");
    std::cout << "Real CUDA support near: " << stats.maskedPanelJobs << " masked panel jobs, "
              << stats.quadratureJobs << " integration chunks, " << stats.batches << " batches, "
              << stats.patches << " patches; live/peak bytes=" << gpu.deviceBytes() << '/'
              << gpu.peakDeviceBytes() << ".\n";
}
void checkQualifiedRegularRules() {
    struct Case {
        double separation, k;
    };
    for (const auto c : std::array<Case, 2>{{{1., 1.}, {.3, 8.}}}) {
        auto mesh = testing::sphere(0);
        const auto source = mesh;
        const auto offset = static_cast<std::uint32_t>(mesh.vertices.size());
        for (auto &triangle : mesh.triangles) {
            triangle.faceId = 1;
            triangle.boundary = BoundaryCondition::driven;
        }
        for (auto point : source.vertices) {
            point.z += c.separation;
            mesh.vertices.push_back(point);
        }
        for (auto triangle : source.triangles) {
            for (auto &vertex : triangle.vertices)
                vertex += offset;
            triangle.faceId = 2;
            triangle.boundary = BoundaryCondition::driven;
            mesh.triangles.push_back(triangle);
        }
        const auto space = galerkin::buildSpace(mesh);
        SupportTreeSettings treeSettings;
        treeSettings.leafCapacity = 64;
        const auto tree = makeSupportTree(space, treeSettings);
        const auto plan =
            makeSupportNearPlan(space, tree, std::array<std::uint32_t, 2>{2, 1}, 1024 * 1024);
        SupportNearOperatorSettings settings;
        settings.quadrature.nearDistanceRatio = 0;
        settings.quadrature.relativeTolerance = 1e-4;
        const auto reference = directReference(space, tree, c.k, settings.quadrature);
        SupportNearOperator gpu(space, tree, plan, c.k, settings);
        const auto counts = gpu.constructionStats().regularRuleJobs;
        require(c.separation > .5 ? counts[0] > 0 : counts[1] > 0 && counts[2] > 0,
                "Qualified regular CUDA rule selection did not exercise its intended degree.");
        const auto errorA = relative(dense(gpu.pressureMatrix()), reference.a);
        const auto errorB = relative(dense(gpu.faceFluxMatrix()), reference.b);
        require(
            errorA < 3e-6 && errorB < 3e-6,
            "Qualified regular near matrix disagrees with independent refined Duffy integration.");
        std::cout << "Qualified near degrees4/6/8=" << counts[0] << '/' << counts[1] << '/'
                  << counts[2] << " A/B=" << errorA << '/' << errorB << '\n';
    }
}
void checkProductionBatch() {
    auto mesh = testing::sphere(2);
    for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
        auto &triangle = mesh.triangles[i];
        triangle.faceId = static_cast<std::uint32_t>(i % 2 + 1);
        triangle.boundary = BoundaryCondition::driven;
    }
    const auto space = galerkin::buildSpace(mesh);
    SupportTreeSettings treeSettings;
    treeSettings.leafCapacity = 32;
    const auto tree = makeSupportTree(space, treeSettings);
    const auto plan =
        makeSupportNearPlan(space, tree, std::array<std::uint32_t, 2>{2, 1}, 1024 * 1024);
    SupportNearOperatorSettings settings;
    // The diagnostic zero-distance policy is explicit and identical in the
    // independent reference; it isolates regular/singular batch scheduling.
    settings.quadrature.nearDistanceRatio = 0;
    constexpr double k = 20;
    const auto reference = directReference(space, tree, k, settings.quadrature);
    SupportNearOperator gpu(space, tree, plan, k, settings);
    const auto stats = gpu.constructionStats();
    require(relative(dense(gpu.pressureMatrix()), reference.a) < 3e-6 &&
                relative(dense(gpu.faceFluxMatrix()), reference.b) < 3e-6,
            "Production near batches changed the direct Galerkin matrix.");
    require(stats.maximumBatchJobs > 256 &&
                stats.maximumBatchGreenSamples <= settings.maximumBatchGreenSamples &&
                stats.maximumKernelMilliseconds > 0,
            "Production near batch gate did not exercise measured bounded larger kernels.");
    auto workLimited = settings;
    workLimited.maximumBatchGreenSamples = 100000;
    SupportNearOperator limited(space, tree, plan, k, workLimited);
    require(limited.constructionStats().maximumBatchGreenSamples <= 100000 &&
                relative(dense(limited.pressureMatrix()), reference.a) < 3e-6,
            "Green sample work cap changed quadrature or failed to split batches.");
    const auto pipelineStats = limited.constructionStats();
    require(pipelineStats.submittedBatches == pipelineStats.batches &&
                pipelineStats.maximumResidentBatches == 2 &&
                pipelineStats.preparationsStartedWhileAssemblyPending > 0,
            "Near host preparation never began with an actually unfinished assembly batch.");
    TestStream borrowed;
    bool observedPendingCancellation = false;
    std::stop_source pendingStop;
    bool cancelledPending = false;
    try {
        SupportNearOperator interrupted(space, tree, plan, k, workLimited, pendingStop.get_token(),
                                        borrowed.value, [&](const auto &progress) {
                                            if (progress.submittedBatches > progress.batches &&
                                                cudaStreamQuery(borrowed.value) ==
                                                    cudaErrorNotReady) {
                                                observedPendingCancellation = true;
                                                pendingStop.request_stop();
                                            }
                                        });
    } catch (const std::runtime_error &error) {
        cancelledPending = std::string(error.what()).find("cancelled") != std::string::npos;
    }
    require(observedPendingCancellation && cancelledPending &&
                cudaStreamQuery(borrowed.value) == cudaSuccess,
            "Pending near cancellation did not drain its borrowed stream before unwinding.");
    {
        // A throwing progress callback is also a constructor failure with live
        // pinned uploads; destruction must drain before any staging memory dies.
        bool callbackThrew = false;
        try {
            SupportNearOperator interrupted(
                space, tree, plan, k, workLimited, {}, borrowed.value, [](const auto &progress) {
                    if (progress.submittedBatches > progress.batches)
                        throw std::runtime_error("Intentional near pending callback failure.");
                });
        } catch (const std::runtime_error &error) {
            callbackThrew = std::string(error.what()).find("Intentional") != std::string::npos;
        }
        require(callbackThrew && cudaStreamQuery(borrowed.value) == cudaSuccess,
                "Near constructor callback failure retained borrowed-stream work.");
    }
    {
        auto invalid = space;
        invalid.panels.front().curls[0].x = std::numeric_limits<double>::quiet_NaN();
        bool attributed = false;
        try {
            SupportNearOperator rejected(invalid, tree, plan, k, workLimited, {}, borrowed.value);
        } catch (const std::runtime_error &error) {
            const std::string message = error.what();
            attributed = message.find("failed in batch 1, status 2") != std::string::npos;
        }
        require(attributed && cudaStreamQuery(borrowed.value) == cudaSuccess,
                "Near non-finite scatter status lost its FIFO batch identity or stream drain.");
    }
    SupportNearOperator recovered(space, tree, plan, k, workLimited, {}, borrowed.value);
    require(
        relative(dense(recovered.pressureMatrix()), reference.a) < 3e-6 &&
            relative(dense(recovered.faceFluxMatrix()), reference.b) < 3e-6,
        "Borrowed stream did not recover after pending cancellation/constructor/status errors.");
    std::cout << "Near two-slot pipeline pinned/device bytes="
              << pipelineStats.pinnedConstructionHostBytes << '/'
              << pipelineStats.deviceConstructionBatchBytes << ", preparation starts while pending="
              << pipelineStats.preparationsStartedWhileAssemblyPending
              << "; FIFO status, pending cancellation, callback failure and borrowed reuse PASS\n";
    std::stop_source stop;
    rejects([&] {
        SupportNearOperator cancelled(space, tree, plan, k, settings, stop.get_token(), nullptr,
                                      [&](const auto &progress) {
                                          if (progress.batches)
                                              stop.request_stop();
                                      });
    });
    std::cout << "Production near batch jobs=" << stats.maximumBatchJobs
              << " Green samples=" << stats.maximumBatchGreenSamples
              << " maximum measured kernel_ms=" << stats.maximumKernelMilliseconds << '\n';
}
} // namespace
int main() {
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess || !n) {
        std::cerr << "CUDA hardware/runtime unavailable; near GPU test did not execute.\n";
        return 77;
    }
    try {
        run();
        checkQualifiedRegularRules();
        checkProductionBatch();
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
