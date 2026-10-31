#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cuda/adaptive_patch_estimator.hpp>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <stop_token>
using namespace hornsculpt;
using namespace hornsculpt::wideband;
namespace {
using C = std::complex<double>;
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
template <class F> void rejects(F action) {
    bool rejected = false;
    try {
        action();
    } catch (const std::exception &) {
        rejected = true;
    }
    require(rejected, "Invalid adaptive estimator request was accepted.");
}
galerkin::Space fixture(double secondGap = .001) {
    SurfaceMesh mesh;
    unsigned component = 0;
    for (const Vec3 origin : std::array<Vec3, 2>{{{0, 0, 0}, {1e6, -2e6, 3e6}}}) {
        const auto n = static_cast<std::uint32_t>(mesh.vertices.size());
        const double gap = component++ ? secondGap : .001;
        for (const Vec3 point : std::array<Vec3, 6>{{{0, 0, 0},
                                                     {.01, 0, 0},
                                                     {0, .008, 0},
                                                     {.0003, .0007, gap},
                                                     {.009, .0007, gap},
                                                     {.0003, .009, gap}}})
            mesh.vertices.push_back({origin.x + point.x, origin.y + point.y, origin.z + point.z});
        mesh.triangles.push_back({{n, n + 1, n + 2}, 1, BoundaryCondition::driven});
        mesh.triangles.push_back({{n + 3, n + 4, n + 5}, 2, BoundaryCondition::driven});
    }
    return galerkin::buildSpace(mesh);
}
Vec3 local(const galerkin::Panel &p, const double *b) {
    return {(p.vertices[1].x - p.vertices[0].x) * b[1] + (p.vertices[2].x - p.vertices[0].x) * b[2],
            (p.vertices[1].y - p.vertices[0].y) * b[1] + (p.vertices[2].y - p.vertices[0].y) * b[2],
            (p.vertices[1].z - p.vertices[0].z) * b[1] +
                (p.vertices[2].z - p.vertices[0].z) * b[2]};
}
std::array<C, 36> moments(const galerkin::Space &space, const AdaptivePatchTask &task,
                          unsigned degree, double k) {
    const auto rule = galerkin::triangleRule(degree);
    std::array<C, 36> result{};
    const auto &t = space.panels[task.testPanel], &s = space.panels[task.sourcePanel];
    const double scale = std::max(t.maximumEdgeLength, s.maximumEdgeLength);
    for (const auto &x : rule)
        for (const auto &y : rule) {
            const auto p = galerkin::mapPairPatch(x, y, task.patch);
            const auto tx = local(t, p.test), sy = local(s, p.trial);
            const Vec3 d{t.vertices[0].x - s.vertices[0].x + tx.x - sy.x,
                         t.vertices[0].y - s.vertices[0].y + tx.y - sy.y,
                         t.vertices[0].z - s.vertices[0].z + tx.z - sy.z};
            const double r = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
            const auto g = std::exp(C(0, -k * r)) / r,
                       gradient = -g * C(1, k * r) * (scale / (r * r));
            const std::array<C, 4> kernels{g, gradient * d.x, gradient * d.y, gradient * d.z};
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j)
                    for (unsigned c = 0; c < 4; ++c)
                        result[c * 9 + 3 * i + j] +=
                            kernels[c] * (p.weight * p.test[i] * p.trial[j]);
        }
    return result;
}
bool samePatch(const galerkin::PairQuadraturePatch &a, const galerkin::PairQuadraturePatch &b) {
    return a.degree == b.degree && a.areaFraction == b.areaFraction &&
           std::equal(a.test, a.test + 9, b.test) && std::equal(a.trial, a.trial + 9, b.trial);
}
void requireSameRule(const galerkin::PairQuadratureRule &a, const galerkin::PairQuadratureRule &b) {
    require(a.maximumDepth == b.maximumDepth && a.patches.size() == b.patches.size() &&
                a.estimatedRelativeError == b.estimatedRelativeError,
            "Independent pair batching changed its accepted adaptive rule/error.");
    for (std::size_t i = 0; i < a.patches.size(); ++i)
        require(samePatch(a.patches[i], b.patches[i]),
                "Independent pair batching reordered or changed accepted patches.");
}
void capacityPlanGate() {
    const auto space = fixture(.05);
    constexpr double k = 36.6;
    galerkin::PairQuadratureSettings settings;
    settings.relativeTolerance = 1e-4;
    settings.maximumSubdivisionDepth = 8;
    struct Contact {
        galerkin::PairQuadraturePatch patch;
        AdaptivePatchMoments estimate;
    };
    std::array<std::vector<Contact>, 4> traces;
    std::array<galerkin::PairQuadratureRule, 4> accepted;
    // This driver observes every request before production commit. A single
    // active pair is the trace baseline; larger capacities interleave pairs.
    for (unsigned capacity : {1u, 512u, 4096u}) {
        galerkin::AdaptiveQuadratureWorkspace workspace(256 * 1024 * 1024);
        AdaptivePatchEstimator estimator(space, k, 4 * capacity);
        std::array<std::unique_ptr<galerkin::AdaptivePairQuadratureState>, 4> states;
        std::array<std::size_t, 4> cursor{};
        for (unsigned i = 0; i < states.size(); ++i) {
            const unsigned t = (i & 1) ? 2 : 0, s = t + 1;
            states[i] = std::make_unique<galerkin::AdaptivePairQuadratureState>(
                space.panels[(i & 2) ? s : t], space.panels[(i & 2) ? t : s], settings, k,
                workspace);
        }
        for (;;) {
            std::vector<AdaptivePatchTask> tasks;
            std::array<std::size_t, 4> first{}, counts{};
            unsigned active{};
            for (unsigned i = 0; i < states.size() && active < capacity; ++i) {
                if (states[i]->done())
                    continue;
                ++active;
                const unsigned t = (i & 1) ? 2 : 0, s = t + 1;
                first[i] = tasks.size();
                for (const auto &patch : states[i]->requests()) {
                    if (capacity != 1) {
                        require(cursor[i] + counts[i] < traces[i].size() &&
                                    samePatch(patch, traces[i][cursor[i] + counts[i]].patch),
                                "Adaptive batch capacity changed per-pair request order.");
                    }
                    tasks.push_back({(i & 2) ? s : t, (i & 2) ? t : s, patch});
                    ++counts[i];
                }
            }
            if (tasks.empty())
                break;
            const auto estimates = estimator.evaluateView(tasks);
            for (unsigned i = 0; i < states.size(); ++i) {
                if (!counts[i])
                    continue;
                for (std::size_t j = 0; j < counts[i]; ++j) {
                    const auto &task = tasks[first[i] + j];
                    const auto &estimate = estimates[first[i] + j];
                    if (capacity == 1) {
                        const auto fine = moments(space, task, task.patch.degree, k),
                                   coarse =
                                       moments(space, task, task.patch.degree == 8 ? 6 : 12, k);
                        double scale{}, error{}, difference{};
                        for (unsigned c = 0; c < 36; ++c) {
                            scale = std::max(scale, std::abs(fine[c]));
                            error = std::max(error, std::abs(fine[c] - coarse[c]));
                            difference =
                                std::max(difference, std::abs(estimate.integral[c] - fine[c]));
                        }
                        require(scale > 0 && difference / scale < 2e-12 &&
                                    std::abs(estimate.absoluteError - error) / scale < 2e-12,
                                "Adaptive request failed its independent 36-moment oracle.");
                        traces[i].push_back({task.patch, estimate});
                    } else {
                        const auto &baseline = traces[i][cursor[i] + j].estimate;
                        require(estimate.integral == baseline.integral &&
                                    estimate.absoluteError == baseline.absoluteError,
                                "Adaptive pair batching changed deterministic patch estimates.");
                    }
                }
                states[i]->commit(estimates.subspan(first[i], counts[i]));
                cursor[i] += counts[i];
            }
        }
        for (unsigned i = 0; i < states.size(); ++i) {
            require(cursor[i] == traces[i].size(), "Adaptive batching omitted patch requests.");
            if (capacity == 1)
                accepted[i] = states[i]->takeResult();
            else
                requireSameRule(states[i]->result(), accepted[i]);
        }
    }
    // Fill actual production queues too, including more than 512 independent
    // pairs. Two difficult pairs exercise refinement; translated easy pairs
    // fill the remaining slots without a redundant large CPU oracle workload.
    constexpr unsigned pairs = 4096;
    for (unsigned capacity : {1u, 512u, 4096u}) {
        AdaptivePairBatch batch(space, k, settings, 256 * 1024 * 1024,
                                std::numeric_limits<std::size_t>::max(), capacity);
        const auto idle = batch.statistics().liveWorkspaceBytes;
        std::vector<bool> seen(pairs);
        const auto kind = [](unsigned id) { return id < 4 ? id : 1u + 2 * (id & 1u); };
        const auto consume = [&](const SupportNearJob &job,
                                 const galerkin::PairQuadratureRule &rule) {
            const unsigned id = job.pressureMask | (unsigned(job.fluxMask) << 9);
            require(id < pairs && !seen[id], "Adaptive production queue lost job identity.");
            seen[id] = true;
            requireSameRule(rule, accepted[kind(id)]);
        };
        for (unsigned id = 0; id < pairs; ++id) {
            const auto i = kind(id), t = (i & 1) ? 2u : 0u, s = t + 1;
            const SupportNearJob job{(i & 2) ? s : t, (i & 2) ? t : s,
                                     static_cast<std::uint16_t>(id & 511u),
                                     static_cast<std::uint8_t>(id >> 9)};
            while (!batch.submit(job))
                batch.advance(consume);
        }
        batch.drain(consume);
        require(std::all_of(seen.begin(), seen.end(), [](bool value) { return value; }) &&
                    batch.statistics().completedPairs == pairs &&
                    batch.statistics().liveWorkspaceBytes == idle &&
                    batch.statistics().peakWorkspaceBytes <= 256 * 1024 * 1024,
                "Adaptive queue changed completion, idle payload or its host budget.");
        std::cout << "Adaptive pair capacity=" << capacity << " requests/accepted rules match, "
                  << "rounds=" << batch.statistics().rounds << '\n';
    }
}
void adaptivePlanGate() {
    const auto space = fixture();
    const double k = 2 * std::numbers::pi * 20000 / 343;
    galerkin::PairQuadratureSettings settings;
    settings.relativeTolerance = 1e-4;
    settings.maximumSubdivisionDepth = 8;
    std::array<galerkin::PairQuadratureRule, 4> reference;
    for (unsigned i = 0; i < 4; ++i) {
        const unsigned t = (i & 1) ? 2 : 0, s = t + 1;
        reference[i] = galerkin::makePairQuadrature(space.panels[(i & 2) ? s : t],
                                                    space.panels[(i & 2) ? t : s], settings, k);
    }
    for (unsigned i = 0; i < 4; ++i)
        std::cout << "CPU adaptive reference " << i << " patches=" << reference[i].patches.size()
                  << " peak=" << reference[i].peakWorkspaceBytes << std::endl;
    AdaptivePairBatch batch(space, k, settings, 256 * 1024 * 1024,
                            std::numeric_limits<std::size_t>::max(), 64);
    const auto idleWorkspaceBytes = batch.statistics().liveWorkspaceBytes;
    unsigned completed = 0;
    auto consume = [&](const SupportNearJob &job, const galerkin::PairQuadratureRule &rule) {
        const unsigned i = job.pressureMask;
        const auto &expected = reference[i];
        require(rule.patches.size() == expected.patches.size() &&
                    rule.maximumDepth == expected.maximumDepth,
                "CUDA adaptive state changed the converged patch count/depth.");
        auto key = [](const galerkin::PairQuadraturePatch &a,
                      const galerkin::PairQuadraturePatch &b) {
            for (unsigned j = 0; j < 9; ++j)
                if (a.test[j] != b.test[j])
                    return a.test[j] < b.test[j];
            for (unsigned j = 0; j < 9; ++j)
                if (a.trial[j] != b.trial[j])
                    return a.trial[j] < b.trial[j];
            return a.degree < b.degree;
        };
        auto actual = rule.patches, oracle = expected.patches;
        std::sort(actual.begin(), actual.end(), key);
        std::sort(oracle.begin(), oracle.end(), key);
        for (std::size_t j = 0; j < actual.size(); ++j) {
            require(actual[j].degree == oracle[j].degree &&
                        actual[j].areaFraction == oracle[j].areaFraction,
                    "CUDA adaptive state changed patch order/measure.");
            for (unsigned c = 0; c < 9; ++c)
                require(actual[j].test[c] == oracle[j].test[c] &&
                            actual[j].trial[c] == oracle[j].trial[c],
                        "CUDA adaptive state changed subdivision decisions.");
        }
        require(std::abs(rule.estimatedRelativeError - expected.estimatedRelativeError) < 1e-11,
                "CUDA adaptive state changed the global integral/error norm.");
        ++completed;
    };
    for (unsigned repeat = 0; repeat < 2; ++repeat) {
        for (unsigned j = 0; j < 64; ++j) {
            const unsigned i = j % 4, t = (i & 1) ? 2 : 0, s = t + 1;
            require(
                batch.submit({(i & 2) ? s : t, (i & 2) ? t : s, static_cast<std::uint16_t>(i), 0}),
                "Batch free slot was lost.");
        }
        require(!batch.submit({0, 1, 0, 0}), "Adaptive batch exceeded its live pair capacity.");
        batch.drain(consume);
        require(batch.statistics().liveWorkspaceBytes == idleWorkspaceBytes,
                "Completed adaptive states retained host workspace claims during reuse.");
    }
    const auto stats = batch.statistics();
    require(completed == 128 && stats.completedPairs == 128 && stats.maximumRoundPatches >= 64 &&
                stats.peakWorkspaceBytes <= 256 * 1024 * 1024,
            "Adaptive CUDA pair states were not batched/budgeted/reused.");
    std::cout << "CUDA adaptive plans pairs=" << completed << " rounds=" << stats.rounds
              << " patches=" << stats.evaluatedPatches << " maxRound=" << stats.maximumRoundPatches
              << " cpu_decision_ms=" << stats.decisionMilliseconds
              << " gpu_ms=" << stats.estimatorKernelMilliseconds
              << " estimator_wall_ms=" << stats.estimatorWallMilliseconds
              << " host_peak_bytes=" << stats.peakWorkspaceBytes << '\n';
    rejects([&] { AdaptivePairBatch tiny(space, k, settings, 1); });
    std::stop_source stop;
    stop.request_stop();
    AdaptivePairBatch stopped(space, k, settings, 1024 * 1024,
                              std::numeric_limits<std::size_t>::max(), 1, nullptr,
                              stop.get_token());
    rejects([&] { stopped.submit({0, 1, 0, 0}); });
}

void run() {
    const auto space = fixture();
    std::vector<AdaptivePatchTask> tasks;
    for (unsigned i = 0; i < 4096; ++i) {
        AdaptivePatchTask task;
        task.testPanel = (i & 1) ? 2 : 0;
        task.sourcePanel = task.testPanel + 1;
        if (i & 2)
            std::swap(task.testPanel, task.sourcePanel);
        task.patch.test[0] = task.patch.test[4] = task.patch.test[8] = 1;
        task.patch.trial[0] = task.patch.trial[4] = task.patch.trial[8] = 1;
        task.patch.areaFraction = 1;
        task.patch.degree = (i & 4) ? 14 : 8;
        if (i & 8) {
            const double child[9]{1, 0, 0, .5, .5, 0, .5, 0, .5};
            std::copy_n(child, 9, task.patch.test);
            task.patch.areaFraction = .25;
        }
        tasks.push_back(task);
    }
    for (const double k : {0., 36.6, 2 * std::numbers::pi * 20000 / 343}) {
        AdaptivePatchEstimator gpu(space, k, tasks.size());
        const auto actual = gpu.evaluate(tasks);
        double maximumMomentError = 0, maximumEstimatorError = 0;
        for (std::size_t i = 0; i < 16; ++i) {
            const auto fine = moments(space, tasks[i], tasks[i].patch.degree, k),
                       coarse = moments(space, tasks[i], tasks[i].patch.degree == 8 ? 6 : 12, k);
            double scale = 0, error = 0, difference = 0;
            for (unsigned j = 0; j < 36; ++j) {
                scale = std::max(scale, std::abs(fine[j]));
                error = std::max(error, std::abs(fine[j] - coarse[j]));
                difference = std::max(difference, std::abs(actual[i].integral[j] - fine[j]));
            }
            maximumMomentError = std::max(maximumMomentError, difference / scale);
            maximumEstimatorError =
                std::max(maximumEstimatorError, std::abs(actual[i].absoluteError - error) / scale);
        }
        for (std::size_t i = 16; i < tasks.size(); ++i)
            require(actual[i].integral == actual[i % 16].integral &&
                        actual[i].absoluteError == actual[i % 16].absoluteError,
                    "Estimator chunking changed repeated independent patch estimates.");
        require(gpu.lastChunkCount() > 1 &&
                    gpu.lastMaximumChunkGreenSamples() <=
                        AdaptivePatchEstimator::maximumChunkGreenSamples &&
                    gpu.hostBytes() >
                        tasks.size() * (sizeof(AdaptivePatchTask) + sizeof(AdaptivePatchMoments)),
                "Large estimator batch exceeded its sample bound or omitted pinned payload.");
        std::cout << "Batched CUDA adaptive moments k=" << k << " patches=" << tasks.size()
                  << " moment/estimator scaled error=" << maximumMomentError << '/'
                  << maximumEstimatorError << " kernel_ms=" << gpu.lastKernelMilliseconds()
                  << " owned_device_bytes=" << gpu.deviceBytes() << '\n';
        require(maximumMomentError < 2e-12 && maximumEstimatorError < 2e-12,
                "CUDA estimator changed the 36-moment norm/scaling.");
        const auto repeated = gpu.evaluate(tasks);
        require(repeated[0].integral == actual[0].integral,
                "Estimator reuse changed a deterministic patch.");
        require(gpu.evaluate({}).empty(), "Empty estimator batch produced data.");
        auto invalid = tasks;
        invalid.push_back(tasks[0]);
        rejects([&] { gpu.evaluate(invalid); });
        invalid = {tasks[0]};
        invalid[0].patch.degree = 4;
        rejects([&] { gpu.evaluate(invalid); });
        invalid[0] = tasks[0];
        invalid[0].patch.test[0] = std::numeric_limits<double>::quiet_NaN();
        rejects([&] { gpu.evaluate(invalid); });
        std::stop_source stop;
        stop.request_stop();
        rejects([&] { gpu.evaluate(tasks, stop.get_token()); });
        rejects([&] {
            AdaptivePatchEstimator rejected(space, k, tasks.size(), gpu.deviceBytes() - 1);
        });
    }
}
void projectedPatchGate() {
    const auto space = fixture(.00001);
    constexpr double k = 293.093;
    std::vector<AdaptivePatchTask> tasks;
    for (unsigned first : {0u, 2u})
        for (unsigned degree : {8u, 14u})
            for (unsigned intervals : {1u, 2u}) {
                galerkin::PairQuadraturePatch patch{};
                for (unsigned i = 0; i < 3; ++i)
                    patch.test[4 * i] = patch.trial[4 * i] = 1;
                patch.areaFraction = 1;
                patch.degree = degree;
                patch.projection = 1;
                patch.projectionIntervals = intervals;
                tasks.push_back({first, first + 1, patch});
            }
    AdaptivePatchEstimator gpu(space, k, tasks.size());
    const auto actual = gpu.evaluate(tasks);
    for (std::size_t i = 0; i < tasks.size(); ++i) {
        const auto expected = galerkin::estimatePairPatchCpu(
            space.panels[tasks[i].testPanel], space.panels[tasks[i].sourcePanel], tasks[i].patch, k);
        double norm = 0, difference = 0;
        for (unsigned j = 0; j < 36; ++j) {
            norm = std::max(norm, std::abs(expected.integral[j]));
            difference = std::max(difference, std::abs(actual[i].integral[j] - expected.integral[j]));
        }
        require(difference < 3e-12 * norm &&
                    std::abs(actual[i].innerError - expected.innerError) < 3e-12 * norm &&
                    std::abs(actual[i].outerError - expected.outerError) < 3e-12 * norm &&
                    std::abs(actual[i].absoluteError - expected.absoluteError) < 6e-12 * norm,
                "CUDA projected patch changed the inner/outer 14/20 and 6/8 or 12/14 estimates.");
    }
    std::cout << "Projected CUDA estimator: local/translated patches, degrees 8/14, intervals 1/2 passed.\n";
}
} // namespace
int main() {
    int devices = 0;
    if (cudaGetDeviceCount(&devices) != cudaSuccess || !devices) {
        std::cerr << "CUDA unavailable; adaptive estimator did not execute.\n";
        return 77;
    }
    try {
        run();
        projectedPatchGate();
        adaptivePlanGate();
        capacityPlanGate();
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
