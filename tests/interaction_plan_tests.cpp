#include "fixture_geometry.hpp"
#include <galerkin/interaction_plan.hpp>
#include <galerkin/near_projection.hpp>
#include <galerkin/regular_qualification.hpp>
#include <galerkin/space.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

using namespace hornsculpt;
using namespace hornsculpt::galerkin;

namespace {
using Pair = std::pair<std::uint32_t, std::uint32_t>;
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

void appendSphere(SurfaceMesh &mesh, Vec3 origin, double scale, unsigned refinement = 0) {
    auto other = hornsculpt::testing::sphere(refinement);
    const auto offset = static_cast<std::uint32_t>(mesh.vertices.size());
    for (const auto point : other.vertices)
        mesh.vertices.push_back(
            {origin.x + point.x * scale, origin.y + point.y * scale, origin.z + point.z * scale});
    for (auto triangle : other.triangles) {
        for (auto &vertex : triangle.vertices)
            vertex += offset;
        mesh.triangles.push_back(triangle);
    }
}

void checkMembership(const Space &space, const InteractionPlan &plan, double ratio, double k = .7,
                     double tolerance = 1e-5) {
    std::set<Pair> touching, near;
    std::vector<RegularPanelBounds> bounds;
    for (const auto &panel : space.panels)
        bounds.push_back(makeRegularPanelBounds(panel));
    for (std::uint32_t test = 0; test < space.panels.size(); ++test)
        for (std::uint32_t trial = 0; trial < space.panels.size(); ++trial) {
            const auto topology = classifyPair(space.panels[test], space.panels[trial]).topology;
            if (topology != PairTopology::disjoint)
                touching.emplace(test, trial);
            else if (isNearPair(space.panels[test], space.panels[trial], ratio) ||
                     qualifiedRegularDegree(bounds[test], bounds[trial], k, tolerance) == 0)
                near.emplace(test, trial);
        }
    std::set<Pair> actualTouching, actualNear;
    Pair previous{};
    bool first = true;
    for (const auto &job : plan.touching) {
        const Pair key{job.test, job.trial};
        require(first || previous < key, "Touching interactions are not unique and sorted");
        first = false;
        previous = key;
        actualTouching.insert(key);
        const auto expected = classifyPair(space.panels[job.test], space.panels[job.trial]);
        require(job.permutation.topology == expected.topology,
                "Touching interaction has an incorrect topology");
        for (unsigned local = 0; local < 3; ++local)
            require(job.permutation.test[local] == expected.test[local] &&
                        job.permutation.trial[local] == expected.trial[local],
                    "Touching interaction has an incorrect node permutation");
    }
    first = true;
    for (const auto &job : plan.near) {
        const Pair key{job.test, job.trial};
        require(first || previous < key, "Near interactions are not unique and sorted");
        first = false;
        previous = key;
        require(job.patchCount > 0 && job.firstPatch + job.patchCount <= plan.patches.size(),
                "Near interaction has an invalid patch range");
        require(job.swapBarycentric == (job.test > job.trial ? 1U : 0U),
                "Reverse interaction does not use the canonical patch direction");
        actualNear.insert(key);
        const auto reverse =
            std::find_if(plan.near.begin(), plan.near.end(), [&](const auto &other) {
                return other.test == job.trial && other.trial == job.test;
            });
        require(reverse != plan.near.end() && reverse->firstPatch == job.firstPatch &&
                    reverse->patchCount == job.patchCount,
                "Reverse near interaction duplicates or loses its patch storage");
    }
    const auto patchNear = actualNear;
    for (std::uint32_t trial = 1; trial < space.panels.size(); ++trial)
        for (std::uint32_t test = 0; test < trial; ++test)
            if (!plan.wholePanelDegrees.empty() &&
                plan.wholePanelDegrees[unorderedPanelPairIndex(test, trial)]) {
                require(!actualNear.contains({test, trial}) && !actualNear.contains({trial, test}),
                        "Full-panel override duplicates an adaptive near job.");
                actualNear.emplace(test, trial);
                actualNear.emplace(trial, test);
            }
    require(actualTouching == touching, "CSR touching plan differs from exhaustive classification");
    require(actualNear == near,
            "Spatial plan misses near pairs or includes bounding false positives");
    require(plan.nearOffsets.size() == space.panels.size() + 1 && plan.nearOffsets.front() == 0 &&
                plan.nearOffsets.back() == plan.near.size() &&
                plan.nearTrials.size() == plan.near.size(),
            "Near-membership CSR dimensions are inconsistent");
    std::set<Pair> csrNear;
    for (std::uint32_t test = 0; test < space.panels.size(); ++test) {
        require(plan.nearOffsets[test] <= plan.nearOffsets[test + 1],
                "Nonmonotone near CSR offsets");
        for (auto index = plan.nearOffsets[test]; index < plan.nearOffsets[test + 1]; ++index) {
            const auto trial = plan.nearTrials[index];
            require(index == plan.nearOffsets[test] || plan.nearTrials[index - 1] < trial,
                    "Near-membership CSR row is not strictly increasing");
            csrNear.emplace(test, trial);
        }
    }
    require(csrNear == patchNear, "GPU-membership CSR differs from adaptive patch jobs");
}

void comparePlans(const InteractionPlan &a, const InteractionPlan &b) {
    require(a.near.size() == b.near.size() && a.touching.size() == b.touching.size() &&
                a.patches.size() == b.patches.size() && a.nearOffsets == b.nearOffsets &&
                a.nearTrials == b.nearTrials && a.wholePanelDegrees == b.wholePanelDegrees &&
                a.maximumDepth == b.maximumDepth &&
                a.maximumEstimatedRelativeError == b.maximumEstimatedRelativeError,
            "Worker count changes the interaction-plan result");
    for (std::size_t index = 0; index < a.near.size(); ++index) {
        const auto &x = a.near[index], &y = b.near[index];
        require(x.test == y.test && x.trial == y.trial && x.firstPatch == y.firstPatch &&
                    x.patchCount == y.patchCount && x.swapBarycentric == y.swapBarycentric,
                "Parallel merge is not deterministic");
    }
    for (std::size_t index = 0; index < a.patches.size(); ++index) {
        const auto &x = a.patches[index], &y = b.patches[index];
        require(x.areaFraction == y.areaFraction && x.degree == y.degree &&
                    x.projection == y.projection && x.projectionIntervals == y.projectionIntervals,
                "Worker count changes patch weights or polynomial accuracy");
        for (unsigned coordinate = 0; coordinate < 9; ++coordinate)
            require(x.test[coordinate] == y.test[coordinate] &&
                        x.trial[coordinate] == y.trial[coordinate],
                    "Worker count changes patch barycentrics");
    }
}

std::array<std::complex<double>, 9> pairMoments(const Panel &test, const Panel &trial,
                                                const NearInteraction &job,
                                                const InteractionPlan &plan, double k) {
    std::array<std::complex<double>, 9> result{};
    for (std::uint32_t index = job.firstPatch; index < job.firstPatch + job.patchCount; ++index) {
        const auto rule = triangleRule(plan.patches[index].degree);
        const auto &patch = plan.patches[index];
        for (unsigned point = 0;
             point < nearPatchPointCount(patch, static_cast<unsigned>(rule.size())); ++point) {
            auto sample = mapNearPatchPoint(rule.data(), static_cast<unsigned>(rule.size()), patch,
                                            point, job.swapBarycentric ? trial : test,
                                            job.swapBarycentric ? test : trial);
            if (job.swapBarycentric)
                for (unsigned local = 0; local < 3; ++local)
                    std::swap(sample.test[local], sample.trial[local]);
            Vec3 delta{};
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j) {
                    const double weight = sample.test[i] * sample.trial[j];
                    delta.x += (test.vertices[i].x - trial.vertices[j].x) * weight;
                    delta.y += (test.vertices[i].y - trial.vertices[j].y) * weight;
                    delta.z += (test.vertices[i].z - trial.vertices[j].z) * weight;
                }
            const double distance = std::hypot(delta.x, delta.y, delta.z);
            const auto green =
                sample.weight * std::exp(std::complex<double>(0, -k * distance)) / distance;
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j)
                    result[3 * i + j] += green * sample.test[i] * sample.trial[j];
        }
    }
    return result;
}

void checkSwappedPatches(const Space &space, const InteractionPlan &plan) {
    for (const auto &job : plan.near) {
        if (job.test > job.trial)
            continue;
        const auto reverse =
            std::find_if(plan.near.begin(), plan.near.end(), [&](const auto &other) {
                return other.test == job.trial && other.trial == job.test;
            });
        const auto forwardMoments =
            pairMoments(space.panels[job.test], space.panels[job.trial], job, plan, .7);
        const auto reverseMoments =
            pairMoments(space.panels[job.trial], space.panels[job.test], *reverse, plan, .7);
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j)
                require(std::abs(forwardMoments[3 * i + j] / reverseMoments[3 * j + i] - 1.0) <
                            2e-12,
                        "Reverse patch metadata violates P1 reciprocity");
    }
}

void membershipAndParallelTests() {
    SurfaceMesh mesh;
    appendSphere(mesh, {0, 0, 0}, 1);
    appendSphere(mesh, {.205, 0, 0}, .001);
    appendSphere(mesh, {5, 0, 0}, .4);
    const auto space = buildSpace(mesh);
    PairQuadratureSettings settings;
    settings.relativeTolerance = 1e-4;
    settings.maximumSubdivisionDepth = 6;
    const auto serial = makeInteractionPlan(space, settings, .7, {}, 1),
               parallel = makeInteractionPlan(space, settings, .7, {}, 4);
    checkMembership(space, serial, settings.nearDistanceRatio, .7, settings.relativeTolerance);
    comparePlans(serial, parallel);
    require(serial.candidatePairs > serial.near.size() / 2,
            "Fixture does not exercise rejection of spatial false positives");
    require(serial.maximumEstimatedRelativeError <= settings.relativeTolerance,
            "Plan assembly silently relaxes quadrature tolerance");
    checkSwappedPatches(space, serial);
    settings.nearDistanceRatio = 0;
    const auto disabled = makeInteractionPlan(space, settings, .7);
    checkMembership(space, disabled, 0, .7, settings.relativeTolerance);
    require(!disabled.near.empty() && !disabled.patches.empty(),
            "Disabling geometric near correction must retain unqualified adaptive pairs");
    std::cout << "Variable-size plan: " << serial.candidatePairs << " candidates, "
              << serial.near.size() / 2 << " unique near pairs, " << serial.patches.size()
              << " shared patches\n";
}

void sparseCandidateTest() {
    SurfaceMesh mesh;
    for (unsigned body = 0; body < 20; ++body)
        appendSphere(mesh, {body * 2.0, 0, 0}, 1);
    const auto space = buildSpace(mesh);
    const auto plan = makeInteractionPlan(space, {}, .7, {}, 4);
    checkMembership(space, plan, 1);
    require(plan.candidatePairs < space.panels.size() * space.panels.size() / 8,
            "Spatial index regressed to enumerating all panel pairs");
}

void targetSectorTest() {
    SurfaceMesh mesh;
    appendSphere(mesh, {0, 0, 0}, 1);
    appendSphere(mesh, {.205, 0, 0}, .001);
    const auto space = buildSpace(mesh);
    PairQuadratureSettings settings;
    settings.relativeTolerance = 1e-4;
    const auto full = makeInteractionPlan(space, settings, .7, {}, 1);
    const std::vector<std::uint32_t> targets{0, 3, 9};
    const auto selected =
        makeInteractionPlan(space, settings, .7, {}, 2, std::size_t{1024} * 1024 * 1024, targets);
    const auto active = [&](unsigned test) {
        return std::find(targets.begin(), targets.end(), test) != targets.end();
    };
    std::set<Pair> expectedTouching, expectedNear, actualTouching, actualNear;
    for (const auto &job : full.touching)
        if (active(job.test))
            expectedTouching.emplace(job.test, job.trial);
    for (const auto &job : full.near)
        if (active(job.test))
            expectedNear.emplace(job.test, job.trial);
    for (const auto &job : selected.touching)
        actualTouching.emplace(job.test, job.trial);
    for (const auto &job : selected.near) {
        actualNear.emplace(job.test, job.trial);
        const auto reference = std::find_if(full.near.begin(), full.near.end(), [&](const auto &v) {
            return v.test == job.test && v.trial == job.trial;
        });
        require(reference != full.near.end(), "Selected target creates an unplanned near pair.");
        const auto a =
            pairMoments(space.panels[job.test], space.panels[job.trial], job, selected, .7);
        const auto b =
            pairMoments(space.panels[job.test], space.panels[job.trial], *reference, full, .7);
        for (unsigned i = 0; i < a.size(); ++i)
            require(a[i] == b[i], "Selected target changes adaptive source integration.");
    }
    require(actualTouching == expectedTouching && actualNear == expectedNear,
            "Selected target plan loses or adds source interactions.");
    for (std::size_t panel = 0; panel < space.panels.size(); ++panel)
        if (!active(static_cast<unsigned>(panel)))
            require(selected.nearOffsets[panel] == selected.nearOffsets[panel + 1],
                    "Inactive target retains near quadrature jobs.");
    for (auto invalid :
         {std::vector<std::uint32_t>{0, 0},
          std::vector<std::uint32_t>{static_cast<std::uint32_t>(space.panels.size())}}) {
        bool rejected{};
        try {
            (void)makeInteractionPlan(space, settings, .7, {}, 1, std::size_t{1024} * 1024 * 1024,
                                      invalid);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Invalid target-panel restriction was accepted.");
    }
}

void coordinateRangeTest() {
    SurfaceMesh mesh;
    // Very distant, individually small planar panels exercise the floating-point
    // margin in spatial cell sizing without requiring a huge panel edge.
    for (const auto [x, offset] :
         std::array<std::array<double, 2>, 3>{{{0, 0}, {1e20, 0}, {1e20, .012}}}) {
        const auto first = static_cast<std::uint32_t>(mesh.vertices.size());
        mesh.vertices.push_back({x, offset, 0});
        mesh.vertices.push_back({x, offset + .008, 0});
        mesh.vertices.push_back({x, offset, .008});
        mesh.triangles.push_back({{first, first + 1, first + 2}, 1});
    }
    const auto space = buildSpace(mesh);
    const auto plan = makeInteractionPlan(space, {}, .7, {}, 2);
    checkMembership(space, plan, 1);
    require(plan.near.size() == 2, "Large coordinate range hides a true near pair");
}

void retainedRulesBudgetTest() {
    SurfaceMesh mesh;
    // Many individually inexpensive rules must fit independently of how rows
    // are distributed among workers. Geometric near selection retains patches
    // even when electrical qualification permits a spatial candidate search.
    constexpr unsigned panelCount = 40;
    for (unsigned panel = 0; panel < panelCount; ++panel) {
        const auto first = static_cast<std::uint32_t>(mesh.vertices.size());
        const double z = .04 * panel;
        mesh.vertices.push_back({0, 0, z});
        mesh.vertices.push_back({.002, 0, z});
        mesh.vertices.push_back({0, .002, z});
        mesh.triangles.push_back({{first, first + 1, first + 2}, panel + 1});
    }
    const auto space = buildSpace(mesh);
    PairQuadratureSettings settings;
    settings.relativeTolerance = 1e-4;
    settings.maximumSubdivisionDepth = 6;
    settings.nearDistanceRatio = 1000;
    constexpr double k = 366.366;
    constexpr std::size_t budget = 1024 * 1024;
    const auto serial = makeInteractionPlan(space, settings, k, {}, 1, budget);
    const auto parallel = makeInteractionPlan(space, settings, k, {}, 8, budget);
    comparePlans(serial, parallel);
    require(parallel.near.size() == panelCount * (panelCount - 1),
            "Memory-bounded planning dropped retained pair rules.");
    require(parallel.maximumEstimatedRelativeError <= .2 * settings.relativeTolerance,
            "Memory-bounded planning relaxed the adaptive stopping tolerance.");
    std::cout << "Retained-rule budget: " << parallel.near.size() / 2 << " pairs, "
              << parallel.patches.size() << " patches, " << budget << " bytes\n";
}

void wholePanelRulesBudgetTest() {
    SurfaceMesh mesh;
    constexpr unsigned panelCount = 40;
    for (unsigned panel = 0; panel < panelCount; ++panel) {
        const auto first = static_cast<std::uint32_t>(mesh.vertices.size());
        const double z = .04 * panel;
        mesh.vertices.push_back({0, 0, z});
        mesh.vertices.push_back({.004, 0, z});
        mesh.vertices.push_back({0, .004, z});
        mesh.triangles.push_back({{first, first + 1, first + 2}, panel + 1});
    }
    const auto space = buildSpace(mesh);
    PairQuadratureSettings settings;
    settings.relativeTolerance = 1e-4;
    settings.maximumSubdivisionDepth = 6;
    constexpr double k = 366.366;
    const auto serial = makeInteractionPlan(space, settings, k, {}, 1, 1024 * 1024);
    const auto parallel = makeInteractionPlan(space, settings, k, {}, 8, 1024 * 1024);
    comparePlans(serial, parallel);
    checkMembership(space, parallel, settings.nearDistanceRatio, k, settings.relativeTolerance);
    require(parallel.wholePanelDegrees.size() == panelCount * (panelCount - 1) / 2 &&
                parallel.near.empty() && parallel.patches.empty(),
            "Unsubdivided rules retain duplicate full-panel patch/job storage.");
    for (unsigned trial = 1; trial < panelCount; ++trial)
        for (unsigned test = 0; test < trial; ++test) {
            const auto rule =
                makePairQuadrature(space.panels[test], space.panels[trial], settings, k, {}, true);
            const auto degree = parallel.wholePanelDegrees[unorderedPanelPairIndex(test, trial)];
            require(rule.patches.size() == 1 && rule.patches.front().degree == degree &&
                        rule.maximumDepth == 0 && !rule.patches.front().projection,
                    "Compact degree differs from independently accepted adaptive rule.");
            require(unorderedPanelPairIndex(test, trial) == unorderedPanelPairIndex(trial, test),
                    "Compact rule is not symmetric under pair reversal.");
        }
    std::cout << "Full-panel rule budget: " << parallel.wholePanelDegrees.size()
              << " bytes for 780 independently accepted pair rules\n";
}

void failureAndCancellationTests() {
    const auto space = buildSpace(hornsculpt::testing::sphere(0));
    auto insufficient = PairQuadratureSettings{};
    bool budgetRejected = false;
    try {
        makeInteractionPlan(space, insufficient, .7, {}, 4, 128);
    } catch (const std::runtime_error &error) {
        budgetRejected = std::string(error.what()).find("budget") != std::string::npos;
    }
    require(budgetRejected, "Dense plan exceeded its hard host payload budget.");
    insufficient.relativeTolerance = 1e-10;
    insufficient.maximumSubdivisionDepth = 0;
    bool rejected = false;
    try {
        makeInteractionPlan(space, insufficient, 400, {}, 4);
    } catch (const std::runtime_error &error) {
        const std::string message(error.what());
        rejected = message.find("subdivision limit") != std::string::npos &&
                   message.find("Near pair panels ") != std::string::npos &&
                   message.find("faces ") != std::string::npos &&
                   message.find("k=400") != std::string::npos &&
                   message.find("gap=") != std::string::npos &&
                   message.find("depth=0") != std::string::npos;
    }
    require(rejected, "Worker quadrature failure was swallowed or silently downgraded");
    std::stop_source stopped;
    stopped.request_stop();
    rejected = false;
    try {
        makeInteractionPlan(space, {}, .7, stopped.get_token());
    } catch (const std::runtime_error &error) {
        rejected = std::string(error.what()).find("cancelled") != std::string::npos;
    }
    require(rejected, "Pre-cancelled interaction plan was computed");

    const auto larger = buildSpace(hornsculpt::testing::sphere(3));
    std::stop_source during;
    std::jthread trigger([&] {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        during.request_stop();
    });
    rejected = false;
    try {
        makeInteractionPlan(larger, {}, .7, during.get_token(), 4);
    } catch (const std::runtime_error &error) {
        rejected = std::string(error.what()).find("cancelled") != std::string::npos;
    }
    require(rejected, "Active interaction planning ignored cancellation");
}
} // namespace

int main() {
    try {
        membershipAndParallelTests();
        sparseCandidateTest();
        targetSectorTest();
        coordinateRangeTest();
        retainedRulesBudgetTest();
        wholePanelRulesBudgetTest();
        failureAndCancellationTests();
        std::cout
            << "Interaction membership, shared patches, workers and cancellation tests passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
