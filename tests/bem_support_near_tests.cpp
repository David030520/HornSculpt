#include "fixture_geometry.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <fmm/support_near_plan.hpp>
#include <galerkin/pair_quadrature.hpp>
#include <iostream>
#include <set>
#include <stdexcept>
#include <stop_token>

using namespace hornsculpt;
using namespace hornsculpt::wideband;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class F> void rejects(F action) {
    bool failed = false;
    try {
        action();
    } catch (const std::exception &) {
        failed = true;
    }
    require(failed, "Invalid near-plan operation was accepted.");
}
SurfaceMesh meshFixture() {
    auto mesh = testing::sphere(1);
    const auto second = testing::sphere(0);
    const auto base = static_cast<std::uint32_t>(mesh.vertices.size());
    for (auto p : second.vertices)
        mesh.vertices.push_back({p.x + 2, p.y, p.z});
    for (auto t : second.triangles) {
        for (auto &v : t.vertices)
            v += base;
        mesh.triangles.push_back(t);
    }
    for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
        mesh.triangles[i].faceId = static_cast<std::uint32_t>(i % 3 + 1);
        mesh.triangles[i].boundary =
            i % 3 == 2 ? BoundaryCondition::wall : BoundaryCondition::driven;
    }
    return mesh;
}
struct Oracle {
    std::vector<unsigned char> pressureNear, fluxNear, pressureSlots, fluxSlots;
    std::vector<std::uint64_t> aSignatures, bSignatures;
    std::vector<double> aMass, bMass;
};
std::uint64_t signature(std::size_t t, std::size_t s, unsigned i, unsigned j) {
    // Integer identities test scatter/add ownership only. These are not
    // acoustic pressures or substitutes for the independent physics oracles.
    return 1000000 * (t + 1) + 1000 * (s + 1) + 10 * i + j + 1;
}
Oracle oracle(const galerkin::Space &space, const SupportTree &tree, const SupportNearPlan &plan) {
    const auto n = space.vertices.size(), p = space.panels.size(), f = plan.activeFaceIds.size();
    Oracle o;
    o.pressureNear.resize(n * n);
    o.fluxNear.resize(n * p);
    o.pressureSlots.resize(p * p * 9);
    o.fluxSlots.resize(p * p * 3);
    o.aSignatures.resize(n * n);
    o.bSignatures.resize(n * f);
    o.aMass.resize(n * n);
    o.bMass.resize(n * f);
    for (auto pair : tree.pressurePressure.nearInteractions) {
        const auto &t = tree.nodes[pair.target], &s = tree.nodes[pair.source];
        for (std::size_t a = t.firstPressure; a < std::size_t(t.firstPressure) + t.pressureCount;
             ++a)
            for (std::size_t b = s.firstPressure;
                 b < std::size_t(s.firstPressure) + s.pressureCount; ++b)
                require(++o.pressureNear[std::size_t(tree.pressureOrder[a]) * n +
                                         tree.pressureOrder[b]] == 1,
                        "Near tree duplicates A ownership.");
    }
    for (auto pair : tree.pressureFlux.nearInteractions) {
        const auto &t = tree.nodes[pair.target], &s = tree.nodes[pair.source];
        for (std::size_t a = t.firstPressure; a < std::size_t(t.firstPressure) + t.pressureCount;
             ++a)
            for (std::size_t b = s.firstFlux; b < std::size_t(s.firstFlux) + s.fluxCount; ++b)
                require(++o.fluxNear[std::size_t(tree.pressureOrder[a]) * p + tree.fluxOrder[b]] ==
                            1,
                        "Near tree duplicates B ownership.");
    }
    // Exhaustive panel Cartesian product and all local basis slots, independent
    // of the production iterator, its support panel lists and its bit masks.
    for (std::size_t t = 0; t < p; ++t)
        for (std::size_t s = 0; s < p; ++s) {
            const auto &tp = space.panels[t], &sp = space.panels[s];
            for (unsigned i = 0; i < 3; ++i) {
                const auto row = tp.nodes[i];
                for (unsigned j = 0; j < 3; ++j)
                    if (o.pressureNear[std::size_t(row) * n + sp.nodes[j]]) {
                        o.pressureSlots[(t * p + s) * 9 + 3 * i + j] = 1;
                        const auto index = std::size_t(row) * n + sp.nodes[j];
                        o.aSignatures[index] += signature(t, s, i, j);
                        if (t == s)
                            o.aMass[index] += tp.area * (i == j ? 2. : 1.) / 12;
                    }
                const auto column = plan.panelToActiveFaceColumn[s];
                if (column != UINT32_MAX && o.fluxNear[std::size_t(row) * p + s]) {
                    o.fluxSlots[(t * p + s) * 3 + i] = 1;
                    o.bSignatures[std::size_t(row) * f + column] += signature(t, s, i, 0);
                    if (t == s)
                        o.bMass[std::size_t(row) * f + column] += tp.area / 3;
                }
            }
        }
    return o;
}
std::size_t slot(const SupportCsr &csr, std::uint32_t row, std::uint32_t column) {
    const auto begin = csr.columns.begin() + csr.rowOffsets[row],
               end = csr.columns.begin() + csr.rowOffsets[std::size_t(row) + 1];
    const auto found = std::lower_bound(begin, end, column);
    require(found != end && *found == column, "Masked panel addend has no CSR destination.");
    return static_cast<std::size_t>(found - csr.columns.begin());
}
void check(const galerkin::Space &space, const SupportTree &tree, const SupportNearPlan &plan,
           std::size_t batchSize) {
    const auto expected = oracle(space, tree, plan);
    const auto n = space.vertices.size(), p = space.panels.size(), f = plan.activeFaceIds.size();
    std::vector<unsigned char> actualA(p * p * 9), actualB(p * p * 3);
    std::vector<std::uint64_t> a(plan.pressure.columns.size()), b(plan.flux.columns.size());
    std::vector<double> massA(a.size()), massB(b.size());
    std::vector<unsigned> pairJobs(p * p);
    std::vector<SupportNearJob> jobs(batchSize);
    SupportNearJobIterator iterator(space, tree, plan);
    std::size_t batches = 0;
    while (!iterator.done()) {
        const auto count = iterator.next(jobs, jobs.size() * sizeof(SupportNearJob));
        ++batches;
        require(count <= batchSize, "Iterator exceeded its explicit job buffer.");
        for (std::size_t k = 0; k < count; ++k) {
            const auto &job = jobs[k];
            const auto &tp = space.panels[job.testPanel], &sp = space.panels[job.sourcePanel];
            const auto pair = std::size_t(job.testPanel) * p + job.sourcePanel;
            require((job.pressureMask & ~0x1ffu) == 0 && (job.fluxMask & ~7u) == 0 &&
                        (job.pressureMask || job.fluxMask),
                    "Invalid pressure/flux local mask.");
            require(++pairJobs[pair] == 1,
                    "Panel pair was integrated repeatedly instead of merging masks.");
            for (unsigned i = 0; i < 3; ++i) {
                for (unsigned j = 0; j < 3; ++j)
                    if (job.pressureMask & (1u << (3 * i + j))) {
                        require(++actualA[pair * 9 + 3 * i + j] == 1,
                                "Pressure panel-local summand was duplicated.");
                        const auto destination = slot(plan.pressure, tp.nodes[i], sp.nodes[j]);
                        a[destination] += signature(job.testPanel, job.sourcePanel, i, j);
                        if (job.testPanel == job.sourcePanel)
                            massA[destination] += tp.massP1P1[3 * i + j];
                    }
                if (job.fluxMask & (1u << i)) {
                    require(++actualB[pair * 3 + i] == 1,
                            "Flux panel-local face addend was duplicated.");
                    const auto column = plan.panelToActiveFaceColumn[job.sourcePanel];
                    require(column < f, "Wall/inactive flux panel became a face RHS column.");
                    const auto destination = slot(plan.flux, tp.nodes[i], column);
                    b[destination] += signature(job.testPanel, job.sourcePanel, i, 0);
                    if (job.testPanel == job.sourcePanel)
                        massB[destination] += tp.massP1DP0[i];
                }
            }
        }
    }
    require(batches > 1 && actualA == expected.pressureSlots && actualB == expected.fluxSlots,
            "Bounded jobs omit legitimate local summands.");
    require(std::count(pairJobs.begin(), pairJobs.end(), 1u) > 1,
            "Fixture did not exercise distinct panel pairs.");
    for (std::size_t row = 0; row < n; ++row) {
        for (const auto *csr : {&plan.pressure, &plan.flux}) {
            const auto first = csr->columns.begin() + csr->rowOffsets[row],
                       last = csr->columns.begin() + csr->rowOffsets[row + 1];
            require(std::is_sorted(first, last) && std::adjacent_find(first, last) == last,
                    "CSR row is not sorted/unique.");
        }
        for (std::size_t column = 0; column < n; ++column) {
            const auto near = expected.pressureNear[row * n + column] != 0;
            const auto begin = plan.pressure.columns.begin() + plan.pressure.rowOffsets[row],
                       end = plan.pressure.columns.begin() + plan.pressure.rowOffsets[row + 1];
            require(std::binary_search(begin, end, static_cast<std::uint32_t>(column)) == near,
                    "A CSR sparsity differs from exact DOF-pair partition.");
            if (near) {
                const auto index = slot(plan.pressure, static_cast<std::uint32_t>(row),
                                        static_cast<std::uint32_t>(column));
                require(a[index] == expected.aSignatures[row * n + column] &&
                            std::abs(massA[index] - expected.aMass[row * n + column]) < 1e-14,
                        "A scatter dropped multiple panel contributions or mass slots.");
            }
        }
        for (std::size_t column = 0; column < f; ++column) {
            const auto expectedValue = expected.bSignatures[row * f + column];
            const auto begin = plan.flux.columns.begin() + plan.flux.rowOffsets[row],
                       end = plan.flux.columns.begin() + plan.flux.rowOffsets[row + 1];
            const bool exists = std::binary_search(begin, end, static_cast<std::uint32_t>(column));
            require(exists == (expectedValue != 0),
                    "B face CSR sparsity differs from exhaustive panel enumeration.");
            if (exists) {
                const auto index = slot(plan.flux, static_cast<std::uint32_t>(row),
                                        static_cast<std::uint32_t>(column));
                require(b[index] == expectedValue &&
                            std::abs(massB[index] - expected.bMass[row * f + column]) < 1e-14,
                        "Face aggregation dropped legitimate multiple-panel addends.");
            }
        }
        // Complex RHS reuse checks face-column order independently of scatter.
        for (unsigned rhs = 0; rhs < 2; ++rhs) {
            std::complex<double> actual{}, reference{};
            for (std::size_t column = 0; column < f; ++column) {
                const std::complex<double> velocity{.2 * (rhs + 1) + column,
                                                    .3 * (column + 1) - rhs};
                reference += static_cast<double>(expected.bSignatures[row * f + column]) * velocity;
                for (auto k = plan.flux.rowOffsets[row]; k < plan.flux.rowOffsets[row + 1]; ++k)
                    if (plan.flux.columns[k] == column)
                        actual += static_cast<double>(b[k]) * velocity;
            }
            require(std::abs(actual - reference) < 1e-12 * std::max(1., std::abs(reference)),
                    "Face-constant complex RHS reuse changed the aggregated action.");
        }
    }
}
void checkAdaptiveStates(const galerkin::Space &space) {
    std::vector<SupportNearJob> inputs;
    galerkin::PairQuadratureSettings settings;
    settings.relativeTolerance = 1e-4;
    for (std::uint32_t t = 0; t < space.panels.size() && inputs.size() < 16; ++t)
        for (std::uint32_t s = 0; s < space.panels.size() && inputs.size() < 16; ++s)
            if (galerkin::classifyPair(space.panels[t], space.panels[s]).topology ==
                    galerkin::PairTopology::disjoint &&
                galerkin::isNearPair(space.panels[t], space.panels[s], settings.nearDistanceRatio))
                inputs.push_back({t, s, 0x1ff, 7});
    require(inputs.size() == 16, "State fixture lacks real adaptive near pairs.");
    galerkin::AdaptiveQuadratureWorkspace workspace(8 * 1024 * 1024);
    std::vector<std::unique_ptr<galerkin::AdaptivePairQuadratureState>> states;
    for (const auto &input : inputs)
        states.push_back(std::make_unique<galerkin::AdaptivePairQuadratureState>(
            space.panels[input.testPanel], space.panels[input.sourcePanel], settings, 20,
            workspace));
    for (std::size_t i = 0; i < inputs.size(); ++i) {
        auto &state = *states[i];
        const auto &job = inputs[i];
        while (!state.done()) {
            std::array<galerkin::PairPatchEstimate, 4> evaluated;
            const auto requests = state.requests();
            for (std::size_t j = 0; j < requests.size(); ++j)
                evaluated[j] = galerkin::estimatePairPatchCpu(
                    space.panels[job.testPanel], space.panels[job.sourcePanel], requests[j], 20);
            state.commit(std::span(evaluated).first(requests.size()));
        }
        const auto expected = galerkin::makePairQuadrature(
            space.panels[job.testPanel], space.panels[job.sourcePanel], settings, 20);
        const auto &actual = state.result();
        require(actual.patches.size() == expected.patches.size() &&
                    actual.estimatedRelativeError == expected.estimatedRelativeError,
                "Externally advanced state changed the complete-pair norm.");
        for (std::size_t j = 0; j < actual.patches.size(); ++j) {
            const auto &a = actual.patches[j], &b = expected.patches[j];
            require(a.degree == b.degree && a.areaFraction == b.areaFraction &&
                        std::equal(std::begin(a.test), std::end(a.test), std::begin(b.test)) &&
                        std::equal(std::begin(a.trial), std::end(a.trial), std::begin(b.trial)),
                    "Externally advanced state changed the subdivision queue.");
        }
    }
    states.clear();
    require(workspace.liveBytes() == 0 && workspace.peakBytes() <= 8 * 1024 * 1024,
            "Shared adaptive state payload leaked/exceeded its hard bound.");
    galerkin::AdaptiveQuadratureWorkspace tiny(0);
    rejects([&] {
        galerkin::AdaptivePairQuadratureState invalid(space.panels[inputs[0].testPanel],
                                                      space.panels[inputs[0].sourcePanel], settings,
                                                      20, tiny);
    });
    std::stop_source stop;
    galerkin::AdaptivePairQuadratureState cancelled(space.panels[inputs[0].testPanel],
                                                    space.panels[inputs[0].sourcePanel], settings,
                                                    20, workspace, stop.get_token());
    stop.request_stop();
    rejects([&] { cancelled.commit({}); });
}

void run() {
    const auto space = galerkin::buildSpace(meshFixture());
    checkAdaptiveStates(space);
    SupportTreeSettings settings;
    settings.leafCapacity = 2;
    const auto tree = makeSupportTree(space, settings);
    const std::array<std::uint32_t, 2> faces{2, 1};
    const auto plan =
        makeSupportNearPlan(space, tree, faces, std::numeric_limits<std::size_t>::max());
    require(plan.activeFaceIds == std::vector<std::uint32_t>({2, 1}),
            "Caller face-column order was changed.");
    for (std::size_t p = 0; p < space.panels.size(); ++p)
        require(plan.panelToActiveFaceColumn[p] == (space.panels[p].faceId == 1   ? 1u
                                                    : space.panels[p].faceId == 2 ? 0u
                                                                                  : UINT32_MAX),
                "Panel-to-active-face mapping is incorrect.");
    for (std::size_t n = 0; n < tree.nodes.size(); ++n) {
        std::set<std::uint32_t> expected;
        const auto &node = tree.nodes[n];
        if (node.leaf())
            for (std::size_t i = node.firstPressure;
                 i < std::size_t(node.firstPressure) + node.pressureCount; ++i) {
                const auto vertex = tree.pressureOrder[i];
                for (auto j = space.incidenceOffsets[vertex];
                     j < space.incidenceOffsets[vertex + 1]; ++j)
                    expected.insert(space.incidences[j].panel);
            }
        const std::vector<std::uint32_t> actual(
            plan.supportPanels.begin() + plan.nodePanelOffsets[n],
            plan.supportPanels.begin() + plan.nodePanelOffsets[n + 1]);
        require(actual == std::vector<std::uint32_t>(expected.begin(), expected.end()),
                "Leaf support panel list is incomplete or not deduplicated.");
    }
    require(plan.supportPanels.size() <= space.incidences.size(),
            "Support panel references exceeded the 3T bound.");
    for (const auto batch : {1u, 7u, 32u})
        check(space, tree, plan, batch);
    const auto bounded = makeSupportNearPlan(space, tree, faces, plan.peakHostBytes);
    require(bounded.peakHostBytes == plan.peakHostBytes &&
                bounded.allocatedBytes() <= bounded.peakHostBytes,
            "Actual buffer-capacity budget is inconsistent.");
    rejects([&] { makeSupportNearPlan(space, tree, faces, plan.peakHostBytes - 1); });
    rejects([&] { makeSupportNearPlan(space, tree, faces, 0); });
    rejects([&] {
        makeSupportNearPlan(space, tree, std::array<std::uint32_t, 2>{1, 1}, plan.peakHostBytes);
    });
    rejects([&] {
        makeSupportNearPlan(space, tree, std::array<std::uint32_t, 1>{3}, plan.peakHostBytes);
    });
    rejects([&] {
        makeSupportNearPlan(space, tree, std::array<std::uint32_t, 1>{99}, plan.peakHostBytes);
    });
    const auto inactive = makeSupportNearPlan(space, tree, {}, plan.peakHostBytes);
    require(inactive.flux.columns.empty() && inactive.pressure.columns == plan.pressure.columns,
            "Inactive face RHS changed pressure sparsity.");
    check(space, tree, inactive, 13);
    std::stop_source stop;
    stop.request_stop();
    rejects([&] { makeSupportNearPlan(space, tree, faces, plan.peakHostBytes, stop.get_token()); });
    SupportNearJobIterator iterator(space, tree, plan);
    SupportNearJobIterator boundedIterator(space, tree, plan, iterator.workspaceBytes());
    require(boundedIterator.workspaceBytes() == iterator.workspaceBytes(),
            "Iterator workspace capacity accounting is inconsistent.");
    rejects(
        [&] { SupportNearJobIterator tooSmall(space, tree, plan, iterator.workspaceBytes() - 1); });
    rejects([&] { iterator.next({}, 0); });
    std::array<SupportNearJob, 2> jobs{};
    jobs[0].testPanel = UINT32_MAX;
    rejects([&] { iterator.next(jobs, sizeof(SupportNearJob)); });
    require(jobs[0].testPanel == UINT32_MAX, "Overbudget batch wrote output.");
    rejects([&] { iterator.next(jobs, sizeof(jobs), stop.get_token()); });
    require(jobs[0].testPanel == UINT32_MAX, "Cancelled batch wrote output.");
    require(iterator.next(jobs, sizeof(jobs)) == 2,
            "Budget rejection advanced or invalidated the iterator.");
    std::cout << "Complete support near CSR/masks: A=" << plan.pressure.columns.size()
              << ", Bface=" << plan.flux.columns.size() << ", owned bytes=" << plan.allocatedBytes()
              << ", peak including scratch=" << plan.peakHostBytes << ".\n";
}
} // namespace
int main() {
    try {
        run();
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
