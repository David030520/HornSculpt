#include "cuda/baffled_solver.hpp"
#include "fixture_geometry.hpp"
#include "galerkin/space.hpp"
#include "planning_frequencies.hpp"
#include <atomic>
#include <hornsculpt/axisymmetric_mesh.hpp>
#include <iostream>
#include <limits>
#include <thread>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class F> void rejects(F function) {
    bool rejected{};
    try {
        function();
    } catch (const std::exception &) {
        rejected = true;
    }
    require(rejected, "An invalid or cancelled plan operation was accepted.");
}
AxisymmetricMesh cylinder() {
    AxisymmetricMeridian m;
    m.points = {{0, -.006}, {.025, -.006}, {.025, 0}, {0, 0}};
    m.segments = {{0, 1, 1, BoundaryCondition::driven},
                  {1, 2, 2, BoundaryCondition::wall},
                  {2, 3, 3, BoundaryCondition::wall}};
    return buildAxisymmetricMesh(m, 8);
}
void keysAndPublication() {
    const auto geometry = cylinder();
    const auto &mesh = geometry.mesh;
    const auto space = galerkin::buildSpace(mesh);
    cuda_detail::BaffledPlanCache cache(mesh);
    galerkin::PairQuadratureSettings q;
    q.nearDistanceRatio = 0;
    const auto acquire = [&](const auto &settings, double k = 4.) {
        return cache.acquire(mesh, space, settings, k, {}, {});
    };
    auto first = acquire(q);
    require(first.plan && !first.reused && first.buildMilliseconds >= 0,
            "Cold host planning was not recorded.");
    auto warm = acquire(q);
    require(warm.reused && warm.plan == first.plan && warm.buildMilliseconds == 0,
            "An exact key failed to reuse the immutable plan.");
    for (unsigned key = 0; key < 5; ++key) {
        (void)acquire(q);
        auto changed = q;
        switch (key) {
        case 0:
            changed.regularOrder = 2;
            break;
        case 1:
            changed.singularOrder = 6;
            break;
        case 2:
            changed.nearDistanceRatio = .001;
            break;
        case 3:
            changed.relativeTolerance *= 2;
            break;
        case 4:
            ++changed.maximumSubdivisionDepth;
            break;
        }
        require(!acquire(changed).reused && acquire(changed).reused,
                "A changed quadrature key reused a mismatched plan.");
    }
    first = acquire(q, 10);
    require(!acquire(q, 5).reused && !acquire(q, 10).reused,
            "Plan reuse assumed an unproved monotonic wavenumber certificate.");
    require(acquire(q, 10).reused, "Identical maximum k was not reused.");
    auto bad = q;
    bad.singularOrder = 21;
    rejects([&] { (void)acquire(bad, 10); });
    require(acquire(q, 10).reused, "Failed planning destroyed the previous valid entry.");
    std::stop_source stop;
    stop.request_stop();
    rejects([&] { (void)cache.acquire(mesh, space, q, 10, {}, {}, stop.get_token()); });
    require(acquire(q, 10).reused, "Cancellation destroyed the previous valid entry.");
    auto otherMesh = mesh;
    rejects([&] { (void)cache.acquire(otherMesh, space, q, 10, {}, {}); });
    require(acquire(q, 10).reused, "A mesh identity mismatch replaced the valid entry.");

    AxisymmetricBem reduction{
        {}, {0, 0, 1}, 8, geometry.pressureDofs, geometry.fluxDofs, geometry.testPanels};
    const auto reduced = [&](const AxisymmetricBem &r) {
        return cache.acquire(mesh, space, q, 10, r, r.testPanels);
    };
    require(!reduced(reduction).reused && reduced(reduction).reused,
            "3D and reduced target plans were confused.");
    auto shiftedOrigin = reduction;
    shiftedOrigin.origin.z = .01;
    require(!reduced(shiftedOrigin).reused,
            "The explicit rotational frame was absent from the key.");
    auto relabelled = shiftedOrigin;
    for (auto &dof : relabelled.pressureDofs) {
        if (dof == 0)
            dof = 1;
        else if (dof == 1)
            dof = 0;
    }
    require(!reduced(relabelled).reused, "Pressure-orbit maps were absent from the key.");
    for (auto &dof : relabelled.fluxDofs) {
        if (dof == 0)
            dof = 1;
        else if (dof == 1)
            dof = 0;
    }
    require(!reduced(relabelled).reused, "Flux-strip maps were absent from the key.");
    std::reverse(relabelled.testPanels.begin(), relabelled.testPanels.end());
    require(!reduced(relabelled).reused && reduced(relabelled).reused,
            "Exact target-row selection was absent from the key.");
}
void cancelledReplacement() {
    const auto mesh = testing::sphere(3);
    const auto space = galerkin::buildSpace(mesh);
    cuda_detail::BaffledPlanCache cache(mesh);
    galerkin::PairQuadratureSettings q;
    q.nearDistanceRatio = 0;
    const auto original = cache.acquire(mesh, space, q, 1, {}, {});
    auto changed = q;
    changed.nearDistanceRatio = 1.3;
    std::stop_source stop;
    std::atomic_bool started{}, rejected{};
    std::jthread planner([&] {
        started.store(true);
        try {
            (void)cache.acquire(mesh, space, changed, 2, {}, {}, stop.get_token());
        } catch (const std::exception &) {
            rejected.store(true);
        }
    });
    while (!started.load())
        std::this_thread::yield();
    stop.request_stop();
    planner.join();
    require(rejected.load(), "An interrupted replacement published a completed entry.");
    const auto restored = cache.acquire(mesh, space, q, 1, {}, {});
    require(restored.reused && restored.plan == original.plan,
            "Interrupted planning replaced the last valid immutable plan.");
}
void planningRaster() {
    const std::vector<double> full{500, 2000, 16000}, low{500}, high{16000};
    const double expected = bem_detail::planningWavenumber(full, {}, 343);
    require(bem_detail::planningWavenumber(low, full, 343) == expected &&
                bem_detail::planningWavenumber(high, full, 343) == expected,
            "Scheduled subset planning differs from the complete sweep policy.");
    for (const auto &bad :
         {std::vector<double>{500, 500}, std::vector<double>{0}, std::vector<double>{-1},
          std::vector<double>{2000}, std::vector<double>{std::numeric_limits<double>::infinity()}})
        rejects([&] { (void)bem_detail::planningWavenumber(low, bad, 343); });
}
} // namespace
int main() {
    try {
        keysAndPublication();
        cancelledReplacement();
        planningRaster();
        std::cout << "Baffled host-plan exact keys, immutable publication, cancellation and raster "
                     "checks passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
