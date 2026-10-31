#pragma once

#include "../galerkin/axisymmetric_quadrature.hpp"
#include "../galerkin/interaction_plan.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <hornsculpt/bem.hpp>
#include <memory>
#include <stdexcept>
#include <tuple>

namespace hornsculpt::cuda_detail {
// One solver owns one immutable mesh. A failed or cancelled replacement leaves
// the previous complete host plan usable for its original numerical key.
class BaffledPlanCache {
  public:
    explicit BaffledPlanCache(const SurfaceMesh &mesh) : mesh_(&mesh) {}
    struct Lookup {
        std::shared_ptr<const galerkin::InteractionPlan> plan;
        double buildMilliseconds{};
        bool reused{};
    };
    struct RingLookup {
        std::shared_ptr<const galerkin::RingInteractionPlan> plan;
        double buildMilliseconds{};
        bool reused{};
    };
    RingLookup acquireRing(const SurfaceMesh &mesh, const galerkin::AxisymmetricSpace &space,
                           const BemNumerics &numerics, const AxisymmetricBem &reduction,
                           std::stop_token stop = {}) {
        if (&mesh != mesh_)
            throw std::logic_error("Baffled ring cache belongs to a different immutable mesh.");
        checkCancelled(stop);
        galerkin::PairQuadratureSettings settings;
        settings.regularOrder = numerics.quadratureOrder;
        settings.singularOrder = numerics.singularOrder;
        settings.nearDistanceRatio = numerics.nearDistanceRatio;
        settings.relativeTolerance = numerics.nearRelativeTolerance;
        settings.maximumSubdivisionDepth = numerics.nearMaximumDepth;
        if (ringEntry_ && sameQuadrature(ringEntry_->settings, settings) &&
            ringEntry_->reduction == reduction)
            return {ringEntry_->plan, 0, true};
        const auto start = std::chrono::steady_clock::now();
        auto plan = std::make_shared<const galerkin::RingInteractionPlan>(
            galerkin::makeRingInteractionPlan(space, numerics, stop));
        const double milliseconds =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                .count();
        auto entry = std::make_unique<RingEntry>(RingEntry{settings, reduction, std::move(plan)});
        checkCancelled(stop);
        ringEntry_ = std::move(entry);
        return {ringEntry_->plan, milliseconds, false};
    }
    Lookup acquire(const SurfaceMesh &mesh, const galerkin::Space &space,
                   const galerkin::PairQuadratureSettings &settings, double maximumWavenumber,
                   const std::optional<AxisymmetricBem> &axisymmetric,
                   std::span<const std::uint32_t> testPanels, std::stop_token stop = {}) {
        if (&mesh != mesh_)
            throw std::logic_error("Baffled plan cache belongs to a different immutable mesh.");
        checkCancelled(stop);
        if (!std::isfinite(maximumWavenumber) || maximumWavenumber < 0)
            throw std::invalid_argument("Baffled plan wavenumber must be finite and nonnegative.");
        if (entry_ && sameQuadrature(entry_->settings, settings) &&
            entry_->axisymmetric == axisymmetric &&
            std::equal(entry_->testPanels.begin(), entry_->testPanels.end(), testPanels.begin(),
                       testPanels.end()) &&
            entry_->maximumWavenumber == maximumWavenumber)
            return {entry_->plan, 0, true};
        const auto start = std::chrono::steady_clock::now();
        auto plan = std::make_shared<const galerkin::InteractionPlan>(
            galerkin::makeInteractionPlan(space, settings, maximumWavenumber, stop, 0, testPanels));
        const double milliseconds =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                .count();
        auto candidate = std::make_unique<Entry>(Entry{
            settings, maximumWavenumber, axisymmetric,
            std::vector<std::uint32_t>(testPanels.begin(), testPanels.end()), std::move(plan)});
        checkCancelled(stop);
        entry_ = std::move(candidate);
        return {entry_->plan, milliseconds, false};
    }

  private:
    static bool sameQuadrature(const galerkin::PairQuadratureSettings &a,
                               const galerkin::PairQuadratureSettings &b) {
        return std::tie(a.regularOrder, a.singularOrder, a.nearDistanceRatio, a.relativeTolerance,
                        a.maximumSubdivisionDepth) ==
               std::tie(b.regularOrder, b.singularOrder, b.nearDistanceRatio, b.relativeTolerance,
                        b.maximumSubdivisionDepth);
    }
    static void checkCancelled(std::stop_token stop) {
        if (stop.stop_requested())
            throw std::runtime_error("Baffled interaction planning cancelled.");
    }
    struct Entry {
        galerkin::PairQuadratureSettings settings;
        double maximumWavenumber{};
        std::optional<AxisymmetricBem> axisymmetric;
        std::vector<std::uint32_t> testPanels;
        std::shared_ptr<const galerkin::InteractionPlan> plan;
    };
    struct RingEntry {
        galerkin::PairQuadratureSettings settings;
        AxisymmetricBem reduction;
        // Geometry and rules only: every frequency certifies its own operator
        // during CUDA assembly, so this cache needs no frequency assumption.
        std::shared_ptr<const galerkin::RingInteractionPlan> plan;
    };
    const SurfaceMesh *mesh_;
    std::unique_ptr<const Entry> entry_;
    std::unique_ptr<const RingEntry> ringEntry_;
};

// The closed cavity mesh has normals out of its air volume. Coordinates are
// local to worldOrigin; the public request retains world-space coordinates.
BemResult solveBaffledCavity(const SurfaceMesh &mesh, Vec3 worldOrigin, const BemRequest &request,
                             BaffledPlanCache &planCache, std::stop_token stop,
                             std::function<void(const BemProgress &)> progress);
} // namespace hornsculpt::cuda_detail
