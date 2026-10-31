#include "acoustic_partition_fixture.hpp"
#include "fixture_geometry.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <hornsculpt/acoustic_partition.hpp>
#include <hornsculpt/bem.hpp>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

using namespace hornsculpt;
namespace {
using Complex = std::complex<double>;
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
struct Fixture {
    SurfaceMesh mesh;
    BemSubdomainModel model;
};
Fixture fixture(unsigned level) {
    Fixture f;
    f.mesh = testing::sphere(level);
    auto outer = testing::sphere(level);
    const auto nv = static_cast<std::uint32_t>(f.mesh.vertices.size());
    const auto nt = static_cast<std::uint32_t>(f.mesh.triangles.size());
    for (const auto v : outer.vertices)
        f.mesh.vertices.push_back({2 * v.x, 2 * v.y, 2 * v.z});
    for (auto t : outer.triangles) {
        for (auto &v : t.vertices)
            v += nv;
        t.faceId = 2;
        t.boundary = BoundaryCondition::wall;
        f.mesh.triangles.push_back(t);
    }
    BemSubdomain shell{1, "shell", BemDomainKind::interior, {}};
    BemSubdomain exterior{2, "exterior", BemDomainKind::exterior, {}};
    BemInterface iface{3, "interface", 1, 2, {}};
    for (std::uint32_t i = 0; i < nt; ++i)
        shell.panels.push_back({i, -1});
    for (std::uint32_t i = nt; i < f.mesh.triangles.size(); ++i) {
        shell.panels.push_back({i, 1});
        exterior.panels.push_back({i, -1});
        iface.triangleIndices.push_back(i);
    }
    f.model.domains = {shell, exterior};
    f.model.interfaces = {iface};
    return f;
}
BemRequest request(const BemSubdomainModel &model) {
    BemRequest r;
    r.subdomains = model;
    r.frequencies = {500};
    r.excitations = {{"unit", {{1, {1, 0}}}}};
    r.observers = {{0, 0, .3}, {1, 0, 0}};
    r.numerics.relativeResidualTolerance = 5e-5;
    return r;
}
Complex analytic(double radius) {
    const Complex i{0, 1};
    const double k = 2 * std::numbers::pi * 500 / 343.;
    return 1.21 * 343. * (i * k * .1) / (1. + i * k * .1) * (.1 / radius) *
           std::exp(-i * k * (radius - .1));
}
const BemSubdomainTrace &trace(const BemExcitationResult &r, std::uint32_t id) {
    auto it = std::find_if(r.subdomainTraces().begin(), r.subdomainTraces().end(),
                           [id](const auto &x) { return x.domainId == id; });
    require(it != r.subdomainTraces().end(), "Missing domain trace");
    return *it;
}
void checkInterface(const BemExcitationResult &r) {
    const auto &a = trace(r, 1), &b = trace(r, 2);
    require(a.interfaceFluxes.size() == 1 && b.interfaceFluxes.size() == 1,
            "Missing P1 interface flux");
    const auto &qa = a.interfaceFluxes.front(), &qb = b.interfaceFluxes.front();
    require(qa.interfaceId == 3 && qb.interfaceId == 3 && !qa.vertexIndices.empty() &&
                qa.vertexIndices.size() == qa.normalDerivative.size() &&
                qb.vertexIndices.size() == qb.normalDerivative.size() &&
                a.vertexIndices.size() == a.pressureTrace.size() &&
                b.vertexIndices.size() == b.pressureTrace.size(),
            "Invalid interface trace dimensions");
    double signal = 0;
    for (std::size_t j = 0; j < qa.vertexIndices.size(); ++j) {
        auto vertex = qa.vertexIndices[j];
        auto q = std::find(qb.vertexIndices.begin(), qb.vertexIndices.end(), vertex);
        auto pa = std::find(a.vertexIndices.begin(), a.vertexIndices.end(), vertex);
        auto pb = std::find(b.vertexIndices.begin(), b.vertexIndices.end(), vertex);
        require(q != qb.vertexIndices.end() && pa != a.vertexIndices.end() &&
                    pb != b.vertexIndices.end(),
                "Missing interface vertex");
        auto x = qa.normalDerivative[j],
             y = qb.normalDerivative[static_cast<std::size_t>(q - qb.vertexIndices.begin())];
        auto p = a.pressureTrace[static_cast<std::size_t>(pa - a.vertexIndices.begin())];
        auto other = b.pressureTrace[static_cast<std::size_t>(pb - b.vertexIndices.begin())];
        require(std::abs(x + y) <= 1e-8 * (1 + std::abs(x)), "Interface flux sign mismatch");
        require(std::abs(p - other) <= 1e-8 * (1 + std::abs(p)),
                "Interface pressure discontinuity");
        signal += std::norm(x);
    }
    require(signal > 1e-16, "Interface flux silently zero");
}
BemExcitationResult solve(const Fixture &f) {
    CudaBemSolver solver(f.mesh, f.model);
    auto out = solver.solve(request(f.model));
    require(out.frequencies.size() == 1 && out.frequencies.front().excitations.size() == 1,
            "Invalid result dimensions");
    auto r = out.frequencies.front().excitations.front();
    require(std::isfinite(r.boundaryResidual) && r.boundaryResidual <= 5e-5 &&
                r.observerPressure.size() == 2,
            "Residual or observer failure");
    checkInterface(r);
    return r;
}
void rejected(const Fixture &f, const char *message) {
    bool ok = false;
    try {
        CudaBemSolver solver(f.mesh, f.model);
        (void)solver.solve(request(f.model));
    } catch (const std::invalid_argument &) {
        ok = true;
    }
    require(ok, message);
}
void invalidModels() {
    const auto original = fixture(0);
    auto bad = original;
    bad.model.domains[0].panels[0].normalSign = 0;
    rejected(bad, "Accepted zero normal sign");
    bad = original;
    bad.model.domains[1].id = 1;
    rejected(bad, "Accepted duplicate domain ID");
    bad = original;
    bad.model.interfaces[0].secondDomainId = 999;
    rejected(bad, "Accepted unknown domain");
    bad = original;
    bad.model.domains[0].panels.pop_back();
    rejected(bad, "Accepted missing ownership");
    bad = original;
    bad.model.interfaces[0].triangleIndices[0] = 999999;
    rejected(bad, "Accepted unknown triangle");
    bad = original;
    bad.model.domains[1].panels[0].normalSign = 1;
    rejected(bad, "Accepted same interface normals");
    bad = original;
    auto third = bad.model.domains[1];
    third.id = 4;
    bad.model.domains.push_back(third);
    bad.model.interfaces.push_back({4, "overlap", 1, 4, bad.model.interfaces[0].triangleIndices});
    rejected(bad, "Accepted multi-domain interface junction");
}
void convergence() {
    std::vector<double> errors;
    for (unsigned level = 0; level <= 2; ++level) {
        auto f = fixture(level);
        auto s = solve(f);
        auto req = request(f.model);
        req.subdomains.reset();
        CudaBemSolver single(testing::sphere(level));
        auto out = single.solve(req);
        const auto &base = out.frequencies.front().excitations.front();
        double worst = 0;
        for (std::size_t i = 0; i < 2; ++i) {
            worst =
                std::max(worst, std::abs(s.observerPressure[i] / analytic(i == 0 ? .3 : 1.) - 1.));
            if (level == 2)
                require(std::abs(s.observerPressure[i] / base.observerPressure[i] - 1.) < .15,
                        "Split/unsplit mismatch");
        }
        errors.push_back(worst);
        std::cout << "Sphere level " << level << ", triangles " << f.mesh.triangles.size()
                  << ", error " << worst << ", residual " << s.boundaryResidual << '\n';
    }
    require(errors[1] < errors[0] && errors[2] < errors[1] && errors[2] < .15,
            "Sphere did not converge");
}
void orientationAndOrder() {
    auto f = fixture(1);
    const auto reference = solve(f);
    for (auto index : f.model.interfaces[0].triangleIndices) {
        std::swap(f.mesh.triangles[index].vertices[1], f.mesh.triangles[index].vertices[2]);
        for (auto &d : f.model.domains)
            for (auto &p : d.panels)
                if (p.triangleIndex == index)
                    p.normalSign = -p.normalSign;
    }
    std::swap(f.model.interfaces[0].firstDomainId, f.model.interfaces[0].secondDomainId);
    std::reverse(f.model.interfaces[0].triangleIndices.begin(),
                 f.model.interfaces[0].triangleIndices.end());
    for (auto &d : f.model.domains)
        std::reverse(d.panels.begin(), d.panels.end());
    std::reverse(f.model.domains.begin(), f.model.domains.end());
    auto changed = solve(f);
    for (std::size_t i = 0; i < reference.observerPressure.size(); ++i)
        require(std::abs(changed.observerPressure[i] / reference.observerPressure[i] - 1.) < 1e-5,
                "Orientation/order changed solution");
}
void cupSplitInvariance() {
    const auto physical = testing::acousticCup();
    MeshingOptions meshing;
    meshing.targetEdgeLength = .02;
    const auto divided = partitionAcousticSurface(
        physical, std::vector<AcousticRegionDefinition>{{7, "Hornluft", {2, 3}}}, meshing);
    BemRequest req;
    req.frequencies = {500};
    req.excitations = {{"Kolben", {{3, {1, 0}}}}};
    req.observers = {{0, 0, .2}, {.2, 0, .2}, {0, 0, .04}};
    CudaBemSolver original(physical);
    const auto baseline = original.solve(req);
    req.subdomains = divided.model;
    CudaBemSolver coupled(divided.mesh, divided.model);
    const auto output = coupled.solve(req);
    const auto &a = baseline.frequencies.front().excitations.front();
    const auto &b = output.frequencies.front().excitations.front();
    double worst = 0;
    for (std::size_t i = 0; i < req.observers.size(); ++i)
        worst = std::max(worst, std::abs(b.observerPressure[i] / a.observerPressure[i] - 1.));
    std::cout << "Cup interface-rim split comparison " << worst << ", residual "
              << b.boundaryResidual << '\n';
    require(worst < .05, "Cup with interface rim changed physical pressure by more than 5%");
    require(b.impedances.size() == 1 && b.impedances.front().faceId == 3,
            "Subdomain physical source impedance missing");
    auto segmented = physical;
    for (auto &t : segmented.triangles)
        if (t.faceId == 2) {
            double z = 0;
            for (auto vertex : t.vertices)
                z += segmented.vertices[vertex].z / 3;
            if (z > .04)
                t.faceId = 6;
        }
    const auto chain = partitionAcousticSurface(
        segmented,
        std::vector<AcousticRegionDefinition>{{7, "Hals", {2, 3}, {}}, {8, "Glocke", {6}, {}}},
        meshing);
    req.subdomains = chain.model;
    CudaBemSolver threeDomains(chain.mesh, chain.model);
    const auto threeResult = threeDomains.solve(req);
    const auto &threePressure =
        threeResult.frequencies.front().excitations.front().observerPressure;
    double chainDifference = 0;
    for (std::size_t i = 0; i < req.observers.size(); ++i)
        chainDifference =
            std::max(chainDifference, std::abs(threePressure[i] / a.observerPressure[i] - 1.));
    std::cout << "Three-domain horn split comparison " << chainDifference << '\n';
    require(chainDifference < .05,
            "Throat/bell/exterior decomposition changed physical pressure by more than 5%");
    const auto volume = std::abs(a.impedances.front().volumeVelocity);
    require(std::abs(b.impedances.front().volumeVelocity - a.impedances.front().volumeVelocity) <
                1e-10 * volume,
            "Subdomain decomposition changed physical driven volume flow");
}
void compareCacheResults(const BemResult &actual, const BemResult &reference, double tolerance) {
    require(actual.frequencies.size() == reference.frequencies.size(), "Cache frequency mismatch");
    const auto compare = [tolerance](Complex a, Complex b) {
        return std::isfinite(std::abs(a)) && std::isfinite(std::abs(b)) &&
               std::abs(a - b) <= tolerance * std::max(std::abs(b), 1e-20);
    };
    for (std::size_t i = 0; i < actual.frequencies.size(); ++i) {
        const auto &a = actual.frequencies[i].excitations.front();
        const auto &b = reference.frequencies[i].excitations.front();
        require(a.observerPressure.size() == b.observerPressure.size() &&
                    a.impedances.size() == b.impedances.size() && !a.impedances.empty(),
                "Cache result dimensions mismatch");
        for (std::size_t j = 0; j < a.observerPressure.size(); ++j)
            require(compare(a.observerPressure[j], b.observerPressure[j]),
                    "Cache changed observer pressure");
        for (std::size_t j = 0; j < a.impedances.size(); ++j)
            require(a.impedances[j].faceId == b.impedances[j].faceId &&
                        compare(a.impedances[j].volumeImpedance, b.impedances[j].volumeImpedance),
                    "Cache changed source impedance");
    }
}
void subdomainPlanCache() {
    const auto physical = testing::acousticCup();
    MeshingOptions meshing;
    meshing.targetEdgeLength = .02;
    const auto divided = partitionAcousticSurface(
        physical, std::vector<AcousticRegionDefinition>{{7, "Hornluft", {2, 3}}}, meshing);
    BemRequest req;
    req.subdomains = divided.model;
    req.frequencies = {350, 500};
    req.interactionPlanFrequencies = req.frequencies;
    req.excitations = {{"Kolben", {{3, {1, 0}}}}};
    req.observers = {{0, 0, .2}, {.2, 0, .2}, {0, 0, .04}};
    CudaBemSolver solver(divided.mesh, divided.model);
    const auto cold = solver.solve(req);
    const auto warm = solver.solve(req);
    require(!cold.frequencies.front().interactionPlanReused &&
                warm.frequencies.front().interactionPlanReused,
            "Subdomain plan cold/warm reuse mismatch");
    require(cold.frequencies.front().assemblyMilliseconds >=
                    cold.frequencies.front().interactionPlanMilliseconds &&
                warm.frequencies.front().interactionPlanMilliseconds == 0,
            "Subdomain planning time is missing from assembly or counted on cache reuse");
    compareCacheResults(warm, cold, 1e-9);

    auto subset = req;
    subset.frequencies = {350};
    const auto cachedSubset = solver.solve(subset);
    require(cachedSubset.frequencies.front().interactionPlanReused,
            "Subset did not reuse full-raster plan");
    auto independentRequest = subset;
    independentRequest.interactionPlanFrequencies = {350};
    CudaBemSolver independent(divided.mesh, divided.model);
    const auto independentSubset = independent.solve(independentRequest);
    compareCacheResults(cachedSubset, independentSubset, 1e-3);

    auto changedDensity = subset;
    changedDensity.medium.density *= 1.1;
    const auto densityResult = solver.solve(changedDensity);
    require(densityResult.frequencies.front().interactionPlanReused,
            "Density change unnecessarily invalidated geometry plan");
    CudaBemSolver densityReference(divided.mesh, divided.model);
    compareCacheResults(densityResult, densityReference.solve(changedDensity), 1e-9);

    auto changedOrder = subset;
    ++changedOrder.numerics.singularOrder;
    const auto orderResult = solver.solve(changedOrder);
    require(!orderResult.frequencies.front().interactionPlanReused,
            "Singular order change did not invalidate plan");
    CudaBemSolver orderReference(divided.mesh, divided.model);
    compareCacheResults(orderResult, orderReference.solve(changedOrder), 1e-9);

    auto largerPlan = changedOrder;
    largerPlan.interactionPlanFrequencies = {350, 500, 800};
    const auto largerResult = solver.solve(largerPlan);
    require(!largerResult.frequencies.front().interactionPlanReused,
            "Larger planning wavenumber did not invalidate plan");
    CudaBemSolver largerReference(divided.mesh, divided.model);
    compareCacheResults(largerResult, largerReference.solve(largerPlan), 1e-9);

    CudaBemSolver cancelled(divided.mesh, divided.model);
    std::stop_source stop;
    bool interrupted = false;
    try {
        (void)cancelled.solve(req, stop.get_token(), [&stop](const BemProgress &progress) {
            if (progress.phase == BemSolvePhase::preparing)
                stop.request_stop();
        });
    } catch (const std::runtime_error &) {
        interrupted = true;
    }
    require(interrupted, "Plan preparation cancellation was ignored");
    const auto afterCancel = cancelled.solve(req);
    require(!afterCancel.frequencies.front().interactionPlanReused,
            "Cancelled plan was installed as complete");
    compareCacheResults(afterCancel, cold, 1e-9);
    // Exercise the compact regular-pair rules at the same electrical panel
    // size as the high-frequency workload, on the identical partition.
    double maximumEdge = 0;
    for (const auto &triangle : divided.mesh.triangles) {
        for (unsigned edge = 0; edge < 3; ++edge) {
            const auto &a = divided.mesh.vertices[triangle.vertices[edge]];
            const auto &b = divided.mesh.vertices[triangle.vertices[(edge + 1) % 3]];
            maximumEdge = std::max(maximumEdge, std::hypot(a.x - b.x, a.y - b.y, a.z - b.z));
        }
    }
    require(maximumEdge > 0 && divided.mesh.triangles.size() <= 10000,
            "Invalid high-electrical-size fixture");
    const double scale = req.medium.soundSpeed / (2 * std::numbers::pi * maximumEdge);
    auto electrical = req;
    electrical.frequencies = {.75 * scale, 2.25 * scale};
    electrical.interactionPlanFrequencies = electrical.frequencies;
    CudaBemSolver electricalSolver(divided.mesh, divided.model);
    const auto electricalCold = electricalSolver.solve(electrical);
    const auto electricalWarm = electricalSolver.solve(electrical);
    require(!electricalCold.frequencies.front().interactionPlanReused &&
                electricalWarm.frequencies.front().interactionPlanReused,
            "High-electrical-size plan cold/warm reuse mismatch");
    const auto checkResiduals = [&electrical](const BemResult &result) {
        for (const auto &frame : result.frequencies)
            for (const auto &excitation : frame.excitations)
                require(std::isfinite(excitation.boundaryResidual) &&
                            excitation.boundaryResidual <=
                                electrical.numerics.relativeResidualTolerance,
                        "High-electrical-size residual exceeds target");
    };
    checkResiduals(electricalCold);
    checkResiduals(electricalWarm);
    compareCacheResults(electricalWarm, electricalCold, 1e-9);
    auto electricalSubset = electrical;
    electricalSubset.frequencies = {electrical.frequencies.front()};
    const auto electricalCached = electricalSolver.solve(electricalSubset);
    require(electricalCached.frequencies.front().interactionPlanReused,
            "High-electrical-size plan was not reused for subset");
    electricalSubset.interactionPlanFrequencies = electricalSubset.frequencies;
    CudaBemSolver electricalReference(divided.mesh, divided.model);
    const auto electricalFresh = electricalReference.solve(electricalSubset);
    checkResiduals(electricalCached);
    checkResiduals(electricalFresh);
    compareCacheResults(electricalCached, electricalFresh, 1e-3);
    auto strict = electrical;
    strict.numerics.nearRelativeTolerance = 1e-6;
    const auto strictResult = electricalSolver.solve(strict);
    require(!strictResult.frequencies.front().interactionPlanReused,
            "A stricter quadrature tolerance reused a looser plan");
    checkResiduals(strictResult);
    compareCacheResults(electricalCold, strictResult, 1e-3);
    std::cout << "high_electrical_max_edge_m=" << maximumEdge
              << " low_hz=" << electrical.frequencies.front()
              << " high_hz=" << electrical.frequencies.back() << " low_kh=0.75 high_kh=2.25\n";
    std::cout << "Subdomain plan cache tests PASS.\n";
}
void benchmark() {
    using Clock = std::chrono::steady_clock;
    const auto physical = testing::acousticCup();
    MeshingOptions meshing;
    meshing.targetEdgeLength = .02;
    const auto partitionStart = Clock::now();
    const auto divided = partitionAcousticSurface(
        physical, std::vector<AcousticRegionDefinition>{{7, "Hornluft", {2, 3}}}, meshing);
    const double partitionMs =
        std::chrono::duration<double, std::milli>(Clock::now() - partitionStart).count();
    require(divided.mesh.triangles.size() <= 10000, "Benchmark mesh exceeds 10000 triangles");
    BemRequest req;
    // This coarse fixture is validated for the 350 to 500 Hz benchmark band.
    // Higher frequencies require a separate mesh-convergence study.
    req.frequencies = {350, 500};
    req.excitations = {{"Kolben", {{3, {1, 0}}}}};
    req.observers = {{0, 0, .2}, {.2, 0, .2}, {0, 0, .04}};
    req.numerics.acceleration = BemAcceleration::dense;
    req.numerics.relativeResidualTolerance = 5e-5;
    CudaBemSolver original(physical);
    CudaBemSolver coupled(divided.mesh, divided.model);
    auto splitRequest = req;
    splitRequest.subdomains = divided.model;
    const auto timedSolve = [](CudaBemSolver &solver, const BemRequest &request) {
        const auto start = Clock::now();
        auto result = solver.solve(request);
        const double milliseconds =
            std::chrono::duration<double, std::milli>(Clock::now() - start).count();
        return std::make_pair(std::move(result), milliseconds);
    };
    const auto coldBaseline = timedSolve(original, req);
    const auto coldCoupled = timedSolve(coupled, splitRequest);
    const auto warmBaseline = timedSolve(original, req);
    const auto warmCoupled = timedSolve(coupled, splitRequest);
    const auto &baseline = warmBaseline.first;
    const auto &output = warmCoupled.first;
    require(baseline.frequencies.size() == req.frequencies.size() &&
                output.frequencies.size() == req.frequencies.size(),
            "Benchmark frequency count mismatch");
    std::cout << std::setprecision(10) << "benchmark=cup_split_dense device=\""
              << CudaBemSolver::deviceDescription() << "\"\n"
              << "physical_triangles=" << physical.triangles.size()
              << " partition_triangles=" << divided.mesh.triangles.size()
              << " physical_pressure_dofs=" << baseline.pressureDofCount
              << " partition_pressure_dofs=" << output.pressureDofCount
              << " partition_ms=" << partitionMs << '\n'
              << "quadrature_order=" << req.numerics.quadratureOrder
              << " singular_order=" << req.numerics.singularOrder
              << " near_distance_ratio=" << req.numerics.nearDistanceRatio
              << " near_relative_tolerance=" << req.numerics.nearRelativeTolerance
              << " near_maximum_depth=" << req.numerics.nearMaximumDepth
              << " residual_target=" << req.numerics.relativeResidualTolerance << '\n'
              << "unsplit_cold_wall_ms=" << coldBaseline.second
              << " split_cold_wall_ms=" << coldCoupled.second
              << " unsplit_warm_wall_ms=" << warmBaseline.second
              << " split_warm_wall_ms=" << warmCoupled.second
              << " unsplit_peak_device_bytes=" << baseline.peakDeviceBytes
              << " split_peak_device_bytes=" << output.peakDeviceBytes << '\n'
              << "timing_scope=each_solve_call_two_frequencies "
                 "order=unsplit_cold,split_cold,unsplit_warm,split_warm "
                 "phase_timings=warm "
                 "geometry_reuse=within_call_not_independently_instrumented\n"
              << "fixture_band_hz=350..500 higher_frequency_mesh_convergence=not_established\n";
    for (std::size_t i = 0; i < req.frequencies.size(); ++i) {
        const auto &aFrame = baseline.frequencies[i], &bFrame = output.frequencies[i];
        require(aFrame.excitations.size() == 1 && bFrame.excitations.size() == 1,
                "Benchmark excitation count mismatch");
        const auto &a = aFrame.excitations.front(), &b = bFrame.excitations.front();
        require(a.observerPressure.size() == req.observers.size() &&
                    b.observerPressure.size() == req.observers.size() && a.impedances.size() == 1 &&
                    b.impedances.size() == 1 && a.impedances.front().faceId == 3 &&
                    b.impedances.front().faceId == 3,
                "Benchmark physical results missing");
        double pressureDifference = 0;
        for (std::size_t j = 0; j < req.observers.size(); ++j) {
            require(std::abs(a.observerPressure[j]) > 0, "Benchmark reference pressure is zero");
            const double difference = std::abs(b.observerPressure[j] / a.observerPressure[j] - 1.);
            require(std::isfinite(difference), "Benchmark pressure is not finite");
            pressureDifference = std::max(pressureDifference, difference);
        }
        require(std::abs(a.impedances.front().volumeImpedance) > 0 &&
                    std::abs(a.impedances.front().volumeVelocity) > 0,
                "Benchmark reference impedance or volume flow is zero");
        const double impedanceDifference = std::abs(
            b.impedances.front().volumeImpedance / a.impedances.front().volumeImpedance - 1.);
        const double flowDifference = std::abs(
            b.impedances.front().volumeVelocity / a.impedances.front().volumeVelocity - 1.);
        std::cout << "frequency_hz=" << aFrame.frequency
                  << " pressure_relative_difference=" << pressureDifference
                  << " impedance_relative_difference=" << impedanceDifference
                  << " volume_flow_relative_difference=" << flowDifference
                  << " unsplit_residual=" << a.boundaryResidual
                  << " split_residual=" << b.boundaryResidual << '\n';
        const auto phases = [](std::string_view name, const BemFrequencyResult &f) {
            std::cout << "model=" << name << " frequency_hz=" << f.frequency
                      << " assembly_ms=" << f.assemblyMilliseconds
                      << " factorization_ms=" << f.factorizationMilliseconds
                      << " triangular_solve_ms=" << f.triangularSolveMilliseconds
                      << " residual_check_ms=" << f.residualCheckMilliseconds
                      << " field_ms=" << f.fieldMilliseconds << '\n';
        };
        phases("unsplit", aFrame);
        phases("split", bFrame);
        require(std::isfinite(a.boundaryResidual) && std::isfinite(b.boundaryResidual) &&
                    a.boundaryResidual <= req.numerics.relativeResidualTolerance &&
                    b.boundaryResidual <= req.numerics.relativeResidualTolerance,
                "Benchmark residual exceeds target");
        // The existing cup split test establishes a 5% pressure discretization tolerance.
        require(pressureDifference < .05 && impedanceDifference < .05 && flowDifference < 1e-10,
                "Benchmark physical split/unsplit agreement failed");
    }
}
} // namespace
int main(int argc, char **argv) {
    try {
        std::cout << "CUDA device: " << CudaBemSolver::deviceDescription() << '\n';
        if (argc == 2 && std::string_view(argv[1]) == "--cache-test") {
            subdomainPlanCache();
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--benchmark") {
            benchmark();
            return 0;
        }
        require(argc == 1, "Usage: subdomain_bem_tests [--benchmark|--cache-test]");
        invalidModels();
        convergence();
        orientationAndOrder();
        cupSplitInvariance();
        std::cout << "Subdomain tests PASS.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Subdomain BEM test failed: " << e.what() << '\n';
        return 1;
    }
}
