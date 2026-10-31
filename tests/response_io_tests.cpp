#include "measured_response.hpp"
#include "result_export.hpp"
#include "simulation_request.hpp"
#include <QCoreApplication>
#include <QFile>
#include <QTemporaryDir>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
using namespace hornsculpt;
namespace {
void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
template <class F> void rejected(F f) {
    try {
        f();
    } catch (const std::invalid_argument &) {
        return;
    }
    throw std::runtime_error("Invalid response was accepted");
}
void observationFrames() {
    const auto nearVector = [](Vec3 actual, Vec3 expected) {
        require(std::hypot(actual.x - expected.x, actual.y - expected.y, actual.z - expected.z) <
                    1e-12,
                "Observation frame changed a world-space direction");
    };
    SimulationSettings settings;
    const std::vector<FaceVelocity> velocities{{5, {.5, .25}}};
    const auto original = makeSimulationRequest(settings, velocities, {1, 2, 3}, {0, 0, 2});
    nearVector(original.farField->origin, {1, 2, 3});
    nearVector(original.observers[0], {1, 2, 4});
    nearVector(original.observers[55], {2, 2, 3});
    nearVector(original.observers[128], {1, 3, 3});
    const auto alongX = makeSimulationRequest(settings, velocities, {1, 2, 3}, {2, 0, 0});
    nearVector(alongX.observers[0], {2, 2, 3});
    nearVector(alongX.observers[55], {1, 3, 3});
    nearVector(alongX.observers[128], {1, 2, 4});

    settings.observationFrame = AcousticFrame{{-.1, .4, -.2}, {0, -2, 0}, {0, 1, 3}};
    const auto request = makeSimulationRequest(settings, velocities, {1, 2, 3}, {2, 0, 0});
    require(request.observers.size() == 147 && request.farField->referenceDistance == 1,
            "Override changed the angular grid or pressure reference distance");
    nearVector(request.farField->origin, {-.1, .4, -.2});
    nearVector(request.observers[0], {-.1, -.6, -.2});
    for (std::size_t index = 0; index < 73; ++index) {
        const double angle = (-180.0 + 5.0 * static_cast<double>(index)) * std::numbers::pi / 180;
        nearVector(request.observers[index + 1],
                   {-.1, .4 - std::cos(angle), -.2 + std::sin(angle)});
        nearVector(request.observers[index + 74],
                   {-.1 - std::sin(angle), .4 - std::cos(angle), -.2});
    }
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    nearVector(simulationFrame(settings, {nan, 0, 0}, {}).origin, {-.1, .4, -.2});
    settings.observationFrame->horizontal = {0, -1, 0};
    rejected([&] { (void)makeSimulationRequest(settings, velocities); });
    settings.observationFrame.reset();
    nearVector(simulationFrame(settings, {1, 2, 3}, {0, 0, 2}).horizontal, {1, 0, 0});
    rejected([&] { (void)makeSimulationRequest(settings, velocities, {}, {}); });
    rejected([&] { (void)makeSimulationRequest(settings, velocities, {nan, 0, 0}); });
}

void coupledRequests() {
    SimulationSettings settings;
    settings.minimumFrequency = settings.maximumFrequency = 300;
    settings.frequencyCount = 1;
    settings.farField = false;
    settings.observationDistance = 2.5;
    settings.couplingEnabled = true;
    const lem::DriverParameters parameters{6.4, .0007, 4.98, .004, .0016, .262, .00544};
    settings.couplingNetwork = lem::makeDriverNetwork(parameters, {{5}, {7}}, {27, 3}, .0004);
    settings.velocityScale = 12;
    settings.phaseDegrees = 90;
    const auto request = makeSimulationRequest(settings, {}, {1, 2, 3}, {0, 0, 2});
    require(request.coupling == settings.couplingNetwork && request.excitations.empty(),
            "A coupled request lost its graph or retained prescribed excitation");
    require(!request.farField && request.frequencies == std::vector<double>{300},
            "A finite-distance coupled request silently became a far-field calculation");
    require(request.observers.size() == 147 && request.observers.front() == Vec3{1, 2, 5.5},
            "The finite observation radius or origin was not used");
    for (const auto observer : request.observers)
        require(std::abs(std::hypot(observer.x - 1, observer.y - 2, observer.z - 3) - 2.5) < 1e-12,
                "A polar observer does not lie on the requested finite radius");
    const auto voltage =
        std::get<lem::EffortSource>(request.coupling->elements.front().component).amplitude;
    require(voltage == lem::Complex{27, 3},
            "The prescribed-velocity multiplier changed the circuit voltage");
    std::get<lem::EffortSource>(settings.couplingNetwork.elements.front().component).amplitude = {
        11, 0};
    require(std::get<lem::EffortSource>(request.coupling->elements.front().component).amplitude ==
                voltage,
            "The simulation request does not own an immutable graph snapshot");
    settings.farField = true;
    const auto far = makeSimulationRequest(settings, {}, {1, 2, 3});
    require(far.farField == FarFieldObservation{{1, 2, 3}, 2.5},
            "The selected far-field reference distance was lost");
    for (const double badDistance : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                                     std::numeric_limits<double>::quiet_NaN()}) {
        settings.observationDistance = badDistance;
        rejected([&] { (void)makeSimulationRequest(settings, {}); });
    }
    settings.observationDistance = 1;
    settings.couplingNetwork.nodes.front().reference = false;
    rejected([&] { (void)makeSimulationRequest(settings, {}); });
    settings.couplingEnabled = false;
    rejected([&] { (void)makeSimulationRequest(settings, {}); });
}

void nearQuadratureRequests() {
    SimulationSettings settings;
    const std::vector<FaceVelocity> velocities{{5, {1, 0}}};
    require(makeSimulationRequest(settings, velocities).numerics.nearMaximumDepth == 6,
            "The desktop near quadrature default changed");
    for (const unsigned depth : {0u, 8u, 10u}) {
        settings.nearMaximumDepth = depth;
        const auto request = makeSimulationRequest(settings, velocities);
        require(request.numerics.nearMaximumDepth == depth &&
                    request.numerics.nearRelativeTolerance == 1e-4,
                "The request lost the subdivision depth or relaxed the near error target");
        settings.nearMaximumDepth = 6;
        require(request.numerics.nearMaximumDepth == depth,
                "Changing settings mutated an existing request snapshot");
    }
    settings.nearMaximumDepth = 11;
    rejected([&] { (void)makeSimulationRequest(settings, velocities); });
}
} // namespace
int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    try {
        observationFrames();
        coupledRequests();
        nearQuadratureRequests();
        QTemporaryDir dir;
        require(dir.isValid(), "Temporary directory unavailable");
        const auto path = dir.filePath("response.txt");
        const auto load = [&](const QByteArray &bytes) {
            QFile f(path);
            require(f.open(QIODevice::WriteOnly), "Cannot create response fixture");
            require(f.write(bytes) == bytes.size(), "Fixture write failed");
            f.close();
            return readMeasuredResponse(path);
        };
        const auto frd = load("# response\n100 82.5 -10\n1000 86 20\n");
        require(frd.samples.size() == 2 && frd.samples[1] == MeasuredPoint{1000, 86},
                "FRD import changed samples");
        const auto csv = load("frequency_Hz;level_dB\n100,5;82,5\n1000;86\n");
        require(csv.samples[0] == MeasuredPoint{100.5, 82.5}, "Decimal comma import failed");
        require(load("frequency_Hz,level_dB\n100,80\n200,81\n").samples.size() == 2,
                "CSV header failed");
        for (const auto *bad : {"200 80\n100 81\n", "100 80\n100 81\n", "0 80\n100 81\n",
                                "100 80\n200 nan\n", "100 80\n200 broken\n", "100 80\n"})
            rejected([&] { (void)load(bad); });
        SimulationSettings s;
        s.minimumFrequency = 100;
        s.maximumFrequency = 1000;
        s.frequencyCount = 3;
        const std::vector<FaceVelocity> velocities{{5, {.5, .25}}};
        auto request = makeSimulationRequest(s, velocities, {1, 2, 3}, {0, 0, 2});
        require(request.frequencies.size() == 3 &&
                    std::abs(request.frequencies[1] - std::sqrt(100000)) < 1e-10,
                "Request logarithmic frequency grid changed");
        require(request.observers.size() == 147 && request.observers[0].z == 4 &&
                    request.excitations[0].faces[0].normalVelocity == velocities[0].normalVelocity,
                "Desktop and study request convention differs");
        s.frequencyCount = 0;
        rejected([&] { (void)makeSimulationRequest(s, velocities); });
        SimulationArchive archive;
        archive.revision = 41;
        archive.request = request;
        archive.request->farField.reset();
        // Serialization fixture only; these values are not an acoustic reference.
        auto result = std::make_shared<BemResult>();
        result->device = "SYNTHETIC CSV TEST - no solver computation";
        BemExcitationResult excitation;
        excitation.name = "test;\"one\"";
        excitation.observerPressure.assign(request.observers.size(), {2, -1});
        excitation.impedances.push_back({5, .001, {}, {}, {100, 50}});
        result->frequencies.push_back({100, {excitation}});
        archive.result = result;
        const auto exported = resultCsv(archive);
        require(exported.contains("revision=41") && exported.contains("current=false") &&
                    exported.contains("observer;0;1;2;4;2;-1;") &&
                    exported.contains("field_mode=finite") &&
                    exported.contains("reference_distance_m=none") &&
                    exported.contains("pressure=Pa_at_observer") &&
                    exported.contains("test;\"\"one\"\"") && exported.contains("impedance;5;"),
                "CSV lost revision, complex pressure, escaping or normalized impedance");
        require(exported.count(";observer;") == static_cast<int>(request.observers.size()) &&
                    exported.count(";impedance;") == 1,
                "CSV did not export all genuine observations and face impedances");
        auto farField = archive;
        farField.request->farField = FarFieldObservation{{1, 2, 3}, 2.5};
        const auto farExport = resultCsv(farField);
        require(farExport.contains("field_mode=far_field") &&
                    farExport.contains("reference_distance_m=2.5") &&
                    farExport.contains("origin_m=1,2,3") &&
                    farExport.contains("world_direction_anchors_m") &&
                    farExport.contains("pressure=Pa_at_reference_distance") &&
                    farExport.contains("far_field;0;1;2;4;2;-1;") &&
                    !farExport.contains(";observer;"),
                "Far-field CSV confuses direction anchors with finite observation points");
        farField.request->farField->referenceDistance = 0;
        rejected([&] { (void)resultCsv(farField); });
        farField.request->farField = FarFieldObservation{farField.request->observers.front(), 1};
        rejected([&] { (void)resultCsv(farField); });
        rejected([&] { (void)resultCsv({}); });
        std::cout << "Measurement import, shared request and exact CSV export passed.\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
