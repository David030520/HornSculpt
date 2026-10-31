#include "fixture_geometry.hpp"
#include "simulation_dock.hpp"

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QMainWindow>
#include <QSettings>
#include <QSpinBox>
#include <QTemporaryDir>
#include <QTimer>
#include <charconv>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <thread>

using namespace hornsculpt;

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

template <class Predicate> void waitUntil(Predicate predicate) {
    QElapsedTimer elapsed;
    elapsed.start();
    while (!predicate()) {
        if (elapsed.elapsed() > 180000)
            throw std::runtime_error("Fast Multipole desktop solve timed out.");
        QCoreApplication::processEvents(QEventLoop::AllEvents, 10);
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
}

void verifyComputedArchive(const SimulationArchive &archive, unsigned expectedWidth) {
    require(archive.current && archive.result && archive.mesh && archive.request,
            "Completed desktop solve did not produce a current, complete archive.");
    require(archive.settings.singularOrder == 8 && archive.request->numerics.singularOrder == 8 &&
                archive.settings.acceleration == BemAcceleration::fastMultipole &&
                archive.settings.fastMultipole.leafSupportCapacity == 32 &&
                archive.settings.fastMultipole.preconditionerBlockWidth == expectedWidth &&
                archive.request->numerics.acceleration == BemAcceleration::fastMultipole &&
                archive.request->numerics.fastMultipole == archive.settings.fastMultipole &&
                archive.result->acceleration == BemAcceleration::fastMultipole,
            "Desktop Fast Multipole selection did not reach request and result metadata.");
    require(archive.result->triangleCount == 2048 && archive.result->frequencies.size() == 1 &&
                archive.result->frequencies.front().frequency == 500,
            "Desktop solve used an unexpected mesh or frequency.");
    const auto &frequency = archive.result->frequencies.front();
    require(
        frequency.excitations.size() == 1 && frequency.fastMultipole.farInteractionCount > 0 &&
            frequency.fastMultipole.pressurePressureFarPairs > 0 &&
            frequency.fastMultipole.pressurePressureNearPairs > 0 &&
            frequency.fastMultipole.pressurePressureNearPairs +
                    frequency.fastMultipole.pressurePressureFarPairs ==
                archive.mesh->vertices.size() * archive.mesh->vertices.size() &&
            frequency.fastMultipole.pressureFluxNearPairs +
                    frequency.fastMultipole.pressureFluxFarPairs ==
                archive.mesh->vertices.size() * archive.mesh->triangles.size() &&
            !frequency.fastMultipole.historicalChebyshevStatistics &&
            frequency.residualCheckMilliseconds > 0,
        "Desktop solve did not execute genuine FMM interactions and a reported residual audit.");
    const auto &excitation = frequency.excitations.front();
    require(
        excitation.iterations > 0 && std::isfinite(excitation.boundaryResidual) &&
            excitation.boundaryResidual <= archive.request->numerics.relativeResidualTolerance &&
            excitation.boundary()->surfacePressure.size() == archive.mesh->triangles.size() &&
            excitation.observerPressure.size() == archive.request->observers.size() &&
            !excitation.observerPressure.empty() &&
            std::isfinite(std::abs(excitation.observerPressure.front())) &&
            std::abs(excitation.observerPressure.front()) > 0 && excitation.impedances.size() == 1,
        "Desktop FMM pressure, impedance, iterations or residual are invalid.");
    if (excitation.residualEvaluation == BemResidualEvaluation::higherOrderFastMultipole)
        require(excitation.fastMultipoleVerificationAccuracySeed <
                        excitation.fastMultipoleComparisonAccuracySeed &&
                    excitation.fastMultipoleComparisonAccuracySeed <
                        excitation.fastMultipoleAccuracySeed &&
                    !excitation.fastMultipoleVerificationLevels.empty() &&
                    !excitation.historicalChebyshevEvidence &&
                    excitation.fastMultipoleOperatorDiscrepancy <=
                        .25 * archive.request->numerics.relativeResidualTolerance,
                "Desktop FMM result lost its approximate audit contract.");
    else
        require(excitation.fastMultipoleVerificationLevels.empty() &&
                    excitation.fastMultipoleComparisonLevels.empty() &&
                    excitation.fastMultipoleOperatorDiscrepancy == 0,
                "Desktop FMM direct audit contains approximate audit orders.");
}
} // namespace

int main(int argc, char **argv) {
    QApplication application(argc, argv);
    QCoreApplication::setOrganizationName("HornSculptTests");
    QCoreApplication::setApplicationName("FastMultipoleDesktop");
    QTemporaryDir directory;
    QSettings::setDefaultFormat(QSettings::IniFormat);
    QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, directory.path());
    try {
        unsigned selectedWidth = 256;
        if (argc != 1) {
            require(argc == 3 && std::string_view(argv[1]) == "--width",
                    "Usage: fmm_gui_tests [--width 256|512|1024|2048]");
            const std::string_view argument(argv[2]);
            const auto parsed =
                std::from_chars(argument.data(), argument.data() + argument.size(), selectedWidth);
            require(parsed.ec == std::errc{} && parsed.ptr == argument.data() + argument.size() &&
                        (selectedWidth == 256 || selectedWidth == 512 || selectedWidth == 1024 ||
                         selectedWidth == 2048),
                    "Invalid isolated desktop preconditioner width.");
        }
        require(directory.isValid(), "Temporary project directory cannot be created.");
        const auto path = directory.filePath("fmm-current.hsproj");
        const auto stalePath = directory.filePath("fmm-stale.hsproj");
        const auto mesh = std::make_shared<SurfaceMesh>(testing::sphere(4));
        const std::vector<FaceVelocity> velocities{{1, {.8, -.3}}};
        constexpr std::uint64_t revision = 17;
        std::complex<double> computedOnAxis;
        unsigned computedIterations{};
        std::size_t computedFarInteractions{};
        {
            QMainWindow window;
            SimulationDock dock(&window);
            dock.setMesh(mesh, revision, velocities);
            dock.setFrequencySweep(500, 500, 1);
            auto *enabled = dock.findChild<QCheckBox *>("fastMultipoleEnabled");
            auto *order = dock.findChild<QSpinBox *>("fastMultipoleLeafCapacity");
            auto *jacobi = dock.findChild<QComboBox *>("fastMultipoleJacobiBlock");
            auto *singular = dock.findChild<QSpinBox *>("singularQuadratureOrder");
            require(enabled && order && jacobi && singular && singular->value() == 4 &&
                        singular->isEnabled() && !enabled->isChecked() &&
                        jacobi->currentData().toUInt() == 32 && !jacobi->isEnabled(),
                    "The explicit Fast Multipole switch is missing or unexpectedly enabled.");
            enabled->setChecked(true);
            order->setValue(32);
            jacobi->setCurrentIndex(jacobi->findData(selectedWidth));
            singular->setValue(8);
            require(jacobi->isEnabled(), "Enabling FMM did not enable the Jacobi-block control.");
            unsigned eventLoopTicks{};
            QTimer heartbeat;
            heartbeat.setInterval(5);
            QObject::connect(&heartbeat, &QTimer::timeout, [&] { ++eventLoopTicks; });
            heartbeat.start();
            dock.requestSolve();
            require(dock.isBusy(), "The Fast Multipole desktop solve did not enter its worker.");
            waitUntil([&] { return !dock.isBusy(); });
            heartbeat.stop();
            if (!dock.lastError().isEmpty())
                throw std::runtime_error(dock.lastError().toStdString());
            require(eventLoopTicks > 0 && dock.hasCurrentResult(),
                    "The event loop stopped processing or the computed result was not published.");
            const auto archive = dock.archive();
            verifyComputedArchive(archive, selectedWidth);
            require(archive.revision == revision, "Computed result lost its snapshot revision.");
            const auto &frequency = archive.result->frequencies.front();
            computedOnAxis = frequency.excitations.front().observerPressure.front();
            computedIterations = frequency.excitations.front().iterations;
            computedFarInteractions = frequency.fastMultipole.farInteractionCount;
            Project project;
            project.revision = project.geometryRevision = project.meshRevision = revision;
            project.preview = mesh;
            project.mesh =
                std::make_shared<MeshingResult>(MeshingResult{*mesh, inspectMesh(*mesh)});
            project.simulation = archive;
            writeProject(path, project);
        }
        auto project = readProject(path);
        verifyComputedArchive(project.simulation, selectedWidth);
        {
            QMainWindow window;
            SimulationDock dock(&window);
            dock.setMesh(std::make_shared<SurfaceMesh>(project.mesh->mesh), revision, velocities);
            dock.restoreArchive(project.simulation);
            require(dock.hasCurrentResult() && !dock.isBusy(),
                    "Reopening the completed FMM project lost its valid result or launched a job.");
            const auto reopened = dock.archive();
            verifyComputedArchive(reopened, selectedWidth);
            const auto &frequency = reopened.result->frequencies.front();
            require(frequency.excitations.front().observerPressure.front() == computedOnAxis &&
                        frequency.excitations.front().iterations == computedIterations &&
                        frequency.fastMultipole.farInteractionCount == computedFarInteractions,
                    "Reopening changed computed pressures or FMM diagnostics.");
            const auto retained = dock.result();
            auto *enabled = dock.findChild<QCheckBox *>("fastMultipoleEnabled");
            auto *order = dock.findChild<QSpinBox *>("fastMultipoleLeafCapacity");
            auto *jacobi = dock.findChild<QComboBox *>("fastMultipoleJacobiBlock");
            auto *singular = dock.findChild<QSpinBox *>("singularQuadratureOrder");
            require(enabled && order && jacobi && singular && singular->value() == 8 &&
                        enabled->isChecked() && order->value() == 32 &&
                        jacobi->currentData().toUInt() == selectedWidth && jacobi->isEnabled(),
                    "Reopening did not restore the selected solver controls.");
            unsigned changes{};
            dock.onChanged = [&] { ++changes; };
            jacobi->setCurrentIndex(jacobi->findData(64u));
            require(
                changes == 1 && !dock.hasCurrentResult() && dock.result() == retained &&
                    dock.archive().settings.fastMultipole.preconditionerBlockWidth == 64 &&
                    dock.archive().request->numerics.fastMultipole.preconditionerBlockWidth ==
                        selectedWidth,
                "Changing Jacobi-block width did not invalidate and retain the solved snapshot.");
            jacobi->setCurrentIndex(jacobi->findData(selectedWidth));
            singular->setValue(10);
            require(changes == 3 && !dock.hasCurrentResult() && dock.result() == retained &&
                        dock.archive().settings.singularOrder == 10 &&
                        dock.archive().request->numerics.singularOrder == 8,
                    "Changing singular order relabeled the solved snapshot or kept it current.");
            singular->setValue(8);
            enabled->setChecked(false);
            require(changes == 5 && !dock.hasCurrentResult() && dock.result() == retained &&
                        !dock.isBusy() && !order->isEnabled() && !jacobi->isEnabled(),
                    "Switching to dense did not retain and invalidate the FMM snapshot.");
            project.simulation = dock.archive();
            require(project.simulation.settings.acceleration == BemAcceleration::dense &&
                        project.simulation.request->numerics.acceleration ==
                            BemAcceleration::fastMultipole &&
                        project.simulation.result->acceleration == BemAcceleration::fastMultipole &&
                        !project.simulation.current,
                    "The stale snapshot was relabeled with the newly selected solver mode.");
            writeProject(stalePath, project);
            const auto stale = readProject(stalePath);
            dock.restoreArchive(stale.simulation);
            require(!dock.hasCurrentResult() && !enabled->isChecked() && !dock.isBusy() &&
                        dock.result()->acceleration == BemAcceleration::fastMultipole &&
                        dock.result()
                                ->frequencies.front()
                                .excitations.front()
                                .observerPressure.front() == computedOnAxis,
                    "Stale FMM results were lost or marked current when reopening in dense mode.");
        }
        std::cout << "Actual CUDA FMM desktop solve, responsive event loop, current/stale project "
                     "save/reopen and backend switching passed. Iterations="
                  << computedIterations << " far_interactions=" << computedFarInteractions
                  << " block_width=" << selectedWidth << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Fast Multipole desktop test: " << error.what() << '\n';
        return 1;
    }
}
