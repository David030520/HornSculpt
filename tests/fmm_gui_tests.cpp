#include "fixture_geometry.hpp"
#include "simulation_dock.hpp"

#include <QApplication>
#include <QCheckBox>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QMainWindow>
#include <QSettings>
#include <QSpinBox>
#include <QTemporaryDir>
#include <QTimer>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
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

void verifyComputedArchive(const SimulationArchive &archive) {
    require(archive.current && archive.result && archive.mesh && archive.request,
            "Completed desktop solve did not produce a current, complete archive.");
    require(archive.settings.acceleration == BemAcceleration::fastMultipole &&
                archive.settings.fastMultipole.interpolationOrder == 4 &&
                archive.request->numerics.acceleration == BemAcceleration::fastMultipole &&
                archive.request->numerics.fastMultipole == archive.settings.fastMultipole &&
                archive.result->acceleration == BemAcceleration::fastMultipole,
            "Desktop Fast Multipole selection did not reach request and result metadata.");
    require(archive.result->triangleCount == 2048 && archive.result->frequencies.size() == 1 &&
                archive.result->frequencies.front().frequency == 500,
            "Desktop solve used an unexpected mesh or frequency.");
    const auto &frequency = archive.result->frequencies.front();
    require(frequency.excitations.size() == 1 && frequency.fastMultipole.farInteractionCount > 0 &&
                frequency.fastMultipole.representedFarPanelPairs > 0 &&
                frequency.fastMultipole.directPanelPairCount > 0 &&
                frequency.fastMultipole.directPanelPairCount +
                        frequency.fastMultipole.representedFarPanelPairs ==
                    2048u * 2048u &&
                frequency.residualCheckMilliseconds > 0,
            "Desktop solve did not execute genuine FMM interactions and direct residual audit.");
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
            auto *order = dock.findChild<QSpinBox *>("fastMultipoleOrder");
            require(enabled && order && !enabled->isChecked(),
                    "The explicit Fast Multipole switch is missing or unexpectedly enabled.");
            enabled->setChecked(true);
            order->setValue(4);
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
            verifyComputedArchive(archive);
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
        verifyComputedArchive(project.simulation);
        {
            QMainWindow window;
            SimulationDock dock(&window);
            dock.setMesh(std::make_shared<SurfaceMesh>(project.mesh->mesh), revision, velocities);
            dock.restoreArchive(project.simulation);
            require(dock.hasCurrentResult() && !dock.isBusy(),
                    "Reopening the completed FMM project lost its valid result or launched a job.");
            const auto reopened = dock.archive();
            verifyComputedArchive(reopened);
            const auto &frequency = reopened.result->frequencies.front();
            require(frequency.excitations.front().observerPressure.front() == computedOnAxis &&
                        frequency.excitations.front().iterations == computedIterations &&
                        frequency.fastMultipole.farInteractionCount == computedFarInteractions,
                    "Reopening changed computed pressures or FMM diagnostics.");
            const auto retained = dock.result();
            auto *enabled = dock.findChild<QCheckBox *>("fastMultipoleEnabled");
            auto *order = dock.findChild<QSpinBox *>("fastMultipoleOrder");
            require(enabled && order && enabled->isChecked() && order->value() == 4,
                    "Reopening did not restore the selected solver controls.");
            unsigned changes{};
            dock.onChanged = [&] { ++changes; };
            enabled->setChecked(false);
            require(changes == 1 && !dock.hasCurrentResult() && dock.result() == retained &&
                        !dock.isBusy() && !order->isEnabled(),
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
                  << computedIterations << " far_interactions=" << computedFarInteractions << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Fast Multipole desktop test: " << error.what() << '\n';
        return 1;
    }
}
