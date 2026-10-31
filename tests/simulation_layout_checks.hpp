#pragma once

#include "simulation_dock.hpp"
#include <QApplication>
#include <QComboBox>
#include <QLabel>
#include <QMainWindow>
#include <QProgressBar>
#include <QPushButton>
#include <QScrollArea>
#include <QScrollBar>
#include <QTabWidget>
#include <stdexcept>

namespace hornsculpt::testing {

// Check the dock's actual layout and navigation without starting a GPU job.
inline void checkSimulationDockLayout() {
    const auto require = [](bool condition, const char *message) {
        if (!condition)
            throw std::runtime_error(message);
    };
    const auto settleLayout = [] {
        for (int i = 0; i < 3; ++i)
            QApplication::processEvents();
    };
    QMainWindow host;
    SimulationDock dock(&host);
    dock.setFloating(true);
    dock.resize(800, 450);
    dock.show();
    settleLayout();
    auto *workspace = dock.findChild<QTabWidget *>("simulationWorkspace");
    auto *results = dock.findChild<QTabWidget *>("acousticResults");
    auto *settings = dock.findChild<QScrollArea *>("simulationSettingsScroll");
    auto *start = dock.findChild<QPushButton *>("startSimulation");
    auto *stop = dock.findChild<QPushButton *>("stopSimulation");
    auto *status = dock.findChild<QLabel *>("simulationStatus");
    auto *progress = dock.findChild<QProgressBar *>("simulationProgress");
    auto *mode = dock.findChild<QComboBox *>("bemGeometryMode");
    auto *meridian = dock.findChild<QWidget *>("axisymmetricControls");
    auto *emptyResults = dock.findChild<QLabel *>("simulationEmptyResults");
    require(workspace && results && settings && start && stop && status && progress && mode &&
                meridian && emptyResults,
            "Simulation workspace or pinned job controls are missing.");
    require(workspace->currentIndex() == 0 && meridian->isHidden(),
            "A new 3D project must start with calculation settings and hide 2D-only controls.");
    const auto pinnedVisible = [&] {
        for (QWidget *control : {static_cast<QWidget *>(start), static_cast<QWidget *>(stop),
                                 static_cast<QWidget *>(status)})
            require(control->isVisible() && dock.rect().contains(QRect(
                                                control->mapTo(&dock, QPoint{}), control->size())),
                    "A calculation action or status has scrolled out of the dock.");
        require(progress->isHidden(), "Idle dock displayed a misleading progress indicator.");
    };
    settings->verticalScrollBar()->setValue(settings->verticalScrollBar()->maximum());
    pinnedVisible();
    mode->setCurrentIndex(1);
    settleLayout();
    require(!meridian->isHidden(), "2D mode did not reveal its meridian settings.");
    mode->setCurrentIndex(0);
    require(meridian->isHidden(), "Returning to 3D left unused meridian settings visible.");
    dock.setFrequencySweep(400, 4000, 17);
    const auto settingsBefore = dock.archive().settings;
    int changes = 0;
    dock.onChanged = [&] { ++changes; };
    dock.revealResultView();
    settleLayout();
    require(workspace->currentIndex() == 1,
            "Opening a result did not select the result workspace.");
    require(emptyResults->isVisible() && !emptyResults->text().isEmpty() && results->isHidden(),
            "A project without results must explain the next step instead of showing empty plots.");
    workspace->setCurrentIndex(0);
    dock.revealResultView();
    settleLayout();
    require(emptyResults->isVisible() && changes == 0 && !dock.isBusy(),
            "Navigating to an empty result page changed the project or started a calculation.");
    // Imported measurements are genuine display inputs even without a BEM result.
    auto measured = dock.archive();
    measured.measurements.push_back(
        {QStringLiteral("Layout-Prüfkurve"), true, 0, {{400, 80}, {4000, 82}}});
    dock.restoreArchive(measured);
    settleLayout();
    require(workspace->currentIndex() == 1 && emptyResults->isHidden() && results->isVisible() &&
                !dock.result() && !dock.isBusy() &&
                dock.archive().measurements == measured.measurements,
            "The empty state concealed an imported measurement or turned it into a BEM result.");
    for (const char *name : {"acousticOnAxis", "acousticPolarH", "acousticIsobarH"}) {
        auto *content = dock.findChild<QWidget *>(name);
        require(content, "Result page content missing.");
        int page = -1;
        for (int i = 0; i < results->count(); ++i)
            if (results->widget(i)->isAncestorOf(content)) { page = i; break; }
        require(page >= 0, "Result page not found.");
        results->setCurrentIndex(page);
        settleLayout();
        auto *scroll = qobject_cast<QScrollArea *>(results->currentWidget());
        require(scroll && scroll->horizontalScrollBar()->maximum() == 0 &&
                    scroll->verticalScrollBar()->maximum() == 0,
                "Standard H/V results do not fit in an 800 by 450 pixel dock.");
        require(dock.width() <= 800 && dock.height() <= 450,
                "Inactive calculation or field controls enlarged the result dock.");
        pinnedVisible();
        const auto shown = [&](const char *name) {
            auto *widget = dock.findChild<QWidget *>(name);
            require(widget != nullptr, "Expected context-specific result controls are missing.");
            return !widget->isHidden();
        };
        const bool polar = QString::fromLatin1(name) == "acousticPolarH";
        const bool response = QString::fromLatin1(name) == "acousticOnAxis";
        require(shown("polarFrequencyControls") == polar &&
                    shown("impedanceControls") == response &&
                    shown("polarNormalizationControls") == !response &&
                    shown("directivityScaleControls") == !response,
                "Result controls do not match the selected diagram type.");
    }
    workspace->setCurrentIndex(0);
    settleLayout();
    pinnedVisible();
    const auto settingsAfter = dock.archive().settings;
    require(changes == 0 && !dock.isBusy() &&
                settingsAfter.minimumFrequency == settingsBefore.minimumFrequency &&
                settingsAfter.maximumFrequency == settingsBefore.maximumFrequency &&
                settingsAfter.frequencyCount == settingsBefore.frequencyCount &&
                settingsAfter.geometryMode == settingsBefore.geometryMode,
            "Navigating settings/results modified the simulation request or started a job.");
    dock.onChanged = {};
    measured.measurements.clear();
    dock.restoreArchive(measured);
    dock.revealResultView();
    settleLayout();
    require(emptyResults->isVisible() && results->isHidden(),
            "Removing all display data did not restore the empty result explanation.");
}

} // namespace hornsculpt::testing
