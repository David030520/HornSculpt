#include "branding.hpp"
#include "main_window.hpp"
#include "qt_localization.hpp"
#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QSettings>
#include <QTemporaryDir>
#include <QTimer>
#include <QTranslator>
#include <cstdio>
#include <exception>

int main(int argc, char **argv) {
    QApplication app(argc, argv);
    QTranslator qtTranslator;
    if (!hornsculpt::installGermanQtTranslations(qtTranslator)) {
        std::fprintf(stderr, "Deutscher Qt-Übersetzungskatalog konnte nicht geladen werden.\n");
        return 1;
    }
    QCoreApplication::setOrganizationName("HornSculpt");
    QCoreApplication::setApplicationName("HornSculpt");
    QCoreApplication::setApplicationVersion(QStringLiteral(HORNSCULPT_VERSION));
    const QDir packageRoot(QCoreApplication::applicationDirPath() + "/..");
    if (packageRoot.exists("portable.ini")) {
        QSettings::setDefaultFormat(QSettings::IniFormat);
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope,
                           packageRoot.filePath("data"));
    }
    QCommandLineParser parser;
    parser.setApplicationDescription("HornSculpt · Geometrie und Akustik");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addOption(
        {"smoke-test", "Prüft Live-Geometrie und Projekt-Roundtrip und beendet die Anwendung."});
    parser.addOption(
        {"mesh-test", "Erzeugt beim Smoke-Test zusätzlich ein geprüftes isotropes Netz."});
    parser.addOption({"solver-test", "Prüft zwei CUDA-Sweeps mit 7 Frequenzen von 400 bis 4000 Hz "
                                     "und vollständigen Ergebnis-Save/Open."});
    parser.addOption(
        {"solve-project",
         "Berechnet das geöffnete Projekt und speichert das Ergebnis unter einem neuen Namen.",
         "output"});
    parser.addOption({"project-test",
                      "Prüft den vollständigen Import einer angegebenen Projektdatei ohne Solve."});
    parser.addOption({"screenshot", "Speichert nach dem Aufbau ein Bild des Fensters.", "path"});
    parser.addOption({"screenshot-view",
                      "Ergebnisansicht: frequency, polar, isobar, angular, phase, field oder lem.",
                      "view"});
    parser.addOption(
        {"view",
         "Öffnet eine Ergebnisansicht: frequency, polar, isobar, angular, phase, field oder lem.",
         "view"});
    parser.addOption(
        {"evaluate-field",
         "Berechnet nach dem Öffnen den Feldschnitt aus gespeicherten BEM-Randwerten. "
         "Ein Screenshot wartet auf das Ergebnis; ohne Screenshot bleibt das Fenster offen."});
    parser.addPositionalArgument("project", "HornSculpt-Projektdatei (*.hsproj).", "[project]");
    parser.process(app);
    if (parser.isSet("project-test") && parser.positionalArguments().isEmpty())
        parser.showHelp(1);
    if (parser.isSet("project-test") &&
        (parser.isSet("smoke-test") || parser.isSet("mesh-test") || parser.isSet("solver-test")))
        parser.showHelp(1);
    const bool solveProject = parser.isSet("solve-project");
    if (solveProject) {
        const auto inputs = parser.positionalArguments();
        const auto output = parser.value("solve-project");
        if (inputs.size() != 1 || output.isEmpty() ||
            QFileInfo(inputs.front()).suffix().compare("hsproj", Qt::CaseInsensitive) != 0 ||
            QFileInfo(output).suffix().compare("hsproj", Qt::CaseInsensitive) != 0 ||
            parser.isSet("smoke-test") || parser.isSet("mesh-test") ||
            parser.isSet("solver-test") || parser.isSet("project-test") ||
            parser.isSet("screenshot") || parser.isSet("screenshot-view") || parser.isSet("view") ||
            parser.isSet("evaluate-field")) {
            std::fprintf(stderr, "--solve-project benötigt genau eine Eingabeprojektdatei und "
                                 "einen neuen .hsproj-Ausgabepfad ohne Testoptionen.\n");
            return 1;
        }
        const QFileInfo inputInfo(inputs.front()), outputInfo(output);
        if (!inputInfo.exists() || !inputInfo.isFile() || outputInfo.exists() ||
            inputInfo.absoluteFilePath().compare(outputInfo.absoluteFilePath(),
                                                 Qt::CaseInsensitive) == 0) {
            std::fprintf(stderr, "Eingabe fehlt oder Ausgabedatei existiert bereits.\n");
            return 1;
        }
    }
    hornsculpt::MainWindow window;
    hornsculpt::installBranding(window);
    if (!parser.positionalArguments().isEmpty()) {
        const auto path = parser.positionalArguments().front();
        const auto suffix = QFileInfo(path).suffix().toLower();
        if (suffix == "step" || suffix == "stp")
            window.openStep(path);
        else
            window.openProject(path);
    }
    window.show();
    if (solveProject) {
        auto *timer = new QTimer(&window);
        timer->setInterval(100);
        auto elapsed = std::make_shared<QElapsedTimer>();
        elapsed->start();
        auto stage = std::make_shared<int>(0);
        auto phaseStart = std::make_shared<qint64>(0);
        const QString output = QFileInfo(parser.value("solve-project")).absoluteFilePath();
        QObject::connect(
            timer, &QTimer::timeout, &window,
            [&app, &window, timer, elapsed, stage, phaseStart, output] {
                try {
                    if (window.sourceImportInProgress())
                        return;
                    if (!window.sourceImportError().isEmpty())
                        throw std::runtime_error(window.sourceImportError().toStdString());
                    if (!window.geometryReady())
                        throw std::runtime_error(
                            "Projektgeometrie ist nicht aktuell oder ungültig.");
                    if (*stage == 0) {
                        std::fprintf(stdout, "Projekt geladen; Revision %llu.\n",
                                     static_cast<unsigned long long>(window.revision()));
                        std::fflush(stdout);
                        if (!window.meshReady()) {
                            window.requestMesh();
                            *phaseStart = elapsed->elapsed();
                            *stage = 1;
                            return;
                        }
                        *stage = 2;
                    }
                    if (*stage == 1) {
                        if (window.meshInProgress())
                            return;
                        if (!window.meshReady()) {
                            const auto error = window.lastMeshError();
                            throw std::runtime_error(error.isEmpty()
                                                         ? "Netzerstellung fehlgeschlagen."
                                                         : error.toStdString());
                        }
                        std::fprintf(stdout, "Netz erstellt in %lld ms.\n",
                                     static_cast<long long>(elapsed->elapsed() - *phaseStart));
                        *stage = 2;
                    }
                    if (*stage == 2) {
                        const auto snapshot = window.projectState();
                        const auto &settings = snapshot.simulation.settings;
                        std::fprintf(stdout,
                                     "Solve: %zu Frequenzen, %.9g bis %.9g Hz, Beschleunigung %d, "
                                     "Geometriemodus %d.\n",
                                     settings.frequencyCount, settings.minimumFrequency,
                                     settings.maximumFrequency,
                                     static_cast<int>(settings.acceleration),
                                     static_cast<int>(settings.geometryMode));
                        std::fflush(stdout);
                        window.requestSolve();
                        if (!window.simulationInProgress()) {
                            const auto error = window.lastSimulationError();
                            throw std::runtime_error(
                                error.isEmpty() ? "Simulation konnte nicht gestartet werden."
                                                : error.toStdString());
                        }
                        *phaseStart = elapsed->elapsed();
                        *stage = 3;
                        return;
                    }
                    if (window.simulationInProgress())
                        return;
                    if (!window.lastSimulationError().isEmpty())
                        throw std::runtime_error(window.lastSimulationError().toStdString());
                    if (!window.simulationReady()) {
                        const auto error = window.lastSimulationError();
                        throw std::runtime_error(
                            error.isEmpty() ? "Simulation lieferte kein aktuelles Ergebnis."
                                            : error.toStdString());
                    }
                    const auto snapshot = window.projectState();
                    if (!snapshot.simulation.current || !snapshot.simulation.request ||
                        !snapshot.simulation.result)
                        throw std::runtime_error("Simulationssnapshot ist unvollständig.");
                    const auto &request = *snapshot.simulation.request;
                    const auto &result = *snapshot.simulation.result;
                    std::fprintf(stdout,
                                 "Ergebnis: %zu/%zu Frequenzen, %zu Dreiecke, %zu "
                                 "Druckfreiheitsgrade; Solve %lld ms.\n",
                                 result.frequencies.size(), request.frequencies.size(),
                                 result.triangleCount, result.pressureDofCount,
                                 static_cast<long long>(elapsed->elapsed() - *phaseStart));
                    for (const auto &f : result.frequencies)
                        std::fprintf(
                            stdout,
                            "%.9g Hz: Aufbau %.3f ms, Lösen %.3f ms, Faktorisierung %.3f ms, "
                            "Rückwärtssubstitution %.3f ms, Residuum %.3f ms, Feld %.3f ms.\n",
                            f.frequency, f.assemblyMilliseconds, f.solveMilliseconds,
                            f.factorizationMilliseconds, f.triangularSolveMilliseconds,
                            f.residualCheckMilliseconds, f.fieldMilliseconds);
                    if (result.frequencies.size() != request.frequencies.size())
                        throw std::runtime_error(
                            "Simulationsergebnis enthält nicht alle angeforderten Frequenzen.");
                    QTemporaryDir temporary(QFileInfo(output).absolutePath() +
                                            "/.hornsculpt-solve-XXXXXX");
                    if (!temporary.isValid())
                        throw std::runtime_error(temporary.errorString().toStdString());
                    const auto temporaryPath = temporary.filePath("result.hsproj");
                    hornsculpt::writeProject(temporaryPath, snapshot);
                    if (!QFile::rename(temporaryPath, output))
                        throw std::runtime_error("Ausgabe konnte nicht angelegt werden; vorhandene "
                                                 "Dateien werden nicht überschrieben.");
                    std::fprintf(stdout, "Projekt gespeichert: %s\n", qPrintable(output));
                    std::fflush(stdout);
                    timer->stop();
                    app.exit(0);
                } catch (const std::exception &e) {
                    std::fprintf(stderr, "Projektberechnung fehlgeschlagen: %s\n", e.what());
                    timer->stop();
                    app.exit(1);
                }
            });
        timer->start();
    }
    if (parser.isSet("smoke-test") || parser.isSet("screenshot") || parser.isSet("solver-test") ||
        parser.isSet("project-test") || parser.isSet("view") || parser.isSet("evaluate-field")) {
        auto *timer = new QTimer(&window);
        timer->setInterval(100);
        auto elapsed = std::make_shared<QElapsedTimer>();
        elapsed->start();
        auto stage = std::make_shared<int>(0);
        auto viewApplied = std::make_shared<bool>(false);
        auto fieldRequested = std::make_shared<bool>(false);
        const bool smoke = parser.isSet("smoke-test"), solve = parser.isSet("solver-test"),
                   mesh = parser.isSet("mesh-test") || solve;
        const bool project = parser.isSet("project-test");
        const bool evaluateField = parser.isSet("evaluate-field");
        const bool verification = smoke || mesh || solve || project || parser.isSet("screenshot");
        const QString screenshot = parser.value("screenshot");
        const QString screenshotView = parser.isSet("screenshot-view")
                                           ? parser.value("screenshot-view")
                                       : parser.isSet("view") ? parser.value("view")
                                       : evaluateField        ? QStringLiteral("field")
                                                              : QString{};
        QObject::connect(
            timer, &QTimer::timeout, &window,
            [&app, &window, timer, elapsed, stage, viewApplied, fieldRequested, smoke, mesh, solve,
             project, screenshot, screenshotView, verification, evaluateField] {
                if (elapsed->elapsed() > 120000) {
                    std::fprintf(stderr, "Desktop verification timed out.\n");
                    app.exit(2);
                    return;
                }
                try {
                    if (window.sourceImportInProgress())
                        return;
                    if (!window.sourceImportError().isEmpty())
                        throw std::runtime_error(window.sourceImportError().toStdString());
                    if (project && *stage == 0) {
                        window.exerciseLoadedProjectRoundTrip();
                        *stage = 6;
                    }
                    if (!project && !window.geometryReady())
                        return;
                    if (smoke && *stage == 0) {
                        window.exerciseProjectRoundTrip();
                        *stage = 1;
                        return;
                    }
                    if (mesh && *stage < 2) {
                        window.requestMesh();
                        *stage = 2;
                        return;
                    }
                    if (mesh && !window.meshReady()) {
                        if (!window.meshInProgress())
                            throw std::runtime_error(window.lastMeshError().toStdString());
                        return;
                    }
                    if (solve && *stage < 3) {
                        window.setSimulationFrequencySweep(400, 4000, 7);
                        window.requestSolve();
                        *stage = 3;
                        return;
                    }
                    if (solve && !window.simulationReady()) {
                        if (!window.simulationInProgress())
                            throw std::runtime_error(window.lastSimulationError().toStdString());
                        return;
                    }
                    if (solve && *stage == 3) {
                        window.requestSolve();
                        *stage = 4;
                        return;
                    }
                    if (solve && *stage == 4) {
                        window.exerciseResultProjectRoundTrip();
                        *stage = 5;
                    }
                    if (!*viewApplied) {
                        if (!screenshotView.isEmpty())
                            window.setSimulationView(screenshotView);
                        *viewApplied = true;
                        // Let the native dock layout settle before starting a
                        // field job or capturing its final presentation.
                        return;
                    }
                    if (evaluateField) {
                        if (!*fieldRequested) {
                            window.requestFieldEvaluation();
                            *fieldRequested = true;
                            elapsed->restart();
                            return;
                        }
                        if (window.fieldEvaluationInProgress())
                            return;
                        if (!window.fieldReady()) {
                            const auto error = window.lastFieldError();
                            throw std::runtime_error(error.isEmpty()
                                                         ? "Feldberechnung lieferte kein Ergebnis."
                                                         : error.toStdString());
                        }
                    }
                    if (!screenshot.isEmpty() && !window.grab().save(screenshot))
                        throw std::runtime_error("Screenshot konnte nicht gespeichert werden.");
                    timer->stop();
                    if (!verification)
                        return;
                    std::fprintf(stdout, "Desktop verification passed; revision %llu.\n",
                                 static_cast<unsigned long long>(window.revision()));
                    app.exit(0);
                } catch (const std::exception &e) {
                    std::fprintf(stderr, "Desktop verification failed: %s\n", e.what());
                    app.exit(1);
                }
            });
        timer->start();
    }
    return app.exec();
}
