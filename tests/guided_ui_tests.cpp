#include "main_window.hpp"
#include "parameter_sweep.hpp"
#include "project_file.hpp"
#include "result_comparison.hpp"
#include "simulation_dock.hpp"
#include "sweep_model.hpp"
#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QCryptographicHash>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QDirIterator>
#include <QElapsedTimer>
#include <QFileDialog>
#include <QMessageBox>
#include <QListWidget>
#include <QProcess>
#include <QLineEdit>
#include <QLabel>
#include <QProgressDialog>
#include <QProgressBar>
#include <QPushButton>
#include <QSpinBox>
#include <QSemaphore>
#include <QStandardPaths>
#include <QTabWidget>
#include <QTableView>
#include <QTemporaryDir>
#include <QThread>
#include <QTimer>
#include <QTest>
#include <QWizard>
#include <cstdio>
#include <numbers>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool value, const char *message) { if (!value) throw std::runtime_error(message); }
void waitUntil(const std::function<bool()> &predicate, int limit = 30000) {
    QElapsedTimer timer; timer.start();
    while (!predicate()) {
        if (timer.elapsed() > limit) throw std::runtime_error("Timed out waiting for GUI action");
        QCoreApplication::processEvents();
        QThread::msleep(5);
    }
}
QByteArray hashFile(const QString &path) {
    QFile file(path); require(file.open(QIODevice::ReadOnly), "Cannot read fixture");
    return QCryptographicHash::hash(file.readAll(), QCryptographicHash::Sha256);
}
void childWrite(const QString &root, const QString &original) {
    MainWindow window(nullptr, {}, {}, root);
    window.show();
    waitUntil([&] { return window.geometryReady() && window.findChild<QTimer *>("recoveryTimer"); });
    window.openProject(original);
    waitUntil([&] { return !window.sourceImportInProgress(); });
    require(window.sourceImportError().isEmpty(), "Original project did not open");
    auto *length = window.findChild<QDoubleSpinBox *>("length");
    require(length, "Length control missing");
    length->setValue(190);
    waitUntil([&] { return window.geometryReady(); });
    auto *timer = window.findChild<QTimer *>("recoveryTimer");
    require(QMetaObject::invokeMethod(timer, "timeout", Qt::DirectConnection), "Cannot request recovery");
    waitUntil([&] { return !QDir(root).entryList({"*.hsproj"}, QDir::Files).isEmpty(); });
    QElapsedTimer settling; settling.start();
    while (settling.elapsed() < 300) { QCoreApplication::processEvents(); QThread::msleep(5); }
    const auto recoveryFile = QDir(root).filePath(QDir(root).entryList({"*.hsproj"}, QDir::Files).front());
    const auto savedAt = QFileInfo(recoveryFile).lastModified();
    require(QMetaObject::invokeMethod(timer, "timeout", Qt::DirectConnection), "Cannot repeat recovery timer");
    settling.restart();
    while (settling.elapsed() < 300) { QCoreApplication::processEvents(); QThread::msleep(5); }
    require(QFileInfo(recoveryFile).lastModified() == savedAt,
            "Unchanged recovery snapshot was unnecessarily rewritten");
    std::puts("RECOVERY_READY"); std::fflush(stdout);
    QCoreApplication::exec();
}
void childRestore(const QString &root, const QString &output) {
    QTimer responder;
    QObject::connect(&responder, &QTimer::timeout, [] {
        if (auto *box = qobject_cast<QMessageBox *>(QApplication::activeModalWidget())) {
            for (auto *button : box->buttons())
                if (button->text() == QStringLiteral("Wiederherstellen")) { button->click(); break; }
        }
    });
    responder.start(20);
    MainWindow window(nullptr, {}, {}, root);
    window.show();
    auto *length = window.findChild<QDoubleSpinBox *>("length");
    waitUntil([&] { return length && length->value() == 190 && !window.sourceImportInProgress() &&
                            window.windowTitle().contains("*"); });
    require(window.sourceImportError().isEmpty(), "Recovery import failed");
    responder.stop();
    QTimer saveResponder;
    QObject::connect(&saveResponder, &QTimer::timeout, [&] {
        if (auto *dialog = qobject_cast<QFileDialog *>(QApplication::activeModalWidget())) {
            dialog->selectFile(output);
            static_cast<QDialog *>(dialog)->accept();
        }
    });
    saveResponder.start(20);
    QAction *save = nullptr;
    for (auto *action : window.findChildren<QAction *>())
        if (action->shortcut() == QKeySequence::Save) save = action;
    require(save, "Save action missing");
    save->trigger();
    waitUntil([&] { return QFileInfo::exists(output) && !window.saveInProgress(); });
    saveResponder.stop();
    require(!window.windowTitle().contains("*"), "Recovered project did not become clean after save");
    window.close();
}
void recoveryUi() {
    QTemporaryDir directory; require(directory.isValid(), "Temporary directory unavailable");
    const auto root = directory.filePath("recovery");
    const auto original = directory.filePath("original.hsproj");
    const auto output = directory.filePath("restored.hsproj");
    writeProject(original, Project{});
    const auto before = hashFile(original);
    QProcess process;
    process.setProcessChannelMode(QProcess::MergedChannels);
    process.start(QCoreApplication::applicationFilePath(), {"--write", root, original});
    require(process.waitForStarted(10000), "Recovery child did not start");
    QByteArray outputText;
    QElapsedTimer timer; timer.start();
    while (!outputText.contains("RECOVERY_READY") && timer.elapsed() < 30000) {
        process.waitForReadyRead(100);
        outputText += process.readAll();
        if (process.state() == QProcess::NotRunning) break;
    }
    require(outputText.contains("RECOVERY_READY"), outputText.constData());
    process.kill(); require(process.waitForFinished(10000), "Recovery child did not terminate");
    process.start(QCoreApplication::applicationFilePath(), {"--restore", root, output});
    require(process.waitForStarted(10000) && process.waitForFinished(45000), "Restore child timed out");
    const auto log = process.readAll();
    require(process.exitStatus() == QProcess::NormalExit && process.exitCode() == 0, log.constData());
    require(hashFile(original) == before, "Recovery changed the original project");
    require(std::abs(readProject(output).parameters.length - .19) < 1e-9,
            "Saved recovery lost edited dimensions");
    std::puts("Real GUI recovery after hard process termination passed.");
}
void guidedCalculation() {
    MainWindow window;
    window.show();
    waitUntil([&] { return window.geometryReady(); });
    auto *dock = dynamic_cast<SimulationDock *>(window.findChild<QDockWidget *>("simulationDock"));
    auto *length = window.findChild<QDoubleSpinBox *>("length");
    auto *meshPreset = window.findChild<QComboBox *>("meshPreset");
    require(dock && length && meshPreset, "Guided baseline controls missing");
    const auto initialSettings = dock->archive().settings;
    const double initialLength = length->value();
    const int initialPreset = meshPreset->currentIndex();
    auto *action = window.findChild<QAction *>("openDraftWizard");
    require(action, "Wizard action missing");
    QTimer::singleShot(0, [&] {
        auto *wizard = window.findChild<QWizard *>("draftWizard");
        if (!wizard) return;
        wizard->findChild<QDoubleSpinBox *>("wizardLength")->setValue(60);
        wizard->button(QWizard::NextButton)->click();
        wizard->findChild<QDoubleSpinBox *>("wizardHorizontalOpening")->setValue(30);
        wizard->findChild<QDoubleSpinBox *>("wizardVerticalOpening")->setValue(30);
        wizard->button(QWizard::NextButton)->click();
        wizard->button(QWizard::NextButton)->click();
        wizard->findChild<QComboBox *>("wizardMeshPreset")->setCurrentIndex(0);
        wizard->findChild<QDoubleSpinBox *>("wizardMinimumFrequency")->setValue(1000);
        wizard->findChild<QDoubleSpinBox *>("wizardMaximumFrequency")->setValue(1000);
        wizard->findChild<QSpinBox *>("wizardFrequencyCount")->setValue(1);
        wizard->button(QWizard::FinishButton)->click();
    });
    action->trigger();
    waitUntil([&] { return window.findChild<QProgressDialog *>("guidedCalculationProgress") != nullptr; });
    waitUntil([&] { return !window.findChild<QProgressDialog *>("guidedCalculationProgress") ||
                            window.simulationReady(); }, 180000);
    require(window.meshReady() && window.simulationReady() && dock && dock->result() &&
                !dock->result()->device.empty(), "Guided calculation did not produce an actual GPU result");
    std::printf("Guided calculation produced a result on %s.\n", dock->result()->device.c_str());
    auto *undo = window.findChild<QAction *>("undoEdit");
    require(undo && undo->isEnabled(), "Guided edit cannot be undone");
    undo->trigger();
    waitUntil([&] { return window.geometryReady(); });
    const auto restored = dock->archive().settings;
    require(length->value() == initialLength && meshPreset->currentIndex() == initialPreset &&
                restored.minimumFrequency == initialSettings.minimumFrequency &&
                restored.maximumFrequency == initialSettings.maximumFrequency &&
                restored.frequencyCount == initialSettings.frequencyCount,
            "One undo did not restore the guided geometry, mesh and frequency settings");
}
void resultFocus() {
    MainWindow window;
    window.show();
    waitUntil([&] { return window.geometryReady(); });
    auto *dock = dynamic_cast<SimulationDock *>(window.findChild<QDockWidget *>("simulationDock"));
    auto *enlarge = window.findChild<QAction *>("enlargeResults");
    auto *tabs = window.findChild<QTabWidget *>("acousticResults");
    require(dock && enlarge && tabs && dock->onFieldChanged, "Result controls missing");
    window.setSimulationView("frequency");
    for (int i = 0; i < tabs->count(); ++i)
        if (tabs->widget(i)->findChild<QWidget *>("cta2034Levels")) tabs->setCurrentIndex(i);
    enlarge->setChecked(true);
    auto plane = std::make_shared<PressureFieldPlane>();
    plane->columns = plane->rows = 2;
    plane->minimumLateral = plane->minimumAxial = -1;
    plane->maximumLateral = plane->maximumAxial = 1;
    plane->frequency = 1000;
    plane->positions = {{-1, 0, -1}, {1, 0, -1}, {-1, 0, 1}, {1, 0, 1}};
    plane->pressure.assign(4, {1, 0});
    plane->valid.assign(4, true);
    dock->onFieldChanged(plane, {});
    require(enlarge->isChecked(), "Field publication stole the enlarged CTA view");
    window.setSimulationView("field");
    require(!enlarge->isChecked(), "Selecting the field did not restore 3D");
    std::puts("CTA focus is preserved across field publication.");
}
void lateGuidedCancel() {
    QTemporaryDir directory;
    require(directory.isValid(), "Temporary directory unavailable");
    const auto projectPath = directory.filePath("replacement.hsproj");
    writeProject(projectPath, Project{});
    MainWindow window;
    window.show();
    waitUntil([&] { return window.geometryReady(); });
    auto *action = window.findChild<QAction *>("openDraftWizard");
    require(action, "Wizard action missing");
    QTimer::singleShot(0, [&] {
        auto *wizard = window.findChild<QWizard *>("draftWizard");
        if (!wizard) return;
        wizard->findChild<QDoubleSpinBox *>("wizardLength")->setValue(60);
        wizard->button(QWizard::NextButton)->click();
        wizard->button(QWizard::NextButton)->click();
        wizard->button(QWizard::NextButton)->click();
        wizard->button(QWizard::FinishButton)->click();
    });
    action->trigger();
    waitUntil([&] { return window.meshInProgress(); });
    auto *progress = window.findChild<QProgressDialog *>("guidedCalculationProgress");
    require(progress, "Guided progress missing");
    auto *timer = progress->findChild<QTimer *>("guidedCalculationTimer");
    require(timer, "Guided timer missing");
    timer->stop(); // Hold the old callback across the project change deterministically.
    auto *cancelMesh = window.findChild<QAction *>("cancelMeshing");
    require(cancelMesh && cancelMesh->isEnabled(), "Guided mesh cancellation missing");
    cancelMesh->trigger();
    waitUntil([&] { return !window.meshInProgress(); });
    window.openProject(projectPath);
    waitUntil([&] { return !window.sourceImportInProgress(); });
    if (!window.sourceImportError().isEmpty())
        std::fprintf(stderr, "Replacement import: %s\n", window.sourceImportError().toUtf8().constData());
    require(window.sourceImportError().isEmpty() && window.geometryReady(), "Replacement project failed");
    window.requestMesh();
    require(window.meshInProgress(), "Replacement mesh did not start");
    require(QMetaObject::invokeMethod(progress, "canceled", Qt::DirectConnection),
            "Could not deliver the old cancellation");
    waitUntil([&] { return !window.meshInProgress(); });
    require(window.meshReady(), "A late guided cancellation stopped the replacement mesh");
    std::puts("Late guided cancellation preserved the replacement mesh.");
}
void rejectedRecovery() {
    QTemporaryDir directory;
    require(directory.isValid(), "Temporary directory unavailable");
    {
        RecoveryStore store(directory.path());
        require(store.save(Project{}, "original.hsproj", {}), "Recovery fixture did not start");
        waitUntil([&] { return !store.busy(); });
    }
    QSemaphore entered, release;
    MainWindow window(nullptr, {}, [&](std::function<void(std::stop_token)> job) {
        return std::jthread([&, job = std::move(job)](std::stop_token stop) {
            entered.release();
            release.acquire();
            job(stop);
        });
    }, directory.path());
    struct ReleaseOnExit { QSemaphore &semaphore; ~ReleaseOnExit() { semaphore.release(); } } guard{release};
    QTimer responder;
    QObject::connect(&responder, &QTimer::timeout, [] {
        if (auto *box = qobject_cast<QMessageBox *>(QApplication::activeModalWidget()))
            for (auto *button : box->buttons())
                if (button->text() == QStringLiteral("Wiederherstellen")) { button->click(); break; }
    });
    responder.start(20);
    window.show();
    waitUntil([&] { return entered.tryAcquire(); });
    responder.stop();
    auto *frequency = window.findChild<QDoubleSpinBox *>("simulationMinimumFrequency");
    require(frequency, "Simulation frequency control missing");
    frequency->setValue(frequency->value() + 1);
    release.release();
    waitUntil([&] { return !window.sourceImportInProgress(); });
    require(window.sourceImportError().contains("Import verworfen"), "Changed import was not rejected");
    auto retained = RecoveryStore::available(directory.path());
    require(retained.size() == 1 && QFileInfo::exists(retained.front().path),
            "Rejected recovery retained its lock or deleted the file");
    RecoveryStore::discard(retained.front());
    std::puts("Rejected recovery released its lock and preserved the file.");
}
void comparisonEpoch() {
    QTemporaryDir directory;
    require(directory.isValid(), "Temporary directory unavailable");
    const auto path = directory.filePath("replacement.hsproj");
    writeProject(path, Project{});
    MainWindow window;
    window.show();
    waitUntil([&] { return window.geometryReady(); });
    auto *action = window.findChild<QAction *>("openResultComparison");
    require(action, "Comparison action missing");
    QString failure;
    QTimer::singleShot(0, [&] {
        try {
            ResultComparisonDialog *dialog = nullptr;
            for (auto *candidate : window.findChildren<QDialog *>())
                if (auto *comparison = dynamic_cast<ResultComparisonDialog *>(candidate)) dialog = comparison;
            require(dialog && dialog->onSnapshotsChanged, "Comparison callback missing");
            window.openProject(path);
            waitUntil([&] { return !window.sourceImportInProgress(); });
            require(window.sourceImportError().isEmpty(), "Replacement project failed");
            dialog->onSnapshotsChanged({{QStringLiteral("Old project"), std::make_shared<const Project>()}});
            dialog->accept();
        } catch (const std::exception &error) {
            failure = QString::fromUtf8(error.what());
            if (auto *dialog = qobject_cast<QDialog *>(QApplication::activeModalWidget())) dialog->reject();
        }
    });
    action->trigger();
    require(failure.isEmpty(), failure.toUtf8().constData());
    bool noOldSnapshots = false;
    QTimer::singleShot(0, [&] {
        for (auto *candidate : window.findChildren<QDialog *>())
            if (auto *dialog = dynamic_cast<ResultComparisonDialog *>(candidate)) {
                auto *list = dialog->findChild<QListWidget *>();
                noOldSnapshots = list && list->count() == 0;
                dialog->accept();
            }
    });
    action->trigger();
    require(noOldSnapshots, "Old comparison dialog overwrote the replacement project");
    std::puts("Old comparison dialog cannot mutate a replacement project.");
}
void studyEdit(bool afterCuda) {
    QTemporaryDir output;
    require(output.isValid(), "Study output directory unavailable");
    MainWindow window;
    window.show();
    waitUntil([&] { return window.geometryReady(); });
    auto *dock = dynamic_cast<SimulationDock *>(window.findChild<QDockWidget *>("simulationDock"));
    auto *action = window.findChild<QAction *>("openParameterSweep");
    require(dock && action, "Study controls missing");
    if (afterCuda) {
        window.requestMesh();
        waitUntil([&] { return !window.meshInProgress(); }, 60000);
        require(window.meshReady(), "Study baseline mesh failed");
        dock->setFrequencySweep(1000, 1000, 1);
        window.requestSolve();
        waitUntil([&] { return !window.simulationInProgress(); }, 120000);
        require(window.simulationReady() && dock->result() && !dock->result()->device.empty(),
                "Study baseline CUDA solve failed");
    } else {
        dock->setRealtimeEnabled(true);
    }
    const auto baseline = dock->result();
    const bool gpuReleasedBefore = dock->externalJobReady();
    QString failure;
    QString analyzedFile;
    QByteArray analyzedHash;
    double analyzedAngle = 0;
    QTimer unexpectedDialog;
    QObject::connect(&unexpectedDialog, &QTimer::timeout, [&] {
        if (auto *message = qobject_cast<QMessageBox *>(QApplication::activeModalWidget())) {
            failure = "Unexpected study message: " + message->text();
            message->reject();
        }
    });
    unexpectedDialog.start(20);
    QTimer::singleShot(0, [&] {
        ParameterSweepDialog *study = nullptr;
        try {
            study = dynamic_cast<ParameterSweepDialog *>(QApplication::activeModalWidget());
            require(study, "Opening an idle study displayed a GPU wait instead of the editor");
            require(afterCuda || dock->realtimeEnabled(), "Opening a study stopped realtime");
            auto *first = study->findChild<QComboBox *>("sweepParameter1");
            auto *second = study->findChild<QCheckBox *>("sweepSecondAxis");
            auto *start = study->findChild<QDoubleSpinBox *>("sweepStart1");
            auto *directory = study->findChild<QLineEdit *>("sweepOutputDirectory");
            auto *startButton = study->findChild<QPushButton *>("sweepStart");
            require(first && second && start && directory && startButton, "Study editor missing controls");
            directory->setText(output.path());
            first->setCurrentIndex(first->findData("horizontalAngle"));
            second->setChecked(true);
            auto *editor = start->findChild<QLineEdit *>();
            require(editor, "Study numeric editor missing");
            editor->setFocus();
            QTest::keyClick(editor, Qt::Key_A, Qt::ControlModifier);
            QTest::keyClicks(editor, "45");
            QTest::keyClick(editor, Qt::Key_Return);
            QTest::keyClick(editor, Qt::Key_Tab);
            QCoreApplication::processEvents();
            require(std::abs(start->value() - 45) < 1e-9, "Keyboard input was not committed");
            require(startButton->isEnabled() && QDir(output.path()).entryList(QDir::Dirs | QDir::NoDotAndDotDot).empty(),
                    "Editing a parameter started a study");
            require(dock->externalJobReady() == gpuReleasedBefore && dock->result() == baseline &&
                        (afterCuda || dock->realtimeEnabled()),
                    "Editing a study claimed GPU resources or changed desktop state");
            if (afterCuda) {
                second->setChecked(false);
                auto *points = study->findChild<QSpinBox *>("sweepPoints1");
                auto *end = study->findChild<QDoubleSpinBox *>("sweepEnd1");
                auto *progress = study->findChild<QProgressBar *>("sweepProgress");
                auto *status = study->findChild<QLabel *>("sweepStatus");
                require(points && end && progress && status, "Study progress controls missing");
                points->setValue(1);
                end->setValue(start->value());
                window.requestSolve();
                require(window.simulationInProgress(), "Desktop CUDA solve did not restart");
                startButton->click();
                require(!startButton->isEnabled() && status->text().contains("Berechnung"),
                        "Active desktop solve was not handed off at study start");
                waitUntil([&] { return startButton->isEnabled(); }, 180000);
                QDirIterator manifests(output.path(), {"study.hsstudy"}, QDir::Files,
                                       QDirIterator::Subdirectories);
                require(manifests.hasNext(), "Completed study manifest missing");
                const auto studyData = readSweepStudy(manifests.next());
                require(studyData.variants.size() == 1 &&
                            studyData.variants.front().state == SweepState::completed,
                        "Real CUDA study did not complete after handoff");
                const auto variant = readProject(QDir(studyData.directory).filePath(
                    studyData.variants.front().projectFile));
                require(variant.simulation.result && baseline &&
                            variant.simulation.current && variant.simulation.resultComplete == true &&
                            variant.simulation.result->device == baseline->device &&
                            !variant.simulation.result->device.empty(),
                        "Study variant did not retain the actual CUDA device result");
                require(progress->value() == 1, "Study progress did not finish");
                window.requestSolve();
                require(window.simulationInProgress(), "Desktop solve stayed locked after study");
                waitUntil([&] { return !window.simulationInProgress(); }, 120000);
                require(window.simulationReady(), "Desktop solve failed after study GPU release");
                analyzedFile = QDir(studyData.directory).filePath(studyData.variants.front().projectFile);
                analyzedHash = hashFile(analyzedFile);
                analyzedAngle = variant.parameters.horizontalAngle;
                auto *table = study->findChild<QTableView *>("sweepVariants");
                auto *analyze = study->findChild<QPushButton *>("sweepAnalyzeVariant");
                require(table && analyze, "Study analysis controls missing");
                table->selectRow(0);
                require(analyze->isEnabled(), "Completed result cannot be selected for analysis");
                unexpectedDialog.stop();
                QTimer::singleShot(0, [&] {
                    if (auto *box = qobject_cast<QMessageBox *>(QApplication::activeModalWidget()))
                        if (auto *discard = box->button(QMessageBox::Discard)) discard->click();
                });
                analyze->click();
            } else {
                study->reject();
            }
        } catch (const std::exception &error) {
            if (failure.isEmpty()) failure = QString::fromUtf8(error.what());
            if (auto *dialog = qobject_cast<QDialog *>(QApplication::activeModalWidget())) dialog->reject();
            if (study && study->isVisible()) study->reject();
        }
    });
    action->trigger();
    if (afterCuda && failure.isEmpty()) {
        waitUntil([&] { return !window.sourceImportInProgress(); }, 60000);
        require(window.sourceImportError().isEmpty() && window.meshReady() &&
                    window.simulationReady() && dock->hasCurrentResult() && !window.simulationInProgress() &&
                    window.windowTitle().contains('*') && hashFile(analyzedFile) == analyzedHash,
                "Study analysis did not load the current result as an unchanged copy");
        auto *angle = window.findChild<QDoubleSpinBox *>("angleH");
        require(angle && std::abs(angle->value() - analyzedAngle * 180 / std::numbers::pi) < .01,
                "Study analysis loaded the wrong geometry parameters");
    }
    if (!afterCuda && failure.isEmpty()) {
        QTimer::singleShot(0, [&] {
            auto *study = dynamic_cast<ParameterSweepDialog *>(QApplication::activeModalWidget());
            if (!study || !dock->realtimeEnabled()) {
                failure = "Reopening the idle study claimed the GPU or stopped realtime";
                if (auto *dialog = qobject_cast<QDialog *>(QApplication::activeModalWidget())) dialog->reject();
                return;
            }
            study->reject();
        });
        action->trigger();
    }
    unexpectedDialog.stop();
    dock->setRealtimeEnabled(false);
    require(failure.isEmpty(), failure.toUtf8().constData());
    std::puts(afterCuda ? "Study editing preserved the completed CUDA cache and result."
                        : "Study editing did not claim GPU or stop realtime.");
}
void studyHandoffCancel() {
    QTemporaryDir output;
    require(output.isValid(), "Study output directory unavailable");
    auto study = std::make_unique<ParameterSweepDialog>(Project{});
    unsigned acquired = 0, cancelled = 0, released = 0;
    std::function<void(bool)> pending;
    study->onAcquireGpu = [&](std::function<void(bool)> done) { ++acquired; pending = std::move(done); };
    study->onCancelGpuAcquire = [&] { ++cancelled; };
    study->onReleaseGpu = [&] { ++released; };
    study->show();
    auto *directory = study->findChild<QLineEdit *>("sweepOutputDirectory");
    auto *start = study->findChild<QPushButton *>("sweepStart");
    auto *cancel = study->findChild<QPushButton *>("sweepCancel");
    require(directory && start && cancel, "Study handoff controls missing");
    directory->setText(output.path());
    start->click();
    require(acquired == 1 && pending && cancel->isEnabled() && !start->isEnabled(),
            "Study did not enter cancellable GPU handoff");
    auto stale = pending;
    cancel->click();
    require(cancelled == 1 && start->isEnabled() && released == 0,
            "Cancelling GPU handoff failed to release its reservation");
    start->click();
    require(acquired == 2 && !start->isEnabled(), "Study did not retry GPU handoff");
    stale(true);
    require(!start->isEnabled() && cancelled == 1 && released == 0,
            "An old GPU handoff callback started a newer study");
    auto late = pending;
    study->reject();
    require(cancelled == 2, "Closing a waiting study left GPU handoff reserved");
    study.reset();
    late(true); // A provider completion may arrive after the dialog is gone.
    require(released == 0, "A cancelled handoff released an unowned worker lease");
    std::puts("Study handoff cancellation, retry and late completion passed.");
}
void studyEpochGuard() {
    QTemporaryDir directory;
    require(directory.isValid(), "Study fixture directory unavailable");
    const auto path = directory.filePath("replacement.hsproj");
    writeProject(path, Project{});
    MainWindow window;
    window.show();
    waitUntil([&] { return window.geometryReady(); });
    auto *action = window.findChild<QAction *>("openParameterSweep");
    auto *dock = dynamic_cast<SimulationDock *>(window.findChild<QDockWidget *>("simulationDock"));
    require(action && dock, "Study or simulation controls missing");
    QString failure;
    QTimer::singleShot(0, [&] {
        try {
            auto *study = dynamic_cast<ParameterSweepDialog *>(QApplication::activeModalWidget());
            require(study, "Study editor is missing");
            study->findChild<QLineEdit *>("sweepOutputDirectory")->setText(directory.path());
            auto *start = study->findChild<QPushButton *>("sweepStart");
            start->click();
            auto *timer = study->findChild<QTimer *>("studyGpuHandoffTimer");
            require(timer && timer->isActive() && !start->isEnabled(), "GPU handoff did not start");
            timer->stop();
            window.openProject(path);
            waitUntil([&] { return !window.sourceImportInProgress(); });
            require(window.sourceImportError().isEmpty(), "Replacement project did not open");
            require(QMetaObject::invokeMethod(timer, "timeout", Qt::DirectConnection),
                    "Old handoff could not be delivered");
            require(start->isEnabled(), "Old handoff continued after project change");
            dock->setRealtimeEnabled(true);
            require(dock->realtimeEnabled(), "Old study left the new project's solver reserved");
            dock->setRealtimeEnabled(false);
            bool rejected = false;
            try { study->onVariantChosen(Project{}); }
            catch (const std::invalid_argument &) { rejected = true; }
            require(rejected, "Old study could overwrite the replacement project");
            study->reject();
        } catch (const std::exception &error) {
            failure = QString::fromUtf8(error.what());
            if (auto *dialog = qobject_cast<QDialog *>(QApplication::activeModalWidget())) dialog->reject();
        }
    });
    action->trigger();
    require(failure.isEmpty(), failure.toUtf8().constData());
    std::puts("Study handoff and adoption reject an obsolete project epoch.");
}
}
int main(int argc, char **argv) {
    QApplication app(argc, argv);
    QStandardPaths::setTestModeEnabled(true);
    QCoreApplication::setOrganizationName("HornSculptGuidedTests");
    QCoreApplication::setApplicationName("GuidedQA");
    try {
        const QStringList args = QCoreApplication::arguments();
        if (args.size() == 4 && args[1] == "--write") childWrite(args[2], args[3]);
        else if (args.size() == 4 && args[1] == "--restore") childRestore(args[2], args[3]);
        else if (args.size() == 2 && args[1] == "--recovery") recoveryUi();
        else if (args.size() == 2 && args[1] == "--calculation") guidedCalculation();
        else if (args.size() == 2 && args[1] == "--focus") resultFocus();
        else if (args.size() == 2 && args[1] == "--late-cancel") lateGuidedCancel();
        else if (args.size() == 2 && args[1] == "--recovery-reject") rejectedRecovery();
        else if (args.size() == 2 && args[1] == "--comparison-epoch") comparisonEpoch();
        else if (args.size() == 2 && args[1] == "--study-idle") studyEdit(false);
        else if (args.size() == 2 && args[1] == "--study-cached") studyEdit(true);
        else if (args.size() == 2 && args[1] == "--study-handoff") studyHandoffCancel();
        else if (args.size() == 2 && args[1] == "--study-epoch") studyEpochGuard();
        else throw std::runtime_error("Unknown test mode");
        return 0;
    } catch (const std::exception &error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
