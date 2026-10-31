#include "step_export_dialog.hpp"

#include <QApplication>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QFile>
#include <QFileInfo>
#include <QLabel>
#include <QPointer>
#include <QPushButton>
#include <QSpinBox>
#include <QTemporaryDir>
#include <QTimer>
#include <chrono>
#include <functional>
#include <hornsculpt/draft_step_export.hpp>
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
template <class Widget> Widget *control(StepExportDialog &dialog, const char *name) {
    auto *widget = dialog.findChild<Widget *>(QString::fromLatin1(name));
    require(widget != nullptr, "Required STEP export control is missing.");
    return widget;
}
void until(const std::function<bool()> &done, int timeout = 60000) {
    QElapsedTimer timer;
    timer.start();
    while (!done()) {
        QCoreApplication::processEvents(QEventLoop::AllEvents, 10);
        require(timer.elapsed() < timeout, "STEP export GUI operation timed out.");
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
}
QByteArray contents(const QString &path) {
    QFile file(path);
    require(file.open(QIODevice::ReadOnly), "Cannot read STEP output.");
    return file.readAll();
}
void save(const QString &path, const QByteArray &bytes) {
    QFile file(path);
    require(file.open(QIODevice::WriteOnly), "Cannot create existing output fixture.");
    require(file.write(bytes) == bytes.size(), "Incomplete existing output fixture.");
}
std::shared_ptr<const DraftGeometry> smallHorn() {
    DraftParameters p;
    p.entry.width = p.entry.height = .02;
    p.length = .03;
    p.wallThickness = .002;
    p.horizontalAngle = p.verticalAngle = .2;
    p.profile.k = p.verticalProfile.k = 0;
    p.profile.terminationScale = p.verticalProfile.terminationScale = 0;
    return std::make_shared<const DraftGeometry>(p);
}
void publicationAndControls(const QTemporaryDir &directory) {
    StepExportDialog dialog(smallHorn(), 73);
    dialog.show();
    auto *kind = control<QComboBox>(dialog, "stepExportKind");
    auto *sections = control<QSpinBox>(dialog, "stepExportSectionCount");
    auto *tolerance = control<QDoubleSpinBox>(dialog, "stepExportTolerance");
    auto *start = control<QPushButton>(dialog, "stepExportStart");
    require(kind->count() == 3 && sections->isHidden() && tolerance->value() == .05,
            "STEP defaults or supported export variants are incorrect.");
    require(control<QLabel>(dialog, "stepExportRevision")->text().contains("73"),
            "Export snapshot revision is missing.");
    kind->setCurrentIndex(kind->findData(static_cast<int>(DraftStepExportKind::sections)));
    require(!sections->isHidden(), "Contour count is unavailable for section export.");
    sections->setValue(3);
    const auto path = directory.filePath("profiles.step");
    save(path, "existing user file");
    unsigned ticks = 0;
    QTimer heartbeat;
    QObject::connect(&heartbeat, &QTimer::timeout, [&] { ++ticks; });
    heartbeat.start(1);
    dialog.startExport(path);
    require(dialog.busy() && !kind->isEnabled() && !tolerance->isEnabled() &&
                !sections->isEnabled() && !start->isEnabled(),
            "An active export permits its immutable settings to change.");
    require(contents(path) == "existing user file",
            "Starting the worker prematurely replaced the destination.");
    until([&] { return !dialog.busy(); });
    require(dialog.lastError().isEmpty(), qPrintable(dialog.lastError()));
    const auto exported = contents(path);
    require(exported.startsWith("ISO-10303-21;") && exported.contains("END-ISO-10303-21;"),
            "Completed GUI export did not publish actual complete STEP bytes.");
    require(ticks > 0 && start->isEnabled(), "STEP construction blocked the event loop.");

    // A failed atomic destination write must report the I/O error and leave all
    // existing files intact. The backend still performs a real contour export.
    dialog.startExport(directory.path());
    until([&] { return !dialog.busy(); });
    require(!dialog.lastError().isEmpty() && contents(path) == exported,
            "Failed publication lost its error or damaged another file.");
}
void cancellationAndLifetime(const QTemporaryDir &directory) {
    const auto path = directory.filePath("cancelled.step");
    const QByteArray previous("keep this existing geometry");
    save(path, previous);
    StepExportDialog dialog(smallHorn(), 91);
    control<QComboBox>(dialog, "stepExportKind")
        ->setCurrentIndex(static_cast<int>(DraftStepExportKind::sections));
    dialog.startExport(path);
    dialog.cancelExport();
    until([&] { return !dialog.busy(); });
    require(dialog.lastError().isEmpty() && contents(path) == previous &&
                control<QLabel>(dialog, "stepExportStatus")->text().contains("Abgebrochen"),
            "Cancelled worker published output or reported a construction error.");

    auto *closing = new StepExportDialog(smallHorn(), 92);
    const QPointer<StepExportDialog> alive(closing);
    control<QComboBox>(*closing, "stepExportKind")
        ->setCurrentIndex(static_cast<int>(DraftStepExportKind::sections));
    closing->startExport(path);
    QElapsedTimer deletion;
    deletion.start();
    delete closing;
    require(alive.isNull() && deletion.elapsed() < 500,
            "Destroying the dialog waited for a background CAD operation.");
    // Let completion callbacks run after the dialog's destruction. They must
    // discard their bytes and never dereference or publish through that owner.
    QElapsedTimer drain;
    drain.start();
    until([&] { return drain.elapsed() >= 1000; }, 3000);
    require(contents(path) == previous, "A destroyed export dialog still published its result.");
}
} // namespace

int main(int argc, char **argv) {
    qputenv("QT_QPA_PLATFORM", "offscreen");
    QApplication application(argc, argv);
    application.setQuitOnLastWindowClosed(false);
    try {
        QTemporaryDir directory;
        require(directory.isValid(), "Cannot create STEP export test directory.");
        publicationAndControls(directory);
        cancellationAndLifetime(directory);
        std::cout << "STEP export GUI: real bytes, atomic publication, responsive cancellation "
                     "and safe dialog destruction passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
