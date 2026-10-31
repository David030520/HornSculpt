#include "mesh_import_dialog.hpp"
#include <QApplication>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QElapsedTimer>
#include <QFile>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QTemporaryDir>
#include <QThread>
#include <functional>
#include <iostream>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
void waitFor(const std::function<bool()> &condition) {
    QElapsedTimer timer;
    timer.start();
    while (!condition() && timer.elapsed() < 5000) {
        QCoreApplication::processEvents();
        QThread::msleep(1);
    }
    require(condition(), "Asynchronous import timed out");
}
} // namespace
int main(int argc, char **argv) {
    QApplication app(argc, argv);
    try {
        QTemporaryDir dir;
        require(dir.isValid(), "Temporary directory failed");
        const auto path = dir.filePath("mesh.msh");
        QFile file(path);
        require(file.open(QIODevice::WriteOnly), "Cannot create fixture");
        file.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n1\n2 7 "
                   "\"Mund\"\n$EndPhysicalNames\n"
                   "$Nodes\n3\n1 0 0 0\n2 1000 0 0\n3 0 1000 0\n$EndNodes\n"
                   "$Elements\n1\n4 2 2 7 5 1 2 3\n$EndElements\n");
        file.close();
        MeshImportDialog dialog(path);
        auto *units = dialog.findChild<QComboBox *>("meshImportUnits");
        auto *groups = dialog.findChild<QListWidget *>("meshImportGroups");
        auto *status = dialog.findChild<QLabel *>("meshImportStatus");
        auto *buttons = dialog.findChild<QDialogButtonBox *>();
        require(units && groups && status && buttons, "Missing import controls");
        require(!buttons->button(QDialogButtonBox::Ok)->isEnabled(), "Unit must be explicit");
        units->setCurrentIndex(1);
        units->setCurrentIndex(2);
        waitFor([&] { return buttons->button(QDialogButtonBox::Ok)->isEnabled(); });
        require(buttons->button(QDialogButtonBox::Ok)->isEnabled() && groups->count() == 1 &&
                    groups->item(0)->text().contains("Mund") && status->text().contains("3"),
                "Missing group or open boundary diagnostics");
        dialog.accept();
        require(dialog.result() == QDialog::Accepted && dialog.selection().has_value(),
                "Import not committed");
        require(dialog.selection()->lengthScaleToMetres == .001 &&
                    dialog.selection()->imported.mesh.vertices[1].x == 1.,
                "Millimetre conversion failed");
        MeshImportDialog bad(dir.filePath("missing.msh"));
        bad.findChild<QComboBox *>("meshImportUnits")->setCurrentIndex(1);
        waitFor([&] {
            return bad.findChild<QLabel *>("meshImportStatus")->text().contains("fehlgeschlagen");
        });
        require(!bad.findChild<QDialogButtonBox *>()->button(QDialogButtonBox::Ok)->isEnabled(),
                "Missing file accepted");
        std::cout << "Mesh import dialog tests PASS.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
