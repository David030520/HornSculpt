#include "acoustic_frame_dialog.hpp"

#include <QApplication>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QLineEdit>
#include <QLocale>
#include <QPushButton>
#include <QTest>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace hornsculpt;

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

template <class Control> Control *control(AcousticFrameDialog &dialog, const char *name) {
    auto *result = dialog.findChild<Control *>(QString::fromLatin1(name));
    require(result != nullptr, "Missing acoustic-frame control.");
    return result;
}

bool near(Vec3 a, Vec3 b) {
    return std::hypot(a.x - b.x, a.y - b.y, a.z - b.z) < 1e-10;
}

void setDirection(AcousticFrameDialog &dialog, const char *prefix, Vec3 vector) {
    const auto name = QString::fromLatin1(prefix);
    auto *preset = dialog.findChild<QComboBox *>(name + QStringLiteral("Preset"));
    require(preset != nullptr, "Missing direction preset.");
    preset->setCurrentIndex(6);
    const std::array<double, 3> values{vector.x, vector.y, vector.z};
    const std::array<const char *, 3> suffixes{"X", "Y", "Z"};
    for (std::size_t i = 0; i < values.size(); ++i) {
        auto *input = dialog.findChild<QDoubleSpinBox *>(name + QString::fromLatin1(suffixes[i]));
        require(input && !input->isReadOnly(), "Free direction cannot be edited.");
        input->setValue(values[i]);
    }
}

void geometryMode() {
    const AcousticFrame source{{.0123, -.0456, .0789}, {0, 1, 0}, {1, 0, 0}};
    AcousticFrameDialog dialog(source);
    require(control<QComboBox>(dialog, "acousticFrameMode")->currentIndex() == 0,
            "A missing override does not follow geometry.");
    require(std::abs(control<QDoubleSpinBox>(dialog, "acousticOriginX")->value() - 12.3) < 1e-12,
            "The frame origin is not displayed in millimetres.");
    require(control<QDoubleSpinBox>(dialog, "acousticOriginX")->isReadOnly() &&
                control<QDoubleSpinBox>(dialog, "acousticForwardY")->isReadOnly() &&
                !control<QComboBox>(dialog, "acousticHorizontalPreset")->isEnabled(),
            "Geometry-derived coordinates can be edited.");
    require(control<QComboBox>(dialog, "acousticForwardPreset")->currentIndex() == 2,
            "The geometry's forward direction was silently replaced by draft +Z.");
    require(control<QLabel>(dialog, "acousticVerticalZ")->text() ==
                dialog.locale().toString(-1., 'f', 6),
            "The vertical direction does not follow forward cross horizontal.");
    dialog.accept();
    require(dialog.result() == QDialog::Accepted && !dialog.observationFrame(),
            "Geometry mode stores a fixed override.");
}

void customAxesAndCancel() {
    const AcousticFrame source{{.01, .02, .03}, {0, 0, 1}, {1, 0, 0}};
    const AcousticFrame initial{{.1, .2, .3}, {-1, 0, 0}, {0, 0, 1}};
    AcousticFrameDialog dialog(source, initial);
    auto *mode = control<QComboBox>(dialog, "acousticFrameMode");
    auto *origin = control<QDoubleSpinBox>(dialog, "acousticOriginX");
    require(mode->currentIndex() == 1 && !origin->isReadOnly(),
            "An existing override is not editable.");
    origin->setValue(42.5);
    control<QComboBox>(dialog, "acousticForwardPreset")->setCurrentIndex(2);
    control<QComboBox>(dialog, "acousticHorizontalPreset")->setCurrentIndex(1);
    require(dialog.observationFrame() == initial, "Editing committed an unaccepted override.");
    mode->setCurrentIndex(0);
    require(origin->value() == 10 && origin->isReadOnly(),
            "Switching to geometry does not show the actual source frame.");
    mode->setCurrentIndex(1);
    require(origin->value() == 42.5 &&
                control<QComboBox>(dialog, "acousticHorizontalPreset")->currentIndex() == 1,
            "Switching modes loses custom origin or roll.");
    dialog.accept();
    const auto result = dialog.observationFrame();
    require(result && near(result->origin, {.0425, .2, .3}) && near(result->forward, {0, 1, 0}) &&
                near(result->horizontal, {-1, 0, 0}) && near(verticalDirection(*result), {0, 0, 1}),
            "Accepted axis presets, roll or millimetre conversion are wrong.");

    AcousticFrameDialog cancelled(source, initial);
    control<QDoubleSpinBox>(cancelled, "acousticOriginY")->setValue(-750);
    control<QPushButton>(cancelled, "acousticFrameCancel")->click();
    require(cancelled.result() == QDialog::Rejected && cancelled.observationFrame() == initial,
            "Cancel changes the original override.");

    AcousticFrameDialog reset(source, initial);
    control<QComboBox>(reset, "acousticFrameMode")->setCurrentIndex(0);
    reset.accept();
    require(!reset.observationFrame(), "Geometry mode cannot clear an existing override.");
}

void validationAndFreeRoll() {
    AcousticFrameDialog dialog({}, AcousticFrame{});
    auto *ok = control<QPushButton>(dialog, "acousticFrameAccept");
    auto *error = control<QLabel>(dialog, "acousticFrameError");
    auto *horizontal = control<QComboBox>(dialog, "acousticHorizontalPreset");
    horizontal->setCurrentIndex(4);
    require(!ok->isEnabled() && error->text().contains(QStringLiteral("parallel")) &&
                horizontal->currentIndex() == 4,
            "Parallel presets are accepted or silently replaced.");
    dialog.accept();
    require(dialog.result() != QDialog::Accepted, "Programmatic acceptance bypasses validation.");
    horizontal->setCurrentIndex(5);
    require(!ok->isEnabled(), "Antiparallel directions are accepted.");
    setDirection(dialog, "acousticForward", {0, 0, 0});
    require(!ok->isEnabled() && error->text().contains(QStringLiteral("null")),
            "A zero forward vector is accepted.");
    setDirection(dialog, "acousticForward", {0, 0, 2});
    setDirection(dialog, "acousticHorizontal", {0, 0, 0});
    require(!ok->isEnabled() && error->text().startsWith(QStringLiteral("H-Richtung")),
            "A zero horizontal vector lacks a specific error.");
    setDirection(dialog, "acousticHorizontal", {2, 0, 1});
    require(ok->isEnabled() && error->text().isEmpty(),
            "A valid oblique horizontal direction is rejected.");
    require(horizontal->currentIndex() == 6,
            "Editing a free vector unexpectedly switches back to a fixed preset.");
    dialog.accept();
    const auto frame = dialog.observationFrame();
    require(frame && near(frame->forward, {0, 0, 1}) && near(frame->horizontal, {1, 0, 0}) &&
                near(verticalDirection(*frame), {0, 1, 0}),
            "Free directions do not use the shared normalized frame contract.");
}

void pendingNumericInput() {
    AcousticFrameDialog dialog({}, AcousticFrame{});
    dialog.setLocale(QLocale(QLocale::German, QLocale::Germany));
    dialog.show();
    QCoreApplication::processEvents();
    auto *origin = control<QDoubleSpinBox>(dialog, "acousticOriginZ");
    auto *editor = origin->findChild<QLineEdit *>();
    require(editor != nullptr, "The numeric editor is inaccessible.");
    editor->setFocus();
    editor->selectAll();
    QTest::keyClicks(editor, "123,25");
    require(origin->value() == 0, "An unfinished coordinate edit committed early.");
    dialog.accept();
    require(dialog.result() == QDialog::Accepted && dialog.observationFrame() &&
                near(dialog.observationFrame()->origin, {0, 0, .12325}),
            "Accept discards a still-focused coordinate entry.");
}
} // namespace

int main(int argc, char **argv) {
    QApplication application(argc, argv);
    try {
        geometryMode();
        customAxesAndCancel();
        validationAndFreeRoll();
        pendingNumericInput();
        if (argc == 3 && QString::fromLocal8Bit(argv[1]) == QStringLiteral("--screenshot")) {
            const auto frame =
                defaultAcousticFrame({.311700818, 0, .186599442}, {-.002264232, -.999997437, 0});
            AcousticFrameDialog dialog(frame, frame);
            dialog.setLocale(QLocale(QLocale::German, QLocale::Germany));
            dialog.show();
            QTest::qWait(100);
            require(dialog.grab().save(QString::fromLocal8Bit(argv[2])),
                    "The acoustic frame dialog screenshot could not be saved.");
        }
        std::cout << "Acoustic frame dialog: geometry, presets, free roll, validation, "
                     "units and transactional editing passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
