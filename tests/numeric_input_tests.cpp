#include "numeric_input.hpp"

#include <QApplication>
#include <QHBoxLayout>
#include <QLineEdit>
#include <QLocale>
#include <QPushButton>
#include <QTest>
#include <QWheelEvent>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace hornsculpt;

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

void type(QAbstractSpinBox &input, const char *text) {
    auto *editor = input.findChild<QLineEdit *>();
    QTest::mouseClick(editor, Qt::LeftButton, Qt::NoModifier,
                      QPoint(editor->width() - 25, editor->height() / 2));
    QCoreApplication::processEvents();
    QTest::keyClicks(editor, text);
}

void decimals(const QLocale &locale) {
    QWidget window;
    QHBoxLayout layout(&window);
    QPushButton other("Focus");
    DecimalInput input;
    input.setLocale(locale);
    input.setRange(.001, 100);
    input.setDecimals(3);
    input.setSuffix(" mm");
    layout.addWidget(&other);
    layout.addWidget(&input);
    window.show();
    window.activateWindow();
    QCoreApplication::processEvents();
    for (const auto text : {"12", "0.25", "0,25"}) {
        input.setValue(8);
        other.setFocus();
        QCoreApplication::processEvents();
        type(input, text);
        require(input.value() == 8, "Typing commits an unfinished numeric entry");
        QTest::keyClick(input.findChild<QLineEdit *>(), Qt::Key_Return);
        const double expected = text[0] == '1' ? 12 : .25;
        if (std::abs(input.value() - expected) >= 1e-12)
            std::cerr << locale.name().toStdString() << " typed " << text << ": "
                      << input.text().toStdString() << "; focused " << input.hasFocus() << '\n';
        require(std::abs(input.value() - expected) < 1e-12,
                "Click and type failed, or decimal point was interpreted as grouping");
    }
    other.setFocus();
    QCoreApplication::processEvents();
    type(input, "20,5");
    QTest::keyClick(input.findChild<QLineEdit *>(), Qt::Key_Tab);
    require(input.value() == 20.5, "Tab does not commit numeric input");
    auto *editor = input.findChild<QLineEdit *>();
    editor->selectAll();
    editor->insert("0.125");
    QTest::keyClick(editor, Qt::Key_Return);
    require(input.value() == .125, "Inserted decimal text is misinterpreted");
    for (const bool focused : {false, true}) {
        const double previous = input.value();
        (focused ? static_cast<QWidget *>(&input) : static_cast<QWidget *>(&other))->setFocus();
        QWheelEvent wheel(QPointF(10, 10), input.mapToGlobal(QPoint(10, 10)), {}, QPoint(0, 120),
                          Qt::NoButton, Qt::NoModifier, Qt::NoScrollPhase, false);
        QApplication::sendEvent(&input, &wheel);
        require(input.value() == previous + input.singleStep() && wheel.isAccepted(),
                "Wheel does not step a decimal field");
    }
    const double previous = input.value();
    QTest::keyClick(&input, Qt::Key_Up);
    require(input.value() == previous + input.singleStep(), "Arrow keys no longer step");
}

void integer() {
    QWidget window;
    QHBoxLayout layout(&window);
    QPushButton other("Focus");
    IntegerInput input;
    input.setRange(1, 100000);
    input.setValue(17);
    layout.addWidget(&other);
    layout.addWidget(&input);
    window.show();
    window.activateWindow();
    other.setFocus();
    QCoreApplication::processEvents();
    type(input, "240");
    QTest::keyClick(input.findChild<QLineEdit *>(), Qt::Key_Return);
    require(input.value() == 240, "Integer click and type failed");
    QWheelEvent wheel(QPointF(10, 10), input.mapToGlobal(QPoint(10, 10)), {}, QPoint(0, -120),
                      Qt::NoButton, Qt::NoModifier, Qt::NoScrollPhase, false);
    QApplication::sendEvent(&input, &wheel);
    require(input.value() == 239 && wheel.isAccepted(), "Wheel does not step an integer field");
}
} // namespace

int main(int argc, char **argv) {
    QApplication application(argc, argv);
    try {
        decimals(QLocale(QLocale::German, QLocale::Germany));
        decimals(QLocale(QLocale::English, QLocale::UnitedStates));
        integer();
        std::cout << "Numeric mouse/keyboard editing, decimal locales and commit tests passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
