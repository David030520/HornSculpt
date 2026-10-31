#include "chart_scaling.hpp"
#include <QCheckBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QLineEdit>
#include <QLocale>
#include <QMessageBox>
#include <QSettings>
#include <QVBoxLayout>
#include <QWidget>
#include <cmath>
#include <stdexcept>

namespace hornsculpt {
namespace {
QLineEdit *numberField(QFormLayout &layout, const QString &label, const QString &name,
                       std::optional<double> number) {
    auto *field = new QLineEdit;
    field->setObjectName(name);
    if (number)
        field->setText(QLocale().toString(*number, 'g', 17));
    layout.addRow(label, field);
    return field;
}
double number(const QLineEdit &field) {
    bool valid = false;
    const double result = QLocale().toDouble(field.text(), &valid);
    if (!valid || !std::isfinite(result))
        throw std::invalid_argument("Bitte eine endliche Zahl eingeben.");
    return result;
}
struct AxisEditor {
    QCheckBox *automatic{};
    QLineEdit *minimum{}, *maximum{}, *step{};
    AxisEditor(QVBoxLayout &parent, const QString &label, const QString &prefix,
               PlotAxisScale scale, bool allowAuto = true) {
        auto *group = new QGroupBox(label);
        auto *form = new QFormLayout(group);
        automatic = new QCheckBox(QObject::tr("Automatisch"), group);
        automatic->setObjectName(prefix + "Automatic");
        automatic->setChecked(!scale.minimum);
        if (allowAuto)
            form->addRow(automatic);
        else
            automatic->hide();
        minimum = numberField(*form, QObject::tr("Minimum"), prefix + "Minimum", scale.minimum);
        maximum = numberField(*form, QObject::tr("Maximum"), prefix + "Maximum", scale.maximum);
        step = numberField(*form, QObject::tr("Schritt"), prefix + "Step", scale.tickStep);
        step->setToolTip(QObject::tr("Teilstrichabstand in Achseneinheiten. Leer: automatisch. "
                                     "Sehr dichte Teilstriche werden beim Zeichnen ausgedünnt."));
        auto enabled = [this] {
            minimum->setEnabled(!automatic->isChecked());
            maximum->setEnabled(!automatic->isChecked());
        };
        QObject::connect(automatic, &QCheckBox::toggled, group, [enabled](bool) { enabled(); });
        enabled();
        parent.addWidget(group);
    }
    PlotAxisScale value() const {
        PlotAxisScale s;
        if (!automatic->isChecked()) {
            s.minimum = number(*minimum);
            s.maximum = number(*maximum);
        }
        if (!step->text().trimmed().isEmpty())
            s.tickStep = number(*step);
        return s;
    }
};
template <class Function> void buttons(QDialog &dialog, QVBoxLayout &layout, Function validate) {
    auto *box = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    layout.addWidget(box);
    QObject::connect(box, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    QObject::connect(box, &QDialogButtonBox::accepted, &dialog, [&dialog, validate] {
        try {
            validate();
            dialog.accept();
        } catch (const std::exception &e) {
            QMessageBox::warning(&dialog, QObject::tr("Skalierung"), QString::fromUtf8(e.what()));
        }
    });
}
QString key(const QWidget &widget, const QString &axis) {
    return "Charts/" + widget.objectName() + "/scales/" + axis;
}
} // namespace
std::optional<PlotAxisScales> editAxisScales(QWidget &owner, PlotAxisScales scales,
                                             const QString &x, const QString &y, bool logarithmic,
                                             bool decibelY) {
    QDialog dialog(&owner);
    dialog.setWindowTitle(QObject::tr("Skalierung"));
    QVBoxLayout layout(&dialog);
    AxisEditor xEditor(layout, x, "scaleX", scales.x), yEditor(layout, y, "scaleY", scales.y);
    buttons(dialog, layout, [&] {
        scales = {xEditor.value(), yEditor.value()};
        validatePlotAxisScale(scales.x, logarithmic);
        validatePlotAxisScale(scales.y);
        if (decibelY && scales.y.minimum && scales.y.tickStep)
            validatePlotDbScale({*scales.y.minimum, *scales.y.maximum, *scales.y.tickStep});
    });
    return dialog.exec() == QDialog::Accepted ? std::optional{scales} : std::nullopt;
}
std::optional<PolarScales> editPolarScales(QWidget &owner, PolarScales scales) {
    QDialog dialog(&owner);
    dialog.setWindowTitle(QObject::tr("Skalierung"));
    QVBoxLayout layout(&dialog);
    AxisEditor db(layout, QObject::tr("Pegel [dB]"), "scaleDb",
                  {scales.decibels.minimum, scales.decibels.maximum, scales.decibels.tickStep},
                  false);
    auto *form = new QFormLayout;
    auto *angle = numberField(*form, QObject::tr("Winkelschritt [°]"), "scaleAngularStep",
                              scales.angularStep);
    layout.addLayout(form);
    buttons(dialog, layout, [&] {
        const auto s = db.value();
        if (!s.tickStep)
            throw std::invalid_argument("Der dB-Schritt darf nicht leer sein.");
        scales.decibels = {*s.minimum, *s.maximum, *s.tickStep};
        validatePlotDbScale(scales.decibels);
        scales.angularStep = number(*angle);
        if (!(scales.angularStep > 0))
            throw std::invalid_argument("Der Winkelschritt muss positiv sein.");
    });
    return dialog.exec() == QDialog::Accepted ? std::optional{scales} : std::nullopt;
}
std::optional<HeatmapScales> editHeatmapScales(QWidget &owner, HeatmapScales scales) {
    QDialog dialog(&owner);
    dialog.setWindowTitle(QObject::tr("Skalierung"));
    QVBoxLayout layout(&dialog);
    AxisEditor x(layout, QObject::tr("Frequenz [Hz]"), "scaleX", scales.axes.x);
    AxisEditor y(layout, QObject::tr("Winkel [°]"), "scaleY", scales.axes.y);
    AxisEditor db(layout, QObject::tr("Pegel [dB]"), "scaleDb",
                  {scales.decibels.minimum, scales.decibels.maximum, scales.decibels.tickStep},
                  false);
    buttons(dialog, layout, [&] {
        scales.axes = {x.value(), y.value()};
        validatePlotAxisScale(scales.axes.x, true);
        validatePlotAxisScale(scales.axes.y);
        const auto s = db.value();
        if (!s.tickStep)
            throw std::invalid_argument("Der dB-Schritt darf nicht leer sein.");
        scales.decibels = {*s.minimum, *s.maximum, *s.tickStep};
        validatePlotDbScale(scales.decibels);
    });
    return dialog.exec() == QDialog::Accepted ? std::optional{scales} : std::nullopt;
}
PlotAxisScale readAxisScale(const QWidget &widget, const QString &axis, PlotAxisScale fallback) {
    if (widget.objectName().isEmpty())
        return fallback;
    QSettings settings;
    settings.beginGroup(key(widget, axis));
    if (!settings.contains("automatic"))
        return fallback;
    PlotAxisScale result;
    if (!settings.value("automatic").toBool()) {
        result.minimum = settings.value("minimum").toDouble();
        result.maximum = settings.value("maximum").toDouble();
    }
    if (settings.contains("step"))
        result.tickStep = settings.value("step").toDouble();
    try {
        validatePlotAxisScale(result);
    } catch (const std::invalid_argument &) {
        return fallback;
    }
    return result;
}
void saveAxisScale(const QWidget &widget, const QString &axis, PlotAxisScale scale) {
    if (widget.objectName().isEmpty())
        return;
    QSettings settings;
    settings.beginGroup(key(widget, axis));
    settings.remove("");
    settings.setValue("automatic", !scale.minimum);
    if (scale.minimum) {
        settings.setValue("minimum", *scale.minimum);
        settings.setValue("maximum", *scale.maximum);
    }
    if (scale.tickStep)
        settings.setValue("step", *scale.tickStep);
}
double readAngularStep(const QWidget &widget, double fallback) {
    if (widget.objectName().isEmpty())
        return fallback;
    const double value = QSettings().value(key(widget, "angularStep"), fallback).toDouble();
    return std::isfinite(value) && value > 0 ? value : fallback;
}
void saveAngularStep(const QWidget &widget, double step) {
    if (!widget.objectName().isEmpty())
        QSettings().setValue(key(widget, "angularStep"), step);
}
} // namespace hornsculpt
