#include "parameter_panel.hpp"
#include "numeric_input.hpp"
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QHeaderView>
#include <QLabel>
#include <QPushButton>
#include <QScrollArea>
#include <QTableWidget>
#include <QTabWidget>
#include <QVBoxLayout>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace hornsculpt {
namespace {
constexpr double degrees = 180.0 / std::numbers::pi;
QGroupBox *section(const QString &title, QVBoxLayout *layout, QFormLayout *&form) {
    auto *group = new QGroupBox(title);
    form = new QFormLayout(group);
    form->setContentsMargins(9, 13, 9, 9);
    form->setSpacing(5);
    form->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
    layout->addWidget(group);
    return group;
}
} // namespace

QDoubleSpinBox *ParameterPanel::createNumber(const QString &key, double value, double min,
                                             double max, int decimals, const QString &suffix,
                                             const QString &explanation) {
    auto *spin = new DecimalInput;
    spin->setObjectName(key);
    spin->setRange(min, max);
    spin->setDecimals(decimals);
    spin->setValue(value);
    spin->setSuffix(suffix);
    spin->setSingleStep(decimals == 0 ? 1.0 : .1);
    spin->setToolTip(explanation);
    spin->setWhatsThis(explanation);
    numbers_[key] = spin;
    connect(spin, &QDoubleSpinBox::valueChanged, this, [this, key] {
        if (!updating_)
            editedNumbers_.insert(key);
        changed();
    });
    return spin;
}
QDoubleSpinBox *ParameterPanel::addNumber(QFormLayout *form, const QString &key,
                                          const QString &label, double value, double min,
                                          double max, int decimals, const QString &suffix,
                                          const QString &explanation) {
    auto *spin = createNumber(key, value, min, max, decimals, suffix, explanation);
    auto *caption = new QLabel(label);
    caption->setBuddy(spin);
    caption->setToolTip(explanation);
    form->addRow(caption, spin);
    return spin;
}

ParameterPanel::ParameterPanel(QWidget *parent) : QWidget(parent) {
    auto *outer = new QVBoxLayout(this);
    outer->setContentsMargins(0, 0, 0, 0);
    geometryStatus_ = new QLabel;
    geometryStatus_->setObjectName("draftGeometryStatus");
    geometryStatus_->setTextFormat(Qt::PlainText);
    geometryStatus_->setWordWrap(true);
    geometryStatus_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    geometryStatus_->hide();
    outer->addWidget(geometryStatus_);
    sections_ = new QTabWidget;
    sections_->setObjectName("draftSections");
    sections_->setAccessibleName(tr("Draft-Parameterbereiche"));
    outer->addWidget(sections_);
    const auto page = [this](const QString &title, const QString &name) {
        auto *scroll = new QScrollArea;
        scroll->setWidgetResizable(true);
        scroll->setFrameShape(QFrame::NoFrame);
        scroll->setObjectName(name);
        auto *content = new QWidget;
        auto *pageLayout = new QVBoxLayout(content);
        pageLayout->setContentsMargins(7, 7, 7, 7);
        pageLayout->setSpacing(7);
        scroll->setWidget(content);
        sections_->addTab(scroll, title);
        return pageLayout;
    };
    auto *layout = page(tr("&Profil"), "draftProfileScroll");
    QFormLayout *form{};
    section(tr("Basis · OS-SE"), layout, form);
    addNumber(form, "entryWidth", tr("Eintrittsbreite"), 25.4, .01, 2000, 3, " mm",
              tr("Volle Breite des Eintrittsquerschnitts entlang X bei z = 0. Einheit: mm."));
    addNumber(form, "entryHeight", tr("Eintrittshöhe"), 25.4, .01, 2000, 3, " mm",
              tr("Volle Höhe des Eintrittsquerschnitts entlang Y bei z = 0. Einheit: mm."));
    addNumber(form, "entryExponent", tr("Eintrittsexponent"), 2, 2, 100, 3, "",
              tr("Superellipsen-Exponent des Eintritts: 2 ist elliptisch, größere Werte sind "
                 "rechteckiger. Unabhängig von Zielquerschnitt und Termination."));
    addNumber(form, "length", tr("Länge L"), 146, 1, 10000, 2, " mm",
              tr("Axiale Länge vom Hals zum Mund in Richtung +Z. Einheit: mm."));
    addNumber(form, "thickness", tr("Wandstärke"), 6, .1, 500, 2, " mm",
              tr("Radiale Materialstärke der geschlossenen Außenwand. Einheit: mm."));
    auto *profile = new QGroupBox(tr("Profil · OS-SE"));
    auto *profileGrid = new QGridLayout(profile);
    profileGrid->setContentsMargins(9, 13, 9, 9);
    profileGrid->setHorizontalSpacing(5);
    profileGrid->setVerticalSpacing(5);
    profileGrid->setColumnStretch(1, 1);
    profileGrid->setColumnStretch(2, 1);
    asymmetricEnabled_ = new QCheckBox(tr("Asymmetrisch · vier Seiten"));
    asymmetricEnabled_->setObjectName("asymmetricEnabled");
    asymmetricEnabled_->setToolTip(
        tr("Eigene OS-SE-Profile für rechts (+X), links (−X), oben (+Y) und unten (−Y). "
           "Bei ausgeschaltetem Modus werden H und V gespiegelt; die gespeicherten "
           "Gegenseiten bleiben erhalten. Eintritt, Länge und Wandstärke gelten gemeinsam."));
    asymmetricEnabled_->setWhatsThis(asymmetricEnabled_->toolTip());
    profileGrid->addWidget(asymmetricEnabled_, 0, 0, 1, 3);
    horizontalProfileLabel_ = new QLabel;
    verticalProfileLabel_ = new QLabel;
    horizontalProfileLabel_->setAlignment(Qt::AlignCenter);
    verticalProfileLabel_->setAlignment(Qt::AlignCenter);
    horizontalProfileLabel_->setToolTip(tr("Rechts: +X, φ = 0°. Symmetrisch: beide H-Seiten."));
    verticalProfileLabel_->setToolTip(tr("Oben: +Y, φ = 90°. Symmetrisch: beide V-Seiten."));
    profileGrid->addWidget(horizontalProfileLabel_, 1, 1);
    profileGrid->addWidget(verticalProfileLabel_, 1, 2);
    auto *opposite = new QGroupBox(tr("Gegenseiten · OS-SE"));
    oppositeProfiles_ = opposite;
    auto *oppositeGrid = new QGridLayout(opposite);
    oppositeGrid->setContentsMargins(9, 13, 9, 9);
    oppositeGrid->setHorizontalSpacing(5);
    oppositeGrid->setVerticalSpacing(5);
    for (int column = 1; column <= 2; ++column) {
        oppositeGrid->setColumnStretch(column, 1);
        auto *label = new QLabel(column == 1 ? tr("Links") : tr("Unten"));
        label->setAlignment(Qt::AlignCenter);
        label->setToolTip(column == 1 ? tr("Links: −X, φ = 180°.") : tr("Unten: −Y, φ = 270°."));
        oppositeGrid->addWidget(label, 0, column);
    }
    int profileRow = 1;
    const auto addProfile = [&](const QString &label, const QString &keyH, const QString &keyV,
                                const QString &keyL, const QString &keyB, double valueH,
                                double valueV, double minimum, double maximum, int decimals,
                                const QString &suffix, const QString &explanation) {
        const auto pair = [&](QGridLayout *grid, int row, bool negative) {
            auto *caption = new QLabel(label);
            caption->setToolTip(explanation);
            grid->addWidget(caption, row, 0);
            for (int column = 1; column <= 2; ++column) {
                const auto direction =
                    negative ? (column == 1 ? tr("Links (−X)") : tr("Unten (−Y)"))
                             : (column == 1 ? tr("H / Rechts (+X)") : tr("V / Oben (+Y)"));
                const auto &key =
                    negative ? (column == 1 ? keyL : keyB) : (column == 1 ? keyH : keyV);
                auto *input =
                    createNumber(key, column == 1 ? valueH : valueV, minimum, maximum, decimals,
                                 suffix, direction + QStringLiteral(": ") + explanation);
                input->setAccessibleName(label + QStringLiteral(" ") + direction);
                grid->addWidget(input, row, column);
                if (column == 1)
                    caption->setBuddy(input);
            }
        };
        pair(profileGrid, profileRow + 1, false);
        pair(oppositeGrid, profileRow, true);
        ++profileRow;
    };
    addProfile(
        tr("Halbwinkel α"), "angleH", "angleV", "angleL", "angleB", 40, 25, 0, 85, 2, " °",
        tr("Halber Öffnungswinkel des OS-Grundprofils in Grad; kein akustischer "
           "Abstrahlwinkel. Bei aktiver Leitkurve bleibt dieser Winkel im Formübergang "
           "wirksam. Die endgültige Kontur trifft trotzdem die Leitkurve."));
    addProfile(tr("Halswinkel α₀"), "throatAngle", "throatAngleV", "throatAngleL", "throatAngleB",
               0, 0, 0, 85, 2, " °",
               tr("Halber Anfangswinkel des OS-SE-Profilabschnitts in Grad; bei aktivem Spalt "
                  "beginnt dieser am Spaltausgang. Ein glatter Spaltanschluss verändert den "
                  "tatsächlichen Wandwinkel in diesem Bereich."));
    addProfile(tr("Formfaktor k"), "k", "kV", "kL", "kB", 1, 1, 0, 100, 3, "",
               tr("OS-Formfaktor: 1 entspricht OS, 0 einem konischen Grundprofil. Dimensionslos."));
    addProfile(tr("Termination s"), "s", "sV", "sL", "sB", .2, .2, 0, 10, 3, "",
               tr("Terminationsstärke; 0 schaltet die Termination ab. Dimensionslos."));
    addProfile(
        tr("Termination q"), "q", "qV", "qL", "qB", .995, .995, .01, 1, 4, "",
        tr("Relative Endposition der Termination. Bei 1 kann die Mundsteigung unendlich werden."));
    addProfile(tr("Termination n"), "terminationN", "terminationNV", "terminationNL",
               "terminationNB", 4, 4, 2, 100, 2, "",
               tr("Terminations-Exponent ≥ 2, unabhängig vom Querschnittsexponenten."));
    auto *copyProfiles = new QPushButton(tr("Rechts → links · Oben → unten"));
    copyProfiles->setObjectName("copyOppositeProfiles");
    copyProfiles->setToolTip(tr("Übernimmt alle sechs Profilwerte der rechten und oberen Seite "
                                "für links und unten. Die bisherigen Gegenseiten werden ersetzt."));
    copyProfiles->setWhatsThis(copyProfiles->toolTip());
    oppositeGrid->addWidget(copyProfiles, profileRow, 0, 1, 3);
    connect(copyProfiles, &QPushButton::clicked, this, [this] {
        updating_ = true;
        const auto copy = [&](const QString &sourceKey, const QString &targetKey, double source,
                              double &target, double units = 1) {
            target = modelNumber(sourceKey, source, units);
            numbers_.at(targetKey)->setValue(target * units);
            editedNumbers_.erase(targetKey);
        };
        copy("angleH", "angleL", model_.horizontalAngle, model_.leftAngle, degrees);
        copy("angleV", "angleB", model_.verticalAngle, model_.bottomAngle, degrees);
        const auto copyProfile = [&](const QString &sourceSuffix, const QString &targetSuffix,
                                     const OsseDirectionalParameters &source,
                                     OsseDirectionalParameters &target) {
            copy("throatAngle" + sourceSuffix, "throatAngle" + targetSuffix, source.throatAngle,
                 target.throatAngle, degrees);
            copy("k" + sourceSuffix, "k" + targetSuffix, source.k, target.k);
            copy("s" + sourceSuffix, "s" + targetSuffix, source.terminationScale,
                 target.terminationScale);
            copy("q" + sourceSuffix, "q" + targetSuffix, source.terminationQ, target.terminationQ);
            copy("terminationN" + sourceSuffix, "terminationN" + targetSuffix,
                 source.terminationExponent, target.terminationExponent);
        };
        copyProfile({}, "L", model_.profile, model_.leftProfile);
        copyProfile("V", "B", model_.verticalProfile, model_.bottomProfile);
        updating_ = false;
        changed();
    });
    layout->addWidget(profile);
    layout->addWidget(opposite);
    guideAngleHint_ = new QPushButton(tr("Leitkurven-Ziel einstellen …"));
    guideAngleHint_->setObjectName("guideAngleHint");
    guideAngleHint_->setToolTip(tr("Die Leitkurve legt einen Zielquerschnitt fest. "
                                   "Die Grundprofilwinkel bleiben im Formübergang wirksam. "
                                   "Hier zu Zielkontur und Position wechseln."));
    guideAngleHint_->setWhatsThis(guideAngleHint_->toolTip());
    connect(guideAngleHint_, &QPushButton::clicked, this, [this] {
        sections_->setCurrentIndex(2);
        guideEnabled_->setFocus();
    });
    layout->insertWidget(1, guideAngleHint_);
    layout->addStretch();
    layout = page(tr("&Form"), "draftShapeScroll");
    section(tr("Querschnitt"), layout, form);
    addNumber(form, "crossExponent", tr("Superellipse n"), 4, 2, 100, 2, "",
              tr("Exponent des Zielquerschnitts. n = 2 ist elliptisch; höhere Werte ergeben "
                 "rechteckigere Konturen. Der Eintritt besitzt einen eigenen Exponenten."));
    guideShapeHint_ = new QPushButton(tr("Zielkontur unter Leitkurve einstellen …"));
    guideShapeHint_->setObjectName("guideShapeHint");
    guideShapeHint_->setToolTip(tr("Bei aktiver Leitkurve ersetzt deren Kontur den direkten "
                                   "Querschnittsexponenten. Der Formverlauf bleibt wirksam."));
    guideShapeHint_->setWhatsThis(guideShapeHint_->toolTip());
    connect(guideShapeHint_, &QPushButton::clicked, this, [this] {
        sections_->setCurrentIndex(2);
        curveType_->setFocus();
    });
    form->addRow(guideShapeHint_);
    transitionMethod_ = new QComboBox;
    transitionMethod_->setObjectName("transitionMethod");
    transitionMethod_->addItems({tr("Bereich · glatt"), tr("Stützpunkte (LUT)")});
    transitionMethod_->setToolTip(
        tr("Steuert den Formübergang über die axiale Länge. Die LUT ordnet jeder relativen "
           "Position einen Formanteil zwischen 0 und 1 zu."));
    transitionMethod_->setWhatsThis(transitionMethod_->toolTip());
    form->addRow(tr("Formübergang"), transitionMethod_);
    rangeSettings_ = new QWidget;
    auto *rangeForm = new QFormLayout(rangeSettings_);
    rangeForm->setContentsMargins(0, 0, 0, 0);
    addNumber(rangeForm, "transitionStart", tr("Beginn"), 0, 0, 100, 1, " % L",
              tr("Beginn des glatten Formübergangs als Prozent der Hornlänge. Vorher gelten "
                 "die H/V-Grundprofile."));
    addNumber(rangeForm, "transitionEnd", tr("Ende"), 65, 0, 100, 1, " % L",
              tr("Ende des glatten Formübergangs als Prozent der Hornlänge. Muss hinter dem Beginn "
                 "liegen."));
    form->addRow(rangeSettings_);
    lutSettings_ = new QWidget;
    auto *lutLayout = new QVBoxLayout(lutSettings_);
    lutLayout->setContentsMargins(0, 0, 0, 0);
    knots_ = new QTableWidget(2, 2);
    knots_->setHorizontalHeaderLabels({tr("Position [% L]"), tr("Form [%]")});
    knots_->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    knots_->verticalHeader()->hide();
    knots_->setMaximumHeight(130);
    knots_->setMinimumHeight(90);
    knots_->setToolTip(tr("Monotone Stützpunkte: Positionen von 0 bis 100 % L, Formanteil 0 bis "
                          "100 %. Endpunkte müssen Hals und Mund abdecken."));
    knots_->setWhatsThis(knots_->toolTip());
    for (int row = 0; row < 2; ++row)
        for (int col = 0; col < 2; ++col)
            knots_->setItem(row, col, new QTableWidgetItem(row == 0 ? "0" : "100"));
    lutLayout->addWidget(knots_);
    auto *buttons = new QHBoxLayout;
    auto *add = new QPushButton(tr("Punkt +"));
    auto *remove = new QPushButton(tr("Punkt −"));
    buttons->addWidget(add);
    buttons->addWidget(remove);
    lutLayout->addLayout(buttons);
    connect(add, &QPushButton::clicked, this, [this] {
        const int row = std::max(1, knots_->currentRow());
        knots_->insertRow(row);
        const double previous = knots_->item(row - 1, 0)->text().toDouble();
        const double next = knots_->item(row + 1, 0)->text().toDouble();
        updating_ = true;
        knots_->setItem(row, 0, new QTableWidgetItem(QString::number((previous + next) * .5)));
        knots_->setItem(
            row, 1,
            new QTableWidgetItem(QString::number((knots_->item(row - 1, 1)->text().toDouble() +
                                                  knots_->item(row + 1, 1)->text().toDouble()) *
                                                 .5)));
        updating_ = false;
        changed();
    });
    connect(remove, &QPushButton::clicked, this, [this] {
        const int row = knots_->currentRow();
        if (row > 0 && row < knots_->rowCount() - 1) {
            knots_->removeRow(row);
            changed();
        }
    });
    connect(knots_, &QTableWidget::cellChanged, this, [this] { changed(); });
    form->addRow(lutSettings_);
    layout->addStretch();
    layout = page(tr("&Leitkurve"), "draftGuideScroll");
    section(tr("Leitkurve"), layout, form);
    guideEnabled_ = new QCheckBox(tr("Profil an Leitkurve ausrichten"));
    guideEnabled_->setObjectName(QStringLiteral("guideEnabled"));
    guideEnabled_->setToolTip(
        tr("Legt einen Zielquerschnitt fest, den die endgültige Hornkontur an Position z "
           "genau trifft. Unter Form wird der Übergang vom Grundprofil zur Leitkurvenform "
           "eingestellt. Bei 100 % L liegt die Leitkurve am Mund."));
    guideEnabled_->setWhatsThis(guideEnabled_->toolTip());
    form->addRow(guideEnabled_);
    guideSettings_ = new QWidget;
    auto *guideForm = new QFormLayout(guideSettings_);
    guideForm->setContentsMargins(0, 0, 0, 0);
    curveType_ = new QComboBox;
    curveType_->setObjectName(QStringLiteral("guideCurveType"));
    curveType_->addItems({tr("Superellipse"), tr("Superformula")});
    guideForm->addRow(tr("Kurve"), curveType_);
    positionMode_ = new QComboBox;
    positionMode_->setObjectName(QStringLiteral("guidePositionMode"));
    positionMode_->addItems({tr("Relativ · % L"), tr("Absolut · mm")});
    guideForm->addRow(tr("Lage"), positionMode_);
    addNumber(guideForm, "guidePosition", tr("Position z"), 100, .001, 100000, 3, " % L",
              tr("Axiale Lage der Leitkurve. Relative Lage in Prozent L oder absolute Lage in mm "
                 "wird ausdrücklich ausgewählt."));
    addNumber(guideForm, "guideWidth", tr("Breite"), 300, .1, 10000, 2, " mm",
              tr("Gesamte Breite der Leitkurve entlang ihrer lokalen X-Achse vor der Rotation. "
                 "Einheit: mm."));
    addNumber(guideForm, "guideAspect", tr("Höhe / Breite"), .7, .01, 100, 3, "",
              tr("Physisches Seitenverhältnis der Leitkurve. Unabhängig von den Formparametern a "
                 "und b."));
    addNumber(guideForm, "guideRotation", tr("Rotation"), 0, -360, 360, 2, " °",
              tr("Rotation der Leitkurve gegen den Uhrzeigersinn um +Z. Einheit: Grad."));
    superellipseSettings_ = new QWidget;
    auto *ellipseForm = new QFormLayout(superellipseSettings_);
    ellipseForm->setContentsMargins(0, 0, 0, 0);
    addNumber(ellipseForm, "guideExponent", tr("Exponent n"), 4, 2, 100, 2, "",
              tr("Superellipsen-Exponent der Leitkurve: 2 ist elliptisch, höhere Werte "
                 "sind rechteckiger. Bei Rotation 0 bleiben die H/V-Achsenmaße gleich; "
                 "die Änderung betrifft die Ecken und ist in der 3D-Ansicht sichtbar."));
    guideForm->addRow(superellipseSettings_);
    superformulaSettings_ = new QWidget;
    auto *sfForm = new QFormLayout(superformulaSettings_);
    sfForm->setContentsMargins(0, 0, 0, 0);
    for (const auto &key : {"a", "b", "m1", "m2", "n1", "n2", "n3"}) {
        const QString name = QString::fromLatin1(key);
        const double value = (name == "a" || name == "b") ? 1 : 4;
        const double minimum = name == "n1" ? -128 : (name.startsWith('m') ? 0 : .01);
        addNumber(sfForm, "sf_" + name, name, value, minimum, 128, 3, "",
                  tr("Dimensionsloser Parameter %1 der Superformula. m1 und m2 steuern die "
                     "Winkelfrequenzen; n1 ≠ 0, a/b/n2/n3 > 0. Nur endliche periodische Konturen. "
                     "Bei a = b, m1 = m2 und n2 = n3 = 2 ist die Grundkurve ein Kreis; n1 und "
                     "gemeinsames m ändern dann die Form nicht.")
                      .arg(name));
    }
    guideForm->addRow(superformulaSettings_);
    guideMetrics_ = new QLabel;
    guideMetrics_->setObjectName("guideMetrics");
    guideMetrics_->setWordWrap(true);
    guideMetrics_->setToolTip(tr("Maße im lokalen Koordinatensystem vor Rotation. "
                                 "Der Zielquerschnitt muss hinter dem Spaltanschluss und "
                                 "nach Beginn des Formverlaufs liegen. Mund-Morphing darf "
                                 "erst an oder hinter der Leitkurve beginnen."));
    guideMetrics_->setWhatsThis(guideMetrics_->toolTip());
    guideForm->addRow(guideMetrics_);
    auto *guideFormLink = new QPushButton(tr("Formverlauf einstellen …"));
    guideFormLink->setObjectName("guideFormLink");
    connect(guideFormLink, &QPushButton::clicked, this, [this] {
        sections_->setCurrentIndex(1);
        transitionMethod_->setFocus();
    });
    guideForm->addRow(guideFormLink);
    form->addRow(guideSettings_);
    layout->addStretch();
    layout = page(tr("&Mund"), "draftMouthScroll");
    section(tr("Mund-Morphing"), layout, form);
    mouthEnabled_ = new QCheckBox(tr("Mundkontur anpassen"));
    mouthEnabled_->setObjectName("mouthMorphEnabled");
    mouthEnabled_->setToolTip(tr("Additive C²-Korrektur zur Zielkontur am Mund. Hals, Leitkurve "
                                 "und Mundtangente bleiben erhalten."));
    mouthEnabled_->setWhatsThis(mouthEnabled_->toolTip());
    form->addRow(mouthEnabled_);
    mouthSettings_ = new QWidget;
    auto *mouthForm = new QFormLayout(mouthSettings_);
    mouthForm->setContentsMargins(0, 0, 0, 0);
    addNumber(mouthForm, "mouthStart", tr("Beginn"), 75, 0, 100, 1, " % L",
              tr("Beginn der C²-Mundkorrektur als Prozent der Hornlänge. Bei aktiver Leitkurve "
                 "muss er an oder hinter deren Position liegen. Keine automatische Verschiebung."));
    addNumber(mouthForm, "mouthEnd", tr("Ende"), 100, 0, 100, 1, " % L",
              tr("Ende der Mundkorrektur als Prozent der Hornlänge; muss hinter dem Beginn "
                 "liegen. Danach gilt die vollständige radiale Korrektur zur Mundkontur."));
    mouthCurveType_ = new QComboBox;
    mouthCurveType_->setObjectName("mouthCurveType");
    mouthCurveType_->addItems({tr("Superellipse"), tr("Superformula")});
    mouthCurveType_->setToolTip(
        tr("Eigenständige Zielkontur am Mund, unabhängig von der Leitkurve."));
    mouthCurveType_->setWhatsThis(mouthCurveType_->toolTip());
    mouthForm->addRow(tr("Zielkontur"), mouthCurveType_);
    addNumber(mouthForm, "mouthWidth", tr("Breite"), 300, .1, 10000, 2, " mm",
              tr("Volle Breite der Mund-Zielkontur entlang ihrer lokalen X-Achse vor Rotation. "
                 "Einheit: mm."));
    addNumber(mouthForm, "mouthAspect", tr("Höhe / Breite"), .7, .01, 100, 3, "",
              tr("Physisches Seitenverhältnis der Mund-Zielkontur, unabhängig von den "
                 "Superformula-Formparametern a und b."));
    addNumber(mouthForm, "mouthRotation", tr("Rotation"), 0, -360, 360, 2, " °",
              tr("Rotation der Mund-Zielkontur gegen den Uhrzeigersinn um +Z. Einheit: Grad."));
    mouthSuperellipseSettings_ = new QWidget;
    auto *mouthEllipseForm = new QFormLayout(mouthSuperellipseSettings_);
    mouthEllipseForm->setContentsMargins(0, 0, 0, 0);
    addNumber(mouthEllipseForm, "mouthExponent", tr("Exponent n"), 4, 2, 100, 2, "",
              tr("Superellipsen-Exponent der Mund-Zielkontur; 2 ist elliptisch. Unabhängig von "
                 "Eintritt, Querschnittsübergang und Termination."));
    mouthForm->addRow(mouthSuperellipseSettings_);
    mouthSuperformulaSettings_ = new QWidget;
    auto *mouthSfForm = new QFormLayout(mouthSuperformulaSettings_);
    mouthSfForm->setContentsMargins(0, 0, 0, 0);
    for (const auto &key : {"a", "b", "m1", "m2", "n1", "n2", "n3"}) {
        const QString name = QString::fromLatin1(key);
        const double value = (name == "a" || name == "b") ? 1 : 4;
        const double minimum = name == "n1" ? -128 : (name.startsWith('m') ? 0 : .01);
        addNumber(
            mouthSfForm, "mouthSf_" + name, name, value, minimum, 128, 3, "",
            tr("Dimensionsloser Superformula-Parameter %1 der Mund-Zielkontur. n1 darf "
               "nicht null sein, a/b/n2/n3 sind positiv. m1/m2 "
               "steuern Winkelfrequenzen, n1/n2/n3 die Form. Die Kontur muss geschlossen sein.")
                .arg(name));
    }
    mouthForm->addRow(mouthSuperformulaSettings_);
    form->addRow(mouthSettings_);
    layout->addStretch();
    layout = page(tr("&Spalt"), "draftSlotScroll");
    section(tr("Diffraktionsspalt"), layout, form);
    slotEnabled_ = new QCheckBox(tr("Spalt vor OS-SE aktivieren"));
    slotEnabled_->setObjectName("slotEnabled");
    slotEnabled_->setToolTip(
        tr("Führt den Eintritt über einen Adapter zum Spalt. Das OS-SE-Profil beginnt am "
           "Spaltausgang; die Gesamtlänge L bleibt erhalten. Ausschalten erhält alle Spaltwerte. "
           "Die akustische Wirkung muss simuliert werden."));
    slotEnabled_->setWhatsThis(slotEnabled_->toolTip());
    form->addRow(slotEnabled_);
    slotSettings_ = new QWidget;
    auto *slotForm = new QFormLayout(slotSettings_);
    slotForm->setContentsMargins(0, 0, 0, 0);
    addNumber(slotForm, "slotWidth", tr("Spaltbreite X"), 12, .01, 2000, 3, " mm",
              tr("Volle lichte Spaltbreite entlang X. Eine schmale X-Abmessung bildet einen "
                 "vertikalen Spalt. Einheit: mm; keine garantierte Abstrahlbreite."));
    addNumber(slotForm, "slotHeight", tr("Spalthöhe Y"), 35, .01, 2000, 3, " mm",
              tr("Volle lichte Spalthöhe entlang Y. Einheit: mm."));
    slotContour_ = new QComboBox;
    slotContour_->setObjectName("slotContour");
    slotContour_->addItems({tr("Superellipse"), tr("Rechteck mit Eckradius")});
    slotContour_->setToolTip(tr("Die Superellipse verwendet einen Formexponenten. Das "
                               "abgerundete Rechteck verwendet echte Kreisbögen mit "
                               "freiem Eckradius; beide Definitionen bleiben getrennt gespeichert."));
    slotContour_->setWhatsThis(slotContour_->toolTip());
    slotForm->addRow(tr("Spaltkontur"), slotContour_);
    addNumber(slotForm, "slotCornerRadius", tr("Eckradius"), 2, 0, 10000, 3, " mm",
              tr("Kreisradius der vier Ecken im Spaltquerschnitt. Einheit: mm; 0 mm ist "
                 "ein scharfes Rechteck. Höchstens die halbe kleinere Spaltabmessung. "
                 "Zu große Werte werden gemeldet, nicht automatisch begrenzt. STEP-Export erfordert "
                 "einen positiven Eckradius; 0 mm ist für Vorschau und Netz zulässig."));
    addNumber(slotForm, "slotExponent", tr("Spaltexponent"), 4, 2, 100, 3, "",
              tr("Superellipsen-Exponent des Spalts: 2 ist elliptisch, höhere Werte sind "
                 "rechteckiger. Unabhängig von Eintritt, Zielquerschnitt und Termination."));
    addNumber(slotForm, "slotAdapterLength", tr("Adapterlänge"), 20, .01, 10000, 3, " mm",
              tr("Axiale Länge des glatten Übergangs vom Eintritt bei z = 0 zum Spalt. "
                 "Einheit: mm; muss positiv sein."));
    addNumber(slotForm, "slotStraightLength", tr("Gerader Spalt"), 3, 0, 10000, 3, " mm",
              tr("Axiale Länge des konstanten Spaltquerschnitts nach dem Adapter. 0 mm lässt "
                 "das OS-SE-Profil unmittelbar am Adapterende beginnen."));
    slotExit_ = new QComboBox;
    slotExit_->setObjectName("slotExit");
    slotExit_->addItems({tr("Axiale Glättung"), tr("Kreisbogen mit Radius")});
    slotExit_->setToolTip(tr("Axiale Glättung steuert die Anschlusslänge. Kreisbogen verwendet "
                            "einen echten Radius im radialen Längsschnitt zum Horn hin. "
                            "Der Eckradius des Querschnitts bleibt davon unabhängig."));
    slotExit_->setWhatsThis(slotExit_->toolTip());
    slotForm->addRow(tr("Spaltausgang"), slotExit_);
    addNumber(slotForm, "slotExitRadius", tr("Ausgangsradius"), 2, 0, 10000, 3, " mm",
              tr("Echter Kreisbogenradius im radialen Längsschnitt vom geraden Spalt "
                 "zum OS-SE-Profil. Einheit: mm; 0 mm erhält den ungerundeten Anschluss. "
                 "Kein allgemeiner 3D-Rollkugelradius. Radius und Profil müssen einen "
                 "tangentialen Anschluss innerhalb der verfügbaren Länge zulassen."));
    addNumber(slotForm, "slotTransitionLength", tr("Glättungslänge"), 5, 0, 10000, 3, " mm",
              tr("Axiale Anschlusslänge zum OS-SE-Profil, kein Verrundungsradius. 0 mm erhält "
                 "die OS-SE-Anfangstangente: ein Knick entsteht bei abweichender Tangente. "
                 "Positive Werte "
                 "schließen C2-stetig an. Formübergang und Mund-Morphing müssen dahinter "
                 "beginnen; die Leitkurve muss dahinter liegen."));
    slotMetrics_ = new QLabel;
    slotMetrics_->setObjectName("slotMetrics");
    slotMetrics_->setWordWrap(true);
    slotMetrics_->setToolTip(
        tr("Lichte Fläche der gewählten Spaltkontur. Die Profil-Restlänge ist L minus "
           "Adapterlänge minus gerader Spalt. Beim Kreisbogen wird der tatsächliche "
           "Anschlussbereich durch die Geometrie bestimmt."));
    slotMetrics_->setWhatsThis(slotMetrics_->toolTip());
    slotForm->addRow(slotMetrics_);
    slotPositionHint_ = new QLabel;
    slotPositionHint_->setObjectName("slotPositionHint");
    slotPositionHint_->setWordWrap(true);
    slotForm->addRow(slotPositionHint_);
    slotFormStart_ = new QPushButton(tr("Formbeginn hinter Spalt setzen"));
    slotFormStart_->setObjectName("slotFormStart");
    slotFormStart_->setToolTip(
        tr("Setzt ausschließlich den Beginn des glatten Formbereichs auf das Ende des "
           "Spaltanschlusses. Das Formende muss bereits dahinter liegen. LUT, Leitkurve und "
           "Mund-Morphing bleiben unverändert und werden in ihren Bereichen eingestellt. "
           "Beim Kreisbogen ist das tatsächlich berechnete Anschlussende erforderlich; "
           "während der Berechnung oder bei ungültigem Anschluss bleibt die Aktion gesperrt."));
    slotFormStart_->setWhatsThis(slotFormStart_->toolTip());
    slotForm->addRow(slotFormStart_);
    form->addRow(slotSettings_);
    layout->addStretch();
    connect(slotContour_, &QComboBox::currentIndexChanged, this, [this] {
        updateEnabled();
        changed();
    });
    connect(slotExit_, &QComboBox::currentIndexChanged, this, [this] {
        updateEnabled();
        changed();
    });
    connect(slotEnabled_, &QCheckBox::toggled, this, [this] {
        updateEnabled();
        changed();
    });
    connect(slotFormStart_, &QPushButton::clicked, this, [this] {
        const auto slot = slotParameters();
        if (slot.exit == SlotExit::circularFillet && !slotTransitionEnd_)
            return;
        const double joinEnd = slot.exit == SlotExit::circularFillet
                                   ? *slotTransitionEnd_
                                   : slot.adapterLength + slot.straightLength + slot.transitionLength;
        const double finish = joinEnd / modelNumber("length", model_.length, 1000);
        updating_ = true;
        model_.transition.startFraction = finish;
        numbers_.at("transitionStart")->setValue(finish * 100);
        editedNumbers_.erase("transitionStart");
        updating_ = false;
        changed();
    });
    connect(asymmetricEnabled_, &QCheckBox::toggled, this, [this] {
        updateEnabled();
        changed();
    });
    connect(guideEnabled_, &QCheckBox::toggled, this, [this] {
        updateEnabled();
        changed();
    });
    connect(curveType_, &QComboBox::currentIndexChanged, this, [this] {
        updateEnabled();
        changed();
    });
    connect(mouthEnabled_, &QCheckBox::toggled, this, [this] {
        updateEnabled();
        changed();
    });
    connect(mouthCurveType_, &QComboBox::currentIndexChanged, this, [this] {
        updateEnabled();
        changed();
    });
    connect(transitionMethod_, &QComboBox::currentIndexChanged, this, [this] {
        updateEnabled();
        changed();
    });
    connect(positionMode_, &QComboBox::currentIndexChanged, this, [this](int index) {
        if (!updating_ && index != previousPositionMode_) {
            updating_ = true;
            const double length = modelNumber("length", model_.length, 1000);
            const double old = modelNumber("guidePosition", model_.guide.position,
                                           previousPositionMode_ == 0 ? 100 : 1000);
            model_.guide.position = index == 1 ? old * length : old / length;
            model_.guide.positionMode =
                index == 1 ? GuidePositionMode::absolute : GuidePositionMode::relative;
            editedNumbers_.erase("guidePosition");
            numbers_.at("guidePosition")
                ->setValue(model_.guide.position * (index == 1 ? 1000 : 100));
            updating_ = false;
        }
        previousPositionMode_ = index;
        numbers_.at("guidePosition")->setSuffix(index == 0 ? " % L" : " mm");
        changed();
    });
    setParameters(DraftParameters{});
    updateEnabled();
}

ParameterPanel::~ParameterPanel() {
    // QWidget owns the inputs, but their callbacks depend on our derived-class members.
    for (auto *child : findChildren<QObject *>())
        disconnect(child, nullptr, this, nullptr);
}

void ParameterPanel::changed() {
    if (!updating_) {
        slotTransitionEnd_.reset();
        updateSlotState();
        updateGuideState();
        if (onChanged)
            onChanged();
    }
}
void ParameterPanel::updateEnabled() {
    const bool asymmetric = asymmetricEnabled_->isChecked();
    oppositeProfiles_->setVisible(asymmetric);
    oppositeProfiles_->setEnabled(asymmetric);
    horizontalProfileLabel_->setText(asymmetric ? tr("Rechts") : tr("H"));
    verticalProfileLabel_->setText(asymmetric ? tr("Oben") : tr("V"));
    guideSettings_->setEnabled(guideEnabled_->isChecked());
    guideSettings_->setVisible(guideEnabled_->isChecked());
    guideAngleHint_->setVisible(guideEnabled_->isChecked());
    superellipseSettings_->setVisible(curveType_->currentIndex() == 0);
    superformulaSettings_->setVisible(curveType_->currentIndex() == 1);
    rangeSettings_->setVisible(transitionMethod_->currentIndex() == 0);
    lutSettings_->setVisible(transitionMethod_->currentIndex() == 1);
    numbers_.at("crossExponent")->setEnabled(!guideEnabled_->isChecked());
    guideShapeHint_->setVisible(guideEnabled_->isChecked());
    mouthSettings_->setEnabled(mouthEnabled_->isChecked());
    mouthSettings_->setVisible(mouthEnabled_->isChecked());
    mouthSuperellipseSettings_->setVisible(mouthCurveType_->currentIndex() == 0);
    mouthSuperformulaSettings_->setVisible(mouthCurveType_->currentIndex() == 1);
    updateSlotState();
    updateGuideState();
}

void ParameterPanel::setGeometryStatus(const QString &message, bool error) {
    geometryStatus_->setText(message);
    geometryStatus_->setStyleSheet(error ? "color: #a72e2b" : "");
    geometryStatus_->setVisible(!message.isEmpty());
}

void ParameterPanel::updateGuideState() {
    if (!guideMetrics_ || updating_)
        return;
    const double width = modelNumber("guideWidth", model_.guide.width, 1000);
    const double height = width * modelNumber("guideAspect", model_.guide.aspectRatio);
    const double position = modelNumber("guidePosition", model_.guide.position,
                                        positionMode_->currentIndex() == 0 ? 100 : 1000);
    const double z = positionMode_->currentIndex() == 0
                         ? position * modelNumber("length", model_.length, 1000)
                         : position;
    guideMetrics_->setText(tr("Zielquerschnitt: %1 × %2 mm bei z = %3 mm")
                               .arg(width * 1000, 0, 'f', 2)
                               .arg(height * 1000, 0, 'f', 2)
                               .arg(z * 1000, 0, 'f', 2));
}

void ParameterPanel::setSlotTransitionEnd(std::optional<double> position) {
    slotTransitionEnd_ = position;
    updateSlotState();
}

DiffractionSlotParameters ParameterPanel::slotParameters() const {
    auto slot = model_.slot;
    slot.enabled = slotEnabled_->isChecked();
    slot.contour = static_cast<SlotContour>(slotContour_->currentIndex());
    slot.exit = static_cast<SlotExit>(slotExit_->currentIndex());
    slot.cornerRadius = modelNumber("slotCornerRadius", slot.cornerRadius, 1000);
    slot.exitRadius = modelNumber("slotExitRadius", slot.exitRadius, 1000);
    slot.aperture.width = modelNumber("slotWidth", slot.aperture.width, 1000);
    slot.aperture.height = modelNumber("slotHeight", slot.aperture.height, 1000);
    slot.aperture.exponent = modelNumber("slotExponent", slot.aperture.exponent);
    slot.adapterLength = modelNumber("slotAdapterLength", slot.adapterLength, 1000);
    slot.straightLength = modelNumber("slotStraightLength", slot.straightLength, 1000);
    slot.transitionLength = modelNumber("slotTransitionLength", slot.transitionLength, 1000);
    return slot;
}

void ParameterPanel::updateSlotState() {
    if (!slotEnabled_ || updating_)
        return;
    const bool active = slotEnabled_->isChecked();
    slotSettings_->setEnabled(active);
    slotSettings_->setVisible(active);
    const bool rounded =
        slotContour_->currentIndex() == static_cast<int>(SlotContour::roundedRectangle);
    const bool circular = slotExit_->currentIndex() == static_cast<int>(SlotExit::circularFillet);
    numbers_.at("slotCornerRadius")->setEnabled(rounded);
    numbers_.at("slotExponent")->setEnabled(!rounded);
    numbers_.at("slotExitRadius")->setEnabled(circular);
    numbers_.at("slotTransitionLength")->setEnabled(!circular);
    const auto slot = slotParameters();
    const double exit = slot.adapterLength + slot.straightLength;
    const double length = modelNumber("length", model_.length, 1000);
    try {
        slotMetrics_->setText(tr("Spaltfläche: %1 mm² · Profil: %2 mm")
                                  .arg(diffractionSlotArea(slot) * 1e6, 0, 'f', 1)
                                  .arg((length - exit) * 1000, 0, 'f', 3));
        slotMetrics_->setToolTip(slotMetrics_->whatsThis());
    } catch (const std::invalid_argument &error) {
        slotMetrics_->setText(tr("Spaltkontur ungültig"));
        slotMetrics_->setToolTip(QString::fromUtf8(error.what()));
    }
    if (circular && !slotTransitionEnd_) {
        slotFormStart_->setEnabled(false);
        slotPositionHint_->setText(
            tr("Spaltausgang: %1 mm · Anschlussende wird aus dem Kreisradius berechnet.")
                .arg(exit * 1000, 0, 'f', 3));
        return;
    }
    const double finish = (circular ? *slotTransitionEnd_ : exit + slot.transitionLength) / length;
    const double start = modelNumber("transitionStart", model_.transition.startFraction, 100);
    const double end = modelNumber("transitionEnd", model_.transition.endFraction, 100);
    const bool range = transitionMethod_->currentIndex() == 0;
    slotFormStart_->setEnabled(active && range && start < finish && finish < end && finish < 1);
    if (finish >= 1) {
        slotPositionHint_->setText(tr("Adapter und Spaltanschluss sind zu lang für L."));
    } else {
        const QString position = tr("Spaltausgang: %1 mm · Formbeginn ab %2 % L.")
                                     .arg(exit * 1000, 0, 'f', 3)
                                     .arg(finish * 100, 0, 'f', 2);
        slotPositionHint_->setText(
            position + (range ? (start < finish ? tr(" Formbeginn anpassen.") : QString{})
                              : tr(" LUT muss bis dahin Formanteil 0 halten.")));
    }
}

double ParameterPanel::modelNumber(const QString &key, double unchangedValue,
                                   double displayUnitsPerSI) const {
    return editedNumbers_.contains(key) ? numbers_.at(key)->value() / displayUnitsPerSI
                                        : unchangedValue;
}

DraftParameters ParameterPanel::parameters() const {
    const auto value = [this](const QString &key, double &field, double units = 1) {
        field = modelNumber(key, field, units);
    };
    DraftParameters p = model_;
    value("entryWidth", p.entry.width, 1000);
    value("entryHeight", p.entry.height, 1000);
    value("entryExponent", p.entry.exponent);
    p.slot = slotParameters();
    value("length", p.length, 1000);
    value("throatAngle", p.profile.throatAngle, degrees);
    value("k", p.profile.k);
    value("s", p.profile.terminationScale);
    value("q", p.profile.terminationQ);
    value("terminationN", p.profile.terminationExponent);
    value("throatAngleV", p.verticalProfile.throatAngle, degrees);
    value("kV", p.verticalProfile.k);
    value("sV", p.verticalProfile.terminationScale);
    value("qV", p.verticalProfile.terminationQ);
    value("terminationNV", p.verticalProfile.terminationExponent);
    value("angleH", p.horizontalAngle, degrees);
    value("angleV", p.verticalAngle, degrees);
    p.asymmetricEnabled = asymmetricEnabled_->isChecked();
    value("angleL", p.leftAngle, degrees);
    value("angleB", p.bottomAngle, degrees);
    const auto oppositeProfile = [&](const QString &suffix, OsseDirectionalParameters &profile) {
        value("throatAngle" + suffix, profile.throatAngle, degrees);
        value("k" + suffix, profile.k);
        value("s" + suffix, profile.terminationScale);
        value("q" + suffix, profile.terminationQ);
        value("terminationN" + suffix, profile.terminationExponent);
    };
    oppositeProfile("L", p.leftProfile);
    oppositeProfile("B", p.bottomProfile);
    value("thickness", p.wallThickness, 1000);
    value("crossExponent", p.crossSectionExponent);
    p.transition.method = transitionMethod_->currentIndex() == 0 ? TransitionMethod::smoothRange
                                                                 : TransitionMethod::lut;
    value("transitionStart", p.transition.startFraction, 100);
    value("transitionEnd", p.transition.endFraction, 100);
    p.transition.knots.clear();
    for (int row = 0; row < knots_->rowCount(); ++row) {
        bool fractionOk = false, blendOk = false;
        const auto knotValue = [&](int column, bool &valid) {
            const auto *item = knots_->item(row, column);
            if (!item)
                return 0.0;
            if (item->data(Qt::UserRole).isValid() &&
                item->text() == item->data(Qt::UserRole + 1).toString()) {
                valid = true;
                return item->data(Qt::UserRole).toDouble();
            }
            return item->text().toDouble(&valid) * .01;
        };
        const double fraction = knotValue(0, fractionOk), blend = knotValue(1, blendOk);
        if (!fractionOk || !blendOk)
            throw std::invalid_argument("LUT: Position und Form müssen Zahlen sein.");
        p.transition.knots.push_back({fraction, blend});
    }
    p.guideEnabled = guideEnabled_->isChecked();
    p.guide.type =
        curveType_->currentIndex() == 0 ? CurveType::superellipse : CurveType::superformula;
    value("guideWidth", p.guide.width, 1000);
    value("guideAspect", p.guide.aspectRatio);
    value("guideRotation", p.guide.rotation, degrees);
    p.guide.positionMode = positionMode_->currentIndex() == 0 ? GuidePositionMode::relative
                                                              : GuidePositionMode::absolute;
    value("guidePosition", p.guide.position, positionMode_->currentIndex() == 0 ? 100 : 1000);
    value("guideExponent", p.guide.superellipseExponent);
    value("sf_a", p.guide.superformula.a);
    value("sf_b", p.guide.superformula.b);
    value("sf_m1", p.guide.superformula.m1);
    value("sf_m2", p.guide.superformula.m2);
    value("sf_n1", p.guide.superformula.n1);
    value("sf_n2", p.guide.superformula.n2);
    value("sf_n3", p.guide.superformula.n3);
    p.mouthMorph.enabled = mouthEnabled_->isChecked();
    value("mouthStart", p.mouthMorph.startFraction, 100);
    value("mouthEnd", p.mouthMorph.endFraction, 100);
    auto &mouth = p.mouthMorph.target;
    mouth.type =
        mouthCurveType_->currentIndex() == 0 ? CurveType::superellipse : CurveType::superformula;
    value("mouthWidth", mouth.width, 1000);
    value("mouthAspect", mouth.aspectRatio);
    value("mouthRotation", mouth.rotation, degrees);
    value("mouthExponent", mouth.superellipseExponent);
    value("mouthSf_a", mouth.superformula.a);
    value("mouthSf_b", mouth.superformula.b);
    value("mouthSf_m1", mouth.superformula.m1);
    value("mouthSf_m2", mouth.superformula.m2);
    value("mouthSf_n1", mouth.superformula.n1);
    value("mouthSf_n2", mouth.superformula.n2);
    value("mouthSf_n3", mouth.superformula.n3);
    return p;
}

void ParameterPanel::setParameters(const DraftParameters &p) {
    updating_ = true;
    const auto set = [this](const QString &key, double value) {
        auto *spin = numbers_.at(key);
        if (!std::isfinite(value) || value < spin->minimum() || value > spin->maximum()) {
            updating_ = false;
            throw std::invalid_argument(
                ("Wert außerhalb des Eingabebereichs: " + key).toStdString());
        }
        spin->setValue(value);
    };
    set("entryWidth", p.entry.width * 1000);
    set("entryHeight", p.entry.height * 1000);
    set("entryExponent", p.entry.exponent);
    slotEnabled_->setChecked(p.slot.enabled);
    slotContour_->setCurrentIndex(static_cast<int>(p.slot.contour));
    slotExit_->setCurrentIndex(static_cast<int>(p.slot.exit));
    set("slotCornerRadius", p.slot.cornerRadius * 1000);
    set("slotExitRadius", p.slot.exitRadius * 1000);
    set("slotWidth", p.slot.aperture.width * 1000);
    set("slotHeight", p.slot.aperture.height * 1000);
    set("slotExponent", p.slot.aperture.exponent);
    set("slotAdapterLength", p.slot.adapterLength * 1000);
    set("slotStraightLength", p.slot.straightLength * 1000);
    set("slotTransitionLength", p.slot.transitionLength * 1000);
    set("length", p.length * 1000);
    set("throatAngle", p.profile.throatAngle * degrees);
    set("k", p.profile.k);
    set("s", p.profile.terminationScale);
    set("q", p.profile.terminationQ);
    set("terminationN", p.profile.terminationExponent);
    set("throatAngleV", p.verticalProfile.throatAngle * degrees);
    set("kV", p.verticalProfile.k);
    set("sV", p.verticalProfile.terminationScale);
    set("qV", p.verticalProfile.terminationQ);
    set("terminationNV", p.verticalProfile.terminationExponent);
    set("angleH", p.horizontalAngle * degrees);
    set("angleV", p.verticalAngle * degrees);
    asymmetricEnabled_->setChecked(p.asymmetricEnabled);
    set("angleL", p.leftAngle * degrees);
    set("angleB", p.bottomAngle * degrees);
    const auto oppositeProfile = [&](const QString &suffix,
                                     const OsseDirectionalParameters &profile) {
        set("throatAngle" + suffix, profile.throatAngle * degrees);
        set("k" + suffix, profile.k);
        set("s" + suffix, profile.terminationScale);
        set("q" + suffix, profile.terminationQ);
        set("terminationN" + suffix, profile.terminationExponent);
    };
    oppositeProfile("L", p.leftProfile);
    oppositeProfile("B", p.bottomProfile);
    set("thickness", p.wallThickness * 1000);
    set("crossExponent", p.crossSectionExponent);
    set("transitionStart", p.transition.startFraction * 100);
    set("transitionEnd", p.transition.endFraction * 100);
    transitionMethod_->setCurrentIndex(p.transition.method == TransitionMethod::smoothRange ? 0
                                                                                            : 1);
    knots_->setRowCount(static_cast<int>(p.transition.knots.size()));
    for (int row = 0; row < knots_->rowCount(); ++row) {
        const auto &knot = p.transition.knots[static_cast<std::size_t>(row)];
        const auto setKnot = [&](int column, double exact) {
            const auto text = QString::number(exact * 100, 'g', 9);
            auto *item = new QTableWidgetItem(text);
            item->setData(Qt::UserRole, exact);
            item->setData(Qt::UserRole + 1, text);
            knots_->setItem(row, column, item);
        };
        setKnot(0, knot.fraction);
        setKnot(1, knot.blend);
    }
    guideEnabled_->setChecked(p.guideEnabled);
    curveType_->setCurrentIndex(p.guide.type == CurveType::superellipse ? 0 : 1);
    positionMode_->setCurrentIndex(p.guide.positionMode == GuidePositionMode::relative ? 0 : 1);
    set("guidePosition",
        p.guide.position * (p.guide.positionMode == GuidePositionMode::relative ? 100 : 1000));
    set("guideWidth", p.guide.width * 1000);
    set("guideAspect", p.guide.aspectRatio);
    set("guideRotation", p.guide.rotation * degrees);
    set("guideExponent", p.guide.superellipseExponent);
    const auto &sf = p.guide.superformula;
    set("sf_a", sf.a);
    set("sf_b", sf.b);
    set("sf_m1", sf.m1);
    set("sf_m2", sf.m2);
    set("sf_n1", sf.n1);
    set("sf_n2", sf.n2);
    set("sf_n3", sf.n3);
    mouthEnabled_->setChecked(p.mouthMorph.enabled);
    set("mouthStart", p.mouthMorph.startFraction * 100);
    set("mouthEnd", p.mouthMorph.endFraction * 100);
    const auto &mouth = p.mouthMorph.target;
    mouthCurveType_->setCurrentIndex(mouth.type == CurveType::superellipse ? 0 : 1);
    set("mouthWidth", mouth.width * 1000);
    set("mouthAspect", mouth.aspectRatio);
    set("mouthRotation", mouth.rotation * degrees);
    set("mouthExponent", mouth.superellipseExponent);
    const auto &mouthSf = mouth.superformula;
    set("mouthSf_a", mouthSf.a);
    set("mouthSf_b", mouthSf.b);
    set("mouthSf_m1", mouthSf.m1);
    set("mouthSf_m2", mouthSf.m2);
    set("mouthSf_n1", mouthSf.n1);
    set("mouthSf_n2", mouthSf.n2);
    set("mouthSf_n3", mouthSf.n3);
    model_ = p;
    slotTransitionEnd_.reset();
    setGeometryStatus({});
    editedNumbers_.clear();
    updating_ = false;
    updateEnabled();
}

} // namespace hornsculpt
