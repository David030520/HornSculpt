#include "draft_wizard.hpp"
#include "numeric_input.hpp"

#include <QAbstractButton>
#include <QCheckBox>
#include <QComboBox>
#include <QFormLayout>
#include <QLabel>
#include <QPainter>
#include <QSpinBox>
#include <QVBoxLayout>
#include <QWizardPage>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <limits>
#include <stdexcept>

namespace hornsculpt {
namespace {
class ProfilePreview final : public QWidget {
  public:
    explicit ProfilePreview(QWidget *parent = nullptr) : QWidget(parent) {
        setObjectName("wizardGeometryPreview");
        setMinimumSize(210, 160);
    }
    void setGeometry(const DraftParameters &parameters) {
        try {
            DraftGeometry geometry(parameters);
            horizontal_.clear();
            vertical_.clear();
            length_ = geometry.length();
            radius_ = 0;
            for (int i = 0; i <= 64; ++i) {
                const double z = length_ * i / 64;
                const double h = geometry.radius(z, 0);
                const double v = geometry.radius(z, std::numbers::pi / 2);
                horizontal_.push_back(h);
                vertical_.push_back(v);
                radius_ = std::max({radius_, h, v});
            }
        } catch (const std::exception &) {
            horizontal_.clear();
            vertical_.clear();
        }
        update();
    }
    void clear() { horizontal_.clear(); vertical_.clear(); update(); }
  protected:
    void paintEvent(QPaintEvent *) override {
        QPainter painter(this);
        painter.fillRect(rect(), palette().base());
        painter.setPen(palette().mid().color());
        const QRectF plot = rect().adjusted(12, 30, -12, -32);
        painter.drawRect(plot);
        painter.drawText(QRectF(8, 4, width() - 16, 22), tr("Geometrische H/V-Kontur"));
        painter.drawText(QRectF(8, height() - 25, width() - 16, 20), tr("Blau: H · Orange: V"));
        if (horizontal_.empty() || radius_ <= 0 || length_ <= 0) return;
        painter.drawLine(QPointF(plot.left(), plot.center().y()), QPointF(plot.right(), plot.center().y()));
        const double scale = std::min(plot.width() / length_, plot.height() / (2 * radius_));
        for (const auto &[values, color] : {std::pair{&horizontal_, QColor(28, 102, 180)},
                                            std::pair{&vertical_, QColor(205, 112, 33)}}) {
            painter.setPen(QPen(color, 1.5));
            for (int sign : {-1, 1}) {
                QPolygonF line;
                for (std::size_t i = 0; i < values->size(); ++i)
                    line << QPointF(plot.left() + scale * length_ * i / (values->size() - 1),
                                    plot.center().y() - sign * scale * (*values)[i]);
                painter.drawPolyline(line);
            }
        }
    }
  private:
    std::vector<double> horizontal_, vertical_;
    double length_{}, radius_{};
};
} // namespace
DraftParameters makeBasicDraft(const BasicDraftSpecification &s) {
    const auto positive = [](double value) { return std::isfinite(value) && value > 0; };
    if (!positive(s.throatDiameterMm) || !positive(s.lengthMm) || !positive(s.wallThicknessMm))
        throw std::invalid_argument("Halsdurchmesser, Länge und Wandstärke müssen positiv sein.");
    for (const auto angle : {s.horizontalOpeningDegrees, s.verticalOpeningDegrees})
        if (!std::isfinite(angle) || angle < 0 || angle >= 180)
            throw std::invalid_argument("Geometrische Vollwinkel müssen zwischen 0° und unter 180° liegen.");
    if (!std::isfinite(s.crossSectionExponent) || s.crossSectionExponent < 2)
        throw std::invalid_argument("Der Querschnittsexponent muss mindestens 2 betragen.");
    DraftParameters p;
    p.entry = {s.throatDiameterMm / 1000, s.throatDiameterMm / 1000, 2};
    p.length = s.lengthMm / 1000;
    p.horizontalAngle = s.horizontalOpeningDegrees * std::numbers::pi / 360;
    p.verticalAngle = s.verticalOpeningDegrees * std::numbers::pi / 360;
    p.leftAngle = p.horizontalAngle;
    p.bottomAngle = p.verticalAngle;
    p.crossSectionExponent = s.crossSectionExponent;
    p.wallThickness = s.wallThicknessMm / 1000;
    // Validate using the same geometry contract as the live Draft pipeline.
    const DraftGeometry checked(p);
    return checked.parameters();
}

DraftWizard::DraftWizard(QWidget *parent) : QWizard(parent) {
    setObjectName("draftWizard");
    setWindowTitle(tr("Horn-Assistent"));
    setWizardStyle(QWizard::ModernStyle);
    setMinimumSize(580, 400);
    setButtonText(QWizard::BackButton, tr("Zurück"));
    setButtonText(QWizard::NextButton, tr("Weiter"));
    setButtonText(QWizard::CancelButton, tr("Abbrechen"));
    setButtonText(QWizard::FinishButton, tr("Draft übernehmen"));
    button(QWizard::NextButton)->setProperty("primaryAction", true);
    button(QWizard::FinishButton)->setProperty("primaryAction", true);
    const auto page = [this](const QString &title, const QString &subtitle) {
        auto *result = new QWizardPage;
        result->setTitle(title);
        result->setSubTitle(subtitle);
        addPage(result);
        return result;
    };
    const auto number = [](QFormLayout *form, const QString &name, const QString &label,
                           double value, double minimum, double maximum, const QString &unit,
                           const QString &help) {
        auto *input = new DecimalInput;
        input->setObjectName(name);
        input->setDecimals(2);
        input->setRange(minimum, maximum);
        input->setSuffix(unit);
        input->setValue(value);
        input->setToolTip(help);
        input->setWhatsThis(help);
        form->addRow(label, input);
        return input;
    };
    auto *entry = page(tr("1 · Anschluss und Bauraum"), tr("Freier Anschlussdurchmesser und verfügbare Abmessungen."));
    auto *extentForm = new QFormLayout(entry);
    diameter_ = number(extentForm, "wizardThroatDiameter", tr("Halsdurchmesser"),
                       25.4, 1, 1000, tr(" mm"), tr("Voller Innendurchmesser bei z = 0. Der Hals bleibt kreisförmig."));
    length_ = number(extentForm, "wizardLength", tr("Hornlänge"), 180, 1, 10000, tr(" mm"),
                     tr("Axiale Länge vom Hals bei z = 0 zum Mund bei z = L."));
    wall_ = number(extentForm, "wizardWallThickness", tr("Wandstärke"), 6, .1, 100, tr(" mm"),
                   tr("Radiale Materialstärke der Hornwand; keine akustische Optimierung."));
    auto *angles = page(tr("2 · Öffnung und Form"), tr("Volle Öffnungswinkel des OS-SE-Grundprofils. Die akustische Abstrahlung wird später berechnet."));
    auto *angleForm = new QFormLayout(angles);
    const auto angleHelp = tr("Voller geometrischer Öffnungswinkel 2α. Intern wird der halbe Winkel α zur +Z-Achse verwendet. Kein garantierter akustischer Abstrahlwinkel.");
    horizontal_ = number(angleForm, "wizardHorizontalOpening", tr("Horizontal · Vollwinkel 2α"),
                         80, 0, 178, tr(" °"), angleHelp);
    vertical_ = number(angleForm, "wizardVerticalOpening", tr("Vertikal · Vollwinkel 2α"),
                       50, 0, 178, tr(" °"), angleHelp);
    contour_ = new QComboBox;
    contour_->setObjectName("wizardContour");
    contour_->addItem(tr("Ellipse · Exponent 2"), 2.0);
    contour_->addItem(tr("Superellipse · Exponent 4"), 4.0);
    contour_->addItem(tr("Superellipse · Exponent 6"), 6.0);
    contour_->setToolTip(tr("H/V-Winkel bestimmen die Abmessungen. Der Formanteil steigt glatt von 0 % bis 65 % der Hornlänge an. Gleiche H/V-Winkel ergeben bei Exponent 2 einen Kreis."));
    contour_->setWhatsThis(contour_->toolTip());
    angleForm->addRow(tr("Querschnitt"), contour_);
    auto *finish = page(tr("3 · Entwurf prüfen"), tr("Die Übernahme ersetzt die Draft-Grundgeometrie in einem rückgängig machbaren Schritt."));
    auto *finishLayout = new QVBoxLayout(finish);
    summary_ = new QLabel;
    summary_->setObjectName("wizardSummary");
    summary_->setWordWrap(true);
    summary_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    finishLayout->addWidget(summary_);
    validation_ = new QLabel;
    validation_->setWordWrap(true);
    validation_->setObjectName("wizardValidation");
    finishLayout->addWidget(validation_);
    finishLayout->addStretch();
    auto *calculation = page(tr("4 · Netz und Berechnung"),
        tr("Wählen Sie eine Netzvorgabe und prüfen Sie Frequenzbereich und Anregung."));
    auto *calculationForm = new QFormLayout(calculation);
    meshPreset_ = new QComboBox;
    meshPreset_->setObjectName("wizardMeshPreset");
    meshPreset_->addItems({tr("Vorschau"), tr("Entwurf"), tr("Referenz")});
    meshPreset_->setCurrentIndex(1);
    meshPreset_->setToolTip(tr("Verwendet die gleichnamige Netzvorgabe; vorhandene Qualitätsgrenzen bleiben erhalten."));
    calculationForm->addRow(tr("Netzvorgabe"), meshPreset_);
    minimumFrequency_ = number(calculationForm, "wizardMinimumFrequency", tr("Frequenz von"),
                               300, 1, 100000, tr(" Hz"), tr("Untere Frequenz des logarithmischen Rasters."));
    maximumFrequency_ = number(calculationForm, "wizardMaximumFrequency", tr("Frequenz bis"),
                               10000, 1, 100000, tr(" Hz"), tr("Obere Frequenz des logarithmischen Rasters."));
    frequencyCount_ = new QSpinBox;
    frequencyCount_->setObjectName("wizardFrequencyCount");
    frequencyCount_->setRange(1, std::numeric_limits<int>::max());
    frequencyCount_->setValue(17);
    calculationForm->addRow(tr("Frequenzpunkte"), frequencyCount_);
    auto *note = new QLabel(tr("Dieser Ablauf verwendet den 3D-BEM-Modus. Der Hals muss als Driven mit gültiger Anregung zugeordnet sein. "
        "Netzgröße und Speicherbedarf hängen von Geometrie und Frequenzbereich ab; "
        "eine verlässliche Laufzeit kann vorab nicht zugesagt werden."));
    note->setWordWrap(true);
    calculationForm->addRow(note);
    startCalculation_ = new QCheckBox(tr("Netz und Berechnung nach Übernahme starten"));
    startCalculation_->setObjectName("wizardStartCalculation");
    startCalculation_->setChecked(true);
    setButtonText(QWizard::FinishButton, tr("Übernehmen und berechnen"));
    connect(startCalculation_, &QCheckBox::toggled, this, [this](bool checked) {
        setButtonText(QWizard::FinishButton, checked ? tr("Übernehmen und berechnen")
                                                   : tr("Draft übernehmen"));
    });
    calculationForm->addRow(startCalculation_);
    calculationValidation_ = new QLabel;
    calculationValidation_->setWordWrap(true);
    calculationForm->addRow(calculationValidation_);
    auto *preview = new ProfilePreview;
    setSideWidget(preview);
    const auto updatePreview = [this, preview] {
        try { preview->setGeometry(parameters()); }
        catch (const std::exception &) { preview->clear(); }
    };
    for (auto *input : {diameter_, length_, wall_, horizontal_, vertical_})
        connect(input, &QDoubleSpinBox::valueChanged, this, [updatePreview] { updatePreview(); });
    connect(contour_, &QComboBox::currentIndexChanged, this, [updatePreview] { updatePreview(); });
    updatePreview();
}

BasicDraftSpecification DraftWizard::specification() const {
    return {diameter_->value(), length_->value(), horizontal_->value(), vertical_->value(),
            contour_->currentData().toDouble(), wall_->value()};
}

DraftParameters DraftWizard::parameters() const { return makeBasicDraft(specification()); }
bool DraftWizard::startCalculation() const { return startCalculation_->isChecked(); }
int DraftWizard::meshPreset() const { return meshPreset_->currentIndex(); }
double DraftWizard::minimumFrequency() const { return minimumFrequency_->value(); }
double DraftWizard::maximumFrequency() const { return maximumFrequency_->value(); }
int DraftWizard::frequencyCount() const { return frequencyCount_->value(); }

void DraftWizard::initializePage(int id) {
    QWizard::initializePage(id);
    if (id != 2)
        return;
    const auto s = specification();
    const auto number = [this](double value) { return locale().toString(value, 'g', 12); };
    summary_->setText(tr("Anschluss Ø %1 mm · Länge %2 mm · Wand %3 mm\n"
                         "Geometrische Vollwinkel H %4° / V %5°\n"
                         "Halbwinkel zur Achse H %6° / V %7°\n"
                         "OS-SE · Querschnittsexponent %8\n\n"
                         "Leitkurve, Mund-Morphing, Asymmetrie und Diffraktionsspalt sind ausgeschaltet. "
                         "Weitere Einstellungen finden Sie anschließend rechts im Parameterdock.")
                          .arg(number(s.throatDiameterMm)).arg(number(s.lengthMm)).arg(number(s.wallThicknessMm))
                          .arg(number(s.horizontalOpeningDegrees)).arg(number(s.verticalOpeningDegrees))
                          .arg(number(s.horizontalOpeningDegrees / 2)).arg(number(s.verticalOpeningDegrees / 2))
                          .arg(number(s.crossSectionExponent)));
    try {
        const DraftGeometry geometry(parameters());
        summary_->setText(summary_->text() + tr("\nMund innen: %1 × %2 mm (H × V)")
            .arg(locale().toString(2000 * geometry.radius(geometry.length(), 0), 'f', 1))
            .arg(locale().toString(2000 * geometry.radius(geometry.length(), std::numbers::pi / 2), 'f', 1)));
        validation_->clear();
    } catch (const std::invalid_argument &error) {
        validation_->setText(tr("Entwurf korrigieren: %1").arg(QString::fromUtf8(error.what())));
    }
}

bool DraftWizard::validateCurrentPage() {
    if (currentId() == 3) {
        if (frequencyCount_->value() > 1 && maximumFrequency_->value() <= minimumFrequency_->value()) {
            calculationValidation_->setText(tr("Die obere Frequenz muss größer als die untere Frequenz sein."));
            return false;
        }
        calculationValidation_->clear();
        return true;
    }
    if (currentId() != 2)
        return QWizard::validateCurrentPage();
    try {
        (void)parameters();
        validation_->clear();
        return true;
    } catch (const std::invalid_argument &error) {
        validation_->setText(tr("Entwurf korrigieren: %1").arg(QString::fromUtf8(error.what())));
        return false;
    }
}
} // namespace hornsculpt
