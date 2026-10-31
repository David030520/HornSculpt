#include "acoustic_field_widget.hpp"
#include "numeric_input.hpp"
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QMessageBox>
#include <QPushButton>
#include <QSaveFile>
#include <QScrollArea>
#include <QSignalBlocker>
#include <QSlider>
#include <QSpinBox>
#include <QTimer>
#include <QToolButton>
#include <QVBoxLayout>
#include <algorithm>
#include <cmath>
#include <hornsculpt/pressure_field_plane.hpp>
#include <new>

namespace hornsculpt {
namespace {
Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3 mul(Vec3 a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}
FieldDisplaySettings cacheKey(FieldDisplaySettings s) {
    s.quantity = 0;
    s.minimumDb = 40;
    s.maximumDb = 110;
    s.phaseDegrees = 0;
    return s;
}
// Only data used by field evaluation and its near-quadrature settings belong here.
bool sameFieldRequest(const std::optional<BemRequest> &a, const std::optional<BemRequest> &b) {
    if (a.has_value() != b.has_value())
        return false;
    if (!a)
        return true;
    if (a->medium.density != b->medium.density || a->medium.soundSpeed != b->medium.soundSpeed ||
        a->numerics.nearRelativeTolerance != b->numerics.nearRelativeTolerance ||
        a->numerics.nearMaximumDepth != b->numerics.nearMaximumDepth ||
        a->axisymmetric != b->axisymmetric || a->infiniteBaffle != b->infiniteBaffle ||
        a->subdomains != b->subdomains || a->symmetry != b->symmetry ||
        a->exteriorBaffle != b->exteriorBaffle || a->observerValidity != b->observerValidity ||
        a->coupling.has_value() != b->coupling.has_value() ||
        (a->coupling && !lem::samePhysics(*a->coupling, *b->coupling)) ||
        a->excitations.size() != b->excitations.size())
        return false;
    for (std::size_t i = 0; i < a->excitations.size(); ++i) {
        const auto &x = a->excitations[i], &y = b->excitations[i];
        if (x.name != y.name || x.faces.size() != y.faces.size())
            return false;
        for (std::size_t j = 0; j < x.faces.size(); ++j)
            if (x.faces[j].faceId != y.faces[j].faceId ||
                x.faces[j].normalVelocity != y.faces[j].normalVelocity)
                return false;
    }
    return true;
}
QDoubleSpinBox *number(const char *name, double low, double high, int decimals,
                       const QString &suffix) {
    auto *spin = new DecimalInput;
    spin->setObjectName(name);
    spin->setRange(low, high);
    spin->setDecimals(decimals);
    spin->setSuffix(suffix);
    return spin;
}
} // namespace
AcousticFieldWidget::AcousticFieldWidget(QWidget *parent) : QWidget(parent) {
    setObjectName("acousticField");
    auto *layout = new QGridLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    plane_ = new QComboBox;
    plane_->setObjectName("fieldPlane");
    plane_->addItems({tr("H-Schnitt"), tr("V-Schnitt")});
    frequency_ = new QComboBox;
    frequency_->setObjectName("fieldFrequency");
    excitation_ = new QComboBox;
    excitation_->setObjectName("fieldExcitation");
    quantity_ = new QComboBox;
    quantity_->setObjectName("fieldQuantity");
    quantity_->addItems({tr("Pegel"), tr("Phase"), tr("Re(p·eⁱᵠ)")});
    evaluate_ = new QPushButton(tr("Feld berechnen"));
    evaluate_->setObjectName("evaluateField");
    evaluate_->setProperty("primaryAction", true);
    cancel_ = new QPushButton(tr("Abbrechen"));
    cancel_->setObjectName("cancelField");
    export_ = new QPushButton(tr("CSV …"));
    export_->setObjectName("exportField");
    auto *selectionGroup = new QWidget;
    auto *selection = new QGridLayout(selectionGroup);
    selection->setContentsMargins(0, 0, 0, 0);
    const auto addField = [](QGridLayout *grid, int row, int column, const QString &text,
                             QWidget *control) {
        auto *label = new QLabel(text);
        label->setBuddy(control);
        grid->addWidget(label, row, column);
        grid->addWidget(control, row, column + 1);
    };
    addField(selection, 0, 0, tr("Schnitt"), plane_);
    addField(selection, 0, 2, tr("Frequenz"), frequency_);
    addField(selection, 1, 0, tr("Anregung"), excitation_);
    addField(selection, 1, 2, tr("Darstellung"), quantity_);
    selection->setColumnStretch(1, 1);
    selection->setColumnStretch(3, 1);
    layout->addWidget(selectionGroup, 0, 0);
    auto *frequencyNavigation = new QHBoxLayout;
    previousFrequency_ = new QToolButton;
    previousFrequency_->setObjectName("fieldPreviousFrequency");
    previousFrequency_->setText(tr("◀"));
    previousFrequency_->setToolTip(tr("Vorherige berechnete Frequenz"));
    frequencySlider_ = new QSlider(Qt::Horizontal);
    frequencySlider_->setObjectName("fieldFrequencySlider");
    frequencySlider_->setAccessibleName(tr("Berechnete Frequenz"));
    frequencySlider_->setToolTip(tr("Nur gespeicherte Frequenzen auswählen"));
    nextFrequency_ = new QToolButton;
    nextFrequency_->setObjectName("fieldNextFrequency");
    nextFrequency_->setText(tr("▶"));
    nextFrequency_->setToolTip(tr("Nächste berechnete Frequenz"));
    frequencyNavigation->addWidget(previousFrequency_);
    frequencyNavigation->addWidget(frequencySlider_, 1);
    frequencyNavigation->addWidget(nextFrequency_);
    layout->addLayout(frequencyNavigation, 1, 0);
    minimumLateral_ = number("fieldMinimumLateral", -100000, 99999, 3, " mm");
    maximumLateral_ = number("fieldMaximumLateral", -99999, 100000, 3, " mm");
    minimumAxial_ = number("fieldMinimumAxial", -100000, 99999, 3, " mm");
    maximumAxial_ = number("fieldMaximumAxial", -99999, 100000, 3, " mm");
    columns_ = new IntegerInput;
    columns_->setObjectName("fieldColumns");
    columns_->setRange(3, 201);
    rows_ = new IntegerInput;
    rows_->setObjectName("fieldRows");
    rows_->setRange(3, 201);
    auto *details = new QWidget;
    auto *detailsLayout = new QVBoxLayout(details);
    detailsLayout->setContentsMargins(0, 0, 0, 0);
    auto *rangeGroup = new QGroupBox(tr("Bereich und Raster"));
    rangeGroup->setProperty("fieldSection", true);
    auto *range = new QGridLayout(rangeGroup);
    addField(range, 0, 0, tr("Quer von"), minimumLateral_);
    addField(range, 0, 2, tr("bis"), maximumLateral_);
    addField(range, 1, 0, tr("Axial von"), minimumAxial_);
    addField(range, 1, 2, tr("bis"), maximumAxial_);
    addField(range, 2, 0, tr("Spalten"), columns_);
    addField(range, 2, 2, tr("Zeilen"), rows_);
    exclusion_ = number("fieldExclusion", .001, 10, 3, " mm");
    addField(range, 3, 0, tr("Randabstand"), exclusion_);
    range->setColumnStretch(1, 1);
    range->setColumnStretch(3, 1);
    detailsLayout->addWidget(rangeGroup);
    minimumDb_ = number("fieldMinimumDb", -300, 299, 1, " dB");
    maximumDb_ = number("fieldMaximumDb", -299, 300, 1, " dB");
    phase_ = number("fieldPhase", -180, 180, 1, " °");
    auto *displayGroup = new QGroupBox(tr("Farbskala und Darstellung"));
    displayGroup->setProperty("fieldSection", true);
    auto *display = new QGridLayout(displayGroup);
    addField(display, 0, 0, tr("Pegel von"), minimumDb_);
    addField(display, 0, 2, tr("bis"), maximumDb_);
    addField(display, 1, 0, tr("Phase φ"), phase_);
    animate_ = new QPushButton(tr("Animation starten"));
    animate_->setObjectName("fieldAnimate");
    animate_->setToolTip(tr("Zeigt den rotierenden Realteil aus dem gespeicherten komplexen Feld; "
                            "keine neue Berechnung. Die eingestellte Phase bleibt erhalten."));
    display->addWidget(animate_, 2, 0, 1, 4);
    animationTimer_ = new QTimer(this);
    animationTimer_->setInterval(40);
    connect(animate_, &QPushButton::clicked, this, [this] {
        if (animationTimer_->isActive()) {
            animationTimer_->stop();
            animate_->setText(tr("Animation starten"));
            publishField();
            return;
        }
        if (!field_)
            return;
        quantity_->setCurrentIndex(2);
        animationPhase_ = settings_.phaseDegrees;
        animationTimer_->start();
        animate_->setText(tr("Animation pausieren"));
    });
    connect(animationTimer_, &QTimer::timeout, this, [this] {
        animationPhase_ = std::fmod(animationPhase_ + 6. + 180., 360.) - 180.;
        publishField();
    });
    display->setColumnStretch(1, 1);
    display->setColumnStretch(3, 1);
    detailsLayout->addWidget(displayGroup);
    detailsLayout->addStretch();
    auto *detailsScroll = new QScrollArea;
    detailsScroll->setObjectName("fieldDetailsScroll");
    detailsScroll->setFrameShape(QFrame::NoFrame);
    detailsScroll->setWidgetResizable(true);
    detailsScroll->setWidget(details);
    auto *detailsDialog = new QDialog(this);
    detailsDialog->setObjectName("fieldDetailsDialog");
    detailsDialog->setWindowTitle(tr("Schallfeld-Einstellungen"));
    detailsDialog->resize(520, 380);
    auto *dialogLayout = new QVBoxLayout(detailsDialog);
    dialogLayout->addWidget(detailsScroll);
    auto *dialogButtons = new QDialogButtonBox(QDialogButtonBox::Close);
    dialogLayout->addWidget(dialogButtons);
    connect(dialogButtons, &QDialogButtonBox::rejected, detailsDialog, &QDialog::close);
    auto *detailsToggle = new QToolButton;
    detailsToggle->setObjectName("fieldDetailsToggle");
    detailsToggle->setText(tr("Optionen"));
    detailsToggle->setAccessibleName(tr("Schallfeld-Einstellungen"));
    detailsToggle->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    detailsToggle->setToolTip(tr("Bereich, Raster und Farbskala bearbeiten."));
    connect(detailsToggle, &QToolButton::clicked, this, [detailsDialog] {
        detailsDialog->show();
        detailsDialog->raise();
    });
    auto *actions = new QHBoxLayout;
    actions->addWidget(evaluate_);
    actions->addWidget(cancel_);
    actions->addWidget(detailsToggle);
    actions->addStretch();
    auto *displayActions = new QHBoxLayout;
    displayActions->addWidget(export_);
    auto *probe = probe_ = new QToolButton;
    probe->setObjectName("fieldProbeMode");
    probe->setText(tr("Messpunkt"));
    probe->setCheckable(true);
    probe->setToolTip(tr("In der 3D-Ansicht einen berechneten Rasterpunkt auswählen; "
                         "kein interpolierter Messwert."));
    displayActions->addWidget(probe);
    connect(probe, &QToolButton::toggled, this, [this](bool enabled) {
        if (onProbeModeChanged)
            onProbeModeChanged(enabled);
    });
    focus_ = new QPushButton(tr("3D einpassen"));
    focus_->setObjectName("focusField3D");
    focus_->setToolTip(tr("Die berechnete Feldebene im 3D-Modell einpassen."));
    focus_->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    displayActions->addWidget(focus_);
    displayActions->addStretch();
    layout->addLayout(actions, 2, 0);
    layout->addLayout(displayActions, 3, 0);
    status_ = new QLabel;
    status_->setObjectName("fieldStatus");
    status_->setWordWrap(true);
    layout->addWidget(status_, 4, 0);
    for (auto *control : {static_cast<QWidget *>(plane_), static_cast<QWidget *>(frequency_),
                          static_cast<QWidget *>(excitation_), static_cast<QWidget *>(quantity_),
                          static_cast<QWidget *>(evaluate_), static_cast<QWidget *>(cancel_),
                          static_cast<QWidget *>(export_), static_cast<QWidget *>(focus_),
                          static_cast<QWidget *>(detailsToggle)}) {
        auto policy = control->sizePolicy();
        policy.setVerticalPolicy(QSizePolicy::Fixed);
        control->setSizePolicy(policy);
    }
    const auto rangeHelp =
        tr("Schnitt relativ zum gespeicherten Beobachtungsrahmen des Ergebnisses, "
           "unabhängig von der Kamera. H/V sind die jeweiligen Querachsen; axial "
           "folgt der 0°-Richtung. Endliche Beobachtungspunkte, kein Fernfeldbild. "
           "Mindestens etwa sechs Punkte je Wellenlänge wählen; enge Spalte "
           "benötigen ein dichteres Raster. Änderungen benötigen Feld berechnen.");
    for (auto *w : {minimumLateral_, maximumLateral_, minimumAxial_, maximumAxial_}) {
        w->setToolTip(rangeHelp);
        w->setWhatsThis(rangeHelp);
    }
    for (auto *w : {columns_, rows_}) {
        w->setToolTip(rangeHelp);
        w->setWhatsThis(rangeHelp);
    }
    phase_->setToolTip(tr("Rotiert den komplexen RMS-Druck für Re(p·exp(iφ)). "
                          "Phasenlage des harmonischen Feldes; kein zusätzlicher Solverlauf. "
                          "Spitzenwerte sind √2-mal so groß."));
    phase_->setWhatsThis(phase_->toolTip());
    exclusion_->setToolTip(
        tr("Punkte innerhalb dieses Abstands zur BEM-Oberfläche bleiben ausgespart. "
           "Keine Grenzwert- oder Singularauswertung. Einheit mm."));
    exclusion_->setWhatsThis(exclusion_->toolTip());
    quantity_->setToolTip(
        tr("Pegel: absoluter RMS-Schalldruck, Bezug 20 µPa, in 3-dB-Stufen ab der Obergrenze. "
           "Phase: zehn zyklische 36°-Sektoren. Re(p·exp(iφ)): zehn Stufen über dem "
           "symmetrischen Referenzdruck in Pa RMS; Null liegt an einer Stufengrenze. "
           "Änderungen der Darstellung benötigen keine neue Feldberechnung. "
           "Ausgesparte Punkte und angrenzende Dreiecke bleiben transparent."));
    quantity_->setWhatsThis(quantity_->toolTip());
    connect(evaluate_, &QPushButton::clicked, this, [this] { requestEvaluation(); });
    connect(cancel_, &QPushButton::clicked, this, [this] {
        cancel();
        status_->setText(tr("Feldberechnung abgebrochen."));
    });
    connect(export_, &QPushButton::clicked, this, [this] { exportSamples(); });
    connect(focus_, &QPushButton::clicked, this, [this] {
        if (onFocusField)
            onFocusField();
    });
    for (auto *w : {minimumLateral_, maximumLateral_, minimumAxial_, maximumAxial_, exclusion_})
        connect(w, &QDoubleSpinBox::valueChanged, this, [this] { controlsChanged(true); });
    for (auto *w : {minimumDb_, maximumDb_, phase_})
        connect(w, &QDoubleSpinBox::valueChanged, this, [this] { controlsChanged(false); });
    for (auto *w : {columns_, rows_})
        connect(w, &QSpinBox::valueChanged, this, [this] { controlsChanged(true); });
    connect(plane_, &QComboBox::currentIndexChanged, this, [this] { controlsChanged(true); });
    connect(quantity_, &QComboBox::currentIndexChanged, this, [this] { controlsChanged(false); });
    connect(previousFrequency_, &QToolButton::clicked, this,
            [this] { frequency_->setCurrentIndex(std::max(0, frequency_->currentIndex() - 1)); });
    connect(nextFrequency_, &QToolButton::clicked, this, [this] {
        frequency_->setCurrentIndex(
            std::min(frequency_->count() - 1, frequency_->currentIndex() + 1));
    });
    connect(frequencySlider_, &QSlider::valueChanged, frequency_, &QComboBox::setCurrentIndex);
    connect(frequency_, &QComboBox::currentIndexChanged, this, [this] {
        cancel();
        clearField();
        settings_.frequencyHz = frequency_->currentData().toDouble();
        {
            const QSignalBlocker blocker(frequencySlider_);
            frequencySlider_->setValue(std::max(0, frequency_->currentIndex()));
        }
        restoreCachedField();
        updateAvailability();
        if (onChanged)
            onChanged();
    });
    connect(excitation_, &QComboBox::currentIndexChanged, this, [this] {
        cancel();
        clearField();
        settings_.excitationIndex = excitation_->currentIndex();
        restoreCachedField();
        updateAvailability();
        if (onChanged)
            onChanged();
    });
    setSettings({});
    updateAvailability();
}
AcousticFieldWidget::~AcousticFieldWidget() {
    for (auto *child : findChildren<QObject *>())
        QObject::disconnect(child, nullptr, this, nullptr);
    onChanged = {};
    onFieldChanged = {};
    onFocusField = {};
    onProbeModeChanged = {};
    cancel();
    // Join while the receiver still exists; interactive cancellation never waits.
    workers_.clear();
}
void AcousticFieldWidget::showEvent(QShowEvent *event) {
    QWidget::showEvent(event);
    for (const auto *name :
         {"fieldPlane", "fieldFrequency", "fieldExcitation", "fieldQuantity", "evaluateField",
          "cancelField", "exportField", "focusField3D", "fieldDetailsToggle"}) {
        auto *control = findChild<QWidget *>(QString::fromLatin1(name));
        control->ensurePolished();
        control->setMinimumHeight(control->minimumSizeHint().height());
    }
    auto *selectionGroup = plane_->parentWidget();
    selectionGroup->setMinimumHeight(selectionGroup->minimumSizeHint().height());
}
void AcousticFieldWidget::hideEvent(QHideEvent *event) {
    if (animationTimer_->isActive()) {
        animationTimer_->stop();
        animate_->setText(tr("Animation starten"));
        publishField();
    }
    if (auto *probe = findChild<QToolButton *>("fieldProbeMode"))
        probe->setChecked(false);
    QWidget::hideEvent(event);
}
void AcousticFieldWidget::cancel() {
    const bool wasBusy = busy_;
    ++jobId_;
    for (auto &worker : workers_)
        worker.thread.request_stop();
    busy_ = false;
    if (wasBusy)
        error_ = tr("Feldberechnung abgebrochen.");
    if (wasBusy)
        clearField();
    if (cancel_)
        updateAvailability();
}
void AcousticFieldWidget::setSettings(const FieldDisplaySettings &s) {
    validateFieldDisplaySettings(s);
    auto sampling = s;
    sampling.quantity = settings_.quantity;
    sampling.minimumDb = settings_.minimumDb;
    sampling.maximumDb = settings_.maximumDb;
    sampling.phaseDegrees = settings_.phaseDegrees;
    if (sampling != settings_) {
        cancel();
        clearField();
    }
    const QSignalBlocker a(plane_), b(quantity_), c(minimumLateral_), d(maximumLateral_),
        e(minimumAxial_), f(maximumAxial_), g(columns_), h(rows_), i(minimumDb_), j(maximumDb_),
        k(phase_), l(exclusion_);
    settings_ = s;
    plane_->setCurrentIndex(s.vertical ? 1 : 0);
    quantity_->setCurrentIndex(s.quantity);
    minimumLateral_->setValue(s.minimumLateral * 1000);
    maximumLateral_->setValue(s.maximumLateral * 1000);
    minimumAxial_->setValue(s.minimumAxial * 1000);
    maximumAxial_->setValue(s.maximumAxial * 1000);
    columns_->setValue(static_cast<int>(s.lateralCount));
    rows_->setValue(static_cast<int>(s.axialCount));
    minimumDb_->setValue(s.minimumDb);
    maximumDb_->setValue(s.maximumDb);
    phase_->setValue(s.phaseDegrees);
    exclusion_->setValue(s.exclusionDistance * 1000);
    publishField();
    phase_->setEnabled(s.quantity == 2);
    minimumDb_->setEnabled(s.quantity == 0);
    maximumDb_->setEnabled(s.quantity != 1);
}
void AcousticFieldWidget::controlsChanged(bool sampling) {
    auto selected = settings_;
    selected.vertical = plane_->currentIndex() == 1;
    selected.quantity = quantity_->currentIndex();
    selected.minimumLateral = minimumLateral_->value() / 1000;
    selected.maximumLateral = maximumLateral_->value() / 1000;
    selected.minimumAxial = minimumAxial_->value() / 1000;
    selected.maximumAxial = maximumAxial_->value() / 1000;
    selected.lateralCount = static_cast<unsigned>(columns_->value());
    selected.axialCount = static_cast<unsigned>(rows_->value());
    selected.minimumDb = minimumDb_->value();
    selected.maximumDb = maximumDb_->value();
    selected.phaseDegrees = phase_->value();
    selected.exclusionDistance = exclusion_->value() / 1000;
    try {
        validateFieldDisplaySettings(selected);
    } catch (const std::invalid_argument &e) {
        status_->setText(QString::fromUtf8(e.what()));
        return;
    }
    if (sampling) {
        cancel();
        clearField();
    }
    settings_ = selected;
    if (selected.quantity != 2 && animationTimer_->isActive()) {
        animationTimer_->stop();
        animate_->setText(tr("Animation starten"));
    }
    if (sampling)
        restoreCachedField();
    publishField();
    phase_->setEnabled(selected.quantity == 2);
    minimumDb_->setEnabled(selected.quantity == 0);
    maximumDb_->setEnabled(selected.quantity != 1);
    updateAvailability();
    if (onChanged)
        onChanged();
}
void AcousticFieldWidget::setSnapshot(std::shared_ptr<const SurfaceMesh> mesh,
                                      std::optional<BemRequest> request,
                                      std::shared_ptr<const BemResult> result, AcousticFrame frame,
                                      std::uint64_t revision, bool current) {
    const auto normalizedFrame = normalizedAcousticFrame(frame);
    const bool changed = mesh_ != mesh || result_ != result || revision_ != revision ||
                         current_ != current || frame_ != normalizedFrame ||
                         !sameFieldRequest(request_, request);
    if (changed) {
        cancel();
        clearField();
        fieldCache_.clear();
    }
    mesh_ = std::move(mesh);
    request_ = std::move(request);
    result_ = std::move(result);
    frame_ = normalizedFrame;
    revision_ = revision;
    current_ = current;
    const auto selected = settings_.frequencyHz;
    const auto selectedExcitation = settings_.excitationIndex;
    const QSignalBlocker a(frequency_), b(excitation_);
    frequency_->clear();
    excitation_->clear();
    if (result_) {
        for (const auto &f : result_->frequencies)
            frequency_->addItem(QString::number(f.frequency, 'g', 7) + " Hz", f.frequency);
        if (!result_->frequencies.empty())
            for (const auto &e : result_->frequencies.front().excitations)
                excitation_->addItem(QString::fromStdString(e.name));
    }
    const int found = selected > 0 ? frequency_->findData(selected) : -1;
    frequency_->setCurrentIndex(found >= 0 ? found : frequency_->count() / 2);
    if (selectedExcitation >= 0 && selectedExcitation < excitation_->count())
        excitation_->setCurrentIndex(selectedExcitation);
    settings_.frequencyHz =
        frequency_->currentIndex() >= 0 ? frequency_->currentData().toDouble() : 0;
    settings_.excitationIndex = excitation_->currentIndex();
    {
        const QSignalBlocker blocker(frequencySlider_);
        frequencySlider_->setRange(0, std::max(0, frequency_->count() - 1));
        frequencySlider_->setValue(std::max(0, frequency_->currentIndex()));
    }
    updateAvailability();
}
bool AcousticFieldWidget::canEvaluate() const noexcept {
    bool ready = mesh_ && request_ && result_ && frequency_->currentIndex() >= 0 &&
                 excitation_->currentIndex() >= 0;
    if (ready) {
        const auto fi = static_cast<std::size_t>(frequency_->currentIndex()),
                   ei = static_cast<std::size_t>(excitation_->currentIndex());
        ready =
            fi < result_->frequencies.size() && ei < result_->frequencies[fi].excitations.size();
        if (ready) {
            const auto &e = result_->frequencies[fi].excitations[ei];
            if (request_->subdomains) {
                ready = e.subdomainTraces().size() == request_->subdomains->domains.size();
                for (const auto &trace : e.subdomainTraces()) {
                    ready = ready && trace.vertexIndices.size() == trace.pressureTrace.size() &&
                            trace.triangleIndices.size() == trace.normalDerivative.size();
                    for (const auto &flux : trace.interfaceFluxes)
                        ready = ready && flux.vertexIndices.size() == flux.normalDerivative.size();
                }
            } else
                ready = e.boundary() &&
                        e.boundary()->pressureTrace.size() == mesh_->vertices.size() &&
                        e.boundary()->normalDerivative.size() == mesh_->triangles.size();
        }
    }
    return ready;
}
bool AcousticFieldWidget::evaluationInProgress() const noexcept {
    return busy_;
}
bool AcousticFieldWidget::fieldReady() const noexcept {
    return !busy_ && error_.isEmpty() && bool(field_);
}
QString AcousticFieldWidget::lastError() const {
    return error_;
}
std::shared_ptr<const PressureFieldPlane> AcousticFieldWidget::fieldPlane() const noexcept {
    return field_;
}
void AcousticFieldWidget::publishField() {
    if (onFieldChanged)
        onFieldChanged(field_,
                       {settings_.quantity, settings_.minimumDb, settings_.maximumDb,
                        animationTimer_->isActive() ? animationPhase_ : settings_.phaseDegrees});
}
void AcousticFieldWidget::clearField() {
    animationTimer_->stop();
    animate_->setText(tr("Animation starten"));
    field_.reset();
    publishField();
}
bool AcousticFieldWidget::restoreCachedField() {
    const auto key = cacheKey(settings_);
    const auto found = std::find_if(fieldCache_.begin(), fieldCache_.end(),
                                    [&](const CachedField &entry) { return entry.key == key; });
    if (found == fieldCache_.end())
        return false;
    fieldCache_.splice(fieldCache_.begin(), fieldCache_, found);
    field_ = fieldCache_.front().plane;
    error_.clear();
    status_->setText(fieldCache_.front().status + tr(" · aus Feldspeicher"));
    publishField();
    return true;
}
void AcousticFieldWidget::updateAvailability() {
    const bool ready = canEvaluate();
    evaluate_->setEnabled(ready && !busy_);
    cancel_->setEnabled(busy_);
    export_->setEnabled(bool(field_) && !busy_);
    focus_->setEnabled(bool(field_) && !busy_);
    if (!field_ || busy_)
        probe_->setChecked(false);
    probe_->setEnabled(bool(field_) && !busy_);
    animate_->setEnabled(bool(field_) && !busy_);
    frequencySlider_->setEnabled(frequency_->count() > 1);
    previousFrequency_->setEnabled(frequency_->currentIndex() > 0);
    nextFrequency_->setEnabled(frequency_->currentIndex() >= 0 &&
                               frequency_->currentIndex() + 1 < frequency_->count());
    if (!ready && !busy_)
        status_->setText(result_ ? tr("Für dieses Ergebnis fehlen gespeicherte BEM-Randwerte. "
                                      "Simulation neu berechnen.")
                                 : tr("Berechnetes BEM-Ergebnis erforderlich."));
    else if (ready && !busy_ && !field_)
        status_->setText(tr("Feldschnitt bereit · Rev. %1%2")
                             .arg(revision_)
                             .arg(current_ ? QString{} : tr(" · veraltet")));
}
void AcousticFieldWidget::requestEvaluation() {
    cancel();
    error_.clear();
    clearField();
    if (!canEvaluate()) {
        error_ =
            result_ ? tr("Für dieses Ergebnis fehlen vollständige gespeicherte BEM-Randwerte. "
                         "Simulation neu berechnen.")
                    : tr("Für die Feldberechnung ist ein gespeichertes BEM-Ergebnis erforderlich.");
        updateAvailability();
        status_->setText(error_);
        return;
    }
    if (restoreCachedField()) {
        updateAvailability();
        return;
    }
    try {
        const auto mesh = mesh_;
        const auto result = result_;
        auto request = *request_;
        const auto frame = frame_;
        const auto settings = settings_;
        const auto fi = static_cast<std::size_t>(frequency_->currentIndex()),
                   ei = static_cast<std::size_t>(excitation_->currentIndex());
        const auto id = ++jobId_, revision = revision_;
        const bool current = current_;
        std::erase_if(workers_, [](const FieldWorker &worker) { return worker.finished->load(); });
        workers_.reserve(workers_.size() + 1);
        auto finished = std::make_shared<std::atomic_bool>(false);
        std::jthread thread([this, mesh, result, request = std::move(request), frame, settings, fi,
                             ei, id, revision, current, finished](std::stop_token stop) {
            struct CompletionFlag {
                std::atomic_bool *finished;
                ~CompletionFlag() {
                    finished->store(true);
                }
            } markFinished{finished.get()};
            try {
                const auto &frequency = result->frequencies.at(fi);
                const auto &solution = frequency.excitations.at(ei);
                const auto lateral =
                    settings.vertical ? verticalDirection(frame) : frame.horizontal;
                std::vector<Vec3> points;
                points.reserve(std::size_t(settings.lateralCount) * settings.axialCount);
                for (unsigned z = 0; z < settings.axialCount; ++z)
                    for (unsigned x = 0; x < settings.lateralCount; ++x) {
                        const double ax = settings.minimumAxial +
                                          (settings.maximumAxial - settings.minimumAxial) * z /
                                              (settings.axialCount - 1);
                        const double lat = settings.minimumLateral +
                                           (settings.maximumLateral - settings.minimumLateral) * x /
                                               (settings.lateralCount - 1);
                        points.push_back(
                            add(frame.origin, add(mul(frame.forward, ax), mul(lateral, lat))));
                    }
                const auto samples = evaluatePressureField(
                    *mesh, request, frequency.frequency, solution, points,
                    {settings.exclusionDistance, request.numerics.nearRelativeTolerance,
                     request.numerics.nearMaximumDepth},
                    stop, [this, id](std::size_t done, std::size_t total) {
                        if (done % 32 == 0)
                            QMetaObject::invokeMethod(
                                this,
                                [this, id, done, total] {
                                    if (id == jobId_)
                                        status_->setText(
                                            tr("Schallfeld %1 / %2 Punkte").arg(done).arg(total));
                                },
                                Qt::QueuedConnection);
                    });
                auto plane = std::make_shared<PressureFieldPlane>();
                plane->columns = settings.lateralCount;
                plane->rows = settings.axialCount;
                plane->frame = frame;
                plane->vertical = settings.vertical;
                plane->minimumLateral = settings.minimumLateral;
                plane->maximumLateral = settings.maximumLateral;
                plane->minimumAxial = settings.minimumAxial;
                plane->maximumAxial = settings.maximumAxial;
                plane->frequency = frequency.frequency;
                plane->revision = revision;
                plane->current = current;
                plane->excitation = QString::fromStdString(solution.name);
                std::size_t unresolved = 0, air = 0;
                for (const auto &sample : samples) {
                    plane->positions.push_back(sample.position);
                    plane->pressure.push_back(sample.pressure);
                    const bool valid = sample.status == FieldPointStatus::air;
                    plane->valid.push_back(valid);
                    air += valid ? 1 : 0;
                    unresolved += sample.status == FieldPointStatus::quadratureNotConverged ? 1 : 0;
                }
                validatePressureFieldPlane(*plane);
                const double largestStep = std::max(
                    (settings.maximumLateral - settings.minimumLateral) /
                        (settings.lateralCount - 1),
                    (settings.maximumAxial - settings.minimumAxial) / (settings.axialCount - 1));
                const bool coarse =
                    largestStep > request.medium.soundSpeed / frequency.frequency / 6;
                QMetaObject::invokeMethod(
                    this,
                    [this, id, key = cacheKey(settings), plane = std::move(plane), air, unresolved,
                     coarse]() mutable {
                        if (id != jobId_)
                            return;
                        busy_ = false;
                        error_.clear();
                        field_ = std::move(plane);
                        publishField();
                        status_->setText(
                            tr("%1 Luftpunkte · %2 Punkte ohne Quadraturkonvergenz%3")
                                .arg(air)
                                .arg(unresolved)
                                .arg(coarse ? tr(" · Raster gröber als λ/6") : QString{}));
                        try {
                            fieldCache_.push_front({key, field_, status_->text()});
                            while (fieldCache_.size() > 8)
                                fieldCache_.pop_back();
                        } catch (const std::bad_alloc &) {
                            status_->setText(status_->text() +
                                             tr(" · Feldspeicher nicht verfügbar"));
                        }
                        updateAvailability();
                    },
                    Qt::QueuedConnection);
            } catch (const std::exception &e) {
                const QString message = QString::fromUtf8(e.what());
                QMetaObject::invokeMethod(
                    this,
                    [this, id, message] {
                        if (id != jobId_)
                            return;
                        busy_ = false;
                        error_ = message;
                        updateAvailability();
                        status_->setText(message);
                    },
                    Qt::QueuedConnection);
            }
        });
        workers_.push_back({std::move(thread), std::move(finished)});
        busy_ = true;
        updateAvailability();
        status_->setText(tr("Schallfeld wird berechnet …"));
    } catch (const std::exception &e) {
        ++jobId_;
        for (auto &worker : workers_)
            worker.thread.request_stop();
        busy_ = false;
        error_ =
            tr("Feldberechnung konnte nicht gestartet werden: %1").arg(QString::fromUtf8(e.what()));
        updateAvailability();
        status_->setText(error_);
    }
}
void AcousticFieldWidget::exportSamples() {
    if (!field_)
        return;
    const auto &grid = *field_;
    const auto path =
        QFileDialog::getSaveFileName(this, tr("Schallfeld exportieren"), {}, tr("CSV (*.csv)"));
    if (path.isEmpty())
        return;
    QByteArray csvBytes(
        "# Complex RMS pressure; exp(+i omega t); coordinates in result observation frame\n");
    csvBytes +=
        "# frequency_Hz=" + QByteArray::number(frequency_->currentData().toDouble(), 'g', 17) +
        '\n';
    csvBytes += "# result_revision=" + QByteArray::number(revision_) + '\n';
    csvBytes += current_ ? "# result_current=1\n" : "# result_current=0\n";
    const auto vectorText = [](Vec3 v) {
        return QByteArray::number(v.x, 'g', 17) + ',' + QByteArray::number(v.y, 'g', 17) + ',' +
               QByteArray::number(v.z, 'g', 17);
    };
    csvBytes += "# origin_world_m=" + vectorText(frame_.origin) + '\n';
    csvBytes += "# axial_direction_world=" + vectorText(frame_.forward) + '\n';
    csvBytes += "# lateral_direction_world=" +
                vectorText(settings_.vertical ? verticalDirection(frame_) : frame_.horizontal) +
                '\n';
    csvBytes +=
        "lateral_m,axial_m,x_world_m,y_world_m,z_world_m,valid,pressure_real_Pa,pressure_imag_Pa\n";
    for (unsigned z = 0; z < grid.rows; ++z)
        for (unsigned x = 0; x < grid.columns; ++x) {
            const auto i = std::size_t(z) * grid.columns + x;
            csvBytes += QByteArray::number(grid.minimumLateral +
                                               (grid.maximumLateral - grid.minimumLateral) * x /
                                                   (grid.columns - 1),
                                           'g', 17) +
                        ',';
            csvBytes +=
                QByteArray::number(grid.minimumAxial + (grid.maximumAxial - grid.minimumAxial) * z /
                                                           (grid.rows - 1),
                                   'g', 17) +
                ',';
            csvBytes += vectorText(grid.positions[i]) + ',';
            csvBytes += grid.valid[i] ? "1," : "0,";
            if (grid.valid[i])
                csvBytes += QByteArray::number(grid.pressure[i].real(), 'g', 17) + ',' +
                            QByteArray::number(grid.pressure[i].imag(), 'g', 17);
            else
                csvBytes += ',';
            csvBytes += '\n';
        }
    QSaveFile file(path);
    if (!file.open(QIODevice::WriteOnly) || file.write(csvBytes) != csvBytes.size() ||
        !file.commit())
        QMessageBox::warning(this, tr("Export fehlgeschlagen"), file.errorString());
}
} // namespace hornsculpt
