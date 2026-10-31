#include "mesh_panel.hpp"
#include "numeric_input.hpp"

#include <QBrush>
#include <QCheckBox>
#include <QColor>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QItemSelectionModel>
#include <QLabel>
#include <QLocale>
#include <QPushButton>
#include <QSignalBlocker>
#include <QTreeWidget>
#include <QVBoxLayout>
#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <stdexcept>

namespace hornsculpt {
namespace {
const std::array<MeshSettings, 3> presets{MeshSettings{.012, {true, .00025, .00035, .6}, {}},
                                          MeshSettings{.008, {true, .0001, .00025, .3}, {}},
                                          MeshSettings{.004, {true, .00005, .0001, .25}, {}}};

bool equalSettings(const MeshSettings &a, const MeshSettings &b) {
    const auto near = [](double x, double y) { return std::abs(x - y) <= 1e-12; };
    return a.adaptive.enabled == b.adaptive.enabled &&
           near(a.targetEdgeLength, b.targetEdgeLength) &&
           near(a.adaptive.minimumEdgeLength, b.adaptive.minimumEdgeLength) &&
           near(a.adaptive.geometricTolerance, b.adaptive.geometricTolerance) &&
           near(a.adaptive.growthRate, b.adaptive.growthRate);
}
} // namespace

MeshPanel::MeshPanel(QWidget *parent) : QWidget(parent) {
    auto *layout = new QVBoxLayout(this);
    layout->setContentsMargins(8, 8, 8, 8);
    auto *form = new QFormLayout;
    form->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
    preset_ = new QComboBox;
    preset_->setObjectName(QStringLiteral("meshPreset"));
    preset_->addItems({tr("Vorschau"), tr("Entwurf"), tr("Referenz"), tr("Benutzerdefiniert")});
    preset_->setToolTip(tr("Setzt das Größenfeld. Die eigenen "
                           "Qualitätsgrenzen bleiben unverändert."));
    form->addRow(tr("Vorgabe"), preset_);
    const auto number = [&](const QString &name, const QString &label, double minimum,
                            double maximum, const QString &suffix, const QString &tooltip) {
        auto *input = new DecimalInput;
        input->setObjectName(name);
        input->setDecimals(3);
        input->setRange(minimum, maximum);
        input->setSuffix(suffix);
        input->setSingleStep(.05);
        input->setToolTip(tooltip);
        input->setWhatsThis(tooltip);
        form->addRow(label, input);
        return input;
    };
    targetEdge_ = number(QStringLiteral("meshTargetEdgeLength"), tr("Zielkante"), .05, 100, " mm",
                         tr("Maximale Zielkante auf großen Flächen. Dasselbe isotrope "
                            "Größenfeld gilt für Draft und FreeCAD."));
    adaptive_ = new QCheckBox(tr("Adaptiv"));
    adaptive_->setObjectName(QStringLiteral("meshAdaptiveSizing"));
    adaptive_->setToolTip(tr("Lokales skalares Größenfeld an kleinen Flächen und "
                             "gekrümmten Details; Elemente bleiben isotrop."));
    form->addRow(adaptive_);
    minimumEdge_ =
        number(QStringLiteral("meshMinimumEdgeLength"), tr("Min. Kante"), .001, 100, " mm",
               tr("Untere Zielkante der lokalen Verfeinerung; darf die "
                  "Zielkante nicht überschreiten. Keine harte Mindestlänge jeder Netzkante: "
                  "Geometrische Detailabweichung, feste Ränder und gültige Topologie können "
                  "kürzere Kanten erfordern."));
    tolerance_ =
        number(QStringLiteral("meshGeometricTolerance"), tr("Detailabweichung"), .001, 100, " mm",
               tr("Genauigkeit der gesamten Draft-Tessellation und Zielwert für "
                  "krümmungsabhängige Netzverfeinerung. STEP/BRep-Netzknoten werden auf "
                  "exakte CAD-Flächen und gemeinsame Randkurven projiziert. Ältere reine "
                  "Polygonquellen bleiben durch ihre Tessellation begrenzt."));
    growth_ = number(QStringLiteral("meshSizeGrowth"), tr("Wachstum"), .01, 1, "",
                     tr("Maximales Wachstum dh/dd der lokalen Zielkante pro Abstand. "
                        "Dimensionslos; kleine Werte ergeben sanftere Übergänge."));
    qualityLimits_ = new QCheckBox(tr("Qualitätsgrenzen"));
    qualityLimits_->setObjectName(QStringLiteral("meshQualityLimitsEnabled"));
    qualityLimits_->setToolTip(
        tr("Optionale Winkel- und Seitenverhältnisgrenzen. Aus: keine Verfeinerung oder "
           "Ablehnung allein wegen dieser Werte. Formschluss und gültige Elemente bleiben "
           "erforderlich."));
    qualityLimits_->setWhatsThis(qualityLimits_->toolTip());
    form->addRow(qualityLimits_);
    minimumAngle_ =
        number(QStringLiteral("meshMinimumAngle"), tr("Min. Winkel"), .001, 59.999, " °",
               tr("Kleinster erlaubter Innenwinkel eines Dreiecks in Grad. "
                  "Nur bei aktiven Qualitätsgrenzen wirksam."));
    minimumAngle_->setSingleStep(1);
    maximumAspect_ =
        number(QStringLiteral("meshMaximumAspect"), tr("Max. Seitenverh."), 1.155, 1000000, "",
               tr("Längste Kante geteilt durch kleinste Höhe. "
                  "Gleichseitig: 1,155. Nur bei aktiven Qualitätsgrenzen wirksam."));
    maximumAspect_->setSingleStep(.1);
    layout->addLayout(form);
    faces_ = new QTreeWidget;
    faces_->setObjectName(QStringLiteral("meshFaces"));
    faces_->setHeaderLabels({tr("Fläche"), tr("h [mm]"), tr("n")});
    faces_->setRootIsDecorated(false);
    faces_->setSelectionMode(QAbstractItemView::ExtendedSelection);
    faces_->setMinimumHeight(100);
    faces_->setToolTip(tr("Lokale Zielkante für ausgewählte Flächen. Strg erweitert die Auswahl; "
                          "ohne Zuweisung gilt das globale Größenfeld."));
    layout->addWidget(faces_, 1);
    auto *localForm = new QFormLayout;
    localTarget_ = new DecimalInput;
    localTarget_->setObjectName(QStringLiteral("meshLocalTarget"));
    localTarget_->setDecimals(3);
    localTarget_->setRange(.001, 100);
    localTarget_->setValue(2);
    localTarget_->setSuffix(QStringLiteral(" mm"));
    localTarget_->setSingleStep(.1);
    localTarget_->setToolTip(
        tr("Lokale isotrope Zielkante der Auswahl. Nur Verfeinerung zwischen "
           "globaler Min. Kante und Zielkante. Übergänge folgen dem eingestellten Wachstum."));
    localForm->addRow(tr("Lokale Zielkante"), localTarget_);
    layout->addLayout(localForm);
    auto *localButtons = new QHBoxLayout;
    assignLocal_ = new QPushButton(tr("Zuweisen"));
    resetLocal_ = new QPushButton(tr("Zurücksetzen"));
    assignLocal_->setObjectName(QStringLiteral("meshAssignLocalTarget"));
    resetLocal_->setObjectName(QStringLiteral("meshResetLocalTarget"));
    resetLocal_->setToolTip(
        tr("Lokale Zuweisung der Auswahl entfernen; danach gilt das globale Größenfeld."));
    localButtons->addWidget(assignLocal_);
    localButtons->addWidget(resetLocal_);
    layout->addLayout(localButtons);
    localStatus_ = new QLabel;
    layout->addWidget(localStatus_);
    reverseFaces_ = new QPushButton(tr("Normalen umkehren"));
    reverseFaces_->setObjectName(QStringLiteral("meshReverseFaces"));
    reverseFaces_->setToolTip(
        tr("Kehrt die Normalen aller ausgewählten Flächen um. Koordinaten und komplexe "
           "Normalgeschwindigkeiten bleiben unverändert. BEM erfordert anschließend ein "
           "konsistent nach außen orientiertes, geschlossenes Netz."));
    reverseFaces_->setWhatsThis(reverseFaces_->toolTip());
    layout->addWidget(reverseFaces_);
    orientationStatus_ = new QLabel;
    orientationStatus_->setObjectName(QStringLiteral("meshOrientationStatus"));
    orientationStatus_->setWordWrap(true);
    layout->addWidget(orientationStatus_);
    connect(reverseFaces_, &QPushButton::clicked, this, [this] {
        const auto ids = selectedFaces();
        if (!ids.empty() && onFaceOrientationChanged)
            onFaceOrientationChanged(ids);
    });
    connect(faces_, &QTreeWidget::itemSelectionChanged, this, [this] {
        selectionChanged();
        if (onFacesSelected)
            onFacesSelected(selectedFaces());
    });
    connect(assignLocal_, &QPushButton::clicked, this, [this] { assignLocalTarget(); });
    connect(resetLocal_, &QPushButton::clicked, this, [this] {
        const auto ids = selectedFaces();
        if (!ids.empty() && onFaceTargetChanged)
            onFaceTargetChanged(ids, {});
    });
    for (auto *input :
         {targetEdge_, minimumEdge_, tolerance_, growth_, minimumAngle_, maximumAspect_})
        connect(input, &QDoubleSpinBox::valueChanged, this, [this] { controlsChanged(); });
    connect(adaptive_, &QCheckBox::toggled, this, [this] { controlsChanged(); });
    connect(qualityLimits_, &QCheckBox::toggled, this, [this] { controlsChanged(); });
    connect(preset_, &QComboBox::currentIndexChanged, this, [this](int index) {
        if (index < 0 || index >= static_cast<int>(presets.size()))
            return;
        auto settings = presets[static_cast<std::size_t>(index)];
        settings.quality = inputSettings().quality;
        setSettings(settings);
        if (onChanged)
            onChanged();
    });
    setSettings(MeshSettings{});
    setOrientationStatus({});
    selectionChanged();
}
MeshSettings MeshPanel::inputSettings() const {
    return {targetEdge_->value() * .001,
            {adaptive_->isChecked(), minimumEdge_->value() * .001, tolerance_->value() * .001,
             growth_->value()},
            {qualityLimits_->isChecked(), minimumAngle_->value(), maximumAspect_->value()}};
}
MeshPanel::~MeshPanel() {
    // Selection teardown must not access sizing maps after their destruction.
    for (auto *child : findChildren<QObject *>())
        disconnect(child, nullptr, this, nullptr);
}

MeshSettings MeshPanel::settings() const {
    const auto result = inputSettings();
    validateMeshSettings(result);
    return result;
}
void MeshPanel::setSettings(const MeshSettings &settings) {
    validateMeshSettings(settings);
    setInputSettings(settings);
}
void MeshPanel::setInputSettings(const MeshSettings &settings) {
    const auto check = [](const QDoubleSpinBox *input, double value) {
        const double scale = std::pow(10.0, input->decimals());
        if (!std::isfinite(value) || value < input->minimum() || value > input->maximum() ||
            std::abs(std::round(value * scale) / scale - value) >
                1e-9 * std::max(1.0, std::abs(value)))
            throw std::invalid_argument(
                "Gespeicherte Netzeinstellung liegt außerhalb des Eingabebereichs.");
    };
    check(targetEdge_, settings.targetEdgeLength * 1000);
    check(minimumEdge_, settings.adaptive.minimumEdgeLength * 1000);
    check(tolerance_, settings.adaptive.geometricTolerance * 1000);
    check(growth_, settings.adaptive.growthRate);
    check(minimumAngle_, settings.quality.minimumAngleDegrees);
    check(maximumAspect_, settings.quality.maximumAspectRatio);
    const QSignalBlocker a(targetEdge_), b(adaptive_), c(minimumEdge_), d(tolerance_), e(growth_),
        f(qualityLimits_), g(minimumAngle_), h(maximumAspect_);
    targetEdge_->setValue(settings.targetEdgeLength * 1000);
    adaptive_->setChecked(settings.adaptive.enabled);
    minimumEdge_->setValue(settings.adaptive.minimumEdgeLength * 1000);
    tolerance_->setValue(settings.adaptive.geometricTolerance * 1000);
    growth_->setValue(settings.adaptive.growthRate);
    qualityLimits_->setChecked(settings.quality.enabled);
    minimumAngle_->setValue(settings.quality.minimumAngleDegrees);
    maximumAspect_->setValue(settings.quality.maximumAspectRatio);
    updateControls();
    updatePreset();
    updateLocalStatus();
}
void MeshPanel::updateControls() {
    for (auto *input : {minimumEdge_, growth_})
        input->setEnabled(adaptive_->isChecked());
    minimumAngle_->setEnabled(qualityLimits_->isChecked());
    maximumAspect_->setEnabled(qualityLimits_->isChecked());
}
void MeshPanel::updatePreset() {
    const auto settings = inputSettings();
    const auto found = std::find_if(presets.begin(), presets.end(), [&](const auto &preset) {
        return equalSettings(settings, preset);
    });
    const QSignalBlocker blocker(preset_);
    preset_->setCurrentIndex(static_cast<int>(std::distance(presets.begin(), found)));
}
void MeshPanel::controlsChanged() {
    updateControls();
    updatePreset();
    updateLocalStatus();
    if (onChanged)
        onChanged();
}

void MeshPanel::setFaces(std::vector<MeshFaceItem> faces, std::map<std::uint32_t, double> targets) {
    if (faces == faceItems_ && targets == faceTargets_)
        return;
    const auto selection = selectedFaces();
    faceItems_ = std::move(faces);
    faceTargets_ = std::move(targets);
    const QSignalBlocker blocker(faces_);
    faces_->clear();
    for (const auto &face : faceItems_) {
        const auto target = faceTargets_.find(face.id);
        auto *item = new QTreeWidgetItem(
            faces_, {face.label, target == faceTargets_.end()
                                     ? tr("global")
                                     : QLocale().toString(target->second * 1000, 'f', 3)});
        item->setData(0, Qt::UserRole, face.id);
        item->setToolTip(0, tr("Flächen-ID %1").arg(face.id));
        item->setSelected(std::binary_search(selection.begin(), selection.end(), face.id));
    }
    faces_->resizeColumnToContents(1);
    setReversedFaces(reversedFaces_);
    selectionChanged();
}
void MeshPanel::setReversedFaces(const std::set<std::uint32_t> &faces) {
    reversedFaces_ = faces;
    const QSignalBlocker blocker(faces_);
    for (int row = 0; row < faces_->topLevelItemCount(); ++row) {
        auto *item = faces_->topLevelItem(row);
        const bool reversed = reversedFaces_.contains(item->data(0, Qt::UserRole).toUInt());
        item->setText(2, reversed ? QStringLiteral("−") : QStringLiteral("+"));
        item->setToolTip(2, reversed ? tr("Gegenüber der Quellgeometrie umgekehrt")
                                     : tr("Orientierung der Quellgeometrie"));
    }
    faces_->resizeColumnToContents(2);
}
void MeshPanel::setOrientationStatus(QString status, bool error) {
    orientationStatus_->setText(std::move(status));
    orientationStatus_->setStyleSheet(error ? QStringLiteral("color: #b1312b") : QString{});
    orientationStatus_->setVisible(!orientationStatus_->text().isEmpty());
}
std::vector<std::uint32_t> MeshPanel::selectedFaces() const {
    std::vector<std::uint32_t> result;
    for (const auto *item : faces_->selectedItems())
        result.push_back(item->data(0, Qt::UserRole).toUInt());
    std::sort(result.begin(), result.end());
    return result;
}
void MeshPanel::setSelectedFaces(const std::vector<std::uint32_t> &faces) {
    if (selectedFaces() == faces)
        return;
    const QSignalBlocker blocker(faces_);
    faces_->clearSelection();
    for (int row = 0; row < faces_->topLevelItemCount(); ++row) {
        auto *item = faces_->topLevelItem(row);
        item->setSelected(std::find(faces.begin(), faces.end(),
                                    item->data(0, Qt::UserRole).toUInt()) != faces.end());
    }
    selectionChanged();
}
void MeshPanel::selectFace(std::uint32_t face, Qt::KeyboardModifiers modifiers) {
    for (int row = 0; row < faces_->topLevelItemCount(); ++row) {
        auto *item = faces_->topLevelItem(row);
        if (item->data(0, Qt::UserRole).toUInt() != face)
            continue;
        const QSignalBlocker blocker(faces_);
        const bool selected = item->isSelected();
        if (!(modifiers & Qt::ControlModifier))
            faces_->clearSelection();
        faces_->setCurrentItem(item, 0, QItemSelectionModel::NoUpdate);
        item->setSelected((modifiers & Qt::ControlModifier) ? !selected : true);
        faces_->scrollToItem(item);
        selectionChanged();
        if (onFacesSelected)
            onFacesSelected(selectedFaces());
        return;
    }
}
void MeshPanel::selectionChanged() {
    const auto selection = selectedFaces();
    localTarget_->setEnabled(!selection.empty());
    assignLocal_->setEnabled(!selection.empty());
    reverseFaces_->setEnabled(!selection.empty());
    resetLocal_->setEnabled(std::any_of(selection.begin(), selection.end(),
                                        [&](auto id) { return faceTargets_.contains(id); }));
    if (!selection.empty() && !localTarget_->hasFocus()) {
        const auto first = faceTargets_.find(selection.front());
        localTarget_->setValue(first == faceTargets_.end() ? targetEdge_->value()
                                                           : first->second * 1000);
    }
    updateLocalStatus();
}
void MeshPanel::updateLocalStatus() {
    const auto settings = inputSettings();
    std::size_t invalid = 0;
    for (int row = 0; row < faces_->topLevelItemCount(); ++row) {
        auto *item = faces_->topLevelItem(row);
        const auto target = faceTargets_.find(item->data(0, Qt::UserRole).toUInt());
        const bool outside =
            target != faceTargets_.end() && (!std::isfinite(target->second) ||
                                             target->second < settings.adaptive.minimumEdgeLength ||
                                             target->second > settings.targetEdgeLength);
        item->setForeground(1, outside ? QBrush(QColor(177, 49, 43)) : QBrush());
        if (outside)
            ++invalid;
    }
    localStatus_->setText(invalid ? tr("%1 lokale Zielkanten außerhalb der Grenzen").arg(invalid)
                          : !faceNotice_.isEmpty()
                              ? faceNotice_
                              : tr("%1 Flächen ausgewählt").arg(selectedFaces().size()));
    localStatus_->setStyleSheet(invalid ? QStringLiteral("color: #b1312b") : QString{});
    localStatus_->setToolTip(tr("Zulässige lokale Zielkante: %1 bis %2 mm. Zuweisungen werden "
                                "nicht automatisch geändert.")
                                 .arg(settings.adaptive.minimumEdgeLength * 1000)
                                 .arg(settings.targetEdgeLength * 1000));
}
void MeshPanel::setFaceNotice(QString notice) {
    faceNotice_ = std::move(notice);
    updateLocalStatus();
}
void MeshPanel::assignLocalTarget() {
    const auto ids = selectedFaces();
    if (ids.empty())
        return;
    try {
        const auto globals = settings();
        const double target = localTarget_->value() * .001;
        if (target < globals.adaptive.minimumEdgeLength || target > globals.targetEdgeLength)
            throw std::invalid_argument(
                "Lokale Zielkante muss zwischen Min. Kante und globaler Zielkante liegen.");
        if (onFaceTargetChanged)
            onFaceTargetChanged(ids, target);
    } catch (const std::invalid_argument &error) {
        localStatus_->setText(QString::fromUtf8(error.what()));
        localStatus_->setStyleSheet(QStringLiteral("color: #b1312b"));
    }
}

} // namespace hornsculpt
