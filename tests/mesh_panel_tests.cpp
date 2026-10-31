#include "mesh_panel.hpp"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QTreeWidget>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace hornsculpt;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
template <class Control> Control *control(MeshPanel &panel, const char *name) {
    auto *result = panel.findChild<Control *>(QString::fromLatin1(name));
    require(result != nullptr, "Missing mesh control.");
    return result;
}
} // namespace
int main(int argc, char **argv) {
    QApplication app(argc, argv);
    try {
        int changes = 0;
        MeshPanel panel;
        require(!panel.findChild<QDoubleSpinBox *>(QStringLiteral("meshEntryAngleTolerance")),
                "The mesh panel still exposes obsolete entry-angle meshing controls.");
        panel.onChanged = [&] { ++changes; };
        auto *enabled = control<QCheckBox>(panel, "meshQualityLimitsEnabled");
        auto *angle = control<QDoubleSpinBox>(panel, "meshMinimumAngle");
        auto *aspect = control<QDoubleSpinBox>(panel, "meshMaximumAspect");
        auto settings = panel.settings();
        require(settings.quality.enabled && angle->value() == 20 && aspect->value() == 4,
                "Default quality settings changed.");
        settings.quality = {false, 12.5, 9};
        panel.setSettings(settings);
        require(changes == 0 && !enabled->isChecked() && !angle->isEnabled() &&
                    !aspect->isEnabled(),
                "Loading disabled limits changed the project or left active controls.");
        enabled->setChecked(true);
        require(changes == 1 && angle->isEnabled() && aspect->isEnabled(),
                "The quality toggle is not connected to the shared mesh settings.");
        angle->setValue(10);
        aspect->setValue(12);
        settings = panel.settings();
        require(changes == 3 && settings.quality.enabled &&
                    settings.quality.minimumAngleDegrees == 10 &&
                    settings.quality.maximumAspectRatio == 12,
                "Numeric quality edits do not reach the mesh settings.");
        control<QComboBox>(panel, "meshPreset")->setCurrentIndex(0);
        settings = panel.settings();
        require(settings.targetEdgeLength == .012 && settings.quality.minimumAngleDegrees == 10 &&
                    settings.quality.maximumAspectRatio == 12,
                "Changing the size preset overwrites the user's quality limits.");
        enabled->setChecked(false);
        require(!panel.settings().quality.enabled && !angle->isEnabled() && !aspect->isEnabled(),
                "Quality limits cannot be disabled after changing presets.");
        control<QCheckBox>(panel, "meshAdaptiveSizing")->setChecked(false);
        auto *tolerance = control<QDoubleSpinBox>(panel, "meshGeometricTolerance");
        require(tolerance->isEnabled(),
                "Global geometry accuracy became unavailable with adaptive sizing disabled.");
        tolerance->setValue(.125);
        require(std::abs(panel.settings().adaptive.geometricTolerance - .000125) < 1e-12,
                "Global geometry accuracy does not reach the mesh settings.");
        panel.setFaces({{1, "Innenwand"}, {5, "Driven"}}, {});
        panel.setSelectedFaces({1, 5});
        std::vector<std::uint32_t> reversed;
        panel.onFaceOrientationChanged = [&](auto ids) { reversed = std::move(ids); };
        auto *reverse = control<QPushButton>(panel, "meshReverseFaces");
        require(reverse->isEnabled(), "Selected faces cannot reverse their normals.");
        reverse->click();
        require(reversed == std::vector<std::uint32_t>{1, 5},
                "Normal reversal does not preserve the complete multi-selection.");
        panel.setReversedFaces({1, 5});
        auto *faces = control<QTreeWidget>(panel, "meshFaces");
        require(faces->topLevelItem(0)->text(2) == QString::fromUtf8("−") &&
                    faces->topLevelItem(1)->text(2) == QString::fromUtf8("−"),
                "Reversed face orientation is not visible in the face table.");
        panel.setSelectedFaces({});
        require(!reverse->isEnabled(), "Normal reversal remains active without selected faces.");
        std::cout << "Shared mesh quality controls and size-preset independence passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
