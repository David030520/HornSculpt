#pragma once

#include "mesh_settings.hpp"
#include <QWidget>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <vector>

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QPushButton;
class QTreeWidget;

namespace hornsculpt {

struct MeshFaceItem {
    std::uint32_t id{};
    QString label;
    bool operator==(const MeshFaceItem &) const = default;
};

class MeshPanel final : public QWidget {
  public:
    explicit MeshPanel(QWidget *parent = nullptr);
    ~MeshPanel() override;
    MeshSettings settings() const;
    void setSettings(const MeshSettings &settings);
    MeshSettings inputSettings() const;
    void setInputSettings(const MeshSettings &settings);
    void setFaces(std::vector<MeshFaceItem> faces, std::map<std::uint32_t, double> targets);
    void setSelectedFaces(const std::vector<std::uint32_t> &faces);
    void selectFace(std::uint32_t face, Qt::KeyboardModifiers modifiers = Qt::NoModifier);
    std::vector<std::uint32_t> selectedFaces() const;
    void setFaceNotice(QString notice);
    void setReversedFaces(const std::set<std::uint32_t> &faces);
    void setOrientationStatus(QString status, bool error = false);
    std::function<void()> onChanged;
    std::function<void(std::vector<std::uint32_t>)> onFacesSelected;
    std::function<void(std::vector<std::uint32_t>, std::optional<double>)> onFaceTargetChanged;
    std::function<void(std::vector<std::uint32_t>)> onFaceOrientationChanged;

  private:
    void controlsChanged();
    void updateControls();
    void updatePreset();
    void selectionChanged();
    void updateLocalStatus();
    void assignLocalTarget();
    QComboBox *preset_{};
    QDoubleSpinBox *targetEdge_{};
    QCheckBox *adaptive_{};
    QDoubleSpinBox *minimumEdge_{};
    QDoubleSpinBox *tolerance_{};
    QDoubleSpinBox *growth_{};
    QCheckBox *qualityLimits_{};
    QDoubleSpinBox *minimumAngle_{};
    QDoubleSpinBox *maximumAspect_{};
    QTreeWidget *faces_{};
    QDoubleSpinBox *localTarget_{};
    QPushButton *assignLocal_{};
    QPushButton *resetLocal_{};
    QLabel *localStatus_{};
    QPushButton *reverseFaces_{};
    QLabel *orientationStatus_{};
    std::vector<MeshFaceItem> faceItems_;
    std::map<std::uint32_t, double> faceTargets_;
    std::set<std::uint32_t> reversedFaces_;
    QString faceNotice_;
};

} // namespace hornsculpt
