#pragma once

#include <QDockWidget>
#include <functional>
#include <hornsculpt/freecad_link.hpp>
#include <map>
#include <optional>
#include <set>
#include <vector>

class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QMainWindow;
class QPushButton;
class QTreeWidget;

namespace hornsculpt {

class FreeCadDock final : public QDockWidget {
  public:
    explicit FreeCadDock(QMainWindow *parent);
    ~FreeCadDock() override;
    void ensureListening();
    bool isSourceAvailable() const;
    bool hasLiveLink() const;
    bool isGeometryPending() const;
    bool setAssignments(const QString &documentId, const QString &sessionId, std::uint64_t revision,
                        const std::vector<CadFaceAssignment> &faces,
                        const std::vector<std::uint32_t> &missing);
    std::shared_ptr<const CadSnapshot> snapshot() const;
    void restoreSnapshot(std::shared_ptr<const CadSnapshot> snapshot);
    void restoreSnapshot(const PreparedCadSnapshot &snapshot);
    void resetProject();
    void selectFace(std::uint32_t faceId, Qt::KeyboardModifiers modifiers = Qt::NoModifier);
    void clearSelection();
    void setSelectedFaces(const std::vector<std::uint32_t> &faces);
    std::vector<std::uint32_t> selectedFaces() const;
    void setFaceTargets(const std::map<std::uint32_t, double> &targets);
    void setReversedFaces(const std::set<std::uint32_t> &faces);
    void setInterfaceFaces(std::set<std::uint32_t> faces, bool active);
    void setBemInterfaceFaces(std::set<std::uint32_t> faces);
    std::function<void(std::shared_ptr<const CadSnapshot>)> onSnapshotChanged;
    std::function<void(bool)> onAvailabilityChanged;
    std::function<void(std::vector<std::uint32_t>)> onFacesSelected;
    std::function<void(QString)> onError;
    std::function<void()> onAssignmentEditStarted;
    std::function<void()> onEditBemRegions;
    std::function<void(std::vector<std::uint32_t>, std::optional<double>)> onFaceTargetChanged;

  private:
    void refresh();
    void selectionChanged();
    void publishAvailability();
    void installWorkbench();
    FreeCadLink link_;
    QLabel *status_{};
    QLabel *document_{};
    QLabel *mapping_{};
    QLabel *selectionLabel_{};
    QPushButton *connect_{};
    QPushButton *install_{};
    QPushButton *apply_{};
    QPushButton *bemRegions_{};
    QPushButton *confirm_{};
    QPushButton *missing_{};
    QTreeWidget *faces_{};
    QComboBox *boundary_{};
    QDoubleSpinBox *velocityReal_{};
    QDoubleSpinBox *velocityImaginary_{};
    bool sourceAvailable_{};
    bool offlineSnapshot_{};
    bool awaitingSnapshot_{};
    QString selectionDocumentId_;
    std::map<std::uint32_t, double> faceTargets_;
    std::set<std::uint32_t> reversedFaces_;
    std::set<std::uint32_t> interfaceFaces_;
    std::set<std::uint32_t> bemInterfaceFaces_;
    bool interfacesActive_{};
};

} // namespace hornsculpt
