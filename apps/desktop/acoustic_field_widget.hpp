#pragma once

#include "field_display_settings.hpp"
#include <QWidget>
#include <hornsculpt/acoustic_frame.hpp>
#include <hornsculpt/field_evaluation.hpp>
#include <hornsculpt/pressure_field_plane.hpp>
#include <atomic>
#include <list>
#include <memory>
#include <optional>
#include <thread>

class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QPushButton;
class QSpinBox;
class QSlider;
class QTimer;
class QToolButton;
namespace hornsculpt {
class AcousticFieldWidget final : public QWidget {
  public:
    explicit AcousticFieldWidget(QWidget *parent = nullptr);
    ~AcousticFieldWidget() override;
    FieldDisplaySettings settings() const noexcept {
        return settings_;
    }
    void setSettings(const FieldDisplaySettings &settings);
    void setSnapshot(std::shared_ptr<const SurfaceMesh> mesh, std::optional<BemRequest> request,
                     std::shared_ptr<const BemResult> result, AcousticFrame frame,
                     std::uint64_t revision, bool current);
    void cancel();
    // Explicit field postprocessing on the stored immutable result snapshot.
    // Unavailable input or worker failure is reported through lastError().
    void requestEvaluation();
    bool evaluationInProgress() const noexcept;
    bool fieldReady() const noexcept;
    QString lastError() const;
    std::shared_ptr<const PressureFieldPlane> fieldPlane() const noexcept;
    std::function<void()> onChanged;
    std::function<void(std::shared_ptr<const PressureFieldPlane>, PressureFieldAppearance)>
        onFieldChanged;
    std::function<void()> onFocusField;
    std::function<void(bool)> onProbeModeChanged;

  protected:
    void showEvent(QShowEvent *event) override;
    void hideEvent(QHideEvent *event) override;

  private:
    void controlsChanged(bool sampling);
    bool canEvaluate() const noexcept;
    void exportSamples();
    void updateAvailability();
    void clearField();
    void publishField();
    bool restoreCachedField();
    FieldDisplaySettings settings_;
    std::shared_ptr<const SurfaceMesh> mesh_;
    std::optional<BemRequest> request_;
    std::shared_ptr<const BemResult> result_;
    AcousticFrame frame_;
    std::uint64_t revision_{}, jobId_{};
    bool current_{}, busy_{};
    QString error_;
    struct FieldWorker {
        std::jthread thread;
        std::shared_ptr<std::atomic_bool> finished;
    };
    std::vector<FieldWorker> workers_;
    QComboBox *plane_{}, *quantity_{}, *frequency_{}, *excitation_{};
    QDoubleSpinBox *minimumLateral_{}, *maximumLateral_{}, *minimumAxial_{}, *maximumAxial_{};
    QDoubleSpinBox *minimumDb_{}, *maximumDb_{}, *phase_{}, *exclusion_{};
    QSpinBox *columns_{}, *rows_{};
    QPushButton *evaluate_{}, *cancel_{}, *export_{};
    QLabel *status_{};
    std::shared_ptr<const PressureFieldPlane> field_;
    struct CachedField {
        FieldDisplaySettings key;
        std::shared_ptr<const PressureFieldPlane> plane;
        QString status;
    };
    std::list<CachedField> fieldCache_;
    QSlider *frequencySlider_{};
    QToolButton *previousFrequency_{}, *nextFrequency_{};
    QToolButton *probe_{};
    QPushButton *animate_{};
    QTimer *animationTimer_{};
    double animationPhase_{};
    QPushButton *focus_{};
};
} // namespace hornsculpt
