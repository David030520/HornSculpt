#pragma once

#include <QDialog>
#include <array>
#include <hornsculpt/acoustic_frame.hpp>
#include <optional>

class QComboBox;
class QDialogButtonBox;
class QLabel;

namespace hornsculpt {

class DecimalInput;

class AcousticFrameDialog final : public QDialog {
  public:
    explicit AcousticFrameDialog(const AcousticFrame &sourceFrame,
                                 std::optional<AcousticFrame> overrideFrame = {},
                                 QWidget *parent = nullptr);

    // Read after Accepted. An empty override follows the current geometry.
    [[nodiscard]] std::optional<AcousticFrame> observationFrame() const;
    void accept() override;

  private:
    struct DirectionControls {
        QComboBox *preset{};
        std::array<DecimalInput *, 3> coordinates{};
    };

    void displayFrame(const AcousticFrame &frame);
    void modeChanged();
    void presetChanged(DirectionControls &controls);
    void coordinatesChanged(DirectionControls &controls);
    void updateReadOnly();
    void validateFrame();
    [[nodiscard]] AcousticFrame inputFrame() const;

    AcousticFrame sourceFrame_;
    AcousticFrame customFrame_;
    std::optional<AcousticFrame> committedFrame_;
    std::optional<AcousticFrame> validatedFrame_;
    QComboBox *mode_{};
    std::array<DecimalInput *, 3> origin_{};
    DirectionControls forward_;
    DirectionControls horizontal_;
    std::array<QLabel *, 3> vertical_{};
    QLabel *error_{};
    QDialogButtonBox *buttons_{};
    bool showingCustom_{};
    bool loading_{};
};

} // namespace hornsculpt
