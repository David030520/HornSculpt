#pragma once

#include <QDialog>
#include <array>
#include <hornsculpt/bem.hpp>
#include <optional>

class QComboBox;
class QDialogButtonBox;
class QLabel;
class QLineEdit;

namespace hornsculpt {
class DecimalInput;

class BaffleDialog final : public QDialog {
  public:
    explicit BaffleDialog(std::optional<InfiniteBaffle> baffle = {}, QWidget *parent = nullptr);
    [[nodiscard]] std::optional<InfiniteBaffle> infiniteBaffle() const;
    void accept() override;

  private:
    void validateInput();
    InfiniteBaffle inputBaffle() const;
    QComboBox *mode_{};
    std::array<DecimalInput *, 3> origin_{};
    std::array<DecimalInput *, 3> normal_{};
    QLineEdit *aperture_{};
    QLabel *error_{};
    QDialogButtonBox *buttons_{};
    std::optional<InfiniteBaffle> committed_;
    SurfaceMesh exteriorDrivenMesh_;
};
} // namespace hornsculpt
