#pragma once

#include "field_display_settings.hpp"
#include <QWidget>
#include <complex>
#include <functional>
#include <vector>

class QCheckBox;
class QDoubleSpinBox;
class QLabel;
namespace hornsculpt {
class PlotWidget;
struct PhaseResponseRow {
    double frequency{};
    std::vector<std::complex<double>> horizontal, vertical;
};
class PhaseResponseWidget final : public QWidget {
  public:
    explicit PhaseResponseWidget(QWidget *parent = nullptr);
    ~PhaseResponseWidget() override;
    PhaseDisplaySettings settings() const noexcept {
        return settings_;
    }
    void setSettings(const PhaseDisplaySettings &settings);
    void setRows(std::vector<double> angles, std::vector<PhaseResponseRow> rows, QString context);
    std::function<void()> onChanged;

  private:
    void refresh();
    void controlsChanged();
    PhaseDisplaySettings settings_;
    std::vector<double> angles_;
    std::vector<PhaseResponseRow> rows_;
    QDoubleSpinBox *angle_{}, *reference_{}, *delay_{};
    QCheckBox *relative_{}, *unwrap_{};
    PlotWidget *horizontal_{}, *vertical_{};
    QLabel *context_{};
};
} // namespace hornsculpt
