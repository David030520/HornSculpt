#pragma once

#include "angular_response_settings.hpp"
#include <QString>
#include <QWidget>
#include <functional>
#include <hornsculpt/acoustic_analysis.hpp>
#include <vector>

class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QLineEdit;

namespace hornsculpt {
class PlotWidget;

struct AngularResponseRow {
    double frequencyHz{};
    std::vector<AcousticsAngularSample> horizontal;
    std::vector<AcousticsAngularSample> vertical;
};

// Owns a presentation snapshot of actual computed rows, never live estimates.
class AngularResponseWidget final : public QWidget {
  public:
    explicit AngularResponseWidget(QWidget *parent = nullptr);
    void setRows(std::vector<AngularResponseRow> rows, QString context = {});
    void setSettings(const AngularResponseSettings &settings);
    const AngularResponseSettings &settings() const noexcept {
        return settings_;
    }
    std::function<void()> onChanged;

  private:
    void selectionChanged();
    void refresh();
    AngularResponseSettings settings_;
    std::vector<AngularResponseRow> rows_;
    QString context_;
    QLineEdit *angles_{};
    QComboBox *normalization_{};
    QDoubleSpinBox *reference_{};
    QLabel *error_{};
    QLabel *coverage_{};
    PlotWidget *horizontal_{};
    PlotWidget *vertical_{};
};
} // namespace hornsculpt
