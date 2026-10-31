#pragma once

#include <QWidget>
#include <functional>
#include <hornsculpt/draft_geometry.hpp>
#include <map>
#include <optional>
#include <set>

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QFormLayout;
class QLabel;
class QTableWidget;
class QTabWidget;
class QPushButton;

namespace hornsculpt {

class ParameterPanel final : public QWidget {
  public:
    explicit ParameterPanel(QWidget *parent = nullptr);
    ~ParameterPanel() override;
    DraftParameters parameters() const;
    void setParameters(const DraftParameters &parameters);
    void setGeometryStatus(const QString &message, bool error = false);
    // The current immutable geometry job supplies the actual circular-join extent.
    void setSlotTransitionEnd(std::optional<double> position);
    std::function<void()> onChanged;

  private:
    QDoubleSpinBox *createNumber(const QString &key, double value, double min, double max,
                                 int decimals, const QString &suffix, const QString &explanation);
    QDoubleSpinBox *addNumber(QFormLayout *form, const QString &key, const QString &label,
                              double value, double min, double max, int decimals,
                              const QString &suffix, const QString &explanation);
    void changed();
    void updateEnabled();
    void updateSlotState();
    DiffractionSlotParameters slotParameters() const;
    void updateGuideState();
    double modelNumber(const QString &key, double unchangedValue,
                       double displayUnitsPerSI = 1) const;
    std::map<QString, QDoubleSpinBox *> numbers_;
    // Display rounding never changes an untouched model parameter.
    DraftParameters model_;
    std::set<QString> editedNumbers_;
    QCheckBox *asymmetricEnabled_{};
    QCheckBox *slotEnabled_{};
    QWidget *slotSettings_{};
    QComboBox *slotContour_{};
    QComboBox *slotExit_{};
    std::optional<double> slotTransitionEnd_;
    QLabel *slotPositionHint_{};
    QLabel *slotMetrics_{};
    QPushButton *slotFormStart_{};
    QTabWidget *sections_{};
    QPushButton *guideAngleHint_{};
    QPushButton *guideShapeHint_{};
    QLabel *geometryStatus_{};
    QLabel *guideMetrics_{};
    QWidget *oppositeProfiles_{};
    QLabel *horizontalProfileLabel_{};
    QLabel *verticalProfileLabel_{};
    QCheckBox *guideEnabled_{};
    QComboBox *curveType_{};
    QComboBox *positionMode_{};
    QComboBox *transitionMethod_{};
    QTableWidget *knots_{};
    QWidget *guideSettings_{};
    QWidget *superellipseSettings_{};
    QWidget *superformulaSettings_{};
    QWidget *rangeSettings_{};
    QWidget *lutSettings_{};
    QCheckBox *mouthEnabled_{};
    QComboBox *mouthCurveType_{};
    QWidget *mouthSettings_{};
    QWidget *mouthSuperellipseSettings_{};
    QWidget *mouthSuperformulaSettings_{};
    bool updating_{};
    int previousPositionMode_{};
};

} // namespace hornsculpt
