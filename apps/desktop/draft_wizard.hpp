#pragma once

#include <QWizard>
#include <hornsculpt/draft_geometry.hpp>

class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QCheckBox;
class QSpinBox;

namespace hornsculpt {

// User-facing lengths are millimetres and opening angles are full angles.
// The conversion to the geometry module's SI/half-angle contract happens here.
struct BasicDraftSpecification {
    double throatDiameterMm{25.4};
    double lengthMm{180};
    double horizontalOpeningDegrees{80};
    double verticalOpeningDegrees{50};
    double crossSectionExponent{2};
    double wallThicknessMm{6};
};

DraftParameters makeBasicDraft(const BasicDraftSpecification &specification);

class DraftWizard final : public QWizard {
  public:
    explicit DraftWizard(QWidget *parent = nullptr);
    DraftParameters parameters() const;
    bool startCalculation() const;
    int meshPreset() const;
    double minimumFrequency() const;
    double maximumFrequency() const;
    int frequencyCount() const;

  protected:
    void initializePage(int id) override;
    bool validateCurrentPage() override;

  private:
    BasicDraftSpecification specification() const;
    QDoubleSpinBox *diameter_{}, *length_{}, *horizontal_{}, *vertical_{}, *wall_{};
    QComboBox *contour_{};
    QComboBox *meshPreset_{};
    QDoubleSpinBox *minimumFrequency_{}, *maximumFrequency_{};
    QSpinBox *frequencyCount_{};
    QCheckBox *startCalculation_{};
    QLabel *summary_{}, *validation_{}, *calculationValidation_{};
};

} // namespace hornsculpt
