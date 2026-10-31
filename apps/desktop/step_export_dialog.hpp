#pragma once

#include <QDialog>
#include <QPointer>
#include <QString>
#include <hornsculpt/draft_geometry.hpp>
#include <memory>

class QComboBox;
class QLabel;
class QProgressBar;
class QPushButton;

namespace hornsculpt {
class DecimalInput;
class IntegerInput;
class StepExportJob;
struct DraftStepExportResult;

// The geometry and revision are an immutable export snapshot. Later Draft
// edits never change a running export or its revision label.
class StepExportDialog final : public QDialog {
  public:
    StepExportDialog(std::shared_ptr<const DraftGeometry> geometry, quint64 revision,
                     QWidget *parent = nullptr);
    ~StepExportDialog() override;

    // The caller has confirmed this destination, including replacement of an
    // existing file. Normal UI use obtains it through the save-file dialog.
    void startExport(const QString &confirmedPath);
    void cancelExport();
    bool busy() const noexcept;
    QString lastError() const;
    void reject() override;

  private:
    friend class StepExportJob;
    void chooseDestination();
    void updateControls();
    void finished(DraftStepExportResult result, const QString &error, bool cancelled);

    const std::shared_ptr<const DraftGeometry> geometry_;
    const quint64 revision_;
    QPointer<StepExportJob> job_;
    QComboBox *kind_{};
    DecimalInput *tolerance_{};
    IntegerInput *sections_{};
    QLabel *sectionsLabel_{}, *status_{}, *destination_{};
    QProgressBar *progress_{};
    QPushButton *export_{}, *cancel_{};
    QString targetPath_, error_;
    bool cancellationRequested_{};
};
} // namespace hornsculpt
