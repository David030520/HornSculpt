#pragma once

#include "project_file.hpp"
#include <QDialog>
#include <functional>
#include <memory>
#include <vector>

namespace hornsculpt {
class ParameterSweepDialog final : public QDialog {
  public:
    explicit ParameterSweepDialog(Project baseline, QWidget *parent = nullptr);
    ~ParameterSweepDialog() override;
    void reject() override;
    void openStudy(const QString &manifestPath);
    void configureConvergence(std::vector<double> frequencies);
    std::function<void(Project)> onVariantChosen;
    // The owner opens the validated completed variant after this dialog exits.
    std::function<void(const QString &)> onVariantAnalysisRequested;
    // Called only after Start has validated the study. The completion runs on
    // the GUI thread when the desktop solver has relinquished GPU resources.
    std::function<void(std::function<void(bool)>)> onAcquireGpu;
    std::function<void()> onCancelGpuAcquire;
    std::function<void()> onReleaseGpu;
    std::function<QString()> gpuWaitMessage;

  private:
    struct Implementation;
    std::unique_ptr<Implementation> implementation_;
};
} // namespace hornsculpt
