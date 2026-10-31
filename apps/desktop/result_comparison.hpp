#pragma once
#include "project_file.hpp"
#include <QDialog>
#include <QImage>
#include <functional>

namespace hornsculpt {
bool canCaptureComparison(const Project &project) noexcept;
class ResultComparisonDialog final : public QDialog {
public:
    explicit ResultComparisonDialog(Project current, QImage currentView = {}, QWidget *parent = nullptr);
    std::function<void(std::vector<NamedProjectSnapshot>)> onSnapshotsChanged;
private:
    struct Implementation;
    std::shared_ptr<Implementation> implementation_;
};
// Exports the given actual project state. An empty image is explicitly noted.
void writeCustomerReport(const QString &path, const Project &project, const QImage &currentView = {},
                         const QString &standName = {});
}
