#pragma once
#include <QDialog>
#include <QString>
#include <hornsculpt/acoustic_partition.hpp>
#include <map>
#include <vector>

class QLabel;
class QListWidget;
class QTableWidget;

namespace hornsculpt {
class AcousticRegionsDialog final : public QDialog {
  public:
    explicit AcousticRegionsDialog(std::vector<AcousticRegionDefinition> regions,
                                   std::map<std::uint32_t, QString> availableFaceLabels,
                                   QWidget *parent = nullptr);
    [[nodiscard]] const std::vector<AcousticRegionDefinition> &regions() const noexcept;
    void accept() override;

  private:
    void addRegion();
    void removeRegion();
    void useSelectedFaces(int column);
    QTableWidget *table_{};
    QListWidget *faces_{};
    QLabel *error_{};
    std::map<std::uint32_t, QString> availableFaces_;
    std::vector<AcousticRegionDefinition> committed_;
    std::uint64_t nextId_{1};
};
} // namespace hornsculpt
