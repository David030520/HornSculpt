#pragma once
#include <QDialog>
#include <QString>
#include <hornsculpt/acoustic_partition.hpp>
#include <map>

class QLabel;
class QListWidget;
class QTreeWidget;
namespace hornsculpt {
class AcousticRegionsDialog final : public QDialog {
  public:
    explicit AcousticRegionsDialog(AcousticModelDefinition model,
                                   std::map<std::uint32_t, QString> availableFaceLabels,
                                   QWidget *parent = nullptr);
    [[nodiscard]] const AcousticModelDefinition &model() const noexcept;
    void selectFaces(const std::vector<std::uint32_t> &faceIds);
    void accept() override;

  private:
    void refresh();
    void addDomain(BemDomainKind kind);
    void addInterface();
    void editCurrent();
    void removeCurrent();
    void assignFaces(bool twoSided);
    void editBaffle();
    void addGroup();
    AcousticModelDefinition working_, committed_;
    std::map<std::uint32_t, QString> availableFaces_;
    QTreeWidget *tree_{};
    QListWidget *faces_{};
    QLabel *error_{};
};
} // namespace hornsculpt
