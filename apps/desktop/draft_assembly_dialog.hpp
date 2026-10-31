#pragma once
#include <QDialog>
#include <hornsculpt/draft_instances.hpp>
class QTableWidget;
class QCheckBox;
namespace hornsculpt {
class DecimalInput;
class DraftAssemblyDialog final : public QDialog {
  public:
    DraftAssemblyDialog(DraftAssembly assembly, bool hasCad, QWidget *parent = nullptr);
    DraftAssembly assembly() const;

  private:
    void addRow(const DraftInstance &instance);
    QTableWidget *table_{};
    QCheckBox *includeCad_{};
    DecimalInput *cadSize_{};
    std::uint32_t nextId_{1};
};
} // namespace hornsculpt
