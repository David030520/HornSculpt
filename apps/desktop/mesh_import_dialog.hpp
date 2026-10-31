#pragma once
#include "gmsh_import_data.hpp"
#include <QDialog>
#include <QString>
#include <cstdint>
#include <optional>
#include <thread>

class QComboBox;
class QDialogButtonBox;
class QLabel;
class QListWidget;

namespace hornsculpt {
class MeshImportDialog final : public QDialog {
  public:
    explicit MeshImportDialog(QString path, QWidget *parent = nullptr);
    ~MeshImportDialog() override;
    [[nodiscard]] const std::optional<MeshImportSelection> &selection() const noexcept;
    void accept() override;

  private:
    void reload();
    void startImport();
    QString path_;
    std::jthread worker_;
    std::uint64_t generation_{};
    QComboBox *units_{};
    QListWidget *groups_{};
    QLabel *status_{};
    QDialogButtonBox *buttons_{};
    std::optional<MeshImportSelection> candidate_;
    std::optional<MeshImportSelection> committed_;
};
} // namespace hornsculpt
