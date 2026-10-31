#pragma once

#include <QWidget>
#include <cstdint>
#include <functional>
#include <hornsculpt/lem.hpp>
#include <hornsculpt/mesh.hpp>
#include <memory>
#include <optional>
#include <vector>

namespace hornsculpt {

struct LemFaceItem {
    std::uint32_t id{};
    QString label;
    BoundaryCondition boundary{BoundaryCondition::driven};
    bool interfaceAllowed{true};
    bool operator==(const LemFaceItem &) const = default;
};

// The editor owns an editable network value. Solver jobs receive their own immutable copy.
class LemEditor final : public QWidget {
  public:
    explicit LemEditor(QWidget *parent = nullptr);
    ~LemEditor() override;

    const lem::Network &network() const;
    void setNetwork(lem::Network network);
    void setFaces(std::vector<LemFaceItem> faces);
    void setSelectedFaces(std::vector<std::uint32_t> faces);
    bool couplingEnabled() const;
    void setCouplingEnabled(bool enabled);
    void setSolution(std::optional<lem::Solution> solution, bool stale = false);

    std::function<void()> onChanged;
    std::function<void(bool)> onEnabledChanged;
    std::function<void(std::vector<std::uint32_t>)> onFacesSelected;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace hornsculpt
