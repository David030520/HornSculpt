#pragma once

#include "mesh_settings.hpp"
#include <QString>
#include <QUndoCommand>
#include <functional>
#include <hornsculpt/acoustic_partition.hpp>
#include <hornsculpt/draft_geometry.hpp>
#include <hornsculpt/draft_instances.hpp>
#include <hornsculpt/freecad_link.hpp>
#include <hornsculpt/lem.hpp>
#include <map>
#include <optional>
#include <set>

namespace hornsculpt {

struct CadAssignmentState {
    QString documentId;
    QString sessionId;
    std::uint64_t revision{};
    std::vector<CadFaceAssignment> faces;
    std::vector<std::uint32_t> missing;
};

struct EditState {
    DraftParameters parameters;
    DraftAssembly assembly;
    MeshSettings meshing;
    std::map<std::uint32_t, double> draftFaceTargets;
    std::map<std::uint32_t, double> cadFaceTargets;
    std::set<std::uint32_t> draftReversedFaceIds;
    std::set<std::uint32_t> cadReversedFaceIds;
    std::optional<CadAssignmentState> cad;
    lem::Network couplingNetwork;
    bool couplingEnabled{};
    std::vector<AcousticRegionDefinition> acousticRegions;
};

bool sameCadGeometry(const std::optional<CadAssignmentState> &a,
                     const std::optional<CadAssignmentState> &b);
bool sameEditState(const EditState &a, const EditState &b);

// The initial edit already happened in the UI; only subsequent undo/redo replays it.
class EditCommand final : public QUndoCommand {
  public:
    EditCommand(EditState before, EditState after, std::function<void(const EditState &)> apply);
    void undo() override;
    void redo() override;

  private:
    EditState before_, after_;
    std::function<void(const EditState &)> apply_;
    bool initial_{true};
};

} // namespace hornsculpt
