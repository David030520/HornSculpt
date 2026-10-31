#include "edit_history.hpp"

#include <QUndoStack>
#include <cstdio>
#include <stdexcept>

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
} // namespace

int main() {
    using namespace hornsculpt;
    try {
        EditState original;
        original.cad =
            CadAssignmentState{"doc",
                               "session",
                               8,
                               {{1, BoundaryCondition::wall, {1, 0}, CadFaceState::newFace}},
                               {9}};
        auto changed = original;
        require(sameEditState(original, changed), "Unchanged edits must compare equal.");
        changed.parameters.verticalProfile.k = 2;
        require(!sameEditState(original, changed), "Vertical profile edits must enter history.");
        for (const auto edit : {0, 1, 2, 3, 4}) {
            changed = original;
            switch (edit) {
            case 0:
                changed.parameters.asymmetricEnabled = true;
                break;
            case 1:
                changed.parameters.leftAngle = .3;
                break;
            case 2:
                changed.parameters.bottomAngle = .2;
                break;
            case 3:
                changed.parameters.leftProfile.k = .8;
                break;
            case 4:
                changed.parameters.bottomProfile.terminationScale = .3;
                break;
            }
            require(!sameEditState(original, changed),
                    "Four-direction Draft controls must enter history, including dormant values.");
        }
        changed = original;
        changed.parameters.entry.height *= 2;
        require(!sameEditState(original, changed), "Entry dimensions must enter history.");
        for (const auto edit : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}) {
            changed = original;
            auto &slot = changed.parameters.slot;
            switch (edit) {
            case 0:
                slot.enabled = true;
                break;
            case 1:
                slot.aperture.width *= 1.1;
                break;
            case 2:
                slot.aperture.height *= 1.1;
                break;
            case 3:
                slot.aperture.exponent += 1;
                break;
            case 4:
                slot.adapterLength *= 1.1;
                break;
            case 5:
                slot.straightLength = 0;
                break;
            case 6:
                slot.transitionLength = 0;
                break;
            case 7:
                slot.contour = SlotContour::roundedRectangle;
                break;
            case 8:
                slot.cornerRadius = .004;
                break;
            case 9:
                slot.exit = SlotExit::circularFillet;
                break;
            case 10:
                slot.exitRadius = .018;
                break;
            }
            require(!sameEditState(original, changed),
                    "Every slot control must enter history, including dormant dimensions.");
        }
        changed = original;
        changed.parameters.transition.knots[1].blend = .9;
        require(!sameEditState(original, changed), "LUT edits must enter history.");
        changed = original;
        changed.parameters.mouthMorph.target.superformula.n3 = 3;
        require(!sameEditState(original, changed), "Independent mouth shape must enter history.");
        changed = original;
        changed.meshing.adaptive.minimumEdgeLength *= 2;
        require(!sameEditState(original, changed), "Mesh controls must enter history.");
        changed = original;
        changed.meshing.quality.enabled = false;
        require(!sameEditState(original, changed), "Mesh quality mode must enter history.");
        changed = original;
        changed.meshing.quality.minimumAngleDegrees = 10;
        require(!sameEditState(original, changed), "Mesh angle constraint must enter history.");
        changed = original;
        changed.meshing.quality.maximumAspectRatio = 100;
        require(!sameEditState(original, changed), "Mesh aspect constraint must enter history.");
        changed = original;
        changed.draftFaceTargets[1] = .002;
        require(!sameEditState(original, changed), "Draft local sizing must enter history.");
        changed = original;
        changed.cadFaceTargets[1] = .002;
        require(!sameEditState(original, changed), "CAD local sizing must enter history.");
        changed = original;
        changed.draftReversedFaceIds = {1, 5};
        require(!sameEditState(original, changed), "Draft orientation must enter history.");
        changed = original;
        changed.cadReversedFaceIds = {1};
        require(!sameEditState(original, changed), "CAD orientation must enter history.");
        changed = original;
        changed.symmetryXYEnabled = true;
        changed.symmetryPlaneTolerance = .000002;
        require(!sameEditState(original, changed), "XY symmetry must enter history.");
        changed = original;
        changed.acousticModel = {
            {{7, "Innen", BemDomainKind::interior, {1}}, {8, "Außen", BemDomainKind::exterior}},
            {{1, "Öffnung", 7, 8, {2}}},
            std::nullopt};
        require(!sameEditState(original, changed), "Adding an acoustic region must enter history.");
        auto regionEdit = changed;
        regionEdit.acousticModel.interfaces.front().faceIds = {3};
        require(!sameEditState(changed, regionEdit),
                "Changing an interface assignment must enter history.");
        regionEdit = changed;
        regionEdit.acousticModel.domains.front().needsAssignmentReview = true;
        require(!sameEditState(changed, regionEdit),
                "Changing region review state must enter history.");
        changed = original;
        changed.cad->faces[0].boundary = BoundaryCondition::driven;
        changed.cad->faces[0].normalVelocity = {.2, -.3};
        changed.cad->faces[0].state = CadFaceState::confirmed;
        changed.cad->missing.clear();
        changed.meshing.quality = {false, 10, 100};
        changed.draftReversedFaceIds = {5};
        changed.cadReversedFaceIds = {1};
        changed.acousticModel = {{{7, "Innen", BemDomainKind::interior, {1}, {}, true},
                                  {8, "Außen", BemDomainKind::exterior}},
                                 {{1, "Öffnung", 7, 8, {2}}},
                                 std::nullopt};
        changed.symmetryXYEnabled = true;
        changed.symmetryPlaneTolerance = .000002;
        changed.parameters.asymmetricEnabled = true;
        changed.parameters.leftAngle = .3;
        changed.parameters.bottomProfile.k = .8;
        changed.parameters.slot.enabled = true;
        changed.parameters.slot.aperture.width = .01;
        changed.parameters.slot.transitionLength = 0;
        changed.parameters.slot.contour = SlotContour::roundedRectangle;
        changed.parameters.slot.cornerRadius = .003;
        changed.parameters.slot.exit = SlotExit::circularFillet;
        changed.parameters.slot.exitRadius = .018;
        require(!sameEditState(original, changed) && sameCadGeometry(original.cad, changed.cad),
                "Assignments change on the same CAD geometry.");
        auto newerCad = changed.cad;
        ++newerCad->revision;
        require(!sameCadGeometry(original.cad, newerCad),
                "A new CAD revision cannot replay old assignments.");
        newerCad = changed.cad;
        newerCad->sessionId = "other-session";
        require(!sameCadGeometry(original.cad, newerCad),
                "CAD session identity must remain bound.");
        EditState actual = changed;
        unsigned calls = 0;
        QUndoStack stack;
        stack.push(new EditCommand(original, changed, [&](const auto &state) {
            actual = state;
            ++calls;
        }));
        require(calls == 0 && sameEditState(actual, changed),
                "Pushing an existing edit must not apply it twice.");
        stack.undo();
        require(calls == 1 && sameEditState(actual, original),
                "Undo must restore the complete small edit state.");
        stack.redo();
        require(calls == 2 && sameEditState(actual, changed),
                "Redo must restore assignments and sizing together.");
        stack.undo();
        auto branch = original;
        branch.parameters.length += .01;
        actual = branch;
        stack.push(new EditCommand(original, branch, [&](const auto &state) { actual = state; }));
        require(!stack.canRedo(), "A new edit must discard the abandoned redo branch.");
        std::puts("Edit history checks passed.");
    } catch (const std::exception &error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
