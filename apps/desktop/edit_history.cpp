#include "edit_history.hpp"

#include <algorithm>
#include <tuple>
#include <utility>

namespace hornsculpt {
namespace {
auto directional(const OsseDirectionalParameters &p) {
    return std::tie(p.throatAngle, p.k, p.terminationScale, p.terminationQ, p.terminationExponent);
}
bool sameParameters(const DraftParameters &a, const DraftParameters &b) {
    const auto base = [](const auto &p) {
        return std::tie(p.entry.width, p.entry.height, p.entry.exponent, p.length,
                        p.horizontalAngle, p.verticalAngle, p.crossSectionExponent, p.guideEnabled,
                        p.wallThickness, p.asymmetricEnabled, p.leftAngle, p.bottomAngle);
    };
    const auto contour = [](const ContourParameters &p) {
        return std::tie(p.type, p.width, p.aspectRatio, p.rotation, p.superellipseExponent,
                        p.superformula.a, p.superformula.b, p.superformula.m1, p.superformula.m2,
                        p.superformula.n1, p.superformula.n2, p.superformula.n3);
    };
    const auto transition = [](const ShapeTransitionParameters &p) {
        return std::tie(p.method, p.startFraction, p.endFraction);
    };
    const auto slot = [](const DiffractionSlotParameters &p) {
        return std::tie(p.enabled, p.aperture.width, p.aperture.height, p.aperture.exponent,
                        p.adapterLength, p.straightLength, p.transitionLength,
                        p.contour, p.cornerRadius, p.exit, p.exitRadius);
    };
    return base(a) == base(b) && slot(a.slot) == slot(b.slot) &&
           directional(a.profile) == directional(b.profile) &&
           directional(a.verticalProfile) == directional(b.verticalProfile) &&
           directional(a.leftProfile) == directional(b.leftProfile) &&
           directional(a.bottomProfile) == directional(b.bottomProfile) &&
           contour(a.guide) == contour(b.guide) &&
           std::tie(a.guide.positionMode, a.guide.position) ==
               std::tie(b.guide.positionMode, b.guide.position) &&
           std::tie(a.mouthMorph.enabled, a.mouthMorph.startFraction, a.mouthMorph.endFraction) ==
               std::tie(b.mouthMorph.enabled, b.mouthMorph.startFraction,
                        b.mouthMorph.endFraction) &&
           contour(a.mouthMorph.target) == contour(b.mouthMorph.target) &&
           transition(a.transition) == transition(b.transition) &&
           a.transition.knots.size() == b.transition.knots.size() &&
           std::equal(a.transition.knots.begin(), a.transition.knots.end(),
                      b.transition.knots.begin(), [](const auto &x, const auto &y) {
                          return x.fraction == y.fraction && x.blend == y.blend;
                      });
}
auto meshSettings(const MeshSettings &s) {
    return std::tie(s.targetEdgeLength, s.adaptive.enabled, s.adaptive.minimumEdgeLength,
                    s.adaptive.geometricTolerance, s.adaptive.growthRate, s.quality.enabled,
                    s.quality.minimumAngleDegrees, s.quality.maximumAspectRatio);
}
} // namespace

bool sameCadGeometry(const std::optional<CadAssignmentState> &a,
                     const std::optional<CadAssignmentState> &b) {
    return a.has_value() == b.has_value() &&
           (!a || std::tie(a->documentId, a->sessionId, a->revision) ==
                      std::tie(b->documentId, b->sessionId, b->revision));
}
bool sameEditState(const EditState &a, const EditState &b) {
    if (a.couplingEnabled != b.couplingEnabled || a.couplingNetwork != b.couplingNetwork)
        return false;
    if (a.assembly != b.assembly || !sameParameters(a.parameters, b.parameters) ||
        meshSettings(a.meshing) != meshSettings(b.meshing) ||
        a.draftFaceTargets != b.draftFaceTargets || a.cadFaceTargets != b.cadFaceTargets ||
        a.draftReversedFaceIds != b.draftReversedFaceIds ||
        a.cadReversedFaceIds != b.cadReversedFaceIds || !sameCadGeometry(a.cad, b.cad))
        return false;
    if (!a.cad)
        return true;
    return a.cad->missing == b.cad->missing && a.cad->faces.size() == b.cad->faces.size() &&
           std::equal(a.cad->faces.begin(), a.cad->faces.end(), b.cad->faces.begin(),
                      [](const auto &x, const auto &y) {
                          return std::tie(x.id, x.boundary, x.normalVelocity, x.state) ==
                                 std::tie(y.id, y.boundary, y.normalVelocity, y.state);
                      });
}
EditCommand::EditCommand(EditState before, EditState after,
                         std::function<void(const EditState &)> apply)
    : before_(std::move(before)), after_(std::move(after)), apply_(std::move(apply)) {}
void EditCommand::undo() {
    apply_(before_);
}
void EditCommand::redo() {
    if (initial_)
        initial_ = false;
    else
        apply_(after_);
}

} // namespace hornsculpt
