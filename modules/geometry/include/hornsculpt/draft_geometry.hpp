#pragma once

#include <hornsculpt/mesh.hpp>
#include <hornsculpt/mouth_morph.hpp>
#include <hornsculpt/osse_profile.hpp>
#include <hornsculpt/radial_curve.hpp>
#include <hornsculpt/shape_transition.hpp>
#include <array>
#include <memory>
#include <optional>
#include <stop_token>

namespace hornsculpt {

enum class DraftFace : std::uint32_t {
    innerWall = 1,
    outerWall = 2,
    mouthRim = 3,
    rear = 4,
    driven = 5,
    adapterWall = 6,
    slotWall = 7,
    upstreamOuterWall = 8
};

struct EntryParameters {
    double width{0.0254};
    double height{0.0254};
    double exponent{2.0};
};

enum class SlotContour { superellipse, roundedRectangle };
enum class SlotExit { axialBlend, circularFillet };

// All dimensions in metres. The aperture is an internal open cross-section;
// the source remains on EntryParameters at z=0. In axialBlend mode,
// transitionLength controls axial smoothing; in circularFillet mode,
// exitRadius is the physical meridian-circle radius. Zero in the active
// exit control retains the one-sided profile tangents. cornerRadius applies
// only to roundedRectangle and is independent of either exit control.
struct DiffractionSlotParameters {
    bool enabled{false};
    EntryParameters aperture{.012, .035, 4};
    double adapterLength{.02};
    double straightLength{.003};
    double transitionLength{.005};
    SlotContour contour{SlotContour::superellipse};
    double cornerRadius{.002};
    SlotExit exit{SlotExit::axialBlend};
    double exitRadius{.002};
};

// Exact area of the aperture's selected contour, in square metres.
double diffractionSlotArea(const DiffractionSlotParameters &slot);

// Circle in a radial (z,r) meridian, not a general 3D rolling-ball fillet.
// start/end are absolute axial positions; the centre is (start, centreRadius).
struct SlotExitFillet {
    double start{}, end{}, centreRadius{}, radius{};
};
struct SlotExitCache;

enum class AxialDerivativeSide { unspecified, before, after };

struct DraftParameters {
    EntryParameters entry;
    DiffractionSlotParameters slot;
    double length{0.18};
    OsseDirectionalParameters profile;
    OsseDirectionalParameters verticalProfile;
    double horizontalAngle{0.6981317007977318};
    double verticalAngle{0.4363323129985824};
    // Positive X/Y use the horizontal/vertical controls. In symmetric mode
    // these also define negative X/Y; inactive opposite-side values are kept.
    bool asymmetricEnabled{false};
    OsseDirectionalParameters leftProfile;
    OsseDirectionalParameters bottomProfile;
    double leftAngle{0.6981317007977318};
    double bottomAngle{0.4363323129985824};
    double crossSectionExponent{2.0};
    bool guideEnabled{false};
    GuideParameters guide;
    ShapeTransitionParameters transition;
    MouthMorphParameters mouthMorph;
    // Radial material thickness; the rear closure is this far behind the throat.
    double wallThickness{0.006};
};

class DraftGeometry {
  public:
    explicit DraftGeometry(DraftParameters parameters = {});
    DraftGeometry(DraftParameters parameters, std::shared_ptr<const LongitudinalProfile> profile);
    DraftGeometry(DraftParameters parameters,
                  std::shared_ptr<const LongitudinalProfile> horizontalProfile,
                  std::shared_ptr<const LongitudinalProfile> verticalProfile);
    // Providers follow +X, +Y, -X, -Y. With a slot, providers start at its
    // aperture with the remaining length; otherwise at the source entry.
    DraftGeometry(DraftParameters parameters,
                  std::shared_ptr<const LongitudinalProfile> rightProfile,
                  std::shared_ptr<const LongitudinalProfile> topProfile,
                  std::shared_ptr<const LongitudinalProfile> leftProfile,
                  std::shared_ptr<const LongitudinalProfile> bottomProfile);
    double radius(double z, double phi) const;
    // An unspecified derivative is NaN at an actual axial tangent jump.
    double radiusDerivative(double z, double phi,
                            AxialDerivativeSide side = AxialDerivativeSide::unspecified) const;
    // Physical segment boundaries in metres, including entry and mouth.
    std::vector<double> axialBreakpoints() const;
    // Meridian-specific tangent stations; an empty result means no circular fillet.
    std::optional<SlotExitFillet> slotExitFillet(double phi) const;
    // Last station whose complete circumference is still in the straight slot.
    // The slot-wall face ends here; the horn-wall face includes the full fillet.
    double slotStraightEnd() const noexcept;
    // End of the actual inlet composition, including its rounded exit.
    double slotTransitionEnd() const noexcept;
    double length() const noexcept;
    double guidePosition() const;
    double guideRadius(double phi) const;
    const DraftParameters &parameters() const noexcept;
    // Display tessellation only. The BEM meshing backend must remesh and validate it.
    // Metres; adaptive chord and triangle-deviation checks cover the full composition.
    SurfaceMesh previewMesh(double geometricTolerance = 0.00025, std::stop_token stop = {}) const;

  private:
    RadialSample crossSection(double z, double phi, double exponent,
                              bool useSlotContour = false, bool withSlope = true,
                              const std::array<double, 2> *cachedAxes = nullptr) const;
    void initializeSlotExit();
    SlotExitFillet solveSlotExitFillet(double phi) const;
    std::optional<RadialSample> filletSample(double z, double phi) const;
    RadialSample guidedSample(double z, double phi, double targetGuideRadius) const;
    RadialSample baseSample(double z, double phi, bool withSlope = true,
                            const std::array<double, 2> *cachedAxes = nullptr) const;
    RadialSample sample(double z, double phi) const;
    double profilePosition(double z) const;
    double profilePositionDerivative(double z) const;
    double profileStart() const noexcept;
    std::optional<RadialSample> upstreamSample(double z, double phi) const;
    DraftParameters parameters_;
    std::shared_ptr<const LongitudinalProfile> horizontalProfile_;
    std::shared_ptr<const LongitudinalProfile> verticalProfile_;
    std::shared_ptr<const LongitudinalProfile> leftProfile_;
    std::shared_ptr<const LongitudinalProfile> bottomProfile_;
    std::optional<PhysicalContour> guideCurve_;
    ShapeTransition transition_;
    std::optional<MouthMorph> mouthMorph_;
    // Thread-safe bounded memoization of pure per-meridian circle solutions.
    std::shared_ptr<SlotExitCache> slotExitCache_;
    double slotFilletStart_{}, slotFilletEnd_{};
};

} // namespace hornsculpt
