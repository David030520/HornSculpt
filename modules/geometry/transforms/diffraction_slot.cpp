#include "../draft_slot.hpp"
#include <hornsculpt/draft_geometry.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <mutex>
#include <numbers>
#include <stdexcept>
#include <string>

namespace hornsculpt {
struct SlotExitCache {
    std::mutex mutex;
    std::map<double, SlotExitFillet> meridians;
};
namespace {
double apertureRadius(const DiffractionSlotParameters &slot, double phi) {
    return slot.contour == SlotContour::roundedRectangle
               ? RoundedRectangle(slot.aperture.width * .5, slot.aperture.height * .5,
                                  slot.cornerRadius)
                     .radius(phi)
               : Superellipse(slot.aperture.width * .5, slot.aperture.height * .5,
                              slot.aperture.exponent)
                     .radius(phi);
}

void validateTransition(const DraftParameters &p, double finish) {
    const auto location = " Anschlussende z=" + std::to_string(finish * 1000) + " mm (" +
                          std::to_string(finish / p.length * 100) + " % L).";
    if (ShapeTransition(p.transition).blend(finish / p.length) != 0)
        throw std::invalid_argument(
            "Spalt: Formuebergang muss hinter Spalt und Ausgangsuebergang beginnen. "
            "Formbeginn oder LUT ausdruecklich anpassen." +
            location);
    if (p.mouthMorph.enabled && p.mouthMorph.startFraction < finish / p.length)
        throw std::invalid_argument("Spalt: Mund-Morphing darf den Halsaufbau nicht ueberlagern." +
                                    location);
    const double guideLimit =
        p.guide.positionMode == GuidePositionMode::relative ? finish / p.length : finish;
    if (p.guideEnabled && p.guide.position <= guideLimit)
        throw std::invalid_argument("Spalt: Leitkurve muss hinter dem Ausgangsuebergang liegen." +
                                    location);
}
} // namespace

double diffractionSlotArea(const DiffractionSlotParameters &slot) {
    return slot.contour == SlotContour::roundedRectangle
               ? RoundedRectangle(slot.aperture.width * .5, slot.aperture.height * .5,
                                  slot.cornerRadius)
                     .area()
               : Superellipse(slot.aperture.width * .5, slot.aperture.height * .5,
                              slot.aperture.exponent)
                     .area();
}

void validateDiffractionSlot(const DraftParameters &p) {
    const auto &s = p.slot;
    if (!s.enabled)
        return;
    const double blendLength = s.exit == SlotExit::axialBlend ? s.transitionLength : 0;
    if (!std::isfinite(s.adapterLength) || s.adapterLength <= 0 ||
        !std::isfinite(s.straightLength) || s.straightLength < 0 || !std::isfinite(blendLength) ||
        blendLength < 0 || !std::isfinite(p.length) || p.length <= 0 ||
        s.adapterLength + s.straightLength + blendLength >= p.length)
        throw std::invalid_argument(
            "Spalt: Adapterlaenge > 0, Spalt-/Uebergangslaenge >= 0 und positive "
            "Restlaenge bis zum Mund erforderlich.");
    if (s.contour != SlotContour::superellipse && s.contour != SlotContour::roundedRectangle)
        throw std::invalid_argument("Spalt: unbekannter Konturmodus.");
    if (s.exit != SlotExit::axialBlend && s.exit != SlotExit::circularFillet)
        throw std::invalid_argument("Spalt: unbekannter Ausgangsmodus.");
    if (s.contour == SlotContour::superellipse &&
        (!std::isfinite(s.aperture.exponent) || s.aperture.exponent < 2))
        throw std::invalid_argument("Spalt: Superellipsen-Exponent muss mindestens 2 sein.");
    (void)diffractionSlotArea(s);
    if (s.exit == SlotExit::circularFillet && (!std::isfinite(s.exitRadius) || s.exitRadius < 0))
        throw std::invalid_argument("Spalt: Ausgangsradius muss endlich und mindestens 0 sein.");
    validateTransition(p, s.adapterLength + s.straightLength + blendLength);
}

double DraftGeometry::profileStart() const noexcept {
    return parameters_.slot.enabled
               ? parameters_.slot.adapterLength + parameters_.slot.straightLength
               : 0;
}

double DraftGeometry::profilePosition(double z) const {
    const double u = z - profileStart();
    const auto &slot = parameters_.slot;
    const double t = slot.enabled && slot.exit == SlotExit::axialBlend ? slot.transitionLength : 0;
    if (t == 0 || u >= t)
        return u;
    const double x = u / t;
    // Quintic axial reparameterization: g(0)=g'(0)=g''(0)=0,
    // g(t)=t, g'(t)=1, g''(t)=0. It is distinct from a physical circular fillet.
    return t * x * x * x * (6 + x * (-8 + 3 * x));
}

double DraftGeometry::profilePositionDerivative(double z) const {
    const double u = z - profileStart();
    const auto &slot = parameters_.slot;
    const double t = slot.enabled && slot.exit == SlotExit::axialBlend ? slot.transitionLength : 0;
    if (t == 0 || u >= t)
        return 1;
    const double x = u / t;
    return x * x * (18 + x * (-32 + 15 * x));
}

SlotExitFillet DraftGeometry::solveSlotExitFillet(double phi) const {
    const auto &s = parameters_.slot;
    const double r = s.exitRadius, start = profileStart(), throat = apertureRadius(s, phi);
    // Position tolerance is a floating-point geometry tolerance in metres,
    // not a meshing tolerance or a licence to shrink the requested radius.
    const double distanceTolerance =
        256 * std::numeric_limits<double>::epsilon() * std::max({length(), r, throat});
    const auto sample = [&](double u) {
        return crossSection(start + u, phi, s.aperture.exponent, true);
    };
    const auto equation = [&](double u) {
        const auto p = sample(u);
        if (!std::isfinite(p.slope) || p.slope < 0)
            throw std::invalid_argument(
                "Spalt-Radius: endliche aufweitende Profiltangente erforderlich.");
        const double norm = std::hypot(1., p.slope);
        // Stable 1-cos(theta), avoiding cancellation at small opening angles.
        const double sag = r * (p.slope / norm) * (p.slope / (norm + 1));
        return p.radius - throat - sag;
    };
    const double limit = std::min(r, (length() - start) * (1 - 1e-12));
    double lo = 0, flo = equation(0), hi = 0;
    bool bracketed = false;
    // A tangent original profile has a trivial root at zero. Search for a
    // strictly nontrivial negative-to-positive crossing instead of accepting it.
    if (flo == 0) {
        lo = std::max(limit * 1e-5, distanceTolerance * 16);
        flo = lo < limit ? equation(lo) : 0;
    }
    constexpr unsigned rootSearchIntervals = 128;
    const double initial = lo;
    for (unsigned i = 1; i <= rootSearchIntervals; ++i) {
        hi = std::lerp(initial, limit, static_cast<double>(i) / rootSearchIntervals);
        const double fhi = equation(hi);
        if (flo < -distanceTolerance && fhi >= 0) {
            bracketed = true;
            break;
        }
        // Keep the last definitely negative point; near-zero cancellation is
        // not evidence of a physical tangent-circle intersection.
        if (fhi < -distanceTolerance) {
            lo = hi;
            flo = fhi;
        }
    }
    if (!bracketed)
        throw std::invalid_argument(
            "Spalt-Radius: kein nichttrivialer tangentialer Kreisbogen bei phi=" +
            std::to_string(phi * 180 / std::numbers::pi) +
            " Grad. Radius und OS-SE-Halswinkel/Profil pruefen; ein bereits tangentialer "
            "Anschluss ist keine Radiusverrundung.");
    double u = (lo + hi) * .5;
    for (unsigned iteration = 0; iteration < 80; ++iteration) {
        u = (lo + hi) * .5;
        const double residual = equation(u);
        // A moving rounded corner may give this equation a much larger
        // derivative than the wall slope. Satisfy both the axial bracket and
        // the attachment residual instead of mistaking incomplete convergence
        // for a discontinuous profile. A genuine jump still fails below.
        if (hi - lo <= distanceTolerance &&
            std::abs(residual) <= distanceTolerance * std::max(1., sample(u).slope))
            break;
        if (u == lo || u == hi)
            break;
        if (residual < 0)
            lo = u;
        else
            hi = u;
    }
    const auto end = sample(u);
    // A sharp contour can make the longitudinal tangent discontinuous when
    // its corner crosses this meridian. A sign change at that jump is not a
    // root: verify positional attachment before accepting the tangent circle.
    const double attachmentTolerance = distanceTolerance * std::max(1., end.slope);
    if (std::abs(equation(u)) > attachmentTolerance)
        throw std::invalid_argument(
            "Spalt-Radius: Tangentialgleichung bei phi=" +
            std::to_string(phi * 180 / std::numbers::pi) +
            " Grad nicht erfuellt (moeglicher Profilsprung). Eckradius und "
            "Halsprofil pruefen; kein geschlossener Kreisbogenanschluss moeglich.");
    const double normalizer = std::hypot(1., end.slope);
    const double localCentre = u - r * end.slope / normalizer;
    if (localCentre > distanceTolerance)
        throw std::invalid_argument(
            "Spalt-Radius: Kreisbogen wuerde den geraden Spalt verlaengern statt die "
            "Ausgangskante abzurunden. Radius und Halsprofil pruefen.");
    if (localCentre < -s.straightLength - distanceTolerance)
        throw std::invalid_argument("Spalt-Radius: gerader Spalt zu kurz; benoetigt mindestens " +
                                    std::to_string(-localCentre * 1000) + " mm bei phi=" +
                                    std::to_string(phi * 180 / std::numbers::pi) + " Grad.");
    if (u <= distanceTolerance || end.slope <= 0)
        throw std::invalid_argument("Spalt-Radius: Kreisbogen hat keine endliche Ausdehnung.");
    return {start + localCentre, start + u, throat + r, r};
}

std::optional<SlotExitFillet> DraftGeometry::slotExitFillet(double phi) const {
    if (!std::isfinite(phi))
        throw std::invalid_argument("Spalt-Radius: endlicher Umfangswinkel erforderlich.");
    if (!slotExitCache_)
        return {};
    double key = std::fmod(phi, 2 * std::numbers::pi);
    if (key < 0)
        key += 2 * std::numbers::pi;
    {
        const std::lock_guard lock(slotExitCache_->mutex);
        if (const auto found = slotExitCache_->meridians.find(key);
            found != slotExitCache_->meridians.end())
            return found->second;
    }
    const auto value = solveSlotExitFillet(key);
    {
        const std::lock_guard lock(slotExitCache_->mutex);
        // Bounded memoization keeps geometry snapshots small during large jobs.
        if (slotExitCache_->meridians.size() < 8192)
            slotExitCache_->meridians.emplace(key, value);
    }
    return value;
}

void DraftGeometry::initializeSlotExit() {
    const auto &s = parameters_.slot;
    slotFilletStart_ = slotFilletEnd_ = profileStart();
    if (!s.enabled || s.exit != SlotExit::circularFillet || s.exitRadius == 0)
        return;
    slotExitCache_ = std::make_shared<SlotExitCache>();
    // Locate circumferential extrema of the tangent stations from dense samples,
    // then refine each detected extremum. Evaluation also solves arbitrary meridians.
    constexpr unsigned directions = 720;
    constexpr double turn = 2 * std::numbers::pi;
    std::array<SlotExitFillet, directions> samples;
    for (unsigned i = 0; i < directions; ++i) {
        samples[i] = *slotExitFillet(turn * i / directions);
        slotFilletStart_ = std::min(slotFilletStart_, samples[i].start);
        slotFilletEnd_ = std::max(slotFilletEnd_, samples[i].end);
    }
    const auto refine = [&](unsigned i, bool maximum) {
        double a = turn * (static_cast<double>(i) - 1) / directions;
        double b = turn * (static_cast<double>(i) + 1) / directions;
        const auto value = [&](double phi) {
            const auto fillet = *slotExitFillet(phi);
            return maximum ? fillet.end : -fillet.start;
        };
        // Golden-section maximization of each isolated sampled extremum.
        constexpr double ratio = .6180339887498948482;
        double c = b - ratio * (b - a), d = a + ratio * (b - a);
        double fc = value(c), fd = value(d);
        for (unsigned k = 0; k < 40; ++k) {
            if (fc < fd) {
                a = c;
                c = d;
                fc = fd;
                d = a + ratio * (b - a);
                fd = value(d);
            } else {
                b = d;
                d = c;
                fd = fc;
                c = b - ratio * (b - a);
                fc = value(c);
            }
        }
        return std::max(fc, fd);
    };
    const double stationTolerance = 256 * std::numeric_limits<double>::epsilon() * length();
    for (unsigned i = 0; i < directions; ++i) {
        const auto &p = samples[(i + directions - 1) % directions];
        const auto &v = samples[i];
        const auto &n = samples[(i + 1) % directions];
        if ((v.start < p.start - stationTolerance && v.start <= n.start) ||
            (v.start <= p.start && v.start < n.start - stationTolerance))
            slotFilletStart_ = std::min(slotFilletStart_, -refine(i, false));
        if ((v.end > p.end + stationTolerance && v.end >= n.end) ||
            (v.end >= p.end && v.end > n.end + stationTolerance))
            slotFilletEnd_ = std::max(slotFilletEnd_, refine(i, true));
    }
    validateTransition(parameters_, slotFilletEnd_);
}

double DraftGeometry::slotStraightEnd() const noexcept {
    return slotExitCache_ ? slotFilletStart_ : profileStart();
}

double DraftGeometry::slotTransitionEnd() const noexcept {
    const auto &s = parameters_.slot;
    if (!s.enabled)
        return 0;
    return s.exit == SlotExit::axialBlend ? profileStart() + s.transitionLength : slotFilletEnd_;
}

std::optional<RadialSample> DraftGeometry::filletSample(double z, double phi) const {
    if (!slotExitCache_ || z < profileStart() - parameters_.slot.exitRadius ||
        z > profileStart() + parameters_.slot.exitRadius)
        return {};
    const auto fillet = *slotExitFillet(phi);
    if (z < fillet.start || z > fillet.end)
        return {};
    const double d = z - fillet.start;
    const double square = (fillet.radius - d) * (fillet.radius + d);
    if (square <= 0)
        throw std::invalid_argument("Spalt-Radius: Kreisbogen erreicht eine senkrechte Tangente.");
    const double height = std::sqrt(square);
    const double sag = d * d / (fillet.radius + height);
    return {RadialSample{fillet.centreRadius - fillet.radius + sag, d / height}};
}

std::optional<RadialSample> DraftGeometry::upstreamSample(double z, double phi) const {
    const auto &s = parameters_.slot;
    if (!s.enabled || z > profileStart())
        return {};
    const double aperture = apertureRadius(s, phi);
    if (z >= s.adapterLength)
        return RadialSample{aperture, 0};
    const auto &e = parameters_.entry;
    const double entry = Superellipse(e.width * .5, e.height * .5, e.exponent).radius(phi);
    const double x = z / s.adapterLength;
    // Convex radial interpolation with zero first/second end derivatives.
    const double blend = x * x * x * (10 + x * (-15 + 6 * x));
    const double derivative = 30 * x * x * (1 - x) * (1 - x) / s.adapterLength;
    return RadialSample{std::lerp(entry, aperture, blend), (aperture - entry) * derivative};
}

std::vector<double> DraftGeometry::axialBreakpoints() const {
    std::vector<double> positions{0, length()};
    if (parameters_.slot.enabled) {
        positions.push_back(parameters_.slot.adapterLength);
        // A circular fillet removes the nominal sharp exit: it is neither a
        // face boundary nor a sampling break. Keep the full fillet envelope.
        positions.push_back(slotStraightEnd());
        positions.push_back(slotTransitionEnd());
    }
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
    return positions;
}

double DraftGeometry::radiusDerivative(double z, double phi, AxialDerivativeSide side) const {
    const auto value = sample(z, phi);
    const auto &s = parameters_.slot;
    const bool rounded =
        s.exit == SlotExit::axialBlend ? s.transitionLength != 0 : s.exitRadius != 0;
    if (!s.enabled || rounded || z != profileStart())
        return value.slope;
    const double before = value.slope;
    const double after = crossSection(z, phi, s.aperture.exponent, true).slope;
    if (side == AxialDerivativeSide::before)
        return before;
    if (side == AxialDerivativeSide::after)
        return after;
    return before == after ? before : std::numeric_limits<double>::quiet_NaN();
}
} // namespace hornsculpt
