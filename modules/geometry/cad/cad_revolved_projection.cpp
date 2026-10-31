#include "cad_revolved_projection.hpp"
#include <Extrema_ExtPC.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Precision.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace hornsculpt::cad_detail {
std::optional<RevolvedProjection> RevolvedProjection::make(const Handle(Geom_Surface) & input,
                                                           std::array<double, 4> bounds) {
    auto surface = input;
    while (const auto trim = Handle(Geom_RectangularTrimmedSurface)::DownCast(surface))
        surface = trim->BasisSurface();
    const auto revolution = Handle(Geom_SurfaceOfRevolution)::DownCast(surface);
    if (revolution.IsNull())
        return {};
    const double period = revolution->UPeriod();
    const double angularRoundoff = 128 * std::numeric_limits<double>::epsilon() * period;
    if (!std::all_of(bounds.begin(), bounds.end(), [](double v) { return std::isfinite(v); }) ||
        bounds[1] - bounds[0] < period - angularRoundoff)
        return {};
    auto basis = revolution->BasisCurve();
    while (const auto trim = Handle(Geom_TrimmedCurve)::DownCast(basis))
        basis = trim->BasisCurve();
    const auto spline = Handle(Geom_BSplineCurve)::DownCast(basis);
    if (spline.IsNull())
        return {};
    RevolvedProjection result;
    result.surface_ = revolution;
    result.curve_ = revolution->BasisCurve();
    result.origin_ = revolution->Axis().Location();
    result.axis_ = gp_Vec(revolution->Axis().Direction());
    result.bounds_ = bounds;
    const auto radial = [&](const gp_Pnt &p) {
        const gp_Vec delta(result.origin_, p);
        return delta - result.axis_ * delta.Dot(result.axis_);
    };
    for (int i = 1; i <= spline->NbPoles(); ++i) {
        const auto r = radial(spline->Pole(i));
        if (r.Magnitude() > result.radialScale_) {
            result.radialScale_ = r.Magnitude();
            result.radial_ = r;
        }
    }
    if (!(result.radialScale_ > 0))
        return {};
    result.radial_.Normalize();
    result.transverse_ = result.axis_.Crossed(result.radial_);
    // Positive B-spline weights preserve the convex hull: the entire meridian
    // lies in this plane and on one side of the axis. This certifies the exact
    // radial reduction without depending on tessellation or mesh tolerances.
    const double planarRoundoff =
        128 * std::numeric_limits<double>::epsilon() * result.radialScale_;
    for (int i = 1; i <= spline->NbPoles(); ++i) {
        const auto r = radial(spline->Pole(i));
        if (std::abs(r.Dot(result.transverse_)) > planarRoundoff ||
            r.Dot(result.radial_) < -planarRoundoff || !(spline->Weight(i) > 0))
            return {};
    }
    return result;
}

SurfaceMinimum RevolvedProjection::minimum(Vec3 approximate, bool requireUnique) const {
    if (!std::isfinite(approximate.x) || !std::isfinite(approximate.y) ||
        !std::isfinite(approximate.z))
        throw std::invalid_argument("CAD-Koordinaten muessen endlich sein.");
    const gp_Pnt query(approximate.x, approximate.y, approximate.z);
    const gp_Vec delta(origin_, query);
    const double axial = delta.Dot(axis_);
    const auto radial = delta - axis_ * axial;
    const double radius = radial.Magnitude();
    const gp_Pnt meridianQuery = origin_.Translated(radial_ * radius + axis_ * axial);
    const GeomAdaptor_Curve curve(curve_, bounds_[2], bounds_[3]);
    TColStd_Array1OfReal spans(1, curve.NbIntervals(GeomAbs_CN) + 1);
    curve.Intervals(spans, GeomAbs_CN);
    std::vector<SurfaceMinimum> candidates;
    const auto append = [&](double v) {
        const auto base = curve_->Value(v);
        const gp_Vec baseDelta(origin_, base);
        const double targetAngle = std::atan2(radial.Dot(transverse_), radial.Dot(radial_));
        const double baseAngle = std::atan2(baseDelta.Dot(transverse_), baseDelta.Dot(radial_));
        const double period = surface_->UPeriod();
        double u = std::fmod(targetAngle - baseAngle - bounds_[0], period);
        if (u < 0)
            u += period;
        u += bounds_[0];
        const auto position = surface_->Value(u, v);
        if (!std::isfinite(u) || !std::isfinite(v) || !std::isfinite(position.X()) ||
            !std::isfinite(position.Y()) || !std::isfinite(position.Z()))
            throw std::runtime_error("CAD-Meridianprojektion liefert nicht endliche Daten.");
        candidates.push_back({position, u, v});
    };
    // Each native polynomial/continuity span contributes its endpoints and
    // stationary minima. A full-box 2D search can miss a short spline minimum.
    for (int span = spans.Lower(); span < spans.Upper(); ++span) {
        const double first = spans(span), last = spans(span + 1);
        append(first);
        append(last);
        const Extrema_ExtPC projection(meridianQuery, curve, first, last, Precision::PConfusion());
        if (!projection.IsDone())
            throw std::runtime_error("CAD-Meridianprojektion findet keine Kurvenextrema.");
        for (int i = 1; i <= projection.NbExt(); ++i)
            if (projection.IsMin(i)) {
                const double v = projection.Point(i).Parameter();
                if (v >= first && v <= last)
                    append(v);
            }
    }
    SurfaceMinimum best;
    double distance = std::numeric_limits<double>::infinity();
    bool ambiguous = false;
    const double positionalRoundoff = 128 * std::numeric_limits<double>::epsilon() * radialScale_;
    for (const auto &candidate : candidates) {
        const double current = query.Distance(candidate.position);
        const double distanceRoundoff =
            128 * std::numeric_limits<double>::epsilon() *
            std::max({radialScale_, current, std::isfinite(distance) ? distance : 0.});
        if (current < distance - distanceRoundoff) {
            best = candidate;
            distance = current;
            ambiguous = false;
        } else if (std::abs(current - distance) <= distanceRoundoff &&
                   candidate.position.Distance(best.position) > positionalRoundoff)
            ambiguous = true;
    }
    const gp_Vec bestDelta(origin_, best.position);
    const double bestRadius = (bestDelta - axis_ * bestDelta.Dot(axis_)).Magnitude();
    if (requireUnique &&
        (ambiguous || (radius <= positionalRoundoff && bestRadius > positionalRoundoff)))
        throw std::invalid_argument("CAD-Projektion ist geometrisch mehrdeutig.");
    if (!std::isfinite(distance))
        throw std::runtime_error("CAD-Meridianprojektion besitzt kein endliches Minimum.");
    return best;
}
} // namespace hornsculpt::cad_detail
