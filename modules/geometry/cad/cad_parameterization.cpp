#include "cad_parameterization.hpp"

#include <BRepAdaptor_Curve.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BndLib_Add2dCurve.hxx>
#include <Bnd_Box2d.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopExp.hxx>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <mutex>
#include <numbers>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace hornsculpt::cad_detail {
namespace {
void finite(double value) {
    if (!std::isfinite(value))
        throw std::runtime_error("Native CAD-Parametrisierung besitzt nicht endliche Daten.");
}
Vec3 value(const gp_XYZ &p) {
    finite(p.X());
    finite(p.Y());
    finite(p.Z());
    return {p.X(), p.Y(), p.Z()};
}
std::optional<std::array<double, 2>> bounds(double first, double last) {
    if (Precision::IsInfinite(first) || Precision::IsInfinite(last))
        return std::nullopt;
    finite(first);
    finite(last);
    if (first >= last)
        throw std::runtime_error("Native CAD-Flaeche besitzt keinen geordneten Parameterbereich.");
    return std::array{first, last};
}
template <class F> auto evaluate(F f) {
    try {
        return f();
    } catch (const Standard_Failure &error) {
        throw std::runtime_error(std::string("Native CAD-Parameterauswertung: ") +
                                 error.GetMessageString());
    }
}
struct CurveUse {
    Handle(Geom2d_Curve) curve;
    double first{}, last{}, tolerance{};
    std::array<double, 2> vertexTolerances{};
    std::uint32_t faceId{};
};
struct CartesianCap {
    Handle(Geom_SurfaceOfRevolution) surface;
    Handle(Geom_SphericalSurface) sphere;
    gp_Pnt origin;
    gp_Vec axis, x, y;
    double first{}, last{}, maximumRadius{}, poleParameter{};
    int radialDirection{};

    gp_Vec radialVector(const gp_Pnt &p) const {
        const gp_Vec delta(origin, p);
        return delta - axis * delta.Dot(axis);
    }
    double radius(double v) const {
        if (!sphere.IsNull())
            return sphere->Radius() * std::cos(v);
        return radialVector(surface->BasisCurve()->Value(v)).Magnitude();
    }
    double radialRoundoff() const {
        return 128 * std::numeric_limits<double>::epsilon() * maximumRadius;
    }
    ParameterPoint coordinates(const gp_Pnt &p) const {
        const gp_Vec delta(origin, p);
        if (radialVector(p).Magnitude() <= radialRoundoff())
            return {};
        return {delta.Dot(x), delta.Dot(y)};
    }
    std::array<double, 2> native(ParameterPoint uv) const {
        const double r = std::hypot(uv.u, uv.v);
        const double tolerance = radialRoundoff();
        if (r > maximumRadius + tolerance)
            throw std::out_of_range("Cartesian CAD cap query exceeds its natural radial domain.");
        if (r <= tolerance)
            return {0, poleParameter};
        if (!sphere.IsNull()) {
            double u = std::atan2(uv.v, uv.u);
            if (u < 0)
                u += 2 * std::numbers::pi;
            // A strictly hemispherical trim is a one-to-one radial graph.
            // Its Cartesian metric remains regular at the spherical pole.
            return {u, std::copysign(std::acos(std::clamp(r / sphere->Radius(), 0.0, 1.0)),
                                     poleParameter)};
        }
        double a = first, b = last;
        for (;;) {
            const double middle = std::midpoint(a, b), current = radius(middle);
            if (std::abs(current - r) <= tolerance || middle == a || middle == b) {
                const auto base = radialVector(surface->BasisCurve()->Value(middle));
                double u = std::atan2(uv.v, uv.u) - std::atan2(base.Dot(y), base.Dot(x));
                const double period = surface->UPeriod();
                u = std::fmod(u, period);
                if (u < 0)
                    u += period;
                return {u, middle};
            }
            if ((current < r) == (radialDirection > 0))
                a = middle;
            else
                b = middle;
        }
    }
};
std::optional<CartesianCap> cartesianCap(const Handle(Geom_Surface) & surface,
                                         const TopoDS_Face &face) {
    if (const auto sphere = Handle(Geom_SphericalSurface)::DownCast(surface)) {
        double u0, u1, v0, v1;
        BRepTools::UVBounds(face, u0, u1, v0, v1);
        constexpr double pole = std::numbers::pi / 2;
        constexpr double angularRoundoff = 256 * std::numeric_limits<double>::epsilon();
        const bool south = std::abs(v0 + pole) <= angularRoundoff && v1 < -angularRoundoff;
        const bool north = std::abs(v1 - pole) <= angularRoundoff && v0 > angularRoundoff;
        if (!south && !north)
            return {};
        CartesianCap cap;
        cap.sphere = sphere;
        cap.origin = sphere->Location();
        cap.axis = gp_Vec(sphere->Position().Direction());
        cap.x = gp_Vec(sphere->Position().XDirection());
        cap.y = gp_Vec(sphere->Position().YDirection());
        cap.first = v0;
        cap.last = v1;
        cap.poleParameter = south ? -pole : pole;
        cap.radialDirection = south ? 1 : -1;
        cap.maximumRadius = cap.radius(south ? v1 : v0);
        return cap;
    }
    const auto revolution = Handle(Geom_SurfaceOfRevolution)::DownCast(surface);
    if (revolution.IsNull())
        return {};
    auto basis = revolution->BasisCurve();
    const double first = basis->FirstParameter(), last = basis->LastParameter();
    while (const auto trimmed = Handle(Geom_TrimmedCurve)::DownCast(basis))
        basis = trimmed->BasisCurve();
    const auto curve = Handle(Geom_BSplineCurve)::DownCast(basis);
    if (curve.IsNull() || curve->IsRational() || curve->IsPeriodic())
        return {};
    CartesianCap cap;
    cap.surface = revolution;
    cap.origin = revolution->Axis().Location();
    cap.axis = gp_Vec(revolution->Axis().Direction());
    cap.first = first;
    cap.last = last;
    const double firstRadius = cap.radius(cap.first), lastRadius = cap.radius(cap.last);
    cap.maximumRadius = std::max(firstRadius, lastRadius);
    // This round-off test identifies an actual collapsed rotational coordinate;
    // it is unrelated to the requested mesh/model approximation budgets.
    const double positionalRoundoff = cap.radialRoundoff();
    if (!(cap.maximumRadius > 0) || std::min(firstRadius, lastRadius) > positionalRoundoff)
        return {};
    cap.poleParameter = firstRadius < lastRadius ? cap.first : cap.last;
    cap.radialDirection = firstRadius < lastRadius ? 1 : -1;
    const auto outer = curve->Value(firstRadius < lastRadius ? cap.last : cap.first);
    cap.x = cap.radialVector(outer).Normalized();
    cap.y = cap.axis.Crossed(cap.x);
    double previous = 0;
    for (int i = 1; i <= curve->NbPoles(); ++i) {
        const auto delta = cap.radialVector(curve->Pole(i));
        const double r = delta.Dot(cap.x);
        // A planar, strictly monotone non-rational B-spline control polygon
        // certifies a one-to-one radial graph, including every polynomial span.
        if (r < -positionalRoundoff || std::abs(delta.Dot(cap.y)) > positionalRoundoff ||
            (i > 1 && cap.radialDirection * (r - previous) <= 0))
            return {};
        previous = r;
    }
    gp_Pnt pole;
    gp_Vec tangent;
    curve->D1(cap.poleParameter, pole, tangent);
    const double radialSlope = tangent.Dot(cap.x);
    // Dimensionless round-off budget for a smooth pole's zero axial/radial
    // slope. A true cone tip does not satisfy this condition.
    constexpr double poleSlopeRoundoff = 128 * std::numeric_limits<double>::epsilon();
    if (!(std::abs(radialSlope) > 0) ||
        std::abs(tangent.Dot(cap.axis) / radialSlope) > poleSlopeRoundoff)
        return {};
    return cap;
}
struct SourceDomainEvidence {
    Bnd_Box2d bounds;
    double modelTolerance{};
};
struct Face {
    TopoDS_Face shape;
    Handle(Geom_Surface) surface;
    Bnd_Box2d sourceBounds;
    std::vector<SourceDomainEvidence> sourceEvidence;
    std::optional<CartesianCap> cap;
};
bool includes(const Bnd_Box2d &box, double u, double v) {
    double u0, v0, u1, v1;
    box.Get(u0, v0, u1, v1);
    const double roundoff = 256 * std::numeric_limits<double>::epsilon() *
                            std::max({std::abs(u0), std::abs(u1), std::abs(v0), std::abs(v1),
                                      std::abs(u), std::abs(v)});
    return u >= u0 - roundoff && u <= u1 + roundoff && v >= v0 - roundoff && v <= v1 + roundoff;
}
void addEvidence(Face &face, const Bnd_Box2d &box, double tolerance) {
    face.sourceBounds.Add(box);
    face.sourceEvidence.push_back({box, tolerance});
}
double surfaceParameter(double parameter, const std::optional<std::array<double, 2>> &range,
                        const std::optional<double> &period) {
    if (!range)
        return parameter;
    const auto [first, last] = *range;
    if (period && (parameter < first || parameter > last)) {
        const double shifted = parameter - first;
        finite(shifted);
        double wrapped = std::fmod(shifted, *period);
        if (wrapped < 0)
            wrapped += *period;
        parameter = first + wrapped;
    }
    const double parameterRoundoff =
        256 * std::numeric_limits<double>::epsilon() *
        std::max({std::abs(first), std::abs(last), last - first, std::abs(parameter)});
    if (parameter < first - parameterRoundoff || parameter > last + parameterRoundoff) {
        std::ostringstream message;
        message << std::setprecision(17) << "Parameter " << parameter << " ausserhalb [" << first
                << ", " << last << "]; Ueberschreitung "
                << std::max(first - parameter, parameter - last) << ".";
        throw std::out_of_range(message.str());
    }
    // Only endpoint arithmetic roundoff is regularized. Native model
    // tolerances never enlarge this mathematical parameter domain.
    return std::clamp(parameter, first, last);
}
class Parameterization : public SurfaceParameterization {
  public:
    SurfaceMesh mesh;
    std::vector<ParameterTriangle> records;
    std::map<std::uint32_t, Face> faces;
    std::vector<TopoDS_Edge> curves;
    std::vector<std::vector<double>> curvePartitions;
    std::vector<CurveUse> uses;
    double geometricTolerance{};
    const SurfaceMesh &sourceMesh() const noexcept override {
        return mesh;
    }
    std::span<const ParameterTriangle> triangles() const noexcept override {
        return records;
    }
    const Face &face(std::uint32_t id) const {
        const auto found = faces.find(id);
        if (found == faces.end())
            throw std::invalid_argument("Unbekannte native CAD-Flaechen-ID.");
        return found->second;
    }
    ParameterDomain domain(std::uint32_t id) const override {
        return evaluate([&] {
            const auto &f = face(id);
            const auto &surface = f.surface;
            ParameterDomain result;
            if (f.cap) {
                result.uBounds = result.vBounds =
                    std::array{-f.cap->maximumRadius, f.cap->maximumRadius};
                result.orientation =
                    (f.shape.Orientation() == TopAbs_REVERSED ? -1 : 1) * -f.cap->radialDirection;
                return result;
            }
            double u0, u1, v0, v1;
            surface->Bounds(u0, u1, v0, v1);
            result.uBounds = bounds(u0, u1);
            result.vBounds = bounds(v0, v1);
            if (!f.sourceBounds.IsVoid()) {
                double su0, sv0, su1, sv1;
                f.sourceBounds.Get(su0, sv0, su1, sv1);
                if (result.uBounds)
                    result.uBounds = std::array{std::min(u0, su0), std::max(u1, su1)};
                if (result.vBounds)
                    result.vBounds = std::array{std::min(v0, sv0), std::max(v1, sv1)};
            }
            if (surface->IsUPeriodic()) {
                finite(surface->UPeriod());
                result.uPeriod = surface->UPeriod();
            }
            if (surface->IsVPeriodic()) {
                finite(surface->VPeriod());
                result.vPeriod = surface->VPeriod();
            }
            result.orientation = f.shape.Orientation() == TopAbs_REVERSED ? -1 : 1;
            return result;
        });
    }
    ParametricSurfaceSample surface(std::uint32_t id, ParameterPoint uv) const override {
        finite(uv.u);
        finite(uv.v);
        return evaluate([&] {
            const auto &sourceFace = face(id);
            if (sourceFace.cap) {
                const auto &cap = *sourceFace.cap;
                const auto native = cap.native(uv);
                gp_Pnt p;
                gp_Vec nativeDu, nativeDv;
                sourceFace.surface->D1(native[0], native[1], p, nativeDu, nativeDv);
                const double r = std::hypot(uv.u, uv.v);
                double slope = 0;
                if (r > cap.radialRoundoff()) {
                    const auto radial = cap.radialVector(p).Normalized();
                    const double dr = nativeDv.Dot(radial);
                    if (dr == 0)
                        throw std::runtime_error(
                            "Cartesian CAD cap lost its certified radial metric.");
                    slope = nativeDv.Dot(cap.axis) / dr;
                }
                const gp_Vec du = cap.x + cap.axis * (r > 0 ? slope * uv.u / r : 0);
                const gp_Vec dv = cap.y + cap.axis * (r > 0 ? slope * uv.v / r : 0);
                return ParametricSurfaceSample{value(p.XYZ()), value(du.XYZ()), value(dv.XYZ())};
            }
            const auto domain = this->domain(id);
            double u, v;
            try {
                u = surfaceParameter(uv.u, domain.uBounds, domain.uPeriod);
                v = surfaceParameter(uv.v, domain.vBounds, domain.vPeriod);
            } catch (const std::out_of_range &error) {
                std::ostringstream message;
                message << std::setprecision(17) << "CAD-Flaeche " << id << ", UV=(" << uv.u << ", "
                        << uv.v << "): " << error.what();
                const auto append = [&](const char *axis, const auto &range) {
                    message << " " << axis << "=";
                    if (range)
                        message << "[" << (*range)[0] << ", " << (*range)[1] << "]";
                    else
                        message << "unbegrenzt";
                };
                append("U", domain.uBounds);
                append("V", domain.vBounds);
                throw std::out_of_range(message.str());
            }
            gp_Pnt p;
            gp_Vec du, dv;
            const auto &f = face(id);
            f.surface->D1(u, v, p, du, dv);
            double u0, u1, v0, v1;
            f.surface->Bounds(u0, u1, v0, v1);
            const bool extendsU = !f.surface->IsUPeriodic() && (u < u0 || u > u1);
            const bool extendsV = !f.surface->IsVPeriodic() && (v < v0 || v > v1);
            if (extendsU || extendsV) {
                // Periodic trim coordinates can lie outside the basis interval
                // without extrapolation. Only clamp genuinely bounded axes.
                const auto onBasis = f.surface->Value(extendsU ? std::clamp(u, u0, u1) : u,
                                                      extendsV ? std::clamp(v, v0, v1) : v);
                const double residual = p.Distance(onBasis);
                finite(residual);
                // A source-defined trim may slightly exceed its spline basis.
                // Preserve its true value and D1, but never accept extrapolation
                // beyond either independent physical accuracy budget.
                double modelTolerance = BRep_Tool::Tolerance(f.shape);
                for (const auto &evidence : f.sourceEvidence)
                    if (includes(evidence.bounds, u, v))
                        modelTolerance = std::max(modelTolerance, evidence.modelTolerance);
                if (residual > modelTolerance || residual > geometricTolerance) {
                    std::ostringstream message;
                    message << std::setprecision(17) << "CAD-Flaeche " << id << ", UV=(" << u
                            << ", " << v << "): Quelldomaenenrest " << residual
                            << " m ueberschreitet Modelltoleranz " << modelTolerance
                            << " m oder Meshgenauigkeit " << geometricTolerance << " m.";
                    throw std::out_of_range(message.str());
                }
            }
            return ParametricSurfaceSample{value(p.XYZ()), value(du.XYZ()), value(dv.XYZ())};
        });
    }
    bool isInsideTrim(std::uint32_t id, ParameterPoint uv) const override {
        finite(uv.u);
        finite(uv.v);
        return evaluate([&] {
            const auto &f = face(id);
            if (f.cap) {
                const auto native = f.cap->native(uv);
                uv = {native[0], native[1]};
            } else {
                const auto domain = this->domain(id);
                uv = {surfaceParameter(uv.u, domain.uBounds, domain.uPeriod),
                      surfaceParameter(uv.v, domain.vBounds, domain.vPeriod)};
            }
            const auto key = std::tuple{id, uv.u, uv.v};
            {
                const std::scoped_lock lock(trimCacheMutex_);
                if (const auto found = trimCache_.find(key); found != trimCache_.end())
                    return found->second;
            }
            const BRepClass_FaceClassifier classifier(f.shape, gp_Pnt2d(uv.u, uv.v),
                                                      Precision::PConfusion(), Standard_True);
            const bool inside = classifier.State() == TopAbs_IN;
            {
                const std::scoped_lock lock(trimCacheMutex_);
                // Exact keys preserve classification at trim boundaries. This
                // snapshot-local cap bounds memory; misses always classify.
                constexpr std::size_t maximumCachedTrimQueries = 65536;
                if (trimCache_.size() < maximumCachedTrimQueries)
                    trimCache_.emplace(key, inside);
            }
            return inside;
        });
    }
    ParametricCurveSample curve(std::uint32_t id, double t) const override {
        finite(t);
        if (id == 0 || id > curves.size())
            throw std::invalid_argument("Unbekannte native CAD-Kurven-ID.");
        return evaluate([&] {
            BRepAdaptor_Curve curve(curves[id - 1]);
            if (t < curve.FirstParameter() || t > curve.LastParameter())
                throw std::out_of_range("CAD-Kurvenparameter liegt ausserhalb des echten Trims.");
            gp_Pnt p;
            gp_Vec tangent;
            curve.D1(t, p, tangent);
            return ParametricCurveSample{value(p.XYZ()), value(tangent.XYZ())};
        });
    }
    std::span<const double> curvePartition(std::uint32_t id) const override {
        if (id == 0 || id > curvePartitions.size())
            throw std::invalid_argument("Unbekannte native CAD-Kurven-ID.");
        return curvePartitions[id - 1];
    }
    CurveUseSample curveUse(std::uint32_t id, double t) const override {
        finite(t);
        if (id == 0 || id > uses.size())
            throw std::invalid_argument("Unbekannte native CAD-Randverwendungs-ID.");
        const auto &use = uses[id - 1];
        if (t < use.first || t > use.last)
            throw std::out_of_range("CAD-Randparameter liegt ausserhalb der echten Verwendung.");
        return evaluate([&] {
            gp_Pnt2d uv;
            gp_Vec2d derivative;
            use.curve->D1(t, uv, derivative);
            finite(uv.X());
            finite(uv.Y());
            finite(derivative.X());
            finite(derivative.Y());
            std::optional<double> vertexTolerance;
            if (t == use.first)
                vertexTolerance = use.vertexTolerances[0];
            else if (t == use.last)
                vertexTolerance = use.vertexTolerances[1];
            const auto &f = face(use.faceId);
            if (f.cap) {
                gp_Pnt p;
                gp_Vec du, dv;
                f.surface->D1(uv.X(), uv.Y(), p, du, dv);
                const auto tangent = du * derivative.X() + dv * derivative.Y();
                return CurveUseSample{f.cap->coordinates(p),
                                      {tangent.Dot(f.cap->x), tangent.Dot(f.cap->y)},
                                      use.tolerance,
                                      vertexTolerance};
            }
            return CurveUseSample{
                {uv.X(), uv.Y()}, {derivative.X(), derivative.Y()}, use.tolerance, vertexTolerance};
        });
    }

  private:
    mutable std::mutex trimCacheMutex_;
    mutable std::map<std::tuple<std::uint32_t, double, double>, bool> trimCache_;
};
std::uint32_t nextId(std::size_t count) {
    if (count >= std::numeric_limits<std::uint32_t>::max() - 1)
        throw std::length_error("Native CAD-Parametrisierung besitzt zu viele Identitaeten.");
    return static_cast<std::uint32_t>(count + 1);
}
} // namespace
struct ParameterizationBuilder::Data : Parameterization {};
ParameterizationBuilder::ParameterizationBuilder(double geometricTolerance)
    : data_(std::make_unique<Data>()) {
    finite(geometricTolerance);
    if (geometricTolerance <= 0)
        throw std::invalid_argument("CAD-Detailgenauigkeit muss positiv sein.");
    data_->geometricTolerance = geometricTolerance;
}
ParameterizationBuilder::~ParameterizationBuilder() = default;
void ParameterizationBuilder::addFace(std::uint32_t id, const TopoDS_Face &face) {
    auto surface = BRep_Tool::Surface(face);
    while (const auto trimmed = Handle(Geom_RectangularTrimmedSurface)::DownCast(surface))
        surface = trimmed->BasisSurface();
    if (!data_->faces.emplace(id, Face{face, surface, {}, {}, cartesianCap(surface, face)}).second)
        throw std::invalid_argument("Doppelte native CAD-Parameterflaechen-ID.");
}
std::uint32_t ParameterizationBuilder::addCurve(const TopoDS_Edge &edge) {
    const auto id = nextId(data_->curves.size());
    const auto partition = evaluate([&] {
        BRepAdaptor_Curve curve(edge);
        // CN requests every spline span, including knots whose continuity is
        // higher than D1. Those knots still bound independent polynomials.
        TColStd_Array1OfReal intervals(1, curve.NbIntervals(GeomAbs_CN) + 1);
        curve.Intervals(intervals, GeomAbs_CN);
        std::vector<double> result;
        for (int i = intervals.Lower(); i <= intervals.Upper(); ++i) {
            finite(intervals(i));
            if (!result.empty() && intervals(i) <= result.back())
                throw std::runtime_error("Native CAD-Kurvenspannen sind nicht streng geordnet.");
            result.push_back(intervals(i));
        }
        if (result.size() < 2 || result.front() != curve.FirstParameter() ||
            result.back() != curve.LastParameter())
            throw std::runtime_error("Native CAD-Kurvenspannen decken den echten Trim nicht ab.");
        return result;
    });
    data_->curves.push_back(edge);
    data_->curvePartitions.push_back(partition);
    return id;
}
std::uint32_t ParameterizationBuilder::addUse(std::uint32_t curveId, const TopoDS_Edge &edge,
                                              const TopoDS_Face &face, double anchorParameter,
                                              ParameterPoint anchorUV) {
    if (curveId == 0 || curveId > data_->curves.size())
        throw std::invalid_argument("Native CAD-Randverwendung ohne Raumkurve.");
    double first, last;
    auto curve = BRep_Tool::CurveOnSurface(edge, face, first, last);
    BRepAdaptor_Curve spatial(data_->curves[curveId - 1]);
    if (curve.IsNull() || first != spatial.FirstParameter() || last != spatial.LastParameter())
        throw std::runtime_error(
            "CAD-Rand besitzt keine gemeinsame native Raum-/Flaechenparameterachse.");
    if (BRep_Tool::IsClosed(edge, face)) {
        // OCCT stores two native pcurves and two tessellation polygons for a
        // periodic seam. Bind the actual polygon's UV side, not a guessed
        // convention between their independent storage orders.
        auto opposite = edge;
        opposite.Reverse();
        double otherFirst, otherLast;
        const auto other = BRep_Tool::CurveOnSurface(opposite, face, otherFirst, otherLast);
        if (other.IsNull() || otherFirst != first || otherLast != last)
            throw std::runtime_error(
                "Periodische CAD-Naht besitzt keine gemeinsame Parameterachse.");
        const gp_Pnt2d anchor(anchorUV.u, anchorUV.v);
        const double a = curve->Value(anchorParameter).Distance(anchor);
        const double b = other->Value(anchorParameter).Distance(anchor);
        finite(a);
        finite(b);
        const double roundoff = 256 * std::numeric_limits<double>::epsilon() *
                                std::max({1.0, std::abs(anchorUV.u), std::abs(anchorUV.v), a, b});
        if (std::abs(a - b) <= roundoff)
            throw std::runtime_error(
                "Periodische CAD-Randverwendung ist in ihren UV-Daten mehrdeutig.");
        if (b < a)
            curve = other;
    }
    // Keep the source parameterization, including its native SameParameter
    // approximation. Consumers measure S(gamma(t))-C(t); no inversion occurs here.
    const double tolerance = std::max(BRep_Tool::Tolerance(edge), BRep_Tool::Tolerance(face));
    finite(tolerance);
    const auto firstVertex = TopExp::FirstVertex(edge, Standard_False);
    const auto lastVertex = TopExp::LastVertex(edge, Standard_False);
    if (firstVertex.IsNull() || lastVertex.IsNull())
        throw std::runtime_error("Native CAD-Randverwendung besitzt keine Trimendpunkte.");
    const std::array vertexTolerances{BRep_Tool::Tolerance(firstVertex),
                                      BRep_Tool::Tolerance(lastVertex)};
    finite(vertexTolerances[0]);
    finite(vertexTolerances[1]);
    const auto id = nextId(data_->uses.size());
    data_->uses.push_back({curve, first, last, tolerance, vertexTolerances});
    for (auto &[faceId, sourceFace] : data_->faces) {
        (void)faceId;
        if (sourceFace.shape.IsSame(face)) {
            data_->uses.back().faceId = faceId;
            // Zero padding: the immutable halo comes from the native trim's
            // conservative curve hull, never from a query or modelling tolerance.
            Bnd_Box2d box;
            BndLib_Add2dCurve::Add(curve, first, last, 0, box);
            addEvidence(sourceFace, box, tolerance);
            for (const auto &[t, vertexTolerance] :
                 {std::pair{first, vertexTolerances[0]}, std::pair{last, vertexTolerances[1]}}) {
                Bnd_Box2d endpoint;
                endpoint.Add(curve->Value(t));
                addEvidence(sourceFace, endpoint, std::max(tolerance, vertexTolerance));
            }
            break;
        }
    }
    if (data_->uses.back().faceId == 0)
        throw std::logic_error("Native CAD curve use has no registered physical face.");
    return id;
}
std::shared_ptr<const SurfaceParameterization>
ParameterizationBuilder::finish(SurfaceMesh mesh, std::vector<ParameterTriangle> triangles) {
    if (mesh.triangles.size() != triangles.size())
        throw std::logic_error("Native CAD-Parameterdreiecke passen nicht zum Quellnetz.");
    for (std::size_t i = 0; i < triangles.size(); ++i) {
        auto &face = data_->faces.at(mesh.triangles[i].faceId);
        for (const auto uv : triangles[i].uv) {
            finite(uv.u);
            finite(uv.v);
            face.sourceBounds.Add(gp_Pnt2d(uv.u, uv.v));
        }
        for (std::size_t corner = 0; corner < 3; ++corner) {
            if (const auto &edge = triangles[i].edges[corner]) {
                const auto &use = data_->uses.at(edge->useId - 1);
                for (const auto &[t, endpointCorner] :
                     {std::pair{edge->t0, corner}, std::pair{edge->t1, (corner + 1) % 3}}) {
                    if (t != use.first && t != use.last)
                        continue;
                    const auto uv = triangles[i].uv[endpointCorner];
                    Bnd_Box2d endpoint;
                    endpoint.Add(gp_Pnt2d(uv.u, uv.v));
                    const auto vertexTolerance = use.vertexTolerances[t == use.first ? 0 : 1];
                    addEvidence(face, endpoint, std::max(use.tolerance, vertexTolerance));
                }
            }
        }
    }
    for (std::size_t i = 0; i < triangles.size(); ++i) {
        const auto &face = data_->faces.at(mesh.triangles[i].faceId);
        if (face.cap)
            for (auto &uv : triangles[i].uv)
                uv = face.cap->coordinates(face.surface->Value(uv.u, uv.v));
    }
    data_->mesh = std::move(mesh);
    data_->records = std::move(triangles);
    return std::shared_ptr<const SurfaceParameterization>(std::move(data_));
}
} // namespace hornsculpt::cad_detail
