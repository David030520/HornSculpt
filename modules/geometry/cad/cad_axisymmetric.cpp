#include <hornsculpt/axisymmetric_geometry.hpp>
#include <hornsculpt/cad_geometry.hpp>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace hornsculpt {
namespace {
Vec3 sub(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
double norm(Vec3 v) {
    return std::hypot(v.x, v.y, v.z);
}
Vec3 vector(gp_Pnt p) {
    return {p.X(), p.Y(), p.Z()};
}
Vec3 vector(gp_Dir p) {
    return {p.X(), p.Y(), p.Z()};
}
bool finite(Vec3 v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Meridianaufbau abgebrochen.");
}
} // namespace

AxisymmetricMeridian makeAxisymmetricMeridian(const CadGeometry &geometry, Vec3 origin, Vec3 axis,
                                              const MeridianSampling &s, std::stop_token stop) {
    cancelled(stop);
    if (!finite(origin) || !finite(axis) || std::abs(norm(axis) - 1) > 1e-12 ||
        !std::isfinite(s.maximumLength) || s.maximumLength <= 0 ||
        !std::isfinite(s.chordTolerance) || s.chordTolerance <= 0 ||
        !std::isfinite(s.rotationalTolerance) || s.rotationalTolerance <= 0)
        throw std::invalid_argument("2D-Rotation: Achse und positive SI-Genauigkeiten prüfen.");
    AxisymmetricMeridian out;
    out.origin = origin;
    out.axis = axis;
    const Vec3 seed = std::abs(axis.x) < .8 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    const auto projection = dot(seed, axis);
    auto radial = Vec3{seed.x - projection * axis.x, seed.y - projection * axis.y,
                       seed.z - projection * axis.z};
    const double radialNorm = norm(radial);
    radial = {radial.x / radialNorm, radial.y / radialNorm, radial.z / radialNorm};
    out.radial = radial;
    const auto transverse = cross(axis, radial);
    const auto local = [&](Vec3 p) {
        const auto d = sub(p, origin);
        return MeridianPoint{dot(d, radial), dot(d, axis)};
    };
    const auto world = [&](MeridianPoint p, double phi) {
        const double c = std::cos(phi), q = std::sin(phi);
        return Vec3{origin.x + p.axial * axis.x + p.radius * (c * radial.x + q * transverse.x),
                    origin.y + p.axial * axis.y + p.radius * (c * radial.y + q * transverse.y),
                    origin.z + p.axial * axis.z + p.radius * (c * radial.z + q * transverse.z)};
    };
    const auto sameAxis = [&](const gp_Ax1 &a) {
        const auto direction = vector(a.Direction());
        const auto d = sub(vector(a.Location()), origin);
        return norm(cross(direction, axis)) < 1e-10 &&
               norm(cross(d, axis)) <= s.rotationalTolerance;
    };
    const auto pointId = [&](MeridianPoint p) {
        if (p.radius < -s.rotationalTolerance)
            throw std::runtime_error("2D-Rotation: negativer Halbschnitt.");
        if (std::abs(p.radius) <= s.rotationalTolerance)
            p.radius = 0;
        for (std::size_t i = 0; i < out.points.size(); ++i)
            if (std::hypot(p.radius - out.points[i].radius, p.axial - out.points[i].axial) <=
                s.rotationalTolerance)
                return static_cast<std::uint32_t>(i);
        if (out.points.size() >= std::numeric_limits<std::uint32_t>::max())
            throw std::length_error("2D-Meridian zu groß.");
        out.points.push_back(p);
        return static_cast<std::uint32_t>(out.points.size() - 1);
    };
    try {
        for (const auto &part : geometry.brepParts()) {
            cancelled(stop);
            std::istringstream stream(part.data);
            BRep_Builder builder;
            TopoDS_Shape shape;
            BRepTools::Read(shape, stream, builder);
            if (shape.IsNull() || stream.bad())
                throw std::runtime_error("2D-Rotation: exakte BRep nicht lesbar.");
            TopTools_IndexedMapOfShape faces;
            TopExp::MapShapes(shape, TopAbs_FACE, faces);
            if (part.faceIds.size() != static_cast<std::size_t>(faces.Extent()))
                throw std::runtime_error("2D-Rotation: Flächenidentitäten fehlen.");
            // Bounds are derived from source geometry and only bound the
            // mathematical cutting sheet; this sheet never enters the mesh.
            double span = 1;
            for (const auto &f : geometry.faces())
                span = std::max(span, 4 * (norm(sub(f.center, origin)) + std::sqrt(f.area)));
            const auto planeNormal = cross(radial, axis);
            const gp_Pln plane(gp_Ax3(gp_Pnt(origin.x, origin.y, origin.z),
                                      gp_Dir(planeNormal.x, planeNormal.y, planeNormal.z),
                                      gp_Dir(radial.x, radial.y, radial.z)));
            const auto sheet = BRepBuilderAPI_MakeFace(plane, 0, span, -span, span).Face();
            for (int f = 1; f <= faces.Extent(); ++f) {
                cancelled(stop);
                const auto face = TopoDS::Face(faces(f));
                const auto id = part.faceIds[static_cast<std::size_t>(f - 1)];
                const BRepAdaptor_Surface surface(face, Standard_True);
                bool supported = false;
                switch (surface.GetType()) {
                case GeomAbs_Plane:
                    supported =
                        norm(cross(vector(surface.Plane().Axis().Direction()), axis)) < 1e-10;
                    break;
                case GeomAbs_Cylinder:
                    supported = sameAxis(surface.Cylinder().Axis());
                    break;
                case GeomAbs_Cone:
                    supported = sameAxis(surface.Cone().Axis());
                    break;
                case GeomAbs_Torus:
                    supported = sameAxis(surface.Torus().Axis());
                    break;
                case GeomAbs_Sphere:
                    supported = norm(cross(sub(vector(surface.Sphere().Location()), origin),
                                           axis)) <= s.rotationalTolerance;
                    break;
                case GeomAbs_SurfaceOfRevolution:
                    supported = sameAxis(surface.AxeOfRevolution());
                    break;
                default:
                    break;
                }
                if (!supported)
                    throw std::invalid_argument(
                        "2D-Rotation: Fläche " + std::to_string(id) +
                        " ist keine nachgewiesene koaxiale Rotationsfläche.");
                BRepAlgoAPI_Section section(face, sheet, Standard_False);
                section.Approximation(Standard_False);
                section.Build();
                if (!section.IsDone())
                    throw std::runtime_error("2D-Rotation: exakter Halbschnitt fehlgeschlagen.");
                bool found = false;
                const auto target =
                    s.faceTargetLengths.contains(id) ? s.faceTargetLengths.at(id) : s.maximumLength;
                if (!std::isfinite(target) || target <= 0)
                    throw std::invalid_argument("2D-Rotation: ungültige lokale Abschnittslänge.");
                for (TopExp_Explorer e(section.Shape(), TopAbs_EDGE); e.More(); e.Next()) {
                    const BRepAdaptor_Curve curve(TopoDS::Edge(e.Current()));
                    const auto first = curve.FirstParameter(), last = curve.LastParameter();
                    if (!std::isfinite(first) || !std::isfinite(last) || last <= first)
                        throw std::runtime_error("2D-Rotation: unbeschränkte Schnittkurve.");
                    std::vector<MeridianPoint> points{local(vector(curve.Value(first)))};
                    std::function<void(double, double, unsigned)> sample = [&](double a, double b,
                                                                               unsigned depth) {
                        cancelled(stop);
                        const auto p = local(vector(curve.Value(a))),
                                   q = local(vector(curve.Value(b)));
                        double deviation = 0;
                        for (double t : {.25, .5, .75}) {
                            const auto x = local(vector(curve.Value(a + (b - a) * t)));
                            deviation =
                                std::max(deviation,
                                         std::hypot(x.radius - ((1 - t) * p.radius + t * q.radius),
                                                    x.axial - ((1 - t) * p.axial + t * q.axial)));
                        }
                        if (std::hypot(q.radius - p.radius, q.axial - p.axial) >
                                std::min(target, s.maximumLength) ||
                            deviation > s.chordTolerance) {
                            if (depth >= 24)
                                throw std::runtime_error(
                                    "2D-Rotation: Meridian erreicht die Genauigkeit nicht.");
                            const double mid = (a + b) * .5;
                            sample(a, mid, depth + 1);
                            sample(mid, b, depth + 1);
                        } else
                            points.push_back(q);
                    };
                    sample(first, last, 0);
                    if (points.size() < 2)
                        continue;
                    for (const auto p : points)
                        for (unsigned angle = 0; angle < 13; ++angle) {
                            const auto rotated = world(p, 2 * std::numbers::pi * angle / 13);
                            if (norm(sub(geometry.projectSurface(id, rotated).position, rotated)) >
                                s.rotationalTolerance)
                                throw std::invalid_argument("2D-Rotation: Trimmung einer Fläche "
                                                            "ist nicht rotationssymmetrisch.");
                        }
                    for (std::size_t i = 1; i < points.size(); ++i) {
                        const auto p = points[i - 1], q = points[i];
                        if (std::hypot(q.radius - p.radius, q.axial - p.axial) <=
                            s.rotationalTolerance)
                            continue;
                        if (p.radius == 0 && q.radius == 0)
                            continue;
                        const MeridianPoint middle{(p.radius + q.radius) * .5,
                                                   (p.axial + q.axial) * .5};
                        const auto normal = geometry.projectSurface(id, world(middle, 0)).normal;
                        const double orientation = (q.axial - p.axial) * dot(normal, radial) -
                                                   (q.radius - p.radius) * dot(normal, axis);
                        auto a = pointId(p), b = pointId(q);
                        if (orientation < 0)
                            std::swap(a, b);
                        if (orientation == 0)
                            throw std::runtime_error(
                                "2D-Rotation: Schnittorientierung unbestimmt.");
                        out.segments.push_back({a, b, id,
                                                s.faceConditions.contains(id)
                                                    ? s.faceConditions.at(id)
                                                    : BoundaryCondition::wall});
                    }
                    found = true;
                }
                if (!found)
                    throw std::invalid_argument(
                        "2D-Rotation: eine Fläche besitzt keinen vollständigen Meridian.");
            }
        }
    } catch (const Standard_Failure &error) {
        throw std::runtime_error(std::string("2D-CAD-Schnitt: ") + error.GetMessageString());
    }
    return out;
}
} // namespace hornsculpt
