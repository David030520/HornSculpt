#include "cad_boundary_corners.hpp"

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace hornsculpt::cad_detail {
namespace {
void finite(const gp_XYZ &v) {
    if (!std::isfinite(v.X()) || !std::isfinite(v.Y()) || !std::isfinite(v.Z()))
        throw std::runtime_error("CAD-Eckabfrage liefert nicht endliche Geometriedaten.");
}
Vec3 value(const gp_Pnt &p) {
    return {p.X(), p.Y(), p.Z()};
}
std::optional<std::array<gp_Pnt, 2>> lineEnds(const TopoDS_Edge &edge, const TopoDS_Face &face) {
    if (BRep_Tool::Degenerated(edge) || BRep_Tool::IsClosed(edge, face))
        return std::nullopt;
    if (edge.Orientation() != TopAbs_FORWARD && edge.Orientation() != TopAbs_REVERSED)
        return std::nullopt;
    BRepAdaptor_Curve curve(edge);
    if (curve.GetType() != GeomAbs_Line)
        return std::nullopt;
    const double first = curve.FirstParameter(), last = curve.LastParameter();
    if (!std::isfinite(first) || !std::isfinite(last))
        throw std::runtime_error("CAD-Eckabfrage besitzt keinen endlichen Linienbereich.");
    if (first >= last || TopExp::FirstVertex(edge, true).IsNull() ||
        TopExp::LastVertex(edge, true).IsNull())
        return std::nullopt;
    std::array<gp_Pnt, 2> points{BRep_Tool::Pnt(TopExp::FirstVertex(edge, true)),
                                 BRep_Tool::Pnt(TopExp::LastVertex(edge, true))};
    const bool reversed = edge.Orientation() == TopAbs_REVERSED;
    const std::array<gp_Pnt, 2> evaluated{curve.Value(reversed ? last : first),
                                          curve.Value(reversed ? first : last)};
    double scale = 0;
    for (const auto &p : {points[0], points[1], evaluated[0], evaluated[1]}) {
        finite(p.XYZ());
        scale = std::max({scale, std::abs(p.X()), std::abs(p.Y()), std::abs(p.Z())});
    }
    const double roundoff = 256 * std::numeric_limits<double>::epsilon() * scale;
    // Native vertex/curve approximation beyond floating-point roundoff is not
    // an exact planar line constraint. Do not expand by CAD model tolerances.
    if (points[0].Distance(evaluated[0]) > roundoff ||
        points[1].Distance(evaluated[1]) > roundoff || points[0].Distance(points[1]) <= roundoff)
        return std::nullopt;
    return points;
}
std::optional<PlanarLineConstraint>
constraint(const TopoDS_Edge &incoming, const TopoDS_Edge &outgoing, const TopoDS_Face &face) {
    BRepAdaptor_Surface surface(face, Standard_True);
    if (surface.GetType() != GeomAbs_Plane)
        return std::nullopt;
    const auto before = lineEnds(incoming, face), after = lineEnds(outgoing, face);
    if (!before || !after)
        return std::nullopt;
    const auto plane = surface.Plane();
    const auto normal = plane.Axis().Direction();
    finite(plane.Location().XYZ());
    finite(normal.XYZ());
    double scale = std::max({std::abs(plane.Location().X()), std::abs(plane.Location().Y()),
                             std::abs(plane.Location().Z())});
    for (const auto &p : {(*before)[0], (*before)[1], (*after)[0], (*after)[1]})
        scale = std::max({scale, std::abs(p.X()), std::abs(p.Y()), std::abs(p.Z())});
    const double epsilon = 256 * std::numeric_limits<double>::epsilon();
    for (const auto &p : {(*before)[0], (*before)[1], (*after)[0], (*after)[1]})
        if (plane.Distance(p) > epsilon * scale)
            return std::nullopt;
    for (auto ray : {gp_Vec((*before)[1], (*before)[0]), gp_Vec((*after)[0], (*after)[1])}) {
        ray.Normalize();
        if (std::abs(ray.Dot(gp_Vec(normal))) > epsilon)
            return std::nullopt;
    }
    return PlanarLineConstraint{{value((*before)[0]), value((*after)[1])},
                                {normal.X(), normal.Y(), normal.Z()}};
}
} // namespace
std::vector<BoundaryCorner> boundaryCorners(const TopoDS_Face &source, std::uint32_t faceId) {
    if (BRepAdaptor_Surface(source, Standard_True).GetType() != GeomAbs_Plane)
        return {};
    // Forward-face wires put the domain on their left, including clockwise
    // hole loops. Solid orientation does not change convex/reflex trim sense.
    auto face = source;
    face.Orientation(TopAbs_FORWARD);
    std::vector<BoundaryCorner> result;
    for (TopExp_Explorer wires(face, TopAbs_WIRE); wires.More(); wires.Next()) {
        const auto wire = TopoDS::Wire(wires.Current());
        std::vector<TopoDS_Edge> edges;
        for (BRepTools_WireExplorer ordered(wire, face); ordered.More(); ordered.Next())
            edges.push_back(ordered.Current());
        std::size_t uses = 0;
        for (TopExp_Explorer all(wire, TopAbs_EDGE); all.More(); all.Next())
            ++uses;
        if (edges.empty() || edges.size() != uses)
            continue; // optional evidence adds no new CAD import requirement
        for (std::size_t i = 0; i < edges.size(); ++i) {
            const auto &before = edges[(i + edges.size() - 1) % edges.size()];
            const auto &after = edges[i];
            const auto vertex = TopExp::FirstVertex(after, true);
            const auto previous = TopExp::LastVertex(before, true);
            if (vertex.IsNull() || previous.IsNull() || !previous.IsSame(vertex))
                continue;
            const auto p = BRep_Tool::Pnt(vertex);
            finite(p.XYZ());
            if (auto evidence = constraint(before, after, face))
                result.push_back({faceId, value(p), *evidence});
        }
    }
    return result;
}
} // namespace hornsculpt::cad_detail
