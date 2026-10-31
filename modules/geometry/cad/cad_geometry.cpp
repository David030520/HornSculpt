#include "cad_boundary_corners.hpp"
#include "cad_parameterization.hpp"
#include "cad_revolved_projection.hpp"
#include <hornsculpt/cad_geometry.hpp>

#include <Adaptor3d_CurveOnSurface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Extrema_ExtPC.hxx>
#include <Extrema_POnCurv.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GProp_GProps.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2d_Curve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <IMeshTools_Parameters.hxx>
#include <Message_ProgressIndicator.hxx>
#include <Poly_PolygonOnTriangulation.hxx>
#include <Precision.hxx>
#include <STEPControl_Reader.hxx>
#include <ShapeProcess.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <gp_Pnt2d.hxx>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace hornsculpt {
namespace {
using cad_detail::SurfaceMinimum;
void checkStop(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("CAD-Tessellation abgebrochen.");
}
void finite(Vec3 p) {
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
        throw std::invalid_argument("CAD-Koordinaten muessen endlich sein.");
}
Vec3 value(const gp_Pnt &p) {
    return {p.X(), p.Y(), p.Z()};
}
gp_Pnt point(Vec3 p) {
    finite(p);
    return {p.x, p.y, p.z};
}

template <class Function> auto cadOperation(const char *name, Function operation) {
    try {
        return operation();
    } catch (const Standard_Failure &error) {
        throw std::runtime_error(std::string(name) + ": " + error.GetMessageString());
    }
}
std::vector<gp_Pnt> nearestCurvePoints(const TopoDS_Edge &edge, Vec3 approximate) {
    const BRepAdaptor_Curve curve(edge);
    const double first = curve.FirstParameter(), last = curve.LastParameter();
    if (!std::isfinite(first) || !std::isfinite(last) || first >= last)
        throw std::runtime_error("CAD-Randkurve besitzt keinen endlichen getrimmten Bereich.");
    // Query the actual bounded curve. The shape-distance operation can miss
    // interior spline extrema and return only a distant topological endpoint.
    const Extrema_ExtPC projection(point(approximate), curve, first, last, Precision::PConfusion());
    // A line has no stationary point in the trim when its orthogonal
    // projection lies outside; its distance is convex, so an endpoint wins.
    // Other failed searches (including infinitely many circle-axis minima)
    // must not silently become a single seam/endpoint solution.
    if (!projection.IsDone() && curve.GetType() != GeomAbs_Line)
        throw std::runtime_error("CAD-Randprojektion findet keine eindeutigen Kurvenextrema.");
    std::vector<gp_Pnt> result{curve.Value(first), curve.Value(last)};
    const int count = projection.IsDone() ? projection.NbExt() : 0;
    for (int i = 1; i <= count; ++i)
        if (projection.IsMin(i)) {
            const double parameter = projection.Point(i).Parameter();
            if (!std::isfinite(parameter))
                throw std::runtime_error(
                    "CAD-Randprojektion liefert nicht endliche Kurvenparameter.");
            // Analytic roots can fall just outside a rounded native trim.
            // Keep the exact trim: its endpoints are already candidates.
            if (parameter >= first && parameter <= last)
                result.push_back(curve.Value(parameter));
        }
    for (const auto &p : result)
        finite(value(p));
    return result;
}
gp_Pnt uniqueNearest(std::span<const gp_Pnt> candidates, Vec3 approximate) {
    const auto origin = point(approximate);
    gp_Pnt best;
    double minimum = std::numeric_limits<double>::infinity();
    bool found = false, ambiguous = false;
    for (const auto &candidate : candidates) {
        const double d = origin.Distance(candidate);
        if (!std::isfinite(d))
            throw std::runtime_error("CAD-Projektion ueberschreitet den Zahlenbereich.");
        const double roundoff = 128 * std::numeric_limits<double>::epsilon() *
                                std::max({1.0, d, found ? minimum : 0.0});
        if (!found || d < minimum - roundoff) {
            best = candidate;
            minimum = d;
            found = true;
            ambiguous = false;
        } else if (std::abs(d - minimum) <= roundoff) {
            const double positionalRoundoff =
                128 * std::numeric_limits<double>::epsilon() *
                std::max({1.0, std::abs(best.X()), std::abs(best.Y()), std::abs(best.Z()),
                          std::abs(candidate.X()), std::abs(candidate.Y()),
                          std::abs(candidate.Z())});
            ambiguous = ambiguous || best.Distance(candidate) > positionalRoundoff;
        }
    }
    if (!found || ambiguous)
        throw std::invalid_argument("CAD-Projektion ist geometrisch mehrdeutig.");
    return best;
}
Vec3 normalAt(const TopoDS_Face &face, double u, double v) {
    BRepAdaptor_Surface surface(face, Standard_True);
    BRepLProp_SLProps properties(surface, u, v, 2, Precision::PConfusion());
    if (!properties.IsNormalDefined())
        throw std::runtime_error("CAD-Flaechennormale ist an dieser Position nicht definiert.");
    auto normal = properties.Normal();
    if (face.Orientation() == TopAbs_REVERSED)
        normal.Reverse();
    return {normal.X(), normal.Y(), normal.Z()};
}
gp_Pnt2d edgeParameters(const TopoDS_Edge &edge, const TopoDS_Face &face, const gp_Pnt &pointOnEdge,
                        double parameter, double tolerance) {
    double first, last;
    const auto curve = BRep_Tool::CurveOnSurface(edge, face, first, last);
    if (curve.IsNull() || !std::isfinite(parameter) || parameter < first || parameter > last)
        throw std::runtime_error("CAD-Randprojektion besitzt keine gueltige Flaechenkurve.");
    const Handle(BRepAdaptor_Surface) surface = new BRepAdaptor_Surface(face, Standard_True);
    const auto initial = curve->Value(parameter);
    const double initialError = surface->Value(initial.X(), initial.Y()).Distance(pointOnEdge);
    const double roundoff = 128 * std::numeric_limits<double>::epsilon() *
                            std::max({1.0, std::abs(pointOnEdge.X()), std::abs(pointOnEdge.Y()),
                                      std::abs(pointOnEdge.Z())});
    if (initialError <= tolerance + roundoff)
        return initial;
    // SameParameter is an approximation, not an exact identity. Invert only
    // this topological pcurve; never search another underlying-surface branch.
    const Handle(Geom2dAdaptor_Curve) curve2d = new Geom2dAdaptor_Curve(curve, first, last);
    const Adaptor3d_CurveOnSurface embedded(curve2d, surface);
    const Extrema_ExtPC projection(pointOnEdge, embedded, first, last, roundoff);
    std::vector<double> parameters;
    std::vector<gp_Pnt> points;
    const auto append = [&](double t) {
        const auto p = embedded.Value(t);
        finite(value(p));
        // A parameter correction must stay local along the curve, including
        // loops whose remote parameters happen to approach the same 3D point.
        const double arc = t == parameter
                               ? 0
                               : GCPnts_AbscissaPoint::Length(embedded, std::min(t, parameter),
                                                              std::max(t, parameter), roundoff);
        if (!std::isfinite(arc) || arc > initialError + tolerance + roundoff)
            return;
        points.push_back(p);
        parameters.push_back(t);
    };
    append(first);
    append(last);
    if (projection.IsDone())
        for (int i = 1; i <= projection.NbExt(); ++i)
            if (projection.IsMin(i))
                append(projection.Point(i).Parameter());
    if (points.empty())
        throw std::runtime_error("CAD-Randkurve und Flaeche sind lokal nicht modellkonsistent.");
    const auto nearest = uniqueNearest(points, value(pointOnEdge));
    for (std::size_t i = 0; i < points.size(); ++i)
        if (points[i].Distance(nearest) <= roundoff)
            return curve->Value(parameters[i]);
    throw std::runtime_error("CAD-Randkurvenparameter sind geometrisch mehrdeutig.");
}
gp_Pnt2d supportedParameters(const BRepExtrema_DistShapeShape &projection, int index,
                             const TopoDS_Face &face) {
    gp_Pnt2d uv;
    const auto support = projection.SupportOnShape2(index);
    double tolerance = BRep_Tool::Tolerance(face);
    switch (projection.SupportTypeShape2(index)) {
    case BRepExtrema_IsInFace: {
        double u, v;
        projection.ParOnFaceS2(index, u, v);
        uv.SetCoord(u, v);
        break;
    }
    case BRepExtrema_IsOnEdge: {
        const auto edge = TopoDS::Edge(support);
        TopTools_IndexedMapOfShape edges;
        TopExp::MapShapes(face, TopAbs_EDGE, edges);
        if (!edges.Contains(edge))
            throw std::runtime_error("CAD-Projektionskante gehoert nicht zur getrimmten Flaeche.");
        double parameter;
        projection.ParOnEdgeS2(index, parameter);
        tolerance = std::max(tolerance, BRep_Tool::Tolerance(edge));
        uv = edgeParameters(edge, face, projection.PointOnShape2(index), parameter, tolerance);
        break;
    }
    case BRepExtrema_IsVertex: {
        const auto vertex = TopoDS::Vertex(support);
        TopTools_IndexedMapOfShape vertices;
        TopExp::MapShapes(face, TopAbs_VERTEX, vertices);
        if (!vertices.Contains(vertex))
            throw std::runtime_error("CAD-Projektionsknoten gehoert nicht zur getrimmten Flaeche.");
        uv = BRep_Tool::Parameters(vertex, face);
        tolerance = std::max(tolerance, BRep_Tool::Tolerance(vertex));
        break;
    }
    default:
        throw std::runtime_error("CAD-Projektion besitzt keinen unterstuetzten Flaechentraeger.");
    }
    if (!std::isfinite(uv.X()) || !std::isfinite(uv.Y()) || !std::isfinite(tolerance) ||
        tolerance <= 0)
        throw std::runtime_error("CAD-Projektion besitzt ungueltige Flaechenparameter.");
    // A native edge/vertex support and its face-specific pcurve are the trim
    // proof at a boundary: the 2D classifier can report OUT even for that exact
    // pcurve. Interior support still requires an IN/ON parameter classification.
    if (projection.SupportTypeShape2(index) == BRepExtrema_IsInFace) {
        const BRepTopAdaptor_FClass2d classifier(face, Precision::PConfusion());
        const auto state = classifier.Perform(uv);
        if (state != TopAbs_IN && state != TopAbs_ON)
            throw std::runtime_error(
                "CAD-Projektionsparameter liegen ausserhalb der getrimmten Flaeche.");
    }
    const auto p = projection.PointOnShape2(index);
    const double roundoff = 128 * std::numeric_limits<double>::epsilon() *
                            std::max({1.0, std::abs(p.X()), std::abs(p.Y()), std::abs(p.Z())});
    const BRepAdaptor_Surface surface(face, Standard_True);
    // Edge pcurves and 3D curves are separate approximations. The returned
    // face projection uses the pcurve's actual surface point; their residual
    // is not a condition that this face point must equal the 3D edge point.
    if (projection.SupportTypeShape2(index) != BRepExtrema_IsOnEdge &&
        surface.Value(uv.X(), uv.Y()).Distance(p) > tolerance + roundoff)
        throw std::runtime_error("CAD-Projektionsparameter rekonstruieren den Randpunkt nicht.");
    return uv;
}
SurfaceProjection trimmedProjection(const TopoDS_Face &face, Vec3 approximate) {
    const auto vertex = BRepBuilderAPI_MakeVertex(point(approximate)).Vertex();
    const BRepExtrema_DistShapeShape projection(vertex, face);
    if (!projection.IsDone() || projection.NbSolution() == 0)
        throw std::runtime_error("CAD-Projektion findet keinen Punkt auf der getrimmten Flaeche.");
    std::vector<gp_Pnt> candidates;
    for (int i = 1; i <= projection.NbSolution(); ++i) {
        const auto p = projection.PointOnShape2(i);
        finite(value(p));
        candidates.push_back(p);
    }
    const auto nearest = uniqueNearest(candidates, approximate);
    const double roundoff =
        128 * std::numeric_limits<double>::epsilon() *
        std::max({1.0, std::abs(nearest.X()), std::abs(nearest.Y()), std::abs(nearest.Z())});
    for (int i = 1; i <= projection.NbSolution(); ++i)
        if (projection.PointOnShape2(i).Distance(nearest) <= roundoff) {
            // Keep the parameters of the actual trimmed-face solution. A new
            // unrestricted inversion can select a remote or OUT branch of a
            // spline and return an unrelated normal for this same boundary point.
            const auto uv = supportedParameters(projection, i, face);
            const BRepAdaptor_Surface surface(face, Standard_True);
            const auto onFace = surface.Value(uv.X(), uv.Y());
            finite(value(onFace));
            return {value(onFace), normalAt(face, uv.X(), uv.Y())};
        }
    throw std::runtime_error("CAD-Projektion verlor ihren eindeutigen Flaechentraeger.");
}
CadFaceMeasures measures(const TopoDS_Face &face, std::uint32_t id) {
    GProp_GProps properties;
    BRepGProp::SurfaceProperties(face, properties);
    CadFaceMeasures result;
    result.faceId = id;
    result.area = properties.Mass();
    result.center = value(properties.CentreOfMass());
    if (!std::isfinite(result.area) || result.area <= 0)
        throw std::invalid_argument("CAD-Flaeche hat keinen positiven Flaecheninhalt.");
    finite(result.center);
    double u0, u1, v0, v1;
    BRepTools::UVBounds(face, u0, u1, v0, v1);
    result.normal = normalAt(face, std::midpoint(u0, u1), std::midpoint(v0, v1));
    result.surfaceType = BRep_Tool::Surface(face)->DynamicType()->Name();
    if (result.surfaceType.starts_with("Geom_"))
        result.surfaceType.erase(0, 5);
    switch (BRepAdaptor_Surface(face).GetType()) {
    case GeomAbs_Cylinder:
        result.surfaceType = "Cylinder";
        break;
    case GeomAbs_Cone:
        result.surfaceType = "Cone";
        break;
    case GeomAbs_Sphere:
        result.surfaceType = "Sphere";
        break;
    case GeomAbs_Torus:
        result.surfaceType = "Torus";
        break;
    default:
        break;
    }
    TopTools_IndexedMapOfShape edges, wires;
    TopExp::MapShapes(face, TopAbs_EDGE, edges);
    TopExp::MapShapes(face, TopAbs_WIRE, wires);
    result.edgeCount = static_cast<std::uint32_t>(edges.Extent());
    result.wireCount = static_cast<std::uint32_t>(wires.Extent());
    return result;
}
void validateShape(const TopoDS_Shape &shape) {
    if (shape.IsNull() || !BRepCheck_Analyzer(shape, Standard_True).IsValid())
        throw std::invalid_argument(
            "CAD-Geometrie ist leer oder topologisch/geometrisch ungueltig.");
    TopTools_IndexedMapOfShape solids, faces;
    TopExp::MapShapes(shape, TopAbs_SOLID, solids);
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    if (solids.IsEmpty() || faces.IsEmpty())
        throw std::invalid_argument("CAD-Import verlangt geschlossene Volumenkoerper.");
    TopTools_MapOfShape covered;
    for (int i = 1; i <= solids.Extent(); ++i) {
        const auto solid = TopoDS::Solid(solids(i));
        GProp_GProps volume;
        BRepGProp::VolumeProperties(solid, volume, Standard_True);
        BRepClass3d_SolidClassifier classifier(solid);
        classifier.PerformInfinitePoint(Precision::Confusion());
        if (!std::isfinite(volume.Mass()) || volume.Mass() <= 0 || classifier.State() != TopAbs_OUT)
            throw std::invalid_argument(
                "CAD-Koerper ist offen, degeneriert oder nach innen orientiert.");
        for (TopExp_Explorer it(solid, TopAbs_FACE); it.More(); it.Next())
            covered.Add(it.Current());
    }
    if (covered.Extent() != faces.Extent())
        throw std::invalid_argument(
            "CAD-Geometrie enthaelt freie Flaechen ausserhalb ihrer Volumenkoerper.");
    TopTools_IndexedMapOfShape edges;
    TopExp::MapShapes(shape, TopAbs_EDGE, edges);
    std::vector<std::pair<unsigned, int>> uses(static_cast<std::size_t>(edges.Extent()));
    for (int f = 1; f <= faces.Extent(); ++f)
        for (TopExp_Explorer it(faces(f), TopAbs_EDGE); it.More(); it.Next()) {
            const auto edge = TopoDS::Edge(it.Current());
            if (BRep_Tool::Degenerated(edge))
                continue;
            auto &use = uses[static_cast<std::size_t>(edges.FindIndex(edge) - 1)];
            ++use.first;
            if (edge.Orientation() == TopAbs_FORWARD)
                ++use.second;
            else if (edge.Orientation() == TopAbs_REVERSED)
                --use.second;
            else
                throw std::invalid_argument(
                    "CAD-Randkante hat keine orientierte aeussere Verwendung.");
        }
    for (int i = 1; i <= edges.Extent(); ++i)
        if (!BRep_Tool::Degenerated(TopoDS::Edge(edges(i))) &&
            uses[static_cast<std::size_t>(i - 1)] != std::pair<unsigned, int>{2, 0})
            throw std::invalid_argument(
                "CAD-Rand ist keine geschlossene orientierte Zweimannigfaltigkeit.");
}
std::string writeBrep(const TopoDS_Shape &shape) {
    std::ostringstream output;
    BRepTools::Write(shape, output, Standard_False, Standard_False, TopTools_FormatVersion_CURRENT);
    if (!output)
        throw std::runtime_error("CAD-BRep kann nicht gespeichert werden.");
    return output.str();
}
class StopProgress final : public Message_ProgressIndicator {
  public:
    explicit StopProgress(std::stop_token stop) : stop_(stop) {}

  protected:
    Standard_Boolean UserBreak() override {
        return stop_.stop_requested();
    }
    void Show(const Message_ProgressScope &, Standard_Boolean) noexcept override {
        // OCCT also calls this presentation hook while destroying progress
        // scopes. Cancellation belongs to UserBreak; throw only after meshing.
    }

  private:
    std::stop_token stop_;
};
struct PreviewNode {
    Vec3 position;
    std::uint32_t parent;
    unsigned priority{};
};
class UnresolvedTessellation final : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};
class PreviewBuilder {
  public:
    PreviewBuilder(double tolerance, std::stop_token stop)
        : stop_(stop), parameterization_(tolerance) {}
    void append(const TopoDS_Shape &source, std::span<const std::uint32_t> ids, double tolerance) {
        checkStop(stop_);
        BRepBuilderAPI_Copy copy(source, Standard_False, Standard_False);
        const auto shape = copy.Shape();
        IMeshTools_Parameters parameters;
        parameters.Deflection = parameters.DeflectionInterior = tolerance;
        parameters.Relative = Standard_False;
        parameters.InParallel = Standard_False;
        parameters.ControlSurfaceDeflection = Standard_True;
        parameters.EnableControlSurfaceDeflectionAllSurfaces = Standard_True;
        Handle(Message_ProgressIndicator) progress = new StopProgress(stop_);
        BRepMesh_IncrementalMesh mesher(shape, parameters, progress->Start());
        checkStop(stop_);
        if (!mesher.IsDone())
            throw std::runtime_error("OpenCASCADE konnte die CAD-Oberflaeche nicht tessellieren.");
        TopTools_IndexedMapOfShape sourceFaces, vertices, edges;
        TopExp::MapShapes(source, TopAbs_FACE, sourceFaces);
        TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
        TopExp::MapShapes(shape, TopAbs_EDGE, edges);
        std::vector<std::uint32_t> vertexNodes(static_cast<std::size_t>(vertices.Extent()), absent);
        std::vector<std::map<double, std::uint32_t>> edgeNodes(
            static_cast<std::size_t>(edges.Extent()));
        std::vector<std::uint32_t> curveIds(static_cast<std::size_t>(edges.Extent()), absent);
        for (int i = 1; i <= edges.Extent(); ++i)
            if (!BRep_Tool::Degenerated(TopoDS::Edge(edges(i))))
                curveIds[static_cast<std::size_t>(i - 1)] =
                    parameterization_.addCurve(TopoDS::Edge(edges(i)));
        for (int i = 1; i <= sourceFaces.Extent(); ++i) {
            checkStop(stop_);
            auto face = TopoDS::Face(copy.ModifiedShape(sourceFaces(i)));
            face.Orientation(sourceFaces(i).Orientation());
            parameterization_.addFace(ids[static_cast<std::size_t>(i - 1)], face);
            appendFace(face, ids[static_cast<std::size_t>(i - 1)], vertices, edges, vertexNodes,
                       edgeNodes, curveIds);
        }
    }
    MeshingSource finish() {
        SurfaceMesh mesh;
        std::vector<ParameterTriangle> records;
        std::map<std::uint32_t, std::uint32_t> indices;
        for (std::size_t index = 0; index < triangles_.size(); ++index) {
            auto triangle = triangles_[index];
            checkStop(stop_);
            for (auto &v : triangle.vertices) {
                const auto root = find(v);
                auto [it, created] =
                    indices.emplace(root, static_cast<std::uint32_t>(mesh.vertices.size()));
                if (created)
                    mesh.vertices.push_back(nodes_[root].position);
                v = it->second;
            }
            if (triangle.vertices[0] == triangle.vertices[1] ||
                triangle.vertices[1] == triangle.vertices[2] ||
                triangle.vertices[2] == triangle.vertices[0])
                continue;
            const auto a = point(mesh.vertices[triangle.vertices[0]]),
                       b = point(mesh.vertices[triangle.vertices[1]]),
                       c = point(mesh.vertices[triangle.vertices[2]]);
            const double area = gp_Vec(a, b).Crossed(gp_Vec(a, c)).Magnitude();
            if (!std::isfinite(area) || area <= 0)
                throw std::runtime_error("CAD-Tessellation enthaelt entartete Dreiecke.");
            mesh.triangles.push_back(triangle);
            records.push_back(parameterTriangles_[index]);
        }
        std::map<std::array<std::uint32_t, 2>, std::pair<unsigned, int>> uses;
        for (const auto &t : mesh.triangles)
            for (unsigned i = 0; i < 3; ++i) {
                const auto a = t.vertices[i], b = t.vertices[(i + 1) % 3];
                auto &use = uses[{std::min(a, b), std::max(a, b)}];
                ++use.first;
                use.second += a < b ? 1 : -1;
            }
        for (const auto &[edge, use] : uses) {
            if (use != std::pair<unsigned, int>{2, 0})
                throw UnresolvedTessellation(
                    "CAD-Tessellation ist an gemeinsamen topologischen Kanten "
                    "nicht wasserdicht orientiert: " +
                    std::to_string(edge[0]) + "-" + std::to_string(edge[1]) +
                    ", Verwendungen=" + std::to_string(use.first) +
                    ", Orientierung=" + std::to_string(use.second) + ".");
        }
        if (mesh.triangles.empty())
            throw std::runtime_error("CAD-Tessellation ist leer.");
        auto parameters = parameterization_.finish(mesh, std::move(records));
        return {std::move(mesh), std::move(parameters)};
    }

  private:
    static constexpr auto absent = std::numeric_limits<std::uint32_t>::max();
    static constexpr std::size_t maximumNodes = 2000000;
    std::stop_token stop_;
    std::vector<PreviewNode> nodes_;
    std::vector<Triangle> triangles_;
    std::vector<ParameterTriangle> parameterTriangles_;
    cad_detail::ParameterizationBuilder parameterization_;
    std::uint32_t find(std::uint32_t n) {
        while (nodes_[n].parent != n) {
            nodes_[n].parent = nodes_[nodes_[n].parent].parent;
            n = nodes_[n].parent;
        }
        return n;
    }
    void unite(std::uint32_t a, std::uint32_t b) {
        a = find(a);
        b = find(b);
        if (a == b)
            return;
        if (nodes_[a].priority < nodes_[b].priority)
            std::swap(a, b);
        nodes_[b].parent = a;
    }
    void canonical(std::uint32_t n, Vec3 p, unsigned priority) {
        n = find(n);
        if (priority > nodes_[n].priority) {
            nodes_[n].position = p;
            nodes_[n].priority = priority;
        }
    }
    void appendFace(const TopoDS_Face &face, std::uint32_t id,
                    const TopTools_IndexedMapOfShape &vertices,
                    const TopTools_IndexedMapOfShape &edges,
                    std::vector<std::uint32_t> &vertexNodes,
                    std::vector<std::map<double, std::uint32_t>> &edgeNodes,
                    std::span<const std::uint32_t> curveIds) {
        TopLoc_Location location;
        const auto triangulation = BRep_Tool::Triangulation(face, location);
        if (triangulation.IsNull() || triangulation->NbTriangles() == 0)
            throw UnresolvedTessellation("CAD-Flaeche besitzt keine Tessellation.");
        if (!triangulation->HasUVNodes())
            throw std::runtime_error("CAD-Tessellation besitzt keine nativen Eckparameter.");
        std::map<std::array<int, 2>, CurveUseInterval> edgeUses;
        if (nodes_.size() + static_cast<std::size_t>(triangulation->NbNodes()) > maximumNodes ||
            triangles_.size() + static_cast<std::size_t>(triangulation->NbTriangles()) >
                2 * maximumNodes)
            throw std::length_error("CAD-Tessellation ueberschreitet das Geometriebudget.");
        const auto first = static_cast<std::uint32_t>(nodes_.size());
        for (int i = 1; i <= triangulation->NbNodes(); ++i) {
            if (i % 256 == 0)
                checkStop(stop_);
            auto p = value(triangulation->Node(i).Transformed(location.Transformation()));
            finite(p);
            auto index = static_cast<std::uint32_t>(nodes_.size());
            nodes_.push_back({p, index, 0});
        }
        for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next()) {
            checkStop(stop_);
            const auto edge = TopoDS::Edge(it.Current());
            const auto polygon = BRep_Tool::PolygonOnTriangulation(edge, triangulation, location);
            if (polygon.IsNull())
                throw std::runtime_error(
                    "CAD-Randkante besitzt kein gemeinsames Tessellationspolygon.");
            const auto edgeIndex = static_cast<std::size_t>(edges.FindIndex(edge) - 1);
            TopoDS_Vertex v0, v1;
            TopExp::Vertices(edge, v0, v1);
            const bool degenerate = BRep_Tool::Degenerated(edge);
            std::optional<BRepAdaptor_Curve> curve;
            if (!degenerate) {
                curve.emplace(edge);
                if (!polygon->HasParameters())
                    throw std::runtime_error("CAD-Randpolygon besitzt keine Kurvenparameter.");
                const auto curveId = curveIds[edgeIndex];
                const int anchor = (polygon->NbNodes() + 1) / 2;
                const auto uv = triangulation->UVNode(polygon->Node(anchor));
                const auto useId = parameterization_.addUse(
                    curveId, edge, face, polygon->Parameter(anchor), {uv.X(), uv.Y()});
                for (int j = 1; j < polygon->NbNodes(); ++j) {
                    const int a = polygon->Node(j), b = polygon->Node(j + 1);
                    if (a == b)
                        continue;
                    const double t0 = polygon->Parameter(j), t1 = polygon->Parameter(j + 1);
                    finite({t0, t1, 0});
                    if (!edgeUses
                             .emplace(std::array{a, b}, CurveUseInterval{curveId, useId, t0, t1})
                             .second ||
                        !edgeUses
                             .emplace(std::array{b, a}, CurveUseInterval{curveId, useId, t1, t0})
                             .second)
                        throw std::runtime_error(
                            "CAD-Quellkante hat mehrdeutige Randverwendungen.");
                }
            }
            for (int j = 1; j <= polygon->NbNodes(); ++j) {
                const auto n = first + static_cast<std::uint32_t>(polygon->Node(j) - 1);
                if (degenerate || j == 1 || j == polygon->NbNodes()) {
                    if (v0.IsNull() || v1.IsNull())
                        throw std::runtime_error(
                            "CAD-Randkante besitzt unvollstaendige Eckpunkte.");
                    const auto p0 = BRep_Tool::Pnt(v0), p1 = BRep_Tool::Pnt(v1);
                    const auto vertex = point(nodes_[n].position).Distance(p0) <=
                                                point(nodes_[n].position).Distance(p1)
                                            ? v0
                                            : v1;
                    const auto vi = static_cast<std::size_t>(vertices.FindIndex(vertex) - 1);
                    canonical(n, value(BRep_Tool::Pnt(vertex)), 2);
                    if (vertexNodes[vi] == absent)
                        vertexNodes[vi] = n;
                    else
                        unite(n, vertexNodes[vi]);
                }
                if (degenerate)
                    continue;
                const double parameter = polygon->Parameter(j),
                             range = curve->LastParameter() - curve->FirstParameter();
                if (!std::isfinite(parameter) || !std::isfinite(range) || range <= 0)
                    throw std::runtime_error(
                        "CAD-Randkurve besitzt keine endliche Parametrisierung.");
                const double normalized = (parameter - curve->FirstParameter()) / range;
                auto &registry = edgeNodes[edgeIndex];
                auto found = registry.lower_bound(normalized);
                constexpr double parameterRoundoff = 256 * std::numeric_limits<double>::epsilon();
                if (found == registry.end() ||
                    std::abs(found->first - normalized) > parameterRoundoff) {
                    if (found != registry.begin() &&
                        std::abs(std::prev(found)->first - normalized) <= parameterRoundoff)
                        found = std::prev(found);
                    else
                        found = registry.emplace(normalized, n).first;
                }
                canonical(n, value(curve->Value(parameter)), 1);
                unite(n, found->second);
            }
        }
        for (int i = 1; i <= triangulation->NbTriangles(); ++i) {
            if (i % 256 == 0)
                checkStop(stop_);
            int a, b, c;
            triangulation->Triangle(i).Get(a, b, c);
            if (face.Orientation() == TopAbs_REVERSED)
                std::swap(b, c);
            ParameterTriangle record;
            const std::array local{a, b, c};
            for (unsigned j = 0; j < 3; ++j) {
                const auto uv = triangulation->UVNode(local[j]);
                finite({uv.X(), uv.Y(), 0});
                record.uv[j] = {uv.X(), uv.Y()};
                const auto found = edgeUses.find({local[j], local[(j + 1) % 3]});
                if (found != edgeUses.end())
                    record.edges[j] = found->second;
            }
            triangles_.push_back({{first + static_cast<std::uint32_t>(a - 1),
                                   first + static_cast<std::uint32_t>(b - 1),
                                   first + static_cast<std::uint32_t>(c - 1)},
                                  id,
                                  BoundaryCondition::wall});
            parameterTriangles_.push_back(record);
        }
    }
};
} // namespace

namespace {
// Each lease owns mutable OCCT search state exclusively. The initialized
// surface grid survives successive points without becoming shared solver state.
class SurfaceProjectorPool {
    struct Entry {
        GeomAPI_ProjectPointOnSurf projector;
        bool inUse{};
    };

  public:
    SurfaceProjectorPool(Handle(Geom_Surface) surface, std::array<double, 4> bounds)
        : surface_(std::move(surface)), bounds_(bounds),
          revolved_(cad_detail::RevolvedProjection::make(surface_, bounds_)) {}

    const std::optional<cad_detail::RevolvedProjection> &revolved() const noexcept {
        return revolved_;
    }

    class Lease {
      public:
        Lease(SurfaceProjectorPool &pool, Entry &entry) : pool_(pool), entry_(entry) {}
        Lease(const Lease &) = delete;
        Lease &operator=(const Lease &) = delete;
        ~Lease() {
            std::lock_guard lock(pool_.mutex_);
            entry_.inUse = false;
        }
        GeomAPI_ProjectPointOnSurf *operator->() {
            return &entry_.projector;
        }

      private:
        SurfaceProjectorPool &pool_;
        Entry &entry_;
    };

    Lease acquire() {
        std::lock_guard lock(mutex_);
        for (auto &entry : entries_)
            if (!entry->inUse) {
                entry->inUse = true;
                return Lease(*this, *entry);
            }
        auto entry = std::make_unique<Entry>();
        entry->projector.Init(surface_, bounds_[0], bounds_[1], bounds_[2], bounds_[3],
                              Precision::PConfusion());
        entry->projector.SetExtremaFlag(Extrema_ExtFlag_MIN);
        entry->inUse = true;
        entries_.push_back(std::move(entry));
        return Lease(*this, *entries_.back());
    }

  private:
    Handle(Geom_Surface) surface_;
    std::array<double, 4> bounds_;
    std::optional<cad_detail::RevolvedProjection> revolved_;
    std::mutex mutex_;
    // Allocation grows only with concurrent leases on this face. Stable entries
    // make release allocation-free even while another caller expands the pool.
    std::vector<std::unique_ptr<Entry>> entries_;
};
std::optional<SurfaceMinimum> surfaceMinimum(SurfaceProjectorPool &pool,
                                             const std::array<double, 4> &bounds, Vec3 approximate,
                                             bool allowParameterRoundoff = false) {
    if (pool.revolved())
        return pool.revolved()->minimum(approximate);
    auto projection = pool.acquire();
    projection->Perform(point(approximate));
    if (projection->IsDone() && projection->NbPoints() > 0) {
        const auto p = projection->NearestPoint();
        const double minimum = projection->LowerDistance();
        const double distanceRoundoff =
            128 * std::numeric_limits<double>::epsilon() * std::max(1.0, minimum);
        const double pointRoundoff =
            128 * std::numeric_limits<double>::epsilon() *
            std::max({1.0, std::abs(p.X()), std::abs(p.Y()), std::abs(p.Z())});
        bool unique = true;
        for (int i = 1; i <= projection->NbPoints(); ++i) {
            finite(value(projection->Point(i)));
            const double distance = projection->Distance(i);
            if (!std::isfinite(distance))
                throw std::runtime_error(
                    "CAD-Flaechenprojektion liefert nicht endliche Abstaende.");
            if (std::abs(distance - minimum) <= distanceRoundoff &&
                projection->Point(i).Distance(p) > pointRoundoff)
                unique = false;
        }
        double u, v;
        projection->LowerDistanceParameters(u, v);
        if (!std::isfinite(minimum) || !std::isfinite(u) || !std::isfinite(v))
            throw std::runtime_error("CAD-Flaechenprojektion liefert nicht endliche Parameter.");
        const double parameterRoundoff =
            allowParameterRoundoff
                ? 128 * std::numeric_limits<double>::epsilon() *
                      std::max({1.0, std::abs(u), std::abs(v), std::abs(bounds[0]),
                                std::abs(bounds[1]), std::abs(bounds[2]), std::abs(bounds[3])})
                : 0;
        if (unique && std::isfinite(minimum) && std::isfinite(u) && std::isfinite(v) &&
            u >= bounds[0] - parameterRoundoff && u <= bounds[1] + parameterRoundoff &&
            v >= bounds[2] - parameterRoundoff && v <= bounds[3] + parameterRoundoff) {
            finite(value(p));
            return SurfaceMinimum{p, u, v};
        }
    }
    return std::nullopt;
}
} // namespace

struct CadGeometry::Data {
    struct Part {
        TopoDS_Shape shape;
        std::vector<std::uint32_t> ids;
    };
    struct Face {
        TopoDS_Face shape;
        TopTools_IndexedMapOfShape edges;
        Handle(Geom_Surface) surface;
        std::array<double, 4> bounds{};
        std::unique_ptr<SurfaceProjectorPool> projectors;
        std::unique_ptr<const BRepTopAdaptor_FClass2d> trimClassifier;
        Handle(Geom_Surface) supportingSurface;
        std::array<double, 4> supportingBounds{};
        std::unique_ptr<SurfaceProjectorPool> supportingProjectors;
    };
    std::vector<Part> parts;
    std::vector<CadBrepPart> artifacts;
    std::vector<CadFaceMeasures> measures;
    std::vector<BoundaryCorner> corners;
    std::map<std::uint32_t, Face> faces;
    void append(TopoDS_Shape shape, std::vector<std::uint32_t> ids) {
        validateShape(shape);
        TopTools_IndexedMapOfShape mapped;
        TopExp::MapShapes(shape, TopAbs_FACE, mapped);
        if (ids.size() != static_cast<std::size_t>(mapped.Extent()))
            throw std::invalid_argument("CAD-Flaechen-IDs passen nicht zum exakten BRep-Artefakt.");
        for (int i = 1; i <= mapped.Extent(); ++i) {
            const auto id = ids[static_cast<std::size_t>(i - 1)];
            if (id == 0 || id == std::numeric_limits<std::uint32_t>::max() || faces.contains(id))
                throw std::invalid_argument(
                    "CAD-Flaechen-IDs muessen eindeutig, positiv und fortsetzbar sein.");
            auto face = TopoDS::Face(mapped(i));
            Face record;
            record.shape = face;
            // The ray classifier can reject interior points of small notched
            // planar trims after SI scaling. Use OCCT's wire classifier as the
            // primary trim strategy; OUT remains OUT, including real holes.
            record.trimClassifier =
                std::make_unique<const BRepTopAdaptor_FClass2d>(face, Precision::PConfusion());
            record.surface = BRep_Tool::Surface(face);
            BRepTools::UVBounds(face, record.bounds[0], record.bounds[1], record.bounds[2],
                                record.bounds[3]);
            record.projectors =
                std::make_unique<SurfaceProjectorPool>(record.surface, record.bounds);
            record.supportingSurface = record.surface;
            while (const auto trimmed =
                       Handle(Geom_RectangularTrimmedSurface)::DownCast(record.supportingSurface))
                record.supportingSurface = trimmed->BasisSurface();
            record.supportingSurface->Bounds(record.supportingBounds[0], record.supportingBounds[1],
                                             record.supportingBounds[2],
                                             record.supportingBounds[3]);
            record.supportingProjectors = std::make_unique<SurfaceProjectorPool>(
                record.supportingSurface, record.supportingBounds);
            TopExp::MapShapes(face, TopAbs_EDGE, record.edges);
            faces.emplace(id, std::move(record));
            measures.push_back(::hornsculpt::measures(face, id));
            const auto evidence = cad_detail::boundaryCorners(face, id);
            corners.insert(corners.end(), evidence.begin(), evidence.end());
        }
        artifacts.push_back({writeBrep(shape), ids});
        parts.push_back({std::move(shape), std::move(ids)});
    }
    const Face &face(std::uint32_t id) const {
        auto found = faces.find(id);
        if (found == faces.end())
            throw std::invalid_argument("Unbekannte exakte CAD-Flaechen-ID.");
        return found->second;
    }
};
CadGeometry::CadGeometry(std::unique_ptr<Data> data) : data_(std::move(data)) {}
CadGeometry::~CadGeometry() = default;
std::shared_ptr<const CadGeometry> CadGeometry::fromBrepParts(std::span<const CadBrepPart> parts,
                                                              double metresPerUnit) {
    return cadOperation("CAD-BRep-Import", [&] {
        if (parts.empty() || !std::isfinite(metresPerUnit) || metresPerUnit <= 0)
            throw std::invalid_argument("CAD-BRep-Einheiten oder Teile sind ungueltig.");
        auto data = std::make_unique<Data>();
        std::size_t bytes = 0;
        for (const auto &part : parts) {
            if (part.data.empty() || part.data.size() > 256ull * 1024 * 1024 - bytes)
                throw std::length_error("CAD-BRep-Artefakte sind leer oder zu gross.");
            bytes += part.data.size();
            std::istringstream stream(part.data);
            TopoDS_Shape shape;
            BRep_Builder builder;
            BRepTools::Read(shape, stream, builder);
            if (shape.IsNull() || stream.bad())
                throw std::invalid_argument("CAD-BRep-Artefakt ist unlesbar.");
            if (metresPerUnit != 1) {
                gp_Trsf scale;
                scale.SetScale(gp_Pnt(0, 0, 0), metresPerUnit);
                shape = BRepBuilderAPI_Transform(shape, scale, Standard_True).Shape();
            }
            data->append(shape, part.faceIds);
        }
        return std::shared_ptr<const CadGeometry>(new CadGeometry(std::move(data)));
    });
}
std::shared_ptr<const CadGeometry> CadGeometry::fromStep(const std::string &path) {
    return cadOperation("STEP-Import", [&] {
        if (path.empty())
            throw std::invalid_argument("STEP-Dateipfad ist leer.");
        STEPControl_Reader reader;
        reader.SetShapeProcessFlags(ShapeProcess::OperationsFlags{});
        if (reader.ReadFile(path.c_str()) != IFSelect_RetDone)
            throw std::invalid_argument("STEP-Datei kann nicht gelesen werden.");
        // OCCT specifies transfer length units in millimetres: 1000 mm = 1 metre.
        reader.SetSystemLengthUnit(1000.0);
        if (reader.TransferRoots() == 0)
            throw std::invalid_argument("STEP-Datei liefert keine CAD-Koerper.");
        auto shape = reader.OneShape();
        TopTools_IndexedMapOfShape faces;
        TopExp::MapShapes(shape, TopAbs_FACE, faces);
        std::vector<std::uint32_t> ids(static_cast<std::size_t>(faces.Extent()));
        std::iota(ids.begin(), ids.end(), 1u);
        auto data = std::make_unique<Data>();
        data->append(shape, std::move(ids));
        return std::shared_ptr<const CadGeometry>(new CadGeometry(std::move(data)));
    });
}
std::span<const CadFaceMeasures> CadGeometry::faces() const noexcept {
    return data_->measures;
}
std::span<const CadBrepPart> CadGeometry::brepParts() const noexcept {
    return data_->artifacts;
}
std::span<const BoundaryCorner> CadGeometry::boundaryCorners() const noexcept {
    return data_->corners;
}
SurfaceMesh CadGeometry::previewMesh(double tolerance, std::stop_token stop) const {
    return meshingSource(tolerance, stop).mesh;
}
MeshingSource CadGeometry::meshingSource(double tolerance, std::stop_token stop) const {
    return cadOperation("CAD-Tessellation", [&] {
        checkStop(stop);
        if (!std::isfinite(tolerance) || tolerance < Precision::Confusion())
            throw std::invalid_argument(
                "CAD-Detailabweichung liegt ausserhalb der endlichen Modellpraezision.");
        // A coarse chord approximation can cross another trim on a thin curved
        // face. OCCT may still report IsDone with an incomplete triangulation.
        // Retry only this discrete boundary failure, with a strictly smaller
        // chord budget. The CAD model and requested accuracy remain unchanged.
        for (double deflection = tolerance;; deflection *= .5) {
            checkStop(stop);
            try {
                PreviewBuilder builder(tolerance, stop);
                for (const auto &part : data_->parts)
                    builder.append(part.shape, part.ids, deflection);
                return builder.finish();
            } catch (const UnresolvedTessellation &) {
                if (deflection * .5 < Precision::Confusion())
                    throw;
            }
        }
    });
}
SurfaceProjection CadGeometry::projectSurface(std::uint32_t id, Vec3 approximate) const {
    return cadOperation("CAD-Flaechenprojektion", [&] {
        const auto &face = data_->face(id);
        if (const auto minimum = surfaceMinimum(*face.projectors, face.bounds, approximate)) {
            const auto state = face.trimClassifier->Perform(gp_Pnt2d(minimum->u, minimum->v));
            // Interior stationary points are admissible trim candidates.
            if (state == TopAbs_IN || state == TopAbs_ON) {
                const SurfaceProjection interior{value(minimum->position),
                                                 normalAt(face.shape, minimum->u, minimum->v)};
                // A stationary surface minimum can be farther than a trim
                // edge when the unconstrained nearest point lies outside UV.
                // Only a zero-distance witness already proves the minimum.
                const double roundoff =
                    128 * std::numeric_limits<double>::epsilon() *
                    std::max({1.0, std::abs(approximate.x), std::abs(approximate.y),
                              std::abs(approximate.z)});
                if (minimum->position.Distance(point(approximate)) <= roundoff)
                    return interior;
                const auto trimmed = trimmedProjection(face.shape, approximate);
                // Both searches already enforce uniqueness within their own
                // candidate sets. Independent numerical approximations of one
                // minimum are not additional geometric roots.
                return point(trimmed.position).Distance(point(approximate)) <
                               minimum->position.Distance(point(approximate))
                           ? trimmed
                           : interior;
            }
        }
        return trimmedProjection(face.shape, approximate);
    });
}
SurfaceProjection CadGeometry::projectSupportingSurface(std::uint32_t id, Vec3 approximate) const {
    return cadOperation("CAD-Traegereflaechenprojektion", [&] {
        const auto &face = data_->face(id);
        const auto context = "CAD-Traegereflaeche " + std::to_string(id);
        const auto minimum =
            surfaceMinimum(*face.supportingProjectors, face.supportingBounds, approximate, true);
        if (!minimum)
            throw SurfaceProjectionDomainError(context + ": kein eindeutiges lokales Minimum.");
        gp_Pnt reconstructed;
        gp_Vec du, dv;
        face.supportingSurface->D1(minimum->u, minimum->v, reconstructed, du, dv);
        finite(value(reconstructed));
        const gp_Vec cross = du.Crossed(dv);
        const double jacobian = cross.Magnitude();
        if (!std::isfinite(jacobian))
            throw std::runtime_error(context + ": nicht endliche Flaechenableitungen.");
        if (!(jacobian > 0))
            throw SurfaceProjectionDomainError(context + ": keine regulaere Tangentialebene.");
        gp_Vec normal = cross / jacobian;
        const gp_Vec residual(minimum->position, point(approximate));
        const double tangent = (residual - normal * residual.Dot(normal)).Magnitude();
        const auto p = value(minimum->position);
        const double roundoff =
            128 * std::numeric_limits<double>::epsilon() *
            std::max({1.0, std::abs(p.x), std::abs(p.y), std::abs(p.z), std::abs(approximate.x),
                      std::abs(approximate.y), std::abs(approximate.z)});
        // MIN is searched in a finite UV box. Its boundary must not masquerade
        // as a regular closest-surface map: require a stationary solution.
        // Propagate the existing parameter precision through the local metric.
        const double stationaryTolerance =
            Precision::PConfusion() * (du.Magnitude() + dv.Magnitude()) + roundoff;
        if (!std::isfinite(tangent) || !std::isfinite(stationaryTolerance))
            throw std::runtime_error(context + ": nicht endliche Projektionsdaten.");
        if (tangent > stationaryTolerance) {
            std::ostringstream message;
            message << std::setprecision(17) << context
                    << ": keine stationaere Projektion im Suchbereich; Punkt [m] (" << approximate.x
                    << ", " << approximate.y << ", " << approximate.z << "), UV (" << minimum->u
                    << ", " << minimum->v << "), Tangentialrest [m] " << tangent << ", Grenze [m] "
                    << stationaryTolerance << ".";
            throw SurfaceProjectionDomainError(message.str());
        }
        if (face.shape.Orientation() == TopAbs_REVERSED)
            normal.Reverse();
        return SurfaceProjection{p, {normal.X(), normal.Y(), normal.Z()}};
    });
}
double CadGeometry::distanceToSurface(std::uint32_t id, Vec3 position) const {
    return cadOperation("CAD-Flaechenabstand", [&] {
        finite(position);
        const auto &face = data_->face(id);
        // A scalar-distance witness may be any minimizing surface point;
        // uniqueness remains mandatory for the separate position/normal map.
        const auto minimum =
            face.projectors->revolved()
                ? std::optional{face.projectors->revolved()->minimum(position, false)}
                : surfaceMinimum(*face.projectors, face.bounds, position);
        std::optional<double> interiorDistance;
        if (minimum) {
            const auto state = face.trimClassifier->Perform(gp_Pnt2d(minimum->u, minimum->v));
            if (state == TopAbs_IN || state == TopAbs_ON)
                interiorDistance = minimum->position.Distance(point(position));
        }
        const auto vertex = BRepBuilderAPI_MakeVertex(point(position)).Vertex();
        const BRepExtrema_DistShapeShape projection(vertex, face.shape);
        if (!projection.IsDone() || projection.NbSolution() == 0 ||
            !std::isfinite(projection.Value()) || projection.Value() < 0)
            throw std::runtime_error("Kein endlicher Abstand zur getrimmten CAD-Flaeche.");
        // Multiple equally near trim points share one distance. No arbitrary
        // support point or normal is selected for this scalar quality query.
        return interiorDistance ? std::min(*interiorDistance, projection.Value())
                                : projection.Value();
    });
}
Vec3 CadGeometry::projectJunction(std::span<const std::uint32_t> ids, Vec3 approximate) const {
    return cadOperation("CAD-Knotenprojektion", [&] {
        finite(approximate);
        if (ids.size() < 3)
            throw std::invalid_argument("Ein CAD-Knoten erfordert mindestens drei Flaechen.");
        const std::set<std::uint32_t> unique(ids.begin(), ids.end());
        if (unique.size() != ids.size())
            throw std::invalid_argument("CAD-Knotenreferenz enthaelt doppelte Flaechen.");
        std::vector<TopTools_IndexedMapOfShape> vertices(ids.size());
        for (std::size_t i = 0; i < ids.size(); ++i)
            TopExp::MapShapes(data_->face(ids[i]).shape, TopAbs_VERTEX, vertices[i]);
        std::optional<Vec3> result;
        for (int i = 1; i <= vertices.front().Extent(); ++i) {
            const auto &candidate = vertices.front()(i);
            if (!std::all_of(vertices.begin() + 1, vertices.end(),
                             [&](const auto &map) { return map.Contains(candidate); }))
                continue;
            const auto vertex = TopoDS::Vertex(candidate);
            const auto position = BRep_Tool::Pnt(vertex);
            const double modelTolerance = BRep_Tool::Tolerance(vertex);
            const auto p = value(position);
            finite(p);
            if (!std::isfinite(modelTolerance) || modelTolerance <= 0)
                throw std::runtime_error(
                    "CAD-Knoten besitzt keine gueltige Modellierungstoleranz.");
            // OCCT vertices represent intersections within their own modelled
            // tolerance. The additional term covers only coordinate roundoff.
            const double roundoff =
                128 * std::numeric_limits<double>::epsilon() *
                std::max({1.0, std::abs(p.x), std::abs(p.y), std::abs(p.z), std::abs(approximate.x),
                          std::abs(approximate.y), std::abs(approximate.z)});
            const double tolerance = modelTolerance + roundoff;
            if (position.Distance(point(approximate)) > tolerance)
                continue;
            for (auto id : ids) {
                const auto &face = data_->face(id);
                // Validate this topological vertex's own face parameters. A
                // separate nearest-point search can select another stationary
                // point and does not test the modelled vertex representation.
                const auto uv = BRep_Tool::Parameters(vertex, face.shape);
                if (!std::isfinite(uv.X()) || !std::isfinite(uv.Y()))
                    throw std::runtime_error(
                        "CAD-Knoten besitzt keine endlichen Flaechenparameter.");
                const auto surface = face.surface->Value(uv.X(), uv.Y());
                finite(value(surface));
                if (position.Distance(surface) > tolerance)
                    throw std::runtime_error(
                        "CAD-Knoten liegt ausserhalb seiner modellierten Flaechentoleranz.");
            }
            if (result)
                throw std::runtime_error(
                    "CAD-Knotenreferenz ist innerhalb der Modellierungstoleranz mehrdeutig.");
            result = p;
        }
        if (!result)
            throw std::invalid_argument(
                "Kein gemeinsamer topologischer CAD-Knoten nahe der Quellposition.");
        return *result;
    });
}
Vec3 CadGeometry::projectBoundary(std::span<const std::uint32_t> ids, Vec3 approximate) const {
    return cadOperation("CAD-Kantenprojektion", [&] {
        finite(approximate);
        if (ids.size() < 2)
            throw std::invalid_argument(
                "Eine CAD-Kante erfordert mindestens zwei verschiedene Flaechen.");
        std::set<std::uint32_t> unique(ids.begin(), ids.end());
        if (unique.size() != ids.size())
            throw std::invalid_argument("CAD-Kantenreferenz enthaelt doppelte Flaechen.");
        std::vector<const Data::Face *> faces;
        for (auto id : ids)
            faces.push_back(&data_->face(id));
        std::vector<gp_Pnt> candidates;
        for (int i = 1; i <= faces.front()->edges.Extent(); ++i) {
            const auto edge = TopoDS::Edge(faces.front()->edges(i));
            if (BRep_Tool::Degenerated(edge))
                continue;
            if (!std::all_of(faces.begin() + 1, faces.end(),
                             [&](auto face) { return face->edges.Contains(edge); }))
                continue;
            const auto points = nearestCurvePoints(edge, approximate);
            candidates.insert(candidates.end(), points.begin(), points.end());
        }
        if (candidates.empty())
            throw std::invalid_argument(
                "Die CAD-Flaechen besitzen keine gemeinsame topologische Kante.");
        return value(uniqueNearest(candidates, approximate));
    });
}
} // namespace hornsculpt
