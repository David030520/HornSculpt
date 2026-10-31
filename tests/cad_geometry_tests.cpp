#include "../modules/meshing/src/surface_charts.hpp"
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <DESTEP_Parameters.hxx>
#include <GProp_GProps.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <STEPControl_Writer.hxx>
#include <StepData_StepModel.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <array>
#include <barrier>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <future>
#include <gp_Ax3.hxx>
#include <gp_Pln.hxx>
#include <hornsculpt/cad_geometry.hpp>
#include <hornsculpt/meshing.hpp>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
void expectNear(double a, double b, double tolerance, const char *message) {
    require(std::isfinite(a) && std::abs(a - b) <= tolerance, message);
}
template <class Function> void rejects(Function f) {
    bool rejected = false;
    try {
        f();
    } catch (const std::exception &) {
        rejected = true;
    }
    require(rejected, "Invalid CAD input/projection was accepted.");
}
double distance(Vec3 a, Vec3 b) {
    return std::hypot(a.x - b.x, a.y - b.y, a.z - b.z);
}
CadBrepPart artifact(const TopoDS_Shape &shape, std::uint32_t firstId) {
    std::ostringstream data;
    BRepTools::Write(shape, data, Standard_False, Standard_False, TopTools_FormatVersion_CURRENT);
    TopTools_IndexedMapOfShape faces;
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    std::vector<std::uint32_t> ids(static_cast<std::size_t>(faces.Extent()));
    std::iota(ids.begin(), ids.end(), firstId);
    return {data.str(), ids};
}
std::shared_ptr<const CadGeometry> geometry(const TopoDS_Shape &shape, double scale = 1,
                                            std::uint32_t firstId = 101) {
    const auto part = artifact(shape, firstId);
    return CadGeometry::fromBrepParts(std::span{&part, 1}, scale);
}
std::uint32_t faceWithNormal(const CadGeometry &source, Vec3 normal) {
    for (const auto &face : source.faces())
        if (distance(face.normal, normal) < 1e-10)
            return face.faceId;
    throw std::runtime_error("Missing expected planar face.");
}
void closedMesh(const CadGeometry &source, const SurfaceMesh &mesh, int euler, double tolerance) {
    require(!mesh.vertices.empty() && !mesh.triangles.empty(), "CAD mesh is empty.");
    std::map<std::array<std::uint32_t, 2>, std::pair<unsigned, int>> edges;
    std::set<std::uint32_t> used, ids;
    for (const auto &face : source.faces())
        ids.insert(face.faceId);
    for (const auto &t : mesh.triangles) {
        require(ids.contains(t.faceId) && t.boundary == BoundaryCondition::wall,
                "CAD source face IDs/Wall tags changed.");
        for (unsigned i = 0; i < 3; ++i) {
            auto a = t.vertices[i], b = t.vertices[(i + 1) % 3];
            require(a < mesh.vertices.size() && b < mesh.vertices.size() && a != b,
                    "Invalid CAD mesh edge.");
            auto &edge = edges[{std::min(a, b), std::max(a, b)}];
            ++edge.first;
            edge.second += a < b ? 1 : -1;
            used.insert(a);
        }
        const auto a = mesh.vertices[t.vertices[0]], b = mesh.vertices[t.vertices[1]],
                   c = mesh.vertices[t.vertices[2]];
        const Vec3 centroid{(a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3, (a.z + b.z + c.z) / 3};
        const auto exact = source.projectSurface(t.faceId, centroid);
        require(distance(centroid, exact.position) <= 1.01 * tolerance,
                "CAD triangle exceeds its geometric deflection.");
        const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
        const Vec3 n{u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
        require(n.x * exact.normal.x + n.y * exact.normal.y + n.z * exact.normal.z > 0,
                "CAD mesh normal points into the solid.");
    }
    for (const auto &[key, use] : edges) {
        (void)key;
        require(use == std::pair<unsigned, int>{2, 0},
                "CAD mesh is not closed and consistently oriented.");
    }
    require(static_cast<long long>(used.size()) - static_cast<long long>(edges.size()) +
                    static_cast<long long>(mesh.triangles.size()) ==
                euler,
            "CAD mesh changed the source genus/component count.");
}
void scaledNotchedPrismProjection() {
    // A small concave planar trim exposes scale-sensitive ray classification.
    // All coordinates define this test; no external CAD document is required.
    const std::array<gp_Pnt, 10> points{{{-.201402311655776, .143092643087061, .23008497937},
                                         {-.0926627189387369, .143092643087061, .23008497937},
                                         {-.0911180398147984, .143820905064740, .230738157970979},
                                         {-.0908418385778138, .143092643087061, .23008497937},
                                         {.0908418385778138, .143092643087061, .23008497937},
                                         {.0911180398147984, .143820905064740, .230738157970979},
                                         {.0926627189387369, .143092643087061, .23008497937},
                                         {.201402311655776, .143092643087061, .23008497937},
                                         {.296, .17282855784, .256755139},
                                         {-.296, .17282855784, .256755139}}};
    const gp_Pnt origin(-6.062414068402007e-18, .1603043827203048, .24552219888527468);
    const gp_Dir normal(0, .6676889841472743, -.7444403404225088);
    const gp_Dir xDirection(0, -.7444403404225088, -.6676889841472743);
    const Vec3 witness{.21563863646359341, .15885395562284071, .24422130995878819};
    for (const double units : {1., 1000.}) {
        BRepBuilderAPI_MakePolygon polygon;
        for (const auto &p : points)
            polygon.Add(gp_Pnt(p.X() * units, p.Y() * units, p.Z() * units));
        polygon.Close();
        const gp_Pln plane(
            gp_Ax3(gp_Pnt(origin.X() * units, origin.Y() * units, origin.Z() * units), normal,
                   xDirection));
        const auto face = BRepBuilderAPI_MakeFace(plane, polygon.Wire()).Face();
        const auto solid = BRepPrimAPI_MakePrism(face, gp_Vec(normal) * (.005 * units)).Shape();
        const auto source = geometry(solid, 1 / units);
        const auto target = std::min_element(
            source->faces().begin(), source->faces().end(), [&](const auto &a, const auto &b) {
                const auto planeDistance = [&](const auto &f) {
                    return std::abs((f.center.y - origin.Y()) * normal.Y() +
                                    (f.center.z - origin.Z()) * normal.Z());
                };
                return planeDistance(a) < planeDistance(b);
            });
        require(target != source->faces().end(), "Notched prism cap is absent.");
        const auto projected = source->projectSurface(target->faceId, witness);
        require(distance(projected.position, witness) < 1e-10,
                "SI trim classification rejected a proven interior point.");
        expectNear(source->distanceToSurface(target->faceId, witness), 0, 1e-10,
                   "SI scalar trim distance changed a proven zero witness.");
        require(source->distanceToSurface(target->faceId, {.4, witness.y, witness.z}) > .1,
                "Trim classification accepted a proven exterior point.");
    }
}

void curvedTrimCompetesWithStationaryMinimum() {
    // The widening graph has a shallow minimum. Just outside its side trim,
    // a remote stationary surface point must not beat a closer boundary point.
    constexpr double length = .08665;
    const auto h = [](double z) { return .025 - .0667 * z + 4 * z * z; };
    const auto y = [](double z) { return .032 - .0667 * z + 2 * z * z; };
    const std::array widths{h(0), h(0) - .0667 * length / 2, h(length)};
    const std::array heights{y(0), y(0) - .0667 * length / 2, y(length)};
    TColgp_Array2OfPnt poles(1, 3, 1, 2);
    for (int u = 1; u <= 3; ++u)
        for (int v = 1; v <= 2; ++v)
            poles.SetValue(u, v,
                           gp_Pnt((v == 1 ? -1 : 1) * widths[std::size_t(u - 1)],
                                  heights[std::size_t(u - 1)], length * (u - 1) / 2));
    TColStd_Array1OfReal uKnots(1, 2), vKnots(1, 2);
    uKnots.SetValue(1, 0);
    uKnots.SetValue(2, length);
    vKnots.SetValue(1, -1);
    vKnots.SetValue(2, 1);
    TColStd_Array1OfInteger uMults(1, 2), vMults(1, 2);
    uMults.Init(3);
    vMults.Init(2);
    Handle(Geom_BSplineSurface) surface =
        new Geom_BSplineSurface(poles, uKnots, vKnots, uMults, vMults, 2, 1, false, false);
    const auto face = BRepBuilderAPI_MakeFace(surface, 1e-10).Face();
    GProp_GProps original;
    BRepGProp::SurfaceProperties(face, original);
    const auto solid = BRepPrimAPI_MakePrism(face, gp_Vec(0, .005, 0)).Shape();
    const auto source = geometry(solid);
    const Vec3 center{original.CentreOfMass().X(), original.CentreOfMass().Y(),
                      original.CentreOfMass().Z()};
    const auto target = std::min_element(
        source->faces().begin(), source->faces().end(), [&](const auto &a, const auto &b) {
            return distance(a.center, center) < distance(b.center, center);
        });
    require(target != source->faces().end(), "Curved prism cap is absent.");
    for (int i = 1; i <= 20; ++i) {
        const double z = .001 * i;
        for (const double offset : {1e-5, 5e-5, 1e-4}) {
            const Vec3 query{h(z) + offset, y(z), z};
            const auto projected = source->projectSurface(target->faceId, query);
            require(distance(projected.position, query) <= offset + 1e-10,
                    "A farther stationary minimum overruled a proven trim-edge witness.");
            require(source->distanceToSurface(target->faceId, query) <= offset + 1e-10,
                    "Scalar surface distance overruled a proven trim-edge witness.");
        }
    }
}

void boxUnitsAndProjection() {
    const auto source = geometry(BRepPrimAPI_MakeBox(20, 30, 40).Shape(), .001);
    require(source->faces().size() == 6, "Box face count changed.");
    double area = 0;
    for (const auto &f : source->faces()) {
        area += f.area;
        require(f.surfaceType == "Plane", "Plane type normalization failed.");
        require(f.edgeCount == 4 && f.wireCount == 1, "Unique box edge/wire counts changed.");
        require(f.faceId >= 101 && f.faceId <= 106, "Explicit artifact IDs changed.");
    }
    expectNear(area, .0052, 1e-12, "BRep millimetres did not convert to SI area.");
    const auto xp = faceWithNormal(*source, {1, 0, 0}), xm = faceWithNormal(*source, {-1, 0, 0}),
               yp = faceWithNormal(*source, {0, 1, 0});
    auto projection = source->projectSurface(xp, {.025, .015, .020});
    expectNear(projection.position.x, .020, 1e-12, "Exact plane projection failed.");
    expectNear(projection.normal.x, 1, 1e-12, "Oriented face normal failed.");
    projection = source->projectSurface(xp, {.025, .04, .06});
    require(distance(projection.position, {.020, .030, .040}) < 1e-11,
            "Projection ignored trimmed face bounds.");
    require(distance(projection.normal, {1, 0, 0}) < 1e-12,
            "Trimmed corner projection lost the face-oriented vertex normal.");
    projection = source->projectSurface(xp, {.025, .04, .018});
    require(distance(projection.position, {.020, .030, .018}) < 1e-11 &&
                distance(projection.normal, {1, 0, 0}) < 1e-12,
            "Trimmed edge projection lost its face-specific support parameters.");
    const std::array<std::uint32_t, 2> edge{xp, yp};
    const auto edgeEndpoint = source->projectBoundary(edge, {.025, .04, .06});
    require(distance(edgeEndpoint, {.020, .030, .040}) < 1e-11,
            "A bounded line projection ignored its endpoint minimum.");
    const auto boundary = source->projectBoundary(edge, {.025, .04, .018});
    require(distance(boundary, {.020, .030, .018}) < 1e-11,
            "Common topological edge projection failed.");
    rejects([&] {
        const std::array<std::uint32_t, 2> ids{xp, xm};
        (void)source->projectBoundary(ids, {});
    });
    rejects([&] {
        const std::array<std::uint32_t, 2> ids{xp, xp};
        (void)source->projectBoundary(ids, {});
    });
    rejects([&] { (void)source->projectSurface(9999, {}); });
    rejects([&] {
        (void)source->projectSurface(xp, {std::numeric_limits<double>::quiet_NaN(), 0, 0});
    });
    const auto restored = CadGeometry::fromBrepParts(source->brepParts(), 1);
    require(restored->faces().size() == source->faces().size(),
            "Canonical BRep roundtrip changed faces.");
    for (std::size_t i = 0; i < source->faces().size(); ++i) {
        const auto &a = source->faces()[i], &b = restored->faces()[i];
        require(a.faceId == b.faceId, "Canonical BRep roundtrip changed artifact IDs.");
        expectNear(a.area, b.area, 1e-14, "Canonical BRep units changed.");
        require(distance(a.center, b.center) < 1e-12, "Canonical BRep centroid changed.");
    }
    closedMesh(*source, source->previewMesh(.00025), 2, .00025);
    std::stop_source stop;
    stop.request_stop();
    rejects([&] { (void)source->previewMesh(.00025, stop.get_token()); });
    rejects([&] { (void)source->previewMesh(0); });
}
void thinCurvedRimTessellation() {
    // Offset seams make coarse chords intersect across an otherwise valid thin
    // annulus. The two cases exercise partial and missing OCCT face meshes.
    constexpr double radius = .076, gap = .0001, phase = .07, tolerance = .00025;
    for (const double shift : {0., .00003}) {
        const auto outer = BRepPrimAPI_MakeCylinder(radius, .01).Shape();
        const gp_Ax2 axis(gp_Pnt(shift, 0, -.01), gp_Dir(0, 0, 1),
                          gp_Dir(std::cos(phase), std::sin(phase), 0));
        const auto inner = BRepPrimAPI_MakeCylinder(axis, radius - gap - shift, .03).Shape();
        const auto source = geometry(BRepAlgoAPI_Cut(outer, inner).Shape());
        const auto original = source->brepParts().front().data;
        const auto input = source->meshingSource(tolerance);
        closedMesh(*source, input.mesh, 0, tolerance);
        const auto quality = inspectMesh(input.mesh);
        require(quality.watertight() && quality.connectedComponents == 1 &&
                    quality.selfIntersectionPairs == 0 && quality.signedVolume > 0,
                "Thin curved CAD rim lost its oriented, embedded boundary.");
        require(source->brepParts().front().data == original && input.parameters,
                "Tessellation refinement changed CAD geometry or lost native parameters.");
        std::stop_source cancelled;
        cancelled.request_stop();
        rejects([&] { (void)source->meshingSource(tolerance, cancelled.get_token()); });
    }
}
void tolerantTopologicalJunction() {
    auto shape = BRepPrimAPI_MakeBox(.02, .03, .04).Shape();
    TopTools_IndexedMapOfShape vertices;
    TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
    const Vec3 displaced{2e-9, 26e-9, 66e-9};
    constexpr double vertexTolerance = 1e-7;
    bool changed = false;
    for (int i = 1; i <= vertices.Extent(); ++i) {
        auto v = TopoDS::Vertex(vertices(i));
        if (BRep_Tool::Pnt(v).Distance(gp_Pnt(0, 0, 0)) > 1e-12)
            continue;
        // Imported BReps may have one topological intersection point which is
        // not exactly on each underlying surface, within its native tolerance.
        BRep_Builder builder;
        builder.UpdateVertex(v, gp_Pnt(displaced.x, displaced.y, displaced.z), vertexTolerance);
        changed = true;
    }
    require(changed, "Missing box junction fixture vertex.");
    const auto source = geometry(shape);
    for (const auto &corner : source->boundaryCorners())
        if (distance(corner.position, displaced) < 1e-12) {
            require(!corner.constraint,
                    "Native model tolerance was used to manufacture a coplanar line constraint.");
        }
    const std::array ids{faceWithNormal(*source, {-1, 0, 0}), faceWithNormal(*source, {0, -1, 0}),
                         faceWithNormal(*source, {0, 0, -1})};
    double largestOffset = 0;
    for (auto id : ids) {
        const double offset = distance(source->projectSurface(id, displaced).position, displaced);
        require(offset <= vertexTolerance, "Fixture exceeds its native vertex tolerance.");
        largestOffset = std::max(largestOffset, offset);
    }
    require(largestOffset > 50e-9, "Fixture no longer exercises nanometre surface offsets.");
    require(distance(source->projectJunction(ids, displaced), displaced) < 1e-14,
            "Valid shared CAD vertex was rejected or moved onto one incident surface.");
    const Vec3 nearby{displaced.x + 2e-8, displaced.y, displaced.z};
    require(distance(source->projectJunction(ids, nearby), displaced) < 1e-14,
            "A local junction query did not recover its common topology vertex.");
    rejects([&] { (void)source->projectJunction(ids, {.001, 0, 0}); });
    rejects([&] {
        const std::array incompatible{ids[0], ids[1], faceWithNormal(*source, {1, 0, 0})};
        (void)source->projectJunction(incompatible, displaced);
    });
    rejects([&] {
        const std::array duplicate{ids[0], ids[1], ids[1]};
        (void)source->projectJunction(duplicate, displaced);
    });
    rejects([&] { (void)source->projectJunction(std::span(ids).first(2), displaced); });
    rejects([&] {
        (void)source->projectJunction(ids, {std::numeric_limits<double>::quiet_NaN(), 0, 0});
    });
}
void cylinderAndTrimmedBore() {
    const auto squarePlate = BRepPrimAPI_MakeBox(.04, .04, .01).Shape();
    const auto squareHole = BRepPrimAPI_MakeBox(gp_Pnt(.015, .015, -.01), .01, .01, .03).Shape();
    const auto squareBore = geometry(BRepAlgoAPI_Cut(squarePlate, squareHole).Shape());
    const auto squareBottom = faceWithNormal(*squareBore, {0, 0, -1});
    const Vec3 equidistant{.02, .02, -.003};
    expectNear(squareBore->distanceToSurface(squareBottom, equidistant), std::hypot(.005, .003),
               1e-11, "Four equally near physical trim edges invalidated their scalar distance.");
    rejects([&] { (void)squareBore->projectSurface(squareBottom, equidistant); });
    const auto cylinder = geometry(BRepPrimAPI_MakeCylinder(.02, .03).Shape());
    const auto coarse = cylinder->previewMesh(.0005), fine = cylinder->previewMesh(.00005);
    closedMesh(*cylinder, coarse, 2, .0005);
    closedMesh(*cylinder, fine, 2, .00005);
    require(fine.triangles.size() > coarse.triangles.size(),
            "CAD deflection did not refine curved faces.");
    const auto outside = BRepPrimAPI_MakeCylinder(.02, .03).Shape();
    const auto cutter =
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -.01), gp_Dir(0, 0, 1)), .004, .05).Shape();
    BRepAlgoAPI_Cut cut(outside, cutter);
    require(cut.IsDone(), "OCC could not build actual bore fixture.");
    const auto bore = geometry(cut.Shape(), 1, 201);
    std::uint32_t inner = 0, bottom = 0;
    for (const auto &f : bore->faces()) {
        if (f.surfaceType == "Cylinder" && f.area < .001)
            inner = f.faceId;
        if (f.surfaceType == "Plane" && f.normal.z < -.9)
            bottom = f.faceId;
    }
    require(inner != 0 && bottom != 0, "Missing cylindrical bore/trimmed cap faces.");
    // The axis has an entire circle of nearest points. Distance is unique even
    // though an interior surface map/normal must remain rejected there.
    expectNear(bore->distanceToSurface(bottom, {0, 0, -.003}), .005, 1e-11,
               "Multiple nearest hole-rim points invalidated their unique distance.");
    expectNear(bore->distanceToSurface(inner, {0, 0, .015}), .004, 1e-11,
               "An ambiguous cylindrical normal invalidated the scalar face distance.");
    expectNear(bore->distanceToSurface(bottom, {.01, 0, -.003}), .003, 1e-11,
               "Interior distance differs from the nearest trimmed plane.");
    rejects([&] { (void)bore->distanceToSurface(9999, {}); });
    rejects([&] {
        (void)bore->distanceToSurface(bottom, {std::numeric_limits<double>::quiet_NaN(), 0, 0});
    });
    const auto wall = bore->projectSurface(inner, {.003, .001, .015});
    expectNear(std::hypot(wall.position.x, wall.position.y), .004, 1e-11,
               "Inner cylindrical surface projection failed.");
    require(wall.normal.x * wall.position.x + wall.normal.y * wall.position.y < 0,
            "Bore normal must point toward the acoustic hole.");
    const auto cap = bore->projectSurface(bottom, {.002, .001, -.003});
    expectNear(std::hypot(cap.position.x, cap.position.y), .004, 1e-11,
               "Trimmed cap projection filled the bore.");
    expectNear(cap.position.z, 0, 1e-11, "Trimmed cap projection left its plane.");
    const auto supportingCap = bore->projectSupportingSurface(bottom, {.002, .001, -.003});
    expectNear(distance(supportingCap.position, {.002, .001, 0}), 0, 1e-12,
               "Supporting cap was clamped to its curved hole trim.");
    require(supportingCap.normal.z < -.99, "Supporting cap changed its oriented normal.");
    const auto exteriorPlane = bore->projectSupportingSurface(bottom, {.04, 0, -.003});
    expectNear(distance(exteriorPlane.position, {.04, 0, 0}), 0, 1e-12,
               "Unbounded supporting plane retained a physical cap clamp.");
    const auto exteriorCylinder = bore->projectSupportingSurface(inner, {.003, .001, .04});
    expectNear(exteriorCylinder.position.z, .04, 1e-12,
               "Unbounded supporting cylinder retained an axial face clamp.");
    rejects([&] { (void)bore->projectSupportingSurface(inner, {0, 0, .015}); });
    rejects([&] { (void)bore->projectSupportingSurface(9999, {}); });
    rejects([&] {
        (void)bore->projectSupportingSurface(bottom,
                                             {std::numeric_limits<double>::quiet_NaN(), 0, 0});
    });
    const std::array<std::uint32_t, 2> ids{inner, bottom};
    const auto edge = bore->projectBoundary(ids, {.005, .001, -.003});
    expectNear(std::hypot(edge.x, edge.y), .004, 1e-11,
               "Hole boundary projection selected the wrong loop.");
    expectNear(edge.z, 0, 1e-11, "Hole boundary projection left the exact circle.");
    closedMesh(*bore, bore->previewMesh(.00015), 0, .00015);
    std::vector<std::future<void>> workers;
    for (unsigned i = 0; i < 4; ++i)
        workers.push_back(std::async(std::launch::async, [&, i] {
            for (unsigned j = 0; j < 32; ++j) {
                const double phi = (i + j) * .21;
                const auto p =
                    bore->projectSurface(inner, {.003 * std::cos(phi), .003 * std::sin(phi), .015});
                expectNear(std::hypot(p.position.x, p.position.y), .004, 1e-11,
                           "Concurrent exact projection changed geometry.");
            }
        }));
    for (auto &worker : workers)
        worker.get();
}
void sphericalCapParameters() {
    for (const double sign : {-1., 1.}) {
        const gp_Ax2 placement(gp_Pnt(.003, .001, -.01), gp_Dir(1, 2, 3));
        constexpr double radius = .02;
        const auto shape =
            BRepPrimAPI_MakeSphere(placement, radius, sign < 0 ? -std::numbers::pi / 2 : .7,
                                   sign < 0 ? -.7 : std::numbers::pi / 2)
                .Shape();
        const auto cad = geometry(shape);
        const auto source = cad->meshingSource(.0002);
        const auto sphere = std::find_if(cad->faces().begin(), cad->faces().end(),
                                         [](const auto &f) { return f.surfaceType == "Sphere"; });
        require(sphere != cad->faces().end(), "Spherical cap regression lacks its sphere face.");
        const auto domain = source.parameters->domain(sphere->faceId);
        require(!domain.uPeriod && !domain.vPeriod,
                "A spherical pole retained its singular periodic coordinate.");
        const auto pole = source.parameters->surface(sphere->faceId, {});
        const auto expected =
            placement.Location().Translated(gp_Vec(placement.Direction()) * (sign * radius));
        expectNear(distance(pole.position, {expected.X(), expected.Y(), expected.Z()}), 0, 1e-12,
                   "Spherical Cartesian pole moved its exact geometry.");
        for (const auto uv : {ParameterPoint{}, ParameterPoint{.002, .003}}) {
            const auto p = source.parameters->surface(sphere->faceId, uv);
            const gp_Vec du(p.du.x, p.du.y, p.du.z), dv(p.dv.x, p.dv.y, p.dv.z);
            require(du.Crossed(dv).Magnitude() > .99,
                    "Spherical Cartesian metric collapses at its pole.");
            const auto h = 1e-7;
            const auto left = source.parameters->surface(sphere->faceId, {uv.u - h, uv.v});
            const auto right = source.parameters->surface(sphere->faceId, {uv.u + h, uv.v});
            const Vec3 derivative{(right.position.x - left.position.x) / (2 * h),
                                  (right.position.y - left.position.y) / (2 * h),
                                  (right.position.z - left.position.z) / (2 * h)};
            expectNear(distance(derivative, p.du), 0, 1e-7,
                       "Spherical Cartesian D1 disagrees with exact positions.");
        }
        MeshingOptions options;
        options.targetEdgeLength = .0035;
        options.adaptive.minimumEdgeLength = .0005;
        options.adaptive.geometricTolerance = .0002;
        options.geometry = cad;
        options.parameterization = source.parameters;
        const auto result = remeshIsotropic(source.mesh, options);
        require(result.quality.watertight() && result.quality.minimumAngleDegrees >= 20 &&
                    result.quality.maximumAspectRatio <= 4,
                "Spherical cap quality mesh failed after pole regularization.");
    }
}
void revolvedSplineNearestMinimum() {
    // This short near-rim minimum was omitted by OCCT's full-box 2D Grad and
    // Tree extrema searches: they returned 223.223 um although an IN surface
    // point was already only 92.667 um away. The meridian minimum is 62.712 um.
    const std::array heights{
        -.000816426954,   -.000522903759973, -.00258731901746, -.00237519465219, -.00378240905996,
        -.00456615778370, -.00488209381670,  -.00507977286729, -.00501878964159, -.00502532847567,
        -.00468412481636, -.00412020665630,  -.00230948687299, -.00167409180989, -.00116230876761,
        -.00129949604632, -.0011123818};
    TColStd_Array1OfReal knots(1, 15);
    TColStd_Array1OfInteger multiplicities(1, 15);
    std::vector<double> expanded;
    for (int i = 0; i < 15; ++i) {
        knots(i + 1) = i / 14.;
        multiplicities(i + 1) = i == 0 || i == 14 ? 4 : 1;
        for (int j = 0; j < multiplicities(i + 1); ++j)
            expanded.push_back(knots(i + 1));
    }
    TColgp_Array1OfPnt poles(1, 17);
    for (std::size_t i = 0; i < heights.size(); ++i) {
        const double greville = (expanded[i + 1] + expanded[i + 2] + expanded[i + 3]) / 3;
        poles(static_cast<int>(i + 1)) = gp_Pnt(.0585 + .012 * greville, 0, heights[i]);
    }
    const Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, knots, multiplicities, 3);
    const gp_Pnt outerBottom(.0705, 0, -.01), innerBottom(.0585, 0, -.01);
    BRepBuilderAPI_MakeWire wire;
    wire.Add(BRepBuilderAPI_MakeEdge(curve).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(poles(17), outerBottom).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(outerBottom, innerBottom).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(innerBottom, poles(1)).Edge());
    const auto shape = BRepPrimAPI_MakeRevol(BRepBuilderAPI_MakeFace(wire.Wire()).Face(),
                                             gp_Ax1(gp_Pnt(), gp_Dir(0, 0, 1)))
                           .Shape();
    const auto source = geometry(shape);
    const auto findMeridianFace = [](const auto &cad) {
        const auto found =
            std::find_if(cad->faces().begin(), cad->faces().end(),
                         [](const auto &f) { return f.surfaceType == "SurfaceOfRevolution"; });
        require(found != cad->faces().end(), "Missing revolved spline regression face.");
        return found->faceId;
    };
    const auto face = findMeridianFace(source);
    const Vec3 query{-.050708470182465465, -.029612689777137602, -.00091525900771126845};
    const auto projection = source->projectSurface(face, query);
    const double error = distance(projection.position, query);
    expectNear(error, .000062712138037644187, 1e-10,
               "Revolved spline projection missed its narrow nearest meridian minimum.");
    expectNear(source->distanceToSurface(face, query), error, 1e-12,
               "Scalar revolved distance lost the nearest narrow meridian minimum.");
    const double axisDistance = source->distanceToSurface(face, {0, 0, query.z});
    require(std::isfinite(axisDistance) && axisDistance > 0,
            "A circle of equally near revolved points invalidated scalar distance.");
    const auto supporting = source->projectSupportingSurface(face, query);
    expectNear(distance(supporting.position, projection.position), 0, 1e-10,
               "Supporting and interior trimmed minima disagree on a regular meridian root.");
    rejects([&] { (void)source->projectSurface(face, {0, 0, query.z}); });
    rejects([&] { (void)source->projectSupportingSurface(face, {0, 0, query.z}); });
    gp_Trsf placement;
    placement.SetRotation(gp_Ax1(gp_Pnt(), gp_Dir(1, 2, 3)), .43);
    placement.SetTranslationPart(gp_Vec(.3, -.2, .6));
    const auto placed = geometry(BRepBuilderAPI_Transform(shape, placement, Standard_True).Shape());
    gp_Pnt movedQuery(query.x, query.y, query.z);
    movedQuery.Transform(placement);
    gp_Pnt movedExpected(projection.position.x, projection.position.y, projection.position.z);
    movedExpected.Transform(placement);
    const auto moved = placed->projectSurface(findMeridianFace(placed),
                                              {movedQuery.X(), movedQuery.Y(), movedQuery.Z()});
    expectNear(distance(moved.position, {movedExpected.X(), movedExpected.Y(), movedExpected.Z()}),
               0, 1e-10, "Rigid placement changed the nearest revolved spline minimum.");
    gp_Pnt movedAxis(0, 0, query.z);
    movedAxis.Transform(placement);
    expectNear(placed->distanceToSurface(findMeridianFace(placed),
                                         {movedAxis.X(), movedAxis.Y(), movedAxis.Z()}),
               axisDistance, 1e-10, "Rigid placement changed an ambiguous-map scalar distance.");
    std::cout << "Revolved spline narrow minimum: " << error << " m\n";
}

void supportingSpherePole() {
    const auto source = geometry(BRepPrimAPI_MakeSphere(.02).Shape());
    require(source->faces().size() == 1, "Unexpected analytic sphere topology.");
    const auto face = source->faces().front().faceId;
    const auto pole = source->projectSupportingSurface(face, {0, 0, .021});
    require(distance(pole.position, {0, 0, .02}) < 1e-12 &&
                distance(pole.normal, {0, 0, 1}) < 1e-12,
            "A coordinate-singular sphere pole lost its regular physical projection.");
    rejects([&] { (void)source->projectSupportingSurface(face, {}); });
}
std::optional<double> cornerAngle(const BoundaryCorner &corner) {
    if (!corner.constraint)
        return std::nullopt;
    const auto &c = *corner.constraint;
    gp_Vec before(gp_Pnt(corner.position.x, corner.position.y, corner.position.z),
                  gp_Pnt(c.boundaryEnds[0].x, c.boundaryEnds[0].y, c.boundaryEnds[0].z));
    gp_Vec after(gp_Pnt(corner.position.x, corner.position.y, corner.position.z),
                 gp_Pnt(c.boundaryEnds[1].x, c.boundaryEnds[1].y, c.boundaryEnds[1].z));
    before.Normalize();
    after.Normalize();
    const auto cross = before.Crossed(after);
    const double minor = std::atan2(cross.Magnitude(), before.Dot(after));
    const gp_Vec normal(c.loopNormal.x, c.loopNormal.y, c.loopNormal.z);
    return normal.Dot(cross) <= 0 ? minor : 2 * std::numbers::pi - minor;
}
void physicalBoundaryCorners() {
    const auto box = geometry(BRepPrimAPI_MakeBox(.02, .03, .04).Shape());
    require(box->boundaryCorners().size() == 24, "Box physical trim corners are missing.");
    for (const auto &corner : box->boundaryCorners()) {
        require(cornerAngle(corner).has_value(), "Regular box corner is unknown.");
        expectNear(*cornerAngle(corner), std::numbers::pi / 2, 1e-12,
                   "A reversed box face changed its convex interior angle.");
    }
    BRepBuilderAPI_MakePolygon polygon;
    for (const auto p : {gp_Pnt(0, 0, 0), gp_Pnt(.02, 0, 0), gp_Pnt(.02, .01, 0),
                         gp_Pnt(.01, .01, 0), gp_Pnt(.01, .02, 0), gp_Pnt(0, .02, 0)})
        polygon.Add(p);
    polygon.Close();
    const auto concave = geometry(
        BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(polygon.Wire()).Face(), gp_Vec(0, 0, .01))
            .Shape());
    std::size_t reflex = 0;
    for (const auto &corner : concave->boundaryCorners())
        if (cornerAngle(corner) && *cornerAngle(corner) > std::numbers::pi) {
            ++reflex;
            expectNear(*cornerAngle(corner), 1.5 * std::numbers::pi, 1e-12,
                       "Reflex physical trim was reported as its convex complement.");
        }
    require(reflex == 2, "Extruded concave polygon lost its two reflex trim corners.");
    const auto hole = geometry(
        BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(.03, .03, .01).Shape(),
                        BRepPrimAPI_MakeBox(gp_Pnt(.01, .01, -.005), .01, .01, .02).Shape())
            .Shape());
    std::size_t holeReflex = 0;
    for (const auto &corner : hole->boundaryCorners())
        if (cornerAngle(corner) && *cornerAngle(corner) > std::numbers::pi) {
            ++holeReflex;
            expectNear(*cornerAngle(corner), 1.5 * std::numbers::pi, 1e-12,
                       "Clockwise inner trim loop was treated as an outer boundary.");
        }
    require(holeReflex == 8, "Through-hole prism lost its eight reflex face corners.");
    BRepBuilderAPI_MakePolygon wedge;
    const double angle = 12 * std::numbers::pi / 180;
    for (const auto p : {gp_Pnt(0, 0, 0), gp_Pnt(.02, 0, 0),
                         gp_Pnt(.02 * std::cos(angle), .02 * std::sin(angle), 0)})
        wedge.Add(p);
    wedge.Close();
    const auto narrow = geometry(
        BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(wedge.Wire()).Face(), gp_Vec(0, 0, .005))
            .Shape());
    std::size_t acuteLineCorners = 0;
    for (const auto &corner : narrow->boundaryCorners())
        if (const auto measured = cornerAngle(corner); measured && *measured < .3) {
            ++acuteLineCorners;
            expectNear(*measured, angle, 1e-12, "Exact planar line wedge angle changed.");
            require(distance(corner.constraint->boundaryEnds[0], corner.position) > .019 &&
                        distance(corner.constraint->boundaryEnds[1], corner.position) > .019,
                    "Constraint endpoints are short tangent samples instead of full trim ends.");
        }
    require(acuteLineCorners == 2,
            "Both cap corners of the narrow linear wedge must be evidenced.");
    for (const auto &source : {geometry(BRepPrimAPI_MakeCylinder(.01, .02).Shape()),
                               geometry(BRepPrimAPI_MakeSphere(.01).Shape())}) {
        require(source->boundaryCorners().empty(),
                "Periodic/nonplanar trims produced unsupported planar-line constraints.");
    }
    const auto path =
        std::filesystem::path(__FILE__).parent_path() / "data" / "spline_boundary_cusp.brep";
    TopoDS_Shape shape;
    BRep_Builder builder;
    require(BRepTools::Read(shape, path.string().c_str(), builder), "Cusp fixture is missing.");
    const auto cusp = geometry(shape);
    const Vec3 position{-.117357879186514, .41604817913815501, -.24851544106888598};
    for (const auto &corner : cusp->boundaryCorners())
        if (distance(corner.position, position) < 1e-10) {
            require(!corner.constraint,
                    "A curved tangent cusp was incorrectly made a hard chord-angle constraint.");
        }
}
void splineCornersAreNotLineConstraints() {
    TColgp_Array1OfPnt poles(1, 3);
    poles(1) = gp_Pnt(0, .02, 0);
    poles(2) = gp_Pnt(0, 0, 0);
    poles(3) = gp_Pnt(.02, 0, 0);
    TColStd_Array1OfReal knots(1, 3);
    TColStd_Array1OfInteger multiplicities(1, 3);
    for (int i = 1; i <= 3; ++i) {
        knots(i) = i - 1;
        multiplicities(i) = i == 2 ? 1 : 2;
    }
    const Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, knots, multiplicities, 1);
    BRepBuilderAPI_MakeWire wire;
    wire.Add(BRepBuilderAPI_MakeEdge(curve, 0, 1).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(curve, 1, 2).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(poles(3), poles(1)).Edge());
    const auto source = geometry(
        BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(wire.Wire()).Face(), gp_Vec(0, 0, .005))
            .Shape());
    for (const auto &corner : source->boundaryCorners())
        if (distance(corner.position, {}) < 1e-12) {
            require(!corner.constraint,
                    "C0 spline trims were heuristically treated as exact analytic lines.");
        }
}
void exactRemeshingIntegration() {
    const std::array<std::shared_ptr<const CadGeometry>, 2> sources{
        geometry(BRepPrimAPI_MakeBox(.012, .014, .018).Shape()),
        geometry(BRepPrimAPI_MakeCylinder(.008, .020).Shape())};
    for (const auto &source : sources) {
        MeshingOptions options;
        options.targetEdgeLength = .004;
        options.adaptive.geometricTolerance = .00015;
        options.geometry = source;
        const auto input = source->meshingSource(.00015);
        options.parameterization = input.parameters;
        const auto result = remeshIsotropic(input.mesh, options);
        require(result.quality.watertight() && result.quality.selfIntersectionPairs == 0 &&
                    result.quality.eulerCharacteristic == 2,
                "Exact CAD remeshing changed closure, genus or embedding.");
        require(result.quality.minimumAngleDegrees >= options.quality.minimumAngleDegrees &&
                    result.quality.maximumAspectRatio <= options.quality.maximumAspectRatio,
                "Exact CAD remeshing did not meet physical angle/aspect limits.");
        for (const auto &triangle : result.mesh.triangles)
            for (const auto index : triangle.vertices) {
                const auto p = result.mesh.vertices[index];
                const auto exact = source->projectSurface(triangle.faceId, p);
                require(distance(exact.position, p) < 1e-10,
                        "Remeshing vertices drifted from exact trimmed CAD geometry.");
            }
    }
}

void thinPlanarBoundaryIntegration() {
    // A retained thin physical face requires many refinement levels even
    // though the requested global size is millimetres. Batch count alone is
    // not evidence that the unchanged hard shape criteria are infeasible.
    constexpr double thickness = 11.5e-6;
    const auto source = geometry(BRepPrimAPI_MakeBox(.03, .02, thickness).Shape());
    MeshingOptions options{.008};
    options.geometry = source;
    const auto input = source->meshingSource(.00025);
    options.parameterization = input.parameters;
    const auto result = remeshIsotropic(input.mesh, options);
    const auto &quality = result.quality;
    require(quality.watertight() && quality.selfIntersectionPairs == 0 &&
                quality.eulerCharacteristic == 2 && quality.connectedComponents == 1 &&
                quality.minimumAngleDegrees >= 20 && quality.maximumAspectRatio <= 4,
            "Retained thin planar CAD faces did not meet unchanged closure/quality gates.");
    expectNear(quality.signedVolume, .03 * .02 * thickness, 1e-16,
               "Thin planar CAD remeshing changed its physical volume.");
    std::cout << "thin planar boundary: " << result.mesh.triangles.size() << " triangles\n";
}

void boundedAndNonsmoothCurves() {
    const auto sector = geometry(BRepPrimAPI_MakeCylinder(.01, .02, 3.14159265358979).Shape());
    std::uint32_t cylinder = 0;
    for (const auto &face : sector->faces())
        if (face.surfaceType == "Cylinder")
            cylinder = face.faceId;
    const std::array circle{cylinder, faceWithNormal(*sector, {0, 0, -1})};
    require(cylinder != 0, "Missing trimmed circular edge fixture.");
    require(distance(sector->projectBoundary(circle, {-.011, 0, 0}), {-.01, 0, 0}) < 1e-11,
            "A circular root outside the rounded trim displaced its true endpoint minimum.");
    rejects([&] { (void)sector->projectBoundary(circle, {0, 0, 0}); });

    TColgp_Array1OfPnt poles(1, 3);
    poles(1) = gp_Pnt(0, .02, 0);
    poles(2) = gp_Pnt(0, 0, 0);
    poles(3) = gp_Pnt(.02, 0, 0);
    TColStd_Array1OfReal knots(1, 3);
    TColStd_Array1OfInteger multiplicities(1, 3);
    for (int i = 1; i <= 3; ++i) {
        knots(i) = i - 1;
        multiplicities(i) = i == 2 ? 1 : 2;
    }
    const Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, knots, multiplicities, 1);
    const auto edge = BRepBuilderAPI_MakeEdge(curve).Edge();
    BRepBuilderAPI_MakeWire wire;
    wire.Add(edge);
    wire.Add(BRepBuilderAPI_MakeEdge(poles(3), poles(1)).Edge());
    const auto shape =
        BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(wire.Wire()).Face(), gp_Vec(0, 0, .005))
            .Shape();
    TopTools_IndexedMapOfShape faces;
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    std::vector<std::uint32_t> adjacent;
    for (int i = 1; i <= faces.Extent(); ++i) {
        TopTools_IndexedMapOfShape edges;
        TopExp::MapShapes(faces(i), TopAbs_EDGE, edges);
        if (edges.Contains(edge))
            adjacent.push_back(static_cast<std::uint32_t>(i));
    }
    require(adjacent.size() == 2, "Missing shared C0 spline edge in the closed prism.");
    const auto kink = geometry(shape, 1, 1);
    require(distance(kink->projectBoundary(adjacent, {-.001, -.001, 0}), {}) < 1e-11,
            "A trimmed spline projection skipped its interior C0 knot minimum.");
    require(distance(kink->projectBoundary(adjacent, {.003, -.001, 0}), {.003, 0, 0}) < 1e-11,
            "A spline projection missed the smooth interval next to its C0 knot.");
}

void finiteSplineBoundaryMetric(const std::shared_ptr<const CadGeometry> &source) {
    const Vec3 originalCenter{-.11271669846169176, .10488884203199959, -.22784981734704488};
    const auto face = std::min_element(
        source->faces().begin(), source->faces().end(), [&](const auto &a, const auto &b) {
            return distance(a.center, originalCenter) < distance(b.center, originalCenter);
        });
    require(distance(face->center, originalCenter) < 1e-10,
            "Missing original finite spline patch in the boundary fixture.");
    detail::SurfaceChart chart;
    chart.faceId = face->faceId;
    chart.geometry = source;
    detail::ChartReferenceTriangle reference;
    reference.uv = {{{0, 0}, {1, 0}, {0, 1}}};
    reference.positions = {{{-.11109415717771753, .10156723238895893, -.2265201316547234},
                            {-.10659322987461201, .099910285658156101, -.226347693540112},
                            {-.107993799039965, .098508800159434001, -.22662017181435601}}};
    reference.sizes = {{.004, .004, .004}};
    const auto normal = source->projectSurface(face->faceId, reference.positions[0]).normal;
    const auto &p = reference.positions;
    const gp_Vec e(p[1].x - p[0].x, p[1].y - p[0].y, p[1].z - p[0].z);
    const gp_Vec f(p[2].x - p[0].x, p[2].y - p[0].y, p[2].z - p[0].z);
    const bool reverseReference = e.Crossed(f).Dot(gp_Vec(normal.x, normal.y, normal.z)) < 0;
    if (reverseReference)
        std::swap(reference.positions[1], reference.positions[2]);
    chart.reference.push_back(reference);
    const detail::ChartPoint edgePoint = reverseReference
                                             ? detail::ChartPoint{.82020214020766302, 0}
                                             : detail::ChartPoint{0, .82020214020766302};
    const auto metric = chart.boundaryMetric(edgePoint);
    require(std::isfinite(metric.du.x) && std::isfinite(metric.dv.x) && metric.targetEdgeLength > 0,
            "Boundary tangent metric incorrectly required a full supporting position map.");
    bool unavailable = false;
    try {
        (void)chart.sample(edgePoint);
    } catch (const SurfaceProjectionDomainError &) {
        unavailable = true;
    }
    require(unavailable, "Finite spline search boundary was accepted as an interior map.");
}

void splineBoundaryProjection() {
    const auto path =
        std::filesystem::path(__FILE__).parent_path() / "data" / "spline_boundary_projection.brep";
    std::ifstream input(path);
    require(input.good(), "Missing closed spline-curve regression fixture.");
    BRep_Builder builder;
    TopoDS_Shape shape;
    BRepTools::Read(shape, input, builder);
    const auto source = geometry(shape, 1, 1);
    require(source->faces().size() == 8,
            "Spline-curve fixture is not the closed eight-face prism.");
    finiteSplineBoundaryMetric(source);
    const std::array ids{3u, 7u};
    const Vec3 point{-.11624951904256421, .11011661426716764, -.23131945739616308};
    const Vec3 expected{-.11625570518745129, .11011530834676361, -.23131782589426733};
    const auto projected = source->projectBoundary(ids, point);
    require(distance(projected, expected) < 1e-10,
            "Trimmed curve projection skipped its interior minimum and jumped to an endpoint.");
    const std::array reversed{7u, 3u};
    require(distance(source->projectBoundary(reversed, point), projected) < 1e-13,
            "Exact common curve depends on face order.");
    for (const Vec3 endpoint : {Vec3{-.11550000000000002, .10962365930961382, -.22994242926252309},
                                Vec3{-.11646228209305091, .11057951430393774, -.23249055358180609}})
        require(distance(source->projectBoundary(ids, endpoint), endpoint) < 1e-10,
                "Trimmed curve projection lost an endpoint minimum.");
}

void splineTopologicalJunction() {
    const auto path =
        std::filesystem::path(__FILE__).parent_path() / "data" / "spline_topological_junction.brep";
    std::ifstream input(path);
    require(input.good(), "Missing closed spline-junction regression fixture.");
    BRep_Builder builder;
    TopoDS_Shape shape;
    BRepTools::Read(shape, input, builder);
    const auto source = geometry(shape, .001, 1);
    require(source->faces().size() == 9,
            "Spline-junction fixture is not the closed nine-face cut.");
    const Vec3 point{-.068188134021083008, .38905774564574302, -.29075482811787301};
    const std::array ids{1u, 6u, 7u};
    require(distance(source->projectJunction(ids, point), point) < 1e-13,
            "Native face-vertex parameters were replaced by a different nearest-surface query.");
    const Vec3 nearby{point.x + 1e-8, point.y, point.z};
    require(distance(source->projectJunction(ids, nearby), point) < 1e-13,
            "Local spline-junction query moved the common topology vertex.");
    rejects([&] { (void)source->projectJunction(ids, {point.x + .001, point.y, point.z}); });
    rejects([&] {
        const std::array unrelated{1u, 2u, 6u};
        (void)source->projectJunction(unrelated, point);
    });
}

void trimmedSplineProjection() {
    const auto cylinder = BRepPrimAPI_MakeCylinder(.02, .03).Shape();
    const auto cutter =
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -.01), gp_Dir(0, 0, 1)), .004, .05).Shape();
    BRepAlgoAPI_Cut cut(cylinder, cutter);
    require(cut.IsDone(), "Could not construct trimmed spline fixture.");
    BRepBuilderAPI_NurbsConvert conversion(cut.Shape(), Standard_True);
    constexpr std::uint32_t firstFaceId = 101;
    const auto source = geometry(conversion.Shape(), 1, firstFaceId);
    std::uint32_t inner = 0, bottom = 0;
    for (const auto &face : source->faces()) {
        require(face.surfaceType == "BSplineSurface", "NURBS fixture retained analytic surfaces.");
        if (face.area < .001)
            inner = face.faceId;
        if (face.normal.z < -.9)
            bottom = face.faceId;
    }
    require(inner && bottom, "Missing trimmed spline fixture faces.");
    // These NURBS surfaces have genuinely finite natural parameter ranges.
    // Missing stationary roots beyond them must never become edge-clamped maps.
    rejects([&] { (void)source->projectSupportingSurface(inner, {.003, .001, .04}); });
    rejects([&] { (void)source->projectSupportingSurface(bottom, {.04, 0, -.003}); });
    TopTools_IndexedMapOfShape nativeFaces;
    TopExp::MapShapes(conversion.Shape(), TopAbs_FACE, nativeFaces);
    const double boundaryTolerance =
        BRep_Tool::Tolerance(TopoDS::Face(nativeFaces(static_cast<int>(inner - firstFaceId + 1))));
    std::vector<std::pair<Vec3, SurfaceProjection>> samples;
    for (unsigned i = 0; i < 32; ++i) {
        const double phi = (i + .25) * 2 * std::numbers::pi / 32;
        for (double z : {-.003, 0., .015, .03, .033}) {
            const Vec3 point{.003 * std::cos(phi), .003 * std::sin(phi), z};
            const auto projected = source->projectSurface(inner, point);
            samples.emplace_back(point, projected);
            // The trimmed-boundary extrema operation resolves the native model
            // tolerance. For a circle, its normal error is bounded by position
            // error / radius. Keep the existing tighter interior reference test.
            const bool interior = z > 0 && z < .03;
            const double positionBudget = interior ? 1e-10 : boundaryTolerance;
            const double normalBudget = interior ? 1e-9 : boundaryTolerance / .004;
            require(distance(projected.position, {.004 * std::cos(phi), .004 * std::sin(phi),
                                                  std::clamp(z, 0., .03)}) < positionBudget,
                    "Spline projection did not attain the nearest trimmed analytic point.");
            require(distance(projected.normal, {-std::cos(phi), -std::sin(phi), 0}) < normalBudget,
                    "Spline interior/boundary projection changed the oriented normal.");
        }
    }
    const auto verify = [&](std::size_t index) {
        const auto &[point, expected] = samples[index];
        const auto actual = source->projectSurface(inner, point);
        require(distance(actual.position, expected.position) < 1e-13 &&
                    distance(actual.normal, expected.normal) < 1e-13,
                "NURBS projection depends on preceding points or concurrent callers.");
    };
    for (std::size_t i = samples.size(); i-- > 0;)
        verify(i);
    std::barrier start(4);
    std::vector<std::future<void>> workers;
    for (unsigned i = 0; i < 4; ++i)
        workers.push_back(std::async(std::launch::async, [&, i] {
            start.arrive_and_wait();
            // An exception must return the exclusive search state before a
            // subsequent valid projection, including trimmed-boundary cases.
            rejects([&] {
                (void)source->projectSurface(inner,
                                             {std::numeric_limits<double>::quiet_NaN(), 0, 0});
            });
            for (std::size_t j = 0; j < samples.size(); ++j)
                verify((j * 17 + i) % samples.size());
        }));
    for (auto &worker : workers)
        worker.get();
    const auto interior = source->projectSurface(bottom, {.013, .001, -.003});
    require(distance(interior.position, {.013, .001, 0}) < 1e-10,
            "Trim classification rejected the spline face interior.");
    const auto hole = source->projectSurface(bottom, {.002, .001, -.003});
    expectNear(std::hypot(hole.position.x, hole.position.y), .004, 1e-10,
               "Spline interior projection ignored its trimmed hole.");
    expectNear(hole.position.z, 0, 1e-10, "Spline trim projection left the cap plane.");
    require(distance(hole.normal, {0, 0, -1}) < 1e-12,
            "Trimmed spline hole used a normal from a different surface location.");
}

void trimmedSplineBoundaryNormal() {
    const auto path =
        std::filesystem::path(__FILE__).parent_path() / "data" / "trimmed_spline_boundary.brep";
    std::ifstream input(path);
    require(input.good(), "Missing trimmed spline boundary regression fixture.");
    BRep_Builder builder;
    TopoDS_Shape shape;
    BRepTools::Read(shape, input, builder);
    require(!shape.IsNull(), "Could not read trimmed spline boundary fixture.");
    const auto source = geometry(shape, 1, 1);
    require(source->faces().size() == 6, "Boundary fixture is not the minimal closed prism.");
    const Vec3 boundaryPoint{-.045261334628069562, .33703602764109225, -.18250201538169269};
    const Vec3 outwardNormal{-.60536021210321656, .75685164860694865, -.24640331897783488};
    const auto projected = source->projectSurface(5, boundaryPoint);
    require(distance(projected.position, boundaryPoint) < 1e-11,
            "Trimmed spline boundary projection moved off its topological edge.");
    require(distance(projected.normal, outwardNormal) < 1e-9,
            "Boundary normal was evaluated at a remote untrimmed surface extremum.");
}

void splineEdgeParameterCorrection() {
    const auto path =
        std::filesystem::path(__FILE__).parent_path() / "data" / "spline_edge_parameter.brep";
    std::ifstream input(path);
    require(input.good(), "Missing spline edge parameter regression fixture.");
    BRep_Builder builder;
    TopoDS_Shape shape;
    BRepTools::Read(shape, input, builder);
    const auto source = geometry(shape, 1, 1);
    require(source->faces().size() == 27, "Spline edge fixture is not the closed 27-face prism.");
    const Vec3 boundaryPoint{-.1109621569458225, .48625926604784769, -.19731814235088041};
    const Vec3 outwardNormal{-.69500275533306066, .65446603873545295, .2977337304058415};
    TopTools_IndexedMapOfShape faces;
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    const auto face = TopoDS::Face(faces(26));
    const gp_Pnt p(boundaryPoint.x, boundaryPoint.y, boundaryPoint.z);
    BRepExtrema_DistShapeShape projection(BRepBuilderAPI_MakeVertex(p).Vertex(), face);
    require(projection.IsDone() && projection.NbSolution() == 1 &&
                projection.SupportTypeShape2(1) == BRepExtrema_IsOnEdge,
            "Spline regression no longer exercises a topological edge projection.");
    const auto edge = TopoDS::Edge(projection.SupportOnShape2(1));
    double parameter, first, last;
    projection.ParOnEdgeS2(1, parameter);
    const auto curve = BRep_Tool::CurveOnSurface(edge, face, first, last);
    const auto uv = curve->Value(parameter);
    const auto direct = BRepAdaptor_Surface(face, Standard_True).Value(uv.X(), uv.Y());
    const double tolerance = std::max(BRep_Tool::Tolerance(face), BRep_Tool::Tolerance(edge));
    require(BRep_Tool::SameParameter(edge) && direct.Distance(p) > tolerance,
            "Spline fixture lost its inaccurate SameParameter approximation.");
    const auto projected = source->projectSurface(26, boundaryPoint);
    const auto onFace =
        BRepAdaptor_Surface(face, Standard_True).Value(.0089264564387950187, .90817324696316182);
    require(distance(projected.position, {onFace.X(), onFace.Y(), onFace.Z()}) < 1e-11 &&
                distance(projected.position, boundaryPoint) < tolerance,
            "Local pcurve inversion did not return the actual trimmed-face point.");
    // The independent 3D curve remains the shared-boundary contract. A face
    // projection must not replace its exact position with a different face's
    // pcurve approximation.
    std::uint32_t adjacent = 0;
    for (int i = 1; i <= faces.Extent(); ++i) {
        TopTools_IndexedMapOfShape edges;
        TopExp::MapShapes(faces(i), TopAbs_EDGE, edges);
        if (i != 26 && edges.Contains(edge))
            adjacent = static_cast<std::uint32_t>(i);
    }
    require(adjacent != 0, "Spline edge fixture lost its adjacent physical face.");
    const auto shared =
        source->projectBoundary(std::array<std::uint32_t, 2>{26, adjacent}, boundaryPoint);
    require(distance(shared, boundaryPoint) < 1e-11,
            "Surface projection changed the separate shared-boundary geometry.");
    require(distance(projected.normal, outwardNormal) < 1e-9,
            "Local pcurve inversion did not use the corresponding surface normal.");
    const Vec3 inconsistentEdgePoint{-.068187871299755828, .49552697964842418, -.2917703551354292};
    const gp_Pnt otherPoint(inconsistentEdgePoint.x, inconsistentEdgePoint.y,
                            inconsistentEdgePoint.z);
    BRepExtrema_DistShapeShape otherProjection(BRepBuilderAPI_MakeVertex(otherPoint).Vertex(),
                                               face);
    require(otherProjection.IsDone() && otherProjection.NbSolution() == 1 &&
                otherProjection.SupportTypeShape2(1) == BRepExtrema_IsOnEdge,
            "Spline residual regression lost its edge support.");
    const auto otherEdge = TopoDS::Edge(otherProjection.SupportOnShape2(1));
    const auto otherCurve = BRep_Tool::CurveOnSurface(otherEdge, face, first, last);
    const auto otherUv = otherCurve->Value(.1250000000042219);
    const auto otherFacePoint =
        BRepAdaptor_Surface(face, Standard_True).Value(otherUv.X(), otherUv.Y());
    const auto otherResult = source->projectSurface(26, inconsistentEdgePoint);
    require(distance(otherResult.position,
                     {otherFacePoint.X(), otherFacePoint.Y(), otherFacePoint.Z()}) < 1e-10,
            "Spline projection returned the 3D curve instead of the actual face point.");
    require(distance(otherResult.position, inconsistentEdgePoint) >
                std::max(BRep_Tool::Tolerance(face), BRep_Tool::Tolerance(otherEdge)),
            "Spline residual regression no longer distinguishes curve and surface geometry.");
}

void invalidSources() {
    const auto box = BRepPrimAPI_MakeBox(.02, .03, .04).Shape();
    auto part = artifact(box, 1);
    rejects([&] { (void)CadGeometry::fromBrepParts({}, 1); });
    rejects([&] { (void)CadGeometry::fromBrepParts(std::span{&part, 1}, 0); });
    auto broken = part;
    broken.faceIds.pop_back();
    rejects([&] { (void)CadGeometry::fromBrepParts(std::span{&broken, 1}, 1); });
    broken = part;
    broken.faceIds[1] = broken.faceIds[0];
    rejects([&] { (void)CadGeometry::fromBrepParts(std::span{&broken, 1}, 1); });
    broken = part;
    broken.data = "not a BRep";
    rejects([&] { (void)CadGeometry::fromBrepParts(std::span{&broken, 1}, 1); });
    rejects([&] { (void)geometry(box.Reversed()); });
    TopExp_Explorer face(box, TopAbs_FACE);
    rejects([&] { (void)geometry(face.Current()); });
}
void stepRoundtrip() {
    const auto directory =
        std::filesystem::temp_directory_path() /
        ("hornsculpt-cad-" +
         std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    require(std::filesystem::create_directory(directory),
            "Cannot create isolated STEP test directory.");
    struct Cleanup {
        std::filesystem::path path;
        ~Cleanup() {
            std::error_code error;
            std::filesystem::remove_all(path, error);
        }
    } cleanup{directory};
    const auto path = (directory / "box-mm.step").string();
    STEPControl_Writer writer;
    writer.Model()->SetLocalLengthUnit(1.0);
    writer.Model()->SetWriteLengthUnit(1.0);
    DESTEP_Parameters parameters;
    parameters.WriteUnit = UnitsMethods_LengthUnit_Millimeter;
    require(writer.Transfer(BRepPrimAPI_MakeBox(20, 30, 40).Shape(), STEPControl_AsIs,
                            parameters) == IFSelect_RetDone,
            "Real STEP fixture transfer failed.");
    require(writer.Write(path.c_str()) == IFSelect_RetDone, "Real STEP fixture write failed.");
    const auto imported = CadGeometry::fromStep(path);
    double area = 0;
    for (const auto &f : imported->faces())
        area += f.area;
    expectNear(area, .0052, 1e-10, "STEP units did not transfer into metres.");
    closedMesh(*imported, imported->previewMesh(.00025), 2, .00025);
    const auto offline = CadGeometry::fromBrepParts(imported->brepParts(), 1);
    expectNear(offline->faces().front().area, imported->faces().front().area, 1e-12,
               "STEP canonical offline BRep changed SI measures.");
    rejects([&] { (void)CadGeometry::fromStep((directory / "missing.step").string()); });
}
int meshStep(const std::string &path, const std::string &targetText,
             MeshQualityLimits limits = {}) {
    std::size_t consumed = 0;
    const double target = std::stod(targetText, &consumed);
    if (consumed != targetText.size() || !std::isfinite(target) || target <= 0)
        throw std::invalid_argument(
            "Target edge length must be a finite positive number in metres.");
    std::cout << std::unitbuf;
    const auto begin = std::chrono::steady_clock::now();
    const auto seconds = [](auto since) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - since).count();
    };
    const auto source = CadGeometry::fromStep(path);
    std::cout << "STEP faces=" << source->faces().size() << " importSeconds=" << seconds(begin)
              << '\n';
    const auto previewBegin = std::chrono::steady_clock::now();
    const auto input = source->meshingSource(.00025);
    const auto &boundary = input.mesh;
    const auto sourceQuality = inspectMesh(boundary);
    std::cout << "Source triangles=" << boundary.triangles.size()
              << " closed=" << sourceQuality.watertight()
              << " Euler=" << sourceQuality.eulerCharacteristic
              << " components=" << sourceQuality.connectedComponents
              << " previewSeconds=" << seconds(previewBegin) << '\n';
    MeshingOptions options;
    options.targetEdgeLength = target;
    options.quality = limits;
    options.adaptive.geometricTolerance = .00025;
    options.geometry = source;
    options.parameterization = input.parameters;
    std::cout << "Limits angle=" << limits.minimumAngleDegrees
              << " aspect=" << limits.maximumAspectRatio << " enabled=" << limits.enabled << '\n';
    const auto meshBegin = std::chrono::steady_clock::now();
    const auto result = [&] {
        try {
            return remeshIsotropic(boundary, options);
        } catch (...) {
            std::cout << "Failed meshSeconds=" << seconds(meshBegin)
                      << " totalSeconds=" << seconds(begin) << '\n';
            throw;
        }
    }();
    const auto &quality = result.quality;
    std::set<std::uint32_t> faces;
    for (const auto &triangle : result.mesh.triangles)
        faces.insert(triangle.faceId);
    std::cout << "Output triangles=" << result.mesh.triangles.size()
              << " minAngle=" << quality.minimumAngleDegrees
              << " aspect=" << quality.maximumAspectRatio << " closed=" << quality.watertight()
              << " Euler=" << quality.eulerCharacteristic
              << " components=" << quality.connectedComponents << " faces=" << faces.size()
              << " selfIntersections=" << quality.selfIntersectionPairs
              << " meshSeconds=" << seconds(meshBegin) << " totalSeconds=" << seconds(begin)
              << '\n';
    require(quality.watertight() && quality.selfIntersectionPairs == 0 &&
                quality.eulerCharacteristic == sourceQuality.eulerCharacteristic &&
                quality.connectedComponents == sourceQuality.connectedComponents &&
                faces.size() == source->faces().size() &&
                (!limits.enabled || (quality.minimumAngleDegrees >= limits.minimumAngleDegrees &&
                                     quality.maximumAspectRatio <= limits.maximumAspectRatio)),
            "Exact STEP remeshing failed physical topology/quality gates.");
    return 0;
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc > 1) {
            if (argc == 2 && std::string(argv[1]) == "--trim-contracts-only") {
                scaledNotchedPrismProjection();
                std::cout << "Scaled notched trim PASS.\n" << std::flush;
                curvedTrimCompetesWithStationaryMinimum();
                std::cout << "Trim boundary versus stationary minimum PASS.\n";
                return 0;
            }
            if (argc == 6 && std::string(argv[1]) == "--mesh-step")
                return meshStep(argv[2], argv[3], {true, std::stod(argv[4]), std::stod(argv[5])});
            if (argc == 5 && std::string(argv[1]) == "--mesh-step" &&
                std::string(argv[4]) == "--no-shape-limits")
                return meshStep(argv[2], argv[3], {false, 20, 4});
            if (argc == 2 && std::string(argv[1]) == "--thin-only") {
                thinPlanarBoundaryIntegration();
                return 0;
            }
            if (argc == 2 && std::string(argv[1]) == "--caps-only") {
                sphericalCapParameters();
                return 0;
            }
            if (argc == 2 && std::string(argv[1]) == "--revolved-only") {
                revolvedSplineNearestMinimum();
                return 0;
            }
            if (argc != 4 || std::string(argv[1]) != "--mesh-step")
                throw std::invalid_argument(
                    "Usage: cad_geometry_tests --mesh-step <path> <target-metres> "
                    "[<minimum-angle-degrees> <maximum-aspect> | --no-shape-limits]");
            return meshStep(argv[2], argv[3]);
        }
        boxUnitsAndProjection();
        scaledNotchedPrismProjection();
        curvedTrimCompetesWithStationaryMinimum();
        tolerantTopologicalJunction();
        splineTopologicalJunction();
        splineBoundaryProjection();
        boundedAndNonsmoothCurves();
        cylinderAndTrimmedBore();
        thinCurvedRimTessellation();
        supportingSpherePole();
        sphericalCapParameters();
        revolvedSplineNearestMinimum();
        physicalBoundaryCorners();
        splineCornersAreNotLineConstraints();
        trimmedSplineProjection();
        trimmedSplineBoundaryNormal();
        splineEdgeParameterCorrection();
        exactRemeshingIntegration();
        thinPlanarBoundaryIntegration();
        invalidSources();
        stepRoundtrip();
        std::cout << "Exact CAD: BRep/STEP SI units, trimmed faces/holes, oriented normals, "
                     "topology-welded meshes, concurrent projection and rejection PASS.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
