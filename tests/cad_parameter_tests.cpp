#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_Surface.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <future>
#include <hornsculpt/cad_geometry.hpp>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class F> void rejects(F f) {
    bool failed = false;
    try {
        f();
    } catch (const std::exception &) {
        failed = true;
    }
    require(failed, "Invalid native parameter query was accepted.");
}
double distance(Vec3 a, Vec3 b) {
    return std::hypot(a.x - b.x, a.y - b.y, a.z - b.z);
}
Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
std::shared_ptr<const CadGeometry> geometry(const TopoDS_Shape &shape) {
    std::ostringstream stream;
    BRepTools::Write(shape, stream, Standard_False, Standard_False, TopTools_FormatVersion_CURRENT);
    TopTools_IndexedMapOfShape faces;
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    std::vector<std::uint32_t> ids(faces.Extent());
    std::iota(ids.begin(), ids.end(), 101u);
    const CadBrepPart part{stream.str(), ids};
    return CadGeometry::fromBrepParts(std::span{&part, 1}, 1);
}
void verify(const MeshingSource &source) {
    require(source.parameters && source.parameters->sourceMesh() == source.mesh,
            "Parameter provider does not own the exact returned source snapshot.");
    auto records = source.parameters->triangles();
    require(records.size() == source.mesh.triangles.size(),
            "Native triangle records lost alignment during weld.");
    std::size_t intervals = 0;
    for (std::size_t i = 0; i < records.size(); ++i) {
        const auto &t = source.mesh.triangles[i];
        const auto &r = records[i];
        auto domain = source.parameters->domain(t.faceId);
        require(domain.orientation == 1 || domain.orientation == -1,
                "Missing native source orientation.");
        ParameterPoint center{};
        for (unsigned j = 0; j < 3; ++j) {
            center.u += r.uv[j].u / 3;
            center.v += r.uv[j].v / 3;
            const auto p = source.parameters->surface(t.faceId, r.uv[j]);
            require(distance(p.position, source.mesh.vertices[t.vertices[j]]) < 2e-12,
                    "Native corner UV was replaced by a welded or projected coordinate.");
            if (auto edge = r.edges[j]) {
                ++intervals;
                require(edge->curveId && edge->useId && edge->t0 != edge->t1,
                        "Invalid native curve-use interval.");
                for (auto [parameter, corner] :
                     {std::pair{edge->t0, j}, std::pair{edge->t1, (j + 1) % 3}}) {
                    auto curve = source.parameters->curve(edge->curveId, parameter);
                    auto use = source.parameters->curveUse(edge->useId, parameter);
                    require(distance(curve.position, source.mesh.vertices[t.vertices[corner]]) <
                                2e-12,
                            "Canonical native curve parameter lost its welded source endpoint.");
                    if (std::hypot(use.uv.u - r.uv[corner].u, use.uv.v - r.uv[corner].v) >= 1e-9)
                        throw std::runtime_error(
                            "Native periodic/face use UV mismatch: face=" +
                            std::to_string(t.faceId) + " curve=" + std::to_string(edge->curveId) +
                            " use=" + std::to_string(edge->useId) +
                            " uv=" + std::to_string(use.uv.u) + "," + std::to_string(use.uv.v) +
                            " corner=" + std::to_string(r.uv[corner].u) + "," +
                            std::to_string(r.uv[corner].v) + " t=" + std::to_string(parameter));
                    require(use.modelToleranceMetres >= 0 &&
                                std::isfinite(use.modelToleranceMetres),
                            "Invalid native modelling tolerance.");
                }
            }
        }
        const auto s = source.parameters->surface(t.faceId, center);
        auto n = cross(
            subtract(source.mesh.vertices[t.vertices[1]], source.mesh.vertices[t.vertices[0]]),
            subtract(source.mesh.vertices[t.vertices[2]], source.mesh.vertices[t.vertices[0]]));
        require(domain.orientation * dot(n, cross(s.du, s.dv)) > 0,
                "Native D1 derivatives were flipped or source orientation was lost.");
    }
    require(intervals > 0, "Source has no explicit native boundary uses.");
}
void derivatives(const MeshingSource &source) {
    auto records = source.parameters->triangles();
    for (std::size_t i = 0; i < std::min<std::size_t>(records.size(), 80); ++i) {
        const auto id = source.mesh.triangles[i].faceId;
        ParameterPoint uv{};
        for (const auto &corner : records[i].uv) {
            uv.u += corner.u / 3;
            uv.v += corner.v / 3;
        }
        const auto s = source.parameters->surface(id, uv);
        const double h = 1e-7;
        const auto u = source.parameters->surface(id, {uv.u + h, uv.v}).position;
        const auto v = source.parameters->surface(id, {uv.u, uv.v + h}).position;
        require(
            distance({(u.x - s.position.x) / h, (u.y - s.position.y) / h, (u.z - s.position.z) / h},
                     s.du) < 1e-6,
            "Native surface du is not dS/du.");
        require(
            distance({(v.x - s.position.x) / h, (v.y - s.position.y) / h, (v.z - s.position.z) / h},
                     s.dv) < 1e-6,
            "Native surface dv is not dS/dv.");
        for (const auto &e : records[i].edges)
            if (e) {
                double t = (e->t0 + e->t1) / 2, step = std::abs(e->t1 - e->t0) * 1e-5;
                auto c = source.parameters->curve(e->curveId, t),
                     p = source.parameters->curve(e->curveId, t + step),
                     m = source.parameters->curve(e->curveId, t - step);
                require(distance({(p.position.x - m.position.x) / (2 * step),
                                  (p.position.y - m.position.y) / (2 * step),
                                  (p.position.z - m.position.z) / (2 * step)},
                                 c.tangent) < 1e-6,
                        "Native spatial curve tangent is not dC/dt.");
                auto use = source.parameters->curveUse(e->useId, t);
                require(!use.vertexToleranceMetres,
                        "A nonendpoint curve use inherited a topological vertex tolerance.");
                auto a = source.parameters->curveUse(e->useId, t + step),
                     b = source.parameters->curveUse(e->useId, t - step);
                require(std::hypot((a.uv.u - b.uv.u) / (2 * step) - use.uvDerivative.u,
                                   (a.uv.v - b.uv.v) / (2 * step) - use.uvDerivative.v) < 1e-6,
                        "Native pcurve tangent is not dgamma/dt.");
            }
    }
}
void nativeSnapshots() {
    for (auto shape :
         {BRepPrimAPI_MakeBox(.02, .03, .04).Shape(), BRepPrimAPI_MakeCylinder(.01, .03).Shape(),
          BRepPrimAPI_MakeSphere(.01).Shape()}) {
        auto cad = geometry(shape);
        auto source = cad->meshingSource(.00025);
        verify(source);
        derivatives(source);
        require(cad->previewMesh(.00025) == source.mesh,
                "previewMesh no longer delegates to the same native source path.");
        const auto immutable = source.parameters->sourceMesh();
        source.mesh.triangles.front().boundary = BoundaryCondition::driven;
        require(source.parameters->sourceMesh() == immutable,
                "Caller BC assignment mutated the immutable provider.");
        cad.reset();
        auto id = immutable.triangles.front().faceId;
        auto uv = source.parameters->triangles().front().uv[0];
        const auto expected = source.parameters->surface(id, uv);
        for (const auto &triangle : source.mesh.triangles) {
            const auto domain = source.parameters->domain(triangle.faceId);
            if (domain.vBounds && !domain.vPeriod) {
                rejects([&] {
                    (void)source.parameters->surface(triangle.faceId,
                                                     {0, (*domain.vBounds)[1] + 1});
                });
                break;
            }
        }
        std::array<std::future<void>, 4> workers;
        for (auto &w : workers)
            w = std::async(std::launch::async, [&] {
                for (unsigned n = 0; n < 128; ++n)
                    require(distance(source.parameters->surface(id, uv).position,
                                     expected.position) == 0,
                            "Concurrent immutable native evaluation changed.");
            });
        for (auto &w : workers)
            w.get();
        rejects([&] { (void)source.parameters->surface(0, {}); });
        rejects([&] { (void)source.parameters->curve(0, 0); });
        rejects([&] { (void)source.parameters->curveUse(0, 0); });
        rejects([&] {
            (void)source.parameters->surface(id, {std::numeric_limits<double>::quiet_NaN(), 0});
        });
    }
}

void explicitTrimEvidence() {
    const auto cad = geometry(BRepPrimAPI_MakeBox(.02, .03, .04).Shape());
    const auto source = cad->meshingSource(.00025);
    const auto &record = source.parameters->triangles().front();
    const auto id = source.mesh.triangles.front().faceId;
    ParameterPoint center{};
    for (const auto uv : record.uv) {
        center.u += uv.u / 3;
        center.v += uv.v / 3;
    }
    require(source.parameters->isInsideTrim(id, center),
            "Interior native point has no strict trim evidence.");
    require(!source.parameters->isInsideTrim(id, record.uv.front()),
            "A physical corner was incorrectly certified strictly inside the trim.");
    const ParameterPoint outside{center.u + 1, center.v + 1};
    const auto sample = source.parameters->surface(id, outside);
    require(std::isfinite(sample.position.x) && !source.parameters->isInsideTrim(id, outside),
            "Supporting-surface evaluation was confused with physical-trim membership.");
    std::array<std::future<void>, 4> workers;
    for (auto &worker : workers)
        worker = std::async(std::launch::async, [&] {
            for (unsigned i = 0; i < 16; ++i)
                require(source.parameters->isInsideTrim(id, center) &&
                            !source.parameters->isInsideTrim(id, outside),
                        "Concurrent trim classification changed its immutable snapshot.");
        });
    for (auto &worker : workers)
        worker.get();
    rejects([&] { (void)source.parameters->isInsideTrim(0, center); });
    rejects([&] {
        (void)source.parameters->isInsideTrim(id, {std::numeric_limits<double>::quiet_NaN(), 0});
    });
}
void periodicUses() {
    const auto cad = geometry(BRepPrimAPI_MakeCylinder(.01, .03).Shape());
    const auto source = cad->meshingSource(.00025);
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::map<std::uint32_t, double>> uses;
    for (std::size_t i = 0; i < source.mesh.triangles.size(); ++i)
        for (const auto &e : source.parameters->triangles()[i].edges)
            if (e)
                uses[{source.mesh.triangles[i].faceId, e->curveId}].emplace(e->useId,
                                                                            (e->t0 + e->t1) / 2);
    std::size_t seams = 0;
    for (auto &[key, list] : uses)
        if (list.size() == 2) {
            ++seams;
            auto first = list.begin(), second = std::next(first);
            const double t = first->second;
            auto a = source.parameters->curveUse(first->first, t),
                 b = source.parameters->curveUse(second->first, t);
            require(std::abs(std::abs(a.uv.u - b.uv.u) - 2 * std::numbers::pi) < 1e-12,
                    "Periodic seam lost its two distinct UV lifts.");
            require(std::abs(a.uv.v - b.uv.v) < 1e-12,
                    "Periodic seam uses disagree on native height.");
            const auto domain = source.parameters->domain(key.first);
            const auto native = source.parameters->surface(key.first, a.uv);
            const auto lifted =
                source.parameters->surface(key.first, {a.uv.u + 2 * std::numbers::pi, a.uv.v});
            require(distance(native.position, lifted.position) < 1e-12 &&
                        distance(native.du, lifted.du) < 1e-12,
                    "A genuine periodic UV lift changed native S/D1.");
            require(domain.uPeriod && std::abs(*domain.uPeriod - 2 * std::numbers::pi) < 1e-12,
                    "Native period is missing.");
        }
    require(seams == 1,
            "Cylinder periodic seam was not recorded as two uses of one spatial curve.");
    for (const auto &[key, list] : uses) {
        const auto domain = source.parameters->domain(key.first);
        if (domain.vBounds && !domain.vPeriod)
            rejects([&] {
                (void)source.parameters->surface(key.first, {0, (*domain.vBounds)[1] + 1});
            });
    }
    auto fine = cad->meshingSource(.00005);
    require(fine.mesh.triangles.size() > source.mesh.triangles.size(),
            "A previous preview cache replaced the requested fresh detail snapshot.");
    std::stop_source stop;
    stop.request_stop();
    rejects([&] { (void)cad->meshingSource(.00025, stop.get_token()); });
}
void shiftedPeriodicSourceDomain() {
    auto shape = BRepPrimAPI_MakeCylinder(.01, .03).Shape();
    // The same cylinder can legitimately carry a source trim over [2pi, 4pi].
    // Shift pcurves only: the surface, solid, spatial edges and tolerances stay fixed.
    for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next()) {
        const auto face = TopoDS::Face(it.Current());
        if (!BRep_Tool::Surface(face)->IsUPeriodic())
            continue;
        TopTools_IndexedMapOfShape edges;
        TopExp::MapShapes(face, TopAbs_EDGE, edges);
        BRep_Builder builder;
        const auto shifted = [&](const TopoDS_Edge &edge) {
            double first, last;
            const auto curve = BRep_Tool::CurveOnSurface(edge, face, first, last);
            auto result = Handle(Geom2d_Curve)::DownCast(curve->Copy());
            result->Translate(gp_Vec2d(2 * std::numbers::pi, 0));
            return result;
        };
        for (int i = 1; i <= edges.Extent(); ++i) {
            const auto edge = TopoDS::Edge(edges(i));
            auto first = shifted(edge);
            if (BRep_Tool::IsClosed(edge, face)) {
                auto second = shifted(TopoDS::Edge(edge.Reversed()));
                if (edge.Orientation() == TopAbs_REVERSED)
                    std::swap(first, second);
                builder.UpdateEdge(edge, first, second, face, BRep_Tool::Tolerance(edge));
            } else {
                builder.UpdateEdge(edge, first, face, BRep_Tool::Tolerance(edge));
            }
        }
        BRepTools::Update(face);
    }
    const auto source = geometry(shape)->meshingSource(.00025);
    verify(source);
    derivatives(source);
    bool beyondBasis = false;
    for (std::size_t i = 0; i < source.mesh.triangles.size(); ++i) {
        const auto id = source.mesh.triangles[i].faceId;
        if (!source.parameters->domain(id).uPeriod)
            continue;
        for (const auto uv : source.parameters->triangles()[i].uv) {
            if (uv.u <= 2 * std::numbers::pi + .1)
                continue;
            beyondBasis = true;
            const auto sample = source.parameters->surface(id, uv);
            const auto canonical =
                source.parameters->surface(id, {uv.u - 2 * std::numbers::pi, uv.v});
            require(distance(sample.position, canonical.position) < 1e-12 &&
                        distance(sample.du, canonical.du) < 1e-12 &&
                        distance(sample.dv, canonical.dv) < 1e-12,
                    "A periodic source-domain lift was mistaken for spline extrapolation.");
        }
    }
    require(beyondBasis, "Synthetic periodic source did not extend past the basis interval.");
}
void endpointTolerances() {
    auto shape = BRepPrimAPI_MakeBox(.02, .03, .04).Shape();
    TopTools_IndexedMapOfShape vertices;
    TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
    BRep_Builder builder;
    builder.UpdateVertex(TopoDS::Vertex(vertices(1)), 1e-6);
    builder.UpdateVertex(TopoDS::Vertex(vertices(2)), 2e-6);
    const auto source = geometry(shape)->meshingSource(.00025);
    std::set<double> seen;
    for (const auto &record : source.parameters->triangles())
        for (const auto &e : record.edges)
            if (e)
                for (double t : {e->t0, e->t1}) {
                    const auto use = source.parameters->curveUse(e->useId, t);
                    if (!use.vertexToleranceMetres)
                        continue;
                    const auto p = source.parameters->curve(e->curveId, t).position;
                    bool matched = false;
                    for (int i = 1; i <= vertices.Extent(); ++i) {
                        auto v = TopoDS::Vertex(vertices(i));
                        auto native = BRep_Tool::Pnt(v);
                        if (distance(p, {native.X(), native.Y(), native.Z()}) < 1e-12) {
                            require(*use.vertexToleranceMetres == BRep_Tool::Tolerance(v),
                                    "Canonical endpoint tolerance was swapped by edge-use "
                                    "orientation.");
                            matched = true;
                        }
                    }
                    require(matched,
                            "Endpoint metadata did not name a real original topological vertex.");
                    if (*use.vertexToleranceMetres > use.modelToleranceMetres)
                        seen.insert(*use.vertexToleranceMetres);
                    auto interior = source.parameters->curveUse(e->useId, (e->t0 + e->t1) / 2);
                    require(!interior.vertexToleranceMetres,
                            "Native vertex tolerance leaked into a curve interior.");
                }
    require(seen.contains(1e-6) && seen.contains(2e-6),
            "Distinct native vertex tolerances were lost or merged into edge/face tolerance.");
}

void rawSplineUse() {
    const auto path =
        std::filesystem::path(__FILE__).parent_path() / "data" / "spline_edge_parameter.brep";
    TopoDS_Shape shape;
    BRep_Builder builder;
    require(BRepTools::Read(shape, path.string().c_str(), builder),
            "Native pcurve fixture missing.");
    auto cad = geometry(shape);
    auto source = cad->meshingSource(.00025);
    require(source.parameters->sourceMesh() == source.mesh,
            "Spline parameter snapshot lost mesh identity.");
    constexpr std::uint32_t face = 126;
    const auto domain = source.parameters->domain(face);
    require(domain.vBounds.has_value(), "Native spline supporting domain is missing.");
    rejects([&] { (void)source.parameters->surface(face, {.5, (*domain.vBounds)[1] + 1e-4}); });
    constexpr double parameter = .75385466356773623;
    const Vec3 expected{-.1109621569458225, .48625926604784769, -.19731814235088041};
    std::set<std::uint32_t> matches;
    for (std::size_t i = 0; i < source.mesh.triangles.size(); ++i)
        if (source.mesh.triangles[i].faceId == face)
            for (const auto &e : source.parameters->triangles()[i].edges)
                if (e && parameter >= std::min(e->t0, e->t1) &&
                    parameter <= std::max(e->t0, e->t1)) {
                    auto curve = source.parameters->curve(e->curveId, parameter);
                    if (distance(curve.position, expected) > 1e-9)
                        continue;
                    matches.insert(e->useId);
                    auto use = source.parameters->curveUse(e->useId, parameter);
                    auto surface = source.parameters->surface(face, use.uv);
                    double residual = distance(surface.position, curve.position);
                    require(residual > 2.3e-7 && residual < 2.5e-7,
                            "Native SameParameter residual was corrected or lost during source "
                            "capture.");
                    require(std::abs(use.modelToleranceMetres - 1e-7) < 1e-15,
                            "Native edge/face tolerance was inflated to mask its pcurve residual.");
                }
    require(matches.size() == 1, "Native spline fixture lost its authoritative face use.");
}

void nativeSourceDomain() {
    // A valid closed box with a native pcurve approximation just beyond its
    // finite spline basis reproduces PH170's source-domain discrepancy.
    auto shape = BRepBuilderAPI_NurbsConvert(BRepPrimAPI_MakeBox(.02, .03, .04).Shape()).Shape();
    TopTools_IndexedMapOfShape faces;
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    const auto face = TopoDS::Face(faces(1));
    const auto surface = BRep_Tool::Surface(face);
    double u0, u1, v0, v1;
    surface->Bounds(u0, u1, v0, v1);
    constexpr double offset = 5e-7;
    bool changed = false;
    for (TopExp_Explorer x(face, TopAbs_EDGE); x.More(); x.Next()) {
        const auto edge = TopoDS::Edge(x.Current());
        double first, last;
        const auto pcurve = BRep_Tool::CurveOnSurface(edge, face, first, last);
        const auto middle = pcurve->Value(std::midpoint(first, last));
        if (std::abs(middle.X() - u1) > 1e-12)
            continue;
        auto shifted = Handle(Geom2d_Curve)::DownCast(pcurve->Copy());
        shifted->Translate(gp_Vec2d(offset, 0));
        BRep_Builder builder;
        builder.UpdateEdge(edge, shifted, face, 1e-5);
        builder.UpdateVertex(TopExp::FirstVertex(edge), 1e-5);
        builder.UpdateVertex(TopExp::LastVertex(edge), 1e-5);
        changed = true;
        break;
    }
    require(changed, "Synthetic native source halo was not constructed.");
    const auto cad = geometry(shape);
    const auto source = cad->meshingSource(.00025);
    const ParameterPoint outside{u1 + offset, std::midpoint(v0, v1)};
    const auto domain = source.parameters->domain(101);
    require(domain.uBounds && (*domain.uBounds)[1] >= outside.u,
            "Native pcurve source bounds were replaced by nominal spline bounds.");
    const auto actual = source.parameters->surface(101, outside);
    gp_Pnt expected;
    gp_Vec du, dv;
    surface->D1(outside.u, outside.v, expected, du, dv);
    require(distance(actual.position, {expected.X(), expected.Y(), expected.Z()}) < 1e-15 &&
                distance(actual.du, {du.X(), du.Y(), du.Z()}) < 1e-15 &&
                distance(actual.dv, {dv.X(), dv.Y(), dv.Z()}) < 1e-15,
            "Source halo position or D1 was clamped to its spline basis.");
    const auto boundary = surface->Value(u1, outside.v);
    require(expected.Distance(boundary) > BRep_Tool::Tolerance(face),
            "Synthetic halo did not require its native edge tolerance beyond face tolerance.");
    rejects([&] { (void)source.parameters->surface(101, {u1 + 2 * offset, outside.v}); });
    // The immutable same source domain does not grant a coarser modelling
    // allowance when a stricter mesh accuracy was explicitly requested.
    const auto strict = cad->meshingSource(1e-7);
    rejects([&] { (void)strict.parameters->surface(101, outside); });
}

void cancellationDuringTessellation() {
    const auto cad = geometry(BRepPrimAPI_MakeSphere(.1).Shape());
    for (const auto delay : {std::chrono::milliseconds(1), std::chrono::milliseconds(5)}) {
        std::stop_source stop;
        std::promise<void> started;
        auto ready = started.get_future();
        auto work = std::async(std::launch::async, [&] {
            started.set_value();
            try {
                (void)cad->meshingSource(1e-7, stop.get_token());
                return false;
            } catch (const std::runtime_error &error) {
                return std::string(error.what()).find("abgebrochen") != std::string::npos;
            }
        });
        ready.wait();
        std::this_thread::sleep_for(delay);
        stop.request_stop();
        require(work.get(), "In-flight native tessellation did not finish with cancellation.");
    }
}

} // namespace
int main() {
    try {
        nativeSnapshots();
        explicitTrimEvidence();
        periodicUses();
        shiftedPeriodicSourceDomain();
        rawSplineUse();
        endpointTolerances();
        nativeSourceDomain();
        cancellationDuringTessellation();
        std::cout << "Native CAD parameter snapshots: UV before weld, seam uses, true D1, "
                     "immutable lifetime and cancellation PASS\n";
        return 0;
    } catch (std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
