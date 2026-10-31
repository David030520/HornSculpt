#include "../modules/geometry/cad/cad_parameterization.hpp"
#include "../modules/meshing/src/frontal_delaunay.hpp"
#include "../modules/meshing/src/surface_charts.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <hornsculpt/cad_geometry.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numbers>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <tuple>

using namespace hornsculpt;
using namespace hornsculpt::detail;
namespace {
using AtlasSegment = std::array<std::uint32_t, 2>;
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
Vec3 sub(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 mul(Vec3 a, double b) {
    return {a.x * b, a.y * b, a.z * b};
}
double length(Vec3 p) {
    return std::hypot(p.x, p.y, p.z);
}
class BoundedPlane final : public SurfaceParameterization {
  public:
    ParameterDomain limits{{std::array{0., 1.}}, {std::array{0., 1.}}, {}, {}, 1};
    mutable unsigned outsideEvaluations{};
    const SurfaceMesh &sourceMesh() const noexcept override {
        return mesh_;
    }
    std::span<const ParameterTriangle> triangles() const noexcept override {
        return {};
    }
    ParameterDomain domain(std::uint32_t) const override {
        return limits;
    }
    ParametricSurfaceSample surface(std::uint32_t, ParameterPoint uv) const override {
        if ((!limits.uPeriod && (uv.u < 0 || uv.u > 1)) ||
            (!limits.vPeriod && (uv.v < 0 || uv.v > 1))) {
            ++outsideEvaluations;
            throw std::out_of_range("Test plane source domain exceeded.");
        }
        return {{uv.u, uv.v, 0}, {1, 0, 0}, {0, 1, 0}};
    }
    ParametricCurveSample curve(std::uint32_t, double) const override {
        throw std::logic_error("Plane candidate test has no physical curve.");
    }
    CurveUseSample curveUse(std::uint32_t, double) const override {
        throw std::logic_error("Plane candidate test has no physical curve use.");
    }
    std::span<const double> curvePartition(std::uint32_t) const override {
        throw std::logic_error("Plane candidate test has no physical curve.");
    }

  private:
    SurfaceMesh mesh_;
};
void optionalCandidateDomain() {
    auto provider = std::make_shared<BoundedPlane>();
    SurfaceChart chart;
    chart.faceId = 1;
    chart.parameterization = provider;
    chart.loops = {{{{{0, 0}, 0}, {{1, 0}, 1}, {{1, 1}, 2}, {{0, 1}, 3}}}};
    const ChartPoint outside{-1e-13, .5};
    require(chart.contains(outside),
            "Regression point no longer exercises polygon boundary tolerance.");
    require(!chart.containsParameterDomain(outside),
            "Polygon boundary tolerance extended the native source domain.");
    require(chart.containsParameterDomain({0, .5}) && chart.containsParameterDomain({1, 1}),
            "Actual source-domain endpoints were rejected.");
    provider->limits.uPeriod = 1;
    require(chart.containsParameterDomain({-1, .5}) && chart.containsParameterDomain({2, .5}),
            "Genuine periodic lifts were rejected as outside candidates.");
    chart.parameterOrientation = -1;
    require(chart.containsParameterDomain({0, -.5}) && !chart.containsParameterDomain({0, .5}),
            "Candidate source-domain checking ignored the reflected native V coordinate.");
    require(provider->outsideEvaluations == 0,
            "Optional candidate-domain checking evaluated an invalid surface point.");
}
class ClusteredSpanParameters final : public SurfaceParameterization {
  public:
    explicit ClusteredSpanParameters(bool reverse) : reverse_(reverse) {
        mesh_.vertices = {{0, 0, 0}, {.016, 0, 0}, {.016, .01, 0}, {0, .01, 0}};
        mesh_.triangles = {{{0, 1, 2}, 1, BoundaryCondition::wall},
                           {{0, 2, 3}, 1, BoundaryCondition::wall}};
        records_.resize(2);
        records_[0].uv = {{{0, 0}, {.016, 0}, {.016, .01}}};
        records_[1].uv = {{{0, 0}, {.016, .01}, {0, .01}}};
        records_[0].edges[0] = CurveUseInterval{1, 1, reverse ? 1. : 0., reverse ? 0. : 1.};
    }
    const SurfaceMesh &sourceMesh() const noexcept override {
        return mesh_;
    }
    std::span<const ParameterTriangle> triangles() const noexcept override {
        return records_;
    }
    ParameterDomain domain(std::uint32_t) const override {
        return {{std::array{0., .016}}, {std::array{0., .01}}, {}, {}, 1};
    }
    ParametricSurfaceSample surface(std::uint32_t, ParameterPoint uv) const override {
        return {{uv.u, uv.v, 0}, {1, 0, 0}, {0, 1, 0}};
    }
    ParametricCurveSample curve(std::uint32_t, double t) const override {
        // An exact straight edge with a narrow cubic parameter span. Its
        // localized speed is invisible to whole/half GL5 panels on [0,1].
        const double s = std::clamp((t - partition_[1]) / (partition_[2] - partition_[1]), 0., 1.);
        const double x = .006 * t + .01 * s * s * (3 - 2 * s);
        const double speed = .006 + .06 * s * (1 - s) / (partition_[2] - partition_[1]);
        return {{reverse_ ? .016 - x : x, 0, 0}, {reverse_ ? -speed : speed, 0, 0}};
    }
    CurveUseSample curveUse(std::uint32_t id, double t) const override {
        const auto value = curve(id, t);
        return {{value.position.x, 0}, {value.tangent.x, 0}, 0, {}};
    }
    std::span<const double> curvePartition(std::uint32_t) const override {
        return partition_;
    }

  private:
    bool reverse_;
    SurfaceMesh mesh_;
    std::vector<ParameterTriangle> records_;
    const std::array<double, 4> partition_{0, .000235, .000576, 1};
};
void clusteredCurveSpans() {
    for (const bool reverse : {false, true}) {
        auto parameters = std::make_shared<ClusteredSpanParameters>(reverse);
        MeshingOptions settings;
        settings.parameterization = parameters;
        settings.targetEdgeLength = .1;
        settings.adaptive.enabled = false;
        auto atlas = makeSurfaceAtlas(parameters->sourceMesh(), settings);
        std::optional<AtlasSegment> selected;
        for (const auto &chart : atlas.charts)
            for (const auto &loop : chart.loops)
                for (std::size_t i = 0; i < loop.vertices.size(); ++i) {
                    const auto a = loop.vertices[i].atlasVertex;
                    const auto b = loop.vertices[(i + 1) % loop.vertices.size()].atlasVertex;
                    const auto p = atlas.vertices[a], q = atlas.vertices[b];
                    if (p.y == 0 && q.y == 0 && std::abs(p.x - q.x) > .01)
                        selected = AtlasSegment{a, b};
                }
        require(selected.has_value(), "Clustered-span fixture lost its long straight edge.");
        const double halfLength =
            (atlas.vertices[(*selected)[0]].x + atlas.vertices[(*selected)[1]].x) / 2;
        const auto before = atlas.vertices.size();
        const std::array<AtlasSegment, 1> segments{*selected};
        require(!refineSurfaceAtlas(atlas, segments).empty(), "Clustered span did not refine.");
        require(atlas.vertices.size() == before + 1, "Curve split inserted unexpected vertices.");
        const auto midpoint = atlas.vertices.back();
        require(std::abs(midpoint.x - halfLength) < 2e-10 && midpoint.y == 0 && midpoint.z == 0,
                "Native curve split missed the physical half length of a narrow span.");
    }
}
class WarpedCylinderParameters final : public SurfaceParameterization {
  public:
    explicit WarpedCylinderParameters(ParameterPoint scale) : scale_(scale) {}
    const SurfaceMesh &sourceMesh() const noexcept override {
        return unused_;
    }
    std::span<const ParameterTriangle> triangles() const noexcept override {
        return {};
    }
    ParameterDomain domain(std::uint32_t) const override {
        return {{std::array{0., scale_.u}},
                {std::array{0., scale_.v * std::numbers::pi / 2}},
                {},
                {},
                1};
    }
    ParametricSurfaceSample surface(std::uint32_t, ParameterPoint uv) const override {
        const double u = uv.u / scale_.u, v = uv.v / scale_.v;
        return {{.04 * std::cos(v), .04 * std::sin(v), .04 * (u + 2 * u * u)},
                {0, 0, .04 * (1 + 4 * u) / scale_.u},
                {-.04 * std::sin(v) / scale_.v, .04 * std::cos(v) / scale_.v, 0}};
    }
    ParametricCurveSample curve(std::uint32_t, double) const override {
        throw std::logic_error("Analytic front fixture has no native curve records.");
    }
    CurveUseSample curveUse(std::uint32_t, double) const override {
        throw std::logic_error("Analytic front fixture has no native curve uses.");
    }
    std::span<const double> curvePartition(std::uint32_t) const override {
        throw std::logic_error("Analytic front fixture has no native curve records.");
    }

  private:
    ParameterPoint scale_;
    SurfaceMesh unused_;
};
void varyingNativeMetric(ParameterPoint scale) {
    auto parameters = std::make_shared<WarpedCylinderParameters>(scale);
    SurfaceChart chart;
    chart.faceId = 1;
    chart.parameterization = parameters;
    std::vector<Vec3> positions;
    SurfaceChartLoop loop;
    const auto uv = [&](double u, double v) { return ChartPoint{u * scale.u, v * scale.v}; };
    const auto add = [&](ChartPoint point) {
        loop.vertices.push_back({point, static_cast<std::uint32_t>(positions.size())});
        positions.push_back(parameters->surface(1, {point.u, point.v}).position);
    };
    for (unsigned i = 0; i < 48; ++i)
        add(uv(i / 48., 0));
    for (unsigned i = 0; i < 40; ++i)
        add(uv(1, i * std::numbers::pi / 80));
    for (unsigned i = 0; i < 48; ++i)
        add(uv(1 - i / 48., std::numbers::pi / 2));
    for (unsigned i = 0; i < 6; ++i)
        add(uv(0, (1 - i / 6.) * std::numbers::pi / 2));
    chart.loops.push_back(std::move(loop));
    const auto reference = [&](std::array<ChartPoint, 3> points) {
        ChartReferenceTriangle triangle;
        triangle.uv = points;
        triangle.sizes = {.008, .008, .008};
        for (unsigned i = 0; i < 3; ++i)
            triangle.positions[i] = parameters->surface(1, {points[i].u, points[i].v}).position;
        chart.reference.push_back(triangle);
    };
    for (unsigned i = 0; i < 32; ++i)
        for (unsigned j = 0; j < 24; ++j) {
            const auto a = uv(i / 32., j * std::numbers::pi / 48);
            const auto b = uv((i + 1) / 32., j * std::numbers::pi / 48);
            const auto c = uv((i + 1) / 32., (j + 1) * std::numbers::pi / 48);
            const auto d = uv(i / 32., (j + 1) * std::numbers::pi / 48);
            reference({a, b, c});
            reference({a, c, d});
        }
    MeshingOptions settings;
    settings.targetEdgeLength = .008;
    for (unsigned pass = 0; pass <= 12; ++pass) {
        try {
            const auto mesh = meshFrontalChart(chart, positions, settings);
            for (const auto &vertex : mesh.vertices)
                if (vertex.atlasVertex != std::numeric_limits<std::uint32_t>::max())
                    require(vertex.position == positions.at(vertex.atlasVertex),
                            "Physical relaxation moved a mandatory boundary/fixed point.");
            for (const auto &triangle : mesh.triangles) {
                std::array<Vec3, 3> p;
                ChartPoint center{};
                for (unsigned i = 0; i < 3; ++i) {
                    const auto &v = mesh.vertices.at(triangle[i]);
                    p[i] = v.position;
                    center.u += v.uv.u / 3;
                    center.v += v.uv.v / 3;
                }
                const auto cross = [](Vec3 a, Vec3 b) {
                    return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                                a.x * b.y - a.y * b.x};
                };
                const auto dot = [](Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; };
                const auto n = cross(sub(p[1], p[0]), sub(p[2], p[0]));
                const auto local = chart.sample(center);
                require(dot(n, cross(local.du, local.dv)) > 0,
                        "Varying-metric native front reversed an output triangle.");
                double longest = 0;
                for (unsigned i = 0; i < 3; ++i) {
                    const auto a = sub(p[(i + 1) % 3], p[i]), b = sub(p[(i + 2) % 3], p[i]);
                    longest = std::max(longest, length(a));
                    require(std::atan2(length(n), dot(a, b)) * 180 / std::numbers::pi >= 20,
                            "Varying-metric native front violated the actual minimum angle.");
                }
                require(longest * longest / length(n) <= 4,
                        "Varying-metric native front violated the actual aspect limit.");
            }
            const auto repeated = meshFrontalChart(chart, positions, settings);
            require(mesh.triangles == repeated.triangles &&
                        mesh.vertices.size() == repeated.vertices.size(),
                    "Native physical optimization changed on identical repeated inputs.");
            for (std::size_t i = 0; i < mesh.vertices.size(); ++i)
                require(mesh.vertices[i].position == repeated.vertices[i].position &&
                            mesh.vertices[i].uv.u == repeated.vertices[i].uv.u &&
                            mesh.vertices[i].uv.v == repeated.vertices[i].uv.v,
                        "Native physical optimization is not exactly repeatable.");
            return;
        } catch (const BoundaryRefinementRequired &required) {
            for (const auto segment : required.segments()) {
                auto &nodes = chart.loops.front().vertices;
                bool found = false;
                for (std::size_t i = 0; i < nodes.size(); ++i) {
                    const auto a = nodes[i], b = nodes[(i + 1) % nodes.size()];
                    if (AtlasSegment{std::min(a.atlasVertex, b.atlasVertex),
                                     std::max(a.atlasVertex, b.atlasVertex)} != segment)
                        continue;
                    const ChartPoint point{std::midpoint(a.uv.u, b.uv.u),
                                           std::midpoint(a.uv.v, b.uv.v)};
                    const auto id = static_cast<std::uint32_t>(positions.size());
                    positions.push_back(parameters->surface(1, {point.u, point.v}).position);
                    nodes.insert(nodes.begin() + i + 1, {point, id});
                    found = true;
                    break;
                }
                require(found, "Native metric test requested a non-current boundary segment.");
            }
        }
    }
    throw std::runtime_error("Native varying-metric front exhausted its shared refinement passes.");
}
std::shared_ptr<const CadGeometry> geometry(const TopoDS_Shape &shape) {
    std::ostringstream bytes;
    BRepTools::Write(shape, bytes, Standard_False, Standard_False, TopTools_FormatVersion_CURRENT);
    TopTools_IndexedMapOfShape faces;
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    std::vector<std::uint32_t> ids(static_cast<std::size_t>(faces.Extent()));
    std::iota(ids.begin(), ids.end(), 1);
    const CadBrepPart part{bytes.str(), ids};
    return CadGeometry::fromBrepParts(std::span{&part, 1}, 1);
}
class ObservedParameters final : public SurfaceParameterization {
  public:
    struct Use {
        std::uint32_t use, curve, face;
        double parameter;
        CurveUseSample sample;
    };
    explicit ObservedParameters(std::shared_ptr<const SurfaceParameterization> source)
        : source_(std::move(source)) {
        for (std::size_t f = 0; f < source_->triangles().size(); ++f)
            for (const auto &edge : source_->triangles()[f].edges)
                if (edge)
                    identity_.emplace(
                        edge->useId,
                        std::pair{edge->curveId, source_->sourceMesh().triangles[f].faceId});
    }
    const SurfaceMesh &sourceMesh() const noexcept override {
        return source_->sourceMesh();
    }
    std::span<const ParameterTriangle> triangles() const noexcept override {
        return source_->triangles();
    }
    ParameterDomain domain(std::uint32_t face) const override {
        return source_->domain(face);
    }
    ParametricSurfaceSample surface(std::uint32_t face, ParameterPoint uv) const override {
        return source_->surface(face, uv);
    }
    ParametricCurveSample curve(std::uint32_t id, double t) const override {
        return source_->curve(id, t);
    }
    std::span<const double> curvePartition(std::uint32_t id) const override {
        return source_->curvePartition(id);
    }
    CurveUseSample curveUse(std::uint32_t id, double t) const override {
        const auto result = source_->curveUse(id, t);
        const auto [curve, face] = identity_.at(id);
        const std::scoped_lock lock(mutex_);
        calls_.push_back({id, curve, face, t, result});
        return result;
    }
    std::vector<Use> calls() const {
        const std::scoped_lock lock(mutex_);
        return calls_;
    }

  private:
    std::shared_ptr<const SurfaceParameterization> source_;
    std::map<std::uint32_t, std::pair<std::uint32_t, std::uint32_t>> identity_;
    mutable std::mutex mutex_;
    mutable std::vector<Use> calls_;
};
MeshingOptions options(const std::shared_ptr<const CadGeometry> &cad,
                       std::shared_ptr<const SurfaceParameterization> parameters) {
    MeshingOptions result;
    result.geometry = cad;
    result.parameterization = std::move(parameters);
    result.targetEdgeLength = .008;
    result.adaptive.enabled = false;
    return result;
}
bool sameModulo(double a, double b, std::optional<double> period) {
    double delta = a - b;
    if (period)
        delta -= std::nearbyint(delta / *period) * *period;
    return std::abs(delta) <=
           512 * std::numeric_limits<double>::epsilon() * std::max({1., std::abs(a), std::abs(b)});
}
void checkActiveBinding(const SurfaceAtlas &atlas, const ObservedParameters &parameters) {
    const auto calls = parameters.calls();
    std::set<std::tuple<double, double, double>> originals;
    for (auto p : parameters.sourceMesh().vertices)
        originals.emplace(p.x, p.y, p.z);
    std::size_t insertedPhysical = 0;
    for (const auto &chart : atlas.charts) {
        require(chart.parameterization != nullptr, "An exact native chart lost its evaluator.");
        const auto domain = parameters.domain(chart.faceId);
        for (const auto &loop : chart.loops)
            for (const auto &node : loop.vertices) {
                const auto p = atlas.vertices.at(node.atlasVertex);
                const ParameterPoint uv{node.uv.u, chart.parameterOrientation * node.uv.v};
                const auto mapped = parameters.surface(chart.faceId, uv).position;
                const bool original = originals.contains({p.x, p.y, p.z});
                bool physical = false, bound = false;
                double model = 0;
                for (const auto &call : calls) {
                    if (call.face != chart.faceId)
                        continue;
                    const auto curve = parameters.curve(call.curve, call.parameter).position;
                    const double fp = 512 * std::numeric_limits<double>::epsilon() *
                                      std::max(length(p), length(curve));
                    if (length(sub(curve, p)) > call.sample.modelToleranceMetres + fp)
                        continue;
                    physical = true;
                    model = std::max(model, call.sample.modelToleranceMetres);
                    bound |= sameModulo(uv.u, call.sample.uv.u, domain.uPeriod) &&
                             sameModulo(uv.v, call.sample.uv.v, domain.vPeriod);
                }
                if (physical && !original) {
                    require(bound,
                            "Inserted native boundary UV is not its actual curve-use parameter.");
                    ++insertedPhysical;
                }
                const double fp = 1024 * std::numeric_limits<double>::epsilon() *
                                  std::max(length(p), length(mapped));
                // Source/topological curve points may carry a native representation
                // residual; inserted artificial cuts must lie directly on S.
                require(
                    length(sub(p, mapped)) <= std::max(model, 1e-12) + fp,
                    "Active native boundary XYZ and S(storedUV) disagree beyond model precision.");
                const auto sample = chart.sample(node.uv);
                require(length(sub(sample.position, mapped)) <= fp,
                        "Native chart evaluation is not direct S at its stored parameters.");
            }
        for (const auto &t : chart.reference) {
            const ChartPoint uv{(t.uv[0].u + t.uv[1].u + t.uv[2].u) / 3,
                                (t.uv[0].v + t.uv[1].v + t.uv[2].v) / 3};
            const auto value = chart.sample(uv);
            const double step = 1e-6;
            const auto du = mul(sub(chart.sample({uv.u + step, uv.v}).position,
                                    chart.sample({uv.u - step, uv.v}).position),
                                .5 / step);
            const auto dv = mul(sub(chart.sample({uv.u, uv.v + step}).position,
                                    chart.sample({uv.u, uv.v - step}).position),
                                .5 / step);
            require(length(sub(du, value.du)) <= 1e-5 * std::max(length(value.du), 1e-8) &&
                        length(sub(dv, value.dv)) <= 1e-5 * std::max(length(value.dv), 1e-8),
                    "Native chart derivatives do not match independent position differences.");
        }
    }
    require(insertedPhysical > 0, "No genuinely inserted native physical sample was checked.");
}
void refineAndCheck(SurfaceAtlas &atlas, const ObservedParameters &parameters) {
    checkActiveBinding(atlas, parameters);
    for (unsigned pass = 0; pass < 3; ++pass) {
        std::set<AtlasSegment> segments;
        for (const auto &chart : atlas.charts)
            for (const auto &loop : chart.loops)
                if (loop.vertices.size() > 2)
                    segments.insert(
                        {std::min(loop.vertices[0].atlasVertex, loop.vertices[1].atlasVertex),
                         std::max(loop.vertices[0].atlasVertex, loop.vertices[1].atlasVertex)});
        const std::vector<AtlasSegment> selected(segments.begin(), segments.end());
        std::cout << "Native refine pass=" << pass << " segments=" << selected.size() << std::endl;
        require(!refineSurfaceAtlas(atlas, selected).empty(),
                "Native split did not change any current charts.");
        checkActiveBinding(atlas, parameters);
    }
}
void nativeCylinderAndHole() {
    for (const auto &cad :
         {geometry(BRepPrimAPI_MakeCylinder(.01, .02).Shape()),
          geometry(
              BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(.03, .03, .008).Shape(),
                              BRepPrimAPI_MakeCylinder(
                                  gp_Ax2(gp_Pnt(.015, .015, -.004), gp_Dir(0, 0, 1)), .004, .016)
                                  .Shape())
                  .Shape())}) {
        std::cout << "Native analytic source faces=" << cad->faces().size() << std::endl;
        const auto source = cad->meshingSource(.00025);
        require(inspectMesh(source.mesh).watertight(), "Native source is not a closed boundary.");
        auto observed = std::make_shared<ObservedParameters>(source.parameters);
        auto assigned = source.mesh;
        const auto assignedFace = assigned.triangles.front().faceId;
        for (auto &triangle : assigned.triangles)
            if (triangle.faceId == assignedFace)
                triangle.boundary = BoundaryCondition::driven;
        auto atlas = makeSurfaceAtlas(assigned, options(cad, observed));
        if (cad->faces().size() == 3) {
            std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<ParameterPoint>> uses;
            for (const auto &chart : atlas.charts)
                for (const auto &loop : chart.loops)
                    for (auto node : loop.vertices)
                        uses[{chart.faceId, node.atlasVertex}].push_back(
                            {node.uv.u, chart.parameterOrientation * node.uv.v});
            bool periodicSeam = false;
            for (const auto &[key, uv] : uses) {
                const auto domain = source.parameters->domain(key.first);
                for (auto a : uv)
                    for (auto b : uv)
                        periodicSeam |= domain.uPeriod &&
                                        std::abs(std::abs(a.u - b.u) - *domain.uPeriod) < 1e-12;
            }
            require(periodicSeam,
                    "The periodic cylinder lost its distinct active native UV seam uses.");
        } else {
            require(std::any_of(atlas.charts.begin(), atlas.charts.end(),
                                [](const auto &chart) { return chart.loops.size() > 1; }),
                    "Native through-hole trim loops were lost.");
        }
        refineAndCheck(atlas, *observed);
        const auto &snapshot = observed->sourceMesh();
        auto wrong = source.mesh;
        wrong.vertices.front().x += .001;
        bool rejected = false;
        try {
            (void)makeSurfaceAtlas(wrong, options(cad, observed));
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Wrong native snapshot positions were silently accepted.");
        require(snapshot.vertices == source.mesh.vertices,
                "A split mutated the immutable provider snapshot.");
    }
}
void nativeSplineBinding() {
    std::ifstream stream(HS_SPLINE_PARAMETER_FIXTURE_PATH, std::ios::binary);
    require(static_cast<bool>(stream), "Missing closed native spline fixture.");
    std::ostringstream bytes;
    bytes << stream.rdbuf();
    const CadBrepPart part{bytes.str(), {1, 2, 3, 4, 5, 6, 7, 8}};
    const auto cad = CadGeometry::fromBrepParts(std::span{&part, 1}, 1);
    std::cout << "Native spline source faces=" << cad->faces().size() << std::endl;
    const auto source = cad->meshingSource(.00025);
    TopoDS_Shape native;
    BRep_Builder builder;
    std::istringstream nativeBytes(bytes.str());
    BRepTools::Read(native, nativeBytes, builder);
    TopTools_IndexedMapOfShape edges;
    TopExp::MapShapes(native, TopAbs_EDGE, edges);
    std::uint32_t curveId = 0;
    std::size_t checkedKnots = 0;
    for (int i = 1; i <= edges.Extent(); ++i) {
        const auto edge = TopoDS::Edge(edges(i));
        if (BRep_Tool::Degenerated(edge))
            continue;
        BRepAdaptor_Curve curve(edge);
        const auto partition = source.parameters->curvePartition(++curveId);
        require(partition.size() >= 2 && partition.front() == curve.FirstParameter() &&
                    partition.back() == curve.LastParameter(),
                "Native curve partition changed its true trim.");
        if (curve.GetType() == GeomAbs_BSplineCurve) {
            const auto spline = curve.BSpline();
            for (int k = 1; k <= spline->NbKnots(); ++k) {
                const double knot = spline->Knot(k);
                if (knot > curve.FirstParameter() && knot < curve.LastParameter()) {
                    ++checkedKnots;
                    require(std::binary_search(partition.begin(), partition.end(), knot),
                            "Native partition omitted a spline polynomial span.");
                }
            }
        }
    }
    require(checkedKnots > 0, "Native spline fixture checked no inner polynomial spans.");
    auto observed = std::make_shared<ObservedParameters>(source.parameters);
    auto atlas = makeSurfaceAtlas(source.mesh, options(cad, observed));
    refineAndCheck(atlas, *observed);
}
void nativeSmoothCap() {
    TColgp_Array1OfPnt poles(1, 4);
    poles(1) = gp_Pnt(0, 0, 0);
    poles(2) = gp_Pnt(.003, 0, 0);
    poles(3) = gp_Pnt(.009, 0, -.002);
    poles(4) = gp_Pnt(.012, 0, -.005);
    TColStd_Array1OfReal knots(1, 2);
    knots(1) = 0;
    knots(2) = 1;
    TColStd_Array1OfInteger multiplicities(1, 2);
    multiplicities(1) = multiplicities(2) = 4;
    const Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, knots, multiplicities, 3);
    BRepBuilderAPI_MakeWire wire;
    wire.Add(BRepBuilderAPI_MakeEdge(curve).Edge());
    const gp_Pnt base(0, 0, -.005);
    wire.Add(BRepBuilderAPI_MakeEdge(poles(4), base).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(base, poles(1)).Edge());
    const auto solid = BRepPrimAPI_MakeRevol(BRepBuilderAPI_MakeFace(wire.Wire()).Face(),
                                             gp_Ax1(gp_Pnt(), gp_Dir(0, 0, 1)))
                           .Shape();
    const auto cad = geometry(solid);
    const auto source = cad->meshingSource(.0001);
    const auto face = std::find_if(cad->faces().begin(), cad->faces().end(), [](const auto &f) {
        return f.surfaceType == "SurfaceOfRevolution";
    });
    require(face != cad->faces().end(), "Smooth native rotational cap fixture is missing.");
    const auto id = face->faceId;
    const auto domain = source.parameters->domain(id);
    require(!domain.uPeriod && !domain.vPeriod,
            "A smooth collapsed rotational pole retained its singular periodic chart.");
    require(source.parameters->isInsideTrim(id, {.001, .002}),
            "An interior Cartesian cap point lost its strict physical-trim evidence.");
    for (const ParameterPoint p :
         {ParameterPoint{}, ParameterPoint{.001, .002}, ParameterPoint{-.005, .003}}) {
        constexpr double step = 1e-7;
        const auto exact = source.parameters->surface(id, p);
        const auto du = mul(sub(source.parameters->surface(id, {p.u + step, p.v}).position,
                                source.parameters->surface(id, {p.u - step, p.v}).position),
                            .5 / step);
        const auto dv = mul(sub(source.parameters->surface(id, {p.u, p.v + step}).position,
                                source.parameters->surface(id, {p.u, p.v - step}).position),
                            .5 / step);
        require(length(sub(du, exact.du)) < 1e-5 && length(sub(dv, exact.dv)) < 1e-5 &&
                    length(exact.du) > 0 && length(exact.dv) > 0,
                "Cartesian native cap D1 disagrees with exact position differences.");
    }
    std::size_t checkedUses = 0;
    for (std::size_t i = 0; i < source.mesh.triangles.size(); ++i)
        if (source.mesh.triangles[i].faceId == id)
            for (const auto &edge : source.parameters->triangles()[i].edges)
                if (edge) {
                    const double t = (edge->t0 + edge->t1) / 2;
                    const double step = std::abs(edge->t1 - edge->t0) * 1e-5;
                    const auto use = source.parameters->curveUse(edge->useId, t);
                    const auto before = source.parameters->curveUse(edge->useId, t - step).uv;
                    const auto after = source.parameters->curveUse(edge->useId, t + step).uv;
                    require(std::hypot((after.u - before.u) / (2 * step) - use.uvDerivative.u,
                                       (after.v - before.v) / (2 * step) - use.uvDerivative.v) <
                                1e-6,
                            "The cap coordinate change lost its exact pcurve derivative.");
                    const auto onSurface = source.parameters->surface(id, use.uv).position;
                    require(length(sub(onSurface,
                                       source.parameters->curve(edge->curveId, t).position)) <
                                use.modelToleranceMetres + 1e-12,
                            "Cartesian cap source and canonical curve bindings disagree.");
                    ++checkedUses;
                }
    require(checkedUses > 0, "No native cap physical curve binding was checked.");
    auto settings = options(cad, source.parameters);
    settings.targetEdgeLength = .003;
    settings.adaptive = {true, .0002, .00005, .3};
    const auto atlas = makeSurfaceAtlas(source.mesh, settings);
    for (const auto &chart : atlas.charts)
        if (chart.faceId == id) {
            require(chart.loops.size() == 1,
                    "A native cap's former parameter seam became a physical boundary.");
            for (const auto &node : chart.loops.front().vertices) {
                const auto p = atlas.vertices[node.atlasVertex];
                require(std::hypot(p.x, p.y) > .011,
                        "An artificial meridian or pole remained constrained inside the cap.");
            }
        }
    const auto mesh = remeshIsotropic(source.mesh, settings);
    require(mesh.quality.watertight() && mesh.quality.selfIntersectionPairs == 0 &&
                mesh.quality.minimumAngleDegrees >= 20 && mesh.quality.maximumAspectRatio <= 4,
            "A smooth native pole violates the common isotropic mesh contract.");
    std::cout << "Native smooth cap: " << mesh.mesh.triangles.size() << " triangles, "
              << mesh.quality.minimumAngleDegrees << " deg / " << mesh.quality.maximumAspectRatio
              << " aspect\n";
}

void nativeAnnularOffsetRemainsPeriodic() {
    // A large placement cannot turn a measurable 5 mm inner radius into an
    // axis pole. Certification uses the physical radial scale, not the offset.
    constexpr double offset = 1e12;
    TColgp_Array1OfPnt poles(1, 4);
    poles(1) = gp_Pnt(offset + .005, 0, 0);
    poles(2) = gp_Pnt(offset + .006, 0, 0);
    poles(3) = gp_Pnt(offset + .009, 0, -.002);
    poles(4) = gp_Pnt(offset + .012, 0, -.005);
    TColStd_Array1OfReal knots(1, 2);
    knots(1) = 0;
    knots(2) = 1;
    TColStd_Array1OfInteger multiplicities(1, 2);
    multiplicities(1) = multiplicities(2) = 4;
    const Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, knots, multiplicities, 3);
    const Handle(Geom_SurfaceOfRevolution) surface =
        new Geom_SurfaceOfRevolution(curve, gp_Ax1(gp_Pnt(offset, 0, 0), gp_Dir(0, 0, 1)));
    const auto face = BRepBuilderAPI_MakeFace(surface, 0, 2 * std::numbers::pi, 0, 1, 1e-7).Face();
    cad_detail::ParameterizationBuilder builder(.0001);
    builder.addFace(1, face);
    const auto provider = builder.finish({}, {});
    require(provider->domain(1).uPeriod.has_value(),
            "World placement incorrectly certified a finite annular radius as a collapsed pole.");
}

void nativeGradedSmoothCap() {
    // Four cubic spans and an increasing outer slope exercise a graded front
    // around a regular Cartesian pole without a project-specific data file.
    constexpr std::array radii{0., .001, .003, .006, .009, .011, .012};
    constexpr std::array heights{
        0., 0., -.0001516794727, -.0009877440049, -.0024940410552, -.0038057183773, -.00488778975};
    constexpr double offset = -.0276622389;
    TColgp_Array1OfPnt poles(1, 7);
    for (int i = 1; i <= 7; ++i)
        poles(i) = gp_Pnt(radii[i - 1], 0, offset + heights[i - 1]);
    TColStd_Array1OfReal knots(1, 5);
    TColStd_Array1OfInteger multiplicities(1, 5);
    for (int i = 1; i <= 5; ++i) {
        knots(i) = .25 * (i - 1);
        multiplicities(i) = i == 1 || i == 5 ? 4 : 1;
    }
    const Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, knots, multiplicities, 3);
    BRepBuilderAPI_MakeWire wire;
    wire.Add(BRepBuilderAPI_MakeEdge(curve).Edge());
    const gp_Pnt base(0, 0, offset + heights.back());
    wire.Add(BRepBuilderAPI_MakeEdge(poles(7), base).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(base, poles(1)).Edge());
    const auto cad = geometry(BRepPrimAPI_MakeRevol(BRepBuilderAPI_MakeFace(wire.Wire()).Face(),
                                                    gp_Ax1(gp_Pnt(), gp_Dir(0, 0, 1)))
                                  .Shape());
    const auto source = cad->meshingSource(.00025);
    auto settings = options(cad, source.parameters);
    settings.targetEdgeLength = .006;
    settings.adaptive = {true, .0006, .0001, .3};
    const auto face = std::find_if(cad->faces().begin(), cad->faces().end(), [](const auto &f) {
        return f.surfaceType == "SurfaceOfRevolution";
    });
    require(face != cad->faces().end(), "Graded smooth cap fixture has no rotational surface.");
    settings.faceTargetEdgeLengths = {{face->faceId, .002}};
    auto atlas = makeSurfaceAtlas(source.mesh, settings);
    for (auto &chart : atlas.charts)
        if (chart.faceId == face->faceId) {
            // All three positions lie on the exact surface, but their nearly
            // collinear UVs make PL derivative reconstruction ill-conditioned.
            // A preview facet cannot override the native differential metric.
            ChartReferenceTriangle thin;
            thin.uv = {{{-.006, 0}, {.006, 0}, {0, 1e-15}}};
            thin.sizes.fill(.002);
            for (int i = 0; i < 3; ++i)
                thin.positions[i] = chart.sample(thin.uv[i]).position;
            chart.reference.push_back(thin);
            const auto refined = meshFrontalChart(chart, atlas.vertices, settings);
            SurfaceMesh patch;
            for (const auto &vertex : refined.vertices)
                patch.vertices.push_back(vertex.position);
            for (const auto triangle : refined.triangles)
                patch.triangles.push_back({triangle, face->faceId});
            const auto quality = inspectMesh(patch);
            require(
                quality.degenerateTriangles == 0 && quality.selfIntersectionPairs == 0 &&
                    quality.minimumAngleDegrees >= 20 && quality.maximumAspectRatio <= 4,
                "A thin preview facet corrupted an exact native surface's differential metric.");
        }
    const auto mesh = remeshIsotropic(source.mesh, settings);
    require(mesh.quality.watertight() && mesh.quality.selfIntersectionPairs == 0 &&
                mesh.quality.minimumAngleDegrees >= 20 && mesh.quality.maximumAspectRatio <= 4,
            "A multi-span smooth cap stranded an underresolved or ill-shaped interior front.");
    std::cout << "Native graded smooth cap: " << mesh.mesh.triangles.size() << " triangles, "
              << mesh.quality.minimumAngleDegrees << " deg / " << mesh.quality.maximumAspectRatio
              << " aspect\n";
}

void nativeTorusDomain() {
    // A torus turns its normal through every octant but has one simple native
    // rectangular UV domain after its two periodic seams are cut. Projection
    // sectors must not create additional boundaries inside that domain.
    const auto cad = geometry(BRepPrimAPI_MakeTorus(.03, .008).Shape());
    const auto source = cad->meshingSource(.0005);
    auto observed = std::make_shared<ObservedParameters>(source.parameters);
    auto atlas = makeSurfaceAtlas(source.mesh, options(cad, observed));
    require(atlas.charts.size() == 1 && atlas.charts.front().loops.size() == 1,
            "A simple native torus UV domain was split by preview normal sectors.");
    require(atlas.charts.front().reference.size() >= source.mesh.triangles.size(),
            "Native chart construction lost curved source triangles.");
    const auto &loop = atlas.charts.front().loops.front().vertices;
    const std::array<AtlasSegment, 1> selected{
        {{std::min(loop[0].atlasVertex, loop[1].atlasVertex),
          std::max(loop[0].atlasVertex, loop[1].atlasVertex)}}};
    require(!refineSurfaceAtlas(atlas, selected).empty(),
            "Native periodic boundary refinement did not insert a sample.");
    refineAndCheck(atlas, *observed);
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 2 && std::string_view(argv[1]) == "--graded-cap-only") {
            nativeGradedSmoothCap();
            return 0;
        }
        optionalCandidateDomain();
        clusteredCurveSpans();
        varyingNativeMetric({1, 1});
        varyingNativeMetric({.001, 17});
        nativeCylinderAndHole();
        nativeTorusDomain();
        nativeSmoothCap();
        nativeGradedSmoothCap();
        nativeAnnularOffsetRemainsPeriodic();
        nativeSplineBinding();
        std::cout << "Native parameter snapshots, active initial/refined pcurve binding, periodic "
                     "uses, holes and exact D1 PASS.\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
