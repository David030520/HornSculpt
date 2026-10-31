#include "../modules/meshing/src/frontal_delaunay.hpp"
#include "../modules/meshing/src/surface_charts.hpp"
#include <hornsculpt/draft_geometry.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>

using namespace hornsculpt;
using namespace hornsculpt::detail;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
template <class F> void rejects(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const std::exception &) {
        rejected = true;
    }
    require(rejected, "An invalid atlas was accepted.");
}
Vec3 sub(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
double norm(Vec3 p) {
    return std::hypot(p.x, p.y, p.z);
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
bool sameUv(ChartPoint a, ChartPoint b) {
    return a.u == b.u && a.v == b.v;
}
bool sameNode(const SurfaceChartVertex &a, const SurfaceChartVertex &b) {
    return a.atlasVertex == b.atlasVertex && sameUv(a.uv, b.uv);
}
bool sameBoundaries(const SurfaceChart &a, const SurfaceChart &b) {
    return a.loops.size() == b.loops.size() &&
           std::equal(a.loops.begin(), a.loops.end(), b.loops.begin(),
                      [](const auto &x, const auto &y) {
                          return x.vertices.size() == y.vertices.size() &&
                                 std::equal(x.vertices.begin(), x.vertices.end(),
                                            y.vertices.begin(), sameNode);
                      });
}
void unchangedChartGeometry(const SurfaceChart &a, const SurfaceChart &b) {
    require(
        a.faceId == b.faceId && a.boundary == b.boundary && a.planar == b.planar &&
            a.planeOrigin == b.planeOrigin && a.planeU == b.planeU && a.planeV == b.planeV &&
            a.reference.size() == b.reference.size() &&
            a.fixedPoints.size() == b.fixedPoints.size() &&
            std::equal(a.fixedPoints.begin(), a.fixedPoints.end(), b.fixedPoints.begin(), sameNode),
        "Shared refinement changed chart tags, parametrization or fixed points.");
    for (std::size_t i = 0; i < a.reference.size(); ++i) {
        const auto &x = a.reference[i], &y = b.reference[i];
        require(x.positions == y.positions && x.sizes == y.sizes &&
                    std::equal(x.uv.begin(), x.uv.end(), y.uv.begin(), sameUv),
                "Shared refinement changed immutable reference geometry or sizing.");
        const ChartPoint centroid{(x.uv[0].u + x.uv[1].u + x.uv[2].u) / 3,
                                  (x.uv[0].v + x.uv[1].v + x.uv[2].v) / 3};
        const auto first = a.sample(centroid), second = b.sample(centroid);
        require(first.position == second.position && first.du == second.du &&
                    first.dv == second.dv && first.targetEdgeLength == second.targetEdgeLength,
                "Shared refinement changed surface evaluation or its physical metric.");
    }
}
auto physicalTriangles(const FrontalChartMesh &mesh) {
    std::vector<std::array<std::array<double, 3>, 3>> triangles;
    for (const auto &triangle : mesh.triangles) {
        std::array<std::array<double, 3>, 3> positions;
        for (unsigned i = 0; i < 3; ++i) {
            const auto p = mesh.vertices.at(triangle[i]).position;
            positions[i] = {p.x, p.y, p.z};
        }
        // Canonical cyclic order preserves winding while ignoring local numbering.
        std::rotate(positions.begin(), std::min_element(positions.begin(), positions.end()),
                    positions.end());
        triangles.push_back(positions);
    }
    std::sort(triangles.begin(), triangles.end());
    return triangles;
}
MeshingOptions settings() {
    MeshingOptions options;
    options.targetEdgeLength = .008;
    options.adaptive.geometricTolerance = .00025;
    return options;
}
SurfaceMesh ring(unsigned segments, bool cylinder) {
    SurfaceMesh mesh;
    for (unsigned layer = 0; layer < 2; ++layer)
        for (unsigned i = 0; i < segments; ++i) {
            const double angle = 2 * std::numbers::pi * i / segments;
            const double radius = cylinder ? .01 : (layer == 0 ? .03 : .01);
            mesh.vertices.push_back(
                {radius * std::cos(angle), radius * std::sin(angle), cylinder ? layer * .06 : 0});
        }
    for (unsigned i = 0; i < segments; ++i) {
        const unsigned j = (i + 1) % segments;
        if (cylinder) {
            mesh.triangles.push_back({{i, j, j + segments}, 7, BoundaryCondition::wall});
            mesh.triangles.push_back({{i, j + segments, i + segments}, 7, BoundaryCondition::wall});
        } else {
            mesh.triangles.push_back({{i, j, i + segments}, 7, BoundaryCondition::wall});
            mesh.triangles.push_back({{j, j + segments, i + segments}, 7, BoundaryCondition::wall});
        }
    }
    return mesh;
}
void referenceAgreement(const SurfaceAtlas &atlas) {
    for (const auto &chart : atlas.charts)
        for (const auto &triangle : chart.reference) {
            ChartPoint uv{};
            Vec3 expected{};
            for (unsigned i = 0; i < 3; ++i) {
                uv.u += triangle.uv[i].u / 3;
                uv.v += triangle.uv[i].v / 3;
                expected.x += triangle.positions[i].x / 3;
                expected.y += triangle.positions[i].y / 3;
                expected.z += triangle.positions[i].z / 3;
            }
            const auto value = chart.sample(uv);
            require(norm(sub(value.position, expected)) < 2e-10,
                    "Chart lifting changed the source geometry.");
            const double g00 = dot(value.du, value.du), g01 = dot(value.du, value.dv),
                         g11 = dot(value.dv, value.dv);
            require(g00 * g11 - g01 * g01 > 0, "The chart metric is singular.");
            require(value.targetEdgeLength > 0 && std::isfinite(value.targetEdgeLength),
                    "Nonfinite chart size.");
        }
    for (const auto &chart : atlas.charts)
        if (!chart.planar) {
            double lowerU = std::numeric_limits<double>::infinity(), lowerV = lowerU;
            double upperU = -lowerU, upperV = -lowerU;
            for (const auto &loop : chart.loops)
                for (const auto &node : loop.vertices) {
                    lowerU = std::min(lowerU, node.uv.u);
                    upperU = std::max(upperU, node.uv.u);
                    lowerV = std::min(lowerV, node.uv.v);
                    upperV = std::max(upperV, node.uv.v);
                }
            unsigned samples = 0;
            for (unsigned i = 0; i < 19; ++i)
                for (unsigned j = 0; j < 19; ++j) {
                    const ChartPoint uv{lowerU + (i + .5) * (upperU - lowerU) / 19,
                                        lowerV + (j + .5) * (upperV - lowerV) / 19};
                    if (chart.contains(uv)) {
                        const auto value = chart.sample(uv);
                        require(std::isfinite(value.position.x),
                                "Resampled chart domain has no valid surface lift.");
                        ++samples;
                    }
                }
            require(samples > 0, "A nonplanar chart has no represented interior.");
        }
}

void sameTagGeometricCreases() {
    SurfaceMesh source;
    source.vertices = {{.1, 0, 0}, {-.1, 0, 0}, {0, .1, 0}, {0, -.1, 0}, {0, 0, .1}, {0, 0, -.1}};
    for (const auto triangle : std::array<std::array<std::uint32_t, 3>, 8>{{{4, 0, 2},
                                                                            {4, 2, 1},
                                                                            {4, 1, 3},
                                                                            {4, 3, 0},
                                                                            {5, 2, 0},
                                                                            {5, 1, 2},
                                                                            {5, 3, 1},
                                                                            {5, 0, 3}}})
        source.triangles.push_back({triangle, 1, BoundaryCondition::wall});
    auto options = settings();
    options.targetEdgeLength = .03;
    options.adaptive.enabled = false;
    const auto atlas = makeSurfaceAtlas(source, options);
    require(atlas.charts.size() == 8,
            "Real geometric creases were bridged inside one physical tag.");
    for (const auto &chart : atlas.charts)
        require(chart.faceId == 1 && chart.planar && chart.reference.size() == 1,
                "Preserving a geometric crease changed physical tags or planes.");
    referenceAgreement(atlas);
}
void planarHole() {
    const auto source = ring(96, false);
    const auto atlas = makeSurfaceAtlas(source, settings());
    require(atlas.charts.size() == 1, "The planar ring was unnecessarily cut into charts.");
    const auto &chart = atlas.charts.front();
    require(chart.planar && chart.loops.size() == 2, "The planar hole boundary was lost.");
    const auto uv = [&](Vec3 p) {
        const auto d = sub(p, chart.planeOrigin);
        return ChartPoint{dot(d, chart.planeU), dot(d, chart.planeV)};
    };
    require(!chart.contains(uv({0, 0, 0})), "A planar hole was filled.");
    require(chart.contains(uv({.02, 0, 0})), "The annulus interior was excluded.");
    std::size_t nodes = 0;
    for (const auto &loop : chart.loops)
        nodes += loop.vertices.size();
    require(nodes < source.vertices.size(),
            "Every source boundary vertex became an output constraint.");
    referenceAgreement(atlas);
    const auto first = chart.sample(uv({.02, 0, 0}));
    require(std::abs(dot(first.du, first.du) - 1) < 1e-12 &&
                std::abs(dot(first.du, first.dv)) < 1e-12,
            "The exact planar chart is not orthonormal.");
}
void cylinderSeam() {
    auto source = ring(32, true);
    source.fixedVertices = {0, 32};
    const auto atlas = makeSurfaceAtlas(source, settings());
    require(atlas.charts.size() == 1 && atlas.charts.front().loops.size() == 1,
            "The cylinder was not cut into one disk chart.");
    require(!atlas.charts.front().planar, "A cylinder was flattened by planar projection.");
    require(atlas.fixedVertices.size() == 2, "Fixed cylinder knots were lost.");
    std::map<std::uint32_t, unsigned> visits;
    for (const auto &loop : atlas.charts.front().loops)
        for (const auto &node : loop.vertices)
            ++visits[node.atlasVertex];
    bool shared = false;
    for (const auto &[id, n] : visits) {
        require(id < atlas.vertices.size(), "Invalid seam node.");
        shared = shared || n > 1;
    }
    require(shared, "The two cylinder seam copies have different global node IDs.");
    referenceAgreement(atlas);
    for (const auto &t : atlas.charts.front().reference) {
        const ChartPoint uv{(t.uv[0].u + t.uv[1].u + t.uv[2].u) / 3,
                            (t.uv[0].v + t.uv[1].v + t.uv[2].v) / 3};
        const auto sample = atlas.charts.front().sample(uv);
        require(std::abs(dot(sample.du, sample.du) - 1) < 1e-8 &&
                    std::abs(dot(sample.dv, sample.dv) - 1) < 1e-8 &&
                    std::abs(dot(sample.du, sample.dv)) < 1e-8,
                "An exactly developable cylinder did not receive an isometric chart.");
    }
    const auto second = makeSurfaceAtlas(source, settings());
    require(second.vertices.size() == atlas.vertices.size(), "The atlas is nondeterministic.");
    for (std::size_t i = 0; i < atlas.vertices.size(); ++i)
        require(norm(sub(atlas.vertices[i], second.vertices[i])) == 0,
                "Boundary registry order is nondeterministic.");
}
void torusTopology(unsigned around = 24, unsigned tube = 12) {
    SurfaceMesh source;
    for (unsigned i = 0; i < around; ++i)
        for (unsigned j = 0; j < tube; ++j) {
            const double u = 2 * std::numbers::pi * i / around, v = 2 * std::numbers::pi * j / tube;
            source.vertices.push_back({(.03 + .008 * std::cos(v)) * std::cos(u),
                                       (.03 + .008 * std::cos(v)) * std::sin(u),
                                       .008 * std::sin(v)});
        }
    for (unsigned i = 0; i < around; ++i)
        for (unsigned j = 0; j < tube; ++j) {
            const unsigned a = i * tube + j, b = ((i + 1) % around) * tube + j,
                           c = ((i + 1) % around) * tube + (j + 1) % tube,
                           d = i * tube + (j + 1) % tube;
            source.triangles.push_back({{a, b, c}, 3, BoundaryCondition::wall});
            source.triangles.push_back({{a, c, d}, 3, BoundaryCondition::wall});
        }
    const auto atlas = makeSurfaceAtlas(source, settings());
    require(atlas.charts.size() > 1 && atlas.charts.size() < source.triangles.size() / 4,
            "The torus was not divided into coarse geometric charts.");
    referenceAgreement(atlas);
    std::size_t triangles = 0;
    std::map<std::array<std::uint32_t, 2>, unsigned> sharedEdges;
    for (const auto &chart : atlas.charts) {
        require(chart.faceId == 3, "A computational sector became a physical face.");
        triangles += chart.reference.size();
        std::set<std::pair<double, double>> vertices;
        std::set<std::array<std::pair<double, double>, 2>> edges;
        for (const auto &t : chart.reference) {
            const ChartPoint uv{(t.uv[0].u + t.uv[1].u + t.uv[2].u) / 3,
                                (t.uv[0].v + t.uv[1].v + t.uv[2].v) / 3};
            const auto sample = chart.sample(uv);
            const double g00 = dot(sample.du, sample.du), g01 = dot(sample.du, sample.dv),
                         g11 = dot(sample.dv, sample.dv);
            const double maximum = (g00 + g11 + std::hypot(g00 - g11, 2 * g01)) / 2;
            require(maximum / std::sqrt(g00 * g11 - g01 * g01) <= 2 + 1e-8,
                    "A normal-cone chart exceeds its distortion bound.");
            for (unsigned i = 0; i < 3; ++i) {
                const std::pair<double, double> a{t.uv[i].u, t.uv[i].v},
                    b{t.uv[(i + 1) % 3].u, t.uv[(i + 1) % 3].v};
                vertices.insert(a);
                edges.insert({std::min(a, b), std::max(a, b)});
            }
        }
        require(static_cast<std::int64_t>(vertices.size()) -
                        static_cast<std::int64_t>(edges.size()) +
                        static_cast<std::int64_t>(chart.reference.size()) ==
                    2 - static_cast<std::int64_t>(chart.loops.size()),
                "A projected torus chart has the wrong Euler characteristic.");
        for (const auto &loop : chart.loops)
            for (std::size_t i = 0; i < loop.vertices.size(); ++i) {
                const auto a = loop.vertices[i].atlasVertex,
                           b = loop.vertices[(i + 1) % loop.vertices.size()].atlasVertex;
                ++sharedEdges[{std::min(a, b), std::max(a, b)}];
            }
    }
    require(triangles >= source.triangles.size(), "Chart partitioning lost reference triangles.");
    for (const auto &[edge, count] : sharedEdges) {
        (void)edge;
        require(count == 2, "Torus computational cuts do not conform.");
    }
}

void shortExtrudedAnnulus() {
    // A generic thin, elliptic tube has two physical rims and no axial crease.
    // Alternate the source diagonals to expose a face-tree cut that winds
    // around the tube instead of opening it along a generator.
    constexpr unsigned count = 64;
    constexpr double height = .003;
    for (const bool reverseDiagonals : {false, true}) {
        SurfaceMesh source;
        for (const double z : {0., .0025, height})
            for (unsigned i = 0; i < count; ++i) {
                const double phi = 2 * std::numbers::pi * i / count + .17;
                source.vertices.push_back({.007 * std::cos(phi), .0175 * std::sin(phi), z});
            }
        for (unsigned layer = 0; layer < 2; ++layer)
            for (unsigned i = 0; i < count; ++i) {
                const unsigned a = layer * count + i, b = layer * count + (i + 1) % count,
                               c = a + count, d = b + count;
                if (((i + layer) % 2 == 0) != reverseDiagonals) {
                    source.triangles.push_back({{a, b, d}, 7, BoundaryCondition::wall});
                    source.triangles.push_back({{a, d, c}, 7, BoundaryCondition::wall});
                } else {
                    source.triangles.push_back({{a, b, c}, 7, BoundaryCondition::wall});
                    source.triangles.push_back({{b, d, c}, 7, BoundaryCondition::wall});
                }
            }
        auto options = settings();
        options.quality.minimumAngleDegrees = 20;
        options.quality.maximumAspectRatio = 4;
        const auto atlas = makeSurfaceAtlas(source, options);
        require(atlas.charts.size() == 1 && atlas.charts.front().loops.size() == 1,
                "A short annulus was not opened into one disk.");
        const auto &chart = atlas.charts.front();
        require(chart.faceId == 7 && chart.reference.size() == source.triangles.size(),
                "Opening the annulus changed its physical face or reference triangles.");
        std::map<std::uint32_t, unsigned> visits;
        const auto &loop = chart.loops.front().vertices;
        for (const auto &node : loop)
            ++visits[node.atlasVertex];
        double seamLength = 0;
        unsigned seamSegments = 0;
        for (std::size_t i = 0; i < loop.size(); ++i) {
            const auto &a = loop[i], &b = loop[(i + 1) % loop.size()];
            if (visits.at(a.atlasVertex) == 2 && visits.at(b.atlasVertex) == 2) {
                const auto p = atlas.vertices[a.atlasVertex], q = atlas.vertices[b.atlasVertex];
                require(std::hypot(p.x - q.x, p.y - q.y) < 1e-12,
                        "The artificial annulus seam winds diagonally around a physical rim.");
                seamLength += norm(sub(p, q));
                ++seamSegments;
            }
            const auto previous = loop[(i + loop.size() - 1) % loop.size()].uv;
            const auto current = a.uv, next = b.uv;
            const double u0 = previous.u - current.u, v0 = previous.v - current.v;
            const double u1 = next.u - current.u, v1 = next.v - current.v;
            const double angle = std::atan2(std::abs(u0 * v1 - v0 * u1), u0 * u1 + v0 * v1);
            require(angle >= std::numbers::pi / 2 - 1e-8,
                    "A short annulus acquired an acute artificial disk corner.");
        }
        require(seamSegments >= 2 && std::abs(seamLength - 2 * height) < 1e-12,
                "The annulus cut is not the shortest shared boundary-to-boundary path.");
        referenceAgreement(atlas);
        const auto meshed = meshFrontalChart(chart, atlas.vertices, options);
        require(!meshed.triangles.empty(), "The short annulus did not produce a quality-checked mesh.");
    }
}
void fixedPatchAndPhysicalSeams() {
    SurfaceMesh source;
    source.vertices = {{0, 0, 0},   {.02, 0, 0},   {.02, .02, 0},   {0, .02, 0},
                       {0, 0, .02}, {.02, 0, .02}, {.02, .02, .02}, {0, .02, .02}};
    source.triangles = {
        {{0, 2, 1}, 1, BoundaryCondition::driven}, {{0, 3, 2}, 1, BoundaryCondition::driven},
        {{4, 5, 6}, 2, BoundaryCondition::wall},   {{4, 6, 7}, 2, BoundaryCondition::wall},
        {{0, 1, 5}, 3, BoundaryCondition::wall},   {{0, 5, 4}, 3, BoundaryCondition::wall},
        {{1, 2, 6}, 4, BoundaryCondition::wall},   {{1, 6, 5}, 4, BoundaryCondition::wall},
        {{2, 3, 7}, 5, BoundaryCondition::wall},   {{2, 7, 6}, 5, BoundaryCondition::wall},
        {{3, 0, 4}, 6, BoundaryCondition::wall},   {{3, 4, 7}, 6, BoundaryCondition::wall}};
    source.fixedVertices = {0, 1, 2, 3};
    const auto atlas = makeSurfaceAtlas(source, settings());
    require(atlas.charts.size() == 5 && atlas.fixedTriangles.size() == 2,
            "The fixed driven patch was triangulated again.");
    require(atlas.fixedVertices.size() == 4, "Fixed positions were lost.");
    for (std::size_t i = 0; i < 4; ++i)
        require(norm(sub(atlas.vertices[atlas.fixedVertices[i]], source.vertices[i])) == 0,
                "A fixed vertex moved.");
    std::map<std::array<std::uint32_t, 2>, unsigned> edges;
    auto record = [&](std::uint32_t a, std::uint32_t b) {
        ++edges[{std::min(a, b), std::max(a, b)}];
    };
    for (const auto &chart : atlas.charts)
        for (const auto &loop : chart.loops)
            for (std::size_t i = 0; i < loop.vertices.size(); ++i)
                record(loop.vertices[i].atlasVertex,
                       loop.vertices[(i + 1) % loop.vertices.size()].atlasVertex);
    for (const auto &t : atlas.fixedTriangles)
        for (unsigned i = 0; i < 3; ++i)
            record(t.vertices[i], t.vertices[(i + 1) % 3]);
    for (const auto &[edge, count] : edges) {
        (void)edge;
        require(count == 2, "Shared physical boundary samples are nonconforming.");
    }
    referenceAgreement(atlas);
}
void invalidAndCancelled() {
    auto source = ring(8, false);
    for (auto &triangle : source.triangles)
        triangle.faceId = 0;
    require(makeSurfaceAtlas(source, settings()).charts.front().faceId == 0,
            "A neutral physical face ID of zero was rejected or changed.");
    source = ring(8, false);
    source.triangles[0].vertices[0] = 999;
    rejects([&] { makeSurfaceAtlas(source, settings()); });
    source = ring(8, false);
    source.fixedVertices = {0, 0};
    rejects([&] { makeSurfaceAtlas(source, settings()); });
    source = ring(8, false);
    source.triangles[0].boundary = BoundaryCondition::driven;
    rejects([&] { makeSurfaceAtlas(source, settings()); });
    source = ring(8, false);
    source.vertices[0].x = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { makeSurfaceAtlas(source, settings()); });
    source = ring(8, false);
    std::stop_source stop;
    stop.request_stop();
    rejects([&] { makeSurfaceAtlas(source, settings(), stop.get_token()); });
}
void curvedHoleReferenceLens() {
    constexpr unsigned segments = 126;
    SurfaceMesh source;
    for (double radius : {.004, .003, .00225})
        for (unsigned i = 0; i < segments; ++i) {
            const double angle = 2 * std::numbers::pi * i / segments;
            source.vertices.push_back(
                {radius * std::cos(angle), radius * std::sin(angle), radius - .00225});
        }
    for (unsigned layer = 0; layer < 2; ++layer)
        for (unsigned i = 0; i < segments; ++i) {
            const unsigned a = layer * segments + i, b = layer * segments + (i + 1) % segments,
                           c = a + segments, d = b + segments;
            source.triangles.push_back({{a, b, c}, 9, BoundaryCondition::wall});
            source.triangles.push_back({{b, d, c}, 9, BoundaryCondition::wall});
        }
    auto options = settings();
    options.targetEdgeLength = .002;
    const auto atlas = makeSurfaceAtlas(source, options);
    require(atlas.charts.size() == 1 && !atlas.charts.front().planar,
            "The curved hole lost its actual nonplanar chart.");
    const auto &chart = atlas.charts.front();
    std::size_t boundary = 0;
    for (const auto &loop : chart.loops)
        boundary += loop.vertices.size();
    require(boundary < segments,
            "Dense source hole tessellation became mandatory output boundary nodes.");
    require(chart.reference.size() > source.triangles.size(),
            "The concave curved-hole domain has no explicit reference extension.");
    for (std::size_t i = source.triangles.size(); i < chart.reference.size(); ++i) {
        const auto &t = chart.reference[i];
        const Vec3 point{(t.positions[0].x + t.positions[1].x + t.positions[2].x) / 3,
                         (t.positions[0].y + t.positions[1].y + t.positions[2].y) / 3,
                         (t.positions[0].z + t.positions[1].z + t.positions[2].z) / 3};
        double distance = 1;
        for (unsigned j = 0; j < segments; ++j) {
            const auto a = source.vertices[2 * segments + j],
                       b = source.vertices[2 * segments + (j + 1) % segments], edge = sub(b, a);
            const double fraction =
                std::clamp(dot(sub(point, a), edge) / dot(edge, edge), 0.0, 1.0);
            distance = std::min(distance,
                                norm(sub(point, {a.x + fraction * edge.x, a.y + fraction * edge.y,
                                                 a.z + fraction * edge.z})));
        }
        require(distance <= options.adaptive.geometricTolerance * (1 + 1e-10),
                "A reference lens exceeded its physical source-curve deviation bound.");
    }
    referenceAgreement(atlas);
}

void sharedCurveRefinement() {
    const auto source = ring(96, false);
    auto atlas = makeSurfaceAtlas(source, settings());
    std::array<std::uint32_t, 2> segment{};
    bool found = false;
    for (const auto &loop : atlas.charts.front().loops)
        for (std::size_t i = 0; i < loop.vertices.size(); ++i) {
            const auto a = loop.vertices[i].atlasVertex,
                       b = loop.vertices[(i + 1) % loop.vertices.size()].atlasVertex;
            const auto p = atlas.vertices[a], q = atlas.vertices[b];
            if (norm(p) > .029 && norm(q) > .029 &&
                .03 - norm({(p.x + q.x) / 2, (p.y + q.y) / 2, 0}) > 1e-5) {
                segment = {a, b};
                found = true;
                break;
            }
        }
    require(found, "The curve fixture has no genuinely curved resampled segment.");
    const auto before = atlas.vertices.size();
    const auto a = atlas.vertices[segment[0]], b = atlas.vertices[segment[1]];
    refineSurfaceAtlas(atlas, std::span{&segment, 1});
    require(atlas.vertices.size() == before + 1,
            "Shared curve refinement did not add exactly one global vertex.");
    const auto point = atlas.vertices.back();
    require(norm(sub(point, {(a.x + b.x) / 2, (a.y + b.y) / 2, 0})) > 1e-5,
            "A curved segment was split at its straight chord midpoint.");
    double distance = 1;
    for (unsigned layer = 0; layer < 2; ++layer)
        for (unsigned i = 0; i < 96; ++i) {
            const auto p = source.vertices[layer * 96 + i],
                       q = source.vertices[layer * 96 + (i + 1) % 96], edge = sub(q, p);
            const double t = std::clamp(dot(sub(point, p), edge) / dot(edge, edge), 0.0, 1.0);
            distance = std::min(
                distance, norm(sub(point, {p.x + t * edge.x, p.y + t * edge.y, p.z + t * edge.z})));
        }
    require(distance < 1e-12, "The new shared vertex left its original source curve.");
    const auto committed = atlas.curves;
    rejects([&] { refineSurfaceAtlas(atlas, std::span{&segment, 1}); });
    require(atlas.vertices.size() == before + 1 && atlas.curves == committed,
            "A stale parent request partially changed the atlas.");
    const std::array<std::array<std::uint32_t, 2>, 2> duplicates{
        {{segment[0], static_cast<std::uint32_t>(before)},
         {static_cast<std::uint32_t>(before), segment[0]}}};
    refineSurfaceAtlas(atlas, duplicates);
    require(atlas.vertices.size() == before + 2,
            "Duplicate curve requests were split more than once.");
    referenceAgreement(atlas);

    auto cylinder = makeSurfaceAtlas(ring(32, true), settings());
    std::map<std::array<std::uint32_t, 2>, unsigned> uses;
    for (const auto &chart : cylinder.charts)
        for (const auto &loop : chart.loops)
            for (std::size_t i = 0; i < loop.vertices.size(); ++i) {
                const auto x = loop.vertices[i].atlasVertex,
                           y = loop.vertices[(i + 1) % loop.vertices.size()].atlasVertex;
                ++uses[{std::min(x, y), std::max(x, y)}];
            }
    found = false;
    for (const auto &[edge, count] : uses)
        if (count == 2) {
            segment = edge;
            found = true;
            break;
        }
    require(found, "The periodic fixture has no shared seam segment.");
    const auto seamId = static_cast<std::uint32_t>(cylinder.vertices.size());
    const auto changed = refineSurfaceAtlas(cylinder, std::span{&segment, 1});
    require(changed == std::vector<std::size_t>{0},
            "Two UV sides of one periodic chart were reported as multiple changed charts.");
    unsigned occurrences = 0;
    for (const auto &chart : cylinder.charts)
        for (const auto &loop : chart.loops)
            for (const auto &vertex : loop.vertices)
                if (vertex.atlasVertex == seamId) {
                    ++occurrences;
                    require(norm(sub(chart.sample(vertex.uv).position, cylinder.vertices[seamId])) <
                                1e-11,
                            "A periodic seam UV copy disagrees with its original curve position.");
                }
    require(occurrences == 2, "A periodic seam split was applied on only one UV side.");
    referenceAgreement(cylinder);
}

void localRefinementPreservesUnchangedCharts() {
    SurfaceMesh source;
    constexpr double width = .016;
    for (unsigned column = 0; column < 4; ++column)
        for (unsigned row = 0; row < 2; ++row)
            source.vertices.push_back({column * width, row * width, 0});
    for (std::uint32_t patch = 0; patch < 3; ++patch) {
        const auto a = 2 * patch, b = a + 2;
        const auto boundary = patch == 0 ? BoundaryCondition::driven : BoundaryCondition::wall;
        source.triangles.push_back({{a, b, b + 1}, patch + 1, boundary});
        source.triangles.push_back({{a, b + 1, a + 1}, patch + 1, boundary});
    }
    source.fixedVertices = {0};
    const auto options = settings();
    auto atlas = makeSurfaceAtlas(source, options);
    require(atlas.charts.size() == 3, "The three physical strip patches were not preserved.");
    std::map<std::uint32_t, std::size_t> index;
    for (std::size_t i = 0; i < atlas.charts.size(); ++i)
        require(index.emplace(atlas.charts[i].faceId, i).second,
                "A planar strip patch was unnecessarily split.");
    const auto before = atlas;
    const auto third = index.at(3);
    const auto originalFront = meshFrontalChart(atlas.charts[third], atlas.vertices, options);
    std::array<std::uint32_t, 2> seam{};
    bool found = false;
    for (const auto &loop : atlas.charts[index.at(1)].loops)
        for (std::size_t i = 0; i < loop.vertices.size(); ++i) {
            const auto a = loop.vertices[i].atlasVertex,
                       b = loop.vertices[(i + 1) % loop.vertices.size()].atlasVertex;
            if (std::abs(atlas.vertices[a].x - width) < 1e-12 &&
                std::abs(atlas.vertices[b].x - width) < 1e-12) {
                seam = {a, b};
                found = true;
                break;
            }
        }
    require(found, "Missing physical seam between the first two patches.");
    const std::array<std::array<std::uint32_t, 2>, 3> requests{{seam, {seam[1], seam[0]}, seam}};
    const auto changed = refineSurfaceAtlas(atlas, requests);
    std::vector<std::size_t> expected{index.at(1), index.at(2)};
    std::sort(expected.begin(), expected.end());
    require(changed == expected,
            "Refinement must report exactly both incident charts, sorted and unique.");
    require(
        atlas.charts.size() == before.charts.size() &&
            atlas.vertices.size() == before.vertices.size() + 1 &&
            std::equal(before.vertices.begin(), before.vertices.end(), atlas.vertices.begin()) &&
            atlas.fixedVertices == before.fixedVertices &&
            atlas.fixedTriangles == before.fixedTriangles,
        "Refinement renumbered existing atlas IDs or changed fixed source geometry.");
    const auto appended = static_cast<std::uint32_t>(before.vertices.size());
    for (std::size_t i = 0; i < atlas.charts.size(); ++i) {
        unchangedChartGeometry(before.charts[i], atlas.charts[i]);
        unsigned copies = 0;
        for (const auto &loop : atlas.charts[i].loops)
            for (const auto &node : loop.vertices)
                if (node.atlasVertex == appended) {
                    ++copies;
                    require(norm(sub(atlas.charts[i].sample(node.uv).position,
                                     atlas.vertices[appended])) < 1e-12,
                            "A shared split has inconsistent UV lifts on its incident charts.");
                }
        require(copies == (i == third ? 0u : 1u),
                "A split was omitted from a neighbor or inserted into an unrelated chart.");
    }
    require(sameBoundaries(before.charts[third], atlas.charts[third]),
            "A local seam request changed the third patch boundary.");
    const auto unchangedFront = meshFrontalChart(atlas.charts[third], atlas.vertices, options);
    require(
        !originalFront.triangles.empty() &&
            physicalTriangles(originalFront) == physicalTriangles(unchangedFront),
        "Appending unrelated shared vertices changed an unchanged chart's physical front mesh.");
    const auto committedVertices = atlas.vertices;
    require(refineSurfaceAtlas(atlas, {}).empty() && atlas.vertices == committedVertices &&
                sameBoundaries(before.charts[third], atlas.charts[third]),
            "An empty refinement request invalidated unchanged chart data.");
}

class CollapsingBoundaryGeometry final : public SurfaceGeometry {
  public:
    std::span<const BoundaryCorner> boundaryCorners() const noexcept override {
        return {};
    }
    SurfaceProjection projectSurface(std::uint32_t, Vec3 p) const override {
        return {p, {0, 0, 1}};
    }
    SurfaceProjection projectSupportingSurface(std::uint32_t face, Vec3 p) const override {
        return projectSurface(face, p);
    }
    Vec3 projectBoundary(std::span<const std::uint32_t> faces, Vec3 p) const override {
        require(faces.size() == 2, "Missing shared physical curve in the collapse fixture.");
        // Simulate a defective nearest-curve solver at a new interior point.
        // Existing samples stay correct; the source itself does not change.
        if (p.x > 0 && p.x < .0004 && std::abs(p.x - p.y) < 1e-15)
            return {};
        return p;
    }
    Vec3 projectJunction(std::span<const std::uint32_t>, Vec3) const override {
        throw std::runtime_error("This planar fixture has no three-face junction.");
    }
};

void collapsedCurveRefinementIsAtomic() {
    SurfaceMesh source;
    source.vertices = {{0, 0, 0}, {.002, 0, 0}, {.002, .002, 0}, {0, .002, 0}};
    source.triangles = {{{0, 1, 2}, 1, BoundaryCondition::driven},
                        {{0, 2, 3}, 2, BoundaryCondition::wall}};
    auto options = settings();
    options.targetEdgeLength = .001;
    options.adaptive.enabled = false;
    options.geometry = std::make_shared<CollapsingBoundaryGeometry>();
    auto atlas = makeSurfaceAtlas(source, options);
    std::array<std::uint32_t, 2> segment{};
    bool found = false;
    for (const auto &chart : atlas.charts)
        for (const auto &loop : chart.loops)
            for (std::size_t i = 0; i < loop.vertices.size(); ++i) {
                const auto a = loop.vertices[i].atlasVertex,
                           b = loop.vertices[(i + 1) % loop.vertices.size()].atlasVertex;
                const auto p = atlas.vertices[a], q = atlas.vertices[b];
                require(p != q, "Collapse fixture starts with an invalid physical boundary.");
                const Vec3 midpoint{(p.x + q.x) / 2, (p.y + q.y) / 2, 0};
                if (std::min(p.x, q.x) == 0 && p.x == p.y && q.x == q.y && midpoint.x > 0 &&
                    midpoint.x < .0004) {
                    segment = {a, b};
                    found = true;
                }
            }
    require(found, "Missing genuinely interior curve point in the collapse fixture.");
    const auto before = atlas;
    std::string failure;
    try {
        (void)refineSurfaceAtlas(atlas, std::span{&segment, 1});
    } catch (const std::exception &error) {
        failure = error.what();
    }
    require(failure.find("collapsed") != std::string::npos,
            "A collapsed physical midpoint was accepted or misdiagnosed.");
    require(atlas.vertices == before.vertices && atlas.curves == before.curves &&
                atlas.fixedVertices == before.fixedVertices &&
                atlas.fixedTriangles == before.fixedTriangles &&
                atlas.charts.size() == before.charts.size(),
            "A collapsed split partially committed new atlas vertices or registry data.");
    for (std::size_t i = 0; i < atlas.charts.size(); ++i) {
        require(sameBoundaries(atlas.charts[i], before.charts[i]),
                "A collapsed split changed an incident boundary loop.");
        unchangedChartGeometry(atlas.charts[i], before.charts[i]);
    }
}

void infeasibleCurveRefinementIsAtomic() {
    SurfaceMesh source;
    source.vertices = {{0, 0, 0}, {.002, 0, 0}, {.002, .002, 0}, {0, .002, 0}};
    source.triangles = {{{0, 1, 2}, 1, BoundaryCondition::driven},
                        {{0, 2, 3}, 2, BoundaryCondition::wall}};
    source.fixedVertices = {0, 1, 2};
    auto atlas = makeSurfaceAtlas(source, settings());
    require(atlas.fixedTriangles.size() == 1, "Missing fixed-triangle refinement fixture.");
    const auto vertices = atlas.vertices.size();
    const auto registry = atlas.curves;
    const std::array<std::uint32_t, 2> fixedEdge{atlas.fixedVertices[0], atlas.fixedVertices[2]};
    rejects([&] { refineSurfaceAtlas(atlas, std::span{&fixedEdge, 1}); });
    require(atlas.vertices.size() == vertices && atlas.curves == registry,
            "An infeasible fixed-edge request modified the atlas.");
    const std::array<std::uint32_t, 2> unknown{0, 999999};
    rejects([&] { refineSurfaceAtlas(atlas, std::span{&unknown, 1}); });
    std::stop_source stop;
    stop.request_stop();
    rejects([&] { refineSurfaceAtlas(atlas, std::span{&fixedEdge, 1}, stop.get_token()); });
    require(atlas.vertices.size() == vertices && atlas.curves == registry,
            "Rejected refinement changed fixed geometry.");
}

void actualDraftCorners() {
    std::vector<DraftParameters> examples(5);
    examples[0].profile.throatAngle = .2;
    examples[0].profile.k = .7;
    examples[1].entry = {.032, .008, 8};
    examples[1].profile.throatAngle = .12;
    examples[1].profile.k = .7;
    examples[1].verticalProfile.throatAngle = .21;
    examples[2].guideEnabled = true;
    examples[2].guide.type = CurveType::superformula;
    examples[2].guide.superformula.n1 = 4;
    examples[2].guide.superformula.n2 = 4;
    examples[2].guide.superformula.n3 = 4;
    examples[3] = examples[2];
    examples[3].guide.rotation = .3;
    examples[4].mouthMorph.enabled = true;
    for (const auto &parameters : examples) {
        const auto source = DraftGeometry(parameters).previewMesh();
        const auto atlas = makeSurfaceAtlas(source, settings());
        const auto triangleKey = [](std::array<Vec3, 3> positions) {
            std::array<std::array<double, 3>, 3> points;
            for (std::size_t i = 0; i < 3; ++i)
                points[i] = {positions[i].x, positions[i].y, positions[i].z};
            std::sort(points.begin(), points.end());
            return points;
        };
        using ReferenceKey = std::pair<std::uint32_t, std::array<std::array<double, 3>, 3>>;
        std::map<ReferenceKey, std::size_t> references;
        for (const auto &chart : atlas.charts)
            for (const auto &triangle : chart.reference)
                ++references[{chart.faceId, triangleKey(triangle.positions)}];
        for (const auto &triangle : atlas.fixedTriangles)
            ++references[{triangle.faceId, triangleKey({atlas.vertices[triangle.vertices[0]],
                                                        atlas.vertices[triangle.vertices[1]],
                                                        atlas.vertices[triangle.vertices[2]]})}];
        for (const auto &triangle : source.triangles) {
            auto &count =
                references[{triangle.faceId, triangleKey({source.vertices[triangle.vertices[0]],
                                                          source.vertices[triangle.vertices[1]],
                                                          source.vertices[triangle.vertices[2]]})}];
            require(count > 0, "Chart conditioning changed or lost original physical geometry.");
            --count;
        }
        for (const auto &chart : atlas.charts)
            for (const auto &loop : chart.loops)
                for (std::size_t i = 0; i < loop.vertices.size(); ++i) {
                    const auto previous =
                        loop.vertices[(i + loop.vertices.size() - 1) % loop.vertices.size()];
                    const auto current = loop.vertices[i],
                               next = loop.vertices[(i + 1) % loop.vertices.size()];
                    const auto a = sub(atlas.vertices[previous.atlasVertex],
                                       atlas.vertices[current.atlasVertex]);
                    const auto b =
                        sub(atlas.vertices[next.atlasVertex], atlas.vertices[current.atlasVertex]);
                    const Vec3 product{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                                       a.x * b.y - a.y * b.x};
                    const double turn =
                        (current.uv.u - previous.uv.u) * (next.uv.v - current.uv.v) -
                        (current.uv.v - previous.uv.v) * (next.uv.u - current.uv.u);
                    double angle = std::atan2(norm(product), dot(a, b)) * 180 / std::numbers::pi;
                    if (turn < 0)
                        angle = 360 - angle;
                    require(angle >= 20, "A real Draft parameter combination contains an "
                                         "infeasible acute chart-domain corner.");
                }
    }
}

void requiredCurveKnotsRespectSpacing() {
    for (bool rectangular : {false, true}) {
        DraftParameters parameters;
        if (rectangular) {
            parameters.entry = {.032, .008, 8};
            parameters.profile.throatAngle = .12;
            parameters.profile.k = .7;
            parameters.verticalProfile.throatAngle = .21;
        }
        const auto source = DraftGeometry(parameters).previewMesh();
        const auto atlas = makeSurfaceAtlas(source, settings());
        for (auto fixed : source.fixedVertices) {
            const auto position = source.vertices[fixed];
            require(std::any_of(atlas.fixedVertices.begin(), atlas.fixedVertices.end(),
                                [&](auto id) {
                                    const auto p = atlas.vertices[id];
                                    return p.x == position.x && p.y == position.y &&
                                           p.z == position.z;
                                }),
                    "Shared-curve resampling moved a fixed source vertex.");
        }
        unsigned checked = 0;
        for (const auto &chart : atlas.charts)
            if (chart.faceId == static_cast<std::uint32_t>(DraftFace::innerWall))
                for (const auto &loop : chart.loops)
                    for (std::size_t i = 0; i < loop.vertices.size(); ++i) {
                        const auto a = loop.vertices[i],
                                   b = loop.vertices[(i + 1) % loop.vertices.size()];
                        const auto p = atlas.vertices[a.atlasVertex],
                                   q = atlas.vertices[b.atlasVertex];
                        // Interior meridians have no fixed entry/mouth constraints.
                        // Adding a required reference knot must redistribute nearby
                        // samples instead of leaving micrometre-scale remnants.
                        if (std::min(p.z, q.z) < .02 ||
                            std::max(p.z, q.z) > parameters.length - .005)
                            continue;
                        const double h = std::min(chart.sample(a.uv).targetEdgeLength,
                                                  chart.sample(b.uv).targetEdgeLength);
                        require(norm(sub(q, p)) >= .05 * h,
                                "A required curve knot created an artificial tiny segment.");
                        ++checked;
                    }
        require(checked > 20, "Missing interior shared-curve spacing coverage.");
        referenceAgreement(atlas);
    }
}
} // namespace
int main() {
    try {
        // The almost straight arc produces a very thin, valid reference lens.
        // A midpoint on its boundary must not fail solely through division by
        // its tiny area; an actually exterior query must still be rejected.
        for (const double scale : {0.001, 1., 1000.}) {
            SurfaceChart chart;
            ChartReferenceTriangle lens;
            lens.uv = {{{0.00062233113789500905, 0.002376103127350788},
                        {0.0012260737817320429, 0.0008230191810444092},
                        {0.0009188962105152246, 0.0016132111085377647}}};
            for (unsigned i = 0; i < 3; ++i) {
                lens.uv[i].u *= scale;
                lens.uv[i].v *= scale;
                lens.positions[i] = {lens.uv[i].u, lens.uv[i].v, 0};
                lens.sizes[i] = .008 * scale;
            }
            chart.reference.push_back(lens);
            const auto a = lens.uv[0], b = lens.uv[1];
            const ChartPoint middle{(a.u + b.u) / 2, (a.v + b.v) / 2};
            const auto value = chart.sample(middle);
            require(norm(sub(value.position, {middle.u, middle.v, 0})) < 1e-17 * scale,
                    "A reference-lens boundary query changed its physical position.");
            const double length = std::hypot(b.u - a.u, b.v - a.v);
            const ChartPoint outside{middle.u + 1e-12 * scale * (b.v - a.v) / length,
                                     middle.v - 1e-12 * scale * (b.u - a.u) / length};
            rejects([&] { (void)chart.sample(outside); });
        }
        planarHole();
        cylinderSeam();
        shortExtrudedAnnulus();
        torusTopology();
        // Coarse torus creases split one physical face into smooth components.
        // Small normal regions must not inspect labels across those cuts.
        torusTopology(8, 6);
        fixedPatchAndPhysicalSeams();
        sameTagGeometricCreases();
        invalidAndCancelled();
        sharedCurveRefinement();
        localRefinementPreservesUnchangedCharts();
        curvedHoleReferenceLens();
        infeasibleCurveRefinementIsAtomic();
        collapsedCurveRefinementIsAtomic();
        actualDraftCorners();
        requiredCurveKnotsRespectSpacing();
        std::cout << "Surface atlas: planar holes, periodic seams, genus, fixed patches, "
                     "metric/lifting, deterministic registry and rejection PASS.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
