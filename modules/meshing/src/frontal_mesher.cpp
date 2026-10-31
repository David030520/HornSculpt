#include "frontal_delaunay.hpp"
#include <hornsculpt/meshing.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>

namespace hornsculpt {
namespace {
void checkStop(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Meshing cancelled.");
}

using Conditions = std::map<std::uint32_t, BoundaryCondition>;
Conditions validateInput(const SurfaceMesh &source, const MeshingOptions &options) {
    if (!std::isfinite(options.targetEdgeLength) || options.targetEdgeLength <= 0 ||
        !std::isfinite(options.quality.minimumAngleDegrees) ||
        options.quality.minimumAngleDegrees <= 0 || options.quality.minimumAngleDegrees >= 60 ||
        !std::isfinite(options.quality.maximumAspectRatio) ||
        options.quality.maximumAspectRatio < 2 / std::sqrt(3.) ||
        !std::isfinite(options.adaptive.minimumEdgeLength) ||
        options.adaptive.minimumEdgeLength <= 0 ||
        options.adaptive.minimumEdgeLength > options.targetEdgeLength ||
        !std::isfinite(options.adaptive.geometricTolerance) ||
        options.adaptive.geometricTolerance <= 0 || !std::isfinite(options.adaptive.growthRate) ||
        options.adaptive.growthRate <= 0 || options.adaptive.growthRate > 1)
        throw std::invalid_argument("Invalid meshing size or quality limits.");
    std::set<std::uint32_t> fixed;
    for (const auto vertex : source.fixedVertices)
        if (vertex >= source.vertices.size() || !fixed.insert(vertex).second)
            throw std::invalid_argument("Fixed vertex indices must be unique and in range.");
    Conditions conditions;
    for (const auto &triangle : source.triangles) {
        const auto [it, inserted] = conditions.emplace(triangle.faceId, triangle.boundary);
        if (!inserted && it->second != triangle.boundary)
            throw std::invalid_argument("A face ID has conflicting boundary conditions.");
    }
    for (const auto &[face, size] : options.faceTargetEdgeLengths)
        if (!conditions.contains(face) || !std::isfinite(size) ||
            size < options.adaptive.minimumEdgeLength || size > options.targetEdgeLength)
            throw std::invalid_argument("Local face targets require an existing face ID and a size "
                                        "within the global bounds.");
    return conditions;
}

using SharedSegment = std::array<std::uint32_t, 2>;
struct ChartMeshingPass {
    struct Failure {
        std::size_t chart{};
        std::optional<double> violationEnergy;
        std::string reason;
    };
    std::vector<SharedSegment> refinements;
    std::vector<Failure> failures;
    std::uint32_t firstRefinementFace{};
    std::string firstRefinementReason;
};
ChartMeshingPass meshCharts(const detail::SurfaceAtlas &atlas, const MeshingOptions &options,
                            std::vector<detail::FrontalChartMesh> &meshes,
                            std::vector<std::uint8_t> &dirty, std::stop_token externalStop) {
    ChartMeshingPass pass;
    std::vector<std::size_t> pending;
    for (std::size_t i = 0; i < dirty.size(); ++i)
        if (dirty[i])
            pending.push_back(i);
    if (pending.empty())
        return pass;
    std::stop_source cancellation;
    std::stop_callback externalCancellation(externalStop, [&] { cancellation.request_stop(); });
    std::atomic<std::size_t> next{};
    std::vector<std::vector<SharedSegment>> refinements(meshes.size());
    std::vector<std::string> refinementReasons(meshes.size());
    std::vector<std::optional<double>> violationEnergies(meshes.size());
    std::exception_ptr failure;
    std::mutex failureMutex;
    const auto count = std::min<std::size_t>(
        pending.size(), std::clamp(std::thread::hardware_concurrency(), 1u, 4u));
    {
        std::vector<std::jthread> workers;
        workers.reserve(count);
        for (std::size_t worker = 0; worker < count; ++worker)
            workers.emplace_back([&] {
                while (!cancellation.stop_requested()) {
                    const auto task = next.fetch_add(1);
                    if (task >= pending.size())
                        return;
                    const auto index = pending[task];
                    try {
                        meshes[index] = detail::meshFrontalChart(
                            atlas.charts[index], atlas.vertices, options, cancellation.get_token());
                        dirty[index] = 0;
                    } catch (const detail::BoundaryRefinementRequired &request) {
                        // Every pending chart finishes before a shared curve changes.
                        // Requests are collected in chart order, independently of worker timing.
                        refinements[index].assign(request.segments().begin(),
                                                  request.segments().end());
                        refinementReasons[index] = request.what();
                        violationEnergies[index] = request.violationEnergy();
                    } catch (...) {
                        // Preserve the chart's actual error; sibling jobs observe cancellation.
                        const std::scoped_lock lock(failureMutex);
                        if (!failure)
                            failure = std::current_exception();
                        cancellation.request_stop();
                    }
                }
            });
    }
    checkStop(externalStop);
    if (failure)
        std::rethrow_exception(failure);
    std::set<SharedSegment> unique;
    for (std::size_t chart = 0; chart < refinements.size(); ++chart) {
        const auto &requests = refinements[chart];
        if (requests.empty() && meshes[chart].triangles.empty())
            throw std::runtime_error("A chart returned neither a mesh nor a refinement request.");
        if (!requests.empty())
            pass.failures.push_back({chart, violationEnergies[chart], refinementReasons[chart]});
        if (unique.empty() && !requests.empty()) {
            pass.firstRefinementFace = atlas.charts[chart].faceId;
            pass.firstRefinementReason = std::move(refinementReasons[chart]);
        }
        for (auto segment : requests) {
            if (segment[1] < segment[0])
                std::swap(segment[0], segment[1]);
            unique.insert(segment);
        }
    }
    pass.refinements.assign(unique.begin(), unique.end());
    return pass;
}

std::uint32_t appendVertex(SurfaceMesh &mesh, Vec3 position) {
    if (mesh.vertices.size() >= std::numeric_limits<std::uint32_t>::max())
        throw std::runtime_error("The generated mesh exceeds its vertex index range.");
    const auto index = static_cast<std::uint32_t>(mesh.vertices.size());
    mesh.vertices.push_back(position);
    return index;
}

void compactVertices(SurfaceMesh &mesh) {
    constexpr auto absent = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> indices(mesh.vertices.size(), absent);
    std::vector<Vec3> vertices;
    for (auto &triangle : mesh.triangles)
        for (auto &vertex : triangle.vertices) {
            if (vertex >= mesh.vertices.size())
                throw std::runtime_error("The surface atlas returned an invalid vertex index.");
            auto &mapped = indices[vertex];
            if (mapped == absent) {
                mapped = static_cast<std::uint32_t>(vertices.size());
                vertices.push_back(mesh.vertices[vertex]);
            }
            vertex = mapped;
        }
    for (auto &vertex : mesh.fixedVertices) {
        if (vertex >= indices.size() || indices[vertex] == absent)
            throw std::runtime_error("Frontal-Delaunay meshing lost a fixed source vertex.");
        vertex = indices[vertex];
    }
    mesh.vertices = std::move(vertices);
}

void validateFixedPositions(const SurfaceMesh &source, const SurfaceMesh &result) {
    using Coordinate = std::tuple<double, double, double>;
    const auto coordinates = [](const SurfaceMesh &mesh) {
        std::multiset<Coordinate> positions;
        for (const auto index : mesh.fixedVertices) {
            const auto p = mesh.vertices.at(index);
            positions.emplace(p.x, p.y, p.z);
        }
        return positions;
    };
    if (coordinates(source) != coordinates(result))
        throw std::runtime_error("Frontal-Delaunay meshing changed fixed source positions.");
}
} // namespace

MeshingResult remeshIsotropic(const SurfaceMesh &source, const MeshingOptions &options,
                              std::stop_token stop) {
    checkStop(stop);
    const auto conditions = validateInput(source, options);
    const auto input = inspectMesh(source);
    if (!input.watertight() || input.signedVolume <= 0)
        throw std::invalid_argument(
            "Meshing requires a closed, consistently outward-oriented nondegenerate boundary.");
    checkStop(stop);
    auto atlas = detail::makeSurfaceAtlas(source, options, stop);
    detail::validateAtlasBoundaryCorners(atlas, options, stop);
    std::vector<detail::FrontalChartMesh> chartMeshes(atlas.charts.size());
    // Shared refinement only appends atlas vertices and changes adjacent loops.
    // Successful unrelated fronts remain valid, including their shared IDs.
    std::vector<std::uint8_t> dirtyCharts(atlas.charts.size(), 1);
    // The number of global batches does not measure work: a real thin trim
    // can need many dyadic levels while unrelated completed charts are cached.
    // Bound the actual added shared vertices and require progress every batch.
    constexpr std::size_t maximumSharedBoundaryInsertions = 100000;
    struct Progress {
        double bestEnergy{std::numeric_limits<double>::infinity()};
        unsigned stagnantPasses{};
    };
    std::vector<Progress> progress(atlas.charts.size());
    std::size_t boundaryInsertions = 0;
    for (;;) {
        auto pass = meshCharts(atlas, options, chartMeshes, dirtyCharts, stop);
        if (pass.refinements.empty()) {
            break;
        }
        for (std::size_t chart = 0; chart < dirtyCharts.size(); ++chart)
            if (!dirtyCharts[chart])
                progress[chart] = {};
        for (const auto &failure : pass.failures) {
            auto &state = progress[failure.chart];
            // Boundary-size preprocessing has no complete chart energy yet.
            if (!failure.violationEnergy) {
                state = {};
                continue;
            }
            if (*failure.violationEnergy < .99 * state.bestEnergy) {
                state.bestEnergy = *failure.violationEnergy;
                state.stagnantPasses = 0;
            } else if (++state.stagnantPasses >= 6) {
                throw std::runtime_error("Frontal-Delaunay shared-boundary refinement stopped "
                                         "without quality progress on face " +
                                         std::to_string(atlas.charts[failure.chart].faceId) +
                                         ". Six retries did not reduce the dimensionless violation "
                                         "energy by one percent; "
                                         "no mesh was accepted. " +
                                         failure.reason);
            }
        }
        if (pass.refinements.size() > maximumSharedBoundaryInsertions - boundaryInsertions)
            throw std::runtime_error(
                "Frontal-Delaunay meshing reached its shared-boundary refinement budget on face " +
                std::to_string(pass.firstRefinementFace) +
                ". "
                "Review the local size and geometric constraints; no mesh was accepted. " +
                pass.firstRefinementReason);
        checkStop(stop);
        const auto previousVertices = atlas.vertices.size();
        const auto changed = detail::refineSurfaceAtlas(atlas, pass.refinements, stop);
        if (changed.empty() || atlas.vertices.size() <= previousVertices)
            throw std::runtime_error("Shared-boundary refinement made no geometric progress.");
        boundaryInsertions += atlas.vertices.size() - previousVertices;
        for (const auto index : changed)
            dirtyCharts.at(index) = 1;
    }
    MeshingResult result;
    result.mesh.vertices = std::move(atlas.vertices);
    result.mesh.triangles = std::move(atlas.fixedTriangles);
    result.mesh.fixedVertices = std::move(atlas.fixedVertices);
    const auto atlasVertexCount = result.mesh.vertices.size();
    constexpr auto interior = std::numeric_limits<std::uint32_t>::max();
    for (std::size_t index = 0; index < chartMeshes.size(); ++index) {
        checkStop(stop);
        const auto &chart = atlas.charts[index];
        const auto &mesh = chartMeshes[index];
        std::vector<std::uint32_t> mapping;
        mapping.reserve(mesh.vertices.size());
        for (const auto &vertex : mesh.vertices) {
            if (vertex.atlasVertex == interior) {
                mapping.push_back(appendVertex(result.mesh, vertex.position));
            } else {
                if (vertex.atlasVertex >= atlasVertexCount)
                    throw std::runtime_error("A chart refers to an unknown shared atlas vertex.");
                mapping.push_back(vertex.atlasVertex);
            }
        }
        for (const auto triangle : mesh.triangles)
            result.mesh.triangles.push_back(
                {{mapping.at(triangle[0]), mapping.at(triangle[1]), mapping.at(triangle[2])},
                 chart.faceId,
                 chart.boundary});
    }
    compactVertices(result.mesh);
    validateFixedPositions(source, result.mesh);
    Conditions outputConditions;
    for (const auto &triangle : result.mesh.triangles) {
        const auto [it, inserted] = outputConditions.emplace(triangle.faceId, triangle.boundary);
        if (!inserted && it->second != triangle.boundary)
            throw std::runtime_error("Frontal-Delaunay meshing mixed boundary conditions.");
    }
    if (outputConditions != conditions)
        throw std::runtime_error("Frontal-Delaunay meshing lost a physical surface patch.");
    checkStop(stop);
    result.quality = inspectMesh(result.mesh);
    checkStop(stop);
    const auto &q = result.quality;
    if (!q.watertight() || q.connectedComponents != input.connectedComponents ||
        q.eulerCharacteristic != input.eulerCharacteristic || q.signedVolume <= 0 ||
        (options.quality.enabled && (q.minimumAngleDegrees < options.quality.minimumAngleDegrees ||
                                     q.maximumAspectRatio > options.quality.maximumAspectRatio))) {
        std::ostringstream message;
        message << "Frontal-Delaunay boundary rejected: minimum angle " << q.minimumAngleDegrees
                << " deg, maximum aspect ratio " << q.maximumAspectRatio << ", open edges "
                << q.boundaryEdges << ", nonmanifold edges " << q.nonManifoldEdges
                << ", inconsistent edges " << q.inconsistentEdges << ", intersecting pairs "
                << q.selfIntersectionPairs << ", nonmanifold vertices " << q.nonManifoldVertices
                << ", degenerate triangles " << q.degenerateTriangles << ", Euler characteristic "
                << q.eulerCharacteristic << " (input " << input.eulerCharacteristic << ").";
        throw std::runtime_error(message.str());
    }
    return result;
}
} // namespace hornsculpt
