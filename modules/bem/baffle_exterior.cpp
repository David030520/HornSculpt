#include "baffle_exterior.hpp"
#include "galerkin/near_projection.hpp"
#include <algorithm>
#include <cmath>
#include <hornsculpt/meshing.hpp>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>

namespace hornsculpt::bem_detail {
namespace {
using Complex = std::complex<double>;
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 sub(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 point(const galerkin::Panel &p, const double *b) {
    return {p.vertices[0].x + b[1] * (p.vertices[1].x - p.vertices[0].x) +
                b[2] * (p.vertices[2].x - p.vertices[0].x),
            p.vertices[0].y + b[1] * (p.vertices[1].y - p.vertices[0].y) +
                b[2] * (p.vertices[2].y - p.vertices[0].y),
            p.vertices[0].z + b[1] * (p.vertices[1].z - p.vertices[0].z) +
                b[2] * (p.vertices[2].z - p.vertices[0].z)};
}
Complex singleLayer(const galerkin::Panel &a, const galerkin::Panel &b,
                    const galerkin::PairQuadratureSettings &settings, double k,
                    std::stop_token stop) {
    const auto plan = galerkin::makePairQuadrature(a, b, settings, k, stop);
    Complex value{};
    const auto add = [&](const galerkin::PairQuadraturePoint &q) {
        auto delta = sub(point(a, q.test), point(b, q.trial));
        double r = std::hypot(delta.x, delta.y, delta.z);
        if (!(r > 0) || !std::isfinite(r))
            throw std::runtime_error("Exterior Rayleigh quadrature is singular.");
        const double weight = a.area * b.area * q.weight / (std::numbers::pi * r);
        value += weight * std::exp(Complex{0, -k * r});
    };
    for (const auto &q : plan.points)
        add(q);
    for (const auto &patch : plan.patches) {
        const auto rule = galerkin::triangleRule(patch.degree == 14 ? 14 : 8);
        for (unsigned i = 0;
             i < galerkin::nearPatchPointCount(patch, static_cast<unsigned>(rule.size())); ++i)
            add(galerkin::mapNearPatchPoint(rule.data(), static_cast<unsigned>(rule.size()), patch,
                                            i, a, b));
    }
    if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
        throw std::runtime_error("Exterior Rayleigh integral is not finite.");
    return value;
}
} // namespace
SurfaceMesh appendBaffleExterior(const SurfaceMesh &cavity, const SurfaceMesh &exterior) {
    auto out = cavity;
    const auto offset = static_cast<std::uint32_t>(out.vertices.size());
    if (exterior.vertices.size() > UINT32_MAX - offset)
        throw std::invalid_argument("Exterior patch index capacity exceeded.");
    out.vertices.insert(out.vertices.end(), exterior.vertices.begin(), exterior.vertices.end());
    for (auto t : exterior.triangles) {
        for (auto &i : t.vertices)
            i += offset;
        out.triangles.push_back(t);
    }
    return out;
}
void validateBaffleExteriorMesh(const SurfaceMesh &cavity, const SurfaceMesh &exterior, Vec3 origin,
                                Vec3 normal) {
    if (exterior.vertices.empty() && exterior.triangles.empty())
        return;
    auto quality = inspectMesh(exterior);
    if (quality.degenerateTriangles || quality.nonManifoldEdges || quality.inconsistentEdges ||
        quality.selfIntersectionPairs || !quality.boundaryEdges)
        throw std::invalid_argument("Exterior baffle sources require consistently oriented, "
                                    "nonoverlapping planar open patches.");
    // inspectMesh's vertex statistic requires a CLOSED link cycle. An open source
    // patch permits one connected link path at each rim vertex as well.
    std::vector<std::vector<std::array<unsigned, 2>>> links(exterior.vertices.size());
    for (const auto &t : exterior.triangles)
        for (unsigned i = 0; i < 3; ++i)
            links[t.vertices[i]].push_back({t.vertices[(i + 1) % 3], t.vertices[(i + 2) % 3]});
    for (const auto &link : links) {
        if (link.empty())
            continue;
        std::map<unsigned, std::vector<unsigned>> graph;
        for (auto edge : link) {
            graph[edge[0]].push_back(edge[1]);
            graph[edge[1]].push_back(edge[0]);
        }
        unsigned endpoints = 0;
        for (const auto &[id, neighbors] : graph) {
            if (neighbors.size() == 1)
                ++endpoints;
            else if (neighbors.size() != 2)
                throw std::invalid_argument("Exterior patch has a nonmanifold vertex.");
        }
        std::set<unsigned> seen;
        std::vector<unsigned> pending{graph.begin()->first};
        while (!pending.empty()) {
            auto id = pending.back();
            pending.pop_back();
            if (!seen.insert(id).second)
                continue;
            for (auto next : graph.at(id))
                pending.push_back(next);
        }
        if ((endpoints != 0 && endpoints != 2) || seen.size() != graph.size())
            throw std::invalid_argument("Exterior patch has a disconnected vertex link.");
    }
    auto space = galerkin::buildSpace(exterior);
    std::set<unsigned> ids;
    for (const auto &t : cavity.triangles)
        ids.insert(t.faceId);
    double extent = 0;
    for (auto p : exterior.vertices)
        extent = std::max(extent, std::hypot(p.x - origin.x, p.y - origin.y, p.z - origin.z));
    double tol = 128 * std::numeric_limits<double>::epsilon() * std::max(extent, 1.);
    for (auto p : exterior.vertices)
        if (std::abs(dot(sub(p, origin), normal)) > tol)
            throw std::invalid_argument("Exterior source vertices must lie on the baffle plane.");
    for (const auto &t : space.panels)
        if (!t.faceId || ids.contains(t.faceId) || t.boundary != BoundaryCondition::driven ||
            dot(t.normal, normal) < 1 - 1e-12)
            throw std::invalid_argument(
                "Exterior sources need unique Driven face IDs and outward plane normals.");
    const auto cavitySpace = cavity.triangles.empty() ? galerkin::Space{} : galerkin::buildSpace(cavity);
    for (const auto &a : space.panels)
        for (const auto &b : cavitySpace.panels)
            if (galerkin::panelDistance(a, b) <= tol)
                throw std::invalid_argument(
                    "Exterior source patches must be disjoint from the cavity boundary.");
}
BaffleExteriorOperator baffleExteriorOperator(const galerkin::Space &cavity,
                                              std::span<const unsigned> apertures,
                                              const galerkin::Space &exterior, const BemNumerics &n,
                                              double k, std::stop_token stop) {
    BaffleExteriorOperator out;
    const auto ne = exterior.panels.size(), na = apertures.size();
    out.cross.resize(ne * na);
    out.self.resize(ne * ne);
    galerkin::PairQuadratureSettings qs{n.quadratureOrder, n.singularOrder, n.nearDistanceRatio,
                                        n.nearRelativeTolerance, n.nearMaximumDepth};
    auto ext = exterior.panels;
    // Separate spaces must never accidentally share topological vertex IDs.
    for (auto &p : ext)
        for (auto &id : p.nodes)
            id += static_cast<unsigned>(cavity.vertices.size());
    for (std::size_t e = 0; e < ne; ++e) {
        if (stop.stop_requested())
            throw std::runtime_error("Exterior Rayleigh assembly cancelled.");
        for (std::size_t a = 0; a < na; ++a)
            out.cross[e * na + a] =
                singleLayer(ext[e], cavity.panels.at(apertures[a]), qs, k, stop) / ext[e].area;
        for (std::size_t j = 0; j <= e; ++j) {
            auto v = singleLayer(ext[e], ext[j], qs, k, stop);
            out.self[e * ne + j] = v / ext[e].area;
            out.self[j * ne + e] = v / ext[j].area;
        }
    }
    return out;
}
} // namespace hornsculpt::bem_detail
