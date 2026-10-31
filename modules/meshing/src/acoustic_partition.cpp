#include "frontal_delaunay.hpp"
#include <hornsculpt/acoustic_partition.hpp>
#include <hornsculpt/symmetry_sector.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace hornsculpt {
namespace {
using Edge = std::pair<std::uint32_t, std::uint32_t>;
using Loop = std::vector<std::uint32_t>;
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Teilgebietsvernetzung abgebrochen.");
}
Vec3 sub(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3 mul(Vec3 a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double norm(Vec3 v) {
    return std::hypot(v.x, v.y, v.z);
}
std::uint32_t index32(std::size_t size) {
    if (size >= std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument("Zu viele Teilgebiets-Netzelemente.");
    return static_cast<std::uint32_t>(size);
}
std::vector<Loop> openingLoops(const SurfaceMesh &source, const BemSubdomain &domain) {
    struct Use {
        unsigned count{};
        Edge direction;
    };
    std::map<Edge, Use> edges;
    for (const auto &panel : domain.panels) {
        auto v = source.triangles.at(panel.triangleIndex).vertices;
        // Work in the inward-air orientation, matching the original physical
        // material faces and the convention for newly constructed caps.
        if (panel.normalSign == 1)
            std::swap(v[1], v[2]);
        for (unsigned j = 0; j < 3; ++j) {
            const auto a = v[j], b = v[(j + 1) % 3];
            auto &use = edges[std::minmax(a, b)];
            ++use.count;
            use.direction = {a, b};
        }
    }
    std::map<std::uint32_t, std::uint32_t> next;
    std::set<std::uint32_t> incoming;
    for (const auto &[edge, use] : edges) {
        (void)edge;
        if (use.count > 2)
            throw std::invalid_argument("Teilgebietsrand ist nicht mannigfaltig.");
        if (use.count != 1)
            continue;
        // A cap must traverse the physical boundary in the opposite direction.
        const auto [b, a] = use.direction;
        if (!next.emplace(a, b).second || !incoming.insert(b).second)
            throw std::invalid_argument("Mehrdeutiger Ã–ffnungsrand oder Dreifachanschluss.");
    }
    std::vector<Loop> loops;
    while (!next.empty()) {
        const auto first = next.begin()->first;
        Loop loop;
        auto vertex = first;
        do {
            const auto found = next.find(vertex);
            if (found == next.end())
                throw std::invalid_argument(
                    "TeilgebietsÃ¶ffnung besitzt keinen geschlossenen Rand.");
            loop.push_back(vertex);
            vertex = found->second;
            next.erase(found);
        } while (vertex != first);
        if (loop.size() < 3)
            throw std::invalid_argument("TeilgebietsÃ¶ffnung besitzt weniger als drei Randpunkte.");
        loops.push_back(std::move(loop));
    }
    return loops;
}
std::vector<Edge> loopKey(const Loop &loop) {
    std::vector<Edge> result;
    for (std::size_t i = 0; i < loop.size(); ++i)
        result.push_back(std::minmax(loop[i], loop[(i + 1) % loop.size()]));
    std::sort(result.begin(), result.end());
    return result;
}
void validateBoundary(const SurfaceMesh &source, const BemSubdomain &domain,
                      const MeshQualityLimits &limits) {
    SurfaceMesh boundary;
    std::map<std::uint32_t, std::uint32_t> vertices;
    for (const auto &panel : domain.panels) {
        auto triangle = source.triangles.at(panel.triangleIndex);
        for (auto &v : triangle.vertices) {
            const auto [it, inserted] = vertices.emplace(v, index32(boundary.vertices.size()));
            if (inserted)
                boundary.vertices.push_back(source.vertices.at(v));
            v = it->second;
        }
        if (panel.normalSign == -1)
            std::swap(triangle.vertices[1], triangle.vertices[2]);
        boundary.triangles.push_back(triangle);
    }
    const auto quality = inspectMesh(boundary);
    if (!quality.watertight())
        throw std::invalid_argument(
            "Gebiet '" + domain.name +
            "': Randabschluss, Orientierung oder Schnittfreiheit ungÃ¼ltig.");
    // Signed volumes must be checked per nested component. Summing them can
    // hide an inverted cavity boundary behind a larger outer shell.
    std::vector<std::size_t> parent(boundary.vertices.size());
    std::iota(parent.begin(), parent.end(), 0);
    const auto root = [&](std::size_t v) {
        while (parent[v] != v) {
            parent[v] = parent[parent[v]];
            v = parent[v];
        }
        return v;
    };
    for (const auto &triangle : boundary.triangles) {
        const auto first = root(triangle.vertices[0]);
        for (unsigned j = 1; j < 3; ++j)
            parent[root(triangle.vertices[j])] = first;
    }
    struct Component {
        std::vector<std::size_t> panels;
        Vec3 lower, upper;
    };
    std::vector<Component> components;
    std::map<std::size_t, std::size_t> lookup;
    for (std::size_t t = 0; t < boundary.triangles.size(); ++t) {
        const auto &triangle = boundary.triangles[t];
        const auto [it, added] = lookup.emplace(root(triangle.vertices[0]), components.size());
        if (added) {
            const auto p = boundary.vertices[triangle.vertices[0]];
            components.push_back({{}, p, p});
        }
        auto &c = components[it->second];
        c.panels.push_back(t);
        for (auto v : triangle.vertices) {
            const auto p = boundary.vertices[v];
            c.lower = {std::min(c.lower.x, p.x), std::min(c.lower.y, p.y),
                       std::min(c.lower.z, p.z)};
            c.upper = {std::max(c.upper.x, p.x), std::max(c.upper.y, p.y),
                       std::max(c.upper.z, p.z)};
        }
    }
    struct Sum {
        double value{}, correction{};
        void add(double x) {
            const double n = value + x;
            correction += std::abs(value) >= std::abs(x) ? (value - n) + x : (x - n) + value;
            value = n;
        }
        double total() const {
            return value + correction;
        }
    };
    const auto contains = [&](const Component &c, Vec3 q) {
        if (q.x < c.lower.x || q.x > c.upper.x || q.y < c.lower.y || q.y > c.upper.y ||
            q.z < c.lower.z || q.z > c.upper.z)
            return false;
        Sum angle;
        for (auto t : c.panels) {
            const auto &triangle = boundary.triangles[t];
            std::array<Vec3, 3> ray{};
            for (unsigned j = 0; j < 3; ++j) {
                ray[j] = sub(boundary.vertices[triangle.vertices[j]], q);
                const double length = norm(ray[j]);
                if (!(length > 0) || !std::isfinite(length))
                    throw std::invalid_argument("Gebietsrandkomponenten berÃ¼hren sich.");
                ray[j] = mul(ray[j], 1 / length);
            }
            angle.add(2 * std::atan2(dot(ray[0], cross(ray[1], ray[2])), 1 + dot(ray[0], ray[1]) +
                                                                             dot(ray[1], ray[2]) +
                                                                             dot(ray[2], ray[0])));
        }
        return std::abs(angle.total()) > 2 * std::numbers::pi;
    };
    std::size_t outerCount = 0;
    for (std::size_t i = 0; i < components.size(); ++i) {
        const auto &component = components[i];
        const auto q = boundary.vertices[boundary.triangles[component.panels.front()].vertices[0]];
        unsigned depth = 0;
        for (std::size_t j = 0; j < components.size(); ++j)
            if (i != j && contains(components[j], q))
                ++depth;
        if (!depth)
            ++outerCount;
        // Deeper shells are disconnected air domains, not part of one region.
        if ((domain.kind == BemDomainKind::exterior && depth != 0) ||
            (domain.kind == BemDomainKind::interior && depth > 1))
            throw std::invalid_argument(
                "Gebiet '" + domain.name +
                "': getrennte Luftgebiete durch verschachtelte Randkomponenten.");
        Sum volume;
        for (auto t : component.panels) {
            const auto &triangle = boundary.triangles[t];
            const auto a = sub(boundary.vertices[triangle.vertices[0]], component.lower);
            const auto b = sub(boundary.vertices[triangle.vertices[1]], component.lower);
            const auto c = sub(boundary.vertices[triangle.vertices[2]], component.lower);
            volume.add(dot(a, cross(b, c)) / 6);
        }
        const bool positive = domain.kind == BemDomainKind::interior && depth == 0;
        if (!std::isfinite(volume.total()) || volume.total() == 0 ||
            (volume.total() > 0) != positive)
            throw std::invalid_argument("Gebiet '" + domain.name +
                                        "': Randkomponente ist falsch orientiert.");
    }
    if (domain.kind == BemDomainKind::interior && outerCount != 1)
        throw std::invalid_argument("Gebiet '" + domain.name +
                                    "': getrennte Luftvolumen benÃ¶tigen eigene Gebiete.");
    if (limits.enabled && (quality.minimumAngleDegrees + 1e-9 < limits.minimumAngleDegrees ||
                           quality.maximumAspectRatio > limits.maximumAspectRatio * (1 + 1e-12)))
        throw std::invalid_argument("Gebiet '" + domain.name +
                                    "': Netz unterschreitet die eingestellten QualitÃ¤tsgrenzen.");
}
void validateSymmetricBoundary(const SurfaceMesh &source, const BemSubdomain &domain,
                               const MeshQualityLimits &limits,
                               const OrthogonalSymmetryPlanes &planes) {
    SurfaceMesh sector;
    std::map<std::uint32_t, std::uint32_t> local;
    for (const auto &panel : domain.panels) {
        auto t = source.triangles.at(panel.triangleIndex);
        for (auto &v : t.vertices) {
            const auto [where, added] = local.emplace(v, index32(sector.vertices.size()));
            if (added)
                sector.vertices.push_back(source.vertices.at(v));
            v = where->second;
        }
        if (panel.normalSign == -1)
            std::swap(t.vertices[1], t.vertices[2]);
        sector.triangles.push_back(t);
    }
    (void)validateSymmetrySector(sector, planes);
    // Project only a private sector copy for orbit validation. This has no
    // more triangles than the input and leaves the model snapshot untouched.
    std::vector<unsigned> seam(sector.vertices.size());
    for (std::size_t i = 0; i < sector.vertices.size(); ++i) {
        auto &p = sector.vertices[i];
        const double d1 = dot(sub(p, planes.origin), planes.firstNormal);
        const double d2 = dot(sub(p, planes.origin), planes.secondNormal);
        if (std::abs(d1) <= planes.planeDistanceTolerance) {
            seam[i] |= 1u;
            p = sub(p, mul(planes.firstNormal, d1));
        }
        if (std::abs(d2) <= planes.planeDistanceTolerance) {
            seam[i] |= 2u;
            p = sub(p, mul(planes.secondNormal, d2));
        }
    }
    const auto sectorReport = validateSymmetrySector(sector, planes);
    const auto quality = inspectMesh(sector);
    if (!sectorReport.firstPlaneEdges || !sectorReport.secondPlaneEdges ||
        sectorReport.intersectionAxisEdges)
        throw std::invalid_argument(
            "Symmetrisches Gebiet benÃ¶tigt eine zusammenhÃ¤ngende Randkomponente mit je einem "
            "Nahtanschluss an beiden Ebenen ohne Randkante auf der Schnittachse.");
    if (limits.enabled && (quality.minimumAngleDegrees + 1e-9 < limits.minimumAngleDegrees ||
                           quality.maximumAspectRatio > limits.maximumAspectRatio * (1 + 1e-12)))
        throw std::invalid_argument(
            "Symmetrisches Gebiet verletzt die eingestellten QualitÃ¤tsgrenzen.");
    struct EdgeUse {
        unsigned count{};
    };
    std::map<Edge, EdgeUse> edges;
    for (const auto &triangle : sector.triangles)
        for (unsigned j = 0; j < 3; ++j)
            ++edges[std::minmax(triangle.vertices[j], triangle.vertices[(j + 1) % 3])].count;
    std::vector<unsigned> axisIncident(sector.vertices.size());
    for (const auto &[edge, use] : edges)
        if (use.count == 1) {
            unsigned common = seam[edge.first] & seam[edge.second];
            if (common != 1u && common != 2u)
                throw std::invalid_argument("Mehrdeutige Symmetrienaht.");
            if (seam[edge.first] == 3u)
                axisIncident[edge.first] |= common;
            if (seam[edge.second] == 3u)
                axisIncident[edge.second] |= common;
        }
    for (std::size_t i = 0; i < seam.size(); ++i)
        if (seam[i] == 3u && axisIncident[i] != 3u)
            throw std::invalid_argument(
                "Schnittachsenknoten benÃ¶tigt NahtanschlÃ¼sse an beide Ebenen.");
    // Components may be nested (e.g. the two sides of an air shell).
    // Evaluate each image on demand; never materialize the fourfold mesh.
    std::vector<std::size_t> parent(sector.vertices.size());
    std::iota(parent.begin(), parent.end(), 0);
    const auto root = [&](std::size_t v) {
        while (parent[v] != v) {
            parent[v] = parent[parent[v]];
            v = parent[v];
        }
        return v;
    };
    for (const auto &t : sector.triangles) {
        const auto a = root(t.vertices[0]);
        parent[root(t.vertices[1])] = a;
        parent[root(t.vertices[2])] = a;
    }
    std::map<std::size_t, std::vector<std::size_t>> components;
    for (std::size_t i = 0; i < sector.triangles.size(); ++i)
        components[root(sector.triangles[i].vertices[0])].push_back(i);
    const auto image = [&](Vec3 p, unsigned mask) {
        if (mask & 1u)
            p = sub(p, mul(planes.firstNormal, 2 * dot(sub(p, planes.origin), planes.firstNormal)));
        if (mask & 2u)
            p = sub(p,
                    mul(planes.secondNormal, 2 * dot(sub(p, planes.origin), planes.secondNormal)));
        return p;
    };
    const auto contains = [&](const std::vector<std::size_t> &panels, Vec3 point) {
        double sum = 0;
        for (unsigned mask = 0; mask < 4; ++mask)
            for (auto index : panels) {
                const auto &t = sector.triangles[index];
                Vec3 a = sub(image(sector.vertices[t.vertices[0]], mask), point);
                Vec3 b = sub(image(sector.vertices[t.vertices[1]], mask), point);
                Vec3 c = sub(image(sector.vertices[t.vertices[2]], mask), point);
                const double la = norm(a), lb = norm(b), lc = norm(c);
                if (!(la > 0 && lb > 0 && lc > 0))
                    throw std::invalid_argument("Symmetrische Randkomponenten berÃ¼hren sich.");
                a = mul(a, 1 / la);
                b = mul(b, 1 / lb);
                c = mul(c, 1 / lc);
                const double angle =
                    2 * std::atan2(dot(a, cross(b, c)), 1 + dot(a, b) + dot(b, c) + dot(c, a));
                sum += (mask == 1 || mask == 2) ? -angle : angle;
            }
        return std::abs(sum) > 2 * std::numbers::pi;
    };
    std::size_t outerCount = 0;
    for (const auto &[key, panels] : components) {
        SurfaceMesh component;
        std::map<std::uint32_t, std::uint32_t> vertices;
        for (auto index : panels) {
            auto t = sector.triangles[index];
            for (auto &v : t.vertices) {
                const auto [it, added] = vertices.emplace(v, index32(component.vertices.size()));
                if (added)
                    component.vertices.push_back(sector.vertices[v]);
                v = it->second;
            }
            component.triangles.push_back(t);
        }
        const auto report = validateSymmetrySector(component, planes);
        if (!report.firstPlaneEdges || !report.secondPlaneEdges || report.intersectionAxisEdges)
            throw std::invalid_argument(
                "Jede symmetrische Randkomponente benÃ¶tigt NahtanschlÃ¼sse an beide Ebenen.");
        const auto point = sector.vertices[sector.triangles[panels.front()].vertices[0]];
        unsigned depth = 0;
        for (const auto &[otherKey, otherPanels] : components)
            if (otherKey != key && contains(otherPanels, point))
                ++depth;
        if (!depth)
            ++outerCount;
        if ((domain.kind == BemDomainKind::exterior && depth != 0) ||
            (domain.kind == BemDomainKind::interior && depth > 1))
            throw std::invalid_argument(
                "Symmetrische Randkomponenten bilden getrennte Luftgebiete.");
        double sum = 0, correction = 0, extent = 0;
        for (auto index : panels) {
            const auto &t = sector.triangles[index];
            const auto a = sub(sector.vertices[t.vertices[0]], planes.origin);
            const auto b = sub(sector.vertices[t.vertices[1]], planes.origin);
            const auto c = sub(sector.vertices[t.vertices[2]], planes.origin);
            extent = std::max({extent, norm(a), norm(b), norm(c)});
            const double term = dot(a, cross(b, c)) / 6;
            const double next = sum + term;
            correction +=
                std::abs(sum) >= std::abs(term) ? (sum - next) + term : (term - next) + sum;
            sum = next;
        }
        const double volume = 4 * (sum + correction);
        const bool positive = domain.kind == BemDomainKind::interior && depth == 0;
        if (!std::isfinite(volume) || !(std::abs(volume) > extent * extent * extent * 1e-14) ||
            ((volume > 0) != positive))
            throw std::invalid_argument("Symmetrische Randkomponente ist falsch orientiert.");
    }
    if (domain.kind == BemDomainKind::interior && outerCount != 1)
        throw std::invalid_argument(
            "Getrennte symmetrische Luftvolumen benÃ¶tigen eigene Gebiete.");
    // After seam projection, image interiors occupy disjoint open half-spaces.
    // All shared-plane geometry belongs to the checked original sector.
}

std::vector<std::uint32_t> appendCap(SurfaceMesh &mesh, const Loop &loop, std::uint32_t faceId,
                                     const MeshingOptions &settings, std::stop_token stop) {
    const auto origin = mesh.vertices.at(loop.front());
    Vec3 area{};
    double diameter{};
    for (std::size_t i = 0; i < loop.size(); ++i) {
        const auto a = sub(mesh.vertices.at(loop[i]), origin);
        const auto b = sub(mesh.vertices.at(loop[(i + 1) % loop.size()]), origin);
        area = add(area, cross(a, b));
        diameter = std::max(diameter, norm(a));
    }
    if (!(norm(area) > diameter * diameter * 1e-14))
        throw std::invalid_argument("InterfaceÃ¶ffnung hat keine eindeutige ebene FlÃ¤che.");
    const auto normal = mul(area, 1 / norm(area));
    auto u = sub(mesh.vertices.at(loop[1]), origin);
    u = mul(u, 1 / norm(u));
    const auto v = cross(normal, u);
    // This is a numerical planarity budget, not the geometric modelling
    // tolerance. Rim points are never projected or moved to flatten an opening.
    const double planarBudget =
        std::max(diameter * 1e-10,
                 256 * std::numeric_limits<double>::epsilon() * std::max(norm(origin), diameter));
    detail::SurfaceChart chart;
    chart.faceId = faceId;
    chart.planar = true;
    chart.planeOrigin = origin;
    chart.planeU = u;
    chart.planeV = v;
    detail::SurfaceChartLoop rim;
    for (const auto id : loop) {
        const auto delta = sub(mesh.vertices.at(id), origin);
        if (std::abs(dot(delta, normal)) > planarBudget)
            throw std::invalid_argument("InterfaceÃ¶ffnung ist nicht eben. CAD-SchnittflÃ¤che oder "
                                        "einen ebenen Ã–ffnungsrand wÃ¤hlen.");
        rim.vertices.push_back({{dot(delta, u), dot(delta, v)}, id});
    }
    chart.loops.push_back(std::move(rim));
    detail::ChartReferenceTriangle reference;
    reference.uv = {detail::ChartPoint{-diameter, -diameter},
                    detail::ChartPoint{3 * diameter, -diameter},
                    detail::ChartPoint{-diameter, 3 * diameter}};
    for (unsigned i = 0; i < 3; ++i) {
        reference.positions[i] =
            add(origin, add(mul(u, reference.uv[i].u), mul(v, reference.uv[i].v)));
        reference.sizes[i] = settings.targetEdgeLength;
    }
    chart.reference.push_back(reference);
    auto options = settings;
    options.geometry.reset();
    options.parameterization.reset();
    options.faceTargetEdgeLengths.clear();
    detail::FrontalChartMesh cap;
    try {
        cap = detail::meshFrontalChart(chart, mesh.vertices, options, stop);
    } catch (const detail::BoundaryRefinementRequired &) {
        // The physical mesh is immutable here. Never split only the cap rim.
        throw std::invalid_argument("Interface benÃ¶tigt feineren gemeinsamen Ã–ffnungsrand. "
                                    "Physisches Netz mit kleineren Zielkanten erneut erzeugen.");
    }
    std::vector<std::uint32_t> mapped;
    for (const auto &vertex : cap.vertices) {
        if (vertex.atlasVertex != std::numeric_limits<std::uint32_t>::max()) {
            mapped.push_back(vertex.atlasVertex);
        } else {
            mapped.push_back(index32(mesh.vertices.size()));
            mesh.vertices.push_back(vertex.position);
        }
    }
    std::vector<std::uint32_t> triangles;
    for (auto vertices : cap.triangles) {
        for (auto &vertex : vertices)
            vertex = mapped.at(vertex);
        triangles.push_back(index32(mesh.triangles.size()));
        mesh.triangles.push_back({vertices, faceId, BoundaryCondition::wall});
    }
    if (triangles.empty())
        throw std::invalid_argument("Interface hat keine Dreiecke erzeugt.");
    return triangles;
}
} // namespace

namespace {
void validatePartition(const SurfaceMesh &mesh, const BemSubdomainModel &model,
                       const OrthogonalSymmetryPlanes *symmetry, const ExteriorBaffle *baffle) {
    if (model.domains.empty())
        throw std::invalid_argument("Keine akustischen Gebiete.");
    // The combined complex has intentional three-face edges, but geometric
    // crossings and duplicate triangles remain invalid independently of that.
    const auto combined = inspectMesh(mesh);
    if (combined.degenerateTriangles || combined.selfIntersectionPairs)
        throw std::invalid_argument(
            "Gebietsnetz enthÃ¤lt entartete oder sich schneidende FlÃ¤chen.");
    std::map<std::uint32_t, const BemSubdomain *> domains;
    std::vector<std::map<std::uint32_t, int>> owners(mesh.triangles.size());
    unsigned exteriorCount{};
    for (const auto &domain : model.domains) {
        if (!domain.id || domain.name.empty() || domain.panels.empty() ||
            !domains.emplace(domain.id, &domain).second ||
            (domain.kind != BemDomainKind::interior && domain.kind != BemDomainKind::exterior))
            throw std::invalid_argument("UngÃ¼ltige oder doppelte akustische GebietsidentitÃ¤t.");
        if (domain.kind == BemDomainKind::exterior)
            ++exteriorCount;
        for (const auto &panel : domain.panels) {
            if (panel.triangleIndex >= mesh.triangles.size() ||
                (panel.normalSign != 1 && panel.normalSign != -1) ||
                !owners[panel.triangleIndex].emplace(domain.id, panel.normalSign).second)
                throw std::invalid_argument(
                    "Gebietsseite besitzt ungÃ¼ltige Indizes, Normalen oder doppelte Panels.");
        }
    }
    if (exteriorCount > 1)
        throw std::invalid_argument("Mehrere unbeschrÃ¤nkte AuÃŸenrÃ¤ume sind im aktuellen "
                                    "Gebietsvertrag nicht zulÃ¤ssig.");
    std::set<std::uint32_t> interfaceIds, interfacePanels, interfaceVertices;
    for (const auto &interface : model.interfaces) {
        if (!interface.id || interface.name.empty() || !interfaceIds.insert(interface.id).second ||
            interface.firstDomainId == interface.secondDomainId ||
            !domains.contains(interface.firstDomainId) ||
            !domains.contains(interface.secondDomainId) || interface.triangleIndices.empty())
            throw std::invalid_argument("UngÃ¼ltige InterfaceidentitÃ¤t oder Gebietsreferenz.");
        std::set<std::uint32_t> vertices;
        int firstNormal{};
        for (const auto triangle : interface.triangleIndices) {
            if (triangle >= mesh.triangles.size() || !interfacePanels.insert(triangle).second)
                throw std::invalid_argument("Interface besitzt doppelte oder fehlende Dreiecke.");
            const auto &owner = owners[triangle];
            if (owner.size() != 2 || !owner.contains(interface.firstDomainId) ||
                !owner.contains(interface.secondDomainId) ||
                owner.at(interface.firstDomainId) + owner.at(interface.secondDomainId) != 0 ||
                mesh.triangles[triangle].boundary != BoundaryCondition::wall)
                throw std::invalid_argument(
                    "Interface benÃ¶tigt genau zwei gegenorientierte Seiten "
                    "ohne zusÃ¤tzliche Anregung.");
            const int sign = owner.at(interface.firstDomainId);
            if (firstNormal && sign != firstNormal)
                throw std::invalid_argument(
                    "Interface-Normalen sind nicht einheitlich orientiert.");
            firstNormal = sign;
            for (const auto vertex : mesh.triangles[triangle].vertices)
                vertices.insert(vertex);
        }
        for (const auto vertex : vertices)
            if (!interfaceVertices.insert(vertex).second)
                throw std::invalid_argument(
                    "Interfaces mit gemeinsamen Randknoten benÃ¶tigen einen "
                    "Mehrfachanschlussvertrag.");
    }
    std::set<std::uint32_t> twoSidedPanels;
    for (const auto &wall : model.twoSidedSurfaces) {
        if (wall.firstDomainId == wall.secondDomainId || !domains.contains(wall.firstDomainId) ||
            !domains.contains(wall.secondDomainId) || wall.triangleIndices.empty())
            throw std::invalid_argument("UngÃ¼ltige zweiseitige Wandzuordnung.");
        for (const auto triangle : wall.triangleIndices) {
            if (triangle >= mesh.triangles.size() || interfacePanels.contains(triangle) ||
                !twoSidedPanels.insert(triangle).second)
                throw std::invalid_argument(
                    "Zweiseitige Fläche besitzt doppelte oder Interface-Dreiecke.");
            const auto &owner = owners[triangle];
            if (owner.size() != 2 || !owner.contains(wall.firstDomainId) ||
                !owner.contains(wall.secondDomainId) ||
                owner.at(wall.firstDomainId) + owner.at(wall.secondDomainId) != 0 ||
                (mesh.triangles[triangle].boundary != BoundaryCondition::wall &&
                 mesh.triangles[triangle].boundary != BoundaryCondition::driven))
                throw std::invalid_argument("Zweiseitige FlÃ¤che benÃ¶tigt zwei gegenorientierte "
                                            "Seiten und Wall oder Driven.");
        }
    }
    for (std::size_t i = 0; i < owners.size(); ++i)
        if (owners[i].empty() || ((interfacePanels.contains(static_cast<std::uint32_t>(i)) ||
                                   twoSidedPanels.contains(static_cast<std::uint32_t>(i)))
                                      ? owners[i].size() != 2
                                      : owners[i].size() != 1))
            throw std::invalid_argument("Netzdreieck fehlt in der Gebietszuordnung oder hat "
                                        "mehrere physische EigentÃ¼mer.");
    const MeshQualityLimits topologyOnly{false, 20, 4};
    for (const auto &domain : model.domains) {
        if (baffle && domain.kind == BemDomainKind::exterior) {
            SurfaceMesh boundary;
            std::map<std::uint32_t, std::uint32_t> local;
            double volume = 0;
            bool allOnBaffle = true;
            for (const auto &panel : domain.panels) {
                auto triangle = mesh.triangles.at(panel.triangleIndex);
                for (auto &v : triangle.vertices) {
                    const auto [it, added] = local.emplace(v, index32(boundary.vertices.size()));
                    if (added) {
                        const auto point = mesh.vertices.at(v);
                        boundary.vertices.push_back(point);
                        allOnBaffle = allOnBaffle &&
                                      std::abs(dot(sub(point, baffle->origin), baffle->normal)) <=
                                          baffle->planeDistanceTolerance;
                    }
                    v = it->second;
                }
                if (panel.normalSign == -1)
                    std::swap(triangle.vertices[1], triangle.vertices[2]);
                const auto a = sub(boundary.vertices[triangle.vertices[0]], baffle->origin);
                const auto b = sub(boundary.vertices[triangle.vertices[1]], baffle->origin);
                const auto c = sub(boundary.vertices[triangle.vertices[2]], baffle->origin);
                volume += dot(a, cross(b, c)) / 6;
                boundary.triangles.push_back(triangle);
            }
            if (symmetry)
                (void)validateSymmetrySector(boundary, *symmetry, baffle);
            else
                validateBaffleBoundary(boundary, *baffle);
            if (!std::isfinite(volume) || (!allOnBaffle && !(volume < 0)))
                throw std::invalid_argument("Außenrand ist falsch orientiert oder besitzt kein "
                                            "Volumen über der Schallwand.");
        } else if (symmetry)
            validateSymmetricBoundary(mesh, domain, topologyOnly, *symmetry);
        else
            validateBoundary(mesh, domain, topologyOnly);
    }
}
} // namespace
void validateAcousticPartition(const SurfaceMesh &mesh, const BemSubdomainModel &model,
                               const ExteriorBaffle *baffle) {
    validatePartition(mesh, model, nullptr, baffle);
}
void validateSymmetryAcousticPartition(const SurfaceMesh &mesh, const BemSubdomainModel &model,
                                       const OrthogonalSymmetryPlanes &planes,
                                       const ExteriorBaffle *baffle) {
    validatePartition(mesh, model, &planes, baffle);
}

AcousticPartition partitionAcousticSurface(const SurfaceMesh &physical,
                                           std::span<const AcousticRegionDefinition> regions,
                                           const MeshingOptions &options, std::stop_token stop,
                                           const OrthogonalSymmetryPlanes *symmetry,
                                           const ExteriorBaffle *baffle) {
    cancelled(stop);
    if (regions.empty())
        throw std::invalid_argument("Keine akustischen Innengebiete definiert.");
    if (!std::isfinite(options.targetEdgeLength) || options.targetEdgeLength <= 0)
        throw std::invalid_argument("UngÃ¼ltige Interface-Zielkante.");
    const bool explicitInterfaces =
        std::any_of(regions.begin(), regions.end(),
                    [](const auto &region) { return !region.interfaceFaceIds.empty(); });
    if (symmetry && !explicitInterfaces)
        throw std::invalid_argument("Symmetrische Gebiete benÃ¶tigen explizite InterfaceflÃ¤chen; "
                                    "automatische Kappen werden nicht erzeugt.");
    const auto sourceQuality = inspectMesh(physical);
    if ((!explicitInterfaces && !sourceQuality.watertight()) || sourceQuality.degenerateTriangles ||
        sourceQuality.selfIntersectionPairs)
        throw std::invalid_argument("Physisches Quellnetz ist nicht geschlossen und schnittfrei; "
                                    "Interfaces reparieren keine ModelllÃ¼cken.");
    AcousticPartition result{physical, {}};
    std::set<std::uint32_t> domainIds, availableFaces;
    std::map<std::uint32_t, std::size_t> assignments;
    std::map<std::uint32_t, std::vector<std::size_t>> interfaceAssignments;
    std::map<std::uint32_t, std::vector<std::size_t>> twoSidedAssignments;
    std::uint32_t lastFace{};
    for (const auto &triangle : physical.triangles) {
        availableFaces.insert(triangle.faceId);
        lastFace = std::max(lastFace, triangle.faceId);
    }
    for (const auto &region : regions) {
        if (region.needsAssignmentReview)
            throw std::invalid_argument(
                "Gebiet '" + region.name +
                "': FlÃ¤chenzuordnung nach Quellen- oder TopologieÃ¤nderung prÃ¼fen.");
        if (!region.id || region.name.empty() ||
            (region.faceIds.empty() && region.interfaceFaceIds.empty() &&
             region.twoSidedFaceIds.empty()) ||
            !domainIds.insert(region.id).second)
            throw std::invalid_argument(
                "Gebiete benÃ¶tigen eindeutige positive IDs, Namen und physische FlÃ¤chen.");
        const auto index = result.model.domains.size();
        result.model.domains.push_back({region.id, region.name, BemDomainKind::interior, {}});
        for (const auto face : region.faceIds)
            if (!availableFaces.contains(face) || !assignments.emplace(face, index).second)
                throw std::invalid_argument("GebietsflÃ¤che fehlt oder ist mehrfach zugeordnet.");
        std::set<std::uint32_t> uniqueSurfaces;
        for (auto face : region.twoSidedFaceIds) {
            if (!availableFaces.contains(face) || !uniqueSurfaces.insert(face).second)
                throw std::invalid_argument(
                    "Zweiseitige WandflÃ¤che fehlt oder ist mehrfach zugeordnet.");
            auto &surfaceOwners = twoSidedAssignments[face];
            surfaceOwners.push_back(index);
            if (surfaceOwners.size() > 2)
                throw std::invalid_argument(
                    "Zweiseitige Fläche besitzt mehr als zwei Gebietsseiten.");
        }
        std::set<std::uint32_t> uniqueInterfaces;
        for (const auto face : region.interfaceFaceIds) {
            if (!availableFaces.contains(face) || !uniqueInterfaces.insert(face).second)
                throw std::invalid_argument(
                    "InterfaceflÃ¤che fehlt oder ist im Gebiet mehrfach zugeordnet.");
            auto &owners = interfaceAssignments[face];
            owners.push_back(index);
            if (owners.size() > 2)
                throw std::invalid_argument(
                    "InterfaceflÃ¤che besitzt mehr als zwei Gebietsseiten.");
        }
    }
    for (const auto &[face, owners] : interfaceAssignments) {
        (void)owners;
        if (assignments.contains(face) || twoSidedAssignments.contains(face))
            throw std::invalid_argument(
                "Eine InterfaceflÃ¤che darf nicht zugleich Wall oder Driven eines Gebiets sein.");
    }
    for (const auto &[face, owners] : twoSidedAssignments) {
        (void)owners;
        if (assignments.contains(face))
            throw std::invalid_argument(
                "Zweiseitige Fläche darf nicht zugleich einseitig zugeordnet sein.");
    }
    std::uint32_t exteriorId = 1;
    while (domainIds.contains(exteriorId)) {
        if (exteriorId == std::numeric_limits<std::uint32_t>::max())
            throw std::invalid_argument("Keine freie Gebiets-ID.");
        ++exteriorId;
    }
    const auto exteriorIndex = result.model.domains.size();
    result.model.domains.push_back({exteriorId, "AuÃŸenraum", BemDomainKind::exterior, {}});
    for (std::size_t i = 0; i < physical.triangles.size(); ++i) {
        if (interfaceAssignments.contains(physical.triangles[i].faceId) ||
            twoSidedAssignments.contains(physical.triangles[i].faceId))
            continue;
        const auto found = assignments.find(physical.triangles[i].faceId);
        const auto domain = found == assignments.end() ? exteriorIndex : found->second;
        result.model.domains[domain].panels.push_back({index32(i), -1});
    }
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> surfacePairs;
    for (std::size_t i = 0; i < physical.triangles.size(); ++i) {
        const auto found = twoSidedAssignments.find(physical.triangles[i].faceId);
        if (found == twoSidedAssignments.end())
            continue;
        if (physical.triangles[i].boundary != BoundaryCondition::wall &&
            physical.triangles[i].boundary != BoundaryCondition::driven)
            throw std::invalid_argument("Zweiseitige FlÃ¤che benÃ¶tigt Wall oder Driven.");
        const auto &surfaceOwners = found->second;
        const auto first = surfaceOwners.front(),
                   second = surfaceOwners.size() == 2 ? surfaceOwners.back() : exteriorIndex;
        const auto pair = std::make_pair(first, second);
        auto [where, inserted] = surfacePairs.emplace(pair, result.model.twoSidedSurfaces.size());
        if (inserted)
            result.model.twoSidedSurfaces.push_back(
                {result.model.domains[first].id, result.model.domains[second].id, {}});
        result.model.twoSidedSurfaces[where->second].triangleIndices.push_back(index32(i));
        result.model.domains[first].panels.push_back({index32(i), -1});
        result.model.domains[second].panels.push_back({index32(i), 1});
    }
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> explicitPairs;
    for (std::size_t i = 0; i < physical.triangles.size(); ++i) {
        const auto found = interfaceAssignments.find(physical.triangles[i].faceId);
        if (found == interfaceAssignments.end())
            continue;
        const auto &owners = found->second;
        const auto first = owners.front(),
                   second = owners.size() == 2 ? owners.back() : exteriorIndex;
        const auto pair = std::make_pair(first, second);
        auto [where, inserted] = explicitPairs.emplace(pair, result.model.interfaces.size());
        if (inserted)
            result.model.interfaces.push_back({index32(result.model.interfaces.size() + 1),
                                               result.model.domains[first].name + " / Interface",
                                               result.model.domains[first].id,
                                               result.model.domains[second].id,
                                               {}});
        result.model.interfaces[where->second].triangleIndices.push_back(index32(i));
        result.model.domains[first].panels.push_back({index32(i), -1});
        // Exterior sides are appended together with generated caps below.
        if (second != exteriorIndex)
            result.model.domains[second].panels.push_back({index32(i), 1});
    }
    std::map<std::vector<Edge>, std::size_t> interfaces;
    for (std::size_t d = 0; d < regions.size(); ++d) {
        cancelled(stop);
        if (symmetry)
            continue;
        const auto loops = openingLoops(physical, result.model.domains[d]);
        if (explicitInterfaces && !loops.empty())
            throw std::invalid_argument(
                "Gebiet mit importierten InterfaceflÃ¤chen ist nicht vollstÃ¤ndig geschlossen. "
                "Fehlende physische FlÃ¤chen oder Interfacezuordnungen korrigieren; keine "
                "automatische LÃ¼ckenreparatur.");
        for (const auto &loop : loops) {
            const auto key = loopKey(loop);
            if (const auto found = interfaces.find(key); found != interfaces.end()) {
                auto &interface = result.model.interfaces[found->second];
                if (interface.secondDomainId != exteriorId)
                    throw std::invalid_argument("Interface hat mehr als zwei Gebietsseiten.");
                interface.secondDomainId = regions[d].id;
                for (const auto t : interface.triangleIndices)
                    result.model.domains[d].panels.push_back({t, 1});
                continue;
            }
            if (lastFace == std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("Keine freie Interface-FlÃ¤chen-ID.");
            const auto index = result.model.interfaces.size();
            auto triangles = appendCap(result.mesh, loop, ++lastFace, options, stop);
            for (const auto t : triangles)
                result.model.domains[d].panels.push_back({t, -1});
            result.model.interfaces.push_back({index32(index + 1),
                                               regions[d].name + " / Ã–ffnung " +
                                                   std::to_string(loops.size() > 1 ? index + 1 : 1),
                                               regions[d].id, exteriorId, std::move(triangles)});
            interfaces.emplace(key, index);
        }
    }
    for (const auto &interface : result.model.interfaces)
        if (interface.secondDomainId == exteriorId)
            for (const auto t : interface.triangleIndices)
                result.model.domains[exteriorIndex].panels.push_back({t, 1});
    if (result.model.interfaces.empty())
        throw std::invalid_argument(
            "Kein Ã–ffnungsrand fÃ¼r ein transparentes Interface gefunden.");
    for (const auto &domain : result.model.domains) {
        cancelled(stop);
        if (baffle && domain.kind == BemDomainKind::exterior)
            continue; // Validated with the complete model below.
        if (symmetry)
            validateSymmetricBoundary(result.mesh, domain, options.quality, *symmetry);
        else
            validateBoundary(result.mesh, domain, options.quality);
    }
    validatePartition(result.mesh, result.model, symmetry, baffle);
    return result;
}
AcousticPartition partitionAcousticSurface(const SurfaceMesh &physical,
                                           const AcousticModelDefinition &definition,
                                           const MeshingOptions &options, std::stop_token stop,
                                           const OrthogonalSymmetryPlanes *symmetry) {
    if (!definition.interfaces.empty()) {
        cancelled(stop);
        AcousticPartition result{physical, {}};
        std::map<std::uint32_t, std::size_t> domainIndex;
        std::map<std::uint32_t, std::size_t> physicalOwners;
        std::map<std::uint32_t, std::vector<std::size_t>> twoSidedOwners;
        std::map<std::uint32_t, std::size_t> interfaceOwners;
        std::set<std::uint32_t> availableFaces, interfaceIds;
        for (const auto &triangle : physical.triangles)
            availableFaces.insert(triangle.faceId);
        std::optional<std::size_t> exteriorIndex;
        for (const auto &domain : definition.domains) {
            const auto index = result.model.domains.size();
            if (!domain.id || domain.name.empty() || domain.needsAssignmentReview ||
                !domainIndex.emplace(domain.id, index).second ||
                (domain.kind != BemDomainKind::interior && domain.kind != BemDomainKind::exterior))
                throw std::invalid_argument("Ungültige oder ungeprüfte Gebietsdefinition.");
            if (domain.kind == BemDomainKind::exterior) {
                if (exteriorIndex)
                    throw std::invalid_argument("Mehrere Außenräume sind nicht zulässig.");
                exteriorIndex = index;
            }
            result.model.domains.push_back({domain.id, domain.name, domain.kind, {}});
            for (auto face : domain.faceIds)
                if (!availableFaces.contains(face) || !physicalOwners.emplace(face, index).second)
                    throw std::invalid_argument(
                        "Gebietsfläche fehlt oder ist mehrfach zugeordnet.");
            std::set<std::uint32_t> unique;
            for (auto face : domain.twoSidedFaceIds) {
                if (!availableFaces.contains(face) || !unique.insert(face).second)
                    throw std::invalid_argument("Zweiseitige Fläche fehlt oder ist doppelt.");
                auto &owners = twoSidedOwners[face];
                owners.push_back(index);
                if (owners.size() > 2)
                    throw std::invalid_argument("Zweiseitige Fläche hat mehr als zwei Seiten.");
            }
        }
        if (!exteriorIndex || result.model.domains.size() < 2)
            throw std::invalid_argument("Innengebiet und genau ein Außenraum erforderlich.");
        for (const auto &link : definition.interfaces) {
            if (!link.id || link.name.empty() || link.needsAssignmentReview ||
                link.faceIds.empty() || !interfaceIds.insert(link.id).second ||
                link.firstDomainId == link.secondDomainId ||
                !domainIndex.contains(link.firstDomainId) ||
                !domainIndex.contains(link.secondDomainId))
                throw std::invalid_argument("Ungültige oder ungeprüfte Interfacezuordnung.");
            const auto index = result.model.interfaces.size();
            result.model.interfaces.push_back(
                {link.id, link.name, link.firstDomainId, link.secondDomainId, {}});
            for (auto face : link.faceIds)
                if (!availableFaces.contains(face) || !interfaceOwners.emplace(face, index).second)
                    throw std::invalid_argument(
                        "Interfacefläche fehlt oder ist mehrfach zugeordnet.");
        }
        for (auto [face, index] : interfaceOwners) {
            (void)index;
            if (physicalOwners.contains(face) || twoSidedOwners.contains(face))
                throw std::invalid_argument("Interfacefläche besitzt eine weitere Flächenrolle.");
        }
        for (const auto &[face, owners] : twoSidedOwners)
            if (physicalOwners.contains(face))
                throw std::invalid_argument(
                    "Zweiseitige Fläche besitzt eine weitere Flächenrolle.");
        std::map<std::pair<std::size_t, std::size_t>, std::size_t> surfacePairs;
        for (std::size_t i = 0; i < physical.triangles.size(); ++i) {
            const auto face = physical.triangles[i].faceId;
            const auto panel = index32(i);
            if (const auto found = interfaceOwners.find(face); found != interfaceOwners.end()) {
                auto &link = result.model.interfaces[found->second];
                link.triangleIndices.push_back(panel);
                result.model.domains[domainIndex.at(link.firstDomainId)].panels.push_back(
                    {panel, -1});
                result.model.domains[domainIndex.at(link.secondDomainId)].panels.push_back(
                    {panel, 1});
            } else if (const auto found = twoSidedOwners.find(face);
                       found != twoSidedOwners.end()) {
                const auto first = found->second.front();
                const auto second =
                    found->second.size() == 2 ? found->second.back() : *exteriorIndex;
                if (first == second)
                    throw std::invalid_argument("Zweiseitige Fläche benötigt zwei Gebiete.");
                const auto [where, added] = surfacePairs.emplace(
                    std::make_pair(first, second), result.model.twoSidedSurfaces.size());
                if (added)
                    result.model.twoSidedSurfaces.push_back(
                        {result.model.domains[first].id, result.model.domains[second].id, {}});
                result.model.twoSidedSurfaces[where->second].triangleIndices.push_back(panel);
                result.model.domains[first].panels.push_back({panel, -1});
                result.model.domains[second].panels.push_back({panel, 1});
            } else {
                const auto ownerIt = physicalOwners.find(face);
                const auto owner =
                    ownerIt == physicalOwners.end() ? *exteriorIndex : ownerIt->second;
                result.model.domains[owner].panels.push_back({panel, -1});
            }
        }
        const auto quality = inspectMesh(physical);
        if (options.quality.enabled &&
            (quality.minimumAngleDegrees + 1e-9 < options.quality.minimumAngleDegrees ||
             quality.maximumAspectRatio > options.quality.maximumAspectRatio * (1 + 1e-12)))
            throw std::invalid_argument("Gebietsnetz verletzt die eingestellten Qualitätsgrenzen.");
        validatePartition(result.mesh, result.model, symmetry,
                          definition.exteriorBaffle ? &*definition.exteriorBaffle : nullptr);
        return result;
    }
    std::vector<AcousticRegionDefinition> regions;
    const AcousticDomainDefinition *exterior = nullptr;
    std::map<std::uint32_t, std::size_t> interiorIndex;
    std::set<std::uint32_t> domainIds;
    for (const auto &domain : definition.domains) {
        if (!domain.id || domain.name.empty() || !domainIds.insert(domain.id).second ||
            domain.needsAssignmentReview)
            throw std::invalid_argument("Ungültige oder ungeprüfte Gebietsdefinition.");
        if (domain.kind == BemDomainKind::exterior) {
            if (exterior || !domain.twoSidedFaceIds.empty())
                throw std::invalid_argument(
                    "Genau ein Außenraum ohne eigene zweiseitige Flächen erforderlich.");
            exterior = &domain;
        } else if (domain.kind == BemDomainKind::interior) {
            interiorIndex.emplace(domain.id, regions.size());
            regions.push_back(
                {domain.id, domain.name, domain.faceIds, {}, false, domain.twoSidedFaceIds});
        } else {
            throw std::invalid_argument("Unbekannte Gebietsart.");
        }
    }
    if (!exterior || regions.empty())
        throw std::invalid_argument("Innengebiet und genau ein Außenraum erforderlich.");
    std::set<std::uint32_t> ids, faces;
    for (const auto &link : definition.interfaces) {
        if (!link.id || link.name.empty() || !ids.insert(link.id).second ||
            link.needsAssignmentReview || link.faceIds.empty() ||
            !interiorIndex.contains(link.firstDomainId) ||
            link.firstDomainId == link.secondDomainId || !domainIds.contains(link.secondDomainId))
            throw std::invalid_argument("Ungültige oder ungeprüfte Interfacezuordnung.");
        const auto first = interiorIndex.at(link.firstDomainId);
        if (link.secondDomainId != exterior->id && (!interiorIndex.contains(link.secondDomainId) ||
                                                    interiorIndex.at(link.secondDomainId) <= first))
            throw std::invalid_argument(
                "Die erste Interfaceseite muss vor der zweiten Innengebietsseite stehen.");
        for (auto face : link.faceIds) {
            if (!faces.insert(face).second)
                throw std::invalid_argument("Interfacefläche ist mehrfach zugeordnet.");
            regions[first].interfaceFaceIds.push_back(face);
            if (link.secondDomainId != exterior->id)
                regions[interiorIndex.at(link.secondDomainId)].interfaceFaceIds.push_back(face);
        }
    }
    auto result =
        partitionAcousticSurface(physical, regions, options, stop, symmetry,
                                 definition.exteriorBaffle ? &*definition.exteriorBaffle : nullptr);
    auto &actualExterior = result.model.domains.back();
    for (auto face : exterior->faceIds) {
        bool found = false;
        for (const auto &panel : actualExterior.panels)
            if (result.mesh.triangles[panel.triangleIndex].faceId == face && !faces.contains(face))
                found = true;
        if (!found)
            throw std::invalid_argument("Außenraumfläche fehlt oder ist anderweitig zugeordnet.");
    }
    const auto oldExteriorId = actualExterior.id;
    actualExterior.id = exterior->id;
    actualExterior.name = exterior->name;
    for (auto &surface : result.model.twoSidedSurfaces) {
        if (surface.firstDomainId == oldExteriorId)
            surface.firstDomainId = exterior->id;
        if (surface.secondDomainId == oldExteriorId)
            surface.secondDomainId = exterior->id;
    }
    if (definition.interfaces.empty()) {
        for (auto &link : result.model.interfaces) {
            if (link.firstDomainId == oldExteriorId)
                link.firstDomainId = exterior->id;
            if (link.secondDomainId == oldExteriorId)
                link.secondDomainId = exterior->id;
        }
    } else {
        result.model.interfaces.clear();
        for (const auto &link : definition.interfaces) {
            BemInterface out{link.id, link.name, link.firstDomainId, link.secondDomainId, {}};
            for (std::size_t i = 0; i < physical.triangles.size(); ++i)
                if (std::find(link.faceIds.begin(), link.faceIds.end(),
                              physical.triangles[i].faceId) != link.faceIds.end())
                    out.triangleIndices.push_back(index32(i));
            result.model.interfaces.push_back(std::move(out));
        }
    }
    validatePartition(result.mesh, result.model, symmetry,
                      definition.exteriorBaffle ? &*definition.exteriorBaffle : nullptr);
    return result;
}
} // namespace hornsculpt
