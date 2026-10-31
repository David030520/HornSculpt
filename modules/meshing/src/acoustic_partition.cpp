#include "frontal_delaunay.hpp"
#include <hornsculpt/acoustic_partition.hpp>

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
            throw std::invalid_argument("Mehrdeutiger Öffnungsrand oder Dreifachanschluss.");
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
                    "Teilgebietsöffnung besitzt keinen geschlossenen Rand.");
            loop.push_back(vertex);
            vertex = found->second;
            next.erase(found);
        } while (vertex != first);
        if (loop.size() < 3)
            throw std::invalid_argument("Teilgebietsöffnung besitzt weniger als drei Randpunkte.");
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
            "': Randabschluss, Orientierung oder Schnittfreiheit ungültig.");
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
                    throw std::invalid_argument("Gebietsrandkomponenten berühren sich.");
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
                                    "': getrennte Luftvolumen benötigen eigene Gebiete.");
    if (limits.enabled && (quality.minimumAngleDegrees + 1e-9 < limits.minimumAngleDegrees ||
                           quality.maximumAspectRatio > limits.maximumAspectRatio * (1 + 1e-12)))
        throw std::invalid_argument("Gebiet '" + domain.name +
                                    "': Netz unterschreitet die eingestellten Qualitätsgrenzen.");
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
        throw std::invalid_argument("Interfaceöffnung hat keine eindeutige ebene Fläche.");
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
            throw std::invalid_argument("Interfaceöffnung ist nicht eben. CAD-Schnittfläche oder "
                                        "einen ebenen Öffnungsrand wählen.");
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
        throw std::invalid_argument("Interface benötigt feineren gemeinsamen Öffnungsrand. "
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

void validateAcousticPartition(const SurfaceMesh &mesh, const BemSubdomainModel &model) {
    if (model.domains.empty())
        throw std::invalid_argument("Keine akustischen Gebiete.");
    // The combined complex has intentional three-face edges, but geometric
    // crossings and duplicate triangles remain invalid independently of that.
    const auto combined = inspectMesh(mesh);
    if (combined.degenerateTriangles || combined.selfIntersectionPairs)
        throw std::invalid_argument("Gebietsnetz enthält entartete oder sich schneidende Flächen.");
    std::map<std::uint32_t, const BemSubdomain *> domains;
    std::vector<std::map<std::uint32_t, int>> owners(mesh.triangles.size());
    unsigned exteriorCount{};
    for (const auto &domain : model.domains) {
        if (!domain.id || domain.name.empty() || domain.panels.empty() ||
            !domains.emplace(domain.id, &domain).second ||
            (domain.kind != BemDomainKind::interior && domain.kind != BemDomainKind::exterior))
            throw std::invalid_argument("Ungültige oder doppelte akustische Gebietsidentität.");
        if (domain.kind == BemDomainKind::exterior)
            ++exteriorCount;
        for (const auto &panel : domain.panels) {
            if (panel.triangleIndex >= mesh.triangles.size() ||
                (panel.normalSign != 1 && panel.normalSign != -1) ||
                !owners[panel.triangleIndex].emplace(domain.id, panel.normalSign).second)
                throw std::invalid_argument(
                    "Gebietsseite besitzt ungültige Indizes, Normalen oder doppelte Panels.");
        }
    }
    if (exteriorCount > 1)
        throw std::invalid_argument(
            "Mehrere unbeschränkte Außenräume sind im aktuellen Gebietsvertrag nicht zulässig.");
    std::set<std::uint32_t> interfaceIds, interfacePanels, interfaceVertices;
    for (const auto &interface : model.interfaces) {
        if (!interface.id || interface.name.empty() || !interfaceIds.insert(interface.id).second ||
            interface.firstDomainId == interface.secondDomainId ||
            !domains.contains(interface.firstDomainId) ||
            !domains.contains(interface.secondDomainId) || interface.triangleIndices.empty())
            throw std::invalid_argument("Ungültige Interfaceidentität oder Gebietsreferenz.");
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
                throw std::invalid_argument("Interface benötigt genau zwei gegenorientierte Seiten "
                                            "ohne zusätzliche Anregung.");
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
                throw std::invalid_argument("Interfaces mit gemeinsamen Randknoten benötigen einen "
                                            "Mehrfachanschlussvertrag.");
    }
    for (std::size_t i = 0; i < owners.size(); ++i)
        if (owners[i].empty() ||
            (interfacePanels.contains(static_cast<std::uint32_t>(i)) ? owners[i].size() != 2
                                                                     : owners[i].size() != 1))
            throw std::invalid_argument(
                "Netzdreieck fehlt in der Gebietszuordnung oder hat mehrere physische Eigentümer.");
    const MeshQualityLimits topologyOnly{false, 20, 4};
    for (const auto &domain : model.domains)
        validateBoundary(mesh, domain, topologyOnly);
}

AcousticPartition partitionAcousticSurface(const SurfaceMesh &physical,
                                           std::span<const AcousticRegionDefinition> regions,
                                           const MeshingOptions &options, std::stop_token stop) {
    cancelled(stop);
    if (regions.empty())
        throw std::invalid_argument("Keine akustischen Innengebiete definiert.");
    if (!std::isfinite(options.targetEdgeLength) || options.targetEdgeLength <= 0)
        throw std::invalid_argument("Ungültige Interface-Zielkante.");
    const bool explicitInterfaces =
        std::any_of(regions.begin(), regions.end(),
                    [](const auto &region) { return !region.interfaceFaceIds.empty(); });
    const auto sourceQuality = inspectMesh(physical);
    if ((!explicitInterfaces && !sourceQuality.watertight()) || sourceQuality.degenerateTriangles ||
        sourceQuality.selfIntersectionPairs)
        throw std::invalid_argument("Physisches Quellnetz ist nicht geschlossen und schnittfrei; "
                                    "Interfaces reparieren keine Modelllücken.");
    AcousticPartition result{physical, {}};
    std::set<std::uint32_t> domainIds, availableFaces;
    std::map<std::uint32_t, std::size_t> assignments;
    std::map<std::uint32_t, std::vector<std::size_t>> interfaceAssignments;
    std::uint32_t lastFace{};
    for (const auto &triangle : physical.triangles) {
        availableFaces.insert(triangle.faceId);
        lastFace = std::max(lastFace, triangle.faceId);
    }
    for (const auto &region : regions) {
        if (region.needsAssignmentReview)
            throw std::invalid_argument(
                "Gebiet '" + region.name +
                "': Flächenzuordnung nach Quellen- oder Topologieänderung prüfen.");
        if (!region.id || region.name.empty() ||
            (region.faceIds.empty() && region.interfaceFaceIds.empty()) ||
            !domainIds.insert(region.id).second)
            throw std::invalid_argument(
                "Gebiete benötigen eindeutige positive IDs, Namen und physische Flächen.");
        const auto index = result.model.domains.size();
        result.model.domains.push_back({region.id, region.name, BemDomainKind::interior, {}});
        for (const auto face : region.faceIds)
            if (!availableFaces.contains(face) || !assignments.emplace(face, index).second)
                throw std::invalid_argument("Gebietsfläche fehlt oder ist mehrfach zugeordnet.");
        std::set<std::uint32_t> uniqueInterfaces;
        for (const auto face : region.interfaceFaceIds) {
            if (!availableFaces.contains(face) || !uniqueInterfaces.insert(face).second)
                throw std::invalid_argument(
                    "Interfacefläche fehlt oder ist im Gebiet mehrfach zugeordnet.");
            auto &owners = interfaceAssignments[face];
            owners.push_back(index);
            if (owners.size() > 2)
                throw std::invalid_argument("Interfacefläche besitzt mehr als zwei Gebietsseiten.");
        }
    }
    for (const auto &[face, owners] : interfaceAssignments) {
        (void)owners;
        if (assignments.contains(face))
            throw std::invalid_argument(
                "Eine Interfacefläche darf nicht zugleich Wall oder Driven eines Gebiets sein.");
    }
    std::uint32_t exteriorId = 1;
    while (domainIds.contains(exteriorId)) {
        if (exteriorId == std::numeric_limits<std::uint32_t>::max())
            throw std::invalid_argument("Keine freie Gebiets-ID.");
        ++exteriorId;
    }
    const auto exteriorIndex = result.model.domains.size();
    result.model.domains.push_back({exteriorId, "Außenraum", BemDomainKind::exterior, {}});
    for (std::size_t i = 0; i < physical.triangles.size(); ++i) {
        if (interfaceAssignments.contains(physical.triangles[i].faceId))
            continue;
        const auto found = assignments.find(physical.triangles[i].faceId);
        const auto domain = found == assignments.end() ? exteriorIndex : found->second;
        result.model.domains[domain].panels.push_back({index32(i), -1});
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
        const auto loops = openingLoops(physical, result.model.domains[d]);
        if (explicitInterfaces && !loops.empty())
            throw std::invalid_argument(
                "Gebiet mit importierten Interfaceflächen ist nicht vollständig geschlossen. "
                "Fehlende physische Flächen oder Interfacezuordnungen korrigieren; keine "
                "automatische Lückenreparatur.");
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
                throw std::invalid_argument("Keine freie Interface-Flächen-ID.");
            const auto index = result.model.interfaces.size();
            auto triangles = appendCap(result.mesh, loop, ++lastFace, options, stop);
            for (const auto t : triangles)
                result.model.domains[d].panels.push_back({t, -1});
            result.model.interfaces.push_back(
                {index32(index + 1),
                 regions[d].name + " / Öffnung " + std::to_string(loops.size() > 1 ? index + 1 : 1),
                 regions[d].id, exteriorId, std::move(triangles)});
            interfaces.emplace(key, index);
        }
    }
    for (const auto &interface : result.model.interfaces)
        if (interface.secondDomainId == exteriorId)
            for (const auto t : interface.triangleIndices)
                result.model.domains[exteriorIndex].panels.push_back({t, 1});
    if (result.model.interfaces.empty())
        throw std::invalid_argument("Kein Öffnungsrand für ein transparentes Interface gefunden.");
    for (const auto &domain : result.model.domains) {
        cancelled(stop);
        validateBoundary(result.mesh, domain, options.quality);
    }
    validateAcousticPartition(result.mesh, result.model);
    return result;
}
} // namespace hornsculpt
