#include <hornsculpt/axisymmetric_mesh.hpp>
#include <hornsculpt/meshing.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <stdexcept>

namespace hornsculpt {
namespace {
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
bool finite(Vec3 v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Rotationsnetz abgebrochen.");
}
} // namespace

AxisymmetricMesh buildAxisymmetricMesh(const AxisymmetricMeridian &m, unsigned sectors,
                                       std::stop_token stop) {
    cancelled(stop);
    if (sectors < 8 || sectors % 4 != 0 || sectors > 4096)
        throw std::invalid_argument(
            "Rotationsnetz: Sektorzahl muss 8..4096 und durch 4 teilbar sein.");
    if (!finite(m.origin) || !finite(m.axis) || !finite(m.radial) ||
        std::abs(dot(m.axis, m.axis) - 1) > 1e-12 ||
        std::abs(dot(m.radial, m.radial) - 1) > 1e-12 || std::abs(dot(m.axis, m.radial)) > 1e-12 ||
        m.points.empty() || m.segments.empty())
        throw std::invalid_argument("Rotationsnetz: Meridian und orthonormale Achsen prüfen.");
    if (m.points.size() > std::numeric_limits<std::uint32_t>::max() / sectors ||
        m.segments.size() > std::numeric_limits<std::uint32_t>::max() / (2 * sectors))
        throw std::invalid_argument("Rotationsnetz ist zu groß.");
    const auto transverse = cross(m.axis, m.radial);
    AxisymmetricMesh out;
    std::vector<std::vector<std::uint32_t>> rings(m.points.size());
    std::vector<unsigned> incoming(m.points.size()), outgoing(m.points.size());
    for (std::size_t p = 0; p < m.points.size(); ++p) {
        cancelled(stop);
        const auto point = m.points[p];
        if (!std::isfinite(point.radius) || !std::isfinite(point.axial) || point.radius < 0)
            throw std::invalid_argument("Rotationsnetz: ungültiger Meridianpunkt.");
        rings[p].resize(sectors);
        const unsigned copies = point.radius == 0 ? 1 : sectors;
        for (unsigned s = 0; s < copies; ++s) {
            const double phi = 2 * std::numbers::pi * s / sectors;
            const double rcos = point.radius * std::cos(phi), rsin = point.radius * std::sin(phi);
            const auto index = static_cast<std::uint32_t>(out.mesh.vertices.size());
            rings[p][s] = index;
            out.mesh.vertices.push_back(
                {m.origin.x + point.axial * m.axis.x + rcos * m.radial.x + rsin * transverse.x,
                 m.origin.y + point.axial * m.axis.y + rcos * m.radial.y + rsin * transverse.y,
                 m.origin.z + point.axial * m.axis.z + rcos * m.radial.z + rsin * transverse.z});
            out.pressureDofs.push_back(static_cast<std::uint32_t>(p));
        }
        if (copies == 1)
            std::fill(rings[p].begin(), rings[p].end(), rings[p][0]);
    }
    std::map<std::uint32_t, BoundaryCondition> conditions;
    for (std::size_t e = 0; e < m.segments.size(); ++e) {
        cancelled(stop);
        const auto edge = m.segments[e];
        if (edge.first >= m.points.size() || edge.second >= m.points.size() ||
            edge.first == edge.second || !edge.faceId ||
            (edge.boundary != BoundaryCondition::wall &&
             edge.boundary != BoundaryCondition::driven))
            throw std::invalid_argument("Rotationsnetz: ungültiger Meridianabschnitt.");
        const auto a = m.points[edge.first], b = m.points[edge.second];
        if ((a.radius == b.radius && a.axial == b.axial) || (a.radius == 0 && b.radius == 0))
            throw std::invalid_argument("Rotationsnetz: entarteter Abschnitt auf der Achse.");
        if (++outgoing[edge.first] > 1 || ++incoming[edge.second] > 1)
            throw std::invalid_argument(
                "Rotationsnetz: verzweigter oder uneinheitlich orientierter Meridian.");
        const auto [it, inserted] = conditions.emplace(edge.faceId, edge.boundary);
        if (!inserted && it->second != edge.boundary)
            throw std::invalid_argument("Rotationsnetz: inkonsistente Flächenrandbedingung.");
        for (unsigned s = 0; s < sectors; ++s) {
            const auto n = (s + 1) % sectors;
            const auto append = [&](std::array<std::uint32_t, 3> tri) {
                if (s == 0)
                    out.testPanels.push_back(static_cast<std::uint32_t>(out.mesh.triangles.size()));
                out.mesh.triangles.push_back({tri, edge.faceId, edge.boundary});
                out.fluxDofs.push_back(static_cast<std::uint32_t>(e));
            };
            if (a.radius == 0)
                append({rings[edge.first][0], rings[edge.second][n], rings[edge.second][s]});
            else if (b.radius == 0)
                append({rings[edge.first][s], rings[edge.first][n], rings[edge.second][0]});
            else {
                append({rings[edge.first][s], rings[edge.first][n], rings[edge.second][s]});
                append({rings[edge.first][n], rings[edge.second][n], rings[edge.second][s]});
            }
        }
    }
    for (std::size_t p = 0; p < m.points.size(); ++p) {
        if (m.points[p].radius == 0) {
            if (incoming[p] + outgoing[p] != 1)
                throw std::invalid_argument("Rotationsnetz: Achspunkt ist kein Meridianende.");
        } else if (incoming[p] != 1 || outgoing[p] != 1)
            throw std::invalid_argument("Rotationsnetz: Meridian besitzt eine offene Kante.");
    }
    const auto quality = inspectMesh(out.mesh);
    if (!quality.watertight() || quality.connectedComponents != 1 || !(quality.signedVolume > 0))
        throw std::invalid_argument(
            "Rotationsnetz: geschlossene, schnittfreie Außenorientierung nicht erfüllt.");
    return out;
}
} // namespace hornsculpt
