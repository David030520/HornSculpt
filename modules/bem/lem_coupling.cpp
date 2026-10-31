#include "lem_coupling.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace hornsculpt::bem_detail {
CouplingBasis makeCouplingBasis(const lem::Network &network, const SurfaceMesh &mesh,
                                std::span<const double> panelAreas, CouplingDomain domain) {
    lem::validateNetwork(network);
    const auto ports = lem::bemPorts(network);
    if (ports.empty() || ports.size() > 1024)
        throw std::invalid_argument("Coupling needs between one and 1024 BEM ports.");
    if (panelAreas.size() != mesh.triangles.size())
        throw std::logic_error("BEM coupling area count does not match the mesh.");
    struct Face {
        double area{};
        BoundaryCondition boundary{};
    };
    std::map<std::uint32_t, Face> faces;
    for (std::size_t i = 0; i < mesh.triangles.size(); ++i) {
        const auto &triangle = mesh.triangles[i];
        const auto [it, inserted] = faces.emplace(triangle.faceId, Face{0, triangle.boundary});
        if (!inserted && it->second.boundary != triangle.boundary)
            throw std::invalid_argument("Conflicting boundary conditions on one face.");
        it->second.area += panelAreas[i];
    }
    CouplingBasis basis;
    std::map<std::uint32_t, int> facePorts;
    std::map<std::uint32_t, double> faceVelocities;
    for (std::size_t j = 0; j < ports.size(); ++j) {
        const auto &port = std::get<lem::BemPort>(ports[j]->component);
        if (port.kind == lem::BemPortKind::interfacePort && domain != CouplingDomain::baffledAir)
            throw std::invalid_argument(
                "BEM interface ports require the baffled air-cavity domain.");
        const auto boundary = port.kind == lem::BemPortKind::interfacePort
                                  ? BoundaryCondition::wall
                                  : BoundaryCondition::driven;
        double area = 0;
        for (const auto id : port.faceIds) {
            const auto found = faces.find(id);
            if (found == faces.end() || found->second.boundary != boundary ||
                !facePorts.emplace(id, static_cast<int>(j)).second)
                throw std::invalid_argument(
                    "Coupled faces must exist, match their port type (Driven surface or Wall "
                    "interface), and belong to exactly one port.");
            const double weight =
                port.velocityWeights.contains(id) ? port.velocityWeights.at(id) : 1;
            area += weight * found->second.area;
        }
        if (!(area > 0) || !std::isfinite(area))
            throw std::invalid_argument("Invalid BEM coupling port area.");
        basis.areas.push_back(area);
        BemExcitation excitation{ports[j]->name, {}};
        for (const auto id : port.faceIds) {
            const double weight =
                port.velocityWeights.contains(id) ? port.velocityWeights.at(id) : 1;
            const double velocity = weight / area;
            if (!(velocity > 0) || !std::isfinite(velocity) ||
                velocity > std::numeric_limits<float>::max())
                throw std::invalid_argument("BEM port modal velocity exceeds storage.");
            excitation.faces.push_back({id, {velocity, 0}});
            faceVelocities.emplace(id, velocity);
        }
        basis.excitations.push_back(std::move(excitation));
    }
    for (const auto &[id, face] : faces)
        if (domain == CouplingDomain::exteriorSolid && face.boundary == BoundaryCondition::driven &&
            !facePorts.contains(id))
            throw std::invalid_argument("Every Driven face needs a BEM port in a coupled request.");
    for (const auto &triangle : mesh.triangles) {
        const auto found = facePorts.find(triangle.faceId);
        basis.panelPorts.push_back(found == facePorts.end() ? -1 : found->second);
        basis.panelVelocities.push_back(
            found == facePorts.end() ? 0 : faceVelocities.at(triangle.faceId));
    }
    return basis;
}

std::vector<FaceVelocity> coupledVelocities(const CouplingBasis &basis,
                                            const lem::Solution &solution) {
    if (solution.ports.size() != basis.excitations.size())
        throw std::logic_error("LEM result port count does not match its BEM basis.");
    std::vector<FaceVelocity> velocities;
    for (std::size_t j = 0; j < solution.ports.size(); ++j)
        for (const auto &face : basis.excitations[j].faces)
            velocities.push_back(
                {face.faceId, face.normalVelocity * solution.ports[j].volumeVelocity});
    return velocities;
}

void validateCouplingClosure(const CouplingBasis &basis, const lem::Solution &solution,
                             std::span<const double> panelAreas,
                             std::span<const std::complex<double>> surfacePressure) {
    if (panelAreas.size() != basis.panelPorts.size() ||
        panelAreas.size() != basis.panelVelocities.size() ||
        surfacePressure.size() != basis.panelPorts.size() ||
        solution.ports.size() != basis.areas.size())
        throw std::logic_error("BEM pressure-closure dimensions are inconsistent.");
    std::vector<std::complex<double>> means(basis.areas.size());
    std::vector<double> magnitudes(basis.areas.size());
    for (std::size_t panel = 0; panel < panelAreas.size(); ++panel) {
        const int port = basis.panelPorts[panel];
        if (port < 0)
            continue;
        const double weight = panelAreas[panel] * basis.panelVelocities[panel];
        means[port] += weight * surfacePressure[panel];
        magnitudes[port] += weight * std::abs(surfacePressure[panel]);
    }
    for (std::size_t port = 0; port < means.size(); ++port) {
        const auto pressure = solution.ports[port].pressure;
        const double scale = std::max(magnitudes[port], std::abs(pressure));
        const double difference = std::abs(means[port] - pressure);
        if (!std::isfinite(scale) || !std::isfinite(difference) ||
            difference > relativePortPressureClosureTolerance * scale)
            throw std::runtime_error(
                "BEM and LEM port pressures fail the feedback consistency check.");
    }
}
} // namespace hornsculpt::bem_detail
