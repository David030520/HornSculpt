#include "simulation_request.hpp"
#include "baffle_exterior.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>

namespace hornsculpt {
namespace {
// Dimensionless round-off budget for persisted normalized direction components.
constexpr double unitNormalTolerance = 1e-12;
Vec3 scaled(Vec3 p, double s) {
    return {p.x * s, p.y * s, p.z * s};
}
Vec3 added(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
} // namespace
long double estimateResultBytes(const ResultMemoryDimensions &d) {
    const long double boundary =
        static_cast<long double>(d.excitations) *
        (sizeof(BemBoundaryData) +
         (2.L * d.triangles + d.vertices + d.exteriorTriangles) * sizeof(std::complex<double>));
    const long double snapshot =
        sizeof(BemFrequencyResult) +
        static_cast<long double>(d.excitations) *
            (sizeof(BemExcitationResult) +
             static_cast<long double>(d.observers) * sizeof(std::complex<double>) +
             static_cast<long double>(d.faces) * (sizeof(FaceImpedance) + sizeof(FaceVelocity)));
    // Domain traces are shared immutably across realtime publications. Every domain
    // may conservatively reference all vertices/panels and interface vertices.
    const long double domainTraces = static_cast<long double>(d.domains) * d.excitations *
                                     (2.L * d.vertices + d.triangles) *
                                     (sizeof(std::complex<double>) + sizeof(std::uint32_t));
    return static_cast<long double>(d.frequencies) *
           (boundary + domainTraces + static_cast<long double>(d.snapshots) * snapshot);
}
OrthogonalSymmetryPlanes inferXYSymmetry(const SurfaceMesh &mesh, double toleranceMetres) {
    if (mesh.vertices.empty() || !std::isfinite(toleranceMetres) || toleranceMetres < 1e-9 ||
        toleranceMetres > .01)
        throw std::invalid_argument(
            "Symmetrie benötigt ein gültiges Netz und eine gültige Toleranz.");
    double minX = std::numeric_limits<double>::infinity(), maxX = -minX;
    double minY = minX, maxY = -minX;
    for (const auto &v : mesh.vertices) {
        if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z))
            throw std::invalid_argument("Symmetrienetz enthält nicht endliche Koordinaten.");
        minX = std::min(minX, v.x);
        maxX = std::max(maxX, v.x);
        minY = std::min(minY, v.y);
        maxY = std::max(maxY, v.y);
    }
    const auto side = [toleranceMetres](double low, double high) {
        const bool positive = high > toleranceMetres, negative = low < -toleranceMetres;
        if (positive == negative)
            throw std::invalid_argument(
                "Das Netz muss je Symmetrieebene eindeutig auf genau einer Seite liegen.");
        return positive ? 1. : -1.;
    };
    return {{0, 0, 0}, {side(minX, maxX), 0, 0}, {0, side(minY, maxY), 0}, toleranceMetres};
}
bool hasBoundInterfacePorts(const lem::Network &network) {
    return std::any_of(network.elements.begin(), network.elements.end(), [](const auto &element) {
        const auto *port = std::get_if<lem::BemPort>(&element.component);
        return port && port->kind == lem::BemPortKind::interfacePort &&
               (port->needsAssignmentReview || !port->faceIds.empty());
    });
}
void validateAxisymmetricSamplingSettings(const AxisymmetricSamplingSettings &s) {
    if (!std::isfinite(s.maximumLength) || s.maximumLength <= 0 ||
        !std::isfinite(s.chordTolerance) || s.chordTolerance <= 0 ||
        !std::isfinite(s.rotationalTolerance) || s.rotationalTolerance <= 0 || s.sectors < 8 ||
        s.sectors > 4096 || s.sectors % 4)
        throw std::invalid_argument(
            "2D-Rotation: Meridianlänge, Sehnenabweichung und Symmetrietoleranz müssen positiv "
            "sein; Sektoren 8 bis 4096, durch 4 teilbar.");
}
void validateInfiniteBaffleSettings(const InfiniteBaffle &baffle) {
    const auto finite = [](Vec3 value) {
        return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
    };
    if (!finite(baffle.origin) || !finite(baffle.normal) ||
        std::abs(std::hypot(baffle.normal.x, baffle.normal.y, baffle.normal.z) - 1) >
            unitNormalTolerance)
        throw std::invalid_argument(
            "Schallwand: Ursprung muss endlich und die Normale ein Einheitsvektor sein.");
    std::set<std::uint32_t> faces;
    for (const auto id : baffle.apertureFaceIds)
        if (!id || !faces.insert(id).second)
            throw std::invalid_argument(
                "Schallwand: Apertur-Flächen müssen positive, eindeutige IDs besitzen.");
    if (faces.empty())
        throw std::invalid_argument("Schallwand: Mindestens eine Apertur-Fläche ist erforderlich.");
    bem_detail::validateBaffleExteriorMesh({}, baffle.exteriorDrivenMesh, baffle.origin,
                                           baffle.normal);
}

AcousticFrame simulationFrame(const SimulationSettings &settings, Vec3 sourceOrigin,
                              Vec3 sourceAxis) {
    return settings.observationFrame ? normalizedAcousticFrame(*settings.observationFrame)
                                     : defaultAcousticFrame(sourceOrigin, sourceAxis);
}

std::vector<std::uint8_t> exteriorObserverValidity(const BemRequest &request) {
    if (!request.exteriorBaffle)
        return {};
    if (!request.farField)
        throw std::invalid_argument("Schallwand benötigt Fernfeldbeobachter.");
    const auto n = request.exteriorBaffle->normal;
    const double nLength = std::hypot(n.x, n.y, n.z);
    if (!std::isfinite(nLength) || std::abs(nLength - 1.) > 1e-8)
        throw std::invalid_argument("Ungültige Schallwandnormale.");
    std::vector<std::uint8_t> mask;
    mask.reserve(request.observers.size());
    for (const auto &point : request.observers) {
        const auto direction = added(point, scaled(request.farField->origin, -1));
        const double length = std::hypot(direction.x, direction.y, direction.z);
        if (!std::isfinite(length) || length <= 0)
            throw std::invalid_argument("Ungültige Fernfeldrichtung.");
        const double projection =
            (direction.x * n.x + direction.y * n.y + direction.z * n.z) / length;
        mask.push_back(projection >= -64 * std::numeric_limits<double>::epsilon() ? 1 : 0);
    }
    return mask;
}
BemRequest makeSimulationRequest(const SimulationSettings &s,
                                 std::span<const FaceVelocity> velocities, Vec3 origin, Vec3 axis) {
    validateAxisymmetricSamplingSettings(s.axisymmetricSampling);
    if (!std::isfinite(s.symmetryPlaneTolerance) || s.symmetryPlaneTolerance < 1e-9 ||
        s.symmetryPlaneTolerance > .01)
        throw std::invalid_argument("Ungültige Symmetrietoleranz.");
    if (s.symmetryXYEnabled &&
        (s.acousticModel.domains.empty() || s.geometryMode != BemGeometryMode::threeDimensional ||
         s.acceleration != BemAcceleration::dense || !s.farField || s.couplingEnabled))
        throw std::invalid_argument("XY-Symmetrie benötigt Teilgebiete, 3D, dichten Solver und "
                                    "Fernfeld ohne LEM-Kopplung.");
    if (s.acousticModel.exteriorBaffle && !s.symmetryXYEnabled)
        throw std::invalid_argument("Die Außenraum-Schallwand benötigt XY-Symmetrie.");
    if ((std::any_of(s.acousticModel.domains.begin(), s.acousticModel.domains.end(),
                     [](const auto &r) { return r.needsAssignmentReview; }) ||
         std::any_of(s.acousticModel.interfaces.begin(), s.acousticModel.interfaces.end(),
                     [](const auto &i) { return i.needsAssignmentReview; })))
        throw std::invalid_argument(
            "Teilgebietsflächen müssen nach dem Geometriewechsel geprüft werden.");
    if (!s.acousticModel.domains.empty() &&
        (s.infiniteBaffle || s.geometryMode != BemGeometryMode::threeDimensional ||
         s.acceleration != BemAcceleration::dense))
        throw std::invalid_argument("Teilgebiete unterstützen derzeit nur 3D mit dichtem Solver "
                                    "ohne unendliche Schallwand.");
    if (s.geometryMode != BemGeometryMode::threeDimensional &&
        s.geometryMode != BemGeometryMode::axisymmetric)
        throw std::invalid_argument("Unbekannter BEM-Geometriemodus.");
    if (s.geometryMode == BemGeometryMode::axisymmetric && s.acceleration != BemAcceleration::dense)
        throw std::invalid_argument("2D-Rotation (m=0) unterstützt den dichten reduzierten Solver; "
                                    "Fast Multipole ausdrücklich ausschalten.");
    if (s.geometryMode == BemGeometryMode::axisymmetric && !s.infiniteBaffle)
        throw std::invalid_argument(
            "2D-Rotation (m=0) benötigt eine Luftkavität mit ebener Apertur in einer "
            "unendlichen Schallwand. Schallwand und Aperturflächen ausdrücklich zuordnen.");
    if (!s.couplingEnabled && hasBoundInterfacePorts(s.couplingNetwork))
        throw std::invalid_argument(
            "Virtuelle BEM↔LEM-Interfaces sind zugeordnet. Kopplung aktivieren; Ausschalten darf "
            "die Rechenöffnung nicht als Materialwand schließen.");
    if (!std::isfinite(s.minimumFrequency) || !std::isfinite(s.maximumFrequency) ||
        s.minimumFrequency <= 0 || s.maximumFrequency < s.minimumFrequency || !s.frequencyCount ||
        (s.frequencyCount > 1 && s.maximumFrequency == s.minimumFrequency) ||
        !std::isfinite(s.medium.density) || s.medium.density <= 0 ||
        !std::isfinite(s.medium.soundSpeed) || s.medium.soundSpeed <= 0 ||
        !std::isfinite(s.velocityScale) || s.velocityScale < 0 || !std::isfinite(s.phaseDegrees) ||
        !std::isfinite(s.observationDistance) || s.observationDistance <= 0 ||
        (!s.couplingEnabled && velocities.empty()))
        throw std::invalid_argument("Ungültige Frequenz-, Medium- oder Antriebseinstellungen.");
    if (s.nearMaximumDepth > 10)
        throw std::invalid_argument(
            "Die maximale Nahquadraturtiefe muss zwischen 0 und 10 liegen.");
    if (s.infiniteBaffle) {
        validateInfiniteBaffleSettings(*s.infiniteBaffle);
        if (s.acceleration != BemAcceleration::dense)
            throw std::invalid_argument(
                "Unendliche Schallwand unterstützt ausschließlich den dichten GPU-Solver.");
        if (!s.farField)
            throw std::invalid_argument(
                "H/V-Auswertung der unendlichen Schallwand unterstützt nur Fernfeld.");
    }
    const auto frame = simulationFrame(s, origin, axis);
    origin = frame.origin;
    axis = frame.forward;
    BemRequest request;
    request.infiniteBaffle = s.infiniteBaffle;
    request.exteriorBaffle = s.acousticModel.exteriorBaffle;
    if (s.farField)
        request.farField = FarFieldObservation{origin, s.observationDistance};
    request.medium = s.medium;
    request.numerics.acceleration = s.acceleration;
    request.numerics.fastMultipole = s.fastMultipole;
    request.numerics.nearMaximumDepth = s.nearMaximumDepth;
    if (s.singularOrder < 2 || s.singularOrder > 20)
        throw std::invalid_argument("Singular quadrature order must be from 2 to 20.");
    request.numerics.singularOrder = s.singularOrder;
    request.frequencies.reserve(s.frequencyCount);
    for (std::size_t i = 0; i < s.frequencyCount; ++i)
        request.frequencies.push_back(
            s.frequencyCount == 1
                ? s.minimumFrequency
                : s.minimumFrequency *
                      std::pow(s.maximumFrequency / s.minimumFrequency,
                               static_cast<double>(i) / static_cast<double>(s.frequencyCount - 1)));
    request.interactionPlanFrequencies = request.frequencies;
    if (s.couplingEnabled) {
        lem::validateNetwork(s.couplingNetwork);
        const auto ports = lem::bemPorts(s.couplingNetwork);
        if (ports.empty())
            throw std::invalid_argument(
                "BEM↔LEM benötigt mindestens einen zugeordneten BEM-Anschluss.");
        const bool interfaces = std::any_of(ports.begin(), ports.end(), [](const auto *element) {
            return std::get<lem::BemPort>(element->component).kind ==
                   lem::BemPortKind::interfacePort;
        });
        if (interfaces && !s.infiniteBaffle)
            throw std::invalid_argument(
                "Virtuelle BEM↔LEM-Interfaces benötigen die unendliche Schallwand.");
        request.coupling = s.couplingNetwork;
    }
    if (!s.couplingEnabled || (s.infiniteBaffle && !velocities.empty())) {
        BemExcitation excitation;
        excitation.name = "Driven";
        const auto factor = std::polar(s.velocityScale, s.phaseDegrees * std::numbers::pi / 180);
        std::set<std::uint32_t> portFaces;
        if (request.coupling)
            for (const auto *element : lem::bemPorts(*request.coupling)) {
                const auto &port = std::get<lem::BemPort>(element->component);
                portFaces.insert(port.faceIds.begin(), port.faceIds.end());
            }
        for (const auto &v : velocities)
            if (!portFaces.contains(v.faceId))
                excitation.faces.push_back({v.faceId, v.normalVelocity * factor});
        if (!excitation.faces.empty())
            request.excitations.push_back(std::move(excitation));
    }
    request.observers.push_back(added(origin, scaled(axis, s.observationDistance)));
    for (const auto transverse : {frame.horizontal, verticalDirection(frame)})
        for (int degrees = -180; degrees <= 180; degrees += 5) {
            const double angle = degrees * std::numbers::pi / 180;
            request.observers.push_back(
                added(origin, added(scaled(axis, s.observationDistance * std::cos(angle)),
                                    scaled(transverse, s.observationDistance * std::sin(angle)))));
        }
    request.observerValidity = exteriorObserverValidity(request);
    return request;
}
} // namespace hornsculpt
