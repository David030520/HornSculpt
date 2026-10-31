#include "axisymmetric_reduction.hpp"
#include "baffle_exterior.hpp"
#include "planning_frequencies.hpp"
#include "project_json.hpp"
#include "simulation_request.hpp"
#include <algorithm>
#include <cmath>
#include <hornsculpt/surface_measures.hpp>
#include <limits>
#include <set>

namespace hornsculpt::project_json {
namespace {
void require(bool value) {
    if (!value)
        invalid(QStringLiteral("Inkonsistentes Simulationsarchiv."));
}
double range(const QJsonValue &v, double low, double high) {
    const double result = number(v);
    require(result >= low && result <= high);
    return result;
}
QString encodeGeometryMode(BemGeometryMode mode) {
    require(mode == BemGeometryMode::threeDimensional || mode == BemGeometryMode::axisymmetric);
    return mode == BemGeometryMode::axisymmetric ? "axisymmetric" : "threeDimensional";
}
BemGeometryMode decodeGeometryMode(const QJsonValue &value) {
    const auto mode = text(value);
    require(mode == "threeDimensional" || mode == "axisymmetric");
    return mode == "axisymmetric" ? BemGeometryMode::axisymmetric
                                  : BemGeometryMode::threeDimensional;
}
QJsonObject encodeAxisymmetricSampling(const AxisymmetricSamplingSettings &s) {
    validateAxisymmetricSamplingSettings(s);
    return {{"maximumLength", s.maximumLength},
            {"chordTolerance", s.chordTolerance},
            {"rotationalTolerance", s.rotationalTolerance},
            {"sectors", double(s.sectors)}};
}
AxisymmetricSamplingSettings decodeAxisymmetricSampling(const QJsonValue &value) {
    const auto o = object(value);
    keys(o, {"maximumLength", "chordTolerance", "rotationalTolerance", "sectors"});
    AxisymmetricSamplingSettings s{positive(o["maximumLength"]), positive(o["chordTolerance"]),
                                   positive(o["rotationalTolerance"]), integer(o["sectors"], 4096)};
    validateAxisymmetricSamplingSettings(s);
    return s;
}
QJsonObject encodeAxisymmetric(const AxisymmetricBem &r) {
    require(r.surfaceContract == AxisymmetricSurfaceContract::triangulatedOrbit ||
            r.surfaceContract == AxisymmetricSurfaceContract::continuousRevolution);
    const auto indices = [](const auto &values) {
        QJsonArray array;
        for (auto value : values)
            array.append(double(value));
        return array;
    };
    return {
        {"origin", encodeVector(r.origin)},
        {"axis", encodeVector(r.axis)},
        {"sectors", double(r.sectors)},
        {"pressureDofs", indices(r.pressureDofs)},
        {"fluxDofs", indices(r.fluxDofs)},
        {"testPanels", indices(r.testPanels)},
        {"surfaceContract", r.surfaceContract == AxisymmetricSurfaceContract::continuousRevolution
                                ? "continuousRevolution"
                                : "triangulatedOrbit"}};
}
AxisymmetricBem decodeAxisymmetric(const QJsonValue &value, const SurfaceMesh &mesh) {
    require(!mesh.vertices.empty() && !mesh.triangles.empty());
    const auto o = object(value);
    auto required = o;
    required.remove("surfaceContract");
    keys(required, {"origin", "axis", "sectors", "pressureDofs", "fluxDofs", "testPanels"});
    AxisymmetricBem r;
    r.surfaceContract = AxisymmetricSurfaceContract::triangulatedOrbit;
    if (o.contains("surfaceContract")) {
        const auto contract = text(o["surfaceContract"]);
        require(contract == "triangulatedOrbit" || contract == "continuousRevolution");
        if (contract == "continuousRevolution")
            r.surfaceContract = AxisymmetricSurfaceContract::continuousRevolution;
    }
    r.origin = vector(o["origin"]);
    r.axis = vector(o["axis"]);
    r.sectors = integer(o["sectors"], 4096);
    for (const auto &id : array(o["pressureDofs"], mesh.vertices.size(), mesh.vertices.size()))
        r.pressureDofs.push_back(integer(id, static_cast<std::uint32_t>(mesh.vertices.size() - 1)));
    for (const auto &id : array(o["fluxDofs"], mesh.triangles.size(), mesh.triangles.size()))
        r.fluxDofs.push_back(integer(id, static_cast<std::uint32_t>(mesh.triangles.size() - 1)));
    for (const auto &id : array(o["testPanels"], mesh.triangles.size(), 1))
        r.testPanels.push_back(integer(id, static_cast<std::uint32_t>(mesh.triangles.size() - 1)));
    bem_detail::validateAxisymmetricReduction(mesh, {}, r);
    return r;
}
QJsonObject encodeMedium(const AcousticMedium &v) {
    return {{"density", v.density}, {"soundSpeed", v.soundSpeed}};
}
AcousticMedium decodeMedium(const QJsonValue &v) {
    const auto o = object(v);
    keys(o, {"density", "soundSpeed"});
    return {positive(o["density"]), positive(o["soundSpeed"])};
}
QJsonObject encodeAngularResponses(const AngularResponseSettings &settings) {
    validateAngularResponseSettings(settings);
    QJsonArray angles;
    for (const double angle : settings.anglesDegrees)
        angles.append(angle);
    return {{"anglesDegrees", angles},
            {"normalizeToAngle", settings.normalizeToAngle},
            {"referenceAngleDegrees", settings.referenceAngleDegrees}};
}
AngularResponseSettings decodeAngularResponses(const QJsonValue &value) {
    const auto o = object(value);
    keys(o, {"anglesDegrees", "normalizeToAngle", "referenceAngleDegrees"});
    AngularResponseSettings settings;
    settings.anglesDegrees.clear();
    for (const auto &angle : array(o["anglesDegrees"], 12, 1))
        settings.anglesDegrees.push_back(number(angle));
    settings.normalizeToAngle = boolean(o["normalizeToAngle"]);
    settings.referenceAngleDegrees = number(o["referenceAngleDegrees"]);
    validateAngularResponseSettings(settings);
    return settings;
}
QJsonObject encodeInfiniteBaffle(const InfiniteBaffle &baffle) {
    validateInfiniteBaffleSettings(baffle);
    QJsonArray faces;
    for (const auto id : baffle.apertureFaceIds)
        faces.append(double(id));
    QJsonObject result{{"origin", encodeVector(baffle.origin)},
            {"normal", encodeVector(baffle.normal)},
            {"apertureFaceIds", faces}};
    if (!baffle.exteriorDrivenMesh.triangles.empty())
        result.insert("exteriorDrivenMesh", encodeMesh(baffle.exteriorDrivenMesh));
    return result;
}
InfiniteBaffle decodeInfiniteBaffle(const QJsonValue &value) {
    const auto o = object(value);
    auto required = o;
    required.remove("exteriorDrivenMesh");
    keys(required, {"origin", "normal", "apertureFaceIds"});
    InfiniteBaffle baffle{vector(o["origin"]), vector(o["normal"]), {}};
    for (const auto &id : array(o["apertureFaceIds"], 2000000, 1))
        baffle.apertureFaceIds.push_back(integer(id));
    if (o.contains("exteriorDrivenMesh"))
        baffle.exteriorDrivenMesh = decodeMesh(o["exteriorDrivenMesh"], false);
    validateInfiniteBaffleSettings(baffle);
    return baffle;
}
QString encodeAcceleration(BemAcceleration acceleration) {
    require(acceleration == BemAcceleration::dense ||
            acceleration == BemAcceleration::fastMultipole);
    return acceleration == BemAcceleration::fastMultipole ? QStringLiteral("fastMultipole")
                                                          : QStringLiteral("dense");
}
BemAcceleration decodeAcceleration(const QJsonValue &v) {
    const auto mode = text(v);
    require(mode == "dense" || mode == "fastMultipole");
    return mode == "fastMultipole" ? BemAcceleration::fastMultipole : BemAcceleration::dense;
}
QJsonObject encodeFastMultipole(const FastMultipoleSettings &s) {
    return {{"interpolationOrder", int(s.interpolationOrder)},
            {"leafPanelCapacity", int(s.leafPanelCapacity)},
            {"separationRatio", s.separationRatio},
            {"maximumIterations", int(s.maximumIterations)},
            {"restartLength", int(s.restartLength)}};
}
FastMultipoleSettings decodeFastMultipole(const QJsonValue &v) {
    const auto o = object(v);
    keys(o, {"interpolationOrder", "leafPanelCapacity", "separationRatio", "maximumIterations",
             "restartLength"});
    FastMultipoleSettings s;
    s.interpolationOrder = integer(o["interpolationOrder"], 8);
    s.leafPanelCapacity = integer(o["leafPanelCapacity"], 256);
    s.separationRatio = range(o["separationRatio"], 1.5, 10);
    s.maximumIterations = integer(o["maximumIterations"], 10000);
    s.restartLength = integer(o["restartLength"], 300);
    require(s.interpolationOrder >= 3 && s.leafPanelCapacity >= 4 && s.maximumIterations >= 1 &&
            s.restartLength >= 2);
    return s;
}
QJsonObject encodeNumerics(const BemNumerics &v) {
    return {{"quadratureOrder", int(v.quadratureOrder)},
            {"singularOrder", int(v.singularOrder)},
            {"nearDistanceRatio", v.nearDistanceRatio},
            {"nearRelativeTolerance", v.nearRelativeTolerance},
            {"nearMaximumDepth", int(v.nearMaximumDepth)},
            {"relativeResidualTolerance", v.relativeResidualTolerance},
            {"maximumRefinementSteps", int(v.maximumRefinementSteps)},
            {"acceleration", encodeAcceleration(v.acceleration)},
            {"fastMultipole", encodeFastMultipole(v.fastMultipole)}};
}
BemNumerics decodeNumerics(const QJsonValue &v) {
    const auto o = object(v);
    auto required = o;
    required.remove("acceleration");
    required.remove("fastMultipole");
    keys(required,
         {"quadratureOrder", "singularOrder", "nearDistanceRatio", "nearRelativeTolerance",
          "nearMaximumDepth", "relativeResidualTolerance", "maximumRefinementSteps"});
    BemNumerics n;
    n.quadratureOrder = integer(o["quadratureOrder"], 4);
    n.singularOrder = integer(o["singularOrder"], 20);
    require((n.quadratureOrder == 2 || n.quadratureOrder == 4) && n.singularOrder >= 2);
    n.nearDistanceRatio = range(o["nearDistanceRatio"], 0, 10);
    n.nearRelativeTolerance = range(o["nearRelativeTolerance"], 1e-10, .01);
    n.nearMaximumDepth = integer(o["nearMaximumDepth"], 10);
    n.relativeResidualTolerance = range(o["relativeResidualTolerance"], 1e-7, .01);
    n.maximumRefinementSteps = integer(o["maximumRefinementSteps"], 10);
    if (o.contains("acceleration"))
        n.acceleration = decodeAcceleration(o["acceleration"]);
    if (o.contains("fastMultipole"))
        n.fastMultipole = decodeFastMultipole(o["fastMultipole"]);
    return n;
}
QJsonArray encodePressure(const std::vector<std::complex<double>> &values) {
    QJsonArray result;
    for (auto value : values)
        result.append(encodeComplex(value));
    return result;
}
QJsonObject encodeFieldDisplay(const FieldDisplaySettings &s) {
    validateFieldDisplaySettings(s);
    return {{"vertical", s.vertical},
            {"minimumLateral", s.minimumLateral},
            {"maximumLateral", s.maximumLateral},
            {"minimumAxial", s.minimumAxial},
            {"maximumAxial", s.maximumAxial},
            {"lateralCount", int(s.lateralCount)},
            {"axialCount", int(s.axialCount)},
            {"quantity", s.quantity},
            {"minimumDb", s.minimumDb},
            {"maximumDb", s.maximumDb},
            {"phaseDegrees", s.phaseDegrees},
            {"exclusionDistance", s.exclusionDistance},
            {"frequencyHz", s.frequencyHz}, {"excitationIndex", s.excitationIndex}};
}
FieldDisplaySettings decodeFieldDisplay(const QJsonValue &value) {
    const auto o = object(value);
    auto required = o;
    required.remove("frequencyHz");
    required.remove("excitationIndex");
    keys(required, {"vertical", "minimumLateral", "maximumLateral", "minimumAxial", "maximumAxial",
             "lateralCount", "axialCount", "quantity", "minimumDb", "maximumDb", "phaseDegrees",
             "exclusionDistance"});
    FieldDisplaySettings s{boolean(o["vertical"]),        number(o["minimumLateral"]),
                           number(o["maximumLateral"]),   number(o["minimumAxial"]),
                           number(o["maximumAxial"]),     integer(o["lateralCount"], 201),
                           integer(o["axialCount"], 201), int(integer(o["quantity"], 2)),
                           number(o["minimumDb"]),        number(o["maximumDb"]),
                           number(o["phaseDegrees"]),     positive(o["exclusionDistance"])};
    if (o.contains("frequencyHz")) s.frequencyHz = number(o["frequencyHz"]);
    if (o.contains("excitationIndex")) {
        const double selection = range(o["excitationIndex"], -1, 10000);
        require(std::floor(selection) == selection);
        s.excitationIndex = static_cast<int>(selection);
    }
    validateFieldDisplaySettings(s);
    return s;
}
QJsonObject encodePhaseDisplay(const PhaseDisplaySettings &s) {
    validatePhaseDisplaySettings(s);
    return {{"angleDegrees", s.angleDegrees},
            {"referenceAngleDegrees", s.referenceAngleDegrees},
            {"delaySeconds", s.delaySeconds},
            {"relative", s.relative},
            {"unwrap", s.unwrap}};
}
PhaseDisplaySettings decodePhaseDisplay(const QJsonValue &value) {
    const auto o = object(value);
    keys(o, {"angleDegrees", "referenceAngleDegrees", "delaySeconds", "relative", "unwrap"});
    PhaseDisplaySettings s{number(o["angleDegrees"]), number(o["referenceAngleDegrees"]),
                           number(o["delaySeconds"]), boolean(o["relative"]), boolean(o["unwrap"])};
    validatePhaseDisplaySettings(s);
    return s;
}
std::vector<std::complex<double>> decodePressure(const QJsonValue &v, std::size_t count) {
    const auto a = array(v, static_cast<qsizetype>(count), static_cast<qsizetype>(count));
    std::vector<std::complex<double>> result;
    result.reserve(count);
    for (const auto &value : a)
        result.push_back(complex(value));
    return result;
}
QJsonArray encodeVelocities(const std::vector<FaceVelocity> &velocities) {
    QJsonArray result;
    for (const auto &v : velocities)
        result.append(QJsonObject{{"faceId", double(v.faceId)},
                                  {"normalVelocity", encodeComplex(v.normalVelocity)}});
    return result;
}
std::vector<FaceVelocity> decodeVelocities(const QJsonValue &value, const SurfaceMesh &mesh,
                                           const std::set<std::uint32_t> &interfaces = {}) {
    auto known = meshFaces(mesh);
    std::set<std::uint32_t> seen;
    std::vector<FaceVelocity> result;
    for (const auto &v : array(value, known.size())) {
        const auto o = object(v);
        keys(o, {"faceId", "normalVelocity"});
        FaceVelocity face{integer(o["faceId"]), complex(o["normalVelocity"])};
        require(known.contains(face.faceId) &&
                (known.at(face.faceId) == BoundaryCondition::driven ||
                 interfaces.contains(face.faceId)) &&
                seen.insert(face.faceId).second);
        result.push_back(face);
    }
    return result;
}
QJsonObject encodeRequest(const BemRequest &r) {
    QJsonArray frequencies, observers, excitations;
    for (auto f : r.frequencies)
        frequencies.append(f);
    for (auto p : r.observers)
        observers.append(encodeVector(p));
    for (const auto &e : r.excitations) {
        QJsonArray faces;
        for (const auto &f : e.faces)
            faces.append(QJsonObject{{"faceId", double(f.faceId)},
                                     {"normalVelocity", encodeComplex(f.normalVelocity)}});
        excitations.append(QJsonObject{{"name", QString::fromStdString(e.name)}, {"faces", faces}});
    }
    QJsonObject result{{"medium", encodeMedium(r.medium)},
                       {"frequencies", frequencies},
                       {"observers", observers},
                       {"excitations", excitations},
                       {"numerics", encodeNumerics(r.numerics)}};
    if (r.farField)
        result.insert("farField",
                      QJsonObject{{"origin", encodeVector(r.farField->origin)},
                                  {"referenceDistance", r.farField->referenceDistance}});
    if (r.coupling)
        result.insert("coupling", encodeLemNetwork(*r.coupling));
    if (r.infiniteBaffle)
        result.insert("infiniteBaffle", encodeInfiniteBaffle(*r.infiniteBaffle));
    if (r.axisymmetric)
        result.insert("axisymmetric", encodeAxisymmetric(*r.axisymmetric));
    if (!r.interactionPlanFrequencies.empty()) {
        QJsonArray raster;
        for (double frequency : r.interactionPlanFrequencies)
            raster.append(frequency);
        result.insert("interactionPlanFrequencies", raster);
    }
    return result;
}
BemRequest decodeRequest(const QJsonValue &requestValue, const SurfaceMesh &mesh) {
    const auto o = object(requestValue);
    auto required = o;
    required.remove("interactionPlanFrequencies");
    required.remove("farField");
    required.remove("coupling");
    required.remove("infiniteBaffle");
    required.remove("axisymmetric");
    keys(required, {"medium", "frequencies", "observers", "excitations", "numerics"});
    BemRequest r;
    r.medium = decodeMedium(o["medium"]);
    r.numerics = decodeNumerics(o["numerics"]);
    if (o.contains("axisymmetric")) {
        r.axisymmetric = decodeAxisymmetric(o["axisymmetric"], mesh);
        require(r.numerics.acceleration == BemAcceleration::dense);
    }
    // The file byte limit bounds materialized arrays; no arbitrary sweep-count cap.
    for (const auto &f : array(o["frequencies"], std::numeric_limits<int>::max(), 1))
        r.frequencies.push_back(positive(f));
    if (o.contains("interactionPlanFrequencies"))
        for (const auto &f :
             array(o["interactionPlanFrequencies"], std::numeric_limits<int>::max(), 1))
            r.interactionPlanFrequencies.push_back(positive(f));
    bem_detail::planningWavenumber(r.frequencies, r.interactionPlanFrequencies,
                                   r.medium.soundSpeed);
    for (const auto &p : array(o["observers"], 100000))
        r.observers.push_back(vector(p));
    if (o.contains("farField")) {
        const auto farField = object(o["farField"]);
        keys(farField, {"origin", "referenceDistance"});
        r.farField = FarFieldObservation{vector(farField["origin"]),
                                         positive(farField["referenceDistance"])};
        for (const auto point : r.observers) {
            const auto origin = r.farField->origin;
            const double distance =
                std::hypot(point.x - origin.x, point.y - origin.y, point.z - origin.z);
            require(std::isfinite(distance) && distance > 0);
        }
    }
    auto known = meshFaces(mesh);
    if (o.contains("infiniteBaffle")) {
        r.infiniteBaffle = decodeInfiniteBaffle(o["infiniteBaffle"]);
        require(r.numerics.acceleration == BemAcceleration::dense);
        for (const auto id : r.infiniteBaffle->apertureFaceIds)
            require(known.contains(id) && known.at(id) == BoundaryCondition::wall);
        const auto &b = *r.infiniteBaffle;
        bem_detail::validateBaffleExteriorMesh(mesh, b.exteriorDrivenMesh, b.origin, b.normal);
        require(!r.axisymmetric || b.exteriorDrivenMesh.triangles.empty());
        const auto externalFaces = meshFaces(b.exteriorDrivenMesh);
        known.insert(externalFaces.begin(), externalFaces.end());
    }
    require(!r.axisymmetric || r.infiniteBaffle.has_value());
    std::set<std::uint32_t> portFaces;
    if (o.contains("coupling")) {
        r.coupling = decodeLemNetwork(o["coupling"]);
        lem::validateNetwork(*r.coupling);
        const auto ports = lem::bemPorts(*r.coupling);
        require(!ports.empty());
        for (const auto *element : ports) {
            const auto &port = std::get<lem::BemPort>(element->component);
            const bool interfacePort = port.kind == lem::BemPortKind::interfacePort;
            require(!interfacePort || r.infiniteBaffle.has_value());
            for (const auto id : port.faceIds) {
                require(known.contains(id) &&
                        known.at(id) ==
                            (interfacePort ? BoundaryCondition::wall : BoundaryCondition::driven) &&
                        portFaces.insert(id).second);
                if (r.infiniteBaffle)
                    require(std::find(r.infiniteBaffle->apertureFaceIds.begin(),
                                      r.infiniteBaffle->apertureFaceIds.end(),
                                      id) == r.infiniteBaffle->apertureFaceIds.end());
            }
        }
        require(!portFaces.empty());
        for (const auto &[id, boundary] : known)
            require(r.infiniteBaffle || boundary != BoundaryCondition::driven ||
                    portFaces.contains(id));
    }
    for (const auto &v : array(o["excitations"], r.coupling ? (r.infiniteBaffle ? 1 : 0) : 1024,
                               r.coupling ? 0 : 1)) {
        const auto entry = object(v);
        keys(entry, {"name", "faces"});
        BemExcitation e;
        e.name = text(entry["name"], true).toStdString();
        std::set<std::uint32_t> assigned;
        bool nonzero = false;
        for (const auto &faceValue :
             array(entry["faces"], static_cast<qsizetype>(known.size()), 1)) {
            const auto f = object(faceValue);
            keys(f, {"faceId", "normalVelocity"});
            FaceVelocity velocity{integer(f["faceId"]), complex(f["normalVelocity"])};
            const auto found = known.find(velocity.faceId);
            require(found != known.end() && found->second == BoundaryCondition::driven &&
                    assigned.insert(velocity.faceId).second &&
                    !portFaces.contains(velocity.faceId));
            nonzero |= std::abs(velocity.normalVelocity) > 0;
            e.faces.push_back(velocity);
        }
        require(nonzero || (r.coupling && r.infiniteBaffle));
        r.excitations.push_back(std::move(e));
    }
    return r;
}
QJsonObject encodeResult(const BemResult &r) {
    QJsonArray frequencies;
    for (const auto &f : r.frequencies) {
        QJsonArray excitations;
        for (const auto &e : f.excitations) {
            require(bool(e.boundary()));
            QJsonArray impedances;
            for (const auto &z : e.impedances)
                impedances.append(
                    QJsonObject{{"faceId", double(z.faceId)},
                                {"area", z.area},
                                {"meanPressure", encodeComplex(z.meanPressure)},
                                {"volumeVelocity", encodeComplex(z.volumeVelocity)},
                                {"volumeImpedance", encodeComplex(z.volumeImpedance)}});
            QJsonObject excitation{{"name", QString::fromStdString(e.name)},
                                   {"surfacePressure", encodePressure(e.boundary()->surfacePressure)},
                                   {"observerPressure", encodePressure(e.observerPressure)},
                                   {"impedances", impedances},
                                   {"refinementSteps", int(e.refinementSteps)},
                                   {"iterations", int(e.iterations)},
                                   {"boundaryResidual", e.boundaryResidual}};
            if (e.lemSolution) {
                excitation.insert("lemSolution", encodeLemSolution(*e.lemSolution));
                excitation.insert("solvedVelocities", encodeVelocities(e.solvedVelocities));
            }
            if (!e.boundary()->exteriorSurfacePressure.empty()) {
                excitation.insert("exteriorSurfacePressure", encodePressure(e.boundary()->exteriorSurfacePressure));
                excitation.insert("solvedVelocities", encodeVelocities(e.solvedVelocities));
            }
            require(e.boundary()->pressureTrace.empty() == e.boundary()->normalDerivative.empty());
            if (!e.boundary()->pressureTrace.empty()) {
                excitation.insert("pressureTrace", encodePressure(e.boundary()->pressureTrace));
                excitation.insert("normalDerivative", encodePressure(e.boundary()->normalDerivative));
            }
            excitations.append(excitation);
        }
        frequencies.append(QJsonObject{
            {"frequency", f.frequency},
            {"excitations", excitations},
            {"assemblyMilliseconds", f.assemblyMilliseconds},
            {"solveMilliseconds", f.solveMilliseconds},
            {"factorizationMilliseconds", f.factorizationMilliseconds},
            {"triangularSolveMilliseconds", f.triangularSolveMilliseconds},
            {"fieldMilliseconds", f.fieldMilliseconds},
            {"residualCheckMilliseconds", f.residualCheckMilliseconds},
            {"interactionPlanMilliseconds", f.interactionPlanMilliseconds},
            {"interactionPlanReused", f.interactionPlanReused},
            {"fastMultipole",
             QJsonObject{{"nodeCount", encodeWide(f.fastMultipole.nodeCount)},
                         {"leafCount", encodeWide(f.fastMultipole.leafCount)},
                         {"farInteractionCount", encodeWide(f.fastMultipole.farInteractionCount)},
                         {"directPanelPairCount", encodeWide(f.fastMultipole.directPanelPairCount)},
                         {"representedFarPanelPairs",
                          encodeWide(f.fastMultipole.representedFarPanelPairs)}}}});
    }
    return {{"device", QString::fromStdString(r.device)},
            {"geometryMode", encodeGeometryMode(r.geometryMode)},
            {"acceleration", encodeAcceleration(r.acceleration)},
            {"triangleCount", double(r.triangleCount)},
            {"pressureDofCount", double(r.pressureDofCount)},
            {"peakDeviceBytes", encodeWide(r.peakDeviceBytes)},
            {"frequencies", frequencies}};
}
struct FaceIntegral {
    double area{};
    std::complex<double> pressure{};
    double absolutePressure{};
};
bool close(std::complex<double> a, std::complex<double> b) {
    if (!std::isfinite(a.real()) || !std::isfinite(a.imag()) || !std::isfinite(b.real()) ||
        !std::isfinite(b.imag()))
        return false;
    return std::abs(a - b) <= 1e-10 * std::max({std::abs(a), std::abs(b), 1e-300});
}
BemResult decodeResult(const QJsonValue &resultValue, const BemRequest &request,
                       const SurfaceMesh &mesh) {
    const auto o = object(resultValue);
    auto requiredResult = o;
    requiredResult.remove("acceleration");
    requiredResult.remove("geometryMode");
    keys(requiredResult,
         {"device", "triangleCount", "pressureDofCount", "peakDeviceBytes", "frequencies"});
    BemResult r;
    if (o.contains("geometryMode"))
        r.geometryMode = decodeGeometryMode(o["geometryMode"]);
    require((r.geometryMode == BemGeometryMode::axisymmetric) == request.axisymmetric.has_value());
    if (o.contains("acceleration"))
        r.acceleration = decodeAcceleration(o["acceleration"]);
    require(r.acceleration == request.numerics.acceleration);
    r.device = text(o["device"], true).toStdString();
    r.triangleCount = integer(o["triangleCount"], 2000000);
    r.pressureDofCount = integer(o["pressureDofCount"], 1000000);
    const auto bytes = wideInteger(o["peakDeviceBytes"]);
    require(bytes <= std::numeric_limits<std::size_t>::max());
    r.peakDeviceBytes = static_cast<std::size_t>(bytes);
    const auto pressureCount =
        request.axisymmetric
            ? bem_detail::validateAxisymmetricReduction(mesh, {}, *request.axisymmetric)
                  .pressureCount
            : mesh.vertices.size();
    require(r.triangleCount == mesh.triangles.size() && r.pressureDofCount == pressureCount);
    const auto measures = surfaceMeasures(mesh, request.axisymmetric);
    const auto exteriorMesh = request.infiniteBaffle ? request.infiniteBaffle->exteriorDrivenMesh : SurfaceMesh{};
    const auto combinedMesh = bem_detail::appendBaffleExterior(mesh, exteriorMesh);
    auto combinedAreas = measures.areas;
    if (!exteriorMesh.triangles.empty()) {
        const auto externalAreas = surfaceMeasures(exteriorMesh).areas;
        combinedAreas.insert(combinedAreas.end(), externalAreas.begin(), externalAreas.end());
    }
    const auto fs = array(o["frequencies"], request.frequencies.size(), request.frequencies.size());
    for (qsizetype fi = 0; fi < fs.size(); ++fi) {
        const auto entry = object(fs[fi]);
        const bool multipole = r.acceleration == BemAcceleration::fastMultipole;
        if (multipole)
            require(entry.contains("fastMultipole") && entry.contains("residualCheckMilliseconds"));
        auto requiredFrequency = entry;
        requiredFrequency.remove("residualCheckMilliseconds");
        requiredFrequency.remove("fastMultipole");
        requiredFrequency.remove("interactionPlanMilliseconds");
        requiredFrequency.remove("interactionPlanReused");
        keys(requiredFrequency,
             {"frequency", "excitations", "assemblyMilliseconds", "solveMilliseconds",
              "factorizationMilliseconds", "triangularSolveMilliseconds", "fieldMilliseconds"});
        BemFrequencyResult f;
        f.frequency = positive(entry["frequency"]);
        require(f.frequency == request.frequencies[fi]);
        f.assemblyMilliseconds = range(entry["assemblyMilliseconds"], 0, 1e15);
        f.solveMilliseconds = range(entry["solveMilliseconds"], 0, 1e15);
        f.factorizationMilliseconds = range(entry["factorizationMilliseconds"], 0, 1e15);
        f.triangularSolveMilliseconds = range(entry["triangularSolveMilliseconds"], 0, 1e15);
        f.fieldMilliseconds = range(entry["fieldMilliseconds"], 0, 1e15);
        if (entry.contains("residualCheckMilliseconds"))
            f.residualCheckMilliseconds = range(entry["residualCheckMilliseconds"], 0, 1e15);
        if (entry.contains("interactionPlanMilliseconds"))
            f.interactionPlanMilliseconds =
                range(entry["interactionPlanMilliseconds"], 0, f.assemblyMilliseconds);
        if (entry.contains("interactionPlanReused"))
            f.interactionPlanReused = boolean(entry["interactionPlanReused"]);
        require(!f.interactionPlanReused || f.interactionPlanMilliseconds == 0);
        if (entry.contains("fastMultipole")) {
            const auto stats = object(entry["fastMultipole"]);
            keys(stats, {"nodeCount", "leafCount", "farInteractionCount", "directPanelPairCount",
                         "representedFarPanelPairs"});
            const auto count = [](const QJsonValue &value) {
                const auto result = wideInteger(value);
                require(result <= std::numeric_limits<std::size_t>::max());
                return static_cast<std::size_t>(result);
            };
            f.fastMultipole = {count(stats["nodeCount"]), count(stats["leafCount"]),
                               count(stats["farInteractionCount"]),
                               count(stats["directPanelPairCount"]),
                               count(stats["representedFarPanelPairs"])};
            const auto panels = mesh.triangles.size();
            require(panels > 0 && panels <= std::numeric_limits<std::size_t>::max() / panels);
            const auto pairs = panels * panels;
            require(f.fastMultipole.leafCount <= f.fastMultipole.nodeCount &&
                    f.fastMultipole.directPanelPairCount <= pairs &&
                    f.fastMultipole.representedFarPanelPairs <= pairs);
            if (multipole)
                // Subtraction verifies the partition without overflowing an addition.
                require(f.fastMultipole.leafCount >= 1 &&
                        f.fastMultipole.directPanelPairCount ==
                            pairs - f.fastMultipole.representedFarPanelPairs);
        }
        const auto excitationCount = request.coupling ? 1 : request.excitations.size();
        const auto es = array(entry["excitations"], excitationCount, excitationCount);
        for (qsizetype ei = 0; ei < es.size(); ++ei) {
            const auto item = object(es[ei]);
            if (multipole)
                require(item.contains("iterations"));
            auto requiredItem = item;
            requiredItem.remove("iterations");
            requiredItem.remove("pressureTrace");
            requiredItem.remove("normalDerivative");
            requiredItem.remove("exteriorSurfacePressure");
            if (!exteriorMesh.triangles.empty())
                requiredItem.remove("solvedVelocities");
            require(item.contains("pressureTrace") == item.contains("normalDerivative"));
            if (request.coupling) {
                require(item.contains("lemSolution") && item.contains("solvedVelocities"));
                requiredItem.remove("lemSolution");
                requiredItem.remove("solvedVelocities");
            }
            keys(requiredItem, {"name", "surfacePressure", "observerPressure", "impedances",
                                "refinementSteps", "boundaryResidual"});
            BemExcitationResult e;
            e.name = text(item["name"], true).toStdString();
            require(e.name == (request.coupling ? "LEM" : request.excitations[ei].name));
            if (request.coupling) {
                e.lemSolution =
                    decodeLemSolution(item["lemSolution"], *request.coupling, f.frequency);
                std::set<std::uint32_t> interfaces;
                for (const auto *element : lem::bemPorts(*request.coupling)) {
                    const auto &port = std::get<lem::BemPort>(element->component);
                    if (port.kind == lem::BemPortKind::interfacePort)
                        interfaces.insert(port.faceIds.begin(), port.faceIds.end());
                }
                e.solvedVelocities = decodeVelocities(item["solvedVelocities"], combinedMesh, interfaces);
                std::map<std::uint32_t, std::complex<double>> velocities;
                for (const auto &v : e.solvedVelocities)
                    velocities.emplace(v.faceId, v.normalVelocity);
                const auto ports = lem::bemPorts(*request.coupling);
                std::set<std::uint32_t> expected;
                for (const auto &excitation : request.excitations)
                    for (const auto &known : excitation.faces) {
                        require(expected.insert(known.faceId).second &&
                                velocities.contains(known.faceId) &&
                                close(velocities.at(known.faceId), known.normalVelocity));
                    }
                for (std::size_t i = 0; i < ports.size(); ++i) {
                    const auto &port = std::get<lem::BemPort>(ports[i]->component);
                    double area = 0;
                    for (std::size_t ti = 0; ti < combinedMesh.triangles.size(); ++ti) {
                        const auto &triangle = combinedMesh.triangles[ti];
                        if (std::find(port.faceIds.begin(), port.faceIds.end(), triangle.faceId) !=
                            port.faceIds.end())
                            area +=
                                combinedAreas[ti] * (port.velocityWeights.contains(triangle.faceId)
                                                          ? port.velocityWeights.at(triangle.faceId)
                                                          : 1.);
                    }
                    require(std::isfinite(area) && area > 0);
                    const auto &solution = e.lemSolution->ports.at(i);
                    require(solution.id == ports[i]->id);
                    for (const auto id : port.faceIds) {
                        expected.insert(id);
                        const auto weight =
                            port.velocityWeights.contains(id) ? port.velocityWeights.at(id) : 1.;
                        require(
                            velocities.contains(id) &&
                            close(velocities.at(id), (weight / area) * solution.volumeVelocity));
                    }
                }
                require(expected.size() == velocities.size());
            }
            BemBoundaryData boundary;
            boundary.surfacePressure = decodePressure(item["surfacePressure"], mesh.triangles.size());
            if (!exteriorMesh.triangles.empty()) {
                boundary.exteriorSurfacePressure = decodePressure(item["exteriorSurfacePressure"], exteriorMesh.triangles.size());
                if (!request.coupling)
                    e.solvedVelocities = decodeVelocities(item["solvedVelocities"], combinedMesh, {});
            } else require(!item.contains("exteriorSurfacePressure"));
            if (item.contains("pressureTrace")) {
                boundary.pressureTrace = decodePressure(item["pressureTrace"], mesh.vertices.size());
                boundary.normalDerivative =
                    decodePressure(item["normalDerivative"], mesh.triangles.size());
                if (request.axisymmetric)
                    bem_detail::validateAxisymmetricTrace(*request.axisymmetric, boundary.pressureTrace,
                                                          boundary.normalDerivative);
                for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
                    const auto &nodes = mesh.triangles[t].vertices;
                    const auto &weights = measures.pressureIntegrals[t];
                    const bool continuous = request.axisymmetric &&
                                            request.axisymmetric->surfaceContract ==
                                                AxisymmetricSurfaceContract::continuousRevolution;
                    const auto mean =
                        continuous
                            ? (weights[0] / measures.areas[t]) * boundary.pressureTrace[nodes[0]] +
                                  (weights[1] / measures.areas[t]) * boundary.pressureTrace[nodes[1]] +
                                  (weights[2] / measures.areas[t]) * boundary.pressureTrace[nodes[2]]
                            : (boundary.pressureTrace[nodes[0]] + boundary.pressureTrace[nodes[1]] +
                               boundary.pressureTrace[nodes[2]]) /
                                  3.;
                    // Continuous cell weights can differ by coordinate roundoff
                    // after centering. Cancellation of nodal pressures must not
                    // turn that roundoff into an apparent relative corruption.
                    constexpr double relativeArchivedPressureMeanTolerance = 1e-10;
                    double scale = std::abs(boundary.surfacePressure[t]);
                    double absoluteMean = 0;
                    for (unsigned i = 0; i < 3; ++i)
                        absoluteMean +=
                            (weights[i] / measures.areas[t]) * std::abs(boundary.pressureTrace[nodes[i]]);
                    scale = std::max(scale, absoluteMean);
                    require(std::isfinite(scale) &&
                            std::abs(mean - boundary.surfacePressure[t]) <=
                                relativeArchivedPressureMeanTolerance * scale);
                }
            }
            e.observerPressure = decodePressure(item["observerPressure"], request.observers.size());
            const unsigned couplingPasses =
                e.lemSolution ? static_cast<unsigned>(e.lemSolution->ports.size()) + 1u : 1u;
            e.refinementSteps = integer(item["refinementSteps"],
                                        request.numerics.maximumRefinementSteps * couplingPasses);
            if (item.contains("iterations"))
                e.iterations =
                    integer(item["iterations"],
                            request.numerics.acceleration == BemAcceleration::fastMultipole
                                ? request.numerics.fastMultipole.maximumIterations *
                                      (request.numerics.maximumRefinementSteps + 1) * couplingPasses
                                : 0);
            e.boundaryResidual =
                range(item["boundaryResidual"], 0, request.numerics.relativeResidualTolerance);
            std::map<std::uint32_t, std::complex<double>> assigned;
            for (const auto &face :
                 request.coupling ? e.solvedVelocities : request.excitations[ei].faces)
                if (std::abs(face.normalVelocity) > 0)
                    assigned.emplace(face.faceId, face.normalVelocity);
            std::map<std::uint32_t, FaceIntegral> integrals;
            for (std::size_t ti = 0; ti < mesh.triangles.size(); ++ti) {
                const auto &triangle = mesh.triangles[ti];
                if (!assigned.contains(triangle.faceId))
                    continue;
                const double area = measures.areas[ti];
                auto &sum = integrals[triangle.faceId];
                sum.area += area;
                sum.pressure += area * boundary.surfacePressure[ti];
                sum.absolutePressure += area * std::abs(boundary.surfacePressure[ti]);
            }
            std::set<std::uint32_t> seen;
            for (std::size_t ti = 0; ti < exteriorMesh.triangles.size(); ++ti) {
                const auto &triangle = exteriorMesh.triangles[ti];
                if (!assigned.contains(triangle.faceId))
                    continue;
                const double area = combinedAreas[mesh.triangles.size() + ti];
                auto &sum = integrals[triangle.faceId];
                sum.area += area;
                sum.pressure += area * boundary.exteriorSurfacePressure[ti];
                sum.absolutePressure += area * std::abs(boundary.exteriorSurfacePressure[ti]);
            }
            for (const auto &v : array(item["impedances"], assigned.size(), assigned.size())) {
                const auto z = object(v);
                keys(z, {"faceId", "area", "meanPressure", "volumeVelocity", "volumeImpedance"});
                FaceImpedance impedance{integer(z["faceId"]), positive(z["area"]),
                                        complex(z["meanPressure"]), complex(z["volumeVelocity"]),
                                        complex(z["volumeImpedance"])};
                require(assigned.contains(impedance.faceId) &&
                        seen.insert(impedance.faceId).second);
                const auto &sum = integrals.at(impedance.faceId);
                const double pressureScale =
                    std::max({std::abs(impedance.meanPressure), std::abs(sum.pressure / sum.area),
                              sum.absolutePressure / sum.area});
                require(std::isfinite(pressureScale) && close(impedance.area, sum.area) &&
                        std::abs(impedance.meanPressure - sum.pressure / sum.area) <=
                            1e-10 * pressureScale &&
                        close(impedance.volumeVelocity, assigned.at(impedance.faceId) * sum.area) &&
                        close(impedance.volumeImpedance,
                              impedance.meanPressure / impedance.volumeVelocity));
                e.impedances.push_back(impedance);
            }
            if (request.coupling) {
                const auto ports = lem::bemPorts(*request.coupling);
                for (std::size_t i = 0; i < ports.size(); ++i) {
                    const auto &port = std::get<lem::BemPort>(ports[i]->component);
                    auto allPressure = boundary.surfacePressure;
                    allPressure.insert(allPressure.end(), boundary.exteriorSurfacePressure.begin(), boundary.exteriorSurfacePressure.end());
                    double area = 0, magnitude = 0;
                    std::complex<double> pressure{};
                    for (std::size_t ti = 0; ti < combinedMesh.triangles.size(); ++ti) {
                        const auto &triangle = combinedMesh.triangles[ti];
                        if (std::find(port.faceIds.begin(), port.faceIds.end(), triangle.faceId) !=
                            port.faceIds.end())
                            area +=
                                combinedAreas[ti] * (port.velocityWeights.contains(triangle.faceId)
                                                          ? port.velocityWeights.at(triangle.faceId)
                                                          : 1.);
                    }
                    require(std::isfinite(area) && area > 0);
                    for (std::size_t ti = 0; ti < combinedMesh.triangles.size(); ++ti)
                        if (std::find(port.faceIds.begin(), port.faceIds.end(),
                                      combinedMesh.triangles[ti].faceId) != port.faceIds.end()) {
                            // Normalize before multiplying by pressure: a global rescaling
                            // of the modal weights must not overflow an unchanged field.
                            const double weight =
                                port.velocityWeights.contains(combinedMesh.triangles[ti].faceId)
                                    ? port.velocityWeights.at(combinedMesh.triangles[ti].faceId)
                                    : 1.;
                            const double a = (weight / area) * combinedAreas[ti];
                            pressure += a * allPressure[ti];
                            magnitude += a * std::abs(allPressure[ti]);
                        }
                    // This is the solver's independent port-pressure closure acceptance limit,
                    // not a pressure-error estimate derived from the Galerkin residual.
                    require(std::abs(e.lemSolution->ports[i].pressure - pressure) <=
                            relativePortPressureClosureTolerance *
                                std::max(magnitude, std::abs(e.lemSolution->ports[i].pressure)));
                }
            }
            e.setBoundary(std::move(boundary));
            f.excitations.push_back(std::move(e));
        }
        r.frequencies.push_back(std::move(f));
    }
    return r;
}
} // namespace

QJsonObject encodeSimulation(const SimulationArchive &a) {
    require(!a.resultComplete.has_value() || bool(a.result));
    require(!a.current || a.resultComplete.value_or(true));
    const auto &s = a.settings;
    QJsonArray measurements;
    for (const auto &curve : a.measurements) {
        validateMeasuredResponse(curve);
        QJsonArray samples;
        for (const auto &p : curve.samples)
            samples.append(QJsonArray{p.frequencyHz, p.decibels});
        measurements.append(QJsonObject{{"name", curve.name},
                                        {"absoluteSpl", curve.absoluteSpl},
                                        {"offsetDb", curve.offsetDb},
                                        {"samples", samples}});
    }
    QJsonObject result{
        {"measurements", measurements},
        {"settings", QJsonObject{{"minimumFrequency", s.minimumFrequency},
                                 {"maximumFrequency", s.maximumFrequency},
                                 {"frequencyCount", double(s.frequencyCount)},
                                 {"medium", encodeMedium(s.medium)},
                                 {"velocityScale", s.velocityScale},
                                 {"phaseDegrees", s.phaseDegrees},
                                 {"acceleration", encodeAcceleration(s.acceleration)},
                                 {"fastMultipole", encodeFastMultipole(s.fastMultipole)},
                                 {"nearMaximumDepth", int(s.nearMaximumDepth)},
                                 {"normalizeToOnAxis", s.normalizeToOnAxis},
                                 {"angularResponses", encodeAngularResponses(s.angularResponses)},
                                 {"fieldDisplay", encodeFieldDisplay(s.fieldDisplay)},
                                 {"phaseDisplay", encodePhaseDisplay(s.phaseDisplay)},
                                 {"normAngle", s.normalizationAngleDegrees},
                                 {"DbScale", QJsonObject{{"minimum", s.directivityMinimumDb},
                                                         {"maximum", s.directivityMaximumDb},
                                                         {"step", s.directivityStepDb}}},
                                 {"beamwidthTarget",
                                  QJsonObject{{"enabled", s.beamwidthTargetEnabled},
                                              {"horizontalDegrees", s.targetHorizontalDegrees},
                                              {"verticalDegrees", s.targetVerticalDegrees},
                                              {"toleranceDegrees", s.targetToleranceDegrees}}}}},
        {"revision", encodeWide(a.revision)},
        {"current", a.current},
        {"solvedMesh", a.mesh ? QJsonValue(encodeMesh(*a.mesh)) : QJsonValue(QJsonValue::Null)},
        {"request",
         a.request ? QJsonValue(encodeRequest(*a.request)) : QJsonValue(QJsonValue::Null)},
        {"result", a.result ? QJsonValue(encodeResult(*a.result)) : QJsonValue(QJsonValue::Null)}};
    if (a.result && a.resultComplete)
        result.insert("resultComplete", *a.resultComplete);
    auto couplingSettings = result["settings"].toObject();
    require(s.ctaExcitationIndex >= -1 && s.ctaExcitationIndex <= 10000);
    couplingSettings.insert("ctaExcitationIndex", s.ctaExcitationIndex);
    couplingSettings.insert("couplingEnabled", s.couplingEnabled);
    couplingSettings.insert("couplingNetwork", encodeLemNetwork(s.couplingNetwork));
    couplingSettings.insert("farField", s.farField);
    couplingSettings.insert("observationDistance", s.observationDistance);
    couplingSettings.insert("geometryMode", encodeGeometryMode(s.geometryMode));
    couplingSettings.insert("axisymmetricSampling",
                            encodeAxisymmetricSampling(s.axisymmetricSampling));
    if (s.infiniteBaffle)
        couplingSettings.insert("infiniteBaffle", encodeInfiniteBaffle(*s.infiniteBaffle));
    result.insert("settings", couplingSettings);
    if (a.axisymmetricSampling)
        result.insert("axisymmetricSampling", encodeAxisymmetricSampling(*a.axisymmetricSampling));
    if (s.observationFrame) {
        const auto frame = normalizedAcousticFrame(*s.observationFrame);
        auto settings = result["settings"].toObject();
        settings.insert("observationFrame",
                        QJsonObject{{"origin", encodeVector(frame.origin)},
                                    {"forward", encodeVector(frame.forward)},
                                    {"horizontal", encodeVector(frame.horizontal)}});
        result.insert("settings", settings);
    }
    return result;
}
SimulationArchive decodeSimulation(const QJsonValue &v) {
    const auto o = object(v);
    auto required = o;
    required.remove("measurements");
    required.remove("axisymmetricSampling");
    required.remove("resultComplete");
    keys(required, {"settings", "revision", "current", "solvedMesh", "request", "result"});
    SimulationArchive a;
    if (o.contains("axisymmetricSampling"))
        a.axisymmetricSampling = decodeAxisymmetricSampling(o["axisymmetricSampling"]);
    if (o.contains("measurements"))
        for (const auto &value : array(o["measurements"], 1000)) {
            const auto curve = object(value);
            keys(curve, {"name", "absoluteSpl", "offsetDb", "samples"});
            MeasuredResponse response{
                text(curve["name"]), boolean(curve["absoluteSpl"]), number(curve["offsetDb"]), {}};
            for (const auto &point : array(curve["samples"], 1000000, 2)) {
                const auto pair = array(point, 2, 2);
                response.samples.push_back({positive(pair[0]), number(pair[1])});
            }
            validateMeasuredResponse(response);
            a.measurements.push_back(std::move(response));
        }
    const bool hasMesh = !o["solvedMesh"].isNull(), hasRequest = !o["request"].isNull(),
               hasResult = !o["result"].isNull();
    require(hasMesh == hasRequest && hasRequest == hasResult);
    const auto settings = object(o["settings"]);
    auto requiredSettings = settings;
    requiredSettings.remove("ctaExcitationIndex");
    requiredSettings.remove("beamwidthTarget");
    requiredSettings.remove("normAngle");
    requiredSettings.remove("DbScale");
    requiredSettings.remove("observationFrame");
    requiredSettings.remove("acceleration");
    requiredSettings.remove("fastMultipole");
    requiredSettings.remove("nearMaximumDepth");
    requiredSettings.remove("couplingEnabled");
    requiredSettings.remove("couplingNetwork");
    requiredSettings.remove("farField");
    requiredSettings.remove("observationDistance");
    requiredSettings.remove("infiniteBaffle");
    requiredSettings.remove("angularResponses");
    requiredSettings.remove("fieldDisplay");
    requiredSettings.remove("phaseDisplay");
    requiredSettings.remove("geometryMode");
    requiredSettings.remove("axisymmetricSampling");
    keys(requiredSettings, {"minimumFrequency", "maximumFrequency", "frequencyCount", "medium",
                            "velocityScale", "phaseDegrees", "normalizeToOnAxis"});
    auto &s = a.settings;
    if (settings.contains("ctaExcitationIndex")) {
        const double selection = range(settings["ctaExcitationIndex"], -1, 10000);
        require(std::floor(selection) == selection);
        s.ctaExcitationIndex = static_cast<int>(selection);
    }
    if (settings.contains("geometryMode"))
        s.geometryMode = decodeGeometryMode(settings["geometryMode"]);
    if (settings.contains("axisymmetricSampling"))
        s.axisymmetricSampling = decodeAxisymmetricSampling(settings["axisymmetricSampling"]);
    if (settings.contains("fieldDisplay"))
        s.fieldDisplay = decodeFieldDisplay(settings["fieldDisplay"]);
    if (settings.contains("phaseDisplay"))
        s.phaseDisplay = decodePhaseDisplay(settings["phaseDisplay"]);
    if (settings.contains("angularResponses"))
        s.angularResponses = decodeAngularResponses(settings["angularResponses"]);
    if (settings.contains("infiniteBaffle"))
        s.infiniteBaffle = decodeInfiniteBaffle(settings["infiniteBaffle"]);
    s.minimumFrequency = range(settings["minimumFrequency"], 1, 100000);
    s.maximumFrequency = range(settings["maximumFrequency"], s.minimumFrequency, 100000);
    s.frequencyCount = integer(settings["frequencyCount"], std::numeric_limits<int>::max());
    require(s.frequencyCount > 0 &&
            (s.frequencyCount == 1 || s.maximumFrequency > s.minimumFrequency));
    s.medium = decodeMedium(settings["medium"]);
    require(s.medium.density >= .01 && s.medium.density <= 100 && s.medium.soundSpeed >= 1 &&
            s.medium.soundSpeed <= 10000);
    s.velocityScale = range(settings["velocityScale"], 0, 100);
    s.phaseDegrees = range(settings["phaseDegrees"], -180, 180);
    if (settings.contains("acceleration"))
        s.acceleration = decodeAcceleration(settings["acceleration"]);
    if (settings.contains("fastMultipole"))
        s.fastMultipole = decodeFastMultipole(settings["fastMultipole"]);
    if (settings.contains("nearMaximumDepth"))
        s.nearMaximumDepth = integer(settings["nearMaximumDepth"], 10);
    s.normalizeToOnAxis = boolean(settings["normalizeToOnAxis"]);
    if (settings.contains("couplingEnabled")) {
        s.couplingEnabled = boolean(settings["couplingEnabled"]);
        s.couplingNetwork = decodeLemNetwork(settings["couplingNetwork"]);
    } else
        require(!settings.contains("couplingNetwork"));
    if (settings.contains("farField")) {
        s.farField = boolean(settings["farField"]);
        s.observationDistance = range(settings["observationDistance"], .001, 100000);
    } else
        require(!settings.contains("observationDistance"));
    if (settings.contains("observationFrame")) {
        const auto frame = object(settings["observationFrame"]);
        keys(frame, {"origin", "forward", "horizontal"});
        s.observationFrame = normalizedAcousticFrame(
            {vector(frame["origin"]), vector(frame["forward"]), vector(frame["horizontal"])});
    }
    if (settings.contains("normAngle"))
        s.normalizationAngleDegrees = range(settings["normAngle"], -180, 180);
    if (settings.contains("DbScale")) {
        const auto scale = object(settings["DbScale"]);
        keys(scale, {"minimum", "maximum", "step"});
        s.directivityMinimumDb = number(scale["minimum"]);
        s.directivityMaximumDb = number(scale["maximum"]);
        s.directivityStepDb = positive(scale["step"]);
        const double span = s.directivityMaximumDb - s.directivityMinimumDb;
        const double ticks = std::ceil(span / s.directivityStepDb) + 1;
        require(std::isfinite(span) && span > 0 && std::isfinite(ticks) && ticks <= 1000);
    }
    if (settings.contains("beamwidthTarget")) {
        const auto target = object(settings["beamwidthTarget"]);
        keys(target, {"enabled", "horizontalDegrees", "verticalDegrees", "toleranceDegrees"});
        s.beamwidthTargetEnabled = boolean(target["enabled"]);
        s.targetHorizontalDegrees = range(target["horizontalDegrees"], 0, 360);
        s.targetVerticalDegrees = range(target["verticalDegrees"], 0, 360);
        s.targetToleranceDegrees = range(target["toleranceDegrees"], 0, 180);
    }
    a.revision = wideInteger(o["revision"]);
    a.current = boolean(o["current"]);
    if (o.contains("resultComplete"))
        a.resultComplete = boolean(o["resultComplete"]);
    // Historic writers only set current after the complete solve finished;
    // a stale result without this field may be complete or a partial raster.
    else if (a.current)
        a.resultComplete = true;
    require(!a.resultComplete.has_value() || hasResult);
    require(!a.current || a.resultComplete.value_or(false));
    require(!a.current || s.couplingEnabled || !hasBoundInterfacePorts(s.couplingNetwork));
    if (!o["solvedMesh"].isNull())
        a.mesh = std::make_shared<SurfaceMesh>(decodeMesh(o["solvedMesh"], true));
    if (!o["request"].isNull()) {
        require(bool(a.mesh));
        a.request = decodeRequest(o["request"], *a.mesh);
    }
    if (!o["result"].isNull()) {
        require(a.mesh && a.request.has_value());
        a.result = std::make_shared<BemResult>(decodeResult(o["result"], *a.request, *a.mesh));
    }
    require(!a.current || bool(a.result));
    require(a.axisymmetricSampling.has_value() == (a.request && a.request->axisymmetric));
    if (a.axisymmetricSampling)
        require(a.axisymmetricSampling->sectors == a.request->axisymmetric->sectors);
    if (a.current) {
        require((s.geometryMode == BemGeometryMode::axisymmetric) ==
                a.request->axisymmetric.has_value());
        if (a.axisymmetricSampling)
            require(*a.axisymmetricSampling == s.axisymmetricSampling);
    }
    require(!a.current || a.request->infiniteBaffle == s.infiniteBaffle);
    if (a.current) {
        require(a.request->coupling.has_value() == s.couplingEnabled);
        if (a.request->coupling)
            require(lem::samePhysics(*a.request->coupling, s.couplingNetwork));
    }
    require(a.result ? a.revision > 0 : a.revision == 0 && !a.current);
    return a;
}
} // namespace hornsculpt::project_json
