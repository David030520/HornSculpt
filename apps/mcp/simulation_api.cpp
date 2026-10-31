#include "simulation_api.hpp"
#include "project_file.hpp"
#include "project_json.hpp"
#include "simulation_request.hpp"
#include <QSet>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace hornsculpt::mcp {
namespace {
namespace json = project_json;
constexpr qsizetype maximumFrequencies = 8192;
constexpr qsizetype maximumObservers = 4096;
constexpr qsizetype maximumExcitations = 16;

[[noreturn]] void invalid(const QString &message) {
    throw std::invalid_argument((QStringLiteral("MCP-Studie: ") + message).toStdString());
}
void require(bool condition, const QString &message) {
    if (!condition)
        invalid(message);
}
void allowedKeys(const QJsonObject &object, std::initializer_list<const char *> allowed) {
    for (auto it = object.begin(); it != object.end(); ++it)
        if (std::none_of(allowed.begin(), allowed.end(), [&](const char *key) {
                return it.key() == QLatin1String(key);
            }))
            invalid(QStringLiteral("Unbekanntes Feld: ") + it.key());
}
double inRange(const QJsonValue &value, double minimum, double maximum) {
    const double result = json::number(value);
    require(result >= minimum && result <= maximum,
            QStringLiteral("Numerikwert außerhalb des unterstützten Bereichs."));
    return result;
}
QString identifier(const QJsonValue &value) {
    const auto result = json::text(value);
    require(!result.trimmed().isEmpty(), QStringLiteral("Leere ID oder Anregungsbezeichnung."));
    return result;
}
bool nonzero(std::complex<double> value) {
    return value.real() != 0 || value.imag() != 0;
}
QJsonArray checkedComplex(std::complex<double> value) {
    require(std::isfinite(value.real()) && std::isfinite(value.imag()),
            QStringLiteral("Nicht endlicher komplexer Ergebniswert."));
    return json::encodeComplex(value);
}
double checkedNonnegative(double value) {
    require(std::isfinite(value) && value >= 0,
            QStringLiteral("Nicht endlicher oder negativer Ergebniswert."));
    return value;
}
double count(std::size_t value) {
    require(value <= 9007199254740991ULL,
            QStringLiteral("Ergebniszähler überschreitet die genaue JSON-Zahlendarstellung."));
    return static_cast<double>(value);
}
QString accelerationName(BemAcceleration value) {
    require(value == BemAcceleration::dense || value == BemAcceleration::fastMultipole,
            QStringLiteral("Unbekannter Beschleunigungsmodus."));
    return value == BemAcceleration::dense ? QStringLiteral("dense")
                                           : QStringLiteral("fastMultipole");
}
QJsonObject encodeNumerics(const BemNumerics &n) {
    return {{"quadratureOrder", int(n.quadratureOrder)},
            {"singularOrder", int(n.singularOrder)},
            {"nearDistanceRatio", n.nearDistanceRatio},
            {"nearRelativeTolerance", n.nearRelativeTolerance},
            {"nearMaximumDepth", int(n.nearMaximumDepth)},
            {"relativeResidualTolerance", n.relativeResidualTolerance},
            {"maximumRefinementSteps", int(n.maximumRefinementSteps)},
            {"acceleration", accelerationName(n.acceleration)},
            {"fastMultipole",
             QJsonObject{{"interpolationOrder", int(n.fastMultipole.interpolationOrder)},
                         {"leafPanelCapacity", int(n.fastMultipole.leafPanelCapacity)},
                         {"separationRatio", n.fastMultipole.separationRatio},
                         {"maximumIterations", int(n.fastMultipole.maximumIterations)},
                         {"restartLength", int(n.fastMultipole.restartLength)}}}};
}
QJsonObject encodeObservation(const BemRequest &request) {
    if (request.farField)
        return {{"mode", "far_field"},
                {"originMetres", json::encodeVector(request.farField->origin)},
                {"referenceDistanceMetres", request.farField->referenceDistance},
                {"coordinateMeaning", "directions_from_origin"},
                {"pressureMeaning", "asymptotic_outgoing_pressure_at_reference_distance"}};
    return {{"mode", "finite"},
            {"originMetres", QJsonValue::Null},
            {"referenceDistanceMetres", QJsonValue::Null},
            {"coordinateMeaning", "physical_world_points"},
            {"pressureMeaning", "finite_distance_pressure"}};
}
const FaceVelocity &referenceVelocity(const std::vector<FaceVelocity> &faces,
                                      std::uint32_t faceId) {
    const auto found = std::find_if(faces.begin(), faces.end(),
                                    [=](const auto &face) { return face.faceId == faceId; });
    require(found != faces.end() && nonzero(found->normalVelocity),
            QStringLiteral("Die Bezugsfläche muss tatsächlich mit v != 0 angeregt werden."));
    return *found;
}
QJsonObject objectSchema(const QJsonObject &properties, const QJsonArray &required) {
    return {{"type", "object"}, {"properties", properties}, {"required", required},
            {"additionalProperties", false}};
}
QJsonObject arraySchema(const QJsonObject &items, int minimum, int maximum) {
    return {{"type", "array"}, {"items", items}, {"minItems", minimum}, {"maxItems", maximum}};
}
QJsonObject inputSchema() {
    const QJsonObject number{{"type", "number"}};
    const QJsonObject positive{{"type", "number"}, {"exclusiveMinimum", 0}};
    const QJsonObject text{{"type", "string"}, {"minLength", 1}, {"maxLength", 4096}};
    const QJsonObject faceId{{"type", "integer"}, {"minimum", 1}, {"maximum", 4294967295.0}};
    auto frequencies = arraySchema(positive, 1, int(maximumFrequencies));
    frequencies.insert("description", "Explicit strictly ascending positive frequencies in Hz.");
    auto observers = arraySchema(
        objectSchema({{"id", text},
                      {"role", QJsonObject{{"type", "string"},
                                           {"enum", QJsonArray{"reference", "validation", "sample"}}}},
                      {"positionMetres", arraySchema(number, 3, 3)}},
                     {"id", "role", "positionMetres"}),
        1, int(maximumObservers));
    observers.insert("description", "Unique IDs; exactly one reference role. World SI coordinates.");
    const auto excitation = objectSchema(
        {{"name", text},
         {"faces", arraySchema(objectSchema({{"faceId", faceId},
                                              {"normalVelocity", arraySchema(number, 2, 2)}},
                                             {"faceId", "normalVelocity"}),
                                1, 65536)}},
        {"name", "faces"});
    const auto finite = objectSchema({{"mode", QJsonObject{{"const", "finite"}}}}, {"mode"});
    const auto far = objectSchema(
        {{"mode", QJsonObject{{"const", "far_field"}}},
         {"originMetres", arraySchema(number, 3, 3)}, {"referenceDistanceMetres", positive}},
        {"mode"});
    return objectSchema(
        {{"frequenciesHz", frequencies}, {"observers", observers},
         {"observation", QJsonObject{{"oneOf", QJsonArray{finite, far}},
                                      {"description", "Default: finite. Far-field defaults: origin [0,0,0], distance 1 m."}}},
         {"excitations", arraySchema(excitation, 1, int(maximumExcitations))},
         {"reference", objectSchema({{"excitation", text}, {"faceId", faceId}},
                                     {"excitation", "faceId"})},
         {"medium", objectSchema({{"densityKgM3", positive},
                                   {"soundSpeedMetresPerSecond", positive}},
                                  {"densityKgM3", "soundSpeedMetresPerSecond"})},
         {"numerics", objectSchema(
                          {{"acceleration", QJsonObject{{"enum", QJsonArray{"dense", "fastMultipole"}}}},
                           {"nearRelativeTolerance", QJsonObject{{"type", "number"}, {"minimum", 1e-10}, {"maximum", .01}}},
                           {"nearMaximumDepth", QJsonObject{{"type", "integer"}, {"minimum", 0}, {"maximum", 10}}},
                           {"relativeResidualTolerance", QJsonObject{{"type", "number"}, {"minimum", 1e-7}, {"maximum", .01}}}},
                          {})}},
        {"frequenciesHz", "observers", "excitations", "reference"});
}
} // namespace

StudyRequest parseStudyRequest(const QJsonObject &input, const Project &project) {
    allowedKeys(input, {"frequenciesHz", "observers", "observation", "excitations", "reference",
                        "medium", "numerics"});
    const auto &settings = project.simulation.settings;
    require(!settings.couplingEnabled && !hasBoundInterfacePorts(settings.couplingNetwork),
            QStringLiteral("Aktive LEM-Kopplung oder zugeordnete Interfaceports werden vom "
                           "Studienadapter noch nicht unterstützt."));
    require(!settings.infiniteBaffle,
            QStringLiteral("Die aktive unendliche Schallwand wird vom Studienadapter noch "
                           "nicht unterstützt."));
    require(settings.geometryMode == BemGeometryMode::threeDimensional,
            QStringLiteral("Der Studienadapter unterstützt ausschließlich vollständige 3D-BEM."));
    require(project.mesh && project.meshRevision == project.revision,
            QStringLiteral("Ein aktuelles geprüftes Netz ist erforderlich."));
    const auto meshFaces = json::meshFaces(project.mesh->mesh);
    StudyRequest result;
    auto &request = result.bem;
    request.medium = settings.medium;
    request.numerics.acceleration = settings.acceleration;
    request.numerics.fastMultipole = settings.fastMultipole;
    request.numerics.nearMaximumDepth = settings.nearMaximumDepth;

    for (const auto &value : json::array(input["frequenciesHz"], maximumFrequencies, 1)) {
        const double frequency = json::positive(value);
        require(request.frequencies.empty() || frequency > request.frequencies.back(),
                QStringLiteral("Frequenzen müssen streng aufsteigend sein."));
        request.frequencies.push_back(frequency);
    }
    request.interactionPlanFrequencies = request.frequencies;
    QSet<QString> observerIds;
    unsigned referenceObservers = 0;
    for (const auto &value : json::array(input["observers"], maximumObservers, 1)) {
        const auto observer = json::object(value);
        json::keys(observer, {"id", "role", "positionMetres"});
        const auto id = identifier(observer["id"]), role = json::text(observer["role"]);
        require(!observerIds.contains(id), QStringLiteral("Doppelte Beobachter-ID: ") + id);
        observerIds.insert(id);
        require(role == "reference" || role == "validation" || role == "sample",
                QStringLiteral("Unbekannte Beobachterrolle."));
        referenceObservers += role == "reference" ? 1U : 0U;
        request.observers.push_back(json::vector(observer["positionMetres"]));
        result.observerMetadata.append(observer);
    }
    require(referenceObservers == 1,
            QStringLiteral("Genau ein Beobachter benötigt die Rolle reference."));
    if (input.contains("observation")) {
        const auto observation = json::object(input["observation"]);
        const auto mode = json::text(observation["mode"]);
        if (mode == "finite")
            json::keys(observation, {"mode"});
        else if (mode == "far_field") {
            allowedKeys(observation, {"mode", "originMetres", "referenceDistanceMetres"});
            FarFieldObservation far;
            if (observation.contains("originMetres"))
                far.origin = json::vector(observation["originMetres"]);
            if (observation.contains("referenceDistanceMetres"))
                far.referenceDistance = json::positive(observation["referenceDistanceMetres"]);
            for (const auto &point : request.observers) {
                const double distance = std::hypot(point.x - far.origin.x, point.y - far.origin.y,
                                                   point.z - far.origin.z);
                require(std::isfinite(distance) && distance > 0,
                        QStringLiteral("Fernfeldrichtung muss endlich und vom Ursprung verschieden sein."));
            }
            request.farField = far;
        } else
            invalid(QStringLiteral("Beobachtungsmodus muss finite oder far_field sein."));
    }
    QSet<QString> excitationNames;
    for (const auto &value : json::array(input["excitations"], maximumExcitations, 1)) {
        const auto inputExcitation = json::object(value);
        json::keys(inputExcitation, {"name", "faces"});
        const auto name = identifier(inputExcitation["name"]);
        require(!excitationNames.contains(name), QStringLiteral("Doppelte Anregung: ") + name);
        excitationNames.insert(name);
        BemExcitation excitation;
        excitation.name = name.toStdString();
        std::set<std::uint32_t> faces;
        bool hasVelocity = false;
        for (const auto &faceValue : json::array(inputExcitation["faces"], 65536, 1)) {
            const auto face = json::object(faceValue);
            json::keys(face, {"faceId", "normalVelocity"});
            const auto id = json::integer(face["faceId"]);
            require(id && meshFaces.contains(id) &&
                        meshFaces.at(id) == BoundaryCondition::driven && faces.insert(id).second,
                    QStringLiteral("Jede Anregung darf vorhandene Driven-Flächen höchstens einmal zuweisen."));
            const auto velocity = json::complex(face["normalVelocity"]);
            hasVelocity = hasVelocity || nonzero(velocity);
            excitation.faces.push_back({id, velocity});
        }
        require(hasVelocity, QStringLiteral("Jede Anregung benötigt mindestens eine Geschwindigkeit != 0."));
        request.excitations.push_back(std::move(excitation));
    }
    const auto reference = json::object(input["reference"]);
    json::keys(reference, {"excitation", "faceId"});
    result.referenceExcitation = identifier(reference["excitation"]);
    result.referenceFaceId = json::integer(reference["faceId"]);
    const auto excitation = std::find_if(request.excitations.begin(), request.excitations.end(),
                                          [&](const auto &e) {
                                              return e.name == result.referenceExcitation.toStdString();
                                          });
    require(excitation != request.excitations.end(), QStringLiteral("Bezugsanregung fehlt."));
    (void)referenceVelocity(excitation->faces, result.referenceFaceId);
    if (input.contains("medium")) {
        const auto medium = json::object(input["medium"]);
        json::keys(medium, {"densityKgM3", "soundSpeedMetresPerSecond"});
        request.medium = {json::positive(medium["densityKgM3"]),
                          json::positive(medium["soundSpeedMetresPerSecond"])};
    }
    require(std::isfinite(request.medium.density) && request.medium.density > 0 &&
                std::isfinite(request.medium.soundSpeed) && request.medium.soundSpeed > 0,
            QStringLiteral("Das Medium benötigt positive endliche Dichte und Schallgeschwindigkeit."));
    if (input.contains("numerics")) {
        const auto numerics = json::object(input["numerics"]);
        allowedKeys(numerics, {"acceleration", "nearRelativeTolerance", "nearMaximumDepth",
                               "relativeResidualTolerance"});
        auto &n = request.numerics;
        if (numerics.contains("acceleration")) {
            const auto mode = json::text(numerics["acceleration"]);
            require(mode == "dense" || mode == "fastMultipole",
                    QStringLiteral("Unbekannter Beschleunigungsmodus."));
            n.acceleration = mode == "dense" ? BemAcceleration::dense : BemAcceleration::fastMultipole;
        }
        // These are the public CUDA Galerkin bounds; the solver remains the final
        // validator of all inherited numerical settings and physical geometry.
        if (numerics.contains("nearRelativeTolerance"))
            n.nearRelativeTolerance = inRange(numerics["nearRelativeTolerance"], 1e-10, .01);
        if (numerics.contains("nearMaximumDepth"))
            n.nearMaximumDepth = json::integer(numerics["nearMaximumDepth"], 10);
        if (numerics.contains("relativeResidualTolerance"))
            n.relativeResidualTolerance = inRange(numerics["relativeResidualTolerance"], 1e-7, .01);
    }
    (void)accelerationName(request.numerics.acceleration);
    require(request.numerics.nearMaximumDepth <= 10,
            QStringLiteral("Maximale Nahquadraturtiefe muss zwischen 0 und 10 liegen."));
    return result;
}

QJsonObject describeStudyCapabilities() {
    return {{"schemaVersion", 1},
            {"inputSchema", inputSchema()},
            {"geometryMode", "threeDimensional"},
            {"sourceModel", "prescribed_normal_velocity"},
            {"observations", QJsonArray{"finite", "far_field"}},
            {"acceleration", QJsonArray{"dense", "fastMultipole"}},
            {"maximumFrequencies", int(maximumFrequencies)},
            {"maximumObservers", int(maximumObservers)},
            {"maximumExcitations", int(maximumExcitations)},
            {"unsupportedActiveProjectPhysics", QJsonArray{"LEM", "infinite_baffle", "axisymmetric_m0"}},
            {"finiteObserverValidation", "CUDA rejects solid and boundary points; near-field quadrature must converge."},
            {"excitationConvention", "Complex RMS m/s into the BEM air domain; omitted Driven faces have zero velocity. Explicit study velocities replace saved GUI amplitude/phase."},
            {"transferDefinition", "p / actual reference-face normal velocity, in Pa/(m/s), for every observer of the selected excitation."},
            {"modalDecomposition", false},
            {"impulseResponse", false},
            {"firDesign", false},
            {"audioRequirements", "Complex frequency samples are not an audio IR. An IR requires a specified uniform DFT raster, sample rate, time window, DC/Nyquist and out-of-band treatment. Log/sparse rasters cannot be treated as FFT bins. A velocity-driven acoustic H is not an electrical loudspeaker response."}};
}

QJsonObject complexResultDocument(const StudyRequest &request, const BemResult &result) {
    const auto &bem = request.bem;
    require(!bem.coupling && !bem.infiniteBaffle && !bem.axisymmetric &&
                result.geometryMode == BemGeometryMode::threeDimensional,
            QStringLiteral("Ergebnisphysik wird vom Studienadapter nicht unterstützt."));
    require(result.frequencies.size() == bem.frequencies.size() && !result.frequencies.empty() &&
                request.observerMetadata.size() == static_cast<qsizetype>(bem.observers.size()) &&
                result.acceleration == bem.numerics.acceleration && !result.device.empty() &&
                result.triangleCount && result.pressureDofCount,
            QStringLiteral("Ergebnis und Studienauftrag passen nicht zusammen."));
    QJsonArray frequencies, requestedExcitations;
    for (const auto &excitation : bem.excitations) {
        QJsonArray faces;
        for (const auto &face : excitation.faces)
            faces.append(QJsonObject{{"faceId", double(face.faceId)},
                                      {"normalVelocity", checkedComplex(face.normalVelocity)}});
        requestedExcitations.append(QJsonObject{{"name", QString::fromStdString(excitation.name)},
                                                 {"faces", faces}});
    }
    QString referenceObserverId;
    for (const auto &metadata : request.observerMetadata)
        if (metadata.toObject()["role"] == "reference")
            referenceObserverId = metadata.toObject()["id"].toString();
    require(!referenceObserverId.isEmpty(), QStringLiteral("Referenzbeobachter fehlt."));
    for (std::size_t f = 0; f < result.frequencies.size(); ++f) {
        const auto &frame = result.frequencies[f];
        require(std::isfinite(frame.frequency) && frame.frequency == bem.frequencies[f] &&
                    frame.excitations.size() == bem.excitations.size(),
                QStringLiteral("Frequenz- oder Anregungszuordnung des Ergebnisses stimmt nicht."));
        QJsonArray excitations;
        QJsonObject selectedTransfer;
        for (std::size_t e = 0; e < frame.excitations.size(); ++e) {
            const auto &excitation = frame.excitations[e];
            const auto &input = bem.excitations[e];
            require(excitation.name == input.name && !excitation.lemSolution &&
                        excitation.observerPressure.size() == bem.observers.size() &&
                        excitation.solvedVelocities.size() == input.faces.size(),
                    QStringLiteral("Komplexe Beobachter- oder Quellendaten sind unvollständig."));
            std::map<std::uint32_t, std::complex<double>> actualVelocities;
            QJsonArray velocities, pressures, impedances;
            for (const auto &face : excitation.solvedVelocities) {
                require(actualVelocities.emplace(face.faceId, face.normalVelocity).second,
                        QStringLiteral("Doppelte tatsächliche Quellfläche."));
                velocities.append(QJsonObject{{"faceId", double(face.faceId)},
                                               {"normalVelocity", checkedComplex(face.normalVelocity)}});
            }
            for (const auto &face : input.faces)
                require(actualVelocities.contains(face.faceId) &&
                            actualVelocities.at(face.faceId) == face.normalVelocity,
                        QStringLiteral("Tatsächliche Quellgeschwindigkeit stimmt nicht mit dem Auftrag überein."));
            for (const auto pressure : excitation.observerPressure)
                pressures.append(checkedComplex(pressure));
            std::set<std::uint32_t> impedanceFaces;
            for (const auto &impedance : excitation.impedances) {
                require(actualVelocities.contains(impedance.faceId) &&
                            nonzero(actualVelocities.at(impedance.faceId)) &&
                            impedanceFaces.insert(impedance.faceId).second &&
                            std::isfinite(impedance.area) && impedance.area > 0 &&
                            nonzero(impedance.volumeVelocity),
                        QStringLiteral("Ungültige Quellimpedanz oder Volumengeschwindigkeit."));
                impedances.append(QJsonObject{
                    {"faceId", double(impedance.faceId)}, {"areaSquareMetres", impedance.area},
                    {"meanPressurePaRms", checkedComplex(impedance.meanPressure)},
                    {"volumeVelocityCubicMetresPerSecondRms", checkedComplex(impedance.volumeVelocity)},
                    {"volumeImpedancePascalSecondsPerCubicMetre", checkedComplex(impedance.volumeImpedance)}});
            }
            for (const auto &face : input.faces)
                require(!nonzero(face.normalVelocity) || impedanceFaces.contains(face.faceId),
                        QStringLiteral("Volumenimpedanz einer angeregten Fläche fehlt."));
            const double residual = checkedNonnegative(excitation.boundaryResidual);
            require(residual <= bem.numerics.relativeResidualTolerance,
                    QStringLiteral("Originalgleichungsresiduum überschreitet die angeforderte Toleranz."));
            excitations.append(QJsonObject{
                {"name", QString::fromStdString(excitation.name)}, {"pressurePaRms", pressures},
                {"actualFaceVelocitiesMetresPerSecondRms", velocities},
                {"sourceImpedances", impedances}, {"relativeBoundaryResidual", residual},
                {"refinementSteps", int(excitation.refinementSteps)},
                {"iterations", int(excitation.iterations)}});
            if (excitation.name == request.referenceExcitation.toStdString()) {
                const auto velocity = referenceVelocity(excitation.solvedVelocities,
                                                         request.referenceFaceId).normalVelocity;
                QJsonArray transfer;
                for (const auto pressure : excitation.observerPressure)
                    transfer.append(checkedComplex(pressure / velocity));
                selectedTransfer = {{"excitation", request.referenceExcitation},
                                    {"referenceFaceId", double(request.referenceFaceId)},
                                    {"referenceVelocity", checkedComplex(velocity)},
                                    {"unit", "Pa/(m/s)"},
                                    {"pressureOverReferenceVelocity", transfer}};
            }
        }
        require(!selectedTransfer.isEmpty(), QStringLiteral("Komplexe Bezugsübertragung fehlt."));
        frequencies.append(QJsonObject{
            {"frequencyHz", frame.frequency}, {"excitations", excitations},
            {"selectedTransfer", selectedTransfer},
            {"timingMilliseconds", QJsonObject{
                {"assembly", checkedNonnegative(frame.assemblyMilliseconds)},
                {"solve", checkedNonnegative(frame.solveMilliseconds)},
                {"factorization", checkedNonnegative(frame.factorizationMilliseconds)},
                {"triangularSolve", checkedNonnegative(frame.triangularSolveMilliseconds)},
                {"field", checkedNonnegative(frame.fieldMilliseconds)},
                {"residualCheck", checkedNonnegative(frame.residualCheckMilliseconds)},
                {"interactionPlan", checkedNonnegative(frame.interactionPlanMilliseconds)}}},
            {"interactionPlanReused", frame.interactionPlanReused},
            {"fastMultipole", QJsonObject{
                {"nodeCount", count(frame.fastMultipole.nodeCount)},
                {"leafCount", count(frame.fastMultipole.leafCount)},
                {"farInteractionCount", count(frame.fastMultipole.farInteractionCount)},
                {"directPanelPairCount", count(frame.fastMultipole.directPanelPairCount)},
                {"representedFarPanelPairs", count(frame.fastMultipole.representedFarPanelPairs)}}}});
    }
    return {{"schema", "hornsculpt.complex-study"}, {"schemaVersion", 1},
            {"conventions", QJsonObject{{"time", "exp(+i*omega*t)"}, {"outgoingWave", "exp(-i*k*r)"},
                                         {"amplitude", "RMS"}, {"complexEncoding", "[real,imaginary]"},
                                         {"lengthUnit", "m"}, {"frequencyUnit", "Hz"},
                                         {"pressureUnit", "Pa"},
                                         {"positiveNormalVelocity", "into_the_BEM_air_domain"}}},
            {"medium", QJsonObject{{"densityKgM3", bem.medium.density},
                                    {"soundSpeedMetresPerSecond", bem.medium.soundSpeed}}},
            {"observation", encodeObservation(bem)}, {"observers", request.observerMetadata},
            {"reference", QJsonObject{{"excitation", request.referenceExcitation},
                                       {"faceId", double(request.referenceFaceId)},
                                       {"observerId", referenceObserverId}}},
            {"requestedExcitations", requestedExcitations},
            {"numerics", encodeNumerics(bem.numerics)},
            {"solver", QJsonObject{{"device", QString::fromStdString(result.device)},
                                    {"triangleCount", count(result.triangleCount)},
                                    {"pressureDofCount", count(result.pressureDofCount)},
                                    {"peakDeviceBytes", count(result.peakDeviceBytes)},
                                    {"acceleration", accelerationName(result.acceleration)},
                                    {"geometryMode", "threeDimensional"},
                                    {"method", "Burton-Miller Galerkin P1/DP0"}}},
            {"frequencies", frequencies},
            {"audioInterpretation", "Finite positive-frequency acoustic samples; no IR, FIR or modal decomposition has been computed."}};
}
} // namespace hornsculpt::mcp
