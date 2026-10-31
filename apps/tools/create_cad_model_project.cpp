#include "cad_snapshot.hpp"
#include "project_file.hpp"
#include "result_export.hpp"
#include "simulation_request.hpp"
#include "step_import.hpp"

#include <QCommandLineParser>
#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <hornsculpt/axisymmetric_geometry.hpp>
#include <hornsculpt/axisymmetric_mesh.hpp>
#include <hornsculpt/acoustic_analysis.hpp>
#include <hornsculpt/cad_geometry.hpp>
#include <hornsculpt/field_evaluation.hpp>
#include <hornsculpt/surface_measures.hpp>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

using namespace hornsculpt;
namespace {
Vec3 vector(const QJsonValue &value) {
    const auto a = value.toArray();
    if (a.size() != 3)
        throw std::invalid_argument("Ein Vektor benötigt drei SI-Komponenten.");
    for (const auto &v : a)
        if (!v.isDouble() || !std::isfinite(v.toDouble()))
            throw std::invalid_argument("Nicht endlicher Modellvektor.");
    return {a[0].toDouble(), a[1].toDouble(), a[2].toDouble()};
}

double distance(Vec3 a, Vec3 b) {
    return std::hypot(a.x - b.x, a.y - b.y, a.z - b.z);
}

std::uint32_t identifier(const QJsonValue &value) {
    const auto id = value.toDouble(-1);
    if (!value.isDouble() || !std::isfinite(id) || id < 1 || id > UINT32_MAX ||
        std::floor(id) != id)
        throw std::invalid_argument("Flächenidentität muss eine positive Ganzzahl sein.");
    return static_cast<std::uint32_t>(id);
}

FaceColor faceColor(const QJsonValue &value) {
    const auto rgb = value.toArray();
    if (rgb.size() != 3)
        throw std::invalid_argument("Flächenfarbe benötigt drei RGB-Kanäle.");
    const auto channel = [](const QJsonValue &v) {
        const auto c = v.toDouble(-1);
        if (!v.isDouble() || !std::isfinite(c) || c < 0 || c > 255 || std::floor(c) != c)
            throw std::invalid_argument("RGB-Kanäle müssen Ganzzahlen zwischen 0 und 255 sein.");
        return static_cast<std::uint8_t>(c);
    };
    return {channel(rgb[0]), channel(rgb[1]), channel(rgb[2])};
}

void saveChecked(const QString &path, const Project &project) {
    writeProject(path, project);
    const auto restored = readProject(path);
    if (restored.simulation.settings.infiniteBaffle != project.simulation.settings.infiniteBaffle ||
        !restored.cad || restored.cad->faces.size() != project.cad->faces.size() ||
        bool(restored.mesh) != bool(project.mesh) ||
        bool(restored.simulation.result) != bool(project.simulation.result) ||
        restored.cadFaceColors != project.cadFaceColors ||
        restored.simulation.settings.couplingEnabled !=
            project.simulation.settings.couplingEnabled ||
        !lem::samePhysics(restored.simulation.settings.couplingNetwork,
                          project.simulation.settings.couplingNetwork))
        throw std::runtime_error("CAD-Projekt-Roundtrip fehlgeschlagen.");
    for (const auto &face : project.cad->faces) {
        const auto found = std::find_if(restored.cad->faces.begin(), restored.cad->faces.end(),
                                        [&](const auto &f) { return f.id == face.id; });
        if (found == restored.cad->faces.end() || found->boundary != face.boundary ||
            found->normalVelocity != face.normalVelocity)
            throw std::runtime_error("CAD-Projekt-Roundtrip verändert Flächenanregungen.");
    }
    std::cout << "Projekt gespeichert und geprüft: " << path.toStdString() << '\n' << std::flush;
}

void reuseMesh(Project &project, const QString &path) {
    const auto previous = readProject(path);
    if (!previous.mesh || previous.meshRevision != previous.geometryRevision || !previous.cad ||
        !project.cad || previous.cad->documentId != project.cad->documentId ||
        previous.cad->brepMetresPerUnit != project.cad->brepMetresPerUnit ||
        previous.cad->brepObjects.size() != project.cad->brepObjects.size() ||
        previous.cad->objects.size() != project.cad->objects.size() ||
        previous.cad->faces.size() != project.cad->faces.size() ||
        previous.cadReversedFaceIds != project.cadReversedFaceIds ||
        previous.cadFaceTargetEdgeLengths != project.cadFaceTargetEdgeLengths)
        throw std::invalid_argument(
            "Netzwiederverwendung benötigt dasselbe bestätigte CAD-Modell.");
    const auto &a = previous.meshing;
    const auto &b = project.meshing;
    if (a.targetEdgeLength != b.targetEdgeLength || a.adaptive.enabled != b.adaptive.enabled ||
        a.adaptive.minimumEdgeLength != b.adaptive.minimumEdgeLength ||
        a.adaptive.geometricTolerance != b.adaptive.geometricTolerance ||
        a.adaptive.growthRate != b.adaptive.growthRate || a.quality.enabled != b.quality.enabled ||
        a.quality.minimumAngleDegrees != b.quality.minimumAngleDegrees ||
        a.quality.maximumAspectRatio != b.quality.maximumAspectRatio)
        throw std::invalid_argument("Netzwiederverwendung: Netzeinstellungen unterscheiden sich.");
    for (std::size_t i = 0; i < previous.cad->brepObjects.size(); ++i) {
        const auto &oldPart = previous.cad->brepObjects[i];
        const auto &newPart = project.cad->brepObjects[i];
        if (oldPart.objectId != newPart.objectId || oldPart.data != newPart.data ||
            oldPart.faceIds != newPart.faceIds)
            throw std::invalid_argument(
                "Netzwiederverwendung: Exaktes BRep oder Flächen-IDs geändert.");
    }
    for (std::size_t i = 0; i < previous.cad->objects.size(); ++i)
        if (previous.cad->objects[i].id != project.cad->objects[i].id ||
            previous.cad->objects[i].placement != project.cad->objects[i].placement)
            throw std::invalid_argument("Netzwiederverwendung: CAD-Platzierung geändert.");
    for (std::size_t i = 0; i < previous.cad->faces.size(); ++i)
        if (previous.cad->faces[i].id != project.cad->faces[i].id ||
            previous.cad->faces[i].boundary != project.cad->faces[i].boundary ||
            previous.cad->faces[i].state != CadFaceState::confirmed)
            throw std::invalid_argument(
                "Netzwiederverwendung: Randbedingung oder Freigabe geändert.");
    project.mesh = previous.mesh;
    project.meshRevision = project.geometryRevision;
    std::cout << "Vorhandenes geprüftes Netz übernommen: " << project.mesh->mesh.triangles.size()
              << " Dreiecke.\n";
}

void alignReferenceModeArea(Project &project, const SurfaceMesh *integrationMesh = nullptr,
                            const std::optional<AxisymmetricBem> &axisymmetric = {}) {
    if (!project.mesh && !integrationMesh)
        return;
    const auto &mesh = integrationMesh ? *integrationMesh : project.mesh->mesh;
    const auto measures = surfaceMeasures(mesh, axisymmetric);
    auto &network = project.simulation.settings.couplingNetwork;
    for (auto &element : network.elements) {
        auto *driver = std::get_if<lem::Driver>(&element.component);
        if (!driver)
            continue;
        for (const auto *candidate : lem::bemPorts(network)) {
            const auto &port = std::get<lem::BemPort>(candidate->component);
            if (port.kind != lem::BemPortKind::drivenSurface ||
                port.terminals.positive != driver->front)
                continue;
            double modeArea = 0;
            for (std::size_t ti = 0; ti < mesh.triangles.size(); ++ti) {
                const auto &triangle = mesh.triangles[ti];
                if (std::find(port.faceIds.begin(), port.faceIds.end(), triangle.faceId) ==
                    port.faceIds.end())
                    continue;
                modeArea += measures.areas[ti] * (port.velocityWeights.contains(triangle.faceId)
                                                      ? port.velocityWeights.at(triangle.faceId)
                                                      : 1.);
            }
            // Match the actual power-conjugate integration surface, including
            // the continuous ring measure in the meridional formulation.
            driver->parameters.pistonArea = modeArea;
            std::cout << "Diskrete Modusfläche des Motors: " << modeArea * 1e4 << " cm².\n";
        }
    }
}

CadSnapshot nativeSnapshot(const QString &path, const CadSnapshot *previous = nullptr) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly) || file.size() > 64 * 1024 * 1024)
        throw std::invalid_argument("Nativer FreeCAD-Snapshot nicht lesbar oder zu groß.");
    const auto bytes = file.readAll();
    const auto root = QJsonDocument::fromJson(bytes).object();
    auto snapshot = cad_detail::decodeSnapshot(bytes, root["sessionId"].toString());
    std::uint32_t nextId = 1;
    if (previous)
        for (const auto &face : previous->faces)
            nextId = std::max(nextId, face.id + 1);
    cad_detail::reconcileSnapshot(snapshot, previous, nextId);
    prepareCadGeometry(snapshot, true);
    return snapshot;
}

Project prepare(const QString &manifestPath, const QString &snapshotPath,
                const QString &driverReferencePath, const QString &interiorReferencePath) {
    QFile file(manifestPath);
    if (!file.open(QIODevice::ReadOnly) || file.size() > 1024 * 1024)
        throw std::runtime_error("Modellmanifest nicht lesbar oder zu groß.");
    QJsonParseError parseError;
    const auto document = QJsonDocument::fromJson(file.readAll(), &parseError);
    if (parseError.error != QJsonParseError::NoError || !document.isObject())
        throw std::invalid_argument("Ungültiges Modellmanifest.");
    const auto root = document.object();
    const auto stepName = root["stepFile"].toString();
    if (stepName.isEmpty())
        throw std::invalid_argument("STEP-Datei fehlt im Modellmanifest.");
    const auto stepPath = QFileInfo(manifestPath).dir().absoluteFilePath(stepName);
    auto cad = snapshotPath.isEmpty()
                   ? std::make_shared<CadSnapshot>(*readStepSnapshot(stepPath, .0001))
                   : std::make_shared<CadSnapshot>(nativeSnapshot(snapshotPath));
    cad->documentName = root["name"].toString(QFileInfo(stepPath).completeBaseName());
    std::optional<InfiniteBaffle> baffle;
    if (root.contains("aperture")) {
        if (!root["aperture"].isObject())
            throw std::invalid_argument("Schallwandapertur muss ein Objekt sein.");
        const auto aperture = root["aperture"].toObject();
        baffle = InfiniteBaffle{};
        baffle->origin = vector(aperture["origin_m"]);
        baffle->normal = vector(aperture["normal"]);
        cad->acousticOrigin = baffle->origin;
        cad->acousticAxis = baffle->normal;
    } else {
        const auto frame = root["acousticFrame"].toObject();
        cad->acousticOrigin = vector(frame["origin_m"]);
        cad->acousticAxis = vector(frame["axis"]);
    }
    std::set<std::uint32_t> matched;
    std::map<std::uint32_t, std::uint32_t> sourceFaces;
    std::map<QString, std::vector<std::uint32_t>> semanticGroups;
    std::set<std::uint32_t> interfaceFaces;
    const auto descriptions = root["faces"].toArray();
    if (descriptions.size() != static_cast<qsizetype>(cad->faces.size()))
        throw std::invalid_argument("Manifest enthält nicht sämtliche STEP-Flächen.");
    Project project;
    project.source = ProjectSource::freecad;
    project.geometryRevision = project.revision;
    project.meshing.targetEdgeLength = .006;
    project.meshing.adaptive = {true, .0006, .0001, .3};
    if (root.contains("meshing")) {
        if (!root["meshing"].isObject())
            throw std::invalid_argument("Netzvorgaben im Manifest müssen ein Objekt sein.");
        const auto settings = root["meshing"].toObject();
        const auto length = [&](const char *key, double fallback) {
            if (!settings.contains(key))
                return fallback;
            const auto value = settings[key];
            if (!value.isDouble() || !std::isfinite(value.toDouble()) || value.toDouble() <= 0)
                throw std::invalid_argument(
                    "Netzlängen im Manifest müssen positiv und endlich sein.");
            return value.toDouble();
        };
        project.meshing.targetEdgeLength = length("targetEdge_m", project.meshing.targetEdgeLength);
        project.meshing.adaptive.minimumEdgeLength =
            length("minimumEdge_m", project.meshing.adaptive.minimumEdgeLength);
        project.meshing.adaptive.geometricTolerance =
            length("geometricTolerance_m", project.meshing.adaptive.geometricTolerance);
    }
    validateMeshSettings(project.meshing);
    for (const auto &entry : descriptions) {
        const auto description = entry.toObject();
        const auto area = description["area_m2"].toDouble(-1);
        const auto sample = vector(description["samplePoint_m"]);
        const auto role = description["role"].toString();
        if (!std::isfinite(area) || area <= 0 ||
            (role != "driven" && role != "wall" && role != "aperture" && role != "interface"))
            throw std::invalid_argument("Ungültige Flächenrolle oder -fläche.");
        std::vector<std::size_t> candidates;
        // A STEP translator may reorder faces. A physical interior sample is
        // independent of face numbering and area-quadrature changes on SI scaling.
        // 0.1 micrometre is an export identity budget, not a scan/model accuracy.
        constexpr double exportIdentityTolerance = 1e-7;
        for (std::size_t i = 0; i < cad->faces.size(); ++i) {
            const auto &f = cad->faces[i];
            if (!matched.contains(f.id) &&
                distance(cad->geometry->projectSurface(f.id, sample).position, sample) <=
                    exportIdentityTolerance)
                candidates.push_back(i);
        }
        if (candidates.size() != 1) {
            throw std::runtime_error("STEP-Flächensignatur ist fehlend oder mehrdeutig: " +
                                     description["label"].toString().toStdString());
        }
        auto &face = cad->faces[candidates.front()];
        matched.insert(face.id);
        face.label = description["label"].toString();
        if (const auto group = description["group"].toString(); !group.isEmpty())
            semanticGroups[group].push_back(face.id);
        face.state = CadFaceState::confirmed;
        face.boundary = role == "driven" ? BoundaryCondition::driven : BoundaryCondition::wall;
        face.normalVelocity = {1, 0};
        if (description.contains("normalVelocity_m_s")) {
            const auto amplitude = description["normalVelocity_m_s"].toArray();
            if (role != "driven" || amplitude.size() != 2 || !amplitude[0].isDouble() ||
                !amplitude[1].isDouble() || !std::isfinite(amplitude[0].toDouble()) ||
                !std::isfinite(amplitude[1].toDouble()))
                throw std::invalid_argument(
                    "Nur Driven-Flächen erhalten komplexe SI-Geschwindigkeit.");
            face.normalVelocity = {amplitude[0].toDouble(), amplitude[1].toDouble()};
        }
        if (description.contains("color_rgb"))
            project.cadFaceColors.emplace(face.id, faceColor(description["color_rgb"]));
        if (description.contains("sourceFace")) {
            const auto source = identifier(description["sourceFace"]);
            if (!sourceFaces.emplace(source, face.id).second)
                throw std::invalid_argument("Doppelte Quellflächenidentität im Modellmanifest.");
        }
        if (role == "interface")
            interfaceFaces.insert(face.id);
        if (role == "aperture") {
            if (!baffle)
                throw std::invalid_argument("Aperturfläche benötigt eine explizite Schallwand.");
            baffle->apertureFaceIds.push_back(face.id);
        }
        const auto localEdge = description["targetEdge_m"].toDouble(0);
        if (localEdge != 0) {
            if (!std::isfinite(localEdge) ||
                localEdge < project.meshing.adaptive.minimumEdgeLength ||
                localEdge > project.meshing.targetEdgeLength)
                throw std::invalid_argument("Lokale Zielkante außerhalb der Netzgrenzen.");
            project.cadFaceTargetEdgeLengths.emplace(face.id, localEdge);
        }
        for (auto &triangle : cad->mesh.triangles)
            if (triangle.faceId == face.id)
                triangle.boundary = face.boundary;
    }
    if (baffle && baffle->apertureFaceIds.empty())
        throw std::invalid_argument("Die Rezessgeometrie benötigt eine Strahlungsapertur.");
    project.cad = cad;
    project.preview = std::make_shared<SurfaceMesh>(cad->mesh);
    project.simulation.settings.infiniteBaffle = baffle;
    project.simulation.settings.minimumFrequency = 1000;
    project.simulation.settings.maximumFrequency = 10000;
    project.simulation.settings.frequencyCount = 7;
    project.simulation.settings.observationFrame =
        defaultAcousticFrame(cad->acousticOrigin, cad->acousticAxis);
    lem::Network network;
    network.nodes.push_back({{1}, lem::Domain::acoustic, true, "Druckreferenz", {0, 180}});
    auto interfaceGroups = root["lemInterfaces"];
    if (!interiorReferencePath.isEmpty()) {
        QFile interiorFile(interiorReferencePath);
        if (!interiorFile.open(QIODevice::ReadOnly) || interiorFile.size() > 1024 * 1024)
            throw std::invalid_argument("LEM-Innenraummodell ist nicht lesbar oder zu groß.");
        QJsonParseError interiorError;
        const auto interior = QJsonDocument::fromJson(interiorFile.readAll(), &interiorError);
        if (interiorError.error != QJsonParseError::NoError || !interior.isObject())
            throw std::invalid_argument("Ungültiges JSON des LEM-Innenraummodells.");
        interfaceGroups = interior.object()["interfaces"];
    }
    if (!interfaceGroups.isUndefined()) {
        if (!interfaceGroups.isArray() || interfaceGroups.toArray().isEmpty())
            throw std::invalid_argument("LEM-Interfaces benötigen eine nicht leere Gruppenliste.");
        std::set<std::uint32_t> bound;
        for (const auto &value : interfaceGroups.toArray()) {
            if (!value.isObject())
                throw std::invalid_argument("LEM-Interfacegruppe muss ein Objekt sein.");
            const auto group = value.toObject();
            const auto sources = group["sourceFaces"].toArray();
            const auto volume = group["volume_m3"].toDouble(-1);
            const auto polarity = group["polarity"].toDouble(1);
            const auto sourceGroups = group["sourceGroups"].toArray();
            if ((sources.isEmpty() == sourceGroups.isEmpty()) || !group["volume_m3"].isDouble() ||
                !std::isfinite(volume) || volume <= 0 || (polarity != 1 && polarity != -1) ||
                (group.contains("polarity") && !group["polarity"].isDouble()))
                throw std::invalid_argument(
                    "LEM-Interface benötigt Flächen, positives SI-Volumen und ±1-Polarität.");
            std::vector<std::uint32_t> faces;
            for (const auto &source : sources) {
                const auto found = sourceFaces.find(identifier(source));
                if (found == sourceFaces.end() || !interfaceFaces.contains(found->second) ||
                    !bound.insert(found->second).second)
                    throw std::invalid_argument(
                        "LEM-Interfacefläche fehlt, hat falsche Rolle oder ist doppelt gebunden.");
                faces.push_back(found->second);
            }
            for (const auto &sourceGroup : sourceGroups) {
                const auto found = semanticGroups.find(sourceGroup.toString());
                if (!sourceGroup.isString() || found == semanticGroups.end())
                    throw std::invalid_argument("LEM-Interface: unbekannte semantische Gruppe.");
                for (const auto id : found->second) {
                    if (!interfaceFaces.contains(id) || !bound.insert(id).second)
                        throw std::invalid_argument(
                            "LEM-Interfacegruppe hat falsche oder doppelte Flächen.");
                    faces.push_back(id);
                }
            }
            const auto node = lem::NodeId{static_cast<std::uint32_t>(network.nodes.size() + 1)};
            const auto name = group["name"].toString("Innenvolumen").toStdString();
            const double x = static_cast<double>(network.nodes.size() - 1) * 260;
            network.nodes.push_back({node, lem::Domain::acoustic, false, name, {x, 0}});
            const lem::Terminals terminals{node, {1}};
            network.elements.push_back({{static_cast<std::uint32_t>(network.elements.size() + 1)},
                                        name + " · Interface",
                                        {x - 60, 90},
                                        lem::BemPort{terminals, faces, static_cast<int>(polarity),
                                                     lem::BemPortKind::interfacePort}});
            auto cavityNode = node;
            if (group.contains("neck")) {
                if (!group["neck"].isObject())
                    throw std::invalid_argument("LEM-Spaltkanal muss ein Parameterobjekt sein.");
                const auto neck = group["neck"].toObject();
                const auto length = [&](const char *key) {
                    const auto parameterValue = neck[key];
                    if (!parameterValue.isDouble() || !std::isfinite(parameterValue.toDouble()) ||
                        parameterValue.toDouble() <= 0)
                        throw std::invalid_argument(std::string("Ungültiger Spaltparameter: ") +
                                                    key);
                    return parameterValue.toDouble();
                };
                cavityNode = lem::NodeId{static_cast<std::uint32_t>(network.nodes.size() + 1)};
                network.nodes.push_back({cavityNode,
                                         lem::Domain::acoustic,
                                         false,
                                         "Druck hinter Spaltkanal",
                                         {x + 180, 0}});
                network.elements.push_back(
                    {{static_cast<std::uint32_t>(network.elements.size() + 1)},
                     "Enger Ringspalt · viskos",
                     {x + 80, 90},
                     lem::AcousticSlit{{node, cavityNode},
                                       length("length_m"),
                                       length("area_m2"),
                                       length("gap_m"),
                                       length("dynamicViscosity_Pa_s")}});
            }
            lem::Component cavity = lem::Cavity{{cavityNode, {1}}, volume};
            if (group.contains("closedSlitCavity")) {
                if (!group["closedSlitCavity"].isObject())
                    throw std::invalid_argument(
                        "Die geschlossene Ringkammer benötigt ein Parameterobjekt.");
                const auto chamber = group["closedSlitCavity"].toObject();
                const auto parameter = [&](const char *key) {
                    const auto parameterValue = chamber[key];
                    if (!parameterValue.isDouble() || !std::isfinite(parameterValue.toDouble()) ||
                        parameterValue.toDouble() <= 0)
                        throw std::invalid_argument(
                            std::string("Ungültiger Ringkammerparameter: ") + key);
                    return parameterValue.toDouble();
                };
                const double depth = parameter("length_m"), area = parameter("area_m2");
                if (std::abs(depth * area - volume) > 1e-10 * volume)
                    throw std::invalid_argument(
                        "Volumen der Ringkammer passt nicht zu Tiefe mal Querschnitt.");
                cavity = lem::ClosedSlitCavity{{cavityNode, {1}},
                                               depth,
                                               area,
                                               parameter("gap_m"),
                                               parameter("dynamicViscosity_Pa_s"),
                                               parameter("heatCapacityRatio"),
                                               parameter("prandtlNumber")};
            }
            network.elements.push_back({{static_cast<std::uint32_t>(network.elements.size() + 1)},
                                        name,
                                        {x + (cavityNode == node ? 60 : 230), 90},
                                        std::move(cavity)});
        }
        if (bound != interfaceFaces)
            throw std::invalid_argument(
                "Nicht sämtliche mathematischen Interfaceflächen sind an LEM gebunden.");
    } else if (!interfaceFaces.empty())
        throw std::invalid_argument("Interfaceflächen erfordern einen expliziten LEM-Anschluss.");
    if (!driverReferencePath.isEmpty()) {
        QFile referenceFile(driverReferencePath);
        if (!referenceFile.open(QIODevice::ReadOnly) || referenceFile.size() > 1024 * 1024)
            throw std::invalid_argument("LEM-Treiberreferenz ist nicht lesbar oder zu groß.");
        QJsonParseError referenceError;
        const auto referenceDocument =
            QJsonDocument::fromJson(referenceFile.readAll(), &referenceError);
        if (referenceError.error != QJsonParseError::NoError || !referenceDocument.isObject())
            throw std::invalid_argument("Ungültiges JSON der LEM-Treiberreferenz.");
        const auto reference = referenceDocument.object();
        const auto parameters = reference["parameters"].toObject();
        const auto number = [&](const char *key, bool allowZero = false) {
            const auto value = parameters[key];
            if (!value.isDouble() || !std::isfinite(value.toDouble()) ||
                (allowZero ? value.toDouble() < 0 : value.toDouble() <= 0))
                throw std::invalid_argument(std::string("Ungültiger LEM-Treiberparameter: ") + key);
            return value.toDouble();
        };
        const auto weights = reference["sourceGroups"].toObject();
        if (weights.isEmpty())
            throw std::invalid_argument("Treiberreferenz benötigt gewichtete Quellgruppen.");
        lem::BemPort port;
        port.kind = lem::BemPortKind::drivenSurface;
        double modeArea = 0;
        for (auto it = weights.begin(); it != weights.end(); ++it) {
            const double weight = it.value().toDouble(-1);
            const auto group = semanticGroups.find(it.key());
            if (!it.value().isDouble() || !std::isfinite(weight) || weight <= 0 ||
                group == semanticGroups.end())
                throw std::invalid_argument(
                    "Treiberreferenz: unbekannte Gruppe oder ungültiges Gewicht.");
            for (const auto id : group->second) {
                const auto face = std::find_if(cad->faces.begin(), cad->faces.end(),
                                               [id](const auto &f) { return f.id == id; });
                if (face == cad->faces.end() || face->boundary != BoundaryCondition::driven)
                    throw std::invalid_argument(
                        "Ein Treibermodus benötigt ausschließlich Driven-Flächen.");
                port.faceIds.push_back(id);
                port.velocityWeights.emplace(id, weight);
                modeArea += weight * face->signature.area;
            }
        }
        // The mechanical velocity is the mode's unit dome velocity. Its force
        // and volume flow are conjugate only with this same weighted area.
        const lem::DriverParameters driverParameters{number("resistance_ohm"),
                                                     number("inductance_H", true),
                                                     number("forceFactor_N_A"),
                                                     number("movingMassInVacuo_kg"),
                                                     number("compliance_m_N"),
                                                     number("mechanicalResistance_N_s_m"),
                                                     modeArea};
        const auto amplitude = reference["voltageRms_V"].toArray();
        if (amplitude.size() != 2 || !amplitude[0].isDouble() || !amplitude[1].isDouble() ||
            !std::isfinite(amplitude[0].toDouble()) || !std::isfinite(amplitude[1].toDouble()))
            throw std::invalid_argument("Treiberquelle benötigt eine komplexe RMS-Spannung.");
        const auto name = reference["name"].toString();
        if (name.isEmpty())
            throw std::invalid_argument("Treiberreferenz benötigt einen sichtbaren Namen.");
        const auto nextNode = static_cast<std::uint32_t>(network.nodes.size() + 1);
        const lem::NodeId front{nextNode}, electricalReference{nextNode + 1},
            electrical{nextNode + 2};
        network.nodes.push_back({front, lem::Domain::acoustic, false, "HF-Frontdruck", {460, 0}});
        network.nodes.push_back({electricalReference,
                                 lem::Domain::electrical,
                                 true,
                                 "Elektrische Referenz",
                                 {700, 180}});
        network.nodes.push_back(
            {electrical, lem::Domain::electrical, false, "HF-Eingang · RMS", {700, 0}});
        port.terminals = {front, {1}};
        const lem::Terminals electricalTerminals{electrical, electricalReference};
        const auto addElement = [&](std::string label, lem::Position position,
                                    lem::Component component) {
            network.elements.push_back({{static_cast<std::uint32_t>(network.elements.size() + 1)},
                                        std::move(label),
                                        position,
                                        std::move(component)});
        };
        addElement("HF · Kalotte 1 / gesamte Sicke 0,5", {420, 90}, std::move(port));
        addElement(name.toStdString(), {580, 90},
                   lem::Driver{electricalTerminals, front, {1}, driverParameters});
        addElement("HF-Spannungsquelle · RMS", {760, 90},
                   lem::EffortSource{electricalTerminals,
                                     {amplitude[0].toDouble(), amplitude[1].toDouble()}});
        std::cout << "LEM-Treiber: " << name.toStdString() << "; Modusfläche " << modeArea * 1e4
                  << " cm²; Kalotte/Sicke 1/0,5; RMS " << amplitude[0].toDouble() << " V\n";
    }
    if (!network.elements.empty()) {
        lem::validateNetwork(network);
        project.simulation.settings.couplingNetwork = std::move(network);
        project.simulation.settings.couplingEnabled = true;
    }
    (void)validateCadSnapshot(*cad);
    return project;
}
} // namespace

int main(int argc, char **argv) {
    QCoreApplication application(argc, argv);
    QCommandLineParser parser;
    parser.setApplicationDescription(
        "STEP-Modell mit semantischen Flächen für freie oder baffelgebundene Rechnung vorbereiten.");
    parser.addHelpOption();
    parser.addOption(QCommandLineOption(
        "manifest", "Modellmanifest mit STEP-Datei und SI-Flächensignaturen.", "json"));
    parser.addOption(QCommandLineOption("input-project", "Gespeichertes CAD-Projekt als Ausgangsstand.", "hsproj"));
    parser.addOption(QCommandLineOption("output", "Zielprojekt.", "hsproj"));
    parser.addOption(QCommandLineOption(
        "snapshot", "Nativer FreeCAD-Snapshot mit persistenten Dokument-/Objektidentitäten.",
        "json"));
    parser.addOption(QCommandLineOption(
        "check-update",
        "Geänderten nativen Snapshot gegen gespeicherte Antriebs-/Aperturidentitäten prüfen.",
        "json"));
    parser.addOption(QCommandLineOption("mesh", "Isotropes Netz aufbauen."));
    parser.addOption(QCommandLineOption(
        "reuse-mesh", "Geprüftes Netz eines exakt identischen CAD-/Netzstands verwenden.",
        "hsproj"));
    parser.addOption(QCommandLineOption(
        "driver-reference", "Explizite LEM-Treiberparameter, RMS-Spannung und Quellgruppen.",
        "json"));
    parser.addOption(QCommandLineOption(
        "interior-reference", "Explizites LEM-Innenraummodell mit optionalem Spaltkanal.", "json"));
    parser.addOption(QCommandLineOption("solve", "Netz und CUDA-Schallwandrechnung ausführen."));
    parser.addOption(
        QCommandLineOption("axisymmetric", "2D-Rotation: Meridianraum m=0 statt 3D-Matrix lösen."));
    parser.addOption(QCommandLineOption(
        "sectors", "Sektoren des gespeicherten Rotationsnetzes [8..4096, durch 4 teilbar].",
        "count", "64"));
    parser.addOption(QCommandLineOption("meridian-step-mm",
                                        "Maximale Meridianabschnittslänge [mm].", "mm", "2"));
    parser.addOption(
        QCommandLineOption("meridian-chord-mm", "Meridian-Sehnenabweichung [mm].", "mm", "0.05"));
    parser.addOption(QCommandLineOption("mesh-scale",
                                        "Skalierung aller globalen/lokalen Zielkanten; Genauigkeit "
                                        "und Qualitätsgrenzen bleiben gleich.",
                                        "factor", "1"));
    parser.addOption(QCommandLineOption(
        "frequencies", "Anzahl logarithmischer Frequenzen im angegebenen Bereich.", "count", "7"));
    parser.addOption(
        QCommandLineOption("minimum-frequency", "Untere Frequenz [Hz].", "Hz", "1000"));
    parser.addOption(
        QCommandLineOption("maximum-frequency", "Obere Frequenz [Hz].", "Hz", "10000"));
    parser.addOption(
        QCommandLineOption("near-depth", "Maximale Nahquadraturtiefe [0–10].", "depth", "6"));
    parser.addOption(QCommandLineOption(
        "export-integration-mesh",
        "2D-Integrationsnetz einschließlich Meridianzuordnungen als SI-JSON exportieren.", "json"));
    parser.addOption(QCommandLineOption("field-csv", "H-Schallfeld aus dem ersten Ergebnis als CSV auswerten.", "csv"));
    parser.addOption(QCommandLineOption("field-lateral-mm", "Halbe Querbreite des Felds [mm].", "mm", "150"));
    parser.addOption(QCommandLineOption("field-axial-start-mm", "Axialer Feldbeginn [mm].", "mm", "-50"));
    parser.addOption(QCommandLineOption("field-axial-end-mm", "Axiales Feldende [mm].", "mm", "350"));
    parser.addOption(QCommandLineOption("field-grid", "Rasterpunkte je Richtung [3–201].", "count", "21"));
    parser.process(application);
    try {
        if (parser.isSet("manifest") == parser.isSet("input-project") || !parser.isSet("output"))
            throw std::invalid_argument("Genau ein --manifest oder --input-project sowie --output sind erforderlich.");
        if (parser.isSet("input-project") &&
            (parser.isSet("snapshot") || parser.isSet("driver-reference") ||
             parser.isSet("interior-reference") || parser.isSet("reuse-mesh")))
            throw std::invalid_argument("Gespeicherte Projekte enthalten bereits ihren vollständigen CAD- und Netzstand.");
        if (parser.isSet("field-csv") && !parser.isSet("solve") && !parser.isSet("input-project"))
            throw std::invalid_argument("--field-csv benötigt eine Simulation oder ein gespeichertes Ergebnis.");
        auto project = parser.isSet("input-project")
            ? readProject(parser.value("input-project"))
            : prepare(parser.value("manifest"), parser.value("snapshot"),
                      parser.value("driver-reference"), parser.value("interior-reference"));
        if (!project.cad || !project.cad->geometry ||
            (parser.isSet("input-project") && (!project.mesh ||
                                                project.meshRevision != project.geometryRevision)))
            throw std::invalid_argument("Gespeicherter CAD- und Qualitätsnetzstand ist nicht aktuell.");
        bool scaleValid = false, countValid = false, minimumValid = false, maximumValid = false,
             depthValid = false;
        const double meshScale = parser.value("mesh-scale").toDouble(&scaleValid);
        const auto frequencyCount = parser.value("frequencies").toUInt(&countValid);
        const auto minimumFrequency = parser.value("minimum-frequency").toDouble(&minimumValid);
        const auto maximumFrequency = parser.value("maximum-frequency").toDouble(&maximumValid);
        const auto nearDepth = parser.value("near-depth").toUInt(&depthValid);
        if (!scaleValid || !std::isfinite(meshScale) || meshScale <= 0 || !countValid ||
            !frequencyCount || !minimumValid || !maximumValid || !std::isfinite(minimumFrequency) ||
            !std::isfinite(maximumFrequency) || minimumFrequency <= 0 ||
            maximumFrequency < minimumFrequency ||
            (frequencyCount > 1 && minimumFrequency == maximumFrequency) || !depthValid ||
            nearDepth > 10)
            throw std::invalid_argument("Ungültige Netzskalierung oder Frequenzfolge.");
        if (parser.isSet("input-project") && (meshScale != 1 || parser.isSet("axisymmetric")))
            throw std::invalid_argument(
                "Gespeichertes Netz darf nicht mit geänderter Netzskalierung oder Geometrieart wiederverwendet werden.");
        project.meshing.targetEdgeLength *= meshScale;
        project.meshing.adaptive.minimumEdgeLength *= meshScale;
        for (auto &[id, length] : project.cadFaceTargetEdgeLengths)
            length *= meshScale;
        validateMeshSettings(project.meshing);
        if (!parser.isSet("input-project") || parser.isSet("solve")) {
            project.simulation.settings.frequencyCount = frequencyCount;
            project.simulation.settings.minimumFrequency = minimumFrequency;
            project.simulation.settings.maximumFrequency = maximumFrequency;
            project.simulation.settings.nearMaximumDepth = nearDepth;
        }
        std::optional<AxisymmetricMesh> rotational;
        std::optional<AxisymmetricBem> rotationalLayout;
        if (parser.isSet("axisymmetric")) {
            if (parser.isSet("reuse-mesh") || parser.isSet("mesh"))
                throw std::invalid_argument("2D-Rotation erzeugt ihr Meridian-Integrationsnetz "
                                            "direkt; --mesh/--reuse-mesh gehören zum 3D-Pfad.");
            bool sectorsValid = false, stepValid = false, chordValid = false;
            auto &settings = project.simulation.settings;
            auto &sampling = settings.axisymmetricSampling;
            sampling.sectors = parser.value("sectors").toUInt(&sectorsValid);
            sampling.maximumLength = parser.value("meridian-step-mm").toDouble(&stepValid) * .001;
            sampling.chordTolerance =
                parser.value("meridian-chord-mm").toDouble(&chordValid) * .001;
            if (!sectorsValid || !stepValid || !chordValid ||
                !std::isfinite(sampling.maximumLength) || sampling.maximumLength <= 0 ||
                !std::isfinite(sampling.chordTolerance) || sampling.chordTolerance <= 0)
                throw std::invalid_argument("2D-Rotation: ungültige Diskretisierung.");
            settings.geometryMode = BemGeometryMode::axisymmetric;
            MeridianSampling geometrySampling;
            geometrySampling.maximumLength = sampling.maximumLength;
            geometrySampling.chordTolerance = sampling.chordTolerance;
            geometrySampling.rotationalTolerance = sampling.rotationalTolerance;
            geometrySampling.faceTargetLengths = project.cadFaceTargetEdgeLengths;
            for (const auto &face : project.cad->faces)
                geometrySampling.faceConditions.emplace(face.id, face.boundary);
            const auto start = std::chrono::steady_clock::now();
            const auto meridian =
                makeAxisymmetricMeridian(*project.cad->geometry, project.cad->acousticOrigin,
                                         project.cad->acousticAxis, geometrySampling);
            rotational = buildAxisymmetricMesh(meridian, sampling.sectors);
            rotationalLayout = AxisymmetricBem{
                project.cad->acousticOrigin, project.cad->acousticAxis, sampling.sectors,
                rotational->pressureDofs,    rotational->fluxDofs,      rotational->testPanels};
            std::cout
                << "2D-Meridian: " << meridian.points.size() << " Druckknoten, "
                << meridian.segments.size() << " Abschnitte; " << sampling.sectors
                << " Sektoren im Ergebnisnetz; Aufbau "
                << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
                << " s.\n"
                << std::flush;
        }
        if (parser.isSet("check-update")) {
            const auto update = nativeSnapshot(parser.value("check-update"), project.cad.get());
            for (const auto &face : project.cad->faces) {
                const bool driven = face.boundary == BoundaryCondition::driven;
                const auto &baffle = project.simulation.settings.infiniteBaffle;
                const bool aperture = baffle &&
                    std::find(baffle->apertureFaceIds.begin(), baffle->apertureFaceIds.end(),
                              face.id) != baffle->apertureFaceIds.end();
                bool interface = false;
                for (const auto *element :
                     lem::bemPorts(project.simulation.settings.couplingNetwork)) {
                    const auto &port = std::get<lem::BemPort>(element->component);
                    interface |= port.kind == lem::BemPortKind::interfacePort &&
                                 std::find(port.faceIds.begin(), port.faceIds.end(), face.id) !=
                                     port.faceIds.end();
                }
                if (!driven && !aperture && !interface)
                    continue;
                const auto found =
                    std::find_if(update.faces.begin(), update.faces.end(),
                                 [&](const auto &candidate) { return candidate.id == face.id; });
                if (found == update.faces.end() || found->objectId != face.objectId ||
                    found->boundary != face.boundary ||
                    found->normalVelocity != face.normalVelocity ||
                    found->state != CadFaceState::confirmed)
                    throw std::runtime_error(
                        "Native Ringänderung verlor Antrieb, Apertur oder LEM-Interfacezuordnung.");
            }
            std::cout << "Native Ringänderung erhält Antriebs-, Apertur- und Interfaceidentität.\n";
        }
        if (parser.isSet("export-integration-mesh")) {
            if (!rotational)
                throw std::invalid_argument("Integrationsnetzexport benötigt --axisymmetric.");
            QJsonArray vertices, triangles, pressure, flux, targets;
            for (const auto p : rotational->mesh.vertices)
                vertices.append(QJsonArray{p.x, p.y, p.z});
            for (const auto &t : rotational->mesh.triangles)
                triangles.append(QJsonArray{
                    static_cast<qint64>(t.vertices[0]), static_cast<qint64>(t.vertices[1]),
                    static_cast<qint64>(t.vertices[2]), static_cast<qint64>(t.faceId),
                    static_cast<int>(t.boundary)});
            for (const auto v : rotational->pressureDofs)
                pressure.append(static_cast<qint64>(v));
            for (const auto v : rotational->fluxDofs)
                flux.append(static_cast<qint64>(v));
            for (const auto v : rotational->testPanels)
                targets.append(static_cast<qint64>(v));
            const auto origin = project.cad->acousticOrigin, axis = project.cad->acousticAxis;
            const auto bytes =
                QJsonDocument(
                    QJsonObject{
                        {"vertices", vertices},
                        {"triangles", triangles},
                        {"origin", QJsonArray{origin.x, origin.y, origin.z}},
                        {"axis", QJsonArray{axis.x, axis.y, axis.z}},
                        {"sectors", static_cast<int>(
                                        project.simulation.settings.axisymmetricSampling.sectors)},
                        {"pressureDofs", pressure},
                        {"fluxDofs", flux},
                        {"testPanels", targets}})
                    .toJson(QJsonDocument::Compact);
            QFile exported(parser.value("export-integration-mesh"));
            if (!exported.open(QIODevice::WriteOnly) || exported.write(bytes) != bytes.size())
                throw std::runtime_error("2D-Integrationsnetzexport fehlgeschlagen.");
        }
        const auto path = QFileInfo(parser.value("output")).absoluteFilePath();
        if (parser.isSet("reuse-mesh"))
            reuseMesh(project, parser.value("reuse-mesh"));
        if (!parser.value("driver-reference").isEmpty())
            alignReferenceModeArea(project, rotational ? &rotational->mesh : nullptr,
                                   rotationalLayout);
        saveChecked(path, project);
        if ((parser.isSet("mesh") || parser.isSet("solve")) && !project.mesh && !rotational) {
            const auto meshStart = std::chrono::steady_clock::now();
            const auto source =
                makeCadMeshingSource(*project.cad, project.meshing.adaptive.geometricTolerance);
            MeshingOptions options;
            options.targetEdgeLength = project.meshing.targetEdgeLength;
            options.adaptive = project.meshing.adaptive;
            options.quality = project.meshing.quality;
            options.faceTargetEdgeLengths = project.cadFaceTargetEdgeLengths;
            options.geometry = project.cad->geometry;
            options.parameterization = source.parameters;
            std::cout << "Isotrope Vernetzung ...\n" << std::flush;
            project.mesh = std::make_shared<MeshingResult>(remeshIsotropic(source.mesh, options));
            project.meshRevision = project.revision;
            if (!parser.value("driver-reference").isEmpty())
                alignReferenceModeArea(project);
            std::cout << "Netz: " << project.mesh->mesh.triangles.size()
                      << " Dreiecke; min. Winkel " << project.mesh->quality.minimumAngleDegrees
                      << "; max. Seitenverhältnis " << project.mesh->quality.maximumAspectRatio
                      << "; Zeit "
                      << std::chrono::duration<double>(std::chrono::steady_clock::now() - meshStart)
                             .count()
                      << " s\n";
            saveChecked(path, project);
        }
        if (parser.isSet("solve")) {
            std::vector<FaceVelocity> velocities;
            for (const auto &face : project.cad->faces)
                if (face.boundary == BoundaryCondition::driven)
                    velocities.push_back({face.id, face.normalVelocity});
            auto request =
                makeSimulationRequest(project.simulation.settings, velocities,
                                      project.cad->acousticOrigin, project.cad->acousticAxis);
            request.axisymmetric = rotationalLayout;
            const auto &solveMesh = rotational ? rotational->mesh : project.mesh->mesh;
            CudaBemSolver solver(solveMesh);
            auto result =
                std::make_shared<BemResult>(solver.solve(request, {}, [](const BemProgress &p) {
                    std::cout << p.frequency << " Hz fertig (" << p.completedFrequencies << '/'
                              << p.totalFrequencies << ")\n"
                              << std::flush;
                }));
            project.simulation.mesh = std::make_shared<SurfaceMesh>(solveMesh);
            if (rotational)
                project.simulation.axisymmetricSampling =
                    project.simulation.settings.axisymmetricSampling;
            project.simulation.request = request;
            project.simulation.result = result;
            project.simulation.revision = project.revision;
            project.simulation.current = true;
            saveChecked(path, project);
            QFile csv(QFileInfo(path).dir().absoluteFilePath(QFileInfo(path).completeBaseName() +
                                                             ".csv"));
            const auto bytes = resultCsv(project.simulation);
            if (!csv.open(QIODevice::WriteOnly) || csv.write(bytes) != bytes.size())
                throw std::runtime_error("Ergebnis-CSV konnte nicht gespeichert werden.");
        }
        if (parser.isSet("field-csv")) {
                const auto &result = project.simulation.result;
                if (!project.simulation.current || !result || !project.simulation.mesh ||
                    !project.simulation.request)
                    throw std::invalid_argument("Kein aktuelles gespeichertes Simulationsergebnis vorhanden.");
                const auto &request = *project.simulation.request;
                const auto &solveMesh = *project.simulation.mesh;
                bool widthOk = false, startOk = false, endOk = false, gridOk = false;
                const auto width = parser.value("field-lateral-mm").toDouble(&widthOk) * .001;
                const auto start = parser.value("field-axial-start-mm").toDouble(&startOk) * .001;
                const auto end = parser.value("field-axial-end-mm").toDouble(&endOk) * .001;
                const auto grid = parser.value("field-grid").toUInt(&gridOk);
                if (!widthOk || !startOk || !endOk || !gridOk || !std::isfinite(width) ||
                    !std::isfinite(start) || !std::isfinite(end) || width <= 0 || end <= start ||
                    grid < 3 || grid > 201 || result->frequencies.empty() ||
                    result->frequencies.front().excitations.empty())
                    throw std::invalid_argument("Ungültiges Schallfeldraster oder fehlendes Ergebnis.");
                const auto frame = *project.simulation.settings.observationFrame;
                std::vector<Vec3> points;
                for (unsigned z = 0; z < grid; ++z)
                    for (unsigned x = 0; x < grid; ++x) {
                        const double ax = start + (end - start) * z / (grid - 1);
                        const double lat = -width + 2 * width * x / (grid - 1);
                        points.push_back({frame.origin.x + ax * frame.forward.x + lat * frame.horizontal.x,
                                          frame.origin.y + ax * frame.forward.y + lat * frame.horizontal.y,
                                          frame.origin.z + ax * frame.forward.z + lat * frame.horizontal.z});
                    }
                const auto &frequency = result->frequencies.front();
                const auto field = evaluatePressureField(
                    solveMesh, request, frequency.frequency, frequency.excitations.front(), points,
                    {.0001, request.numerics.nearRelativeTolerance, request.numerics.nearMaximumDepth});
                QByteArray data("# SI world coordinates; complex RMS pressure; exp(+i omega t)\n");
                data += "# frequency_Hz=" + QByteArray::number(frequency.frequency, 'g', 17) + '\n';
                data += "x_m,y_m,z_m,status,pressure_real_Pa,pressure_imag_Pa\n";
                std::size_t valid = 0, unresolved = 0;
                double minDb = std::numeric_limits<double>::infinity();
                double maxDb = -std::numeric_limits<double>::infinity();
                for (const auto &sample : field) {
                    data += QByteArray::number(sample.position.x, 'g', 17) + ',' +
                            QByteArray::number(sample.position.y, 'g', 17) + ',' +
                            QByteArray::number(sample.position.z, 'g', 17) + ',';
                    const bool air = sample.status == FieldPointStatus::air;
                    valid += air;
                    unresolved += sample.status == FieldPointStatus::quadratureNotConverged;
                    const char *status = air ? "air" : sample.status == FieldPointStatus::boundary
                        ? "boundary" : sample.status == FieldPointStatus::excluded
                        ? "excluded" : "quadrature_not_converged";
                    data += QByteArray(status) + ',';
                    if (air)
                    {
                        data += QByteArray::number(sample.pressure.real(), 'g', 17) + ',' +
                                QByteArray::number(sample.pressure.imag(), 'g', 17);
                        const double db = soundPressureLevelDecibels(sample.pressure);
                        if (std::isfinite(db)) {
                            minDb = std::min(minDb, db);
                            maxDb = std::max(maxDb, db);
                        }
                    }
                    else
                        data += ',';
                    data += '\n';
                }
                QFile output(parser.value("field-csv"));
                if (!output.open(QIODevice::WriteOnly) || output.write(data) != data.size())
                    throw std::runtime_error("Schallfeld-CSV konnte nicht gespeichert werden.");
                auto &display = project.simulation.settings.fieldDisplay;
                display.minimumLateral = -width;
                display.maximumLateral = width;
                display.minimumAxial = start;
                display.maximumAxial = end;
                display.lateralCount = grid;
                display.axialCount = grid;
                if (std::isfinite(minDb) && std::isfinite(maxDb)) {
                    display.minimumDb = std::max(-300., std::floor(minDb / 5) * 5);
                    display.maximumDb = std::min(300., std::ceil(maxDb / 5) * 5);
                    if (display.maximumDb <= display.minimumDb)
                        display.maximumDb = std::min(300., display.minimumDb + 5);
                    if (display.maximumDb <= display.minimumDb)
                        display.minimumDb = display.maximumDb - 5;
                }
                validateFieldDisplaySettings(display);
                saveChecked(path, project);
                std::cout << "Schallfeld " << frequency.frequency << " Hz: " << valid
                          << " Luftpunkte, " << unresolved << " ohne Quadraturkonvergenz.\n";
        }
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "CAD-Modellprojekt: " << e.what() << '\n';
        return 1;
    }
}
