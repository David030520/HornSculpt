#include "galerkin/exterior_domain.hpp"
#include "gmsh_source.hpp"
#include "project_json.hpp"
#include "simulation_request.hpp"
#include <QFile>
#include <QJsonDocument>
#include <QSaveFile>
#include <QSet>
#include <algorithm>
#include <cmath>
#include <hornsculpt/mesh_orientation.hpp>
#include <limits>
#include <stdexcept>
#include <vector>

namespace hornsculpt {
namespace {
using namespace project_json;
constexpr qint64 maximumFileBytes = 256 * 1024 * 1024;

void uniqueKeys(const QByteArray &bytes) {
    // Qt keeps the last duplicate key; reject ambiguous input before schema decoding.
    std::vector<QSet<QString>> stack;
    for (qsizetype i = 0; i < bytes.size(); ++i) {
        const char c = bytes[i];
        if (c == '{' || c == '[') {
            stack.emplace_back();
            if (stack.size() > 256)
                invalid(QStringLiteral("Zu tiefe JSON-Verschachtelung."));
        } else if (c == '}' || c == ']') {
            if (!stack.empty())
                stack.pop_back();
        } else if (c == '"') {
            const qsizetype begin = i++;
            for (; i < bytes.size(); ++i) {
                if (bytes[i] == '\\')
                    ++i;
                else if (bytes[i] == '"')
                    break;
            }
            qsizetype next = i + 1;
            while (next < bytes.size() && (bytes[next] == ' ' || bytes[next] == '\t' ||
                                           bytes[next] == '\n' || bytes[next] == '\r'))
                ++next;
            if (next < bytes.size() && bytes[next] == ':' && !stack.empty()) {
                const auto encoded = QByteArray("[") + bytes.mid(begin, i - begin + 1) + ']';
                const auto key = QJsonDocument::fromJson(encoded).array()[0].toString();
                if (stack.back().contains(key))
                    invalid(QStringLiteral("Doppelter JSON-Schlüssel."));
                stack.back().insert(key);
            }
        }
    }
}

AdaptiveSizing decodeAdaptive(const QJsonValue &value) {
    const auto sizing = object(value);
    keys(sizing, {"enabled", "minimumEdgeLength", "geometricTolerance", "growthRate"});
    return {boolean(sizing["enabled"]), positive(sizing["minimumEdgeLength"]),
            positive(sizing["geometricTolerance"]), positive(sizing["growthRate"])};
}

std::map<std::uint32_t, double> decodeFaceSizes(const QJsonValue &value) {
    std::map<std::uint32_t, double> result;
    for (const auto &item : array(value, 100000)) {
        const auto pair = array(item, 2, 2);
        const auto id = integer(pair[0]);
        if (!id || !result.emplace(id, positive(pair[1])).second)
            invalid(QStringLiteral("Doppelte oder ungültige lokale Netzfläche."));
    }
    return result;
}

QJsonArray encodeFaceSizes(const std::map<std::uint32_t, double> &sizes) {
    QJsonArray result;
    for (const auto &[id, size] : sizes)
        result.append(QJsonArray{static_cast<double>(id), size});
    return result;
}

FaceColors decodeFaceColors(const QJsonValue &value) {
    FaceColors result;
    for (const auto &item : array(value, 100000)) {
        const auto pair = array(item, 2, 2);
        const auto id = integer(pair[0]);
        const auto rgb = array(pair[1], 3, 3);
        const FaceColor color{static_cast<std::uint8_t>(integer(rgb[0], 255)),
                              static_cast<std::uint8_t>(integer(rgb[1], 255)),
                              static_cast<std::uint8_t>(integer(rgb[2], 255))};
        if (!id || !result.emplace(id, color).second)
            invalid(QStringLiteral("Doppelte oder ungültige Farbfläche."));
    }
    return result;
}
QJsonArray encodeFaceColors(const FaceColors &colors) {
    QJsonArray result;
    for (const auto &[id, color] : colors)
        result.append(QJsonArray{double(id), QJsonArray{color.red, color.green, color.blue}});
    return result;
}
void validateFaceColors(const Project &project) {
    for (const auto &[id, color] : project.draftFaceColors) {
        (void)color;
        if (!id || id > static_cast<std::uint32_t>(DraftFace::upstreamOuterWall))
            invalid(QStringLiteral("Unbekannte Draft-Fläche für Farbzuweisung."));
    }
    const auto faces =
        project.cad ? meshFaces(project.cad->mesh) : std::map<std::uint32_t, BoundaryCondition>{};
    for (const auto &[id, color] : project.cadFaceColors) {
        (void)color;
        if (!faces.contains(id))
            invalid(QStringLiteral("Unbekannte CAD-Fläche für Farbzuweisung."));
    }
}

std::set<std::uint32_t> decodeOrientations(const QJsonValue &value) {
    std::set<std::uint32_t> result;
    for (const auto &item : array(value, 100000)) {
        const auto id = integer(item);
        if (!id || !result.insert(id).second)
            invalid(QStringLiteral("Doppelte oder ungültige Orientierungsfläche."));
    }
    return result;
}

QJsonArray encodeOrientations(const std::set<std::uint32_t> &ids) {
    QJsonArray result;
    for (const auto id : ids)
        result.append(static_cast<double>(id));
    return result;
}

void validateOrientations(const Project &project) {
    for (const auto id : project.draftReversedFaceIds)
        if (!id || id > static_cast<std::uint32_t>(DraftFace::upstreamOuterWall))
            invalid(QStringLiteral("Unbekannte Draft-Fläche für Orientierungsumkehr."));
    const auto faces =
        project.cad ? meshFaces(project.cad->mesh) : std::map<std::uint32_t, BoundaryCondition>{};
    for (const auto id : project.cadReversedFaceIds)
        if (!faces.contains(id))
            invalid(QStringLiteral("Unbekannte CAD-Fläche für Orientierungsumkehr."));
}

void validateFaceSizes(const Project &p) {
    for (const auto *sizes : {&p.draftFaceTargetEdgeLengths, &p.cadFaceTargetEdgeLengths})
        for (const auto &[id, size] : *sizes)
            if (!id || !std::isfinite(size) || size < p.meshing.adaptive.minimumEdgeLength ||
                size > (sizes == &p.cadFaceTargetEdgeLengths && p.assembly.includeCad
                            ? p.assembly.cadTargetEdgeLength
                            : p.meshing.targetEdgeLength))
                invalid(QStringLiteral("Lokale Zielkante außerhalb der Netzeinstellungen."));
    for (const auto &[id, size] : p.draftFaceTargetEdgeLengths) {
        (void)size;
        if (id > static_cast<std::uint32_t>(DraftFace::upstreamOuterWall))
            invalid(QStringLiteral("Unbekannte Draft-Fläche für lokale Vernetzung."));
    }
    const auto cadFaces =
        p.cad ? meshFaces(p.cad->mesh) : std::map<std::uint32_t, BoundaryCondition>{};
    for (const auto &[id, size] : p.cadFaceTargetEdgeLengths) {
        (void)size;
        if (!cadFaces.contains(id))
            invalid(QStringLiteral("Unbekannte CAD-Fläche für lokale Vernetzung."));
    }
}

Project decode(const QJsonObject &root, bool allowComparisons = true) {
    if (text(root["format"]) != "HornSculpt" || text(root["units"]) != "SI")
        invalid(QStringLiteral("Unbekanntes Format oder Einheit."));
    const auto version = integer(root["version"], 18);
    auto structural = root;
    if (version >= 15) {
        if (!structural.contains("comparisonSnapshots"))
            invalid(QStringLiteral("Vergleichsstände fehlen."));
        if (!allowComparisons && !array(structural["comparisonSnapshots"], 64).isEmpty())
            invalid(QStringLiteral("Verschachtelte Vergleichsstände sind unzulässig."));
        structural.remove("comparisonSnapshots");
    }
    if (version >= 13) {
        if (!structural.contains("slot"))
            invalid(QStringLiteral("Spaltparameter fehlen."));
        structural.remove("slot");
    }
    structural.remove("gmshMetadata");
    structural.remove("assembly");
    if (version >= 11 && !structural.contains("asymmetry"))
        invalid(QStringLiteral("Asymmetrische Draft-Parameter fehlen."));
    structural.remove("asymmetry");
    structural.remove("draftFaceColors");
    structural.remove("cadFaceColors");
    if (version == 1) {
        keys(structural, {"format", "version", "units", "profileType", "profile", "horizontalAngle",
                          "verticalAngle", "crossSectionExponent", "wallThickness",
                          "targetEdgeLength", "guideEnabled", "transition", "guide"});
    } else if (version == 2 || version == 3) {
        auto common = structural;
        if (version == 3) {
            if (!common.contains("cadAdaptiveSizing"))
                invalid(QStringLiteral("Adaptive CAD-Netzeinstellungen fehlen."));
            common.remove("cadAdaptiveSizing");
        }
        keys(common, {"format",
                      "version",
                      "units",
                      "profileType",
                      "profile",
                      "horizontalAngle",
                      "verticalAngle",
                      "crossSectionExponent",
                      "wallThickness",
                      "targetEdgeLength",
                      "guideEnabled",
                      "transition",
                      "guide",
                      "source",
                      "revision",
                      "geometryRevision",
                      "meshRevision",
                      "preview",
                      "mesh",
                      "cad",
                      "cadTargetEdgeLength",
                      "simulation"});
    } else if (version >= 4) {
        auto common = structural;
        if (version >= 9) {
            if (!common.contains("draftReversedFaces") || !common.contains("cadReversedFaces"))
                invalid(QStringLiteral("Flächenorientierungen fehlen."));
            common.remove("draftReversedFaces");
            common.remove("cadReversedFaces");
        }
        keys(common, {"format",
                      "version",
                      "units",
                      "profileType",
                      "profile",
                      "length",
                      "entry",
                      "horizontalAngle",
                      "verticalAngle",
                      "crossSectionExponent",
                      "wallThickness",
                      "guideEnabled",
                      "transition",
                      "guide",
                      "source",
                      "revision",
                      "geometryRevision",
                      "meshRevision",
                      "preview",
                      "mesh",
                      "cad",
                      "meshing",
                      "draftFaceSizes",
                      "cadFaceSizes",
                      "mouthMorph",
                      "simulation"});
    } else
        invalid(QStringLiteral("Unbekannte Formatversion."));
    Project project;
    decodeDraft(root, project);
    if (root.contains("draftFaceColors"))
        project.draftFaceColors = decodeFaceColors(root["draftFaceColors"]);
    if (root.contains("cadFaceColors"))
        project.cadFaceColors = decodeFaceColors(root["cadFaceColors"]);
    if (version == 1) {
        validateMeshSettings(project.meshing);
        validateFaceColors(project);
        return project;
    }
    const auto source = text(root["source"]);
    if (version >= 12)
        project.assembly = decodeAssembly(root["assembly"]);
    if (source == "draft")
        project.source = ProjectSource::draft;
    else if (source == "freecad")
        project.source = ProjectSource::freecad;
    else
        invalid(QStringLiteral("Unbekannte Geometriequelle."));
    project.revision = wideInteger(root["revision"]);
    project.geometryRevision = wideInteger(root["geometryRevision"]);
    project.meshRevision = wideInteger(root["meshRevision"]);
    if (!project.revision || project.revision == std::numeric_limits<std::uint64_t>::max() ||
        project.geometryRevision > project.revision || project.meshRevision > project.revision)
        invalid(QStringLiteral("Ungültige Projektrevisionen."));
    if (version < 4) {
        MeshSettings cadMeshing;
        cadMeshing.targetEdgeLength = positive(root["cadTargetEdgeLength"]);
        if (version == 3)
            cadMeshing.adaptive = decodeAdaptive(root["cadAdaptiveSizing"]);
        validateMeshSettings(cadMeshing);
        // Old projects stored separate settings; the active source supplies the common set.
        if (project.source == ProjectSource::freecad)
            project.meshing = cadMeshing;
    } else {
        const auto settings = object(root["meshing"]);
        auto common = settings;
        if (version >= 5) {
            const auto quality = object(settings["quality"]);
            keys(quality, {"enabled", "minimumAngleDegrees", "maximumAspectRatio"});
            project.meshing.quality = {boolean(quality["enabled"]),
                                       positive(quality["minimumAngleDegrees"]),
                                       positive(quality["maximumAspectRatio"])};
            common.remove("quality");
        }
        if (version < 6) {
            const double discardedAngle = positive(settings["entryAngleToleranceDegrees"]);
            if (discardedAngle < .05 || discardedAngle > 5)
                invalid(QStringLiteral("Ungültige frühere Eintrittswinkel-Netzeinstellung."));
            common.remove("entryAngleToleranceDegrees");
        }
        keys(common, {"targetEdgeLength", "adaptiveSizing"});
        project.meshing.targetEdgeLength = positive(settings["targetEdgeLength"]);
        project.meshing.adaptive = decodeAdaptive(settings["adaptiveSizing"]);
        project.draftFaceTargetEdgeLengths = decodeFaceSizes(root["draftFaceSizes"]);
        if (version < 6) {
            auto &sizes = project.draftFaceTargetEdgeLengths;
            if (const auto collar = sizes.find(6); collar != sizes.end()) {
                const auto [inner, inserted] = sizes.emplace(1, collar->second);
                if (!inserted)
                    inner->second = std::min(inner->second, collar->second);
                sizes.erase(collar);
            }
        }
        project.cadFaceTargetEdgeLengths = decodeFaceSizes(root["cadFaceSizes"]);
    }
    validateMeshSettings(project.meshing);
    if (!root["preview"].isNull())
        project.preview = std::make_shared<SurfaceMesh>(decodeMesh(root["preview"], false));
    if (!root["mesh"].isNull()) {
        auto mesh = std::make_shared<MeshingResult>();
        mesh->mesh = decodeMesh(root["mesh"], !root.contains("gmshMetadata"));
        mesh->quality = inspectMesh(mesh->mesh);
        project.mesh = std::move(mesh);
    }
    if (!root["cad"].isNull())
        project.cad = std::make_shared<CadSnapshot>(decodeCad(root["cad"], version >= 7));
    if (root.contains("gmshMetadata")) {
        if (!project.cad || project.source != ProjectSource::freecad ||
            project.cad->sessionId != "gmsh-import-offline")
            invalid(QStringLiteral("Gmsh-Metadaten benötigen eine importierte Polygonquelle."));
        project.gmsh = decodeGmshMetadata(object(root["gmshMetadata"]), project.cad->mesh);
    } else if (project.cad && project.cad->sessionId == "gmsh-import-offline")
        invalid(QStringLiteral("Gmsh-Metadaten fehlen."));
    if (isPlacedAssembly(project.assembly) && project.source != ProjectSource::draft)
        invalid(QStringLiteral("Eine Draft-Baugruppe benötigt den Draft-Arbeitsbereich."));
    if (project.assembly.includeCad && !project.cad)
        invalid(QStringLiteral("CAD-Geometrie der Draft-Baugruppe fehlt."));
    validateFaceSizes(project);
    validateFaceColors(project);
    if (version >= 9) {
        project.draftReversedFaceIds = decodeOrientations(root["draftReversedFaces"]);
        project.cadReversedFaceIds = decodeOrientations(root["cadReversedFaces"]);
    }
    validateOrientations(project);
    if (bool(project.preview) != (project.geometryRevision != 0) ||
        bool(project.mesh) != (project.meshRevision != 0))
        invalid(QStringLiteral("Netz und Revision passen nicht zusammen."));
    project.simulation = decodeSimulation(root["simulation"], version);
    if (version < 17 && (std::any_of(
                             project.simulation.settings.acousticModel.domains.begin(),
                             project.simulation.settings.acousticModel.domains.end(),
                             [](const auto &r) { return !r.twoSidedFaceIds.empty(); }) ||
                         (project.simulation.request && project.simulation.request->subdomains &&
                          !project.simulation.request->subdomains->twoSidedSurfaces.empty()) ||
                         project.simulation.settings.symmetryXYEnabled ||
                         project.simulation.settings.symmetryPlaneTolerance != 1e-6 ||
                         (project.simulation.request && project.simulation.request->symmetry)))
        invalid(QStringLiteral("XY-Symmetrie und beidseitige Flächen benötigen Projektformat 17."));
    if (project.gmsh && project.mesh && project.meshRevision == project.revision) {
        const auto orientedSource =
            applyFaceOrientation(project.cad->mesh, project.cadReversedFaceIds);
        if (orientedSource != project.mesh->mesh)
            invalid(QStringLiteral("Geprüftes Gmsh-Netz stimmt nicht mit dem Importnetz überein."));
        MeshingOptions options;
        options.targetEdgeLength = project.meshing.targetEdgeLength;
        options.adaptive = project.meshing.adaptive;
        options.quality = project.meshing.quality;
        options.faceTargetEdgeLengths = project.cadFaceTargetEdgeLengths;
        const auto &settings = project.simulation.settings;
        const auto symmetry = settings.symmetryXYEnabled
                                  ? std::optional{inferXYSymmetry(project.mesh->mesh,
                                                                  settings.symmetryPlaneTolerance)}
                                  : std::nullopt;
        (void)validateGmshMesh(project.mesh->mesh, settings.acousticModel, options, {},
                               symmetry ? &*symmetry : nullptr);
    }
    const auto &a = project.simulation;
    if (a.revision > project.revision)
        invalid(QStringLiteral("Simulationsrevision liegt in der Zukunft."));
    const bool reduced = a.request && a.request->axisymmetric;
    if (a.current &&
        (project.geometryRevision != project.revision || a.revision != project.revision ||
         (!reduced &&
          (!project.mesh || project.meshRevision != project.revision ||
           !sameMesh(project.mesh->mesh, *(a.request->subdomains ? a.sourceMesh : a.mesh))))))
        invalid(QStringLiteral("Aktuelles Ergebnis gehört nicht zum aktuellen Netz."));
    const auto &activeOrientations = project.source == ProjectSource::draft
                                         ? project.draftReversedFaceIds
                                         : project.cadReversedFaceIds;
    if (!activeOrientations.empty() &&
        !(project.gmsh && !project.simulation.settings.acousticModel.domains.empty()) &&
        ((project.mesh && project.meshRevision == project.revision) || (a.current && reduced))) {
        if (!project.preview || project.geometryRevision != project.revision ||
            !inspectMesh(*project.preview).watertight())
            invalid(QStringLiteral("Orientierungskonflikt einer aktuellen Randfläche."));
        galerkin::validateExteriorComponents(*project.preview);
    }
    if (project.source == ProjectSource::freecad) {
        if ((project.preview && project.geometryRevision == project.revision) || a.current) {
            if (!project.cad)
                invalid(QStringLiteral("Aktuelle CAD-Geometrie fehlt."));
            if (project.preview && project.geometryRevision == project.revision &&
                !sameMesh(*project.preview,
                          applyFaceOrientation(project.cad->mesh, project.cadReversedFaceIds)))
                invalid(QStringLiteral("Aktuelle CAD-Vorschau passt nicht zum Snapshot."));
        }
        if (project.cad && project.mesh && project.meshRevision == project.revision) {
            const auto cadFaces = meshFaces(project.cad->mesh);
            const auto solvedFaces = meshFaces(project.mesh->mesh);
            if (cadFaces != solvedFaces)
                invalid(QStringLiteral("CAD-Flächenzuordnung passt nicht zum Netz."));
        }
        if (a.current && project.cad->hasUnresolvedAssignments())
            invalid(
                QStringLiteral("CAD-Flächen eines aktuellen Ergebnisses sind nicht bestätigt."));
        if (a.current && reduced) {
            const auto sourceAxis = defaultAcousticFrame({}, project.cad->acousticAxis).forward;
            const auto resultAxis = a.request->axisymmetric->axis;
            if (meshFaces(project.cad->mesh) != meshFaces(*a.mesh) ||
                a.request->axisymmetric->origin != project.cad->acousticOrigin ||
                std::hypot(sourceAxis.x - resultAxis.x, sourceAxis.y - resultAxis.y,
                           sourceAxis.z - resultAxis.z) > 1e-12)
                invalid(QStringLiteral(
                    "2D-Integrationsnetz oder Rotationsachse passt nicht zum CAD-Snapshot."));
        }
    }
    if (version >= 15 && allowComparisons) {
        std::set<QString> names;
        for (const auto &value : array(root["comparisonSnapshots"], 64)) {
            const auto entry = object(value);
            keys(entry, {"name", "project"});
            const auto rawName = text(entry["name"]);
            const auto name = rawName.trimmed();
            if (name.isEmpty() || name != rawName || name.size() > 128 ||
                !names.insert(name).second)
                invalid(QStringLiteral("Vergleichsname ist leer, zu lang oder doppelt."));
            auto snapshot = decode(object(entry["project"]), false);
            if (!snapshot.simulation.result || !snapshot.simulation.request ||
                !snapshot.simulation.mesh)
                invalid(QStringLiteral("Vergleichsstand benötigt ein berechnetes Ergebnis."));
            project.comparisonSnapshots.push_back(
                {name, std::make_shared<const Project>(std::move(snapshot))});
        }
    }
    return project;
}

QJsonObject encode(const Project &p) {
    auto root = encodeDraft(p);
    root.insert("format", "HornSculpt");
    root.insert("version", 18);
    root.insert("assembly", encodeAssembly(p.assembly));
    root.insert("units", "SI");
    switch (p.source) {
    case ProjectSource::draft:
        root.insert("source", "draft");
        break;
    case ProjectSource::freecad:
        root.insert("source", "freecad");
        break;
    default:
        invalid(QStringLiteral("Unbekannte Geometriequelle."));
    }
    root.insert("revision", encodeWide(p.revision));
    root.insert("geometryRevision", encodeWide(p.geometryRevision));
    root.insert("meshRevision", encodeWide(p.meshRevision));
    root.insert("preview",
                p.preview ? QJsonValue(encodeMesh(*p.preview)) : QJsonValue(QJsonValue::Null));
    root.insert("mesh",
                p.mesh ? QJsonValue(encodeMesh(p.mesh->mesh)) : QJsonValue(QJsonValue::Null));
    root.insert("cad", p.cad ? QJsonValue(encodeCad(*p.cad)) : QJsonValue(QJsonValue::Null));
    if (p.gmsh)
        root.insert("gmshMetadata", encodeGmshMetadata(*p.gmsh));
    const auto &a = p.meshing.adaptive;
    root.insert(
        "meshing",
        QJsonObject{{"targetEdgeLength", p.meshing.targetEdgeLength},
                    {"adaptiveSizing", QJsonObject{{"enabled", a.enabled},
                                                   {"minimumEdgeLength", a.minimumEdgeLength},
                                                   {"geometricTolerance", a.geometricTolerance},
                                                   {"growthRate", a.growthRate}}},
                    {"quality",
                     QJsonObject{{"enabled", p.meshing.quality.enabled},
                                 {"minimumAngleDegrees", p.meshing.quality.minimumAngleDegrees},
                                 {"maximumAspectRatio", p.meshing.quality.maximumAspectRatio}}}});
    root.insert("draftFaceSizes", encodeFaceSizes(p.draftFaceTargetEdgeLengths));
    root.insert("cadFaceSizes", encodeFaceSizes(p.cadFaceTargetEdgeLengths));
    root.insert("draftFaceColors", encodeFaceColors(p.draftFaceColors));
    root.insert("cadFaceColors", encodeFaceColors(p.cadFaceColors));
    root.insert("draftReversedFaces", encodeOrientations(p.draftReversedFaceIds));
    root.insert("cadReversedFaces", encodeOrientations(p.cadReversedFaceIds));
    root.insert("simulation", encodeSimulation(p.simulation));
    QJsonArray snapshots;
    if (p.comparisonSnapshots.size() > 64)
        invalid(QStringLiteral("Zu viele Vergleichsstände."));
    std::set<QString> names;
    for (const auto &entry : p.comparisonSnapshots) {
        const auto name = entry.name.trimmed();
        if (name.isEmpty() || name != entry.name || name.size() > 128 ||
            !names.insert(name).second || !entry.project ||
            !entry.project->comparisonSnapshots.empty() || !entry.project->simulation.result ||
            !entry.project->simulation.request || !entry.project->simulation.mesh)
            invalid(QStringLiteral("Ungültiger, doppelter oder verschachtelter Vergleichsstand."));
        snapshots.append(QJsonObject{{"name", name}, {"project", encode(*entry.project)}});
    }
    root.insert("comparisonSnapshots", snapshots);
    return root;
}
} // namespace

void writeProject(const QString &path, const Project &project) {
    const auto root = encode(project);
    (void)decode(root);
    const auto bytes = QJsonDocument(root).toJson(QJsonDocument::Compact);
    if (bytes.size() > maximumFileBytes)
        project_json::invalid(QStringLiteral("Datei ist größer als 256 MiB."));
    QSaveFile file(path);
    if (!file.open(QIODevice::WriteOnly) || file.write(bytes) != bytes.size() || !file.commit())
        throw std::runtime_error(file.errorString().toStdString());
}

Project readProject(const QString &path) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly))
        throw std::runtime_error(file.errorString().toStdString());
    if (file.size() > maximumFileBytes)
        project_json::invalid(QStringLiteral("Datei ist größer als 256 MiB."));
    const auto bytes = file.read(maximumFileBytes + 1);
    if (bytes.size() > maximumFileBytes)
        project_json::invalid(QStringLiteral("Datei ist größer als 256 MiB."));
    if (file.error() != QFileDevice::NoError)
        throw std::runtime_error(file.errorString().toStdString());
    QJsonParseError error;
    const auto document = QJsonDocument::fromJson(bytes, &error);
    if (error.error != QJsonParseError::NoError || !document.isObject())
        project_json::invalid(QStringLiteral("JSON-Objekt erwartet."));
    uniqueKeys(bytes);
    return decode(document.object());
}

FaceColors retainedCadFaceColors(const FaceColors &colors, const CadSnapshot *previous,
                                 const CadSnapshot *current) {
    if (!previous || !current || previous->documentId != current->documentId)
        return {};
    std::map<std::uint32_t, const CadFace *> before, after;
    for (const auto &face : previous->faces)
        before.emplace(face.id, &face);
    for (const auto &face : current->faces)
        after.emplace(face.id, &face);
    const bool sameGeometry =
        previous->sessionId == current->sessionId && previous->revision == current->revision;
    auto retained = colors;
    std::erase_if(retained, [&](const auto &entry) {
        const auto old = before.find(entry.first), next = after.find(entry.first);
        return old == before.end() || next == after.end() ||
               old->second->objectId != next->second->objectId ||
               next->second->state == CadFaceState::needsReview ||
               (!sameGeometry && next->second->state == CadFaceState::newFace);
    });
    return retained;
}

AcousticModelDefinition reconciledCadAcousticModel(AcousticModelDefinition model,
                                                   const CadSnapshot *previous,
                                                   const CadSnapshot *current) {
    const bool sameDocument = previous && current && previous->documentId == current->documentId;
    const bool sameGeometry = sameDocument && previous->sessionId == current->sessionId &&
                              previous->revision == current->revision;
    const auto needsReview = [&](std::uint32_t id) {
        if (!sameDocument)
            return true;
        const auto find = [id](const auto &faces) {
            return std::find_if(faces.begin(), faces.end(),
                                [id](const auto &face) { return face.id == id; });
        };
        const auto before = find(previous->faces), after = find(current->faces);
        return before == previous->faces.end() || after == current->faces.end() ||
               before->objectId != after->objectId || after->state == CadFaceState::needsReview ||
               (!sameGeometry && after->state == CadFaceState::newFace);
    };
    for (auto &domain : model.domains)
        domain.needsAssignmentReview |=
            std::any_of(domain.faceIds.begin(), domain.faceIds.end(), needsReview) ||
            std::any_of(domain.twoSidedFaceIds.begin(), domain.twoSidedFaceIds.end(), needsReview);
    for (auto &interface : model.interfaces)
        interface.needsAssignmentReview |=
            std::any_of(interface.faceIds.begin(), interface.faceIds.end(), needsReview);
    return model;
}
lem::Network reconciledCadBemAssignments(lem::Network network, const CadSnapshot *previous,
                                         const CadSnapshot *current) {
    const bool sameDocument = previous && current && previous->documentId == current->documentId;
    const bool sameGeometry = sameDocument && previous->sessionId == current->sessionId &&
                              previous->revision == current->revision;
    for (auto &element : network.elements)
        if (auto *port = std::get_if<lem::BemPort>(&element.component)) {
            const auto count = port->faceIds.size();
            std::erase_if(port->faceIds, [&](const auto id) {
                if (!sameDocument)
                    return true;
                const auto find = [id](const auto &faces) {
                    return std::find_if(faces.begin(), faces.end(),
                                        [id](const auto &face) { return face.id == id; });
                };
                const auto before = find(previous->faces), after = find(current->faces);
                const auto boundary = port->kind == lem::BemPortKind::interfacePort
                                          ? BoundaryCondition::wall
                                          : BoundaryCondition::driven;
                return before == previous->faces.end() || after == current->faces.end() ||
                       before->objectId != after->objectId || after->boundary != boundary ||
                       after->state == CadFaceState::needsReview ||
                       (!sameGeometry && after->state == CadFaceState::newFace);
            });
            port->needsAssignmentReview |= count != port->faceIds.size();
            std::erase_if(port->velocityWeights, [&](const auto &entry) {
                return std::find(port->faceIds.begin(), port->faceIds.end(), entry.first) ==
                       port->faceIds.end();
            });
        }
    return network;
}
} // namespace hornsculpt
