#include "axisymmetric_reduction.hpp"
#include "baffle_exterior.hpp"
#include "planning_frequencies.hpp"
#include "project_json.hpp"
#include "simulation_request.hpp"
#include <algorithm>
#include <cmath>
#include <hornsculpt/acoustic_dof_map.hpp>
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
std::vector<AcousticRegionDefinition> decodeAcousticRegions(const QJsonValue &value) {
    std::vector<AcousticRegionDefinition> result;
    std::set<std::uint32_t> ids, faces;
    std::map<std::uint32_t, unsigned> interfaces, twoSidedSurfaces;
    for (const auto &item : array(value, 10000)) {
        const auto o = object(item);
        auto required = o;
        required.remove("interfaceFaceIds");
        required.remove("needsAssignmentReview");
        required.remove("twoSidedFaceIds");
        keys(required, {"id", "name", "faceIds"});
        AcousticRegionDefinition region;
        region.id = integer(o["id"], std::numeric_limits<std::uint32_t>::max());
        region.name = text(o["name"]).toStdString();
        if (o.contains("needsAssignmentReview"))
            region.needsAssignmentReview = boolean(o["needsAssignmentReview"]);
        require(region.id != 0 && !region.name.empty() && ids.insert(region.id).second);
        for (const auto &f : array(o["faceIds"], 1000000)) {
            auto id = integer(f, std::numeric_limits<std::uint32_t>::max());
            require(id != 0 && faces.insert(id).second);
            region.faceIds.push_back(id);
        }
        if (o.contains("interfaceFaceIds")) {
            std::set<std::uint32_t> local;
            for (const auto &f : array(o["interfaceFaceIds"], 1000000)) {
                auto id = integer(f, std::numeric_limits<std::uint32_t>::max());
                require(id != 0 && local.insert(id).second && ++interfaces[id] <= 2);
                region.interfaceFaceIds.push_back(id);
            }
        }
        if (o.contains("twoSidedFaceIds")) {
            std::set<std::uint32_t> local;
            for (const auto &f : array(o["twoSidedFaceIds"], 1000000)) {
                const auto id = integer(f, std::numeric_limits<std::uint32_t>::max());
                require(id != 0 && local.insert(id).second && ++twoSidedSurfaces[id] <= 2);
                region.twoSidedFaceIds.push_back(id);
            }
        }
        require(!region.faceIds.empty() || !region.interfaceFaceIds.empty() ||
                !region.twoSidedFaceIds.empty());
        result.push_back(std::move(region));
    }
    for (const auto &[id, count] : interfaces) {
        (void)count;
        require(!faces.contains(id));
    }
    for (const auto &[id, count] : twoSidedSurfaces) {
        (void)count;
        require(!faces.contains(id) && !interfaces.contains(id));
    }
    return result;
}
AcousticModelDefinition
migrateAcousticRegions(const std::vector<AcousticRegionDefinition> &regions) {
    AcousticModelDefinition model;
    if (regions.empty())
        return model;
    std::set<std::uint32_t> used;
    for (const auto &r : regions) {
        used.insert(r.id);
        model.domains.push_back({r.id, r.name, BemDomainKind::interior, r.faceIds,
                                 r.twoSidedFaceIds, r.needsAssignmentReview});
    }
    std::uint32_t exterior = 1;
    while (used.contains(exterior)) {
        require(exterior < std::numeric_limits<std::uint32_t>::max());
        ++exterior;
    }
    model.domains.push_back({exterior, "Außenraum", BemDomainKind::exterior, {}, {}, false});
    std::map<std::uint32_t, std::vector<std::size_t>> owners;
    for (std::size_t i = 0; i < regions.size(); ++i)
        for (auto face : regions[i].interfaceFaceIds)
            owners[face].push_back(i);
    std::map<std::pair<std::size_t, std::size_t>, std::vector<std::uint32_t>> grouped;
    for (const auto &[face, indices] : owners) {
        require(!indices.empty() && indices.size() <= 2);
        grouped[{indices[0], indices.size() == 2 ? indices[1] : regions.size()}].push_back(face);
    }
    for (const auto &[pair, faces] : grouped) {
        const auto id = static_cast<std::uint32_t>(model.interfaces.size() + 1);
        model.interfaces.push_back(
            {id, regions[pair.first].name + " / Interface", regions[pair.first].id,
             pair.second == regions.size() ? exterior : regions[pair.second].id, faces,
             regions[pair.first].needsAssignmentReview ||
                 (pair.second != regions.size() && regions[pair.second].needsAssignmentReview)});
    }
    return model;
}
QJsonArray encodeElementGroups(const std::vector<AcousticElementGroup> &groups) {
    QJsonArray result;
    for (const auto &group : groups) {
        QJsonArray faces;
        for (auto id : group.faceIds)
            faces.append(double(id));
        result.append(QJsonObject{{"id", double(group.id)},
                                  {"name", QString::fromStdString(group.name)},
                                  {"faceIds", faces}});
    }
    return result;
}
std::vector<AcousticElementGroup> decodeElementGroups(const QJsonValue &value,
                                                      const std::vector<std::uint32_t> &allowed) {
    std::vector<AcousticElementGroup> result;
    std::set<std::uint32_t> ids, used;
    const std::set<std::uint32_t> source(allowed.begin(), allowed.end());
    for (const auto &item : array(value, 100000)) {
        const auto o = object(item);
        keys(o, {"id", "name", "faceIds"});
        AcousticElementGroup group;
        group.id = integer(o["id"], std::numeric_limits<std::uint32_t>::max());
        group.name = text(o["name"]).toStdString();
        require(group.id && !group.name.empty() && ids.insert(group.id).second);
        for (const auto &face : array(o["faceIds"], 1000000, 1)) {
            auto id = integer(face, std::numeric_limits<std::uint32_t>::max());
            require(source.contains(id) && used.insert(id).second);
            group.faceIds.push_back(id);
        }
        result.push_back(std::move(group));
    }
    return result;
}
QJsonObject encodeAcousticModel(const AcousticModelDefinition &model) {
    QJsonArray domains, interfaces;
    for (const auto &d : model.domains) {
        QJsonArray faces, twoSided;
        for (auto id : d.faceIds)
            faces.append(double(id));
        for (auto id : d.twoSidedFaceIds)
            twoSided.append(double(id));
        domains.append(
            QJsonObject{{"id", double(d.id)},
                        {"name", QString::fromStdString(d.name)},
                        {"kind", d.kind == BemDomainKind::interior ? "interior" : "exterior"},
                        {"faceIds", faces},
                        {"twoSidedFaceIds", twoSided},
                        {"needsAssignmentReview", d.needsAssignmentReview},
                        {"groups", encodeElementGroups(d.groups)}});
    }
    for (const auto &i : model.interfaces) {
        QJsonArray faces;
        for (auto id : i.faceIds)
            faces.append(double(id));
        interfaces.append(QJsonObject{{"id", double(i.id)},
                                      {"name", QString::fromStdString(i.name)},
                                      {"firstDomainId", double(i.firstDomainId)},
                                      {"secondDomainId", double(i.secondDomainId)},
                                      {"faceIds", faces},
                                      {"needsAssignmentReview", i.needsAssignmentReview},
                                      {"groups", encodeElementGroups(i.groups)}});
    }
    QJsonValue baffle(QJsonValue::Null);
    if (model.exteriorBaffle) {
        const auto &b = *model.exteriorBaffle;
        baffle = QJsonObject{{"origin", encodeVector(b.origin)},
                             {"normal", encodeVector(b.normal)},
                             {"planeDistanceTolerance", b.planeDistanceTolerance}};
    }
    return {{"domains", domains}, {"interfaces", interfaces}, {"exteriorBaffle", baffle}};
}
AcousticModelDefinition decodeAcousticModel(const QJsonValue &value) {
    const auto o = object(value);
    keys(o, {"domains", "interfaces", "exteriorBaffle"});
    AcousticModelDefinition result;
    std::set<std::uint32_t> domainIds, interfaceIds, physicalFaces, interfaceFaces;
    std::map<std::uint32_t, unsigned> twoSidedFaces;
    unsigned exteriorCount = 0;
    auto readFaces = [](const QJsonValue &v) {
        std::vector<std::uint32_t> faces;
        std::set<std::uint32_t> seen;
        for (const auto &item : array(v, 1000000)) {
            const auto id = integer(item, std::numeric_limits<std::uint32_t>::max());
            require(id && seen.insert(id).second);
            faces.push_back(id);
        }
        return faces;
    };
    for (const auto &item : array(o["domains"], 10000)) {
        const auto d = object(item);
        auto required = d;
        required.remove("groups");
        keys(required,
             {"id", "name", "kind", "faceIds", "twoSidedFaceIds", "needsAssignmentReview"});
        AcousticDomainDefinition domain;
        domain.id = integer(d["id"], std::numeric_limits<std::uint32_t>::max());
        domain.name = text(d["name"]).toStdString();
        const auto kind = text(d["kind"]);
        require(kind == "interior" || kind == "exterior");
        domain.kind = kind == "interior" ? BemDomainKind::interior : BemDomainKind::exterior;
        exteriorCount += domain.kind == BemDomainKind::exterior;
        require(domain.id && !domain.name.empty() && domainIds.insert(domain.id).second);
        domain.faceIds = readFaces(d["faceIds"]);
        domain.twoSidedFaceIds = readFaces(d["twoSidedFaceIds"]);
        domain.needsAssignmentReview = boolean(d["needsAssignmentReview"]);
        if (d.contains("groups")) {
            auto allowed = domain.faceIds;
            allowed.insert(allowed.end(), domain.twoSidedFaceIds.begin(),
                           domain.twoSidedFaceIds.end());
            domain.groups = decodeElementGroups(d["groups"], allowed);
        }
        for (auto id : domain.faceIds)
            require(physicalFaces.insert(id).second);
        for (auto id : domain.twoSidedFaceIds)
            require(++twoSidedFaces[id] <= 2);
        result.domains.push_back(std::move(domain));
    }
    require(result.domains.empty() ? exteriorCount == 0 : exteriorCount == 1);
    for (const auto &item : array(o["interfaces"], 10000)) {
        const auto i = object(item);
        auto required = i;
        required.remove("groups");
        keys(required,
             {"id", "name", "firstDomainId", "secondDomainId", "faceIds", "needsAssignmentReview"});
        AcousticInterfaceDefinition interface;
        interface.id = integer(i["id"], std::numeric_limits<std::uint32_t>::max());
        interface.name = text(i["name"]).toStdString();
        interface.firstDomainId =
            integer(i["firstDomainId"], std::numeric_limits<std::uint32_t>::max());
        interface.secondDomainId =
            integer(i["secondDomainId"], std::numeric_limits<std::uint32_t>::max());
        interface.faceIds = readFaces(i["faceIds"]);
        interface.needsAssignmentReview = boolean(i["needsAssignmentReview"]);
        if (i.contains("groups"))
            interface.groups = decodeElementGroups(i["groups"], interface.faceIds);
        require(interface.id && !interface.name.empty() &&
                interfaceIds.insert(interface.id).second &&
                domainIds.contains(interface.firstDomainId) &&
                domainIds.contains(interface.secondDomainId) &&
                interface.firstDomainId != interface.secondDomainId && !interface.faceIds.empty());
        for (auto id : interface.faceIds)
            require(!physicalFaces.contains(id) && interfaceFaces.insert(id).second);
        result.interfaces.push_back(std::move(interface));
    }
    if (!o["exteriorBaffle"].isNull()) {
        const auto b = object(o["exteriorBaffle"]);
        keys(b, {"origin", "normal", "planeDistanceTolerance"});
        ExteriorBaffle baffle{vector(b["origin"]), vector(b["normal"]),
                              range(b["planeDistanceTolerance"], 1e-9, .01)};
        const double length = std::hypot(baffle.normal.x, baffle.normal.y, baffle.normal.z);
        require(std::isfinite(length) && std::abs(length - 1.) < 1e-8 && exteriorCount == 1);
        result.exteriorBaffle = baffle;
    }
    for (const auto &d : result.domains)
        for (auto id : d.twoSidedFaceIds)
            require(!physicalFaces.contains(id) && !interfaceFaces.contains(id));
    return result;
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
QJsonValue encodeGraphEvidence(const std::optional<BemFastMultipoleGraphEvidence> &value) {
    if (!value)
        return QJsonValue::Null;
    const auto &e = *value;
    return QJsonObject{
        {"captureAttempts", encodeWide(e.captureAttempts)},
        {"instantiations", encodeWide(e.instantiations)},
        {"launches", encodeWide(e.launches)},
        {"reuses", encodeWide(e.reuses)},
        {"executedBatches", encodeWide(e.executedBatches)},
        {"canonicalChunks", encodeWide(e.canonicalChunks)},
        {"canonicalBatches", encodeWide(e.canonicalBatches)},
        {"invalidatedExecutables", encodeWide(e.invalidatedExecutables)},
        {"budgetBypasses", encodeWide(e.budgetBypasses)},
        {"admissionClosures", encodeWide(e.admissionClosures)},
        {"maximumCachedExecutables", encodeWide(e.maximumCachedExecutables)},
        {"maximumCachedNodes", encodeWide(e.maximumCachedNodes)},
        {"peakOwnedHostMetadataBytes", encodeWide(e.peakOwnedHostMetadataBytes)},
        {"captureMilliseconds", e.captureMilliseconds},
        {"instantiationMilliseconds", e.instantiationMilliseconds},
        {"memoryObservationAvailable", e.memoryObservationAvailable},
        {"minimumObservedFreeDeviceBytes", encodeWide(e.minimumObservedFreeDeviceBytes)},
        {"totalDeviceBytes", encodeWide(e.totalDeviceBytes)},
        {"requiredFreeReserveBytes", encodeWide(e.requiredFreeReserveBytes)}};
}
BemFastMultipoleGraphEvidence decodeGraphEvidence(const QJsonValue &value) {
    const auto o = object(value);
    keys(o, {"captureAttempts", "instantiations", "launches", "reuses", "executedBatches",
             "canonicalChunks", "canonicalBatches", "invalidatedExecutables", "budgetBypasses",
             "admissionClosures", "maximumCachedExecutables", "maximumCachedNodes",
             "peakOwnedHostMetadataBytes", "captureMilliseconds", "instantiationMilliseconds",
             "memoryObservationAvailable", "minimumObservedFreeDeviceBytes", "totalDeviceBytes",
             "requiredFreeReserveBytes"});
    BemFastMultipoleGraphEvidence e;
    e.captureAttempts = wideInteger(o["captureAttempts"]);
    e.instantiations = wideInteger(o["instantiations"]);
    e.launches = wideInteger(o["launches"]);
    e.reuses = wideInteger(o["reuses"]);
    e.executedBatches = wideInteger(o["executedBatches"]);
    e.canonicalChunks = wideInteger(o["canonicalChunks"]);
    e.canonicalBatches = wideInteger(o["canonicalBatches"]);
    e.invalidatedExecutables = wideInteger(o["invalidatedExecutables"]);
    e.budgetBypasses = wideInteger(o["budgetBypasses"]);
    e.admissionClosures = wideInteger(o["admissionClosures"]);
    e.maximumCachedExecutables = wideInteger(o["maximumCachedExecutables"]);
    e.maximumCachedNodes = wideInteger(o["maximumCachedNodes"]);
    e.peakOwnedHostMetadataBytes = wideInteger(o["peakOwnedHostMetadataBytes"]);
    e.minimumObservedFreeDeviceBytes = wideInteger(o["minimumObservedFreeDeviceBytes"]);
    e.totalDeviceBytes = wideInteger(o["totalDeviceBytes"]);
    e.requiredFreeReserveBytes = wideInteger(o["requiredFreeReserveBytes"]);
    e.captureMilliseconds = range(o["captureMilliseconds"], 0, 1e15);
    e.instantiationMilliseconds = range(o["instantiationMilliseconds"], 0, 1e15);
    require(o["memoryObservationAvailable"].isBool());
    e.memoryObservationAvailable = o["memoryObservationAvailable"].toBool();
    require(e.instantiations <= e.captureAttempts && e.reuses <= e.launches &&
            e.launches <= e.executedBatches && e.canonicalChunks <= e.canonicalBatches &&
            e.invalidatedExecutables <= e.instantiations &&
            e.maximumCachedExecutables <= e.instantiations &&
            e.maximumCachedExecutables <= e.maximumCachedNodes);
    if (e.memoryObservationAvailable)
        require(e.totalDeviceBytes > 0 && e.minimumObservedFreeDeviceBytes <= e.totalDeviceBytes &&
                e.requiredFreeReserveBytes <= e.totalDeviceBytes);
    else
        require(e.minimumObservedFreeDeviceBytes == 0 && e.totalDeviceBytes == 0 &&
                e.requiredFreeReserveBytes == 0);
    return e;
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
    return {{"leaf_support_capacity", int(s.leafSupportCapacity)},
            {"preconditioner_block_width", int(s.preconditionerBlockWidth)},
            {"separationRatio", s.separationRatio},
            {"maximumIterations", int(s.maximumIterations)},
            {"restartLength", int(s.restartLength)}};
}
FastMultipoleSettings decodeFastMultipole(const QJsonValue &v) {
    const auto o = object(v);
    FastMultipoleSettings s;
    if (o.contains("leaf_support_capacity")) {
        auto required = o;
        required.remove("preconditioner_block_width");
        keys(required,
             {"leaf_support_capacity", "separationRatio", "maximumIterations", "restartLength"});
        require(!o.contains("leafPanelCapacity") && !o.contains("leaf_panel_capacity") &&
                !o.contains("interpolationOrder"));
        s.leafSupportCapacity = integer(o["leaf_support_capacity"], 256);
        if (o.contains("preconditioner_block_width"))
            s.preconditionerBlockWidth = integer(o["preconditioner_block_width"], 2048);
    } else {
        // Historical panel hierarchy settings migrate only the neutral leaf
        // capacity. Chebyshev p has no meaning as a spherical/angular order.
        require(o.contains("leafPanelCapacity") != o.contains("leaf_panel_capacity"));
        if (o.contains("leafPanelCapacity"))
            keys(o, {"interpolationOrder", "leafPanelCapacity", "separationRatio",
                     "maximumIterations", "restartLength"});
        else
            keys(o, {"interpolationOrder", "leaf_panel_capacity", "separationRatio",
                     "maximumIterations", "restartLength"});
        s.leafSupportCapacity = integer(
            o[o.contains("leafPanelCapacity") ? "leafPanelCapacity" : "leaf_panel_capacity"], 256);
        const auto historicalOrder = integer(o["interpolationOrder"], 8);
        require(historicalOrder >= 3);
    }
    s.separationRatio = range(o["separationRatio"], 1.5, 10);
    s.maximumIterations = integer(o["maximumIterations"], 10000);
    s.restartLength = integer(o["restartLength"], 300);
    require(s.leafSupportCapacity >= 4 && s.maximumIterations >= 1 && s.restartLength >= 2);
    require(s.preconditionerBlockWidth == 32 || s.preconditionerBlockWidth == 64 ||
            s.preconditionerBlockWidth == 128 || s.preconditionerBlockWidth == 256 ||
            s.preconditionerBlockWidth == 512 || s.preconditionerBlockWidth == 1024 ||
            s.preconditionerBlockWidth == 2048);
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
            {"frequencyHz", s.frequencyHz},
            {"excitationIndex", s.excitationIndex}};
}
FieldDisplaySettings decodeFieldDisplay(const QJsonValue &value) {
    const auto o = object(value);
    auto required = o;
    required.remove("frequencyHz");
    required.remove("excitationIndex");
    keys(required, {"vertical", "minimumLateral", "maximumLateral", "minimumAxial", "maximumAxial",
                    "lateralCount", "axialCount", "quantity", "minimumDb", "maximumDb",
                    "phaseDegrees", "exclusionDistance"});
    FieldDisplaySettings s{boolean(o["vertical"]),        number(o["minimumLateral"]),
                           number(o["maximumLateral"]),   number(o["minimumAxial"]),
                           number(o["maximumAxial"]),     integer(o["lateralCount"], 201),
                           integer(o["axialCount"], 201), int(integer(o["quantity"], 2)),
                           number(o["minimumDb"]),        number(o["maximumDb"]),
                           number(o["phaseDegrees"]),     positive(o["exclusionDistance"])};
    if (o.contains("frequencyHz"))
        s.frequencyHz = number(o["frequencyHz"]);
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
QJsonArray encodeIndices(const std::vector<std::uint32_t> &values) {
    QJsonArray result;
    for (auto value : values)
        result.append(double(value));
    return result;
}
std::vector<std::uint32_t> decodeIndices(const QJsonValue &value, std::size_t limit) {
    std::vector<std::uint32_t> result;
    require(limit > 0);
    std::set<std::uint32_t> seen;
    for (const auto &v : array(value, limit)) {
        const auto index = integer(v, static_cast<std::uint32_t>(limit - 1));
        require(seen.insert(index).second);
        result.push_back(index);
    }
    return result;
}
QJsonObject encodeSubdomains(const BemSubdomainModel &model) {
    QJsonArray domains, interfaces, twoSidedSurfaces;
    for (const auto &d : model.domains) {
        QJsonArray panels;
        for (const auto &p : d.panels)
            panels.append(
                QJsonObject{{"triangle", double(p.triangleIndex)}, {"sign", p.normalSign}});
        domains.append(QJsonObject{{"id", double(d.id)},
                                   {"name", QString::fromStdString(d.name)},
                                   {"exterior", d.kind == BemDomainKind::exterior},
                                   {"panels", panels}});
    }
    for (const auto &i : model.interfaces)
        interfaces.append(QJsonObject{{"id", double(i.id)},
                                      {"name", QString::fromStdString(i.name)},
                                      {"first", double(i.firstDomainId)},
                                      {"second", double(i.secondDomainId)},
                                      {"triangles", encodeIndices(i.triangleIndices)}});
    for (const auto &wall : model.twoSidedSurfaces)
        twoSidedSurfaces.append(QJsonObject{{"first", double(wall.firstDomainId)},
                                            {"second", double(wall.secondDomainId)},
                                            {"triangles", encodeIndices(wall.triangleIndices)}});
    return {
        {"domains", domains}, {"interfaces", interfaces}, {"twoSidedSurfaces", twoSidedSurfaces}};
}
BemSubdomainModel decodeSubdomains(const QJsonValue &value, const SurfaceMesh &mesh) {
    const auto o = object(value);
    auto required = o;
    required.remove("twoSidedSurfaces");
    keys(required, {"domains", "interfaces"});
    BemSubdomainModel result;
    std::set<std::uint32_t> domains, interfaces;
    for (const auto &v : array(o["domains"], 10000, 1)) {
        const auto d = object(v);
        keys(d, {"id", "name", "exterior", "panels"});
        BemSubdomain domain;
        domain.id = integer(d["id"], UINT32_MAX);
        domain.name = text(d["name"]).toStdString();
        require(domains.insert(domain.id).second && !domain.name.empty());
        domain.kind = boolean(d["exterior"]) ? BemDomainKind::exterior : BemDomainKind::interior;
        std::set<std::uint32_t> seen;
        for (const auto &p : array(d["panels"], mesh.triangles.size(), 1)) {
            const auto panel = object(p);
            keys(panel, {"triangle", "sign"});
            const auto index =
                integer(panel["triangle"], static_cast<std::uint32_t>(mesh.triangles.size() - 1));
            const double sign = number(panel["sign"]);
            require((sign == -1 || sign == 1) && seen.insert(index).second);
            domain.panels.push_back({index, static_cast<int>(sign)});
        }
        result.domains.push_back(std::move(domain));
    }
    for (const auto &v : array(o["interfaces"], 10000)) {
        const auto i = object(v);
        keys(i, {"id", "name", "first", "second", "triangles"});
        BemInterface face;
        face.id = integer(i["id"], UINT32_MAX);
        face.name = text(i["name"]).toStdString();
        face.firstDomainId = integer(i["first"], UINT32_MAX);
        face.secondDomainId = integer(i["second"], UINT32_MAX);
        require(interfaces.insert(face.id).second && !face.name.empty() &&
                face.firstDomainId != face.secondDomainId && domains.contains(face.firstDomainId) &&
                domains.contains(face.secondDomainId));
        face.triangleIndices = decodeIndices(i["triangles"], mesh.triangles.size());
        require(!face.triangleIndices.empty());
        result.interfaces.push_back(std::move(face));
    }
    if (o.contains("twoSidedSurfaces"))
        for (const auto &v : array(o["twoSidedSurfaces"], 10000)) {
            const auto w = object(v);
            keys(w, {"first", "second", "triangles"});
            BemTwoSidedSurface wall;
            wall.firstDomainId = integer(w["first"], UINT32_MAX);
            wall.secondDomainId = integer(w["second"], UINT32_MAX);
            require(wall.firstDomainId != wall.secondDomainId &&
                    domains.contains(wall.firstDomainId) && domains.contains(wall.secondDomainId));
            wall.triangleIndices = decodeIndices(w["triangles"], mesh.triangles.size());
            require(!wall.triangleIndices.empty());
            result.twoSidedSurfaces.push_back(std::move(wall));
        }
    return result;
}
QJsonObject encodeSymmetry(const OrthogonalSymmetryPlanes &s) {
    return {{"origin", encodeVector(s.origin)},
            {"firstNormal", encodeVector(s.firstNormal)},
            {"secondNormal", encodeVector(s.secondNormal)},
            {"planeDistanceTolerance", s.planeDistanceTolerance}};
}
OrthogonalSymmetryPlanes decodeSymmetry(const QJsonValue &value) {
    const auto o = object(value);
    keys(o, {"origin", "firstNormal", "secondNormal", "planeDistanceTolerance"});
    OrthogonalSymmetryPlanes s{vector(o["origin"]), vector(o["firstNormal"]),
                               vector(o["secondNormal"]),
                               range(o["planeDistanceTolerance"], 1e-9, .01)};
    require(s.origin == Vec3{} &&
            (s.firstNormal == Vec3{1, 0, 0} || s.firstNormal == Vec3{-1, 0, 0}) &&
            (s.secondNormal == Vec3{0, 1, 0} || s.secondNormal == Vec3{0, -1, 0}));
    return s;
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
    require(r.observerValidity == exteriorObserverValidity(r));
    if (!r.observerValidity.empty()) {
        QJsonArray validity;
        for (auto value : r.observerValidity)
            validity.append(int(value));
        result.insert("observerValidity", validity);
    }
    if (r.symmetry)
        result.insert("symmetry", encodeSymmetry(*r.symmetry));
    if (r.exteriorBaffle) {
        const auto &b = *r.exteriorBaffle;
        result.insert("exteriorBaffle",
                      QJsonObject{{"origin", encodeVector(b.origin)},
                                  {"normal", encodeVector(b.normal)},
                                  {"planeDistanceTolerance", b.planeDistanceTolerance}});
    }
    if (r.subdomains)
        result.insert("subdomains", encodeSubdomains(*r.subdomains));
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
    required.remove("subdomains");
    required.remove("symmetry");
    required.remove("exteriorBaffle");
    required.remove("observerValidity");
    keys(required, {"medium", "frequencies", "observers", "excitations", "numerics"});
    BemRequest r;
    if (o.contains("symmetry"))
        r.symmetry = decodeSymmetry(o["symmetry"]);
    if (o.contains("exteriorBaffle")) {
        const auto b = object(o["exteriorBaffle"]);
        keys(b, {"origin", "normal", "planeDistanceTolerance"});
        r.exteriorBaffle = ExteriorBaffle{vector(b["origin"]), vector(b["normal"]),
                                          range(b["planeDistanceTolerance"], 1e-9, .01)};
    }
    if (o.contains("subdomains")) {
        r.subdomains = decodeSubdomains(o["subdomains"], mesh);
        if (r.symmetry)
            validateSymmetryAcousticPartition(mesh, *r.subdomains, *r.symmetry,
                                              r.exteriorBaffle ? &*r.exteriorBaffle : nullptr);
        else
            validateAcousticPartition(mesh, *r.subdomains);
    }
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
        for (const auto &faceValue :
             array(entry["faces"], static_cast<qsizetype>(known.size()), 1)) {
            const auto f = object(faceValue);
            keys(f, {"faceId", "normalVelocity"});
            FaceVelocity velocity{integer(f["faceId"]), complex(f["normalVelocity"])};
            const auto found = known.find(velocity.faceId);
            require(found != known.end() && found->second == BoundaryCondition::driven &&
                    assigned.insert(velocity.faceId).second &&
                    !portFaces.contains(velocity.faceId));
            e.faces.push_back(velocity);
        }
        r.excitations.push_back(std::move(e));
    }
    if (o.contains("observerValidity")) {
        for (const auto &value :
             array(o["observerValidity"], r.observers.size(), r.observers.size()))
            r.observerValidity.push_back(static_cast<std::uint8_t>(integer(value, 1)));
    }
    require(r.observerValidity == exteriorObserverValidity(r));
    return r;
}
QJsonArray encodeLevelEvidence(const std::vector<BemFastMultipoleLevelEvidence> &levels) {
    QJsonArray result;
    for (const auto &level : levels)
        result.append(QJsonObject{
            {"level", int(level.level)},
            {"representation", level.representation == BemFastMultipoleRepresentation::planeWave
                                   ? "planeWave"
                                   : "regularSpherical"},
            {"momentBandlimit", int(level.momentBandlimit)},
            {"translationOrder", int(level.translationOrder)},
            {"quadratureOrder", int(level.quadratureOrder)}});
    return result;
}
std::vector<BemFastMultipoleLevelEvidence> decodeLevelEvidence(const QJsonValue &value) {
    std::vector<BemFastMultipoleLevelEvidence> result;
    for (const auto &entry : array(value, 31)) {
        const auto item = object(entry);
        keys(item,
             {"level", "representation", "momentBandlimit", "translationOrder", "quadratureOrder"});
        BemFastMultipoleLevelEvidence level;
        level.level = integer(item["level"], 30);
        const auto representation = text(item["representation"]);
        require(representation == "planeWave" || representation == "regularSpherical");
        level.representation = representation == "planeWave"
                                   ? BemFastMultipoleRepresentation::planeWave
                                   : BemFastMultipoleRepresentation::regularSpherical;
        level.momentBandlimit = integer(item["momentBandlimit"], 512);
        level.translationOrder = integer(item["translationOrder"], 1024);
        level.quadratureOrder = integer(item["quadratureOrder"], 1024);
        require(level.translationOrder >= 2 * level.momentBandlimit &&
                level.quadratureOrder >= level.translationOrder &&
                (result.empty() || level.level > result.back().level));
        result.push_back(level);
    }
    return result;
}
void encodeAccuracyEvidence(QJsonObject &item, const BemExcitationResult &e) {
    item.insert("quadratureEvidence",
                QJsonObject{{"farComparisonDegree", int(e.farQuadratureComparisonDegree)},
                            {"farVerificationDegree", int(e.farQuadratureVerificationDegree)},
                            {"farRefinements", int(e.farQuadratureRefinements)},
                            {"farDiscrepancy", e.farQuadratureDiscrepancy},
                            {"farPreparationMilliseconds", e.farQuadraturePreparationMilliseconds},
                            {"farPressureMilliseconds", e.farQuadraturePressureMilliseconds},
                            {"farFluxMilliseconds", e.farQuadratureFluxMilliseconds},
                            {"fieldComparisonDegree", int(e.fieldQuadratureComparisonDegree)},
                            {"fieldVerificationDegree", int(e.fieldQuadratureVerificationDegree)},
                            {"fieldDiscrepancy", e.fieldQuadratureDiscrepancy}});
    item.insert("fastMultipoleAccuracySeed", e.fastMultipoleAccuracySeed);
    item.insert("fastMultipoleComparisonAccuracySeed", e.fastMultipoleComparisonAccuracySeed);
    item.insert("fastMultipoleVerificationAccuracySeed", e.fastMultipoleVerificationAccuracySeed);
    item.insert("fastMultipoleAccuracyRefinements", int(e.fastMultipoleAccuracyRefinements));
    item.insert("fastMultipoleComparisonLevels",
                encodeLevelEvidence(e.fastMultipoleComparisonLevels));
    item.insert("fastMultipoleVerificationLevels",
                encodeLevelEvidence(e.fastMultipoleVerificationLevels));
    item.insert("historicalChebyshevEvidence",
                e.historicalChebyshevEvidence
                    ? QJsonValue(QJsonObject{
                          {"comparisonOrder", int(e.historicalChebyshevEvidence->comparisonOrder)},
                          {"auditOrder", int(e.historicalChebyshevEvidence->auditOrder)}})
                    : QJsonValue::Null);
}
void decodeQuadratureEvidence(BemExcitationResult &e, const QJsonObject &item, unsigned version,
                              bool multipole, bool hasFar, double tolerance) {
    // Earlier archives did not perform this independent rule check.
    if (!item.contains("quadratureEvidence")) {
        require(version < 16);
        return;
    }
    require(version >= 12);
    const auto q = object(item["quadratureEvidence"]);
    keys(q, {"farComparisonDegree", "farVerificationDegree", "farRefinements", "farDiscrepancy",
             "farPreparationMilliseconds", "farPressureMilliseconds", "farFluxMilliseconds",
             "fieldComparisonDegree", "fieldVerificationDegree", "fieldDiscrepancy"});
    e.farQuadratureComparisonDegree = integer(q["farComparisonDegree"], 14);
    e.farQuadratureVerificationDegree = integer(q["farVerificationDegree"], 14);
    e.farQuadratureRefinements = integer(q["farRefinements"], 3);
    e.farQuadratureDiscrepancy = range(q["farDiscrepancy"], 0, 1e15);
    e.farQuadraturePreparationMilliseconds = range(q["farPreparationMilliseconds"], 0, 1e15);
    e.farQuadraturePressureMilliseconds = range(q["farPressureMilliseconds"], 0, 1e15);
    e.farQuadratureFluxMilliseconds = range(q["farFluxMilliseconds"], 0, 1e15);
    e.fieldQuadratureComparisonDegree = integer(q["fieldComparisonDegree"], 14);
    e.fieldQuadratureVerificationDegree = integer(q["fieldVerificationDegree"], 14);
    e.fieldQuadratureDiscrepancy = range(q["fieldDiscrepancy"], 0, 1e15);
    const auto validPair = [](unsigned a, unsigned b) {
        return (a == 4 && b == 6) || (a == 6 && b == 8) || (a == 8 && b == 12) ||
               (a == 12 && b == 14);
    };
    if (e.farQuadratureComparisonDegree || e.farQuadratureVerificationDegree)
        require(multipole && hasFar && !e.historicalChebyshevEvidence &&
                validPair(e.farQuadratureComparisonDegree, e.farQuadratureVerificationDegree) &&
                e.farQuadratureRefinements <= (e.farQuadratureComparisonDegree == 4   ? 0u
                                               : e.farQuadratureComparisonDegree == 6 ? 1u
                                               : e.farQuadratureComparisonDegree == 8 ? 2u
                                                                                      : 3u) &&
                e.farQuadratureDiscrepancy <= .025 * tolerance);
    else
        require(e.farQuadratureRefinements == 0 && e.farQuadratureDiscrepancy == 0 &&
                e.farQuadraturePreparationMilliseconds == 0 &&
                e.farQuadraturePressureMilliseconds == 0 && e.farQuadratureFluxMilliseconds == 0);
    if (e.fieldQuadratureComparisonDegree || e.fieldQuadratureVerificationDegree)
        require(!e.historicalChebyshevEvidence &&
                validPair(e.fieldQuadratureComparisonDegree, e.fieldQuadratureVerificationDegree) &&
                e.fieldQuadratureDiscrepancy <= .025 * tolerance);
    else
        require(e.fieldQuadratureDiscrepancy == 0);
}
void decodeAccuracyEvidence(BemExcitationResult &e, const QJsonObject &item, unsigned version,
                            bool multipole, bool hasFar, double tolerance) {
    if (item.contains("residualEvaluation")) {
        const auto evaluation = text(item["residualEvaluation"]);
        require(evaluation == "directGalerkin" || evaluation == "higherOrderFastMultipole");
        e.residualEvaluation = evaluation == "directGalerkin"
                                   ? BemResidualEvaluation::directGalerkin
                                   : BemResidualEvaluation::higherOrderFastMultipole;
        e.fastMultipoleKrylovResidual = range(item["fastMultipoleKrylovResidual"], 0, 1e15);
        e.fastMultipoleOperatorDiscrepancy =
            range(item["fastMultipoleOperatorDiscrepancy"], 0, 1e15);
    } else {
        require(version < 16 && !item.contains("fastMultipoleKrylovResidual") &&
                !item.contains("fastMultipoleOperatorDiscrepancy") &&
                !item.contains("fastMultipoleComparisonOrder") &&
                !item.contains("fastMultipoleAuditOrder"));
    }
    if (item.contains("fastMultipoleKrylovTarget")) {
        require(version >= 12);
        e.fastMultipoleKrylovTarget = range(item["fastMultipoleKrylovTarget"], 0, 1);
        require(multipole || e.fastMultipoleKrylovTarget == 0);
    }
    if (version < 16 && !item.contains("fastMultipoleAccuracySeed")) {
        for (const auto *name :
             {"fastMultipoleAccuracySeed", "fastMultipoleComparisonAccuracySeed",
              "fastMultipoleVerificationAccuracySeed", "fastMultipoleAccuracyRefinements",
              "fastMultipoleComparisonLevels", "fastMultipoleVerificationLevels",
              "historicalChebyshevEvidence"})
            require(!item.contains(name));
        BemHistoricalChebyshevEvidence historical;
        if (item.contains("fastMultipoleComparisonOrder") ||
            item.contains("fastMultipoleAuditOrder")) {
            historical.comparisonOrder = integer(item["fastMultipoleComparisonOrder"], 9);
            historical.auditOrder = integer(item["fastMultipoleAuditOrder"], 10);
        }
        if (multipole)
            e.historicalChebyshevEvidence = historical;
    } else {
        require(!item.contains("fastMultipoleComparisonOrder") &&
                !item.contains("fastMultipoleAuditOrder"));
        e.fastMultipoleAccuracySeed = range(item["fastMultipoleAccuracySeed"], 0, 1);
        e.fastMultipoleComparisonAccuracySeed =
            range(item["fastMultipoleComparisonAccuracySeed"], 0, 1);
        e.fastMultipoleVerificationAccuracySeed =
            range(item["fastMultipoleVerificationAccuracySeed"], 0, 1);
        e.fastMultipoleAccuracyRefinements = integer(item["fastMultipoleAccuracyRefinements"], 32);
        e.fastMultipoleComparisonLevels =
            decodeLevelEvidence(item["fastMultipoleComparisonLevels"]);
        e.fastMultipoleVerificationLevels =
            decodeLevelEvidence(item["fastMultipoleVerificationLevels"]);
        require(item.contains("historicalChebyshevEvidence"));
        if (!item["historicalChebyshevEvidence"].isNull()) {
            const auto historical = object(item["historicalChebyshevEvidence"]);
            keys(historical, {"comparisonOrder", "auditOrder"});
            e.historicalChebyshevEvidence = BemHistoricalChebyshevEvidence{
                integer(historical["comparisonOrder"], 9), integer(historical["auditOrder"], 10)};
        }
    }
    const bool approximate =
        e.residualEvaluation == BemResidualEvaluation::higherOrderFastMultipole;
    if (e.historicalChebyshevEvidence) {
        const auto &old = *e.historicalChebyshevEvidence;
        require(multipole && e.fastMultipoleKrylovTarget == 0 && e.fastMultipoleAccuracySeed == 0 &&
                e.fastMultipoleComparisonAccuracySeed == 0 &&
                e.fastMultipoleVerificationAccuracySeed == 0 &&
                e.fastMultipoleAccuracyRefinements == 0 &&
                e.fastMultipoleComparisonLevels.empty() &&
                e.fastMultipoleVerificationLevels.empty());
        require(approximate
                    ? (old.comparisonOrder >= 3 && old.auditOrder == old.comparisonOrder + 1)
                    : (old.comparisonOrder == 0 && old.auditOrder == 0));
    } else if (multipole) {
        require(e.fastMultipoleAccuracySeed > 0 && e.fastMultipoleAccuracySeed < 1);
        if (approximate) {
            require(e.fastMultipoleComparisonAccuracySeed > 0 &&
                    e.fastMultipoleComparisonAccuracySeed <= e.fastMultipoleAccuracySeed &&
                    e.fastMultipoleVerificationAccuracySeed > 0 &&
                    e.fastMultipoleVerificationAccuracySeed <
                        e.fastMultipoleComparisonAccuracySeed &&
                    e.fastMultipoleComparisonLevels.size() ==
                        e.fastMultipoleVerificationLevels.size() &&
                    (!hasFar || !e.fastMultipoleComparisonLevels.empty()));
            for (std::size_t i = 0; i < e.fastMultipoleComparisonLevels.size(); ++i) {
                const auto &a = e.fastMultipoleComparisonLevels[i],
                           &b = e.fastMultipoleVerificationLevels[i];
                require(a.level == b.level && b.momentBandlimit >= a.momentBandlimit &&
                        b.translationOrder >= a.translationOrder &&
                        b.quadratureOrder >= a.quadratureOrder);
            }
        } else
            require(e.fastMultipoleComparisonAccuracySeed == 0 &&
                    e.fastMultipoleVerificationAccuracySeed == 0 &&
                    e.fastMultipoleComparisonLevels.empty() &&
                    e.fastMultipoleVerificationLevels.empty());
    } else
        require(e.fastMultipoleAccuracySeed == 0 && e.fastMultipoleComparisonAccuracySeed == 0 &&
                e.fastMultipoleVerificationAccuracySeed == 0 &&
                e.fastMultipoleAccuracyRefinements == 0 &&
                e.fastMultipoleComparisonLevels.empty() &&
                e.fastMultipoleVerificationLevels.empty());
    if (approximate)
        require(multipole && e.fastMultipoleOperatorDiscrepancy <=
                                 fastMultipoleOrderComparisonBudgetFraction * tolerance);
    else
        require((multipole || e.fastMultipoleKrylovResidual == 0) &&
                e.fastMultipoleOperatorDiscrepancy == 0);
}
QJsonArray encodeSubdomainTraces(const std::vector<BemSubdomainTrace> &traces) {
    QJsonArray result;
    for (const auto &trace : traces) {
        require(trace.vertexIndices.size() == trace.pressureTrace.size() &&
                trace.triangleIndices.size() == trace.normalDerivative.size());
        QJsonArray fluxes;
        for (const auto &flux : trace.interfaceFluxes) {
            require(flux.vertexIndices.size() == flux.normalDerivative.size());
            fluxes.append(QJsonObject{{"interface", double(flux.interfaceId)},
                                      {"vertices", encodeIndices(flux.vertexIndices)},
                                      {"normalDerivative", encodePressure(flux.normalDerivative)}});
        }
        result.append(QJsonObject{{"domain", double(trace.domainId)},
                                  {"vertices", encodeIndices(trace.vertexIndices)},
                                  {"pressure", encodePressure(trace.pressureTrace)},
                                  {"triangles", encodeIndices(trace.triangleIndices)},
                                  {"normalDerivative", encodePressure(trace.normalDerivative)},
                                  {"interfaceFluxes", fluxes}});
    }
    return result;
}
std::vector<BemSubdomainTrace> decodeSubdomainTraces(const QJsonValue &value,
                                                     const SurfaceMesh &mesh,
                                                     const BemSubdomainModel &model) {
    std::vector<BemSubdomainTrace> result;
    std::set<std::uint32_t> domains;
    for (const auto &v : array(value, model.domains.size(), model.domains.size())) {
        const auto t = object(v);
        keys(t, {"domain", "vertices", "pressure", "triangles", "normalDerivative",
                 "interfaceFluxes"});
        BemSubdomainTrace trace;
        trace.domainId = integer(t["domain"], UINT32_MAX);
        require(domains.insert(trace.domainId).second &&
                std::any_of(model.domains.begin(), model.domains.end(),
                            [&](const auto &d) { return d.id == trace.domainId; }));
        trace.vertexIndices = decodeIndices(t["vertices"], mesh.vertices.size());
        trace.pressureTrace = decodePressure(t["pressure"], trace.vertexIndices.size());
        trace.triangleIndices = decodeIndices(t["triangles"], mesh.triangles.size());
        trace.normalDerivative =
            decodePressure(t["normalDerivative"], trace.triangleIndices.size());
        std::set<std::uint32_t> interfaces;
        for (const auto &f : array(t["interfaceFluxes"], model.interfaces.size())) {
            const auto o = object(f);
            keys(o, {"interface", "vertices", "normalDerivative"});
            BemSubdomainTrace::InterfaceFlux flux;
            flux.interfaceId = integer(o["interface"], UINT32_MAX);
            require(
                interfaces.insert(flux.interfaceId).second &&
                std::any_of(model.interfaces.begin(), model.interfaces.end(), [&](const auto &i) {
                    return i.id == flux.interfaceId && (i.firstDomainId == trace.domainId ||
                                                        i.secondDomainId == trace.domainId);
                }));
            flux.vertexIndices = decodeIndices(o["vertices"], mesh.vertices.size());
            flux.normalDerivative =
                decodePressure(o["normalDerivative"], flux.vertexIndices.size());
            trace.interfaceFluxes.push_back(std::move(flux));
        }
        const auto domain = std::find_if(model.domains.begin(), model.domains.end(),
                                         [&](const auto &d) { return d.id == trace.domainId; });
        std::set<std::uint32_t> expectedTriangles, expectedVertices;
        for (const auto &panel : domain->panels) {
            expectedTriangles.insert(panel.triangleIndex);
            for (auto vertex : mesh.triangles[panel.triangleIndex].vertices)
                expectedVertices.insert(vertex);
        }
        require(std::set<std::uint32_t>(trace.triangleIndices.begin(),
                                        trace.triangleIndices.end()) == expectedTriangles &&
                std::set<std::uint32_t>(trace.vertexIndices.begin(), trace.vertexIndices.end()) ==
                    expectedVertices);
        for (const auto &interface : model.interfaces) {
            if (interface.firstDomainId != trace.domainId &&
                interface.secondDomainId != trace.domainId)
                continue;
            const auto flux =
                std::find_if(trace.interfaceFluxes.begin(), trace.interfaceFluxes.end(),
                             [&](const auto &f) { return f.interfaceId == interface.id; });
            require(flux != trace.interfaceFluxes.end());
            std::set<std::uint32_t> vertices;
            for (auto triangle : interface.triangleIndices)
                for (auto vertex : mesh.triangles[triangle].vertices)
                    vertices.insert(vertex);
            require(std::set<std::uint32_t>(flux->vertexIndices.begin(),
                                            flux->vertexIndices.end()) == vertices);
        }
        result.push_back(std::move(trace));
    }
    // Both sides share the same P1 pressure and opposite outward P1 flux.
    constexpr double relativeInterfaceRoundoff = 1e-11;
    const auto near = [](std::complex<double> a, std::complex<double> b) {
        const double scale = std::max({std::abs(a), std::abs(b), 1e-300});
        return std::abs(a - b) <= relativeInterfaceRoundoff * scale;
    };
    for (const auto &interface : model.interfaces) {
        const auto findTrace = [&](std::uint32_t id) -> const BemSubdomainTrace & {
            return *std::find_if(result.begin(), result.end(),
                                 [id](const auto &t) { return t.domainId == id; });
        };
        const auto &a = findTrace(interface.firstDomainId),
                   &b = findTrace(interface.secondDomainId);
        const auto findFlux =
            [&](const BemSubdomainTrace &trace) -> const BemSubdomainTrace::InterfaceFlux & {
            return *std::find_if(trace.interfaceFluxes.begin(), trace.interfaceFluxes.end(),
                                 [&](const auto &f) { return f.interfaceId == interface.id; });
        };
        const auto &fa = findFlux(a), &fb = findFlux(b);
        std::map<std::uint32_t, std::complex<double>> pa, pb, qa, qb;
        for (std::size_t i = 0; i < a.vertexIndices.size(); ++i)
            pa.emplace(a.vertexIndices[i], a.pressureTrace[i]);
        for (std::size_t i = 0; i < b.vertexIndices.size(); ++i)
            pb.emplace(b.vertexIndices[i], b.pressureTrace[i]);
        for (std::size_t i = 0; i < fa.vertexIndices.size(); ++i)
            qa.emplace(fa.vertexIndices[i], fa.normalDerivative[i]);
        for (std::size_t i = 0; i < fb.vertexIndices.size(); ++i)
            qb.emplace(fb.vertexIndices[i], fb.normalDerivative[i]);
        for (const auto &[vertex, flux] : qa)
            require(pa.contains(vertex) && pb.contains(vertex) && qb.contains(vertex) &&
                    near(pa.at(vertex), pb.at(vertex)) && near(flux, -qb.at(vertex)));
    }
    return result;
}
QJsonObject encodeResult(const BemResult &r) {
    QJsonArray frequencies;
    for (const auto &f : r.frequencies) {
        QJsonArray excitations;
        for (const auto &e : f.excitations) {
            require(bool(e.boundary()) != !e.subdomainTraces().empty());
            QJsonArray impedances;
            for (const auto &z : e.impedances)
                impedances.append(
                    QJsonObject{{"faceId", double(z.faceId)},
                                {"area", z.area},
                                {"meanPressure", encodeComplex(z.meanPressure)},
                                {"volumeVelocity", encodeComplex(z.volumeVelocity)},
                                {"volumeImpedance", encodeComplex(z.volumeImpedance)}});
            QJsonObject excitation{
                {"name", QString::fromStdString(e.name)},
                {"observerPressure", encodePressure(e.observerPressure)},
                {"impedances", impedances},
                {"refinementSteps", int(e.refinementSteps)},
                {"iterations", int(e.iterations)},
                {"residualEvaluation", e.residualEvaluation == BemResidualEvaluation::directGalerkin
                                           ? "directGalerkin"
                                           : "higherOrderFastMultipole"},
                {"fastMultipoleKrylovResidual", e.fastMultipoleKrylovResidual},
                {"fastMultipoleKrylovTarget", e.fastMultipoleKrylovTarget},
                {"fastMultipoleOperatorDiscrepancy", e.fastMultipoleOperatorDiscrepancy},
                {"boundaryResidual", e.boundaryResidual}};
            if (e.boundary())
                excitation.insert("surfacePressure", encodePressure(e.boundary()->surfacePressure));
            if (!e.subdomainTraces().empty())
                excitation.insert("subdomainTraces", encodeSubdomainTraces(e.subdomainTraces()));
            encodeAccuracyEvidence(excitation, e);
            if (e.lemSolution) {
                excitation.insert("lemSolution", encodeLemSolution(*e.lemSolution));
                excitation.insert("solvedVelocities", encodeVelocities(e.solvedVelocities));
            }
            if (e.boundary() && !e.boundary()->exteriorSurfacePressure.empty()) {
                excitation.insert("exteriorSurfacePressure",
                                  encodePressure(e.boundary()->exteriorSurfacePressure));
                excitation.insert("solvedVelocities", encodeVelocities(e.solvedVelocities));
            }
            if (e.boundary())
                require(e.boundary()->pressureTrace.empty() ==
                        e.boundary()->normalDerivative.empty());
            if (e.boundary() && !e.boundary()->pressureTrace.empty()) {
                excitation.insert("pressureTrace", encodePressure(e.boundary()->pressureTrace));
                excitation.insert("normalDerivative",
                                  encodePressure(e.boundary()->normalDerivative));
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
             QJsonObject{
                 {"nodeCount", encodeWide(f.fastMultipole.nodeCount)},
                 {"leafCount", encodeWide(f.fastMultipole.leafCount)},
                 {"farInteractionCount", encodeWide(f.fastMultipole.farInteractionCount)},
                 {"pressurePressureNearPairs",
                  encodeWide(f.fastMultipole.pressurePressureNearPairs)},
                 {"pressurePressureFarPairs", encodeWide(f.fastMultipole.pressurePressureFarPairs)},
                 {"pressureFluxNearPairs", encodeWide(f.fastMultipole.pressureFluxNearPairs)},
                 {"pressureFluxFarPairs", encodeWide(f.fastMultipole.pressureFluxFarPairs)},
                 {"graphExecution", encodeGraphEvidence(f.fastMultipole.graphExecution)},
                 {"historicalChebyshevStatistics",
                  f.fastMultipole.historicalChebyshevStatistics
                      ? QJsonValue(
                            QJsonObject{{"directPanelPairCount",
                                         encodeWide(f.fastMultipole.historicalChebyshevStatistics
                                                        ->directPanelPairCount)},
                                        {"representedFarPanelPairs",
                                         encodeWide(f.fastMultipole.historicalChebyshevStatistics
                                                        ->representedFarPanelPairs)}})
                      : QJsonValue::Null}}}});
    }
    return {{"device", QString::fromStdString(r.device)},
            {"geometryMode", encodeGeometryMode(r.geometryMode)},
            {"acceleration", encodeAcceleration(r.acceleration)},
            {"triangleCount", double(r.triangleCount)},
            {"pressureDofCount", double(r.pressureDofCount)},
            {"peakDeviceBytes", encodeWide(r.peakDeviceBytes)},
            {"liveDeviceBytes",
             r.liveDeviceBytes ? QJsonValue(encodeWide(r.liveDeviceBytes)) : QJsonValue::Null},
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
                       const SurfaceMesh &mesh, unsigned version) {
    const auto o = object(resultValue);
    auto requiredResult = o;
    requiredResult.remove("acceleration");
    requiredResult.remove("geometryMode");
    if (version >= 16)
        require(o.contains("liveDeviceBytes"));
    requiredResult.remove("liveDeviceBytes");
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
    if (o.contains("liveDeviceBytes") && !o["liveDeviceBytes"].isNull()) {
        const auto live = wideInteger(o["liveDeviceBytes"]);
        require(live <= bytes);
        r.liveDeviceBytes = static_cast<std::size_t>(live);
    }
    const auto pressureCount =
        request.axisymmetric
            ? bem_detail::validateAxisymmetricReduction(mesh, {}, *request.axisymmetric)
                  .pressureCount
            : (request.subdomains ? makeAcousticPressureDofMap(mesh, *request.subdomains).count
                                  : mesh.vertices.size());
    require(r.triangleCount == mesh.triangles.size() && r.pressureDofCount == pressureCount);
    const auto measures = surfaceMeasures(mesh, request.axisymmetric);
    const auto exteriorMesh =
        request.infiniteBaffle ? request.infiniteBaffle->exteriorDrivenMesh : SurfaceMesh{};
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
            auto stats = object(entry["fastMultipole"]);
            if (stats.contains("graphExecution")) {
                if (!stats["graphExecution"].isNull()) {
                    require(version >= 12 && multipole);
                    f.fastMultipole.graphExecution = decodeGraphEvidence(stats["graphExecution"]);
                }
                stats.remove("graphExecution");
            }
            const auto count = [](const QJsonValue &value) {
                const auto result = wideInteger(value);
                require(result <= std::numeric_limits<std::size_t>::max());
                return static_cast<std::size_t>(result);
            };
            f.fastMultipole.nodeCount = count(stats["nodeCount"]);
            f.fastMultipole.leafCount = count(stats["leafCount"]);
            f.fastMultipole.farInteractionCount = count(stats["farInteractionCount"]);
            const auto panels = mesh.triangles.size();
            require(panels > 0 && panels <= std::numeric_limits<std::size_t>::max() / panels);
            const auto pairs = panels * panels;
            require(f.fastMultipole.leafCount <= f.fastMultipole.nodeCount &&
                    (!multipole || f.fastMultipole.leafCount >= 1));
            if (version < 16 && stats.contains("directPanelPairCount")) {
                keys(stats, {"nodeCount", "leafCount", "farInteractionCount",
                             "directPanelPairCount", "representedFarPanelPairs"});
                const auto directPairs = count(stats["directPanelPairCount"]);
                const auto farPairs = count(stats["representedFarPanelPairs"]);
                if (multipole)
                    f.fastMultipole.historicalChebyshevStatistics =
                        BemHistoricalChebyshevStatistics{directPairs, farPairs};
                else
                    require(directPairs == 0 && farPairs == 0);
            } else {
                keys(stats,
                     {"nodeCount", "leafCount", "farInteractionCount", "pressurePressureNearPairs",
                      "pressurePressureFarPairs", "pressureFluxNearPairs", "pressureFluxFarPairs",
                      "historicalChebyshevStatistics"});
                f.fastMultipole.pressurePressureNearPairs =
                    wideInteger(stats["pressurePressureNearPairs"]);
                f.fastMultipole.pressurePressureFarPairs =
                    wideInteger(stats["pressurePressureFarPairs"]);
                f.fastMultipole.pressureFluxNearPairs = wideInteger(stats["pressureFluxNearPairs"]);
                f.fastMultipole.pressureFluxFarPairs = wideInteger(stats["pressureFluxFarPairs"]);
                if (!stats["historicalChebyshevStatistics"].isNull()) {
                    const auto old = object(stats["historicalChebyshevStatistics"]);
                    keys(old, {"directPanelPairCount", "representedFarPanelPairs"});
                    f.fastMultipole.historicalChebyshevStatistics =
                        BemHistoricalChebyshevStatistics{count(old["directPanelPairCount"]),
                                                         count(old["representedFarPanelPairs"])};
                }
            }
            const auto &s = f.fastMultipole;
            if (s.graphExecution && s.farInteractionCount == 0)
                require(*s.graphExecution == BemFastMultipoleGraphEvidence{});
            if (s.historicalChebyshevStatistics) {
                require(!s.graphExecution);
                require(multipole && s.pressurePressureNearPairs == 0 &&
                        s.pressurePressureFarPairs == 0 && s.pressureFluxNearPairs == 0 &&
                        s.pressureFluxFarPairs == 0);
                const auto &old = *s.historicalChebyshevStatistics;
                require(old.representedFarPanelPairs <= pairs &&
                        old.directPanelPairCount == pairs - old.representedFarPanelPairs);
            } else if (multipole) {
                const auto nodes = std::uint64_t(mesh.vertices.size()),
                           pressurePairs = nodes * nodes, fluxPairs = nodes * panels;
                require(s.pressurePressureFarPairs <= pressurePairs &&
                        s.pressurePressureNearPairs == pressurePairs - s.pressurePressureFarPairs &&
                        s.pressureFluxFarPairs <= fluxPairs &&
                        s.pressureFluxNearPairs == fluxPairs - s.pressureFluxFarPairs);
            }
        }
        const auto excitationCount = request.coupling ? 1 : request.excitations.size();
        const auto es = array(entry["excitations"], excitationCount, excitationCount);
        for (qsizetype ei = 0; ei < es.size(); ++ei) {
            const auto item = object(es[ei]);
            if (version >= 16)
                require(item.contains("residualEvaluation"));
            if (multipole)
                require(item.contains("iterations"));
            auto requiredItem = item;
            requiredItem.remove("subdomainTraces");
            if (request.subdomains)
                requiredItem.remove("surfacePressure");
            requiredItem.remove("iterations");
            requiredItem.remove("pressureTrace");
            requiredItem.remove("normalDerivative");
            requiredItem.remove("exteriorSurfacePressure");
            if (!exteriorMesh.triangles.empty())
                requiredItem.remove("solvedVelocities");
            require(item.contains("pressureTrace") == item.contains("normalDerivative"));
            for (const auto *field :
                 {"residualEvaluation", "fastMultipoleKrylovResidual", "fastMultipoleKrylovTarget",
                  "fastMultipoleOperatorDiscrepancy", "fastMultipoleComparisonOrder",
                  "fastMultipoleAuditOrder", "fastMultipoleAccuracySeed",
                  "fastMultipoleComparisonAccuracySeed", "fastMultipoleVerificationAccuracySeed",
                  "fastMultipoleAccuracyRefinements", "fastMultipoleComparisonLevels",
                  "fastMultipoleVerificationLevels", "historicalChebyshevEvidence",
                  "quadratureEvidence"})
                requiredItem.remove(field);
            if (request.coupling) {
                require(item.contains("lemSolution") && item.contains("solvedVelocities"));
                requiredItem.remove("lemSolution");
                requiredItem.remove("solvedVelocities");
            }
            if (request.subdomains)
                keys(requiredItem, {"name", "observerPressure", "impedances", "refinementSteps",
                                    "boundaryResidual"});
            else
                keys(requiredItem, {"name", "surfacePressure", "observerPressure", "impedances",
                                    "refinementSteps", "boundaryResidual"});
            BemExcitationResult e;
            require(item.contains("subdomainTraces") == request.subdomains.has_value());
            if (request.subdomains)
                e.setSubdomainTraces(
                    decodeSubdomainTraces(item["subdomainTraces"], mesh, *request.subdomains));
            decodeAccuracyEvidence(e, item, version, multipole,
                                   f.fastMultipole.farInteractionCount > 0,
                                   request.numerics.relativeResidualTolerance);
            decodeQuadratureEvidence(e, item, version, multipole,
                                     f.fastMultipole.farInteractionCount > 0,
                                     request.numerics.relativeResidualTolerance);
            require(bool(e.historicalChebyshevEvidence) ==
                    bool(f.fastMultipole.historicalChebyshevStatistics));
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
                e.solvedVelocities =
                    decodeVelocities(item["solvedVelocities"], combinedMesh, interfaces);
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
            if (!request.subdomains)
                boundary.surfacePressure =
                    decodePressure(item["surfacePressure"], mesh.triangles.size());
            if (!exteriorMesh.triangles.empty()) {
                boundary.exteriorSurfacePressure =
                    decodePressure(item["exteriorSurfacePressure"], exteriorMesh.triangles.size());
                if (!request.coupling)
                    e.solvedVelocities =
                        decodeVelocities(item["solvedVelocities"], combinedMesh, {});
            } else
                require(!item.contains("exteriorSurfacePressure"));
            if (item.contains("pressureTrace")) {
                boundary.pressureTrace =
                    decodePressure(item["pressureTrace"], mesh.vertices.size());
                boundary.normalDerivative =
                    decodePressure(item["normalDerivative"], mesh.triangles.size());
                if (request.axisymmetric)
                    bem_detail::validateAxisymmetricTrace(
                        *request.axisymmetric, boundary.pressureTrace, boundary.normalDerivative);
                for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
                    const auto &nodes = mesh.triangles[t].vertices;
                    const auto &weights = measures.pressureIntegrals[t];
                    const bool continuous = request.axisymmetric &&
                                            request.axisymmetric->surfaceContract ==
                                                AxisymmetricSurfaceContract::continuousRevolution;
                    const auto mean =
                        continuous
                            ? (weights[0] / measures.areas[t]) * boundary.pressureTrace[nodes[0]] +
                                  (weights[1] / measures.areas[t]) *
                                      boundary.pressureTrace[nodes[1]] +
                                  (weights[2] / measures.areas[t]) *
                                      boundary.pressureTrace[nodes[2]]
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
                        absoluteMean += (weights[i] / measures.areas[t]) *
                                        std::abs(boundary.pressureTrace[nodes[i]]);
                    scale = std::max(scale, absoluteMean);
                    require(std::isfinite(scale) &&
                            std::abs(mean - boundary.surfacePressure[t]) <=
                                relativeArchivedPressureMeanTolerance * scale);
                }
            }
            e.observerPressure = decodePressure(item["observerPressure"], request.observers.size());
            e.observerValidity = request.observerValidity;
            const unsigned couplingPasses =
                e.lemSolution ? static_cast<unsigned>(e.lemSolution->ports.size()) + 1u : 1u;
            std::uint64_t auditPasses = 1;
            if (e.historicalChebyshevEvidence) {
                const auto &old = *e.historicalChebyshevEvidence;
                // Historical Chebyshev runs used p>=3. Keep a conservative
                // counter bound without retaining p as a live solver setting.
                auditPasses = old.comparisonOrder ? old.comparisonOrder - 3u : 7u;
            } else if (multipole)
                // Earlier quadrature attempts can each have all five seed stages;
                // only the final attempt's refinement count is stored.
                auditPasses = 2u * (std::uint64_t(e.fastMultipoleAccuracyRefinements) + 1 +
                                    5u * e.farQuadratureRefinements);
            const auto quadraturePasses = std::uint64_t(e.farQuadratureRefinements) + 1;
            const auto maximumCorrections =
                std::uint64_t(request.numerics.maximumRefinementSteps) * auditPasses;
            const auto bound = [](std::uint64_t value) {
                return static_cast<std::uint32_t>(std::min<std::uint64_t>(value, UINT32_MAX));
            };
            e.refinementSteps =
                integer(item["refinementSteps"], bound(maximumCorrections * couplingPasses));
            if (item.contains("iterations"))
                e.iterations = integer(
                    item["iterations"],
                    multipole
                        ? bound(std::uint64_t(request.numerics.fastMultipole.maximumIterations) *
                                (maximumCorrections + quadraturePasses) * couplingPasses)
                        : 0);
            e.boundaryResidual =
                range(item["boundaryResidual"], 0, request.numerics.relativeResidualTolerance);
            if (!request.subdomains) {
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
                    keys(z,
                         {"faceId", "area", "meanPressure", "volumeVelocity", "volumeImpedance"});
                    FaceImpedance impedance{
                        integer(z["faceId"]), positive(z["area"]), complex(z["meanPressure"]),
                        complex(z["volumeVelocity"]), complex(z["volumeImpedance"])};
                    require(assigned.contains(impedance.faceId) &&
                            seen.insert(impedance.faceId).second);
                    const auto &sum = integrals.at(impedance.faceId);
                    const double pressureScale = std::max({std::abs(impedance.meanPressure),
                                                           std::abs(sum.pressure / sum.area),
                                                           sum.absolutePressure / sum.area});
                    require(
                        std::isfinite(pressureScale) && close(impedance.area, sum.area) &&
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
                        allPressure.insert(allPressure.end(),
                                           boundary.exteriorSurfacePressure.begin(),
                                           boundary.exteriorSurfacePressure.end());
                        double area = 0, magnitude = 0;
                        std::complex<double> pressure{};
                        for (std::size_t ti = 0; ti < combinedMesh.triangles.size(); ++ti) {
                            const auto &triangle = combinedMesh.triangles[ti];
                            if (std::find(port.faceIds.begin(), port.faceIds.end(),
                                          triangle.faceId) != port.faceIds.end())
                                area += combinedAreas[ti] *
                                        (port.velocityWeights.contains(triangle.faceId)
                                             ? port.velocityWeights.at(triangle.faceId)
                                             : 1.);
                        }
                        require(std::isfinite(area) && area > 0);
                        for (std::size_t ti = 0; ti < combinedMesh.triangles.size(); ++ti)
                            if (std::find(port.faceIds.begin(), port.faceIds.end(),
                                          combinedMesh.triangles[ti].faceId) !=
                                port.faceIds.end()) {
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
                        require(
                            std::abs(e.lemSolution->ports[i].pressure - pressure) <=
                            relativePortPressureClosureTolerance *
                                std::max(magnitude, std::abs(e.lemSolution->ports[i].pressure)));
                    }
                }
            } else {
                if (!request.coupling)
                    e.solvedVelocities = request.excitations[ei].faces;
                std::map<std::uint32_t, std::complex<double>> assigned;
                for (const auto &face : e.solvedVelocities)
                    if (std::abs(face.normalVelocity) > 0)
                        assigned.emplace(face.faceId, face.normalVelocity);
                std::set<std::uint32_t> interfaceTriangles;
                for (const auto &interface : request.subdomains->interfaces)
                    interfaceTriangles.insert(interface.triangleIndices.begin(),
                                              interface.triangleIndices.end());
                std::map<std::uint32_t, FaceIntegral> integrals;
                std::map<std::uint32_t, std::complex<double>> panelMeans;
                std::map<std::pair<std::uint32_t, std::uint32_t>, std::complex<double>> sideMeans;
                std::map<std::uint32_t, std::vector<std::uint32_t>> owners;
                for (const auto &trace : e.subdomainTraces()) {
                    std::map<std::uint32_t, std::complex<double>> pressure;
                    for (std::size_t i = 0; i < trace.vertexIndices.size(); ++i)
                        pressure.emplace(trace.vertexIndices[i], trace.pressureTrace[i]);
                    for (auto index : trace.triangleIndices) {
                        if (interfaceTriangles.contains(index))
                            continue;
                        const auto &triangle = mesh.triangles[index];
                        const auto mean =
                            (pressure.at(triangle.vertices[0]) + pressure.at(triangle.vertices[1]) +
                             pressure.at(triangle.vertices[2])) /
                            3.;
                        require(
                            sideMeans.emplace(std::make_pair(index, trace.domainId), mean).second);
                        owners[index].push_back(trace.domainId);
                    }
                }
                std::map<std::uint32_t, std::pair<std::uint32_t, std::uint32_t>> twoSided;
                for (const auto &surface : request.subdomains->twoSidedSurfaces)
                    for (auto index : surface.triangleIndices)
                        require(twoSided
                                    .emplace(index, std::make_pair(surface.firstDomainId,
                                                                   surface.secondDomainId))
                                    .second);
                const double symmetryFactor = request.symmetry ? 4. : 1.;
                for (const auto &[index, domainIds] : owners) {
                    const auto &triangle = mesh.triangles[index];
                    std::complex<double> mean{};
                    if (twoSided.contains(index)) {
                        const auto [first, second] = twoSided.at(index);
                        require(domainIds.size() == 2 && sideMeans.contains({index, first}) &&
                                sideMeans.contains({index, second}));
                        mean = sideMeans.at({index, first}) - sideMeans.at({index, second});
                    } else {
                        require(domainIds.size() == 1);
                        mean = sideMeans.at({index, domainIds.front()});
                        panelMeans.emplace(index, mean);
                    }
                    if (!assigned.contains(triangle.faceId))
                        continue;
                    auto &sum = integrals[triangle.faceId];
                    const auto area = symmetryFactor * measures.areas[index];
                    sum.area += area;
                    sum.pressure += area * mean;
                    sum.absolutePressure += area * std::abs(mean);
                }
                if (request.coupling) {
                    const auto ports = lem::bemPorts(*request.coupling);
                    for (std::size_t i = 0; i < ports.size(); ++i) {
                        const auto &port = std::get<lem::BemPort>(ports[i]->component);
                        double area = 0, magnitude = 0;
                        std::complex<double> pressure{};
                        for (const auto &[index, mean] : panelMeans) {
                            const auto face = mesh.triangles[index].faceId;
                            if (std::find(port.faceIds.begin(), port.faceIds.end(), face) ==
                                port.faceIds.end())
                                continue;
                            const double weight = port.velocityWeights.contains(face)
                                                      ? port.velocityWeights.at(face)
                                                      : 1.;
                            const double weightedArea = measures.areas[index] * weight;
                            area += weightedArea;
                            pressure += weightedArea * mean;
                            magnitude += weightedArea * std::abs(mean);
                        }
                        require(std::isfinite(area) && area > 0);
                        pressure /= area;
                        magnitude /= area;
                        require(
                            std::abs(e.lemSolution->ports[i].pressure - pressure) <=
                            relativePortPressureClosureTolerance *
                                std::max(magnitude, std::abs(e.lemSolution->ports[i].pressure)));
                    }
                }
                std::set<std::uint32_t> seen;
                for (const auto &v : array(item["impedances"], assigned.size(), assigned.size())) {
                    const auto z = object(v);
                    keys(z,
                         {"faceId", "area", "meanPressure", "volumeVelocity", "volumeImpedance"});
                    FaceImpedance impedance{
                        integer(z["faceId"]), positive(z["area"]), complex(z["meanPressure"]),
                        complex(z["volumeVelocity"]), complex(z["volumeImpedance"])};
                    require(assigned.contains(impedance.faceId) &&
                            seen.insert(impedance.faceId).second &&
                            integrals.contains(impedance.faceId));
                    const auto &sum = integrals.at(impedance.faceId);
                    const double scale = std::max({std::abs(impedance.meanPressure),
                                                   std::abs(sum.pressure / sum.area),
                                                   sum.absolutePressure / sum.area});
                    require(
                        std::isfinite(scale) && close(impedance.area, sum.area) &&
                        std::abs(impedance.meanPressure - sum.pressure / sum.area) <=
                            1e-10 * scale &&
                        close(impedance.volumeVelocity, assigned.at(impedance.faceId) * sum.area) &&
                        close(impedance.volumeImpedance,
                              impedance.meanPressure / impedance.volumeVelocity));
                    e.impedances.push_back(impedance);
                }
            }
            if (!request.subdomains)
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
    (void)decodeAcousticModel(encodeAcousticModel(s.acousticModel));
    (void)decodeAcousticModel(encodeAcousticModel(a.resultModel));
    if (a.result && a.request)
        for (const auto &frequency : a.result->frequencies)
            for (const auto &excitation : frequency.excitations)
                require(excitation.observerValidity == a.request->observerValidity);
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
                                 {"singularOrder", int(s.singularOrder)},
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
    if (a.request && a.request->subdomains) {
        require(bool(a.sourceMesh) && !a.resultModel.domains.empty());
        result.insert("sourceMesh", encodeMesh(*a.sourceMesh));
        result.insert("resultModel", encodeAcousticModel(a.resultModel));
    } else
        require(!a.sourceMesh && a.resultModel.domains.empty());
    if (a.result && a.resultComplete)
        result.insert("resultComplete", *a.resultComplete);
    auto couplingSettings = result["settings"].toObject();
    require(s.ctaExcitationIndex >= -1 && s.ctaExcitationIndex <= 10000);
    couplingSettings.insert("acousticModel", encodeAcousticModel(s.acousticModel));
    require(std::isfinite(s.symmetryPlaneTolerance) && s.symmetryPlaneTolerance >= 1e-9 &&
            s.symmetryPlaneTolerance <= .01);
    couplingSettings.insert("symmetryXY",
                            QJsonObject{{"enabled", s.symmetryXYEnabled},
                                        {"planeDistanceTolerance", s.symmetryPlaneTolerance}});
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
SimulationArchive decodeSimulation(const QJsonValue &v, unsigned formatVersion) {
    const auto o = object(v);
    auto required = o;
    required.remove("measurements");
    required.remove("axisymmetricSampling");
    required.remove("resultComplete");
    required.remove("sourceMesh");
    required.remove("resultRegions");
    required.remove("resultModel");
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
    requiredSettings.remove("acousticRegions");
    requiredSettings.remove("acousticModel");
    requiredSettings.remove("symmetryXY");
    requiredSettings.remove("ctaExcitationIndex");
    requiredSettings.remove("beamwidthTarget");
    requiredSettings.remove("normAngle");
    requiredSettings.remove("DbScale");
    requiredSettings.remove("observationFrame");
    requiredSettings.remove("acceleration");
    requiredSettings.remove("fastMultipole");
    requiredSettings.remove("nearMaximumDepth");
    requiredSettings.remove("singularOrder");
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
    if (settings.contains("symmetryXY")) {
        const auto symmetry = object(settings["symmetryXY"]);
        keys(symmetry, {"enabled", "planeDistanceTolerance"});
        s.symmetryXYEnabled = boolean(symmetry["enabled"]);
        s.symmetryPlaneTolerance = range(symmetry["planeDistanceTolerance"], 1e-9, .01);
    }
    require(!(settings.contains("acousticRegions") && settings.contains("acousticModel")));
    require(formatVersion < 18 ? !settings.contains("acousticModel")
                               : settings.contains("acousticModel"));
    if (settings.contains("acousticRegions"))
        s.acousticModel =
            migrateAcousticRegions(decodeAcousticRegions(settings["acousticRegions"]));
    if (settings.contains("acousticModel"))
        s.acousticModel = decodeAcousticModel(settings["acousticModel"]);
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
    if (settings.contains("singularOrder")) {
        s.singularOrder = integer(settings["singularOrder"], 20);
        require(s.singularOrder >= 2);
    }
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
    const bool partitioned = !o["request"].isNull() && object(o["request"]).contains("subdomains");
    require(o.contains("sourceMesh") == partitioned);
    require(!(o.contains("resultRegions") && o.contains("resultModel")));
    require((formatVersion < 18 ? o.contains("resultRegions") : o.contains("resultModel")) ==
            partitioned);
    if (partitioned) {
        a.sourceMesh = std::make_shared<SurfaceMesh>(decodeMesh(o["sourceMesh"], false));
        a.resultModel = formatVersion < 18
                            ? migrateAcousticRegions(decodeAcousticRegions(o["resultRegions"]))
                            : decodeAcousticModel(o["resultModel"]);
        require(!a.resultModel.domains.empty());
    }
    if (!o["solvedMesh"].isNull())
        a.mesh = std::make_shared<SurfaceMesh>(decodeMesh(o["solvedMesh"], !partitioned));
    if (!o["request"].isNull()) {
        require(bool(a.mesh));
        a.request = decodeRequest(o["request"], *a.mesh);
        if (partitioned) {
            require(a.sourceMesh->vertices.size() <= a.mesh->vertices.size() &&
                    a.sourceMesh->triangles.size() <= a.mesh->triangles.size());
            require(std::equal(a.sourceMesh->vertices.begin(), a.sourceMesh->vertices.end(),
                               a.mesh->vertices.begin()) &&
                    std::equal(a.sourceMesh->triangles.begin(), a.sourceMesh->triangles.end(),
                               a.mesh->triangles.begin()));
        }
    }
    if (!o["result"].isNull()) {
        require(a.mesh && a.request.has_value());
        a.result = std::make_shared<BemResult>(
            decodeResult(o["result"], *a.request, *a.mesh, formatVersion));
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
    require(!a.current || a.request->exteriorBaffle == s.acousticModel.exteriorBaffle);
    require(!a.current || !partitioned || a.resultModel == s.acousticModel);
    if (a.current) {
        require(a.request->symmetry.has_value() == s.symmetryXYEnabled);
        if (a.request->symmetry)
            require(a.request->symmetry->planeDistanceTolerance == s.symmetryPlaneTolerance);
    }
    require(!a.current ||
            (std::none_of(s.acousticModel.domains.begin(), s.acousticModel.domains.end(),
                          [](const auto &r) { return r.needsAssignmentReview; }) &&
             std::none_of(s.acousticModel.interfaces.begin(), s.acousticModel.interfaces.end(),
                          [](const auto &i) { return i.needsAssignmentReview; })));
    if (a.current) {
        require(a.request->coupling.has_value() == s.couplingEnabled);
        if (a.request->coupling)
            require(lem::samePhysics(*a.request->coupling, s.couplingNetwork));
    }
    require(a.result ? a.revision > 0 : a.revision == 0 && !a.current);
    return a;
}
} // namespace hornsculpt::project_json
