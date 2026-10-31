#include "project_json.hpp"
#include "galerkin/exterior_domain.hpp"
#include <QSet>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hornsculpt::project_json {
[[noreturn]] void invalid(const QString &message) {
    throw std::invalid_argument(
        (QStringLiteral("Ungültige Projektdatei: ") + message).toStdString());
}
void keys(const QJsonObject &value, std::initializer_list<const char *> expected) {
    if (value.size() != static_cast<qsizetype>(expected.size()))
        invalid(QStringLiteral("Unbekannte oder fehlende Felder."));
    for (const auto key : expected)
        if (!value.contains(QLatin1String(key)))
            invalid(QStringLiteral("Feld fehlt: ") + QLatin1String(key));
}
QJsonObject object(const QJsonValue &value) {
    if (!value.isObject())
        invalid(QStringLiteral("Objekt erwartet."));
    return value.toObject();
}
QJsonArray array(const QJsonValue &value, qsizetype maximum, qsizetype minimum) {
    if (!value.isArray())
        invalid(QStringLiteral("Liste erwartet."));
    const auto result = value.toArray();
    if (result.size() < minimum || result.size() > maximum)
        invalid(QStringLiteral("Listenlänge außerhalb des erlaubten Bereichs."));
    return result;
}
double number(const QJsonValue &value) {
    if (!value.isDouble() || !std::isfinite(value.toDouble()))
        invalid(QStringLiteral("Endliche Zahl erwartet."));
    return value.toDouble();
}
double positive(const QJsonValue &value) {
    const double result = number(value);
    if (result <= 0)
        invalid(QStringLiteral("Positive Zahl erwartet."));
    return result;
}
bool boolean(const QJsonValue &value) {
    if (!value.isBool())
        invalid(QStringLiteral("Boolescher Wert erwartet."));
    return value.toBool();
}
QString text(const QJsonValue &value, bool allowEmpty) {
    if (!value.isString() || value.toString().size() > 4096 ||
        (!allowEmpty && value.toString().isEmpty()))
        invalid(QStringLiteral("Ungültiger Text."));
    return value.toString();
}
std::uint32_t integer(const QJsonValue &value, std::uint32_t maximum) {
    const double result = number(value);
    if (result < 0 || result > maximum || result != std::floor(result))
        invalid(QStringLiteral("Ungültige Ganzzahl oder ID."));
    return static_cast<std::uint32_t>(result);
}
std::uint64_t wideInteger(const QJsonValue &value) {
    const auto string = text(value);
    if (string.size() > 20 || (string.size() > 1 && string[0] == QChar('0')))
        invalid(QStringLiteral("Ungültige 64-Bit-Dezimalzahl."));
    for (const auto character : string)
        if (character < QChar('0') || character > QChar('9'))
            invalid(QStringLiteral("Ungültige 64-Bit-Dezimalzahl."));
    bool ok = false;
    const auto result = string.toULongLong(&ok);
    if (!ok)
        invalid(QStringLiteral("64-Bit-Zahl außerhalb des Bereichs."));
    return result;
}
QString encodeWide(std::uint64_t value) {
    return QString::number(static_cast<qulonglong>(value));
}
Vec3 vector(const QJsonValue &value) {
    const auto values = array(value, 3, 3);
    return {number(values[0]), number(values[1]), number(values[2])};
}
QJsonArray encodeVector(Vec3 value) {
    return {value.x, value.y, value.z};
}
std::complex<double> complex(const QJsonValue &value) {
    const auto values = array(value, 2, 2);
    return {number(values[0]), number(values[1])};
}
QJsonArray encodeComplex(std::complex<double> value) {
    return {value.real(), value.imag()};
}
bool sameMesh(const SurfaceMesh &a, const SurfaceMesh &b) {
    return a == b;
}
std::map<std::uint32_t, BoundaryCondition> meshFaces(const SurfaceMesh &mesh) {
    std::map<std::uint32_t, BoundaryCondition> result;
    for (const auto &triangle : mesh.triangles) {
        if (!triangle.faceId || (triangle.boundary != BoundaryCondition::wall &&
                                 triangle.boundary != BoundaryCondition::driven))
            invalid(QStringLiteral("Ungültige Netzfläche oder Randbedingung."));
        const auto [it, inserted] = result.emplace(triangle.faceId, triangle.boundary);
        if (!inserted && it->second != triangle.boundary)
            invalid(QStringLiteral("Widersprüchliche Randbedingungen einer Fläche."));
    }
    return result;
}
QJsonObject encodeMesh(const SurfaceMesh &mesh) {
    QJsonArray vertices, triangles;
    for (const auto vertex : mesh.vertices)
        vertices.append(encodeVector(vertex));
    for (const auto &triangle : mesh.triangles)
        triangles.append(QJsonArray{
            static_cast<double>(triangle.vertices[0]), static_cast<double>(triangle.vertices[1]),
            static_cast<double>(triangle.vertices[2]), static_cast<double>(triangle.faceId),
            static_cast<int>(triangle.boundary)});
    QJsonObject result{{"vertices", vertices}, {"triangles", triangles}};
    if (!mesh.fixedVertices.empty()) {
        QJsonArray fixed;
        for (const auto vertex : mesh.fixedVertices)
            fixed.append(static_cast<double>(vertex));
        result.insert("fixedVertices", fixed);
    }
    return result;
}
SurfaceMesh decodeMesh(const QJsonValue &value, bool closed) {
    const auto mesh = object(value);
    auto required = mesh;
    required.remove("fixedVertices");
    keys(required, {"vertices", "triangles"});
    const auto vertices = array(mesh["vertices"], 1000000, 3);
    const auto triangles = array(mesh["triangles"], 2000000, 1);
    SurfaceMesh result;
    result.vertices.reserve(static_cast<std::size_t>(vertices.size()));
    result.triangles.reserve(static_cast<std::size_t>(triangles.size()));
    for (const auto &vertex : vertices)
        result.vertices.push_back(vector(vertex));
    for (const auto &triangle : triangles) {
        const auto data = array(triangle, 5, 5);
        result.triangles.push_back(
            {{integer(data[0], static_cast<std::uint32_t>(vertices.size() - 1)),
              integer(data[1], static_cast<std::uint32_t>(vertices.size() - 1)),
              integer(data[2], static_cast<std::uint32_t>(vertices.size() - 1))},
             integer(data[3]),
             static_cast<BoundaryCondition>(integer(data[4], 1))});
    }
    if (mesh.contains("fixedVertices")) {
        QSet<std::uint32_t> seen;
        for (const auto &indexValue : array(mesh["fixedVertices"], vertices.size())) {
            const auto vertex =
                integer(indexValue, static_cast<std::uint32_t>(vertices.size() - 1));
            if (seen.contains(vertex))
                invalid(QStringLiteral("Doppelt definierter fester Netzknoten."));
            seen.insert(vertex);
            result.fixedVertices.push_back(vertex);
        }
    }
    (void)meshFaces(result);
    const auto quality = inspectMesh(result);
    if (quality.degenerateTriangles || (closed && !quality.watertight()))
        invalid(QStringLiteral("Ungültige oder offene BEM-Randfläche."));
    if (closed) {
        galerkin::validateExteriorComponents(result);
    }
    return result;
}
namespace {
QJsonObject encodeDirection(const OsseDirectionalParameters &p) {
    return {{"throatAngle", p.throatAngle},
            {"k", p.k},
            {"terminationScale", p.terminationScale},
            {"terminationQ", p.terminationQ},
            {"terminationExponent", p.terminationExponent}};
}
OsseDirectionalParameters decodeDirection(const QJsonValue &value) {
    const auto p = object(value);
    keys(p, {"throatAngle", "k", "terminationScale", "terminationQ", "terminationExponent"});
    return {number(p["throatAngle"]), number(p["k"]), number(p["terminationScale"]),
            number(p["terminationQ"]), number(p["terminationExponent"])};
}
QJsonObject encodeContour(const ContourParameters &c) {
    const auto &sf = c.superformula;
    return {{"type", static_cast<int>(c.type)},
            {"width", c.width},
            {"aspectRatio", c.aspectRatio},
            {"rotation", c.rotation},
            {"superellipseExponent", c.superellipseExponent},
            {"superformula", QJsonObject{{"a", sf.a},
                                         {"b", sf.b},
                                         {"m1", sf.m1},
                                         {"m2", sf.m2},
                                         {"n1", sf.n1},
                                         {"n2", sf.n2},
                                         {"n3", sf.n3}}}};
}
ContourParameters decodeContour(const QJsonObject &c) {
    keys(c, {"type", "width", "aspectRatio", "rotation", "superellipseExponent", "superformula"});
    ContourParameters result;
    result.type = static_cast<CurveType>(integer(c["type"], 1));
    result.width = positive(c["width"]);
    result.aspectRatio = positive(c["aspectRatio"]);
    result.rotation = number(c["rotation"]);
    result.superellipseExponent = number(c["superellipseExponent"]);
    const auto sf = object(c["superformula"]);
    keys(sf, {"a", "b", "m1", "m2", "n1", "n2", "n3"});
    result.superformula = {number(sf["a"]),  number(sf["b"]),  number(sf["m1"]), number(sf["m2"]),
                           number(sf["n1"]), number(sf["n2"]), number(sf["n3"])};
    (void)Superellipse(1, 1, result.superellipseExponent);
    (void)Superformula(result.superformula);
    (void)PhysicalContour(result);
    return result;
}
} // namespace
QJsonObject encodeDraft(const Project &project) {
    const auto &p = project.parameters;
    const auto &g = p.guide;
    auto guide = encodeContour(g);
    guide.insert("positionMode", static_cast<int>(g.positionMode));
    guide.insert("position", g.position);
    QJsonArray knots;
    for (const auto &knot : p.transition.knots)
        knots.append(QJsonArray{knot.fraction, knot.blend});
    auto profile = encodeDirection(p.profile);
    profile.insert("vertical", encodeDirection(p.verticalProfile));
    return {{"profileType", "OS-SE"},
            {"profile", profile},
            {"horizontalAngle", p.horizontalAngle},
            {"length", p.length},
            {"entry", QJsonObject{{"width", p.entry.width},
                                  {"height", p.entry.height},
                                  {"exponent", p.entry.exponent}}},
            {"slot", QJsonObject{{"enabled", p.slot.enabled},
                                 {"aperture", QJsonObject{{"width", p.slot.aperture.width},
                                                          {"height", p.slot.aperture.height},
                                                          {"exponent", p.slot.aperture.exponent}}},
                                 {"adapterLength", p.slot.adapterLength},
                                 {"straightLength", p.slot.straightLength},
                                 {"transitionLength", p.slot.transitionLength},
                                 {"contour", static_cast<int>(p.slot.contour)},
                                 {"cornerRadius", p.slot.cornerRadius},
                                 {"exit", static_cast<int>(p.slot.exit)},
                                 {"exitRadius", p.slot.exitRadius}}},
            {"verticalAngle", p.verticalAngle},
            {"asymmetry", QJsonObject{{"enabled", p.asymmetricEnabled},
                                      {"leftAngle", p.leftAngle},
                                      {"bottomAngle", p.bottomAngle},
                                      {"leftProfile", encodeDirection(p.leftProfile)},
                                      {"bottomProfile", encodeDirection(p.bottomProfile)}}},
            {"crossSectionExponent", p.crossSectionExponent},
            {"wallThickness", p.wallThickness},
            {"guideEnabled", p.guideEnabled},
            {"transition", QJsonObject{{"method", static_cast<int>(p.transition.method)},
                                       {"startFraction", p.transition.startFraction},
                                       {"endFraction", p.transition.endFraction},
                                       {"knots", knots}}},
            {"guide", guide},
            {"mouthMorph", QJsonObject{{"enabled", p.mouthMorph.enabled},
                                       {"startFraction", p.mouthMorph.startFraction},
                                       {"endFraction", p.mouthMorph.endFraction},
                                       {"target", encodeContour(p.mouthMorph.target)}}}};
}
void decodeDraft(const QJsonObject &root, Project &project) {
    if (text(root["profileType"]) != "OS-SE")
        invalid(QStringLiteral("Unbekanntes Profil."));
    auto &p = project.parameters;
    const auto os = object(root["profile"]);
    auto horizontal = os;
    horizontal.remove("vertical");
    const auto version = integer(root["version"], 15);
    const bool oldFormat = version < 4;
    if (oldFormat) {
        keys(horizontal, {"throatRadius", "length", "throatAngle", "k", "terminationScale",
                          "terminationQ", "terminationExponent"});
        p.entry = {2 * positive(os["throatRadius"]), 2 * positive(os["throatRadius"]), 2};
        p.length = positive(os["length"]);
        project.meshing.targetEdgeLength = positive(root["targetEdgeLength"]);
    } else {
        keys(horizontal,
             {"throatAngle", "k", "terminationScale", "terminationQ", "terminationExponent"});
        if (!os.contains("vertical"))
            invalid(QStringLiteral("Vertikale Profilparameter fehlen."));
        const auto entry = object(root["entry"]);
        keys(entry, {"width", "height", "exponent"});
        p.entry = {positive(entry["width"]), positive(entry["height"]),
                   positive(entry["exponent"])};
        p.length = positive(root["length"]);
    }
    p.slot = {};
    if (version >= 13) {
        const auto slot = object(root["slot"]);
        if (version >= 14)
            keys(slot, {"enabled", "aperture", "adapterLength", "straightLength",
                        "transitionLength", "contour", "cornerRadius", "exit", "exitRadius"});
        else
            keys(slot, {"enabled", "aperture", "adapterLength", "straightLength",
                        "transitionLength"});
        const auto aperture = object(slot["aperture"]);
        keys(aperture, {"width", "height", "exponent"});
        p.slot.enabled = boolean(slot["enabled"]);
        p.slot.aperture = {positive(aperture["width"]), positive(aperture["height"]),
                           positive(aperture["exponent"])};
        (void)Superellipse(p.slot.aperture.width / 2, p.slot.aperture.height / 2,
                           p.slot.aperture.exponent);
        p.slot.adapterLength = positive(slot["adapterLength"]);
        p.slot.straightLength = number(slot["straightLength"]);
        p.slot.transitionLength = number(slot["transitionLength"]);
        if (p.slot.straightLength < 0 || p.slot.transitionLength < 0)
            invalid(QStringLiteral("Spalt- und Anschlusslänge dürfen nicht negativ sein."));
        if (version >= 14) {
            p.slot.contour = static_cast<SlotContour>(integer(slot["contour"], 1));
            p.slot.cornerRadius = number(slot["cornerRadius"]);
            p.slot.exit = static_cast<SlotExit>(integer(slot["exit"], 1));
            p.slot.exitRadius = number(slot["exitRadius"]);
            if (p.slot.cornerRadius < 0 || p.slot.exitRadius < 0)
                invalid(QStringLiteral("Spaltradien dürfen nicht negativ sein."));
        }
    }
    p.profile.throatAngle = number(os["throatAngle"]);
    p.profile.k = number(os["k"]);
    p.profile.terminationScale = number(os["terminationScale"]);
    p.profile.terminationQ = number(os["terminationQ"]);
    p.profile.terminationExponent = number(os["terminationExponent"]);
    p.verticalProfile = static_cast<const OsseDirectionalParameters &>(p.profile);
    if (os.contains("vertical")) {
        p.verticalProfile = decodeDirection(os["vertical"]);
    }
    p.horizontalAngle = number(root["horizontalAngle"]);
    p.verticalAngle = number(root["verticalAngle"]);
    p.leftProfile = p.profile;
    p.bottomProfile = p.verticalProfile;
    p.leftAngle = p.horizontalAngle;
    p.bottomAngle = p.verticalAngle;
    p.asymmetricEnabled = false;
    if (root.contains("asymmetry")) {
        const auto asymmetry = object(root["asymmetry"]);
        keys(asymmetry, {"enabled", "leftAngle", "bottomAngle", "leftProfile", "bottomProfile"});
        p.asymmetricEnabled = boolean(asymmetry["enabled"]);
        p.leftAngle = number(asymmetry["leftAngle"]);
        p.bottomAngle = number(asymmetry["bottomAngle"]);
        p.leftProfile = decodeDirection(asymmetry["leftProfile"]);
        p.bottomProfile = decodeDirection(asymmetry["bottomProfile"]);
        // Dormant controls must remain valid when the user re-enables the mode.
        for (const auto vertical : {false, true}) {
            OsseParameters direction;
            static_cast<OsseDirectionalParameters &>(direction) =
                vertical ? p.bottomProfile : p.leftProfile;
            direction.throatRadius = .5 * (vertical ? p.entry.height : p.entry.width);
            direction.length = p.length;
            const OsseProfile profile(direction);
            const double radius = profile.radius(p.length, vertical ? p.bottomAngle : p.leftAngle);
            if (!std::isfinite(radius) || radius <= 0)
                invalid(QStringLiteral("Ungültiges gespeichertes Seitenprofil."));
        }
    }
    p.crossSectionExponent = number(root["crossSectionExponent"]);
    p.wallThickness = number(root["wallThickness"]);
    p.guideEnabled = boolean(root["guideEnabled"]);
    const auto transition = object(root["transition"]);
    keys(transition, {"method", "startFraction", "endFraction", "knots"});
    p.transition.method = static_cast<TransitionMethod>(integer(transition["method"], 1));
    p.transition.startFraction = number(transition["startFraction"]);
    p.transition.endFraction = number(transition["endFraction"]);
    p.transition.knots.clear();
    for (const auto &knot : array(transition["knots"], 1000, 2)) {
        const auto values = array(knot, 2, 2);
        p.transition.knots.push_back({number(values[0]), number(values[1])});
    }
    const auto g = object(root["guide"]);
    keys(g, {"type", "width", "aspectRatio", "rotation", "positionMode", "position",
             "superellipseExponent", "superformula"});
    auto contour = g;
    contour.remove("positionMode");
    contour.remove("position");
    static_cast<ContourParameters &>(p.guide) = decodeContour(contour);
    p.guide.positionMode = static_cast<GuidePositionMode>(integer(g["positionMode"], 1));
    p.guide.position = number(g["position"]);
    if (!oldFormat) {
        const auto morph = object(root["mouthMorph"]);
        keys(morph, {"enabled", "startFraction", "endFraction", "target"});
        p.mouthMorph.enabled = boolean(morph["enabled"]);
        p.mouthMorph.startFraction = number(morph["startFraction"]);
        p.mouthMorph.endFraction = number(morph["endFraction"]);
        p.mouthMorph.target = decodeContour(object(morph["target"]));
        (void)MouthMorph(p.mouthMorph);
    }
    auto transitionParameters = p.transition;
    transitionParameters.method = TransitionMethod::smoothRange;
    (void)ShapeTransition(transitionParameters);
    transitionParameters.method = TransitionMethod::lut;
    (void)ShapeTransition(transitionParameters);
    if (p.guide.width <= 0 || p.guide.aspectRatio <= 0 || p.guide.position <= 0 ||
        (p.guide.positionMode == GuidePositionMode::relative && p.guide.position > 1))
        invalid(QStringLiteral("Ungültige Leitkurvenparameter."));
    (void)DraftGeometry(p);
}
} // namespace hornsculpt::project_json
