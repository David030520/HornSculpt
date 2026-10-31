#include "sweep_model.hpp"
#include "project_json.hpp"
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <QStringList>
#include <QUuid>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace hornsculpt {
namespace {
struct Binding {
    SweepParameter parameter;
    double (*get)(const DraftParameters &);
    void (*set)(DraftParameters &, double);
    bool asymmetricOnly{};
    QString symmetricLabel;
};
#define SWEEP_BIND(key, label, unit, scale, member)                                                \
    {{QStringLiteral(key), QStringLiteral(label), QStringLiteral(unit), scale},                    \
     [](const DraftParameters &p) { return p.member; },                                            \
     [](DraftParameters &p, double value) { p.member = value; }}
#define SWEEP_DIRECTION(key, label, symmetric, unit, scale, member, only)                           \
    {{QStringLiteral(key), QStringLiteral(label), QStringLiteral(unit), scale},                    \
     [](const DraftParameters &p) { return p.member; },                                            \
     [](DraftParameters &p, double value) { p.member = value; }, only, QStringLiteral(symmetric)}
const std::vector<Binding> &bindings() {
    static const std::vector<Binding> result{
        {{"meshTargetEdgeLength", "Zielkante", "mm", 1000.}, nullptr, nullptr},
        SWEEP_BIND("length", "Länge", "mm", 1000., length),
        SWEEP_BIND("entry.width", "Eintritt Breite", "mm", 1000., entry.width),
        SWEEP_BIND("entry.height", "Eintritt Höhe", "mm", 1000., entry.height),
        SWEEP_BIND("entry.exponent", "Eintritt Exponent", "", 1., entry.exponent),
        SWEEP_BIND("slot.aperture.width", "Spalt Breite X", "mm", 1000., slot.aperture.width),
        SWEEP_BIND("slot.aperture.height", "Spalt Höhe Y", "mm", 1000., slot.aperture.height),
        SWEEP_BIND("slot.aperture.exponent", "Spalt Exponent", "", 1., slot.aperture.exponent),
        SWEEP_BIND("slot.cornerRadius", "Spalt Eckradius", "mm", 1000., slot.cornerRadius),
        SWEEP_BIND("slot.exitRadius", "Spalt Ausgangsradius", "mm", 1000., slot.exitRadius),
        SWEEP_BIND("slot.adapterLength", "Spalt Adapterlänge", "mm", 1000., slot.adapterLength),
        SWEEP_BIND("slot.straightLength", "Spalt gerade Länge", "mm", 1000., slot.straightLength),
        SWEEP_BIND("slot.transitionLength", "Spalt Anschlusslänge", "mm", 1000.,
                   slot.transitionLength),
        SWEEP_DIRECTION("horizontalAngle", "Öffnungswinkel H / rechts", "Öffnungswinkel H", "°",
                        180. / 3.14159265358979323846, horizontalAngle, false),
        SWEEP_DIRECTION("verticalAngle", "Öffnungswinkel V / oben", "Öffnungswinkel V", "°",
                        180. / 3.14159265358979323846, verticalAngle, false),
        SWEEP_DIRECTION("leftAngle", "Öffnungswinkel links", "", "°", 180. / 3.14159265358979323846,
                        leftAngle, true),
        SWEEP_DIRECTION("bottomAngle", "Öffnungswinkel unten", "", "°", 180. / 3.14159265358979323846,
                        bottomAngle, true),
        SWEEP_DIRECTION("profile.throatAngle", "Halswinkel H / rechts", "Halswinkel H", "°",
                        180. / 3.14159265358979323846, profile.throatAngle, false),
        SWEEP_DIRECTION("verticalProfile.throatAngle", "Halswinkel V / oben", "Halswinkel V", "°",
                        180. / 3.14159265358979323846, verticalProfile.throatAngle, false),
        SWEEP_DIRECTION("leftProfile.throatAngle", "Halswinkel links", "", "°",
                        180. / 3.14159265358979323846, leftProfile.throatAngle, true),
        SWEEP_DIRECTION("bottomProfile.throatAngle", "Halswinkel unten", "", "°",
                        180. / 3.14159265358979323846, bottomProfile.throatAngle, true),
        SWEEP_DIRECTION("profile.k", "k H / rechts", "k H", "", 1., profile.k, false),
        SWEEP_DIRECTION("verticalProfile.k", "k V / oben", "k V", "", 1., verticalProfile.k, false),
        SWEEP_DIRECTION("leftProfile.k", "k links", "", "", 1., leftProfile.k, true),
        SWEEP_DIRECTION("bottomProfile.k", "k unten", "", "", 1., bottomProfile.k, true),
        SWEEP_DIRECTION("profile.terminationScale", "s H / rechts", "s H", "", 1.,
                        profile.terminationScale, false),
        SWEEP_DIRECTION("verticalProfile.terminationScale", "s V / oben", "s V", "", 1.,
                        verticalProfile.terminationScale, false),
        SWEEP_DIRECTION("leftProfile.terminationScale", "s links", "", "", 1.,
                        leftProfile.terminationScale, true),
        SWEEP_DIRECTION("bottomProfile.terminationScale", "s unten", "", "", 1.,
                        bottomProfile.terminationScale, true),
        SWEEP_DIRECTION("profile.terminationQ", "q H / rechts", "q H", "", 1.,
                        profile.terminationQ, false),
        SWEEP_DIRECTION("verticalProfile.terminationQ", "q V / oben", "q V", "", 1.,
                        verticalProfile.terminationQ, false),
        SWEEP_DIRECTION("leftProfile.terminationQ", "q links", "", "", 1.,
                        leftProfile.terminationQ, true),
        SWEEP_DIRECTION("bottomProfile.terminationQ", "q unten", "", "", 1.,
                        bottomProfile.terminationQ, true),
        SWEEP_DIRECTION("profile.terminationExponent", "n H / rechts", "n H", "", 1.,
                        profile.terminationExponent, false),
        SWEEP_DIRECTION("verticalProfile.terminationExponent", "n V / oben", "n V", "", 1.,
                        verticalProfile.terminationExponent, false),
        SWEEP_DIRECTION("leftProfile.terminationExponent", "n links", "", "", 1.,
                        leftProfile.terminationExponent, true),
        SWEEP_DIRECTION("bottomProfile.terminationExponent", "n unten", "", "", 1.,
                        bottomProfile.terminationExponent, true),
        SWEEP_BIND("crossSectionExponent", "Querschnitt Exponent", "", 1., crossSectionExponent),
        SWEEP_BIND("wallThickness", "Wandstärke radial", "mm", 1000., wallThickness),
        SWEEP_BIND("guide.width", "Leitkurve Breite", "mm", 1000., guide.width),
        SWEEP_BIND("guide.aspectRatio", "Leitkurve Seitenverhältnis", "", 1., guide.aspectRatio),
        SWEEP_BIND("guide.rotation", "Leitkurve Rotation", "°", 180. / 3.14159265358979323846,
                   guide.rotation),
        SWEEP_BIND("guide.superellipseExponent", "Leitkurve Superellipse n", "", 1.,
                   guide.superellipseExponent),
        SWEEP_BIND("guide.superformula.a", "Superformula a", "", 1., guide.superformula.a),
        SWEEP_BIND("guide.superformula.b", "Superformula b", "", 1., guide.superformula.b),
        SWEEP_BIND("guide.superformula.m1", "Superformula m1", "", 1., guide.superformula.m1),
        SWEEP_BIND("guide.superformula.m2", "Superformula m2", "", 1., guide.superformula.m2),
        SWEEP_BIND("guide.superformula.n1", "Superformula n1", "", 1., guide.superformula.n1),
        SWEEP_BIND("guide.superformula.n2", "Superformula n2", "", 1., guide.superformula.n2),
        SWEEP_BIND("guide.superformula.n3", "Superformula n3", "", 1., guide.superformula.n3),
        SWEEP_BIND("transition.startFraction", "Übergang Beginn", "%", 100.,
                   transition.startFraction),
        SWEEP_BIND("transition.endFraction", "Übergang Ende", "%", 100., transition.endFraction)};
    return result;
}
#undef SWEEP_BIND
#undef SWEEP_DIRECTION
const Binding &binding(const QString &key) {
    const auto &all = bindings();
    const auto found =
        std::find_if(all.begin(), all.end(), [&](const auto &b) { return b.parameter.key == key; });
    if (found == all.end())
        throw std::invalid_argument("Unbekannter Studienparameter.");
    return *found;
}
SweepParameter presentation(const Binding &entry, const DraftParameters &baseline) {
    auto parameter = entry.parameter;
    if (!baseline.asymmetricEnabled && !entry.symmetricLabel.isEmpty())
        parameter.label = entry.symmetricLabel;
    return parameter;
}
[[noreturn]] void invalid() {
    throw std::invalid_argument("Ungültige HornSculpt-Studie.");
}
void keys(const QJsonObject &object, std::initializer_list<const char *> expected) {
    if (object.size() != static_cast<qsizetype>(expected.size()))
        invalid();
    for (const auto key : expected)
        if (!object.contains(QLatin1String(key)))
            invalid();
}
QJsonObject object(const QJsonValue &value) {
    if (!value.isObject())
        invalid();
    return value.toObject();
}
double number(const QJsonValue &value) {
    if (!value.isDouble() || !std::isfinite(value.toDouble()))
        invalid();
    return value.toDouble();
}
QString text(const QJsonValue &value) {
    if (!value.isString() || value.toString().size() > 4096)
        invalid();
    return value.toString();
}
QString stateName(SweepState state) {
    switch (state) {
    case SweepState::pending:
        return "pending";
    case SweepState::running:
        return "running";
    case SweepState::completed:
        return "completed";
    case SweepState::failed:
        return "failed";
    case SweepState::cancelled:
        return "cancelled";
    }
    invalid();
}
SweepState stateValue(const QString &text) {
    for (const auto state : {SweepState::pending, SweepState::running, SweepState::completed,
                             SweepState::failed, SweepState::cancelled})
        if (stateName(state) == text)
            return state;
    invalid();
}
QString fileName(std::size_t index) {
    return QString("variant-%1.hsproj").arg(static_cast<qulonglong>(index + 1), 6, 10, QChar('0'));
}
Project inputsOnly(Project project) {
    project.comparisonSnapshots.clear();
    project.revision = 1;
    project.geometryRevision = project.meshRevision = 0;
    project.preview.reset();
    project.mesh.reset();
    project.simulation.mesh.reset();
    project.simulation.request.reset();
    project.simulation.result.reset();
    project.simulation.resultComplete.reset();
    project.simulation.axisymmetricSampling.reset();
    project.simulation.revision = 0;
    project.simulation.current = false;
    return project;
}
void validateAxesForBaseline(const Project &baseline, std::span<const SweepAxis> axes) {
    validateSweepMeshSettings(baseline);
    for (const auto &axis : axes) {
        if (binding(axis.parameter).asymmetricOnly && !baseline.parameters.asymmetricEnabled)
            throw std::invalid_argument(
                "Linke und untere Profilparameter erfordern den asymmetrischen Draft-Modus.");
        if (axis.parameter.startsWith("guide.") && !baseline.parameters.guideEnabled)
            throw std::invalid_argument("Die Leitkurve ist im Ausgangsprojekt deaktiviert.");
        if (axis.parameter.startsWith("slot.") && !baseline.parameters.slot.enabled)
            throw std::invalid_argument("Der Diffraktionsspalt ist im Ausgangsprojekt deaktiviert.");
        if (axis.parameter == "slot.cornerRadius" &&
            baseline.parameters.slot.contour != SlotContour::roundedRectangle)
            throw std::invalid_argument(
                "Der Spalt-Eckradius erfordert die Kontur Rechteck mit Eckradius.");
        if (axis.parameter == "slot.aperture.exponent" &&
            baseline.parameters.slot.contour != SlotContour::superellipse)
            throw std::invalid_argument("Der Spaltexponent erfordert die Kontur Superellipse.");
        if (axis.parameter == "slot.exitRadius" &&
            baseline.parameters.slot.exit != SlotExit::circularFillet)
            throw std::invalid_argument("Der Spalt-Ausgangsradius erfordert den Kreisbogen-Anschluss.");
        if (axis.parameter == "slot.transitionLength" &&
            baseline.parameters.slot.exit != SlotExit::axialBlend)
            throw std::invalid_argument("Die Anschlusslänge erfordert den axialen Glättungsmodus.");
        if (axis.parameter.startsWith("transition.") &&
            baseline.parameters.transition.method != TransitionMethod::smoothRange)
            throw std::invalid_argument("Übergangsbereiche erfordern den Prozentmodus statt LUT.");
        if (axis.parameter != "meshTargetEdgeLength")
            continue;
        auto variant = baseline;
        variant.meshing.targetEdgeLength = std::min(axis.start, axis.end);
        validateSweepMeshSettings(variant);
        variant.meshing.targetEdgeLength = std::max(axis.start, axis.end);
        validateSweepMeshSettings(variant);
    }
}
void validateConvergence(const SweepStudy &study) {
    if (study.frequencyOverride.empty())
        return;
    if (study.axes.size() != 1 || study.axes.front().parameter != "meshTargetEdgeLength" ||
        study.axes.front().count < 2 || !(study.axes.front().end > 0) ||
        !(study.axes.front().end < study.axes.front().start))
        throw std::invalid_argument(
            "Netzkonvergenz benötigt mindestens zwei Zielkanten von grob nach fein.");
}
void checkGridMemory(std::size_t count) {
    const long double bytes =
        static_cast<long double>(count) * (sizeof(SweepVariant) + 2 * sizeof(double) + 1024);
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        count > std::vector<SweepVariant>{}.max_size() ||
        bytes > std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("Das Studienraster ist zu groß.");
#ifdef _WIN32
    MEMORYSTATUSEX memory{};
    memory.dwLength = sizeof(memory);
    if (GlobalMemoryStatusEx(&memory) && bytes > memory.ullAvailPhys)
        throw std::invalid_argument("Das Studienraster benötigt mehr verfügbaren Arbeitsspeicher.");
#endif
}
} // namespace

std::span<const SweepParameter> sweepParameters() {
    static const auto parameters = [] {
        std::vector<SweepParameter> result;
        for (const auto &b : bindings())
            result.push_back(b.parameter);
        return result;
    }();
    return parameters;
}
std::vector<SweepParameter> sweepParameters(const DraftParameters &baseline) {
    std::vector<SweepParameter> result;
    for (const auto &entry : bindings()) {
        if (baseline.asymmetricEnabled || !entry.asymmetricOnly)
            result.push_back(presentation(entry, baseline));
    }
    return result;
}
double sweepParameterValue(const DraftParameters &parameters, const QString &key) {
    const auto &b = binding(key);
    if (!b.get)
        throw std::invalid_argument("Die Zielkante gehört zu den Netzeinstellungen.");
    return b.get(parameters);
}
void setSweepParameterValue(DraftParameters &parameters, const QString &key, double value) {
    if (!std::isfinite(value))
        invalid();
    const auto &b = binding(key);
    if (!b.set)
        throw std::invalid_argument("Die Zielkante gehört zu den Netzeinstellungen.");
    b.set(parameters, value);
}
void validateSweepMeshSettings(const Project &project) {
    validateMeshSettings(project.meshing);
    QStringList conflicts;
    const auto check = [&](const auto &targets, const QString &mode) {
        for (const auto &[face, size] : targets)
            if (!face || !std::isfinite(size) ||
                size < project.meshing.adaptive.minimumEdgeLength ||
                size > project.meshing.targetEdgeLength)
                conflicts << mode + " Fläche " + QString::number(face);
    };
    check(project.draftFaceTargetEdgeLengths, "Draft");
    check(project.cadFaceTargetEdgeLengths, "FreeCAD");
    if (!conflicts.empty())
        throw std::invalid_argument(
            ("Lokale Zielkante außerhalb Mindest-/Globalkante: " + conflicts.join(", ") + ".")
                .toStdString());
}
std::vector<SweepVariant> makeSweepGrid(std::span<const SweepAxis> axes) {
    if (axes.empty() || axes.size() > 2 ||
        (axes.size() == 2 && axes[0].parameter == axes[1].parameter))
        throw std::invalid_argument("Eine oder zwei verschiedene Parameterachsen erforderlich.");
    std::size_t total = 1;
    for (const auto &axis : axes) {
        (void)binding(axis.parameter);
        if (!axis.count || !std::isfinite(axis.start) || !std::isfinite(axis.end) ||
            (axis.count == 1 && axis.start != axis.end) ||
            axis.count > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
            total > std::numeric_limits<std::size_t>::max() / axis.count)
            throw std::invalid_argument("Ungültiger Parameterbereich oder Punktzahl.");
        total *= axis.count;
    }
    checkGridMemory(total);
    std::vector<SweepVariant> variants;
    variants.reserve(total);
    for (std::size_t i = 0; i < total; ++i) {
        SweepVariant variant;
        auto remaining = i;
        for (const auto &axis : axes) {
            const auto coordinate = remaining % axis.count;
            remaining /= axis.count;
            const double value = axis.count == 1
                                     ? axis.start
                                     : std::lerp(axis.start, axis.end,
                                                 static_cast<double>(coordinate) /
                                                     static_cast<double>(axis.count - 1));
            variant.values.push_back(value);
        }
        variant.projectFile = fileName(i);
        variants.push_back(std::move(variant));
    }
    return variants;
}
Project sweepVariantInputs(const SweepStudy &study, std::size_t index) {
    if (index >= study.variants.size() || study.variants[index].values.size() != study.axes.size())
        invalid();
    validateAxesForBaseline(study.baseline, study.axes);
    validateConvergence(study);
    auto result = inputsOnly(study.baseline);
    for (std::size_t a = 0; a < study.axes.size(); ++a) {
        if (study.axes[a].parameter == "meshTargetEdgeLength")
            result.meshing.targetEdgeLength = study.variants[index].values[a];
        else {
            if (result.source != ProjectSource::draft)
                invalid();
            setSweepParameterValue(result.parameters, study.axes[a].parameter,
                                   study.variants[index].values[a]);
        }
    }
    validateSweepMeshSettings(result);
    return result;
}
bool sweepVariantMatches(const SweepStudy &study, std::size_t index, const Project &project) {
    auto expected = sweepVariantInputs(study, index);
    if (!study.frequencyOverride.empty()) {
        expected.simulation.settings.minimumFrequency = study.frequencyOverride.front();
        expected.simulation.settings.maximumFrequency = study.frequencyOverride.back();
        expected.simulation.settings.frequencyCount = study.frequencyOverride.size();
    }
    const auto actual = inputsOnly(project);
    const auto &a = actual.meshing;
    const auto &b = expected.meshing;
    return actual.source == expected.source &&
           project_json::encodeDraft(actual) == project_json::encodeDraft(expected) &&
           project_json::encodeAssembly(actual.assembly) ==
               project_json::encodeAssembly(expected.assembly) &&
           bool(actual.cad) == bool(expected.cad) &&
           (!actual.cad ||
            project_json::encodeCad(*actual.cad) == project_json::encodeCad(*expected.cad)) &&
           a.targetEdgeLength == b.targetEdgeLength && a.adaptive.enabled == b.adaptive.enabled &&
           a.adaptive.minimumEdgeLength == b.adaptive.minimumEdgeLength &&
           a.adaptive.geometricTolerance == b.adaptive.geometricTolerance &&
           a.adaptive.growthRate == b.adaptive.growthRate &&
           a.quality.enabled == b.quality.enabled &&
           a.quality.minimumAngleDegrees == b.quality.minimumAngleDegrees &&
           a.quality.maximumAspectRatio == b.quality.maximumAspectRatio &&
           actual.draftFaceTargetEdgeLengths == expected.draftFaceTargetEdgeLengths &&
           actual.cadFaceTargetEdgeLengths == expected.cadFaceTargetEdgeLengths &&
           actual.draftReversedFaceIds == expected.draftReversedFaceIds &&
           actual.cadReversedFaceIds == expected.cadReversedFaceIds &&
           project_json::encodeSimulation(actual.simulation)["settings"] ==
               project_json::encodeSimulation(expected.simulation)["settings"];
}
QString sweepManifestPath(const SweepStudy &study) {
    return QDir(study.directory).filePath("study.hsstudy");
}
SweepStudy createSweepStudy(const QString &parentDirectory, const Project &baseline,
                            std::vector<SweepAxis> axes, std::vector<double> frequencies) {
    if (!QFileInfo(parentDirectory).isDir())
        throw std::invalid_argument("Ein vorhandenes Ausgabeverzeichnis ist erforderlich.");
    if (baseline.source == ProjectSource::freecad &&
        (!baseline.cad || std::any_of(axes.begin(), axes.end(), [](const auto &a) {
            return a.parameter != "meshTargetEdgeLength";
        })))
        throw std::invalid_argument(
            "CAD-Studien unterstützen die Netzkonvergenz desselben Snapshots.");
    SweepStudy study;
    study.baseline = inputsOnly(baseline);
    study.axes = std::move(axes);
    study.variants = makeSweepGrid(study.axes);
    validateAxesForBaseline(study.baseline, study.axes);
    study.frequencyOverride = std::move(frequencies);
    validateConvergence(study);
    double previous = 0;
    for (const double frequency : study.frequencyOverride) {
        if (!std::isfinite(frequency) || frequency <= previous)
            invalid();
        previous = frequency;
    }
    const auto name = "Studie-" + QDateTime::currentDateTime().toString("yyyyMMdd-HHmmss") + "-" +
                      QUuid::createUuid().toString(QUuid::WithoutBraces);
    QDir parent(parentDirectory);
    if (!parent.mkdir(name))
        throw std::runtime_error("Studienverzeichnis konnte nicht angelegt werden.");
    study.directory = parent.absoluteFilePath(name);
    writeProject(QDir(study.directory).filePath("baseline.hsproj"), study.baseline);
    writeSweepStudy(study);
    return study;
}
void writeSweepStudy(const SweepStudy &study) {
    validateAxesForBaseline(study.baseline, study.axes);
    validateConvergence(study);
    const auto expected = makeSweepGrid(study.axes);
    if (expected.size() != study.variants.size() || !QFileInfo(study.directory).isDir())
        invalid();
    QJsonArray axes, variants;
    for (const auto &axis : study.axes)
        axes.append(QJsonObject{{"parameter", axis.parameter},
                                {"start", axis.start},
                                {"end", axis.end},
                                {"count", static_cast<double>(axis.count)}});
    for (std::size_t i = 0; i < study.variants.size(); ++i) {
        const auto &v = study.variants[i];
        if (v.values != expected[i].values || v.projectFile != fileName(i) || v.error.size() > 4096)
            invalid();
        QJsonArray values;
        for (const double value : v.values)
            values.append(value);
        variants.append(QJsonObject{{"values", values},
                                    {"state", stateName(v.state)},
                                    {"project", v.projectFile},
                                    {"error", v.error}});
    }
    QJsonArray frequencies;
    double previous = 0;
    for (const double frequency : study.frequencyOverride) {
        if (!std::isfinite(frequency) || frequency <= previous)
            invalid();
        frequencies.append(frequency);
        previous = frequency;
    }
    const auto data = QJsonDocument(QJsonObject{{"format", "HornSculptStudy"},
                                                {"version", 1},
                                                {"units", "SI"},
                                                {"baseline", "baseline.hsproj"},
                                                {"axes", axes},
                                                {"variants", variants},
                                                {"frequencies", frequencies}})
                          .toJson();
    QSaveFile file(sweepManifestPath(study));
    if (!file.open(QIODevice::WriteOnly) || file.write(data) != data.size() || !file.commit())
        throw std::runtime_error("Studienmanifest konnte nicht atomar gespeichert werden.");
}
SweepStudy readSweepStudy(const QString &manifestPath) {
    QFile file(manifestPath);
    if (!file.open(QIODevice::ReadOnly) || file.size() > 256 * 1024 * 1024)
        invalid();
    QJsonParseError error;
    const auto document = QJsonDocument::fromJson(file.readAll(), &error);
    if (error.error != QJsonParseError::NoError || !document.isObject())
        invalid();
    const auto root = document.object();
    keys(root, {"format", "version", "units", "baseline", "axes", "variants", "frequencies"});
    if (text(root["format"]) != "HornSculptStudy" || number(root["version"]) != 1 ||
        text(root["units"]) != "SI" || text(root["baseline"]) != "baseline.hsproj" ||
        !root["axes"].isArray() || !root["variants"].isArray())
        invalid();
    SweepStudy study;
    study.directory = QFileInfo(manifestPath).absolutePath();
    if (!root["frequencies"].isArray())
        invalid();
    double previous = 0;
    for (const auto value : root["frequencies"].toArray()) {
        const double frequency = number(value);
        if (frequency <= previous)
            invalid();
        study.frequencyOverride.push_back(frequency);
        previous = frequency;
    }
    for (const auto value : root["axes"].toArray()) {
        const auto a = object(value);
        keys(a, {"parameter", "start", "end", "count"});
        const double count = number(a["count"]);
        if (count < 1 || count != std::floor(count) || count > std::numeric_limits<int>::max())
            invalid();
        study.axes.push_back({text(a["parameter"]), number(a["start"]), number(a["end"]),
                              static_cast<std::size_t>(count)});
    }
    study.variants = makeSweepGrid(study.axes);
    const auto stored = root["variants"].toArray();
    if (stored.size() != static_cast<qsizetype>(study.variants.size()))
        invalid();
    for (qsizetype i = 0; i < stored.size(); ++i) {
        const auto v = object(stored[i]);
        keys(v, {"values", "state", "project", "error"});
        auto &variant = study.variants[static_cast<std::size_t>(i)];
        if (!v["values"].isArray() ||
            v["values"].toArray().size() != static_cast<qsizetype>(study.axes.size()) ||
            text(v["project"]) != variant.projectFile)
            invalid();
        const auto values = v["values"].toArray();
        for (qsizetype a = 0; a < values.size(); ++a)
            if (number(values[a]) != variant.values[static_cast<std::size_t>(a)])
                invalid();
        variant.state = stateValue(text(v["state"]));
        variant.error = text(v["error"]);
        if (variant.state == SweepState::completed &&
            !QFileInfo(QDir(study.directory).filePath(variant.projectFile)).isFile())
            invalid();
    }
    study.baseline = readProject(QDir(study.directory).filePath("baseline.hsproj"));
    if (study.baseline.source == ProjectSource::freecad &&
        (!study.baseline.cad ||
         std::any_of(study.axes.begin(), study.axes.end(),
                     [](const auto &a) { return a.parameter != "meshTargetEdgeLength"; })))
        invalid();
    validateAxesForBaseline(study.baseline, study.axes);
    validateConvergence(study);
    return study;
}
QString sweepVariantLabel(const SweepStudy &study, std::size_t index) {
    if (index >= study.variants.size() || study.variants[index].values.size() != study.axes.size())
        invalid();
    QStringList parts;
    for (std::size_t a = 0; a < study.axes.size(); ++a) {
        const auto p = presentation(binding(study.axes[a].parameter), study.baseline.parameters);
        parts << QString("%1=%2%3")
                     .arg(p.label)
                     .arg(study.variants[index].values[a] * p.displayScale, 0, 'g', 7)
                     .arg(p.unit.isEmpty() ? QString{} : " " + p.unit);
    }
    return parts.join("; ");
}
} // namespace hornsculpt
