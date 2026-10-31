#include "project_json.hpp"

#include <hornsculpt/lem.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <type_traits>
#include <utility>

namespace hornsculpt::project_json {
namespace {
constexpr qsizetype maximumArray = std::numeric_limits<int>::max();

void require(bool value) {
    if (!value)
        invalid(QStringLiteral("Inkonsistente LEM-Daten."));
}
std::uint32_t identity(const QJsonValue &value) {
    const auto result = integer(value);
    require(result != 0);
    return result;
}
QJsonArray encodePosition(lem::Position position) {
    return {position.x, position.y};
}
lem::Position decodePosition(const QJsonValue &value) {
    const auto values = array(value, 2, 2);
    return {number(values[0]), number(values[1])};
}
QJsonArray encodeTerminals(lem::Terminals terminals) {
    return {double(terminals.positive.value), double(terminals.negative.value)};
}
lem::Terminals decodeTerminals(const QJsonValue &value) {
    const auto values = array(value, 2, 2);
    return {{identity(values[0])}, {identity(values[1])}};
}
QString encodeDomain(lem::Domain domain) {
    switch (domain) {
    case lem::Domain::electrical:
        return QStringLiteral("electrical");
    case lem::Domain::mechanical:
        return QStringLiteral("mechanical");
    case lem::Domain::acoustic:
        return QStringLiteral("acoustic");
    }
    invalid(QStringLiteral("Unbekannte LEM-Domäne."));
}
lem::Domain decodeDomain(const QJsonValue &value) {
    const auto name = text(value);
    if (name == QStringLiteral("electrical"))
        return lem::Domain::electrical;
    if (name == QStringLiteral("mechanical"))
        return lem::Domain::mechanical;
    if (name == QStringLiteral("acoustic"))
        return lem::Domain::acoustic;
    invalid(QStringLiteral("Unbekannte LEM-Domäne."));
}
QJsonObject encodeParameters(const lem::DriverParameters &p) {
    return {{"resistance", p.resistance},   {"inductance", p.inductance},
            {"forceFactor", p.forceFactor}, {"movingMassInVacuo", p.movingMass},
            {"compliance", p.compliance},   {"mechanicalResistance", p.mechanicalResistance},
            {"pistonArea", p.pistonArea}};
}
lem::DriverParameters decodeParameters(const QJsonValue &value) {
    const auto o = object(value);
    keys(o, {"resistance", "inductance", "forceFactor", "movingMassInVacuo", "compliance",
             "mechanicalResistance", "pistonArea"});
    return {number(o["resistance"]),  number(o["inductance"]),
            number(o["forceFactor"]), number(o["movingMassInVacuo"]),
            number(o["compliance"]),  number(o["mechanicalResistance"]),
            number(o["pistonArea"])};
}
QJsonObject encodeComponent(const lem::Component &component) {
    return std::visit(
        [](const auto &part) -> QJsonObject {
            using Part = std::decay_t<decltype(part)>;
            if constexpr (std::is_same_v<Part, lem::Driver>) {
                return {{"type", "driver"},
                        {"electrical", encodeTerminals(part.electrical)},
                        {"front", double(part.front.value)},
                        {"rear", double(part.rear.value)},
                        {"parameters", encodeParameters(part.parameters)}};
            } else {
                QJsonObject result{{"terminals", encodeTerminals(part.terminals)}};
                if constexpr (std::is_same_v<Part, lem::Resistance>) {
                    result.insert("type", "resistance");
                    result.insert("resistance", part.resistance);
                } else if constexpr (std::is_same_v<Part, lem::Inertance>) {
                    result.insert("type", "inertance");
                    result.insert("inertance", part.inertance);
                } else if constexpr (std::is_same_v<Part, lem::Compliance>) {
                    result.insert("type", "compliance");
                    result.insert("compliance", part.compliance);
                } else if constexpr (std::is_same_v<Part, lem::Cavity>) {
                    result.insert("type", "cavity");
                    result.insert("volume", part.volume);
                } else if constexpr (std::is_same_v<Part, lem::AcousticSlit>) {
                    result.insert("type", "acousticSlit");
                    result.insert("length", part.length);
                    result.insert("area", part.area);
                    result.insert("gap", part.gap);
                    result.insert("dynamicViscosity", part.dynamicViscosity);
                } else if constexpr (std::is_same_v<Part, lem::ClosedSlitCavity>) {
                    result.insert("type", "closedSlitCavity");
                    result.insert("length", part.length);
                    result.insert("area", part.area);
                    result.insert("gap", part.gap);
                    result.insert("dynamicViscosity", part.dynamicViscosity);
                    result.insert("heatCapacityRatio", part.heatCapacityRatio);
                    result.insert("prandtlNumber", part.prandtlNumber);
                } else if constexpr (std::is_same_v<Part, lem::BemPort>) {
                    result.insert("type", "bemPort");
                    require(part.kind == lem::BemPortKind::drivenSurface ||
                            part.kind == lem::BemPortKind::interfacePort);
                    result.insert("kind", part.kind == lem::BemPortKind::interfacePort
                                              ? "interface"
                                              : "drivenSurface");
                    result.insert("needsAssignmentReview", part.needsAssignmentReview);
                    result.insert("polarity", part.polarity);
                    QJsonArray faces;
                    for (const auto face : part.faceIds)
                        faces.append(double(face));
                    result.insert("faceIds", faces);
                    QJsonArray weights;
                    require(part.velocityWeights.empty() ||
                            part.kind == lem::BemPortKind::drivenSurface);
                    for (const auto &[id, weight] : part.velocityWeights) {
                        require(std::find(part.faceIds.begin(), part.faceIds.end(), id) !=
                                    part.faceIds.end() &&
                                std::isfinite(weight) && weight > 0);
                        weights.append(QJsonArray{double(id), weight});
                    }
                    result.insert("velocityWeights", weights);
                } else {
                    result.insert("type", std::is_same_v<Part, lem::EffortSource> ? "effortSource"
                                                                                  : "flowSource");
                    result.insert("amplitude", encodeComplex(part.amplitude));
                }
                return result;
            }
        },
        component);
}
lem::Component decodeComponent(const QJsonValue &value) {
    const auto o = object(value);
    const auto type = text(o["type"]);
    if (type == QStringLiteral("driver")) {
        keys(o, {"type", "electrical", "front", "rear", "parameters"});
        return lem::Driver{decodeTerminals(o["electrical"]),
                           {identity(o["front"])},
                           {identity(o["rear"])},
                           decodeParameters(o["parameters"])};
    }
    const auto terminals = decodeTerminals(o["terminals"]);
    if (type == QStringLiteral("resistance")) {
        keys(o, {"type", "terminals", "resistance"});
        return lem::Resistance{terminals, number(o["resistance"])};
    }
    if (type == QStringLiteral("inertance")) {
        keys(o, {"type", "terminals", "inertance"});
        return lem::Inertance{terminals, number(o["inertance"])};
    }
    if (type == QStringLiteral("compliance")) {
        keys(o, {"type", "terminals", "compliance"});
        return lem::Compliance{terminals, number(o["compliance"])};
    }
    if (type == QStringLiteral("cavity")) {
        keys(o, {"type", "terminals", "volume"});
        return lem::Cavity{terminals, number(o["volume"])};
    }
    if (type == QStringLiteral("acousticSlit")) {
        keys(o, {"type", "terminals", "length", "area", "gap", "dynamicViscosity"});
        return lem::AcousticSlit{terminals, number(o["length"]), number(o["area"]),
                                 number(o["gap"]), number(o["dynamicViscosity"])};
    }
    if (type == QStringLiteral("closedSlitCavity")) {
        keys(o, {"type", "terminals", "length", "area", "gap", "dynamicViscosity",
                 "heatCapacityRatio", "prandtlNumber"});
        return lem::ClosedSlitCavity{terminals,
                                     number(o["length"]),
                                     number(o["area"]),
                                     number(o["gap"]),
                                     number(o["dynamicViscosity"]),
                                     number(o["heatCapacityRatio"]),
                                     number(o["prandtlNumber"])};
    }
    if (type == QStringLiteral("bemPort")) {
        auto required = o;
        required.remove("kind");
        required.remove("needsAssignmentReview");
        required.remove("velocityWeights");
        keys(required, {"type", "terminals", "faceIds", "polarity"});
        const auto polarity = number(o["polarity"]);
        require(polarity == 1 || polarity == -1);
        lem::BemPort result{terminals, {}, static_cast<int>(polarity)};
        if (o.contains("needsAssignmentReview"))
            result.needsAssignmentReview = boolean(o["needsAssignmentReview"]);
        if (o.contains("kind")) {
            const auto kind = text(o["kind"]);
            require(kind == "drivenSurface" || kind == "interface");
            result.kind = kind == "interface" ? lem::BemPortKind::interfacePort
                                              : lem::BemPortKind::drivenSurface;
        }
        for (const auto &face : array(o["faceIds"], maximumArray))
            result.faceIds.push_back(integer(face));
        if (o.contains("velocityWeights"))
            for (const auto &weightValue : array(o["velocityWeights"], maximumArray)) {
                const auto pair = array(weightValue, 2, 2);
                const auto id = identity(pair[0]);
                require(result.kind == lem::BemPortKind::drivenSurface &&
                        std::find(result.faceIds.begin(), result.faceIds.end(), id) !=
                            result.faceIds.end() &&
                        result.velocityWeights.emplace(id, positive(pair[1])).second);
            }
        return result;
    }
    if (type == QStringLiteral("effortSource") || type == QStringLiteral("flowSource")) {
        keys(o, {"type", "terminals", "amplitude"});
        const auto amplitude = complex(o["amplitude"]);
        if (type == QStringLiteral("effortSource"))
            return lem::EffortSource{terminals, amplitude};
        return lem::FlowSource{terminals, amplitude};
    }
    invalid(QStringLiteral("Unbekanntes LEM-Bauteil."));
}
bool close(lem::Complex a, lem::Complex b) {
    return std::isfinite(a.real()) && std::isfinite(a.imag()) && std::isfinite(b.real()) &&
           std::isfinite(b.imag()) &&
           std::abs(a - b) <= 1e-9 * std::max({std::abs(a), std::abs(b), 1e-300});
}
double power(lem::Complex effort, lem::Complex flow) {
    return (effort * std::conj(flow)).real();
}
} // namespace

QJsonObject encodeLemNetwork(const lem::Network &network) {
    QJsonArray nodes, elements;
    for (const auto &node : network.nodes)
        nodes.append(QJsonObject{{"id", double(node.id.value)},
                                 {"domain", encodeDomain(node.domain)},
                                 {"reference", node.reference},
                                 {"name", QString::fromStdString(node.name)},
                                 {"position", encodePosition(node.position)}});
    for (const auto &element : network.elements)
        elements.append(QJsonObject{{"id", double(element.id.value)},
                                    {"name", QString::fromStdString(element.name)},
                                    {"position", encodePosition(element.position)},
                                    {"component", encodeComponent(element.component)}});
    return {{"nodes", nodes}, {"elements", elements}};
}
lem::Network decodeLemNetwork(const QJsonValue &value) {
    const auto o = object(value);
    keys(o, {"nodes", "elements"});
    lem::Network result;
    std::map<std::uint32_t, lem::Domain> domains;
    std::set<std::uint32_t> elementIds;
    for (const auto &valueNode : array(o["nodes"], maximumArray)) {
        const auto n = object(valueNode);
        keys(n, {"id", "domain", "reference", "name", "position"});
        lem::Node node{{identity(n["id"])},
                       decodeDomain(n["domain"]),
                       boolean(n["reference"]),
                       text(n["name"], true).toStdString(),
                       decodePosition(n["position"])};
        require(domains.emplace(node.id.value, node.domain).second);
        result.nodes.push_back(std::move(node));
    }
    const auto domain = [&](lem::NodeId id) {
        require(domains.contains(id.value));
        return domains.at(id.value);
    };
    for (const auto &valueElement : array(o["elements"], maximumArray)) {
        const auto e = object(valueElement);
        keys(e, {"id", "name", "position", "component"});
        lem::Element element{{identity(e["id"])},
                             text(e["name"], true).toStdString(),
                             decodePosition(e["position"]),
                             decodeComponent(e["component"])};
        require(elementIds.insert(element.id.value).second);
        std::visit(
            [&](const auto &part) {
                using Part = std::decay_t<decltype(part)>;
                if constexpr (std::is_same_v<Part, lem::Driver>) {
                    require(domain(part.electrical.positive) == lem::Domain::electrical &&
                            domain(part.electrical.negative) == lem::Domain::electrical &&
                            domain(part.front) == lem::Domain::acoustic &&
                            domain(part.rear) == lem::Domain::acoustic);
                } else {
                    require(domain(part.terminals.positive) == domain(part.terminals.negative));
                    if constexpr (std::is_same_v<Part, lem::Cavity> ||
                                  std::is_same_v<Part, lem::BemPort> ||
                                  std::is_same_v<Part, lem::AcousticSlit> ||
                                  std::is_same_v<Part, lem::ClosedSlitCavity>)
                        require(domain(part.terminals.positive) == lem::Domain::acoustic);
                }
            },
            element.component);
        result.elements.push_back(std::move(element));
    }
    // An unfinished schematic must remain saveable: references, passive values,
    // port mappings, and solvability are checked when a simulation is requested.
    return result;
}

QJsonObject encodeLemSolution(const lem::Solution &solution) {
    QJsonArray nodes, elements, drivers, ports;
    for (const auto &node : solution.nodes)
        nodes.append(
            QJsonObject{{"id", double(node.id.value)}, {"effort", encodeComplex(node.effort)}});
    for (const auto &element : solution.elements)
        elements.append(QJsonObject{{"id", double(element.id.value)},
                                    {"effort", encodeComplex(element.effort)},
                                    {"flow", encodeComplex(element.flow)},
                                    {"activePower", element.activePower}});
    for (const auto &driver : solution.drivers)
        drivers.append(QJsonObject{{"id", double(driver.id.value)},
                                   {"current", encodeComplex(driver.current)},
                                   {"velocity", encodeComplex(driver.velocity)},
                                   {"displacement", encodeComplex(driver.displacement)},
                                   {"electricalPower", driver.electricalPower},
                                   {"copperPower", driver.copperPower},
                                   {"mechanicalPower", driver.mechanicalPower},
                                   {"acousticPower", driver.acousticPower}});
    for (const auto &port : solution.ports)
        ports.append(QJsonObject{{"id", double(port.id.value)},
                                 {"pressure", encodeComplex(port.pressure)},
                                 {"volumeVelocity", encodeComplex(port.volumeVelocity)},
                                 {"activePower", port.activePower}});
    return {{"frequency", solution.frequency},
            {"relativeResidual", solution.relativeResidual},
            {"nodes", nodes},
            {"elements", elements},
            {"drivers", drivers},
            {"ports", ports}};
}

lem::Solution decodeLemSolution(const QJsonValue &value, const lem::Network &network,
                                double frequency) {
    lem::validateNetwork(network);
    const auto o = object(value);
    keys(o, {"frequency", "relativeResidual", "nodes", "elements", "drivers", "ports"});
    lem::Solution result;
    result.frequency = positive(o["frequency"]);
    require(result.frequency == frequency);
    result.relativeResidual = number(o["relativeResidual"]);
    require(result.relativeResidual >= 0 &&
            result.relativeResidual <= lem::maximumRelativeResidual);
    std::map<std::uint32_t, const lem::Node *> knownNodes;
    std::map<std::uint32_t, const lem::Element *> knownElements, knownDrivers, knownPorts;
    for (const auto &node : network.nodes)
        knownNodes.emplace(node.id.value, &node);
    for (const auto &element : network.elements) {
        knownElements.emplace(element.id.value, &element);
        if (std::holds_alternative<lem::Driver>(element.component))
            knownDrivers.emplace(element.id.value, &element);
        if (std::holds_alternative<lem::BemPort>(element.component))
            knownPorts.emplace(element.id.value, &element);
    }
    std::map<std::uint32_t, lem::Complex> efforts;
    std::map<std::uint32_t, lem::ElementResult> elementResults;
    for (const auto &valueNode : array(o["nodes"], knownNodes.size(), knownNodes.size())) {
        const auto n = object(valueNode);
        keys(n, {"id", "effort"});
        lem::NodeResult node{{identity(n["id"])}, complex(n["effort"])};
        require(knownNodes.contains(node.id.value) &&
                efforts.emplace(node.id.value, node.effort).second);
        require(!knownNodes.at(node.id.value)->reference || node.effort == lem::Complex{});
        result.nodes.push_back(node);
    }
    for (const auto &valueElement :
         array(o["elements"], knownElements.size(), knownElements.size())) {
        const auto e = object(valueElement);
        keys(e, {"id", "effort", "flow", "activePower"});
        lem::ElementResult element{{identity(e["id"])},
                                   complex(e["effort"]),
                                   complex(e["flow"]),
                                   number(e["activePower"])};
        require(knownElements.contains(element.id.value) &&
                elementResults.emplace(element.id.value, element).second);
        const auto &component = knownElements.at(element.id.value)->component;
        std::visit(
            [&](const auto &part) {
                using Part = std::decay_t<decltype(part)>;
                if constexpr (std::is_same_v<Part, lem::Driver>)
                    require(close(element.effort, efforts.at(part.electrical.positive.value) -
                                                      efforts.at(part.electrical.negative.value)));
                else {
                    require(close(element.effort, efforts.at(part.terminals.positive.value) -
                                                      efforts.at(part.terminals.negative.value)));
                    if constexpr (std::is_same_v<Part, lem::FlowSource>)
                        require(close(element.flow, part.amplitude));
                    if constexpr (std::is_same_v<Part, lem::EffortSource>) {
                        // A zero differential source can sit at a large common-mode
                        // effort. Test its equation with the original node scale.
                        const double scale = std::abs(efforts.at(part.terminals.positive.value)) +
                                             std::abs(efforts.at(part.terminals.negative.value)) +
                                             std::abs(part.amplitude);
                        require(std::isfinite(scale) && std::abs(element.effort - part.amplitude) <=
                                                            1e-9 * std::max(scale, 1e-300));
                    }
                }
            },
            component);
        require(close(element.activePower, power(element.effort, element.flow)));
        result.elements.push_back(element);
    }
    std::set<std::uint32_t> seenDrivers, seenPorts;
    const double omega = 2 * std::numbers::pi * frequency;
    require(std::isfinite(omega) && omega > 0);
    for (const auto &valueDriver : array(o["drivers"], knownDrivers.size(), knownDrivers.size())) {
        const auto d = object(valueDriver);
        keys(d, {"id", "current", "velocity", "displacement", "electricalPower", "copperPower",
                 "mechanicalPower", "acousticPower"});
        lem::DriverResult driver{{identity(d["id"])},          complex(d["current"]),
                                 complex(d["velocity"]),       complex(d["displacement"]),
                                 number(d["electricalPower"]), number(d["copperPower"]),
                                 number(d["mechanicalPower"]), number(d["acousticPower"])};
        require(knownDrivers.contains(driver.id.value) &&
                seenDrivers.insert(driver.id.value).second);
        const auto &element = elementResults.at(driver.id.value);
        const auto &part = std::get<lem::Driver>(knownDrivers.at(driver.id.value)->component);
        const auto &p = part.parameters;
        require(
            close(driver.current, element.flow) &&
            close(driver.displacement,
                  {driver.velocity.imag() / omega, -driver.velocity.real() / omega}) &&
            close(driver.electricalPower, element.activePower) &&
            close(driver.copperPower, std::norm(driver.current) * p.resistance) &&
            close(driver.mechanicalPower, std::norm(driver.velocity) * p.mechanicalResistance) &&
            close(driver.acousticPower,
                  power(efforts.at(part.front.value) - efforts.at(part.rear.value),
                        p.pistonArea * driver.velocity)));
        result.drivers.push_back(driver);
    }
    for (const auto &valuePort : array(o["ports"], knownPorts.size(), knownPorts.size())) {
        const auto p = object(valuePort);
        keys(p, {"id", "pressure", "volumeVelocity", "activePower"});
        lem::PortResult port{{identity(p["id"])},
                             complex(p["pressure"]),
                             complex(p["volumeVelocity"]),
                             number(p["activePower"])};
        require(knownPorts.contains(port.id.value) && seenPorts.insert(port.id.value).second);
        const auto &element = elementResults.at(port.id.value);
        const double sign =
            std::get<lem::BemPort>(knownPorts.at(port.id.value)->component).polarity;
        require(close(port.pressure, sign * element.effort) &&
                close(port.volumeVelocity, sign * element.flow) &&
                close(port.activePower, element.activePower));
        result.ports.push_back(port);
    }
    return result;
}
} // namespace hornsculpt::project_json
