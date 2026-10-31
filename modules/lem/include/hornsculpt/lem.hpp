#pragma once

#include <complex>
#include <cstdint>
#include <map>
#include <stop_token>
#include <string>
#include <variant>
#include <vector>

namespace hornsculpt::lem {

using Complex = std::complex<double>;

struct NodeId {
    std::uint32_t value{};
    bool operator==(const NodeId &) const = default;
};
struct ElementId {
    std::uint32_t value{};
    bool operator==(const ElementId &) const = default;
};
enum class Domain : std::uint8_t { electrical, mechanical, acoustic };
struct Position {
    double x{};
    double y{};
    bool operator==(const Position &) const = default;
};
struct Node {
    NodeId id;
    Domain domain{Domain::electrical};
    bool reference{};
    std::string name;
    Position position;
    bool operator==(const Node &) const = default;
};
struct Terminals {
    NodeId positive;
    NodeId negative;
    bool operator==(const Terminals &) const = default;
};

// Effort/flow are V/A, N/(m/s), and Pa/(m^3/s), respectively. All phasors
// are RMS with exp(+i omega t). Passive flow runs positive -> negative.
struct Resistance {
    Terminals terminals;
    double resistance{}; // Ohm, N s/m, or Pa s/m^3 according to the node domain.
    bool operator==(const Resistance &) const = default;
};
struct Inertance {
    Terminals terminals;
    double inertance{}; // H, kg, or Pa s^2/m^3 according to the node domain.
    bool operator==(const Inertance &) const = default;
};
struct Compliance {
    Terminals terminals;
    double compliance{}; // F, m/N, or m^3/Pa according to the node domain.
    bool operator==(const Compliance &) const = default;
};
struct EffortSource {
    Terminals terminals;
    Complex amplitude{1.0, 0.0}; // V, N, or Pa; positive-terminal effort rise.
    bool operator==(const EffortSource &) const = default;
};
struct FlowSource {
    Terminals terminals;
    Complex amplitude{1.0, 0.0}; // A, m/s, or m^3/s; positive -> negative flow.
    bool operator==(const FlowSource &) const = default;
};
struct Cavity {
    Terminals terminals;
    double volume{}; // m^3; acoustic compliance is V/(rho c^2), no thermal model.
    bool operator==(const Cavity &) const = default;
};
struct AcousticSlit {
    Terminals terminals;
    double length{}; // m; only the neck length outside the explicit BEM domain.
    double area{};   // m², total flow cross section of the parallel-plate slit.
    double gap{};    // m, full separation between the two stationary no-slip plates.
    double dynamicViscosity{1.84e-5}; // Pa s.
    // Lumped unsteady-Stokes series impedance with exp(+i omega t).
    // No end correction, thermal model, or distributed longitudinal waves.
    bool operator==(const AcousticSlit &) const = default;
};
struct ClosedSlitCavity {
    Terminals terminals;
    double length{}; // m; distributed 1D slit, rigidly closed at its remote end.
    double area{};   // m², constant cross section; the contained volume is area*length.
    double gap{};    // m, full plate separation.
    double dynamicViscosity{1.84e-5}; // Pa s.
    double heatCapacityRatio{1.4};
    double prandtlNumber{.71};
    // Input impedance of the thermoviscous plane mode; no azimuthal modes or
    // entrance/exit end correction. The terminals carry input pressure and Q.
    bool operator==(const ClosedSlitCavity &) const = default;
};

struct DriverParameters {
    double resistance{};  // Re, Ohm.
    double inductance{};  // Le, H; zero is allowed.
    double forceFactor{}; // Bl, N/A.
    double movingMass{};  // Mmd, kg, in vacuo: external acoustic inertia comes from the load.
    double compliance{};  // Cms, m/N.
    double mechanicalResistance{}; // Intrinsic Rms, N s/m, excluding external acoustic loss.
    double pistonArea{};           // Sd, m^2.
    bool operator==(const DriverParameters &) const = default;
};
struct Driver {
    Terminals electrical;
    NodeId front;
    NodeId rear;
    DriverParameters parameters;
    // Positive diaphragm velocity injects Sd*v into front and withdraws it
    // from rear. The pressure difference opposes positive diaphragm motion.
    bool operator==(const Driver &) const = default;
};
enum class BemPortKind : std::uint8_t { drivenSurface, interfacePort };
struct BemPort {
    Terminals terminals; // Acoustic positive node and acoustic reference node.
    std::vector<std::uint32_t> faceIds;
    // Q_BEM = polarity*Q_network. Physical Q enters the BEM air domain:
    // out of solid for drivenSurface, opposite the air mesh normal for interface.
    int polarity{1};
    BemPortKind kind{BemPortKind::drivenSurface};
    bool needsAssignmentReview{};
    // Positive dimensionless normal-velocity shape on drivenSurface faces.
    // Missing face IDs have weight 1. Interfaces remain spatially uniform.
    // Q=sum(A_f*v_f), S_w=sum(A_f*w_f), v_f=(w_f/S_w)*Q. The power-dual
    // port pressure is sum(w_f*integral_f(p))/S_w, not the unweighted mean.
    std::map<std::uint32_t, double> velocityWeights;
    bool operator==(const BemPort &) const = default;
};
using Component = std::variant<Resistance, Inertance, Compliance, EffortSource, FlowSource, Cavity,
                               Driver, BemPort, AcousticSlit, ClosedSlitCavity>;
struct Element {
    ElementId id;
    std::string name;
    Position position;
    Component component{Resistance{}};
    bool operator==(const Element &) const = default;
};
struct Network {
    std::vector<Node> nodes;
    std::vector<Element> elements;
    bool operator==(const Network &) const = default;
};
struct Medium {
    double density{1.21};     // kg/m^3.
    double soundSpeed{343.0}; // m/s.
};
struct AcousticLoad {
    // Row-major Z_ij = power-dual modal pressure at port i / volume velocity at
    // physical port j, Pa s/m^3. Ports follow BemPort order in Network::elements.
    // Includes all mutual terms. Physical flow enters the BEM air domain.
    std::vector<Complex> impedance;
    // Optional pressure generated by independent prescribed BEM excitations.
    // Empty means zero. The same physical port ordering is used.
    std::vector<Complex> blockedPressure;
};
struct NodeResult {
    NodeId id;
    Complex effort;
};
struct ElementResult {
    ElementId id;
    Complex effort;
    Complex flow; // Positive -> negative; driver uses electrical current.
    double
        activePower{}; // Terminal input W (electrical for Driver); delivering sources are negative.
};
struct DriverResult {
    ElementId id;
    Complex current;          // A RMS, electrical positive -> negative.
    Complex velocity;         // m/s RMS, toward the front acoustic node.
    Complex displacement;     // m RMS; multiply magnitude by sqrt(2) for peak.
    double electricalPower{}; // W into the electrical terminals.
    double copperPower{};     // W, |I|^2 Re.
    double mechanicalPower{}; // W dissipated in Rms.
    double acousticPower{};   // W into front and rear acoustic loads together.
};
struct PortResult {
    ElementId id;
    Complex pressure;       // Power-dual physical modal BEM pressure, Pa RMS.
    Complex volumeVelocity; // Physical Q entering BEM air, m^3/s RMS.
    double activePower{};   // W into the BEM domain.
};
// Dimensionless componentwise backward error bound for the original mixed-SI
// equations. Allows FP64 roundoff without accepting an unresolved equation.
inline constexpr double maximumRelativeResidual = 1e-10;
struct Solution {
    double frequency{};
    std::vector<NodeResult> nodes;
    std::vector<ElementResult> elements;
    std::vector<DriverResult> drivers;
    std::vector<PortResult> ports;
    // Maximum componentwise backward error of the unscaled MNA equations.
    double relativeResidual{};
};

// Invalid values, mixed domains, duplicate identities/faces, missing references,
// disconnected terminal groups, or malformed components throw invalid_argument.
// Singular/inconsistent circuits throw runtime_error during solve, never receive
// an implicit shunt, grounding connection, or alternate simulation path.
void validateNetwork(const Network &network);
// Simulation identity excludes display names and schematic positions. Ordering
// is retained because port/element result arrays follow network order. This
// comparison also accepts unfinished editor networks; it does not validate them.
bool samePhysics(const Network &first, const Network &second);
std::vector<const Element *> bemPorts(const Network &network);
Solution solve(const Network &network, double frequency, const AcousticLoad &load = {},
               const Medium &medium = {}, std::stop_token stop = {});

// Creates real voltage-driven transducers, separate sealed rear cavities, and
// one BEM port per face group. It contains no manufacturer presets or results.
Network makeDriverNetwork(const DriverParameters &parameters,
                          const std::vector<std::vector<std::uint32_t>> &faceGroups,
                          Complex voltageRms, double rearVolume);

} // namespace hornsculpt::lem
