#include <hornsculpt/lem.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace hornsculpt::lem {
namespace {

constexpr std::size_t noIndex = std::numeric_limits<std::size_t>::max();
constexpr Complex imaginary{0.0, 1.0};
bool finite(Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}
void require(bool condition, const char *message) {
    if (!condition)
        throw std::invalid_argument(message);
}
void positive(double value, const char *message) {
    require(std::isfinite(value) && value > 0.0, message);
}
void nonnegative(double value, const char *message) {
    require(std::isfinite(value) && value >= 0.0, message);
}
void checkStop(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("LEM-Berechnung abgebrochen.");
}
void validateParameters(const DriverParameters &p) {
    positive(p.resistance, "Treiber: Re muss endlich und positiv sein [Ohm].");
    nonnegative(p.inductance, "Treiber: Le muss endlich und nichtnegativ sein [H].");
    positive(p.forceFactor, "Treiber: Bl muss endlich und positiv sein [N/A].");
    positive(p.movingMass,
             "Treiber: Mmd ohne externe Luftlast muss endlich und positiv sein [kg].");
    positive(p.compliance, "Treiber: Cms muss endlich und positiv sein [m/N].");
    nonnegative(p.mechanicalResistance, "Treiber: Rms muss endlich und nichtnegativ sein [Ns/m].");
    positive(p.pistonArea, "Treiber: Sd muss endlich und positiv sein [m^2].");
}

class Topology {
  public:
    explicit Topology(const Network &network) : network_(network), links_(network.nodes.size()) {
        require(!network.nodes.empty() && !network.elements.empty(), "Das LEM-Netzwerk ist leer.");
        for (std::size_t i = 0; i < network.nodes.size(); ++i) {
            const auto &node = network.nodes[i];
            require(node.id.value != 0, "LEM-Knoten benötigen eine positive Identität.");
            require(indices_.emplace(node.id.value, i).second, "Doppelte LEM-Knotenidentität.");
            require(node.domain == Domain::electrical || node.domain == Domain::mechanical ||
                        node.domain == Domain::acoustic,
                    "Unbekannte LEM-Knotendomäne.");
            require(std::isfinite(node.position.x) && std::isfinite(node.position.y),
                    "Die LEM-Knotenposition ist nicht endlich.");
        }
        std::unordered_set<std::uint32_t> elementIds;
        std::unordered_set<std::uint32_t> faces;
        for (const auto &element : network.elements) {
            require(element.id.value != 0 && elementIds.insert(element.id.value).second,
                    "Ungültige oder doppelte LEM-Bauteilidentität.");
            require(std::isfinite(element.position.x) && std::isfinite(element.position.y),
                    "Die LEM-Bauteilposition ist nicht endlich.");
            std::visit(
                [&](const auto &part) {
                    using Part = std::decay_t<decltype(part)>;
                    if constexpr (std::is_same_v<Part, Driver>) {
                        validateParameters(part.parameters);
                        terminals(part.electrical, Domain::electrical);
                        node(part.front, Domain::acoustic);
                        node(part.rear, Domain::acoustic);
                        link(part.front, part.rear);
                    } else {
                        terminals(part.terminals);
                        if constexpr (std::is_same_v<Part, Resistance>)
                            positive(part.resistance,
                                     "Ein LEM-Widerstand muss endlich und positiv sein.");
                        else if constexpr (std::is_same_v<Part, Inertance>)
                            positive(part.inertance,
                                     "Eine LEM-Induktivität/Masse muss endlich und positiv sein.");
                        else if constexpr (std::is_same_v<Part, Compliance>)
                            positive(
                                part.compliance,
                                "Eine LEM-Kapazität/Nachgiebigkeit muss endlich und positiv sein.");
                        else if constexpr (std::is_same_v<Part, EffortSource> ||
                                           std::is_same_v<Part, FlowSource>) {
                            require(finite(part.amplitude),
                                    "Die komplexe LEM-Anregung ist nicht endlich.");
                            if constexpr (std::is_same_v<Part, EffortSource>)
                                require(part.terminals.positive != part.terminals.negative,
                                        "Eine ideale LEM-Potentialquelle darf nicht "
                                        "kurzgeschlossen sein.");
                        } else if constexpr (std::is_same_v<Part, Cavity>) {
                            terminals(part.terminals, Domain::acoustic);
                            positive(part.volume,
                                     "Das Kammervolumen muss endlich und positiv sein [m^3].");
                        } else if constexpr (std::is_same_v<Part, AcousticSlit> ||
                                             std::is_same_v<Part, ClosedSlitCavity>) {
                            terminals(part.terminals, Domain::acoustic);
                            positive(
                                part.length,
                                "Die akustische Spaltlänge muss endlich und positiv sein [m].");
                            positive(
                                part.area,
                                "Die akustische Spaltfläche muss endlich und positiv sein [m²].");
                            positive(part.gap, "Die volle akustische Spalthöhe muss endlich und "
                                               "positiv sein [m].");
                            positive(
                                part.dynamicViscosity,
                                "Die dynamische Viskosität muss endlich und positiv sein [Pa s].");
                            if constexpr (std::is_same_v<Part, ClosedSlitCavity>) {
                                require(std::isfinite(part.heatCapacityRatio) &&
                                            part.heatCapacityRatio >= 1,
                                        "Das Verhältnis der Wärmekapazitäten muss endlich und "
                                        "mindestens 1 sein.");
                                positive(part.prandtlNumber,
                                         "Die Prandtl-Zahl muss endlich und positiv sein.");
                            }
                        } else if constexpr (std::is_same_v<Part, BemPort>) {
                            terminals(part.terminals, Domain::acoustic);
                            require(!part.needsAssignmentReview,
                                    "Die BEM-Flächenzuordnung muss nach einem Geometriewechsel "
                                    "bestätigt werden.");
                            require(at(part.terminals.negative).reference,
                                    "Der negative BEM-Anschluss muss akustische Referenz sein.");
                            require(!at(part.terminals.positive).reference,
                                    "Der positive BEM-Anschluss darf nicht auf Referenz liegen.");
                            require(part.polarity == 1 || part.polarity == -1,
                                    "Die BEM-Anschlusspolarität muss +1 oder -1 sein.");
                            require(part.kind == BemPortKind::drivenSurface ||
                                        part.kind == BemPortKind::interfacePort,
                                    "Unbekannter Typ des BEM-Anschlusses.");
                            require(!part.faceIds.empty(),
                                    "Ein BEM-Anschluss benötigt mindestens eine Fläche.");
                            for (const auto face : part.faceIds)
                                require(
                                    faces.insert(face).second,
                                    "Eine BEM-Fläche ist mehrfach einem LEM-Anschluss zugeordnet.");
                            require(
                                part.kind == BemPortKind::drivenSurface ||
                                    part.velocityWeights.empty(),
                                "Virtuelle BEM-Interfaces besitzen eine uniforme Geschwindigkeit.");
                            for (const auto &[face, weight] : part.velocityWeights) {
                                require(std::find(part.faceIds.begin(), part.faceIds.end(), face) !=
                                            part.faceIds.end(),
                                        "Ein BEM-Geschwindigkeitsgewicht verweist auf eine fremde "
                                        "Fläche.");
                                positive(weight, "BEM-Geschwindigkeitsgewichte müssen endlich und "
                                                 "positiv sein.");
                            }
                        }
                    }
                },
                element.component);
        }
        // Absolute effort needs an explicit reference in every connected
        // same-domain component. No hidden electrical/acoustic ground is added.
        std::vector<bool> reached(network.nodes.size(), false);
        std::vector<std::size_t> pending;
        for (std::size_t i = 0; i < network.nodes.size(); ++i)
            if (network.nodes[i].reference) {
                reached[i] = true;
                pending.push_back(i);
            }
        while (!pending.empty()) {
            const auto current = pending.back();
            pending.pop_back();
            for (const auto adjacent : links_[current])
                if (!reached[adjacent]) {
                    reached[adjacent] = true;
                    pending.push_back(adjacent);
                }
        }
        require(std::all_of(reached.begin(), reached.end(), [](bool value) { return value; }),
                "Ein LEM-Knoten besitzt keine Verbindung zu einer Referenz seiner Domäne.");
    }

    const Node &at(NodeId id) const {
        const auto found = indices_.find(id.value);
        require(found != indices_.end(), "Ein LEM-Anschluss verweist auf einen fehlenden Knoten.");
        return network_.nodes[found->second];
    }
    std::size_t index(NodeId id) const {
        (void)at(id);
        return indices_.at(id.value);
    }

  private:
    void node(NodeId id, Domain domain) const {
        require(at(id).domain == domain,
                "Ein LEM-Anschluss ist mit der falschen Domäne verbunden.");
    }
    void terminals(Terminals terminals) {
        require(at(terminals.positive).domain == at(terminals.negative).domain,
                "Ein LEM-Bauteil verbindet unterschiedliche Domänen.");
        link(terminals.positive, terminals.negative);
    }
    void terminals(Terminals terminals, Domain domain) {
        node(terminals.positive, domain);
        node(terminals.negative, domain);
        link(terminals.positive, terminals.negative);
    }
    void link(NodeId a, NodeId b) {
        const auto i = index(a), j = index(b);
        links_[i].push_back(j);
        links_[j].push_back(i);
    }
    const Network &network_;
    std::unordered_map<std::uint32_t, std::size_t> indices_;
    std::vector<std::vector<std::size_t>> links_;
};

class DenseSystem {
  public:
    explicit DenseSystem(std::size_t size) : size_(size), rhs(size) {
        require(size > 0 && size <= std::vector<Complex>().max_size() / size,
                "Die Größe des LEM-Gleichungssystems ist ungültig.");
        matrix.resize(size * size);
    }
    Complex &at(std::size_t row, std::size_t column) {
        return matrix[row * size_ + column];
    }
    std::vector<Complex> solve(std::stop_token stop, double &backwardError) const {
        checkStop(stop);
        std::vector<Complex> lu(matrix);
        std::vector<double> rowScale(size_), columnScale(size_);
        std::vector<std::size_t> pivots(size_);
        for (std::size_t row = 0; row < size_; ++row) {
            double largest = 0;
            for (std::size_t column = 0; column < size_; ++column) {
                if (!finite(lu[row * size_ + column]))
                    throw std::runtime_error(
                        "Die LEM-Gleichungen enthalten nichtendliche Koeffizienten.");
                largest = std::max(largest, std::abs(lu[row * size_ + column]));
            }
            if (largest == 0 || !finite(rhs[row]))
                throw std::runtime_error(
                    "Das LEM-Netzwerk enthält eine unbestimmte oder ungültige Gleichung.");
            rowScale[row] = 1 / largest;
            if (!std::isfinite(rowScale[row]))
                throw std::runtime_error(
                    "Die LEM-Zeilenskalierung überschreitet den FP64-Zahlenbereich.");
            for (std::size_t column = 0; column < size_; ++column)
                lu[row * size_ + column] *= rowScale[row];
        }
        for (std::size_t column = 0; column < size_; ++column) {
            double largest = 0;
            for (std::size_t row = 0; row < size_; ++row)
                largest = std::max(largest, std::abs(lu[row * size_ + column]));
            if (largest == 0)
                throw std::runtime_error(
                    "Das LEM-Netzwerk besitzt einen unbestimmten Freiheitsgrad.");
            columnScale[column] = 1 / largest;
            if (!std::isfinite(columnScale[column]))
                throw std::runtime_error(
                    "Die LEM-Spaltenskalierung überschreitet den FP64-Zahlenbereich.");
            for (std::size_t row = 0; row < size_; ++row)
                lu[row * size_ + column] *= columnScale[column];
        }
        const double pivotFloor =
            32 * std::numeric_limits<double>::epsilon() * static_cast<double>(size_);
        for (std::size_t column = 0; column < size_; ++column) {
            checkStop(stop);
            std::size_t pivot = column;
            for (std::size_t row = column + 1; row < size_; ++row)
                if (std::abs(lu[row * size_ + column]) > std::abs(lu[pivot * size_ + column]))
                    pivot = row;
            if (!finite(lu[pivot * size_ + column]) ||
                std::abs(lu[pivot * size_ + column]) <= pivotFloor)
                throw std::runtime_error(
                    "Das LEM-Netzwerk ist singulär oder numerisch nicht eindeutig lösbar.");
            pivots[column] = pivot;
            if (pivot != column)
                for (std::size_t j = 0; j < size_; ++j)
                    std::swap(lu[pivot * size_ + j], lu[column * size_ + j]);
            for (std::size_t row = column + 1; row < size_; ++row) {
                lu[row * size_ + column] /= lu[column * size_ + column];
                for (std::size_t j = column + 1; j < size_; ++j)
                    lu[row * size_ + j] -= lu[row * size_ + column] * lu[column * size_ + j];
            }
        }
        const auto apply = [&](const std::vector<Complex> &right) {
            std::vector<Complex> result(size_);
            for (std::size_t i = 0; i < size_; ++i)
                result[i] = right[i] * rowScale[i];
            for (std::size_t i = 0; i < size_; ++i)
                if (pivots[i] != i)
                    std::swap(result[i], result[pivots[i]]);
            for (std::size_t i = 0; i < size_; ++i)
                for (std::size_t j = 0; j < i; ++j)
                    result[i] -= lu[i * size_ + j] * result[j];
            for (std::size_t i = size_; i-- > 0;) {
                for (std::size_t j = i + 1; j < size_; ++j)
                    result[i] -= lu[i * size_ + j] * result[j];
                result[i] /= lu[i * size_ + i];
            }
            for (std::size_t i = 0; i < size_; ++i) {
                result[i] *= columnScale[i];
                if (!finite(result[i]))
                    throw std::runtime_error("Die LEM-Lösung enthält nichtendliche Werte.");
            }
            return result;
        };
        auto solution = apply(rhs);
        for (unsigned correction = 0; correction <= 2; ++correction) {
            checkStop(stop);
            std::vector<Complex> residual(size_);
            backwardError = 0;
            for (std::size_t row = 0; row < size_; ++row) {
                std::complex<long double> remainder{rhs[row].real(), rhs[row].imag()};
                long double scale = std::abs(remainder);
                for (std::size_t column = 0; column < size_; ++column) {
                    const auto a = matrix[row * size_ + column], x = solution[column];
                    const std::complex<long double> al{a.real(), a.imag()}, xl{x.real(), x.imag()};
                    remainder -= al * xl;
                    scale += std::abs(al) * std::abs(xl);
                }
                residual[row] = {static_cast<double>(remainder.real()),
                                 static_cast<double>(remainder.imag())};
                if (!finite(residual[row]) || !std::isfinite(scale))
                    throw std::runtime_error(
                        "Der LEM-Gleichungsrest überschreitet den FP64-Zahlenbereich.");
                if (scale > 0)
                    backwardError =
                        std::max(backwardError, static_cast<double>(std::abs(remainder) / scale));
            }
            if (std::isfinite(backwardError) && backwardError <= maximumRelativeResidual)
                return solution;
            if (correction != 2) {
                const auto delta = apply(residual);
                for (std::size_t i = 0; i < size_; ++i)
                    solution[i] += delta[i];
            }
        }
        throw std::runtime_error("Die LEM-Lösung erfüllt die Gleichungsrestgrenze nicht.");
    }
    std::size_t size_;
    std::vector<Complex> matrix;
    std::vector<Complex> rhs;
};

Complex stokesSeries(double gap, double viscosity, double omega, double density,
                     double logLengthArea) {
    // Fully developed no-slip parallel-plate flow, positive time convention.
    // Christensen, COMSOL (2024), "Nonrational Transfer Function: Microacoustics":
    // https://www.comsol.com/blogs?p=463571 . Stinson (1991), Eqs. (47),(83)-(86),
    // https://doi.org/10.1121/1.400379 . Full gap h=2*l_y, chi²=i*omega*rho*h²/(4*mu).
    // Logs avoid overflowing intermediate dimensional products when the final
    // impedance and dimensionless shape remain representable in FP64.
    const double logGap = std::log(gap), logOmega = std::log(omega), logDensity = std::log(density),
                 logViscosity = std::log(viscosity),
                 logChi = logGap - std::log(2.0) + .5 * (logOmega + logDensity - logViscosity);
    Complex result;
    if (logChi < std::log(.1)) {
        const Complex z{0, std::exp(2 * logChi)};
        const double resistance =
            std::exp(std::log(12.0) + logViscosity + logLengthArea - 2 * logGap);
        // chi²/[3*(1-tanh(chi)/chi)] expanded in z=chi². At |chi|<=0.1
        // the first omitted term is below 4.1e-14 relative; direct subtraction
        // instead loses all low-frequency resistance when chi is tiny.
        result = resistance * (1.0 + z * (2.0 / 5 + z * (-1.0 / 525 + z * (2.0 / 23625))));
    } else {
        const double inverseComponent = std::exp(-logChi) / std::sqrt(2.0);
        const Complex inverseChi{inverseComponent, -inverseComponent};
        Complex hyperbolic{1, 0};
        // With Re(chi)>20, |tanh(chi)-1|<2*exp(-40), below FP64 roundoff.
        if (logChi <= std::log(20 * std::sqrt(2.0))) {
            const double component = std::exp(logChi) / std::sqrt(2.0);
            hyperbolic = std::tanh(Complex{component, component});
        }
        const double omegaMass = std::exp(logOmega + logDensity + logLengthArea);
        result = Complex{0, omegaMass} / (1.0 - hyperbolic * inverseChi);
    }
    if (!finite(result) || !(std::abs(result) > 0) || result.real() < 0 || result.imag() < 0)
        throw std::runtime_error(
            "Die passive akustische Spaltimpedanz ist im FP64-Zahlenbereich nicht darstellbar.");
    return result;
}
Complex slitImpedance(const AcousticSlit &slit, double omega, double density) {
    return stokesSeries(slit.gap, slit.dynamicViscosity, omega, density,
                        std::log(slit.length) - std::log(slit.area));
}
Complex tanhRatio(double logChi) {
    if (logChi < std::log(.1)) {
        const Complex z{0, std::exp(2 * logChi)};
        return 1.0 +
               z * (-1.0 / 3 +
                    z * (2.0 / 15 + z * (-17.0 / 315 + z * (62.0 / 2835 - z * (1382.0 / 155925)))));
    }
    const double inverseComponent = std::exp(-logChi) / std::sqrt(2.0);
    const Complex inverse{inverseComponent, -inverseComponent};
    if (logChi > std::log(20 * std::sqrt(2.0)))
        return inverse;
    const double component = std::exp(logChi) / std::sqrt(2.0);
    return std::tanh(Complex{component, component}) * inverse;
}
Complex closedSlitImpedance(const ClosedSlitCavity &slit, double omega, const Medium &medium) {
    // Stinson (1991), Eqs.(65),(66),(83)-(87); COMSOL equivalent-fluid LRF theory:
    // https://doc.comsol.com/6.4/doc/com.comsol.help.aco/aco_ug_pressure.05.186.html
    // rho_eff=rho/(1-Fv), beta_eff=(1+(gamma-1)*Ft)/(rho*c²).
    // The closed transmission line has Zin=-i*Zc*cot(k*L), exp(+i omega t).
    const auto seriesDensity =
        stokesSeries(slit.gap, slit.dynamicViscosity, omega, medium.density, 0);
    const auto density = seriesDensity / Complex{0, omega};
    const double thermalLogChi =
        std::log(slit.gap) - std::log(2.0) +
        .5 * (std::log(omega) + std::log(medium.density) - std::log(slit.dynamicViscosity) +
              std::log(slit.prandtlNumber));
    const Complex beta = (1.0 + (slit.heatCapacityRatio - 1) * tanhRatio(thermalLogChi)) /
                         (medium.density * medium.soundSpeed * medium.soundSpeed);
    const auto k = omega * std::sqrt(density * beta);
    // Derive Zc from the SAME k branch, rather than selecting an independent root.
    const auto characteristic = k / (omega * beta * slit.area);
    const auto x = k * slit.length;
    if (!finite(density) || !finite(beta) || !finite(k) || !finite(characteristic) || !finite(x) ||
        !(k.real() >= 0) || !(k.imag() <= 0) || !(characteristic.real() > 0))
        throw std::runtime_error(
            "Die passive geschlossene Spaltleitung besitzt keinen darstellbaren Wellenzweig.");
    Complex result;
    if (std::abs(x) < .1) {
        const auto z = x * x;
        // x*cot(x), avoiding the singular cotangent and its cancellation with k.
        const auto factor = 1.0 + z * (-1.0 / 3 + z * (-1.0 / 45 + z * (-2.0 / 945 - z / 4725.0)));
        result = factor / (Complex{0, omega} * beta * slit.area * slit.length);
    } else if (x.imag() < -20) {
        // cot(x)->i in the lower half-plane, error below FP64 at this bound.
        result = characteristic;
    } else {
        result = Complex{0, -1} * characteristic / std::tan(x);
    }
    if (!finite(result) || !(std::abs(result) > 0) || result.real() < 0)
        throw std::runtime_error(
            "Die passive Eingangsimpedanz der geschlossenen Spaltleitung ist nicht darstellbar.");
    return result;
}
Complex impedance(const Component &component, double omega, const Medium &medium) {
    return std::visit(
        [&](const auto &part) -> Complex {
            using Part = std::decay_t<decltype(part)>;
            if constexpr (std::is_same_v<Part, Resistance>)
                return part.resistance;
            else if constexpr (std::is_same_v<Part, Inertance>)
                return imaginary * omega * part.inertance;
            else if constexpr (std::is_same_v<Part, Compliance>)
                return {0, -1.0 / (omega * part.compliance)};
            else if constexpr (std::is_same_v<Part, Cavity>)
                return {0, -medium.density * medium.soundSpeed * medium.soundSpeed /
                               (omega * part.volume)};
            else if constexpr (std::is_same_v<Part, AcousticSlit>)
                return slitImpedance(part, omega, medium.density);
            else if constexpr (std::is_same_v<Part, ClosedSlitCavity>)
                return closedSlitImpedance(part, omega, medium);
            else
                return 0.0;
        },
        component);
}
double activePower(Complex effort, Complex flow) {
    return (effort * std::conj(flow)).real();
}
} // namespace

void validateNetwork(const Network &network) {
    (void)Topology(network);
}

bool samePhysics(const Network &first, const Network &second) {
    if (first.nodes.size() != second.nodes.size() ||
        first.elements.size() != second.elements.size())
        return false;
    for (std::size_t i = 0; i < first.nodes.size(); ++i) {
        const auto &a = first.nodes[i];
        const auto &b = second.nodes[i];
        if (a.id != b.id || a.domain != b.domain || a.reference != b.reference)
            return false;
    }
    for (std::size_t i = 0; i < first.elements.size(); ++i) {
        const auto &a = first.elements[i];
        const auto &b = second.elements[i];
        if (a.id != b.id || a.component != b.component)
            return false;
    }
    return true;
}

std::vector<const Element *> bemPorts(const Network &network) {
    std::vector<const Element *> result;
    for (const auto &element : network.elements)
        if (std::holds_alternative<BemPort>(element.component))
            result.push_back(&element);
    return result;
}

Solution solve(const Network &network, double frequency, const AcousticLoad &load,
               const Medium &medium, std::stop_token stop) {
    checkStop(stop);
    positive(frequency, "Die LEM-Frequenz muss endlich und positiv sein [Hz].");
    positive(medium.density, "Die Luftdichte muss endlich und positiv sein [kg/m^3].");
    positive(medium.soundSpeed, "Die Schallgeschwindigkeit muss endlich und positiv sein [m/s].");
    const Topology topology(network);
    const auto ports = bemPorts(network);
    const auto portCount = ports.size();
    require(portCount == 0 || portCount <= load.impedance.max_size() / portCount,
            "Zu viele akustische LEM-Anschlüsse.");
    require(load.impedance.size() == portCount * portCount,
            "Die BEM-Impedanzmatrix passt nicht zur Anzahl der LEM-Anschlüsse.");
    require(load.blockedPressure.empty() || load.blockedPressure.size() == portCount,
            "Der BEM-Sperrdruck passt nicht zur Anzahl der LEM-Anschlüsse.");
    require(std::all_of(load.impedance.begin(), load.impedance.end(), finite) &&
                std::all_of(load.blockedPressure.begin(), load.blockedPressure.end(), finite),
            "Die akustische BEM-Last enthält nichtendliche Werte.");
    const double omega = 2 * std::numbers::pi * frequency;
    require(std::isfinite(omega), "Die LEM-Kreisfrequenz ist nicht endlich.");
    std::vector<std::size_t> nodeDofs(network.nodes.size(), noIndex),
        branches(network.elements.size(), noIndex);
    std::size_t size = 0;
    for (std::size_t i = 0; i < network.nodes.size(); ++i)
        if (!network.nodes[i].reference)
            nodeDofs[i] = size++;
    std::unordered_map<std::uint32_t, std::size_t> portBranches;
    for (std::size_t i = 0; i < network.elements.size(); ++i) {
        const auto &component = network.elements[i].component;
        if (!std::holds_alternative<FlowSource>(component)) {
            branches[i] = size;
            size += std::holds_alternative<Driver>(component) ? 2 : 1;
            if (std::holds_alternative<BemPort>(component))
                portBranches.emplace(network.elements[i].id.value, branches[i]);
        }
    }
    DenseSystem system(size);
    const auto nodeDof = [&](NodeId id) { return nodeDofs[topology.index(id)]; };
    const auto stampFlow = [&](NodeId id, std::size_t branch, double coefficient) {
        const auto row = nodeDof(id);
        if (row != noIndex)
            system.at(row, branch) += coefficient;
    };
    const auto stampEffort = [&](std::size_t row, NodeId id, double coefficient) {
        const auto column = nodeDof(id);
        if (column != noIndex)
            system.at(row, column) += coefficient;
    };
    std::size_t portIndex = 0;
    for (std::size_t index = 0; index < network.elements.size(); ++index) {
        checkStop(stop);
        const auto &element = network.elements[index];
        const auto branch = branches[index];
        std::visit(
            [&](const auto &part) {
                using Part = std::decay_t<decltype(part)>;
                if constexpr (std::is_same_v<Part, Driver>) {
                    // Impedance analogy: V = Ze*I + Bl*v and Bl*I = Zm*v + Sd*(pf-pr).
                    // The opposing pressure force and back EMF conserve terminal power.
                    // COMSOL, "Modeling Speaker Drivers - Lumped Methods", Methods 2-3:
                    // https://www.comsol.com/support/learning-center/course/modeling-speaker-drivers-in-comsol-multiphysics-202/modeling-speaker-drivers-lumped-methods-88401
                    const auto &p = part.parameters;
                    const auto velocity = branch + 1;
                    stampFlow(part.electrical.positive, branch, 1);
                    stampFlow(part.electrical.negative, branch, -1);
                    stampFlow(part.front, velocity, -p.pistonArea);
                    stampFlow(part.rear, velocity, p.pistonArea);
                    stampEffort(branch, part.electrical.positive, 1);
                    stampEffort(branch, part.electrical.negative, -1);
                    system.at(branch, branch) = -(p.resistance + imaginary * omega * p.inductance);
                    system.at(branch, velocity) = -p.forceFactor;
                    system.at(velocity, branch) = p.forceFactor;
                    system.at(velocity, velocity) =
                        -(p.mechanicalResistance + imaginary * omega * p.movingMass +
                          Complex{0, -1.0 / (omega * p.compliance)});
                    stampEffort(velocity, part.front, -p.pistonArea);
                    stampEffort(velocity, part.rear, p.pistonArea);
                } else if constexpr (std::is_same_v<Part, FlowSource>) {
                    const auto positiveDof = nodeDof(part.terminals.positive);
                    const auto negativeDof = nodeDof(part.terminals.negative);
                    if (positiveDof != noIndex)
                        system.rhs[positiveDof] -= part.amplitude;
                    if (negativeDof != noIndex)
                        system.rhs[negativeDof] += part.amplitude;
                } else {
                    stampFlow(part.terminals.positive, branch, 1);
                    stampFlow(part.terminals.negative, branch, -1);
                    stampEffort(branch, part.terminals.positive, 1);
                    stampEffort(branch, part.terminals.negative, -1);
                    if constexpr (std::is_same_v<Part, BemPort>) {
                        for (std::size_t j = 0; j < portCount; ++j) {
                            const auto &other = std::get<BemPort>(ports[j]->component);
                            system.at(branch, portBranches.at(ports[j]->id.value)) -=
                                static_cast<double>(part.polarity * other.polarity) *
                                load.impedance[portIndex * portCount + j];
                        }
                        if (!load.blockedPressure.empty())
                            system.rhs[branch] = static_cast<double>(part.polarity) *
                                                 load.blockedPressure[portIndex];
                        ++portIndex;
                    } else if constexpr (std::is_same_v<Part, EffortSource>)
                        system.rhs[branch] = part.amplitude;
                    else
                        system.at(branch, branch) -= impedance(element.component, omega, medium);
                }
            },
            element.component);
    }
    Solution result;
    result.frequency = frequency;
    const auto unknowns = system.solve(stop, result.relativeResidual);
    const auto effort = [&](NodeId node) -> Complex {
        const auto dof = nodeDof(node);
        return dof == noIndex ? Complex{} : unknowns[dof];
    };
    for (const auto &node : network.nodes)
        result.nodes.push_back({node.id, effort(node.id)});
    for (std::size_t index = 0; index < network.elements.size(); ++index) {
        const auto &element = network.elements[index];
        const auto branch = branches[index];
        ElementResult item;
        item.id = element.id;
        std::visit(
            [&](const auto &part) {
                using Part = std::decay_t<decltype(part)>;
                if constexpr (std::is_same_v<Part, Driver>) {
                    item.effort =
                        effort(part.electrical.positive) - effort(part.electrical.negative);
                    item.flow = unknowns[branch];
                    DriverResult driver;
                    driver.id = element.id;
                    driver.current = item.flow;
                    driver.velocity = unknowns[branch + 1];
                    driver.displacement = {driver.velocity.imag() / omega,
                                           -driver.velocity.real() / omega};
                    driver.electricalPower = activePower(item.effort, driver.current);
                    driver.copperPower = std::norm(driver.current) * part.parameters.resistance;
                    driver.mechanicalPower =
                        std::norm(driver.velocity) * part.parameters.mechanicalResistance;
                    driver.acousticPower =
                        activePower(effort(part.front) - effort(part.rear),
                                    part.parameters.pistonArea * driver.velocity);
                    if (!finite(driver.displacement) || !std::isfinite(driver.copperPower) ||
                        !std::isfinite(driver.mechanicalPower) ||
                        !std::isfinite(driver.acousticPower))
                        throw std::runtime_error("Die LEM-Treibergrößen sind nicht endlich.");
                    result.drivers.push_back(driver);
                } else {
                    item.effort = effort(part.terminals.positive) - effort(part.terminals.negative);
                    if constexpr (std::is_same_v<Part, FlowSource>)
                        item.flow = part.amplitude;
                    else
                        item.flow = unknowns[branch];
                    if constexpr (std::is_same_v<Part, BemPort>) {
                        const auto sign = static_cast<double>(part.polarity);
                        result.ports.push_back({element.id, sign * item.effort, sign * item.flow,
                                                activePower(item.effort, item.flow)});
                    }
                }
                item.activePower = activePower(item.effort, item.flow);
            },
            element.component);
        if (!finite(item.effort) || !finite(item.flow) || !std::isfinite(item.activePower))
            throw std::runtime_error("Die LEM-Ergebnisgrößen sind nicht endlich.");
        result.elements.push_back(item);
    }
    return result;
}

Network makeDriverNetwork(const DriverParameters &parameters,
                          const std::vector<std::vector<std::uint32_t>> &faceGroups,
                          Complex voltageRms, double rearVolume) {
    require(!faceGroups.empty(), "Für das Treibernetzwerk fehlt eine BEM-Flächengruppe.");
    require(faceGroups.size() < (std::numeric_limits<std::uint32_t>::max() - 4U) / 3U,
            "Zu viele Treiber für das LEM-Netzwerk.");
    Network result;
    result.nodes = {{{1}, Domain::electrical, true, "Masse", {0, 160}},
                    {{2}, Domain::acoustic, true, "Umgebung", {720, 160}},
                    {{3}, Domain::electrical, false, "Verstärker", {160, 0}}};
    result.elements.push_back({{1}, "Spannung", {0, 0}, EffortSource{{{3}, {1}}, voltageRms}});
    for (std::size_t i = 0; i < faceGroups.size(); ++i) {
        const auto count = static_cast<std::uint32_t>(i);
        const NodeId front{4 + 2 * count}, rear{5 + 2 * count};
        const double y = 220 * static_cast<double>(i);
        result.nodes.push_back(
            {front, Domain::acoustic, false, "Vorn " + std::to_string(i + 1), {520, y}});
        result.nodes.push_back(
            {rear, Domain::acoustic, false, "Hinten " + std::to_string(i + 1), {360, y + 100}});
        result.elements.push_back({{2 + 3 * count},
                                   "Treiber " + std::to_string(i + 1),
                                   {320, y},
                                   Driver{{{3}, {1}}, front, rear, parameters}});
        result.elements.push_back({{3 + 3 * count},
                                   "Rückkammer " + std::to_string(i + 1),
                                   {520, y + 100},
                                   Cavity{{rear, {2}}, rearVolume}});
        result.elements.push_back({{4 + 3 * count},
                                   "BEM " + std::to_string(i + 1),
                                   {680, y},
                                   BemPort{{front, {2}}, faceGroups[i], 1}});
    }
    validateNetwork(result);
    return result;
}

} // namespace hornsculpt::lem
