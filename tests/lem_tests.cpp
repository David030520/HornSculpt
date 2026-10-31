#include <hornsculpt/lem.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <stop_token>

using namespace hornsculpt::lem;
namespace {
constexpr Complex imaginary{0, 1};
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
void near(Complex actual, Complex expected, const char *message) {
    require(std::isfinite(actual.real()) && std::isfinite(actual.imag()) &&
                std::abs(actual - expected) <= 2e-10 * std::max(1.0, std::abs(expected)),
            message);
}
template <class Exception = std::invalid_argument, class Function> void rejects(Function function) {
    bool rejected = false;
    try {
        function();
    } catch (const Exception &) {
        rejected = true;
    }
    require(rejected, "An invalid or underdetermined network was accepted.");
}
const ElementResult &element(const Solution &solution, std::uint32_t id) {
    const auto found = std::find_if(solution.elements.begin(), solution.elements.end(),
                                    [=](const auto &item) { return item.id.value == id; });
    require(found != solution.elements.end(), "An element result is missing.");
    return *found;
}
Complex node(const Solution &solution, std::uint32_t id) {
    const auto found = std::find_if(solution.nodes.begin(), solution.nodes.end(),
                                    [=](const auto &item) { return item.id.value == id; });
    require(found != solution.nodes.end(), "A node result is missing.");
    return found->effort;
}
DriverParameters driverParameters() {
    const double mass = 0.004;
    const double resonance = 2 * std::numbers::pi * 62;
    return {6.4, 0.0007, 4.98, mass, 1 / (resonance * resonance * mass), 0.262, 0.00544};
}

void seriesCircuits() {
    // Electrical and mechanical effort/flow analogies have the same series law.
    for (const auto domain : {Domain::electrical, Domain::mechanical, Domain::acoustic}) {
        Network network;
        network.nodes = {{{1}, domain, true, {}, {}},
                         {{2}, domain, false, {}, {}},
                         {{3}, domain, false, {}, {}},
                         {{4}, domain, false, {}, {}}};
        network.elements = {{{1}, "source", {}, EffortSource{{{2}, {1}}, {2.5, -0.75}}},
                            {{2}, "R", {}, Resistance{{{2}, {3}}, 12}},
                            {{3}, "L", {}, Inertance{{{3}, {4}}, 0.017}},
                            {{4}, "C", {}, Compliance{{{4}, {1}}, 4e-5}}};
        for (const double frequency : {20.0, 300.0, 20000.0}) {
            const double omega = 2 * std::numbers::pi * frequency;
            const Complex impedance =
                12.0 + imaginary * omega * 0.017 + 1.0 / (imaginary * omega * 4e-5);
            const Complex current = Complex{2.5, -0.75} / impedance;
            const auto result = solve(network, frequency);
            near(element(result, 1).flow, -current, "The source current has the wrong direction.");
            near(element(result, 2).flow, current,
                 "The RLC current differs from the series reference.");
            near(element(result, 3).flow, current, "A series inertance violates continuity.");
            near(element(result, 4).flow, current, "A series compliance violates continuity.");
            near(element(result, 2).activePower, 12 * std::norm(current),
                 "The resistor RMS power is wrong.");
            near(element(result, 3).activePower, 0, "An ideal inertance dissipates power.");
            near(element(result, 4).activePower, 0, "An ideal compliance dissipates power.");
            near(element(result, 1).activePower + element(result, 2).activePower, 0,
                 "The source and load power do not balance.");
            require(result.relativeResidual < 1e-12, "The simple circuit residual is excessive.");
        }
    }
}

void cavityAndFlowSource() {
    Network network;
    network.nodes = {{{1}, Domain::acoustic, true, {}, {}}, {{2}, Domain::acoustic, false, {}, {}}};
    network.elements = {{{1}, "flow", {}, FlowSource{{{1}, {2}}, {0.004, -0.001}}},
                        {{2}, "leak", {}, Resistance{{{2}, {1}}, 15000}},
                        {{3}, "air", {}, Cavity{{{2}, {1}}, 0.0012}}};
    for (const Medium medium : {Medium{1.21, 343}, Medium{1.05, 365}}) {
        const double omega = 2 * std::numbers::pi * 300;
        const double compliance = 0.0012 / (medium.density * medium.soundSpeed * medium.soundSpeed);
        const Complex pressure =
            Complex{0.004, -0.001} / (1.0 / 15000 + imaginary * omega * compliance);
        const auto result = solve(network, 300, {}, medium);
        near(node(result, 2), pressure, "The cavity was not recalculated from the current medium.");
        near(element(result, 3).flow, imaginary * omega * compliance * pressure,
             "The cavity flow has the wrong sign or units.");
        near(element(result, 1).activePower + element(result, 2).activePower, 0,
             "The flow-source power balance is wrong.");
    }
}
Network slitNetwork(AcousticSlit slit) {
    slit.terminals = {{2}, {1}};
    Network network;
    network.nodes = {{{1}, Domain::acoustic, true, {}, {}}, {{2}, Domain::acoustic, false, {}, {}}};
    network.elements = {{{1}, "Flow", {}, FlowSource{{{1}, {2}}, {1, 0}}},
                        {{2}, "Parallel-plate slit", {}, slit}};
    return network;
}
Complex measuredSlitImpedance(const AcousticSlit &slit, double frequency, Medium medium = {}) {
    const auto solution = solve(slitNetwork(slit), frequency, {}, medium);
    const auto &branch = element(solution, 2);
    const auto z = branch.effort / branch.flow;
    require(z.real() > 0 && z.imag() > 0 && solution.relativeResidual < 1e-12,
            "Unsteady-Stokes slit is nonpassive or violates its MNA equations.");
    require(std::abs(element(solution, 1).activePower + branch.activePower) <
                1e-12 * branch.activePower,
            "Slit RMS input power is not balanced by its driving source.");
    return z;
}
void acousticSlit() {
    const AcousticSlit slit{{{2}, {1}}, .006, 8e-5, 1.5e-4, 1.84e-5};
    const Medium medium{1.21, 343};
    const double resistance =
                     12 * slit.dynamicViscosity * slit.length / (slit.area * slit.gap * slit.gap),
                 mass = medium.density * slit.length / slit.area;
    for (double frequency : {1e-12, 1e-6, .001, 1.0, 100.0, 1e4, 1e8, 1e12}) {
        const auto z = measuredSlitImpedance(slit, frequency, medium);
        require(std::arg(z) > 0 && std::arg(z) < std::numbers::pi / 2,
                "Passive slit pressure must lead its volume flow by between zero and pi/2.");
        if (frequency <= 1e-6) {
            require(std::abs(z.real() / resistance - 1) < 1e-12,
                    "Low-frequency slit resistance differs from planar Poiseuille flow.");
            require(std::abs(z.imag() / (2 * std::numbers::pi * frequency) / (1.2 * mass) - 1) <
                        2e-10,
                    "Low-frequency slit inertance is not six-fifths of plug-flow mass.");
        }
        if (frequency == 1e12)
            require(std::abs(z.imag() / (2 * std::numbers::pi * frequency) / mass - 1) < 3e-5 &&
                        z.real() / z.imag() < 3e-5,
                    "High-frequency slit did not approach its inviscid inertance.");
    }
    // Exercise both sides of the cancellation-free series and the large-chi
    // tanh saturation against the analytic expression where it is well scaled.
    for (double chiMagnitude : {.099999, .100001, 1.0, 10.0, 28.2, 28.4, 40.0}) {
        const double omega = chiMagnitude * chiMagnitude * 4 * slit.dynamicViscosity /
                             (medium.density * slit.gap * slit.gap);
        const auto chi =
            (slit.gap / 2) * std::sqrt(imaginary * omega * medium.density / slit.dynamicViscosity);
        const auto expected = imaginary * omega * mass / (1.0 - std::tanh(chi) / chi);
        const auto actual = measuredSlitImpedance(slit, omega / (2 * std::numbers::pi), medium);
        require(std::abs(actual / expected - 1.0) < 2e-11,
                "Slit transition or viscosity phase differs from the analytical Stokes solution.");
    }
    auto viscous = slit;
    viscous.dynamicViscosity *= 2;
    const auto low = measuredSlitImpedance(slit, 1e-6, medium),
               twice = measuredSlitImpedance(viscous, 1e-6, medium),
               dense = measuredSlitImpedance(slit, 1e-6, {2 * medium.density, 300});
    require(std::abs(twice.real() / low.real() - 2) < 1e-12 &&
                std::abs(twice.imag() / low.imag() - 1) < 1e-10 &&
                std::abs(dense.real() / low.real() - 1) < 1e-12 &&
                std::abs(dense.imag() / low.imag() - 2) < 1e-10,
            "Viscosity/density do not independently set the correct slit resistance/inertia.");

    Network loaded;
    loaded.nodes = {{{1}, Domain::acoustic, true, {}, {}},
                    {{2}, Domain::acoustic, false, {}, {}},
                    {{3}, Domain::acoustic, false, {}, {}}};
    auto series = slit;
    series.terminals = {{2}, {3}};
    loaded.elements = {
        {{1}, "Interface", {}, BemPort{{{2}, {1}}, {42}, 1, BemPortKind::interfacePort}},
        {{2}, "Slit", {}, series},
        {{3}, "Volume", {}, Cavity{{{3}, {1}}, 5e-6}}};
    const AcousticLoad front{{{2e4, 5e3}}, {{100, 30}}};
    for (double frequency : {200.0, 2500.0, 8000.0}) {
        const auto zs = measuredSlitImpedance(slit, frequency, medium);
        const auto zc = medium.density * medium.soundSpeed * medium.soundSpeed /
                        (imaginary * (2 * std::numbers::pi * frequency) * 5e-6);
        const auto expected = -front.blockedPressure[0] / (front.impedance[0] + zs + zc);
        const auto solution = solve(loaded, frequency, front, medium);
        require(std::abs(solution.ports[0].volumeVelocity / expected - 1.0) < 1e-10,
                "Series slit/cavity does not feed pressure back through the BEM load.");
        near(solution.ports[0].volumeVelocity, -element(solution, 2).flow,
             "BEM and slit volume flows are discontinuous.");
        near(element(solution, 2).flow, element(solution, 3).flow,
             "Slit and cavity volume flows are discontinuous.");
        near(solution.ports[0].pressure, -(zs + zc) * expected,
             "Interface pressure does not include slit pressure drop before cavity pressure.");
        require(std::abs(solution.ports[0].activePower + element(solution, 2).activePower) <
                    1e-10 * element(solution, 2).activePower,
                "Passive slit loss and power extracted from BEM do not balance.");
    }
    const auto valid = slitNetwork(slit);
    for (int parameter = 0; parameter < 4; ++parameter) {
        for (double value : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()}) {
            auto bad = valid;
            auto &part = std::get<AcousticSlit>(bad.elements[1].component);
            double *field = parameter == 0   ? &part.length
                            : parameter == 1 ? &part.area
                            : parameter == 2 ? &part.gap
                                             : &part.dynamicViscosity;
            *field = value;
            rejects([&] { validateNetwork(bad); });
        }
        auto edited = valid;
        auto &part = std::get<AcousticSlit>(edited.elements[1].component);
        if (parameter == 0)
            part.length *= 2;
        if (parameter == 1)
            part.area *= 2;
        if (parameter == 2)
            part.gap *= 2;
        if (parameter == 3)
            part.dynamicViscosity *= 2;
        require(!samePhysics(valid, edited),
                "Changing a slit parameter does not invalidate the solve.");
    }
    auto wrongDomain = valid;
    for (auto &n : wrongDomain.nodes)
        n.domain = Domain::electrical;
    rejects([&] { validateNetwork(wrongDomain); });
}
struct SlitLineReference {
    Complex beta, k, characteristic;
};
SlitLineReference lineReference(const ClosedSlitCavity &line, double frequency,
                                Medium medium = {}) {
    const double omega = 2 * std::numbers::pi * frequency;
    const auto chi =
        (line.gap / 2) * std::sqrt(imaginary * omega * medium.density / line.dynamicViscosity);
    const auto thermal = chi * std::sqrt(line.prandtlNumber);
    const auto density = medium.density / (1.0 - std::tanh(chi) / chi);
    const auto beta = (1.0 + (line.heatCapacityRatio - 1) * std::tanh(thermal) / thermal) /
                      (medium.density * medium.soundSpeed * medium.soundSpeed);
    const auto k = omega * std::sqrt(density * beta);
    return {beta, k, std::sqrt(density / beta) / line.area};
}
Complex measuredClosedImpedance(ClosedSlitCavity line, double frequency, Medium medium = {}) {
    line.terminals = {{2}, {1}};
    Network network;
    network.nodes = {{{1}, Domain::acoustic, true, {}, {}}, {{2}, Domain::acoustic, false, {}, {}}};
    network.elements = {{{1}, "Flow", {}, FlowSource{{{1}, {2}}, {1, 0}}},
                        {{2}, "Closed line", {}, line}};
    const auto solution = solve(network, frequency, {}, medium);
    const auto &branch = element(solution, 2);
    const auto z = branch.effort / branch.flow;
    require(std::isfinite(z.real()) && std::isfinite(z.imag()) && z.real() > 0 &&
                solution.relativeResidual < 1e-12,
            "Closed thermoviscous slit is nonpassive or has an invalid MNA residual.");
    require(std::abs(element(solution, 1).activePower + branch.activePower) <
                1e-10 * branch.activePower,
            "Closed slit and its source violate RMS power conservation.");
    return z;
}
void closedSlitCavity() {
    const Medium medium{1.21, 343};
    const ClosedSlitCavity line{{{2}, {1}}, .03, 5e-6 / .03, .00175, 1.84e-5, 1.4, .71};
    const double lowFrequency = 1e-6, omegaLow = 2 * std::numbers::pi * lowFrequency,
                 volume = line.area * line.length,
                 isothermal = line.heatCapacityRatio * volume /
                              (medium.density * medium.soundSpeed * medium.soundSpeed);
    const auto low = measuredClosedImpedance(line, lowFrequency, medium);
    require(std::abs((-1 / (omegaLow * low.imag())) / isothermal - 1) < 1e-9,
            "Static closed slit does not approach its isothermal compliance.");
    auto adiabatic = line;
    adiabatic.heatCapacityRatio = 1;
    const auto adiabaticLow = measuredClosedImpedance(adiabatic, lowFrequency, medium);
    require(std::abs(adiabaticLow.imag() / low.imag() - line.heatCapacityRatio) < 1e-9,
            "Thermal and adiabatic static compliances have the wrong ratio.");
    for (unsigned i = 0; i < 96; ++i) {
        const double frequency = std::pow(10.0, -3.0 + 10.0 * i / 95);
        (void)measuredClosedImpedance(line, frequency, medium);
    }
    for (double frequency : {10.0, 100.0, 2858.0, 5716.0, 16000.0}) {
        const auto reference = lineReference(line, frequency, medium);
        require(reference.k.real() > 0 && reference.k.imag() < 0,
                "Reference line picked a growing wave for positive time convention.");
        const auto expected =
            -imaginary * reference.characteristic / std::tan(reference.k * line.length);
        const auto actual = measuredClosedImpedance(line, frequency, medium);
        require(std::abs(actual / expected - 1.0) < 1e-10,
                "Thermoviscous closed-line impedance differs from its analytical transfer matrix.");
    }
    auto lossless = line;
    lossless.dynamicViscosity = 1e-14;
    lossless.heatCapacityRatio = 1;
    const double quarter = medium.soundSpeed / (4 * line.length),
                 z0 = medium.density * medium.soundSpeed / line.area;
    const auto quarterZ = measuredClosedImpedance(lossless, quarter, medium),
               halfZ = measuredClosedImpedance(lossless, 2 * quarter, medium),
               below = measuredClosedImpedance(lossless, .95 * quarter, medium),
               above = measuredClosedImpedance(lossless, 1.05 * quarter, medium);
    require(std::abs(quarterZ) / z0 < 1e-4 && std::abs(halfZ) / z0 > 1e4 && below.imag() < 0 &&
                above.imag() > 0,
            "Closed line does not approach its quarter-wave zero and half-wave pole.");
    auto small = line;
    const auto reference = lineReference(line, 1000, medium);
    for (double argument : {.099999, .100001}) {
        small.length = argument / std::abs(reference.k);
        const auto expected =
            -imaginary * reference.characteristic / std::tan(reference.k * small.length);
        require(std::abs(measuredClosedImpedance(small, 1000, medium) / expected - 1.0) < 2e-11,
                "Small-kL closed-line expansion is discontinuous or inaccurate.");
    }
    auto damped = line;
    damped.gap = 1e-6;
    const auto damping = lineReference(damped, 1e4, medium);
    for (double attenuation : {19.9999, 20.0001, 50.0}) {
        damped.length = attenuation / (-damping.k.imag());
        const auto expected = attenuation > 20 ? damping.characteristic
                                               : -imaginary * damping.characteristic /
                                                     std::tan(damping.k * damped.length);
        require(
            std::abs(measuredClosedImpedance(damped, 1e4, medium) / expected - 1.0) < 1e-8,
            "Highly attenuated closed line does not stably approach its characteristic impedance.");
    }
    Network valid;
    valid.nodes = {{{1}, Domain::acoustic, true, {}, {}}, {{2}, Domain::acoustic, false, {}, {}}};
    valid.elements = {{{1}, "Line", {}, line}};
    for (int parameter = 0; parameter < 6; ++parameter) {
        for (double value : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()}) {
            auto bad = valid;
            auto &part = std::get<ClosedSlitCavity>(bad.elements[0].component);
            double *field = parameter == 0   ? &part.length
                            : parameter == 1 ? &part.area
                            : parameter == 2 ? &part.gap
                            : parameter == 3 ? &part.dynamicViscosity
                            : parameter == 4 ? &part.heatCapacityRatio
                                             : &part.prandtlNumber;
            *field = value;
            rejects([&] { validateNetwork(bad); });
        }
    }
    auto bad = valid;
    std::get<ClosedSlitCavity>(bad.elements[0].component).heatCapacityRatio = .9;
    rejects([&] { validateNetwork(bad); });
    auto edited = valid;
    std::get<ClosedSlitCavity>(edited.elements[0].component).length *= 2;
    require(!samePhysics(valid, edited),
            "Closed-line geometry edit does not invalidate its solution.");
}

struct DriverReference {
    Complex current;
    Complex velocity;
    Complex pressure;
};
DriverReference analyticDriver(const DriverParameters &parameters, Complex voltage,
                               Complex radiationImpedance, Complex blockedPressure,
                               double frequency, double volume, Medium medium) {
    const double omega = 2 * std::numbers::pi * frequency;
    const auto &p = parameters;
    const Complex ze = p.resistance + imaginary * omega * p.inductance;
    const Complex rearImpedance =
        medium.density * medium.soundSpeed * medium.soundSpeed / (imaginary * omega * volume);
    const Complex zm = p.mechanicalResistance + imaginary * omega * p.movingMass +
                       1.0 / (imaginary * omega * p.compliance) +
                       p.pistonArea * p.pistonArea * (radiationImpedance + rearImpedance);
    const Complex velocity = (p.forceFactor * voltage - ze * p.pistonArea * blockedPressure) /
                             (ze * zm + p.forceFactor * p.forceFactor);
    const Complex current = (voltage - p.forceFactor * velocity) / ze;
    return {current, velocity, radiationImpedance * p.pistonArea * velocity + blockedPressure};
}

void driverPressureFeedback() {
    const auto parameters = driverParameters();
    const Medium medium{1.2, 343};
    const Complex voltage{27, 4};
    constexpr double rearVolume = 0.000344;
    auto network = makeDriverNetwork(parameters, {{17, 19}}, voltage, rearVolume);
    for (const double frequency : {60.0, 300.0, 1500.0}) {
        for (const Complex radiation : {Complex{12000, 4000}, Complex{80000, -12000}}) {
            const Complex blocked{80, -30};
            const auto expected = analyticDriver(parameters, voltage, radiation, blocked, frequency,
                                                 rearVolume, medium);
            const auto result = solve(network, frequency, {{radiation}, {blocked}}, medium);
            require(result.drivers.size() == 1 && result.ports.size() == 1,
                    "Driver or port output is missing.");
            const auto &driver = result.drivers.front();
            near(driver.current, expected.current,
                 "The back EMF does not affect electrical current correctly.");
            near(driver.velocity, expected.velocity,
                 "The driver does not respond to acoustic loading.");
            near(driver.displacement,
                 expected.velocity / (imaginary * 2.0 * std::numbers::pi * frequency),
                 "Displacement does not follow the positive-time RMS convention.");
            near(result.ports.front().pressure, expected.pressure,
                 "The pressure feedback is inconsistent.");
            near(result.ports.front().volumeVelocity, parameters.pistonArea * expected.velocity,
                 "Piston velocity and BEM volume velocity are inconsistent.");
            near(driver.electricalPower,
                 driver.copperPower + driver.mechanicalPower + driver.acousticPower,
                 "The electromechanical power balance is not closed.");
            near(driver.acousticPower, result.ports.front().activePower,
                 "The sealed rear cavity absorbed real power.");
        }
    }
    const auto initial = solve(network, 300, {{{12000, 4000}}, {}}, medium);
    const auto loaded = solve(network, 300, {{{80000, -12000}}, {}}, medium);
    require(std::abs(initial.drivers[0].current - loaded.drivers[0].current) > 0.05,
            "Changing the acoustic load was applied only after the electrical solve.");
}

void receivingDriverAndPolarity() {
    const auto parameters = driverParameters();
    constexpr double volume = 0.000344;
    const Medium medium{1.2, 343};
    auto network = makeDriverNetwork(parameters, {{42}}, {}, volume);
    const AcousticLoad load{{{15000, 3000}}, {{100, 40}}};
    auto received = solve(network, 300, load, medium);
    const auto expected = analyticDriver(parameters, {}, load.impedance[0], load.blockedPressure[0],
                                         300, volume, medium);
    near(received.drivers[0].velocity, expected.velocity,
         "Incident pressure failed to move the diaphragm.");
    require(received.drivers[0].copperPower > 0 && received.drivers[0].acousticPower < 0,
            "A shorted receiving driver has the wrong acoustic power sign.");
    std::get<BemPort>(network.elements.back().component).polarity = -1;
    const auto reversed = solve(network, 300, load, medium);
    near(reversed.drivers[0].velocity, -received.drivers[0].velocity,
         "Reversing the port polarity did not reverse the pressure force.");
    near(reversed.ports[0].volumeVelocity, received.ports[0].volumeVelocity,
         "A port sign change broke the physical flow convention.");
    near(reversed.ports[0].pressure, received.ports[0].pressure,
         "A port sign change broke the physical pressure convention.");
}

void interfaceSubdomains() {
    const Medium medium{1.18, 347};
    // Two distinct passive rear subdomains, coupled by a full front-domain Z.
    // Physical Q>0 enters the front BEM. Rear admittances therefore obey
    // Q_BEM=-Y_rear*p and p=p0+Z_front*Q_BEM, including the blocked source.
    Network network;
    network.nodes = {{{1}, Domain::acoustic, true, {}, {}},
                     {{2}, Domain::acoustic, false, {}, {}},
                     {{3}, Domain::acoustic, false, {}, {}}};
    constexpr double volumes[]{5e-6, 8e-6}, resistances[]{2e7, 3e7};
    for (std::uint32_t i = 0; i < 2; ++i) {
        const Terminals terminals{{i + 2}, {1}};
        network.elements.push_back({{3 * i + 1}, "Rear cavity", {}, Cavity{terminals, volumes[i]}});
        network.elements.push_back(
            {{3 * i + 2}, "Loss", {}, Resistance{terminals, resistances[i]}});
        network.elements.push_back(
            {{3 * i + 3},
             "Interface",
             {},
             BemPort{terminals, {20 + i}, i ? -1 : 1, BemPortKind::interfacePort}});
    }
    for (const double frequency : {150.0, 1700.0, 9000.0}) {
        const double omega = 2 * std::numbers::pi * frequency;
        const AcousticLoad load{
            {{1.2e6, omega * 24}, {2e5, omega * 5}, {2e5, omega * 5}, {1.7e6, omega * 31}},
            {{140, -30}, {-20, 85}}};
        const Complex y0 =
            1 / resistances[0] + imaginary * omega * volumes[0] /
                                     (medium.density * medium.soundSpeed * medium.soundSpeed);
        const Complex y1 =
            1 / resistances[1] + imaginary * omega * volumes[1] /
                                     (medium.density * medium.soundSpeed * medium.soundSpeed);
        const Complex a = 1.0 + y0 * load.impedance[0], b = y0 * load.impedance[1],
                      c = y1 * load.impedance[2], d = 1.0 + y1 * load.impedance[3];
        const Complex rhs0 = -y0 * load.blockedPressure[0], rhs1 = -y1 * load.blockedPressure[1];
        const Complex q[]{(d * rhs0 - b * rhs1) / (a * d - b * c),
                          (a * rhs1 - c * rhs0) / (a * d - b * c)};
        const auto solution = solve(network, frequency, load, medium);
        double passivePower{};
        for (std::size_t i = 0; i < 2; ++i) {
            const auto p = load.blockedPressure[i] + load.impedance[2 * i] * q[0] +
                           load.impedance[2 * i + 1] * q[1];
            near(solution.ports[i].volumeVelocity, q[i],
                 "Subdomain volume flow violates feedback.");
            near(solution.ports[i].pressure, p, "Blocked interface pressure or phase is lost.");
            near(p, -q[i] / (i ? y1 : y0), "Front/rear interface flow continuity is broken.");
            passivePower += solution.ports[i].activePower +
                            element(solution, static_cast<std::uint32_t>(3 * i + 2)).activePower;
        }
        near(passivePower, 0, "Interface and rear-subdomain active powers do not balance.");
        const auto radiatedPower =
            std::real((load.impedance[0] * q[0] + load.impedance[1] * q[1]) * std::conj(q[0]) +
                      (load.impedance[2] * q[0] + load.impedance[3] * q[1]) * std::conj(q[1]));
        require(radiatedPower >= 0 && solution.relativeResidual < 1e-12,
                "Passive acoustic multiport or original MNA residual is invalid.");
    }

    Network sealed;
    sealed.nodes = {{{1}, Domain::acoustic, true, {}, {}}, {{2}, Domain::acoustic, false, {}, {}}};
    sealed.elements = {
        {{1}, "Tiny rear cavity", {}, Cavity{{{2}, {1}}, 1e-12}},
        {{2}, "Interface", {}, BemPort{{{2}, {1}}, {42}, 1, BemPortKind::interfacePort}}};
    const AcousticLoad blocked{{{1e6, 2e5}}, {{170, -80}}};
    const auto tiny = solve(sealed, 700, blocked, medium);
    require(std::abs(tiny.ports[0].pressure / blocked.blockedPressure[0] - 1.0) < 1e-7,
            "Vanishing rear volume did not approach the blocked interface limit.");
    const auto compliance = 1e-12 / (medium.density * medium.soundSpeed * medium.soundSpeed);
    const auto expectedFlow =
        -imaginary * (2 * std::numbers::pi * 700) * compliance * tiny.ports[0].pressure;
    require(std::abs(tiny.ports[0].volumeVelocity / expectedFlow - 1.0) < 1e-10,
            "Sealed rear volume has the wrong positive-time compliance sign.");
}

void coupledDrivers() {
    const auto p = driverParameters();
    constexpr double volume = 0.000344;
    constexpr double frequency = 300;
    const Medium medium{1.2, 343};
    const Complex voltage{27, 2};
    auto network = makeDriverNetwork(p, {{101}, {202}}, voltage, volume);
    const AcousticLoad load{{{12000, 3000}, {5000, -1000}, {5000, -1000}, {18000, 4000}}, {}};
    const auto result = solve(network, frequency, load, medium);
    const double omega = 2 * std::numbers::pi * frequency;
    const Complex ze = p.resistance + imaginary * omega * p.inductance;
    const Complex source = p.forceFactor * voltage / ze;
    const Complex internal = p.mechanicalResistance + imaginary * omega * p.movingMass +
                             1.0 / (imaginary * omega * p.compliance) +
                             p.forceFactor * p.forceFactor / ze +
                             p.pistonArea * p.pistonArea * medium.density * medium.soundSpeed *
                                 medium.soundSpeed / (imaginary * omega * volume);
    const double areaSquared = p.pistonArea * p.pistonArea;
    const Complex a = internal + areaSquared * load.impedance[0];
    const Complex b = areaSquared * load.impedance[1];
    const Complex c = areaSquared * load.impedance[2];
    const Complex d = internal + areaSquared * load.impedance[3];
    const Complex determinant = a * d - b * c;
    const Complex v0 = source * (d - b) / determinant;
    const Complex v1 = source * (a - c) / determinant;
    near(result.drivers[0].velocity, v0, "The first diaphragm ignores mutual acoustic impedance.");
    near(result.drivers[1].velocity, v1, "The second diaphragm ignores mutual acoustic impedance.");
    near(result.ports[0].pressure, p.pistonArea * (load.impedance[0] * v0 + load.impedance[1] * v1),
         "The first multiport pressure is incorrect.");
    near(result.ports[1].pressure, p.pistonArea * (load.impedance[2] * v0 + load.impedance[3] * v1),
         "The second multiport pressure is incorrect.");
    near(element(result, 1).flow, -result.drivers[0].current - result.drivers[1].current,
         "Parallel drivers do not share the source current correctly.");
    double power = element(result, 1).activePower;
    for (const auto &driver : result.drivers)
        power += driver.copperPower + driver.mechanicalPower;
    for (const auto &port : result.ports)
        power += port.activePower;
    near(power, 0, "The coupled multi-driver network violates energy conservation.");
    auto uncoupled = load;
    uncoupled.impedance[1] = uncoupled.impedance[2] = {};
    const auto withoutMutual = solve(network, frequency, uncoupled, medium);
    require(std::abs(withoutMutual.drivers[0].velocity - result.drivers[0].velocity) > 0.01,
            "The mutual-impedance test did not exercise a material interaction.");
}

void physicsIdentity() {
    const auto original = makeDriverNetwork(driverParameters(), {{11}}, {10, 2}, .0004);
    auto edited = original;
    for (auto &item : edited.nodes) {
        item.name += " renamed";
        item.position = {item.position.x + 1, item.position.y - 2};
    }
    for (auto &item : edited.elements) {
        item.name += " renamed";
        item.position = {item.position.x - 3, item.position.y + 4};
    }
    require(samePhysics(original, edited) && samePhysics(edited, original),
            "Editor names or layout changed physical network identity.");
    const auto changed = [&](auto edit) {
        auto candidate = original;
        edit(candidate);
        require(!samePhysics(original, candidate) && !samePhysics(candidate, original),
                "A physical or ordered network identity change was ignored.");
    };
    changed([](Network &n) { n.nodes.front().reference = !n.nodes.front().reference; });
    changed([](Network &n) { n.nodes.front().domain = Domain::mechanical; });
    changed([](Network &n) { n.nodes.front().id = {123}; });
    changed([](Network &n) { n.elements.front().id = {456}; });
    changed(
        [](Network &n) { std::get<EffortSource>(n.elements.front().component).amplitude *= 2; });
    changed(
        [](Network &n) { std::get<Driver>(n.elements[1].component).parameters.resistance += 1; });
    changed([](Network &n) { std::get<Driver>(n.elements[1].component).front = {234}; });
    changed([](Network &n) { std::get<BemPort>(n.elements.back().component).faceIds = {12}; });
    changed([](Network &n) { std::get<BemPort>(n.elements.back().component).polarity = -1; });
    changed([](Network &n) {
        std::get<BemPort>(n.elements.back().component).kind = BemPortKind::interfacePort;
    });
    changed([](Network &n) {
        std::get<BemPort>(n.elements.back().component).needsAssignmentReview = true;
    });
    changed([](Network &n) {
        std::get<BemPort>(n.elements.back().component).velocityWeights = {{11, .5}};
    });
    changed([](Network &n) { n.elements.front().component = FlowSource{{{3}, {1}}, {10, 2}}; });
    changed([](Network &n) { std::swap(n.nodes[0], n.nodes[1]); });
    changed([](Network &n) { std::swap(n.elements[0], n.elements[1]); });
    changed([](Network &n) { n.nodes.pop_back(); });
    changed([](Network &n) { n.elements.pop_back(); });
    require(samePhysics({}, {}), "Empty editor networks have unstable simulation identity.");
}

void invalidNetworks() {
    const auto valid = makeDriverNetwork(driverParameters(), {{11}}, {10, 0}, 0.0004);
    const AcousticLoad load{{{12000, 0}}, {}};
    for (const double frequency : {0.0, -300.0, std::numeric_limits<double>::infinity()})
        rejects([&] { (void)solve(valid, frequency, load); });
    rejects([&] { (void)solve(valid, 300); });
    rejects([&] { (void)solve(valid, 300, {{Complex{1, 0}}, {{1, 0}, {2, 0}}}); });
    rejects(
        [&] { (void)solve(valid, 300, {{{std::numeric_limits<double>::quiet_NaN(), 0}}, {}}); });
    rejects([&] { (void)solve(valid, 300, load, {0, 343}); });
    auto bad = valid;
    bad.nodes[0].domain = Domain::acoustic;
    rejects([&] { validateNetwork(bad); });
    bad = valid;
    bad.nodes[0].reference = false;
    rejects([&] { validateNetwork(bad); });
    bad = valid;
    bad.nodes[1].id = bad.nodes[0].id;
    rejects([&] { validateNetwork(bad); });
    bad = valid;
    std::get<Driver>(bad.elements[1].component).front = {999};
    rejects([&] { validateNetwork(bad); });
    bad = valid;
    std::get<Driver>(bad.elements[1].component).parameters.compliance = 0;
    rejects([&] { validateNetwork(bad); });
    bad = valid;
    std::get<BemPort>(bad.elements.back().component).polarity = 0;
    rejects([&] { validateNetwork(bad); });
    bad = valid;
    std::get<BemPort>(bad.elements.back().component).kind = static_cast<BemPortKind>(255);
    rejects([&] { validateNetwork(bad); });
    bad = valid;
    std::get<BemPort>(bad.elements.back().component).needsAssignmentReview = true;
    rejects([&] { validateNetwork(bad); });
    bad = valid;
    std::get<BemPort>(bad.elements.back().component).faceIds = {11, 11};
    rejects([&] { validateNetwork(bad); });
    for (double weight : {0.0, -.5, std::numeric_limits<double>::infinity(),
                          std::numeric_limits<double>::quiet_NaN()}) {
        bad = valid;
        std::get<BemPort>(bad.elements.back().component).velocityWeights = {{11, weight}};
        rejects([&] { validateNetwork(bad); });
    }
    bad = valid;
    std::get<BemPort>(bad.elements.back().component).velocityWeights = {{12, .5}};
    rejects([&] { validateNetwork(bad); });
    bad = valid;
    auto &weightedInterface = std::get<BemPort>(bad.elements.back().component);
    weightedInterface.kind = BemPortKind::interfacePort;
    weightedInterface.velocityWeights = {{11, .5}};
    rejects([&] { validateNetwork(bad); });
    bad = valid;
    std::get<BemPort>(bad.elements.back().component).velocityWeights = {{11, .5}};
    validateNetwork(bad);
    rejects([&] { (void)makeDriverNetwork(driverParameters(), {{11}, {11}}, {10, 0}, 0.0004); });
    bad = valid;
    auto duplicateSource = bad.elements.front();
    duplicateSource.id = {999};
    bad.elements.push_back(duplicateSource);
    rejects<std::runtime_error>([&] { (void)solve(bad, 300, load); });
    std::stop_source cancellation;
    cancellation.request_stop();
    rejects<std::runtime_error>(
        [&] { (void)solve(valid, 300, load, {}, cancellation.get_token()); });
}
} // namespace

int main() {
    try {
        seriesCircuits();
        cavityAndFlowSource();
        acousticSlit();
        closedSlitCavity();
        driverPressureFeedback();
        receivingDriverAndPolarity();
        interfaceSubdomains();
        coupledDrivers();
        physicsIdentity();
        invalidNetworks();
        std::cout << "LEM circuit, driver feedback, mutual impedance, power and validation tests "
                     "passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
