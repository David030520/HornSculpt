#include <cmath>
#include <hornsculpt/driver_model.hpp>
#include <numbers>
#include <stdexcept>

namespace hornsculpt::lem {
namespace {
void positive(double value, const char *name) {
    if (!std::isfinite(value) || value <= 0)
        throw std::invalid_argument(std::string(name) + " muss endlich und positiv sein.");
}
void nonnegative(double value, const char *name) {
    if (!std::isfinite(value) || value < 0)
        throw std::invalid_argument(std::string(name) + " muss endlich und nicht negativ sein.");
}
void validate(const DriverParameters &p) {
    positive(p.resistance, "Re");
    nonnegative(p.inductance, "Le");
    positive(p.forceFactor, "Bl");
    positive(p.movingMass, "Mmd");
    positive(p.compliance, "Cms");
    positive(p.mechanicalResistance, "Rms");
    positive(p.pistonArea, "Sd");
}
} // namespace
DriverParameters driverParameters(const DriverDefinition &driver, double density,
                                  double soundSpeed) {
    if (driver.entryMode == DriverEntryMode::electromechanical) {
        validate(driver.electromechanical);
        return driver.electromechanical;
    }
    if (driver.entryMode != DriverEntryMode::thieleSmall)
        throw std::invalid_argument("Unbekannte Treiberparametrisierung.");
    positive(density, "Luftdichte");
    positive(soundSpeed, "Schallgeschwindigkeit");
    const auto &p = driver.thieleSmall;
    positive(p.resistance, "Re");
    nonnegative(p.inductance, "Le");
    positive(p.pistonArea, "Sd");
    positive(p.resonanceFrequency, "Fs");
    positive(p.mechanicalQ, "Qms");
    positive(p.electricalQ, "Qes");
    positive(p.equivalentVolume, "Vas");
    if (!p.externalAirMass)
        throw std::invalid_argument("Externe Luftmasse explizit angeben.");
    nonnegative(*p.externalAirMass, "Externe Luftmasse (explizit angeben)");
    const double omega = 2 * std::numbers::pi * p.resonanceFrequency;
    const double compliance =
        p.equivalentVolume / (density * soundSpeed * soundSpeed * p.pistonArea * p.pistonArea);
    const double loadedMass = 1 / (omega * omega * compliance);
    DriverParameters result{p.resistance,
                            p.inductance,
                            std::sqrt(omega * loadedMass * p.resistance / p.electricalQ),
                            loadedMass - *p.externalAirMass,
                            compliance,
                            omega * loadedMass / p.mechanicalQ,
                            p.pistonArea};
    validate(result);
    return result;
}
DriverDefinition scanSpeakD2608ReferenceDriver(std::uint32_t id) {
    DriverDefinition driver;
    driver.id = id;
    driver.name = "Scan-Speak D2608/913000 (Referenznäherung)";
    driver.entryMode = DriverEntryMode::electromechanical;
    const double radius = std::sqrt(0.0007 / std::numbers::pi);
    // Low-ka, front-only baffled circular piston: Ma=8*rho*a^3/3.
    const double frontAirMass = 8 * 1.21 * radius * radius * radius / 3;
    driver.electromechanical = {5.6, 0.00004, 2.6, 0.00018 - frontAirMass, 0.00029, 1.54, 0.0007};
    driver.provenance = "https://www.scan-speak.dk/datasheet/pdf/d2608-913000.pdf; gerundete "
                        "elektromechanische Herstellerwerte.";
    driver.assumptions =
        "Konzeptmodell; angenommene Luftlastkorrektur. Effektives Kleinsignalmodell mit "
        "integriertem Rückraum; keine zusätzliche Rückkammer. "
        "Mms=0,18 g minus angenommene frontseitige Kolben-Luftmasse 0,010732 g (rho=1,21 kg/m³, "
        "low-ka). "
        "Keine gemessene trockene Masse; Cms/Rms enthalten unbekannte interne Lastanteile. "
        "Gerundete Daten erfüllen T/S-Beziehungen nur näherungsweise; keine "
        "Herstellerfrequenzgang-Parität.";
    return driver;
}
} // namespace hornsculpt::lem
