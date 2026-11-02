#pragma once
#include <cstdint>
#include <hornsculpt/lem.hpp>
#include <optional>
#include <string>

namespace hornsculpt::lem {
enum class DriverKind : std::uint8_t { tweeter, woofer, compressionDriver };
enum class DriverEntryMode : std::uint8_t { thieleSmall, electromechanical };
struct ThieleSmallParameters {
    double resistance{};         // Ohm.
    double inductance{};         // H.
    double pistonArea{};         // m^2.
    double resonanceFrequency{}; // Hz, free-air resonance.
    double mechanicalQ{};
    double electricalQ{};
    double equivalentVolume{}; // m^3.
    // Explicit external air inertia already included in Mms, kg. Zero is a
    // deliberate approximation; an empty value means the correction has not been specified.
    std::optional<double> externalAirMass;
    bool operator==(const ThieleSmallParameters &) const = default;
};
struct DriverDefinition {
    std::uint32_t id{};
    std::string name;
    DriverKind kind{DriverKind::tweeter};
    DriverEntryMode entryMode{DriverEntryMode::thieleSmall};
    ThieleSmallParameters thieleSmall;
    DriverParameters electromechanical;
    std::string provenance;
    std::string assumptions;
    bool operator==(const DriverDefinition &) const = default;
};
// Converts free-air T/S definitions using Vas=rho*c^2*Sd^2*Cms and the
// resonance/Q definitions. The solver receives dry Mmd, never loaded Mms.
// Incomplete library entries are editable; conversion rejects them explicitly.
DriverParameters driverParameters(const DriverDefinition &driver, double density = 1.204,
                                  double soundSpeed = 343.0);
// Effective reference approximation with integral rear loading retained; not a measured dry driver.
DriverDefinition scanSpeakD2608ReferenceDriver(std::uint32_t id);
} // namespace hornsculpt::lem
