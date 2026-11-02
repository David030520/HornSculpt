#include <algorithm>
#include <cmath>
#include <hornsculpt/driver_model.hpp>
#include <iostream>
#include <numbers>
#include <stdexcept>
using namespace hornsculpt::lem;
namespace {
void require(bool v, const char *m) {
    if (!v)
        throw std::runtime_error(m);
}
void near(double a, double b) {
    require(std::abs(a - b) <= 1e-11 * std::max(std::abs(b), 1e-10),
            "Driver conversion differs from physical reference.");
}
template <class F> void rejects(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "Invalid driver accepted.");
}
} // namespace
int main() {
    try {
        constexpr double rho = 1.204, c = 343, sd = .013, mms = .015, cms = .0006, re = 5.8,
                         bl = 7.2, rms = 1.1, air = .001;
        const double omega = 1 / std::sqrt(mms * cms);
        DriverDefinition d;
        d.thieleSmall = {re,
                         .0003,
                         sd,
                         omega / (2 * std::numbers::pi),
                         omega * mms / rms,
                         omega * mms * re / (bl * bl),
                         rho * c * c * sd * sd * cms,
                         air};
        const auto p = driverParameters(d);
        near(p.resistance, re);
        near(p.inductance, .0003);
        near(p.pistonArea, sd);
        near(p.forceFactor, bl);
        near(p.movingMass, mms - air);
        near(p.compliance, cms);
        near(p.mechanicalResistance, rms);
        d.thieleSmall.externalAirMass.reset();
        rejects([&] { (void)driverParameters(d); });
        d.thieleSmall.externalAirMass = mms * 2;
        rejects([&] { (void)driverParameters(d); });
        d.thieleSmall.externalAirMass = 0;
        near(driverParameters(d).movingMass, mms);
        rejects([&] { (void)driverParameters(d, 0, c); });
        d.entryMode = DriverEntryMode::electromechanical;
        d.electromechanical = p;
        require(driverParameters(d) == p, "Direct driver changed.");
        d.electromechanical.movingMass = 0;
        rejects([&] { (void)driverParameters(d); });
        const auto reference = scanSpeakD2608ReferenceDriver(9);
        near(driverParameters(reference).movingMass, .00016926810651601092);
        require(!reference.provenance.empty() && !reference.assumptions.empty(),
                "Reference provenance missing.");
        std::cout << "Driver model tests passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
