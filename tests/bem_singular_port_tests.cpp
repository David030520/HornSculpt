#include <cuda_runtime_api.h>
#include <hornsculpt/bem.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace hornsculpt;
using Complex = std::complex<double>;
constexpr double relativePortTolerance = 5e-5;
constexpr double phaseToleranceDegrees = .02;

void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
bool finite(Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}
double relative(Complex actual, Complex expected) {
    require(finite(actual) && finite(expected) && std::abs(expected) > 0,
            "Nonfinite or zero singular-port comparison value.");
    return std::abs(actual - expected) / std::abs(expected);
}
struct PortReference {
    std::string excitation;
    std::uint32_t face{};
    double area{};
    Complex volumeVelocity, meanPressure, volumeImpedance;
};
std::vector<PortReference> readReference() {
    std::ifstream input(HS_SINGULAR_PORT_REFERENCE_PATH);
    require(input.good(), "Cannot read the pinned singular-port reference CSV.");
    std::vector<PortReference> result;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#')
            continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        fields.imbue(std::locale::classic());
        PortReference value;
        double qr{}, qi{}, pr{}, pi{}, zr{}, zi{};
        require(bool(fields >> value.excitation >> value.face >> value.area >> qr >> qi >> pr >>
                     pi >> zr >> zi),
                "Malformed singular-port reference row.");
        fields >> std::ws;
        require(fields.eof(), "Unexpected trailing singular-port reference fields.");
        value.volumeVelocity = {qr, qi};
        value.meanPressure = {pr, pi};
        value.volumeImpedance = {zr, zi};
        require(std::isfinite(value.area) && value.area > 0 &&
                    relative(value.meanPressure / value.volumeVelocity, value.volumeImpedance) <
                        1e-12,
                "Reference port does not satisfy its physical p/Q convention.");
        result.push_back(value);
    }
    require(input.eof() && result.size() == 4, "Expected exactly four physical reference ports.");
    return result;
}
SurfaceMesh tetrahedron() {
    SurfaceMesh mesh;
    // Metres; every pair of this convex closed tetrahedron's panels touches.
    // These outward windings and face identities match the independent reference.
    mesh.vertices = {{0, 0, 0}, {.2, 0, 0}, {.03, .14, 0}, {.02, .025, .09}};
    mesh.triangles = {{{0, 2, 1}, 50, BoundaryCondition::driven},
                      {{0, 1, 3}, 51, BoundaryCondition::driven},
                      {{0, 3, 2}, 52, BoundaryCondition::wall},
                      {{1, 2, 3}, 53, BoundaryCondition::wall}};
    return mesh;
}
void checkSingularPorts() {
    const auto reference = readReference();
    BemRequest request;
    request.medium = {1.21, 343};
    request.frequencies = {400};
    // RMS outward normal velocities, exp(+i omega t), outgoing exp(-ikr).
    // Fixed +/-1 m/s on unequal physical faces; no aggregate-flow impedance.
    request.excitations = {BemExcitation{"uniform", {{50, {1, 0}}, {51, {1, 0}}}},
                           BemExcitation{"antisymmetric", {{50, {1, 0}}, {51, {-1, 0}}}}};
    request.numerics.acceleration = BemAcceleration::fastMultipole;
    request.numerics.quadratureOrder = 4;
    request.numerics.singularOrder = 8;
    request.numerics.nearDistanceRatio = 1;
    request.numerics.nearRelativeTolerance = 1e-4;
    request.numerics.nearMaximumDepth = 6;
    request.numerics.relativeResidualTolerance = 5e-5;
    request.numerics.fastMultipole.leafSupportCapacity = 32;
    request.numerics.fastMultipole.preconditionerBlockWidth = 128;
    request.numerics.fastMultipole.separationRatio = 1.5;
    // Ports use exact P1 surface means. No observer-field reference is claimed.
    CudaBemSolver solver(tetrahedron());
    const auto result = solver.solve(request);
    require(result.acceleration == BemAcceleration::fastMultipole && result.triangleCount == 4 &&
                result.pressureDofCount == 4 && result.frequencies.size() == 1,
            "The public singular-port solve changed its physical scope.");
    const auto &frequency = result.frequencies.front();
    require(frequency.frequency == 400 && frequency.excitations.size() == 2 &&
                frequency.fastMultipole.farInteractionCount == 0,
            "The touching-panel fixture must retain both RHS and contain no far interactions.");
    for (std::size_t rhs = 0; rhs < request.excitations.size(); ++rhs) {
        const auto &expectedExcitation = request.excitations[rhs];
        const auto &actual = frequency.excitations[rhs];
        require(actual.name == expectedExcitation.name && actual.impedances.size() == 2 &&
                    actual.solvedVelocities.size() == 2,
                "The public solve lost a physical excitation or port.");
        require(actual.residualEvaluation == BemResidualEvaluation::directGalerkin &&
                    std::isfinite(actual.boundaryResidual) &&
                    actual.boundaryResidual <= request.numerics.relativeResidualTolerance,
                "The singular-port solve failed its direct original-equation residual.");
        for (const auto &velocity : expectedExcitation.faces) {
            require(std::count_if(actual.solvedVelocities.begin(), actual.solvedVelocities.end(),
                                  [&](const auto &v) { return v.faceId == velocity.faceId; }) == 1,
                    "Solved physical face velocity is missing or duplicated.");
            const auto solved =
                std::find_if(actual.solvedVelocities.begin(), actual.solvedVelocities.end(),
                             [&](const auto &v) { return v.faceId == velocity.faceId; });
            require(std::abs(solved->normalVelocity - velocity.normalVelocity) < 1e-12,
                    "The fixed physical normal velocity changed.");
            require(std::count_if(actual.impedances.begin(), actual.impedances.end(),
                                  [&](const auto &p) { return p.faceId == velocity.faceId; }) == 1,
                    "Actual physical port is missing or duplicated.");
            const auto port =
                std::find_if(actual.impedances.begin(), actual.impedances.end(),
                             [&](const auto &p) { return p.faceId == velocity.faceId; });
            const auto matches = [&](const auto &p) {
                return p.excitation == expectedExcitation.name && p.face == velocity.faceId;
            };
            require(std::count_if(reference.begin(), reference.end(), matches) == 1,
                    "Reference physical port is missing or duplicated.");
            const auto expected = std::find_if(reference.begin(), reference.end(), matches);
            require(std::isfinite(port->area) &&
                        std::abs(port->area / expected->area - 1) < 1e-12 &&
                        relative(port->volumeVelocity, expected->volumeVelocity) < 1e-12 &&
                        relative(expected->volumeVelocity,
                                 expected->area * velocity.normalVelocity) < 1e-12,
                    "Port area or volume velocity differs from the independent physical fixture.");
            require(relative(port->meanPressure / port->volumeVelocity, port->volumeImpedance) <
                        1e-12,
                    "The physical port impedance differs from p/Q.");
            const double pressureError = relative(port->meanPressure, expected->meanPressure),
                         impedanceError =
                             relative(port->volumeImpedance, expected->volumeImpedance),
                         phaseError =
                             std::abs(std::arg(port->volumeImpedance / expected->volumeImpedance)) *
                             180 / std::numbers::pi;
            std::cout << actual.name << " face=" << velocity.faceId << " pressure=" << pressureError
                      << " impedance=" << impedanceError << " phase_deg=" << phaseError
                      << " direct_residual=" << actual.boundaryResidual << '\n';
            require(pressureError <= relativePortTolerance &&
                        impedanceError <= relativePortTolerance &&
                        phaseError <= phaseToleranceDegrees,
                    "Singular-order8 physical port fails the refined CF64 reference.");
        }
    }
}
} // namespace

int main() {
    int count{};
    if (cudaGetDeviceCount(&count) != cudaSuccess || !count) {
        std::cout << "No CUDA device; singular-port numerical test skipped.\n";
        return 77;
    }
    try {
        std::cout << std::setprecision(12);
        checkSingularPorts();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
