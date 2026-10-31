#include "fixture_geometry.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <hornsculpt/bem.hpp>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace hornsculpt;

namespace {
using Complex = std::complex<double>;
constexpr Complex imaginary{0, 1};
constexpr double radius = .1;
constexpr double soundSpeed = 343;
constexpr double density = 1.21;
constexpr Complex mixture{.3, .4};

// Gates are fixed before the first GPU run. They include the faceted-sphere
// geometry and DP0 approximation of translation velocity, as well as BEM error.
constexpr double coarseErrorLimit = .08;
constexpr double fineErrorLimit = .03;
constexpr double refinementFactor = .7;
constexpr double dipoleNullLimit = 2e-4;
constexpr double superpositionLimit = 5e-4;

void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}

double length(Vec3 value) {
    return std::sqrt(value.x * value.x + value.y * value.y + value.z * value.z);
}

Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

// Independent closed forms of outgoing spherical Hankel functions for exp(+iwt).
Complex hankel(unsigned mode, double z) {
    const auto phase = std::exp(-imaginary * z);
    if (mode == 0)
        return imaginary * phase / z;
    return (-1.0 / z + imaginary / (z * z)) * phase;
}

Complex hankelDerivative(unsigned mode, double z) {
    const auto phase = std::exp(-imaginary * z);
    if (mode == 0)
        return (1.0 / z - imaginary / (z * z)) * phase;
    return (imaginary / z + 2.0 / (z * z) - 2.0 * imaginary / (z * z * z)) * phase;
}

Complex analyticPressure(unsigned mode, double frequency, Vec3 observer) {
    const double omega = 2 * std::numbers::pi * frequency;
    const double k = omega / soundSpeed;
    const double r = length(observer);
    const double angular = mode == 0 ? 1.0 : observer.z / r;
    return (-imaginary * density * omega / k) * hankel(mode, k * r) /
           hankelDerivative(mode, k * radius) * angular;
}

void checkAnalyticOracle() {
    for (const double z :
         std::array{.4, std::numbers::pi, 4.493409457909064, 2 * std::numbers::pi}) {
        for (unsigned mode = 0; mode < 2; ++mode) {
            const double step = 1e-5;
            const auto difference = (hankel(mode, z + step) - hankel(mode, z - step)) / (2 * step);
            require(std::abs(difference / hankelDerivative(mode, z) - 1.0) < 2e-9,
                    "Independent Hankel derivative is inconsistent");
        }
        const double k = z / radius;
        const double frequency = k * soundSpeed / (2 * std::numbers::pi);
        const double r = .3;
        const auto monopole = density * soundSpeed * (imaginary * z) / (1.0 + imaginary * z) *
                              (radius / r) * std::exp(-imaginary * k * (r - radius));
        require(std::abs(analyticPressure(0, frequency, {0, 0, r}) / monopole - 1.0) < 1e-13,
                "Hankel monopole disagrees with the independent pulsating-sphere formula");
        require(analyticPressure(1, frequency, {r, 0, 0}) == Complex{},
                "Analytic dipole does not vanish in its equatorial plane");
    }
}

BemRequest makeRequest(SurfaceMesh &mesh) {
    BemRequest request;
    request.medium = {density, soundSpeed};
    request.observers = {{0, 0, .3}, {0, 0, -.3}, {0, 0, 1}, {0, 0, -1}, {.3, 0, 0}, {.3, .4, .5}};
    for (const double ka : std::array{std::numbers::pi, 4.493409457909064, 2 * std::numbers::pi})
        for (const double detuning : std::array{-.002, 0.0, .002})
            request.frequencies.push_back(ka * (1 + detuning) * soundSpeed /
                                          (2 * std::numbers::pi * radius));
    request.excitations = {{"radial", {}}, {"translation-z", {}}, {"mixed", {}}};
    for (std::size_t panel = 0; panel < mesh.triangles.size(); ++panel) {
        auto &triangle = mesh.triangles[panel];
        triangle.faceId = static_cast<std::uint32_t>(panel + 1);
        triangle.boundary = BoundaryCondition::driven;
        const auto a = mesh.vertices[triangle.vertices[0]];
        const auto b = mesh.vertices[triangle.vertices[1]];
        const auto c = mesh.vertices[triangle.vertices[2]];
        const auto normal = cross(subtract(b, a), subtract(c, a));
        const double normalZ = normal.z / length(normal);
        request.excitations[0].faces.push_back({triangle.faceId, 1.0});
        // Rigid translation in +Z has a spatially varying normal velocity.
        request.excitations[1].faces.push_back({triangle.faceId, normalZ});
        request.excitations[2].faces.push_back({triangle.faceId, 1.0 + mixture * normalZ});
    }
    return request;
}

struct ErrorSummary {
    std::array<double, 2> squaredFrequencyErrors{};
    std::array<double, 2> maximumErrors{};
    double worstNull{};
    double worstSuperposition{};
};

ErrorSummary checkResult(const BemResult &result, const BemRequest &request, double errorLimit) {
    require(result.frequencies.size() == request.frequencies.size(), "Incomplete resonance sweep");
    require(result.pressureDofCount == result.triangleCount / 2 + 2,
            "Incorrect P1 pressure space on the closed spherical mesh");
    ErrorSummary summary;
    for (std::size_t fi = 0; fi < request.frequencies.size(); ++fi) {
        const auto &frequency = result.frequencies[fi];
        require(frequency.frequency == request.frequencies[fi], "Resonance frequency mismatch");
        require(frequency.excitations.size() == 3, "Independent excitation result is missing");
        for (const auto &excitation : frequency.excitations) {
            require(excitation.observerPressure.size() == request.observers.size(),
                    "Incomplete spatial pressure result");
            require(std::isfinite(excitation.boundaryResidual) &&
                        excitation.boundaryResidual <= request.numerics.relativeResidualTolerance,
                    "Original Burton-Miller equation residual exceeds its tolerance");
            for (const auto value : excitation.observerPressure)
                require(std::isfinite(value.real()) && std::isfinite(value.imag()),
                        "Nonfinite resonant pressure");
        }
        for (unsigned mode = 0; mode < 2; ++mode) {
            double errorSquared = 0, referenceSquared = 0;
            for (std::size_t observer = 0; observer < request.observers.size(); ++observer) {
                const auto expected =
                    analyticPressure(mode, frequency.frequency, request.observers[observer]);
                errorSquared +=
                    std::norm(frequency.excitations[mode].observerPressure[observer] - expected);
                referenceSquared += std::norm(expected);
            }
            const double error = std::sqrt(errorSquared / referenceSquared);
            summary.squaredFrequencyErrors[mode] += error * error;
            summary.maximumErrors[mode] = std::max(summary.maximumErrors[mode], error);
            std::cout << result.triangleCount << " triangles, ka="
                      << 2 * std::numbers::pi * frequency.frequency * radius / soundSpeed
                      << ", mode=" << mode << ", complex spatial L2=" << error
                      << ", residual=" << frequency.excitations[mode].boundaryResidual << '\n';
            require(error <= errorLimit,
                    "Analytic mode " + std::to_string(mode) + " fails near an interior resonance");
        }

        const auto &radial = frequency.excitations[0].observerPressure;
        const auto &dipole = frequency.excitations[1].observerPressure;
        const auto &mixed = frequency.excitations[2].observerPressure;
        // A relative comparison to the zero analytic value would be undefined.
        const double axialScale =
            std::abs(analyticPressure(1, frequency.frequency, request.observers.front()));
        const double nullError = std::abs(dipole[4]) / axialScale;
        summary.worstNull = std::max(summary.worstNull, nullError);
        require(nullError <= dipoleNullLimit, "Translation dipole has a spurious equatorial field");
        double differenceSquared = 0, scaleSquared = 0;
        for (std::size_t observer = 0; observer < request.observers.size(); ++observer) {
            differenceSquared +=
                std::norm(mixed[observer] - radial[observer] - mixture * dipole[observer]);
            scaleSquared +=
                std::pow(std::abs(radial[observer]) + std::abs(mixture * dipole[observer]), 2);
        }
        const double superpositionError = std::sqrt(differenceSquared / scaleSquared);
        summary.worstSuperposition = std::max(summary.worstSuperposition, superpositionError);
        require(superpositionError <= superpositionLimit,
                "Spatially independent multiple right-hand sides fail complex superposition");
    }
    return summary;
}

ErrorSummary solveMesh(unsigned refinement, double errorLimit) {
    auto mesh = hornsculpt::testing::sphere(refinement);
    const auto request = makeRequest(mesh);
    CudaBemSolver solver(std::move(mesh));
    const auto result = solver.solve(request);
    const auto errors = checkResult(result, request, errorLimit);
    std::cout << result.triangleCount << " triangles: worst dipole-null=" << errors.worstNull
              << ", multiple-RHS superposition=" << errors.worstSuperposition << '\n';
    return errors;
}
} // namespace

int main() {
    try {
        std::cout << std::setprecision(10);
        checkAnalyticOracle();
        std::cout << "CUDA device: " << CudaBemSolver::deviceDescription() << '\n';
        const auto coarse = solveMesh(3, coarseErrorLimit);
        const auto fine = solveMesh(4, fineErrorLimit);
        for (unsigned mode = 0; mode < 2; ++mode) {
            const double ratio =
                std::sqrt(fine.squaredFrequencyErrors[mode] / coarse.squaredFrequencyErrors[mode]);
            std::cout << "Mode " << mode << ": resonance-sweep refinement error ratio " << ratio
                      << ", worst fine error " << fine.maximumErrors[mode] << '\n';
            require(ratio <= refinementFactor,
                    "Refinement fails to reduce independent resonant-mode pressure error");
        }
        std::cout << "Burton-Miller resonance, translation dipole and independent multiple-RHS "
                     "tests passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
