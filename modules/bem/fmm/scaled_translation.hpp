#pragma once

#include "wideband_transfer.hpp"
#include <array>
#include <stop_token>

namespace hornsculpt::wideband {
struct ScalarWave {
    Complex value{};
    std::array<Complex, 3> gradient{};
};
struct ScaledWaveBasis {
    std::vector<ScalarWave> waves;
};
// r^l*Y_lm is a polynomial, so both the regular basis and its gradient have
// their correct finite limits at the expansion center. Normalized radial
// factors tend to one as k->0. Shared numerical limits are declared in
// wideband_transfer.hpp and checked by both the policy and these primitives.
ScaledWaveBasis scaledRegularWaves(unsigned order, double wavenumber, double radius, Vec3 point);
ScaledWaveBasis scaledOutgoingWaves(unsigned order, double wavenumber, double radius, Vec3 point);

struct ScaledMultipole {
    unsigned order{};
    double radius{}; // Basis normalization scale; physical support bounds are separate.
    // M_lm=integral density*(r/a)^l*Jhat_l(kr)*conj(Y_lm) dS.
    // A double-layer source uses the normal derivative of this regular basis.
    std::vector<Complex> coefficients;
};
struct ScaledLocal {
    unsigned order{};
    double radius{}; // Basis normalization scale, not a Cauchy probe radius.
    // u(x)=sum L_lm*(r/a)^l*Jhat_l(kr)*Y_lm(x/r).
    std::vector<Complex> coefficients;
};
ScalarWave evaluateScaledMultipole(const ScaledMultipole &source, double wavenumber, Vec3 point);
ScalarWave evaluateScaledLocal(const ScaledLocal &local, double wavenumber, Vec3 point);
std::vector<Complex> scaledMultipoleToSpectrum(const ScaledMultipole &source, double wavenumber,
                                               const SphereRule &rule);
ScaledLocal spectrumToScaledLocal(std::span<const Complex> incoming, const SphereRule &rule,
                                  unsigned order, double wavenumber, double radius);

enum class ScaledAxialKind { multipoleToParent, multipoleToLocal, localToChild };
struct ScaledAxialPlan {
    ScaledAxialKind kind{};
    unsigned sourceOrder{}, targetOrder{}, polarCount{};
    // Block m has rows l=m..targetOrder and columns n=m..sourceOrder,
    // row-major. The same block applies to +m and -m. No cross-m entries exist.
    std::vector<std::size_t> blockOffsets;
    std::vector<Complex> coefficients;
    double projectionRadius{};
};
// Axial displacement is positive +Z. Cauchy traces remove division by zeros of
// individual regular Bessel functions. One-time 1D-polar plan construction uses
// O(p^4) reference work; stored blocks and application require O(p^3).
// Normalization radii do not define integration surfaces or physical supports.
// An explicit M2L probe must be positive and smaller than axialDistance; the
// default is min(targetRadius, axialDistance/2). The finite outgoing basis is
// analytic away from its center, even inside its normalization sphere. Actual
// source/target support separation and local truncation convergence remain the
// caller's independent geometric and physical accuracy contracts. M2M defaults
// to a 1.01*sourceRadius Cauchy surface, and L2L to targetRadius.
ScaledAxialPlan makeScaledAxialPlan(ScaledAxialKind kind, unsigned sourceOrder,
                                    unsigned targetOrder, double wavenumber, double sourceRadius,
                                    double targetRadius, double axialDistance, unsigned polarCount,
                                    double projectionRadius = 0, std::stop_token stop = {});
std::vector<Complex> applyScaledAxialPlan(const ScaledAxialPlan &plan,
                                          std::span<const Complex> source);

} // namespace hornsculpt::wideband
