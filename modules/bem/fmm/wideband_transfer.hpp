#pragma once

#include "../galerkin/space.hpp"
#include <complex>
#include <span>
#include <vector>

namespace hornsculpt::wideband {

using Complex = std::complex<double>;
inline constexpr unsigned maximumScaledSphericalOrder = 256;
inline constexpr double maximumScaledElectricalRadius = 128;

struct SphereSample {
    Vec3 direction{};
    double weight{};
};
struct SphereRule {
    unsigned order{};
    std::vector<SphereSample> samples;
};

struct PolarSample {
    double cosine{}, weight{};
};
// Separable Fourier/associated-Legendre angular transform. Recurrence arrays
// use l*(l+1)/2+m, 0<=m<=l; no dense source-sample/target-sample matrix is stored.
// A GPU implementation can generate one active m-block in scratch and use an
// azimuth FFT. The CPU path below uses the same plan with explicit Fourier sums.
struct AngularTransferPlan {
    unsigned sourceAzimuthCount{}, targetAzimuthCount{}, bandlimit{};
    std::vector<PolarSample> sourcePolar, targetPolar;
    std::vector<double> recurrenceFirst, recurrenceSecond, diagonalScale;
};

AngularTransferPlan makeAngularTransferPlan(const SphereRule &source, const SphereRule &target,
                                            unsigned bandlimit);
std::vector<Complex> resampleSpectrum(std::span<const Complex> values,
                                      const AngularTransferPlan &plan);

// Gauss-Legendre in cos(theta), equispaced azimuth. Products of spherical
// harmonics through order L integrate exactly (up to floating-point rounding).
SphereRule makeSphereRule(unsigned order);

// Orthonormal complex Y_lm with Condon-Shortley phase. Index l*l+l+m,
// -l <= m <= l; directions must be finite unit vectors.
std::vector<Complex> sphericalHarmonics(unsigned order, Vec3 direction);

struct RegularRadialTable {
    std::vector<double> values, derivatives;
};
// j_l(x) and d/dx j_l(x), including x=0. Small arguments use a convergent
// power series; the turning/evanescent regime uses normalized Miller recursion.
// This is the regular radial basis for the low-frequency representation.
RegularRadialTable regularRadialTable(unsigned order, double argument);

enum class LevelRepresentation { inactive, regularSpherical, planeWave };
struct WidebandLevelPlan {
    LevelRepresentation representation{LevelRepresentation::inactive};
    double activeRadius{}, electricalRadius{}, accuracySeed{};
    unsigned signatureOrder{}, translationOrder{}, quadratureOrder{};
    std::size_t sphericalCoefficients{}, angularSamples{};
};
// Initial monopole bandwidth estimate of Gumerov/Duraiswami (2009), Eq. 18,
// using only the radius of active complete supports. This seed is not a bound
// on the Burton-Miller matrix or a solver acceptance criterion. Dipoles,
// derivatives, quadrature and final angular convergence need independent gates.
// The full 2p translation/quadrature rule is budgeted, not a p-grid shortcut.
// convergenceSeparation is b/a (center to closest remote support point / own
// support radius), NOT the public MAC multiplying the sum of both radii.
WidebandLevelPlan makeWidebandLevelPlan(double maximumActiveSupportRadius, double wavenumber,
                                        double convergenceSeparation, double accuracySeed);

// Outgoing convention: F(d) = 4*pi*sum_lm i^l M_lm Y_lm(d).
// Incoming convention: u(r) = sum_lm L_lm j_l(kr) Y_lm(r/|r|), where
// L_lm = 4*pi*(-i)^l integral a(d)*conj(Y_lm(d)) dOmega and
// u(r) = integral a(d)*exp(-ik*d.r) dOmega.
std::vector<Complex> multipoleToSpectrum(std::span<const Complex> multipoles, unsigned order,
                                         const SphereRule &rule);
std::vector<Complex> spectrumToLocal(std::span<const Complex> incoming, const SphereRule &rule,
                                     unsigned order);

// Bandlimited angular interpolation. Incoming down-transfer is the quadrature-
// weighted Hermitian adjoint of outgoing up-transfer. The opposite phase
// signs follow exp(+ik*d.(y-c)) for source signatures and exp(-ik*d.(x-c))
// for local fields, including unequal parent/child angular resolutions.
std::vector<Complex> resampleSpectrum(std::span<const Complex> values, const SphereRule &source,
                                      const SphereRule &target, unsigned bandlimit);
std::vector<Complex> outgoingToParent(std::span<const Complex> child, const SphereRule &childRule,
                                      const SphereRule &parentRule, unsigned childBandlimit,
                                      double wavenumber, Vec3 childMinusParent);
std::vector<Complex> incomingToChild(std::span<const Complex> parent, const SphereRule &parentRule,
                                     const SphereRule &childRule, unsigned childBandlimit,
                                     double wavenumber, Vec3 childMinusParent);

// Diagonal outgoing Helmholtz translation, exp(-ikr)/(4*pi*r):
// (-ik)/(16*pi^2) sum_l (2l+1)*(-i)^l*h_l^(2)(kR)*P_l(d.Rhat).
// R = targetCenter-sourceCenter. The caller must establish geometrical
// separation and convergence in expansionOrder; this formula alone is not an
// error certificate and is not a stable low-frequency translation.
std::vector<Complex> diagonalTranslation(const SphereRule &rule, unsigned expansionOrder,
                                         double wavenumber, Vec3 targetMinusSource);

} // namespace hornsculpt::wideband
