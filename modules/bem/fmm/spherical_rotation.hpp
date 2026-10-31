#pragma once
#include "wideband_transfer.hpp"
#include <cstddef>

namespace hornsculpt::wideband {
// LF RCR device contract only. HF angular signatures and their quadrature
// orders have separate limits and must not be restricted by this bound.
inline constexpr unsigned maximumDeviceScaledOrder = 128;
struct SphericalRotationPlan {
    unsigned order{};
    // Integer-l Wigner d(theta), row/column m=-l..l, row-major per degree.
    std::vector<std::size_t> degreeOffsets;
    std::vector<double> smallD;
    // exp(+i*m*phi), m=-order..order. Q=Rz(phi)*Ry(theta), Q*ez=axis.
    std::vector<Complex> azimuthPhases;
};
// Maximal-spin Clebsch-Gordan coupling builds each next half-integer degree
// from a unitary spin-1/2 rotation; no factorial/Racah alternating sum occurs.
// O(p^3) work/storage, with O(p^2) temporary state. This supplies the rotations
// in the rotation-coaxial-rotation factorization, not a dense p^4 translation.
SphericalRotationPlan makeSphericalRotationPlan(unsigned order, Vec3 axis);
// Pullback f(Q*r): c'=d(theta)^T*diag(exp(+im*phi))*c. The inverse restores c.
std::vector<Complex> rotateToAxis(const SphericalRotationPlan &plan,
                                  std::span<const Complex> values);
std::vector<Complex> rotateFromAxis(const SphericalRotationPlan &plan,
                                    std::span<const Complex> values);
} // namespace hornsculpt::wideband
