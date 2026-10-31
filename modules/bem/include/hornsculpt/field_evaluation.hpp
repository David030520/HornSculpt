#pragma once

#include <hornsculpt/bem.hpp>
#include <stop_token>

namespace hornsculpt {

enum class FieldPointStatus { air, boundary, excluded, quadratureNotConverged };
struct PressureFieldSample {
    Vec3 position;
    std::complex<double> pressure;
    FieldPointStatus status{FieldPointStatus::excluded};
};
struct PressureFieldSettings {
    // A sampling exclusion distance, not a geometry modification. Surface
    // limits require singular trace evaluation and are deliberately excluded.
    double minimumBoundaryDistance{0.0001}; // m
    double relativeQuadratureTolerance{1e-4};
    unsigned maximumQuadratureDepth{8};
};

// Finite-distance Helmholtz representation from actual P1 pressure and DP0
// outward normal derivative traces. Subdomain interfaces use their stored P1
// derivatives and domain-specific outward air normals. exp(+i wt), outgoing exp(-ikr). Exterior:
// Dp-Sq; interior air cavity: Sq-Dp. Baffled exterior: -2 S_aperture q.
// Solid regions and the non-radiating half-space outside the cavity are masked;
// no extrapolation from observer/far-field samples is performed. A stop request
// throws std::runtime_error; quadrature failures remain explicitly masked.
std::vector<PressureFieldSample>
evaluatePressureField(const SurfaceMesh &mesh, const BemRequest &request, double frequency,
                      const BemExcitationResult &solution, const std::vector<Vec3> &points,
                      PressureFieldSettings settings = {}, std::stop_token stop = {},
                      std::function<void(std::size_t, std::size_t)> progress = {});

} // namespace hornsculpt
