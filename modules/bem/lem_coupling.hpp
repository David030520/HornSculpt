#pragma once

#include <hornsculpt/bem.hpp>
#include <span>

namespace hornsculpt::bem_detail {
struct CouplingBasis {
    std::vector<BemExcitation> excitations;
    std::vector<double> areas; // Modal areas S_w=sum(w_f*A_f); uniform ports use w=1.
    // -1 denotes a panel outside the network ports.
    std::vector<int> panelPorts;
    // Normal velocity for physical Q=1 at each panel, w_f/S_w [1/m²].
    // This is also the pressure projection density, preserving boundary power.
    std::vector<double> panelVelocities;
};

enum class CouplingDomain { exteriorSolid, baffledAir };
// Unit physical flow entering BEM air per port. Polarity belongs to the LEM,
// so neither basis velocities nor radiation impedance include that sign.
CouplingBasis makeCouplingBasis(const lem::Network &network, const SurfaceMesh &mesh,
                                std::span<const double> panelAreas,
                                CouplingDomain domain = CouplingDomain::exteriorSolid);
std::vector<FaceVelocity> coupledVelocities(const CouplingBasis &basis,
                                            const lem::Solution &solution);
// Verifies the actual reconstructed pressure, including any final BEM refinement,
// against the pressure fed back through the network. Spatial cancellation does
// not artificially tighten the budget: the scale includes the mean |p|.
void validateCouplingClosure(const CouplingBasis &basis, const lem::Solution &solution,
                             std::span<const double> panelAreas,
                             std::span<const std::complex<double>> surfacePressure);
} // namespace hornsculpt::bem_detail
