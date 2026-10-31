#pragma once

#include "panel_integrals.cuh"
#include <cstdint>
#include <cuda_runtime.h>
#include <hornsculpt/bem.hpp>
#include <memory>

namespace hornsculpt::cuda_detail {
// Global row and pressure-column indices for one oriented domain panel.
// fluxColumns contain the common P1 interface-flow unknowns, or -1 for
// a physical panel. fluxSign converts that unknown to this domain's q/k.
struct SubdomainPanel {
    Panel panel;
    int pressureColumns[3];
    int rows[3];
    int fluxColumns[3];
    double fluxSign;
};
struct SubdomainPlan;
class SubdomainPlanCache {
  public:
    explicit SubdomainPlanCache(const SurfaceMesh &mesh) : mesh_(&mesh) {}

  private:
    friend BemResult solveSubdomains(const SurfaceMesh &, Vec3, const BemRequest &,
                                     SubdomainPlanCache &, std::stop_token,
                                     std::function<void(const BemProgress &)>, bool);
    const SurfaceMesh *mesh_;
    std::shared_ptr<const SubdomainPlan> plan_;
};
void assembleSubdomainMass(const SubdomainPanel *panels, std::size_t panelCount,
                           std::size_t fluxStride, unsigned dimension, double waveNumber,
                           unsigned excitationCount, const cuDoubleComplex *knownFlux,
                           const double *rowScale, cuDoubleComplex *matrix, cuDoubleComplex *rhs,
                           int *error, cudaStream_t stream);
BemResult solveSubdomains(const SurfaceMesh &mesh, Vec3 worldOrigin, const BemRequest &request,
                          SubdomainPlanCache &planCache, std::stop_token stop,
                          std::function<void(const BemProgress &)> progress,
                          bool allowZero = false);
} // namespace hornsculpt::cuda_detail
