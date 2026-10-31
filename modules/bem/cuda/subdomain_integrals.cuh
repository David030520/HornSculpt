#pragma once

#include "panel_integrals.cuh"
#include <cfloat>

namespace hornsculpt::cuda_detail {
// Outward-air interior trace, including the unbounded exterior. Interface q
// uses P1; physical-panel normal flux is DP0.
struct SubdomainMoments {
    DoubleComplex pressure[9]{}, flux[9]{}, rhs[3]{}, green{};
};
__device__ inline void accumulateSubdomainMoment(SubdomainMoments &m, const Panel &x,
                                                 const Panel &y,
                                                 const galerkin::PairQuadraturePoint &p, double k,
                                                 bool coincident, bool interfaceFlux,
                                                 bool qualifiedMixed, int *error) {
    const auto delta = pairDelta(x, y, p.test, p.trial);
    DoubleComplex g{}, a{}, b{};
    bool mixed = false;
    if (qualifiedMixed && !coincident) {
        const float dx = static_cast<float>(delta.x), dy = static_cast<float>(delta.y),
                    dz = static_cast<float>(delta.z);
        const float radiusSquared = dx * dx + dy * dy + dz * dz;
        const float kf = static_cast<float>(k), inverseK = static_cast<float>(1.0 / k);
        const double physicalWeight = p.weight * 4 * x.area * y.area;
        const float weight = static_cast<float>(physicalWeight);
        // The shared mixed kernel assumes a normal finite FP32 squared radius,
        // representable scale and finite phase. Unsafe samples remain FP64.
        if (isfinite(radiusSquared) && radiusSquared >= FLT_MIN && isfinite(kf) && kf >= FLT_MIN &&
            isfinite(inverseK) && inverseK >= FLT_MIN && isfinite(weight) &&
            fabsf(weight) >= FLT_MIN) {
            const float inverseRadius = rsqrtf(radiusSquared);
            const float radius = radiusSquared * inverseRadius;
            const float phase = kf * radius;
            const float scale = weight * (0.07957747154594766788f * inverseRadius);
            if (isfinite(inverseRadius) && isfinite(phase) && isfinite(scale) &&
                fabsf(scale) >= FLT_MIN) {
                const auto sample =
                    mixedKernelSample(delta, x, y, weight, k, inverseK, false, error);
                if (isfinite(sample.green.x) && isfinite(sample.green.y)) {
                    const double r = length(delta);
                    const double inverseRadiusSquared = 1.0 / (r * r);
                    const double source = dot(delta, y.normal) * inverseRadiusSquared;
                    const double target = dot(delta, x.normal) * inverseRadiusSquared;
                    const double nn = dot(x.normal, y.normal);
                    const DoubleComplex pressureFactor{source, k * (r * source + nn)};
                    const DoubleComplex fluxFactor{1 + r * target, -target / k};
                    const auto green = promoted(sample.green);
                    const auto pressure = product(green, pressureFactor);
                    const auto flux = product(green, fluxFactor);
                    if (isfinite(r) && r > 0 && isfinite(pressure.real) &&
                        isfinite(pressure.imag) && isfinite(flux.real) && isfinite(flux.imag)) {
                        g = green;
                        a = pressure;
                        b = flux;
                        mixed = true;
                    }
                }
            }
        }
    }
    if (!mixed) {
        const double r = length(delta);
        if (!(r > 0) || !isfinite(r)) {
            atomicExch(error, 1);
            return;
        }
        double sine, cosine;
        sincos(k * r, &sine, &cosine);
        const double factor = p.weight * 4 * x.area * y.area / (12.56637061435917295385 * r);
        g = {cosine * factor, -sine * factor};
        const DoubleComplex radial{g.real - g.imag * k * r, g.imag + g.real * k * r};
        const auto d = times(radial, coincident ? 0 : dot(delta, y.normal) / (r * r));
        const auto kp = times(radial, coincident ? 0 : -dot(delta, x.normal) / (r * r));
        const double nn = dot(x.normal, y.normal);
        a = {d.real - k * nn * g.imag, d.imag + k * nn * g.real};
        b = {g.real - kp.imag / k, g.imag + kp.real / k};
    }
    m.green = plus(m.green, g);
    for (unsigned i = 0; i < 3; ++i) {
        if (!interfaceFlux)
            m.rhs[i] = plus(m.rhs[i], times(b, p.test[i]));
        for (unsigned j = 0; j < 3; ++j) {
            const auto basis = p.test[i] * p.trial[j];
            m.pressure[3 * i + j] = plus(m.pressure[3 * i + j], times(a, basis));
            if (interfaceFlux)
                m.flux[3 * i + j] = plus(m.flux[3 * i + j], times(b, basis));
        }
    }
}
} // namespace hornsculpt::cuda_detail
