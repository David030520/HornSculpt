#pragma once

#include "axisymmetric_space.hpp"
#include "ring_kernels.hpp"
#include <array>
#include <stop_token>

#if defined(__CUDACC__)
#define HS_MERIDIAN_HD __host__ __device__
#else
#define HS_MERIDIAN_HD
#endif

namespace hornsculpt::galerkin {
struct RingSegment {
    double radius{}, axial{}, dr{}, dz{}, length{}, nr{}, nz{};
    std::uint32_t first{}, second{};
};
struct RingRulePoint {
    double position{}, weight{}, logFirst{}, logLast{};
};
struct RingPairPatch {
    double firstU{}, lastU{1}, firstV{}, lastV{1};
    // 0=Cartesian, 1/2=self Duffy, 3/4=common-endpoint Duffy,
    // 5=disjoint projected sinh. Both Duffy sectors are always included.
    unsigned mapping{}, depth{};
};
struct RingPairPoint {
    RingGeometry geometry;
    double test[2]{}, source[2]{};
    double weight{}, logCorrection{}, logFirstWeight{}, logLastWeight{};
    unsigned status{};
};
struct RingPairMoments {
    // Four kernels times two test and two source P1 functions, component-major.
    RingComplex value[16]{};
    double angularError[16]{};
};
struct RingPairRange {
    std::uint32_t first{}, count{};
};
struct RingInteractionPlan {
    std::vector<RingSegment> segments;
    // Ordered test/source index, including self. Ranges are geometric only;
    // actual embedded quadrature/error control takes place at each frequency.
    std::vector<RingPairRange> pairs;
    std::vector<RingPairPatch> roots;
    std::vector<RingRulePoint> regular, coarse, fine;
    unsigned regularOrder{}, coarseOrder{}, fineOrder{};
};

std::vector<RingRulePoint> ringRule(unsigned order);
RingInteractionPlan makeRingInteractionPlan(const AxisymmetricSpace &space,
                                            const BemNumerics &numerics, std::stop_token stop = {});

namespace ring_pair_detail {
HS_MERIDIAN_HD inline double dot(double ar, double az, double br, double bz) {
    return ar * br + az * bz;
}
HS_MERIDIAN_HD inline double clamped(double x) {
    return x < 0 ? 0 : (x > 1 ? 1 : x);
}
struct SinhPoint {
    double position{}, jacobian{1};
    unsigned status{};
};
HS_MERIDIAN_HD inline SinhPoint projectedSinh(double xr, double xz, double dr, double dz,
                                              double length, double parameter) {
    const double center = dot(xr, xz, dr, dz) / (length * length);
    const double closest = clamped(center);
    const double distance = ::hypot(xr - closest * dr, xz - closest * dz);
    const double width = distance / length;
    if (!(width > 0))
        return {0, 0, 1};
    const double low = ::asinh(-center / width), high = ::asinh((1 - center) / width);
    const double transformed = low + (high - low) * parameter;
    return {center + width * ::sinh(transformed), width * (high - low) * ::cosh(transformed), 0};
}
HS_MERIDIAN_HD inline unsigned sharedEnd(const RingSegment &a, const RingSegment &b) {
    return a.first == b.first || a.first == b.second ? 0 : 1;
}
} // namespace ring_pair_detail

HS_MERIDIAN_HD inline RingPairPoint mapRingPairPoint(const RingSegment &a, const RingSegment &b,
                                                     const RingPairPatch &patch, RingRulePoint gu,
                                                     RingRulePoint gv) {
    using namespace ring_pair_detail;
    RingPairPoint point;
    const double du = patch.lastU - patch.firstU, dv = patch.lastV - patch.firstV;
    const double u = patch.firstU + du * gu.position, rawV = patch.firstV + dv * gv.position;
    double v = rawV, x = u, y = rawV, jacobian = 1;
    bool logU = false, logV = false;
    const bool self = patch.mapping == 1 || patch.mapping == 2;
    const bool shared = patch.mapping == 3 || patch.mapping == 4;
    const bool flipSelf = self && a.radius != 0 && a.radius + a.dr == 0;
    unsigned endA = 0, endB = 0;
    double outAr = a.dr, outAz = a.dz, outBr = b.dr, outBz = b.dz;
    double distanceA = 0, distanceB = 0;
    if (shared) {
        endA = sharedEnd(a, b);
        endB = sharedEnd(b, a);
        if (endA) {
            outAr = -outAr;
            outAz = -outAz;
        }
        if (endB) {
            outBr = -outBr;
            outBz = -outBz;
        }
        const auto transformed = patch.mapping == 3
                                     ? projectedSinh(outAr, outAz, outBr, outBz, b.length, rawV)
                                     : projectedSinh(outBr, outBz, outAr, outAz, a.length, rawV);
        if (transformed.status) {
            point.status = transformed.status;
            return point;
        }
        v = transformed.position;
        jacobian = u * transformed.jacobian;
        distanceA = patch.mapping == 3 ? u : u * v;
        distanceB = patch.mapping == 3 ? u * v : u;
        x = endA ? 1 - distanceA : distanceA;
        y = endB ? 1 - distanceB : distanceB;
        logU = patch.firstU == 0;
    } else if (self) {
        x = patch.mapping == 1 ? u : u * v;
        y = patch.mapping == 1 ? u * v : u;
        if (flipSelf) {
            x = 1 - x;
            y = 1 - y;
        }
        jacobian = u;
        const bool axisAnchor = a.radius == 0 || flipSelf;
        logU = !axisAnchor && patch.firstU == 0;
        logV = patch.lastV == 1;
    } else if (patch.mapping == 5) {
        const double xr = a.radius - b.radius + x * a.dr;
        const double xz = a.axial - b.axial + x * a.dz;
        const auto transformed = projectedSinh(xr, xz, b.dr, b.dz, b.length, rawV);
        if (transformed.status) {
            point.status = transformed.status;
            return point;
        }
        y = transformed.position;
        jacobian = transformed.jacobian;
    }
    point.test[0] = 1 - x;
    point.test[1] = x;
    point.source[0] = 1 - y;
    point.source[1] = y;
    auto &g = point.geometry;
    g.testRadius = a.radius + x * a.dr;
    g.sourceRadius = b.radius + y * b.dr;
    g.testRadialNormal = a.nr;
    g.testAxialNormal = a.nz;
    g.sourceRadialNormal = b.nr;
    g.sourceAxialNormal = b.nz;
    if (self) {
        g.deltaRadius = (x - y) * a.dr;
        g.deltaZ = (x - y) * a.dz;
        if (a.radius == 0) {
            g.testRadius = x * a.dr;
            g.sourceRadius = y * a.dr;
        } else if (flipSelf) {
            const double axisX = patch.mapping == 1 ? u : u * v;
            const double axisY = patch.mapping == 1 ? u * v : u;
            g.testRadius = -axisX * a.dr;
            g.sourceRadius = -axisY * a.dr;
            g.deltaRadius = -(axisX - axisY) * a.dr;
            g.deltaZ = -(axisX - axisY) * a.dz;
        }
        g.testNormalDistance = g.sourceNormalDistance = 0;
    } else if (shared) {
        g.deltaRadius = distanceA * outAr - distanceB * outBr;
        g.deltaZ = distanceA * outAz - distanceB * outBz;
        g.testNormalDistance = -distanceB * dot(a.nr, a.nz, outBr, outBz);
        g.sourceNormalDistance = -distanceA * dot(b.nr, b.nz, outAr, outAz);
    } else {
        g.deltaRadius = a.radius - b.radius + x * a.dr - y * b.dr;
        g.deltaZ = a.axial - b.axial + x * a.dz - y * b.dz;
        g.testNormalDistance = dot(a.nr, a.nz, a.radius - b.radius, a.axial - b.axial) -
                               y * dot(a.nr, a.nz, b.dr, b.dz);
        g.sourceNormalDistance = dot(b.nr, b.nz, b.radius - a.radius, b.axial - a.axial) -
                                 x * dot(b.nr, b.nz, a.dr, a.dz);
    }
    const double common = 2 * ring_detail::pi * g.testRadius * g.sourceRadius * a.length *
                          b.length * du * dv * jacobian;
    point.weight = common * gu.weight * gv.weight;
    if (logU) {
        point.logCorrection += ::log(u) - ::log(du);
        point.logFirstWeight = common * gu.logFirst * gv.weight;
    }
    if (logV) {
        point.logCorrection += ::log(1 - v) - ::log(dv);
        point.logLastWeight = common * gu.weight * gv.logLast;
    }
    if (!ring_detail::finite(point.weight))
        point.status = 2;
    return point;
}

HS_MERIDIAN_HD inline unsigned accumulateRingPairPoint(RingPairMoments &moments,
                                                       const RingPairPoint &point, double k,
                                                       double angularTolerance) {
    if (point.status)
        return point.status;
    const auto kernels = ringKernels(point.geometry, k, angularTolerance);
    if (kernels.status)
        return kernels.status;
    RingStaticSplit split;
    const bool correction = point.logFirstWeight != 0 || point.logLastWeight != 0;
    if (correction) {
        split = ringStaticSplit(point.geometry);
        if (split.status)
            return split.status;
    }
    for (unsigned c = 0; c < 4; ++c) {
        auto value = ringScale(kernels.value[c], point.weight);
        if (correction)
            value.real += split.logarithmic[c] * (point.weight * point.logCorrection -
                                                  point.logFirstWeight - point.logLastWeight);
        for (unsigned i = 0; i < 2; ++i)
            for (unsigned j = 0; j < 2; ++j) {
                const unsigned index = 4 * c + 2 * i + j;
                const double basis = point.test[i] * point.source[j];
                moments.value[index] = ringAdd(moments.value[index], ringScale(value, basis));
                moments.angularError[index] +=
                    ::fabs(point.weight * basis) * kernels.angularAbsoluteError[c];
            }
    }
    return 0;
}

HS_MERIDIAN_HD inline RingPairMoments ringSubtract(RingPairMoments a, const RingPairMoments &b) {
    for (unsigned i = 0; i < 16; ++i) {
        a.value[i].real -= b.value[i].real;
        a.value[i].imag -= b.value[i].imag;
        a.angularError[i] += b.angularError[i];
    }
    return a;
}
HS_MERIDIAN_HD inline void ringAccumulate(RingPairMoments &a, const RingPairMoments &b) {
    for (unsigned i = 0; i < 16; ++i) {
        a.value[i] = ringAdd(a.value[i], b.value[i]);
        a.angularError[i] += b.angularError[i];
    }
}

// Operator-dimensional control: S/diameter, D, K', diameter*H. In particular
// the Maue curl coefficient is 1/(L_test L_source), never a tangent dot product.
HS_MERIDIAN_HD inline double ringOperatorNorm(const RingPairMoments &m, const RingSegment &a,
                                              const RingSegment &b, double k, bool errors = false) {
    const double size = a.length > b.length ? a.length : b.length;
    double cosine = 0;
    RingComplex totalCosine{};
    for (unsigned j = 0; j < 4; ++j) {
        cosine += m.angularError[4 + j];
        totalCosine = ringAdd(totalCosine, m.value[4 + j]);
    }
    double norm = 0;
    for (unsigned i = 0; i < 2; ++i)
        for (unsigned j = 0; j < 2; ++j) {
            const unsigned q = 2 * i + j;
            const double ci = (i ? -1.0 : 1.0) / a.length;
            const double cj = (j ? -1.0 : 1.0) / b.length;
            if (errors) {
                norm = ring_detail::maximum(norm, m.angularError[q] / size);
                norm = ring_detail::maximum(norm, m.angularError[8 + q]);
                norm = ring_detail::maximum(norm, m.angularError[12 + q]);
                const double h =
                    ::fabs(ci * cj) * cosine + k * k *
                                                   (::fabs(a.nr * b.nr) * m.angularError[4 + q] +
                                                    ::fabs(a.nz * b.nz) * m.angularError[q]);
                norm = ring_detail::maximum(norm, size * h);
            } else {
                norm = ring_detail::maximum(norm, ringMagnitude(m.value[q]) / size);
                norm = ring_detail::maximum(norm, ringMagnitude(m.value[8 + q]));
                norm = ring_detail::maximum(norm, ringMagnitude(m.value[12 + q]));
                const auto h = ringAdd(ringScale(totalCosine, ci * cj),
                                       ringAdd(ringScale(m.value[4 + q], -k * k * a.nr * b.nr),
                                               ringScale(m.value[q], -k * k * a.nz * b.nz)));
                norm = ring_detail::maximum(norm, size * ringMagnitude(h));
            }
        }
    return norm;
}

struct RingPairIntegral {
    RingPairMoments moments;
    double estimatedRelativeError{};
    unsigned patches{}, maximumDepth{}, status{};
};
// Independent host execution of the same dimensional quadrature contract,
// used for CPU/CUDA parity and numerical diagnostic tests.
RingPairIntegral integrateRingPair(const RingSegment &test, const RingSegment &source,
                                   const std::vector<RingPairPatch> &roots,
                                   const std::vector<RingRulePoint> &coarse,
                                   const std::vector<RingRulePoint> &fine, double k,
                                   double tolerance, unsigned maximumDepth,
                                   std::stop_token stop = {});
} // namespace hornsculpt::galerkin

#undef HS_MERIDIAN_HD
