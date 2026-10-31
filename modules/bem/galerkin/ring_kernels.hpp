#pragma once

#include <cfloat>
#include <cmath>

#if defined(__CUDACC__)
#define HS_RING_HD __host__ __device__
#else
#define HS_RING_HD
#endif

namespace hornsculpt::galerkin {
// Positive-time convention, G=exp(-ikR)/(4*pi*R). The angular integral is
// unnormalised: g0=int G dpsi, g1=int G cos(psi) dpsi over a full revolution.
struct RingComplex {
    double real{}, imag{};
};
HS_RING_HD inline RingComplex ringAdd(RingComplex a, RingComplex b) {
    return {a.real + b.real, a.imag + b.imag};
}
HS_RING_HD inline RingComplex ringScale(RingComplex a, double b) {
    return {a.real * b, a.imag * b};
}
HS_RING_HD inline double ringMagnitude(RingComplex a) {
    return ::hypot(a.real, a.imag);
}
inline constexpr unsigned ringGreen = 0, ringCosine = 1, ringTargetNormal = 2, ringSourceNormal = 3;
struct RingGeometry {
    double testRadius{}, sourceRadius{}, deltaRadius{}, deltaZ{};
    double testRadialNormal{}, testAxialNormal{}, sourceRadialNormal{}, sourceAxialNormal{};
    // Stable meridional n_x.(x-y), n_y.(y-x). For a common straight interval
    // both are exactly zero; do not subtract large derivative components later.
    double testNormalDistance{}, sourceNormalDistance{};
};
struct RingKernel {
    RingComplex value[4]{};
    double angularError{};
    double angularAbsoluteError[4]{};
    unsigned evaluations{};
    // 0=success, 1=invalid/coincident geometry, 2=nonfinite, 3=angular limit.
    unsigned status{};
};
struct RingStaticSplit {
    // Static kernel = logarithmic[i]*log(4/sqrt(q)) + regular[i]. The smooth
    // cutoff in q only redistributes these two EXACTLY added terms.
    double logarithmic[4]{}, regular[4]{}, q{}, logFactor{};
    unsigned status{};
};

namespace ring_detail {
inline constexpr double pi = 3.141592653589793238462643383279502884;
HS_RING_HD inline bool finite(double x) {
    return x == x && ::fabs(x) <= DBL_MAX;
}
struct Elliptic {
    double first{}, second{};
    unsigned status{};
};
// AGM / complementary expansions, NIST DLMF 19.8 and 19.12. The input is
// q=1-m computed geometrically, never the rounded subtraction 1-m near m=1.
HS_RING_HD inline Elliptic completeElliptic(double q) {
    if (!(q > 0 && q <= 1))
        return {0, 0, 1};
    if (q < 1e-8) {
        const double l = ::log(4.0) - .5 * ::log(q), q2 = q * q;
        return {l + .25 * q * (l - 1) + (9.0 / 64) * q2 * (l - 7.0 / 6),
                1 + .5 * q * (l - .5) + (3.0 / 16) * q2 * (l - 13.0 / 12), 0};
    }
    double a = 1, b = ::sqrt(q), sum = .5 * (1 - q), power = 1;
    for (unsigned iteration = 0; iteration < 32; ++iteration) {
        const double c = .5 * (a - b);
        sum += power * c * c;
        const double next = .5 * (a + b);
        b = ::sqrt(a * b);
        a = next;
        power *= 2;
        if (::fabs(c) <= 2e-16 * a) {
            const double k = pi / (2 * a);
            return {k, k * (1 - sum), 0};
        }
    }
    return {0, 0, 2};
}

HS_RING_HD inline double maximum(double a, double b) {
    return a > b ? a : b;
}
HS_RING_HD inline double scale(const RingGeometry &g) {
    return ::sqrt(g.deltaRadius * g.deltaRadius + g.deltaZ * g.deltaZ +
                  4 * g.testRadius * g.sourceRadius);
}

HS_RING_HD inline void staticKernels(const RingGeometry &g, double *values, unsigned &status) {
    const double r = g.testRadius, t = g.sourceRadius, z = g.deltaZ;
    const double a = r * r + t * t + z * z, b = 2 * r * t;
    const double d2 = g.deltaRadius * g.deltaRadius + z * z, total2 = d2 + 4 * r * t;
    if (!(r >= 0 && t >= 0 && d2 > 0 && finite(total2))) {
        status = 1;
        return;
    }
    const double d = ::sqrt(total2), m = 4 * r * t / total2, q = d2 / total2;
    if (m < .01) {
        // Direct even-power angular averages avoid loss of g1 and normal
        // derivatives in K/E differences near either symmetry-axis endpoint.
        double moments[17]{};
        moments[0] = 1;
        for (unsigned j = 2; j <= 16; j += 2)
            moments[j] = moments[j - 2] * (j - 1.0) / j;
        const double e = b / a;
        double coefficient = 1, power = 1, s0 = 0, s1 = 0, d0 = 0, d1 = 0;
        for (unsigned j = 0; j <= 14; ++j) {
            s0 += coefficient * power * moments[j];
            s1 += coefficient * power * moments[j + 1];
            d0 += (2 * j + 1) * coefficient * power * moments[j];
            d1 += (2 * j + 1) * coefficient * power * moments[j + 1];
            coefficient *= (2 * j + 1.0) / (2 * j + 2.0);
            power *= e;
        }
        const double root = ::sqrt(a), f = -.5 / (a * root);
        values[0] = .5 * s0 / root;
        values[1] = .5 * s1 / root;
        values[2] = f * (g.testNormalDistance * d0 + g.testRadialNormal * t * (d0 - d1));
        values[3] = f * (g.sourceNormalDistance * d0 + g.sourceRadialNormal * r * (d0 - d1));
        return;
    }
    const auto ell = completeElliptic(q);
    if (ell.status) {
        status = ell.status;
        return;
    }
    const double k = ell.first, e = ell.second;
    values[0] = k / (pi * d);
    values[1] = ((2 - m) * k - 2 * e) / (pi * d * m);
    const double testSum = g.testNormalDistance + 2 * t * g.testRadialNormal;
    const double sourceSum = g.sourceNormalDistance + 2 * r * g.sourceRadialNormal;
    const double denominator = pi * d * total2 * m;
    values[2] = ((k - e / q) * g.testNormalDistance + (e - k) * testSum) / denominator;
    values[3] = ((k - e / q) * g.sourceNormalDistance + (e - k) * sourceSum) / denominator;
}

struct RulePoint {
    double position{}, weight{};
};
HS_RING_HD inline RulePoint gauss(unsigned order, unsigned index) {
    const bool mirror = index >= order / 2;
    const unsigned i = mirror ? order - index - 1 : index;
    RulePoint p{};
    if (order == 8) {
        switch (i) {
        case 0:
            p = {.019855071751231856, .05061426814518839};
            break;
        case 1:
            p = {.10166676129318658, .11119051722668723};
            break;
        case 2:
            p = {.2372337950418355, .1568533229389437};
            break;
        default:
            p = {.4082826787521751, .181341891689181};
            break;
        }
    } else {
        switch (i) {
        case 0:
            p = {.005299532504175031, .013576229705875926};
            break;
        case 1:
            p = {.0277124884633837, .031126761969323888};
            break;
        case 2:
            p = {.06718439880608412, .04757925584124645};
            break;
        case 3:
            p = {.1222977958224985, .06231448562776697};
            break;
        case 4:
            p = {.19106187779867811, .07479799440828841};
            break;
        case 5:
            p = {.2709916111713863, .08457825969750121};
            break;
        case 6:
            p = {.35919822461037054, .0913017075224618};
            break;
        default:
            p = {.4524937450811813, .09472530522753424};
            break;
        }
    }
    if (mirror)
        p.position = 1 - p.position;
    return p;
}

HS_RING_HD inline void remainderSample(const RingGeometry &g, double k, double psi,
                                       RingComplex *values) {
    const double h = ::sin(.5 * psi), oneMinusCosine = 2 * h * h;
    const double d2 = g.deltaRadius * g.deltaRadius + g.deltaZ * g.deltaZ;
    const double radius = ::sqrt(d2 + 2 * g.testRadius * g.sourceRadius * oneMinusCosine);
    const double x = k * radius, x2 = x * x;
    RingComplex potential{}, derivative{};
    if (::fabs(x) < 1e-3) {
        potential = {x2 * (-.5 + x2 * (1.0 / 24 - x2 / 720)) / radius,
                     x * (-1 + x2 * (1.0 / 6 - x2 / 120)) / radius};
        derivative = {x2 * (.5 + x2 * (-.125 + x2 / 144)),
                      x * x2 * (-1.0 / 3 + x2 * (1.0 / 30 - x2 / 840))};
    } else {
        const double s = ::sin(x), c = ::cos(x);
        potential = {(c - 1) / radius, -s / radius};
        derivative = {c + x * s - 1, x * c - s};
    }
    values[0] = potential;
    values[1] = ringScale(potential, 1 - oneMinusCosine);
    const double testDistance =
        g.testNormalDistance + g.testRadialNormal * g.sourceRadius * oneMinusCosine;
    const double sourceDistance =
        g.sourceNormalDistance + g.sourceRadialNormal * g.testRadius * oneMinusCosine;
    values[2] = ringScale(derivative, -testDistance / (radius * radius * radius));
    values[3] = ringScale(derivative, -sourceDistance / (radius * radius * radius));
}
struct AngularTransform {
    double width{}, stretch{};
    bool clustered{};
};
HS_RING_HD inline AngularTransform angularTransform(const RingGeometry &g) {
    AngularTransform transform;
    const double radialScale = ::sqrt(g.testRadius) * ::sqrt(g.sourceRadius);
    transform.width = ::hypot(g.deltaRadius, g.deltaZ) / radialScale;
    transform.clustered = transform.width < pi;
    if (transform.clustered)
        transform.stretch = ::asinh(pi / transform.width);
    return transform;
}
HS_RING_HD inline void remainderIntegral(const RingGeometry &g, const AngularTransform &transform,
                                         double k, double first, double last, unsigned order,
                                         RingComplex *values) {
    for (unsigned j = 0; j < order; ++j) {
        const auto p = gauss(order, j);
        RingComplex point[4]{};
        const double coordinate = first + (last - first) * p.position;
        const double psi = transform.clustered
                               ? transform.width * ::sinh(transform.stretch * coordinate)
                               : pi * coordinate;
        const double jacobian = transform.clustered ? transform.width * transform.stretch *
                                                          ::cosh(transform.stretch * coordinate)
                                                    : pi;
        remainderSample(g, k, psi, point);
        // Evenness gives 2/(4*pi), including both signs of the azimuth.
        const double weight = p.weight * (last - first) * jacobian / (2 * pi);
        for (unsigned c = 0; c < 4; ++c)
            values[c] = ringAdd(values[c], ringScale(point[c], weight));
    }
}
} // namespace ring_detail

// Exact Laplace split. A C3 cutoff keeps the logarithmic coefficient confined
// to q<1/4; it is not a geometry or accuracy threshold. The added-back regular
// part is defined by subtraction, so the identity also holds in the transition.
HS_RING_HD inline RingStaticSplit ringStaticSplit(const RingGeometry &g) {
    using namespace ring_detail;
    RingStaticSplit result;
    double full[4]{};
    staticKernels(g, full, result.status);
    if (result.status)
        return result;
    const double d2 = g.deltaRadius * g.deltaRadius + g.deltaZ * g.deltaZ;
    const double total = scale(g), total2 = total * total;
    result.q = d2 / total2;
    result.logFactor = ::log(4.0) - .5 * ::log(result.q);
    for (unsigned i = 0; i < 4; ++i)
        result.regular[i] = full[i];
    if (result.q >= .25)
        return result;
    double a = 1, cOverQ = 0, coefficient = 1, power = 1;
    for (unsigned n = 1; n <= 48; ++n) {
        coefficient *= ((n - .5) / n) * ((n - .5) / n);
        cOverQ += coefficient * (2.0 * n / (2 * n - 1)) * power;
        power *= result.q;
        const double term = coefficient * power;
        a += term;
        if (::fabs(term) < 1e-18)
            break;
    }
    double cutoff = 1;
    if (result.q > .125) {
        const double t = (result.q - .125) / .125, t2 = t * t;
        cutoff -= t2 * t2 * (35 + t * (-84 + t * (70 - 20 * t)));
    }
    const double m = 4 * g.testRadius * g.sourceRadius / total2;
    const double c = result.q * cOverQ, divisor = pi * total * total2 * m;
    const double testSum = g.testNormalDistance + 2 * g.sourceRadius * g.testRadialNormal;
    const double sourceSum = g.sourceNormalDistance + 2 * g.testRadius * g.sourceRadialNormal;
    result.logarithmic[0] = cutoff * a / (pi * total);
    result.logarithmic[1] = cutoff * ((2 - m) * a - 2 * c) / (pi * total * m);
    result.logarithmic[2] =
        cutoff * ((a - cOverQ) * g.testNormalDistance + (c - a) * testSum) / divisor;
    result.logarithmic[3] =
        cutoff * ((a - cOverQ) * g.sourceNormalDistance + (c - a) * sourceSum) / divisor;
    for (unsigned i = 0; i < 4; ++i)
        result.regular[i] = full[i] - result.logarithmic[i] * result.logFactor;
    return result;
}

// All four kernels share one adaptive angular rule. Its absolute error scale
// includes the static part and scales normal derivatives by pair diameter.
// Failure is explicit; there is no fixed-order or accuracy fallback.
HS_RING_HD inline RingKernel ringKernels(const RingGeometry &g, double k,
                                         double relativeTolerance = 1e-7) {
    using namespace ring_detail;
    RingKernel result;
    double statics[4]{};
    staticKernels(g, statics, result.status);
    if (result.status || !(k >= 0 && relativeTolerance > 0) || !finite(k) ||
        !finite(relativeTolerance)) {
        if (!result.status)
            result.status = 1;
        return result;
    }
    const double size = scale(g);
    double reference[4]{};
    for (unsigned c = 0; c < 4; ++c) {
        result.value[c].real = statics[c];
        reference[c] = maximum(::fabs(statics[c]), c < 2 ? 1 / size : 1 / (size * size));
    }
    if (k == 0)
        return result;
    if (g.testRadius == 0 || g.sourceRadius == 0) {
        // On the axis the entire azimuth has identical distance. Avoid angular
        // cancellation in g1 and unnecessary quadrature in every axis field.
        const double distance = ::sqrt(g.deltaRadius * g.deltaRadius + g.deltaZ * g.deltaZ);
        const double phase = k * distance, c = ::cos(phase), s = ::sin(phase);
        result.value[0] = {.5 * c / distance, -.5 * s / distance};
        result.value[1] = {};
        const RingComplex radial{-.5 * (c + phase * s) / (distance * distance * distance),
                                 -.5 * (phase * c - s) / (distance * distance * distance)};
        result.value[2] =
            ringScale(radial, g.testNormalDistance + g.testRadialNormal * g.sourceRadius);
        result.value[3] =
            ringScale(radial, g.sourceNormalDistance + g.sourceRadialNormal * g.testRadius);
        return result;
    }
    const auto transform = angularTransform(g);
    struct Interval {
        double first, last;
        unsigned depth;
    };
    Interval stack[18]{};
    unsigned count = 1;
    stack[0] = {0, 1, 0};
    while (count) {
        const auto interval = stack[--count];
        RingComplex low[4]{}, high[4]{};
        remainderIntegral(g, transform, k, interval.first, interval.last, 8, low);
        remainderIntegral(g, transform, k, interval.first, interval.last, 16, high);
        result.evaluations += 24;
        double error = 0;
        for (unsigned c = 0; c < 4; ++c) {
            const RingComplex difference{high[c].real - low[c].real, high[c].imag - low[c].imag};
            error = maximum(error, ringMagnitude(difference) / reference[c]);
        }
        if (!finite(error)) {
            result.status = 2;
            return result;
        }
        const double budget = relativeTolerance * (interval.last - interval.first);
        if (error <= budget) {
            result.angularError += error;
            for (unsigned c = 0; c < 4; ++c) {
                result.angularAbsoluteError[c] +=
                    ringMagnitude({high[c].real - low[c].real, high[c].imag - low[c].imag});
                result.value[c] = ringAdd(result.value[c], high[c]);
            }
        } else {
            if (interval.depth >= 16 || count + 2 > 18 || result.evaluations > 100000) {
                result.status = 3;
                return result;
            }
            const double mid = .5 * (interval.first + interval.last);
            stack[count++] = {mid, interval.last, interval.depth + 1};
            stack[count++] = {interval.first, mid, interval.depth + 1};
        }
    }
    for (unsigned c = 0; c < 4; ++c)
        if (!finite(result.value[c].real) || !finite(result.value[c].imag))
            result.status = 2;
    return result;
}
} // namespace hornsculpt::galerkin

#undef HS_RING_HD
