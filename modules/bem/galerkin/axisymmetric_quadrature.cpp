#include "axisymmetric_quadrature.hpp"
#include "pair_quadrature.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace hornsculpt::galerkin {
namespace {
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Meridian integration cancelled.");
}
double pointDistance(double xr, double xz, const RingSegment &segment) {
    const double dr = xr - segment.radius, dz = xz - segment.axial;
    const double t = std::clamp(
        (dr * segment.dr + dz * segment.dz) / (segment.length * segment.length), 0.0, 1.0);
    return std::hypot(dr - t * segment.dr, dz - t * segment.dz);
}
double segmentDistance(const RingSegment &a, const RingSegment &b) {
    const double determinant = a.dr * b.dz - a.dz * b.dr;
    const double xr = b.radius - a.radius, xz = b.axial - a.axial;
    if (determinant != 0) {
        const double u = (xr * b.dz - xz * b.dr) / determinant;
        const double v = (xr * a.dz - xz * a.dr) / determinant;
        if (u >= 0 && u <= 1 && v >= 0 && v <= 1)
            return 0;
    }
    return std::min(
        {pointDistance(a.radius, a.axial, b), pointDistance(a.radius + a.dr, a.axial + a.dz, b),
         pointDistance(b.radius, b.axial, a), pointDistance(b.radius + b.dr, b.axial + b.dz, a)});
}
std::vector<RingPairPatch> pairRoots(const RingSegment &a, const RingSegment &b, double nearRatio) {
    if (a.first == b.first && a.second == b.second)
        return {{0, 1, 0, 1, 1, 0}, {0, 1, 0, 1, 2, 0}};
    if (a.first == b.first || a.first == b.second || a.second == b.first || a.second == b.second)
        return {{0, 1, 0, 1, 3, 0}, {0, 1, 0, 1, 4, 0}};
    const double distance = segmentDistance(a, b);
    if (!(distance > 0))
        throw std::invalid_argument("Nonadjacent meridian intervals intersect.");
    if (nearRatio == 0 || distance >= nearRatio * std::max(a.length, b.length))
        return {};
    std::vector<double> cuts{0, 1};
    const double initial =
        ((a.radius - b.radius) * b.dr + (a.axial - b.axial) * b.dz) / (b.length * b.length);
    const double slope = (a.dr * b.dr + a.dz * b.dz) / (b.length * b.length);
    if (slope != 0)
        for (double boundary : {0.0, 1.0}) {
            const double u = (boundary - initial) / slope;
            if (u > 0 && u < 1)
                cuts.push_back(u);
        }
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    std::vector<RingPairPatch> roots;
    for (std::size_t i = 1; i < cuts.size(); ++i)
        roots.push_back({cuts[i - 1], cuts[i], 0, 1, 5, 0});
    return roots;
}
RingPairMoments patchMoments(const RingSegment &a, const RingSegment &b, const RingPairPatch &patch,
                             const std::vector<RingRulePoint> &rule, double k,
                             double angularTolerance, unsigned &status, std::stop_token stop) {
    RingPairMoments result;
    for (const auto &u : rule) {
        cancelled(stop);
        for (const auto &v : rule) {
            status = accumulateRingPairPoint(result, mapRingPairPoint(a, b, patch, u, v), k,
                                             angularTolerance);
            if (status) {
                const auto point = mapRingPairPoint(a, b, patch, u, v);
                const auto diagnostic = ringKernels(point.geometry, k, angularTolerance);
                std::ostringstream message;
                message << std::setprecision(17) << "Ring pair kernel failure status=" << status
                        << ", mapping=" << patch.mapping << ", r=" << point.geometry.testRadius
                        << ", sourceR=" << point.geometry.sourceRadius
                        << ", dr=" << point.geometry.deltaRadius << ", dz=" << point.geometry.deltaZ
                        << ", nr=" << point.geometry.testRadialNormal
                        << ", nz=" << point.geometry.testAxialNormal
                        << ", sourceNr=" << point.geometry.sourceRadialNormal
                        << ", sourceNz=" << point.geometry.sourceAxialNormal
                        << ", nDistance=" << point.geometry.testNormalDistance
                        << ", sourceNDistance=" << point.geometry.sourceNormalDistance
                        << ", k=" << k << ", psiEvaluations=" << diagnostic.evaluations;
                throw std::runtime_error(message.str());
            }
        }
    }
    return result;
}
} // namespace

std::vector<RingRulePoint> ringRule(unsigned order) {
    const auto gauss = gaussLegendre01(order);
    std::vector<RingRulePoint> result;
    result.reserve(order);
    // Lagrange cardinal functions expanded in shifted Legendre polynomials.
    // int P_l(2x-1)log(x)dx=(-1)^(l+1)/(l(l+1)), l>0; l=0 gives-1.
    // Thus logarithmic weights integrate every polynomial of degree n-1 exactly.
    for (const auto &g : gauss) {
        const double x = 2 * g.position - 1;
        double first = -1, last = -1, p0 = 1, p1 = x;
        for (unsigned l = 1; l < order; ++l) {
            const double p = l == 1 ? p1 : ((2 * l - 1) * x * p1 - (l - 1) * p0) / l;
            if (l > 1) {
                p0 = p1;
                p1 = p;
            }
            const double moment = 1.0 / (l * (l + 1.0));
            first += (2 * l + 1) * p * (l % 2 ? moment : -moment);
            last -= (2 * l + 1) * p * moment;
        }
        result.push_back({g.position, g.weight, g.weight * first, g.weight * last});
    }
    return result;
}

RingInteractionPlan makeRingInteractionPlan(const AxisymmetricSpace &space,
                                            const BemNumerics &numerics, std::stop_token stop) {
    if (space.elements.empty() ||
        (numerics.quadratureOrder != 2 && numerics.quadratureOrder != 4) ||
        numerics.singularOrder < 2 || numerics.singularOrder > 20 ||
        !std::isfinite(numerics.nearDistanceRatio) || numerics.nearDistanceRatio < 0 ||
        !std::isfinite(numerics.nearRelativeTolerance) || numerics.nearRelativeTolerance <= 0 ||
        numerics.nearMaximumDepth > 16)
        throw std::invalid_argument("Invalid meridian quadrature settings.");
    RingInteractionPlan plan;
    plan.regularOrder = numerics.quadratureOrder;
    plan.coarseOrder = std::clamp(numerics.singularOrder, 4u, 16u);
    plan.fineOrder = plan.coarseOrder + 4;
    plan.regular = ringRule(plan.regularOrder);
    plan.coarse = ringRule(plan.coarseOrder);
    plan.fine = ringRule(plan.fineOrder);
    for (const auto &element : space.elements) {
        cancelled(stop);
        const auto a = space.meridian.points[element.first],
                   b = space.meridian.points[element.second];
        plan.segments.push_back({a.radius, a.axial, b.radius - a.radius, b.axial - a.axial,
                                 element.length, element.radialNormal, element.axialNormal,
                                 element.first, element.second});
    }
    const auto count = plan.segments.size();
    if (count > std::numeric_limits<std::uint32_t>::max() / count)
        throw std::length_error("Meridian pair index capacity exceeded.");
    plan.pairs.reserve(count * count);
    for (std::size_t i = 0; i < count; ++i)
        for (std::size_t j = 0; j < count; ++j) {
            cancelled(stop);
            const auto roots =
                pairRoots(plan.segments[i], plan.segments[j], numerics.nearDistanceRatio);
            if (plan.roots.size() + roots.size() > std::numeric_limits<std::uint32_t>::max())
                throw std::length_error("Meridian quadrature root capacity exceeded.");
            plan.pairs.push_back({static_cast<std::uint32_t>(plan.roots.size()),
                                  static_cast<std::uint32_t>(roots.size())});
            plan.roots.insert(plan.roots.end(), roots.begin(), roots.end());
        }
    return plan;
}

RingPairIntegral integrateRingPair(const RingSegment &a, const RingSegment &b,
                                   const std::vector<RingPairPatch> &roots,
                                   const std::vector<RingRulePoint> &coarse,
                                   const std::vector<RingRulePoint> &fine, double k,
                                   double tolerance, unsigned maximumDepth, std::stop_token stop) {
    RingPairIntegral output;
    if (roots.empty() || coarse.empty() || fine.empty() || !(tolerance > 0) || maximumDepth > 16) {
        output.status = 1;
        return output;
    }
    double angularTolerance = std::min(1e-8, tolerance * .005);
    RingPairMoments reference;
    double parameterArea = 0;
    for (const auto &root : roots) {
        const auto value = patchMoments(a, b, root, fine, k, angularTolerance, output.status, stop);
        if (output.status)
            return output;
        ringAccumulate(reference, value);
        parameterArea += (root.lastU - root.firstU) * (root.lastV - root.firstV);
    }
    const double referenceNorm = ringOperatorNorm(reference, a, b, k);
    if (!(referenceNorm > 0)) {
        output.status = 2;
        return output;
    }
    double localFactor = 1;
    for (unsigned attempt = 0; attempt < 3; ++attempt) {
        output.moments = {};
        output.patches = output.maximumDepth = 0;
        double absoluteError = 0, angularErrorSum = 0, meridianErrorSum = 0;
        auto pending = roots;
        unsigned work = 0;
        while (!pending.empty()) {
            cancelled(stop);
            const auto patch = pending.back();
            pending.pop_back();
            const auto low =
                patchMoments(a, b, patch, coarse, k, angularTolerance, output.status, stop);
            if (output.status)
                return output;
            const auto high =
                patchMoments(a, b, patch, fine, k, angularTolerance, output.status, stop);
            if (output.status)
                return output;
            const auto difference = ringSubtract(high, low);
            const double meridianError = ringOperatorNorm(difference, a, b, k);
            const double angularError = ringOperatorNorm(difference, a, b, k, true);
            const double error = meridianError + angularError;
            const double fraction =
                (patch.lastU - patch.firstU) * (patch.lastV - patch.firstV) / parameterArea;
            const double budget = tolerance * .1 * referenceNorm * fraction * localFactor;
            // Angular error cannot be repaired by meridian subdivision. Retain
            // it in the global sum and tighten its OWN rule only if necessary.
            // Likewise an at-limit local contribution can fit the global budget.
            if (meridianError <= budget || patch.depth == maximumDepth) {
                ringAccumulate(output.moments, high);
                absoluteError += error;
                angularErrorSum += angularError;
                meridianErrorSum += meridianError;
                ++output.patches;
                output.maximumDepth = std::max(output.maximumDepth, patch.depth);
            } else {
                if (++work > 50000) {
                    std::ostringstream message;
                    message << std::setprecision(17)
                            << "Meridian pair error limit: mapping=" << patch.mapping
                            << ", depth=" << patch.depth << ", u=[" << patch.firstU << ','
                            << patch.lastU << "], v=[" << patch.firstV << ',' << patch.lastV << "]"
                            << ", quadratureError=" << ringOperatorNorm(difference, a, b, k)
                            << ", angularError=" << ringOperatorNorm(difference, a, b, k, true)
                            << ", budget=" << budget << ", globalReference=" << referenceNorm;
                    throw std::runtime_error(message.str());
                }
                const double u = .5 * (patch.firstU + patch.lastU),
                             v = .5 * (patch.firstV + patch.lastV);
                for (unsigned i = 0; i < 4; ++i)
                    pending.push_back({i & 1 ? u : patch.firstU, i & 1 ? patch.lastU : u,
                                       i & 2 ? v : patch.firstV, i & 2 ? patch.lastV : v,
                                       patch.mapping, patch.depth + 1});
            }
        }
        output.estimatedRelativeError = absoluteError / ringOperatorNorm(output.moments, a, b, k);
        if (output.estimatedRelativeError <= tolerance * .2)
            return output;
        if (angularErrorSum > meridianErrorSum)
            angularTolerance *= .1;
        else
            localFactor *= .25;
    }
    output.status = 5;
    return output;
}
} // namespace hornsculpt::galerkin
