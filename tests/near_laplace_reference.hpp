#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <galerkin/pair_quadrature.hpp>
#include <galerkin/space.hpp>
#include <vector>

namespace hornsculpt::testing {
namespace laplace_reference {
inline Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline Vec3 sub(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline Vec3 mul(Vec3 a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}
inline double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(Vec3 a) {
    return std::hypot(a.x, a.y, a.z);
}
struct Inner {
    std::array<double, 3> potential{};
    std::array<Vec3, 3> gradient{};
};

// Independent semi-analytic reference: integrate 1/r and its gradient over
// the entire inner triangle exactly using the planar divergence theorem.
// Only the outer triangle is uniformly integrated; neither the production
// projection fan, sinh variables nor its adaptive estimate are used here.
inline Inner inner(Vec3 x, const galerkin::Panel &trial) {
    const auto relative = sub(x, trial.vertices[0]);
    const double h = dot(relative, trial.normal);
    std::array<Vec3, 3> radius{};
    for (unsigned j = 0; j < 3; ++j)
        radius[j] = sub(x, trial.vertices[j]);
    const double r0 = norm(radius[0]), r1 = norm(radius[1]), r2 = norm(radius[2]);
    const double angle =
        2 * std::atan2(dot(radius[0], cross(radius[1], radius[2])),
                       r0 * r1 * r2 + dot(radius[0], radius[1]) * r2 +
                           dot(radius[1], radius[2]) * r0 + dot(radius[2], radius[0]) * r1);
    double scalar{};
    Vec3 firstMoment{}, boundary{};
    std::array<Vec3, 3> weightedBoundary{};
    for (unsigned edge = 0; edge < 3; ++edge) {
        const auto vector = sub(trial.vertices[(edge + 1) % 3], trial.vertices[edge]);
        const double length = norm(vector);
        const auto tangent = mul(vector, 1 / length), outward = cross(tangent, trial.normal);
        const auto a = add(sub(trial.vertices[edge], x), mul(trial.normal, h));
        const double distance = dot(a, outward), ta = dot(a, tangent), tb = ta + length;
        const double width = std::hypot(distance, h), ra = std::hypot(ta, width),
                     rb = std::hypot(tb, width);
        const double logarithm = std::asinh(tb / width) - std::asinh(ta / width);
        const double primitive = .5 * (tb * rb - ta * ra + width * width * logarithm);
        const auto edgeMoment = add(mul(tangent, rb - ra), mul(outward, distance * logarithm));
        scalar += distance * logarithm;
        firstMoment = add(firstMoment, mul(outward, primitive));
        boundary = add(boundary, mul(outward, logarithm));
        for (unsigned j = 0; j < 3; ++j)
            weightedBoundary[j] =
                add(weightedBoundary[j], mul(outward, dot(trial.gradients[j], edgeMoment)));
    }
    scalar -= h * angle;
    Inner result;
    for (unsigned j = 0; j < 3; ++j) {
        const double phi = (j == 0 ? 1.0 : 0.0) + dot(trial.gradients[j], relative);
        result.potential[j] = phi * scalar + dot(trial.gradients[j], firstMoment);
        result.gradient[j] =
            sub(sub(mul(trial.gradients[j], scalar), mul(boundary, phi)), weightedBoundary[j]);
        result.gradient[j] =
            sub(result.gradient[j],
                mul(trial.normal, phi * angle - h * dot(trial.gradients[j], boundary)));
    }
    return result;
}
} // namespace laplace_reference

inline std::array<double, 36> semiAnalyticLaplacePair(const galerkin::Panel &test,
                                                      const galerkin::Panel &trial,
                                                      unsigned subdivisions, unsigned order = 16) {
    using namespace laplace_reference;
    using Bary = std::array<double, 3>;
    using Triangle = std::array<Bary, 3>;
    std::vector<Triangle> leaves{{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}};
    double fraction = 1;
    const auto midpoint = [](Bary a, Bary b) {
        return Bary{.5 * (a[0] + b[0]), .5 * (a[1] + b[1]), .5 * (a[2] + b[2])};
    };
    for (unsigned depth = 0; depth < subdivisions; ++depth) {
        std::vector<Triangle> children;
        children.reserve(4 * leaves.size());
        for (const auto &t : leaves) {
            const auto ab = midpoint(t[0], t[1]), bc = midpoint(t[1], t[2]),
                       ca = midpoint(t[2], t[0]);
            children.push_back({t[0], ab, ca});
            children.push_back({ab, t[1], bc});
            children.push_back({ca, bc, t[2]});
            children.push_back({ab, bc, ca});
        }
        leaves = std::move(children);
        fraction *= .25;
    }
    const auto gauss = galerkin::gaussLegendre01(order);
    const double scale = std::max(test.maximumEdgeLength, trial.maximumEdgeLength);
    std::array<double, 36> result{};
    for (const auto &leaf : leaves)
        for (const auto &u : gauss)
            for (const auto &v : gauss) {
                const double local[3]{1 - u.position, u.position * (1 - v.position),
                                      u.position * v.position};
                Bary phi{};
                for (unsigned i = 0; i < 3; ++i)
                    for (unsigned j = 0; j < 3; ++j)
                        phi[i] += local[j] * leaf[j][i];
                const auto x = add(test.vertices[0],
                                   add(mul(sub(test.vertices[1], test.vertices[0]), phi[1]),
                                       mul(sub(test.vertices[2], test.vertices[0]), phi[2])));
                const auto source = inner(x, trial);
                const double weight =
                    fraction * u.weight * v.weight * u.position / (2 * trial.area);
                for (unsigned i = 0; i < 3; ++i)
                    for (unsigned j = 0; j < 3; ++j) {
                        const unsigned index = 3 * i + j;
                        result[index] += weight * phi[i] * source.potential[j];
                        result[9 + index] += weight * phi[i] * scale * source.gradient[j].x;
                        result[18 + index] += weight * phi[i] * scale * source.gradient[j].y;
                        result[27 + index] += weight * phi[i] * scale * source.gradient[j].z;
                    }
            }
    return result;
}
} // namespace hornsculpt::testing
