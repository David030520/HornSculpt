#pragma once

#include <galerkin/near_projection.hpp>
#include <galerkin/space.hpp>

namespace hornsculpt::testing {
// Unfactored point formula retained ONLY as a numerical reference. It checks
// that caching the outer/edge/ray geometry moves neither quadrature nodes nor
// signed weights. The separate analytic Laplace reference checks accuracy.
inline galerkin::PairQuadraturePoint
directProjectionPoint(const galerkin::TriangleQuadraturePoint &outer, unsigned edge,
                      galerkin::GaussPoint angular, galerkin::GaussPoint radial,
                      const galerkin::PairQuadraturePatch &patch, const galerkin::Panel &test,
                      const galerkin::Panel &trial) {
    using namespace galerkin::projection_detail;
    galerkin::PairQuadraturePoint point{};
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            point.test[i] += outer.barycentric[j] * patch.test[3 * j + i];
    const auto e1 = subtract(trial.vertices[1], trial.vertices[0]);
    const auto e2 = subtract(trial.vertices[2], trial.vertices[0]);
    const auto tx1 = subtract(test.vertices[1], test.vertices[0]);
    const auto tx2 = subtract(test.vertices[2], test.vertices[0]);
    auto relative = subtract(test.vertices[0], trial.vertices[0]);
    relative.x += tx1.x * point.test[1] + tx2.x * point.test[2];
    relative.y += tx1.y * point.test[1] + tx2.y * point.test[2];
    relative.z += tx1.z * point.test[1] + tx2.z * point.test[2];
    const auto g1 = cross(e2, trial.normal), g2 = cross(trial.normal, e1);
    double projection[3]{};
    projection[1] = dot(relative, g1) / (2 * trial.area);
    projection[2] = dot(relative, g2) / (2 * trial.area);
    projection[0] = 1 - projection[1] - projection[2];
    const double height = dot(relative, trial.normal), h = ::fabs(height);
    const unsigned next = (edge + 1) % 3, opposite = (edge + 2) % 3;
    const auto a = subtract(trial.vertices[edge], trial.vertices[0]);
    const auto direction = subtract(trial.vertices[next], trial.vertices[edge]);
    const Vec3 offset{a.x - relative.x + height * trial.normal.x,
                      a.y - relative.y + height * trial.normal.y,
                      a.z - relative.z + height * trial.normal.z};
    const double lengthSquared = dot(direction, direction);
    const double center = -dot(offset, direction) / lengthSquared;
    const Vec3 perpendicular{offset.x + center * direction.x, offset.y + center * direction.y,
                             offset.z + center * direction.z};
    const double width = ::sqrt((dot(perpendicular, perpendicular) + h * h) / lengthSquared);
    const double lower = ::asinh(-center / width), upper = ::asinh((1 - center) / width);
    const double parameter = lower + (upper - lower) * angular.position;
    const double u = center + width * ::sinh(parameter);
    const double du = width * (upper - lower) * ::cosh(parameter);
    const Vec3 ray{offset.x + u * direction.x, offset.y + u * direction.y,
                   offset.z + u * direction.z};
    const double stretch = ::asinh(::sqrt(dot(ray, ray)) / h);
    const double denominator = ::sinh(stretch);
    const double t =
        stretch == 0 ? radial.position : ::sinh(stretch * radial.position) / denominator;
    const double dt = stretch == 0 ? 1 : stretch * ::cosh(stretch * radial.position) / denominator;
    for (unsigned i = 0; i < 3; ++i) {
        const double boundary = i == edge ? 1 - u : (i == next ? u : 0);
        point.trial[i] = projection[i] + t * (boundary - projection[i]);
    }
    point.weight = outer.weight * patch.areaFraction * projection[opposite] * t * dt * du *
                   angular.weight * radial.weight;
    return point;
}
} // namespace hornsculpt::testing
