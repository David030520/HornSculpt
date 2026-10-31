#pragma once

#include "pair_quadrature.hpp"
#include <cmath>
#include <hornsculpt/mesh.hpp>

#if defined(__CUDACC__)
#define HS_NEAR_HD __host__ __device__
#else
#define HS_NEAR_HD
#endif

namespace hornsculpt::galerkin {
namespace projection_detail {
HS_NEAR_HD inline Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
HS_NEAR_HD inline double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
HS_NEAR_HD inline Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
// Gauss--Legendre rules on [0,1], shared by CPU planning and CUDA evaluation.
// The tests independently compare these literals with generated Gauss roots.
HS_NEAR_HD inline GaussPoint gauss(unsigned order, unsigned index) {
    const bool mirror = index >= order / 2;
    const unsigned i = mirror ? order - 1 - index : index;
    GaussPoint point{};
    if (order == 8) {
        switch (i) {
        case 0:
            point = {.019855071751231856, .05061426814518839};
            break;
        case 1:
            point = {.10166676129318658, .11119051722668723};
            break;
        case 2:
            point = {.2372337950418355, .1568533229389437};
            break;
        default:
            point = {.4082826787521751, .181341891689181};
            break;
        }
    } else if (order == 14) {
        switch (i) {
        case 0:
            point = {.006858095651593843, .01755973016587453};
            break;
        case 1:
            point = {.035782558168213185, .04007904357988019};
            break;
        case 2:
            point = {.08639934246511749, .06075928534395156};
            break;
        case 3:
            point = {.15635354759415726, .0786015835790968};
            break;
        case 4:
            point = {.24237568182092295, .09276919873896892};
            break;
        case 5:
            point = {.3404438155360551, .10259923186064775};
            break;
        default:
            point = {.44597252564632817, .10763192673157888};
            break;
        }
    } else {
        switch (i) {
        case 0:
            point = {.0034357004074525577, .008807003569575289};
            break;
        case 1:
            point = {.018014036361043095, .020300714900193525};
            break;
        case 2:
            point = {.04388278587433703, .03133602416705452};
            break;
        case 3:
            point = {.08044151408889055, .04163837078835237};
            break;
        case 4:
            point = {.1268340467699246, .05096505990862024};
            break;
        case 5:
            point = {.1819731596367425, .05909726598075916};
            break;
        case 6:
            point = {.24456649902458644, .0658443192245883};
            break;
        case 7:
            point = {.3131469556422902, .07104805465919095};
            break;
        case 8:
            point = {.38610707442917747, .07458649323630191};
            break;
        default:
            point = {.46173673943325133, .07637669356536299};
            break;
        }
    }
    if (mirror)
        point.position = 1 - point.position;
    return point;
}
} // namespace projection_detail

// Exact oriented fan decomposition of the trial triangle about the orthogonal
// projection of x. The projection may lie outside the triangle: signed fans
// then cancel, without clipping or changing either physical panel. Both sinh
// maps cluster according to geometric distance, never an artificial gap.
// Related near-integral transformations: Johnston & Elliott,
// doi:10.1002/nme.1208 and doi:10.1002/nme.1816; the formula below is derived
// directly from y=p+t*((1-u)*(a-p)+u*(b-p)) and its signed area Jacobian.
// These immutable states expose the exact same point map at its natural loop
// levels. CPU planning reuses them across inner samples; device evaluation can
// still map an arbitrary individual index without storing a point cloud.
struct ProjectionOuter {
    double test[3]{}, projection[3]{};
    Vec3 relative{};
    double height{}, h{}, weight{};
};
struct ProjectionEdge {
    Vec3 offset{}, direction{};
    double center{}, width{}, lower{}, upper{};
    unsigned edge{}, next{}, opposite{};
};
struct ProjectionRay {
    double boundary[3]{};
    double stretch{}, denominator{}, fanWeight{}, du{}, angularWeight{};
};

template <class PanelType>
HS_NEAR_HD inline ProjectionOuter
makeProjectionOuter(const TriangleQuadraturePoint &outer, const PairQuadraturePatch &patch,
                    const PanelType &test, const PanelType &trial) {
    using namespace projection_detail;
    ProjectionOuter point{};
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            point.test[i] += outer.barycentric[j] * patch.test[3 * j + i];
    const auto e1 = subtract(trial.vertices[1], trial.vertices[0]);
    const auto e2 = subtract(trial.vertices[2], trial.vertices[0]);
    const auto tx1 = subtract(test.vertices[1], test.vertices[0]);
    const auto tx2 = subtract(test.vertices[2], test.vertices[0]);
    point.relative = subtract(test.vertices[0], trial.vertices[0]);
    point.relative.x += tx1.x * point.test[1] + tx2.x * point.test[2];
    point.relative.y += tx1.y * point.test[1] + tx2.y * point.test[2];
    point.relative.z += tx1.z * point.test[1] + tx2.z * point.test[2];
    const auto g1 = cross(e2, trial.normal), g2 = cross(trial.normal, e1);
    point.projection[1] = dot(point.relative, g1) / (2 * trial.area);
    point.projection[2] = dot(point.relative, g2) / (2 * trial.area);
    point.projection[0] = 1 - point.projection[1] - point.projection[2];
    point.height = dot(point.relative, trial.normal);
    point.h = ::fabs(point.height);
    point.weight = outer.weight * patch.areaFraction;
    return point;
}

template <class PanelType>
HS_NEAR_HD inline ProjectionEdge makeProjectionEdge(const ProjectionOuter &outer, unsigned edge,
                                                    const PanelType &trial) {
    using namespace projection_detail;
    ProjectionEdge value{};
    value.edge = edge;
    value.next = (edge + 1) % 3;
    value.opposite = (edge + 2) % 3;
    const auto a = subtract(trial.vertices[edge], trial.vertices[0]);
    value.direction = subtract(trial.vertices[value.next], trial.vertices[edge]);
    value.offset = {a.x - outer.relative.x + outer.height * trial.normal.x,
                    a.y - outer.relative.y + outer.height * trial.normal.y,
                    a.z - outer.relative.z + outer.height * trial.normal.z};
    const double lengthSquared = dot(value.direction, value.direction);
    value.center = -dot(value.offset, value.direction) / lengthSquared;
    const Vec3 perpendicular{value.offset.x + value.center * value.direction.x,
                             value.offset.y + value.center * value.direction.y,
                             value.offset.z + value.center * value.direction.z};
    value.width = ::sqrt((dot(perpendicular, perpendicular) + outer.h * outer.h) / lengthSquared);
    value.lower = ::asinh(-value.center / value.width);
    value.upper = ::asinh((1 - value.center) / value.width);
    return value;
}

HS_NEAR_HD inline ProjectionRay makeProjectionRay(const ProjectionOuter &outer,
                                                  const ProjectionEdge &edge, GaussPoint angular) {
    using namespace projection_detail;
    ProjectionRay value{};
    const double parameter = edge.lower + (edge.upper - edge.lower) * angular.position;
    const double u = edge.center + edge.width * ::sinh(parameter);
    value.du = edge.width * (edge.upper - edge.lower) * ::cosh(parameter);
    const Vec3 ray{edge.offset.x + u * edge.direction.x, edge.offset.y + u * edge.direction.y,
                   edge.offset.z + u * edge.direction.z};
    value.stretch = ::asinh(::sqrt(dot(ray, ray)) / outer.h);
    value.denominator = ::sinh(value.stretch);
    for (unsigned i = 0; i < 3; ++i)
        value.boundary[i] = i == edge.edge ? 1 - u : (i == edge.next ? u : 0);
    value.fanWeight = outer.weight * outer.projection[edge.opposite];
    value.angularWeight = angular.weight;
    return value;
}

HS_NEAR_HD inline PairQuadraturePoint
mapProjectionRay(const ProjectionOuter &outer, const ProjectionRay &ray, GaussPoint radial) {
    PairQuadraturePoint point{};
    const double stretch = ray.stretch, denominator = ray.denominator;
    const double t =
        stretch == 0 ? radial.position : ::sinh(stretch * radial.position) / denominator;
    const double dt = stretch == 0 ? 1 : stretch * ::cosh(stretch * radial.position) / denominator;
    for (unsigned i = 0; i < 3; ++i) {
        point.test[i] = outer.test[i];
        point.trial[i] = outer.projection[i] + t * (ray.boundary[i] - outer.projection[i]);
    }
    point.weight = ray.fanWeight * t * dt * ray.du * ray.angularWeight * radial.weight;
    return point;
}

template <class PanelType>
HS_NEAR_HD inline PairQuadraturePoint
mapProjectionPatch(const TriangleQuadraturePoint &outer, unsigned edge, GaussPoint angular,
                   GaussPoint radial, const PairQuadraturePatch &patch, const PanelType &test,
                   const PanelType &trial) {
    const auto x = makeProjectionOuter(outer, patch, test, trial);
    return mapProjectionRay(x, makeProjectionRay(x, makeProjectionEdge(x, edge, trial), angular),
                            radial);
}

HS_NEAR_HD inline unsigned nearPatchPointCount(const PairQuadraturePatch &patch,
                                               unsigned outerCount) {
    const unsigned order = 20 * patch.projectionIntervals;
    return patch.projection ? outerCount * 3 * order * order : outerCount * outerCount;
}

template <class PanelType>
HS_NEAR_HD inline PairQuadraturePoint
mapNearPatchPoint(const TriangleQuadraturePoint *rule, unsigned count,
                  const PairQuadraturePatch &patch, unsigned index, const PanelType &test,
                  const PanelType &trial) {
    if (!patch.projection)
        return mapPairPatch(rule[index / count], rule[index % count], patch);
    const unsigned order = 20 * patch.projectionIntervals;
    const unsigned radial = index % order;
    index /= order;
    const unsigned angular = index % order;
    index /= order;
    const unsigned edge = index % 3, outer = index / 3;
    auto u = projection_detail::gauss(20, angular % 20);
    auto v = projection_detail::gauss(20, radial % 20);
    u.position = (angular / 20 + u.position) / patch.projectionIntervals;
    v.position = (radial / 20 + v.position) / patch.projectionIntervals;
    u.weight /= patch.projectionIntervals;
    v.weight /= patch.projectionIntervals;
    return mapProjectionPatch(rule[outer], edge, u, v, patch, test, trial);
}
} // namespace hornsculpt::galerkin

#undef HS_NEAR_HD
