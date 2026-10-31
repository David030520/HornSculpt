#include "axisymmetric_field.hpp"
#include "axisymmetric_reduction.hpp"
#include "galerkin/axisymmetric_quadrature.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>

namespace hornsculpt::bem_detail {
namespace {
using Complex = std::complex<double>;
Vec3 sub(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 mul(Vec3 a, double b) {
    return {a.x * b, a.y * b, a.z * b};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
double norm(Vec3 a) {
    return std::hypot(a.x, a.y, a.z);
}
bool finite(Vec3 p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
Complex complex(galerkin::RingComplex p) {
    return {p.real, p.imag};
}
struct Point {
    double radius, axial;
};
double boundaryDistance(Point p, const galerkin::RingSegment &s) {
    const double t = std::clamp(((p.radius - s.radius) * s.dr + (p.axial - s.axial) * s.dz) /
                                    (s.length * s.length),
                                0., 1.);
    return std::hypot(p.radius - s.radius - t * s.dr, p.axial - s.axial - t * s.dz);
}
bool inside(Point p, const std::vector<galerkin::RingSegment> &segments) {
    int winding = 0;
    for (const auto &s : segments) {
        const double lastZ = s.axial + s.dz;
        const double side = s.dr * (p.axial - s.axial) - s.dz * (p.radius - s.radius);
        if (s.axial <= p.axial) {
            if (lastZ > p.axial && side > 0)
                ++winding;
        } else if (lastZ <= p.axial && side < 0)
            --winding;
    }
    // An open pole-to-pole generator is closed on r=0 for winding ONLY.
    // Its axis edge cannot cross the rightward ray for r>=0, including r=0;
    // it is deliberately absent from physical boundary-distance calculations.
    return winding != 0;
}
struct Integral {
    Complex value{};
    double magnitude{}, angularError{};
    unsigned status{};
};
Integral integrate(const galerkin::RingSegment &s, Point p, double first, double last,
                   const std::vector<galerkin::RingRulePoint> &rule, double k, Complex p0,
                   Complex p1, Complex q, bool rayleigh, bool sinhMap, double angularTolerance,
                   std::stop_token stop) {
    Integral result;
    for (const auto &node : rule) {
        if (stop.stop_requested())
            throw std::runtime_error("Meridian field evaluation cancelled.");
        const double u = first + (last - first) * node.position;
        auto mapped = galerkin::ring_pair_detail::SinhPoint{u, 1, 0};
        if (sinhMap)
            mapped = galerkin::ring_pair_detail::projectedSinh(
                p.radius - s.radius, p.axial - s.axial, s.dr, s.dz, s.length, u);
        if (mapped.status) {
            result.status = mapped.status;
            return result;
        }
        const double radius = s.radius + mapped.position * s.dr;
        galerkin::RingGeometry geometry{p.radius,
                                        radius,
                                        p.radius - radius,
                                        p.axial - s.axial - mapped.position * s.dz,
                                        0,
                                        0,
                                        s.nr,
                                        s.nz,
                                        0,
                                        s.nr * (s.radius - p.radius) + s.nz * (s.axial - p.axial)};
        const auto kernel = galerkin::ringKernels(geometry, k, angularTolerance);
        if (kernel.status) {
            result.status = kernel.status;
            return result;
        }
        const auto pressure = (1 - mapped.position) * p0 + mapped.position * p1;
        const auto single = complex(kernel.value[0]) * q;
        const auto dual = complex(kernel.value[3]) * pressure;
        const double weight = (last - first) * node.weight * mapped.jacobian * s.length * radius;
        result.value += weight * (rayleigh ? -2. * single : single - dual);
        result.magnitude +=
            weight * (rayleigh ? 2 * std::abs(single) : std::abs(single) + std::abs(dual));
        result.angularError +=
            weight * (rayleigh ? 2 * std::abs(q) * kernel.angularAbsoluteError[0]
                               : std::abs(q) * kernel.angularAbsoluteError[0] +
                                     std::abs(pressure) * kernel.angularAbsoluteError[3]);
    }
    return result;
}
Integral adaptive(const galerkin::RingSegment &s, Point p, double k, Complex p0, Complex p1,
                  Complex q, bool rayleigh, const std::vector<galerkin::RingRulePoint> &coarse,
                  const std::vector<galerkin::RingRulePoint> &fine, PressureFieldSettings settings,
                  std::stop_token stop) {
    const bool sinhMap = boundaryDistance(p, s) < s.length;
    double angularTolerance = std::min(1e-8, settings.relativeQuadratureTolerance * .005),
           localFactor = 1;
    const auto reference =
        integrate(s, p, 0, 1, fine, k, p0, p1, q, rayleigh, sinhMap, angularTolerance, stop);
    if (reference.status)
        return reference;
    if (reference.magnitude == 0)
        return reference;
    for (unsigned attempt = 0; attempt < 3; ++attempt) {
        struct Patch {
            double first, last;
            unsigned depth;
        };
        std::vector<Patch> pending{{0, 1, 0}};
        Integral result;
        double spatialError = 0;
        unsigned work = 0;
        while (!pending.empty()) {
            const auto patch = pending.back();
            pending.pop_back();
            const auto low = integrate(s, p, patch.first, patch.last, coarse, k, p0, p1, q,
                                       rayleigh, sinhMap, angularTolerance, stop);
            const auto high = integrate(s, p, patch.first, patch.last, fine, k, p0, p1, q, rayleigh,
                                        sinhMap, angularTolerance, stop);
            if (low.status || high.status) {
                result.status = low.status ? low.status : high.status;
                return result;
            }
            const double error = std::abs(high.value - low.value);
            const double budget = settings.relativeQuadratureTolerance * .1 * reference.magnitude *
                                  (patch.last - patch.first) * localFactor;
            if (error <= budget || patch.depth == settings.maximumQuadratureDepth) {
                result.value += high.value;
                result.magnitude += high.magnitude;
                result.angularError += high.angularError + low.angularError;
                spatialError += error;
            } else {
                if (++work > 50000) {
                    result.status = 4;
                    return result;
                }
                const double mid = .5 * (patch.first + patch.last);
                pending.push_back({mid, patch.last, patch.depth + 1});
                pending.push_back({patch.first, mid, patch.depth + 1});
            }
        }
        if (spatialError + result.angularError <=
            settings.relativeQuadratureTolerance * .2 * result.magnitude)
            return result;
        if (result.angularError > spatialError)
            angularTolerance *= .1;
        else
            localFactor *= .25;
    }
    Integral failed;
    failed.status = 5;
    return failed;
}
} // namespace

std::vector<PressureFieldSample>
evaluateAxisymmetricPressureField(const SurfaceMesh &mesh, const BemRequest &request,
                                  double frequency, const BemExcitationResult &solution,
                                  const std::vector<Vec3> &points, PressureFieldSettings settings,
                                  std::stop_token stop,
                                  std::function<void(std::size_t, std::size_t)> progress) {
    if (!request.axisymmetric || !request.infiniteBaffle ||
        request.axisymmetric->surfaceContract != AxisymmetricSurfaceContract::continuousRevolution)
        throw std::invalid_argument(
            "Continuous meridian field requires its recorded surface contract and baffle.");
    validateAxisymmetricTrace(*request.axisymmetric, solution.boundary()->pressureTrace,
                              solution.boundary()->normalDerivative);
    const auto space = galerkin::buildAxisymmetricSpace(mesh, {}, *request.axisymmetric);
    // Only segment/rule construction is needed; avoid constructing an O(M^2)
    // interaction plan for a post-processed field raster.
    std::vector<galerkin::RingSegment> segments;
    for (const auto &e : space.elements) {
        const auto a = space.meridian.points[e.first], b = space.meridian.points[e.second];
        segments.push_back({a.radius, a.axial, b.radius - a.radius, b.axial - a.axial, e.length,
                            e.radialNormal, e.axialNormal, e.first, e.second});
    }
    std::vector<Complex> pressure(space.meridian.points.size()), derivative(space.elements.size());
    for (std::size_t v = 0; v < mesh.vertices.size(); ++v)
        pressure[request.axisymmetric->pressureDofs[v]] = solution.boundary()->pressureTrace[v];
    for (std::size_t t = 0; t < mesh.triangles.size(); ++t)
        derivative[request.axisymmetric->fluxDofs[t]] = solution.boundary()->normalDerivative[t];
    std::set<std::uint32_t> aperture(request.infiniteBaffle->apertureFaceIds.begin(),
                                     request.infiniteBaffle->apertureFaceIds.end());
    const auto normal = request.infiniteBaffle->normal;
    if (!finite(normal) || std::abs(norm(normal) - 1) > 1e-12 ||
        !finite(request.infiniteBaffle->origin) || aperture.empty())
        throw std::invalid_argument("Invalid continuous meridian field baffle.");
    if (dot(normal, space.meridian.axis) < 1 - 1e-12)
        throw std::invalid_argument(
            "Continuous field baffle normal must follow the rotation axis.");
    const auto coarse = galerkin::ringRule(8), fine = galerkin::ringRule(16);
    const double k = 2 * std::numbers::pi * frequency / request.medium.soundSpeed;
    std::vector<PressureFieldSample> output;
    output.reserve(points.size());
    for (const auto world : points) {
        if (stop.stop_requested())
            throw std::runtime_error("Meridian field evaluation cancelled.");
        if (!finite(world))
            throw std::invalid_argument("Non-finite meridian field point.");
        const auto relative = sub(world, space.meridian.origin);
        const double z = dot(relative, space.meridian.axis);
        const Point point{norm(sub(relative, mul(space.meridian.axis, z))), z};
        const double height = dot(sub(world, request.infiniteBaffle->origin), normal);
        PressureFieldSample sample{world, {}, FieldPointStatus::excluded};
        bool boundary = std::abs(height) <= settings.minimumBoundaryDistance;
        if (!boundary)
            for (const auto &s : segments)
                if (boundaryDistance(point, s) <= settings.minimumBoundaryDistance) {
                    boundary = true;
                    break;
                }
        const bool rayleigh = height > settings.minimumBoundaryDistance;
        if (boundary)
            sample.status = FieldPointStatus::boundary;
        else if (rayleigh || inside(point, segments)) {
            sample.status = FieldPointStatus::air;
            for (std::size_t t = 0; t < segments.size(); ++t) {
                if (rayleigh && !aperture.contains(space.elements[t].faceId))
                    continue;
                const auto &s = segments[t];
                if (derivative[t] == Complex{} && (rayleigh || (pressure[s.first] == Complex{} &&
                                                                pressure[s.second] == Complex{})))
                    continue;
                const auto integral =
                    adaptive(s, point, k, pressure[s.first], pressure[s.second], derivative[t],
                             rayleigh, coarse, fine, settings, stop);
                if (integral.status) {
                    sample.status = FieldPointStatus::quadratureNotConverged;
                    sample.pressure = {};
                    break;
                }
                sample.pressure += integral.value;
            }
            if (!std::isfinite(sample.pressure.real()) || !std::isfinite(sample.pressure.imag())) {
                sample.status = FieldPointStatus::quadratureNotConverged;
                sample.pressure = {};
            }
        }
        output.push_back(sample);
        if (progress)
            progress(output.size(), points.size());
    }
    return output;
}
} // namespace hornsculpt::bem_detail
