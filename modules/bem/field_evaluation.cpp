#include "axisymmetric_field.hpp"
#include "galerkin/pair_quadrature.hpp"
#include "galerkin/space.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <hornsculpt/field_evaluation.hpp>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>

namespace hornsculpt {
namespace {
using Complex = std::complex<double>;
Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3 sub(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 mul(Vec3 a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double norm(Vec3 a) {
    return std::hypot(a.x, a.y, a.z);
}
bool finite(Vec3 a) {
    return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z);
}
bool finite(Complex a) {
    return std::isfinite(a.real()) && std::isfinite(a.imag());
}

// Closest-point regions of a triangle (Ericson, Real-Time Collision Detection,
// section 5.1.5). The mesh contract excludes degenerate triangles.
double distanceToPanel(Vec3 p, const galerkin::Panel &t) {
    const auto a = t.vertices[0], ab = sub(t.vertices[1], a), ac = sub(t.vertices[2], a);
    const auto ap = sub(p, a);
    const double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0)
        return norm(ap);
    const auto bp = sub(p, t.vertices[1]);
    const double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3)
        return norm(bp);
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0)
        return norm(sub(ap, mul(ab, d1 / (d1 - d3))));
    const auto cp = sub(p, t.vertices[2]);
    const double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6)
        return norm(cp);
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0)
        return norm(sub(ap, mul(ac, d2 / (d2 - d6))));
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0 && d4 >= d3 && d5 >= d6) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return norm(sub(bp, mul(sub(t.vertices[2], t.vertices[1]), w)));
    }
    return std::abs(dot(ap, t.normal));
}
double solidAngle(Vec3 p, const galerkin::Panel &t) {
    const auto a = sub(t.vertices[0], p), b = sub(t.vertices[1], p), c = sub(t.vertices[2], p);
    const double la = norm(a), lb = norm(b), lc = norm(c);
    return 2 * std::atan2(dot(a, cross(b, c)),
                          la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb);
}
struct Patch {
    std::array<Vec3, 3> bary{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    double areaFraction{1};
};
struct Integral {
    Complex value{};
    double magnitude{};
};
Integral integrate(const galerkin::Panel &panel, const Patch &patch, Vec3 observer, double k,
                   const std::array<Complex, 3> &pressure, const std::array<Complex, 3> &derivative,
                   bool rayleigh, const std::vector<galerkin::TriangleQuadraturePoint> &rule) {
    Integral result;
    for (const auto &point : rule) {
        Vec3 b{};
        for (unsigned i = 0; i < 3; ++i)
            b = add(b, mul(patch.bary[i], point.barycentric[i]));
        const auto y = add(add(mul(panel.vertices[0], b.x), mul(panel.vertices[1], b.y)),
                           mul(panel.vertices[2], b.z));
        const auto displacement = sub(observer, y);
        const double r = norm(displacement);
        if (!(r > 0))
            throw std::runtime_error("Field observation lies on a quadrature point.");
        const Complex g = std::polar(1 / (4 * std::numbers::pi * r), -k * r);
        const Complex single =
            g * (b.x * derivative[0] + b.y * derivative[1] + b.z * derivative[2]);
        const Complex doubleLayer = g * Complex{1 / r, k} * (dot(displacement, panel.normal) / r) *
                                    (b.x * pressure[0] + b.y * pressure[1] + b.z * pressure[2]);
        // triangleRule weights integrate the reference triangle (area 1/2).
        const double weight = 2 * panel.area * patch.areaFraction * point.weight;
        result.value += weight * (rayleigh ? -2. * single : doubleLayer - single);
        result.magnitude +=
            weight * (rayleigh ? 2 * std::abs(single) : std::abs(doubleLayer) + std::abs(single));
    }
    return result;
}
std::array<Patch, 4> children(const Patch &p) {
    const auto a = p.bary[0], b = p.bary[1], c = p.bary[2];
    const auto ab = mul(add(a, b), .5), bc = mul(add(b, c), .5), ca = mul(add(c, a), .5);
    const double area = .25 * p.areaFraction;
    return {{{{a, ab, ca}, area}, {{ab, b, bc}, area}, {{ca, bc, c}, area}, {{ab, bc, ca}, area}}};
}
bool adaptive(const galerkin::Panel &panel, const Patch &patch, Vec3 observer, double k,
              const std::array<Complex, 3> &p, const std::array<Complex, 3> &q, bool rayleigh,
              const std::vector<galerkin::TriangleQuadraturePoint> &coarseRule,
              const std::vector<galerkin::TriangleQuadraturePoint> &fineRule,
              const PressureFieldSettings &settings, unsigned depth, Complex &value,
              std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Field evaluation cancelled.");
    const auto coarse = integrate(panel, patch, observer, k, p, q, rayleigh, coarseRule);
    const auto fine = integrate(panel, patch, observer, k, p, q, rayleigh, fineRule);
    const double phaseSpan = k * panel.maximumEdgeLength * std::sqrt(patch.areaFraction);
    if (phaseSpan <= std::numbers::pi &&
        std::abs(fine.value - coarse.value) <=
            settings.relativeQuadratureTolerance * std::max(fine.magnitude, 1e-300)) {
        value += fine.value;
        return true;
    }
    if (depth >= settings.maximumQuadratureDepth)
        return false;
    for (const auto &child : children(patch))
        if (!adaptive(panel, child, observer, k, p, q, rayleigh, coarseRule, fineRule, settings,
                      depth + 1, value, stop))
            return false;
    return true;
}
bool adaptive(const galerkin::Panel &panel, const Patch &patch, Vec3 observer, double k,
              const std::array<Complex, 3> &p, Complex q, bool rayleigh,
              const std::vector<galerkin::TriangleQuadraturePoint> &coarseRule,
              const std::vector<galerkin::TriangleQuadraturePoint> &fineRule,
              const PressureFieldSettings &settings, unsigned depth, Complex &value,
              std::stop_token stop) {
    return adaptive(panel, patch, observer, k, p, std::array<Complex, 3>{q, q, q}, rayleigh,
                    coarseRule, fineRule, settings, depth, value, stop);
}
struct FieldDomain {
    BemDomainKind kind;
    galerkin::Space space;
    std::vector<std::array<Complex, 3>> pressure;
    std::vector<std::array<Complex, 3>> flux;
    std::vector<bool> interfacePanels;
};
std::vector<PressureFieldSample>
evaluateSubdomainField(const SurfaceMesh &mesh, const BemRequest &request, double frequency,
                       const BemExcitationResult &solution, const std::vector<Vec3> &points,
                       const PressureFieldSettings &settings, std::stop_token stop,
                       const std::function<void(std::size_t, std::size_t)> &progress) {
    std::vector<FieldDomain> domains;
    for (const auto &definition : request.subdomains->domains) {
        const auto found =
            std::find_if(solution.subdomainTraces().begin(), solution.subdomainTraces().end(),
                         [&](const auto &x) { return x.domainId == definition.id; });
        if (found == solution.subdomainTraces().end() ||
            found->vertexIndices.size() != found->pressureTrace.size() ||
            found->triangleIndices.size() != found->normalDerivative.size())
            throw std::invalid_argument("Incomplete subdomain field traces.");
        std::map<std::uint32_t, Complex> pressures, means;
        for (std::size_t i = 0; i < found->vertexIndices.size(); ++i)
            if (!finite(found->pressureTrace[i]) ||
                !pressures.emplace(found->vertexIndices[i], found->pressureTrace[i]).second)
                throw std::invalid_argument("Invalid subdomain pressure trace.");
        for (std::size_t i = 0; i < found->triangleIndices.size(); ++i)
            if (!finite(found->normalDerivative[i]) ||
                !means.emplace(found->triangleIndices[i], found->normalDerivative[i]).second)
                throw std::invalid_argument("Invalid subdomain flux trace.");
        std::map<std::uint32_t, std::map<std::uint32_t, Complex>> interfaceFlux;
        for (const auto &trace : found->interfaceFluxes) {
            if (trace.vertexIndices.size() != trace.normalDerivative.size() ||
                interfaceFlux.contains(trace.interfaceId))
                throw std::invalid_argument("Invalid P1 interface trace.");
            auto &entries = interfaceFlux[trace.interfaceId];
            for (std::size_t i = 0; i < trace.vertexIndices.size(); ++i)
                if (!finite(trace.normalDerivative[i]) ||
                    !entries.emplace(trace.vertexIndices[i], trace.normalDerivative[i]).second)
                    throw std::invalid_argument("Invalid P1 interface flux.");
        }
        std::map<std::uint32_t, std::uint32_t> interfaceOf;
        for (const auto &interface : request.subdomains->interfaces)
            if (interface.firstDomainId == definition.id ||
                interface.secondDomainId == definition.id)
                for (auto t : interface.triangleIndices)
                    if (!interfaceOf.emplace(t, interface.id).second)
                        throw std::invalid_argument("Overlapping field interfaces.");
        SurfaceMesh local;
        std::map<std::uint32_t, std::uint32_t> ids;
        FieldDomain output;
        output.kind = definition.kind;
        for (const auto &entry : definition.panels) {
            if (entry.triangleIndex >= mesh.triangles.size() ||
                (entry.normalSign != 1 && entry.normalSign != -1) ||
                !means.contains(entry.triangleIndex))
                throw std::invalid_argument("Invalid subdomain field panel.");
            auto triangle = mesh.triangles[entry.triangleIndex];
            if (entry.normalSign < 0)
                std::swap(triangle.vertices[1], triangle.vertices[2]);
            std::array<Complex, 3> p{}, q{};
            for (unsigned j = 0; j < 3; ++j) {
                const auto global = triangle.vertices[j];
                if (global >= mesh.vertices.size() || !pressures.contains(global))
                    throw std::invalid_argument("Missing field pressure node.");
                p[j] = pressures.at(global);
                q[j] = means.at(entry.triangleIndex);
                if (interfaceOf.contains(entry.triangleIndex)) {
                    const auto face = interfaceOf.at(entry.triangleIndex);
                    if (!interfaceFlux.contains(face) || !interfaceFlux.at(face).contains(global))
                        throw std::invalid_argument("Missing field interface flux node.");
                    q[j] = interfaceFlux.at(face).at(global);
                }
                const auto [it, added] =
                    ids.emplace(global, static_cast<std::uint32_t>(ids.size()));
                if (added)
                    local.vertices.push_back(mesh.vertices[global]);
                triangle.vertices[j] = it->second;
            }
            local.triangles.push_back(triangle);
            output.pressure.push_back(p);
            output.flux.push_back(q);
            output.interfacePanels.push_back(interfaceOf.contains(entry.triangleIndex));
        }
        output.space = galerkin::buildSpace(local);
        domains.push_back(std::move(output));
    }
    const auto coarse = galerkin::triangleRule(4), fine = galerkin::triangleRule(8);
    const double k = 2 * std::numbers::pi * frequency / request.medium.soundSpeed;
    std::vector<PressureFieldSample> output;
    output.reserve(points.size());
    for (const auto point : points) {
        if (stop.stop_requested())
            throw std::runtime_error("Field evaluation cancelled.");
        if (!finite(point))
            throw std::invalid_argument("Invalid pressure field point.");
        PressureFieldSample sample{point, {}, FieldPointStatus::excluded};
        const FieldDomain *selected = nullptr;
        bool boundary = false;
        std::optional<Complex> interfacePressure;
        // Physical surfaces take priority, including interface rims. The
        // interface itself is transparent; only roundoff-close points use
        // its P1 trace, all other points follow adaptive volume integration.
        for (const auto &domain : domains) {
            for (std::size_t t = 0; t < domain.space.panels.size(); ++t) {
                const auto &panel = domain.space.panels[t];
                if (!domain.interfacePanels[t]) {
                    if (distanceToPanel(point, panel) <= settings.minimumBoundaryDistance)
                        boundary = true;
                    continue;
                }
                const double scale =
                    std::max({norm(point), norm(panel.vertices[0]), panel.maximumEdgeLength,
                              std::numeric_limits<double>::min()});
                const double roundoffDistance =
                    128 * std::numeric_limits<double>::epsilon() * scale;
                if (distanceToPanel(point, panel) > roundoffDistance)
                    continue;
                const auto a = panel.vertices[0], u = sub(panel.vertices[1], a),
                           v = sub(panel.vertices[2], a), w = sub(point, a);
                const double uu = dot(u, u), uv = dot(u, v), vv = dot(v, v), wu = dot(w, u),
                             wv = dot(w, v);
                const double determinant = uu * vv - uv * uv;
                if (!(determinant > 0))
                    throw std::invalid_argument("Degenerate interface field panel.");
                const double b = (wu * vv - wv * uv) / determinant,
                             c = (wv * uu - wu * uv) / determinant;
                const auto value = (1 - b - c) * domain.pressure[t][0] + b * domain.pressure[t][1] +
                                   c * domain.pressure[t][2];
                if (!finite(value))
                    throw std::invalid_argument("Non-finite interface pressure.");
                if (!interfacePressure)
                    interfacePressure = value;
            }
        }
        if (!boundary && !interfacePressure)
            for (const auto &domain : domains) {
                double angle = 0;
                for (const auto &panel : domain.space.panels)
                    angle += solidAngle(point, panel);
                const bool contains = domain.kind == BemDomainKind::interior
                                          ? std::abs(angle) > 2 * std::numbers::pi
                                          : std::abs(angle) < 2 * std::numbers::pi;
                if (contains) {
                    if (selected)
                        throw std::invalid_argument("Ambiguous subdomain field point.");
                    selected = &domain;
                }
            }
        if (boundary)
            sample.status = FieldPointStatus::boundary;
        else if (interfacePressure) {
            sample.status = FieldPointStatus::air;
            sample.pressure = *interfacePressure;
        } else if (selected) {
            sample.status = FieldPointStatus::air;
            for (std::size_t t = 0; t < selected->space.panels.size(); ++t)
                if (!adaptive(selected->space.panels[t], {}, point, k, selected->pressure[t],
                              selected->flux[t], false, coarse, fine, settings, 0, sample.pressure,
                              stop)) {
                    sample.status = FieldPointStatus::quadratureNotConverged;
                    sample.pressure = {};
                    break;
                }
            // All domain normals point out of air: representation is S q - D p.
            sample.pressure = -sample.pressure;
            if (!finite(sample.pressure)) {
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
} // namespace

std::vector<PressureFieldSample>
evaluatePressureField(const SurfaceMesh &mesh, const BemRequest &request, double frequency,
                      const BemExcitationResult &solution, const std::vector<Vec3> &points,
                      PressureFieldSettings settings, std::stop_token stop,
                      std::function<void(std::size_t, std::size_t)> progress) {
    if (!std::isfinite(frequency) || frequency <= 0 || !std::isfinite(request.medium.soundSpeed) ||
        request.medium.soundSpeed <= 0 || !std::isfinite(settings.minimumBoundaryDistance) ||
        settings.minimumBoundaryDistance <= 0 ||
        !std::isfinite(settings.relativeQuadratureTolerance) ||
        settings.relativeQuadratureTolerance <= 0 || settings.relativeQuadratureTolerance > .01 ||
        settings.maximumQuadratureDepth > 12)
        throw std::invalid_argument("Invalid pressure field settings.");
    if (request.symmetry)
        throw std::invalid_argument("Finite field evaluation with symmetry is not yet supported.");
    if (request.subdomains)
        return evaluateSubdomainField(mesh, request, frequency, solution, points, settings, stop,
                                      progress);
    if (!solution.boundary() || solution.boundary()->pressureTrace.size() != mesh.vertices.size() ||
        solution.boundary()->normalDerivative.size() != mesh.triangles.size())
        throw std::invalid_argument(
            "Complete BEM pressure and derivative traces are required. Recompute this result.");
    for (const auto &p : solution.boundary()->pressureTrace)
        if (!finite(p))
            throw std::invalid_argument("Non-finite pressure trace.");
    for (const auto &q : solution.boundary()->normalDerivative)
        if (!finite(q))
            throw std::invalid_argument("Non-finite derivative trace.");
    if (request.axisymmetric &&
        request.axisymmetric->surfaceContract == AxisymmetricSurfaceContract::continuousRevolution)
        return bem_detail::evaluateAxisymmetricPressureField(
            mesh, request, frequency, solution, points, settings, stop, std::move(progress));
    const auto space = galerkin::buildSpace(mesh);
    galerkin::Space exterior;
    std::map<std::uint32_t, Complex> exteriorVelocities;
    if (request.infiniteBaffle && !request.infiniteBaffle->exteriorDrivenMesh.triangles.empty()) {
        exterior = galerkin::buildSpace(request.infiniteBaffle->exteriorDrivenMesh);
        auto velocities = solution.solvedVelocities;
        if (velocities.empty() && !request.coupling) {
            const auto found = std::find_if(request.excitations.begin(), request.excitations.end(),
                                            [&](const auto &e) { return e.name == solution.name; });
            if (found != request.excitations.end())
                velocities = found->faces;
        }
        for (const auto &v : velocities)
            exteriorVelocities.emplace(v.faceId, v.normalVelocity);
        if (exteriorVelocities.empty())
            throw std::invalid_argument(
                "Exterior source velocities are required for field evaluation.");
    }
    const auto coarse = galerkin::triangleRule(4), fine = galerkin::triangleRule(8);
    const double k = 2 * std::numbers::pi * frequency / request.medium.soundSpeed;
    std::set<std::uint32_t> aperture;
    Vec3 normal{};
    if (request.infiniteBaffle) {
        const auto &b = *request.infiniteBaffle;
        if (!finite(b.origin) || !finite(b.normal) || norm(b.normal) <= 0 ||
            b.apertureFaceIds.empty())
            throw std::invalid_argument("Invalid field baffle.");
        normal = mul(b.normal, 1 / norm(b.normal));
        aperture.insert(b.apertureFaceIds.begin(), b.apertureFaceIds.end());
    }
    std::vector<PressureFieldSample> result;
    result.reserve(points.size());
    for (const auto point : points) {
        if (stop.stop_requested())
            throw std::runtime_error("Field evaluation cancelled.");
        if (!finite(point))
            throw std::invalid_argument("Invalid pressure field point.");
        PressureFieldSample sample{point, {}, FieldPointStatus::excluded};
        const double height =
            request.infiniteBaffle ? dot(sub(point, request.infiniteBaffle->origin), normal) : 0;
        const bool rayleigh = request.infiniteBaffle && height > settings.minimumBoundaryDistance;
        bool boundary = false;
        double angle = 0;
        for (const auto &panel : space.panels) {
            if (distanceToPanel(point, panel) <= settings.minimumBoundaryDistance) {
                boundary = true;
                break;
            }
            if (!rayleigh)
                angle += solidAngle(point, panel);
        }
        const bool inside = std::abs(angle) > 2 * std::numbers::pi;
        if (boundary ||
            (request.infiniteBaffle && std::abs(height) <= settings.minimumBoundaryDistance))
            sample.status = FieldPointStatus::boundary;
        else if (rayleigh || (request.infiniteBaffle ? inside : !inside)) {
            sample.status = FieldPointStatus::air;
            for (std::size_t t = 0; t < space.panels.size(); ++t) {
                const auto &panel = space.panels[t];
                if (rayleigh && !aperture.contains(panel.faceId))
                    continue;
                const std::array<Complex, 3> p{solution.boundary()->pressureTrace[panel.nodes[0]],
                                               solution.boundary()->pressureTrace[panel.nodes[1]],
                                               solution.boundary()->pressureTrace[panel.nodes[2]]};
                if (solution.boundary()->normalDerivative[t] == Complex{} &&
                    (rayleigh || (p[0] == Complex{} && p[1] == Complex{} && p[2] == Complex{})))
                    continue;
                if (!adaptive(panel, {}, point, k, p, solution.boundary()->normalDerivative[t],
                              rayleigh, coarse, fine, settings, 0, sample.pressure, stop)) {
                    sample.status = FieldPointStatus::quadratureNotConverged;
                    sample.pressure = {};
                    break;
                }
            }
            if (request.infiniteBaffle && !rayleigh)
                sample.pressure = -sample.pressure;
            if (rayleigh && sample.status == FieldPointStatus::air)
                for (const auto &panel : exterior.panels) {
                    const auto found = exteriorVelocities.find(panel.faceId);
                    if (found == exteriorVelocities.end() || found->second == Complex{})
                        continue;
                    const auto q =
                        Complex{0, -2 * std::numbers::pi * frequency * request.medium.density} *
                        found->second;
                    if (!adaptive(panel, {}, point, k, {}, q, true, coarse, fine, settings, 0,
                                  sample.pressure, stop)) {
                        sample.status = FieldPointStatus::quadratureNotConverged;
                        sample.pressure = {};
                        break;
                    }
                }
            if (!finite(sample.pressure)) {
                sample.status = FieldPointStatus::quadratureNotConverged;
                sample.pressure = {};
            }
        }
        result.push_back(sample);
        if (progress)
            progress(result.size(), points.size());
    }
    return result;
}
} // namespace hornsculpt
