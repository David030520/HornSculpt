#include "draft_slot.hpp"
#include <hornsculpt/draft_geometry.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <unordered_map>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace hornsculpt {
namespace {
constexpr double fullTurn = 2.0 * std::numbers::pi;
void checkTessellationStop(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Draft-Tessellation abgebrochen.");
}
bool validAngle(double angle) {
    return std::isfinite(angle) && angle >= 0.0 && angle < std::numbers::pi / 2.0;
}
OsseParameters profileParameters(const DraftParameters &parameters, bool vertical,
                                 bool opposite = false) {
    validateDiffractionSlot(parameters);
    OsseParameters result;
    static_cast<OsseDirectionalParameters &>(result) =
        opposite && parameters.asymmetricEnabled
            ? (vertical ? parameters.bottomProfile : parameters.leftProfile)
            : (vertical ? parameters.verticalProfile : parameters.profile);
    const auto &entry = parameters.slot.enabled ? parameters.slot.aperture : parameters.entry;
    result.throatRadius = .5 * (vertical ? entry.height : entry.width);
    result.length =
        parameters.length - (parameters.slot.enabled
                                 ? parameters.slot.adapterLength + parameters.slot.straightLength
                                 : 0);
    return result;
}
double axisComponent(double component) {
    return std::abs(component) < 8.0 * std::numeric_limits<double>::epsilon() ? 0.0
                                                                              : std::abs(component);
}

Vec3 difference(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double magnitude(Vec3 p) {
    return std::hypot(p.x, p.y, p.z);
}
double segmentDistance(Vec3 p, Vec3 a, Vec3 b) {
    const auto edge = difference(b, a), delta = difference(p, a);
    const double denominator = dot(edge, edge);
    if (!std::isfinite(denominator))
        throw std::invalid_argument("Draft-Tessellation ueberschreitet den Zahlenbereich.");
    const double t = denominator > 0 ? std::clamp(dot(delta, edge) / denominator, 0.0, 1.0) : 0;
    return magnitude({delta.x - t * edge.x, delta.y - t * edge.y, delta.z - t * edge.z});
}
double triangleDistance(Vec3 p, Vec3 a, Vec3 b, Vec3 c) {
    const auto n = cross(difference(b, a), difference(c, a));
    const double area = magnitude(n);
    if (!std::isfinite(area))
        throw std::invalid_argument("Draft-Tessellation ueberschreitet den Zahlenbereich.");
    if (area > 0 && dot(cross(difference(b, a), difference(p, a)), n) >= 0 &&
        dot(cross(difference(c, b), difference(p, b)), n) >= 0 &&
        dot(cross(difference(a, c), difference(p, c)), n) >= 0)
        return std::abs(dot(difference(p, a), n)) / area;
    return std::min({segmentDistance(p, a, b), segmentDistance(p, b, c), segmentDistance(p, c, a)});
}
Vec3 surfacePoint(const DraftGeometry &geometry, double z, double phi, bool outer,
                  const std::function<double(double, double)> &evaluate) {
    const double r = evaluate(z, phi) + (outer ? geometry.parameters().wallThickness : 0);
    const Vec3 result{r * std::cos(phi), r * std::sin(phi), z};
    if (!std::isfinite(result.x) || !std::isfinite(result.y) || !std::isfinite(result.z))
        throw std::invalid_argument("Draft-Tessellation erzeugt nicht endliche Koordinaten.");
    return result;
}
struct SurfaceGrid {
    std::vector<double> axial, angular;
};
SurfaceGrid initialGrid(const DraftGeometry &geometry) {
    SurfaceGrid grid;
    grid.axial = geometry.axialBreakpoints();
    for (unsigned i = 0; i <= 16; ++i)
        grid.axial.push_back(geometry.length() * i / 16.0);
    for (unsigned i = 0; i <= 64; ++i)
        grid.angular.push_back(fullTurn * i / 64.0);
    const auto &p = geometry.parameters();
    if (p.slot.enabled && p.slot.contour == SlotContour::roundedRectangle) {
        const double a = .5 * p.slot.aperture.width, b = .5 * p.slot.aperture.height;
        for (const double angle : {std::atan2(b - p.slot.cornerRadius, a),
                                   std::atan2(b, a - p.slot.cornerRadius)}) {
            grid.angular.push_back(angle);
            grid.angular.push_back(std::numbers::pi - angle);
            grid.angular.push_back(std::numbers::pi + angle);
            grid.angular.push_back(fullTurn - angle);
        }
        std::sort(grid.angular.begin(), grid.angular.end());
        grid.angular.erase(std::unique(grid.angular.begin(), grid.angular.end()), grid.angular.end());
    }
    const auto axialFeature = [&](double fraction) {
        grid.axial.push_back(geometry.length() * fraction);
    };
    if (p.transition.method == TransitionMethod::lut)
        for (const auto &knot : p.transition.knots)
            axialFeature(knot.fraction);
    else {
        axialFeature(p.transition.startFraction);
        axialFeature(p.transition.endFraction);
    }
    if (p.guideEnabled)
        grid.axial.push_back(geometry.guidePosition());
    if (p.mouthMorph.enabled) {
        axialFeature(p.mouthMorph.startFraction);
        axialFeature(p.mouthMorph.endFraction);
    }
    std::sort(grid.axial.begin(), grid.axial.end());
    grid.axial.erase(std::unique(grid.axial.begin(), grid.axial.end()), grid.axial.end());
    return grid;
}
struct RadiusKeyHash {
    std::size_t operator()(const std::pair<double, double> &key) const noexcept {
        const auto first = std::hash<double>{}(key.first);
        const auto second = std::hash<double>{}(key.second);
        return first ^ (second + 0x9e3779b9U + (first << 6) + (first >> 2));
    }
};
using RadiusCache = std::unordered_map<std::pair<double, double>, double, RadiusKeyHash>;
Vec3 cachedSurfacePoint(const DraftGeometry &geometry, double z, double phi, bool outer,
                        RadiusCache &cache,
                        const std::function<double(double, double)> *evaluate = nullptr) {
    auto [entry, created] = cache.try_emplace({z, phi});
    if (created)
        entry->second = evaluate ? (*evaluate)(z, phi) : geometry.radius(z, phi);
    const double r = entry->second + (outer ? geometry.parameters().wallThickness : 0);
    const Vec3 p{r * std::cos(phi), r * std::sin(phi), z};
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
        throw std::invalid_argument("Draft-Tessellation erzeugt nicht endliche Koordinaten.");
    return p;
}
std::pair<bool, bool> refineCell(const DraftGeometry &geometry, double z0, double z1, double phi0,
                                 double phi1, bool outer, double tolerance, RadiusCache &cache,
                                 const std::function<double(double, double)> &evaluate) {
    // A cell revisits the same quadrature points for edge and interior checks.
    // Retain their coordinates locally while preserving the shared radius cache.
    std::array<Vec3, 49> points{};
    std::array<bool, 49> computed{};
    const auto index = [](double value) {
        if (value == 0.0) return 0;
        if (value == .25) return 1;
        if (value == 1.0 / 3) return 2;
        if (value == .5) return 3;
        if (value == 2.0 / 3) return 4;
        if (value == .75) return 5;
        if (value == 1.0) return 6;
        throw std::logic_error("Unbekannter Tessellations-Prüfpunkt.");
    };
    const auto point = [&](double u, double v) {
        const auto at = static_cast<std::size_t>(index(u) * 7 + index(v));
        if (!computed[at]) {
            points[at] = cachedSurfacePoint(geometry, z0 + u * (z1 - z0),
                                            phi0 + v * (phi1 - phi0), outer, cache, &evaluate);
            computed[at] = true;
        }
        return points[at];
    };
    const auto a = point(0, 0), b = point(0, 1), c = point(1, 0), d = point(1, 1);
    double axialError = 0, angularError = 0, surfaceError = 0;
    for (const double fraction : {.25, .5, .75}) {
        for (const double other : {0.0, .5, 1.0}) {
            axialError = std::max(axialError, segmentDistance(point(fraction, other),
                                                              point(0, other), point(1, other)));
            angularError =
                std::max(angularError,
                         segmentDistance(point(other, fraction), point(other, 0), point(other, 1)));
        }
        for (const double other : {.25, .5, .75}) {
            const auto p = point(fraction, other);
            surfaceError =
                std::max(surfaceError, fraction + other <= 1 ? triangleDistance(p, a, b, c)
                                                             : triangleDistance(p, b, d, c));
            // A curved diagonal belongs to both emitted triangles. Its nearest
            // point can lie inside only one of them, so checking that half alone
            // underestimates the other triangle's final surface-deviation test.
            if (fraction + other == 1)
                surfaceError = std::max(
                    {surfaceError, triangleDistance(p, a, b, c), triangleDistance(p, b, d, c)});
        }
    }
    surfaceError = std::max({surfaceError, triangleDistance(point(1.0 / 3, 1.0 / 3), a, b, c),
                             triangleDistance(point(2.0 / 3, 2.0 / 3), b, d, c)});
    // Share the error budget between curved edges and triangle interiors.
    bool axial = axialError > .5 * tolerance, angular = angularError > .5 * tolerance;
    if (surfaceError > .75 * tolerance && !axial && !angular)
        axial = angular = true;
    return {axial, angular};
}
void splitIntervals(std::vector<double> &nodes, const std::vector<bool> &refine) {
    std::vector<double> next;
    next.reserve(2 * nodes.size());
    for (std::size_t i = 0; i + 1 < nodes.size(); ++i) {
        next.push_back(nodes[i]);
        if (refine[i]) {
            const double midpoint = std::midpoint(nodes[i], nodes[i + 1]);
            if (!(midpoint > nodes[i] && midpoint < nodes[i + 1]))
                throw std::invalid_argument("Draft-Tessellation erreicht die Zahlengenauigkeit.");
            next.push_back(midpoint);
        }
    }
    next.push_back(nodes.back());
    nodes = std::move(next);
}
SurfaceGrid sampleSurface(const DraftGeometry &geometry, double tolerance, std::stop_token stop,
                          const std::function<double(double, double)> &evaluate) {
    checkTessellationStop(stop);
    if (!std::isfinite(tolerance) || tolerance <= 0)
        throw std::invalid_argument(
            "Die geometrische Detailabweichung muss positiv und endlich sein.");
    auto grid = initialGrid(geometry);
    for (unsigned pass = 0; pass <= 24; ++pass) {
        checkTessellationStop(stop);
        // Includes both walls and closures; reject before allocating a huge mesh.
        if (grid.axial.size() > 8192 || grid.angular.size() > 8192 ||
            grid.axial.size() * grid.angular.size() > 500000)
            throw std::invalid_argument("Draft-Tessellation ueberschreitet das Geometriebudget.");
        std::vector<bool> axial(grid.axial.size() - 1), angular(grid.angular.size() - 1);
        bool changed = false;
        for (std::size_t i = 0; i < axial.size(); ++i) {
            checkTessellationStop(stop);
            // Reuse expensive composed-profile evaluations across adjacent cells
            // and both walls. A row bounds the cache independently of mesh size.
            RadiusCache cache;
            for (std::size_t j = 0; j < angular.size(); ++j)
                for (bool outer : {false, true}) {
                    const auto [u, v] =
                        refineCell(geometry, grid.axial[i], grid.axial[i + 1], grid.angular[j],
                                   grid.angular[j + 1], outer, tolerance, cache, evaluate);
                    axial[i] = axial[i] || u;
                    angular[j] = angular[j] || v;
                    changed = changed || u || v;
                }
        }
        if (!changed)
            return grid;
        if (pass == 24)
            throw std::invalid_argument(
                "Draft-Tessellation erreicht die geometrische Detailabweichung nicht.");
        splitIntervals(grid.axial, axial);
        splitIntervals(grid.angular, angular);
    }
    throw std::logic_error("Unerreichbarer Tessellationszustand.");
}
struct SampledRing {
    double z{};
    bool outer{};
    std::vector<Vec3> reference;
    std::vector<bool> kept;
    std::vector<std::uint32_t> output;

    std::vector<std::size_t> indices() const {
        std::vector<std::size_t> result;
        for (std::size_t i = 0; i < kept.size(); ++i)
            if (kept[i])
                result.push_back(i);
        return result;
    }
};
SampledRing sampleRing(const DraftGeometry &geometry, double z, bool outer, const SurfaceGrid &grid,
                       double tolerance, std::stop_token stop,
                       const std::function<double(double, double)> &evaluate) {
    checkTessellationStop(stop);
    SampledRing ring;
    ring.z = z;
    ring.outer = outer;
    for (std::size_t i = 0; i + 1 < grid.angular.size(); ++i) {
        auto p = surfacePoint(geometry, std::max(0.0, z), grid.angular[i], outer, evaluate);
        p.z = z;
        ring.reference.push_back(p);
    }
    ring.reference.push_back(ring.reference.front());
    ring.kept.resize(grid.angular.size());
    std::vector<std::size_t> anchors{0};
    for (unsigned quadrant = 1; quadrant < 4; ++quadrant)
        anchors.push_back(static_cast<std::size_t>(
            std::lower_bound(grid.angular.begin(), grid.angular.end(), quadrant * fullTurn / 4) -
            grid.angular.begin()));
    anchors.push_back(grid.angular.size() - 1);
    std::vector<std::pair<std::size_t, std::size_t>> intervals;
    for (std::size_t i = 1; i < anchors.size(); ++i) {
        ring.kept[anchors[i - 1]] = ring.kept[anchors[i]] = true;
        intervals.emplace_back(anchors[i - 1], anchors[i]);
    }
    while (!intervals.empty()) {
        checkTessellationStop(stop);
        const auto [first, last] = intervals.back();
        intervals.pop_back();
        double error = 0;
        std::size_t worst = first;
        for (std::size_t i = first + 1; i < last; ++i) {
            const double distance =
                segmentDistance(ring.reference[i], ring.reference[first], ring.reference[last]);
            if (distance > error) {
                error = distance;
                worst = i;
            }
        }
        // The accurate reference curve already consumes at most half the
        // tolerance. Local simplification reserves another quarter for its chord.
        if (error > .25 * tolerance) {
            ring.kept[worst] = true;
            intervals.emplace_back(first, worst);
            intervals.emplace_back(worst, last);
        }
    }
    return ring;
}
struct RingVertex {
    const SampledRing *ring;
    std::size_t index;
};
using RingTriangle = std::array<RingVertex, 3>;
template <class Emit>
void zipRings(const SampledRing &first, const SampledRing &next, const std::vector<double> &angles,
              Emit emit) {
    const auto a = first.indices(), b = next.indices();
    std::size_t i = 0, j = 0;
    while (i + 1 < a.size() || j + 1 < b.size()) {
        const RingVertex p{&first, a[i]}, q{&next, b[j]};
        if (i + 1 < a.size() && (j + 1 == b.size() || angles[a[i + 1]] <= angles[b[j + 1]])) {
            emit(RingTriangle{p, RingVertex{&first, a[i + 1]}, q});
            ++i;
        } else {
            emit(RingTriangle{p, RingVertex{&next, b[j + 1]}, q});
            ++j;
        }
    }
}
void restoreArc(SampledRing &ring, const std::vector<double> &angles, double first, double last) {
    auto begin = std::lower_bound(angles.begin(), angles.end(), first);
    for (auto it = begin; it != angles.end() && *it <= last; ++it) {
        const auto i = static_cast<std::size_t>(it - angles.begin());
        ring.kept[i] = true;
    }
}
double triangleSurfaceError(const DraftGeometry &geometry, const std::vector<double> &angles,
                            const RingTriangle &triangle, RadiusCache &cache,
                            const std::function<double(double, double)> &evaluate) {
    std::array<Vec3, 3> points;
    for (std::size_t i = 0; i < 3; ++i)
        points[i] = triangle[i].ring->reference[triangle[i].index];
    const auto distance = [&](std::array<double, 3> weights) {
        double z = 0, phi = 0;
        for (std::size_t i = 0; i < 3; ++i) {
            z += weights[i] * triangle[i].ring->z;
            phi += weights[i] * angles[triangle[i].index];
        }
        auto p = cachedSurfacePoint(geometry, std::clamp(z, 0.0, geometry.length()), phi,
                                    triangle[0].ring->outer, cache, &evaluate);
        p.z = z;
        return triangleDistance(p, points[0], points[1], points[2]);
    };
    double error = distance({1.0 / 3, 1.0 / 3, 1.0 / 3});
    for (std::size_t i = 0; i < 3; ++i) {
        for (double t : {.25, .5, .75}) {
            std::array<double, 3> weights{};
            weights[i] = 1 - t;
            weights[(i + 1) % 3] = t;
            error = std::max(error, distance(weights));
        }
        std::array<double, 3> weights{.25, .25, .25};
        weights[i] = .5;
        error = std::max(error, distance(weights));
    }
    return error;
}
bool validateStrip(const DraftGeometry &geometry, SampledRing &first, SampledRing &next,
                   const std::vector<double> &angles, double tolerance, std::stop_token stop,
                   const std::function<double(double, double)> &evaluate) {
    checkTessellationStop(stop);
    RadiusCache cache;
    struct Refinement {
        double first, last, error;
    };
    std::vector<Refinement> required;
    zipRings(first, next, angles, [&](const RingTriangle &triangle) {
        const double error = triangleSurfaceError(geometry, angles, triangle, cache, evaluate);
        if (error <= tolerance)
            return;
        double lo = fullTurn, hi = 0;
        for (const auto vertex : triangle) {
            lo = std::min(lo, angles[vertex.index]);
            hi = std::max(hi, angles[vertex.index]);
        }
        required.push_back({lo, hi, error});
    });
    // Diagnose against the immutable strip, then apply all refinements together.
    // An overlapping earlier request must not make a later stale cell look final.
    for (const auto &request : required) {
        bool refinable = false;
        for (auto it = std::lower_bound(angles.begin(), angles.end(), request.first);
             it != angles.end() && *it <= request.last; ++it) {
            const auto i = static_cast<std::size_t>(it - angles.begin());
            refinable = refinable || !first.kept[i] || !next.kept[i];
        }
        if (!refinable)
            throw std::invalid_argument(
                "Draft-Tessellation erreicht die Flaechenabweichung nicht: z=" +
                std::to_string(first.z) + ".." + std::to_string(next.z) +
                " m, phi=" + std::to_string(request.first) + ".." + std::to_string(request.last) +
                " rad, Abweichung=" + std::to_string(request.error) + " m.");
    }
    for (const auto &request : required) {
        restoreArc(first, angles, request.first, request.last);
        restoreArc(next, angles, request.first, request.last);
    }
    return !required.empty();
}
void conformSurface(const DraftGeometry &geometry, std::vector<SampledRing> &inner,
                    std::vector<SampledRing> &outer, const std::vector<double> &angles,
                    double tolerance, std::stop_token stop,
                    const std::function<double(double, double)> &evaluate) {
    const auto availableKnots = angles.size() * (inner.size() + outer.size());
    for (std::size_t pass = 0; pass <= availableKnots; ++pass) {
        checkTessellationStop(stop);
        // Both boundaries of the planar lip share one angular partition; this
        // preserves a positive annular strip even for very thin material.
        for (std::size_t i = 0; i < angles.size(); ++i)
            inner.back().kept[i] = outer.back().kept[i] =
                inner.back().kept[i] || outer.back().kept[i];
        bool changed = false;
        for (std::size_t i = 1; i < inner.size(); ++i)
            changed =
                validateStrip(geometry, inner[i - 1], inner[i], angles, tolerance, stop, evaluate) || changed;
        for (std::size_t i = 1; i < outer.size(); ++i)
            changed =
                validateStrip(geometry, outer[i - 1], outer[i], angles, tolerance, stop, evaluate) || changed;
        if (!changed)
            return;
    }
    throw std::runtime_error("Draft-Tessellation ueberschreitet das Konformitaetsbudget.");
}
} // namespace

DraftGeometry::DraftGeometry(DraftParameters parameters)
    : DraftGeometry(parameters, std::make_shared<OsseProfile>(profileParameters(parameters, false)),
                    std::make_shared<OsseProfile>(profileParameters(parameters, true)),
                    std::make_shared<OsseProfile>(profileParameters(parameters, false, true)),
                    std::make_shared<OsseProfile>(profileParameters(parameters, true, true))) {}

DraftGeometry::DraftGeometry(DraftParameters parameters,
                             std::shared_ptr<const LongitudinalProfile> profile)
    : DraftGeometry(std::move(parameters), profile, profile) {}

DraftGeometry::DraftGeometry(DraftParameters parameters,
                             std::shared_ptr<const LongitudinalProfile> horizontalProfile,
                             std::shared_ptr<const LongitudinalProfile> verticalProfile)
    : DraftGeometry(std::move(parameters), horizontalProfile, verticalProfile, horizontalProfile,
                    verticalProfile) {}

DraftGeometry::DraftGeometry(DraftParameters parameters,
                             std::shared_ptr<const LongitudinalProfile> rightProfile,
                             std::shared_ptr<const LongitudinalProfile> topProfile,
                             std::shared_ptr<const LongitudinalProfile> leftProfile,
                             std::shared_ptr<const LongitudinalProfile> bottomProfile)
    : parameters_(std::move(parameters)), horizontalProfile_(std::move(rightProfile)),
      verticalProfile_(std::move(topProfile)), leftProfile_(std::move(leftProfile)),
      bottomProfile_(std::move(bottomProfile)), transition_(parameters_.transition) {
    const auto &p = parameters_;
    validateDiffractionSlot(p);
    if (!std::isfinite(p.length) || p.length <= 0.0 || !std::isfinite(p.crossSectionExponent) ||
        p.crossSectionExponent < 2.0 || !std::isfinite(p.wallThickness) || p.wallThickness <= 0.0)
        throw std::invalid_argument("Draft: Profil, Winkel, Querschnitt oder Wanddicke ungueltig.");
    (void)Superellipse(.5 * p.entry.width, .5 * p.entry.height, p.entry.exponent);
    const LongitudinalProfile *profiles[]{horizontalProfile_.get(), verticalProfile_.get(),
                                          leftProfile_.get(), bottomProfile_.get()};
    const double angles[]{p.horizontalAngle, p.verticalAngle, p.leftAngle, p.bottomAngle};
    for (unsigned i = 0; i < (p.asymmetricEnabled ? 4u : 2u); ++i) {
        const auto *profile = profiles[i];
        if (!profile || profile->length() != p.length - profileStart() || !validAngle(angles[i]))
            throw std::invalid_argument("Draft: Richtungsprofil, Laenge oder Winkel ungueltig.");
        const double throat = profile->radius(0.0, angles[i]);
        const auto &entry = p.slot.enabled ? p.slot.aperture : p.entry;
        const double entryRadius = .5 * (i % 2 == 0 ? entry.width : entry.height);
        if (!std::isfinite(throat) || throat <= 0.0 ||
            std::abs(throat - entryRadius) > 1e-12 * throat)
            throw std::invalid_argument(
                "Draft: Profile muessen zu den Eintrittshalbachsen des Profilabschnitts passen.");
    }
    initializeSlotExit();
    if (p.guideEnabled) {
        const auto &guide = p.guide;
        if (!std::isfinite(guide.width) || guide.width <= 0.0 ||
            !std::isfinite(guide.aspectRatio) || guide.aspectRatio <= 0.0 ||
            !std::isfinite(guide.rotation) || !std::isfinite(guide.position) ||
            guidePosition() <= 0.0 || guidePosition() > length())
            throw std::invalid_argument("Leitkurve: Groesse, Rotation oder Position ungueltig.");
        try {
            guideCurve_.emplace(guide);
        } catch (const std::invalid_argument &error) {
            throw std::invalid_argument(std::string("Leitkurve: ") + error.what());
        }
        const double blend = transition_.blend(guidePosition() / length());
        if (blend <= 1e-12)
            throw std::invalid_argument("Leitkurve liegt vor Beginn des Formuebergangs.");
    }
    if (p.mouthMorph.enabled) {
        if (p.guideEnabled && p.mouthMorph.startFraction < guidePosition() / length())
            throw std::invalid_argument("Mundmorphing muss an oder hinter der Leitkurve beginnen.");
        try {
            mouthMorph_.emplace(p.mouthMorph);
        } catch (const std::invalid_argument &error) {
            throw std::invalid_argument(std::string("Mund-Morphing: ") + error.what());
        }
    }
    // Positive polar radii and z-monotone rings rule out radial self-intersections.
    // Sample the complete composition, including guide inversion and blend overshoot.
    constexpr std::size_t validationAngles = 720;
    const auto breaks = axialBreakpoints();
    for (std::size_t i = 0; i < validationAngles; ++i) {
        const double phi = fullTurn * static_cast<double>(i) / validationAngles;
        for (std::size_t segment = 1; segment < breaks.size(); ++segment) {
            const auto stations =
                std::max(16u, static_cast<unsigned>(std::ceil(
                                  64 * (breaks[segment] - breaks[segment - 1]) / length())));
            for (unsigned j = 0; j <= stations; ++j) {
                const double r = radius(std::lerp(breaks[segment - 1], breaks[segment],
                                                  static_cast<double>(j) / stations),
                                        phi);
                if (!std::isfinite(r) || r <= 0.0)
                    throw std::invalid_argument("Draft erzeugt keinen positiven endlichen Radius.");
            }
        }
    }
}

double DraftGeometry::radius(double z, double phi) const {
    // Radius-only evaluation dominates adaptive display tessellation. The
    // untransformed path shares the same composition without computing slopes.
    if (!parameters_.slot.enabled && !parameters_.guideEnabled && !mouthMorph_)
        return baseSample(z, phi, false).radius;
    return sample(z, phi).radius;
}

RadialSample DraftGeometry::crossSection(double z, double phi, double exponent,
                                         bool useSlotContour, bool withSlope,
                                         const std::array<double, 2> *cachedAxes) const {
    const double localZ = profilePosition(z);
    const bool left = parameters_.asymmetricEnabled && std::cos(phi) < 0;
    const bool bottom = parameters_.asymmetricEnabled && std::sin(phi) < 0;
    const auto &horizontal = left ? *leftProfile_ : *horizontalProfile_;
    const auto &vertical = bottom ? *bottomProfile_ : *verticalProfile_;
    const double horizontalAngle = left ? parameters_.leftAngle : parameters_.horizontalAngle;
    const double verticalAngle = bottom ? parameters_.bottomAngle : parameters_.verticalAngle;
    const double rh = cachedAxes ? (*cachedAxes)[0] : horizontal.radius(localZ, horizontalAngle);
    const double rv = cachedAxes ? (*cachedAxes)[1] : vertical.radius(localZ, verticalAngle);
    if (useSlotContour && parameters_.slot.enabled &&
        parameters_.slot.contour == SlotContour::roundedRectangle) {
        const auto &slot = parameters_.slot;
        const double a = slot.aperture.width * .5, b = slot.aperture.height * .5;
        const RoundedRectangle contour(a, b, slot.cornerRadius);
        // Affine expansion preserves the exact slot contour at the interface.
        // Beyond the straight slot its circular corners become elliptic arcs.
        const double qx = std::cos(phi) * a / rh, qy = std::sin(phi) * b / rv;
        const double norm = std::hypot(qx, qy), theta = std::atan2(qy, qx);
        const double dx = -qx * horizontal.slope(localZ, horizontalAngle) / rh;
        const double dy = -qy * vertical.slope(localZ, verticalAngle) / rv;
        const double axialDerivative = profilePositionDerivative(z);
        const double value = contour.radius(theta);
        const double angularRate = (qx * dy - qy * dx) / (norm * norm);
        const double relativeRate = (qx * dx + qy * dy) / (norm * norm);
        return {value / norm, (contour.angularDerivative(theta) * angularRate -
                              value * relativeRate) / norm * axialDerivative};
    }
    // Adjacent quadrant superellipses share the cardinal radius and tangent.
    // Exponents >= 2 make the disappearing half-axis contribution C1 at a seam.
    const double r = Superellipse(rh, rv, exponent).radius(phi);
    if (!withSlope)
        return {r, 0.0};
    const double axialDerivative = profilePositionDerivative(z);
    const double x = axisComponent(std::cos(phi)) / rh;
    const double y = axisComponent(std::sin(phi)) / rv;
    const double largest = std::max(x, y);
    const double xp = std::pow(x / largest, exponent), yp = std::pow(y / largest, exponent);
    double slope = 0.0;
    if (xp != 0.0)
        slope += r * xp / (xp + yp) * horizontal.slope(localZ, horizontalAngle) / rh;
    if (yp != 0.0)
        slope += r * yp / (xp + yp) * vertical.slope(localZ, verticalAngle) / rv;
    return {r, slope * axialDerivative};
}

RadialSample DraftGeometry::sample(double z, double phi) const {
    const auto base = baseSample(z, phi);
    if (mouthMorph_)
        return mouthMorph_->apply(base, baseSample(length(), phi).radius, z / length(), length(),
                                  phi);
    return base;
}

RadialSample DraftGeometry::guidedSample(double z, double phi, double targetGuideRadius) const {
    const double guideZ = profilePosition(guidePosition()), localZ = profilePosition(z);
    const double cosine = axisComponent(std::cos(phi)), sine = axisComponent(std::sin(phi));
    const double weights[2]{cosine * cosine, sine * sine};
    const LongitudinalProfile *profiles[]{
        parameters_.asymmetricEnabled && std::cos(phi) < 0 ? leftProfile_.get()
                                                           : horizontalProfile_.get(),
        parameters_.asymmetricEnabled && std::sin(phi) < 0 ? bottomProfile_.get()
                                                           : verticalProfile_.get()};
    double zeroAngle[2]{}, baseline = 0;
    for (unsigned i = 0; i < 2; ++i)
        if (weights[i] > 0) {
            zeroAngle[i] = profiles[i]->radius(guideZ, 0);
            baseline += weights[i] * zeroAngle[i];
        }
    const double increment = targetGuideRadius - baseline;
    const double roundoff = 64 * std::numeric_limits<double>::epsilon() *
                            std::max(std::abs(targetGuideRadius), std::abs(baseline));
    if (increment < -roundoff)
        throw std::invalid_argument("Leitkurve liegt unter dem erreichbaren Richtungsprofil.");
    RadialSample result{};
    // Apply one radial opening increment above each provider's zero-angle
    // baseline. Requiring both directional profiles to reach the same absolute
    // polar radius incorrectly rejects valid guides with unequal entry axes.
    // The weighted guide radii still sum exactly to targetGuideRadius.
    for (unsigned i = 0; i < 2; ++i)
        if (weights[i] > 0) {
            const double angle =
                profiles[i]->angleThrough(guideZ, zeroAngle[i] + std::max(0., increment));
            result.radius += weights[i] * profiles[i]->radius(localZ, angle);
            result.slope +=
                weights[i] * profiles[i]->slope(localZ, angle) * profilePositionDerivative(z);
        }
    return result;
}

RadialSample DraftGeometry::baseSample(double z, double phi, bool withSlope,
                                       const std::array<double, 2> *cachedAxes) const {
    if (!std::isfinite(z) || z < 0.0 || z > length() || !std::isfinite(phi))
        throw std::invalid_argument("Draft: Position ausserhalb der Geometrie.");
    if (const auto fillet = filletSample(z, phi))
        return *fillet;
    if (const auto upstream = upstreamSample(z, phi))
        return *upstream;
    const double exponent =
        parameters_.slot.enabled ? parameters_.slot.aperture.exponent : parameters_.entry.exponent;
    const auto base = crossSection(z, phi, exponent, true, withSlope, cachedAxes);
    const double blend = transition_.blend(z / length());
    if (blend == 0.0)
        return base;
    RadialSample target{};
    if (parameters_.guideEnabled) {
        const double guideZ = guidePosition();
        const double guideBlend = transition_.blend(guideZ / length());
        const double baseGuide = crossSection(guideZ, phi, exponent, true).radius;
        const double targetGuide = baseGuide + (guideRadius(phi) - baseGuide) / guideBlend;
        target = guidedSample(z, phi, targetGuide);
    } else {
        target = crossSection(z, phi, parameters_.crossSectionExponent, false, withSlope, cachedAxes);
    }
    if (!withSlope)
        return {base.radius + blend * (target.radius - base.radius), 0.0};
    const double blendSlope = transition_.derivative(z / length()) / length();
    double slope = blendSlope == 0.0 ? 0.0 : blendSlope * (target.radius - base.radius);
    if (blend != 1.0)
        slope += (1.0 - blend) * base.slope;
    slope += blend * target.slope;
    return {base.radius + blend * (target.radius - base.radius), slope};
}

double DraftGeometry::length() const noexcept {
    return parameters_.length;
}

double DraftGeometry::guidePosition() const {
    return parameters_.guide.positionMode == GuidePositionMode::relative
               ? parameters_.guide.position * length()
               : parameters_.guide.position;
}

double DraftGeometry::guideRadius(double phi) const {
    if (!parameters_.guideEnabled || !guideCurve_)
        throw std::logic_error("Keine Leitkurve aktiv.");
    return guideCurve_->radius(phi);
}

const DraftParameters &DraftGeometry::parameters() const noexcept {
    return parameters_;
}

SurfaceMesh DraftGeometry::previewMesh(double geometricTolerance, std::stop_token stop) const {
    // Axial profile radii are invariant across angular samples. Keep this
    // cache local to this immutable preview job; transformed cases use the
    // full evaluator because their axial composition depends on angle.
    const bool cacheAxes = !parameters_.slot.enabled && !parameters_.guideEnabled && !mouthMorph_;
    std::unordered_map<double, std::array<double, 4>> axialRadii;
    const std::function<double(double, double)> evaluate = [&](double z, double phi) {
        if (!cacheAxes)
            return radius(z, phi);
        auto [entry, inserted] = axialRadii.try_emplace(z);
        if (inserted) {
            const double localZ = profilePosition(z);
            auto &r = entry->second;
            r[0] = horizontalProfile_->radius(localZ, parameters_.horizontalAngle);
            r[1] = verticalProfile_->radius(localZ, parameters_.verticalAngle);
            r[2] = parameters_.asymmetricEnabled
                       ? leftProfile_->radius(localZ, parameters_.leftAngle) : r[0];
            r[3] = parameters_.asymmetricEnabled
                       ? bottomProfile_->radius(localZ, parameters_.bottomAngle) : r[1];
        }
        const auto &r = entry->second;
        const std::array<double, 2> axes{
            r[parameters_.asymmetricEnabled && std::cos(phi) < 0 ? 2 : 0],
            r[parameters_.asymmetricEnabled && std::sin(phi) < 0 ? 3 : 1]};
        return baseSample(z, phi, false, &axes).radius;
    };
    const auto grid = sampleSurface(*this, geometricTolerance, stop, evaluate);
    std::vector<SampledRing> inner, outer;
    for (const auto z : grid.axial)
        inner.push_back(sampleRing(*this, z, false, grid, geometricTolerance, stop, evaluate));
    outer.push_back(
        sampleRing(*this, -parameters_.wallThickness, true, grid, geometricTolerance, stop, evaluate));
    for (const auto z : grid.axial)
        outer.push_back(sampleRing(*this, z, true, grid, geometricTolerance, stop, evaluate));
    conformSurface(*this, inner, outer, grid.angular, geometricTolerance, stop, evaluate);
    SurfaceMesh mesh;
    const auto emitRing = [&](SampledRing &ring) {
        checkTessellationStop(stop);
        ring.output.resize(grid.angular.size());
        for (auto i : ring.indices()) {
            if (i + 1 == grid.angular.size())
                continue;
            ring.output[i] = static_cast<std::uint32_t>(mesh.vertices.size());
            mesh.vertices.push_back(ring.reference[i]);
        }
        ring.output.back() = ring.output.front();
    };
    for (auto &ring : inner)
        emitRing(ring);
    for (auto &ring : outer)
        emitRing(ring);
    const auto triangle = [&](std::array<std::uint32_t, 3> vertices, DraftFace face) {
        mesh.triangles.push_back(
            {vertices, static_cast<std::uint32_t>(face),
             face == DraftFace::driven ? BoundaryCondition::driven : BoundaryCondition::wall});
    };
    const auto join = [&](const SampledRing &first, const SampledRing &next, DraftFace face,
                          bool reverse) {
        checkTessellationStop(stop);
        zipRings(first, next, grid.angular, [&](const RingTriangle &points) {
            std::array<std::uint32_t, 3> vertices;
            for (std::size_t i = 0; i < 3; ++i)
                vertices[i] = points[i].ring->output[points[i].index];
            if (reverse)
                std::swap(vertices[1], vertices[2]);
            triangle(vertices, face);
        });
    };
    for (auto i : inner.front().indices())
        if (i + 1 < grid.angular.size())
            mesh.fixedVertices.push_back(inner.front().output[i]);
    for (std::size_t i = 1; i < inner.size(); ++i) {
        const auto face = !parameters_.slot.enabled || inner[i].z > slotStraightEnd()
                              ? DraftFace::innerWall
                          : inner[i].z <= parameters_.slot.adapterLength ? DraftFace::adapterWall
                                                                         : DraftFace::slotWall;
        join(inner[i - 1], inner[i], face, true);
    }
    for (std::size_t i = 1; i < outer.size(); ++i)
        join(outer[i - 1], outer[i],
             parameters_.slot.enabled && outer[i].z <= slotStraightEnd() ? DraftFace::upstreamOuterWall
                                                                      : DraftFace::outerWall,
             false);
    join(inner.back(), outer.back(), DraftFace::mouthRim, true);
    const auto cap = [&](const SampledRing &ring, DraftFace face, bool reverse) {
        const auto center = static_cast<std::uint32_t>(mesh.vertices.size());
        mesh.vertices.push_back({0, 0, ring.z});
        const auto indices = ring.indices();
        for (std::size_t i = 1; i < indices.size(); ++i) {
            if (i % 256 == 0)
                checkTessellationStop(stop);
            auto a = ring.output[indices[i - 1]], b = ring.output[indices[i]];
            if (reverse)
                std::swap(a, b);
            triangle({center, a, b}, face);
        }
    };
    cap(inner.front(), DraftFace::driven, false);
    cap(outer.front(), DraftFace::rear, true);
    return mesh;
}

} // namespace hornsculpt
