#include <hornsculpt/axisymmetric_geometry.hpp>
#include <hornsculpt/draft_geometry.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numbers>
#include <stdexcept>

namespace hornsculpt {
AxisymmetricMeridian makeAxisymmetricMeridian(const DraftGeometry &geometry,
                                              const MeridianSampling &s, std::stop_token stop) {
    if (!std::isfinite(s.maximumLength) || s.maximumLength <= 0 ||
        !std::isfinite(s.chordTolerance) || s.chordTolerance <= 0 ||
        !std::isfinite(s.rotationalTolerance) || s.rotationalTolerance <= 0)
        throw std::invalid_argument("2D-Draft: positive SI-Genauigkeiten erforderlich.");
    const auto cancel = [&] {
        if (stop.stop_requested())
            throw std::runtime_error("2D-Draft abgebrochen.");
    };
    cancel();
    const auto &p = geometry.parameters();
    // A circular entry and circular direct cross-sections are necessary even
    // when an enabled transition happens to mask them at one sampled station.
    if (std::abs(p.entry.width - p.entry.height) > s.rotationalTolerance || p.entry.exponent != 2 ||
        p.crossSectionExponent != 2 ||
        (p.slot.enabled &&
         (std::abs(p.slot.aperture.width - p.slot.aperture.height) > s.rotationalTolerance ||
          (p.slot.contour == SlotContour::superellipse
               ? p.slot.aperture.exponent != 2
               : std::abs(p.slot.cornerRadius - .5 * p.slot.aperture.width) >
                     s.rotationalTolerance))))
        throw std::invalid_argument("2D-Rotation verlangt kreisförmige Draft-Querschnitte.");
    const auto radius = [&](double z) {
        cancel();
        const double r = geometry.radius(z, 0);
        for (unsigned j = 1; j < 127; ++j)
            if (std::abs(geometry.radius(z, 2 * std::numbers::pi * j / 127) - r) >
                s.rotationalTolerance)
                throw std::invalid_argument(
                    "2D-Rotation: Draft-Profile, Leitkurve oder Morphing sind richtungsabhängig.");
        return r;
    };
    // Validate the entire final composition, including transition and morph.
    for (unsigned i = 0; i <= 64; ++i)
        (void)radius(geometry.length() * i / 64);
    const auto breaks = geometry.axialBreakpoints();
    for (std::size_t i = 1; i < breaks.size(); ++i)
        for (unsigned j = 0; j <= 16; ++j)
            (void)radius(std::lerp(breaks[i - 1], breaks[i], j / 16.0));
    AxisymmetricMeridian out;
    const double thickness = p.wallThickness;
    out.points.push_back({0, -thickness});
    const auto append = [&](MeridianPoint point, DraftFace face) {
        const auto first = static_cast<std::uint32_t>(out.points.size() - 1);
        out.points.push_back(point);
        const auto id = static_cast<std::uint32_t>(face);
        const auto boundary = s.faceConditions.contains(id) ? s.faceConditions.at(id)
                              : face == DraftFace::driven   ? BoundaryCondition::driven
                                                            : BoundaryCondition::wall;
        out.segments.push_back({first, first + 1, id, boundary});
    };
    const auto sampleCurve = [&](auto evaluate, DraftFace face) {
        const auto id = static_cast<std::uint32_t>(face);
        const double target = s.faceTargetLengths.contains(id)
                                  ? std::min(s.maximumLength, s.faceTargetLengths.at(id))
                                  : s.maximumLength;
        if (!std::isfinite(target) || target <= 0)
            throw std::invalid_argument("2D-Draft: lokale Länge ungültig.");
        std::function<void(double, double, unsigned)> split = [&](double a, double b,
                                                                  unsigned depth) {
            cancel();
            const auto u = evaluate(a), v = evaluate(b);
            double error = 0;
            for (double fraction : {.25, .5, .75}) {
                const auto x = evaluate(a + (b - a) * fraction);
                error = std::max(
                    error, std::hypot(x.radius - ((1 - fraction) * u.radius + fraction * v.radius),
                                      x.axial - ((1 - fraction) * u.axial + fraction * v.axial)));
            }
            if (std::hypot(v.radius - u.radius, v.axial - u.axial) > target ||
                error > s.chordTolerance) {
                if (depth >= 24)
                    throw std::runtime_error("2D-Draft: Abschnitt erreicht Genauigkeit nicht.");
                const double mid = (a + b) * .5;
                split(a, mid, depth + 1);
                split(mid, b, depth + 1);
            } else
                append(v, face);
        };
        split(0, 1, 0);
    };
    const double r0 = radius(0), rL = radius(geometry.length());
    sampleCurve([&](double t) { return MeridianPoint{t * (r0 + thickness), -thickness}; },
                DraftFace::rear);
    sampleCurve([&](double t) { return MeridianPoint{r0 + thickness, (t - 1) * thickness}; },
                p.slot.enabled ? DraftFace::upstreamOuterWall : DraftFace::outerWall);
    const double slotEnd = geometry.slotStraightEnd();
    for (std::size_t i = 1; i < breaks.size(); ++i)
        sampleCurve(
            [&](double t) {
                const double z = std::lerp(breaks[i - 1], breaks[i], t);
                return MeridianPoint{radius(z) + thickness, z};
            },
            p.slot.enabled && breaks[i] <= slotEnd ? DraftFace::upstreamOuterWall
                                                   : DraftFace::outerWall);
    sampleCurve(
        [&](double t) { return MeridianPoint{rL + (1 - t) * thickness, geometry.length()}; },
        DraftFace::mouthRim);
    for (std::size_t i = breaks.size() - 1; i > 0; --i)
        sampleCurve(
            [&](double t) {
                const double z = std::lerp(breaks[i], breaks[i - 1], t);
                return MeridianPoint{radius(z), z};
            },
            !p.slot.enabled || breaks[i] > slotEnd ? DraftFace::innerWall
            : breaks[i] <= p.slot.adapterLength    ? DraftFace::adapterWall
                                                   : DraftFace::slotWall);
    sampleCurve([&](double t) { return MeridianPoint{(1 - t) * r0, 0}; }, DraftFace::driven);
    return out;
}
} // namespace hornsculpt
