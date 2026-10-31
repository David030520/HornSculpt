#include <algorithm>
#include <cmath>
#include <hornsculpt/level_palette.hpp>
#include <hornsculpt/pressure_field_plane.hpp>
#include <limits>
#include <stdexcept>

namespace hornsculpt {
void validatePressureFieldAppearance(const PressureFieldAppearance &s) {
    if (s.quantity < 0 || s.quantity > 2 || !std::isfinite(s.minimumDb) ||
        !std::isfinite(s.maximumDb) || s.maximumDb <= s.minimumDb || s.minimumDb < -300 ||
        s.maximumDb > 300 || !std::isfinite(s.phaseDegrees) || std::abs(s.phaseDegrees) > 180)
        throw std::invalid_argument("Ungültige Schallfeld-Farbskala.");
}
void validatePressureFieldPlane(const PressureFieldPlane &p) {
    const auto count = std::size_t(p.columns) * p.rows;
    if (p.columns < 2 || p.columns > 201 || p.rows < 2 || p.rows > 201 ||
        p.positions.size() != count || p.pressure.size() != count || p.valid.size() != count ||
        !std::isfinite(p.frequency) || p.frequency <= 0 || !std::isfinite(p.minimumLateral) ||
        !std::isfinite(p.maximumLateral) || !std::isfinite(p.minimumAxial) ||
        !std::isfinite(p.maximumAxial) || p.maximumLateral <= p.minimumLateral ||
        p.maximumAxial <= p.minimumAxial)
        throw std::invalid_argument("Ungültige Schallfeldebene.");
    const auto frame = normalizedAcousticFrame(p.frame);
    const auto lateral = p.vertical ? verticalDirection(frame) : frame.horizontal;
    const double planeExtent =
        std::max(p.maximumLateral - p.minimumLateral, p.maximumAxial - p.minimumAxial);
    if (!std::isfinite(planeExtent))
        throw std::invalid_argument("Nicht endliche Schallfeldausdehnung.");
    for (std::size_t i = 0; i < count; ++i) {
        const auto v = p.positions[i];
        if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z) ||
            (p.valid[i] &&
             (!std::isfinite(p.pressure[i].real()) || !std::isfinite(p.pressure[i].imag()))))
            throw std::invalid_argument("Nicht endlicher Schallfeldpunkt.");
        const auto row = i / p.columns, column = i % p.columns;
        const double axial =
            p.minimumAxial + (p.maximumAxial - p.minimumAxial) * double(row) / (p.rows - 1);
        const double across = p.minimumLateral + (p.maximumLateral - p.minimumLateral) *
                                                     double(column) / (p.columns - 1);
        const Vec3 expected{frame.origin.x + frame.forward.x * axial + lateral.x * across,
                            frame.origin.y + frame.forward.y * axial + lateral.y * across,
                            frame.origin.z + frame.forward.z * axial + lateral.z * across};
        if (!std::isfinite(expected.x) || !std::isfinite(expected.y) || !std::isfinite(expected.z))
            throw std::invalid_argument("Nicht endlicher Schallfeldrahmen.");
        const double worldScale =
            std::max({std::abs(v.x), std::abs(v.y), std::abs(v.z), std::abs(expected.x),
                      std::abs(expected.y), std::abs(expected.z), planeExtent});
        // SI allowance for a handful of multiply/add roundoffs, including a
        // translated CAD origin. This checks the presentation identity only;
        // it neither moves samples nor changes any acoustic tolerance.
        const double coordinateRoundoffMetres =
            16 * std::numeric_limits<double>::epsilon() * worldScale + 1e-11 * planeExtent;
        if (std::max({std::abs(v.x - expected.x), std::abs(v.y - expected.y),
                      std::abs(v.z - expected.z)}) > coordinateRoundoffMetres)
            throw std::invalid_argument(
                "Schallfeld-Weltpunkte stimmen nicht mit Rahmen und Raster überein.");
    }
}
std::vector<std::uint32_t> pressureFieldTriangleIndices(const PressureFieldPlane &p) {
    validatePressureFieldPlane(p);
    std::vector<std::uint32_t> indices;
    for (unsigned row = 0; row + 1 < p.rows; ++row)
        for (unsigned column = 0; column + 1 < p.columns; ++column) {
            const auto a = row * p.columns + column, b = a + 1, c = a + p.columns, d = c + 1;
            for (const auto triangle :
                 {std::array<std::uint32_t, 3>{a, b, d}, std::array<std::uint32_t, 3>{a, d, c}})
                if (p.valid[triangle[0]] && p.valid[triangle[1]] && p.valid[triangle[2]])
                    indices.insert(indices.end(), triangle.begin(), triangle.end());
        }
    return indices;
}
std::array<double, 2> pressureFieldRange(PressureFieldAppearance s) {
    validatePressureFieldAppearance(s);
    if (s.quantity == 1)
        return {-180, 180};
    if (s.quantity == 2) {
        const double maximum = 20e-6 * std::pow(10., s.maximumDb / 20);
        return {-maximum, maximum};
    }
    return {s.minimumDb, s.maximumDb};
}
QColor pressureFieldColor(double value, PressureFieldAppearance s) {
    const auto range = pressureFieldRange(s);
    // A zero-amplitude SPL is -infinity and clips to the floor; a phase or
    // real-pressure infinity has no physical colour.
    if (std::isnan(value) || (s.quantity != 0 && !std::isfinite(value)))
        return {};
    if (s.quantity == 1) {
        // Sector zero is centered on ±180°, so the wrap is inside a band.
        double wrapped = std::fmod(value + 198., 360.);
        if (wrapped < 0.)
            wrapped += 360.;
        const double sector = std::floor(wrapped / 36.);
        return steppedLevelColor(360. - (sector + .5) * 36., {0., 360., 36.});
    }
    if (s.quantity == 2)
        return steppedLevelColor(value / range[1], {-1., 1., .2});
    return steppedLevelColor(value, {range[0], range[1], 3.});
}
} // namespace hornsculpt
