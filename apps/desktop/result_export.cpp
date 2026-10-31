#include "result_export.hpp"
#include "simulation_request.hpp"
#include <QIODevice>
#include <QLocale>
#include <QTextStream>
#include <algorithm>
#include <cmath>
#include <hornsculpt/acoustic_analysis.hpp>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace hornsculpt {
namespace {
Vec3 scaled(Vec3 a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}
Vec3 added(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
std::vector<double> anglesDegrees(std::size_t count) {
    std::vector<double> values;
    for (std::size_t i = 0; i < count; ++i)
        values.push_back((count == 73 ? -180. : -90.) + 5. * i);
    return values;
}
} // namespace
std::size_t angularSampleCount(const BemRequest &request) {
    if (request.observers.size() != 75 && request.observers.size() != 147)
        return 0;
    const auto &points = request.observers;
    const std::size_t count = points.size() == 147 ? 73 : 37;
    const auto subtract = [](Vec3 a, Vec3 b) { return added(a, scaled(b, -1)); };
    const auto length = [](Vec3 a) { return std::hypot(a.x, a.y, a.z); };
    double coordinateScale = 1;
    for (const auto &point : points) {
        if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z))
            return 0;
        coordinateScale =
            std::max({coordinateScale, std::abs(point.x), std::abs(point.y), std::abs(point.z)});
    }
    const auto center = count == 73 ? scaled(added(points[0], points[1]), .5)
                                    : scaled(added(points[1], points[count]), .5);
    const double measurementDistance = length(subtract(points[0], center));
    if (!std::isfinite(measurementDistance) || measurementDistance <= 0)
        return 0;
    // Relative geometric consistency plus translated-coordinate roundoff.
    const double tolerance = std::min(
        1e-6 * measurementDistance,
        1e-8 * measurementDistance + 64 * std::numeric_limits<double>::epsilon() * coordinateScale);
    if (request.farField && length(subtract(center, request.farField->origin)) > tolerance)
        return 0;
    const double directionTolerance = tolerance / measurementDistance;
    const std::size_t positive90 = count == 73 ? 54 : 36;
    const auto axis = subtract(points[0], center);
    const auto horizontal = subtract(points[1 + positive90], center);
    const auto vertical = subtract(points[1 + count + positive90], center);
    for (const auto direction : {axis, horizontal, vertical})
        if (!std::isfinite(length(direction)) || length(direction) <= 0 ||
            std::abs(length(direction) - measurementDistance) > tolerance)
            return 0;
    const auto a = scaled(axis, 1 / length(axis));
    const auto h = scaled(horizontal, 1 / length(horizontal));
    const auto v = scaled(vertical, 1 / length(vertical));
    if (std::abs(dot(a, h)) > directionTolerance || std::abs(dot(a, v)) > directionTolerance ||
        std::abs(dot(h, v)) > directionTolerance ||
        dot(cross(a, h), v) < 1 - 4 * directionTolerance)
        return 0;
    const auto angles = anglesDegrees(count);
    for (std::size_t plane = 0; plane < 2; ++plane)
        for (std::size_t i = 0; i < count; ++i) {
            const double angle = angles[i] * std::numbers::pi / 180;
            const auto expected =
                added(center, scaled(added(scaled(a, std::cos(angle)),
                                           scaled(plane == 0 ? h : v, std::sin(angle))),
                                     measurementDistance));
            if (length(subtract(points[1 + plane * count + i], expected)) > tolerance)
                return 0;
        }
    return count;
}
QByteArray resultCsv(const SimulationArchive &a) {
    if (!a.result || !a.request || a.result->frequencies.empty())
        throw std::invalid_argument("Keine berechneten Ergebnisse vorhanden.");
    if (a.request->farField) {
        const auto &farField = *a.request->farField;
        if (!std::isfinite(farField.origin.x) || !std::isfinite(farField.origin.y) ||
            !std::isfinite(farField.origin.z) || !std::isfinite(farField.referenceDistance) ||
            farField.referenceDistance <= 0)
            throw std::invalid_argument("Ungültige Fernfeld-Referenz im Rechenauftrag.");
        for (const auto point : a.request->observers) {
            const double direction =
                std::hypot(point.x - farField.origin.x, point.y - farField.origin.y,
                           point.z - farField.origin.z);
            if (!std::isfinite(direction) || direction <= 0)
                throw std::invalid_argument(
                    "Fernfeld-Richtungsanker muss vom Ursprung verschieden sein.");
        }
    }
    QByteArray bytes;
    QTextStream out(&bytes, QIODevice::WriteOnly);
    out.setLocale(QLocale::c());
    out.setRealNumberPrecision(17);
    out << "# HornSculpt; revision=" << a.revision << "; current=" << (a.current ? "true" : "false")
        << '\n';
    out << "# RMS; exp(+i*omega*t); pressure=Pa; impedance=Z*S/(rho*c); exact computed frequencies "
           "only\n";
    const bool continuous =
        a.request->axisymmetric && a.request->axisymmetric->surfaceContract ==
                                       AxisymmetricSurfaceContract::continuousRevolution;
    if (a.request->axisymmetric && !continuous &&
        a.request->axisymmetric->surfaceContract != AxisymmetricSurfaceContract::triangulatedOrbit)
        throw std::invalid_argument("Unbekannter Rotationsflächenvertrag im Ergebnis.");
    out << "# geometry_mode="
        << (a.result->geometryMode == BemGeometryMode::axisymmetric ? "axisymmetric_m0"
                                                                    : "three_dimensional")
        << "; pressure_dofs=" << static_cast<qulonglong>(a.result->pressureDofCount)
        << (continuous ? "; carrier_triangles=" : "; integration_triangles=")
        << static_cast<qulonglong>(a.result->triangleCount) << "; residual_system="
        << (a.result->geometryMode == BemGeometryMode::axisymmetric
                ? "reduced_m0_Burton_Miller"
                : "three_dimensional_Burton_Miller")
        << '\n';
    out << "# surface_contract="
        << (continuous                ? "continuousRevolution"
            : a.request->axisymmetric ? "triangulatedOrbit"
                                      : "cartesianTriangles")
        << "; surface_pressure=area_weighted_P1_mean\n";
    if (a.request->farField) {
        const auto &farField = *a.request->farField;
        out << "# field_mode=far_field; reference_distance_m=" << farField.referenceDistance
            << "; origin_m=" << farField.origin.x << ',' << farField.origin.y << ','
            << farField.origin.z
            << "; observer_coordinates=world_direction_anchors_m; "
               "pressure=Pa_at_reference_distance\n";
    } else
        out << "# field_mode=finite; observer_coordinates=world_points_m; "
               "reference_distance_m=none; "
               "pressure=Pa_at_observer\n";
    for (const auto &f : a.result->frequencies)
        out << "# interaction_plan; frequency_Hz=" << f.frequency
            << "; cpu_build_ms=" << f.interactionPlanMilliseconds
            << "; host_cache_hit=" << (f.interactionPlanReused ? "true" : "false") << '\n';
    out << "frequency_Hz;excitation;kind;id;x_m;y_m;z_m;real;imag;magnitude;phase_deg;unit\n";
    for (const auto &f : a.result->frequencies)
        for (const auto &e : f.excitations) {
            auto name = QString::fromStdString(e.name);
            name.replace('"', "\"\"");
            const auto row = [&](const char *kind, std::size_t id, Vec3 point,
                                 std::complex<double> value, const char *unit) {
                out << f.frequency << ";\"" << name << "\";" << kind << ';'
                    << static_cast<qulonglong>(id) << ';' << point.x << ';' << point.y << ';'
                    << point.z << ';' << value.real() << ';' << value.imag() << ';'
                    << std::abs(value) << ';' << std::arg(value) * 180 / std::numbers::pi << ';'
                    << unit << '\n';
            };
            if (e.observerPressure.size() != a.request->observers.size())
                throw std::invalid_argument("Messpunkte und Ergebnisdaten passen nicht zusammen.");
            for (std::size_t i = 0; i < e.observerPressure.size(); ++i)
                if (observerValid(*a.request, i))
                    row(a.request->farField ? "far_field" : "observer", i, a.request->observers[i],
                        e.observerPressure[i], "Pa");
            for (const auto &z : e.impedances)
                row("impedance", z.faceId, {},
                    normalizedVolumeImpedance(z.volumeImpedance, z.area, a.request->medium.density,
                                              a.request->medium.soundSpeed),
                    "1");
            if (e.lemSolution) {
                for (const auto &driver : e.lemSolution->drivers) {
                    const auto id = driver.id.value;
                    row("driver_current_rms", id, {}, driver.current, "A");
                    row("driver_velocity_rms", id, {}, driver.velocity, "m/s");
                    row("driver_displacement_rms", id, {}, driver.displacement, "m");
                    row("driver_electrical_power", id, {}, {driver.electricalPower, 0}, "W");
                    row("driver_copper_power", id, {}, {driver.copperPower, 0}, "W");
                    row("driver_mechanical_loss", id, {}, {driver.mechanicalPower, 0}, "W");
                    row("driver_acoustic_power", id, {}, {driver.acousticPower, 0}, "W");
                    const auto electrical =
                        std::find_if(e.lemSolution->elements.begin(), e.lemSolution->elements.end(),
                                     [&](const auto &element) { return element.id == driver.id; });
                    if (electrical != e.lemSolution->elements.end() && std::abs(driver.current) > 0)
                        row("driver_electrical_impedance", id, {},
                            electrical->effort / driver.current, "ohm");
                }
                for (const auto &port : e.lemSolution->ports) {
                    row("port_pressure", port.id.value, {}, port.pressure, "Pa");
                    row("port_volume_velocity", port.id.value, {}, port.volumeVelocity, "m3/s");
                    row("port_active_power", port.id.value, {}, {port.activePower, 0}, "W");
                }
            }
        }
    out.flush();
    if (out.status() != QTextStream::Ok)
        throw std::runtime_error("CSV-Ergebnisdaten konnten nicht erzeugt werden.");
    return bytes;
}

std::vector<CtaResultRow> ctaResults(const BemRequest &request, const BemResult &result,
                                     std::size_t excitation) {
    if (!request.farField)
        throw std::invalid_argument("CTA benötigt ein Fernfeldergebnis.");
    const std::size_t count = angularSampleCount(request);
    if (!count)
        throw std::invalid_argument("CTA benötigt ein gültiges H/V-Winkelraster.");
    const auto observers = request.observers.size();
    const int first = count == 73 ? -180 : -90;
    std::vector<CtaResultRow> rows;
    for (const auto &f : result.frequencies) {
        if (excitation >= f.excitations.size()) {
            rows.push_back({f.frequency, {}});
            continue;
        }
        const auto &p = f.excitations[excitation].observerPressure;
        if (p.size() != observers) {
            rows.push_back({f.frequency, {}});
            continue;
        }
        std::vector<AcousticsAngularSample> h, v;
        const auto level = [](std::complex<double> value) {
            return std::isfinite(value.real()) && std::isfinite(value.imag())
                       ? soundPressureLevelDecibels(value)
                       : std::numeric_limits<double>::quiet_NaN();
        };
        for (std::size_t i = 0; i < count; ++i) {
            const int angle = first + 5 * static_cast<int>(i);
            const auto hi = angle == 0 ? 0 : 1 + i, vi = angle == 0 ? 0 : 1 + count + i;
            const auto missing = std::numeric_limits<double>::quiet_NaN();
            h.push_back({double(angle), observerValid(request, hi) ? level(p[hi]) : missing});
            v.push_back({double(angle), observerValid(request, vi) ? level(p[vi]) : missing});
        }
        rows.push_back({f.frequency, cta2034Cuts(h, v)});
    }
    return rows;
}
QByteArray ctaCsv(const SimulationArchive &a, std::size_t excitation) {
    if (!a.request || !a.result)
        throw std::invalid_argument("CTA-Ergebnis fehlt.");
    const auto rows = ctaResults(*a.request, *a.result, excitation);
    if (rows.empty())
        throw std::invalid_argument("Keine CTA-Daten für diese Anregung.");
    QByteArray bytes;
    QTextStream out(&bytes, QIODevice::WriteOnly);
    out.setLocale(QLocale::c());
    out.setRealNumberPrecision(17);
    const auto frequencyWithExcitation =
        std::find_if(a.result->frequencies.begin(), a.result->frequencies.end(),
                     [excitation](const auto &f) { return excitation < f.excitations.size(); });
    if (frequencyWithExcitation == a.result->frequencies.end())
        throw std::invalid_argument("CTA-Anregung fehlt im Ergebnis.");
    auto name = QString::fromStdString(frequencyWithExcitation->excitations[excitation].name);
    name.replace('\n', ' ');
    name.replace('\r', ' ');
    out << "# HornSculpt CTA-2034 simulation; revision=" << a.revision
        << "; current=" << (a.current ? "true" : "false") << "; completeness="
        << (!a.resultComplete   ? "unknown"
            : *a.resultComplete ? "complete"
                                : "partial")
        << '\n';
    out << "# excitation=" << name << "; reference_pressure_Pa_RMS=2e-5; far_field_reference_m="
        << a.request->farField->referenceDistance << '\n';
    out << "# sound_power=70_direction_approximation; predicted_in_room=estimate; missing=empty\n";
    out << "frequency_Hz;on_axis_dB_SPL;listening_window_dB_SPL;early_reflections_dB_SPL;"
           "sound_power_dB_SPL;predicted_in_room_dB_SPL;early_reflections_DI_dB;sound_power_DI_"
           "dB\n";
    const auto value = [&](const std::optional<double> &v) {
        if (!v || std::isnan(*v))
            return;
        if (*v == -std::numeric_limits<double>::infinity())
            out << "-inf";
        else if (std::isfinite(*v))
            out << *v;
    };
    for (const auto &row : rows) {
        const auto &c = row.curves;
        out << row.frequency << ';';
        value(c.onAxisDb);
        out << ';';
        value(c.listeningWindowDb);
        out << ';';
        value(c.earlyReflectionsDb);
        out << ';';
        value(c.soundPowerDb);
        out << ';';
        value(c.predictedInRoomDb);
        out << ';';
        value(c.earlyReflectionsDiDb);
        out << ';';
        value(c.soundPowerDiDb);
        out << '\n';
    }
    out.flush();
    if (out.status() != QTextStream::Ok)
        throw std::runtime_error("CTA-CSV konnte nicht erzeugt werden.");
    return bytes;
}
} // namespace hornsculpt
