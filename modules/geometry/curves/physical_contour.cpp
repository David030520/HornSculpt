#include <hornsculpt/physical_contour.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace hornsculpt {

PhysicalContour::PhysicalContour(const ContourParameters &parameters)
    : rotation_(parameters.rotation) {
    if (!std::isfinite(parameters.width) || parameters.width <= 0.0 ||
        !std::isfinite(parameters.aspectRatio) || parameters.aspectRatio <= 0.0 ||
        !std::isfinite(rotation_))
        throw std::invalid_argument("Kontur: Breite, Seitenverhaeltnis oder Rotation ungueltig.");
    double extentX = 2.0, extentY = 2.0;
    if (parameters.type == CurveType::superellipse) {
        curve_ = std::make_shared<Superellipse>(1.0, 1.0, parameters.superellipseExponent);
    } else if (parameters.type == CurveType::superformula) {
        curve_ = std::make_shared<Superformula>(parameters.superformula);
        double minX = std::numeric_limits<double>::max(), maxX = -minX;
        double minY = minX, maxY = -minX;
        constexpr std::size_t samples = 8192;
        for (std::size_t i = 0; i < samples; ++i) {
            const double phi = 2.0 * std::numbers::pi * static_cast<double>(i) / samples;
            const double r = curve_->radius(phi), x = r * std::cos(phi), y = r * std::sin(phi);
            minX = std::min(minX, x);
            maxX = std::max(maxX, x);
            minY = std::min(minY, y);
            maxY = std::max(maxY, y);
        }
        extentX = maxX - minX;
        extentY = maxY - minY;
    } else {
        throw std::invalid_argument("Kontur: unbekannter Kurventyp.");
    }
    scaleX_ = parameters.width / extentX;
    scaleY_ = parameters.width * parameters.aspectRatio / extentY;
    if (!std::isfinite(scaleX_) || !std::isfinite(scaleY_) || scaleX_ <= 0.0 || scaleY_ <= 0.0)
        throw std::invalid_argument("Kontur besitzt keine endliche positive Breite oder Hoehe.");
}

double PhysicalContour::radius(double phi) const {
    if (!std::isfinite(phi))
        throw std::invalid_argument("Konturwinkel muss endlich sein.");
    const double rotated = phi - rotation_;
    const double x = std::cos(rotated) / scaleX_, y = std::sin(rotated) / scaleY_;
    const double r = curve_->radius(std::atan2(y, x)) / std::hypot(x, y);
    if (!std::isfinite(r) || r <= 0.0)
        throw std::invalid_argument("Kontur erzeugt keinen endlichen positiven Radius.");
    return r;
}

} // namespace hornsculpt
