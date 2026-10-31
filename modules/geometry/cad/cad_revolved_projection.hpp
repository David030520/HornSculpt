#pragma once

#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <array>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <hornsculpt/mesh.hpp>
#include <optional>

namespace hornsculpt::cad_detail {
struct SurfaceMinimum {
    gp_Pnt position;
    double u{}, v{};
};

// A full revolution of a certified planar, one-sided spline meridian admits
// an exact closest-point reduction to that curve's radius/axial plane.
class RevolvedProjection {
  public:
    static std::optional<RevolvedProjection> make(const Handle(Geom_Surface) & surface,
                                                  std::array<double, 4> bounds);
    // Scalar-distance callers may accept any minimizer; mappings require one.
    SurfaceMinimum minimum(Vec3 approximate, bool requireUnique = true) const;

  private:
    Handle(Geom_SurfaceOfRevolution) surface_;
    Handle(Geom_Curve) curve_;
    gp_Pnt origin_;
    gp_Vec axis_, radial_, transverse_;
    std::array<double, 4> bounds_{};
    double radialScale_{};
};
} // namespace hornsculpt::cad_detail
