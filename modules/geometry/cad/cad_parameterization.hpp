#pragma once

#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <hornsculpt/surface_parameters.hpp>
#include <memory>
#include <vector>

namespace hornsculpt::cad_detail {
class ParameterizationBuilder {
  public:
    explicit ParameterizationBuilder(double geometricTolerance);
    ~ParameterizationBuilder();
    void addFace(std::uint32_t id, const TopoDS_Face &face);
    std::uint32_t addCurve(const TopoDS_Edge &edge);
    std::uint32_t addUse(std::uint32_t curveId, const TopoDS_Edge &edge, const TopoDS_Face &face,
                         double anchorParameter, ParameterPoint anchorUV);
    std::shared_ptr<const SurfaceParameterization> finish(SurfaceMesh mesh,
                                                          std::vector<ParameterTriangle> triangles);

  private:
    struct Data;
    std::unique_ptr<Data> data_;
};
} // namespace hornsculpt::cad_detail
