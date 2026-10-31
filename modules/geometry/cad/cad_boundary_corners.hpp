#pragma once

#include <TopoDS_Face.hxx>
#include <hornsculpt/surface_geometry.hpp>
#include <vector>

namespace hornsculpt::cad_detail {
std::vector<BoundaryCorner> boundaryCorners(const TopoDS_Face &face, std::uint32_t faceId);
}
