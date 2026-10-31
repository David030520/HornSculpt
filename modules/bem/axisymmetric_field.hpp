#pragma once
#include <hornsculpt/field_evaluation.hpp>

namespace hornsculpt::bem_detail {
std::vector<PressureFieldSample>
evaluateAxisymmetricPressureField(const SurfaceMesh &mesh, const BemRequest &request,
                                  double frequency, const BemExcitationResult &solution,
                                  const std::vector<Vec3> &points, PressureFieldSettings settings,
                                  std::stop_token stop,
                                  std::function<void(std::size_t, std::size_t)> progress);
} // namespace hornsculpt::bem_detail
