#pragma once

#include <map>

namespace hornsculpt::live {

using Trace = std::map<double, double>;

// Display estimates only: preserve baseline detail between adjacent computed
// anchors using a correction linear in log(f). Nonfinite values retain gaps.
Trace estimate(const Trace &baseline, const Trace &computed);

} // namespace hornsculpt::live
