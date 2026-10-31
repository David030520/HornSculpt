#pragma once

#include <QString>
#include <vector>

namespace hornsculpt {
// Plain Hz values; units belong to the axis title or the caller's readout.
QString frequencyLabel(double frequency);
std::vector<double> frequencyTicks(double minimum, double maximum, double widthPixels);
} // namespace hornsculpt
