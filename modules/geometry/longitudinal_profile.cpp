#include <hornsculpt/longitudinal_profile.hpp>

#include <algorithm>

namespace hornsculpt {

double LongitudinalProfile::slope(double z, double angle) const {
    const double step = length() * 1e-5;
    if (z < step)
        return (-3.0 * radius(z, angle) + 4.0 * radius(z + step, angle) -
                radius(z + 2.0 * step, angle)) /
               (2.0 * step);
    if (z > length() - step)
        return (3.0 * radius(z, angle) - 4.0 * radius(z - step, angle) +
                radius(z - 2.0 * step, angle)) /
               (2.0 * step);
    return (radius(z + step, angle) - radius(z - step, angle)) / (2.0 * step);
}

} // namespace hornsculpt
