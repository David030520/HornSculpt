#include "spherical_rotation.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace hornsculpt::wideband {
namespace {
constexpr unsigned maximumOrder = maximumDeviceScaledOrder;
void validate(const SphericalRotationPlan &plan, std::span<const Complex> values) {
    if (plan.order > maximumOrder ||
        values.size() != static_cast<std::size_t>(plan.order + 1) * (plan.order + 1) ||
        plan.degreeOffsets.size() != plan.order + 2 || plan.degreeOffsets.front() != 0 ||
        plan.degreeOffsets.back() != plan.smallD.size() ||
        plan.azimuthPhases.size() != 2 * plan.order + 1)
        throw std::invalid_argument("Spherical rotation dimensions are inconsistent.");
    for (unsigned l = 0; l <= plan.order; ++l)
        if (plan.degreeOffsets[l] > plan.degreeOffsets[l + 1] ||
            plan.degreeOffsets[l + 1] - plan.degreeOffsets[l] !=
                static_cast<std::size_t>(2 * l + 1) * (2 * l + 1))
            throw std::invalid_argument("Spherical rotation degree block is invalid.");
    for (const auto value : values)
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::invalid_argument("Spherical rotation source is non-finite.");
    for (const auto value : plan.smallD)
        if (!std::isfinite(value))
            throw std::invalid_argument("Spherical rotation coefficient is non-finite.");
    for (const auto value : plan.azimuthPhases)
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::invalid_argument("Spherical rotation azimuth phase is non-finite.");
}
} // namespace
SphericalRotationPlan makeSphericalRotationPlan(unsigned order, Vec3 axis) {
    if (order > maximumOrder)
        throw std::invalid_argument("Spherical rotation order exceeds the supported contract.");
    const double radius = std::hypot(axis.x, axis.y, axis.z);
    if (!std::isfinite(radius) || radius <= 0)
        throw std::invalid_argument("Spherical rotation axis must be finite and nonzero.");
    const double theta = std::atan2(std::hypot(axis.x, axis.y), axis.z),
                 phi = std::atan2(axis.y, axis.x), c = std::cos(theta / 2), s = std::sin(theta / 2);
    SphericalRotationPlan result;
    result.order = order;
    result.degreeOffsets = {0, 1};
    result.smallD = {1};
    for (int m = -static_cast<int>(order); m <= static_cast<int>(order); ++m)
        result.azimuthPhases.emplace_back(std::cos(m * phi), std::sin(m * phi));
    std::vector<double> previous{1};
    // With q=2j, the two maximal-spin coupling factors for row i are
    // sqrt(i/q) and sqrt((q-i)/q). Their columns form an isometry. Thus
    // D^j=C^T*(D^(j-1/2) tensor d^(1/2))*C cannot amplify prior matrix error
    // in exact arithmetic. Every elementary coefficient has magnitude <=1.
    for (unsigned q = 1; q <= 2 * order; ++q) {
        const auto dimension = static_cast<std::size_t>(q) + 1;
        std::vector<double> current(dimension * dimension), up(dimension), down(dimension);
        for (unsigned i = 0; i <= q; ++i) {
            up[i] = std::sqrt(static_cast<double>(i) / q);
            down[i] = std::sqrt(static_cast<double>(q - i) / q);
        }
        for (unsigned i = 0; i <= q; ++i)
            for (unsigned j = 0; j <= q; ++j) {
                double value = 0;
                if (i && j)
                    value +=
                        up[i] * up[j] * c * previous[(static_cast<std::size_t>(i) - 1) * q + j - 1];
                if (i && j < q)
                    value -=
                        up[i] * down[j] * s * previous[(static_cast<std::size_t>(i) - 1) * q + j];
                if (i < q && j)
                    value +=
                        down[i] * up[j] * s * previous[static_cast<std::size_t>(i) * q + j - 1];
                if (i < q && j < q)
                    value += down[i] * down[j] * c * previous[static_cast<std::size_t>(i) * q + j];
                current[static_cast<std::size_t>(i) * dimension + j] = value;
            }
        if (!(q % 2)) {
            result.smallD.insert(result.smallD.end(), current.begin(), current.end());
            result.degreeOffsets.push_back(result.smallD.size());
        }
        previous = std::move(current);
    }
    return result;
}
std::vector<Complex> rotateToAxis(const SphericalRotationPlan &plan,
                                  std::span<const Complex> values) {
    validate(plan, values);
    std::vector<Complex> result(values.size());
    for (unsigned l = 0; l <= plan.order; ++l) {
        const auto count = 2 * l + 1;
        const auto offset = static_cast<std::size_t>(l) * l;
        for (unsigned target = 0; target < count; ++target)
            for (unsigned source = 0; source < count; ++source)
                result[offset + target] +=
                    plan.smallD[plan.degreeOffsets[l] + static_cast<std::size_t>(source) * count +
                                target] *
                    plan.azimuthPhases[plan.order + source - l] * values[offset + source];
    }
    return result;
}
std::vector<Complex> rotateFromAxis(const SphericalRotationPlan &plan,
                                    std::span<const Complex> values) {
    validate(plan, values);
    std::vector<Complex> result(values.size());
    for (unsigned l = 0; l <= plan.order; ++l) {
        const auto count = 2 * l + 1;
        const auto offset = static_cast<std::size_t>(l) * l;
        for (unsigned target = 0; target < count; ++target) {
            for (unsigned source = 0; source < count; ++source)
                result[offset + target] +=
                    plan.smallD[plan.degreeOffsets[l] + static_cast<std::size_t>(target) * count +
                                source] *
                    values[offset + source];
            result[offset + target] *= std::conj(plan.azimuthPhases[plan.order + target - l]);
        }
    }
    return result;
}
} // namespace hornsculpt::wideband
