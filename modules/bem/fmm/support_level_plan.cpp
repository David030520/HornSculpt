#include "support_level_plan.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hornsculpt::wideband {
namespace {
void requireSphericalDomain(const WidebandLevelPlan &plan) {
    if (plan.electricalRadius > maximumScaledElectricalRadius ||
        plan.signatureOrder > maximumScaledSphericalOrder)
        throw std::runtime_error(
            "Conditioned support LF plan exceeds the supported spherical domain.");
}
} // namespace
double enforceSupportDiagonalStability(WidebandLevelPlan &plan, double k, double distance) {
    if (plan.representation != LevelRepresentation::planeWave || distance == 0)
        return 0;
    if (!std::isfinite(distance) || distance <= 0 || !std::isfinite(k) || k <= 0)
        throw std::invalid_argument(
            "Support diagonal conditioning requires a positive physical distance.");
    const double z = k * distance;
    if (!std::isfinite(z) || !(z > 0))
        throw std::invalid_argument("Support diagonal electrical distance is invalid.");
    Complex previous{std::sin(z) / z, std::cos(z) / z}, current = previous * Complex{1 / z, 1};
    double condition = z * std::abs(previous);
    for (unsigned l = 1; l <= plan.translationOrder; ++l) {
        condition += z * (2 * l + 1) * std::abs(current);
        if (!std::isfinite(condition))
            break;
        const auto next = (static_cast<double>(2 * l + 1) / z) * current - previous;
        previous = current;
        current = next;
    }
    if (!std::isfinite(condition) ||
        std::numeric_limits<double>::epsilon() * condition > .01 * plan.accuracySeed) {
        requireSphericalDomain(plan);
        plan.representation = LevelRepresentation::regularSpherical;
    }
    return condition;
}
SupportExpansionPlan makeSupportExpansionPlan(const SupportTree &tree,
                                              std::span<const std::uint32_t> faceColumns, double k,
                                              double seed, std::stop_token stop) {
    const auto cancelled = [&] {
        if (stop.stop_requested())
            throw std::runtime_error("BEM calculation cancelled.");
    };
    cancelled();
    if (faceColumns.size() != tree.fluxLeaf.size() || !std::isfinite(k) || k <= 0 ||
        !std::isfinite(seed) || seed <= 0 || seed >= 1)
        throw std::invalid_argument("Support expansion geometry or accuracy seed is invalid.");
    constexpr double infinity = std::numeric_limits<double>::infinity();
    SupportExpansionPlan result;
    result.levels.resize(tree.levels.size());
    result.activeChannels.resize(tree.nodes.size());
    std::vector<bool> hasFlux(tree.nodes.size());
    std::vector<std::array<double, 3>> delta(tree.nodes.size(), {infinity, infinity, infinity});
    for (std::size_t panel = 0; panel < faceColumns.size(); ++panel) {
        if (faceColumns[panel] == UINT32_MAX)
            continue;
        auto node = tree.fluxLeaf[panel];
        while (node != noSupportNode) {
            if (node >= tree.nodes.size())
                throw std::invalid_argument("Support expansion flux leaf is invalid.");
            if (hasFlux[node])
                break;
            hasFlux[node] = true;
            node = tree.nodes[node].parent;
        }
    }
    const auto pair = [&](SupportInteraction interaction, bool flux) {
        cancelled();
        if (interaction.source >= tree.nodes.size() || interaction.target >= tree.nodes.size())
            throw std::invalid_argument("Support expansion far endpoint is invalid.");
        if (flux && !hasFlux[interaction.source])
            return;
        const auto &source = tree.nodes[interaction.source],
                   &target = tree.nodes[interaction.target];
        const double rs = flux ? source.fluxRadius : source.pressureRadius,
                     rt = target.pressureRadius,
                     distance = std::hypot(target.center.x - source.center.x,
                                           target.center.y - source.center.y,
                                           target.center.z - source.center.z);
        if (source.level != target.level || source.level >= tree.levels.size() || !(rs > 0) ||
            !(rt > 0) || !std::isfinite(distance) || !(distance > rs + rt))
            throw std::invalid_argument(
                "Support expansion requires separated equal-level far balls.");
        const auto channel = flux ? 1u : 0u;
        delta[interaction.source][channel] =
            std::min(delta[interaction.source][channel], (distance - rt) / rs);
        delta[interaction.target][2] = std::min(delta[interaction.target][2], (distance - rs) / rt);
        auto &minimum = result.levels[source.level].minimumFarDistance;
        minimum = minimum > 0 ? std::min(minimum, distance) : distance;
    };
    for (const auto interaction : tree.pressurePressure.farInteractions)
        pair(interaction, false);
    for (const auto interaction : tree.pressureFlux.farInteractions)
        pair(interaction, true);
    for (std::size_t l = 0; l < tree.levels.size(); ++l) {
        cancelled();
        auto &level = result.levels[l];
        level.sourceSeparation = level.targetSeparation = infinity;
        for (const auto id : tree.levels[l].nodes) {
            if (id >= tree.nodes.size() || tree.nodes[id].level != l)
                throw std::invalid_argument("Support expansion level inventory is invalid.");
            const auto &node = tree.nodes[id];
            if (node.parent != noSupportNode) {
                if (node.parent >= tree.nodes.size() || tree.nodes[node.parent].level >= l)
                    throw std::invalid_argument("Support expansion parent ordering is invalid.");
                const auto &parent = tree.nodes[node.parent];
                const double shift =
                    std::hypot(node.center.x - parent.center.x, node.center.y - parent.center.y,
                               node.center.z - parent.center.z);
                const auto nested = [&](double childRadius, double parentRadius) {
                    return shift + childRadius <=
                           parentRadius * (1 + 32 * std::numeric_limits<double>::epsilon());
                };
                if ((node.pressureCount && !nested(node.pressureRadius, parent.pressureRadius)) ||
                    (node.fluxCount && hasFlux[id] && !nested(node.fluxRadius, parent.fluxRadius)))
                    throw std::invalid_argument(
                        "Support expansion requires nested physical support balls.");
                for (unsigned channel = 0; channel < 3; ++channel)
                    delta[id][channel] = std::min(delta[id][channel], delta[node.parent][channel]);
            }
            if (node.pressureCount) {
                if (std::isfinite(delta[id][0])) {
                    result.activeChannels[id] |= pressureOutgoing;
                    level.sourceSeparation = std::min(level.sourceSeparation, delta[id][0]);
                }
                if (std::isfinite(delta[id][2])) {
                    result.activeChannels[id] |= pressureIncoming;
                    level.targetSeparation = std::min(level.targetSeparation, delta[id][2]);
                }
                if (result.activeChannels[id] & (pressureOutgoing | pressureIncoming))
                    level.maximumActiveRadius =
                        std::max(level.maximumActiveRadius, node.pressureRadius);
            }
            if (node.fluxCount && hasFlux[id] && std::isfinite(delta[id][1])) {
                result.activeChannels[id] |= fluxOutgoing;
                level.sourceSeparation = std::min(level.sourceSeparation, delta[id][1]);
                level.maximumActiveRadius = std::max(level.maximumActiveRadius, node.fluxRadius);
            }
        }
    }
    // Different far-root branches can make the raw level maxima nonmonotone.
    // Common level representations use a conservative monotone envelope; no
    // inactive ancestors are activated merely to enforce this envelope.
    double radius = 0;
    for (auto i = result.levels.rbegin(); i != result.levels.rend(); ++i)
        if (i->maximumActiveRadius > 0) {
            radius = std::max(radius, i->maximumActiveRadius);
            i->expansion.activeRadius = radius;
        }
    double policyDelta = infinity;
    bool spherical = false;
    for (auto &level : result.levels) {
        if (!(level.maximumActiveRadius > 0)) {
            level.sourceSeparation = level.targetSeparation = 0;
            continue;
        }
        policyDelta = std::min({policyDelta, level.sourceSeparation, level.targetSeparation});
        if (!std::isfinite(policyDelta) || !(policyDelta > 1))
            throw std::invalid_argument("Support expansion has no finite convergence geometry.");
        level.policySeparation = policyDelta;
        level.expansion = makeWidebandLevelPlan(level.expansion.activeRadius, k, policyDelta, seed);
        level.diagonalConditionIndicator =
            enforceSupportDiagonalStability(level.expansion, k, level.minimumFarDistance);
        if (spherical) {
            requireSphericalDomain(level.expansion);
            level.expansion.representation = LevelRepresentation::regularSpherical;
        }
        spherical |= level.expansion.representation == LevelRepresentation::regularSpherical;
        if (!std::isfinite(level.sourceSeparation))
            level.sourceSeparation = 0;
        if (!std::isfinite(level.targetSeparation))
            level.targetSeparation = 0;
    }
    return result;
}
} // namespace hornsculpt::wideband
