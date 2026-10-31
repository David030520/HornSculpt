#pragma once
#include <cstdint>
#include <hornsculpt/acoustic_domains.hpp>
#include <hornsculpt/mesh.hpp>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace hornsculpt {
// A pressure unknown belongs to a (domain, vertex) pair. Only vertices on a
// declared transparent interface identify two domains; physical Wall sides do
// not. The geometry/ownership validator remains a separate prerequisite.
struct AcousticPressureDofMap {
    std::vector<std::map<std::uint32_t, std::uint32_t>> domains;
    std::uint32_t count{};
};
inline AcousticPressureDofMap makeAcousticPressureDofMap(const SurfaceMesh &mesh,
                                                         const BemSubdomainModel &model) {
    AcousticPressureDofMap result;
    result.domains.resize(model.domains.size());
    std::map<std::uint32_t, std::size_t> domainIndex;
    std::vector<unsigned> parent;
    for (std::size_t d = 0; d < model.domains.size(); ++d) {
        if (!domainIndex.emplace(model.domains[d].id, d).second)
            throw std::invalid_argument("Duplicate acoustic domain ID.");
        for (const auto &side : model.domains[d].panels) {
            if (side.triangleIndex >= mesh.triangles.size())
                throw std::invalid_argument("Invalid acoustic domain panel index.");
            for (auto vertex : mesh.triangles[side.triangleIndex].vertices) {
                if (vertex >= mesh.vertices.size())
                    throw std::invalid_argument("Invalid acoustic vertex index.");
                if (!result.domains[d].contains(vertex)) {
                    if (parent.size() >= std::numeric_limits<std::uint32_t>::max())
                        throw std::invalid_argument("Too many acoustic pressure nodes.");
                    const auto index = static_cast<std::uint32_t>(parent.size());
                    result.domains[d].emplace(vertex, index);
                    parent.push_back(index);
                }
            }
        }
    }
    const auto findRoot = [&](unsigned node) {
        while (parent[node] != node) {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        return node;
    };
    for (const auto &interface : model.interfaces) {
        if (!domainIndex.contains(interface.firstDomainId) ||
            !domainIndex.contains(interface.secondDomainId) ||
            interface.firstDomainId == interface.secondDomainId)
            throw std::invalid_argument("Invalid acoustic interface domains.");
        const auto a = domainIndex.at(interface.firstDomainId),
                   b = domainIndex.at(interface.secondDomainId);
        for (auto triangle : interface.triangleIndices) {
            if (triangle >= mesh.triangles.size())
                throw std::invalid_argument("Invalid acoustic interface panel.");
            for (auto vertex : mesh.triangles[triangle].vertices) {
                if (!result.domains[a].contains(vertex) || !result.domains[b].contains(vertex))
                    throw std::invalid_argument("Interface vertex is absent from a domain.");
                const auto x = findRoot(result.domains[a].at(vertex)),
                           y = findRoot(result.domains[b].at(vertex));
                if (x != y)
                    parent[y] = x;
            }
        }
    }
    std::map<unsigned, std::uint32_t> compact;
    for (auto &domain : result.domains)
        for (auto &[vertex, index] : domain) {
            (void)vertex;
            const auto [it, inserted] =
                compact.emplace(findRoot(index), static_cast<std::uint32_t>(compact.size()));
            (void)inserted;
            index = it->second;
        }
    result.count = static_cast<std::uint32_t>(compact.size());
    return result;
}
} // namespace hornsculpt
