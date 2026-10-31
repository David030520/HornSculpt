#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace hornsculpt {

enum class BemDomainKind { interior, exterior };

// normalSign multiplies the mesh winding normal to obtain the normal pointing
// out of this domain's air. It must be exactly +1 or -1.
struct BemDomainPanel {
    std::uint32_t triangleIndex{};
    int normalSign{1};
    bool operator==(const BemDomainPanel &) const = default;
};
struct BemSubdomain {
    std::uint32_t id{};
    std::string name;
    BemDomainKind kind{BemDomainKind::interior};
    std::vector<BemDomainPanel> panels;
    bool operator==(const BemSubdomain &) const = default;
};
struct BemInterface {
    std::uint32_t id{};
    std::string name;
    std::uint32_t firstDomainId{};
    std::uint32_t secondDomainId{};
    std::vector<std::uint32_t> triangleIndices;
    bool operator==(const BemInterface &) const = default;
};
struct BemSubdomainModel {
    std::vector<BemSubdomain> domains;
    std::vector<BemInterface> interfaces;
    bool operator==(const BemSubdomainModel &) const = default;
};

} // namespace hornsculpt
