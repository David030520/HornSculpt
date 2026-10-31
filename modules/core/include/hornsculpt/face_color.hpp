#pragma once

#include <cstdint>
#include <map>

namespace hornsculpt {

// Display-only RGB surface color in [0,255], keyed by the owner's face identity.
struct FaceColor {
    std::uint8_t red{}, green{}, blue{};
    bool operator==(const FaceColor &) const = default;
};

using FaceColors = std::map<std::uint32_t, FaceColor>;

} // namespace hornsculpt
