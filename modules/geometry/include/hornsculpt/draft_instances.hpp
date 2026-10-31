#pragma once

#include <array>
#include <hornsculpt/mesh.hpp>
#include <string>
#include <vector>

namespace hornsculpt {

// Instances share one Draft definition. The rigid frame maps Draft XYZ into
// world coordinates. A mirrored instance reflects local X before placement.
// mouthAnchored retains the mouth location when the shared length changes.
struct DraftInstance {
    std::uint32_t id{1};
    std::string name{"Horn 1"};
    Vec3 position{};
    std::array<Vec3, 3> axes{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    bool mirrorHorizontal{};
    bool mouthAnchored{};
    bool operator==(const DraftInstance &) const = default;
};

struct DraftAssembly {
    std::vector<DraftInstance> instances{DraftInstance{}};
    bool includeCad{};
    double cadTargetEdgeLength{.2}; // m, independent scalar CAD size field.
    bool operator==(const DraftAssembly &) const = default;
};

// IDs survive reordering and removal. Draft face IDs occupy eight slots per
// instance. CAD IDs use a disjoint namespace in the composed solver mesh.
inline constexpr std::uint32_t assemblyCadFaceOffset = 0x40000000;
void validateDraftAssembly(const DraftAssembly &assembly);
bool isPlacedAssembly(const DraftAssembly &assembly);
std::uint32_t instanceFaceId(std::uint32_t instance, std::uint32_t localFace);
Vec3 instancePoint(const DraftInstance &instance, Vec3 local, double hornLength);
SurfaceMesh composeDraftInstances(const SurfaceMesh &draft, const DraftAssembly &assembly,
                                  double hornLength, const SurfaceMesh *cad = nullptr);
} // namespace hornsculpt
