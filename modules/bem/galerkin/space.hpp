#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <hornsculpt/mesh.hpp>
#include <vector>

namespace hornsculpt::galerkin {

struct Panel {
    std::array<std::uint32_t, 3> nodes{};
    std::array<Vec3, 3> vertices{};
    Vec3 centroid;
    Vec3 normal;
    double area{};
    double maximumEdgeLength{};
    // Piecewise affine scalar gradients; scalar surface curl is n cross gradient.
    std::array<Vec3, 3> gradients{};
    std::array<Vec3, 3> curls{};
    // Row-major local P1/P1 mass; mixed mass integrates each P1 against panel DP0=1.
    std::array<double, 9> massP1P1{};
    std::array<double, 3> massP1DP0{};
    std::uint32_t faceId{};
    BoundaryCondition boundary{BoundaryCondition::wall};
};

struct BasisIncidence {
    std::uint32_t panel{};
    std::uint8_t local{};
};

struct Space {
    // Pressure DOF i corresponds exactly to mesh vertex i; flux DOF t to triangle t.
    std::vector<Vec3> vertices;
    std::vector<Panel> panels;
    std::vector<std::size_t> incidenceOffsets;
    std::vector<BasisIncidence> incidences;
    std::vector<double> massDiagonal;
    std::vector<double> lumpedMass;
    double surfaceArea{};
};

// Frequency-independent cache. Rejects invalid coordinates, indices, degenerate
// panels and orphan vertices. Closed oriented topology is checked by the caller.
Space buildSpace(const SurfaceMesh &mesh);

} // namespace hornsculpt::galerkin
