#include "fixture_geometry.hpp"
#include "galerkin/exterior_domain.hpp"
#include <hornsculpt/mesh_orientation.hpp>
#include <hornsculpt/meshing.hpp>
#include <iostream>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class Function> void rejects(Function function) {
    bool rejected = false;
    try {
        function();
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "Invalid orientation was accepted.");
}
} // namespace

int main() {
    try {
        auto source = testing::sphere(0);
        source.fixedVertices = {0, 1};
        std::set<std::uint32_t> all;
        for (std::size_t i = 0; i < source.triangles.size(); ++i) {
            auto &triangle = source.triangles[i];
            triangle.faceId = static_cast<std::uint32_t>(i + 1);
            triangle.boundary = i == 0 ? BoundaryCondition::driven : BoundaryCondition::wall;
            all.insert(triangle.faceId);
        }
        const auto original = source;
        require(applyFaceOrientation(source, {}) == source, "Empty mask changed the source.");
        const auto selected = applyFaceOrientation(source, {1});
        require(source == original && selected.vertices == source.vertices &&
                    selected.fixedVertices == source.fixedVertices,
                "Orientation changed coordinates, fixed vertices or the immutable source.");
        for (std::size_t i = 0; i < source.triangles.size(); ++i) {
            const auto &a = source.triangles[i], &b = selected.triangles[i];
            require(a.faceId == b.faceId && a.boundary == b.boundary,
                    "Orientation changed physical tags or boundary conditions.");
            auto expected = a.vertices;
            if (a.faceId == 1)
                std::swap(expected[1], expected[2]);
            require(b.vertices == expected, "Wrong physical face was reversed.");
        }
        require(inspectMesh(selected).inconsistentEdges == 3 &&
                    !inspectMesh(selected).watertight(),
                "Single-face orientation conflict was hidden.");
        require(applyFaceOrientation(selected, {1}) == original,
                "Reversing twice did not restore the exact source.");
        const auto inward = applyFaceOrientation(source, all);
        require(inspectMesh(inward).watertight() && inspectMesh(inward).signedVolume < 0,
                "Complete reversal did not preserve topology and reverse volume.");
        rejects([&] { galerkin::validateExteriorComponents(inward); });
        rejects([&] { (void)applyFaceOrientation(source, {9999}); });
        require(source == original, "Rejected reversal modified the source.");
        galerkin::validateExteriorComponents(source);
        std::cout << "Physical face orientation, reversible masks and exterior gates PASS.\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
