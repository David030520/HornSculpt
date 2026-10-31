#pragma once

#include <cstdint>
#include <filesystem>
#include <hornsculpt/mesh.hpp>
#include <istream>
#include <stop_token>
#include <string>
#include <vector>

namespace hornsculpt {
struct GmshImportOptions {
    // Mandatory conversion from file coordinates to SI metres.
    double lengthScaleToMetres{};
    std::stop_token stop{};
};
struct GmshPhysicalGroup {
    int dimension{};
    std::uint32_t tag{};
    std::string name;
};
struct GmshTriangleMetadata {
    std::uint64_t elementTag{};
    std::uint32_t entityTag{};
    std::vector<std::uint32_t> physicalTags;
};
struct GmshTopologyDiagnostics {
    std::size_t skippedPoints{};
    std::size_t skippedLines{};
    std::size_t boundaryEdges{};
    std::size_t nonManifoldEdges{};
    std::size_t inconsistentEdges{};
    std::size_t degenerateTriangles{};
};
struct GmshSurfaceImport {
    SurfaceMesh mesh;
    std::vector<std::uint64_t> nodeTags; // Original node tags in mesh vertex order.
    std::vector<GmshPhysicalGroup> physicalGroups;
    std::vector<GmshTriangleMetadata> triangleMetadata;
    GmshTopologyDiagnostics diagnostics;
};
// ASCII MSH 2.2 and 4.1; linear triangles and planar convex quadrilaterals.
// Points and lines are checked and counted as auxiliary elements. Reject unsupported
// sections and elements, preserve triangle winding; do not repair topology.
GmshSurfaceImport importGmsh(std::istream &input, GmshImportOptions options);
GmshSurfaceImport importGmshFile(const std::filesystem::path &path, GmshImportOptions options);
} // namespace hornsculpt
