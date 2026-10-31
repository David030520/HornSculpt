#pragma once
#include "gmsh_import_data.hpp"
#include <QJsonObject>
#include <hornsculpt/acoustic_partition.hpp>
#include <hornsculpt/freecad_link.hpp>
#include <memory>

namespace hornsculpt {
std::shared_ptr<const CadSnapshot> makeGmshCadSnapshot(const MeshImportSelection &selection,
                                                       const QString &filename);
MeshingResult validateGmshMesh(const SurfaceMesh &mesh,
                               std::span<const AcousticRegionDefinition> regions,
                               const MeshingOptions &options, std::stop_token stop = {});
QJsonObject encodeGmshMetadata(const MeshImportSelection &selection);
MeshImportSelection decodeGmshMetadata(const QJsonObject &json, const SurfaceMesh &mesh);
} // namespace hornsculpt
