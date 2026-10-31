#pragma once

#include <QJsonArray>
#include <QJsonObject>
#include <QString>
#include <cstdint>
#include <hornsculpt/bem.hpp>

namespace hornsculpt {
struct Project;
namespace mcp {

struct StudyRequest {
    BemRequest bem;
    QJsonArray observerMetadata;
    QString referenceExcitation;
    std::uint32_t referenceFaceId{};
};

// Requires the current, qualified project mesh. The CUDA solver additionally
// rejects finite observers in solid/on the boundary and unconverged near fields.
StudyRequest parseStudyRequest(const QJsonObject &input, const Project &project);
QJsonObject describeStudyCapabilities();

// Preserves complex RMS quantities. H is p/v of the explicitly selected physical
// Driven face, not an electrical transfer function or an impulse response.
QJsonObject complexResultDocument(const StudyRequest &request, const BemResult &result);

} // namespace mcp
} // namespace hornsculpt
