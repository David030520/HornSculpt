#pragma once
#include "project_file.hpp"
#include <QByteArray>
#include <hornsculpt/acoustic_analysis.hpp>
namespace hornsculpt {
QByteArray resultCsv(const SimulationArchive &archive);
std::size_t angularSampleCount(const BemRequest &request);
inline std::size_t angularSampleCount(const std::optional<BemRequest> &request) {
    return request ? angularSampleCount(*request) : 0;
}
struct CtaResultRow { double frequency; Cta2034Cuts curves; };
std::vector<CtaResultRow> ctaResults(const BemRequest &request, const BemResult &result, std::size_t excitation);
QByteArray ctaCsv(const SimulationArchive &archive, std::size_t excitation);
}
