#pragma once

#include "project_file.hpp"
#include <QString>
#include <array>
#include <span>
#include <vector>

namespace hornsculpt {
struct SweepParameter {
    QString key, label, unit;
    double displayScale{1};
};
std::span<const SweepParameter> sweepParameters();
// The source's asymmetry determines which directional axes are meaningful.
std::vector<SweepParameter> sweepParameters(const DraftParameters &baseline);
double sweepParameterValue(const DraftParameters &parameters, const QString &key);
void setSweepParameterValue(DraftParameters &parameters, const QString &key, double value);

struct SweepAxis {
    QString parameter;
    double start{}, end{}; // SI units; angles in radians.
    std::size_t count{1};
};
enum class SweepState { pending, running, completed, failed, cancelled };
struct SweepVariant {
    std::vector<double> values;
    SweepState state{SweepState::pending};
    QString projectFile, error;
};
struct SweepStudy {
    QString directory;
    Project baseline;
    std::vector<SweepAxis> axes;
    std::vector<SweepVariant> variants;
    std::vector<double> frequencyOverride;
};

// Shared global bounds apply to both source-specific local target maps.
void validateSweepMeshSettings(const Project &project);
std::vector<SweepVariant> makeSweepGrid(std::span<const SweepAxis> axes);
Project sweepVariantInputs(const SweepStudy &study, std::size_t index);
// Checks immutable study inputs, not the independently computed mesh or result.
bool sweepVariantMatches(const SweepStudy &study, std::size_t index, const Project &project);
SweepStudy createSweepStudy(const QString &parentDirectory, const Project &baseline,
                            std::vector<SweepAxis> axes, std::vector<double> frequencies = {});
void writeSweepStudy(const SweepStudy &study);
SweepStudy readSweepStudy(const QString &manifestPath);
QString sweepManifestPath(const SweepStudy &study);
QString sweepVariantLabel(const SweepStudy &study, std::size_t index);
} // namespace hornsculpt
