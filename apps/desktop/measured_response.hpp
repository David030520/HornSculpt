#pragma once

#include <QString>
#include <vector>

namespace hornsculpt {
struct MeasuredPoint {
    double frequencyHz{};
    double decibels{};
    bool operator==(const MeasuredPoint &) const = default;
};
struct MeasuredResponse {
    QString name;
    bool absoluteSpl{true};
    double offsetDb{}; // Explicit display alignment; imported samples remain unchanged.
    std::vector<MeasuredPoint> samples;
    bool operator==(const MeasuredResponse &) const = default;
};
void validateMeasuredResponse(const MeasuredResponse &response);
MeasuredResponse readMeasuredResponse(const QString &path);
} // namespace hornsculpt
