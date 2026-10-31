#include "measured_response.hpp"
#include <QFile>
#include <QFileInfo>
#include <QRegularExpression>
#include <cmath>
#include <stdexcept>

namespace hornsculpt {
void validateMeasuredResponse(const MeasuredResponse &r) {
    if (r.name.isEmpty() || r.name.size() > 4096 || !std::isfinite(r.offsetDb) ||
        std::abs(r.offsetDb) > 300 || r.samples.size() < 2 || r.samples.size() > 1000000)
        throw std::invalid_argument("Ungültige Messkurve.");
    double previous = 0;
    for (const auto &p : r.samples) {
        if (!std::isfinite(p.frequencyHz) || p.frequencyHz <= previous ||
            !std::isfinite(p.decibels))
            throw std::invalid_argument(
                "Messfrequenzen müssen positiv und streng aufsteigend sein.");
        previous = p.frequencyHz;
    }
}
MeasuredResponse readMeasuredResponse(const QString &path) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly))
        throw std::runtime_error(file.errorString().toStdString());
    if (file.size() > 32 * 1024 * 1024)
        throw std::invalid_argument("Messdatei ist größer als 32 MiB.");
    MeasuredResponse result;
    result.name = QFileInfo(path).completeBaseName();
    qsizetype lineNumber = 0;
    for (auto line : QString::fromUtf8(file.readAll()).split('\n')) {
        ++lineNumber;
        line = line.trimmed();
        if (line.isEmpty() || line.startsWith('#') || line.startsWith('*') || line.startsWith("//"))
            continue;
        QStringList columns;
        if (line.contains(';')) {
            columns = line.split(';');
            for (auto &column : columns)
                column.replace(',', '.');
        } else if (line.contains(','))
            columns = line.split(',');
        else
            columns = line.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
        if (columns.size() < 2 || columns.size() > 3)
            throw std::invalid_argument(
                QString("Messdatei, Zeile %1: erwartet Hz, dB, optional Phase.")
                    .arg(lineNumber)
                    .toStdString());
        bool fOk = false, dbOk = false;
        const double frequency = columns[0].trimmed().toDouble(&fOk);
        const double level = columns[1].trimmed().toDouble(&dbOk);
        if (!fOk || !dbOk) {
            // One explicit unit header is accepted; arbitrary malformed rows are errors.
            const auto header = columns[0].trimmed().toLower();
            if (result.samples.empty() &&
                (header == "frequency_hz" || header == "frequency" || header == "frequenz" ||
                 header == "frequenz_hz" || header == "hz") &&
                columns[1].contains("db", Qt::CaseInsensitive))
                continue;
            throw std::invalid_argument(
                QString("Messdatei, Zeile %1: ungültige Zahl.").arg(lineNumber).toStdString());
        }
        if (columns.size() == 3) {
            bool phaseOk = false;
            const double phase = columns[2].trimmed().toDouble(&phaseOk);
            if (!phaseOk || !std::isfinite(phase))
                throw std::invalid_argument("Ungültige Phase in der Messdatei.");
        }
        result.samples.push_back({frequency, level});
    }
    if (file.error() != QFileDevice::NoError)
        throw std::runtime_error(file.errorString().toStdString());
    validateMeasuredResponse(result);
    return result;
}
} // namespace hornsculpt
