#pragma once

#include "project_file.hpp"
#include <QJsonArray>
#include <QJsonObject>
#include <QJsonValue>
#include <complex>
#include <cstdint>
#include <initializer_list>
#include <map>

namespace hornsculpt::project_json {
[[noreturn]] void invalid(const QString &message);
void keys(const QJsonObject &object, std::initializer_list<const char *> expected);
QJsonObject object(const QJsonValue &value);
QJsonArray array(const QJsonValue &value, qsizetype maximum, qsizetype minimum = 0);
double number(const QJsonValue &value);
double positive(const QJsonValue &value);
bool boolean(const QJsonValue &value);
QString text(const QJsonValue &value, bool allowEmpty = false);
std::uint32_t integer(const QJsonValue &value, std::uint32_t maximum = UINT32_MAX);
std::uint64_t wideInteger(const QJsonValue &value);
QString encodeWide(std::uint64_t value);
Vec3 vector(const QJsonValue &value);
QJsonArray encodeVector(Vec3 value);
std::complex<double> complex(const QJsonValue &value);
QJsonArray encodeComplex(std::complex<double> value);
bool sameMesh(const SurfaceMesh &a, const SurfaceMesh &b);
std::map<std::uint32_t, BoundaryCondition> meshFaces(const SurfaceMesh &mesh);
QJsonObject encodeMesh(const SurfaceMesh &mesh);
SurfaceMesh decodeMesh(const QJsonValue &value, bool closed);
QJsonObject encodeDraft(const Project &project);
QJsonObject encodeAssembly(const DraftAssembly &assembly);
DraftAssembly decodeAssembly(const QJsonValue &value);
void decodeDraft(const QJsonObject &root, Project &project);
QJsonObject encodeCad(const CadSnapshot &snapshot);
CadSnapshot decodeCad(const QJsonValue &value, bool exactGeometryFields);
QJsonObject encodeSimulation(const SimulationArchive &archive);
SimulationArchive decodeSimulation(const QJsonValue &value, unsigned formatVersion);
QJsonObject encodeLemNetwork(const lem::Network &network);
lem::Network decodeLemNetwork(const QJsonValue &value);
QJsonObject encodeLemSolution(const lem::Solution &solution);
lem::Solution decodeLemSolution(const QJsonValue &value, const lem::Network &network,
                                double frequency);
} // namespace hornsculpt::project_json
