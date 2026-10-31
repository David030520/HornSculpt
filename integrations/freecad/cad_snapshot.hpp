#pragma once

#include <QByteArray>
#include <hornsculpt/freecad_link.hpp>

namespace hornsculpt::cad_detail {

void validatePlacement(const std::array<double, 16> &placement);
CadSnapshot decodeSnapshot(const QByteArray &payload, const QString &sessionId);
void reconcileSnapshot(CadSnapshot &current, const CadSnapshot *previous,
                       std::uint32_t &nextFaceId);
QByteArray frame(const QByteArray &payload);

} // namespace hornsculpt::cad_detail
