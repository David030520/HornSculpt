#pragma once

#include <QJsonArray>
#include <QJsonObject>
#include <QJsonValue>
#include <QString>
#include <functional>
#include <optional>

namespace hornsculpt::mcp {

struct ProtocolHandlers {
    std::function<QJsonArray()> tools;
    std::function<QJsonObject(const QString &, const QJsonObject &)> call;
    std::function<QJsonArray()> resources;
    std::function<QJsonObject(const QString &)> readResource;
};

// MCP 2026-07-28 per-request metadata and the explicit 2025-11-25 legacy handshake.
// The callbacks own application state; this class owns only protocol negotiation.
class ProtocolServer {
  public:
    ProtocolServer(QString version, ProtocolHandlers handlers);
    [[nodiscard]] std::optional<QJsonObject> handle(const QJsonValue &message);

  private:
    enum class LegacyState { Uninitialized, AwaitingNotification, Ready };
    QString version_;
    ProtocolHandlers handlers_;
    LegacyState legacyState_{LegacyState::Uninitialized};

    [[nodiscard]] QJsonObject serverInfo() const;
    [[nodiscard]] QJsonObject capabilities() const;
    [[nodiscard]] QJsonObject result(const QJsonValue &id, QJsonObject body, bool modern) const;
};

} // namespace hornsculpt::mcp
