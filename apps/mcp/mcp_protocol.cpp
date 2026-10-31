#include "mcp_protocol.hpp"

#include <QJsonDocument>
#include <algorithm>
#include <cmath>
#include <exception>
#include <utility>

namespace hornsculpt::mcp {
namespace {
const QString modernVersion = QStringLiteral("2026-07-28");
const QString legacyVersion = QStringLiteral("2025-11-25");
const QString versionKey = QStringLiteral("io.modelcontextprotocol/protocolVersion");
const QString clientInfoKey = QStringLiteral("io.modelcontextprotocol/clientInfo");
const QString clientCapabilitiesKey = QStringLiteral("io.modelcontextprotocol/clientCapabilities");
const QString serverInfoKey = QStringLiteral("io.modelcontextprotocol/serverInfo");

QJsonObject error(const QJsonValue &id, int code, const QString &message,
                  const QJsonValue &data = QJsonValue::Undefined) {
    QJsonObject details{{"code", code}, {"message", message}};
    if (!data.isUndefined())
        details.insert("data", data);
    return {{"jsonrpc", "2.0"}, {"id", id}, {"error", details}};
}

bool validId(const QJsonValue &id) {
    return id.isString() || (id.isDouble() && std::isfinite(id.toDouble()) &&
                             std::floor(id.toDouble()) == id.toDouble());
}

bool validImplementation(const QJsonValue &value) {
    const auto object = value.toObject();
    return value.isObject() && object.value("name").isString() &&
           object.value("version").isString();
}

bool includes(const QJsonArray &items, const QString &key, const QString &value) {
    return std::any_of(items.begin(), items.end(), [&](const QJsonValue &item) {
        return item.toObject().value(key).toString() == value;
    });
}

QJsonArray textContent(const QString &text) {
    return {QJsonObject{{"type", "text"}, {"text", text}}};
}
} // namespace

ProtocolServer::ProtocolServer(QString version, ProtocolHandlers handlers)
    : version_(std::move(version)), handlers_(std::move(handlers)) {}

QJsonObject ProtocolServer::serverInfo() const {
    return {{"name", "HornSculpt"}, {"version", version_}};
}

QJsonObject ProtocolServer::capabilities() const {
    QJsonObject available;
    if (handlers_.tools && handlers_.call)
        available.insert("tools", QJsonObject{});
    if (handlers_.resources && handlers_.readResource)
        available.insert("resources", QJsonObject{});
    return available;
}

QJsonObject ProtocolServer::result(const QJsonValue &id, QJsonObject body, bool modern) const {
    if (modern) {
        body.insert("resultType", "complete");
        auto metadata = body.value("_meta").toObject();
        metadata.insert(serverInfoKey, serverInfo());
        body.insert("_meta", metadata);
    }
    return {{"jsonrpc", "2.0"}, {"id", id}, {"result", body}};
}

std::optional<QJsonObject> ProtocolServer::handle(const QJsonValue &message) {
    if (!message.isObject())
        return error(QJsonValue::Null, -32600, "Expected one JSON-RPC request object");
    const auto object = message.toObject();
    const auto rawId = object.value("id");
    const bool notification = rawId.isUndefined();
    const auto id = validId(rawId) ? rawId : QJsonValue(QJsonValue::Null);
    const auto fail = [&](int code, const QString &description) -> std::optional<QJsonObject> {
        if (notification && object.value("method").isString())
            return std::nullopt;
        return error(id, code, description);
    };
    if (object.value("jsonrpc") != QJsonValue("2.0") || !object.value("method").isString() ||
        (!notification && !validId(rawId)) || object.contains("result") || object.contains("error"))
        return fail(-32600, "Invalid JSON-RPC 2.0 request");
    if (object.contains("params") && !object.value("params").isObject())
        return fail(-32602, "params must be an object");

    const auto method = object.value("method").toString();
    const auto params = object.value("params").toObject();
    if (notification) {
        // Request methods sent without an ID must never execute side effects.
        if (method == "notifications/initialized" &&
            legacyState_ == LegacyState::AwaitingNotification)
            legacyState_ = LegacyState::Ready;
        return std::nullopt;
    }
    if (params.contains("_meta") && !params.value("_meta").isObject())
        return fail(-32602, "_meta must be an object");
    const auto metadata = params.value("_meta").toObject();
    const bool explicitVersion = metadata.contains(versionKey);

    if (method == "initialize" && !explicitVersion) {
        if (legacyState_ != LegacyState::Uninitialized)
            return fail(-32600, "The legacy connection is already initialized");
        if (!params.value("protocolVersion").isString() ||
            params.value("protocolVersion").toString().isEmpty() ||
            !params.value("capabilities").isObject() ||
            !validImplementation(params.value("clientInfo")))
            return fail(-32602, "initialize requires protocolVersion, capabilities and clientInfo");
        legacyState_ = LegacyState::AwaitingNotification;
        // A legacy client may propose another version; it must decide whether to accept this one.
        return result(id, {{"protocolVersion", legacyVersion},
                           {"capabilities", capabilities()},
                           {"serverInfo", serverInfo()}},
                      false);
    }

    const bool modern = explicitVersion || legacyState_ == LegacyState::Uninitialized ||
                        method == "server/discover";
    // Legacy lifecycle explicitly permits ping before initialization.
    if (method == "ping" && !explicitVersion && !metadata.contains(clientCapabilitiesKey))
        return result(id, {}, false);
    if (modern) {
        if (!metadata.value(versionKey).isString() ||
            !metadata.value(clientCapabilitiesKey).isObject() ||
            (metadata.contains(clientInfoKey) &&
             !validImplementation(metadata.value(clientInfoKey))))
            return fail(-32602, "Each modern request requires protocolVersion and "
                                "clientCapabilities in params._meta");
        const auto requested = metadata.value(versionKey).toString();
        if (requested != modernVersion)
            return error(id, -32022, "Unsupported protocol version",
                         QJsonObject{{"supported", QJsonArray{modernVersion, legacyVersion}},
                                     {"requested", requested}});
    } else if (legacyState_ != LegacyState::Ready) {
        return fail(-32600, "Send notifications/initialized before using the server");
    }

    if (params.contains("task"))
        return fail(-32602, "MCP Tasks are not supported by this server");
    if (method == "ping")
        return result(id, {}, modern);
    if (method == "server/discover")
        return result(id, {{"supportedVersions", QJsonArray{modernVersion, legacyVersion}},
                           {"capabilities", capabilities()}},
                      true);
    try {
        if (method == "tools/list" && handlers_.tools && handlers_.call) {
            if (params.contains("cursor"))
                return fail(-32602, "Invalid cursor: tools/list returns a single complete page");
            return result(id, {{"tools", handlers_.tools()}}, modern);
        }
        if (method == "resources/list" && handlers_.resources && handlers_.readResource) {
            if (params.contains("cursor"))
                return fail(-32602, "Invalid cursor: resources/list returns a single complete page");
            return result(id, {{"resources", handlers_.resources()}}, modern);
        }
        if (method == "tools/call" && handlers_.tools && handlers_.call) {
            const auto name = params.value("name");
            if (!name.isString() || name.toString().isEmpty() ||
                (params.contains("arguments") && !params.value("arguments").isObject()))
                return fail(-32602, "tools/call requires a name and optional arguments object");
            if (!includes(handlers_.tools(), "name", name.toString()))
                return fail(-32602, "Unknown tool: " + name.toString());
            try {
                const auto value = handlers_.call(name.toString(), params.value("arguments").toObject());
                return result(id, {{"content", textContent(QString::fromUtf8(
                                                 QJsonDocument(value).toJson(QJsonDocument::Compact)))},
                                   {"structuredContent", value},
                                   {"isError", false}},
                              modern);
            } catch (const std::exception &exception) {
                return result(id, {{"content", textContent(QString::fromUtf8(exception.what()))},
                                   {"isError", true}},
                              modern);
            }
        }
        if (method == "resources/read" && handlers_.resources && handlers_.readResource) {
            const auto uri = params.value("uri");
            if (!uri.isString() || uri.toString().isEmpty())
                return fail(-32602, "resources/read requires a nonempty uri");
            if (!includes(handlers_.resources(), "uri", uri.toString()))
                return error(id, modern ? -32602 : -32002, "Resource not found",
                             QJsonObject{{"uri", uri}});
            return result(id, handlers_.readResource(uri.toString()), modern);
        }
    } catch (const std::exception &exception) {
        return fail(-32603, QString::fromUtf8(exception.what()));
    }
    return fail(-32601, "Method not found: " + method);
}

} // namespace hornsculpt::mcp
