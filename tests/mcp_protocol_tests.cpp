#include "mcp_protocol.hpp"

#include <QCoreApplication>
#include <QJsonDocument>
#include <iostream>
#include <stdexcept>

using namespace hornsculpt::mcp;

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}

// Synthetic protocol fixtures deliberately do not stand in for geometry or BEM tests.
struct ProtocolFixture {
    int calls{};
    int reads{};
    ProtocolHandlers handlers() {
        return {
            [] {
                return QJsonArray{QJsonObject{
                    {"name", "fixture.echo"},
                    {"inputSchema", QJsonObject{{"type", "object"}}}}};
            },
            [this](const QString &, const QJsonObject &arguments) {
                ++calls;
                if (arguments.value("throw").toBool())
                    throw std::runtime_error("Fixture execution failed");
                return arguments;
            },
            [] {
                return QJsonArray{
                    QJsonObject{{"uri", "fixture://status"}, {"name", "Protocol fixture"}}};
            },
            [this](const QString &uri) {
                ++reads;
                return QJsonObject{
                    {"contents", QJsonArray{QJsonObject{{"uri", uri},
                                                       {"mimeType", "application/json"},
                                                       {"text", "{\"fixture\":true}"}}}}};
            }};
    }
};

QJsonObject request(const QString &method, QJsonObject params = {}, bool modern = true,
                    const QJsonValue &id = QJsonValue(1)) {
    if (modern)
        params.insert("_meta", QJsonObject{
                                   {"io.modelcontextprotocol/protocolVersion", "2026-07-28"},
                                   {"io.modelcontextprotocol/clientCapabilities", QJsonObject{}}});
    return {{"jsonrpc", "2.0"}, {"id", id}, {"method", method}, {"params", params}};
}

QJsonObject answer(ProtocolServer &server, const QJsonValue &message) {
    const auto response = server.handle(message);
    require(response.has_value(), "Request must receive an answer");
    require(response->value("jsonrpc") == QJsonValue("2.0"), "JSON-RPC version is preserved");
    return *response;
}

void expectError(const QJsonObject &response, int code) {
    require(!response.contains("result"), "Protocol failure must not have a result");
    require(response.value("error").toObject().value("code").toInt() == code,
            "Unexpected JSON-RPC error code");
}

void modernRequests() {
    ProtocolFixture fixture;
    ProtocolServer server("test-version", fixture.handlers());
    auto response = answer(server, request("tools/list"));
    const auto result = response.value("result").toObject();
    require(result.value("resultType") == QJsonValue("complete"), "Modern complete result required");
    require(result.value("tools").toArray().size() == 1, "Direct modern tool discovery works");
    require(result.value("_meta")
                .toObject()
                .value("io.modelcontextprotocol/serverInfo")
                .toObject()
                .value("version") == QJsonValue("test-version"),
            "Modern results identify the server");

    response = answer(server, request("server/discover", {}, true, "discovery-id"));
    require(response.value("id") == QJsonValue("discovery-id"), "String request IDs preserved");
    const auto discovery = response.value("result").toObject();
    require(discovery.value("supportedVersions").toArray() ==
                QJsonArray{"2026-07-28", "2025-11-25"},
            "Only implemented revisions advertised");
    const auto capabilities = discovery.value("capabilities").toObject();
    require(capabilities.contains("tools") && capabilities.contains("resources") &&
                capabilities.size() == 2,
            "No unimplemented capabilities advertised");
    require(capabilities.value("resources").toObject().isEmpty(),
            "Subscriptions and listChanged are not advertised");
    require(answer(server, request("ping")).value("result").toObject().value("resultType") ==
                QJsonValue("complete"),
            "Modern ping returns a complete result");

    ProtocolServer emptyServer("test", {});
    require(answer(emptyServer, request("server/discover"))
                .value("result")
                .toObject()
                .value("capabilities")
                .toObject()
                .isEmpty(),
            "Missing handlers are not advertised");
    expectError(answer(emptyServer, request("tools/list")), -32601);
}

void legacyLifecycle() {
    ProtocolFixture fixture;
    ProtocolServer server("test", fixture.handlers());
    require(answer(server, request("ping", {}, false)).contains("result"),
            "Legacy ping works before initialization");
    expectError(answer(server, request("tools/list", {}, false)), -32602);
    auto initialize = request("initialize",
                              {{"protocolVersion", "2025-11-25"},
                               {"capabilities", QJsonObject{}},
                               {"clientInfo", QJsonObject{{"name", "fixture"}, {"version", "1"}}}},
                              false);
    const auto response = answer(server, initialize).value("result").toObject();
    require(response.value("protocolVersion") == QJsonValue("2025-11-25"),
            "Legacy revision negotiated explicitly");
    require(!response.contains("resultType") && response.contains("serverInfo"),
            "Legacy initialization uses its own schema");
    expectError(answer(server, request("tools/list", {}, false)), -32600);
    expectError(answer(server, initialize), -32600);
    require(!server.handle(QJsonObject{{"jsonrpc", "2.0"},
                                      {"method", "notifications/initialized"}}),
            "Initialized is a notification");
    require(answer(server, request("tools/list", {}, false)).contains("result"),
            "Legacy requests work after handshake");
    require(answer(server, request("server/discover")).contains("result"),
            "Modern requests remain independent of legacy session state");
    expectError(answer(server, request("resources/read", {{"uri", "fixture://missing"}}, false)),
                -32002);

    ProtocolServer negotiation("test", fixture.handlers());
    auto params = initialize.value("params").toObject();
    params.insert("protocolVersion", "1900-01-01");
    initialize.insert("params", params);
    require(answer(negotiation, initialize)
                .value("result")
                .toObject()
                .value("protocolVersion") == QJsonValue("2025-11-25"),
            "Legacy negotiation offers the implemented legacy revision");
}

void validationAndNotifications() {
    ProtocolFixture fixture;
    ProtocolServer server("test", fixture.handlers());
    expectError(answer(server, QJsonArray{request("ping")}), -32600);
    expectError(answer(server, QJsonValue(7)), -32600);
    expectError(answer(server, request("tools/list", {}, true, QJsonValue::Null)), -32600);
    expectError(answer(server, request("tools/list", {}, true, 1.5)), -32600);
    auto malformed = request("tools/list");
    malformed.insert("params", QJsonArray{});
    expectError(answer(server, malformed), -32602);
    malformed = request("tools/list", {{"_meta", QJsonObject{}}}, false);
    expectError(answer(server, malformed), -32602);
    auto params = QJsonObject{{"_meta", QJsonObject{
                                           {"io.modelcontextprotocol/protocolVersion", "1900-01-01"},
                                           {"io.modelcontextprotocol/clientCapabilities", QJsonObject{}}}};
    const auto unsupported = answer(server, request("tools/list", params, false));
    expectError(unsupported, -32022);
    require(unsupported.value("error").toObject().value("data").toObject().value("requested") ==
                QJsonValue("1900-01-01"),
            "Version error provides the unsupported request version");
    params = {{"_meta", QJsonObject{{"io.modelcontextprotocol/protocolVersion", "2026-07-28"}}}};
    expectError(answer(server, request("tools/list", params, false)), -32602);
    expectError(answer(server, request("missing/method")), -32601);
    expectError(answer(server, request("tools/list", {{"cursor", "unissued"}})), -32602);
    expectError(answer(server, request("resources/list", {{"cursor", 17}})), -32602);
    expectError(answer(server, request("tools/call", {{"name", "fixture.echo"},
                                                      {"task", QJsonObject{}}})),
                -32602);

    auto notification = request("tools/call", {{"name", "fixture.echo"}});
    notification.remove("id");
    require(!server.handle(notification), "Notifications never receive responses");
    require(fixture.calls == 0, "Request methods sent as notifications never mutate the application");
    notification.insert("method", "notifications/unknown");
    require(!server.handle(notification), "Unknown notifications are ignored");
    notification.insert("params", QJsonArray{});
    require(!server.handle(notification), "Malformed notification does not receive an error response");

    const auto largeId = QJsonDocument::fromJson(
                             R"({"jsonrpc":"2.0","id":9007199254740993,"method":"ping"})")
                             .object();
    require(answer(server, largeId).value("id") == largeId.value("id"),
            "Integer IDs are not converted through double precision");
}

void toolsAndResources() {
    ProtocolFixture fixture;
    ProtocolServer server("test", fixture.handlers());
    const QJsonObject arguments{{"label", QString::fromUtf8("Hörner\nSpalt")},
                                {"value", 17.5}};
    const auto result = answer(server, request("tools/call", {{"name", "fixture.echo"},
                                                              {"arguments", arguments}}))
                            .value("result")
                            .toObject();
    require(result.value("structuredContent").toObject() == arguments,
            "Structured tool data is preserved");
    const auto text = result.value("content").toArray().first().toObject().value("text").toString();
    require(QJsonDocument::fromJson(text.toUtf8()).object() == arguments,
            "Text and structured content represent the same data");
    require(!result.value("isError").toBool() && fixture.calls == 1, "Tool execution succeeds once");
    const auto failed = answer(server, request("tools/call", {{"name", "fixture.echo"},
                                                              {"arguments", QJsonObject{{"throw", true}}}}));
    require(!failed.contains("error") &&
                failed.value("result").toObject().value("isError").toBool(),
            "Tool failures are model-readable execution errors, not JSON-RPC failures");
    expectError(answer(server, request("tools/call", {{"name", "unknown"}})), -32602);
    expectError(answer(server, request("tools/call", {{"name", "fixture.echo"},
                                                      {"arguments", QJsonArray{}}})),
                -32602);
    expectError(answer(server, request("tools/call")), -32602);
    require(fixture.calls == 2, "Malformed and unknown calls do not reach the engine");
    require(answer(server, request("resources/list"))
                .value("result").toObject().value("resources").toArray().size() == 1,
            "Resource listing uses the handler");
    require(answer(server, request("resources/read", {{"uri", "fixture://status"}}))
                .value("result").toObject().value("contents").toArray().size() == 1,
            "Resource reading uses the handler");
    expectError(answer(server, request("resources/read", {{"uri", "fixture://missing"}})), -32602);
    expectError(answer(server, request("resources/read", {{"uri", 17}})), -32602);
    require(fixture.reads == 1, "Invalid resource requests do not reach the engine");
    auto failingHandlers = fixture.handlers();
    failingHandlers.resources = []() -> QJsonArray { throw std::runtime_error("Fixture list failure"); };
    ProtocolServer failing("test", failingHandlers);
    expectError(answer(failing, request("resources/list")), -32603);
}
} // namespace

int main(int argc, char **argv) {
    QCoreApplication application(argc, argv);
    try {
        modernRequests();
        legacyLifecycle();
        validationAndNotifications();
        toolsAndResources();
        std::cout << "MCP protocol fixture tests passed\n";
        return 0;
    } catch (const std::exception &exception) {
        std::cerr << "MCP protocol test failure: " << exception.what() << '\n';
        return 1;
    }
}
