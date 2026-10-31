#include "stdio_transport.hpp"

#include "mcp_protocol.hpp"

#include <QByteArray>
#include <QJsonDocument>
#include <QJsonParseError>
#include <QStringConverter>
#include <cstdio>
#include <memory>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#else
#include <unistd.h>
#endif

namespace hornsculpt::mcp {
namespace {
// Bounds both protocol framing memory and single-response size. Large engineering
// data is retrieved through bounded application queries, never an unbounded line.
constexpr qsizetype maximumMessageBytes = 8 * 1024 * 1024;
using File = std::unique_ptr<FILE, decltype(&std::fclose)>;

File isolateOutput() {
    std::fflush(stdout);
#ifdef _WIN32
    const auto output = _dup(_fileno(stdout));
    if (output < 0)
        return File(nullptr, &std::fclose);
    if (_dup2(_fileno(stderr), _fileno(stdout)) != 0) {
        _close(output);
        return File(nullptr, &std::fclose);
    }
    if (_setmode(_fileno(stdin), _O_BINARY) < 0 || _setmode(output, _O_BINARY) < 0) {
        _close(output);
        return File(nullptr, &std::fclose);
    }
    auto *stream = _fdopen(output, "wb");
    if (!stream)
        _close(output);
#else
    const auto output = dup(fileno(stdout));
    if (output < 0)
        return File(nullptr, &std::fclose);
    if (dup2(fileno(stderr), fileno(stdout)) < 0) {
        close(output);
        return File(nullptr, &std::fclose);
    }
    auto *stream = fdopen(output, "wb");
    if (!stream)
        close(output);
#endif
    return File(stream, &std::fclose);
}

QJsonObject framingError(int code, const QString &message) {
    return {{"jsonrpc", "2.0"},
            {"id", QJsonValue::Null},
            {"error", QJsonObject{{"code", code}, {"message", message}}}};
}

bool writeMessage(FILE *output, QJsonObject message) {
    auto bytes = QJsonDocument(message).toJson(QJsonDocument::Compact);
    if (bytes.size() > maximumMessageBytes) {
        const auto id = message.value("id");
        message = framingError(-32603, "Response exceeds the 8 MiB protocol limit");
        message.insert("id", id);
        bytes = QJsonDocument(message).toJson(QJsonDocument::Compact);
    }
    bytes.append('\n');
    return std::fwrite(bytes.constData(), 1, static_cast<size_t>(bytes.size()), output) ==
               static_cast<size_t>(bytes.size()) &&
           std::fflush(output) == 0;
}

bool dispatchLine(ProtocolServer &server, FILE *output, QByteArray line) {
    if (line.endsWith('\r'))
        line.chop(1);
    QStringDecoder decoder(QStringDecoder::Utf8);
    const QString decoded = decoder.decode(line);
    if (decoder.hasError() || decoded.toUtf8() != line)
        return writeMessage(output, framingError(-32700, "Message is not valid UTF-8"));
    QJsonParseError parseError;
    const auto document = QJsonDocument::fromJson(line, &parseError);
    if (parseError.error != QJsonParseError::NoError)
        return writeMessage(output, framingError(-32700, parseError.errorString()));
    const auto response = server.handle(document.isObject() ? QJsonValue(document.object())
                                                           : QJsonValue(document.array()));
    return !response || writeMessage(output, *response);
}
} // namespace

int serveStdio(ProtocolServer &server) {
    auto output = isolateOutput();
    if (!output) {
        std::fputs("HornSculpt MCP: cannot isolate standard streams\n", stderr);
        return 1;
    }
    QByteArray line;
    bool oversized = false;
    for (;;) {
        const auto character = std::fgetc(stdin);
        if (character == EOF) {
            if (std::ferror(stdin)) {
                std::fputs("HornSculpt MCP: standard input read failed\n", stderr);
                return 1;
            }
            // EOF is shutdown, but a complete last request without its newline can
            // still be processed. A partial JSON document receives a parse error.
            if (oversized)
                return writeMessage(output.get(), framingError(-32600, "Message exceeds 8 MiB"))
                           ? 0
                           : 1;
            return line.isEmpty() || dispatchLine(server, output.get(), line) ? 0 : 1;
        }
        if (character == '\n') {
            const bool written = oversized
                                     ? writeMessage(output.get(),
                                                    framingError(-32600, "Message exceeds 8 MiB"))
                                     : dispatchLine(server, output.get(), line);
            if (!written)
                return 1;
            line.clear();
            oversized = false;
        } else if (!oversized) {
            if (line.size() == maximumMessageBytes) {
                line.clear();
                oversized = true;
            } else {
                line.append(static_cast<char>(character));
            }
        }
    }
}
} // namespace hornsculpt::mcp
