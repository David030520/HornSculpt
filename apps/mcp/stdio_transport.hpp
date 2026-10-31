#pragma once

namespace hornsculpt::mcp {
class ProtocolServer;

// Owns newline-delimited UTF-8 I/O until EOF. Responses use a duplicate of the
// original stdout; process stdout stays redirected to stderr, including after return.
// Call before application libraries can write to stdout. Returns 0 on clean EOF.
int serveStdio(ProtocolServer &server);
} // namespace hornsculpt::mcp
