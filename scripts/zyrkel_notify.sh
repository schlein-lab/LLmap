#!/usr/bin/env bash
# Send a Zyrkel notification by invoking the local claude CLI with the
# zyrkel MCP available. Falls back silently if claude is missing.
#
# Usage: zyrkel_notify.sh "message text"

set -uo pipefail

MSG="${1:-LLmap notification (no message)}"

if ! command -v claude >/dev/null 2>&1; then
    echo "[zyrkel_notify] claude CLI missing, message dropped: $MSG" >&2
    exit 1
fi

# Short prompt — claude calls mcp__zyrkel__notify with the message.
PROMPT="Bitte rufe das tool mcp__zyrkel__notify mit folgender message auf, sonst nichts: ${MSG}"

# Strict timeout — notify is best-effort, must never block driver
timeout 60 claude --print --output-format text "$PROMPT" >/dev/null 2>&1 || {
    echo "[zyrkel_notify] claude invocation failed, message dropped: $MSG" >&2
    exit 1
}
exit 0
