#!/usr/bin/env bash
# timing_harness.sh — wrap a mapper invocation in /usr/bin/time and emit timing.json
#
# Usage:
#   timing_harness.sh <output_dir> -- <mapper command> [args...]
#
# The double-dash is mandatory and separates harness arguments from the
# mapped command line.  The harness writes two files into <output_dir>:
#
#   timing.raw        — verbose output of /usr/bin/time (-v style)
#   timing.json       — normalised JSON with keys:
#                         wall_sec, user_sec, sys_sec, max_rss_kb, cpu_pct,
#                         exit_status, command
#
# The harness is idempotent — re-running overwrites timing.* in place.
# stdout/stderr of the wrapped command pass through untouched so the caller
# can still capture SAM/BAM via redirection.

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "usage: $0 <output_dir> -- <command...>" >&2
    exit 64
fi

OUT_DIR="$1"; shift
if [[ "$1" != "--" ]]; then
    echo "ERROR: expected '--' after output_dir, got '$1'" >&2
    exit 64
fi
shift

mkdir -p "$OUT_DIR"
RAW="$OUT_DIR/timing.raw"
JSON="$OUT_DIR/timing.json"

# Locate /usr/bin/time — GNU time, not the bash builtin.  On some distros
# it lives under /usr/local/bin or is provided by the 'time' package.
TIME_BIN=""
for cand in /usr/bin/time /usr/local/bin/time; do
    if [[ -x "$cand" ]]; then
        TIME_BIN="$cand"
        break
    fi
done
if [[ -z "$TIME_BIN" ]]; then
    echo "ERROR: /usr/bin/time not found — install package 'time'" >&2
    exit 69
fi

# Run the wrapped command.  Field order: wall, user, sys, max_rss_kb, cpu_pct, exit
# %e wall (sec), %U user (sec), %S sys (sec), %M max RSS (kB), %P CPU %
FMT='WALL=%e USER=%U SYS=%S MAXRSS_KB=%M CPU_PCT=%P EXIT=%x'

set +e
"$TIME_BIN" -f "$FMT" -o "$RAW" -- "$@"
EXIT_CODE=$?
set -e

# Parse the single-line format.  We tolerate locale-specific decimal commas
# (de_DE.UTF-8 puts %e as "31,86") by normalising to dots.
LINE=$(tail -n 1 "$RAW" | tr ',' '.')

extract() {
    local key="$1"
    # shellcheck disable=SC2001
    echo "$LINE" | sed -n "s/.*${key}=\([^ ]*\).*/\1/p"
}

WALL=$(extract WALL)
USER=$(extract USER)
SYS=$(extract SYS)
MAXRSS=$(extract MAXRSS_KB)
CPU=$(extract CPU_PCT)
CPU="${CPU%\%}"      # strip trailing percent sign if present

# Fallback defaults so the JSON is always well-formed.
: "${WALL:=0}"
: "${USER:=0}"
: "${SYS:=0}"
: "${MAXRSS:=0}"
: "${CPU:=0}"

# Build the JSON manually — no jq dependency.  Strings are escaped via a
# python one-liner if python3 is available, else we trust the inputs (which
# are numeric except for COMMAND).
CMD_STR="$*"
if command -v python3 >/dev/null 2>&1; then
    python3 - "$JSON" "$WALL" "$USER" "$SYS" "$MAXRSS" "$CPU" "$EXIT_CODE" "$CMD_STR" <<'PY'
import json, sys
out, wall, user, sysc, rss, cpu, exit_code, cmd = sys.argv[1:9]
def f(x):
    try: return float(x)
    except: return 0.0
def i(x):
    try: return int(float(x))
    except: return 0
payload = {
    "wall_sec":    f(wall),
    "user_sec":    f(user),
    "sys_sec":     f(sysc),
    "max_rss_kb":  i(rss),
    "cpu_pct":     f(cpu),
    "exit_status": i(exit_code),
    "command":     cmd,
}
with open(out, "w", encoding="utf-8") as fh:
    json.dump(payload, fh, indent=2, ensure_ascii=False)
    fh.write("\n")
PY
else
    # Last-resort plain-bash emitter — assumes no double quotes in CMD_STR.
    cat > "$JSON" <<EOF
{
  "wall_sec":    ${WALL},
  "user_sec":    ${USER},
  "sys_sec":     ${SYS},
  "max_rss_kb":  ${MAXRSS},
  "cpu_pct":     ${CPU},
  "exit_status": ${EXIT_CODE},
  "command":     "${CMD_STR//\"/\\\"}"
}
EOF
fi

exit "$EXIT_CODE"
