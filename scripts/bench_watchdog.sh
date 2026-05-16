#!/usr/bin/env bash
# Bench watchdog — re-launches run_all_species.sh if no instance is running.
# Add to cron: */5 * * * * /home/<user>/llmap-local/scripts/bench_watchdog.sh

set -uo pipefail

LLMAP_HOME="/home/<user>/llmap-local"
LOG="/tmp/full_bench_v2.log"

log() { printf '[%s] %s\n' "$(date -Iseconds)" "$*" >> "$LLMAP_HOME/bench_watchdog.log"; }

# Already running?
if pgrep -x "run_all_speci" >/dev/null 2>&1 || pgrep -af "bash.*run_all_species\.sh" | grep -v watchdog | grep -v "bash -c" >/dev/null 2>&1; then
    log "bench healthy (pgrep hit)"
    exit 0
fi

# Need llmap binary
if [ ! -x "$LLMAP_HOME/build/src/llmap" ]; then
    log "llmap binary missing — building"
    cmake --build "$LLMAP_HOME/build" -j 8 >> "$LLMAP_HOME/bench_watchdog.log" 2>&1 || {
        log "build failed"; exit 1; }
fi

# Re-launch
log "no bench running — relaunching"
cd "$LLMAP_HOME"
nohup setsid bash benchmarks/scripts/run_all_species.sh < /dev/null >> "$LOG" 2>&1 &
disown
log "relaunched bench PID=$!"
