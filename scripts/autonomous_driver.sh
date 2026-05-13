#!/usr/bin/env bash
# LLmap autonomous driver — invoked every 15 min by cron.
#
# Responsibilities:
#   1. Refuse to overlap with a previous iteration (flock).
#   2. Stop after the 96h horizon.
#   3. Probe Hummel-2; alert via Zyrkel on UP↔DOWN transitions.
#   4. git pull, spawn `claude` CLI to continue Phase 0-9 work, git push.
#   5. Log every iteration.

set -uo pipefail

# --- environment: cron does not inherit shell env, source secrets ourselves
if [[ -f "${HOME}/.env" ]]; then
    set -a
    # shellcheck disable=SC1090
    . "${HOME}/.env"
    set +a
fi
# Normalize token name (the .env uses GITHUB_TOKEN_SCHLEIN_LAB,
# git_askpass.sh and other scripts expect SCHLEIN_LAB_TOKEN)
: "${SCHLEIN_LAB_TOKEN:=${GITHUB_TOKEN_SCHLEIN_LAB:-}}"
export SCHLEIN_LAB_TOKEN
export GIT_ASKPASS="${HOME}/llmap-local/scripts/git_askpass.sh"
export GIT_TERMINAL_PROMPT=0

LLMAP_HOME="${LLMAP_HOME:-${HOME}/llmap-local}"
LOG="${LLMAP_HOME}/autonomous_run.log"
STATE="${LLMAP_HOME}/STATE.md"
DEADLINE_FILE="${LLMAP_HOME}/.deadline_epoch"
HUMMEL_STATUS_FILE="${LLMAP_HOME}/.hummel_status"
LOCK="${LLMAP_HOME}/.driver.lock"
ITER_FILE="${LLMAP_HOME}/.iteration_count"

ts()  { date +"%Y-%m-%dT%H:%M:%S%z"; }
log() { mkdir -p "$(dirname "$LOG")"; echo "[$(ts)] $*" | tee -a "$LOG" >&2; }

cd "$LLMAP_HOME" || { log "FATAL: $LLMAP_HOME missing"; exit 1; }

# --- mutual exclusion via flock ------------------------------------------
exec 200>"$LOCK"
if ! flock -n 200; then
    log "previous iteration still running, skipping"
    exit 0
fi

# --- 96h deadline guard --------------------------------------------------
if [[ ! -f "$DEADLINE_FILE" ]]; then
    echo $(( $(date +%s) + 96 * 3600 )) > "$DEADLINE_FILE"
    log "deadline set: $(date -d @$(cat "$DEADLINE_FILE"))"
fi
now=$(date +%s)
deadline=$(cat "$DEADLINE_FILE")
if (( now > deadline )); then
    log "96h deadline reached — stopping autonomous loop"
    # Send Zyrkel notification on stop (once)
    if [[ ! -f "${LLMAP_HOME}/.deadline_notified" ]]; then
        "${LLMAP_HOME}/scripts/zyrkel_notify.sh" \
            "🏁 LLmap autonomous run complete (96h deadline reached). Iterations: $(cat "$ITER_FILE" 2>/dev/null || echo 0). Check schlein-lab/LLmap repo for progress." \
            || true
        touch "${LLMAP_HOME}/.deadline_notified"
    fi
    exit 0
fi

# --- iteration counter ---------------------------------------------------
iter=$(( $(cat "$ITER_FILE" 2>/dev/null || echo 0) + 1 ))
echo "$iter" > "$ITER_FILE"
log "=== iteration $iter ==="

# --- Hummel health check + transition alerts -----------------------------
prev_status=$(cat "$HUMMEL_STATUS_FILE" 2>/dev/null || echo "UNKNOWN")
if "${LLMAP_HOME}/scripts/hummel_health.sh"; then
    cur_status="UP"
else
    cur_status="DOWN"
fi
echo "$cur_status" > "$HUMMEL_STATUS_FILE"

if [[ "$cur_status" != "$prev_status" ]]; then
    log "Hummel status transition: $prev_status -> $cur_status"
    # Suppress notification on first run when prev was UNKNOWN — initial probe
    # is not a real transition.
    if [[ "$prev_status" != "UNKNOWN" ]]; then
        if [[ "$cur_status" == "DOWN" ]]; then
            "${LLMAP_HOME}/scripts/zyrkel_notify.sh" \
                "🚨 Hummel-2 unreachable. VPN renewal needed. LLmap autonomous loop paused on Hummel-dependent tasks." \
                || true
        else
            "${LLMAP_HOME}/scripts/zyrkel_notify.sh" \
                "✅ Hummel-2 reachable again. LLmap autonomous loop resuming." \
                || true
        fi
    fi
fi

if [[ "$cur_status" == "DOWN" ]]; then
    log "Hummel down — skipping Hummel-dependent work this iteration"
    # Still allow local-only Claude work (docs, host-side code) to proceed
fi

# --- git sync (stash any uncommitted to survive concurrent edits) --------
stashed=0
if [[ -n "$(git status --porcelain)" ]]; then
    log "uncommitted changes present, stashing before pull"
    git stash push -u -m "driver auto-stash iteration $iter" >> "$LOG" 2>&1 && stashed=1
fi
log "git pull --rebase origin main"
git pull --rebase origin main >> "$LOG" 2>&1 || log "WARN: git pull failed"
if (( stashed )); then
    log "restoring auto-stash"
    git stash pop >> "$LOG" 2>&1 || log "WARN: git stash pop failed (manual inspection needed)"
fi

# --- spawn claude continuation -------------------------------------------
PROMPT_FILE="${LLMAP_HOME}/scripts/continuation_prompt.md"
if [[ ! -f "$PROMPT_FILE" ]]; then
    log "FATAL: $PROMPT_FILE missing"
    exit 1
fi

log "spawning claude continuation (subprocess, non-interactive)"
# --print: non-interactive
# --output-format text: clean log lines
# Timeout: 13 minutes — leaves 2 minutes headroom before next cron cycle
timeout 780 claude --print --output-format text \
    --add-dir "$LLMAP_HOME" \
    "$(cat "$PROMPT_FILE")" \
    >> "$LOG" 2>&1
claude_rc=$?
log "claude subprocess exit code: $claude_rc"

# --- commit + push if Claude made changes --------------------------------
if [[ -n "$(git status --porcelain)" ]]; then
    log "Claude made uncommitted changes — committing fallback safety net"
    git add -A
    git commit -m "autonomous: iteration $iter fallback commit" >> "$LOG" 2>&1 || true
fi

# Push whatever is on main
git push origin main >> "$LOG" 2>&1 || log "WARN: git push failed"

# --- stuck-iteration detection -------------------------------------------
LAST_COMMIT_FILE="${LLMAP_HOME}/.last_commit_sha"
cur_sha=$(git rev-parse HEAD)
prev_sha=$(cat "$LAST_COMMIT_FILE" 2>/dev/null || echo "none")
echo "$cur_sha" > "$LAST_COMMIT_FILE"

STUCK_FILE="${LLMAP_HOME}/.stuck_count"
if [[ "$cur_sha" == "$prev_sha" ]]; then
    stuck=$(( $(cat "$STUCK_FILE" 2>/dev/null || echo 0) + 1 ))
    echo "$stuck" > "$STUCK_FILE"
    log "no new commits this iteration (stuck count: $stuck)"
    if (( stuck >= 3 )) && [[ ! -f "${LLMAP_HOME}/.stuck_notified_$stuck" ]]; then
        "${LLMAP_HOME}/scripts/zyrkel_notify.sh" \
            "⚠️ LLmap autonomous loop stuck — $stuck consecutive iterations without commits. Iteration $iter. Check STATE.md and logs." \
            || true
        touch "${LLMAP_HOME}/.stuck_notified_$stuck"
    fi
else
    echo 0 > "$STUCK_FILE"
    rm -f "${LLMAP_HOME}"/.stuck_notified_* 2>/dev/null
fi

log "=== iteration $iter complete ==="
