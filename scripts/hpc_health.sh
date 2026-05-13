#!/usr/bin/env bash
# Probe the HPC cluster login node reachability.
# Exit 0 if reachable, 1 if not.

set -uo pipefail

HOST="${HUMMEL_HOST:-<hpc-login>}"
TIMEOUT="${HUMMEL_TIMEOUT:-5}"

if timeout "$TIMEOUT" ssh -o ConnectTimeout=3 -o BatchMode=yes -o StrictHostKeyChecking=accept-new \
        "$HOST" "uptime" >/dev/null 2>&1; then
    exit 0
fi
exit 1
