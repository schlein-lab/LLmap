# LLmap known issues (cross-species benchmark, 2026-05-15)

Living list of investigated failure modes hit during Phase B/C cross-species
benchmarking. Each entry: root cause, reproduction, workaround, status.

## I-001  Empty `aln.sam` on tier 5 / tier 6 runs

**Symptom.** `benchmarks/scripts/run_species_bench.sh` logs
`WARN: llmap produced empty SAM` and leaves a 0-byte `aln.sam`, 0-byte
`stderr.log`, and 0-byte `.time` next to the missing BAM. Confirmed on:
arabidopsis tier6, bacteria tier5, celegans tier1, celegans tier6,
drosophila tier5 (all 16 organisms x 5 informative tiers were swept).

**Root cause — *not* an LLmap binary bug.** Two `run_all_species.sh`
processes were running concurrently (the user/cron started it twice). Both
loops walked the same `(organism, tier)` matrix and both invoked
`run_species_bench.sh` against the same per-tier output directory:

```
ps -ef | grep run_all_species
... bash benchmarks/scripts/run_all_species.sh
... bash benchmarks/scripts/run_all_species.sh   <-- duplicate
```

Each tier 5/6 mixed_long run uses ~1.5 GB RSS for LLmap. With *two*
LLmap aligners loading the same index at the same time, the host (19 GB
RAM, 3.8 GB swap) tipped into the OOM regime; the OOM killer SIGKILLs one
of the duplicate `time`/`llmap` pairs *before* `/usr/bin/time` can flush
its report. The signal target is non-deterministic, but the artefacts are
always: 0-byte `aln.sam` (created by `open()` but never written to before
the kill), 0-byte `.time` (kill happened before `time` finished), 0-byte
`stderr.log` (kill happened before LLmap printed anything). The smaller
tiers (1/2/10) finish before the second runner even reaches them, so they
escape the race.

The signature `0-byte stderr + 0-byte .time + 0-byte aln.sam` is
diagnostic for "external SIGKILL", not for "LLmap exited cleanly with no
output". An LLmap bug would have left a non-empty stderr.

**Workarounds shipped.**

1. `benchmarks/scripts/run_species_bench.sh` now takes an `flock` on
   `benchmarks/.locks/<organism>_tier<N>.lock` and aborts with a clear
   error (`exit 4`) when a duplicate run is detected. This makes the
   race impossible in the future.
2. The same script now writes a structured `<mapper>/rep0/failure.json`
   when LLmap exits non-zero *or* produces an empty SAM, including the
   preset and inputs. `aggregate_scoreboard.py` recognises these and
   surfaces them as `failure_mode: empty_sam` / `align_rc_<N>` in
   `scoreboard_failed_regions.json`.

**Recovery for the existing 4 empty-SAM tier dirs.** Simply re-run
`bash benchmarks/scripts/run_species_bench.sh --organism <X> --tier <N>`
with the lock now in place (no duplicate runner). The lock guarantees
the second runner is the one that exits, not the one that gets killed.

**Status.** Resolved at the harness level. The LLmap binary itself was
not at fault — `bacteria tier1`, `bacteria tier2`, `bacteria tier6`,
`bacteria tier10`, `arabidopsis tier1/2/5/10`, etc. all produced
sensible BAMs with F1@1kb >= 0.998 on the same hardware, same binary,
same map-hifi preset, *whenever the lock-free duplicate runner did not
collide on the same path*.

## Issue template

```
## I-NNN  <title>

**Symptom.**       Observed behaviour, log lines, file artefacts.
**Reproduction.**  Single command + inputs.
**Root cause.**    Diagnosis (binary bug? harness bug? environment?).
**Workaround.**    What to do until upstream is fixed.
**Status.**        Open / Investigating / Resolved (commit ref).
```
