# LLmap Autonomous Build State

This file is the source of truth for autonomous-driver continuation. The driver script reads `current_task` to decide what to do next, and writes back progress + next task after each iteration.

---

## Run metadata

| Field | Value |
|---|---|
| Run start | 2026-05-13 (autonomous bootstrap) |
| Run horizon | 96h |
| Driver cadence | every 15 min |
| Hummel-2 status | required for heavy jobs |
| Local-box status | required for driver + Claude CLI |
| Last successful iteration | 2 |
| Total iterations | 2 |

---

## Phase progress

- [x] **Phase -1: Bootstrap** — repo + scaffolding + driver
- [ ] **Phase 0: Foundations** — CMake, AlignmentRecord, BucketPyramid, synthetic data
  - [x] Phase 0.1: AlignmentRecord type + tests (8 tests pass)
  - [x] Phase 0.2: BucketPyramid data structure + serialization (11 tests pass)
- [ ] Phase 1: Foundation Model Integration
- [ ] Phase 2: Stage 1 Self-Interference
- [ ] Phase 3: Stage 2 Reference WaveCollapse
- [ ] Phase 4: Classical Path + WFA2
- [ ] **Phase 5: KILL-SWITCH VALIDATION** ★
- [ ] Phase 6: Dual Output (BAM + Parquet)
- [ ] Phase 7: Claude Agent Integration
- [ ] Phase 8: Performance Optimization
- [ ] Phase 9: Single-Cell + Paralog Production

---

## Current task

```
phase: 0
task: wave_state_sparse_csr
substep: 3/4
last_action: implemented BucketPyramid with L0/L1/L2 hierarchy, serialization, 11 tests pass
next_action: implement src/core/wave_state.{h,cpp} with sparse CSR format for Read×Bucket probability matrix; CPU implementation with GPU stub (d_ vectors placeholder); write tests/unit/test_wave_state.cpp
acceptance: WaveState builds, stores per-read bucket probabilities in CSR format; accessors work; collapse/dropout flag management works; tests pass
```

---

## Next task queue (FIFO)

1. ~~Phase 0.1: AlignmentRecord type + tests~~ ✅ done
2. ~~Phase 0.2: BucketPyramid data structure + serialization~~ ✅ done
3. Phase 0.3: WaveState sparse CSR + GPU stub ← CURRENT
4. Phase 0.4: Synthetic IGH-locus generator + simulator wrappers (pbsim3, ART)
5. Phase 1.1: ONNX Runtime CUDA-EP wired
6. Phase 1.2: Caduceus-Ph distilled embedder
7. Phase 1.3: Bucket embedder via Evo-1.5B
8. ... (continues per LLmap_SPEC.md)

---

## Blockers / open questions

- _none yet_

---

## Hummel-2 workspace location

```
~/llmap/                          # working copy clone
/beegfs/u/bbg6775/llmap/          # build + data (when needed)
```

---

## Driver invocations log

| Iteration | Timestamp | Hummel-up? | Action | Outcome |
|---|---|---|---|---|
| 0 | 2026-05-13 | yes | bootstrap | repo created, scaffold committed |
| 1 | 2026-05-13 | n/a | verify alignment_record | confirmed 8 tests pass; advanced to phase 0.2 |
| 2 | 2026-05-13 | n/a | implement bucket_pyramid | L0/L1/L2 hierarchy + serialization; 11 tests pass |

---

## How the autonomous driver reads this file

The driver script `scripts/autonomous_driver.sh` runs every 15 min via cron. Each invocation:

1. Pings Hummel — if down, sends Zyrkel alert and exits cleanly (will retry next cycle)
2. `cd ~/llmap-local && git pull origin main`
3. Reads `STATE.md` → extracts `current_task` block
4. Spawns `claude -p "$(cat scripts/continuation_prompt.md)"` — Claude CLI subprocess with full context
5. Claude continues the work: writes code, builds, tests, commits, pushes
6. Driver verifies Claude updated `STATE.md` with progress
7. Logs to `autonomous_run.log` and `driver invocations log` table above

If Claude fails to advance for 3 consecutive iterations, driver pages via Zyrkel.
