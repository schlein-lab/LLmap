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
| Last successful iteration | 12 |
| Total iterations | 12 |

---

## Phase progress

- [x] **Phase -1: Bootstrap** — repo + scaffolding + driver
- [x] **Phase 0: Foundations** — CMake, AlignmentRecord, BucketPyramid, synthetic data
  - [x] Phase 0.1: AlignmentRecord type + tests (8 tests pass)
  - [x] Phase 0.2: BucketPyramid data structure + serialization (11 tests pass)
  - [x] Phase 0.3: WaveState sparse CSR + GPU stub (21 tests pass)
  - [x] Phase 0.4: Synthetic IGH-locus generator + simulator wrappers (56 tests pass)
- [x] **Phase 1: Foundation Model Integration**
  - [x] Phase 1.1: ONNX Runtime CUDA-EP wired (75 tests pass)
  - [x] Phase 1.2: Test ONNX model + verified embedder inference (85 tests pass)
  - [x] Phase 1.3: BucketEmbedder + benchmark (110 tests pass)
- [ ] Phase 2: Stage 1 Self-Interference
  - [x] Phase 2.1: FAISS-GPU IndexIVFFlat wrapper (140 tests pass)
  - [x] Phase 2.2: Sparse k-NN extraction → similarity graph (180 tests pass)
  - [x] Phase 2.3: Leiden community detection (214 tests pass)
  - [x] Phase 2.4: Self-WaveCollapse intra-cluster EM (249 tests pass)
  - [x] Phase 2.5: Cluster representative selection (288 tests pass)
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
phase: 2
task: self_interference_stage1
substep: 6/6
last_action: ClusterRepSelector for medoid-based representative selection; 39 new tests; 288 total pass
next_action: implement llmap allpair CLI command
acceptance: llmap allpair reads FASTQ, clusters, outputs Parquet with cluster assignments
```

---

## Next task queue (FIFO)

1. ~~Phase 0.1: AlignmentRecord type + tests~~ ✅ done
2. ~~Phase 0.2: BucketPyramid data structure + serialization~~ ✅ done
3. ~~Phase 0.3: WaveState sparse CSR + GPU stub~~ ✅ done
4. ~~Phase 0.4: Synthetic IGH-locus generator + simulator wrappers~~ ✅ done
5. ~~Phase 1.1: ONNX Runtime CUDA-EP wired~~ ✅ done
6. ~~Phase 1.2: Test ONNX model + verified embedder inference~~ ✅ done
7. ~~Phase 1.3: BucketEmbedder + benchmark~~ ✅ done
8. ~~Phase 2.1: FAISS-GPU IndexIVFFlat wrapper~~ ✅ done
9. ~~Phase 2.2: Sparse k-NN extraction → similarity graph~~ ✅ done
10. ~~Phase 2.3: Leiden community detection~~ ✅ done
11. ~~Phase 2.4: Self-WaveCollapse intra-cluster EM~~ ✅ done
12. ~~Phase 2.5: Cluster representative selection~~ ✅ done
13. Phase 2.6: `llmap allpair` CLI command ← CURRENT
14. ... (continues per LLmap_SPEC.md)

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
| 3 | 2026-05-13 | n/a | implement wave_state | sparse CSR format + collapse/level mgmt; 21 tests pass |
| 4 | 2026-05-13 | n/a | implement synthetic_data_generator | IGH-locus generator + pbsim3/ART wrappers; 56 tests pass |
| 5 | 2026-05-13 | n/a | wire ONNX Runtime + FoundationEmbedder | CMake finds ONNX RT; FoundationEmbedder CPU/CUDA/TRT; 75 tests pass |
| 6 | 2026-05-13 | n/a | test ONNX model + embedder verification | Python model generator; real inference tests; 85 tests pass |
| 7 | 2026-05-13 | n/a | BucketEmbedder + throughput benchmark | bucket_embedder.{h,cpp}; 25 new tests; bench_embedder_throughput; 110 tests pass |
| 8 | 2026-05-13 | n/a | FAISS wrapper for ANN search | faiss_wrapper.{h,cpp}; CMake FAISS detection; 30 new tests; 140 total pass |
| 9 | 2026-05-13 | n/a | SimilarityGraph for sparse k-NN | similarity_graph.{h,cpp}; CSR format; 40 new tests; 180 total pass |
| 10 | 2026-05-13 | n/a | Leiden community detection | leiden_clustering.{h,cpp}; modularity-based partitioning; 34 new tests; 214 total pass |
| 11 | 2026-05-13 | n/a | Self-WaveCollapse intra-cluster EM | self_wavecollapse.{h,cpp}; EM-based read refinement within clusters; 35 new tests; 249 total pass |
| 12 | 2026-05-13 | n/a | Cluster representative selection | cluster_rep.{h,cpp}; medoid-based representative selection; 39 new tests; 288 total pass |

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
