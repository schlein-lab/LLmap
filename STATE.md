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
| Last successful iteration | 29 |
| Total iterations | 29 |

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
- [x] **Phase 2: Stage 1 Self-Interference** ✓
  - [x] Phase 2.1: FAISS-GPU IndexIVFFlat wrapper (140 tests pass)
  - [x] Phase 2.2: Sparse k-NN extraction → similarity graph (180 tests pass)
  - [x] Phase 2.3: Leiden community detection (214 tests pass)
  - [x] Phase 2.4: Self-WaveCollapse intra-cluster EM (249 tests pass)
  - [x] Phase 2.5: Cluster representative selection (288 tests pass)
  - [x] Phase 2.6: `llmap allpair` CLI + tests (348 tests pass)
- [x] **Phase 3: Stage 2 Reference WaveCollapse** ✓
  - [x] Phase 3.1: Reference index structure (375 tests pass)
  - [x] Phase 3.2: EM iteration kernel (CPU fallback) (400 tests pass)
  - [x] Phase 3.3: Collapse check + dropout (426 tests pass)
  - [x] Phase 3.4: Refinement (coarse→fine expansion) (443 tests pass)
  - [x] Phase 3.5: Member propagation (469 tests pass)
  - [x] Phase 3.6: Stage 2 pipeline orchestrator (493 tests pass)
- [x] **Phase 4: Classical Path + WFA2** ✓
  - [x] Phase 4.1: Minimizer index structure (529 tests pass)
  - [x] Phase 4.2: Chain extraction (552 tests pass)
  - [x] Phase 4.3: WFA2-lib FFI for gap-affine extension (587 tests pass)
  - [x] Phase 4.4: CPU-only fallback end-to-end (610 tests pass)
- [ ] **Phase 5: KILL-SWITCH VALIDATION** ★
  - [x] Phase 5.1: Kill-switch validation framework (631 tests pass)
  - [x] Phase 5.2: End-to-end synthetic validation (654 tests pass)
  - [x] Phase 5.3: Real reference integration infrastructure (684 tests pass)
  - [x] Phase 5.4: `validate-real` CLI + SLURM script (684 tests pass)
- [ ] Phase 6: Dual Output (BAM + Parquet)
- [ ] Phase 7: Claude Agent Integration
- [ ] Phase 8: Performance Optimization
- [ ] Phase 9: Single-Cell + Paralog Production

---

## Current task

```
phase: 5
task: kill_switch_validation
substep: 5/6
last_action: Phase 5.4 — added `llmap validate-real` CLI command + scripts/slurm_validate_real.sh; CLI parses --reference, --reads, --baseline, --truth, --k, --w, --min-identity, --gpu flags and calls RunRealReferenceValidation; SLURM script loads cuda/gcc modules and submits GPU job; 684 tests pass
next_action: Phase 5.5 — Submit Hummel SLURM job for GPU validation
  - SSH to Hummel and sync repo
  - Build LLmap with CUDA enabled on Hummel
  - Download hg38 chr14 IGH locus reference (~1MB)
  - Generate synthetic reads with ground truth
  - Run minimap2 baseline alignment
  - Submit SLURM job for GPU-enabled LLmap validation
  - Mark as awaiting_slurm_job_<id> if submitted
acceptance: GPU validation job completes with recall ≥ 99.5% of minimap2
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
13. ~~Phase 2.6: `llmap allpair` CLI command~~ ✅ done
14. ~~Phase 3.1: Reference index structure~~ ✅ done
15. ~~Phase 3.2: EM iteration kernel (CPU fallback)~~ ✅ done
16. ~~Phase 3.3: Collapse check + dropout~~ ✅ done
17. ~~Phase 3.4: Refinement (coarse→fine expansion)~~ ✅ done
18. ~~Phase 3.5: Member propagation~~ ✅ done
19. ~~Phase 3.6: Stage 2 pipeline orchestrator~~ ✅ done
20. ~~Phase 3 refactor: bucket_embedder.cpp split~~ ✅ done
21. ~~Phase 4.1: Minimizer index structure~~ ✅ done
22. ~~Phase 4.2: Chain extraction~~ ✅ done
23. ~~Phase 4.3: WFA2-lib FFI for gap-affine extension~~ ✅ done
24. ~~Phase 4.4: CPU-only fallback end-to-end~~ ✅ done
25. ~~Phase 5.1: Kill-switch validation framework~~ ✅ done
26. ~~Phase 5.2: End-to-end synthetic validation~~ ✅ done
27. ~~Phase 5.3: Real reference integration infrastructure~~ ✅ done
28. ~~Phase 5.4 refactor: real_reference.cpp split~~ ✅ done
29. ~~Phase 5.4 validate-real CLI + SLURM script~~ ✅ done
30. Phase 5.5: Submit Hummel SLURM job for GPU validation ← NEXT
31. ... (continues per LLmap_SPEC.md)

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
| 13 | 2026-05-13 | n/a | llmap allpair CLI tests + bugfix | test_llmap_cli.cpp (14 tests); fixed FastqReader::HasMore() peek EOF; 348 total pass |
| 14 | 2026-05-13 | n/a | Phase 3.1 ReferenceIndex structure | reference_index.{h,cpp}; Builder pattern; save/load serialization; spatial bucket lookup; 27 new tests; 375 total pass |
| 15 | 2026-05-13 | n/a | refactor faiss_wrapper + Phase 3.2 em_iterator | split faiss_wrapper.cpp (665 LOC) → 3 files; em_iterator.{h,cpp} for EM step CPU fallback; 25 new tests; 400 total pass; monolith count 10→9 |
| 16 | 2026-05-13 | n/a | refactor leiden_clustering + Phase 3.3 collapse_check | split leiden_clustering.cpp (624 LOC) → 3 files; collapse_check.{h,cpp} for convergence/dropout; 26 new tests; 426 total pass; monolith count 9→8 |
| 17 | 2026-05-13 | n/a | refactor allpair_pipeline + Phase 3.4 refinement | split allpair_pipeline.cpp (575 LOC) → 3 files; refinement.{h,cpp} for coarse→fine expansion; ChildIndex for parent→child mapping; 17 new tests; 443 total pass; monolith count 8→7 |
| 18 | 2026-05-13 | n/a | refactor cluster_rep + Phase 3.5 member_propagation | split cluster_rep.cpp (573 LOC) → 3 files; member_propagation.{h,cpp} for propagating rep positions to members; fixed allpair_pipeline build error; 26 new tests; 469 total pass; monolith count 7→6 |
| 19 | 2026-05-13 | n/a | refactor reference_index + Phase 3.6 stage2_pipeline | split reference_index.cpp (543 LOC) → 3 files; stage2_pipeline.{h,cpp} orchestrator wiring all Stage 2 components; 24 new tests; 493 total pass; monolith count 6→5; Phase 3 complete |
| 20 | 2026-05-13 | n/a | refactor bucket_embedder | split bucket_embedder.cpp (509 LOC) → 3 files (core, tokenize, embed); internal impl header for PIMPL; 493 tests pass; monolith count 6→5 |
| 21 | 2026-05-13 | n/a | Phase 4.1 minimizer_index + refactor | minimizer_index.{h,cpp} with k-mer extraction, Builder, query, serialization; split self_wavecollapse.cpp (495 LOC) → 3 files; split minimizer_index → 3 files; 36 new tests; 529 total pass; monolith count 5→4 |
| 22 | 2026-05-13 | n/a | Phase 4.2 chain extraction + refactor | chain.{h,cpp} with colinear DP chaining, anchor scoring, 0.9× best filter; split similarity_graph.cpp (488 LOC) → 3 files; 23 new tests; 552 total pass; monolith count 4→3 |
| 23 | 2026-05-13 | n/a | Phase 4.3 WFA2 aligner + refactor | wfa2_aligner.{h,cpp} with Gotoh gap-affine DP fallback; PIMPL for native WFA2-lib; Align/ExtendLeft/ExtendRight/AlignBatch APIs; split foundation_embedder.cpp (420 LOC) → 3 files; 35 new tests; 587 total pass; monolith count 3→2 |
| 24 | 2026-05-13 | n/a | Phase 4.4 classical_pipeline + refactor | classical_pipeline.{h,cpp} orchestrator wiring minimizer_index + chain + wfa2_aligner; seed-chain-extend end-to-end; split stage2_pipeline.cpp (415 LOC) → 3 files; 23 new tests; 610 total pass; monolith count 2→1; Phase 4 complete |
| 25 | 2026-05-13 | n/a | Phase 5.1 killswitch validation + refactor | killswitch.{h,cpp} validation framework; KillSwitchValidator class for lossless/position/origin checks; baseline comparison; split igh_locus_generator.cpp (410 LOC) → 2 files; 21 new tests; 631 total pass; monolith count 1→0 |
| 26 | 2026-05-13 | n/a | Phase 5.2 end-to-end synthetic validation | end_to_end.cpp with ClassicalAlignment→AlignmentRecord conversion; EndToEndConfig presets; RunEndToEndValidation orchestrator; position/timing/verdict metrics; 23 new tests; 654 total pass |
| 27 | 2026-05-13 | yes | Phase 5.3 real reference infrastructure | real_reference.{h,cpp} with RealReferenceConfig/RealGroundTruth/RealReferenceResult; BED/BAM parsing; SLURM job management; fasta_reader.{h,cpp}; SLURM template script; 30 new tests; 684 total pass |
| 28 | 2026-05-13 | n/a | Phase 5.4 refactor real_reference | split real_reference.cpp (611 LOC) → 3 files (parse, validate, slurm); internal header for shared utils; monolith count 1→0; 684 tests pass |
| 29 | 2026-05-13 | yes | Phase 5.4 validate-real CLI + SLURM | added `llmap validate-real` CLI command with arg parsing; scripts/slurm_validate_real.sh for GPU job submission; 684 tests pass |

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
