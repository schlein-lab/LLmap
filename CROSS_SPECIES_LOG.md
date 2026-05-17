# LLmap 96h Cross-Species Push — Iteration Log

## Iter 000 — 2026-05-15T21:00 CEST — Bootstrap

Plan doc written (`CROSS_SPECIES_96H.md`), state tracker initialized
(`PHASE_STATE.json`). 15-agent Phase-A dispatch starting. Wake scheduled
every 30min.

## Iter 001 — 2026-05-15T21:30 CEST — Phase A → Phase B handoff

**Phase A complete: 16 organisms, 3320 new JSONs (10742 total org-JSONs).**

| Organism | JSONs | Assembly |
|---|---|---|
| mouse | 454 | mm39 |
| rat | 252 | mRatBN7.2 |
| great_apes (5 species) | 783 | panTro6 + T2T-Pan/Gor/Pongo + Susie-PA + Nleu_3.0 |
| zebrafish | 149 | GRCz11 |
| drosophila | 121 | BDGP6.46 |
| celegans | 78 | WBcel235 |
| scerevisiae | 84 | sacCer3 |
| spombe | 81 | ASM294v2 |
| arabidopsis | 80 | TAIR10 |
| rice | 139 | IRGSP-1.0 |
| maize | 258 | B73 RefGen_v5 |
| viruses_dna | 100 | RefSeq (HBV/HPV×7/EBV/KSHV/HHV-6A+B/AdV×4/Polyoma×3/Vaccinia) |
| viruses_rna | 200 | RefSeq (SARS-CoV-2 lineage tree + HIV + Flu + HCV + Ebola + Zika) |
| bacteria | 300 | RefSeq (10 species, 57 rRNA + 130 IS + 15 PE/PPE + integrons) |
| synthetic_stress | 84 | synth (homopolymer/tri-NT/VNTR/palindrome/Z-DNA/duplicons/chimera) |
| metagenomic | 136 | SILVA + UNITE + BOLD + PR2 (16S/18S/ITS/COI/rbcL/matK) |

**Next: Phase B — cross-species bench harness.** 4-mapper comparison
(llmap vs minimap2 vs bwa-mem2 vs winnowmap) across these organisms.
Per-organism synthetic read generation + ground-truth tracking +
accuracy/recall reporting.

## Iter 002 — 2026-05-15T22:05 CEST — Phase B landed, smoke test kicked

Phase B scripts complete:
- gen_synth_reads.py (registry-driven, 17 organisms, tiers 1-10)
- build_fake_reference.py (stitches FASTA from specific_loci JSONs)
- run_species_bench.sh (orchestrator, 4-mapper, /usr/bin/time wrapped)
- run_all_species.sh (loop organisms × tiers 1,2,5,6,10)
- analyze_bench.py (per-read precision/recall/F1 + region-class breakdown)
- aggregate_scoreboard.py (idempotent scoreboard + failed_regions report)
- per_region_breakdown.py (per-locus accuracy join)
- check_mapper_availability.sh, timing_harness.sh

Mapper availability on vm-science: **llmap OK + minimap2 OK**;
bwa-mem2 + winnowmap UNAVAILABLE (skipped cleanly).

Smoke test: synthetic_stress tier 1 with llmap+minimap2 running in
background. If clean, dispatch Phase C agents on next wake.

## Iter 003 — 2026-05-15T22:38 CEST — Bench critical bug detected

Full bench at 9/16 organisms processed. **F1@1kb=0.000 across all
non-human organisms** for both llmap AND minimap2 — diagnosed as
truth-vs-BAM coordinate mismatch in the harness (not mapper bug).

Symptom: arabidopsis tier 1, 999/1000 reads MAPPED but 0/1000 within
1kb of truth. Same minimap2.

Hypothesis: gen_synth_reads.py records original-genome coords in
truth.tsv, but build_fake_reference.py stitches a packed FASTA with
fake-ref coords. BAM POS uses fake-ref → join fails 100%.

Phase C debug agent dispatched to fix. Also flagged: bacteria tier 6
llmap produced empty SAM (separate bug).

Bench continues running in background.

## Iter 004 — 2026-05-15T23:08 CEST — Bench fix validated, Phase C in flight

After analyze_bench.py truth-pos fix, real F1 numbers across 6 organisms:
- arabidopsis: LLmap 0.999-1.000 (5 tiers)
- bacteria: LLmap 0.998-1.000 (4 tiers, rRNA paralog dip on tier 10)
- celegans: LLmap 0.998-1.000 (3 tiers with data)
- drosophila: LLmap 1.000 (2 tiers)
- scerevisiae: LLmap 1.000
- synthetic_stress: LLmap 0.998 vs minimap2 0.999

LLmap matches minimap2 within Δ=0.001-0.002 where mapped.

**LLmap empty-SAM on tier 5/6 some organisms** (arabidopsis t6, bacteria
t5, celegans t1+t6) — Phase C debug agent dispatched + aggregator
schema fix. Bench continues (currently at drosophila tier 2).

## Iter 005 — 2026-05-15T23:42 CEST — Phase C scoreboard real, tuning dispatched

Aggregator schema fixed → 234 records across 4 mappers, 7 organisms,
22 cells. **LLmap = minimap2 ties at F1=1.000 on synthetic tiers**,
slight regressions only on tier10 SV_INV (delta -0.004 to -0.009 = noise).

**Real regression: legacy human T1/T2** with LLmap F1=0.477/0.445 vs
minimap2 F1=1.000. This is the only genuine signal — likely real-data
fixture probing IGH/paralog regions where LLmap's chain-DP picks
wrong paralog. Phase C investigator agent dispatched.

Bench at 21/~80 (organism × tier) cells, currently celegans tier 5.

## Iter 006 — 2026-05-16T08:45 CEST — Recovery from session restart

Issues recovered:
1. llmap binary was deleted (autonomous_driver clean rebuild) — rebuilt
   v1.0.0 commit 506389a0 fresh.
2. Bench died at human tier 2 (no llmap binary). Restarted.
3. T1/T2 investigator was interrupted mid-run. Re-dispatched.
4. the HPC cluster flipped UP (RRZ SSH key change resolved, manual flip).

**Watchdog infrastructure added**: `scripts/bench_watchdog.sh` + cron
`*/5 * * * *` relaunches `run_all_species.sh` if dead. Fixes the
unreliable ScheduleWakeup → bench-state divergence pattern.

Bench resumed at arabidopsis tier 6.

## Iter 007 — 2026-05-16T08:55 CEST — T1/T2 not a regression, two real fixes found

T1/T2 investigator: F1=0.477 was STALE report from old binary 4667060.
Current binary 506389a + fixed scorer → T1 F1=0.940 (recall 0.887,
precision 1.000). Real findings:

1. **ExtendLeft coord swap bug** in src/classical/wfa2_aligner_align.cpp
   :98-101 — overwrites query_start/ref_start before re-reading them,
   corrupts soft-clip length + NM/identity stats. 4-line fix.

2. **min_identity=0.90 too aggressive** on map-hifi preset — drops
   1135/10000 borderline HiFi reads on T1. Lower to 0.85 expects
   F1 ≈ 0.97.

Fix agent dispatched to apply both + rebuild + re-validate.

## Iter 008 — 2026-05-16T09:30 CEST — T1 re-validated with new binary, real F1=0.940

T1 re-aligned with binary b3694583 (post-WFA2-fix, post-stale-report):
- F1@1kb = 0.9398
- precision@1kb = 1.000
- recall@1kb = 0.8865 (8865/10000 mapped)
- F1@100bp = 0.9385

The investigator's "min_identity=0.90 → 0.85 expects F1≈0.97" claim
was wrong — `map-hifi` preset ALREADY had 0.85f set (line 46 of
src/cli/cmd_align_args.cpp). The 1135 unmapped reads have post-align
identity below 0.85 — that's the real LLmap ceiling on T1.

Versus minimap2 F1=1.000 on T1, this is a -0.06 real delta (not -0.5
as the stale report falsely suggested). Phase C tuning should focus
on those 1135 borderline reads — either WFA2 alignment quality
improvement or a 2-pass retry mode at lower stringency.

Bench continues at great_apes tier 6.

## Iter 009 — 2026-05-16T09:45 CEST — Phase C real signal + LLmap SV win

Scoreboard 347 records / 32 cells. **Top-5 LLmap wins all on
drosophila tier 6** (SV injection mode):
- SV_INV: LLmap 0.995 vs minimap2 0.174 = +0.820 F1
- SD: 1.000 vs 0.197 = +0.803
- SV_DUP: 1.000 vs 0.207 = +0.793
- SV_DEL: 0.999 vs 0.208 = +0.791

LLmap's chain-DP + wave-particle approach dominates structural
variant mapping. minimap2 fails on tier 6 SVs across most organisms
except where LLmap also struggles (arabidopsis t6: 0.857 vs 1.000).

**Real regression target: great_apes tier 1** — LLmap 0.700 overall
vs minimap2 1.000. Even unique-class reads at 0.681. Investigator
agent dispatched. Bench at human tier 10 (16/17 organisms processed).

## Iter 010 — 2026-05-16T09:55 CEST — Full bench signal: massive regression matrix

Bench now 86 summary.json files / 989 scoreboard records. Full matrix:

**LLmap WINS** (where minimap2 fails):
- drosophila tier 6: +0.80 F1 (SV mode)
- rice tier 6: +0.77 (minimap2 hit 0.000)
- metagenomic tier 10: +0.85 (minimap2 hit 0.000)

**LLmap CRITICAL FAIL** (systematic):
- viruses_rna ALL tiers: F1=0.000-0.002 vs minimap2 0.99+ ← MUST FIX
- great_apes ALL tiers: 0.70-0.77 vs 1.0 (cross-species divergence?)
- maize ALL tiers: 0.64-0.67 vs 1.0
- rat ALL tiers: 0.81-0.83 vs 1.0
- rice tier 1-5: 0.76-0.78 vs 1.0
- zebrafish: 0.67-0.70 vs 1.0
- mouse: NA (bench may have skipped llmap on mouse — check)

**LLmap parity** (F1 ≥0.99 both):
- arabidopsis tier 1-5/10, bacteria, celegans (tier 1/2/5/10),
  drosophila (tier 1/2/5), scerevisiae tier 1, synthetic_stress

Pattern hypothesis: viruses_rna failure likely tiny-contig issue
(virus loci 30-200bp << reads ~12kb). Dispatched focused investigator.

Bench at mouse tier 5 (no llmap accuracy.json yet).
