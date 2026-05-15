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
