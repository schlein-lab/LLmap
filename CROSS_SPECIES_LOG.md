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
