# LLmap 96h Cross-Species Push

**Start**: 2026-05-15 ~21:00 CEST
**End**:   2026-05-19 ~21:00 CEST
**Driver**: ScheduleWakeup every 30min + system cron every 15min
**Goal**: Best mapping tool in the world across every species, with curated
SD/low-complexity/specific-loci knowledge per species, validated against
minimap2/bwa-mem2/winnowmap, deployed to losslessmap.com.

---

## Phase A — Cross-species knowledge expansion (0-24h)

Parallel-curated `knowledge/organisms/<org>/specific_loci/` + `regions/`
JSONs covering SDs, centromeres, telomeres, low-complexity, and species-
specific genomic disorders (where applicable). 10-15 parallel agents.

| Agent | Organism | Reference | Focus |
|---|---|---|---|
| A1 | mouse | mm39/GRCm39 | per-chr SDs (Bailey-mouse-SD), centromeres, MHC H2, Igh κ λ, olfactory, KRAB ZFP cluster |
| A2 | rat   | mRatBN7.2  | per-chr SDs, MHC RT1, Igh-rat |
| A3 | chimp | T2T-Pan | great-ape specific HSDs, AmpliconY-Pan, NPHP/IGH paralog families |
| A4 | gorilla | T2T-Gor | GORSAT, gorilla-specific SDs, GAGE-Y |
| A5 | orangutan + bonobo | Susie-PA + T2T-bonobo | LCR16a sister, orang Y |
| A6 | zebrafish | GRCz11 | hox cluster duplications, MHC zebrafish, vasa paralogs |
| A7 | drosophila | BDGP6.46 | heterochromatin Y-loops, NAHR @ 4 chromosome dot, X/Y NORs |
| A8 | C.elegans | WBcel235 | sperm-specific MSP duplicons, X-cluster Plg-1 |
| A9 | yeast (cerevisiae + pombe) | sacCer3 + ASM294v2 | rDNA repeats, Ty retrotransposon LTRs, telomeric Y' Y" |
| A10 | arabidopsis + rice + maize | TAIR10 + IRGSP-1.0 + B73 RefGen_v5 | NLR clusters, centromere/satDNA, transposable elements |
| A11 | viruses-DNA | HBV/HPV/EBV/KSHV/AdV refs | hairpins, inverted repeats, integration hotspots |
| A12 | viruses-RNA | SARS-CoV-2 + HIV-1 + Influenza H1-H9 | conserved RNA structures, recombination breakpoints |
| A13 | bacteria | E.coli K-12/MTB H37Rv/S.aureus N315 | rrn operons, IS-elements, phage islands, ICEs |
| A14 | synthetic stress | LLmap-synth-bench | extreme GC, low-complexity, tandems, palindromes |
| A15 | metagenomic panels | SILVA 16S + UNITE ITS + RefSoil | rRNA highly-conserved + species-discriminating regions |

Each agent writes `knowledge/organisms/<org>/specific_loci/**/*.json` with
the standard schema + mapping_hints, and `knowledge/organisms/<org>/regions/`
classifier rules. Goal: 100-500 JSONs per organism.

## Phase B — Cross-species benchmark harness (24-48h)

Wire each curated organism into `benchmarks/datasets/cache/<org>_t1..t10`
with:
- Per-organism synthetic reads (long + short, accurate + error-prone)
- Truth tables (ground-truth source position per read)
- 4-mapper comparison: llmap / minimap2 / bwa-mem2 / winnowmap
- Per-read accuracy + region-aware accuracy (recovery on SDs/telomeres/centromeres)
- BAM + RGN audit tags

`benchmarks/scripts/run_species_bench.sh <org> <tier>` produces:
- `benchmarks/reports/<org>/<mapper>/rep<N>/{alignments.bam,probabilities.csv,stats.json}`

## Phase C — Iterate (48-72h)

- For each species: identify failing regions, add classifier rules, tune
  per-region mapping_hints, re-run bench.
- LLM checkpoint hooks at AmbiguousChain for cross-species PSV resolution.
- Wave-particle adaptive λ(complexity) tuned per organism.
- Online learning of priors during streaming run (stochastic Bayes update).

## Phase D — Deploy to losslessmap.com (72-96h)

- Static benchmark site under `web/losslessmap.com/`
- Per-species pages: scoreboard, per-region heatmaps, per-mapper failures
- Aggregate cross-species scoreboard (rank by accuracy * speed * RAM)
- GitHub Pages or Cloudflare Pages deploy from `gh-pages` branch
- Daily snapshot tag (e.g. `bench-2026-05-16`, `bench-2026-05-17`, ...)

---

## State tracking

`PHASE_STATE.json` updated each wakeup with:
- current_phase: A/B/C/D
- completed_organisms: [list]
- failed_regions: [list]
- last_benchmark_score: dict
- pending_actions: [list]

## Iteration log

Per-wake entry appended to `CROSS_SPECIES_LOG.md`:
```
## Iter NNN — 2026-MM-DDTHH:MM CEST
- Phase: ...
- Done: ...
- Next: ...
```
