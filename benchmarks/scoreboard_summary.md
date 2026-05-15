# LLmap cross-species scoreboard

Aggregated over **234** result rows across **7** organisms, **23** (organism, tier) cells, **4** mappers (bwa-mem2, llmap, minimap2, winnowmap).

## Head-to-head (LLmap vs other mappers)

- **vs minimap2** [all regions (overall+by_class)]: LLmap wins **0**, ties **59**, loses **26** of **85** cells.
- **vs minimap2** [region_class=SD]: LLmap wins **0**, ties **6**, loses **0** of **6** cells.
- **vs minimap2** [region_class=centromere]: LLmap wins **0**, ties **10**, loses **0** of **10** cells.
- **vs minimap2** [region_class=telomere]: LLmap wins **0**, ties **11**, loses **0** of **11** cells.
- **vs minimap2** [region_class=low_complexity]: LLmap wins **0**, ties **10**, loses **1** of **11** cells.
- **vs minimap2** [region_class=SV_DEL]: LLmap wins **0**, ties **0**, loses **4** of **4** cells.
- **vs minimap2** [region_class=SV_DUP]: LLmap wins **0**, ties **0**, loses **4** of **4** cells.
- **vs minimap2** [region_class=SV_INV]: LLmap wins **0**, ties **0**, loses **4** of **4** cells.

## Overall F1 by organism x tier (region_class=all)

| organism | tier | bwa-mem2 | llmap | minimap2 | winnowmap |
|---|---|---:|---:|---:|---:|
| arabidopsis | tier1 | NA | 0.999 | 1.000 | NA |
| arabidopsis | tier10 | NA | 0.999 | 1.000 | NA |
| arabidopsis | tier2 | NA | 1.000 | 1.000 | NA |
| arabidopsis | tier5 | NA | 1.000 | 1.000 | NA |
| arabidopsis | tier6 | NA | NA | 1.000 | NA |
| bacteria | tier1 | NA | 1.000 | 1.000 | NA |
| bacteria | tier10 | NA | 0.998 | 1.000 | NA |
| bacteria | tier2 | NA | 1.000 | 1.000 | NA |
| bacteria | tier5 | NA | 1.000 | 1.000 | NA |
| bacteria | tier6 | NA | 0.999 | 1.000 | NA |
| celegans | tier1 | NA | 1.000 | 1.000 | NA |
| celegans | tier10 | NA | 0.998 | 1.000 | NA |
| celegans | tier2 | NA | 1.000 | 1.000 | NA |
| celegans | tier5 | NA | 1.000 | 1.000 | NA |
| celegans | tier6 | NA | NA | 1.000 | NA |
| drosophila | tier1 | NA | 1.000 | 1.000 | NA |
| drosophila | tier2 | NA | 1.000 | 1.000 | NA |
| drosophila | tier5 | NA | NA | 1.000 | NA |
| drosophila | tier6 | NA | NA | 0.197 | NA |
| homo_sapiens | T1 | NA | 0.477 | 1.000 | NA |
| homo_sapiens | T2 | NA | 0.445 | 1.000 | NA |
| scerevisiae | tier1 | NA | 1.000 | 1.000 | NA |
| synthetic_stress | tier1 | NA | 0.998 | 0.999 | NA |

## Top-5 LLmap wins (largest F1 gain over best other mapper)

| organism | tier | region_class | best_other | LLmap F1 | other F1 | delta |
|---|---|---|---|---:|---:|---:|
| synthetic_stress | tier1 | synth_stress | minimap2 | 0.857 | 0.857 | +0.000 |
| scerevisiae | tier1 | unique | minimap2 | 1.000 | 1.000 | +0.000 |
| scerevisiae | tier1 | telomere | minimap2 | 1.000 | 1.000 | +0.000 |
| scerevisiae | tier1 | synth_stress | minimap2 | 1.000 | 1.000 | +0.000 |
| scerevisiae | tier1 | all | minimap2 | 1.000 | 1.000 | +0.000 |

## Top-5 LLmap regressions (Phase C tuning targets)

| organism | tier | region_class | best_other | LLmap F1 | other F1 | delta |
|---|---|---|---|---:|---:|---:|
| homo_sapiens | T2 | all | minimap2 | 0.445 | 1.000 | -0.555 |
| homo_sapiens | T1 | all | minimap2 | 0.477 | 1.000 | -0.523 |
| celegans | tier10 | SV_INV | minimap2 | 0.991 | 1.000 | -0.009 |
| bacteria | tier10 | SV_INV | minimap2 | 0.993 | 1.000 | -0.007 |
| arabidopsis | tier10 | SV_INV | minimap2 | 0.996 | 1.000 | -0.004 |
