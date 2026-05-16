# LLmap cross-species scoreboard

Aggregated over **347** result rows across **10** organisms, **34** (organism, tier) cells, **4** mappers (bwa-mem2, llmap, minimap2, winnowmap).

## Head-to-head (LLmap vs other mappers)

- **vs minimap2** [all regions (overall+by_class)]: LLmap wins **5**, ties **64**, loses **66** of **135** cells.
- **vs minimap2** [region_class=SD]: LLmap wins **1**, ties **7**, loses **6** of **14** cells.
- **vs minimap2** [region_class=centromere]: LLmap wins **0**, ties **11**, loses **5** of **16** cells.
- **vs minimap2** [region_class=telomere]: LLmap wins **0**, ties **12**, loses **5** of **17** cells.
- **vs minimap2** [region_class=low_complexity]: LLmap wins **0**, ties **11**, loses **2** of **13** cells.
- **vs minimap2** [region_class=SV_DEL]: LLmap wins **1**, ties **0**, loses **8** of **9** cells.
- **vs minimap2** [region_class=SV_DUP]: LLmap wins **1**, ties **0**, loses **8** of **9** cells.
- **vs minimap2** [region_class=SV_INV]: LLmap wins **1**, ties **0**, loses **8** of **9** cells.

## Overall F1 by organism x tier (region_class=all)

| organism | tier | bwa-mem2 | llmap | minimap2 | winnowmap |
|---|---|---:|---:|---:|---:|
| arabidopsis | tier1 | NA | 0.999 | 1.000 | NA |
| arabidopsis | tier10 | NA | 0.999 | 1.000 | NA |
| arabidopsis | tier2 | NA | 1.000 | 1.000 | NA |
| arabidopsis | tier5 | NA | 1.000 | 1.000 | NA |
| arabidopsis | tier6 | NA | 0.857 | 1.000 | NA |
| bacteria | tier1 | NA | 1.000 | 1.000 | NA |
| bacteria | tier10 | NA | 0.998 | 1.000 | NA |
| bacteria | tier2 | NA | 1.000 | 1.000 | NA |
| bacteria | tier5 | NA | 1.000 | 1.000 | NA |
| bacteria | tier6 | NA | 0.999 | 1.000 | NA |
| celegans | tier1 | NA | 1.000 | 1.000 | NA |
| celegans | tier10 | NA | 0.998 | 1.000 | NA |
| celegans | tier2 | NA | 1.000 | 1.000 | NA |
| celegans | tier5 | NA | 1.000 | 1.000 | NA |
| celegans | tier6 | NA | 0.892 | 1.000 | NA |
| drosophila | tier1 | NA | 1.000 | 1.000 | NA |
| drosophila | tier10 | NA | 0.861 | 1.000 | NA |
| drosophila | tier2 | NA | 1.000 | 1.000 | NA |
| drosophila | tier5 | NA | 1.000 | 1.000 | NA |
| drosophila | tier6 | NA | 0.999 | 0.197 | NA |
| great_apes | tier1 | NA | 0.700 | 1.000 | NA |
| great_apes | tier10 | NA | NA | 1.000 | NA |
| great_apes | tier2 | NA | 0.737 | 1.000 | NA |
| great_apes | tier5 | NA | 0.731 | 1.000 | NA |
| great_apes | tier6 | NA | 0.766 | 1.000 | NA |
| homo_sapiens | T1 | NA | 0.477 | 1.000 | NA |
| homo_sapiens | T2 | NA | 0.445 | 1.000 | NA |
| scerevisiae | tier1 | NA | 1.000 | 1.000 | NA |
| synth_t1 | tier1 | NA | 0.940 | NA | NA |
| synthetic_stress | tier1 | NA | 0.998 | 0.999 | NA |

## Top-5 LLmap wins (largest F1 gain over best other mapper)

| organism | tier | region_class | best_other | LLmap F1 | other F1 | delta |
|---|---|---|---|---:|---:|---:|
| drosophila | tier6 | SV_INV | minimap2 | 0.995 | 0.174 | +0.820 |
| drosophila | tier6 | SD | minimap2 | 1.000 | 0.197 | +0.803 |
| drosophila | tier6 | all | minimap2 | 0.999 | 0.197 | +0.802 |
| drosophila | tier6 | SV_DUP | minimap2 | 1.000 | 0.207 | +0.793 |
| drosophila | tier6 | SV_DEL | minimap2 | 0.999 | 0.208 | +0.791 |

## Top-5 LLmap regressions (Phase C tuning targets)

| organism | tier | region_class | best_other | LLmap F1 | other F1 | delta |
|---|---|---|---|---:|---:|---:|
| homo_sapiens | T2 | all | minimap2 | 0.445 | 1.000 | -0.555 |
| homo_sapiens | T1 | all | minimap2 | 0.477 | 1.000 | -0.523 |
| great_apes | tier1 | telomere | minimap2 | 0.529 | 1.000 | -0.471 |
| great_apes | tier1 | unique | minimap2 | 0.681 | 1.000 | -0.319 |
| great_apes | tier1 | SD | minimap2 | 0.692 | 1.000 | -0.308 |
