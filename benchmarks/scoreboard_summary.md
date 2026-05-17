# LLmap cross-species scoreboard

Aggregated over **989** result rows across **19** organisms, **88** (organism, tier) cells, **4** mappers (bwa-mem2, llmap, minimap2, winnowmap).

## Head-to-head (LLmap vs other mappers)

- **vs minimap2** [all regions (overall+by_class)]: LLmap wins **15**, ties **67**, loses **301** of **383** cells.
- **vs minimap2** [region_class=SD]: LLmap wins **2**, ties **7**, loses **26** of **35** cells.
- **vs minimap2** [region_class=centromere]: LLmap wins **0**, ties **11**, loses **22** of **33** cells.
- **vs minimap2** [region_class=telomere]: LLmap wins **0**, ties **15**, loses **34** of **49** cells.
- **vs minimap2** [region_class=low_complexity]: LLmap wins **0**, ties **11**, loses **15** of **26** cells.
- **vs minimap2** [region_class=SV_DEL]: LLmap wins **3**, ties **0**, loses **27** of **30** cells.
- **vs minimap2** [region_class=SV_DUP]: LLmap wins **3**, ties **0**, loses **27** of **30** cells.
- **vs minimap2** [region_class=SV_INV]: LLmap wins **3**, ties **0**, loses **27** of **30** cells.

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
| great_apes | tier10 | NA | 0.769 | 1.000 | NA |
| great_apes | tier2 | NA | 0.737 | 1.000 | NA |
| great_apes | tier5 | NA | 0.731 | 1.000 | NA |
| great_apes | tier6 | NA | 0.766 | 1.000 | NA |
| homo_sapiens | T1 | NA | 0.477 | 1.000 | NA |
| homo_sapiens | T2 | NA | 0.445 | 1.000 | NA |
| maize | tier1 | NA | 0.647 | 1.000 | NA |
| maize | tier10 | NA | 0.672 | 1.000 | NA |
| maize | tier2 | NA | 0.643 | 1.000 | NA |
| maize | tier5 | NA | 0.639 | 1.000 | NA |
| maize | tier6 | NA | 0.673 | 1.000 | NA |
| metagenomic | tier1 | NA | 0.919 | 1.000 | NA |
| metagenomic | tier10 | NA | 0.851 | 0.000 | NA |
| metagenomic | tier2 | NA | 0.907 | 0.999 | NA |
| metagenomic | tier5 | NA | 0.903 | 0.999 | NA |
| metagenomic | tier6 | NA | 0.849 | 0.999 | NA |
| mouse | tier1 | NA | NA | 1.000 | NA |
| mouse | tier10 | NA | NA | 1.000 | NA |
| mouse | tier2 | NA | NA | 1.000 | NA |
| mouse | tier5 | NA | NA | 1.000 | NA |
| mouse | tier6 | NA | NA | 1.000 | NA |
| rat | tier1 | NA | 0.817 | 1.000 | NA |
| rat | tier10 | NA | 0.831 | 1.000 | NA |
| rat | tier2 | NA | 0.814 | 1.000 | NA |
| rat | tier5 | NA | 0.821 | 1.000 | NA |
| rat | tier6 | NA | 0.830 | 1.000 | NA |
| rice | tier1 | NA | 0.759 | 1.000 | NA |
| rice | tier10 | NA | 0.777 | 1.000 | NA |
| rice | tier2 | NA | 0.768 | 1.000 | NA |
| rice | tier5 | NA | 0.770 | 1.000 | NA |
| rice | tier6 | NA | 0.772 | 0.000 | NA |
| scerevisiae | tier1 | NA | 1.000 | 1.000 | NA |
| scerevisiae | tier10 | NA | 0.898 | 1.000 | NA |
| scerevisiae | tier2 | NA | 0.927 | 1.000 | NA |
| scerevisiae | tier5 | NA | 0.929 | 1.000 | NA |
| scerevisiae | tier6 | NA | 0.903 | 1.000 | NA |
| spombe | tier1 | NA | 0.925 | 1.000 | NA |
| spombe | tier10 | NA | 0.899 | 1.000 | NA |
| spombe | tier2 | NA | 0.934 | 1.000 | NA |
| spombe | tier5 | NA | 0.928 | 1.000 | NA |
| spombe | tier6 | NA | 0.904 | 1.000 | NA |
| synth_t1 | tier1 | NA | 0.940 | NA | NA |
| synthetic_stress | tier1 | NA | 0.998 | 0.999 | NA |
| synthetic_stress | tier10 | NA | 0.913 | 0.999 | NA |
| synthetic_stress | tier2 | NA | 0.942 | 0.999 | NA |
| synthetic_stress | tier5 | NA | 0.943 | 1.000 | NA |
| synthetic_stress | tier6 | NA | 0.911 | 0.999 | NA |
| viruses_dna | tier1 | NA | 0.958 | 0.999 | NA |
| viruses_dna | tier10 | NA | 0.921 | 0.999 | NA |
| viruses_dna | tier2 | NA | 0.957 | 1.000 | NA |
| viruses_dna | tier5 | NA | 0.960 | 1.000 | NA |
| viruses_dna | tier6 | NA | 0.924 | 1.000 | NA |
| viruses_rna | tier1 | NA | 0.002 | 0.998 | NA |
| viruses_rna | tier10 | NA | 0.001 | 0.995 | NA |
| viruses_rna | tier2 | NA | 0.000 | 0.996 | NA |
| viruses_rna | tier5 | NA | 0.001 | 0.997 | NA |
| viruses_rna | tier6 | NA | 0.002 | 0.995 | NA |
| zebrafish | tier1 | NA | 0.668 | 1.000 | NA |
| zebrafish | tier10 | NA | 0.703 | 1.000 | NA |
| zebrafish | tier2 | NA | 0.680 | 1.000 | NA |
| zebrafish | tier5 | NA | 0.678 | 1.000 | NA |
| zebrafish | tier6 | NA | 0.701 | 1.000 | NA |

## Top-5 LLmap wins (largest F1 gain over best other mapper)

| organism | tier | region_class | best_other | LLmap F1 | other F1 | delta |
|---|---|---|---|---:|---:|---:|
| metagenomic | tier10 | unique | minimap2 | 0.902 | 0.000 | +0.902 |
| rice | tier6 | SV_DUP | minimap2 | 0.864 | 0.000 | +0.864 |
| rice | tier6 | SV_DEL | minimap2 | 0.859 | 0.000 | +0.859 |
| metagenomic | tier10 | SV_DEL | minimap2 | 0.853 | 0.000 | +0.853 |
| metagenomic | tier10 | all | minimap2 | 0.851 | 0.000 | +0.851 |

## Top-5 LLmap regressions (Phase C tuning targets)

| organism | tier | region_class | best_other | LLmap F1 | other F1 | delta |
|---|---|---|---|---:|---:|---:|
| viruses_rna | tier1 | all | minimap2 | 0.002 | 0.998 | -0.996 |
| viruses_rna | tier1 | unique | minimap2 | 0.002 | 0.998 | -0.996 |
| viruses_rna | tier2 | all | minimap2 | 0.000 | 0.996 | -0.996 |
| viruses_rna | tier2 | unique | minimap2 | 0.000 | 0.996 | -0.996 |
| viruses_rna | tier5 | all | minimap2 | 0.001 | 0.997 | -0.996 |
