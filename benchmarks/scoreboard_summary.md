# LLmap cross-species scoreboard

Aggregated over **16** result rows.

## Head-to-head

- LLmap wins vs minimap2 on 0/2 (organism,tier) pairs in all regions.

## Top-5 LLmap improvements (largest F1 gain over best other)

| organism | tier | region_class | best_other | LLmap F1 | other F1 | delta |
|---|---|---|---|---:|---:|---:|
| homo_sapiens | T1 | all | minimap2 | 0.477 | 1.000 | -0.523 |
| homo_sapiens | T2 | all | minimap2 | 0.445 | 1.000 | -0.555 |

## Top-5 LLmap regressions (Phase C tuning targets)

| organism | tier | region_class | best_other | LLmap F1 | other F1 | delta |
|---|---|---|---|---:|---:|---:|
| homo_sapiens | T2 | all | minimap2 | 0.445 | 1.000 | -0.555 |
| homo_sapiens | T1 | all | minimap2 | 0.477 | 1.000 | -0.523 |
