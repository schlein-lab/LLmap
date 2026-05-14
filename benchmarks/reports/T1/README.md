# Benchmark Results: Task T1

Generated: 2026-05-14 10:42:51

## Task Description

Synthetic-truth WGS (1M simulated HiFi reads)

## Summary Table

| Tool | Replicates | Mapping Rate | Recall | Precision | F1 | Wallclock (s) | Peak RSS (MB) |
|------|------------|--------------|--------|-----------|-----|---------------|---------------|
| llmap | 1 | 46.2% | 38.7% | 62.3% | 47.7% | 31.9 | 305 |
| minimap2 | 1 | 91.8% | 100.0% | 100.0% | 100.0% | 9.1 | 252 |

## Detailed Statistics

### llmap

- **Replicates**: 1
- **Primary mapped**: 5.0K +/- 0.00
- **Mapping rate**: 46.2 +/- 0.0%
- **Recall**: 38.7 +/- 0.0%
- **Precision**: 62.3 +/- 0.0%
- **F1**: 47.7 +/- 0.0%
- **Wallclock**: 31.86 +/- 0.00 seconds
- **Peak RSS**: 305.17 +/- 0.00 MB

### minimap2

- **Replicates**: 1
- **Primary mapped**: 10.0K +/- 0.00
- **Mapping rate**: 91.8 +/- 0.0%
- **Recall**: 100.0 +/- 0.0%
- **Precision**: 100.0 +/- 0.0%
- **F1**: 100.0 +/- 0.0%
- **Wallclock**: 9.06 +/- 0.00 seconds
- **Peak RSS**: 252.11 +/- 0.00 MB
