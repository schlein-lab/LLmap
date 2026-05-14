# Benchmark Results: Task T2

Generated: 2026-05-14 10:42:51

## Task Description

Synthetic paralog stress (500k reads from duplicated regions)

## Summary Table

| Tool | Replicates | Mapping Rate | Recall | Precision | F1 | Wallclock (s) | Peak RSS (MB) |
|------|------------|--------------|--------|-----------|-----|---------------|---------------|
| llmap | 1 | 45.0% | 36.0% | 58.5% | 44.5% | 16.2 | 69 |
| minimap2 | 1 | 91.9% | 100.0% | 100.0% | 100.0% | 5.6 | 106 |

## Detailed Statistics

### llmap

- **Replicates**: 1
- **Primary mapped**: 2.4K +/- 0.00
- **Mapping rate**: 45.0 +/- 0.0%
- **Recall**: 36.0 +/- 0.0%
- **Precision**: 58.5 +/- 0.0%
- **F1**: 44.5 +/- 0.0%
- **Wallclock**: 16.19 +/- 0.00 seconds
- **Peak RSS**: 68.80 +/- 0.00 MB

### minimap2

- **Replicates**: 1
- **Primary mapped**: 5.0K +/- 0.00
- **Mapping rate**: 91.9 +/- 0.0%
- **Recall**: 100.0 +/- 0.0%
- **Precision**: 100.0 +/- 0.0%
- **F1**: 100.0 +/- 0.0%
- **Wallclock**: 5.63 +/- 0.00 seconds
- **Peak RSS**: 105.71 +/- 0.00 MB
