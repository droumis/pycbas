# Rat CBAS Validation Report

## Summary

| | pycbas | Paper (Kastner et al.) |
|---|---|---|
| Rats | 85 (46 control, 39 lesion) | 85 (46 control, 39 lesion) |
| Max seq length | 6 | 6 |
| Criterion | 800 | 800 |
| Resamples | 10000 | 10,000 |
| Sequences evaluated | 16,500 | 24,342 |
| Significant | 178 (1.1%) | 409 (1.7%) |
| Control > Lesion | 91 | not separately reported |
| Lesion > Control | 87 | not separately reported |
| k (k-FWER) | 9 | not reported |
| Runtime | 8.2s | not reported |

## Timing Profile

| Stage | Time (s) | % Total |
|---|---|---|
| build_count_matrix | 0.22 | 2.7% |
| compute_test_stats | 0.00 | 0.1% |
| bootstrap | 3.69 | 45.2% |
| k_fwer | 4.25 | 52.1% |
| **TOTAL** | **8.17** | |

## Figures

### Manhattan Plot
![Manhattan Plot](figures/manhattan.png)

### Significant Sequences by Direction
![Direction Counts](figures/direction_counts.png)

### Null Distribution vs Observed
![Null vs Observed](figures/null_vs_observed.png)

### Sequence Space
![Sequence Space](figures/sequence_space.png)

### g-value Distribution
![g-value Distribution](figures/gvalue_dist.png)

## Top Significant Sequences

| Sequence | Direction | ζ-value | Decoded (arm, * = rewarded) |
|---|---|---|---|
| 4-3-8-7-8 | control>lesion | 0.0002 | 5 4 3* 2* 3* |
| 4-3-8-7-8-9 | control>lesion | 0.0003 | 5 4 3* 2* 3* 4* |
| 0-1-8-9-8-3 | control>lesion | 0.0003 | 1 2 3* 4* 3* 4 |
| 7-8-9-4-3 | control>lesion | 0.0003 | 2* 3* 4* 5 4 |
| 8-7-8-9-4-3 | control>lesion | 0.0003 | 3* 2* 3* 4* 5 4 |
| 9-4-3-8-7-8 | control>lesion | 0.0003 | 4* 5 4 3* 2* 3* |
| 7-8-9-4-3-8 | control>lesion | 0.0003 | 2* 3* 4* 5 4 3* |
| 7-3-8-3 | lesion>control | 0.0003 | 2* 4 3* 4 |
| 3-8-7-8-9-4 | control>lesion | 0.0003 | 4 3* 2* 3* 4* 5 |
| 8-7-8-1-0-1 | control>lesion | 0.0003 | 3* 2* 3* 2 1 2 |
| 7-8-1-0-1 | control>lesion | 0.0003 | 2* 3* 2 1 2 |
| 0-1 | control>lesion | 0.0004 | 1 2 |
| 8-1-8-1-8-1 | lesion>control | 0.0004 | 3* 2 3* 2 3* 2 |
| 8-9-4-3 | control>lesion | 0.0004 | 3* 4* 5 4 |
| 8-7-3-8-3 | lesion>control | 0.0004 | 3* 2* 4 3* 4 |
| 9-4-3 | control>lesion | 0.0005 | 4* 5 4 |
| 8-3-8-7-3-8 | lesion>control | 0.0005 | 3* 4 3* 2* 4 3* |
| 3-8-7-8-9 | control>lesion | 0.0008 | 4 3* 2* 3* 4* |
| 7-3-8 | lesion>control | 0.0008 | 2* 4 3* |
| 8-7-3-8 | lesion>control | 0.0008 | 3* 2* 4 3* |
| 1-8-9-8-3-4 | control>lesion | 0.0008 | 2 3* 4* 3* 4 5 |
| 8-4-8 | lesion>control | 0.0012 | 3* 5 3* |
| 8-3-8-7-8-9 | control>lesion | 0.0013 | 3* 4 3* 2* 3* 4* |
| 0-1-8-9-8 | control>lesion | 0.0015 | 1 2 3* 4* 3* |
| 7-3-8-7-8 | lesion>control | 0.0015 | 2* 4 3* 2* 3* |
