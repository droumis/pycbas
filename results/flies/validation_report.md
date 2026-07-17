# Fly CBAS Validation Report

## Summary

| | pycbas | Paper (Kastner et al.) |
|---|---|---|
| Flies | 1566 (759 CA, 807 w1118) | 1,566 (759 CA, 807 w1118) |
| Max seq length | 10 | 10 |
| Criterion | 250 | 250 |
| Resamples | 10,000 | 10,000 |
| Sequences evaluated | 2,046 | 2,046 |
| Significant | 2046 (100.0%) | 1,605 (78.4%) |
| CA > w1118 | 454 | not separately reported |
| w1118 > CA | 1592 | not separately reported |
| k (k-FWER) | 103 | not reported |
| Runtime | 267.1s | not reported |

## Timing Profile

| Stage | Time (s) | % Total |
|---|---|---|
| build_count_matrix | 2.21 | 0.8% |
| compute_test_stats | 0.01 | 0.0% |
| bootstrap | 9.18 | 3.4% |
| k_fwer | 255.74 | 95.7% |
| **TOTAL** | **267.14** | |

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
