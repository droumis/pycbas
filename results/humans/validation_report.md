# Human CBAS Validation Report (Correlative Mode)

## Summary

| | pycbas | Paper (Kastner et al.) |
|---|---|---|
| Subjects | 1413 | 1,413 |
| Max seq length | 4 | 4 |
| Criterion | 400 | 400 |
| Resamples | 10,000 | 10,000 |
| Sequences evaluated | 408 | 408 |
| Significant | 69 (16.9%) | 31 (7.6%) |
| Positive correlation (↑ CBIT → ↑ usage) | 51 | not separately reported |
| Negative correlation (↑ CBIT → ↓ usage) | 18 | not separately reported |
| k (k-FWER) | 4 | not reported |
| Runtime | 4.1s | not reported |

## Notes

- **Mode:** Correlative — tests Pearson correlation between each sequence's usage
  count across subjects and each subject's CBIT score (a compulsivity measure).
- **Interpretation:** Positive correlation means higher CBIT (more compulsive)
  subjects use that sequence more. Negative means less.

## Timing Profile

| Stage | Time (s) | % Total |
|---|---|---|
| build_count_matrix | 0.72 | 17.8% |
| compute_test_stats | 0.01 | 0.2% |
| bootstrap | 2.84 | 69.6% |
| k_fwer | 0.51 | 12.4% |
| **TOTAL** | **4.08** | |

## Figures

### Manhattan Plot
![Manhattan Plot](figures/manhattan.png)

### Significant Sequences by Correlation Direction
![Direction Counts](figures/direction_counts.png)

### Null Distribution vs Observed
![Null vs Observed](figures/null_vs_observed.png)

### Sequence Space
![Sequence Space](figures/sequence_space.png)

### g-value Distribution
![g-value Distribution](figures/gvalue_dist.png)
