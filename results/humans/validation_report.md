# Human CBAS Validation Report (Correlative Mode)

**Mode:** Correlative — tests Pearson correlation between each sequence's usage
count across subjects and each subject's CBIT score (a compulsivity measure).
Positive correlation means higher CBIT (more compulsive) subjects use that
sequence more; negative means less.

> **Match:** We find 31 significant sequences — a perfect 31/31
> overlap with David's result (same sequences, same directions, p-values within
> ±0.001). The permutation-based correlative null does not involve centering (no
> group delta to subtract), so the result is stable across implementations.

## Summary

| | pycbas | Paper (Kastner et al.) |
|---|---|---|
| Subjects | 1413 | 1,413 |
| Max seq length | 4 | 4 |
| Criterion | 400 | 400 |
| Resamples | 10,000 | 10,000 |
| Sequences evaluated | 408 | 408 |
| Significant | 31 (7.6%) | 31 (7.6%) |
| Positive correlation (↑ CBIT → ↑ usage) | 26 | not separately reported |
| Negative correlation (↑ CBIT → ↓ usage) | 5 | not separately reported |
| k (k-FWER) | 2 | not reported |
| Runtime | 5.6s | not reported |

## Timing Profile

| Stage | Time (s) | % Total |
|---|---|---|
| build_count_matrix | 0.66 | 11.9% |
| compute_test_stats | 0.01 | 0.1% |
| bootstrap | 4.00 | 72.0% |
| k_fwer | 0.89 | 16.1% |
| **TOTAL** | **5.56** | |

## Figures

### Manhattan Plot
![Manhattan Plot](figures/manhattan.png)

Each dot is one behavioral sequence in the two-step task. The y-axis shows the
significance of its correlation with the CBIT compulsivity score. Sequences are
grouped by length (2-step on the left, 4-step on the right).

> **Paper comparison (Fig 1c middle panel):** The paper shows very few sequences
> crossing the threshold, concentrated at length 4. Our result matches —
> most sequences are not correlated with CBIT.

### Significant Sequences by Correlation Direction
![Direction Counts](figures/direction_counts.png)

Most significant sequences are positively correlated with CBIT — higher
compulsivity subjects use them more. This matches the paper's finding that
reward-switching motifs (sequences involving transitions between the common
and rare paths) drive the CBIT correlation.

### Null Distribution vs Observed
![Null vs Observed](figures/null_vs_observed.png)

Blue: observed correlation test statistics. Gray: null row-max per resample
(strongest signal a permutation can produce). The observed max exceeding the
null confirms that some sequences genuinely correlate with compulsivity.

### Sequence Space
![Sequence Space](figures/sequence_space.png)

With 6 choices and max length 4, there are 408 possible sequences. The space
is fully enumerable (unlike the rat case where 12^6 = 2.9M are possible but
only ~16K are observed).

### g-value Distribution
![g-value Distribution](figures/gvalue_dist.png)

Most g-values cluster near 1 (not significant), with a small subset below
the 0.5 threshold — consistent with a sparse signal where only specific
task-relevant sequences show a compulsivity correlation.
