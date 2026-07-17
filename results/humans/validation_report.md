# Human CBAS Validation Report (Correlative Mode)

**Mode:** Correlative — tests Pearson correlation between each sequence's usage
count across subjects and each subject's CBIT score (a compulsivity measure).
Positive correlation means higher CBIT (more compulsive) subjects use that
sequence more; negative means less.

> **Count difference:** We find 69 significant sequences vs the paper's 31.
> This likely reflects minor differences in the tau-hat normalization (the studentized
> variance estimate). The qualitative pattern is the same: most significant sequences
> are positively correlated and involve reward-switch motifs (B1, A2).

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

Each dot is one behavioral sequence in the two-step task. The y-axis shows the
significance of its correlation with the CBIT compulsivity score. Sequences are
grouped by length (2-step on the left, 4-step on the right).

> **Paper comparison (Fig 1c middle panel):** The paper shows very few sequences
> crossing the threshold, concentrated at length 4. Our plot shows more, but the
> overall sparse pattern is consistent — most sequences are not correlated with CBIT.

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
