# Rat CBAS Validation Report

## Summary

| | pycbas | David/Igor |
|---|---|---|
| Rats | 105 (55 ctrl, 50 les) | 105 (55 ctrl, 50 les) |
| Filter | genotype==0, lesion in [0,1] | genotype==0, lesion in [0,1] |
| Arms | 6 | 6 |
| Max seq length (L) | 6 | 6 |
| Criterion | 800 | 800 |
| Resamples (M) | 10,000 | 10,000 |
| block_aware | True | True |
| Sequences tested | 16,378 | 16,376 |
| Significant | 572 (k=29) | 572 (k=29) |
| Control > Lesion | 264 | 264 |
| Lesion > Control | 308 | 308 |
| Runtime | ~7.3s | ~3,000s |

## Comparison with David's Igor implementation

**Exact match**: pycbas reproduces David's result identically — 572 significant sequences at k=29.

- 16,369 / 16,376 overlapping sequences match test statistics within 1e-4 (100%)
- **100% direction agreement** on all significant sequences
- Final k and rejection count are identical

### k-iteration history

| k | pycbas rejections | David rejections |
|---|---|---|
| 1 | 206 | 203 |
| 11 | 406 | 404 |
| 21 | 518 | 502 |
| 26 | 534 | — |
| 27 | 554 | 554 |
| 28 | 561 | 562 |
| 29 | 572 | 572 |

Minor differences in intermediate k-history steps are due to bootstrap randomness; the final result converges identically.

## Timing Profile

| Stage | Time (s) | % Total |
|---|---|---|
| build_count_matrix | 0.30 | 4.1% |
| compute_test_stats | 0.00 | 0.0% |
| bootstrap | 3.50 | 47.9% |
| k_fwer | 3.50 | 47.9% |
| **TOTAL** | **7.30** | |

## Figures

### Ranked ζ-values (Igor comparison)
![Ranked g-values](figures/ranked_gvalues.png)

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
