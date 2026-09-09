# Rat CBAS Validation Report

## Validation against the Igor reference implementation

Test statistics were compared sequence-by-sequence against the Igor reference output (`ratTestStats.txt`). Of 16,376 overlapping sequences, 100.0000% match within 1e-6, with a maximum absolute difference of 2.66e-07. This run found 572 significant sequences with k=29.

## Results

| Parameter | Value |
|---|---|
| Cohort | all_published (experiments 0-3, genotype 0, lesion known) |
| Subjects | 105 (55 control, 50 lesion) |
| Max sequence length | 6 |
| Criterion | 800 |
| block_aware | True |
| Resamples | 10,000 |
| Sequences evaluated | 16,378 |
| Significant | 572 (3.5%) |
| Control > Lesion | 264 |
| Lesion > Control | 308 |
| k (k-FWER) | 29 |
| Runtime | 12.3s |

## Timing Profile

| Stage | Time (s) | % Total |
|---|---|---|
| build_count_matrix | 0.27 | 2.2% |
| compute_test_stats | 0.01 | 0.1% |
| bootstrap | 4.97 | 40.4% |
| k_fwer | 7.06 | 57.4% |
| **TOTAL** | **12.31** | |

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

### Ranked ζ-values (Igor comparison)
![Ranked ζ-values](figures/ranked_gvalues.png)

## Top Significant Sequences

| Sequence | Direction | ζ-value | Decoded (arm, * = rewarded) |
|---|---|---|---|
| 3 | lesion>control | 0.0001 | 4 |
| 9 | control>lesion | 0.0001 | 4* |
| 8-9 | control>lesion | 0.0001 | 3* 4* |
| 7-8 | control>lesion | 0.0001 | 2* 3* |
| 9-8 | control>lesion | 0.0001 | 4* 3* |
| 8-7-8 | control>lesion | 0.0001 | 3* 2* 3* |
| 0-1 | control>lesion | 0.0001 | 1 2 |
| 7-8-1-8-9 | control>lesion | 0.0001 | 2* 3* 2 3* 4* |
| 8-7-8-1-8-9 | control>lesion | 0.0001 | 3* 2* 3* 2 3* 4* |
| 0-1-8 | control>lesion | 0.0001 | 1 2 3* |
| 7-0-1 | control>lesion | 0.0001 | 2* 1 2 |
| 4-3 | control>lesion | 0.0001 | 5 4 |
| 8-7-0-1 | control>lesion | 0.0001 | 3* 2* 1 2 |
| 0-3 | lesion>control | 0.0001 | 1 4 |
| 7-3 | lesion>control | 0.0001 | 2* 4 |
| 0-1-8-9 | control>lesion | 0.0001 | 1 2 3* 4* |
| 9-4 | control>lesion | 0.0001 | 4* 5 |
| 8-9-4 | control>lesion | 0.0001 | 3* 4* 5 |
| 5-4 | control>lesion | 0.0001 | 6 5 |
| 8-7-3 | lesion>control | 0.0001 | 3* 2* 4 |
| 4-3-8 | control>lesion | 0.0001 | 5 4 3* |
| 7-3-8 | lesion>control | 0.0001 | 2* 4 3* |
| 7-0-1-8 | control>lesion | 0.0001 | 2* 1 2 3* |
| 8-7-3-8 | lesion>control | 0.0001 | 3* 2* 4 3* |
| 8-7-0-1-8 | control>lesion | 0.0001 | 3* 2* 1 2 3* |
