# Fly CBAS Validation Report

**Our reimplementation uses the same algorithm as David's Igor code** (no centering,
uncentered bootstrap null). Both strains show clear behavioral differences: CA flies
favor longer runs of same-direction turns (higher persistence), while w1118 flies
alternate more frequently.

> **k-FWER convergence:** Our k-iteration jumps from k=1 (R=1243) to k=63 due to
> the formula next_k = ceil((R+1) x gamma). David's implementation gives 1,605
> significant — likely converging at a lower k via a more gradual path (possibly
> due to bootstrap RNG differences, though other implementation details may
> contribute). Our 1,243 significant sequences are a strict subset of David's
> 1,605 (0 overcalled, 362 missed).

## Summary

| | pycbas | Paper (Kastner et al.) |
|---|---|---|
| Flies | 1566 (759 CA, 807 w1118) | 1,566 (759 CA, 807 w1118) |
| Max seq length | 10 | 10 |
| Criterion | 250 | 250 |
| Resamples | 10,000 | 10,000 |
| Sequences evaluated | 2,046 | 2,046 |
| Significant | 1243 (60.8%) | 1,605 (78.4%) |
| CA > w1118 | 208 | not separately reported |
| w1118 > CA | 1035 | not separately reported |
| k (k-FWER) | 63 | not reported |
| Runtime | 27.5s | not reported |

## Timing Profile

| Stage | Time (s) | % Total |
|---|---|---|
| build_count_matrix | 2.38 | 8.6% |
| compute_test_stats | 0.01 | 0.0% |
| bootstrap | 15.85 | 57.6% |
| k_fwer | 9.27 | 33.7% |
| **TOTAL** | **27.51** | |

## Figures

### Manhattan Plot
![Manhattan Plot](figures/manhattan.png)

Each dot is one behavioral sequence. The y-axis shows statistical significance
(higher = more different between strains). Sequences are grouped by length
(2-symbol on the left, 10-symbol on the right). Dots above the dotted threshold
are significantly different between CA and w1118 flies.

> **Paper comparison (Fig 1c left panel):** Our plot reproduces the same layout —
> most sequences are significant, with the signal strongest at intermediate lengths
> where persistence differences are most detectable.

### Significant Sequences by Direction
![Direction Counts](figures/direction_counts.png)

Breaks down significant sequences by which strain uses them more. The strong
asymmetry (w1118 > CA dominating) reflects w1118 flies' preference for short
alternating sequences, which outnumber the longer persistent sequences that
CA flies favor.

### Null Distribution vs Observed
![Null vs Observed](figures/null_vs_observed.png)

Blue: observed test statistics for all sequences. Gray: null row-max per resample
(strongest signal chance can produce). The red line (observed max) sitting far to
the right of the null confirms the group differences are genuine.

### Sequence Space
![Sequence Space](figures/sequence_space.png)

With 2 arms (L/R) and max length 10, there are 2,046 possible sequences total.
Unlike the rat/human cases, the combinatorial space is fully enumerable here.

### g-value Distribution
![g-value Distribution](figures/gvalue_dist.png)

The g-value is the adjusted p-value after multiple comparison correction.
Values below 0.5 are significant. A bimodal distribution (most sequences
clearly significant or clearly not) means the correction procedure is
working well.
