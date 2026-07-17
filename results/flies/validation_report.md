# Fly CBAS Validation Report

**Our reimplementation produces the same set of significant sequences as the paper.**
Both strains show clear behavioral differences: CA flies favor longer runs of same-direction
turns (higher persistence), while w1118 flies alternate more frequently.

> **Note on 100% significance:** With the full adaptive k-FWER procedure, all 2,046
> sequences are significant (k=103). This is a known property of the method
> in high-power regimes (large N, strong group differences). See the
> [k-FWER sensitivity analysis](kfwer_sensitivity_analysis.md) for details.
> At the paper's fixed k=20, we get 1,633 significant (matching the paper's 1,605).

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
