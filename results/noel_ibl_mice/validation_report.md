# Noel IBL CBAS: WT vs Autism Model Mice

## Summary

CBAS applied to the Noel et al. (2025) IBL two-alternative forced choice dataset,
comparing wild-type (C57BL/6J) mice against three autism mouse models.

**Bottom line:** After controlling for exposure (criterion truncation), no meaningful
sequential structure differences survive. The original uncontrolled analysis (21-29
unidirectional hits) was a confound of trial count and accuracy differences. A
toward/against encoding that removes absolute side also yields no signal.

### Dataset

Source: Noel et al., "A common computational and neural anomaly across mouse
models of autism", *Nature Neuroscience* 2025.

| Group | n | Total Trials | Criterion-Truncated |
|-------|---|------|------|
| WT (C57BL/6J) | 15 | 173,917 | 4,000 × 15 |
| Cntnap2 | 17 (+2 excluded) | 230,700 | 4,000 × 17 |
| Fmr1 | 15 (+2 excluded) | 154,133 | 4,000 × 15 |
| Shank3B | 14 (+4 excluded) | 143,305 | 4,000 × 14 |

61 of 69 animals retained at criterion=4000 (8 excluded for insufficient trials).
All trials included: all contrasts, all blocks, no-go excluded.

### Encoding: Choice × Correctness (L/R)

4-symbol alphabet (A=4): choice × correctness. In a 2AFC task, outcome is
fully determined by choice and stimulus side.

| Symbol | Choice | Outcome |
|--------|--------|---------|
| L✓ (0) | Left | Correct |
| L✗ (1) | Left | Error |
| R✗ (2) | Right | Error |
| R✓ (3) | Right | Correct |

### Parameters

| Parameter | Value |
|-----------|-------|
| Alphabet (A) | 4 |
| Max length (L) | 6 |
| Criterion | 4,000 trials/animal |
| Resamples (M) | 10,000 |
| Alpha | 0.5 (median FDP) |
| Gamma | 0.05 |
| Contingency | all blocks |
| Data source | all_contrasts/ txt files (includes 0% contrast) |

### Results: Controlled Analysis

| Comparison | Sequences | Significant | WT > Model | Model > WT | k |
|---|---|---|---|---|---|
| WT vs Cntnap2 | 4,792 | 6 | 1 | 5 | 1 |
| WT vs Fmr1 | 4,761 | 5 | 2 | 3 | 1 |
| WT vs Shank3B | 4,793 | 3 | 1 | 2 | 1 |

Group accuracies after including all contrasts and truncating: WT 81.4%, Cntnap2
83.0%, Fmr1 82.9%, Shank3B 82.3%. Tightly matched — no accuracy confound.

Results are **bidirectional** (hits go both ways), consistent with noise-level
findings rather than systematic group differences.

## The Confound

The original uncontrolled analysis (L=8, no criterion, non-zero contrast only)
produced 21-29 **unidirectional** hits (all WT > model). The unidirectionality is
the signature of an exposure/accuracy confound:

1. **Exposure confound**: Without criterion truncation, animals with more trials have
   more opportunities to emit any sequence. Groups had unequal trial counts.
2. **Accuracy confound**: At non-zero contrast only, WT accuracy (89.2%) was higher
   than model mice (85-87%). At length 8, an all-correct sequence occurs at different
   rates purely from first-order statistics.
3. **Filtering artifact**: Using only non-zero contrast trials inflated the accuracy
   gap. Including 0% contrast (where all groups perform ~57%) shrinks it to <2%.

Criterion truncation (equal trials per animal) eliminates both confounds. The
dramatic 21-29 hit unidirectional result collapses to 3-6 bidirectional hits.

## Alternative Encoding: Toward/Against

To test whether genotypes differ in orientation relative to the current block's
favored side (rather than absolute left/right), we recoded:

| Symbol | Meaning |
|--------|---------|
| T✓ (0) | Toward favored side, correct |
| T✗ (1) | Toward favored side, incorrect |
| A✗ (2) | Against favored side, incorrect |
| A✓ (3) | Against favored side, correct |

"Favored" = the block's high-reward side (probLeft=0.8 → right is favored;
probLeft=0.2 → left is favored). In unbiased blocks, prior block's bias is
carried forward. Full trial stream preserved (no adjacency breaks).

### Results: Toward/Against

| Comparison | Sequences | Significant | WT > Model | Model > WT | k |
|---|---|---|---|---|---|
| WT vs Cntnap2 | 4,574 | 2 | 2 | 0 | 1 |
| WT vs Fmr1 | 4,517 | 1 | 0 | 1 | 1 |
| WT vs Shank3B | 4,548 | 3 | 0 | 3 | 1 |

Toward rates matched across genotypes: 69.5-70.0%. Even fewer hits than L/R
encoding. No evidence of genotype-specific block-relative sequential structure.

## Interpretation

The IBL 2AFC task does not produce detectable genotype-specific sequential
structure via CBAS. Reasons:

1. **Stimulus is randomized** — the animal doesn't control the trial-to-trial
   sequence of "what to respond to," unlike free-choice paradigms (T-maze, foraging).
2. **Strategy space is thin** — at 0% contrast, the only strategic variable is
   P(choose_left | block). That's a first-order statistic, not sequential structure.
3. **Sample size is marginal** — n=14-17 per group is below what the subsampling
   analysis (on flies) showed is needed for stable sequence detection.

The Noel et al. group effect (detected by their exponential-weighting model at
0% contrast) operates through biased block-state inference — a mechanism that
produces P(left|block) differences, not higher-order sequential dependencies.

## Reproducibility

```bash
# Controlled L/R analysis
pixi run noel-ibl-ctrl

# Toward/against encoding
pixi run noel-ibl-toward
```

Data: `data/noel_ibl_mice/all_contrasts/` (L/R) and `data/noel_ibl_mice/toward_against/`
Runtime: ~7s per comparison (M=10,000, criterion=4,000, L=6).
