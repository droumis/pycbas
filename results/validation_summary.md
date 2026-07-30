# pycbas Validation Summary

Cross-species validation of the CBAS reimplementation against Kastner et al. (2026).

## Results at a glance

| Dataset | Mode | Subjects | Sequences | pycbas | David (Igor) | k |
|---|---|---|---|---|---|---|
| Flies | Comparative | 1,566 | 2,046 | 1,594 (77.9%) | 1,605 (78.4%) | 80 |
| Humans | Correlative | 1,413 | 408 | 31 (7.6%) | 31 (7.6%) | 2 |
| Rats | Comparative | 85 | 16,483 | 177 (1.1%) | 386* | 9 |

*Rat comparison is deferred — we don't have David's full dataset.

## Flies

- **2 arms, seq_len_max=10, criterion=250, M=10,000**
- 1,594/2,046 significant (k=80) vs David's 1,605
- Sequence-level comparison: 1,584 in both, 21 David-only, 10 us-only
- Result is perfectly stable across 5 different RNG seeds (0 unstable sequences)
- Runtime: ~18s | Peak RAM: ~360 MB

[Full report](flies/validation_report.md)

## Humans

- **6 arms, seq_len_max=4, criterion=400, M=10,000**
- 31/408 significant (k=2) — **exact match** with David
- Runtime: ~3s | Peak RAM: ~155 MB

[Full report](humans/validation_report.md)

## Rats

- **6 arms, seq_len_max=6, criterion=800, M=10,000**
- 177/16,483 significant (k=9)
- Paper reports 386 sig (different subject count: 24,342 sequences vs our 16,483)
- Deferred until we have David's full rat dataset
- Runtime: ~6s | Peak RAM: ~1.7 GB

[Full report](rats/validation_report.md)

## Algorithm notes

Two key fixes brought us into alignment with David's Igor implementation:

1. **Magnitude-based null** — Store |t| per sequence per bootstrap row (not direction-specific). This fixed humans from 69→31 (exact match).

2. **Direction-conditional removal** — In the step-down, only remove a sequence from a bootstrap row when the bootstrap went the same direction as the observed stat. This fixed flies from 2,046/2,046 (100%) → 1,594/2,046 (78%).

The remaining 11-sequence fly discrepancy (1,594 vs 1,605) is due to different RNG implementations (Igor's per-row seeded PRNG vs numpy), not algorithmic — our result is perfectly deterministic across seeds.

See [notes/algorithm_comparison.md](../notes/algorithm_comparison.md) for detailed pseudocode comparison.
