# pycbas Validation Summary

Cross-species validation of the CBAS reimplementation against Kastner et al. (2026).

## Results at a glance

| Dataset | Mode | Subjects | Sequences | pycbas | David (Igor) | k |
|---|---|---|---|---|---|---|
| Flies | Comparative | 1,566 | 2,046 | 1,605 (78.4%) | 1,605 (78.4%) | 81 |
| Humans | Correlative | 1,413 | 408 | 31 (7.6%) | 31 (7.6%) | 2 |
| Rats | Comparative | 105 | 16,378 | 572 (3.5%) | 572 (3.5%) | 29 |

## Flies

- **2 arms, seq_len_max=10, criterion=250, M=10,000**
- 1,605/2,046 significant (k=81) — **exact match** with David
- Test statistics match to 1.2e-06; rank ordering matches 2044/2046 (2 ties at floating point precision)
- Runtime: ~21s | Peak RAM: ~560 MB

[Full report](flies/validation_report.md)

## Humans

- **6 arms, seq_len_max=4, criterion=400, M=10,000**
- 31/408 significant (k=2) — **exact match** with David
- Runtime: ~3s | Peak RAM: ~155 MB

[Full report](humans/validation_report.md)

## Rats

- **6 arms, seq_len_max=6, criterion=800, M=10,000, block_aware=True**
- 572/16,378 significant (k=29) — **exact match** with David
- 105 subjects (55 control, 50 lesion), all_published cohort (genotype==0, lesion known)
- Test statistics match David's Igor within 1e-6 on all 16,376 overlapping sequences
- All 572 significant sequences identical between pycbas and Igor
- Runtime: ~10s | Peak RAM: ~3.6 GB

[Full report](rats/validation_report.md)

## Algorithm notes

Three key fixes brought us into exact alignment with David's Igor implementation:

1. **Magnitude-based null** — Store |t| per sequence per bootstrap row (not direction-specific). This fixed humans from 69→31 (exact match).

2. **Direction-conditional removal** — In the step-down, only remove a sequence from a bootstrap row when the bootstrap went the same direction as the observed stat. This fixed flies from 2,046/2,046 (100%) → 1,594/2,046 (78%).

3. **Criterion boundary (inclusive)** — Igor checks `start_position <= criterion` (inclusive, 0-based), giving 251 counting windows per subject. We were using `stream[:criterion]` which gave only 250 elements (250−L+1 windows). This missed L counting windows per subject per sequence length and caused the 11-sequence discrepancy: 1,594 → 1,605 (exact match).

See [notes/algorithm_comparison.md](../notes/algorithm_comparison.md) for detailed pseudocode comparison.
