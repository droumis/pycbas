# pycbas Validation Summary

Cross-species validation of the CBAS reimplementation against Kastner et al. (2026).

## Results at a glance

| Dataset | Mode | Subjects | Sequences | Significant | k |
|---|---|---|---|---|---|
| Flies | Comparative | 1566 | 2,046 | 1243 (60.8%) | 63 |
| Humans | Correlative | 1413 | 408 | 31 (7.6%) | 2 |
| Rats | Comparative | 85 | 16,483 | 111 (0.7%) | 6 |

## Flies

- **2 arms, seq_len_max=10, criterion=250, M=10,000**
- 1243/2046 significant (k=63)
- Runtime: 27.5s

[Full report](flies/validation_report.md)
## Humans

- **6 arms, seq_len_max=4, criterion=400, M=10,000**
- 31/408 significant (k=2)
- Runtime: 5.6s

[Full report](humans/validation_report.md)
## Rats

- **6 arms, seq_len_max=6, criterion=800, M=10,000**
- 111/16483 significant (k=6)
- Runtime: 24.9s

[Full report](rats/validation_report.md)

## Notes

- **Bootstrap null:** Default is uncentered (no centering), matching David's Igor implementation. Centering per Clarke et al. 2020 is available via `CBASParams(centering=True)` but produces a more liberal (less conservative) result.
- **Fly k-FWER:** Our iteration converges at k=63 giving 1243 significant — a strict subset of David's 1,605 (0 overcalled, 362 missed). The convergence path differs (our k jumps 1→63; David's lands at 1,605 via what we estimate is a more gradual path, though other implementation differences may also contribute).
- **Human:** Perfect match (31/408 = paper's 31).
- **Rat:** Paper reports 409/24,342 sig. Our different sequence count (16483 vs 24,342) reflects subject subset differences at longer lengths.
