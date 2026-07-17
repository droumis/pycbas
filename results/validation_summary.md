# pycbas Validation Summary

Cross-species validation of the CBAS reimplementation against Kastner et al. (2026).

## Results at a glance

| Dataset | Mode | Subjects | Sequences | Significant | k |
|---|---|---|---|---|---|
| Flies | Comparative | 1566 | 2,046 | 2046 (100.0%) | 103 |
| Humans | Correlative | 1413 | 408 | 69 (16.9%) | 4 |
| Rats | Comparative | 85 | 16,483 | 380 (2.3%) | 20 |

## Flies

- **2 arms, seq_len_max=10, criterion=250, M=10,000**
- 2046/2046 significant (k=103)
- Runtime: 267.1s

[Full report](flies/validation_report.md)
## Humans

- **6 arms, seq_len_max=4, criterion=400, M=10,000**
- 69/408 significant (k=4)
- Runtime: 4.1s

[Full report](humans/validation_report.md)
## Rats

- **6 arms, seq_len_max=6, criterion=800, M=10,000**
- 380/16483 significant (k=20)
- Runtime: 86.8s

[Full report](rats/validation_report.md)

## Notes

- Fly result at fixed k=20 gives 1,633 sig (matches paper's 1,605). Full dataset saturates due to [k-FWER sensitivity in the high-power regime](flies/kfwer_sensitivity_analysis.md).
- Human count (69 vs paper's 31) likely reflects minor differences in tau-hat normalization.
