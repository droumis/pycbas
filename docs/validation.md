# Validation

pycbas is validated against the original Igor Pro implementation by David Kastner.

## Results

| Dataset | Mode | Subjects | Sequences | pycbas | Igor | k |
|---|---|---|---|---|---|---|
| Flies | Comparative | 1,566 | 2,046 | 1,605 (78.4%) | 1,605 (78.4%) | 81 |
| Humans | Correlative | 1,413 | 408 | 31 (7.6%) | 31 (7.6%) | 2 |
| Rats | Comparative | 85 | 16,500 | 178 (1.1%) | 386* | 9 |

*Rat comparison uses a subset of the full dataset (85 of the subjects used in the paper). Not directly comparable.

## Fly validation

Binary left/right choice task (spontaneous alternation), comparing Canton-S (CA) and w1118 strains. Parameters: 2 arms, L=10, criterion=250, M=10,000.

Significance counts match exactly (1,605/2,046). Test statistics match to a maximum absolute difference of 1.2e-06. Rank ordering of test statistics matches for 2044 of 2046 sequences. The 2 mismatches are at adjacent ranks (899 and 900) where both sequences have t values of approximately 4.168849, differing by 1.5e-07 (floating point tie-breaking).

## Human validation

Two-step decision task (correlative mode with CBIT compulsivity score). Parameters: 6 arms, L=4, criterion=400, M=10,000.

Exact match: 31/408 significant sequences with k=2.

## Key implementation details that affect correctness

Three behaviors in the Igor code were non-obvious and required careful matching:

**1. Magnitude-based null with direction tracking.** The bootstrap stores |t| per sequence per resample (not direction-specific values). Direction is tracked separately for the removal step. Getting this wrong produced 69 significant sequences for humans instead of 31.

**2. Direction-conditional removal in step-down.** When a sequence is removed from the active set during step-down, it is only removed from bootstrap rows where the bootstrap went the same direction as the observed statistic. If the bootstrap went the other direction, that row keeps the magnitude active. Getting this wrong produced 2,046/2,046 significant for flies (everything significant).

**3. Criterion boundary (inclusive start position).** Igor checks `start_position <= criterion` (0-based inclusive), giving `criterion + 1` valid start positions. A boundary error here (using strict less-than, giving `criterion` start positions) caused an 11-sequence discrepancy in flies (1,594 instead of 1,605).

