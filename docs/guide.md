# User Guide

## What CBAS does

CBAS tests whether behavioral sequences (e.g. LLRRL, choice patterns in a maze) occur at different rates between experimental groups or correlate with a continuous measure. It evaluates all subsequences up to a specified length, controls for multiple comparisons using Romano-Wolf step-down, and controls false discovery proportion via iterative k-FWER.

The output is a set of significant sequences with adjusted p-values (called g-values or zeta-values in the paper).

## Modes

**Comparative** tests whether sequences differ between two groups (e.g. control vs lesion). The test statistic is a studentized difference in usage rates.

**Correlative** tests whether sequences correlate with a continuous covariate (e.g. a compulsivity score). The test statistic is a studentized Pearson correlation.

## Data format

Each subject's data is a comma-separated text file with four columns per row:

```
session, choice, reward, contingency
0,3,1,2
0,2,0,2
0,5,1,2
1,1,0,2
```

- `session` is an integer session ID (unused by the algorithm but required in the format)
- `choice` is the arm/symbol chosen (0-indexed, must be < num_arms)
- `reward` is 0 or 1
- `contingency` is the block type (used to filter trials)

For binary tasks like the fly spontaneous alternation, there are 2 arms and no reward encoding. For tasks like the rat 8-arm maze, there are 6 arms with reward encoding (doubling the effective alphabet to 12 symbols).

## Running a comparative analysis

```python
from pycbas import CBASParams, load_subject_data, run_cbas_comparative

# Load data
files_group0 = [...]  # paths to control subject files
files_group1 = [...]  # paths to experimental subject files
subjects_data = [load_subject_data(f) for f in files_group0 + files_group1]
group_labels = [0] * len(files_group0) + [1] * len(files_group1)

# Configure
params = CBASParams(
    num_arms=6,         # number of choice symbols
    seq_len_max=6,      # test all lengths 1..6
    criterion=800,      # use first 800 trials per subject
    resample_number=10000,
)

# Run
result = run_cbas_comparative(
    subjects_data, group_labels, params,
    contingency=2,       # filter to block type 2
    encode_reward=True,  # symbol = choice + reward * num_arms
)

print(f"{result.n_significant} significant sequences (k={result.k_final})")
```

## Running a correlative analysis

```python
from pycbas import CBASParams, load_subject_data, run_cbas_correlative
import numpy as np

subjects_data = [load_subject_data(f) for f in subject_files]
scores = np.array([...])  # one continuous value per subject

params = CBASParams(
    num_arms=6,
    seq_len_max=4,
    criterion=400,
    resample_number=10000,
)

result = run_cbas_correlative(subjects_data, scores, params)
```

## Working with results

`run_cbas_comparative` and `run_cbas_correlative` return a `CBASResult`:

```python
result.sequences          # list of tuples, e.g. [(0,1), (1,0,1), ...]
result.n_significant      # count of significant sequences
result.k_final            # converged k value
result.significant_mask   # boolean array, one per sequence
result.test_stats         # array of shape (2 * n_sequences,)
result.g_values           # adjusted p-values, same shape as test_stats
```

The test_stats and g_values arrays use a paired layout. For sequence index `i`:
- `test_stats[i*2]` is the positive-direction statistic (group 0 > group 1)
- `test_stats[i*2 + 1]` is the negative-direction statistic (group 1 > group 0)
- Exactly one of the pair is NaN (the direction that wasn't observed)

To find which direction a significant sequence went:

```python
for i, seq in enumerate(result.sequences):
    if not result.significant_mask[i]:
        continue
    pos_g = result.g_values[i * 2]
    neg_g = result.g_values[i * 2 + 1]
    if not np.isnan(pos_g) and pos_g < 0.5:
        print(f"{seq} group0 > group1, g={pos_g:.4f}")
    elif not np.isnan(neg_g) and neg_g < 0.5:
        print(f"{seq} group1 > group0, g={neg_g:.4f}")
```

## Parameters

| Parameter | Default | Effect |
|-----------|---------|--------|
| `num_arms` | 6 | Base alphabet size. With `encode_reward=True`, effective alphabet is `num_arms * 2`. |
| `seq_len_max` | 6 | Maximum sequence length L. All lengths 1 through L are tested. |
| `criterion` | 800 | How many trials per subject to use for counting. Sequences starting at positions 0 through `criterion` (inclusive) are counted. |
| `resample_number` | 10,000 | Number of bootstrap resamples M. More gives tighter p-values but costs linearly in time and memory. |
| `alpha` | 0.5 | Significance threshold for FDP control. The paper uses 0.5 (median FDP). |
| `gamma` | 0.05 | FDP tolerance. Fraction of rejections allowed to be false. |
| `centering` | False | Whether to center the bootstrap null by subtracting the observed delta. False matches the Igor implementation. True is slightly more liberal. |

## Choosing parameters

**num_arms** is determined by your task. Binary choice = 2. Six-arm maze = 6.

**encode_reward** should be True when reward is informative and not fully determined by the choice. Set False for deterministic tasks where the choice symbol already encodes the outcome (e.g. a two-alternative forced choice where left always means stimulus A).

**seq_len_max** controls how many sequences exist in the hypothesis space. The total is `sum(A^l for l in 1..L)` where A is the effective alphabet. For 12 symbols and L=6, that's about 2 million sequences. In practice, only observed sequences are tested, which is much smaller.

**criterion** should match the minimum usable trial count across your subjects. If some subjects have only 300 trials, set criterion to 300 or less.

**resample_number** of 10,000 is standard. For exploratory work, 1,000 is faster with coarser p-values.

## Memory and chunked mode

The main memory cost is the bootstrap null matrix (M rows by number of valid test statistics). For large hypothesis spaces (e.g. rats with 16,500 sequences), this can reach several GB.

The `chunked=True` option (default in `run_cbas_comparative`) generates the bootstrap in row-chunks directly into the sorted submatrix, avoiding the full intermediate allocation. This uses about 40% less peak memory at the cost of about 30% more time.

```python
# Lower memory (default)
result = run_cbas_comparative(subjects_data, labels, params, chunked=True)

# Faster, more memory
result = run_cbas_comparative(subjects_data, labels, params, chunked=False)
```

Use `estimate_resources` to check before running:

```python
from pycbas import estimate_resources, print_resource_estimate

est = estimate_resources(num_arms=12, seq_len_max=8, n_observed=5000)
print_resource_estimate(est)
```

## Using individual pipeline stages

For custom workflows, you can call each stage separately:

```python
from pycbas import (
    CBASParams, load_subject_data, build_count_matrix,
    compute_test_stats, bootstrap_test_stats, find_k_fwer,
)

params = CBASParams(num_arms=2, seq_len_max=10, criterion=250, resample_number=10000)

# Build count matrix
sequences, count_matrix = build_count_matrix(subjects_data, params, contingency=1)

# Compute test statistics
group_indices = [np.where(labels == 0)[0], np.where(labels == 1)[0]]
test_stats = compute_test_stats(count_matrix, group_indices)

# Generate bootstrap null
null_matrix, null_directions = bootstrap_test_stats(count_matrix, group_indices, params)

# Run step-down with k-FWER
g_values, k_final = find_k_fwer(
    test_stats, null_matrix, alpha=0.5, gamma=0.05,
    null_directions=null_directions,
)
```

This is useful when you want to inspect intermediate results, use a custom bootstrap, or run the step-down with different alpha/gamma without recomputing the null.
