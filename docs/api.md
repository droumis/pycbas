# API Reference

## Pipeline functions

### `run_cbas_comparative`

```python
run_cbas_comparative(subjects_data, group_labels, params=None,
                     contingency=2, encode_reward=True, chunked=True,
                     block_aware=False)
```

Run the full comparative CBAS pipeline from raw data to significant sequences.

**Arguments**

- `subjects_data` (list of ndarray) - One array per subject from `load_subject_data`.
- `group_labels` (array-like of int) - 0 or 1 per subject indicating group membership.
- `params` (CBASParams, optional) - Analysis parameters. Uses defaults if None.
- `contingency` (int or None) - Trial condition to filter on. None uses all trials.
- `encode_reward` (bool) - If True, symbol = choice + reward * num_arms (doubles alphabet). Set False for tasks where choice already encodes the outcome.
- `chunked` (bool) - If True, use the memory-efficient chunked pipeline.
- `block_aware` (bool) - If True, sequences cannot span block/session boundaries. Enable for multi-session experiments.

**Returns** `CBASResult`

---

### `run_cbas_correlative`

```python
run_cbas_correlative(subjects_data, covariate, params=None,
                     contingency=2, encode_reward=True, block_aware=False)
```

Run the full correlative CBAS pipeline.

**Arguments**

- `subjects_data` (list of ndarray) - One array per subject from `load_subject_data`.
- `covariate` (array-like of float) - One continuous value per subject (e.g. a behavioral score).
- `params` (CBASParams, optional) - Analysis parameters.
- `contingency` (int or None) - Trial condition to filter on. None uses all trials.
- `encode_reward` (bool) - If True, symbol = choice + reward * num_arms (doubles alphabet).
- `block_aware` (bool) - If True, sequences cannot span block/session boundaries.

**Returns** `CBASResult`

---

## Data classes

### `CBASParams`

```python
CBASParams(num_arms=6, seq_len_max=6, criterion=800, criterion_order=0,
           resample_number=10000, alpha=0.5, gamma=0.05, centering=False)
```

All analysis parameters. See the [User Guide](guide.md#parameters) for descriptions.

---

### `run_cbas_multicontingency`

```python
run_cbas_multicontingency(records, group_labels, params=None, blocks=None,
                          encode_reward=True, chunked=True)
```

Comparative CBAS across several contingencies, each counted as its own set of
hypotheses. `records` comes from `load_cohort_with_contingencies`. Returns a
`CBASResult` whose `sequences` entries are `(block, sequence_tuple)` pairs rather than
bare tuples. Everything downstream of the count matrix is the ordinary comparative
path. See [multiple contingencies](guide.md#multiple-contingencies).

---

### `load_cohort_with_contingencies`

```python
load_cohort_with_contingencies(directory, allow_mid_session_change=False)
```

Load every `an*.txt` subject file in a directory, ordered numerically, plus the
`anInfo.txt` table if present. Returns `(records, info)`. Raises if the file count and
the info row count disagree.

---

### `shared_contingency_blocks`

```python
shared_contingency_blocks(records)
```

Sorted contingency block indices common to every subject, excluding exploration.
Raises if a block index denotes different arms for different subjects, or if some
subject did not run a contingency the others did.

---

### `build_multicontingency_count_matrix`

```python
build_multicontingency_count_matrix(records, params, blocks=None, encode_reward=True)
```

The count matrix behind `run_cbas_multicontingency`. Returns
`(sequences, count_matrix)` with columns grouped contiguously by contingency, so a
single contingency can be sliced out.

---

### `criterion_trial`

```python
criterion_trial(reward_blocks, order, count)
```

Trial index at which `count` units of the given order have accumulated, or `inf` if the
subject never gets there. `reward_blocks` is a list of per-session 0/1 arrays.
`reached_criterion(value)` tests the result.

---

### `subject_criteria`

```python
subject_criteria(subjects_data, params, contingency=2, block_aware=False)
```

Per-subject criterion trial index, as a float array so that `inf` survives. `inf` marks
a subject that never reached a higher-order criterion, which means it is not truncated
and contributes every window it has. Pass the same `contingency` and `block_aware` used
for the run, since the criterion is expressed in the same coordinates as the
enumeration. Only interesting when `params.criterion_order` is nonzero; see
[higher-order criteria](guide.md#higher-order-criteria).

---

### `CBASResult`

Returned by the pipeline functions.

**Attributes**

- `sequences` (list of tuple) - All evaluated sequences, sorted by total frequency descending.
- `test_stats` (ndarray, shape 2S) - Observed test statistics. Paired layout: indices `[i*2]` and `[i*2+1]` are positive and negative directions for sequence i.
- `g_values` (ndarray, shape 2S) - Adjusted p-values from the converged step-down. Same layout as test_stats.
- `k_final` (int) - Converged k value from k-FWER iteration.
- `significant_mask` (ndarray of bool, shape S) - True for sequences significant in either direction.

**Properties**

- `n_significant` (int) - Count of significant sequences.

---

## I/O functions

### `load_subject_data`

```python
load_subject_data(filepath)
```

Load a single subject's data file. Expects comma-separated rows with four columns: session, choice, reward, contingency.

**Returns** ndarray of shape (n_trials, 4), dtype int32.

---

### `extract_choice_stream`

```python
extract_choice_stream(subject_data, contingency=2, num_arms=6, encode_reward=True)
```

Extract the choice stream from a subject's data array, filtering by contingency and optionally encoding reward into the symbol.

**Returns** 1D ndarray of integer symbols.

---

### `extract_choice_streams_by_block`

```python
extract_choice_streams_by_block(subject_data, contingency=2, num_arms=6, encode_reward=True)
```

Extract choice streams split by block/session boundaries. Used internally when `block_aware=True`.

**Returns** list of 1D ndarrays, one per block.

---

### `enumerate_sequences`

```python
enumerate_sequences(choice_stream, seq_len, criterion)
```

Count all subsequences of a given length with start position <= criterion.

**Returns** dict mapping sequence tuple to count.

---

### `enumerate_sequences_block_aware`

```python
enumerate_sequences_block_aware(block_streams, seq_len, criterion)
```

Count sequences within blocks, never crossing block boundaries. The criterion is a
maximum global start position in the concatenated stream, inclusive; positions that
cannot start a sequence within their own block still advance that global position. So
it is not the same as a cap on how many sequences are counted.

**Returns** dict mapping sequence tuple to count.

---

## Core computation

### `build_count_matrix`

```python
build_count_matrix(subjects_data, params, contingency=2, encode_reward=True,
                   block_aware=False)
```

Build the full sequence count matrix across all subjects and all sequence lengths 1 through `params.seq_len_max`.

**Arguments**

- `subjects_data` (list of ndarray) - One array per subject.
- `params` (CBASParams) - Analysis parameters.
- `contingency` (int or None) - Trial condition to filter on. None uses all trials.
- `encode_reward` (bool) - If True, symbol = choice + reward * num_arms.
- `block_aware` (bool) - If True, sequences cannot span block/session boundaries.

**Returns** `(sequences, count_matrix)` where sequences is a list of tuples and count_matrix is ndarray of shape (n_subjects, n_sequences).

---

### `compute_test_stats`

```python
compute_test_stats(count_matrix, group_indices)
```

Compute studentized two-sample test statistics for all sequences using two one-tailed tests.

**Arguments**

- `count_matrix` (ndarray, shape N x S)
- `group_indices` (list of two arrays) - Indices into count_matrix rows for each group.

**Returns** ndarray of shape (2S,). NaN where the test stat is not in that direction.

---

### `compute_test_stats_correlative`

```python
compute_test_stats_correlative(count_matrix, covariate)
```

Compute studentized correlation test statistics using the robust tau estimator.

**Returns** ndarray of shape (2S,).

---

## Bootstrap

### `bootstrap_test_stats`

```python
bootstrap_test_stats(count_matrix, group_indices, params, rng=None)
```

Generate the bootstrap null by resampling subjects from the pooled sample ignoring group labels.

**Returns** `(null_matrix, null_directions)` where null_matrix is (M, S) float64 magnitudes and null_directions is (M, S) int8 (0=positive, 1=negative, -1=undefined).

---

### `bootstrap_test_stats_correlative`

```python
bootstrap_test_stats_correlative(count_matrix, covariate, params, rng=None)
```

Generate the bootstrap null for correlative mode by permuting the covariate.

**Returns** `(null_matrix, null_directions)` same format as above.

---

## Step-down and k-FWER

### `romano_wolf_stepdown`

```python
romano_wolf_stepdown(test_stats, null_matrix, null_directions=None, k=1)
```

Apply the Romano-Wolf step-down procedure at a fixed k.

**Returns** ndarray of shape (2S,) adjusted p-values. NaN where test_stats is NaN.

---

### `find_k_fwer`

```python
find_k_fwer(test_stats, null_matrix, alpha=0.5, gamma=0.05, null_directions=None)
```

Run iterative k-FWER to convergence and return the final adjusted p-values.

**Returns** `(g_values, k_final)`.

---

### `find_k_fwer_chunked`

```python
find_k_fwer_chunked(test_stats, count_matrix, group_indices, params,
                    chunk_size=500, rng=None, return_history=False, tie_rtol=0.0)
```

Memory-efficient variant that generates bootstrap directly into the null submatrix in
row-chunks. Produces results identical to `find_k_fwer(..., null_directions=...)`; it
always applies direction-conditional removal, so it does not match a call to
`find_k_fwer` that omits the directions.

**Returns** `(g_values, k_final)`.

---

### `tie_rtol` on the step-down functions

`romano_wolf_stepdown`, `find_k_fwer`, `find_k_fwer_k1` and `find_k_fwer_chunked` all
accept `tie_rtol`, a relative slack on the `null >= observed` comparison. It defaults
to `0.0`, which is correct for an integer count matrix: there the observed and
bootstrap statistics agree bitwise, so the comparison is exact and tolerating anything
would be arbitrary.

A non-integer matrix, such as one normalised to rates, has no such guarantee, and a
strict comparison then discards the resamples that represent the observed value, which
adds false positives. `compute_test_stats` warns once when it sees such a matrix.
Passing the integer count matrix is the better fix where the normalising denominator is
common to all subjects, since the statistic is scale-invariant. Otherwise a derived
bound is available as `pycbas._moments.tie_rtol_for(matrix)`; it is private and
experimental rather than part of the supported surface.

---

### `find_k_fwer_k1`

```python
find_k_fwer_k1(test_stats, null_matrix, alpha=0.5, gamma=0.05,
               null_directions=None, tie_rtol=0.0)
```

Conservative variant that always uses k=1 (standard FWER). Useful for comparison and debugging.

**Returns** `(g_values, k_final)` where k_final is what k *would* be from the iteration formula.

---

## Resource estimation

### `estimate_resources`

```python
estimate_resources(num_arms, seq_len_max, n_subjects=None, n_observed=None,
                   resample_number=10000, encode_reward=True)
```

Estimate memory and time requirements before running an analysis.

**Returns** dict with keys: `alphabet`, `seq_len_max`, `total_sequences`, `observed_sequences`, `resample_number`, `n_subjects`, `memory_full_null_gb`, `memory_chunked_gb`, `est_time_seconds`, `recommendation`.

---

### `print_resource_estimate`

```python
print_resource_estimate(est)
```

Pretty-print the output of `estimate_resources`.
