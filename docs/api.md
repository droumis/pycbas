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
- `group_labels` (array-like of int) - 0 or 1 per subject indicating group membership,
  in the same order as `subjects_data`. Raises `ValueError` if the count does not match
  the number of subjects, if any label is outside {0, 1}, or if either group is empty:
  each of those would otherwise analyse a cohort you did not pass, or produce all-NaN
  statistics, and return a complete-looking result either way. Exclude a subject by
  dropping it from `subjects_data`, not by labelling it into neither group.
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
- `covariate` (array-like of float) - One continuous value per subject (e.g. a behavioral
  score), in the same order as `subjects_data`. Raises `ValueError` on a length mismatch.
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
path, including the `group_labels` validation described under
`run_cbas_comparative`, counted against `records`. See
[multiple contingencies](guide.md#multiple-contingencies).

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
single contingency can be sliced out. Row `i` is `records[i]`, and within a contingency
the columns follow the same rule as `build_count_matrix`; see
[Row and column order](#row-and-column-order). Each `sequences` entry is a
`(block, sequence)` pair rather than a bare tuple.

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

### `decode_symbol`

```python
decode_symbol(sym, num_arms=6, encode_reward=True)
```

Invert `extract_choice_stream`'s encoding.

**Returns** `(choice, rewarded)`, choice 0-based and rewarded a bool, or None when
`encode_reward` is False, since the outcome is then not recoverable from the symbol.
Raises `ValueError` if the symbol is out of range for `num_arms`, because decoding with
the wrong `num_arms` otherwise reports a plausible wrong arm.

---

### `decode_sequence`

```python
decode_sequence(entry, num_arms=6, encode_reward=True, join=" ")
```

Readable label for one entry of a `sequences` list, in the published convention: arms
numbered from 1, a trailing `*` for a rewarded choice. With `num_arms=6` the symbol 2
reads as `3` and the symbol 8 as `3*`. A `(block, sequence)` entry from the
multi-contingency pipeline is prefixed with its block.

**Returns** str.

---

### `split_sequence_entry`

```python
split_sequence_entry(entry)
```

`(block, symbols)` for one entry of a `sequences` list, with block None for a bare
sequence. Use this rather than `len(entry)` to get a sequence's length: a
multi-contingency entry is `(block, (symbols...))`, so `len` is 2 for every hypothesis
whatever its actual length.

---

### `enumerate_sequences`

```python
enumerate_sequences(choice_stream, seq_len, criterion)
```

Count all subsequences of a given length with start position <= criterion.

**Returns** dict mapping sequence tuple to count.

---

### `record_criteria`

```python
record_criteria(records, params, blocks=None)
```

Criterion trial index per subject per contingency, `inf` where a subject never reached
it. The multi-contingency counterpart of `subject_criteria`, which takes per-subject
arrays and so cannot be used on `SubjectRecord`s.

The criterion applies within each contingency, so a subject can reach it in one and fall
short in another, and the result is a matrix rather than a vector. Worth checking before
a run: a subject that falls short is not truncated and contributes every window it has,
so the weakest subjects contribute the most data, and counting more contingencies means
more chances to fall short.

**Returns** `(criteria, blocks)` where `criteria` has shape `(n_subjects, n_blocks)`.

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

#### Row and column order

**Row `i` is `subjects_data[i]`.** The order you pass in is the order you get back;
subjects are never sorted or grouped, and this function is not given the group labels,
so group structure cannot affect it. Group membership is applied later as indices into
these rows, so the two groups need not be contiguous. Align group labels, covariates
and any per-subject metadata to the matrix positionally.

Note that the caller decides that order, and a cohort loader may impose one of its own:
`load_cohort_with_contingencies` orders subjects by the digits in their filenames, and a
caller that sorts its subjects by group before building the matrix gets group-blocked
rows. Nothing in the returned value records which order was used.

**Column `j` is `sequences[j]`.** Columns are ordered by total count summed over every
subject, descending, with ties broken by sequence length and then by sequence value.
They are *not* in order of first appearance, and the order depends on the cohort:
because it is driven by cohort-wide totals, adding or removing one subject can move
most columns. Only sequences observed in at least one subject get a column at all.

So keep `sequences` with the matrix and index through it. A bare column position is not
meaningful across two runs, even two runs on nearly the same cohort. `decode_sequence`
turns an entry into a readable label.

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
