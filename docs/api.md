# API Reference

## Pipeline functions

### `run_cbas_comparative`

```python
run_cbas_comparative(cohort, group_labels, params=None,
                     contingency=2, encode_reward=True, chunked=True,
                     block_aware=False)
```

Run the full comparative CBAS pipeline from raw data to significant sequences.

**Arguments**

- `cohort` (Cohort) - Subjects from `load_cohort`, or built directly.
- `group_labels` - how to group the cohort, in one of three forms. A
  `{subject_id: 0/1}` mapping, or the name of a `meta` column to derive one from, ties
  each label to a subject by name: the cohort's order then cannot affect the grouping,
  and neither can the count matrix's. A list or array is read in cohort order instead,
  which is the one form where the labels can be in the wrong order, and that is not
  detectable from the data. Prefer a mapping or a column name when the labels and the
  subject files were assembled separately. Raises
  `ValueError` if a subject has no label, if the count does not match the cohort, if any
  label is outside {0, 1}, or if either group is empty: each of those would otherwise
  analyse a cohort you did not pass, or produce all-NaN statistics, and return a
  complete-looking result either way. Exclude a subject with `Cohort.filter`, not by
  labelling it into neither group.
- `params` (CBASParams, optional) - Analysis parameters. Uses defaults if None.
- `contingency` (int or None) - Trial condition to filter on. None uses all trials.
- `encode_reward` (bool) - If True, symbol = choice + reward * num_arms (doubles alphabet). Set False for tasks where choice already encodes the outcome.
- `chunked` (bool) - If True, use the memory-efficient chunked pipeline.
- `block_aware` (bool) - If True, sequences cannot span block/session boundaries. Enable for multi-session experiments.

**Returns** `CBASResult`

---

### `run_cbas_correlative`

```python
run_cbas_correlative(cohort, covariate, params=None,
                     contingency=2, encode_reward=True, block_aware=False)
```

Run the full correlative CBAS pipeline.

**Arguments**

- `cohort` (Cohort) - Subjects from `load_cohort`, or built directly.
- `covariate` - one continuous value per subject, in the same three forms as
  `group_labels` above and with the same trade-off: a `{subject_id: value}` mapping or a
  `meta` column name is matched by name, while a list or array is read in cohort order.
  Raises `ValueError` on a length mismatch or a missing subject.
- `params` (CBASParams, optional) - Analysis parameters.
- `contingency` (int or None) - Trial condition to filter on. None uses all trials.
- `encode_reward` (bool) - If True, symbol = choice + reward * num_arms (doubles alphabet).
- `block_aware` (bool) - If True, sequences cannot span block/session boundaries.

**Returns** `CBASResult`

---

## Cohorts and subjects

Everything per-subject is keyed by subject id rather than by list position, so a
reordering or a filter cannot pair a subject with another's group, covariate or counts.

### `Subject`

```python
Subject(id, trials, meta=None, blocks=None)
```

One subject's trials, under the identity it was loaded with. `trials` is an
(n_trials, 4) array of session, choice, reward, condition, exposed also as the
`.session`, `.choice`, `.reward` and `.condition` properties. `condition` is the file's
fourth column for single-contingency data and the contingency block index for
multi-contingency data. Both are selected by value, so one filter path serves both, but
what a given value means depends on the format it came from; `has_blocks` distinguishes
them.

`meta` holds whatever the cohort table supplied, such as `lesion`, `sex` or `genotype`.
Group membership is deliberately not a field: the same subject is control in one
comparison and, say, male in another, so a label belongs to an analysis rather than to
the data.

`blocks` is present only for multi-contingency subjects; `has_blocks` reports it, and
the block methods raise without it rather than treating a raw condition column as a
validated block structure.

---

### `Cohort`

```python
Cohort(subjects)
```

An ordered collection of `Subject`, indexed by position or slice, with `len` and
iteration. Ids must be unique, since everything else is resolved through them, and they
are coerced to strings: cohort tables often number their animals, and a lookup that
worked for `"200"` but not `200` would be the worst of both.

Indexing is deliberately not overloaded to take an id as well. Ids are often numbers, so
`cohort[0]` and `cohort["0"]` would be two different questions in one syntax; `by_id`
asks the second one.

**Attributes and methods**

- `ids` - subject ids, in cohort order.
- `by_id(subject_id)` - the subject with that id, whether given as a string or a number.
- `filter(predicate=None, **meta_equals)` - the subjects passing a predicate and
  matching every `meta` value given, e.g. `cohort.filter(genotype="WT")`. A value may be
  a set or list to keep several.
- `reorder(order)` - a cohort in the given order of positions, identity following each
  subject.
- `labels_from(key, coder=None)` - 0/1 group labels from a `meta` column. `coder` maps a
  raw value to 0 or 1, or to None when it cannot be used; the default understands `0`/`1`
  and words like `control`, `sham`, `wt`, `lesion`, `ko`, `mutant`. Unusable subjects are
  named in a `ValueError` rather than dropped, since dropping them here would quietly
  change the cohort you assembled. Remove them with `filter` first.
- `covariate_from(key)` - a float covariate from a `meta` column.
- `meta_values(key)` - one `meta` value per subject, in cohort order.
- `has_blocks` - whether every subject carries contingency blocks.

---

### `CountMatrix`

```python
CountMatrix(counts, subject_ids, sequences)
```

Sequence counts with the ids of its rows and the labels of its columns, so the
correspondence can be checked rather than assumed. `counts[i, j]` is how often subject
`subject_ids[i]` used sequence `sequences[j]`.

**Attributes and methods**

- `shape` - `counts.shape`.
- `row(subject_id)` - one subject's counts, by id.
- `reorder(order)` - rows in the given order, ids following them.
- `select_block(block)` - the columns of one contingency block. Raises on a
  single-contingency matrix rather than selecting nothing.
- `column_labels(num_arms=6, encode_reward=True, join=" ")` - readable labels for every
  column, in the published convention.

---

### `load_cohort`

```python
load_cohort(source, pattern="*.txt", meta=None, ids=None)
```

A cohort from a directory or an explicit list of files. Files from a directory are
ordered by the digits in their names, so `an2` precedes `an10`; an explicit list keeps
the order given. Ids default to filename stems. `meta` may be a `{id: dict}` mapping or
a list in file order.

---

### `load_subject`

```python
load_subject(filepath, id=None, meta=None)
```

One subject from the single-contingency text format, with the filename stem as its
default id.

---

### `resolve_labels`

```python
resolve_labels(labels, ids, cohort_ids=None, what="labels")
```

Values for `ids`, from a `{id: value}` mapping or a list in cohort order. A list is
converted to a mapping first, so the result follows `ids` even when those are not in
cohort order. This is what the pipelines use to align labels and covariates to matrix
rows; `what` names the thing being resolved so the error messages read correctly.

---

### `default_group_coder`

```python
default_group_coder(value)
```

0, 1, or None for a raw group value from a cohort table. None means "not usable", which
is distinct from either group: a blank or unrecognised value is not evidence of
membership, so it is not guessed at.

Exact words are matched first, then a substring, so `Hippocampal Lesion` reads as 1. That
last step reads words rather than sentences, so a table spelling out `no lesion evident`
is also read as 1; pass an explicit `coder` to `labels_from` when a table phrases its
groups that way.

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
run_cbas_multicontingency(cohort, group_labels, params=None, blocks=None,
                          encode_reward=True, chunked=True)
```

Comparative CBAS across several contingencies, each counted as its own set of
hypotheses. `cohort` comes from `load_cohort_with_contingencies`. Returns a
`CBASResult` whose `sequences` entries are `(block, sequence_tuple)` pairs rather than
bare tuples. Everything downstream of the count matrix is the ordinary comparative
path, including the `group_labels` validation described under
`run_cbas_comparative`, counted against the cohort. See
[multiple contingencies](guide.md#multiple-contingencies).

---

### `load_cohort_with_contingencies`

```python
load_cohort_with_contingencies(directory, allow_mid_session_change=False)
```

Load every `an*.txt` subject file in a directory, ordered numerically, plus the
`anInfo.txt` table if present. Returns a `Cohort` whose subjects carry the info
table's columns as `meta`, so a grouping or filter can be named rather than assembled
alongside. Raises if the file count and the info row count disagree.

---

### `load_subject_with_contingencies`

```python
load_subject_with_contingencies(filepath, allow_mid_session_change=False, id=None,
                                meta=None)
```

One subject from the multi-contingency text format. Returns a `Subject` whose `condition`
column is the contingency block index and whose `blocks` describe each block's arms, with
the filename stem as the default id. Raises by default when a contingency changes partway
through a session; `allow_mid_session_change=True` permits it.

---

### `load_cohort_info`

```python
load_cohort_info(filepath)
```

Parse an `anInfo.txt` cohort table into a list of dicts, one per subject, in file order.
`load_cohort_with_contingencies` calls this and attaches each row to its subject as
`meta`, so you only need it directly to inspect the table itself.

---

### `shared_contingency_blocks`

```python
shared_contingency_blocks(cohort)
```

Sorted contingency block indices common to every subject, excluding exploration.
Raises if a block index denotes different arms for different subjects, or if some
subject did not run a contingency the others did.

---

### `build_multicontingency_count_matrix`

```python
build_multicontingency_count_matrix(cohort, params, blocks=None, encode_reward=True)
```

The count matrix behind `run_cbas_multicontingency`. Returns
a `CountMatrix` with columns grouped contiguously by contingency, so `select_block`
can take one. Row `i` is `cohort[i]`, and within a contingency
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
subject_criteria(cohort, params, contingency=2, block_aware=False)
```

Criterion trial index per subject, as `{subject_id: trial}` so a shortfall can be
reported by name. The values are floats so that `inf` survives; `inf` marks a subject
that never reached a higher-order criterion, which means it is not truncated and
contributes every window it has. Pass the same `contingency` and `block_aware` used
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

### `contingency_criteria`

```python
contingency_criteria(cohort, params, blocks=None)
```

Criterion trial index per subject per contingency, `inf` where unreached. Returns
`({subject_id: {block: trial}}, blocks_used)`. Separate from `subject_criteria` because
it answers a different question: the criterion applies within each contingency, so a
subject can reach it in one and fall short in another.

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
build_count_matrix(cohort, params, contingency=2, encode_reward=True,
                   block_aware=False)
```

Build the full sequence count matrix across all subjects and all sequence lengths 1 through `params.seq_len_max`.

**Arguments**

- `cohort` (Cohort) - Subjects from `load_cohort`, or built directly. A bare list of
  arrays raises `TypeError`, since it carries no subject identity.
- `params` (CBASParams) - Analysis parameters.
- `contingency` (int or None) - Trial condition to filter on. None uses all trials.
- `encode_reward` (bool) - If True, symbol = choice + reward * num_arms.
- `block_aware` (bool) - If True, sequences cannot span block/session boundaries.

**Returns** a `CountMatrix` with `.counts` of shape (n_subjects, n_sequences),
`.subject_ids` naming the rows, and `.sequences` labelling the columns.

#### Row and column order

**Row `i` is `cohort[i]`, named by `subject_ids[i]`.** The order you pass in is the order
you get back; subjects are never sorted or grouped, and this function is not given the
group labels, so group structure cannot affect it.

The ids are why nothing per-subject has to be aligned positionally. Group labels and
covariates are resolved through them, so reordering a cohort or a matrix relabels
nothing, and `CountMatrix.reorder` moves the rows and their ids together. Anything else
you align to the rows should be keyed by id too.

The caller still decides the order, and a loader may impose one of its own:
`load_cohort` and `load_cohort_with_contingencies` order subjects by the digits in their
filenames, and a caller that sorts by group before building gets group-blocked rows.

**Column `j` is `sequences[j]`.** Columns are ordered by total count summed over every
subject, descending, with ties broken by sequence length and then by sequence value.
They are *not* in order of first appearance, and the order depends on the cohort:
because it is driven by cohort-wide totals, adding or removing one subject can move
most columns. Only sequences observed in at least one subject get a column at all.

`sequences` travels with the matrix, so index through it. A bare column position is not
meaningful across two runs, even two runs on nearly the same cohort.
`CountMatrix.column_labels` renders them all readably, and `decode_sequence` does one.

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
romano_wolf_stepdown(test_stats, null_matrix, null_directions=None, k=1,
                     tie_rtol=0.0)
```

Apply the Romano-Wolf step-down procedure at a fixed k.

**Returns** ndarray of shape (2S,) adjusted p-values. NaN where test_stats is NaN.

---

### `find_k_fwer`

```python
find_k_fwer(test_stats, null_matrix, alpha=0.5, gamma=0.05, null_directions=None,
            return_history=False, tie_rtol=0.0)
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
                   resample_number=10000, encode_reward=True, n_hypothesis_sets=1)
```

Estimate memory and time requirements before running an analysis.

**Returns** dict with keys: `alphabet`, `seq_len_max`, `total_sequences`, `observed_sequences`, `resample_number`, `n_subjects`, `memory_full_null_gb`, `memory_chunked_gb`, `est_time_seconds`, `recommendation`.

---

### `print_resource_estimate`

```python
print_resource_estimate(est)
```

Pretty-print the output of `estimate_resources`.
