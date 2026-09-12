# User Guide

## What CBAS does

CBAS tests whether behavioral sequences (e.g. LLRRL, choice patterns in a maze) occur at different rates between experimental groups or correlate with a continuous measure. It evaluates all subsequences up to a specified length, controls for multiple comparisons using Romano-Wolf step-down, and controls false discovery proportion via iterative k-FWER.

The output is a set of significant sequences with adjusted p-values (called g-values or zeta-values in the paper).

## Modes

**Comparative** tests whether sequences differ between two groups (e.g. control vs lesion). The test statistic is a studentized difference in usage rates.

**Correlative** tests whether sequences correlate with a continuous covariate (e.g. a compulsivity score). The test statistic is a studentized Pearson correlation.

## What the algorithm needs

At its core, CBAS needs two things per subject:

1. **A choice stream** — a 1D array of integers representing the sequence of choices (or symbols) that subject made. Each integer is in the range `[0, alphabet_size)`.
2. **A grouping variable** — either a binary group label (comparative mode) or a continuous score (correlative mode).

Everything else (file format, reward encoding, session filtering) is about getting your raw experimental data into that shape.

## Data format

### Option A: bring your own arrays

If you already have choice streams as numpy arrays (e.g. from your own preprocessing pipeline), you can skip `load_subject_data` entirely and pass data directly to `build_count_matrix`. Each element in `subjects_data` should be a 2D array with columns `(session, choice, reward, contingency)`. If you don't use session filtering or reward encoding, you can fill those columns with zeros:

```python
import numpy as np
from pycbas import CBASParams, run_cbas_comparative

# Suppose you have a list of 1D choice arrays
choice_streams = [np.array([0, 1, 1, 0, 1, ...]), ...]  # one per subject

# Wrap each into the expected 4-column format
subjects_data = []
for stream in choice_streams:
    n = len(stream)
    arr = np.zeros((n, 4), dtype=np.int32)
    arr[:, 1] = stream       # column 1 = choice
    arr[:, 3] = 1            # column 3 = contingency (set to match your filter)
    subjects_data.append(arr)

group_labels = [0] * 20 + [1] * 20  # e.g. 20 per group

params = CBASParams(
    num_arms=2,           # size of your choice alphabet
    seq_len_max=4,        # max pattern length to test
    criterion=200,        # use first 200 choices per subject
    resample_number=10000,
)

result = run_cbas_comparative(
    subjects_data, group_labels, params,
    contingency=1,        # must match the value you put in column 3
    encode_reward=False,  # no reward encoding since column 2 is zeros
)
```

For correlative mode, replace `group_labels` with a numpy array of continuous scores (one per subject) and call `run_cbas_correlative`.

### Option B: the CSV loader

If your data is stored as one CSV per subject with columns `session, choice, reward, contingency`, you can use the built-in loader:

```
0,3,1,2
0,2,0,2
0,5,1,2
1,1,0,2
```

- `session` — integer session ID (used only if you want to filter by session)
- `choice` — the arm/symbol chosen (0-indexed integer, must be < `num_arms`)
- `reward` — 0 or 1 (used when `encode_reward=True` to double the alphabet)
- `contingency` — trial condition integer (used to filter trials by condition)

```python
from pycbas import load_subject_data

subjects_data = [load_subject_data(f) for f in my_file_list]
```

### Reward encoding

When `encode_reward=True`, each trial's symbol becomes `choice + reward * num_arms`. This doubles the effective alphabet. Use this when the same physical choice can lead to different outcomes and that distinction matters (e.g. a 6-arm bandit where getting reward vs not changes the behavioral meaning of the choice). Set `encode_reward=False` for deterministic tasks where reward is fully predicted by the choice (like a two-alternative forced choice task).

### How many symbols?

Set `num_arms` to the number of distinct choices your task offers. With `encode_reward=True`, the effective alphabet becomes `num_arms * 2`. Examples:

- Binary maze (left/right, no reward): `num_arms=2`, `encode_reward=False` → 2 symbols
- 6-arm bandit with reward: `num_arms=6`, `encode_reward=True` → 12 symbols
- 4-symbol task (choice × side): `num_arms=4`, `encode_reward=False` → 4 symbols

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
    contingency=2,       # filter to trial condition 2
    encode_reward=True,  # symbol = choice + reward * num_arms
)

print(f"{result.n_significant} significant sequences (k={result.k_final})")
```

## Running a correlative analysis

Use correlative mode when each subject has a continuous measure (a clinical score, age, reaction time, performance metric, etc.) and you want to find sequences whose usage tracks with that measure across subjects.

```python
from pycbas import CBASParams, load_subject_data, run_cbas_correlative
import numpy as np

subjects_data = [load_subject_data(f) for f in subject_files]

# One score per subject, in the same order as subjects_data.
# This can be any continuous measure: clinical scores, ages, performance, etc.
scores = np.array([72.1, 58.3, 85.0, ...])  # length must equal len(subjects_data)

params = CBASParams(
    num_arms=6,
    seq_len_max=4,
    criterion=400,
    resample_number=10000,
)

result = run_cbas_correlative(subjects_data, scores, params)
```

The `scores` array is the covariate. CBAS will test, for every sequence in the count matrix, whether that sequence's usage (across subjects) correlates with these scores. The order must match `subjects_data` — `scores[i]` is the score for `subjects_data[i]`.

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
| `criterion` | 800 | How much of each subject's stream to count. With the default `criterion_order=0` this is a trial index: sequences starting at positions 0 through `criterion` (inclusive) are counted. |
| `criterion_order` | 0 | What `criterion` counts. 0 is a trial index. 1 stops each subject at its `criterion`-th reward. Higher orders stop at the `criterion`-th run of that many consecutive rewarded choices. See [Higher-order criteria](#higher-order-criteria). |
| `resample_number` | 10,000 | Number of bootstrap resamples M. More gives tighter p-values but costs linearly in time and memory. |
| `alpha` | 0.5 | Significance threshold for FDP control. The paper uses 0.5 (median FDP). |
| `gamma` | 0.05 | FDP tolerance. Fraction of rejections allowed to be false. |
| `centering` | False | Whether to center the bootstrap null by subtracting the observed delta. False matches the Igor implementation. True is slightly more liberal. |

Pipeline functions (`run_cbas_comparative`, `run_cbas_correlative`) also accept:

| Parameter | Default | Description |
|-----------|---------|--------|
| `contingency` | 2 | Filter trials by contingency column value, or None for all trials. |
| `encode_reward` | True | Encode reward into symbols (doubles alphabet). Set False for deterministic-outcome tasks. |
| `block_aware` | False | When True, sequences cannot span block/session boundaries. Enable for multi-session experiments where concatenating across sessions would create artificial adjacencies. |

## Choosing parameters

**num_arms** is determined by your task. Binary choice = 2. Six-arm maze = 6.

**encode_reward** should be True when reward is informative and not fully determined by the choice. Set False for deterministic tasks where the choice symbol already encodes the outcome (e.g. a two-alternative forced choice where left always means stimulus A).

**seq_len_max** controls how many sequences exist in the hypothesis space. The total is `sum(A^l for l in 1..L)` where A is the effective alphabet. For 12 symbols and L=6, that's about 3.3 million sequences. In practice, only observed sequences are tested, which is much smaller.

**criterion** should match the minimum usable trial count across your subjects. If some subjects have only 300 trials, set criterion to 300 or less.

**resample_number** of 10,000 is standard. For exploratory work, 1,000 is faster with coarser p-values.

## Higher-order criteria

By default the criterion is a fixed number of trials, which gives every subject the
same amount of data. An alternative is to stop each subject once it has *achieved*
some amount, matching subjects on performance rather than on exposure. This is
optional and off by default.

```python
params = CBASParams(
    num_arms=6,
    seq_len_max=6,
    criterion=100,
    criterion_order=4,   # stop at the 100th run of 4 consecutive rewarded choices
)
```

| `criterion_order` | `criterion` counts | Stops each subject at |
|---|---|---|
| 0 (default) | trials | trial number `criterion` |
| 1 | rewards | its `criterion`-th reward |
| k | runs of k consecutive rewarded choices | its `criterion`-th such run |

Runs are counted in overlapping windows, so seven consecutive rewarded trials contain
four runs of length four. Runs never span a session boundary when `block_aware=True`.

Three consequences worth understanding before using it.

**Subjects contribute unequal amounts of data.** Reaching the same achievement takes
different numbers of trials, so subjects no longer contribute equally. Since CBAS
averages raw counts rather than rates, a subject that took twice as long has roughly
twice the counts. If time-to-criterion differs between your groups, that difference
alone shifts every sequence.

For multi-contingency data use `record_criteria`, which returns one criterion per
subject per contingency, since the criterion applies within each one and a subject can
reach it in one contingency and fall short in another.

**Subjects that never reach the criterion are not truncated.** They contribute every
trial they have, so the weakest subjects contribute the most data. Check for this
before running:

```python
from pycbas import subject_criteria
import numpy as np

criteria = subject_criteria(subjects_data, params, contingency=2, block_aware=True)
reached = criteria[np.isfinite(criteria)]
print(f"{criteria.size - reached.size} of {criteria.size} never reached the criterion")
if reached.size:
    print(f"trials used: {reached.min():.0f} to {reached.max():.0f}")
```

The interactive app shows the same information automatically once a higher order is
selected.

**The achievable count is set by your weakest subject.** If you want every subject to
reach the criterion, `criterion` cannot exceed what the worst performer manages.

Reward information comes from the reward column, so higher orders work regardless of
the `encode_reward` setting.

## Multiple contingencies

When subjects run several task contingencies in sequence, each one can be analysed as
its own set of hypotheses. The same arm sequence under two contingencies then becomes
two columns rather than one, and the multiplicity correction runs over all of them
jointly. With 100 sequences in the first contingency and 200 in the second, you are
testing 300.

This needs a data format that records which contingency each trial belongs to, which
the single-contingency loader does not carry:

```python
from pycbas import (CBASParams, load_cohort_with_contingencies,
                    run_cbas_multicontingency)

records, info = load_cohort_with_contingencies("data/my_cohort")
labels = ...   # 0/1 per subject, aligned with records

params = CBASParams(num_arms=6, seq_len_max=4, criterion=100, criterion_order=4)
result = run_cbas_multicontingency(records, labels, params)

# sequences are (contingency_block, sequence) pairs
for (block, seq), significant in zip(result.sequences, result.significant_mask):
    if significant:
        print(block, seq)
```

Four things behave differently from the single-contingency path.

**A contingency is a block of sessions, not a set of arms.** If a task returns to an
earlier configuration later in training, that recurrence is a separate contingency with
its own columns. Identity comes from the order the blocks appear, so keying on the arms
would silently merge the two. Exploration or pretraining phases are excluded.

**The criterion applies within each contingency.** Each is its own learning episode, so
a subject gets a separate cutoff per contingency, and the trial index is local to it.

**Counting is always session-respecting.** Sequences never span a session boundary
within a contingency.

**Every subject must run every contingency.** A subject missing one has no data for
those columns, and filling zero would assert it never produced those sequences rather
than that it was not measured. Representing that honestly needs missing-value support
in the statistic and the bootstrap, which is not implemented, so the mismatch raises.
Restrict to the shared blocks explicitly if you need to:

```python
from pycbas import shared_contingency_blocks

print(shared_contingency_blocks(records))          # e.g. [1, 2, 3, 4, 5, 6]
result = run_cbas_multicontingency(records, labels, params, blocks=[1, 2])
```

**Watch the hypothesis space.** Counting several contingencies multiplies it, and the
step-down's peak allocation grows with it directly. The null it holds is
`resample_number` by hypotheses, at nine bytes an entry: eight for the magnitude and
one for the direction. At `M=10,000` that is roughly 0.9 GB per ten thousand
hypotheses, so a space in the tens of thousands is a multi-gigabyte run, and raising
`seq_len_max` by one can multiply it several times over.

Measure your own space rather than extrapolating from someone else's cohort, because
the count depends on how much of the sequence space the subjects actually visit:
`build_multicontingency_count_matrix` returns the sequence list, and
`estimate_resources` turns a hypothesis count into a memory and time estimate.

Note that `chunked=True`, the default, does **not** remove this cost. It avoids the
intermediate full-width matrix, roughly halving peak memory, but still allocates
`M × n_valid` for the sorted null and its directions. Reducing `resample_number` or
`seq_len_max` is what actually lowers the floor. Pass the multi-contingency hypothesis
count to `estimate_resources` as `n_observed` rather than relying on a
single-contingency estimate.

## Memory and chunked mode

The main memory cost is the bootstrap null matrix (M rows by number of valid test statistics). For large hypothesis spaces (e.g. rats with 16,378 sequences), this can reach several GB.

The `chunked=True` option (default in `run_cbas_comparative`) generates the bootstrap in row-chunks directly into the sorted submatrix, avoiding the full intermediate allocation. This uses roughly half the peak memory at the cost of about 30% more time.

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

## Working with the count matrix

`build_count_matrix` produces the matrix every later stage reads, so it is where to start for anything CBAS does not do itself: a different statistic, an embedding, an EM fit.

```python
from pathlib import Path
from pycbas import CBASParams, build_count_matrix, decode_sequence, load_subject_data

# One pass over the files, so the ids, the arrays and the labels cannot fall out of step.
subject_ids, subjects_data, group_labels = [], [], []
for path in sorted(Path("data/cohort").glob("*.txt")):
    subject_ids.append(path.stem)
    subjects_data.append(load_subject_data(path))
    group_labels.append(0 if "control" in path.stem else 1)   # your own group rule here

params = CBASParams(num_arms=6, seq_len_max=3, criterion=400)
sequences, counts = build_count_matrix(subjects_data, params, contingency=2)

counts.shape        # (len(subjects_data), len(sequences))
```

**Row `i` is `subjects_data[i]`.** Subjects are never sorted or grouped on the way in, and `build_count_matrix` is not given the group labels, so nothing about the groups can move a row. Membership is applied afterwards, as row indices.

**Column `j` is `sequences[j]`.** Columns are ordered by each sequence's total count over the whole cohort, not by order of appearance, and that order changes when the cohort changes. Keep `sequences` with the matrix and look columns up through it rather than by position. The full rule is in the [API reference](api.md#row-and-column-order).

Print both to see the mapping on your own data:

```python
for i, subject_id in enumerate(subject_ids):
    print(f"row {i}  {subject_id}  group {group_labels[i]}")

for j in range(5):
    print(f"col {j}  {sequences[j]}  {decode_sequence(sequences[j], num_arms=6)}")
```

`decode_sequence` writes a symbol tuple in the published convention, so `(8,)` with six arms reads as `3*`, meaning arm 3, rewarded.

### Reordering rows

Nothing in the matrix records which order was used, so permute every per-subject array in the same step:

```python
import numpy as np

order = np.argsort(group_labels, kind="stable")   # group 0 first, then group 1
counts = counts[order]
subject_ids = [subject_ids[i] for i in order]
group_labels = [group_labels[i] for i in order]
```

The pipeline rejects a label vector whose length does not match the cohort, holds values outside `{0, 1}`, or leaves a group empty. It cannot detect a permuted vector of the right length, so assemble the cohort in one pass as above.

### One contingency at a time

[Multiple contingencies](#multiple-contingencies) covers that format and its loader. Its count matrix labels every column with a `(block, sequence)` pair and keeps the blocks contiguous, so one contingency slices out, with `records` from that section:

```python
from pycbas import build_multicontingency_count_matrix, split_sequence_entry

sequences, counts = build_multicontingency_count_matrix(records, params)

# sequences[j] is a pair here, e.g. (1, (8,)), so select on its block.
in_block_1 = [j for j, entry in enumerate(sequences)
              if split_sequence_entry(entry)[0] == 1]
counts[:, in_block_1]
```

`split_sequence_entry` also gives a sequence's length, which `len(entry)` does not: `len` is 2 for every pair, whatever the sequence length. Given a bare entry from `build_count_matrix` it reports a block of `None`, so a block test against those selects nothing rather than raising.

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
