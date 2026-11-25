"""Where each counted sequence occurred, subject by subject.

The count matrix says how often each subject used a sequence. This module says
where. Every row of the occurrence table is one window the count matrix counted, so
for any subject and sequence the number of rows equals that cell of the matrix.
Nothing is added to the analysis and nothing is summarised. Medians, normalisation
and how to treat subjects that rarely or never use a sequence are analytic choices,
and they are left to the caller.

Positions are reported in the coordinate the criterion uses. That is a 0-based index
into the subject's counted stream, which holds only the trials of the selected
contingency, counted across sessions. In multi-contingency data it restarts in each
contingency. `row` maps a position back to `Subject.trials`.

The windows are found by a scan written separately from `enumerate_sequences`, and
every call checks the two against each other through the count matrix. A
disagreement raises rather than returning positions that do not match the analysis.
"""

import warnings

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from .cohort import CHOICE, CONDITION, REWARD, SESSION, CountMatrix, as_cohort
from .contingency import build_multicontingency_count_matrix, shared_contingency_blocks
from .core import build_count_matrix, reward_blocks
from .criterion import as_enumeration_cutoff, criterion_trial
from .io import extract_choice_stream, extract_choice_streams_by_block
from .params import CBASResult, _counts_digest

__all__ = ["sequence_occurrences"]

_MISMATCH = (
    "the recount does not match {what}. The cohort, the parameters or the counting "
    "settings differ from the ones the analysis used. Check encode_reward, "
    "block_aware, contingency or blocks, criterion and criterion_order, and whether "
    "the cohort was filtered or reordered since.")


def sequence_occurrences(cohort, params, sequences, *, contingency=2,
                         block_aware=False, encode_reward=True, blocks=None,
                         counted=None):
    """Start positions of every counted occurrence of the given sequences.

    Pass the same cohort, parameters and counting settings as the analysis, and pass
    its result as `counted` so the recount can be checked against it.

    Args:
        cohort: Cohort the analysis ran on.
        params: CBASParams the analysis used.
        sequences: entries in the form of `result.sequences`, bare tuples for a
            single-contingency analysis or `(block, tuple)` pairs for a
            multi-contingency one. For the significant sequences of a result,
            `[s for s, m in zip(result.sequences, result.significant_mask) if m]`.
        contingency: condition value the analysis filtered on, or None for all
            trials. Single-contingency only.
        block_aware: whether the analysis kept sequences within sessions.
            Single-contingency only. Multi-contingency counting always does.
        encode_reward: whether the analysis folded reward into the symbols.
        blocks: contingency blocks the analysis included. Multi-contingency only,
            default all shared by every subject.
        counted: the analysis's CBASResult, or a CountMatrix. The recount must match
            it exactly. A result is matched by its `counts_digest`. A result without
            one, from an older version, is matched by its sequence list, which is
            not exact, and a warning says so.

    Returns:
        `(occurrences, subjects)`, two dicts of equal-length numpy arrays, so
        `pandas.DataFrame(occurrences)` gives a table.

        `occurrences` has one row per counted window, ordered by subject, then by
        the order of `sequences`, then by position.

        - `subject_id`
        - `block`, the contingency block, or the `contingency` value for a
          single-contingency analysis
        - `sequence`, the entry as passed
        - `start`, the 0-based position of the window's first trial in the counted
          stream. Its last trial is `start + len(sequence) - 1`, and without
          `block_aware` a window can span two sessions.
        - `row`, the index of that first trial in `Subject.trials`
        - `session`, the session of that first trial

        `subjects` has one row per subject and block, including subjects that never
        used any of the sequences.

        - `subject_id`, `block`
        - `n_trials`, the length of the counted stream
        - `reached`, whether the subject reached the criterion
        - `cutoff`, the last position at which a window may start. It is the
          criterion trial when reached and `n_trials` otherwise, so positions are
          bounded by how quickly each subject reached the criterion.

    Raises:
        ValueError: if the recount does not match `counted`, or a sequence is
            malformed, repeated, longer than `params.seq_len_max`, or was never
            counted under these settings.
    """
    cohort = as_cohort(cohort)
    entries, multi = _normalise(sequences)
    if not multi and blocks is not None:
        raise ValueError("blocks applies to multi-contingency sequences, which are "
                         "(block, sequence) pairs")
    for entry in entries:
        length = len(entry[1] if multi else entry)
        if not 1 <= length <= params.seq_len_max:
            raise ValueError(f"sequence {entry!r} has length {length}, outside 1 to "
                             f"seq_len_max={params.seq_len_max}")

    if multi:
        matrix = build_multicontingency_count_matrix(
            cohort, params, blocks=blocks, encode_reward=encode_reward)
        used_blocks = sorted(blocks) if blocks is not None else \
            shared_contingency_blocks(cohort)
    else:
        matrix = build_count_matrix(cohort, params, contingency=contingency,
                                    encode_reward=encode_reward,
                                    block_aware=block_aware)
        used_blocks = [contingency]
    _verify(counted, matrix)

    column = {entry: j for j, entry in enumerate(matrix.sequences)}
    missing = [entry for entry in entries if entry not in column]
    if missing:
        raise ValueError(
            f"{len(missing)} sequence(s) were never counted under these settings, "
            f"e.g. {missing[:3]!r}. A wrong encode_reward or num_arms is the usual "
            f"cause.")

    targets = {}
    for entry in entries:
        block, seq = entry if multi else (contingency, entry)
        targets.setdefault(block, []).append(seq)

    occ = {"subject_id": [], "block": [], "sequence": [], "start": [], "row": [],
           "session": []}
    subj = {"subject_id": [], "block": [], "n_trials": [], "reached": [],
            "cutoff": []}
    order = getattr(params, "criterion_order", 0)
    matrix_row = {sid: i for i, sid in enumerate(matrix.subject_ids)}

    for subject in cohort:
        found = {}
        for block in used_blocks:
            if multi:
                symbols, rows, rewards = _block_segments(
                    subject, block, params.num_arms, encode_reward)
            else:
                symbols, rows, rewards = _single_segments(
                    subject, contingency, block_aware, params.num_arms, encode_reward)
            _check_rows(subject, symbols, rows, params.num_arms, encode_reward)

            n_trials = sum(len(s) for s in symbols)
            criterion = criterion_trial(rewards, order, params.criterion)
            cutoff = as_enumeration_cutoff(criterion, n_trials)
            subj["subject_id"].append(subject.id)
            subj["block"].append(block)
            subj["n_trials"].append(n_trials)
            subj["reached"].append(bool(np.isfinite(criterion)))
            subj["cutoff"].append(cutoff)

            for seq, hits in _scan(symbols, rows, targets.get(block, []),
                                   cutoff).items():
                found[(block, seq) if multi else seq] = hits

        counts = matrix.counts[matrix_row[subject.id]]
        for entry in entries:
            starts, hit_rows = found[entry]
            if len(starts) != counts[column[entry]]:
                raise RuntimeError(
                    f"found {len(starts)} occurrences of {entry!r} for subject "
                    f"{subject.id!r} where the count matrix has "
                    f"{counts[column[entry]]:g}. This is a bug in pycbas.")
            n = len(starts)
            occ["subject_id"].extend([subject.id] * n)
            occ["block"].extend([entry[0] if multi else contingency] * n)
            occ["sequence"].extend([entry] * n)
            occ["start"].append(starts)
            occ["row"].append(hit_rows)
            occ["session"].append(subject.trials[hit_rows, SESSION])

    return _columns(occ, ("start", "row", "session")), \
        _columns(subj, ("n_trials", "cutoff"), bools=("reached",))


def _normalise(sequences):
    """Entries as hashable tuples of builtin ints, and whether they are paired."""
    entries, kinds = [], set()
    for entry in sequences:
        paired = len(entry) == 2 and isinstance(entry[1], (tuple, list, np.ndarray))
        kinds.add(paired)
        if paired:
            entries.append((int(entry[0]), tuple(int(s) for s in entry[1])))
        else:
            entries.append(tuple(int(s) for s in entry))
    if not entries:
        raise ValueError("sequences is empty")
    if len(kinds) > 1:
        raise ValueError("sequences mixes bare sequences with (block, sequence) pairs")
    if len(set(entries)) != len(entries):
        raise ValueError("sequences contains repeated entries")
    return entries, kinds.pop()


def _verify(counted, matrix):
    if counted is None:
        return
    if isinstance(counted, CountMatrix):
        same = (list(counted.subject_ids) == list(matrix.subject_ids)
                and list(counted.sequences) == list(matrix.sequences)
                and np.array_equal(counted.counts, matrix.counts))
        what = "the given CountMatrix"
    elif isinstance(counted, CBASResult):
        if counted.counts_digest is not None:
            same = counted.counts_digest == _counts_digest(
                matrix.counts, matrix.subject_ids, matrix.sequences)
        else:
            warnings.warn(
                "this CBASResult has no counts_digest, so the recount is checked "
                "against its sequence list only. That check can pass when the "
                "settings differ, especially with few arms and short sequences. "
                "Re-run the analysis with this version for an exact check.",
                UserWarning, stacklevel=3)
            same = list(counted.sequences) == list(matrix.sequences)
        what = "the given CBASResult"
    else:
        raise TypeError(f"counted must be a CBASResult or a CountMatrix, got "
                        f"{type(counted).__name__}")
    if not same:
        raise ValueError(_MISMATCH.format(what=what))


def _single_segments(subject, contingency, block_aware, num_arms, encode_reward):
    """Symbol and row segments as `build_count_matrix` sees them, with rewards."""
    trials = subject.trials
    if contingency is None:
        mask = np.ones(len(trials), dtype=bool)
    else:
        mask = trials[:, CONDITION] == contingency
    if block_aware:
        symbols = extract_choice_streams_by_block(trials, contingency, num_arms,
                                                  encode_reward=encode_reward)
        rows = _rows_by_session(trials, mask)
    else:
        symbols = [extract_choice_stream(trials, contingency, num_arms,
                                         encode_reward=encode_reward)]
        rows = [np.flatnonzero(mask)]
    return symbols, rows, reward_blocks(trials, contingency, block_aware)


def _block_segments(subject, block, num_arms, encode_reward):
    """Segments as `build_multicontingency_count_matrix` sees one block."""
    mask = subject.trials[:, CONDITION] == block
    return (subject.symbol_blocks_for(block, num_arms, encode_reward),
            _rows_by_session(subject.trials, mask),
            subject.reward_blocks_for(block))


def _rows_by_session(trials, mask):
    """Row indices per session, in the order `np.unique` gives the counting code."""
    sessions = trials[:, SESSION]
    return [np.flatnonzero(mask & (sessions == s)) for s in np.unique(sessions[mask])]


def _check_rows(subject, symbols, rows, num_arms, encode_reward):
    """The row map must reproduce the stream the counting code used."""
    trials = subject.trials
    for sym, idx in zip(symbols, rows, strict=True):
        rebuilt = trials[idx, CHOICE]
        if encode_reward:
            rebuilt = rebuilt + trials[idx, REWARD] * num_arms
        if not np.array_equal(rebuilt, sym):
            raise RuntimeError(f"row map for subject {subject.id!r} does not match "
                               f"its counted stream. This is a bug in pycbas.")


def _scan(symbols, rows, targets, cutoff):
    """Counted start positions and rows of each target, by the enumeration rules.

    A window starting at global position p is counted when p <= cutoff and the
    window lies inside one segment. Global positions run across segments, and every
    segment adds its full length, including segments too short for the window.
    """
    by_length = {}
    for seq in targets:
        by_length.setdefault(len(seq), []).append(seq)
    starts = {seq: [] for seq in targets}
    hit_rows = {seq: [] for seq in targets}
    offset = 0
    for sym, idx in zip(symbols, rows):
        if offset > cutoff:
            break
        for length, seqs in by_length.items():
            if len(sym) < length:
                continue
            windows = sliding_window_view(np.asarray(sym), length)
            for seq in seqs:
                hits = np.flatnonzero((windows == np.asarray(seq)).all(axis=1))
                hits = hits[offset + hits <= cutoff]
                starts[seq].append(offset + hits)
                hit_rows[seq].append(idx[hits])
        offset += len(sym)
    return {seq: (_concat(starts[seq]), _concat(hit_rows[seq])) for seq in targets}


def _concat(parts):
    return np.concatenate(parts).astype(np.int64) if parts else \
        np.empty(0, dtype=np.int64)


def _columns(table, int_arrays, bools=()):
    """Lists of values, or of int arrays for `int_arrays`, as numpy columns."""
    out = {}
    for key, values in table.items():
        if key in int_arrays and values and isinstance(values[0], np.ndarray):
            out[key] = _concat(values)
        elif key in int_arrays:
            out[key] = np.asarray(values, dtype=np.int64)
        elif key in bools:
            out[key] = np.asarray(values, dtype=bool)
        else:
            column = np.empty(len(values), dtype=object)
            column[:] = values
            out[key] = column
    return out
