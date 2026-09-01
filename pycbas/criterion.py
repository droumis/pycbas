"""Higher-order criteria: where to stop counting for each subject.

The criterion decides how much of each subject's stream enters the count matrix.
`enumerate_sequences` and `enumerate_sequences_block_aware` take it as a maximum
start position, so a criterion is always ultimately a trial index.

Orders, following Kastner's Igor implementation:

    order 0   a fixed trial index, the only mode pycbas supported historically
    order 1   the trial where the subject has earned `count` rewards
    order k   the trial where the subject has produced `count` runs of k
              consecutive rewarded choices

Order 1 is the k=1 case of order k, since a single rewarded trial is a run of
length one, so orders 1 and above share one code path.

Two properties are easy to get wrong and are both pinned by tests against the
reference Igor output:

Runs are counted in **overlapping** windows. Seven rewarded trials in a row
contain four distinct runs of length four, not one.

The returned index counts **trials**, not windows, and it counts across session
boundaries. A session contributes `len(session)` to the running position even
though only `len(session) - k + 1` of those positions can start a run. Counting
windows instead shifts every value by a few percent, which is small enough to
look plausible and wrong everywhere.

A subject who never reaches `count` gets `inf`, which the enumeration reads as
"no cutoff", so that subject contributes every window it has. Kastner confirmed
this is current intended behaviour, with the caveat that a fixed criterion is
only strictly appropriate when every subject in every group reaches it. Callers
that want those subjects excluded should test for `inf` themselves;
`reached_criterion` is provided for that.
"""

import numpy as np

__all__ = [
    "perfect_run_starts",
    "criterion_trial",
    "reached_criterion",
    "criterion_trials_by_subject",
]


def perfect_run_starts(reward_blocks, order):
    """Trial indices that start a run of `order` consecutive rewarded trials.

    Args:
        reward_blocks: list of per-session 0/1 reward arrays, in temporal order.
            Runs never span a session boundary.
        order: run length, >= 1.

    Returns:
        Sorted array of trial indices, where the index counts every trial in the
        concatenation of blocks, not just those that can start a run.
    """
    if order < 1:
        raise ValueError(f"order must be >= 1 for run counting, got {order}")

    starts = []
    position = 0
    for block in reward_blocks:
        rewards = np.asarray(block).astype(bool)
        n = len(rewards)
        if n >= order:
            # rolling AND over a window of `order`: a run is perfect when the
            # number of rewarded trials in the window equals the window length
            cumulative = np.concatenate([[0], np.cumsum(rewards)])
            window_sums = cumulative[order:] - cumulative[:-order]
            local = np.flatnonzero(window_sums == order)
            starts.append(local + position)
        position += n

    if not starts:
        return np.empty(0, dtype=np.int64)
    return np.concatenate(starts).astype(np.int64)


def criterion_trial(reward_blocks, order, count):
    """Trial index at which `count` units of the given order have accumulated.

    Args:
        reward_blocks: list of per-session 0/1 reward arrays, in temporal order.
        order: 0 for a fixed trial index, >= 1 for runs of that many rewarded
            trials.
        count: for order 0 the trial index itself, otherwise how many runs are
            required.

    Returns:
        A trial index, or `inf` when the subject never accumulates `count`.
        `inf` means the enumeration applies no cutoff.
    """
    if count < 1:
        raise ValueError(f"count must be >= 1, got {count}")

    if order == 0:
        total_trials = sum(len(b) for b in reward_blocks)
        # Matching getZerothOrder: a subject with fewer trials than the criterion
        # gets no cutoff. That is equivalent to the historical pycbas behaviour,
        # where min(criterion, len - seq_len) already admitted every window.
        return count if total_trials >= count else np.inf

    starts = perfect_run_starts(reward_blocks, order)
    if len(starts) < count:
        return np.inf
    return int(starts[count - 1])


def reached_criterion(criterion):
    """Whether a criterion value represents a subject that actually reached it."""
    return bool(np.isfinite(criterion))


def criterion_trials_by_subject(reward_blocks_by_subject, order, count):
    """Apply `criterion_trial` across subjects.

    Args:
        reward_blocks_by_subject: sequence of per-subject block lists.
        order, count: as for `criterion_trial`.

    Returns:
        float array of length n_subjects, holding trial indices and `inf`.
        Float rather than int so that `inf` survives.
    """
    return np.array(
        [criterion_trial(blocks, order, count) for blocks in reward_blocks_by_subject],
        dtype=np.float64,
    )
