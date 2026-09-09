"""Count matrix construction and test statistics."""

import warnings

import numpy as np
from .io import (extract_choice_stream, extract_choice_streams_by_block,
                 enumerate_sequences, enumerate_sequences_block_aware)
from .criterion import criterion_trial, as_enumeration_cutoff
from ._moments import sigma_from_sums

#: One warning per session is enough: this is a property of the caller's data,
#: not of any single call.
_WARNED_NON_INTEGER = False


def reward_blocks(subj_data, contingency=2, block_aware=False):
    """Per-subject reward arrays in the same coordinate system as enumeration.

    The criterion index has to mean the same thing as the enumeration's start
    position, so the block structure must match. With `block_aware` the stream is
    split by session, exactly as `extract_choice_streams_by_block` splits it, and
    a run of rewarded trials cannot span a session. Without it the stream is one
    block, matching the pooled enumeration, where sequences may span sessions.
    """
    if contingency is None:
        data = subj_data
    else:
        data = subj_data[subj_data[:, 3] == contingency]

    rewards = data[:, 2]
    if not block_aware:
        return [rewards]
    sessions = data[:, 0]
    return [rewards[sessions == s] for s in np.unique(sessions)]


def subject_criteria(subjects_data, params, contingency=2, block_aware=False):
    """Per-subject criterion trial index, with `inf` where a subject fell short.

    Exposed because the shortfall is worth reporting rather than absorbing. A
    subject with an infinite criterion contributes every window it has, so with a
    higher-order criterion the weakest subjects contribute the most data, and the
    shortfall is often uneven across groups. Check it on your own data rather than
    assuming it is negligible.

    Returns:
        float array of length n_subjects.
    """
    order = getattr(params, "criterion_order", 0)
    return np.array([
        criterion_trial(reward_blocks(d, contingency, block_aware),
                        order, params.criterion)
        for d in subjects_data
    ], dtype=np.float64)


def build_count_matrix(subjects_data, params, contingency=2, encode_reward=True,
                       block_aware=False):
    """Build the full sequence count matrix.

    Args:
        subjects_data: list of subject data arrays (from load_subject_data)
        params: CBASParams instance
        contingency: block type to filter on, or None for all trials
        encode_reward: if True, encode symbol + reward*num_arms. Set False for
            tasks where outcome is deterministic from the symbol (e.g., 2AFC).
        block_aware: if True, sequences cannot span block/session boundaries.
            Matches Igor's counting for multi-session experiments.

    Returns:
        sequences: list of all unique sequence tuples (sorted by total frequency descending)
        count_matrix: ndarray of shape (n_subjects, n_sequences) with usage counts
    """
    n_subjects = len(subjects_data)
    order = getattr(params, "criterion_order", 0)
    all_seq_counts = []
    for subj_data in subjects_data:
        subj_counts = {}
        if block_aware:
            block_streams = extract_choice_streams_by_block(
                subj_data, contingency, params.num_arms, encode_reward=encode_reward)
            n_trials = sum(len(b) for b in block_streams)
            cutoff = as_enumeration_cutoff(
                criterion_trial(reward_blocks(subj_data, contingency, True),
                                order, params.criterion),
                n_trials)
            for seq_len in range(1, params.seq_len_max + 1):
                seq_counts = enumerate_sequences_block_aware(
                    block_streams, seq_len, cutoff)
                subj_counts.update(seq_counts)
        else:
            stream = extract_choice_stream(subj_data, contingency, params.num_arms,
                                           encode_reward=encode_reward)
            cutoff = as_enumeration_cutoff(
                criterion_trial(reward_blocks(subj_data, contingency, False),
                                order, params.criterion),
                len(stream))
            for seq_len in range(1, params.seq_len_max + 1):
                seq_counts = enumerate_sequences(stream, seq_len, cutoff)
                subj_counts.update(seq_counts)
        all_seq_counts.append(subj_counts)

    all_sequences = set()
    for sc in all_seq_counts:
        all_sequences.update(sc.keys())

    seq_totals = {}
    for seq in all_sequences:
        seq_totals[seq] = sum(sc.get(seq, 0) for sc in all_seq_counts)

    sequences = sorted(seq_totals.keys(), key=lambda s: (-seq_totals[s], len(s), s))

    seq_to_idx = {s: i for i, s in enumerate(sequences)}
    count_matrix = np.zeros((n_subjects, len(sequences)), dtype=np.float64)
    for subj_idx, sc in enumerate(all_seq_counts):
        for seq, count in sc.items():
            count_matrix[subj_idx, seq_to_idx[seq]] = count

    return sequences, count_matrix


def compute_test_stats(count_matrix, group_indices):
    """Compute studentized two-sample test statistics for all sequences.

    Uses two one-tailed tests per sequence (type III error handling).
    Returns array of shape (n_sequences * 2,) where:
      - even indices: positive direction (group0 > group1)
      - odd indices: negative direction (group1 > group0)
    NaN where the test stat is not in that direction or is undefined.
    """
    grp0 = group_indices[0]
    grp1 = group_indices[1]

    counts0 = np.ascontiguousarray(count_matrix[grp0], dtype=np.float64)
    counts1 = np.ascontiguousarray(count_matrix[grp1], dtype=np.float64)

    # Warned here rather than in the step-down because this is the one place every
    # path passes through that can see the matrix. `find_k_fwer` only ever receives
    # the null, so it cannot detect this for itself.
    global _WARNED_NON_INTEGER
    if not _WARNED_NON_INTEGER:
        finite = count_matrix[np.isfinite(count_matrix)]
        if finite.size and not np.all(finite == np.rint(finite)):
            _WARNED_NON_INTEGER = True
            warnings.warn(
                "count matrix is not integer-valued. Observed and bootstrap "
                "statistics are then not guaranteed to agree bitwise, so the "
                "step-down's `null >= observed` can discard the resamples that "
                "represent the observed value, which adds false positives. Pass "
                "tie_rtol=pycbas._moments.tie_rtol_for(matrix) to the step-down, "
                "or pass the integer count matrix instead if the normalising "
                "denominator is common to every subject.",
                RuntimeWarning, stacklevel=2)

    n0 = len(grp0)
    n1 = len(grp1)

    # Built from raw sums, and combined by `_moments.sigma_from_sums`, so that
    # this agrees bitwise with the bootstrap for any sequence whose resampled
    # multiset has the same sums. See pycbas/_moments.py for why that matters:
    # the step-down's `>=` is an equality test on a large block of resamples
    # whenever a statistic is an exact small rational, which is common for rare
    # sequences. Do not rewrite this as `.mean()` and `.std()`.
    sum0 = counts0.sum(axis=0)
    sum1 = counts1.sum(axis=0)
    sq0 = (counts0 * counts0).sum(axis=0)
    sq1 = (counts1 * counts1).sum(axis=0)

    delta = sum0 / n0 - sum1 / n1
    sigma = sigma_from_sums(sum0, sq0, float(n0), sum1, sq1, float(n1))

    n_seq = count_matrix.shape[1]
    stats = np.full(n_seq * 2, np.nan)

    valid = (sigma > 0) & (delta != 0)
    safe_sigma = np.where(sigma > 0, sigma, 1.0)
    t_vals = np.where(valid, delta / safe_sigma, np.nan)

    pos_mask = valid & (delta > 0)
    neg_mask = valid & (delta < 0)
    stats[0::2] = np.where(pos_mask, t_vals, np.nan)
    stats[1::2] = np.where(neg_mask, -t_vals, np.nan)

    return stats


def compute_test_stats_correlative(count_matrix, covariate):
    """Compute studentized correlation test statistics (eq. 2-4 in paper).

    For each sequence, computes the studentized Pearson correlation between
    that sequence's usage counts across subjects and the covariate (e.g. CBIT).

    Uses two one-tailed tests: positive correlation and negative correlation.
    Returns array of shape (n_sequences * 2,) where:
      - even indices: positive correlation (rho > 0)
      - odd indices: negative correlation (rho < 0)
    """
    n = count_matrix.shape[0]
    n_seq = count_matrix.shape[1]
    Y = np.asarray(covariate, dtype=np.float64)
    Y_bar = Y.mean()
    Y_dev = Y - Y_bar
    ss_Y = np.sum(Y_dev ** 2)

    stats = np.full(n_seq * 2, np.nan)

    for s in range(n_seq):
        X = count_matrix[:, s]
        X_bar = X.mean()
        X_dev = X - X_bar
        ss_X = np.sum(X_dev ** 2)

        if ss_X == 0 or ss_Y == 0:
            continue

        rho = (np.sum(X * Y) - n * X_bar * Y_bar) / np.sqrt(ss_X * ss_Y)

        # (x_dev * y_dev) ** 2, not (x_dev ** 2) * (y_dev ** 2): the two are not
        # bitwise equal, and `_bootstrap_correlative_parallel` uses the former.
        # The correlative null is a permutation null, so permuting only within a
        # tied block of the covariate gives a mathematically identical statistic;
        # the step-down's `>=` then decides on rounding. See pycbas/_moments.py.
        tau_num = np.sqrt(np.sum((X_dev * Y_dev) ** 2) / n)
        tau_den = np.sqrt(ss_X / n) * np.sqrt(ss_Y / n)
        tau = tau_num / tau_den

        if tau == 0:
            continue

        t_val = np.sqrt(n) * rho / tau

        if rho > 0:
            stats[s * 2] = t_val
        elif rho < 0:
            stats[s * 2 + 1] = -t_val

    return stats
