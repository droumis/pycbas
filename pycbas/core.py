"""Count matrix construction and test statistics."""

import numpy as np
from .io import (extract_choice_stream, extract_choice_streams_by_block,
                 enumerate_sequences, enumerate_sequences_block_aware)


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
    all_seq_counts = []
    for subj_data in subjects_data:
        subj_counts = {}
        if block_aware:
            block_streams = extract_choice_streams_by_block(
                subj_data, contingency, params.num_arms, encode_reward=encode_reward)
            for seq_len in range(1, params.seq_len_max + 1):
                seq_counts = enumerate_sequences_block_aware(
                    block_streams, seq_len, params.criterion)
                subj_counts.update(seq_counts)
        else:
            stream = extract_choice_stream(subj_data, contingency, params.num_arms,
                                           encode_reward=encode_reward)
            for seq_len in range(1, params.seq_len_max + 1):
                seq_counts = enumerate_sequences(stream, seq_len, params.criterion)
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

    counts0 = count_matrix[grp0]
    counts1 = count_matrix[grp1]

    n0 = len(grp0)
    n1 = len(grp1)

    mean0 = counts0.mean(axis=0)
    mean1 = counts1.mean(axis=0)
    sem0 = counts0.std(axis=0, ddof=1) / np.sqrt(n0)
    sem1 = counts1.std(axis=0, ddof=1) / np.sqrt(n1)

    delta = mean0 - mean1
    sigma = np.sqrt(sem0**2 + sem1**2)

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

        tau_num = np.sqrt(np.sum(X_dev ** 2 * Y_dev ** 2) / n)
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
