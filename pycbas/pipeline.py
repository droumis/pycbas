"""High-level CBAS pipeline functions."""

import numpy as np
from .params import CBASParams, CBASResult
from .core import build_count_matrix, compute_test_stats, compute_test_stats_correlative
from .bootstrap import bootstrap_test_stats, bootstrap_test_stats_correlative
from .stepdown import find_k_fwer, find_k_fwer_chunked


def run_cbas_multicontingency(records, group_labels, params=None, blocks=None,
                              encode_reward=True, chunked=True):
    """Comparative CBAS across several contingencies, each counted separately.

    Sequences from different contingencies are distinct hypotheses, so the same
    arm sequence under two contingencies contributes two columns, and the
    multiplicity correction runs over all of them jointly. Everything downstream
    of the count matrix is the ordinary comparative path, unchanged.

    Args:
        records: list of SubjectRecord from `load_subject_data_with_contingencies`
        group_labels: array of 0/1 indicating group membership
        params: CBASParams instance; the criterion applies within each contingency
        blocks: contingency blocks to include, default all shared by every subject
        encode_reward: encode reward into symbols
        chunked: use the memory-efficient chunked pipeline. Strongly advised here,
            since counting several contingencies multiplies the hypothesis space.

    Returns:
        CBASResult, whose `sequences` entries are (block, sequence_tuple) pairs
        rather than bare tuples.
    """
    from .contingency import build_multicontingency_count_matrix

    if params is None:
        params = CBASParams()

    group_labels = np.asarray(group_labels)
    group_indices = [
        np.where(group_labels == 0)[0],
        np.where(group_labels == 1)[0],
    ]

    sequences, count_matrix = build_multicontingency_count_matrix(
        records, params, blocks=blocks, encode_reward=encode_reward)

    return _finish_comparative(sequences, count_matrix, group_indices, params,
                               chunked)


def run_cbas_comparative(subjects_data, group_labels, params=None,
                         contingency=2, encode_reward=True, chunked=True,
                         block_aware=False):
    """Run the full comparative CBAS pipeline.

    Args:
        subjects_data: list of subject data arrays (from load_subject_data)
        group_labels: array of 0/1 indicating group membership
        params: CBASParams instance
        contingency: block type to filter on, or None for all trials
        encode_reward: if True, encode symbol + reward*num_arms. Set False for 2AFC.
        chunked: if True (default), use memory-efficient chunked pipeline
        block_aware: if True, sequences cannot span block/session boundaries.

    Returns:
        CBASResult
    """
    if params is None:
        params = CBASParams()

    group_labels = np.asarray(group_labels)
    group_indices = [
        np.where(group_labels == 0)[0],
        np.where(group_labels == 1)[0],
    ]

    sequences, count_matrix = build_count_matrix(subjects_data, params,
                                                 contingency=contingency,
                                                 encode_reward=encode_reward,
                                                 block_aware=block_aware)
    return _finish_comparative(sequences, count_matrix, group_indices, params,
                               chunked)


def _finish_comparative(sequences, count_matrix, group_indices, params, chunked):
    """Statistic, null, step-down and k-FWER for a comparative count matrix.

    Shared by the single- and multi-contingency entry points so both run through
    exactly the same validated numerical path, differing only in how the count
    matrix was built and how its columns are labelled.
    """
    test_stats = compute_test_stats(count_matrix, group_indices)

    if chunked:
        g_values, k_final, k_history = find_k_fwer_chunked(
            test_stats, count_matrix, group_indices, params, return_history=True)
    else:
        null_matrix, null_directions = bootstrap_test_stats(count_matrix, group_indices, params)
        g_values, k_final, k_history = find_k_fwer(
            test_stats, null_matrix, params.alpha, params.gamma,
            null_directions=null_directions, return_history=True)

    significant = np.zeros(len(sequences), dtype=bool)
    for i in range(len(sequences)):
        pos_p = g_values[i * 2]
        neg_p = g_values[i * 2 + 1]
        if (not np.isnan(pos_p) and pos_p < params.alpha) or \
           (not np.isnan(neg_p) and neg_p < params.alpha):
            significant[i] = True

    return CBASResult(
        sequences=sequences,
        test_stats=test_stats,
        g_values=g_values,
        k_final=k_final,
        significant_mask=significant,
        k_history=k_history,
    )


def run_cbas_correlative(subjects_data, covariate, params=None,
                         contingency=2, encode_reward=True, block_aware=False):
    """Run the full correlative CBAS pipeline.

    Args:
        subjects_data: list of subject data arrays (from load_subject_data)
        covariate: array of continuous values (e.g. CBIT scores), one per subject
        params: CBASParams instance
        contingency: block type to filter on, or None for all trials
        encode_reward: if True, encode symbol + reward*num_arms. Set False for 2AFC.
        block_aware: if True, sequences cannot span block/session boundaries.

    Returns:
        CBASResult
    """
    if params is None:
        params = CBASParams()

    covariate = np.asarray(covariate, dtype=np.float64)
    sequences, count_matrix = build_count_matrix(subjects_data, params,
                                                 contingency=contingency,
                                                 encode_reward=encode_reward,
                                                 block_aware=block_aware)
    test_stats = compute_test_stats_correlative(count_matrix, covariate)
    null_matrix, null_directions = bootstrap_test_stats_correlative(count_matrix, covariate, params)
    g_values, k_final, k_history = find_k_fwer(
        test_stats, null_matrix, params.alpha, params.gamma,
        null_directions=null_directions, return_history=True)

    significant = np.zeros(len(sequences), dtype=bool)
    for i in range(len(sequences)):
        pos_p = g_values[i * 2]
        neg_p = g_values[i * 2 + 1]
        if (not np.isnan(pos_p) and pos_p < params.alpha) or \
           (not np.isnan(neg_p) and neg_p < params.alpha):
            significant[i] = True

    return CBASResult(
        sequences=sequences,
        test_stats=test_stats,
        g_values=g_values,
        k_final=k_final,
        significant_mask=significant,
        k_history=k_history,
    )
