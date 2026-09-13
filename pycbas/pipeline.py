"""High-level CBAS pipeline functions."""

import numpy as np
from .params import CBASParams, CBASResult
from .core import build_count_matrix, compute_test_stats, compute_test_stats_correlative
from .bootstrap import bootstrap_test_stats, bootstrap_test_stats_correlative
from .stepdown import find_k_fwer, find_k_fwer_chunked


def _group_indices(group_labels, cohort, row_ids):
    """Row indices for each group, resolved through the ids of those rows.

    Resolving by id rather than by position is what makes the row order of the count
    matrix irrelevant here: whatever order the builder produced, each row is grouped
    by the label belonging to *that* subject.

    Every one of these used to produce a complete-looking CBASResult from a cohort
    that was not the one the caller passed, so they fail here instead:

    - too few labels: `np.where` simply found fewer rows, so an analysis of five of
      six subjects was indistinguishable from a finished one
    - a label outside {0, 1}: that subject is in neither group, so it is dropped
      just as silently, which is what a cohort coded 1/2 instead of 0/1 looks like
    - an empty group: the studentized statistic divides by that group's size, so
      every value is NaN, reported only as a RuntimeWarning

    Callers with subjects to exclude should drop them from the cohort, with
    `Cohort.filter`, rather than labelling them into neither group.
    """
    from .cohort import resolve_labels
    if isinstance(group_labels, str):
        labels = resolve_labels(dict(zip(cohort.ids,
                                         cohort.labels_from(group_labels).tolist())),
                                row_ids)
    else:
        labels = resolve_labels(group_labels, row_ids, cohort_ids=cohort.ids)

    n_subjects = len(row_ids)
    groups = [np.where(labels == 0)[0], np.where(labels == 1)[0]]
    assigned = len(groups[0]) + len(groups[1])
    if assigned != n_subjects:
        stray = sorted({v.item() for v in np.unique(labels)} - {0, 1})
        raise ValueError(
            f"group labels must be 0 or 1, but {n_subjects - assigned} of "
            f"{n_subjects} subjects are labelled {stray}, which puts them in "
            f"neither group; drop them from the cohort instead")
    for group, index in zip(groups, (0, 1)):
        if len(group) == 0:
            raise ValueError(
                f"group {index} has no subjects; a comparative analysis needs both")
    return groups


def _check_covariate(covariate, cohort, row_ids):
    """The correlative counterpart of `_group_indices`, also resolved by id.

    A short covariate raised from a broadcast deep in the statistic, naming array
    shapes rather than the mistake; a long one was truncated by the permutation
    indices without comment.
    """
    from .cohort import resolve_labels
    if isinstance(covariate, str):
        values = resolve_labels(dict(zip(cohort.ids,
                                         cohort.covariate_from(covariate).tolist())),
                                row_ids, what="covariate values")
    else:
        values = resolve_labels(covariate, row_ids, cohort_ids=cohort.ids,
                                what="covariate values")
    return np.asarray(values, dtype=np.float64)


def run_cbas_multicontingency(cohort, group_labels, params=None, blocks=None,
                              encode_reward=True, chunked=True):
    """Comparative CBAS across several contingencies, each counted separately.

    Sequences from different contingencies are distinct hypotheses, so the same
    arm sequence under two contingencies contributes two columns, and the
    multiplicity correction runs over all of them jointly. Everything downstream
    of the count matrix is the ordinary comparative path, unchanged.

    Args:
        cohort: Cohort from `load_cohort_with_contingencies`
        group_labels: how to group the cohort. A `{subject_id: 0/1}` mapping, a
            sequence in cohort order, or the name of a `meta` column to derive it
            from. Every subject must be 0 or 1 and neither group may be empty.
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

    matrix = build_multicontingency_count_matrix(
        cohort, params, blocks=blocks, encode_reward=encode_reward)
    group_indices = _group_indices(group_labels, cohort, matrix.subject_ids)

    return _finish_comparative(matrix.sequences, matrix.counts, group_indices,
                               params, chunked)


def run_cbas_comparative(cohort, group_labels, params=None,
                         contingency=2, encode_reward=True, chunked=True,
                         block_aware=False):
    """Run the full comparative CBAS pipeline.

    Args:
        cohort: Cohort of subjects (from `load_cohort`)
        group_labels: how to group the cohort. A `{subject_id: 0/1}` mapping, a
            sequence in cohort order, or the name of a `meta` column to derive it
            from. Every subject must be 0 or 1 and neither group may be empty.
        params: CBASParams instance
        contingency: condition value to filter on, or None for all trials
        encode_reward: if True, encode symbol + reward*num_arms. Set False for 2AFC.
        chunked: if True (default), use memory-efficient chunked pipeline
        block_aware: if True, sequences cannot span block/session boundaries.

    Returns:
        CBASResult
    """
    if params is None:
        params = CBASParams()

    matrix = build_count_matrix(cohort, params, contingency=contingency,
                                encode_reward=encode_reward,
                                block_aware=block_aware)
    group_indices = _group_indices(group_labels, cohort, matrix.subject_ids)
    return _finish_comparative(matrix.sequences, matrix.counts, group_indices,
                               params, chunked)


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


def run_cbas_correlative(cohort, covariate, params=None,
                         contingency=2, encode_reward=True, block_aware=False):
    """Run the full correlative CBAS pipeline.

    Args:
        cohort: Cohort of subjects (from `load_cohort`)
        covariate: continuous value per subject: a `{subject_id: value}` mapping, a
            sequence in cohort order, or the name of a `meta` column.
        params: CBASParams instance
        contingency: block type to filter on, or None for all trials
        encode_reward: if True, encode symbol + reward*num_arms. Set False for 2AFC.
        block_aware: if True, sequences cannot span block/session boundaries.

    Returns:
        CBASResult
    """
    if params is None:
        params = CBASParams()

    matrix = build_count_matrix(cohort, params, contingency=contingency,
                                encode_reward=encode_reward,
                                block_aware=block_aware)
    sequences, count_matrix = matrix.sequences, matrix.counts
    covariate = _check_covariate(covariate, cohort, matrix.subject_ids)
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
