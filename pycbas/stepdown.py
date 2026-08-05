"""Romano-Wolf step-down and k-FWER iteration."""

import numpy as np
from ._numba import njit, prange
from .bootstrap import _bootstrap_chunk_into


@njit(cache=True, parallel=True)
def _stepdown_core(sorted_stats, null_sub, k, max_pval):
    """Numba-accelerated inner loop of Romano-Wolf step-down.

    Args:
        sorted_stats: (n_valid,) test stats sorted descending
        null_sub: (M, n_valid) null matrix columns, NaN replaced with -inf,
                  columns ordered to match sorted_stats
        k: k-FWER parameter
        max_pval: stop computing once p-value reaches this (1.0 for full, alpha for counting)

    Returns:
        p_values: (n_valid,) adjusted p-values with monotonicity enforced
    """
    n_valid = len(sorted_stats)
    M = null_sub.shape[0]
    p_values = np.empty(n_valid)
    active = np.ones(n_valid, dtype=np.bool_)
    last_p = 0.0

    for step in range(n_valid):
        n_active = int(active.sum())
        if n_active <= k - 1:
            comparison_val = np.full(M, -np.inf)
        else:
            comparison_val = np.empty(M)
            for m in prange(M):
                buf = np.empty(k)
                buf[:] = -np.inf
                for col in range(n_valid):
                    if not active[col]:
                        continue
                    val = null_sub[m, col]
                    if val > buf[k - 1]:
                        buf[k - 1] = val
                        for j in range(k - 2, -1, -1):
                            if buf[j + 1] > buf[j]:
                                buf[j], buf[j + 1] = buf[j + 1], buf[j]
                            else:
                                break
                comparison_val[m] = buf[k - 1]

        current_stat = sorted_stats[step]
        count = 0
        for m in range(M):
            if comparison_val[m] >= current_stat:
                count += 1
        p_val = (count + 1) / (M + 1)

        if p_val < last_p:
            p_val = last_p
        last_p = p_val
        p_values[step] = p_val

        active[step] = False

        if p_val >= max_pval:
            for remaining in range(step + 1, n_valid):
                p_values[remaining] = max_pval
            break

    return p_values


@njit(cache=True, parallel=True)
def _stepdown_core_directional(sorted_stats, null_sub, dir_sub, obs_directions, k, max_pval):
    """Step-down with direction-conditional removal (matches David's Igor).

    In David's implementation, a sequence is only removed from a bootstrap row
    if the bootstrap went the SAME direction as the observed stat being processed.
    If the bootstrap went the other direction, that sequence's magnitude stays
    active for subsequent k-th largest computations.

    Optimization: caches per-row comparison values. After each step, only rows
    where the direction matched (and the active set changed) are recomputed.
    Rows where the removed column wasn't in their top-k are skipped entirely.
    """
    n_valid = len(sorted_stats)
    M = null_sub.shape[0]
    p_values = np.empty(n_valid)
    row_active = np.ones((M, n_valid), dtype=np.bool_)
    comparison_val = np.empty(M)
    last_p = 0.0

    for m in prange(M):
        buf = np.empty(k)
        buf[:] = -np.inf
        n_active_m = 0
        for col in range(n_valid):
            if not row_active[m, col]:
                continue
            n_active_m += 1
            val = null_sub[m, col]
            if val > buf[k - 1]:
                buf[k - 1] = val
                for j in range(k - 2, -1, -1):
                    if buf[j + 1] > buf[j]:
                        buf[j], buf[j + 1] = buf[j + 1], buf[j]
                    else:
                        break
        if n_active_m <= k - 1:
            comparison_val[m] = -np.inf
        else:
            comparison_val[m] = buf[k - 1]

    for step in range(n_valid):
        current_stat = sorted_stats[step]
        count = 0
        for m in range(M):
            if comparison_val[m] >= current_stat:
                count += 1
        p_val = (count + 1) / (M + 1)

        if p_val < last_p:
            p_val = last_p
        last_p = p_val
        p_values[step] = p_val

        if p_val >= max_pval:
            for remaining in range(step + 1, n_valid):
                p_values[remaining] = max_pval
            break

        obs_dir = obs_directions[step]
        for m in prange(M):
            if dir_sub[m, step] == obs_dir:
                row_active[m, step] = False
                if null_sub[m, step] >= comparison_val[m]:
                    buf = np.empty(k)
                    buf[:] = -np.inf
                    n_active_m = 0
                    for col in range(n_valid):
                        if not row_active[m, col]:
                            continue
                        n_active_m += 1
                        val = null_sub[m, col]
                        if val > buf[k - 1]:
                            buf[k - 1] = val
                            for j in range(k - 2, -1, -1):
                                if buf[j + 1] > buf[j]:
                                    buf[j], buf[j + 1] = buf[j + 1], buf[j]
                                else:
                                    break
                    if n_active_m <= k - 1:
                        comparison_val[m] = -np.inf
                    else:
                        comparison_val[m] = buf[k - 1]

    return p_values


def _prepare_null_sub(test_stats, null_matrix, null_directions=None):
    """Prepare sorted stats and null submatrix (shared across k-iterations).

    test_stats: shape (2S,) -- one-sided stats, one direction per sequence non-NaN
    null_matrix: shape (M, S) -- magnitude per sequence per resample
    null_directions: shape (M, S) -- direction per value (0=pos, 1=neg, -1=none)

    Maps each valid observed stat to its sequence's column in the null matrix.
    Returns direction info needed for direction-conditional removal.
    """
    valid_mask = ~np.isnan(test_stats)
    valid_indices = np.where(valid_mask)[0]
    valid_stats = test_stats[valid_indices]

    sort_order = np.argsort(-valid_stats)
    sorted_indices = valid_indices[sort_order]
    sorted_stats = valid_stats[sort_order]

    seq_indices = sorted_indices // 2
    null_sub = null_matrix[:, seq_indices].copy()
    null_sub[np.isnan(null_sub)] = -np.inf

    obs_directions = np.array(sorted_indices % 2, dtype=np.int8)

    dir_sub = None
    if null_directions is not None:
        dir_sub = null_directions[:, seq_indices].copy()

    return sorted_stats, sorted_indices, null_sub, obs_directions, dir_sub


def romano_wolf_stepdown(test_stats, null_matrix, null_directions=None, k=1):
    """Apply Romano-Wolf step-down procedure with k-FWER.

    Args:
        test_stats: shape (2S,) observed test statistics (NaN for unused directions)
        null_matrix: shape (M, S) bootstrap null magnitudes
        null_directions: shape (M, S) bootstrap directions (enables directional removal)
        k: the k for k-FWER (k-th largest value per row instead of max)

    Returns:
        p_values: shape (2S,) adjusted p-values (NaN where test stat is NaN)
    """
    sorted_stats, sorted_indices, null_sub, obs_directions, dir_sub = \
        _prepare_null_sub(test_stats, null_matrix, null_directions)

    if dir_sub is not None:
        step_p_values = _stepdown_core_directional(
            sorted_stats, null_sub, dir_sub, obs_directions, k, 1.0)
    else:
        step_p_values = _stepdown_core(sorted_stats, null_sub, k, 1.0)

    p_values = np.full_like(test_stats, np.nan)
    for i in range(len(sorted_indices)):
        p_values[sorted_indices[i]] = step_p_values[i]

    return p_values


def _count_rejections_directional(sorted_stats, null_sub, dir_sub, obs_directions, k, alpha):
    """Fast rejection count with directional removal."""
    step_p_values = _stepdown_core_directional(
        sorted_stats, null_sub, dir_sub, obs_directions, k, alpha)
    return int(np.sum(step_p_values < alpha))


def _count_rejections(sorted_stats, null_sub, k, alpha):
    """Fast rejection count -- stops as soon as p >= alpha."""
    step_p_values = _stepdown_core(sorted_stats, null_sub, k, alpha)
    return int(np.sum(step_p_values < alpha))


def find_k_fwer(test_stats, null_matrix, alpha=0.5, gamma=0.05, null_directions=None,
                return_history=False):
    """Compute adjusted p-values with FDP control via iterative k-FWER.

    Iterates: run step-down at current k -> count rejections -> update k ->
    repeat until convergence. Reports p-values from the CONVERGED k (full
    step-down, not early-stopped). This matches David's Igor implementation.

    When null_directions is provided, uses direction-conditional removal in the
    step-down (matching David's behavior: only remove a sequence from a bootstrap
    row if the bootstrap went the same direction as the observed stat).

    Convergence: stop when rejections < (k / gamma) - 1.

    Returns (g_values, k_final) where g_values are the adjusted p-values
    from the converged-k step-down, and k_final is the converged k.
    """
    sorted_stats, sorted_indices, null_sub, obs_directions, dir_sub = \
        _prepare_null_sub(test_stats, null_matrix, null_directions)

    k = 1
    k_history = []
    for _ in range(100):
        if dir_sub is not None:
            rejections = _count_rejections_directional(
                sorted_stats, null_sub, dir_sub, obs_directions, k, alpha)
        else:
            rejections = _count_rejections(sorted_stats, null_sub, k, alpha)
        k_history.append({"k": k, "rejections": int(rejections)})
        if rejections < (k / gamma) - 1:
            break
        new_k = int(np.ceil((rejections + 1) * gamma))
        k = k + 1 if new_k == k else new_k

    if dir_sub is not None:
        step_p_values = _stepdown_core_directional(
            sorted_stats, null_sub, dir_sub, obs_directions, k, 1.0)
    else:
        step_p_values = _stepdown_core(sorted_stats, null_sub, k, 1.0)

    p_values = np.full_like(test_stats, np.nan)
    for i in range(len(sorted_indices)):
        p_values[sorted_indices[i]] = step_p_values[i]

    if return_history:
        return p_values, k, k_history
    return p_values, k


def find_k_fwer_k1(test_stats, null_matrix, alpha=0.5, gamma=0.05, null_directions=None):
    """Conservative variant: always report k=1 p-values.

    This is the standard Romano-Wolf (max-based) step-down without the
    iterative k relaxation. More conservative than the full k-FWER procedure.
    """
    sorted_stats, sorted_indices, null_sub, obs_directions, dir_sub = \
        _prepare_null_sub(test_stats, null_matrix, null_directions)

    if dir_sub is not None:
        step_p_values = _stepdown_core_directional(
            sorted_stats, null_sub, dir_sub, obs_directions, 1, 1.0)
    else:
        step_p_values = _stepdown_core(sorted_stats, null_sub, 1, 1.0)

    p_values = np.full_like(test_stats, np.nan)
    for i in range(len(sorted_indices)):
        p_values[sorted_indices[i]] = step_p_values[i]

    rejections = int(np.sum(step_p_values < alpha))
    k = max(1, int(np.ceil((rejections + 1) * gamma)))

    return p_values, k


def find_k_fwer_chunked(test_stats, count_matrix, group_indices, params,
                        chunk_size=500, rng=None, return_history=False):
    """Memory-efficient CBAS: generates bootstrap directly into null_sub in chunks.

    Instead of allocating the full null matrix (M x 2S) and then extracting the
    valid/sorted subset, this generates bootstrap stats directly into the
    sorted null_sub array (M x n_valid) in row-chunks. This halves peak memory
    by avoiding the intermediate full-width matrix.

    Args:
        test_stats: shape (2S,) observed test statistics
        count_matrix: shape (N, S) count matrix
        group_indices: [grp0_indices, grp1_indices]
        params: CBASParams
        chunk_size: bootstrap rows generated per chunk (default 500)
        rng: numpy random Generator (default: seeded at 42)
        return_history: if True, return k_history as third element

    Returns (g_values, k_final) or (g_values, k_final, k_history).
    """
    if rng is None:
        rng = np.random.default_rng(2)

    grp0 = group_indices[0]
    grp1 = group_indices[1]
    n0 = len(grp0)
    n1 = len(grp1)
    n_total = n0 + n1
    n_seq = count_matrix.shape[1]
    M = params.resample_number

    if params.centering:
        obs_delta = count_matrix[grp0].mean(axis=0) - count_matrix[grp1].mean(axis=0)
    else:
        obs_delta = np.zeros(n_seq)

    boot_indices_0 = rng.integers(0, n_total, size=(M, n0))
    boot_indices_1 = rng.integers(0, n_total, size=(M, n1))

    valid_mask = ~np.isnan(test_stats)
    valid_indices = np.where(valid_mask)[0]
    valid_stats = test_stats[valid_indices]
    sort_order = np.argsort(-valid_stats)
    sorted_indices = valid_indices[sort_order]
    sorted_stats = valid_stats[sort_order].astype(np.float64)
    n_valid = len(sorted_stats)

    sorted_col_indices = np.ascontiguousarray(sorted_indices, dtype=np.int64)
    count_matrix_f = np.ascontiguousarray(count_matrix, dtype=np.float64)
    obs_delta_f = np.ascontiguousarray(obs_delta, dtype=np.float64)

    obs_directions = np.array(sorted_indices % 2, dtype=np.int8)

    null_sub = np.full((M, n_valid), -np.inf, dtype=np.float64)
    dir_sub = np.full((M, n_valid), np.int8(-1), dtype=np.int8)

    for chunk_start in range(0, M, chunk_size):
        chunk_end = min(chunk_start + chunk_size, M)
        _bootstrap_chunk_into(
            count_matrix_f, boot_indices_0, boot_indices_1,
            n0, n1, obs_delta_f, sorted_col_indices,
            chunk_start, chunk_end, null_sub, dir_sub
        )

    k = 1
    alpha = params.alpha
    gamma = params.gamma
    k_history = []
    for _ in range(100):
        step_p_values = _stepdown_core_directional(
            sorted_stats, null_sub, dir_sub, obs_directions, k, alpha)
        rejections = int(np.sum(step_p_values < alpha))
        k_history.append({"k": k, "rejections": rejections})
        if rejections < (k / gamma) - 1:
            break
        new_k = int(np.ceil((rejections + 1) * gamma))
        k = k + 1 if new_k == k else new_k

    step_p_values = _stepdown_core_directional(
        sorted_stats, null_sub, dir_sub, obs_directions, k, 1.0)

    p_values = np.full_like(test_stats, np.nan)
    for i in range(n_valid):
        p_values[sorted_indices[i]] = step_p_values[i]

    if return_history:
        return p_values, k, k_history
    return p_values, k
