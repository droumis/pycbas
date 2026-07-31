"""Bootstrap resampling for null distribution generation."""

import numpy as np
from ._numba import njit, prange


@njit(cache=True, parallel=True)
def _bootstrap_parallel(count_matrix, boot_indices_0, boot_indices_1, n0, n1, n_seq, M,
                        obs_delta):
    """Numba-parallelized bootstrap computation.

    Stores |t| magnitude per sequence per resample (one column per sequence).
    Also stores which direction each bootstrap value went (0=positive, 1=negative)
    for direction-conditional removal in step-down.
    """
    null_stats = np.full((M, n_seq), np.nan)
    null_directions = np.full((M, n_seq), np.int8(-1))

    for m in prange(M):
        for s in range(n_seq):
            sum0 = 0.0
            sum1 = 0.0
            for i in range(n0):
                sum0 += count_matrix[boot_indices_0[m, i], s]
            for i in range(n1):
                sum1 += count_matrix[boot_indices_1[m, i], s]
            mean0 = sum0 / n0
            mean1 = sum1 / n1

            var0 = 0.0
            var1 = 0.0
            for i in range(n0):
                diff = count_matrix[boot_indices_0[m, i], s] - mean0
                var0 += diff * diff
            for i in range(n1):
                diff = count_matrix[boot_indices_1[m, i], s] - mean1
                var1 += diff * diff

            sem0 = np.sqrt(var0 / (n0 * (n0 - 1)))
            sem1 = np.sqrt(var1 / (n1 * (n1 - 1)))
            sigma = np.sqrt(sem0 * sem0 + sem1 * sem1)

            if sigma > 0.0:
                delta_centered = (mean0 - mean1) - obs_delta[s]
                if delta_centered > 0.0:
                    null_stats[m, s] = delta_centered / sigma
                    null_directions[m, s] = 0
                elif delta_centered < 0.0:
                    null_stats[m, s] = -delta_centered / sigma
                    null_directions[m, s] = 1

    return null_stats, null_directions


def bootstrap_test_stats(count_matrix, group_indices, params, rng=None):
    """Generate bootstrap null distribution by resampling ignoring group labels.

    When params.centering=True, uses Clarke et al. (2020) eq 5:
        t*_s,m = (d*_s,m - d_s) / s*_s,m
    When params.centering=False (default, matches David's Igor implementation):
        t*_s,m = d*_s,m / s*_s,m

    Returns:
        null_stats: (M, S) magnitude per sequence per resample
        null_directions: (M, S) direction per value (0=pos, 1=neg, -1=no value)
    """
    if rng is None:
        rng = np.random.default_rng(42)

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

    count_matrix_f = np.ascontiguousarray(count_matrix, dtype=np.float64)
    obs_delta_f = np.ascontiguousarray(obs_delta, dtype=np.float64)
    null_stats, null_directions = _bootstrap_parallel(
        count_matrix_f, boot_indices_0, boot_indices_1, n0, n1, n_seq, M,
        obs_delta_f
    )

    return null_stats, null_directions


@njit(cache=True, parallel=True)
def _bootstrap_correlative_parallel(count_matrix, perm_indices, covariate, n, n_seq, M):
    """Numba-parallelized correlative bootstrap (permute covariate).

    Stores |t| magnitude per sequence per resample and direction.
    """
    null_stats = np.full((M, n_seq), np.nan)
    null_directions = np.full((M, n_seq), np.int8(-1))

    for m in prange(M):
        Y_bar = 0.0
        for i in range(n):
            Y_bar += covariate[perm_indices[m, i]]
        Y_bar /= n

        ss_Y = 0.0
        for i in range(n):
            d = covariate[perm_indices[m, i]] - Y_bar
            ss_Y += d * d

        if ss_Y == 0.0:
            continue

        for s in range(n_seq):
            X_bar = 0.0
            for i in range(n):
                X_bar += count_matrix[i, s]
            X_bar /= n

            ss_X = 0.0
            for i in range(n):
                d = count_matrix[i, s] - X_bar
                ss_X += d * d

            if ss_X == 0.0:
                continue

            sum_XY = 0.0
            tau_num_sq = 0.0
            for i in range(n):
                x_dev = count_matrix[i, s] - X_bar
                y_dev = covariate[perm_indices[m, i]] - Y_bar
                sum_XY += count_matrix[i, s] * covariate[perm_indices[m, i]]
                tau_num_sq += (x_dev * y_dev) ** 2

            rho = (sum_XY - n * X_bar * Y_bar) / np.sqrt(ss_X * ss_Y)

            tau_num = np.sqrt(tau_num_sq / n)
            tau_den = np.sqrt(ss_X / n) * np.sqrt(ss_Y / n)
            tau = tau_num / tau_den

            if tau == 0.0:
                continue

            t_val = np.sqrt(n) * rho / tau
            null_stats[m, s] = abs(t_val)
            if t_val > 0.0:
                null_directions[m, s] = 0
            elif t_val < 0.0:
                null_directions[m, s] = 1

    return null_stats, null_directions


def bootstrap_test_stats_correlative(count_matrix, covariate, params, rng=None):
    """Generate bootstrap null for correlative mode by permuting covariate.

    Returns:
        null_stats: (M, S) magnitude per sequence per resample
        null_directions: (M, S) direction per value (0=pos, 1=neg, -1=none)
    """
    if rng is None:
        rng = np.random.default_rng(42)

    n = count_matrix.shape[0]
    n_seq = count_matrix.shape[1]
    M = params.resample_number

    perm_indices = np.empty((M, n), dtype=np.int64)
    for m in range(M):
        perm_indices[m] = rng.permutation(n)

    count_matrix_f = np.ascontiguousarray(count_matrix, dtype=np.float64)
    covariate_f = np.ascontiguousarray(covariate, dtype=np.float64)

    null_stats, null_directions = _bootstrap_correlative_parallel(
        count_matrix_f, perm_indices, covariate_f, n, n_seq, M
    )

    return null_stats, null_directions


@njit(cache=True, parallel=True)
def _bootstrap_chunk_into(count_matrix, boot_indices_0, boot_indices_1,
                          n0, n1, obs_delta, sorted_col_indices,
                          chunk_start, chunk_end, out, out_dir):
    """Generate bootstrap null stats directly into a pre-allocated output array.

    Writes rows [chunk_start:chunk_end] of `out` (shape M x n_valid).
    Columns follow sorted_col_indices order (positions in 2S test_stats array).
    Stores |t| magnitude and direction (0=pos, 1=neg, -1=none).
    """
    C = chunk_end - chunk_start
    n_valid = len(sorted_col_indices)

    for cm in prange(C):
        m = chunk_start + cm
        for ci in range(n_valid):
            col_2s = sorted_col_indices[ci]
            seq_idx = col_2s // 2

            sum0 = 0.0
            sum1 = 0.0
            for i in range(n0):
                sum0 += count_matrix[boot_indices_0[m, i], seq_idx]
            for i in range(n1):
                sum1 += count_matrix[boot_indices_1[m, i], seq_idx]
            mean0 = sum0 / n0
            mean1 = sum1 / n1

            var0 = 0.0
            var1 = 0.0
            for i in range(n0):
                diff = count_matrix[boot_indices_0[m, i], seq_idx] - mean0
                var0 += diff * diff
            for i in range(n1):
                diff = count_matrix[boot_indices_1[m, i], seq_idx] - mean1
                var1 += diff * diff

            sem0 = np.sqrt(var0 / (n0 * (n0 - 1)))
            sem1 = np.sqrt(var1 / (n1 * (n1 - 1)))
            sigma = np.sqrt(sem0 * sem0 + sem1 * sem1)

            val = -np.inf
            d = np.int8(-1)
            if sigma > 0.0:
                delta_centered = (mean0 - mean1) - obs_delta[seq_idx]
                if delta_centered > 0.0:
                    val = delta_centered / sigma
                    d = np.int8(0)
                elif delta_centered < 0.0:
                    val = -delta_centered / sigma
                    d = np.int8(1)
            out[m, ci] = val
            out_dir[m, ci] = d
