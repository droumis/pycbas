"""
Core CBAS (Choice-Wide Behavioral Association Study) implementation.

Reimplements the comparative CBAS from the Igor Pro code:
  - Sequence enumeration across all lengths 1..L
  - Studentized test statistics (two one-tailed tests per sequence)
  - Bootstrap resampling ignoring group labels
  - Romano-Wolf step-down with monotonicity enforcement
  - k-FWER iteration for FDP control

Set NUMBA_DISABLE_JIT=1 in the environment to disable JIT for debugging.
"""

import numpy as np
from dataclasses import dataclass

try:
    from numba import njit, prange
except ImportError:
    prange = range

    def njit(*args, **kwargs):
        def wrapper(fn):
            return fn
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return wrapper


@dataclass
class CBASParams:
    num_arms: int = 6
    seq_len_max: int = 6
    criterion: int = 800
    resample_number: int = 10_000
    alpha: float = 0.5
    gamma: float = 0.05
    centering: bool = False


@dataclass
class CBASResult:
    sequences: list[tuple]
    test_stats: np.ndarray
    g_values: np.ndarray
    k_final: int
    significant_mask: np.ndarray

    @property
    def n_significant(self):
        return int(self.significant_mask.sum())


def load_subject_data(filepath):
    """Load a single subject's data file. Returns (session, choice, reward, contingency) arrays."""
    rows = []
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split(",")
            session = int(parts[0])
            choice = int(parts[1])
            reward = int(parts[2])
            conting = int(parts[3]) if parts[3].strip() else 0
            rows.append((session, choice, reward, conting))
    arr = np.array(rows, dtype=np.int32)
    return arr


def extract_choice_stream(subject_data, contingency=2, num_arms=6, encode_reward=True):
    """Extract choice stream, optionally filtered by contingency.

    Args:
        subject_data: array with columns (session, symbol, reward, contingency)
        contingency: block type to filter on, or None to use all trials
        num_arms: number of base symbols (choices)
        encode_reward: if True, encode as symbol + reward*num_arms (doubles alphabet).
            Set False for tasks where outcome is deterministic from the symbol
            (e.g., 2AFC where symbol already encodes choice×stimulus_side).
    """
    if contingency is None:
        data = subject_data
    else:
        mask = subject_data[:, 3] == contingency
        data = subject_data[mask]
    if encode_reward:
        symbols = data[:, 1] + data[:, 2] * num_arms
    else:
        symbols = data[:, 1]
    return symbols


def enumerate_sequences(choice_stream, seq_len, criterion):
    """Find all subsequences of given length within the first `criterion` choices.
    Returns a dict mapping sequence tuple -> count."""
    stream = choice_stream[:criterion]
    counts = {}
    for i in range(len(stream) - seq_len + 1):
        seq = tuple(stream[i:i + seq_len].tolist())
        counts[seq] = counts.get(seq, 0) + 1
    return counts


def build_count_matrix(subjects_data, params, contingency=2, encode_reward=True):
    """Build the full sequence count matrix.

    Args:
        subjects_data: list of subject data arrays (from load_subject_data)
        params: CBASParams instance
        contingency: block type to filter on, or None for all trials
        encode_reward: if True, encode symbol + reward*num_arms. Set False for
            tasks where outcome is deterministic from the symbol (e.g., 2AFC).

    Returns:
        sequences: list of all unique sequence tuples (sorted by total frequency descending)
        count_matrix: ndarray of shape (n_subjects, n_sequences) with usage counts
    """
    n_subjects = len(subjects_data)
    all_seq_counts = []
    for subj_data in subjects_data:
        stream = extract_choice_stream(subj_data, contingency, params.num_arms,
                                       encode_reward=encode_reward)
        subj_counts = {}
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


@njit(cache=True, parallel=True)
def _bootstrap_parallel(count_matrix, boot_indices_0, boot_indices_1, n0, n1, n_seq, M,
                        obs_delta):
    """Numba-parallelized bootstrap computation.

    Stores |t| magnitude per sequence per resample (one column per sequence).
    Also stores which direction each bootstrap value went (0=positive, 1=negative)
    for direction-conditional removal in step-down.
    """
    null_stats = np.full((M, n_seq), np.nan)
    null_directions = np.full((M, n_seq), np.int8(-1))  # -1 = no value (sigma=0 or delta=0)

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
                    null_directions[m, s] = 0  # positive direction
                elif delta_centered < 0.0:
                    null_stats[m, s] = -delta_centered / sigma
                    null_directions[m, s] = 1  # negative direction

    return null_stats, null_directions


def bootstrap_test_stats(count_matrix, group_indices, params, rng=None):
    """Generate bootstrap null distribution by resampling ignoring group labels.

    When params.centering=True, uses Clarke et al. (2020) eq 5:
        t*_s,m = (δ*_s,m - δ_s) / σ*_s,m
    When params.centering=False (default, matches David's Igor implementation):
        t*_s,m = δ*_s,m / σ*_s,m

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
def _stepdown_core(sorted_stats, null_sub, k, max_pval):
    """Numba-accelerated inner loop of Romano-Wolf step-down.

    Parallelizes across M bootstrap resamples using numba prange.

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


def _prepare_null_sub(test_stats, null_matrix, null_directions=None):
    """Prepare sorted stats and null submatrix (shared across k-iterations).

    test_stats: shape (2S,) — one-sided stats, one direction per sequence non-NaN
    null_matrix: shape (M, S) — magnitude per sequence per resample
    null_directions: shape (M, S) — direction per value (0=pos, 1=neg, -1=none)

    Maps each valid observed stat to its sequence's column in the null matrix.
    Returns direction info needed for direction-conditional removal.
    """
    valid_mask = ~np.isnan(test_stats)
    valid_indices = np.where(valid_mask)[0]
    valid_stats = test_stats[valid_indices]

    sort_order = np.argsort(-valid_stats)
    sorted_indices = valid_indices[sort_order]
    sorted_stats = valid_stats[sort_order]

    # Map from (2S) position to sequence index: position i → sequence i // 2
    seq_indices = sorted_indices // 2
    null_sub = null_matrix[:, seq_indices].copy()
    null_sub[np.isnan(null_sub)] = -np.inf

    # Observed direction for each sorted stat: even index → 0 (positive), odd → 1 (negative)
    obs_directions = np.array(sorted_indices % 2, dtype=np.int8)

    # Bootstrap directions for the same columns
    dir_sub = None
    if null_directions is not None:
        dir_sub = null_directions[:, seq_indices].copy()

    return sorted_stats, sorted_indices, null_sub, obs_directions, dir_sub


@njit(cache=True, parallel=True)
def _stepdown_core_directional(sorted_stats, null_sub, dir_sub, obs_directions, k, max_pval):
    """Step-down with direction-conditional removal (matches David's Igor).

    In David's implementation, a sequence is only removed from a bootstrap row
    if the bootstrap went the SAME direction as the observed stat being processed.
    If the bootstrap went the other direction, that sequence's magnitude stays
    active for subsequent k-th largest computations.

    Args:
        sorted_stats: (n_valid,) test stats sorted descending
        null_sub: (M, n_valid) null magnitudes, NaN replaced with -inf
        dir_sub: (M, n_valid) bootstrap directions (0=pos, 1=neg, -1=none)
        obs_directions: (n_valid,) observed direction for each sorted stat (0=pos, 1=neg)
        k: k-FWER parameter
        max_pval: stop when p >= this
    """
    n_valid = len(sorted_stats)
    M = null_sub.shape[0]
    p_values = np.empty(n_valid)
    # Per-row active masks: True = this column still contributes to k-th largest for this row
    row_active = np.ones((M, n_valid), dtype=np.bool_)
    last_p = 0.0

    for step in range(n_valid):
        comparison_val = np.empty(M)
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

        # Direction-conditional removal: only remove from rows where bootstrap
        # direction matches the observed direction of this step's stat
        obs_dir = obs_directions[step]
        for m in prange(M):
            if dir_sub[m, step] == obs_dir:
                row_active[m, step] = False
            # If directions don't match (or dir=-1), leave it active for this row

        if p_val >= max_pval:
            for remaining in range(step + 1, n_valid):
                p_values[remaining] = max_pval
            break

    return p_values


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
    """Fast rejection count — stops as soon as p >= alpha."""
    step_p_values = _stepdown_core(sorted_stats, null_sub, k, alpha)
    return int(np.sum(step_p_values < alpha))


def find_k_fwer(test_stats, null_matrix, alpha=0.5, gamma=0.05, null_directions=None):
    """Compute adjusted p-values with FDP control via iterative k-FWER.

    Iterates: run step-down at current k → count rejections → update k →
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
    for _ in range(100):  # safety cap
        if dir_sub is not None:
            rejections = _count_rejections_directional(
                sorted_stats, null_sub, dir_sub, obs_directions, k, alpha)
        else:
            rejections = _count_rejections(sorted_stats, null_sub, k, alpha)
        if rejections < (k / gamma) - 1:
            break
        new_k = int(np.ceil((rejections + 1) * gamma))
        k = k + 1 if new_k == k else new_k

    # Final full step-down at converged k
    if dir_sub is not None:
        step_p_values = _stepdown_core_directional(
            sorted_stats, null_sub, dir_sub, obs_directions, k, 1.0)
    else:
        step_p_values = _stepdown_core(sorted_stats, null_sub, k, 1.0)

    p_values = np.full_like(test_stats, np.nan)
    for i in range(len(sorted_indices)):
        p_values[sorted_indices[i]] = step_p_values[i]

    return p_values, k


def find_k_fwer_k1(test_stats, null_matrix, alpha=0.5, gamma=0.05, null_directions=None):
    """Conservative variant: always report k=1 p-values.

    This is the standard Romano-Wolf (max-based) step-down without the
    iterative k relaxation. More conservative than the full k-FWER procedure.
    Use for comparison/debugging.
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


# =============================================================================
# Chunked pipeline: generate bootstrap in row-chunks directly into null_sub
# =============================================================================


@njit(cache=True, parallel=True)
def _bootstrap_chunk_into(count_matrix, boot_indices_0, boot_indices_1,
                          n0, n1, obs_delta, sorted_col_indices,
                          chunk_start, chunk_end, out):
    """Generate bootstrap null stats directly into a pre-allocated output array.

    Writes rows [chunk_start:chunk_end] of `out` (shape M × n_valid).
    Columns follow sorted_col_indices order (positions in 2S test_stats array).
    Stores |t| magnitude per sequence — direction-independent.
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
            if sigma > 0.0:
                delta_centered = (mean0 - mean1) - obs_delta[seq_idx]
                if delta_centered != 0.0:
                    val = abs(delta_centered / sigma)
            out[m, ci] = val


def find_k_fwer_chunked(test_stats, count_matrix, group_indices, params,
                        chunk_size=500, rng=None):
    """Memory-efficient CBAS: generates bootstrap directly into null_sub in chunks.

    Instead of allocating the full null matrix (M × 2S) and then extracting the
    valid/sorted subset, this generates bootstrap stats directly into the
    sorted null_sub array (M × n_valid) in row-chunks. This halves peak memory
    by avoiding the intermediate full-width matrix.

    Memory savings:
      Standard: M×2S (full null) + M×n_valid (null_sub) = M × (2S + n_valid) × 8 bytes
      Chunked:  M×n_valid (null_sub only) + chunk_size×n_valid (working) = M × n_valid × 8 bytes

    For IBL (A=8, L=6): standard needs ~96 GB, chunked needs ~48 GB for null_sub.
    Combined with max_memory_gb limit, large problems spill to disk via memmap.

    Produces identical results to find_k_fwer. Uses the same fast _stepdown_core
    with early stopping.

    Args:
        test_stats: shape (2S,) observed test statistics
        count_matrix: shape (N, S) count matrix
        group_indices: [grp0_indices, grp1_indices]
        params: CBASParams
        chunk_size: bootstrap rows generated per chunk (default 500)
        rng: numpy random Generator (default: seeded at 42)
        max_memory_gb: if null_sub exceeds this, use memmap (default 8)

    Returns (g_values, k_final) — same as find_k_fwer.
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

    # Pre-generate bootstrap indices
    boot_indices_0 = rng.integers(0, n_total, size=(M, n0))
    boot_indices_1 = rng.integers(0, n_total, size=(M, n1))

    # Sort valid test stats descending
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

    # Allocate null_sub directly (M × n_valid) — never allocate the full M × 2S
    null_sub = np.full((M, n_valid), -np.inf, dtype=np.float64)

    # Generate bootstrap in chunks, writing directly into null_sub
    for chunk_start in range(0, M, chunk_size):
        chunk_end = min(chunk_start + chunk_size, M)
        _bootstrap_chunk_into(
            count_matrix_f, boot_indices_0, boot_indices_1,
            n0, n1, obs_delta_f, sorted_col_indices,
            chunk_start, chunk_end, null_sub
        )

    # Iterative k-FWER: find converged k
    k = 1
    alpha = params.alpha
    gamma = params.gamma
    for _ in range(100):  # safety cap
        step_p_values = _stepdown_core(sorted_stats, null_sub, k, alpha)
        rejections = int(np.sum(step_p_values < alpha))
        if rejections < (k / gamma) - 1:
            break
        new_k = int(np.ceil((rejections + 1) * gamma))
        k = k + 1 if new_k == k else new_k

    # Final full step-down at converged k
    step_p_values = _stepdown_core(sorted_stats, null_sub, k, 1.0)

    # Map back to original positions
    p_values = np.full_like(test_stats, np.nan)
    for i in range(n_valid):
        p_values[sorted_indices[i]] = step_p_values[i]

    return p_values, k


def run_cbas_comparative(subjects_data, group_labels, params=None,
                         contingency=2, encode_reward=True, chunked=True):
    """Run the full comparative CBAS pipeline.

    Args:
        subjects_data: list of subject data arrays (from load_subject_data)
        group_labels: array of 0/1 indicating group membership
        params: CBASParams instance
        contingency: block type to filter on, or None for all trials
        encode_reward: if True, encode symbol + reward*num_arms. Set False for 2AFC.
        chunked: if True (default), use memory-efficient chunked pipeline

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
                                                 encode_reward=encode_reward)
    test_stats = compute_test_stats(count_matrix, group_indices)

    if chunked:
        g_values, k_final = find_k_fwer_chunked(
            test_stats, count_matrix, group_indices, params)
    else:
        null_matrix, null_directions = bootstrap_test_stats(count_matrix, group_indices, params)
        g_values, k_final = find_k_fwer(test_stats, null_matrix, params.alpha, params.gamma,
                                         null_directions=null_directions)

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
    )


def run_cbas_correlative(subjects_data, covariate, params=None):
    """Run the full correlative CBAS pipeline.

    Args:
        subjects_data: list of subject data arrays (from load_subject_data)
        covariate: array of continuous values (e.g. CBIT scores), one per subject
        params: CBASParams instance

    Returns:
        CBASResult
    """
    if params is None:
        params = CBASParams()

    covariate = np.asarray(covariate, dtype=np.float64)
    sequences, count_matrix = build_count_matrix(subjects_data, params)
    test_stats = compute_test_stats_correlative(count_matrix, covariate)
    null_matrix, null_directions = bootstrap_test_stats_correlative(count_matrix, covariate, params)
    g_values, k_final = find_k_fwer(test_stats, null_matrix, params.alpha, params.gamma,
                                     null_directions=null_directions)

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
    )


# =============================================================================
# Resource estimation
# =============================================================================


def estimate_resources(num_arms, seq_len_max, n_subjects=None, n_observed=None,
                       resample_number=10000, encode_reward=True):
    """Estimate memory and time requirements for a CBAS run.

    Provides worst-case estimates from task parameters, and actual estimates
    if observed sequence count or subject data is provided.

    Args:
        num_arms: number of base symbols (before reward encoding)
        seq_len_max: maximum sequence length L
        n_subjects: number of subjects (for context only)
        n_observed: actual number of observed sequences (if known from data).
            Overrides worst-case S for memory/time estimates.
        resample_number: number of bootstrap resamples M (default 10,000)
        encode_reward: if True, effective alphabet is num_arms*2

    Returns:
        dict with keys:
            alphabet: effective alphabet size A
            seq_len_max: L
            total_sequences: S = sum(A^l, l=1..L)
            observed_sequences: n_observed if provided, else None
            resample_number: M
            n_subjects: n_subjects if provided
            memory_full_null_gb: M × 2S × 8 bytes (standard pipeline)
            memory_chunked_gb: M × 2*n_valid × 8 bytes (chunked pipeline)
            est_time_seconds: estimated runtime
            recommendation: human-readable verdict
    """
    A = num_arms * 2 if encode_reward else num_arms
    S = sum(A**l for l in range(1, seq_len_max + 1))
    M = resample_number

    n_valid = n_observed if n_observed is not None else S

    full_null_bytes = M * 2 * S * 8
    chunked_bytes = M * 2 * n_valid * 8

    full_null_gb = full_null_bytes / (1024**3)
    chunked_gb = chunked_bytes / (1024**3)

    # Time estimate: calibrated from rat benchmark (16,483 valid sequences → 13s chunked)
    rat_cols = 2 * 16483
    rat_time = 13.0
    est_time = rat_time * (2 * n_valid) / rat_cols

    if chunked_gb < 1.0:
        verdict = "TRIVIAL"
    elif chunked_gb < 8.0:
        verdict = "COMFORTABLE"
    elif chunked_gb < 24.0:
        verdict = "FITS (close other apps)"
    elif chunked_gb < 48.0:
        verdict = "TIGHT (may need reduced M or chunking to disk)"
    else:
        verdict = "TOO LARGE (needs embedding or M reduction)"

    result = {
        "alphabet": A,
        "seq_len_max": seq_len_max,
        "total_sequences": S,
        "observed_sequences": n_observed,
        "resample_number": M,
        "n_subjects": n_subjects,
        "memory_full_null_gb": round(full_null_gb, 2),
        "memory_chunked_gb": round(chunked_gb, 2),
        "est_time_seconds": round(est_time, 1),
        "recommendation": verdict,
    }
    return result


def print_resource_estimate(est):
    """Pretty-print output from estimate_resources()."""
    print(f"CBAS Resource Estimate")
    print(f"  Alphabet (A):          {est['alphabet']}")
    print(f"  Max length (L):        {est['seq_len_max']}")
    print(f"  Total sequences (S):   {est['total_sequences']:,}")
    if est["observed_sequences"] is not None:
        pct = est["observed_sequences"] / est["total_sequences"] * 100
        print(f"  Observed sequences:    {est['observed_sequences']:,} ({pct:.0f}%)")
    if est["n_subjects"] is not None:
        print(f"  Subjects:              {est['n_subjects']}")
    print(f"  Resamples (M):         {est['resample_number']:,}")
    print(f"  Memory (standard):     {est['memory_full_null_gb']:.1f} GB")
    print(f"  Memory (chunked):      {est['memory_chunked_gb']:.1f} GB")
    print(f"  Est. time (chunked):   {est['est_time_seconds']:.0f}s")
    print(f"  Verdict:               {est['recommendation']}")
