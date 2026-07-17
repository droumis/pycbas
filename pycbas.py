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


def extract_choice_stream(subject_data, contingency=2, num_arms=6):
    """Extract choice stream for a given contingency. Encodes as choice + reward*num_arms."""
    mask = subject_data[:, 3] == contingency
    data = subject_data[mask]
    symbols = data[:, 1] + data[:, 2] * num_arms
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


def build_count_matrix(subjects_data, params, contingency=2):
    """Build the full sequence count matrix.

    Returns:
        sequences: list of all unique sequence tuples (sorted by total frequency descending)
        count_matrix: ndarray of shape (n_subjects, n_sequences) with usage counts
    """
    n_subjects = len(subjects_data)
    all_seq_counts = []
    for subj_data in subjects_data:
        stream = extract_choice_stream(subj_data, contingency, params.num_arms)
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
def _bootstrap_parallel(count_matrix, boot_indices_0, boot_indices_1, n0, n1, n_seq, M):
    """Numba-parallelized bootstrap computation."""
    null_stats = np.full((M, n_seq * 2), np.nan)

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

            delta = mean0 - mean1
            if sigma > 0.0 and delta != 0.0:
                t_val = delta / sigma
                if delta > 0.0:
                    null_stats[m, s * 2] = t_val
                else:
                    null_stats[m, s * 2 + 1] = -t_val

    return null_stats


def bootstrap_test_stats(count_matrix, group_indices, params, rng=None):
    """Generate bootstrap null distribution by resampling ignoring group labels.

    Returns array of shape (resample_number, n_sequences * 2).
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

    # Pre-generate all bootstrap indices
    boot_indices_0 = rng.integers(0, n_total, size=(M, n0))
    boot_indices_1 = rng.integers(0, n_total, size=(M, n1))

    count_matrix_f = np.ascontiguousarray(count_matrix, dtype=np.float64)
    null_stats = _bootstrap_parallel(
        count_matrix_f, boot_indices_0, boot_indices_1, n0, n1, n_seq, M
    )

    return null_stats


@njit(cache=True, parallel=True)
def _bootstrap_correlative_parallel(count_matrix, perm_indices, covariate, n, n_seq, M):
    """Numba-parallelized correlative bootstrap (permute covariate)."""
    null_stats = np.full((M, n_seq * 2), np.nan)

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

            if rho > 0.0:
                null_stats[m, s * 2] = t_val
            elif rho < 0.0:
                null_stats[m, s * 2 + 1] = -t_val

    return null_stats


def bootstrap_test_stats_correlative(count_matrix, covariate, params, rng=None):
    """Generate bootstrap null for correlative mode by permuting covariate.

    Returns array of shape (resample_number, n_sequences * 2).
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

    null_stats = _bootstrap_correlative_parallel(
        count_matrix_f, perm_indices, covariate_f, n, n_seq, M
    )

    return null_stats


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


def _prepare_null_sub(test_stats, null_matrix):
    """Prepare sorted stats and null submatrix (shared across k-iterations)."""
    valid_mask = ~np.isnan(test_stats)
    valid_indices = np.where(valid_mask)[0]
    valid_stats = test_stats[valid_indices]

    sort_order = np.argsort(-valid_stats)
    sorted_indices = valid_indices[sort_order]
    sorted_stats = valid_stats[sort_order]

    null_sub = null_matrix[:, sorted_indices].copy()
    null_sub[np.isnan(null_sub)] = -np.inf

    return sorted_stats, sorted_indices, null_sub


def romano_wolf_stepdown(test_stats, null_matrix, k=1):
    """Apply Romano-Wolf step-down procedure with k-FWER.

    Args:
        test_stats: shape (2S,) observed test statistics (NaN for unused directions)
        null_matrix: shape (M, 2S) bootstrap null statistics
        k: the k for k-FWER (k-th largest value per row instead of max)

    Returns:
        p_values: shape (2S,) adjusted p-values (NaN where test stat is NaN)
    """
    sorted_stats, sorted_indices, null_sub = _prepare_null_sub(test_stats, null_matrix)
    step_p_values = _stepdown_core(sorted_stats, null_sub, k, 1.0)

    p_values = np.full_like(test_stats, np.nan)
    for i in range(len(sorted_indices)):
        p_values[sorted_indices[i]] = step_p_values[i]

    return p_values


def _count_rejections(sorted_stats, null_sub, k, alpha):
    """Fast rejection count — stops as soon as p >= alpha."""
    step_p_values = _stepdown_core(sorted_stats, null_sub, k, alpha)
    return int(np.sum(step_p_values < alpha))


def find_k_fwer(test_stats, null_matrix, alpha=0.5, gamma=0.05):
    """Iterate to find k for FDP control.

    Uses early-exit step-down during k-iteration (only counts rejections),
    then does one final full pass at the converged k.

    Returns (g_values, k_final) where g_values are the adjusted p-values
    under the final k, and k_final is the converged k.
    """
    sorted_stats, sorted_indices, null_sub = _prepare_null_sub(test_stats, null_matrix)

    k = 1
    prev_rejections = None

    while True:
        rejections = _count_rejections(sorted_stats, null_sub, k, alpha)

        if prev_rejections is not None and rejections < prev_rejections:
            break

        required_k = int(np.ceil((rejections + 1) * gamma))
        if required_k <= k:
            break

        prev_rejections = rejections
        k = required_k

    # Final full pass at converged k
    step_p_values = _stepdown_core(sorted_stats, null_sub, k, 1.0)

    p_values = np.full_like(test_stats, np.nan)
    for i in range(len(sorted_indices)):
        p_values[sorted_indices[i]] = step_p_values[i]

    return p_values, k


def run_cbas_comparative(subjects_data, group_labels, params=None):
    """Run the full comparative CBAS pipeline.

    Args:
        subjects_data: list of subject data arrays (from load_subject_data)
        group_labels: array of 0/1 indicating group membership
        params: CBASParams instance

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

    sequences, count_matrix = build_count_matrix(subjects_data, params)
    test_stats = compute_test_stats(count_matrix, group_indices)
    null_matrix = bootstrap_test_stats(count_matrix, group_indices, params)
    g_values, k_final = find_k_fwer(test_stats, null_matrix, params.alpha, params.gamma)

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
    null_matrix = bootstrap_test_stats_correlative(count_matrix, covariate, params)
    g_values, k_final = find_k_fwer(test_stats, null_matrix, params.alpha, params.gamma)

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
