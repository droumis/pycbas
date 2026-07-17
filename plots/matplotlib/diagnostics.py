import numpy as np
import matplotlib.pyplot as plt


def null_vs_observed(test_stats, null_row_maxes, ax=None):
    """Overlaid histograms of observed test statistics and null row maxima.

    Parameters
    ----------
    test_stats : array
        Observed test statistics (may contain NaNs).
    null_row_maxes : array
        Null distribution row maxima (may contain NaNs).
    ax : matplotlib Axes or None
        If None, creates a new figure.

    Returns
    -------
    fig, ax
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 3.5))
    else:
        fig = ax.get_figure()

    valid_observed = test_stats[~np.isnan(test_stats)]
    valid_null = null_row_maxes[~np.isnan(null_row_maxes)]

    if len(valid_observed) > 0:
        ax.hist(valid_observed, bins=80, density=True, alpha=0.6, color="steelblue",
                label="Observed test statistics", edgecolor="white", linewidth=0.3)
        obs_max = np.nanmax(valid_observed)
        ax.axvline(obs_max, color="red", linewidth=2,
                   label=f"Observed max = {obs_max:.2f}")

    if len(valid_null) > 0:
        ax.hist(valid_null, bins=50, density=True, alpha=0.6, color="gray",
                label="Null row-max (per resample)", edgecolor="white", linewidth=0.3)

    ax.set_xlabel("Test statistic")
    ax.set_ylabel("Density")
    ax.set_title("Null Distribution vs Observed")
    ax.legend(fontsize=8)

    fig.tight_layout()
    return fig, ax


def gvalue_distribution(g_values, n_seq, alpha=0.5, ax=None):
    """Histogram of best-direction g-values with significance threshold.

    Parameters
    ----------
    g_values : array of shape [2*n_seq]
        Interleaved positive/negative g-values.
    n_seq : int
        Number of sequences.
    alpha : float
        Threshold drawn as vertical line.
    ax : matplotlib Axes or None
        If None, creates a new figure.

    Returns
    -------
    fig, ax
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))
    else:
        fig = ax.get_figure()

    all_g = []
    for i in range(n_seq):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        best_g = min(pos_g if not np.isnan(pos_g) else 1.0,
                     neg_g if not np.isnan(neg_g) else 1.0)
        all_g.append(best_g)
    all_g = np.array(all_g)

    ax.hist(all_g[all_g < 1.0], bins=50, color="steelblue", alpha=0.7, edgecolor="white")
    ax.axvline(alpha, color="red", linestyle="--", linewidth=1.5,
               label=f"Threshold (g = {alpha})")
    n_below = int((all_g < alpha).sum())
    ax.annotate(f"{n_below} significant\n(g < {alpha})",
                xy=(0.25, 0.85), xycoords="axes fraction",
                fontsize=11, ha="center", color="steelblue", fontweight="bold")
    ax.set_xlabel("g-value (best direction per sequence)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of g-values")
    ax.legend()

    fig.tight_layout()
    return fig, ax


def sequence_space(seq_lengths, num_arms=None, ax=None):
    """Bar chart of unique sequences per length.

    Parameters
    ----------
    seq_lengths : array
        Sequence length for each sequence.
    num_arms : int or None
        If provided, overlays the theoretical maximum (num_arms^length) as a line.
    ax : matplotlib Axes or None
        If None, creates a new figure.

    Returns
    -------
    fig, ax
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 3))
    else:
        fig = ax.get_figure()

    unique_lens = sorted(set(seq_lengths))
    counts = [int(np.sum(seq_lengths == slen)) for slen in unique_lens]

    ax.bar(unique_lens, counts, color="#4488cc")
    ax.set_xlabel("Sequence length")
    ax.set_ylabel("# unique sequences")
    ax.set_title("Sequence Space by Length")

    if num_arms is not None:
        theoretical_max = [num_arms ** slen for slen in unique_lens]
        ax.plot(unique_lens, theoretical_max, "k--", alpha=0.6, label="Theoretical max")
        ax.legend(fontsize=8)

    fig.tight_layout()
    return fig, ax
