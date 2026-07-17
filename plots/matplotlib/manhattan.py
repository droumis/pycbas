import numpy as np
import matplotlib.pyplot as plt


FLY_COLORS = {
    1: "#00e5ff", 2: "#00aaff", 3: "#0044cc", 4: "#88dd00",
    5: "#44bb00", 6: "#008800", 7: "#ff6600", 8: "#cc3300",
    9: "#990000", 10: "#660066",
}


def manhattan_plot(g_values, seq_lengths, alpha=0.5, colors=None, title=None, ax=None):
    """Manhattan-style plot of -log10(best g-value) grouped by sequence length.

    Parameters
    ----------
    g_values : array of shape [2*n_seq]
        Interleaved positive/negative g-values: g_values[i*2] = positive,
        g_values[i*2+1] = negative direction.
    seq_lengths : array of length n_seq
        Sequence length for each sequence.
    alpha : float
        Significance threshold drawn as horizontal line at -log10(alpha).
    colors : dict or None
        Mapping from sequence length to hex color string. If None, uses FLY_COLORS.
    title : str or None
        Plot title. If None, no title is set.
    ax : matplotlib Axes or None
        If None, creates a new figure.

    Returns
    -------
    fig, ax
    """
    if colors is None:
        colors = FLY_COLORS

    n_seq = len(seq_lengths)
    neg_log_g = np.full(n_seq, np.nan)

    for i in range(n_seq):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        best_g = np.nan
        if not np.isnan(pos_g) and not np.isnan(neg_g):
            best_g = min(pos_g, neg_g)
        elif not np.isnan(pos_g):
            best_g = pos_g
        elif not np.isnan(neg_g):
            best_g = neg_g
        if not np.isnan(best_g) and best_g > 0:
            neg_log_g[i] = -np.log10(best_g)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))
    else:
        fig = ax.get_figure()

    unique_lens = sorted(set(seq_lengths))
    x_pos = np.zeros(n_seq)
    band_width = 1.0
    gap = 0.3

    for band_idx, slen in enumerate(unique_lens):
        mask = seq_lengths == slen
        indices = np.where(mask)[0]
        n_in_band = len(indices)
        if n_in_band > 1:
            positions = np.logspace(0, np.log10(n_in_band), n_in_band)
            positions = (positions - positions.min()) / (positions.max() - positions.min())
        else:
            positions = np.array([0.5])
        band_start = band_idx * (band_width + gap)
        for j, idx in enumerate(indices):
            x_pos[idx] = band_start + positions[j] * band_width

    valid = ~np.isnan(neg_log_g)
    for slen in unique_lens:
        mask = (seq_lengths == slen) & valid
        c = colors.get(slen, "#999999")
        ax.scatter(x_pos[mask], neg_log_g[mask], s=20, alpha=0.7, c=c,
                   edgecolors="black", linewidths=0.2, label=f"len={slen}")

    threshold = -np.log10(alpha)
    ax.axhline(threshold, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_ylabel(r"-log$_{10}$($\zeta$)")
    ax.set_xlabel("Sequence length")

    if title is not None:
        ax.set_title(title)

    ax.legend(loc="upper right", fontsize=7, ncol=2)

    xtick_pos = [i * (band_width + gap) + band_width / 2 for i in range(len(unique_lens))]
    ax.set_xticks(xtick_pos)
    ax.set_xticklabels(unique_lens)

    fig.tight_layout()
    return fig, ax
