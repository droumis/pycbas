import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.cm import ScalarMappable


def _paper_colors(n_lengths):
    """Generate colors matching the paper's cyan→blue→green gradient."""
    if n_lengths <= 4:
        palette = ["#00e5ff", "#0099ff", "#0044dd", "#00aa44"]
    elif n_lengths <= 6:
        palette = ["#00e5ff", "#00bbff", "#0066dd", "#00cc44", "#009922", "#006600"]
    else:
        palette = [
            "#00e5ff", "#00ccff", "#0099ff", "#0055dd",
            "#00bb33", "#00aa00", "#228800",
            "#ccaa00", "#dd6600", "#cc0000",
        ]
    return palette[:n_lengths]


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
        Mapping from sequence length to hex color string. If None, auto-generates
        from a paper-matched gradient.
    title : str or None
        Plot title. If None, no title is set.
    ax : matplotlib Axes or None
        If None, creates a new figure.

    Returns
    -------
    fig, ax
    """
    n_seq = len(seq_lengths)
    unique_lens = sorted(set(seq_lengths))
    n_lengths = len(unique_lens)

    if colors is None:
        palette = _paper_colors(n_lengths)
        colors = {slen: palette[i] for i, slen in enumerate(unique_lens)}

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

    # Paper-style: single continuous x-axis, sequences ordered by length then
    # by index within each length. X position = cumulative rank (1-based).
    x_pos = np.zeros(n_seq)
    rank = 1
    for slen in unique_lens:
        mask = seq_lengths == slen
        indices = np.where(mask)[0]
        for j, idx in enumerate(indices):
            x_pos[idx] = rank
            rank += 1

    valid = ~np.isnan(neg_log_g)
    for slen in unique_lens:
        mask = (seq_lengths == slen) & valid
        c = colors.get(slen, "#999999")
        ax.scatter(x_pos[mask], neg_log_g[mask], s=20, alpha=0.8, c=c,
                   edgecolors="black", linewidths=0.3)

    threshold = -np.log10(alpha)
    ax.axhline(threshold, color="black", linestyle=":", linewidth=1.0, alpha=0.7)
    ax.set_ylabel(r"$-\log_{10}(\zeta)$", fontsize=11)
    ax.set_xlabel("Sequence", fontsize=11, fontweight="bold")
    ax.set_ylim(-0.1, None)
    ax.set_xscale("log")
    ax.set_xlim(0.8, n_seq * 1.1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if title is not None:
        ax.set_title(title, fontsize=11)

    # Vertical colorbar on the right
    color_list = [colors[slen] for slen in unique_lens]
    cmap = ListedColormap(color_list)
    boundaries = np.arange(0.5, n_lengths + 1.5)
    norm = BoundaryNorm(boundaries, cmap.N)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", fraction=0.03, pad=0.02,
                        ticks=np.arange(1, n_lengths + 1))
    cbar.ax.set_yticklabels([str(s) for s in unique_lens], fontsize=8)
    cbar.ax.set_ylabel("sequence length", fontsize=9, rotation=270, labelpad=12)

    fig.tight_layout()
    return fig, ax
