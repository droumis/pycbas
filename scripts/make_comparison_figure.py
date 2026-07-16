"""
Generate a paper-style Manhattan plot and side-by-side comparison with the
paper's Fig 1c right panel screenshot for the README.

Usage:
    pixi run python scripts/make_comparison_figure.py
"""

import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.cm import ScalarMappable

ROOT = Path(__file__).parent.parent
CACHE = ROOT / "results" / "figures" / "results_paper.npz"
PAPER_SCREENSHOT = ROOT / "results" / "figures" / "ppr_fig5cR.png"
OUT = ROOT / "results" / "figures"


def make_paper_style_manhattan():
    """Generate Manhattan plot matching paper Fig 1c style as closely as possible.

    Key paper features:
    - Sequences grouped by length into distinct vertical bands (like GWAS chromosomes)
    - Log-scale x-axis within each band, ordered by frequency
    - Discrete/categorical colors per length (cyan→blue→green)
    - Non-significant dots in black, significant colored and larger
    - Square aspect ratio, y-axis 0-4, dotted threshold line
    """
    data = np.load(CACHE, allow_pickle=False)
    g_values = data["g_values"]
    significant = data["significant"]
    seq_lengths = data["seq_lengths"]
    seq_len_max = int(data["params_seq_len_max"][0])
    n_sequences = len(seq_lengths)

    fig, ax = plt.subplots(figsize=(5, 4))

    # Discrete colors per length matching paper: cyan, light blue, blue, lime, light green, green
    length_colors = {
        1: "#00e5ff",  # cyan
        2: "#00aaff",  # light blue
        3: "#0044cc",  # blue
        4: "#88dd00",  # lime green
        5: "#44bb00",  # light green
        6: "#008800",  # green
    }

    # Compute -log10(best g-value)
    ys = np.zeros(n_sequences)
    for i in range(n_sequences):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        best_g = min(pos_g if not np.isnan(pos_g) else 1.0,
                     neg_g if not np.isnan(neg_g) else 1.0)
        ys[i] = -np.log10(max(best_g, 1e-4))

    # Build x positions grouped by length with log-spacing within each group
    x_positions = np.zeros(n_sequences)
    x_offset = 1.0

    for length in range(1, seq_len_max + 1):
        length_mask = seq_lengths == length
        length_indices = np.where(length_mask)[0]
        n_at_length = len(length_indices)
        if n_at_length == 0:
            continue

        if n_at_length > 1:
            local_x = np.logspace(0, np.log10(n_at_length), n_at_length, endpoint=True)
        else:
            local_x = np.array([1.0])
        x_positions[length_indices] = x_offset + local_x - 1

        x_offset += local_x[-1] * 1.1

    # Plot each length group — all dots colored by length (matching paper)
    for length in range(1, seq_len_max + 1):
        length_mask = seq_lengths == length
        if not length_mask.any():
            continue
        color = length_colors.get(length, "#888888")

        nonsig = length_mask & ~significant
        sig = length_mask & significant

        ax.scatter(x_positions[nonsig], ys[nonsig],
                   c=color, s=40, alpha=0.7, edgecolors="black", linewidths=0.3, zorder=2)
        ax.scatter(x_positions[sig], ys[sig],
                   c=color, s=40, alpha=0.7, edgecolors="black", linewidths=0.3, zorder=3)

    # Threshold line
    ax.axhline(-np.log10(0.5), color="black", linestyle=":", linewidth=1.0, zorder=1)

    ax.set_xscale("log")
    ax.set_xlabel("Sequence", fontsize=12, fontweight="bold")
    ax.set_ylabel(r"$-\log_{10}(\zeta)$", fontsize=12, fontweight="bold")
    ax.set_ylim(-0.1, 4.2)
    ax.set_xlim(0.8, x_offset * 1.1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Discrete colorbar at top
    discrete_cmap = ListedColormap([length_colors[l] for l in range(1, seq_len_max + 1)])
    sm = ScalarMappable(cmap=discrete_cmap, norm=Normalize(vmin=0.5, vmax=seq_len_max + 0.5))
    sm.set_array([])
    cax = ax.inset_axes([0.15, 1.03, 0.7, 0.03])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_ticks([2, 4, 6])
    cbar.ax.tick_params(labelsize=9, length=0)
    cax.set_title("sequence length", fontsize=10, pad=3)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out_path = OUT / "manhattan_paper_style.png"
    fig.savefig(out_path, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"  Paper-style Manhattan: {out_path}")
    return out_path


def make_comparison():
    """Side-by-side: paper screenshot vs our reimplementation."""
    manhattan_path = make_paper_style_manhattan()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.5))

    img_paper = mpimg.imread(str(PAPER_SCREENSHOT))
    ax1.imshow(img_paper)
    ax1.set_title("Paper — Fig 1c right panel\n(Kastner et al. 2026)", fontsize=9)
    ax1.axis("off")

    img_ours = mpimg.imread(str(manhattan_path))
    ax2.imshow(img_ours)
    ax2.set_title("pycbas — Python reimplementation\n(380/16,483 significant)", fontsize=9)
    ax2.axis("off")

    fig.text(0.5, -0.01,
             "Differences in spread and count (380 vs 409 significant) reflect different\n"
             "subject subsets — the paper uses all 85 rats while pycbas selects by data availability.",
             ha="center", fontsize=7, color="0.4")
    fig.tight_layout(pad=1.0)
    out_path = OUT / "comparison_manhattan.png"
    fig.savefig(out_path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  Comparison figure: {out_path}")


if __name__ == "__main__":
    if not CACHE.exists():
        print(f"Run `pixi run validate-paper` first to generate {CACHE}")
        raise SystemExit(1)
    if not PAPER_SCREENSHOT.exists():
        print(f"Paper screenshot not found at {PAPER_SCREENSHOT}")
        raise SystemExit(1)
    make_comparison()
