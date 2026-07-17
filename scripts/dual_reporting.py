"""
Dual-reporting demonstration: adaptive-k discovery set + k=1 ranking.

Produces a figure showing:
  Left: Manhattan at adaptive k (everything pegged at ceiling)
  Right: Manhattan at k=1 (nice spread, ranking preserved)

This is the concrete implementation of the suggestion for David.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from pycbas import (
    CBASParams,
    load_subject_data,
    build_count_matrix,
    compute_test_stats,
    bootstrap_test_stats,
    find_k_fwer,
    _prepare_null_sub,
    _stepdown_core,
)

ROOT_DIR = Path(__file__).parent.parent
DATA_DIR = ROOT_DIR / "data" / "flies"
FIG_DIR = ROOT_DIR / "results" / "flies" / "figures"


def load_flies():
    info_path = DATA_DIR / "flyInfo.txt"
    info = {}
    with open(info_path) as f:
        for line in f:
            parts = line.strip().split(",")
            info[int(parts[0])] = int(parts[1])

    subjects_data = []
    group_labels = []
    for fly_id in sorted(info.keys()):
        fpath = DATA_DIR / f"fly{fly_id}.txt"
        if fpath.exists():
            subjects_data.append(load_subject_data(fpath))
            group_labels.append(info[fly_id])

    return subjects_data, np.array(group_labels)


def compute_gvalues_at_fixed_k(test_stats, null_matrix, k, alpha=0.5):
    """Run step-down at a fixed k and return g-values (adjusted p-values)."""
    sorted_stats, sorted_indices, null_sub = _prepare_null_sub(test_stats, null_matrix)
    g_values_sorted = _stepdown_core(sorted_stats, null_sub, k, alpha)

    # Unsort back to original order
    g_values = np.full(len(test_stats), np.nan)
    for i, idx in enumerate(sorted_indices):
        g_values[idx] = g_values_sorted[i]

    return g_values


def plot_manhattan(ax, g_values, seq_lengths, title, seq_len_max=10, threshold=0.5):
    """Plot a manhattan-style plot from g-values."""
    n_seq = len(seq_lengths)
    xs = []
    ys = []
    colors = []

    for i in range(n_seq):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]

        # Take the more significant direction
        if np.isnan(pos_g) and np.isnan(neg_g):
            continue
        g = np.nanmin([pos_g, neg_g])

        neg_log_g = -np.log10(g) if g > 0 else 4.0
        neg_log_g = min(neg_log_g, 4.0)

        seq_len = int(seq_lengths[i])
        color_val = (seq_len - 1) / (seq_len_max - 1)

        xs.append(i + 1)
        ys.append(neg_log_g)
        colors.append(plt.cm.turbo(color_val))

    ax.scatter(xs, ys, c=colors, s=12, alpha=0.7, edgecolors="none")
    ax.axhline(-np.log10(threshold), color="black", linestyle=":", linewidth=0.8,
               label=f"ζ = {threshold}")
    ax.set_xlabel("Sequence")
    ax.set_ylabel("-log₁₀(ζ)")
    ax.set_title(title)
    ax.set_ylim(-0.1, 4.3)
    ax.set_xscale("log")

    n_sig = sum(1 for y in ys if y > -np.log10(threshold))
    ax.text(0.98, 0.02, f"{n_sig}/{n_seq} significant",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    return n_sig


def main():
    print("=" * 60)
    print("Dual Reporting: adaptive k vs fixed k=1")
    print("=" * 60)

    subjects_data, group_labels = load_flies()
    print(f"Loaded {len(subjects_data)} flies")

    params = CBASParams(num_arms=2, seq_len_max=10, criterion=250, resample_number=1000)
    group_indices = [
        np.where(group_labels == 0)[0],
        np.where(group_labels == 1)[0],
    ]

    print("Building count matrix...")
    sequences, count_matrix = build_count_matrix(subjects_data, params, contingency=1)
    seq_lengths = np.array([len(s) for s in sequences])
    n_seq = len(sequences)
    print(f"  {n_seq} sequences")

    print("Computing test statistics...")
    test_stats = compute_test_stats(count_matrix, group_indices)

    print("Bootstrapping (M=1000)...")
    null_matrix = bootstrap_test_stats(count_matrix, group_indices, params)

    # Adaptive k
    print("Running adaptive k-FWER...")
    g_values_adaptive, k_final = find_k_fwer(test_stats, null_matrix, params.alpha, params.gamma)
    print(f"  Adaptive k={k_final}")

    # Fixed k=1
    print("Running step-down at fixed k=1...")
    g_values_k1 = compute_gvalues_at_fixed_k(test_stats, null_matrix, k=1)

    # Fixed k=5
    print("Running step-down at fixed k=5...")
    g_values_k5 = compute_gvalues_at_fixed_k(test_stats, null_matrix, k=5)

    # Fixed k=20
    print("Running step-down at fixed k=20...")
    g_values_k20 = compute_gvalues_at_fixed_k(test_stats, null_matrix, k=20)

    # Make the figure
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    n_sig = plot_manhattan(axes[0, 0], g_values_adaptive, seq_lengths,
                           f"Adaptive k={k_final} (current method)\nAll sequences significant — no ranking",
                           seq_len_max=10)
    print(f"  Adaptive k={k_final}: {n_sig} significant")

    n_sig = plot_manhattan(axes[0, 1], g_values_k1, seq_lengths,
                           "Fixed k=1 (FWER)\nRanking preserved — spread of ζ values",
                           seq_len_max=10)
    print(f"  Fixed k=1: {n_sig} significant")

    n_sig = plot_manhattan(axes[1, 0], g_values_k5, seq_lengths,
                           "Fixed k=5\nMore discoveries, still informative ranking",
                           seq_len_max=10)
    print(f"  Fixed k=5: {n_sig} significant")

    n_sig = plot_manhattan(axes[1, 1], g_values_k20, seq_lengths,
                           "Fixed k=20 (≈ paper's result)\nBalances discovery and ranking",
                           seq_len_max=10)
    print(f"  Fixed k=20: {n_sig} significant")

    fig.suptitle("Fly CBAS: Dual reporting preserves Manhattan ranking\n"
                 "Suggestion: report adaptive-k discovery set + k=1 ζ for ranking",
                 fontsize=12, fontweight="bold")

    plt.tight_layout(rect=[0, 0.05, 1, 0.94])

    # Add categorical colorbar below all panels
    cbar_ax = fig.add_axes([0.25, 0.015, 0.5, 0.02])
    from matplotlib.colors import BoundaryNorm
    boundaries = np.arange(0.5, 11.5, 1)
    norm = BoundaryNorm(boundaries, ncolors=256)
    sm = plt.cm.ScalarMappable(cmap="turbo", norm=norm)
    cb = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal",
                      label="Sequence length")
    cb.set_ticks(range(1, 11))
    cb.set_ticklabels([str(i) for i in range(1, 11)])
    out_path = FIG_DIR / "dual_reporting_manhattan.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nFigure saved: {out_path}")


if __name__ == "__main__":
    main()
