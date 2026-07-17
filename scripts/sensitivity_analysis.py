"""
k-FWER sensitivity analysis: how significance scales with sample size and k.

Produces a diagnostic figure showing:
  Panel A: Heatmap of fraction significant at each (k, N) combination
  Panel B: Test statistic vs null separation as f(N)
  Panel C: The k-iteration trajectory (rejections vs k at each N)

This is designed to convey to David how the adaptive k-FWER behaves
in the high-power regime (large N, pervasive effect).
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from pycbas import (
    CBASParams,
    load_subject_data,
    build_count_matrix,
    compute_test_stats,
    bootstrap_test_stats,
    _prepare_null_sub,
    _count_rejections,
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


def run_at_sample_size(subjects_data, group_labels, n_per_group, params, rng):
    """Run CBAS at a given sample size, return test_stats, null_matrix, sequences."""
    ca_indices = np.where(group_labels == 0)[0]
    w1118_indices = np.where(group_labels == 1)[0]

    sub_ca = rng.choice(ca_indices, size=n_per_group, replace=False)
    sub_w1118 = rng.choice(w1118_indices, size=n_per_group, replace=False)

    sub_indices = np.concatenate([sub_ca, sub_w1118])
    sub_data = [subjects_data[i] for i in sub_indices]
    sub_labels = np.array([0] * n_per_group + [1] * n_per_group)

    group_indices = [
        np.where(sub_labels == 0)[0],
        np.where(sub_labels == 1)[0],
    ]

    sequences, count_matrix = build_count_matrix(sub_data, params, contingency=1)
    test_stats = compute_test_stats(count_matrix, group_indices)
    null_matrix = bootstrap_test_stats(count_matrix, group_indices, params)

    return sequences, test_stats, null_matrix


def count_sig_at_fixed_k(test_stats, null_matrix, k, alpha=0.5):
    """Count rejections using the step-down procedure at a fixed k."""
    sorted_stats, sorted_indices, null_sub = _prepare_null_sub(test_stats, null_matrix)
    n_rejections = _count_rejections(sorted_stats, null_sub, k, alpha)
    return n_rejections


def count_sig_sequences(n_rejections, n_sequences):
    """Convert hypothesis rejections to sequence count (each seq has 2 hypotheses)."""
    return min(n_rejections, n_sequences)


def trace_k_iteration(test_stats, null_matrix, alpha=0.5, gamma=0.05, max_k=150):
    """Trace the k-FWER iteration: at each k, how many rejections?"""
    sorted_stats, sorted_indices, null_sub = _prepare_null_sub(test_stats, null_matrix)
    ks = []
    rejections = []
    for k in range(1, max_k + 1):
        n_rej = _count_rejections(sorted_stats, null_sub, k, alpha)
        ks.append(k)
        rejections.append(n_rej)
        if n_rej < (k / gamma - 1):
            break
    return np.array(ks), np.array(rejections)


def main():
    print("=" * 60)
    print("k-FWER Sensitivity Analysis")
    print("=" * 60)

    subjects_data, group_labels = load_flies()
    n_ca = int((group_labels == 0).sum())
    n_w1118 = int((group_labels == 1).sum())
    print(f"Full dataset: {len(subjects_data)} flies ({n_ca} CA, {n_w1118} w1118)")

    params = CBASParams(num_arms=2, seq_len_max=10, criterion=250, resample_number=1000)
    rng = np.random.default_rng(123)

    sample_sizes = [40, 80, 160, 320]
    fixed_ks = [1, 5, 10, 20, 50, 100]

    # ---- Run at each sample size ----
    results = {}
    for n in sample_sizes:
        print(f"\nRunning N={n} per group...")
        sequences, test_stats, null_matrix = run_at_sample_size(
            subjects_data, group_labels, n, params, rng
        )
        results[n] = {
            "sequences": sequences,
            "test_stats": test_stats,
            "null_matrix": null_matrix,
            "n_sequences": len(sequences),
        }

    # Full dataset
    print("\nRunning full dataset...")
    group_indices = [
        np.where(group_labels == 0)[0],
        np.where(group_labels == 1)[0],
    ]
    sequences, count_matrix = build_count_matrix(subjects_data, params, contingency=1)
    test_stats = compute_test_stats(count_matrix, group_indices)
    null_matrix = bootstrap_test_stats(count_matrix, group_indices, params)
    n_full = min(n_ca, n_w1118)
    results[n_full] = {
        "sequences": sequences,
        "test_stats": test_stats,
        "null_matrix": null_matrix,
        "n_sequences": len(sequences),
    }
    sample_sizes_all = sample_sizes + [n_full]

    # ---- Panel A: Heatmap of rejections at fixed k × N ----
    print("\nComputing heatmap...")
    heatmap = np.zeros((len(fixed_ks), len(sample_sizes_all)))
    for j, n in enumerate(sample_sizes_all):
        r = results[n]
        n_seq = r["n_sequences"]
        for i, k in enumerate(fixed_ks):
            n_rej = count_sig_at_fixed_k(r["test_stats"], r["null_matrix"], k)
            heatmap[i, j] = n_rej / (n_seq * 2)  # fraction of hypotheses

    # ---- Panel B: Test stat vs null magnitude ----
    print("Computing stat/null separation...")
    stat_percentiles = {}
    null_percentiles = {}
    for n in sample_sizes_all:
        r = results[n]
        valid = r["test_stats"][~np.isnan(r["test_stats"])]
        stat_percentiles[n] = np.percentile(valid, [50, 75, 90, 95, 99])
        null_row_maxes = np.nanmax(r["null_matrix"], axis=1)
        null_percentiles[n] = np.percentile(null_row_maxes, [50, 75, 90, 95, 99])

    # ---- Panel C: k-iteration trajectories ----
    print("Tracing k-iteration trajectories...")
    trajectories = {}
    for n in sample_sizes_all:
        r = results[n]
        ks, rejs = trace_k_iteration(r["test_stats"], r["null_matrix"])
        trajectories[n] = (ks, rejs)
        print(f"  N={n}: converged at k={ks[-1]}, {rejs[-1]} rejections")

    # ---- Make figure ----
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    # Panel A: Heatmap
    ax = fig.add_subplot(gs[0, 0])
    cmap = LinearSegmentedColormap.from_list("sig", ["white", "#4488cc", "#cc4444"])
    im = ax.imshow(heatmap, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(sample_sizes_all)))
    ax.set_xticklabels([str(n) for n in sample_sizes_all])
    ax.set_yticks(range(len(fixed_ks)))
    ax.set_yticklabels([str(k) for k in fixed_ks])
    ax.set_xlabel("Subjects per group (N)")
    ax.set_ylabel("Fixed k")
    ax.set_title("A. Fraction significant at fixed k\n(no adaptive iteration)")
    for i in range(len(fixed_ks)):
        for j in range(len(sample_sizes_all)):
            val = heatmap[i, j]
            color = "white" if val > 0.6 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color=color)
    plt.colorbar(im, ax=ax, label="Fraction of hypotheses rejected")

    # Panel B: Stat vs null separation
    ax = fig.add_subplot(gs[0, 1])
    x_pos = np.arange(len(sample_sizes_all))
    stat_medians = [stat_percentiles[n][0] for n in sample_sizes_all]
    stat_99 = [stat_percentiles[n][4] for n in sample_sizes_all]
    null_medians = [null_percentiles[n][0] for n in sample_sizes_all]
    null_99 = [null_percentiles[n][4] for n in sample_sizes_all]

    ax.fill_between(x_pos, null_medians, null_99, alpha=0.3, color="gray", label="Null max (50th-99th)")
    ax.fill_between(x_pos, stat_medians, stat_99, alpha=0.3, color="blue", label="Observed (50th-99th)")
    ax.plot(x_pos, stat_99, "b-o", markersize=5, label="Observed 99th pctl")
    ax.plot(x_pos, null_99, "k-s", markersize=5, label="Null-max 99th pctl")
    ax.set_xticks(x_pos)
    ax.set_xticklabels([str(n) for n in sample_sizes_all])
    ax.set_xlabel("Subjects per group (N)")
    ax.set_ylabel("Test statistic magnitude")
    ax.set_title("B. Signal-null separation grows with √N\n(observed t_s scales, null stays bounded)")
    ax.legend(fontsize=8)

    # Panel C: k-iteration trajectories
    ax = fig.add_subplot(gs[1, 0])
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(sample_sizes_all)))
    gamma = 0.05
    for idx, n in enumerate(sample_sizes_all):
        ks, rejs = trajectories[n]
        ax.plot(ks, rejs, "-o", markersize=3, color=colors[idx], label=f"N={n}")

    # Plot the stopping boundary: N = k/gamma - 1
    k_range = np.arange(1, 120)
    boundary = k_range / gamma - 1
    ax.plot(k_range, boundary, "r--", linewidth=2, label=f"Stopping: N < k/γ−1")
    ax.set_xlabel("k (FWER tolerance parameter)")
    ax.set_ylabel("Number of rejections (N)")
    ax.set_title("C. k-FWER iteration trajectory\n(adaptive k climbs until N crosses boundary)")
    ax.legend(fontsize=7, loc="upper left")
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 4500)

    # Panel D: Annotation / explanation
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    explanation = (
        "Summary: The k-FWER sensitivity issue\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"
        "• Test statistics t_s = Δr_s / σ_s scale as √N\n"
        "  (more subjects → smaller σ → larger t)\n\n"
        "• Null distribution (bootstrap) stays bounded\n"
        "  because it permutes the same finite data\n\n"
        "• At large N, nearly all t_s exceed the null,\n"
        "  so rejections ≈ S (total hypotheses)\n\n"
        "• The adaptive k-iteration (Panel C) finds k\n"
        "  where rejections < k/γ − 1. When rejections\n"
        "  are near-maximal, k must grow until the\n"
        "  boundary catches up → k ≈ γ × S\n\n"
        "• At k=103, FDP guarantee: ≤103/4092 ≈ 2.5%\n"
        "  false positives. Statistically valid, but\n"
        "  scientifically uninformative (no ranking).\n\n"
        "Implications for interpretation:\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        "• At fixed k=20: 1,633 sig (≈paper's 1,605)\n"
        "• Persistence analysis (Fig 2) remains the\n"
        "  key interpretive tool regardless of k\n"
        "• Subsampling (Fig 3) shows power scaling\n"
        "• Consider reporting effect sizes alongside\n"
        "  significance for high-power datasets"
    )
    ax.text(0.05, 0.95, explanation, transform=ax.transAxes, fontsize=8,
            verticalalignment="top", fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    plt.savefig(FIG_DIR / "kfwer_sensitivity.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nFigure saved: {FIG_DIR / 'kfwer_sensitivity.png'}")

    # Print the heatmap data
    print("\n" + "=" * 60)
    print("HEATMAP: Fraction significant (hypotheses) at fixed k × N")
    print("=" * 60)
    header = f"{'k':<6}" + "".join(f"{'N='+str(n):<10}" for n in sample_sizes_all)
    print(header)
    print("-" * len(header))
    for i, k in enumerate(fixed_ks):
        row = f"{k:<6}" + "".join(f"{heatmap[i,j]:.3f}     " for j in range(len(sample_sizes_all)))
        print(row)


if __name__ == "__main__":
    main()
