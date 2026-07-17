"""
Persistence analysis at different k values.

Shows how the interpretive conclusion (which persistence levels differ
between strains) depends on k. If the pattern is robust across k,
the interpretation is solid. If it changes, the choice of k matters
for scientific conclusions.
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
    """Run step-down at a fixed k and return g-values."""
    sorted_stats, sorted_indices, null_sub = _prepare_null_sub(test_stats, null_matrix)
    g_values_sorted = _stepdown_core(sorted_stats, null_sub, k, alpha)

    g_values = np.full(len(test_stats), np.nan)
    for i, idx in enumerate(sorted_indices):
        g_values[idx] = g_values_sorted[i]

    return g_values


def max_persistence(seq):
    """Max consecutive identical symbols in a sequence (tuple of ints)."""
    if len(seq) <= 1:
        return 1
    max_run = 1
    current_run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1
    return max_run


def persistence_breakdown(sequences, g_values, threshold=0.5):
    """Get persistence distribution of significant sequences, split by direction."""
    n_seq = len(sequences)
    ca_gt_w1118 = {}  # persistence -> count where CA > w1118 (positive direction)
    w1118_gt_ca = {}  # persistence -> count where w1118 > CA (negative direction)

    for i in range(n_seq):
        pos_g = g_values[i * 2]   # CA > w1118
        neg_g = g_values[i * 2 + 1]  # w1118 > CA

        persist = max_persistence(sequences[i])

        if not np.isnan(pos_g) and pos_g < threshold:
            ca_gt_w1118[persist] = ca_gt_w1118.get(persist, 0) + 1
        if not np.isnan(neg_g) and neg_g < threshold:
            w1118_gt_ca[persist] = w1118_gt_ca.get(persist, 0) + 1

    return ca_gt_w1118, w1118_gt_ca


def main():
    print("=" * 60)
    print("Persistence Analysis at Different k Values")
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
    n_seq = len(sequences)
    print(f"  {n_seq} sequences")

    print("Computing test statistics...")
    test_stats = compute_test_stats(count_matrix, group_indices)

    print("Bootstrapping (M=1000)...")
    null_matrix = bootstrap_test_stats(count_matrix, group_indices, params)

    # Compute g-values at different k
    k_values = [1, 5, 20, 103]
    k_labels = ["k=1 (FWER)\n1,245 sig", "k=5\n1,467 sig",
                "k=20 (≈paper)\n1,632 sig", "k=103 (adaptive)\n2,046 sig"]

    all_gvalues = {}
    for k in k_values:
        print(f"  Step-down at k={k}...")
        if k == 103:
            g, _ = find_k_fwer(test_stats, null_matrix, params.alpha, params.gamma)
            all_gvalues[k] = g
        else:
            all_gvalues[k] = compute_gvalues_at_fixed_k(test_stats, null_matrix, k)

    # Also compute persistence for ALL sequences (reference)
    all_persist = [max_persistence(s) for s in sequences]
    persist_range = range(1, 11)

    # Make figure
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    for idx, (k, label) in enumerate(zip(k_values, k_labels)):
        ax = axes[idx // 2, idx % 2]
        ca_counts, w1118_counts = persistence_breakdown(sequences, all_gvalues[k])

        ca_vals = [ca_counts.get(p, 0) for p in persist_range]
        w1118_vals = [w1118_counts.get(p, 0) for p in persist_range]

        x = np.arange(len(persist_range))
        width = 0.35
        bars1 = ax.bar(x - width / 2, ca_vals, width, color="#cc5522", label="CA > w1118")
        bars2 = ax.bar(x + width / 2, w1118_vals, width, color="#2266cc", label="w1118 > CA")

        ax.set_xticks(x)
        ax.set_xticklabels([str(p) for p in persist_range])
        ax.set_xlabel("Max same turns in a row (persistence)")
        ax.set_ylabel("Significant sequences")
        ax.set_title(label)
        ax.legend(fontsize=8)

        # Annotate peak
        if w1118_vals:
            peak_w = np.argmax(w1118_vals)
            ax.annotate(f"{w1118_vals[peak_w]}", (x[peak_w] + width/2, w1118_vals[peak_w]),
                       ha="center", va="bottom", fontsize=7, color="#2266cc")
        if ca_vals:
            peak_c = np.argmax(ca_vals)
            ax.annotate(f"{ca_vals[peak_c]}", (x[peak_c] - width/2, ca_vals[peak_c]),
                       ha="center", va="bottom", fontsize=7, color="#cc5522")

    fig.suptitle("Does the persistence interpretation depend on k?\n"
                 "(cf. Paper Fig 2: w1118 = low persistence, CA = high persistence)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    out_path = FIG_DIR / "persistence_by_k.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nFigure saved: {out_path}")

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY: Ratio of w1118>CA to CA>w1118 at low vs high persistence")
    print("=" * 60)
    for k in k_values:
        ca_counts, w1118_counts = persistence_breakdown(sequences, all_gvalues[k])
        low_persist_w = sum(w1118_counts.get(p, 0) for p in [1, 2, 3])
        low_persist_c = sum(ca_counts.get(p, 0) for p in [1, 2, 3])
        high_persist_w = sum(w1118_counts.get(p, 0) for p in [5, 6, 7, 8, 9, 10])
        high_persist_c = sum(ca_counts.get(p, 0) for p in [5, 6, 7, 8, 9, 10])
        print(f"  k={k:>3}: low persist w1118:CA = {low_persist_w}:{low_persist_c} "
              f"({low_persist_w/(low_persist_c+1e-9):.1f}x) | "
              f"high persist CA:w1118 = {high_persist_c}:{high_persist_w} "
              f"({high_persist_c/(high_persist_w+1e-9):.1f}x)")


if __name__ == "__main__":
    main()
