"""
Subsample analysis for fly data, replicating Figure 3 from the paper.

Runs CBAS at multiple sample sizes to show scaling behavior.
Uses M=1000 (paper uses this for repeats).

Usage:
    pixi run flies-subsample
"""

import sys
import time
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from pycbas import (
    CBASParams,
    load_subject_data,
    build_count_matrix,
    compute_test_stats,
    bootstrap_test_stats,
    find_k_fwer,
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


def run_subsample(subjects_data, group_labels, n_per_group, params, rng):
    """Run CBAS on a random subsample of n_per_group from each group."""
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
    g_values, k_final = find_k_fwer(test_stats, null_matrix, params.alpha, params.gamma)

    n_sig = 0
    for i in range(len(sequences)):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        if (not np.isnan(pos_g) and pos_g < params.alpha) or \
           (not np.isnan(neg_g) and neg_g < params.alpha):
            n_sig += 1

    return n_sig, len(sequences), k_final


def main():
    print("=" * 60)
    print("Fly CBAS Subsampling Analysis (cf. Paper Fig 3)")
    print("=" * 60)

    subjects_data, group_labels = load_flies()
    n_ca = int((group_labels == 0).sum())
    n_w1118 = int((group_labels == 1).sum())
    print(f"\nFull dataset: {len(subjects_data)} flies ({n_ca} CA, {n_w1118} w1118)")

    params = CBASParams(num_arms=2, seq_len_max=10, criterion=250, resample_number=1000)
    sample_sizes = [40, 80, 160, 320]
    n_repeats = 10
    rng = np.random.default_rng(42)

    results = {n: [] for n in sample_sizes}

    for n_per_group in sample_sizes:
        print(f"\n--- {n_per_group} flies per group ({n_repeats} repeats) ---")
        for rep in range(n_repeats):
            t0 = time.perf_counter()
            n_sig, n_seq, k_final = run_subsample(
                subjects_data, group_labels, n_per_group, params, rng
            )
            elapsed = time.perf_counter() - t0
            results[n_per_group].append((n_sig, n_seq, k_final))
            print(f"  Rep {rep+1}: {n_sig}/{n_seq} sig (k={k_final}) [{elapsed:.1f}s]")

    # Also run full dataset at M=1000 for reference
    print(f"\n--- Full dataset (M=1000) ---")
    group_indices = [
        np.where(group_labels == 0)[0],
        np.where(group_labels == 1)[0],
    ]
    t0 = time.perf_counter()
    sequences, count_matrix = build_count_matrix(subjects_data, params, contingency=1)
    test_stats = compute_test_stats(count_matrix, group_indices)
    null_matrix = bootstrap_test_stats(count_matrix, group_indices, params)
    g_values, k_final = find_k_fwer(test_stats, null_matrix, params.alpha, params.gamma)
    n_sig_full = 0
    for i in range(len(sequences)):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        if (not np.isnan(pos_g) and pos_g < params.alpha) or \
           (not np.isnan(neg_g) and neg_g < params.alpha):
            n_sig_full += 1
    elapsed = time.perf_counter() - t0
    print(f"  {n_sig_full}/{len(sequences)} sig (k={k_final}) [{elapsed:.1f}s]")

    # Make figure
    import matplotlib.pyplot as plt

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Panel a: number of significant sequences vs sample size
    ax = axes[0]
    positions = []
    sig_counts = []
    for n_per_group in sample_sizes:
        sigs = [r[0] for r in results[n_per_group]]
        positions.append(n_per_group)
        sig_counts.append(sigs)

    bp = ax.boxplot(sig_counts, positions=sample_sizes, widths=15, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("#4488cc")
        patch.set_alpha(0.5)
    ax.axhline(n_sig_full, color="black", linestyle="--", linewidth=1,
               label=f"Full dataset ({n_sig_full})")
    ax.set_xlabel("Flies per group")
    ax.set_ylabel("Significant sequences")
    ax.set_title("CBAS Scaling with Sample Size")
    ax.legend()

    # Panel b: fraction significant
    ax = axes[1]
    frac_sig = []
    for n_per_group in sample_sizes:
        fracs = [r[0] / r[1] if r[1] > 0 else 0 for r in results[n_per_group]]
        frac_sig.append(fracs)

    bp = ax.boxplot(frac_sig, positions=sample_sizes, widths=15, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("#44cc88")
        patch.set_alpha(0.5)
    ax.axhline(n_sig_full / len(sequences), color="black", linestyle="--", linewidth=1,
               label=f"Full dataset ({n_sig_full/len(sequences)*100:.0f}%)")
    ax.set_xlabel("Flies per group")
    ax.set_ylabel("Fraction significant")
    ax.set_title("Fraction of Sequences Significant")
    ax.set_ylim(0, 1.05)
    ax.legend()

    plt.tight_layout()
    plt.savefig(FIG_DIR / "subsample_scaling.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nFigure saved: {FIG_DIR / 'subsample_scaling.png'}")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"{'Sample size':<15} {'Median sig':<12} {'Median frac':<12} {'Median k'}")
    print("-" * 55)
    for n_per_group in sample_sizes:
        sigs = [r[0] for r in results[n_per_group]]
        fracs = [r[0] / r[1] if r[1] > 0 else 0 for r in results[n_per_group]]
        ks = [r[2] for r in results[n_per_group]]
        print(f"{n_per_group:<15} {np.median(sigs):<12.0f} {np.median(fracs):<12.2f} {np.median(ks):.0f}")
    print(f"{'Full (759/807)':<15} {n_sig_full:<12} {n_sig_full/len(sequences):<12.2f} {k_final}")


if __name__ == "__main__":
    main()
