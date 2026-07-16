"""
Benchmark the CBAS pipeline stages.

Reports timing for each stage so you can see where computation time is spent.
Uses paper-matched parameters by default (seq_len_max=6, M=10,000, 85 subjects).

Usage:
    pixi run benchmark
    pixi run python scripts/benchmark.py --quick   # reduced params for fast check
"""

import argparse
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
    _stepdown_core,
    _prepare_null_sub,
    _count_rejections,
)

DATA_DIR = Path(__file__).parent.parent / "igor_cbas" / "data"


def load_rats(n_ctrl_max=None, n_les_max=None):
    ctrl_data, les_data = [], []
    for f in sorted(DATA_DIR.glob("*.txt")):
        name = f.stem
        if "Control" in name:
            ctrl_data.append(load_subject_data(f))
        elif "Lesion" in name:
            les_data.append(load_subject_data(f))
    if n_ctrl_max:
        ctrl_data = ctrl_data[:n_ctrl_max]
    if n_les_max:
        les_data = les_data[:n_les_max]
    subjects_data = ctrl_data + les_data
    group_labels = np.array([0] * len(ctrl_data) + [1] * len(les_data))
    return subjects_data, group_labels


def run_benchmark(quick=False):
    if quick:
        n_ctrl_max, n_les_max = None, None
        params = CBASParams(seq_len_max=4, criterion=800, resample_number=1000)
        label = "Quick (seq_len_max=4, M=1,000, all subjects)"
    else:
        n_ctrl_max, n_les_max = 46, 39
        params = CBASParams(seq_len_max=6, criterion=800, resample_number=10000)
        label = "Paper params (seq_len_max=6, M=10,000, 85 subjects)"

    print("=" * 60)
    print(f"CBAS Benchmark — {label}")
    print("=" * 60)

    subjects_data, group_labels = load_rats(n_ctrl_max, n_les_max)
    n_subjects = len(subjects_data)
    n_ctrl = int((group_labels == 0).sum())
    n_les = int((group_labels == 1).sum())
    print(f"\nData: {n_subjects} subjects ({n_ctrl} ctrl, {n_les} les)")
    print(f"Params: seq_len_max={params.seq_len_max}, criterion={params.criterion}, "
          f"M={params.resample_number}")

    group_indices = [
        np.where(group_labels == 0)[0],
        np.where(group_labels == 1)[0],
    ]

    timings = {}

    # --- Stage 1: Count matrix ---
    t0 = time.perf_counter()
    sequences, count_matrix = build_count_matrix(subjects_data, params)
    timings["build_count_matrix"] = time.perf_counter() - t0
    n_seq = count_matrix.shape[1]
    print(f"\n  Count matrix: {n_subjects} x {n_seq} ({timings['build_count_matrix']:.2f}s)")

    # --- Stage 2: Test statistics ---
    t0 = time.perf_counter()
    test_stats = compute_test_stats(count_matrix, group_indices)
    timings["compute_test_stats"] = time.perf_counter() - t0
    n_valid = int(np.sum(~np.isnan(test_stats)))
    print(f"  Test stats: {n_valid} valid of {len(test_stats)} ({timings['compute_test_stats']:.3f}s)")

    # --- Stage 3: Bootstrap ---
    t0 = time.perf_counter()
    null_matrix = bootstrap_test_stats(count_matrix, group_indices, params)
    timings["bootstrap"] = time.perf_counter() - t0
    print(f"  Bootstrap ({params.resample_number} resamples): {timings['bootstrap']:.2f}s")

    # --- Stage 4: k-FWER ---
    # Warm up numba (exclude compilation from timing)
    sorted_stats, sorted_indices, null_sub = _prepare_null_sub(test_stats, null_matrix)
    _ = _stepdown_core(sorted_stats[:10], null_sub[:, :10].copy(), 1, 0.5)

    t0 = time.perf_counter()
    g_values, k_final = find_k_fwer(test_stats, null_matrix, params.alpha, params.gamma)
    timings["k_fwer"] = time.perf_counter() - t0
    n_sig = int(np.sum(g_values < 0.5))
    print(f"  k-FWER (converged k={k_final}): {timings['k_fwer']:.2f}s")

    timings["total"] = sum(timings.values())

    # --- Detailed k-FWER breakdown ---
    print(f"\n  k-FWER iteration detail:")
    sorted_stats, sorted_indices, null_sub = _prepare_null_sub(test_stats, null_matrix)
    k = 1
    prev_rejections = None
    iteration = 0
    while True:
        t0 = time.perf_counter()
        rejections = _count_rejections(sorted_stats, null_sub, k, params.alpha)
        elapsed = time.perf_counter() - t0
        iteration += 1
        print(f"    k={k:>2d}: {rejections} rejections ({elapsed:.2f}s)")

        if prev_rejections is not None and rejections < prev_rejections:
            break
        required_k = int(np.ceil((rejections + 1) * params.gamma))
        if required_k <= k:
            break
        prev_rejections = rejections
        k = required_k

    t0 = time.perf_counter()
    pvals = _stepdown_core(sorted_stats, null_sub.copy(), k, 1.0)
    elapsed = time.perf_counter() - t0
    n_steps = int(np.sum(pvals < 1.0))
    print(f"    Final pass (max_pval=1.0): {elapsed:.2f}s ({n_steps} steps before p=1)")

    # --- Summary table ---
    print("\n" + "=" * 60)
    print("TIMING SUMMARY")
    print("=" * 60)
    print(f"{'Stage':<25s} {'Time (s)':>10s} {'% Total':>10s}")
    print("-" * 47)
    for stage in ["build_count_matrix", "compute_test_stats", "bootstrap", "k_fwer"]:
        t = timings[stage]
        pct = t / timings["total"] * 100
        print(f"  {stage:<23s} {t:>9.2f}  {pct:>8.1f}%")
    print("-" * 47)
    print(f"  {'TOTAL':<23s} {timings['total']:>9.2f}")
    print()
    print(f"Results: {n_sig} significant sequences (k={k_final})")
    print()


if __name__ == "__main__":
    if not DATA_DIR.exists():
        print(f"Data not found at {DATA_DIR}")
        print("Run: git clone https://github.com/dbkastner/CBAS.git igor_cbas")
        raise SystemExit(1)

    parser = argparse.ArgumentParser(description="Benchmark CBAS pipeline")
    parser.add_argument("--quick", action="store_true",
                        help="Use reduced parameters for a fast check (~2s)")
    args = parser.parse_args()
    run_benchmark(quick=args.quick)
