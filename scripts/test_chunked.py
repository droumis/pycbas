"""
Verify that find_k_fwer_chunked produces identical results to find_k_fwer,
then benchmark memory and time.
"""

import time
import numpy as np
from pathlib import Path

from pycbas import (
    CBASParams,
    load_subject_data,
    build_count_matrix,
    compute_test_stats,
    bootstrap_test_stats,
    find_k_fwer,
    find_k_fwer_chunked,
)

ROOT_DIR = Path(__file__).parent.parent
DATA_DIR = ROOT_DIR / "igor_cbas" / "data"


def load_rats(n_ctrl_max=46, n_les_max=39):
    ctrl_data, les_data = [], []
    for f in sorted(DATA_DIR.glob("*.txt")):
        name = f.stem
        if "Control" in name:
            ctrl_data.append(load_subject_data(f))
        elif "Lesion" in name:
            les_data.append(load_subject_data(f))
    ctrl_data = ctrl_data[:n_ctrl_max]
    les_data = les_data[:n_les_max]
    subjects_data = ctrl_data + les_data
    group_labels = np.array([0] * len(ctrl_data) + [1] * len(les_data))
    return subjects_data, group_labels


def main():
    print("=" * 60)
    print("TEST: Chunked vs standard pipeline")
    print("=" * 60)

    subjects_data, group_labels = load_rats()
    params = CBASParams(num_arms=6, seq_len_max=6, criterion=800, resample_number=10000)

    group_indices = [
        np.where(group_labels == 0)[0],
        np.where(group_labels == 1)[0],
    ]

    print(f"\nBuilding count matrix...")
    t0 = time.perf_counter()
    sequences, count_matrix = build_count_matrix(subjects_data, params)
    n_seq = len(sequences)
    print(f"  {n_seq:,} sequences, {time.perf_counter() - t0:.1f}s")

    print(f"\nComputing test stats...")
    test_stats = compute_test_stats(count_matrix, group_indices)
    n_valid = int(np.sum(~np.isnan(test_stats)))
    print(f"  {n_valid:,} valid test stats")

    null_sub_mem = 10000 * n_valid * 8 / 1024**3
    full_null_mem = 10000 * n_seq * 2 * 8 / 1024**3
    print(f"\n  Full null matrix would be: {full_null_mem:.2f} GB")
    print(f"  null_sub (sorted valid only): {null_sub_mem:.2f} GB")
    print(f"  Standard peak (both): {full_null_mem + null_sub_mem:.2f} GB")
    print(f"  Chunked peak (null_sub only): {null_sub_mem:.2f} GB")

    # --- Standard approach ---
    print(f"\n--- Standard approach ---")
    t0 = time.perf_counter()
    null_matrix, _ = bootstrap_test_stats(count_matrix, group_indices, params)
    t_boot = time.perf_counter() - t0
    print(f"  Bootstrap: {t_boot:.1f}s")

    t0 = time.perf_counter()
    g_values_std, k_std = find_k_fwer(test_stats, null_matrix, params.alpha, params.gamma)
    t_stepdown = time.perf_counter() - t0
    n_sig_std = sum(1 for i in range(n_seq) if
                    (not np.isnan(g_values_std[i*2]) and g_values_std[i*2] < 0.5) or
                    (not np.isnan(g_values_std[i*2+1]) and g_values_std[i*2+1] < 0.5))
    print(f"  Step-down: {t_stepdown:.1f}s")
    print(f"  Result: {n_sig_std} significant, k={k_std}")
    print(f"  Total: {t_boot + t_stepdown:.1f}s")

    del null_matrix

    # --- Chunked approach ---
    print(f"\n--- Chunked approach (chunk_size=500) ---")
    t0 = time.perf_counter()
    g_values_chk, k_chk = find_k_fwer_chunked(
        test_stats, count_matrix, group_indices, params, chunk_size=500
    )
    t_chunked = time.perf_counter() - t0
    n_sig_chk = sum(1 for i in range(n_seq) if
                    (not np.isnan(g_values_chk[i*2]) and g_values_chk[i*2] < 0.5) or
                    (not np.isnan(g_values_chk[i*2+1]) and g_values_chk[i*2+1] < 0.5))
    print(f"  Total: {t_chunked:.1f}s")
    print(f"  Result: {n_sig_chk} significant, k={k_chk}")

    # Verify match
    both_valid = ~np.isnan(g_values_std) & ~np.isnan(g_values_chk)
    n_both = int(both_valid.sum())
    if np.allclose(g_values_std[both_valid], g_values_chk[both_valid], rtol=1e-10):
        print(f"  MATCH: all {n_both:,} valid p-values identical")
    else:
        diffs = np.abs(g_values_std[both_valid] - g_values_chk[both_valid])
        max_diff = np.max(diffs)
        n_differ = int(np.sum(diffs > 1e-10))
        print(f"  MISMATCH: {n_differ}/{n_both} differ, max diff = {max_diff:.2e}")

    print("\n" + "=" * 60)


if __name__ == "__main__":
    main()
