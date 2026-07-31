"""
Run CBAS on the Noel et al. 2025 IBL autism model mice.

Compares WT (C57BL/6J) to each of three autism models: Cntnap2, Fmr1, Shank3B.
4-symbol encoding (choice × stimulus_side), reward is deterministic (2AFC).

Usage:
    pixi run noel-ibl                  # all 3 comparisons, L=6
    pixi run noel-ibl --seq-len 8      # L=8 (~5 GB, ~30s)
    pixi run noel-ibl --model Cntnap2  # single comparison
    pixi run noel-ibl --quick          # reduced M for fast check
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
    find_k_fwer_chunked,
    estimate_resources,
    print_resource_estimate,
)
from results_io import save_results_json, compute_significance_summary

ROOT_DIR = Path(__file__).parent.parent
DATA_DIR = ROOT_DIR / "data" / "noel_ibl_mice"
RESULTS_DIR = ROOT_DIR / "results" / "noel_ibl_mice"
FIG_DIR = RESULTS_DIR / "figures"


def load_noel_ibl():
    """Load Noel IBL data, returning subjects grouped by genotype."""
    info_path = DATA_DIR / "ibl_mice_info.txt"
    subjects_by_genotype = {"WT": [], "Cntnap2": [], "Fmr1": [], "Shank3B": []}

    with open(info_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split(",")
            if len(parts) < 3:
                continue
            animal_id, genotype = parts[0], parts[1]
            filepath = DATA_DIR / f"{animal_id}.txt"
            if filepath.exists():
                subjects_by_genotype[genotype].append(load_subject_data(filepath))

    return subjects_by_genotype


def run_comparison(wt_subjects, model_subjects, model_name, params):
    """Run one WT vs model comparison."""
    subjects_data = wt_subjects + model_subjects
    group_labels = np.array([0] * len(wt_subjects) + [1] * len(model_subjects))
    group_indices = [np.where(group_labels == 0)[0], np.where(group_labels == 1)[0]]

    timings = {}

    t0 = time.perf_counter()
    sequences, count_matrix = build_count_matrix(
        subjects_data, params, contingency=None, encode_reward=False
    )
    timings["build_count_matrix"] = time.perf_counter() - t0
    n_seq = len(sequences)

    t0 = time.perf_counter()
    test_stats = compute_test_stats(count_matrix, group_indices)
    timings["compute_test_stats"] = time.perf_counter() - t0
    n_valid = int(np.sum(~np.isnan(test_stats)))

    t0 = time.perf_counter()
    g_values, k_final = find_k_fwer_chunked(
        test_stats, count_matrix, group_indices, params, chunk_size=500
    )
    timings["bootstrap_and_stepdown"] = time.perf_counter() - t0
    timings["total"] = sum(timings.values())

    sig_summary = compute_significance_summary(g_values, n_seq, params.alpha)
    n_sig = sig_summary["n_significant"]

    print(f"  [{timings['build_count_matrix']:.1f}s] Count matrix: "
          f"{len(subjects_data)} × {n_seq} ({n_valid} valid)")
    print(f"  [{timings['bootstrap_and_stepdown']:.1f}s] Bootstrap + step-down: k={k_final}")
    print(f"  Result: {n_sig}/{n_seq} significant ({n_sig/n_seq*100:.1f}%)")
    print(f"    WT > {model_name}: {sig_summary['n_positive']}, "
          f"{model_name} > WT: {sig_summary['n_negative']}")
    print(f"  Total: {timings['total']:.1f}s")

    return {
        "sequences": sequences,
        "test_stats": test_stats,
        "g_values": g_values,
        "k_final": k_final,
        "n_seq": n_seq,
        "n_valid": n_valid,
        "sig_summary": sig_summary,
        "timings": timings,
        "n_wt": len(wt_subjects),
        "n_model": len(model_subjects),
    }


def make_figures(results_dict, model_name, params):
    """Generate manhattan plot for a comparison."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT_DIR / "plots" / "matplotlib"))
    from manhattan import manhattan_plot

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    g_values = results_dict["g_values"]
    seq_lengths = np.array([len(s) for s in results_dict["sequences"]])
    n_seq = results_dict["n_seq"]

    fig, ax = manhattan_plot(
        g_values, seq_lengths, alpha=params.alpha,
        title=f"Noel IBL CBAS: WT vs {model_name} (L={params.seq_len_max})"
    )
    fig.savefig(FIG_DIR / f"manhattan_wt_vs_{model_name.lower()}.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)


def run_analysis(models=None, seq_len=6, quick=False):
    print("=" * 60)
    print("CBAS — Noel et al. 2025 IBL Autism Models")
    print("=" * 60)

    subjects_by_genotype = load_noel_ibl()
    for g, subjs in subjects_by_genotype.items():
        n_trials = sum(len(s) for s in subjs)
        print(f"  {g:8s}: {len(subjs):2d} mice, {n_trials:>7,} trials")

    if models is None:
        models = ["Cntnap2", "Fmr1", "Shank3B"]

    M = 1000 if quick else 10000
    params = CBASParams(
        num_arms=4, seq_len_max=seq_len, criterion=99999,
        resample_number=M
    )
    print(f"\nParams: A=4 (no reward encoding), L={seq_len}, "
          f"criterion=all, M={M:,}")

    # Show resource estimate for first comparison
    est = estimate_resources(
        num_arms=4, seq_len_max=seq_len, encode_reward=False,
        resample_number=M
    )
    print(f"Worst-case: {est['total_sequences']:,} sequences, "
          f"{est['memory_chunked_gb']:.1f} GB (chunked)")
    print()

    wt_subjects = subjects_by_genotype["WT"]
    all_results = {}

    for model_name in models:
        model_subjects = subjects_by_genotype[model_name]
        print(f"\n--- WT ({len(wt_subjects)}) vs {model_name} ({len(model_subjects)}) ---")

        result = run_comparison(wt_subjects, model_subjects, model_name, params)
        all_results[model_name] = result

        # Save results JSON
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        results_json = {
            "dataset": "noel_ibl_mice",
            "mode": "comparative",
            "comparison": f"WT_vs_{model_name}",
            "groups": {"WT": result["n_wt"], model_name: result["n_model"]},
            "params": {
                "num_arms": params.num_arms,
                "seq_len_max": params.seq_len_max,
                "encode_reward": False,
                "resample_number": params.resample_number,
                "alpha": params.alpha,
                "gamma": params.gamma,
            },
            "results": {
                "n_subjects": result["n_wt"] + result["n_model"],
                "n_sequences": result["n_seq"],
                "n_valid": result["n_valid"],
                "n_significant": result["sig_summary"]["n_significant"],
                "n_positive": result["sig_summary"]["n_positive"],
                "n_negative": result["sig_summary"]["n_negative"],
                "fraction_significant": result["sig_summary"]["fraction_significant"],
                "k_final": result["k_final"],
            },
            "timing": result["timings"],
            "labels": {
                "positive_direction": f"WT > {model_name}",
                "negative_direction": f"{model_name} > WT",
            },
        }
        json_path = RESULTS_DIR / f"results_wt_vs_{model_name.lower()}.json"
        save_results_json(json_path, results_json)

        try:
            make_figures(result, model_name, params)
        except Exception as e:
            print(f"  (figures skipped: {e})")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"{'Comparison':<20} {'Sequences':>10} {'Significant':>12} {'k':>4} {'Time':>8}")
    print("-" * 60)
    for model_name, result in all_results.items():
        n_sig = result["sig_summary"]["n_significant"]
        print(f"WT vs {model_name:<13} {result['n_seq']:>10,} "
              f"{n_sig:>8} ({n_sig/result['n_seq']*100:.1f}%) "
              f"{result['k_final']:>4} {result['timings']['total']:>7.1f}s")


def main():
    parser = argparse.ArgumentParser(description="Run CBAS on Noel IBL mice")
    parser.add_argument("--quick", action="store_true",
                        help="Reduced M=1000 for fast check")
    parser.add_argument("--seq-len", type=int, default=6,
                        help="Max sequence length (default 6)")
    parser.add_argument("--model", type=str, choices=["Cntnap2", "Fmr1", "Shank3B"],
                        help="Run single comparison (default: all three)")
    args = parser.parse_args()

    models = [args.model] if args.model else None
    run_analysis(models=models, seq_len=args.seq_len, quick=args.quick)


if __name__ == "__main__":
    main()
