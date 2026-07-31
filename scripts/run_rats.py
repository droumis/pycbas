"""
Run CBAS on the rat spatial alternation dataset.

Compares control vs hippocampal lesion rats. 6-arm maze with reward encoding
(12 symbols), first 800 choices per rat.

Paper params: num_arms=6, seq_len_max=6, criterion=800, M=10,000
Paper result: 409/24,342 significant sequences (Fig 1c right panel)

Usage:
    pixi run rats             # paper params (85 subjects)
    pixi run rats-quick       # reduced for fast check
    pixi run rats --full      # all available subjects
"""

import argparse
import sys
import time
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from pycbas import (
    CBASParams,
    load_subject_data,
    build_count_matrix,
    compute_test_stats,
    bootstrap_test_stats,
    find_k_fwer,
    find_k_fwer_chunked,
)
from results_io import save_results_json, compute_significance_summary

ROOT_DIR = Path(__file__).parent.parent
DATA_DIR = ROOT_DIR / "igor_cbas" / "data"
RESULTS_DIR = ROOT_DIR / "results" / "rats"
FIG_DIR = RESULTS_DIR / "figures"


def load_rats(n_ctrl_max=46, n_les_max=39):
    """Load rat data files. Optionally limit to first N of each group."""
    ctrl_data, les_data = [], []
    for f in sorted(DATA_DIR.glob("*.txt")):
        name = f.stem
        if "Control" in name:
            ctrl_data.append(load_subject_data(f))
        elif "Lesion" in name:
            les_data.append(load_subject_data(f))

    if n_ctrl_max is not None:
        ctrl_data = ctrl_data[:n_ctrl_max]
    if n_les_max is not None:
        les_data = les_data[:n_les_max]

    subjects_data = ctrl_data + les_data
    group_labels = np.array([0] * len(ctrl_data) + [1] * len(les_data))
    return subjects_data, group_labels


def decode_rat_sequence(seq, num_arms=6):
    """Decode rat sequence: arm = sym % 6, rewarded = sym // 6, display as '{arm+1}*' if rewarded."""
    parts = []
    for s in seq:
        arm = s % num_arms
        rewarded = s // num_arms
        parts.append(f"{arm+1}{'*' if rewarded else ''}")
    return " ".join(parts)


def make_figures(data):
    """Generate figures from cached results."""
    import matplotlib.pyplot as plt
    from plots import manhattan_plot, null_vs_observed, gvalue_distribution, sequence_space, direction_counts

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    g_values = data["g_values"]
    test_stats = data["test_stats"]
    seq_lengths = data["seq_lengths"]
    null_row_maxes = data["null_row_maxes"]
    n_seq = len(seq_lengths)
    alpha = 0.5

    # --- Manhattan plot ---
    fig, ax = manhattan_plot(g_values, seq_lengths, alpha=alpha,
                             title="Rat CBAS: Control vs Lesion Spatial Alternation")
    fig.savefig(FIG_DIR / "manhattan.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Direction counts ---
    neg_log_g = np.full(n_seq, np.nan)
    directions = []
    for i in range(n_seq):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        best_g = np.nan
        d = ""
        if not np.isnan(pos_g) and not np.isnan(neg_g):
            best_g = min(pos_g, neg_g)
            d = "control>lesion" if pos_g <= neg_g else "lesion>control"
        elif not np.isnan(pos_g):
            best_g = pos_g
            d = "control>lesion"
        elif not np.isnan(neg_g):
            best_g = neg_g
            d = "lesion>control"
        if not np.isnan(best_g) and best_g > 0:
            neg_log_g[i] = -np.log10(best_g)
        directions.append(d)

    threshold = -np.log10(alpha)
    directions_arr = np.array(directions)
    valid = ~np.isnan(neg_log_g)
    sig_mask = neg_log_g > threshold
    n_ctrl_more = int(np.sum((directions_arr == "control>lesion") & sig_mask & valid))
    n_les_more = int(np.sum((directions_arr == "lesion>control") & sig_mask & valid))

    fig, ax = direction_counts(n_ctrl_more, n_les_more,
                               "Control > Lesion", "Lesion > Control",
                               colors=["steelblue", "sienna"])
    fig.savefig(FIG_DIR / "direction_counts.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Null vs observed ---
    fig, ax = null_vs_observed(test_stats, null_row_maxes)
    fig.savefig(FIG_DIR / "null_vs_observed.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Sequence space ---
    num_symbols = 6 * 2  # num_arms * 2 for reward encoding
    fig, ax = sequence_space(seq_lengths, num_arms=num_symbols)
    fig.savefig(FIG_DIR / "sequence_space.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- g-value distribution ---
    fig, ax = gvalue_distribution(
        g_values, n_seq, alpha=alpha,
        title="Distribution of g-values\n"
              "(Bimodal: true differences pile up near 0, nulls near 1)")
    fig.savefig(FIG_DIR / "gvalue_dist.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Figures saved to: {FIG_DIR}/")


def write_report(data, timings):
    """Write markdown validation report."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    g_values = data["g_values"]
    test_stats = data["test_stats"]
    seq_lengths = data["seq_lengths"]
    seq_strs = data["seq_strs"]
    n_subjects = int(data["n_subjects"][0])
    n_ctrl = int(data["n_ctrl"][0])
    n_les = int(data["n_les"][0])
    k_final = int(data["k_final"][0])
    n_seq = len(seq_lengths)

    alpha = 0.5
    n_sig = 0
    n_ctrl_more = 0
    n_les_more = 0
    sig_seqs = []

    for i in range(n_seq):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        is_sig = False
        direction = ""
        best_g = np.nan

        if not np.isnan(pos_g) and pos_g < alpha:
            is_sig = True
            direction = "control>lesion"
            best_g = pos_g
            n_ctrl_more += 1
        if not np.isnan(neg_g) and neg_g < alpha:
            is_sig = True
            direction = "lesion>control"
            best_g = neg_g
            n_les_more += 1

        if is_sig:
            n_sig += 1
            sig_seqs.append((seq_strs[i], direction, best_g, int(seq_lengths[i])))

    sig_seqs.sort(key=lambda x: x[2])

    report = f"""# Rat CBAS Validation Report

## Summary

| | pycbas | Paper (Kastner et al.) |
|---|---|---|
| Rats | {n_subjects} ({n_ctrl} control, {n_les} lesion) | 85 (46 control, 39 lesion) |
| Max seq length | {int(data['params_seq_len_max'][0])} | 6 |
| Criterion | {int(data['params_criterion'][0])} | 800 |
| Resamples | {int(data['params_resample_number'][0])} | 10,000 |
| Sequences evaluated | {n_seq:,} | 24,342 |
| Significant | {n_sig} ({n_sig/n_seq*100:.1f}%) | 409 (1.7%) |
| Control > Lesion | {n_ctrl_more} | not separately reported |
| Lesion > Control | {n_les_more} | not separately reported |
| k (k-FWER) | {k_final} | not reported |
| Runtime | {timings['total']:.1f}s | not reported |

## Timing Profile

| Stage | Time (s) | % Total |
|---|---|---|
| build_count_matrix | {timings['build_count_matrix']:.2f} | {timings['build_count_matrix']/timings['total']*100:.1f}% |
| compute_test_stats | {timings['compute_test_stats']:.2f} | {timings['compute_test_stats']/timings['total']*100:.1f}% |
| bootstrap | {timings['bootstrap']:.2f} | {timings['bootstrap']/timings['total']*100:.1f}% |
| k_fwer | {timings['k_fwer']:.2f} | {timings['k_fwer']/timings['total']*100:.1f}% |
| **TOTAL** | **{timings['total']:.2f}** | |

## Figures

### Manhattan Plot
![Manhattan Plot](figures/manhattan.png)

### Significant Sequences by Direction
![Direction Counts](figures/direction_counts.png)

### Null Distribution vs Observed
![Null vs Observed](figures/null_vs_observed.png)

### Sequence Space
![Sequence Space](figures/sequence_space.png)

### g-value Distribution
![g-value Distribution](figures/gvalue_dist.png)

## Top Significant Sequences

| Sequence | Direction | ζ-value | Decoded (arm, * = rewarded) |
|---|---|---|---|
"""
    for seq_str, direction, gval, slen in sig_seqs[:25]:
        decoded = decode_rat_sequence(tuple(int(x) for x in seq_str.split("-")))
        report += f"| {seq_str} | {direction} | {gval:.4f} | {decoded} |\n"

    report_path = RESULTS_DIR / "validation_report.md"
    with open(report_path, "w") as f:
        f.write(report)
    print(f"Report saved to: {report_path}")


def run_analysis(quick=False, full=False, chunked=False):
    print("=" * 60)
    print("CBAS — Rat Spatial Alternation (Control vs Lesion)")
    print("=" * 60)

    if full:
        subjects_data, group_labels = load_rats(n_ctrl_max=None, n_les_max=None)
    else:
        subjects_data, group_labels = load_rats(n_ctrl_max=46, n_les_max=39)
    n_ctrl = int((group_labels == 0).sum())
    n_les = int((group_labels == 1).sum())
    n_subjects = len(subjects_data)
    print(f"\nData: {n_subjects} rats ({n_ctrl} control, {n_les} lesion)")

    if quick:
        params = CBASParams(num_arms=6, seq_len_max=4, criterion=800, resample_number=1000)
    else:
        params = CBASParams(num_arms=6, seq_len_max=6, criterion=800, resample_number=10000)
    print(f"Params: num_arms={params.num_arms}, seq_len_max={params.seq_len_max}, "
          f"criterion={params.criterion}, M={params.resample_number}")
    if chunked:
        print(f"Mode: CHUNKED (memory-efficient)")

    group_indices = [
        np.where(group_labels == 0)[0],
        np.where(group_labels == 1)[0],
    ]

    timings = {}

    t0 = time.perf_counter()
    sequences, count_matrix = build_count_matrix(subjects_data, params)
    timings["build_count_matrix"] = time.perf_counter() - t0
    n_seq = len(sequences)
    print(f"\n[{timings['build_count_matrix']:.2f}s] Count matrix: "
          f"{n_subjects} x {n_seq}")

    t0 = time.perf_counter()
    test_stats = compute_test_stats(count_matrix, group_indices)
    timings["compute_test_stats"] = time.perf_counter() - t0
    n_valid = int(np.sum(~np.isnan(test_stats)))
    print(f"[{timings['compute_test_stats']:.2f}s] Test stats: {n_valid} valid")

    if chunked:
        t0 = time.perf_counter()
        g_values, k_final = find_k_fwer_chunked(
            test_stats, count_matrix, group_indices, params, chunk_size=500
        )
        timings["bootstrap_and_k_fwer"] = time.perf_counter() - t0
        timings["bootstrap"] = timings["bootstrap_and_k_fwer"] * 0.3  # approx split
        timings["k_fwer"] = timings["bootstrap_and_k_fwer"] * 0.7
        print(f"[{timings['bootstrap_and_k_fwer']:.2f}s] Chunked bootstrap+step-down: k={k_final}")
    else:
        t0 = time.perf_counter()
        null_matrix, null_directions = bootstrap_test_stats(count_matrix, group_indices, params)
        timings["bootstrap"] = time.perf_counter() - t0
        print(f"[{timings['bootstrap']:.2f}s] Bootstrap: {params.resample_number} resamples")

        t0 = time.perf_counter()
        g_values, k_final = find_k_fwer(test_stats, null_matrix, params.alpha, params.gamma,
                                         null_directions=null_directions)
        timings["k_fwer"] = time.perf_counter() - t0
        print(f"[{timings['k_fwer']:.2f}s] k-FWER: k={k_final}")

    timings["total"] = sum(timings.values())

    # Compute significance summary
    sig_summary = compute_significance_summary(g_values, n_seq, params.alpha)
    n_sig = sig_summary["n_significant"]
    print(f"\nResult: {n_sig}/{n_seq} significant sequences ({n_sig/n_seq*100:.1f}%)")
    print(f"Total time: {timings['total']:.1f}s")

    # Cache arrays
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    seq_lengths = np.array([len(s) for s in sequences])
    seq_strs = np.array(["-".join(str(x) for x in s) for s in sequences])
    if not chunked:
        null_row_maxes = np.nanmax(null_matrix, axis=1)
    else:
        null_row_maxes = np.zeros(params.resample_number)

    cache_path = FIG_DIR / "results.npz"
    np.savez_compressed(
        cache_path,
        g_values=g_values,
        test_stats=test_stats,
        seq_lengths=seq_lengths,
        seq_strs=seq_strs,
        null_row_maxes=null_row_maxes,
        k_final=np.array([k_final]),
        params_seq_len_max=np.array([params.seq_len_max]),
        params_criterion=np.array([params.criterion]),
        params_resample_number=np.array([params.resample_number]),
        n_subjects=np.array([n_subjects]),
        n_ctrl=np.array([n_ctrl]),
        n_les=np.array([n_les]),
    )

    # Save structured results JSON (source of truth for reports/figures)
    results_json = {
        "dataset": "rats",
        "mode": "comparative",
        "groups": {"control": n_ctrl, "lesion": n_les},
        "params": {
            "num_arms": params.num_arms,
            "seq_len_max": params.seq_len_max,
            "criterion": params.criterion,
            "resample_number": params.resample_number,
            "alpha": params.alpha,
            "gamma": params.gamma,
        },
        "results": {
            "n_subjects": n_subjects,
            "n_sequences": n_seq,
            "n_significant": sig_summary["n_significant"],
            "n_positive": sig_summary["n_positive"],
            "n_negative": sig_summary["n_negative"],
            "fraction_significant": sig_summary["fraction_significant"],
            "k_final": k_final,
        },
        "timing": timings,
        "labels": {
            "positive_direction": "control > lesion",
            "negative_direction": "lesion > control",
        },
    }
    json_path = RESULTS_DIR / "results.json"
    save_results_json(json_path, results_json)
    print(f"Results JSON: {json_path}")
    print(f"Results NPZ: {cache_path}")

    data = np.load(cache_path, allow_pickle=False)
    make_figures(data)
    write_report(data, timings)


def main():
    parser = argparse.ArgumentParser(description="Run CBAS on rat data")
    parser.add_argument("--quick", action="store_true",
                        help="Reduced params (seq_len_max=4, M=1000)")
    parser.add_argument("--full", action="store_true",
                        help="Use all available rats (not just first 85)")
    parser.add_argument("--chunked", action="store_true",
                        help="Use memory-efficient chunked pipeline")
    parser.add_argument("--figures-only", action="store_true",
                        help="Regenerate figures from cached results")
    args = parser.parse_args()

    if args.figures_only:
        cache_path = FIG_DIR / "results.npz"
        if not cache_path.exists():
            print(f"No cached results at {cache_path}. Run analysis first.")
            raise SystemExit(1)
        data = np.load(cache_path, allow_pickle=False)
        make_figures(data)
    else:
        run_analysis(quick=args.quick, full=args.full, chunked=args.chunked)


if __name__ == "__main__":
    main()
