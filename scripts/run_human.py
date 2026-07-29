"""
Run CBAS on the human two-step task dataset (correlative mode).

Correlates sequence usage with CBIT scores (compulsivity measure).
6 choice values × 2 reward states = 12 possible symbols per position.

Paper params: num_arms=6, seq_len_max=4, criterion=400, M=10,000
Paper result: 31/408 significant sequences (Fig 1c middle panel)

Usage:
    pixi run human             # paper params
    pixi run human-quick       # reduced for fast check
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
    compute_test_stats_correlative,
    bootstrap_test_stats_correlative,
    find_k_fwer,
)
from results_io import save_results_json, compute_significance_summary

ROOT_DIR = Path(__file__).parent.parent
DATA_DIR = ROOT_DIR / "data" / "humans"
RESULTS_DIR = ROOT_DIR / "results" / "humans"
FIG_DIR = RESULTS_DIR / "figures"


def load_humans():
    """Load human data and CBIT scores from info file."""
    info_path = DATA_DIR / "humanInfo.txt"
    info = {}
    with open(info_path) as f:
        for line in f:
            parts = line.strip().split(",")
            info[int(parts[0])] = float(parts[1])

    subjects_data = []
    covariate = []
    for subj_id in sorted(info.keys()):
        fpath = DATA_DIR / f"subject{subj_id}.txt"
        if fpath.exists():
            subjects_data.append(load_subject_data(fpath))
            covariate.append(info[subj_id])

    return subjects_data, np.array(covariate)


def decode_human_symbol(sym, num_arms=6):
    """Decode human symbol into choice description."""
    choice = sym % num_arms
    rewarded = sym // num_arms
    choice_names = ["L1", "R1", "L2", "R2", "NC1", "NC2"]
    name = choice_names[choice] if choice < len(choice_names) else f"c{choice}"
    if rewarded:
        name = name.upper()
    return name


def decode_human_sequence(seq, num_arms=6):
    return " ".join(decode_human_symbol(s, num_arms) for s in seq)


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
                             title="Human CBAS: Two-Step Task × CBIT (Correlative)")
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
            d = "pos_corr" if pos_g <= neg_g else "neg_corr"
        elif not np.isnan(pos_g):
            best_g = pos_g
            d = "pos_corr"
        elif not np.isnan(neg_g):
            best_g = neg_g
            d = "neg_corr"
        if not np.isnan(best_g) and best_g > 0:
            neg_log_g[i] = -np.log10(best_g)
        directions.append(d)

    threshold = -np.log10(alpha)
    directions_arr = np.array(directions)
    valid = ~np.isnan(neg_log_g)
    sig_mask = neg_log_g > threshold
    n_pos = int(np.sum((directions_arr == "pos_corr") & sig_mask & valid))
    n_neg = int(np.sum((directions_arr == "neg_corr") & sig_mask & valid))

    fig, ax = direction_counts(
        n_pos, n_neg,
        "Positive corr\n(↑ CBIT → ↑ usage)",
        "Negative corr\n(↑ CBIT → ↓ usage)",
        colors=["#cc3300", "#0066cc"],
        title="Significant Sequences by Correlation Direction",
        text_offset=0.5)
    fig.savefig(FIG_DIR / "direction_counts.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Null vs observed ---
    fig, ax = null_vs_observed(test_stats, null_row_maxes)
    fig.savefig(FIG_DIR / "null_vs_observed.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Sequence space ---
    fig, ax = sequence_space(seq_lengths, num_arms=6)
    fig.savefig(FIG_DIR / "sequence_space.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- g-value distribution ---
    fig, ax = gvalue_distribution(g_values, n_seq, alpha=alpha)
    fig.savefig(FIG_DIR / "gvalue_dist.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Figures saved to: {FIG_DIR}/")


def write_report(data, timings):
    """Write markdown validation report."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    g_values = data["g_values"]
    seq_lengths = data["seq_lengths"]
    seq_strs = data["seq_strs"]
    n_subjects = int(data["n_subjects"][0])
    k_final = int(data["k_final"][0])
    n_seq = len(seq_lengths)

    alpha = 0.5
    n_sig = 0
    n_pos = 0
    n_neg = 0
    sig_seqs = []

    for i in range(n_seq):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        is_sig = False
        direction = ""
        best_g = np.nan

        if not np.isnan(pos_g) and pos_g < alpha:
            is_sig = True
            direction = "pos_corr"
            best_g = pos_g
            n_pos += 1
        if not np.isnan(neg_g) and neg_g < alpha:
            is_sig = True
            direction = "neg_corr"
            best_g = neg_g
            n_neg += 1

        if is_sig:
            n_sig += 1
            sig_seqs.append((seq_strs[i], direction, best_g, int(seq_lengths[i])))

    sig_seqs.sort(key=lambda x: x[2])

    report = f"""# Human CBAS Validation Report (Correlative Mode)

## Summary

| | pycbas | Paper (Kastner et al.) |
|---|---|---|
| Subjects | {n_subjects} | 1,413 |
| Max seq length | {int(data['params_seq_len_max'][0])} | 4 |
| Criterion | {int(data['params_criterion'][0])} | 400 |
| Resamples | {int(data['params_resample_number'][0])} | 10,000 |
| Sequences evaluated | {n_seq:,} | 408 |
| Significant | {n_sig} ({n_sig/n_seq*100:.1f}%) | 31 (7.6%) |
| Positive correlation | {n_pos} | not separately reported |
| Negative correlation | {n_neg} | not separately reported |
| k (k-FWER) | {k_final} | not reported |
| Runtime | {timings['total']:.1f}s | not reported |

## Notes

- **Mode:** Correlative — tests Pearson correlation between each sequence's usage
  count across subjects and each subject's CBIT score (a compulsivity measure).
- **Symbol encoding:** choice + reward × 6. Choices: 0=L1, 1=R1, 2=L2, 3=R2,
  4=no-choice-stage1, 5=no-choice-stage2. UPPERCASE = rewarded.
- **Interpretation:** Positive correlation means higher CBIT (more compulsive)
  subjects use that sequence more. Negative means less.

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

### Significant Sequences by Correlation Direction
![Direction Counts](figures/direction_counts.png)

### Null Distribution vs Observed
![Null vs Observed](figures/null_vs_observed.png)

### Sequence Space
![Sequence Space](figures/sequence_space.png)

### g-value Distribution
![g-value Distribution](figures/gvalue_dist.png)

## Top Significant Sequences

| Sequence | Direction | ζ-value | Decoded |
|---|---|---|---|
"""
    for seq_str, direction, gval, slen in sig_seqs[:25]:
        decoded = decode_human_sequence(tuple(int(x) for x in seq_str.split("-")))
        dir_label = "+" if direction == "pos_corr" else "−"
        report += f"| {seq_str} | {dir_label} | {gval:.4f} | {decoded} |\n"

    report_path = RESULTS_DIR / "validation_report.md"
    with open(report_path, "w") as f:
        f.write(report)
    print(f"Report saved to: {report_path}")


def run_analysis(quick=False):
    print("=" * 60)
    print("CBAS — Human Two-Step Task (Correlative with CBIT)")
    print("=" * 60)

    subjects_data, covariate = load_humans()
    n_subjects = len(subjects_data)
    print(f"\nData: {n_subjects} subjects")
    print(f"CBIT scores: mean={covariate.mean():.3f}, "
          f"std={covariate.std():.3f}, range=[{covariate.min():.2f}, {covariate.max():.2f}]")

    if quick:
        params = CBASParams(num_arms=6, seq_len_max=2, criterion=400, resample_number=1000)
    else:
        params = CBASParams(num_arms=6, seq_len_max=4, criterion=400, resample_number=10000)
    print(f"Params: num_arms={params.num_arms}, seq_len_max={params.seq_len_max}, "
          f"criterion={params.criterion}, M={params.resample_number}")

    timings = {}

    t0 = time.perf_counter()
    sequences, count_matrix = build_count_matrix(subjects_data, params, contingency=1)
    timings["build_count_matrix"] = time.perf_counter() - t0
    n_seq = len(sequences)
    print(f"\n[{timings['build_count_matrix']:.2f}s] Count matrix: "
          f"{n_subjects} x {n_seq}")

    t0 = time.perf_counter()
    test_stats = compute_test_stats_correlative(count_matrix, covariate)
    timings["compute_test_stats"] = time.perf_counter() - t0
    n_valid = int(np.sum(~np.isnan(test_stats)))
    print(f"[{timings['compute_test_stats']:.2f}s] Test stats: {n_valid} valid")

    t0 = time.perf_counter()
    null_matrix, null_directions = bootstrap_test_stats_correlative(count_matrix, covariate, params)
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
    null_row_maxes = np.nanmax(null_matrix, axis=1)

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
    )

    # Save structured results JSON
    results_json = {
        "dataset": "humans",
        "mode": "correlative",
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
            "positive_direction": "Positive correlation (↑ CBIT → ↑ usage)",
            "negative_direction": "Negative correlation (↑ CBIT → ↓ usage)",
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
    parser = argparse.ArgumentParser(description="Run CBAS on human two-step task data")
    parser.add_argument("--quick", action="store_true",
                        help="Reduced params (seq_len_max=2, M=1000)")
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
        run_analysis(quick=args.quick)


if __name__ == "__main__":
    main()
