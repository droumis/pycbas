"""
CBAS on Noel et al. 2025 IBL mice — corrected analysis.

Restricts to 0%-contrast trials in biased blocks, truncated to a common
per-animal criterion. This isolates prior-driven decisions, removes the
sensory-accuracy confound, and eliminates trial-count bias.

Also runs a virtual-twin surrogate control: for each animal, generates a
synthetic choice stream matched on marginal accuracy and side bias, then
runs CBAS on the surrogates. If the real result reproduces under surrogates,
there's nothing beyond first-order statistics.

Usage:
    pixi run noel-ibl-zero           # default L=6, criterion=600
    pixi run noel-ibl-zero --seq-len 4
    pixi run noel-ibl-zero --criterion 400
    pixi run noel-ibl-zero --skip-surrogates
"""

import argparse
import sys
import time
import json
import numpy as np
import scipy.io as sio
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent))
from pycbas import (
    CBASParams,
    build_count_matrix,
    compute_test_stats,
    find_k_fwer_chunked,
    enumerate_sequences,
)
from results_io import save_results_json, compute_significance_summary

ROOT_DIR = Path(__file__).parent.parent
DATA_DIR = ROOT_DIR / "data" / "noel_ibl_mice"
RESULTS_DIR = ROOT_DIR / "results" / "noel_ibl_mice_zero_contrast"
FIG_DIR = RESULTS_DIR / "figures"

SYM = {0: "L✓", 1: "L✗", 2: "R✗", 3: "R✓"}
IS_CORRECT = {0: True, 1: False, 2: False, 3: True}


def get_genotype(animal_id):
    if animal_id.startswith("NYU"):
        return "WT"
    elif animal_id.startswith("CSP"):
        return "Cntnap2"
    elif animal_id.startswith("FMR"):
        return "Fmr1"
    elif animal_id.startswith("SH"):
        return "Shank3B"
    return "Unknown"


def load_zero_contrast_data(criterion=600):
    """Load 0%-contrast trials from biased blocks, truncated to criterion.

    At 0% contrast there is no physical stimulus. The 'correct' side is
    determined by the block prior (probLeft=0.8 → left is rewarded, etc).

    Encoding: choice × rewarded_side → 4 symbols:
      0 = chose left,  reward side = left  (correct)
      1 = chose left,  reward side = right (error)
      2 = chose right, reward side = left  (error)
      3 = chose right, reward side = right (correct)

    Returns:
        subjects_by_genotype: dict mapping genotype → list of subject arrays
            Each array has columns (0, symbol, reward, block_type) matching
            the format expected by build_count_matrix (session=0, all one block).
        stats: dict with per-animal trial counts and accuracy
        block_info: dict mapping animal_id → per-trial block sequence and
            per-block choice probabilities (for surrogate generation)
    """
    mat = sio.loadmat(DATA_DIR / "raw" / "summary_behavior.mat")
    animal_arr = np.array([str(a[0]) for a in mat["master_animal"].ravel()])
    choice = mat["master_choice"].ravel()       # 1=right, -1=left
    contrast = mat["master_contrast"].ravel()   # signed contrast
    feedback = mat["master_feedback"].ravel()   # 1=correct, -1=incorrect
    probLeft = mat["master_probLeft"].ravel()   # 0.2, 0.5, 0.8

    # Filter: 0% contrast, biased blocks, responded
    mask = (contrast == 0) & (probLeft != 0.5) & (choice != 0)

    # Determine reward side from block: probLeft=0.8 → left rewarded,
    # probLeft=0.2 → right rewarded
    reward_side_is_left = probLeft == 0.8  # shape: full array

    subjects_by_genotype = {"WT": [], "Cntnap2": [], "Fmr1": [], "Shank3B": []}
    stats = {"animals": {}, "excluded": []}
    block_info = {}

    for aid in sorted(np.unique(animal_arr)):
        animal_mask = (animal_arr == aid) & mask
        n_available = animal_mask.sum()

        if n_available < criterion:
            stats["excluded"].append((aid, get_genotype(aid), int(n_available)))
            continue

        # Get this animal's trials
        idx = np.where(animal_mask)[0][:criterion]  # truncate to criterion

        ch = choice[idx]          # 1=right, -1=left
        fb = feedback[idx]        # 1=correct, -1=incorrect
        rew_left = reward_side_is_left[idx]

        # Encode: choice × reward_side → 4 symbols
        chose_left = ch == -1
        symbols = np.zeros(criterion, dtype=np.int32)
        symbols[chose_left & rew_left] = 0       # L✓
        symbols[chose_left & ~rew_left] = 1      # L✗
        symbols[~chose_left & ~rew_left] = 3     # R✓
        symbols[~chose_left & rew_left] = 2      # R✗

        # Pack into the format expected by build_count_matrix:
        # columns: session, symbol, reward, contingency
        # We use session=0, contingency=0 (will use contingency=None)
        reward_col = (fb == 1).astype(np.int32)
        arr = np.column_stack([
            np.zeros(criterion, dtype=np.int32),
            symbols,
            reward_col,
            np.zeros(criterion, dtype=np.int32),
        ])

        genotype = get_genotype(aid)
        subjects_by_genotype[genotype].append(arr)

        accuracy = (fb == 1).mean()
        left_bias = chose_left.mean()
        stats["animals"][aid] = {
            "genotype": genotype,
            "n_available": int(n_available),
            "n_used": criterion,
            "accuracy": float(accuracy),
            "left_bias": float(left_bias),
        }

        # Store block info for surrogate generation
        p_left_left_block = (chose_left[rew_left].mean()
                             if rew_left.any() else 0.5)
        p_left_right_block = (chose_left[~rew_left].mean()
                              if (~rew_left).any() else 0.5)
        block_info[aid] = {
            "genotype": genotype,
            "block_seq": rew_left.astype(np.int32),  # 1=left-biased, 0=right-biased
            "p_left_left_block": float(p_left_left_block),
            "p_left_right_block": float(p_left_right_block),
        }

    return subjects_by_genotype, stats, block_info


def generate_surrogates(block_info, subjects_by_genotype, criterion, rng):
    """Block-aware virtual-twin surrogate.

    Preserves:
      - The exact block-transition sequence each animal experienced
      - Each animal's P(choose_left) separately for left-biased and right-biased blocks

    Destroys:
      - All trial-to-trial sequential dependence within blocks
      - Any history-dependent or engagement-state structure

    This is a stronger null than simple iid-marginal surrogates because
    block transitions generate non-trivial sequential structure (symbol
    distribution shifts at block boundaries) even with iid choices within blocks.
    """
    surrogate_by_genotype = {"WT": [], "Cntnap2": [], "Fmr1": [], "Shank3B": []}

    for genotype, subjects in subjects_by_genotype.items():
        animal_ids = [aid for aid, info in block_info.items()
                      if info["genotype"] == genotype]

        for aid in animal_ids:
            info = block_info[aid]
            block_seq = info["block_seq"]
            p_left_L = info["p_left_left_block"]
            p_left_R = info["p_left_right_block"]

            # Generate iid choices per trial, conditioned on block type
            p_left_per_trial = np.where(block_seq == 1, p_left_L, p_left_R)
            chose_left = rng.random(criterion) < p_left_per_trial
            reward_side_left = block_seq.astype(bool)

            # Encode
            symbols = np.zeros(criterion, dtype=np.int32)
            symbols[chose_left & reward_side_left] = 0       # L✓
            symbols[chose_left & ~reward_side_left] = 1      # L✗
            symbols[~chose_left & ~reward_side_left] = 3     # R✓
            symbols[~chose_left & reward_side_left] = 2      # R✗

            is_correct = (symbols == 0) | (symbols == 3)
            arr = np.column_stack([
                np.zeros(criterion, dtype=np.int32),
                symbols,
                is_correct.astype(np.int32),
                np.zeros(criterion, dtype=np.int32),
            ])
            surrogate_by_genotype[genotype].append(arr)

    return surrogate_by_genotype


def run_comparison(grp0_subjects, grp1_subjects, grp0_name, grp1_name, params):
    """Run one CBAS comparison. Returns results dict."""
    subjects_data = grp0_subjects + grp1_subjects
    group_labels = np.array([0] * len(grp0_subjects) + [1] * len(grp1_subjects))
    group_indices = [np.where(group_labels == 0)[0], np.where(group_labels == 1)[0]]

    t0 = time.perf_counter()
    sequences, count_matrix = build_count_matrix(
        subjects_data, params, contingency=None, encode_reward=False
    )
    test_stats = compute_test_stats(count_matrix, group_indices)
    g_values, k_final = find_k_fwer_chunked(
        test_stats, count_matrix, group_indices, params, chunk_size=500
    )
    elapsed = time.perf_counter() - t0

    n_seq = len(sequences)
    sig_summary = compute_significance_summary(g_values, n_seq, params.alpha)

    # Collect significant sequences with details
    sig_seqs = []
    for i in range(n_seq):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        if not np.isnan(pos_g) and pos_g < params.alpha:
            sig_seqs.append({
                "seq": sequences[i],
                "direction": f"{grp0_name}>{grp1_name}",
                "g_value": float(pos_g),
            })
        if not np.isnan(neg_g) and neg_g < params.alpha:
            sig_seqs.append({
                "seq": sequences[i],
                "direction": f"{grp1_name}>{grp0_name}",
                "g_value": float(neg_g),
            })
    sig_seqs.sort(key=lambda x: x["g_value"])

    return {
        "sequences": sequences,
        "g_values": g_values,
        "test_stats": test_stats,
        "k_final": k_final,
        "n_seq": n_seq,
        "sig_summary": sig_summary,
        "sig_seqs": sig_seqs,
        "elapsed": elapsed,
        "n_grp0": len(grp0_subjects),
        "n_grp1": len(grp1_subjects),
    }


def run_analysis(seq_len=6, criterion=600, quick=False, skip_surrogates=False,
                 n_surrogate_runs=5):
    print("=" * 70)
    print("CBAS — Noel IBL: 0%% Contrast, Biased Blocks, Criterion-Truncated")
    print("=" * 70)

    # Load data
    print(f"\nLoading 0%%-contrast biased-block trials (criterion={criterion})...")
    subjects_by_genotype, stats, block_info = load_zero_contrast_data(criterion=criterion)

    for g in ["WT", "Cntnap2", "Fmr1", "Shank3B"]:
        n_mice = len(subjects_by_genotype[g])
        if n_mice > 0:
            accs = [stats["animals"][a]["accuracy"]
                    for a, s in stats["animals"].items() if s["genotype"] == g]
            biases = [stats["animals"][a]["left_bias"]
                      for a, s in stats["animals"].items() if s["genotype"] == g]
            print(f"  {g:8s}: {n_mice:2d} mice × {criterion} trials, "
                  f"acc={np.mean(accs)*100:.1f}%, left_bias={np.mean(biases)*100:.1f}%")

    if stats["excluded"]:
        print(f"\n  Excluded ({len(stats['excluded'])} animals < {criterion} trials):")
        for aid, geno, n in stats["excluded"]:
            print(f"    {aid} ({geno}): {n} trials")

    M = 1000 if quick else 10000
    params = CBASParams(
        num_arms=4, seq_len_max=seq_len, criterion=criterion,
        resample_number=M
    )
    print(f"\nParams: A=4, L={seq_len}, criterion={criterion}, M={M:,}")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # WT vs each model
    comparisons = [
        ("WT", "Cntnap2"),
        ("WT", "Fmr1"),
        ("WT", "Shank3B"),
    ]

    all_results = {}
    for grp0_name, grp1_name in comparisons:
        grp0 = subjects_by_genotype[grp0_name]
        grp1 = subjects_by_genotype[grp1_name]
        if not grp0 or not grp1:
            print(f"\n--- {grp0_name} vs {grp1_name}: SKIPPED (empty group) ---")
            continue

        print(f"\n--- {grp0_name} ({len(grp0)}) vs {grp1_name} ({len(grp1)}) ---")
        result = run_comparison(grp0, grp1, grp0_name, grp1_name, params)
        all_results[f"{grp0_name}_vs_{grp1_name}"] = result

        n_sig = result["sig_summary"]["n_significant"]
        n_pos = result["sig_summary"]["n_positive"]
        n_neg = result["sig_summary"]["n_negative"]
        print(f"  Sequences: {result['n_seq']:,}, Significant: {n_sig}, k={result['k_final']}")
        print(f"  {grp0_name} > {grp1_name}: {n_pos}, {grp1_name} > {grp0_name}: {n_neg}")
        print(f"  Time: {result['elapsed']:.1f}s")

        if result["sig_seqs"]:
            print(f"  Significant sequences:")
            for s in result["sig_seqs"][:10]:
                decoded = " ".join(SYM[x] for x in s["seq"])
                n_err = sum(1 for x in s["seq"] if not IS_CORRECT[x])
                print(f"    g={s['g_value']:.4f} {s['direction']:20s} "
                      f"len={len(s['seq'])} err={n_err} {decoded}")

        # Save JSON
        comp_key = f"{grp0_name.lower()}_vs_{grp1_name.lower()}"
        results_json = {
            "dataset": "noel_ibl_mice_zero_contrast",
            "mode": "comparative",
            "comparison": f"{grp0_name}_vs_{grp1_name}",
            "filter": "0% contrast, biased blocks only",
            "groups": {grp0_name: len(grp0), grp1_name: len(grp1)},
            "params": {
                "num_arms": params.num_arms,
                "seq_len_max": params.seq_len_max,
                "criterion": criterion,
                "encode_reward": False,
                "resample_number": params.resample_number,
                "alpha": params.alpha,
                "gamma": params.gamma,
            },
            "results": {
                "n_subjects": len(grp0) + len(grp1),
                "n_sequences": result["n_seq"],
                "n_significant": result["sig_summary"]["n_significant"],
                "n_positive": result["sig_summary"]["n_positive"],
                "n_negative": result["sig_summary"]["n_negative"],
                "fraction_significant": result["sig_summary"]["fraction_significant"],
                "k_final": result["k_final"],
            },
            "timing": {"total": result["elapsed"]},
        }
        save_results_json(RESULTS_DIR / f"results_{comp_key}.json", results_json)

    # Cross-model comparisons
    cross_comparisons = [
        ("Fmr1", "Shank3B"),
        ("Cntnap2", "Fmr1"),
        ("Cntnap2", "Shank3B"),
    ]

    print("\n" + "=" * 70)
    print("CROSS-MODEL COMPARISONS")
    print("=" * 70)

    for grp0_name, grp1_name in cross_comparisons:
        grp0 = subjects_by_genotype[grp0_name]
        grp1 = subjects_by_genotype[grp1_name]
        if not grp0 or not grp1:
            continue

        print(f"\n--- {grp0_name} ({len(grp0)}) vs {grp1_name} ({len(grp1)}) ---")
        result = run_comparison(grp0, grp1, grp0_name, grp1_name, params)
        comp_key = f"{grp0_name.lower()}_vs_{grp1_name.lower()}"
        all_results[comp_key] = result

        n_sig = result["sig_summary"]["n_significant"]
        n_pos = result["sig_summary"]["n_positive"]
        n_neg = result["sig_summary"]["n_negative"]
        print(f"  Sequences: {result['n_seq']:,}, Significant: {n_sig}, k={result['k_final']}")
        print(f"  {grp0_name} > {grp1_name}: {n_pos}, {grp1_name} > {grp0_name}: {n_neg}")

        if result["sig_seqs"]:
            for s in result["sig_seqs"][:5]:
                decoded = " ".join(SYM[x] for x in s["seq"])
                print(f"    g={s['g_value']:.4f} {s['direction']:20s} {decoded}")

        results_json = {
            "dataset": "noel_ibl_mice_zero_contrast",
            "mode": "comparative",
            "comparison": f"{grp0_name}_vs_{grp1_name}",
            "filter": "0% contrast, biased blocks only",
            "groups": {grp0_name: len(grp0), grp1_name: len(grp1)},
            "params": {
                "num_arms": params.num_arms,
                "seq_len_max": params.seq_len_max,
                "criterion": criterion,
                "encode_reward": False,
                "resample_number": params.resample_number,
            },
            "results": {
                "n_sequences": result["n_seq"],
                "n_significant": result["sig_summary"]["n_significant"],
                "n_positive": result["sig_summary"]["n_positive"],
                "n_negative": result["sig_summary"]["n_negative"],
                "k_final": result["k_final"],
            },
            "timing": {"total": result["elapsed"]},
        }
        save_results_json(RESULTS_DIR / f"results_{comp_key}.json", results_json)

    # Summary table
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'Comparison':<25} {'Seqs':>8} {'Sig':>6} {'Pos':>5} {'Neg':>5} {'k':>3} {'Time':>6}")
    print("-" * 70)
    for key, r in all_results.items():
        n_sig = r["sig_summary"]["n_significant"]
        print(f"{key:<25} {r['n_seq']:>8,} {n_sig:>6} "
              f"{r['sig_summary']['n_positive']:>5} "
              f"{r['sig_summary']['n_negative']:>5} "
              f"{r['k_final']:>3} {r['elapsed']:>5.1f}s")

    # Virtual-twin surrogate control
    if skip_surrogates:
        print("\n(Surrogates skipped)")
        return all_results, stats

    print("\n" + "=" * 70)
    print(f"VIRTUAL-TWIN SURROGATE CONTROL ({n_surrogate_runs} runs)")
    print("=" * 70)
    print("Block-aware surrogates: preserve block-transition structure,")
    print("iid choices within each block at animal's block-specific rate.")
    print("Tests whether real result exceeds what block structure + marginals produce.\n")

    surrogate_sig_counts = defaultdict(list)
    rng = np.random.default_rng(42)

    for run_idx in range(n_surrogate_runs):
        surrogate_by_genotype = generate_surrogates(
            block_info, subjects_by_genotype, criterion, rng
        )

        for grp0_name, grp1_name in comparisons:
            grp0 = surrogate_by_genotype[grp0_name]
            grp1 = surrogate_by_genotype[grp1_name]
            if not grp0 or not grp1:
                continue

            result = run_comparison(grp0, grp1, grp0_name, grp1_name, params)
            key = f"{grp0_name}_vs_{grp1_name}"
            surrogate_sig_counts[key].append(result["sig_summary"]["n_significant"])

        print(f"  Run {run_idx + 1}/{n_surrogate_runs} complete")

    # Report surrogate results
    print(f"\n{'Comparison':<25} {'Real sig':>10} {'Surrogate sig (mean±sd)':>25}")
    print("-" * 65)
    for grp0_name, grp1_name in comparisons:
        key = f"{grp0_name}_vs_{grp1_name}"
        if key in all_results and key in surrogate_sig_counts:
            real_sig = all_results[key]["sig_summary"]["n_significant"]
            surr_counts = surrogate_sig_counts[key]
            surr_mean = np.mean(surr_counts)
            surr_std = np.std(surr_counts)
            surr_max = max(surr_counts)
            exceeds = "*** EXCEEDS SURROGATES" if real_sig > surr_max else ""
            print(f"{key:<25} {real_sig:>10} "
                  f"{surr_mean:>10.1f} ± {surr_std:.1f} (max={surr_max}) {exceeds}")

    # Save surrogate results
    surrogate_json = {
        "n_runs": n_surrogate_runs,
        "method": "block-aware iid surrogates",
        "description": "Preserves each animal's block-transition sequence and per-block "
                       "P(choose_left). Choices are iid within each block. Destroys "
                       "trial-to-trial sequential dependence while preserving block "
                       "structure and marginal statistics.",
        "results": {},
    }
    for key, counts in surrogate_sig_counts.items():
        real_sig = all_results[key]["sig_summary"]["n_significant"] if key in all_results else 0
        surrogate_json["results"][key] = {
            "real_n_significant": real_sig,
            "surrogate_n_significant": counts,
            "surrogate_mean": float(np.mean(counts)),
            "surrogate_std": float(np.std(counts)),
            "surrogate_max": int(max(counts)),
            "real_exceeds_all_surrogates": bool(real_sig > max(counts)),
        }
    save_results_json(RESULTS_DIR / "surrogate_control.json", surrogate_json)

    # Generate figures
    try:
        make_figures(all_results, params)
    except Exception as e:
        print(f"\n(Figures skipped: {e})")

    return all_results, stats


def make_figures(all_results, params):
    """Generate manhattan plots."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT_DIR / "plots" / "matplotlib"))
    from manhattan import manhattan_plot

    for key, result in all_results.items():
        g_values = result["g_values"]
        seq_lengths = np.array([len(s) for s in result["sequences"]])
        n_sig = result["sig_summary"]["n_significant"]

        fig, ax = manhattan_plot(
            g_values, seq_lengths, alpha=params.alpha,
            title=f"Noel IBL 0%% contrast: {key} ({n_sig} sig / {result['n_seq']:,} seq)"
        )
        fig.savefig(FIG_DIR / f"manhattan_{key.lower()}.png",
                    dpi=150, bbox_inches="tight")
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="CBAS on Noel IBL — 0%% contrast, biased blocks, criterion-truncated"
    )
    parser.add_argument("--seq-len", type=int, default=6,
                        help="Max sequence length (default 6)")
    parser.add_argument("--criterion", type=int, default=600,
                        help="Per-animal trial truncation (default 600)")
    parser.add_argument("--quick", action="store_true",
                        help="Reduced M=1000 for fast check")
    parser.add_argument("--skip-surrogates", action="store_true",
                        help="Skip virtual-twin surrogate control")
    parser.add_argument("--surrogate-runs", type=int, default=5,
                        help="Number of surrogate runs (default 5)")
    args = parser.parse_args()

    run_analysis(
        seq_len=args.seq_len,
        criterion=args.criterion,
        quick=args.quick,
        skip_surrogates=args.skip_surrogates,
        n_surrogate_runs=args.surrogate_runs,
    )


if __name__ == "__main__":
    main()
