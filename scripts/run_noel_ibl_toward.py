"""
CBAS on Noel et al. 2025 IBL mice — toward/against encoding.

Recodes choice as toward vs. against the currently favored (rewarded) side,
rather than left vs. right. This removes absolute side from the encoding,
asking whether genotypes differ in how they orient relative to the block.

Symbols: T✓ (toward, correct), T✗ (toward, incorrect),
         A✗ (against, incorrect), A✓ (against, correct)

Usage:
    pixi run noel-ibl-toward
    pixi run noel-ibl-toward-quick
"""

import argparse
import sys
import time
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent))
from pycbas import (
    CBASParams,
    build_count_matrix,
    compute_test_stats,
    find_k_fwer_chunked,
)
from results_io import save_results_json, compute_significance_summary

ROOT_DIR = Path(__file__).parent.parent
DATA_DIR = ROOT_DIR / "data" / "noel_ibl_mice" / "toward_against"
RESULTS_DIR = ROOT_DIR / "results" / "noel_ibl_mice_toward"
FIG_DIR = RESULTS_DIR / "figures"

SYM = {0: "T✓", 1: "T✗", 2: "A✗", 3: "A✓"}


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


def load_data(criterion=4000):
    subjects_by_genotype = {"WT": [], "Cntnap2": [], "Fmr1": [], "Shank3B": []}
    stats = {"animals": {}, "excluded": []}

    for txt_path in sorted(DATA_DIR.glob("*.txt")):
        aid = txt_path.stem
        data = np.loadtxt(txt_path, delimiter=",", dtype=np.int32)
        n_available = len(data)

        if n_available < criterion:
            stats["excluded"].append((aid, get_genotype(aid), n_available))
            continue

        data = data[:criterion]
        genotype = get_genotype(aid)
        subjects_by_genotype[genotype].append(data)

        symbols = data[:, 1]
        toward_rate = ((symbols == 0) | (symbols == 1)).mean()
        accuracy = ((symbols == 0) | (symbols == 3)).mean()
        stats["animals"][aid] = {
            "genotype": genotype,
            "n_available": n_available,
            "n_used": criterion,
            "toward_rate": float(toward_rate),
            "accuracy": float(accuracy),
        }

    return subjects_by_genotype, stats


def run_comparison(grp0_subjects, grp1_subjects, grp0_name, grp1_name, params):
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

    sig_seqs = []
    for i in range(n_seq):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        if not np.isnan(pos_g) and pos_g < params.alpha:
            sig_seqs.append({"seq": sequences[i], "direction": f"{grp0_name}>{grp1_name}", "g_value": float(pos_g)})
        if not np.isnan(neg_g) and neg_g < params.alpha:
            sig_seqs.append({"seq": sequences[i], "direction": f"{grp1_name}>{grp0_name}", "g_value": float(neg_g)})
    sig_seqs.sort(key=lambda x: x["g_value"])

    return {
        "sequences": sequences,
        "g_values": g_values,
        "k_final": k_final,
        "n_seq": n_seq,
        "sig_summary": sig_summary,
        "sig_seqs": sig_seqs,
        "elapsed": elapsed,
    }


def run_analysis(seq_len=6, criterion=4000, quick=False):
    print("=" * 70)
    print("CBAS — Noel IBL: Toward/Against Encoding, Criterion-Truncated")
    print("=" * 70)

    subjects_by_genotype, stats = load_data(criterion=criterion)

    for g in ["WT", "Cntnap2", "Fmr1", "Shank3B"]:
        n_mice = len(subjects_by_genotype[g])
        if n_mice > 0:
            trs = [s["toward_rate"] for s in stats["animals"].values() if s["genotype"] == g]
            accs = [s["accuracy"] for s in stats["animals"].values() if s["genotype"] == g]
            print(f"  {g:8s}: {n_mice:2d} mice × {criterion} trials, "
                  f"toward={np.mean(trs)*100:.1f}% ± {np.std(trs)*100:.1f}%, "
                  f"acc={np.mean(accs)*100:.1f}%")

    if stats["excluded"]:
        print(f"\n  Excluded ({len(stats['excluded'])} animals < {criterion} trials):")
        for aid, geno, n in stats["excluded"]:
            print(f"    {aid} ({geno}): {n}")

    M = 1000 if quick else 10000
    params = CBASParams(
        num_arms=4, seq_len_max=seq_len, criterion=criterion, resample_number=M
    )
    print(f"\nParams: A=4, L={seq_len}, criterion={criterion}, M={M:,}")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    comparisons = [("WT", "Cntnap2"), ("WT", "Fmr1"), ("WT", "Shank3B")]
    all_results = {}

    for grp0_name, grp1_name in comparisons:
        grp0 = subjects_by_genotype[grp0_name]
        grp1 = subjects_by_genotype[grp1_name]
        if not grp0 or not grp1:
            continue

        print(f"\n--- {grp0_name} ({len(grp0)}) vs {grp1_name} ({len(grp1)}) ---")
        result = run_comparison(grp0, grp1, grp0_name, grp1_name, params)
        key = f"{grp0_name}_vs_{grp1_name}"
        all_results[key] = result

        n_sig = result["sig_summary"]["n_significant"]
        n_pos = result["sig_summary"]["n_positive"]
        n_neg = result["sig_summary"]["n_negative"]
        print(f"  Sequences: {result['n_seq']:,}, Significant: {n_sig}, k={result['k_final']}")
        print(f"  {grp0_name} > {grp1_name}: {n_pos}, {grp1_name} > {grp0_name}: {n_neg}")
        print(f"  Time: {result['elapsed']:.1f}s")

        if result["sig_seqs"]:
            print(f"  Top significant:")
            for s in result["sig_seqs"][:8]:
                decoded = " ".join(SYM[x] for x in s["seq"])
                print(f"    g={s['g_value']:.4f} {s['direction']:20s} {decoded}")

        results_json = {
            "dataset": "noel_ibl_mice_toward",
            "comparison": f"{grp0_name}_vs_{grp1_name}",
            "encoding": "toward/against favored side x correctness",
            "groups": {grp0_name: len(grp0), grp1_name: len(grp1)},
            "params": {"num_arms": 4, "seq_len_max": seq_len, "criterion": criterion,
                       "resample_number": M},
            "results": {
                "n_sequences": result["n_seq"],
                "n_significant": n_sig, "n_positive": n_pos, "n_negative": n_neg,
                "k_final": result["k_final"],
            },
            "timing": {"total": result["elapsed"]},
        }
        save_results_json(RESULTS_DIR / f"results_{key.lower()}.json", results_json)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'Comparison':<25} {'Seqs':>8} {'Sig':>5} {'Pos':>5} {'Neg':>5} {'k':>3} {'Time':>6}")
    print("-" * 60)
    for key, r in all_results.items():
        n_sig = r["sig_summary"]["n_significant"]
        print(f"{key:<25} {r['n_seq']:>8,} {n_sig:>5} "
              f"{r['sig_summary']['n_positive']:>5} "
              f"{r['sig_summary']['n_negative']:>5} "
              f"{r['k_final']:>3} {r['elapsed']:>5.1f}s")

    return all_results, stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seq-len", type=int, default=6)
    parser.add_argument("--criterion", type=int, default=4000)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    run_analysis(seq_len=args.seq_len, criterion=args.criterion, quick=args.quick)


if __name__ == "__main__":
    main()
