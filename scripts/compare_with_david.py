"""
Compare our significant sequences with David's ground-truth lists.

Parses David's exported sig sequence files (notes/flyCBASsigSeq.txt,
notes/humanCBASsigSeq.txt) and compares against our results.
"""

import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

ROOT_DIR = Path(__file__).parent.parent
NOTES_DIR = ROOT_DIR / "notes"


def parse_david_file(path, max_seq_len):
    """Parse David's significant sequence file.

    Format: seq_symbol_1, seq_symbol_2, ..., [empty padding], direction, p_value
    Fixed-width: max_seq_len columns for sequence, then direction, then p-value.
    """
    sequences = []
    with open(path) as f:
        for line in f:
            line = line.strip().replace("\r", "")
            if not line:
                continue
            parts = line.split(",")
            seq_parts = parts[:max_seq_len]
            seq = tuple(int(x) for x in seq_parts if x.strip() != "")
            direction = int(parts[max_seq_len])
            pvalue = float(parts[max_seq_len + 1])
            sequences.append({"seq": seq, "direction": direction, "pvalue": pvalue})
    return sequences


def load_our_results(species):
    """Load our g-values and sequences from cached NPZ."""
    npz_path = ROOT_DIR / "results" / species / "figures" / "results.npz"
    data = np.load(npz_path, allow_pickle=False)
    g_values = data["g_values"]
    seq_strs = data["seq_strs"]
    n_seq = len(seq_strs)

    results = {}
    for i in range(n_seq):
        seq = tuple(int(x) for x in seq_strs[i].split("-"))
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        results[seq] = {"pos_g": pos_g, "neg_g": neg_g, "index": i}
    return results


def compare_flies():
    print("=" * 70)
    print("FLY COMPARISON")
    print("=" * 70)

    david = parse_david_file(NOTES_DIR / "flyCBASsigSeq.txt", max_seq_len=10)
    ours = load_our_results("flies")

    print(f"David: {len(david)} significant sequences")
    print(f"Ours: {len(ours)} total sequences evaluated")

    alpha = 0.5
    our_sig = set()
    for seq, vals in ours.items():
        pos_g, neg_g = vals["pos_g"], vals["neg_g"]
        if (not np.isnan(pos_g) and pos_g < alpha) or (not np.isnan(neg_g) and neg_g < alpha):
            our_sig.add(seq)

    david_sig = set(d["seq"] for d in david)

    print(f"\nOurs significant: {len(our_sig)}")
    print(f"David significant: {len(david_sig)}")

    both = our_sig & david_sig
    only_ours = our_sig - david_sig
    only_david = david_sig - our_sig

    print(f"\nBoth: {len(both)}")
    print(f"Only ours (overcalled): {len(only_ours)}")
    print(f"Only David (missed): {len(only_david)}")

    if only_david:
        print(f"\n--- Sequences David found significant but we MISSED ---")
        for d in david:
            if d["seq"] in only_david:
                seq = d["seq"]
                our_vals = ours.get(seq)
                if our_vals:
                    print(f"  {seq}  david_p={d['pvalue']:.6f}  "
                          f"our_pos_g={our_vals['pos_g']:.6f}  our_neg_g={our_vals['neg_g']:.6f}")
                else:
                    print(f"  {seq}  david_p={d['pvalue']:.6f}  NOT IN OUR SEQUENCE SET")

    # p-value comparison for shared sequences
    print(f"\n--- P-value comparison (David vs ours, first 20 by David's p-value) ---")
    david_sorted = sorted(david, key=lambda x: x["pvalue"])
    print(f"{'Sequence':<30} {'Dir':<5} {'David_p':<12} {'Our_g':<12} {'Match?'}")
    for d in david_sorted[:20]:
        seq = d["seq"]
        our_vals = ours.get(seq)
        if our_vals:
            # direction 0 = CA>w1118 (positive), 1 = CA<w1118 (negative)
            our_g = our_vals["neg_g"] if d["direction"] == 1 else our_vals["pos_g"]
            match = "~" if abs(our_g - d["pvalue"]) < 0.05 else "X"
            print(f"  {str(seq):<28} {d['direction']:<5} {d['pvalue']:<12.6f} {our_g:<12.6f} {match}")

    # Direction agreement check
    print(f"\n--- Direction agreement (for sequences both call significant) ---")
    dir_agree = 0
    dir_disagree = 0
    for d in david:
        seq = d["seq"]
        if seq in both:
            our_vals = ours[seq]
            pos_g, neg_g = our_vals["pos_g"], our_vals["neg_g"]
            our_dir = 0 if (np.isnan(neg_g) or pos_g <= neg_g) else 1
            if our_dir == d["direction"]:
                dir_agree += 1
            else:
                dir_disagree += 1
    print(f"  Agree: {dir_agree}, Disagree: {dir_disagree}")

    # Overcalled sequences: what are their g-values?
    if only_ours:
        overcalled_gs = []
        for seq in only_ours:
            vals = ours[seq]
            best_g = min(
                vals["pos_g"] if not np.isnan(vals["pos_g"]) else 1.0,
                vals["neg_g"] if not np.isnan(vals["neg_g"]) else 1.0,
            )
            overcalled_gs.append(best_g)
        overcalled_gs = np.array(overcalled_gs)
        print(f"\n--- Overcalled sequences g-value distribution ---")
        print(f"  Count: {len(overcalled_gs)}")
        print(f"  Min: {overcalled_gs.min():.6f}")
        print(f"  Max: {overcalled_gs.max():.6f}")
        print(f"  Mean: {overcalled_gs.mean():.6f}")
        print(f"  Median: {np.median(overcalled_gs):.6f}")
        pctiles = [10, 25, 50, 75, 90]
        for p in pctiles:
            print(f"  P{p}: {np.percentile(overcalled_gs, p):.6f}")


def compare_humans():
    print("\n" + "=" * 70)
    print("HUMAN COMPARISON")
    print("=" * 70)

    david = parse_david_file(NOTES_DIR / "humanCBASsigSeq.txt", max_seq_len=4)
    ours = load_our_results("humans")

    print(f"David: {len(david)} significant sequences")
    print(f"Ours: {len(ours)} total sequences evaluated")

    alpha = 0.5
    our_sig = set()
    for seq, vals in ours.items():
        pos_g, neg_g = vals["pos_g"], vals["neg_g"]
        if (not np.isnan(pos_g) and pos_g < alpha) or (not np.isnan(neg_g) and neg_g < alpha):
            our_sig.add(seq)

    david_sig = set(d["seq"] for d in david)

    print(f"\nOurs significant: {len(our_sig)}")
    print(f"David significant: {len(david_sig)}")

    both = our_sig & david_sig
    only_ours = our_sig - david_sig
    only_david = david_sig - our_sig

    print(f"\nBoth: {len(both)}")
    print(f"Only ours (overcalled): {len(only_ours)}")
    print(f"Only David (missed): {len(only_david)}")

    if only_david:
        print(f"\n--- Sequences David found significant but we MISSED ---")
        for d in david:
            if d["seq"] in only_david:
                seq = d["seq"]
                our_vals = ours.get(seq)
                if our_vals:
                    print(f"  {seq}  david_p={d['pvalue']:.6f}  "
                          f"our_pos_g={our_vals['pos_g']:.6f}  our_neg_g={our_vals['neg_g']:.6f}")
                else:
                    print(f"  {seq}  david_p={d['pvalue']:.6f}  NOT IN OUR SEQUENCE SET")

    # p-value comparison
    print(f"\n--- P-value comparison (all of David's sequences) ---")
    david_sorted = sorted(david, key=lambda x: x["pvalue"])
    print(f"{'Sequence':<20} {'Dir':<5} {'David_p':<12} {'Our_g':<12} {'Diff':<10}")
    for d in david_sorted:
        seq = d["seq"]
        our_vals = ours.get(seq)
        if our_vals:
            # direction 0 = pos_corr, 1 = neg_corr
            our_g = our_vals["neg_g"] if d["direction"] == 1 else our_vals["pos_g"]
            diff = our_g - d["pvalue"]
            print(f"  {str(seq):<18} {d['direction']:<5} {d['pvalue']:<12.6f} {our_g:<12.6f} {diff:<+10.6f}")
        else:
            print(f"  {str(seq):<18} {d['direction']:<5} {d['pvalue']:<12.6f} {'N/A':<12}")

    # Direction agreement
    print(f"\n--- Direction agreement (for sequences both call significant) ---")
    dir_agree = 0
    dir_disagree = 0
    for d in david:
        seq = d["seq"]
        if seq in both:
            our_vals = ours[seq]
            pos_g, neg_g = our_vals["pos_g"], our_vals["neg_g"]
            our_dir = 0 if (np.isnan(neg_g) or pos_g <= neg_g) else 1
            if our_dir == d["direction"]:
                dir_agree += 1
            else:
                dir_disagree += 1
    print(f"  Agree: {dir_agree}, Disagree: {dir_disagree}")

    # Overcalled
    if only_ours:
        overcalled_gs = []
        for seq in only_ours:
            vals = ours[seq]
            best_g = min(
                vals["pos_g"] if not np.isnan(vals["pos_g"]) else 1.0,
                vals["neg_g"] if not np.isnan(vals["neg_g"]) else 1.0,
            )
            overcalled_gs.append(best_g)
        overcalled_gs = np.array(overcalled_gs)
        print(f"\n--- Overcalled sequences g-value distribution ---")
        print(f"  Count: {len(overcalled_gs)}")
        print(f"  Min: {overcalled_gs.min():.6f}")
        print(f"  Max: {overcalled_gs.max():.6f}")
        print(f"  Mean: {overcalled_gs.mean():.6f}")
        print(f"  Median: {np.median(overcalled_gs):.6f}")


if __name__ == "__main__":
    compare_flies()
    compare_humans()
