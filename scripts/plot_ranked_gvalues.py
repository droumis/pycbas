"""
Generate ranked ζ-value plots for all species.

For species with David's (Igor) ground truth results, shows agreement/disagreement.
For species without, shows significant vs not significant.

Outputs: results/<species>/figures/ranked_gvalues.png

Usage:
    pixi run python scripts/plot_ranked_gvalues.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent))
from compare_with_david import parse_david_file, load_our_results

ROOT_DIR = Path(__file__).parent.parent
NOTES_DIR = ROOT_DIR / "notes"

SPECIES_CONFIG = [
    {
        "name": "flies",
        "david_file": "flyCBASsigSeq.txt",
        "max_seq_len": 10,
        "title": "Fly CBAS: CA vs w1118",
    },
    {
        "name": "humans",
        "david_file": "humanCBASsigSeq.txt",
        "max_seq_len": 4,
        "title": "Human CBAS: Correlative with CBIT",
    },
    {
        "name": "rats",
        "david_file": "ratSigSeq.txt",
        "max_seq_len": 6,
        "title": "Rat CBAS: Control vs Lesion (block_aware=True)",
    },
]


def get_best_g(vals):
    pos_g, neg_g = vals["pos_g"], vals["neg_g"]
    if not np.isnan(pos_g) and not np.isnan(neg_g):
        return min(pos_g, neg_g)
    elif not np.isnan(pos_g):
        return pos_g
    elif not np.isnan(neg_g):
        return neg_g
    return np.nan


def plot_with_comparison(species_cfg, alpha=0.5):
    """Ranked ζ-values with Igor comparison overlay."""
    david = parse_david_file(
        NOTES_DIR / species_cfg["david_file"],
        max_seq_len=species_cfg["max_seq_len"],
    )
    ours = load_our_results(species_cfg["name"])

    david_sig = set(d["seq"] for d in david)
    all_seqs = sorted(ours.keys(), key=lambda s: ours[s]["index"])
    n_seq = len(all_seqs)

    g_vals = np.full(n_seq, np.nan)
    categories = np.full(n_seq, "", dtype=object)

    our_sig = set()
    for seq in all_seqs:
        vals = ours[seq]
        i = vals["index"]
        best_g = get_best_g(vals)
        g_vals[i] = best_g
        if not np.isnan(best_g) and best_g < alpha:
            our_sig.add(seq)

    both = our_sig & david_sig
    only_ours = our_sig - david_sig
    only_david = david_sig - our_sig
    neither = set(all_seqs) - our_sig - david_sig

    for seq in all_seqs:
        i = ours[seq]["index"]
        if seq in both:
            categories[i] = "both"
        elif seq in only_ours:
            categories[i] = "us_only"
        elif seq in only_david:
            categories[i] = "david_only"
        else:
            categories[i] = "neither"

    valid = ~np.isnan(g_vals)
    sort_idx = np.argsort(g_vals[valid])
    ranked_g = g_vals[valid][sort_idx]
    ranked_cat = categories[valid][sort_idx]
    ranks = np.arange(len(ranked_g))

    overlap_pct = len(both) / len(david_sig) * 100 if david_sig else 0
    style = {
        "both": {"color": "#2ca02c", "marker": ".", "s": 12,
                 "label": f"Both significant ({len(both)}, {overlap_pct:.1f}% of Igor)"},
        "david_only": {"color": "#d94801", "marker": "x", "s": 24,
                       "label": f"Igor only ({len(only_david)})"},
        "us_only": {"color": "#6a3d9a", "marker": "x", "s": 24,
                    "label": f"pycbas only ({len(only_ours)})"},
        "neither": {"color": "#bdbdbd", "marker": ".", "s": 8,
                    "label": f"Neither ({len(neither)})"},
    }

    fig, ax_rank = plt.subplots(figsize=(8, 5))

    # Left panel: ranked g-values
    for cat in ["neither", "both", "us_only", "david_only"]:
        mask = ranked_cat == cat
        if not np.any(mask):
            continue
        ax_rank.scatter(
            ranks[mask], ranked_g[mask],
            c=style[cat]["color"],
            marker=style[cat]["marker"],
            s=style[cat]["s"],
            label=style[cat]["label"],
            zorder=3 if cat in ("david_only", "us_only") else 2,
            edgecolors="none" if style[cat]["marker"] == "." else style[cat]["color"],
            linewidths=1.0,
        )

    ax_rank.axhline(alpha, color="#6baed6", linestyle="--", linewidth=1.2, label=f"α = {alpha}")
    ax_rank.set_xlabel("Sequence rank (by ζ-value)")
    ax_rank.set_ylabel("ζ-value (adjusted p-value)")
    ax_rank.set_title(f"{species_cfg['title']} — Ranked ζ-values")
    ax_rank.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax_rank.set_ylim(-0.02, 1.02)
    ax_rank.set_xlim(-20, len(ranked_g) + 20)

    # Inset zoomed to the threshold region
    near_mask = (ranked_g > alpha - 0.06) & (ranked_g < alpha + 0.06)
    if np.any(near_mask):
        near_ranks = ranks[near_mask]
        ax_inset = ax_rank.inset_axes([0.38, 0.35, 0.42, 0.42])
        for cat in ["neither", "both", "us_only", "david_only"]:
            mask = (ranked_cat == cat) & near_mask
            if not np.any(mask):
                continue
            ax_inset.scatter(
                ranks[mask], ranked_g[mask],
                c=style[cat]["color"],
                marker=style[cat]["marker"],
                s=style[cat]["s"] * 3,
                zorder=3 if cat in ("david_only", "us_only") else 2,
                edgecolors="none" if style[cat]["marker"] == "." else style[cat]["color"],
                linewidths=1.5,
            )
        ax_inset.axhline(alpha, color="#6baed6", linestyle="--", linewidth=1.0)
        ax_inset.set_xlim(near_ranks.min() - 10, near_ranks.max() + 10)
        ax_inset.set_ylim(alpha - 0.055, alpha + 0.055)
        ax_inset.set_xlabel("Rank", fontsize=7)
        ax_inset.set_ylabel("ζ-value", fontsize=7)
        ax_inset.tick_params(labelsize=7)
        ax_inset.set_title("Threshold region", fontsize=8)
        ax_rank.indicate_inset_zoom(ax_inset, edgecolor="0.5", linewidth=0.8)

    fig.tight_layout()
    out_dir = ROOT_DIR / "results" / species_cfg["name"] / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "ranked_gvalues.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {out_path}")


def plot_without_comparison(species_cfg, alpha=0.5):
    """Ranked ζ-values without Igor comparison (significant vs not)."""
    ours = load_our_results(species_cfg["name"])
    all_seqs = sorted(ours.keys(), key=lambda s: ours[s]["index"])
    n_seq = len(all_seqs)

    g_vals = np.full(n_seq, np.nan)
    sig_mask = np.zeros(n_seq, dtype=bool)

    for seq in all_seqs:
        vals = ours[seq]
        i = vals["index"]
        best_g = get_best_g(vals)
        g_vals[i] = best_g
        if not np.isnan(best_g) and best_g < alpha:
            sig_mask[i] = True

    valid = ~np.isnan(g_vals)
    sort_idx = np.argsort(g_vals[valid])
    ranked_g = g_vals[valid][sort_idx]
    ranked_sig = sig_mask[valid][sort_idx]
    ranks = np.arange(len(ranked_g))

    n_sig = int(ranked_sig.sum())
    n_not = int((~ranked_sig).sum())

    fig, ax_rank = plt.subplots(figsize=(8, 5))

    ax_rank.scatter(ranks[~ranked_sig], ranked_g[~ranked_sig],
                    c="#bdbdbd", marker=".", s=8, label=f"Not significant ({n_not})", zorder=2)
    ax_rank.scatter(ranks[ranked_sig], ranked_g[ranked_sig],
                    c="#2ca02c", marker=".", s=12, label=f"Significant ({n_sig})", zorder=2)

    ax_rank.axhline(alpha, color="#6baed6", linestyle="--", linewidth=1.2, label=f"α = {alpha}")
    ax_rank.set_xlabel("Sequence rank (by ζ-value)")
    ax_rank.set_ylabel("ζ-value (adjusted p-value)")
    ax_rank.set_title(f"{species_cfg['title']} — Ranked ζ-values")
    ax_rank.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax_rank.set_ylim(-0.02, 1.02)
    ax_rank.set_xlim(-20, len(ranked_g) + 20)

    fig.tight_layout()
    out_dir = ROOT_DIR / "results" / species_cfg["name"] / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "ranked_gvalues.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {out_path}")


def main():
    for cfg in SPECIES_CONFIG:
        npz = ROOT_DIR / "results" / cfg["name"] / "figures" / "results.npz"
        if not npz.exists():
            print(f"  SKIP {cfg['name']}: no results.npz")
            continue
        if cfg["david_file"] and (NOTES_DIR / cfg["david_file"]).exists():
            plot_with_comparison(cfg)
        else:
            plot_without_comparison(cfg)


if __name__ == "__main__":
    main()
