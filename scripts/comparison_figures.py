"""
Generate side-by-side comparison figures: paper panels vs our results.

Uses the user-provided paper sub-figure screenshots:
  results/figures/ppr_fig5cL-flies.png  (fly manhattan from paper)
  results/figures/ppr_fig5cL-humans.png (human manhattan from paper)
  results/figures/ppr_fig5cR.png        (rat manhattan from paper)
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread
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
PPR_FIG_DIR = ROOT_DIR / "results" / "figures"
FLY_FIG_DIR = ROOT_DIR / "results" / "flies" / "figures"
HUMAN_FIG_DIR = ROOT_DIR / "results" / "humans" / "figures"


def load_fly_results():
    """Load precomputed fly results or rerun."""
    npz_path = ROOT_DIR / "results" / "flies" / "figures" / "results.npz"
    if npz_path.exists():
        data = np.load(npz_path, allow_pickle=True)
        return data
    return None


def plot_our_manhattan(ax, results_npz, title, seq_len_max=10):
    """Plot our manhattan plot on the given axes."""
    g_values = results_npz["g_values"]
    seq_lengths = results_npz["seq_lengths"]

    n_seq = len(seq_lengths)
    for i in range(n_seq):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        g = np.nanmin([pos_g, neg_g]) if not (np.isnan(pos_g) and np.isnan(neg_g)) else np.nan

        if np.isnan(g):
            continue

        neg_log_g = -np.log10(g) if g > 0 else 4.0
        neg_log_g = min(neg_log_g, 4.0)

        seq_len = int(seq_lengths[i])

        color_val = (seq_len - 1) / (seq_len_max - 1)
        color = plt.cm.turbo(color_val)

        ax.scatter(i + 1, neg_log_g, c=[color], s=15, alpha=0.7, edgecolors="none")

    ax.axhline(-np.log10(0.5), color="black", linestyle=":", linewidth=0.5)
    ax.set_xlabel("Sequence")
    ax.set_ylabel("-log₁₀(ζ)")
    ax.set_title(title)
    ax.set_ylim(-0.1, 4.3)
    ax.set_xscale("log")

    sm = plt.cm.ScalarMappable(cmap="turbo", norm=plt.Normalize(1, seq_len_max))
    plt.colorbar(sm, ax=ax, label="Sequence length", shrink=0.8)


def main():
    print("Generating comparison figures...")

    # ---- Fly comparison ----
    paper_fly_path = PPR_FIG_DIR / "ppr_fig5cL-flies.png"
    if paper_fly_path.exists():
        fly_results = load_fly_results()
        if fly_results is not None:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))

            paper_img = imread(str(paper_fly_path))
            axes[0].imshow(paper_img)
            axes[0].axis("off")
            axes[0].set_title("Paper Fig 5c (Flies): CA vs w1118\n1,605 significant sequences")

            plot_our_manhattan(axes[1], fly_results,
                             "pycbas (Flies): CA vs w1118\n2,046 significant (all; adaptive k=103)")

            plt.tight_layout()
            FLY_FIG_DIR.mkdir(parents=True, exist_ok=True)
            plt.savefig(FLY_FIG_DIR / "comparison_manhattan.png", dpi=150, bbox_inches="tight")
            plt.close()
            print(f"  Saved: {FLY_FIG_DIR / 'comparison_manhattan.png'}")
        else:
            print("  Skipping fly comparison (no results.npz)")
    else:
        print(f"  Skipping fly comparison (no paper figure at {paper_fly_path})")

    # ---- Human comparison ----
    paper_human_path = PPR_FIG_DIR / "ppr_fig5cL-humans.png"
    human_npz_path = ROOT_DIR / "results" / "humans" / "figures" / "results.npz"
    if paper_human_path.exists() and human_npz_path.exists():
        human_results = np.load(human_npz_path, allow_pickle=True)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        paper_img = imread(str(paper_human_path))
        axes[0].imshow(paper_img)
        axes[0].axis("off")
        axes[0].set_title("Paper Fig 5c (Humans): CBIT correlation\n31 significant sequences")

        plot_our_manhattan(axes[1], human_results,
                         "pycbas (Humans): CBIT correlation\n(correlative mode)",
                         seq_len_max=4)

        plt.tight_layout()
        HUMAN_FIG_DIR.mkdir(parents=True, exist_ok=True)
        plt.savefig(HUMAN_FIG_DIR / "comparison_manhattan.png", dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {HUMAN_FIG_DIR / 'comparison_manhattan.png'}")
    else:
        print(f"  Skipping human comparison (missing files)")

    print("\nDone.")


if __name__ == "__main__":
    main()
