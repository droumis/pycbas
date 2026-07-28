"""
Generate side-by-side comparison figures: paper screenshots vs our Manhattan plots.

Usage:
    python scripts/make_comparison_figures.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pathlib import Path

ROOT = Path(__file__).parent.parent
PAPER_DIR = ROOT / "results" / "figures"
OUT_DIR = ROOT / "results" / "figures"


SPECIES = [
    {
        "name": "flies",
        "paper_img": PAPER_DIR / "ppr_fig5cL-flies.png",
        "our_img": ROOT / "results" / "flies" / "figures" / "manhattan.png",
        "paper_title": "Paper — Fig 1c left\n(Kastner et al. 2026)",
        "our_title": "pycbas — 1,243/2,046 significant (k=63)",
        "out_name": "comparison_flies.png",
    },
    {
        "name": "humans",
        "paper_img": PAPER_DIR / "ppr_fig5cL-humans.png",
        "our_img": ROOT / "results" / "humans" / "figures" / "manhattan.png",
        "paper_title": "Paper — Fig 1c middle\n(Kastner et al. 2026)",
        "our_title": "pycbas — 31/408 significant (k=2)",
        "out_name": "comparison_humans.png",
    },
    {
        "name": "rats",
        "paper_img": PAPER_DIR / "ppr_fig5cR.png",
        "our_img": ROOT / "results" / "rats" / "figures" / "manhattan.png",
        "paper_title": "Paper — Fig 1c right\n(Kastner et al. 2026)",
        "our_title": "pycbas — 111/16,483 significant (k=6)",
        "out_name": "comparison_rats.png",
    },
]


def make_comparison(spec):
    if not spec["paper_img"].exists():
        print(f"  SKIP {spec['name']}: paper screenshot not found at {spec['paper_img']}")
        return
    if not spec["our_img"].exists():
        print(f"  SKIP {spec['name']}: our manhattan not found at {spec['our_img']}")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    img_paper = mpimg.imread(str(spec["paper_img"]))
    ax1.imshow(img_paper)
    ax1.set_title(spec["paper_title"], fontsize=11)
    ax1.axis("off")

    img_ours = mpimg.imread(str(spec["our_img"]))
    ax2.imshow(img_ours)
    ax2.set_title(spec["our_title"], fontsize=11)
    ax2.axis("off")

    fig.tight_layout(pad=1.5)
    out_path = OUT_DIR / spec["out_name"]
    fig.savefig(out_path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  {out_path}")


if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for spec in SPECIES:
        make_comparison(spec)
    print("\nDone.")
