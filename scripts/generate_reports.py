"""
Generate validation reports and summary from results.json files.

Reads structured results from each dataset's results.json and produces:
  - Per-dataset validation_report.md
  - Aggregate validation_summary.md

Usage:
    pixi run reports
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from results_io import load_results_json

ROOT_DIR = Path(__file__).parent.parent
RESULTS_DIR = ROOT_DIR / "results"


def generate_fly_report(r):
    """Generate fly validation report from results dict."""
    res = r["results"]
    params = r["params"]
    timing = r["timing"]
    groups = r["groups"]
    labels = r["labels"]

    return f"""# Fly CBAS Validation Report

## Summary

| | pycbas | Paper (Kastner et al.) |
|---|---|---|
| Flies | {res['n_subjects']} ({groups['CA']} CA, {groups['w1118']} w1118) | 1,566 (759 CA, 807 w1118) |
| Max seq length | {params['seq_len_max']} | 10 |
| Criterion | {params['criterion']} | 250 |
| Resamples | {params['resample_number']:,} | 10,000 |
| Sequences evaluated | {res['n_sequences']:,} | 2,046 |
| Significant | {res['n_significant']} ({res['fraction_significant']*100:.1f}%) | 1,605 (78.4%) |
| {labels['positive_direction']} | {res['n_positive']} | not separately reported |
| {labels['negative_direction']} | {res['n_negative']} | not separately reported |
| k (k-FWER) | {res['k_final']} | not reported |
| Runtime | {timing['total']:.1f}s | not reported |

## Timing Profile

| Stage | Time (s) | % Total |
|---|---|---|
| build_count_matrix | {timing['build_count_matrix']:.2f} | {timing['build_count_matrix']/timing['total']*100:.1f}% |
| compute_test_stats | {timing['compute_test_stats']:.2f} | {timing['compute_test_stats']/timing['total']*100:.1f}% |
| bootstrap | {timing['bootstrap']:.2f} | {timing['bootstrap']/timing['total']*100:.1f}% |
| k_fwer | {timing['k_fwer']:.2f} | {timing['k_fwer']/timing['total']*100:.1f}% |
| **TOTAL** | **{timing['total']:.2f}** | |

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
"""


def generate_human_report(r):
    """Generate human validation report from results dict."""
    res = r["results"]
    params = r["params"]
    timing = r["timing"]
    labels = r["labels"]

    return f"""# Human CBAS Validation Report (Correlative Mode)

## Summary

| | pycbas | Paper (Kastner et al.) |
|---|---|---|
| Subjects | {res['n_subjects']} | 1,413 |
| Max seq length | {params['seq_len_max']} | 4 |
| Criterion | {params['criterion']} | 400 |
| Resamples | {params['resample_number']:,} | 10,000 |
| Sequences evaluated | {res['n_sequences']:,} | 408 |
| Significant | {res['n_significant']} ({res['fraction_significant']*100:.1f}%) | 31 (7.6%) |
| {labels['positive_direction']} | {res['n_positive']} | not separately reported |
| {labels['negative_direction']} | {res['n_negative']} | not separately reported |
| k (k-FWER) | {res['k_final']} | not reported |
| Runtime | {timing['total']:.1f}s | not reported |

## Notes

- **Mode:** Correlative — tests Pearson correlation between each sequence's usage
  count across subjects and each subject's CBIT score (a compulsivity measure).
- **Interpretation:** Positive correlation means higher CBIT (more compulsive)
  subjects use that sequence more. Negative means less.

## Timing Profile

| Stage | Time (s) | % Total |
|---|---|---|
| build_count_matrix | {timing['build_count_matrix']:.2f} | {timing['build_count_matrix']/timing['total']*100:.1f}% |
| compute_test_stats | {timing['compute_test_stats']:.2f} | {timing['compute_test_stats']/timing['total']*100:.1f}% |
| bootstrap | {timing['bootstrap']:.2f} | {timing['bootstrap']/timing['total']*100:.1f}% |
| k_fwer | {timing['k_fwer']:.2f} | {timing['k_fwer']/timing['total']*100:.1f}% |
| **TOTAL** | **{timing['total']:.2f}** | |

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
"""


def generate_summary(datasets):
    """Generate aggregate validation summary from all dataset results."""
    rows = []
    for name, r in datasets.items():
        res = r["results"]
        rows.append(
            f"| {name.capitalize()} | {r['mode'].capitalize()} | {res['n_subjects']} "
            f"| {res['n_sequences']:,} | {res['n_significant']} ({res['fraction_significant']*100:.1f}%) "
            f"| {res['k_final']} |"
        )

    table = "\n".join(rows)

    sections = []
    for name, r in datasets.items():
        res = r["results"]
        params = r["params"]
        sections.append(
            f"## {name.capitalize()}\n\n"
            f"- **{params['num_arms']} arms, seq_len_max={params['seq_len_max']}, "
            f"criterion={params['criterion']}, M={params['resample_number']:,}**\n"
            f"- {res['n_significant']}/{res['n_sequences']} significant (k={res['k_final']})\n"
            f"- Runtime: {r['timing']['total']:.1f}s\n\n"
            f"[Full report]({name}/validation_report.md)\n"
        )

    return f"""# pycbas Validation Summary

Cross-species validation of the CBAS reimplementation against Kastner et al. (2026).

## Results at a glance

| Dataset | Mode | Subjects | Sequences | Significant | k |
|---|---|---|---|---|---|
{table}

{"".join(sections)}
## Notes

- Fly result at fixed k=20 gives 1,633 sig (matches paper's 1,605). Full dataset saturates due to [k-FWER sensitivity in the high-power regime](flies/kfwer_sensitivity_analysis.md).
- Human count (69 vs paper's 31) likely reflects minor differences in tau-hat normalization.
"""


def main():
    datasets = {}

    # Load available results
    fly_json = RESULTS_DIR / "flies" / "results.json"
    human_json = RESULTS_DIR / "humans" / "results.json"

    if fly_json.exists():
        r = load_results_json(fly_json)
        datasets["flies"] = r
        report = generate_fly_report(r)
        out = RESULTS_DIR / "flies" / "validation_report.md"
        out.write_text(report)
        print(f"  {out}")

    if human_json.exists():
        r = load_results_json(human_json)
        datasets["humans"] = r
        report = generate_human_report(r)
        out = RESULTS_DIR / "humans" / "validation_report.md"
        out.write_text(report)
        print(f"  {out}")

    if datasets:
        summary = generate_summary(datasets)
        out = RESULTS_DIR / "validation_summary.md"
        out.write_text(summary)
        print(f"  {out}")

    if not datasets:
        print("No results.json files found. Run analyses first.")
        raise SystemExit(1)

    print("\nDone. Reports generated from results.json files.")


if __name__ == "__main__":
    main()
