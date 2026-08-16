# pycbas

[![docs](https://img.shields.io/github/actions/workflow/status/droumis/pycbas/docs.yml?style=flat-square&branch=main&label=docs&logo=materialformkdocs&logoColor=white)](https://github.com/droumis/pycbas/actions/workflows/docs.yml)
[![pypi-version](https://img.shields.io/pypi/v/pycbas.svg?logo=pypi&logoColor=white&style=flat-square)](https://pypi.org/project/pycbas)
[![python-version](https://img.shields.io/badge/python-3.11%2B-blue?logo=python&logoColor=white&style=flat-square)](https://pypi.org/project/pycbas)
[![license](https://img.shields.io/github/license/droumis/pycbas?style=flat-square&color=blue&logo=github&logoColor=white)](https://github.com/droumis/pycbas/blob/main/LICENSE)

Python implementation of the [CBAS algorithm](https://github.com/dbkastner/CBAS) (Choice-Wide Behavioral Association Study) for identifying behavioral sequences that differ significantly between experimental groups or correlate with a continuous measure.

Uses Romano-Wolf step-down for multiple comparison correction and k-FWER iteration for false discovery proportion control.

**Reference:** Kastner et al., "Choice-Wide Behavioral Association Study" [(2026 preprint)](https://www.biorxiv.org/content/10.1101/2024.02.26.582115v4)

## How it works

A sliding window walks each subject's choice stream and counts every subsequence up to length `seq_len_max`. Each unique sequence becomes one column of a subject-by-sequence count matrix, and one hypothesis test.

![Sliding window counting subsequences in a choice stream](https://raw.githubusercontent.com/droumis/pycbas/main/docs/img/concept-sequences.gif)

Testing thousands of sequences needs multiple-comparison correction, but Bonferroni's single fixed threshold is far too strict here. Romano-Wolf step-down instead recomputes the threshold from the bootstrap null after every rejection, over only the sequences that remain. The bar drops as strong effects are peeled off, so moderate effects can still clear it:

![Step-down procedure lowering the threshold after each rejection](https://raw.githubusercontent.com/droumis/pycbas/main/docs/img/concept-stepdown.gif)

k-FWER iteration then relaxes "no false positives" to "at most k", raising k until the false discovery proportion is bounded by `gamma`.

Both animations are interactive in the [walkthrough](https://droumis.github.io/pycbas/walkthrough/), which builds up the whole algorithm step by step.

## Installation

```bash
pip install pycbas
```

For the interactive GUI:

```bash
pip install 'pycbas[gui]'
pycbas gui
```

Load data, confirm the auto-detected mode, set parameters, run, and explore results — no code required.

<img src="https://raw.githubusercontent.com/droumis/pycbas/main/docs/img/gui-overview.png" alt="pyCBAS GUI running the human dataset" width="760">

See the [GUI documentation](https://droumis.github.io/pycbas/app/) for details.

### Development install

We recommend installing in a dedicated environment (conda, mamba, or pixi) rather than your base environment.

```bash
git clone https://github.com/droumis/pycbas.git
cd pycbas

# option 1: pixi (handles everything)
pixi install

# option 2: conda/mamba + pip
conda create -n pycbas python=3.11
conda activate pycbas
pip install -e '.[dev]'
```

## Quick start

### Comparative mode (group differences)

```python
from pycbas import CBASParams, load_subject_data, run_cbas_comparative

subjects_data = [load_subject_data(f) for f in data_files]
group_labels = [0, 0, 0, 1, 1, 1]

params = CBASParams(
    num_arms=6,
    seq_len_max=6,
    criterion=800,
    resample_number=10000,
)

result = run_cbas_comparative(subjects_data, group_labels, params)
print(f"{result.n_significant} significant sequences (k={result.k_final})")
```

### Correlative mode (continuous covariate)

```python
from pycbas import run_cbas_correlative

result = run_cbas_correlative(subjects_data, cbit_scores, params)
```

### Resource estimation

```python
from pycbas import estimate_resources, print_resource_estimate

est = estimate_resources(num_arms=12, seq_len_max=8, n_observed=5000)
print_resource_estimate(est)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `num_arms` | 6 | Number of base symbols (choices) |
| `seq_len_max` | 6 | Maximum sequence length L |
| `criterion` | 800 | Number of trials used per subject |
| `resample_number` | 10,000 | Bootstrap resamples M |
| `alpha` | 0.5 | Significance threshold for FDP control |
| `gamma` | 0.05 | FDP tolerance |
| `centering` | False | Center bootstrap null (False matches Igor) |
| `block_aware` | False | Prevent sequences from spanning block/session boundaries |

## Performance

| Dataset | Subjects | Sequences | Time | Peak RAM |
|---|---|---|---|---|
| Flies (2-arm, L=10) | 1,566 | 2,046 | ~21s | ~560 MB |
| Humans (6-arm, L=4) | 1,413 | 408 | ~3s | ~155 MB |
| Rats (6-arm, L=6) | 105 | 16,378 | ~11s | ~3.6 GB |

Timings on Apple M-series. The chunked pipeline (`chunked=True`, default) trades ~30% more time for ~40% less memory. Bootstrap and step-down are parallelized via numba. Set `NUMBA_DISABLE_JIT=1` to disable for debugging.

## Validation

Exact match with the original Igor implementation on flies (1,605/2,046, k=81) and humans (31/408, k=2). Test statistics match to floating-point precision. Rats (105 subjects, `block_aware=True`): 572/16,378 significant (k=29), exact match with David's Igor implementation. Test statistics agree within 1e-6 on all 16,376 overlapping sequences.

See [results/validation_summary.md](results/validation_summary.md) for details, or per-dataset reports:
- [Flies](results/flies/validation_report.md)
- [Humans](results/humans/validation_report.md)
- [Rats](results/rats/validation_report.md)

## Documentation

Full docs at **[droumis.github.io/pycbas](https://droumis.github.io/pycbas/)**

- [User Guide](https://droumis.github.io/pycbas/guide/) - data format, parameter selection, working with results
- [Algorithm](https://droumis.github.io/pycbas/algorithm/) - the step-down and k-FWER procedure in detail
- [API Reference](https://droumis.github.io/pycbas/api/) - all public functions and classes

## Development

```bash
pixi install          # set up environment
pixi run test         # run tests
pixi run flies        # run fly analysis (paper params)
pixi run human        # run human analysis
pixi run rats         # run rat analysis
```

## License

MIT
