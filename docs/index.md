# pyCBAS

Python implementation of the [CBAS algorithm](https://github.com/dbkastner/CBAS) (Choice-Wide Behavioral Association Study) for identifying behavioral sequences that differ significantly between experimental groups or correlate with a continuous measure.

Uses Romano-Wolf step-down for multiple comparison correction and k-FWER iteration for false discovery proportion control.

## Interactive GUI

A no-code interface for running CBAS analyses. Load data, auto-detect parameters and mode, run the pipeline, and explore results visually.

![pyCBAS GUI](img/gui-results.png)

[Launch the GUI](app.md){ .md-button .md-button--primary }

## Install

We recommend installing in a dedicated environment (conda, mamba, or pixi) rather than your base environment.

```bash
git clone https://github.com/droumis/pycbas.git
cd pycbas

# option 1: pixi (recommended)
pixi install

# option 2: conda/mamba + pip
conda create -n pycbas python=3.11
conda activate pycbas
pip install -e ".[dev]"
```

To run the GUI:

```bash
# pixi
pixi run gui

# pip
pip install -e ".[gui]"
pycbas gui
```

## Minimal example

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

## Performance

| Dataset | Subjects | Sequences | Time | Peak RAM |
|---|---|---|---|---|
| Flies (2-arm, L=10) | 1,566 | 2,046 | ~21s | ~560 MB |
| Humans (6-arm, L=4) | 1,413 | 408 | ~3s | ~155 MB |
| Rats (6-arm, L=6) | 111 | 19,013 | ~8s | ~4.1 GB |

Timings on Apple M-series. Bootstrap and step-down are parallelized via numba.

## Reference

Kastner et al., "Choice-Wide Behavioral Association Study" [(2026 preprint)](https://www.biorxiv.org/content/10.1101/2024.02.26.582115v4)
