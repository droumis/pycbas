# CBAS — Choice-Wide Behavioral Association Study

Python reimplementation of the core CBAS algorithm originally written in Igor Pro by [David Kastner](https://github.com/dkastner).

CBAS identifies behavioral sequences that differ significantly between experimental groups (comparative mode) or correlate with a continuous measure (correlative mode). It uses Romano-Wolf step-down for multiple comparison correction and k-FWER iteration for false discovery proportion control.

**Reference:** Kastner et al., "Choice-Wide Behavioral Association Study" (2026 preprint)

## Setup

Requires [pixi](https://pixi.sh):

```bash
pixi install
```

## Usage

```python
from pycbas import CBASParams, load_subject_data, run_cbas_comparative

# Load your data (CSV: session, choice, reward, contingency)
subjects_data = [load_subject_data(f) for f in data_files]
group_labels = [0, 0, 0, 1, 1, 1]  # group assignment per subject

params = CBASParams(
    num_arms=6,          # number of discrete choices
    seq_len_max=6,       # max sequence length to enumerate
    criterion=800,       # number of choices per subject to analyze
    resample_number=10000,  # bootstrap resamples
)

result = run_cbas_comparative(subjects_data, group_labels, params)
print(f"{result.n_significant} significant sequences (k={result.k_final})")
```

## Validation

Run the fast validation (~7s, reduced parameters):

```bash
pixi run validate
```

Run with paper-matched parameters (~8 min, seq_len_max=6, M=10,000, 85 subjects):

```bash
pixi run validate-paper
```

Regenerate figures from cached results (no recomputation):

```bash
pixi run figures          # from default run
pixi run figures-paper    # from paper-params run
```

## Tests

```bash
pixi run test    # 19 unit + integration tests, ~1s
```

## Performance

The step-down procedure is accelerated with numba JIT (~10x speedup). Set `NUMBA_DISABLE_JIT=1` in the environment to disable for debugging.
