# Interactive App

## Try it in your browser (no install needed)

<div style="margin: 16px 0; padding: 16px; background: #fff8e1; border-left: 4px solid #f9a825; border-radius: 0 6px 6px 0;">
<strong>Browser demo</strong> runs entirely in your browser using WebAssembly. It works for small datasets (< 50 subjects, short sequences) but is slower than a local install since it can't use numba acceleration.
</div>

<a href="../demo/" class="md-button md-button--primary" target="_blank">Launch browser demo</a>

## Install locally for full performance

For real analyses with larger datasets, install and run locally. This gives you full numba acceleration and access to all your system's memory.

### Option 1: pipx (recommended for most users)

```bash
pipx install pycbas[gui]
cbas gui
```

This installs CBAS in an isolated environment and gives you a `cbas` command. The GUI opens in your browser at `localhost:5007`.

To use a different port:

```bash
cbas gui --port 5008
```

### Option 2: pip into an existing environment

```bash
pip install pycbas[gui]
cbas gui
```

### Option 3: pixi (for developers)

```bash
git clone https://github.com/droumis/pycbas.git
cd pycbas
pixi run install
pixi run gui
```

## What the app does

The GUI walks you through the full CBAS pipeline in five steps:

1. **Choose mode** — comparative (two groups) or correlative (continuous score)
2. **Load data** — upload a spreadsheet or try the built-in demo data
3. **Configure parameters** — set alphabet size, sequence length, resamples, with live resource estimates
4. **Run analysis** — executes the full pipeline (counting, test stats, bootstrap, step-down, k-FWER)
5. **View results** — Manhattan plot, significant sequences table, CSV export

The resource estimator shows expected memory usage and runtime before you commit to a run, so you know whether your configuration is feasible on your hardware.
