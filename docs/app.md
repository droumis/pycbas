# Interactive App

A no-code interface for running CBAS analyses. Upload your data, configure parameters, run the analysis, and explore results visually.

## Install and launch

### Option 1: pixi (recommended)

```bash
git clone https://github.com/droumis/pycbas.git
cd pycbas
pixi run install
pixi run gui
```

To use a different port:

```bash
pixi run -- cbas gui --port 5008
```

### Option 2: pip from the repo

```bash
git clone https://github.com/droumis/pycbas.git
cd pycbas
pip install -e ".[gui]"
cbas gui
```

### Once published to PyPI

```bash
pipx install "pycbas[gui]"
cbas gui
```

## What the app does

The GUI walks you through the full CBAS pipeline in five steps:

1. **Choose mode** — comparative (two groups) or correlative (continuous score)
2. **Load data** — upload a spreadsheet or try the built-in demo data
3. **Configure parameters** — set alphabet size, sequence length, resamples, with live resource estimates showing expected memory and runtime
4. **Run analysis** — executes the full pipeline (counting, test stats, bootstrap, step-down, k-FWER)
5. **View results** — Manhattan plot, significant sequences table, CSV export
