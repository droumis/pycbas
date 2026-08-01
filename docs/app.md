# Interactive App

A no-code interface for running CBAS analyses. Upload your data, configure parameters, run the analysis, and explore results visually.

## Install and launch

### Option 1: pipx (recommended)

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
3. **Configure parameters** — set alphabet size, sequence length, resamples, with live resource estimates showing expected memory and runtime
4. **Run analysis** — executes the full pipeline (counting, test stats, bootstrap, step-down, k-FWER)
5. **View results** — Manhattan plot, significant sequences table, CSV export
