"""
CBAS Interactive GUI

A no-code interface for running Choice-Wide Behavioral Association Studies.
Launch with: panel serve app.py --show
Or install: pipx install pycbas[gui] && cbas gui
"""

import io
import sys
import time
import threading
from pathlib import Path

import numpy as np
import pandas as pd
import panel as pn
import param
import holoviews as hv
from bokeh.models import HoverTool

BROWSER_MODE = hasattr(sys, "_emscripte_info") or "pyodide" in sys.modules

hv.extension("bokeh")
pn.extension(sizing_mode="stretch_width", notifications=True)


# =============================================================================
# State and logic
# =============================================================================

class CBASApp(param.Parameterized):
    # --- Mode ---
    mode = param.Selector(
        default="Comparative",
        objects=["Comparative", "Correlative"],
        doc="Analysis mode",
    )

    # --- Data ---
    subjects_data = param.List(default=[], doc="Loaded subject arrays")
    group_labels = param.Array(default=np.array([]), doc="Group labels (comparative)")
    covariate = param.Array(default=np.array([]), doc="Covariate scores (correlative)")
    n_subjects = param.Integer(default=0)
    data_loaded = param.Boolean(default=False)

    # --- Parameters ---
    num_arms = param.Integer(default=2, bounds=(2, 20), doc="Number of choice symbols")
    seq_len_max = param.Integer(default=4, bounds=(2, 10), doc="Max sequence length")
    criterion = param.Integer(default=200, bounds=(10, 5000), doc="Trials per subject to use")
    resample_number = param.Integer(default=10000, bounds=(100, 50000), doc="Bootstrap resamples")
    encode_reward = param.Boolean(default=False, doc="Encode reward into symbols")
    contingency = param.Integer(default=1, bounds=(0, 10), doc="Block type filter")

    # --- Run state ---
    running = param.Boolean(default=False)
    progress_value = param.Integer(default=0)
    status_text = param.String(default="")
    result = param.Parameter(default=None)

    def load_choice_streams(self, choice_data, labels_or_scores):
        """Load from arrays directly (choice streams + labels/scores)."""
        self.subjects_data = []
        for stream in choice_data:
            stream = np.asarray(stream, dtype=np.int32)
            n = len(stream)
            arr = np.zeros((n, 4), dtype=np.int32)
            arr[:, 1] = stream
            arr[:, 3] = self.contingency
            self.subjects_data.append(arr)

        if self.mode == "Comparative":
            self.group_labels = np.asarray(labels_or_scores, dtype=np.int32)
        else:
            self.covariate = np.asarray(labels_or_scores, dtype=np.float64)

        self.n_subjects = len(self.subjects_data)
        self.data_loaded = True

    def load_csv_files(self, file_contents_list, labels_or_scores):
        """Load from uploaded CSV file contents."""
        from pycbas import load_subject_data
        self.subjects_data = []
        for content in file_contents_list:
            rows = []
            for line in content.strip().split("\n"):
                parts = line.strip().split(",")
                if len(parts) >= 4:
                    rows.append([int(p) if p.strip() else 0 for p in parts[:4]])
            if rows:
                self.subjects_data.append(np.array(rows, dtype=np.int32))

        if self.mode == "Comparative":
            self.group_labels = np.asarray(labels_or_scores, dtype=np.int32)
        else:
            self.covariate = np.asarray(labels_or_scores, dtype=np.float64)

        self.n_subjects = len(self.subjects_data)
        self.data_loaded = True

    def run_analysis(self):
        """Run the full CBAS pipeline."""
        from pycbas import (
            CBASParams, run_cbas_comparative, run_cbas_correlative,
        )
        params = CBASParams(
            num_arms=self.num_arms,
            seq_len_max=self.seq_len_max,
            criterion=self.criterion,
            resample_number=self.resample_number,
        )
        if self.mode == "Comparative":
            self.result = run_cbas_comparative(
                self.subjects_data, self.group_labels, params,
                contingency=self.contingency,
                encode_reward=self.encode_reward,
            )
        else:
            self.result = run_cbas_correlative(
                self.subjects_data, self.covariate, params,
            )

    def get_resource_estimate(self):
        from pycbas import estimate_resources
        return estimate_resources(
            num_arms=self.num_arms,
            seq_len_max=self.seq_len_max,
            resample_number=self.resample_number,
            encode_reward=self.encode_reward,
        )


app_state = CBASApp()


# =============================================================================
# UI Components
# =============================================================================

# --- Header ---
header = pn.pane.Markdown(
    """# CBAS
### Choice-Wide Behavioral Association Study

Find behavioral sequences that differ between groups or correlate with a continuous measure.
""",
    styles={"margin-bottom": "10px"},
)

# --- Step indicators ---
def make_step_indicator(number, title, active=False):
    color = "#4361ee" if active else "#ccc"
    return pn.pane.HTML(f"""
        <div style="display:flex; align-items:center; gap:8px; margin:4px 0;">
            <div style="width:28px; height:28px; border-radius:50%; background:{color};
                        color:white; display:flex; align-items:center; justify-content:center;
                        font-size:14px; font-weight:600;">{number}</div>
            <span style="font-size:14px; color:{'#333' if active else '#999'}; font-weight:{'600' if active else '400'};">{title}</span>
        </div>
    """)


# --- Mode selection ---
mode_selector = pn.widgets.RadioButtonGroup(
    name="Analysis Mode",
    options=["Comparative", "Correlative"],
    value="Comparative",
    button_type="primary",
    button_style="outline",
)

mode_explanation = pn.pane.Markdown("", styles={"font-size": "13px", "color": "#666"})

def update_mode_explanation(event=None):
    value = event.new if event else mode_selector.value
    app_state.mode = value
    if value == "Comparative":
        mode_explanation.object = (
            "**Comparative mode** tests whether sequences appear at different rates "
            "between two groups of subjects (e.g. control vs. experimental). "
            "You will provide group labels (0 or 1) for each subject."
        )
    else:
        mode_explanation.object = (
            "**Correlative mode** tests whether sequence usage correlates with a "
            "continuous measure across subjects (e.g. a clinical score, age, or "
            "performance metric). You will provide one numeric score per subject."
        )

mode_selector.param.watch(update_mode_explanation, "value")
update_mode_explanation()


# --- Data upload ---
data_upload_type = pn.widgets.RadioButtonGroup(
    name="Data Format",
    options=["CSV per subject", "Single spreadsheet"],
    value="Single spreadsheet",
    button_type="default",
    button_style="outline",
)

# Multi-file upload for CSVs
csv_file_input = pn.widgets.FileInput(
    accept=".csv,.txt",
    multiple=True,
    name="Upload subject files",
)

# Single spreadsheet upload
spreadsheet_input = pn.widgets.FileInput(
    accept=".csv,.tsv,.xlsx",
    multiple=False,
    name="Upload data spreadsheet",
)

spreadsheet_format_help = pn.pane.Markdown("""
**Expected format:** A table where each row is one trial and columns include:

| Column | Description |
|--------|-------------|
| `subject` | Subject identifier (string or number) |
| `choice` | Integer choice/symbol (0-indexed) |

**For comparative mode**, also include:

| Column | Description |
|--------|-------------|
| `group` | Group label (0 or 1) |

**For correlative mode**, also include:

| Column | Description |
|--------|-------------|
| `score` | Continuous measure per subject |

Optionally include `reward` (0 or 1) and `contingency` (block type integer) columns.
Trials should be in order within each subject.

**Example (comparative):**
```
subject,choice,group
s01,0,0
s01,1,0
s01,0,0
s02,1,1
s02,0,1
```
""", styles={"font-size": "13px"})

csv_format_help = pn.pane.Markdown("""
**Expected format:** One CSV file per subject, each with columns:
```
choice
0
1
0
1
```
Or with optional columns: `session,choice,reward,contingency`

You will also need to upload a separate labels/scores file (one value per line,
matching the order of your subject files).
""", styles={"font-size": "13px"})

labels_file_input = pn.widgets.FileInput(
    accept=".csv,.txt",
    multiple=False,
    name="Upload labels/scores file",
)

# Data loading status
data_status = pn.pane.Alert("No data loaded yet.", alert_type="light")

def parse_spreadsheet(event):
    if spreadsheet_input.value is None:
        return
    try:
        content = spreadsheet_input.value
        filename = spreadsheet_input.filename or ""

        if filename.endswith(".xlsx"):
            df = pd.read_excel(io.BytesIO(content))
        elif filename.endswith(".tsv"):
            df = pd.read_csv(io.BytesIO(content), sep="\t")
        else:
            df = pd.read_csv(io.BytesIO(content))

        df.columns = df.columns.str.strip().str.lower()

        if "subject" not in df.columns or "choice" not in df.columns:
            data_status.object = "**Error:** spreadsheet must have 'subject' and 'choice' columns."
            data_status.alert_type = "danger"
            return

        subjects = df["subject"].unique()
        choice_streams = []
        for subj in subjects:
            subj_df = df[df["subject"] == subj]
            stream = subj_df["choice"].values.astype(np.int32)
            choice_streams.append(stream)

        if app_state.mode == "Comparative":
            if "group" not in df.columns:
                data_status.object = "**Error:** comparative mode requires a 'group' column."
                data_status.alert_type = "danger"
                return
            labels = df.groupby("subject")["group"].first().loc[subjects].values
            app_state.load_choice_streams(choice_streams, labels)
        else:
            if "score" not in df.columns:
                data_status.object = "**Error:** correlative mode requires a 'score' column."
                data_status.alert_type = "danger"
                return
            scores = df.groupby("subject")["score"].first().loc[subjects].values
            app_state.load_choice_streams(choice_streams, scores)

        max_choice = max(s.max() for s in choice_streams)
        suggested_arms = int(max_choice) + 1
        if suggested_arms > app_state.num_arms:
            app_state.num_arms = suggested_arms
            num_arms_widget.value = suggested_arms

        min_len = min(len(s) for s in choice_streams)
        if app_state.criterion > min_len:
            app_state.criterion = min_len
            criterion_widget.value = min_len

        data_status.object = (
            f"**Loaded {app_state.n_subjects} subjects** from spreadsheet. "
            f"Choices range 0-{max_choice}, "
            f"min trials per subject: {min_len}."
        )
        data_status.alert_type = "success"
        update_resource_estimate()

    except Exception as e:
        data_status.object = f"**Error parsing file:** {str(e)}"
        data_status.alert_type = "danger"

spreadsheet_input.param.watch(parse_spreadsheet, "value")


def parse_csv_files(event):
    if csv_file_input.value is None or labels_file_input.value is None:
        if csv_file_input.value is not None:
            data_status.object = "Subject files received. Now upload the labels/scores file."
            data_status.alert_type = "warning"
        return

    try:
        file_contents = []
        if isinstance(csv_file_input.value, list):
            for content in csv_file_input.value:
                file_contents.append(content.decode("utf-8") if isinstance(content, bytes) else content)
        else:
            file_contents.append(
                csv_file_input.value.decode("utf-8")
                if isinstance(csv_file_input.value, bytes)
                else csv_file_input.value
            )

        labels_content = (
            labels_file_input.value.decode("utf-8")
            if isinstance(labels_file_input.value, bytes)
            else labels_file_input.value
        )
        labels_values = []
        for line in labels_content.strip().split("\n"):
            line = line.strip()
            if line:
                labels_values.append(float(line.split(",")[-1]))

        app_state.load_csv_files(file_contents, labels_values)
        data_status.object = f"**Loaded {app_state.n_subjects} subjects** from CSV files."
        data_status.alert_type = "success"
        update_resource_estimate()

    except Exception as e:
        data_status.object = f"**Error:** {str(e)}"
        data_status.alert_type = "danger"

csv_file_input.param.watch(parse_csv_files, "value")
labels_file_input.param.watch(parse_csv_files, "value")


# --- Parameters ---
num_arms_widget = pn.widgets.IntInput(
    name="Number of choice symbols",
    value=2, start=2, end=20, step=1,
)
seq_len_max_widget = pn.widgets.IntInput(
    name="Max sequence length (L)",
    value=4, start=2, end=10, step=1,
)
criterion_widget = pn.widgets.IntInput(
    name="Trials per subject (criterion)",
    value=200, start=10, end=5000, step=10,
)
resample_widget = pn.widgets.IntInput(
    name="Bootstrap resamples (M)",
    value=10000, start=100, end=50000, step=100,
)
encode_reward_widget = pn.widgets.Checkbox(
    name="Encode reward into symbols (doubles alphabet)",
    value=False,
)
contingency_widget = pn.widgets.IntInput(
    name="Contingency filter (block type)",
    value=1, start=0, end=10, step=1,
)

def sync_params(*events):
    for e in events:
        if e.obj is num_arms_widget:
            app_state.num_arms = e.new
        elif e.obj is seq_len_max_widget:
            app_state.seq_len_max = e.new
        elif e.obj is criterion_widget:
            app_state.criterion = e.new
        elif e.obj is resample_widget:
            app_state.resample_number = e.new
        elif e.obj is encode_reward_widget:
            app_state.encode_reward = e.new
        elif e.obj is contingency_widget:
            app_state.contingency = e.new
    update_resource_estimate()

for w in [num_arms_widget, seq_len_max_widget, criterion_widget,
          resample_widget, encode_reward_widget, contingency_widget]:
    w.param.watch(sync_params, "value")

resource_estimate_pane = pn.pane.Alert("", alert_type="light")

def update_resource_estimate():
    est = app_state.get_resource_estimate()
    alphabet = est["alphabet"]
    total_seq = est["total_sequences"]
    mem_gb = est["memory_chunked_gb"]
    est_time = est["est_time_seconds"]
    verdict = est["recommendation"]

    color_map = {
        "TRIVIAL": "success",
        "COMFORTABLE": "success",
        "FITS (close other apps)": "warning",
        "TIGHT (may need reduced M or chunking to disk)": "warning",
        "TOO LARGE (needs embedding or M reduction)": "danger",
    }
    alert_type = color_map.get(verdict, "light")

    browser_note = ""
    if BROWSER_MODE and est_time > 30:
        browser_note = (
            "\n\n**Warning:** this configuration may be too large for browser mode. "
            "Consider reducing resamples (M) or max sequence length, "
            "or install locally for full performance."
        )
        alert_type = "danger"

    resource_estimate_pane.object = (
        f"**Effective alphabet:** {alphabet} symbols | "
        f"**Hypothesis space:** {total_seq:,} sequences\n\n"
        f"**Est. memory:** {mem_gb:.1f} GB | "
        f"**Est. time:** {est_time:.0f}s | "
        f"**Verdict:** {verdict}"
        f"{browser_note}"
    )
    resource_estimate_pane.alert_type = alert_type

update_resource_estimate()


# --- Run button and progress ---
run_button = pn.widgets.Button(
    name="Run CBAS Analysis",
    button_type="primary",
    icon="play",
    disabled=True,
    width=250,
    height=45,
)
progress_bar = pn.indicators.Progress(
    name="Progress", value=0, max=100, active=False,
    bar_color="primary", width=400,
)
run_status = pn.pane.Markdown("", styles={"font-size": "13px"})

def check_ready(*events):
    run_button.disabled = not app_state.data_loaded

app_state.param.watch(check_ready, "data_loaded")

def on_run_click(event):
    if not app_state.data_loaded:
        pn.state.notifications.error("Load data first.")
        return

    run_button.disabled = True
    run_button.name = "Running..."
    progress_bar.active = True
    progress_bar.value = -1
    run_status.object = "Starting analysis..."

    if BROWSER_MODE:
        try:
            t0 = time.perf_counter()
            app_state.run_analysis()
            elapsed = time.perf_counter() - t0
            finish_run(elapsed)
        except Exception as e:
            fail_run(str(e))
    else:
        def do_run():
            try:
                t0 = time.perf_counter()
                app_state.run_analysis()
                elapsed = time.perf_counter() - t0
                pn.state.execute(lambda: finish_run(elapsed))
            except Exception as e:
                pn.state.execute(lambda: fail_run(str(e)))

        thread = threading.Thread(target=do_run, daemon=True)
        thread.start()

def finish_run(elapsed):
    progress_bar.active = False
    progress_bar.value = 100
    run_button.disabled = False
    run_button.name = "Run CBAS Analysis"
    result = app_state.result
    run_status.object = (
        f"**Done in {elapsed:.1f}s.** "
        f"Found **{result.n_significant}** significant sequences "
        f"(k={result.k_final}) out of {len(result.sequences):,} tested."
    )
    results_tabs.objects = build_results_tabs()

def fail_run(error_msg):
    progress_bar.active = False
    progress_bar.value = 0
    run_button.disabled = False
    run_button.name = "Run CBAS Analysis"
    run_status.object = f"**Error:** {error_msg}"
    pn.state.notifications.error(f"Analysis failed: {error_msg}")

run_button.on_click(on_run_click)


# --- Results ---
results_tabs = pn.Column()

def build_results_tabs():
    result = app_state.result
    if result is None:
        return []

    tabs = pn.Tabs(dynamic=True)

    # --- Summary tab ---
    n_seq = len(result.sequences)
    n_sig = result.n_significant
    summary_md = f"""
## Results Summary

| Metric | Value |
|--------|-------|
| Subjects | {app_state.n_subjects} |
| Sequences tested | {n_seq:,} |
| Significant sequences | {n_sig} ({n_sig/n_seq*100:.1f}%) |
| k (k-FWER) | {result.k_final} |
| Mode | {app_state.mode} |
"""
    tabs.append(("Summary", pn.pane.Markdown(summary_md)))

    # --- Manhattan plot tab ---
    manhattan = make_manhattan_plot(result)
    tabs.append(("Manhattan Plot", manhattan))

    # --- Significant sequences table ---
    sig_table = make_sig_table(result)
    tabs.append(("Significant Sequences", sig_table))

    # --- Download ---
    download_section = make_download_section(result)
    tabs.append(("Export", download_section))

    return [tabs]


def make_manhattan_plot(result):
    n_seq = len(result.sequences)
    g_values = result.g_values
    seq_lengths = np.array([len(s) for s in result.sequences])

    points_data = []
    for i in range(n_seq):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        best_g = np.nan
        direction = ""
        if not np.isnan(pos_g) and not np.isnan(neg_g):
            if pos_g <= neg_g:
                best_g, direction = pos_g, "positive"
            else:
                best_g, direction = neg_g, "negative"
        elif not np.isnan(pos_g):
            best_g, direction = pos_g, "positive"
        elif not np.isnan(neg_g):
            best_g, direction = neg_g, "negative"

        if not np.isnan(best_g) and best_g > 0:
            neg_log = -np.log10(best_g)
            seq_str = "".join(str(x) for x in result.sequences[i])
            points_data.append({
                "x": i,
                "neg_log_g": neg_log,
                "length": int(seq_lengths[i]),
                "direction": direction,
                "significant": bool(result.significant_mask[i]),
                "sequence": seq_str,
                "g_value": best_g,
            })

    if not points_data:
        return pn.pane.Markdown("No valid g-values to plot.")

    df = pd.DataFrame(points_data)
    threshold = -np.log10(0.5)

    sig_df = df[df["significant"]]
    nonsig_df = df[~df["significant"]]

    hover = HoverTool(tooltips=[
        ("Sequence", "@sequence"),
        ("g-value", "@g_value{0.0000}"),
        ("-log10(g)", "@neg_log_g{0.00}"),
        ("Direction", "@direction"),
        ("Length", "@length"),
    ])

    scatter_nonsig = hv.Points(
        nonsig_df, kdims=["x", "neg_log_g"],
        vdims=["sequence", "g_value", "direction", "length"],
    ).opts(color="#ccc", size=4, alpha=0.5, tools=[hover])

    scatter_sig = hv.Points(
        sig_df, kdims=["x", "neg_log_g"],
        vdims=["sequence", "g_value", "direction", "length"],
    ).opts(color="direction", cmap={"positive": "#e76f51", "negative": "#4361ee"},
           size=7, alpha=0.8, tools=[hover])

    hline = hv.HLine(threshold).opts(
        color="red", line_dash="dashed", line_width=1,
    )

    plot = (scatter_nonsig * scatter_sig * hline).opts(
        width=700, height=400,
        xlabel="Sequence index", ylabel="-log10(g-value)",
        title="Manhattan Plot",
        show_legend=True,
    )
    return pn.pane.HoloViews(plot, sizing_mode="stretch_width")


def make_sig_table(result):
    if result.n_significant == 0:
        return pn.pane.Markdown("No significant sequences found.")

    rows = []
    for i, seq in enumerate(result.sequences):
        if not result.significant_mask[i]:
            continue
        pos_g = result.g_values[i * 2]
        neg_g = result.g_values[i * 2 + 1]
        if not np.isnan(pos_g) and pos_g < 0.5:
            direction = "positive" if app_state.mode == "Correlative" else "group0 > group1"
            g_val = pos_g
        elif not np.isnan(neg_g) and neg_g < 0.5:
            direction = "negative" if app_state.mode == "Correlative" else "group1 > group0"
            g_val = neg_g
        else:
            continue

        rows.append({
            "Sequence": " → ".join(str(x) for x in seq),
            "Length": len(seq),
            "Direction": direction,
            "g-value": round(g_val, 6),
        })

    df = pd.DataFrame(rows).sort_values("g-value")
    return pn.widgets.Tabulator(
        df, show_index=False, sizing_mode="stretch_width",
        page_size=25, pagination="remote",
        frozen_columns=["Sequence"],
    )


def make_download_section(result):
    # Build a downloadable CSV of all results
    rows = []
    for i, seq in enumerate(result.sequences):
        pos_g = result.g_values[i * 2]
        neg_g = result.g_values[i * 2 + 1]
        pos_t = result.test_stats[i * 2]
        neg_t = result.test_stats[i * 2 + 1]
        rows.append({
            "sequence": "-".join(str(x) for x in seq),
            "length": len(seq),
            "t_positive": pos_t if not np.isnan(pos_t) else "",
            "t_negative": neg_t if not np.isnan(neg_t) else "",
            "g_positive": pos_g if not np.isnan(pos_g) else "",
            "g_negative": neg_g if not np.isnan(neg_g) else "",
            "significant": result.significant_mask[i],
        })

    df = pd.DataFrame(rows)
    csv_buffer = io.StringIO()
    df.to_csv(csv_buffer, index=False)

    download_btn = pn.widgets.FileDownload(
        callback=lambda: io.StringIO(csv_buffer.getvalue()),
        filename="cbas_results.csv",
        button_type="success",
        label="Download full results (CSV)",
    )

    sig_df = df[df["significant"] == True]
    sig_buffer = io.StringIO()
    sig_df.to_csv(sig_buffer, index=False)

    download_sig_btn = pn.widgets.FileDownload(
        callback=lambda: io.StringIO(sig_buffer.getvalue()),
        filename="cbas_significant.csv",
        button_type="primary",
        label="Download significant only (CSV)",
    )

    return pn.Column(
        pn.pane.Markdown(f"**{len(df)} total sequences**, {len(sig_df)} significant."),
        pn.Row(download_btn, download_sig_btn),
        pn.pane.Markdown(
            "The CSV contains all sequences with their test statistics, "
            "g-values (adjusted p-values), and significance status."
        ),
    )


# --- Demo data generator ---
def generate_demo_data(event):
    """Generate synthetic data for testing."""
    np.random.seed(42)
    n_per_group = 25
    n_trials = 200

    if app_state.mode == "Comparative":
        choice_streams = []
        # Group 0: slight bias toward alternation (0→1→0→1)
        for _ in range(n_per_group):
            stream = np.zeros(n_trials, dtype=np.int32)
            for t in range(1, n_trials):
                if np.random.rand() < 0.65:
                    stream[t] = 1 - stream[t - 1]
                else:
                    stream[t] = stream[t - 1]
            choice_streams.append(stream)
        # Group 1: slight bias toward perseveration
        for _ in range(n_per_group):
            stream = np.zeros(n_trials, dtype=np.int32)
            for t in range(1, n_trials):
                if np.random.rand() < 0.35:
                    stream[t] = 1 - stream[t - 1]
                else:
                    stream[t] = stream[t - 1]
            choice_streams.append(stream)

        labels = [0] * n_per_group + [1] * n_per_group
        app_state.load_choice_streams(choice_streams, labels)
    else:
        n_subjects = 80
        scores = np.random.normal(50, 15, n_subjects)
        choice_streams = []
        for score in scores:
            stream = np.zeros(n_trials, dtype=np.int32)
            alt_prob = 0.3 + 0.008 * (score - 50)
            alt_prob = np.clip(alt_prob, 0.1, 0.9)
            for t in range(1, n_trials):
                if np.random.rand() < alt_prob:
                    stream[t] = 1 - stream[t - 1]
                else:
                    stream[t] = stream[t - 1]
            choice_streams.append(stream)
        app_state.load_choice_streams(choice_streams, scores)

    data_status.object = (
        f"**Demo data loaded:** {app_state.n_subjects} subjects, "
        f"{n_trials} trials each, {app_state.mode.lower()} mode. "
        "This synthetic dataset has an embedded alternation bias to demonstrate the pipeline."
    )
    data_status.alert_type = "info"
    update_resource_estimate()

demo_button = pn.widgets.Button(
    name="Load demo data",
    button_type="light",
    icon="flask",
    width=160,
)
demo_button.on_click(generate_demo_data)


# =============================================================================
# Layout
# =============================================================================

def data_upload_panel():
    """Build the data upload section based on selected format."""
    if data_upload_type.value == "Single spreadsheet":
        return pn.Column(
            spreadsheet_format_help,
            spreadsheet_input,
        )
    else:
        return pn.Column(
            csv_format_help,
            csv_file_input,
            pn.pane.Markdown("**Labels/scores file** (one value per line, same order as subject files):"),
            labels_file_input,
        )

data_upload_content = pn.Column(data_upload_panel())

def on_upload_type_change(event):
    data_upload_content.objects = [data_upload_panel()]

data_upload_type.param.watch(on_upload_type_change, "value")


# Assemble the sidebar
sidebar = pn.Column(
    pn.pane.Markdown("### Steps", styles={"margin-top": "0"}),
    make_step_indicator(1, "Choose mode"),
    make_step_indicator(2, "Load data"),
    make_step_indicator(3, "Set parameters"),
    make_step_indicator(4, "Run analysis"),
    make_step_indicator(5, "View results"),
    pn.layout.Divider(),
    pn.pane.Markdown(
        "**Tip:** Start with demo data to see how it works before loading your own.",
        styles={"font-size": "12px", "color": "#888"},
    ),
    width=220,
    styles={"padding": "20px 16px", "background": "#f8f9fa", "border-radius": "8px"},
)


browser_banner = pn.pane.Alert(
    "**Browser demo mode.** Running without numba acceleration. "
    "Suitable for small datasets (< 50 subjects, short sequences). "
    "For full performance, install locally: `pipx install pycbas[gui]` then run `cbas gui`.",
    alert_type="warning",
) if BROWSER_MODE else None

# Main content
main_content = pn.Column(
    header,
    *([browser_banner] if browser_banner else []),
    pn.layout.Divider(),

    # Step 1: Mode
    pn.pane.Markdown("## 1. Choose analysis mode"),
    mode_selector,
    mode_explanation,
    pn.layout.Divider(),

    # Step 2: Data
    pn.pane.Markdown("## 2. Load your data"),
    pn.Row(data_upload_type, demo_button),
    data_upload_content,
    data_status,
    pn.layout.Divider(),

    # Step 3: Parameters
    pn.pane.Markdown("## 3. Configure parameters"),
    pn.Row(
        pn.Column(num_arms_widget, seq_len_max_widget, encode_reward_widget, width=300),
        pn.Column(criterion_widget, resample_widget, contingency_widget, width=300),
    ),
    resource_estimate_pane,
    pn.layout.Divider(),

    # Step 4: Run
    pn.pane.Markdown("## 4. Run analysis"),
    pn.Row(run_button, progress_bar),
    run_status,
    pn.layout.Divider(),

    # Step 5: Results
    pn.pane.Markdown("## 5. Results"),
    results_tabs,

    sizing_mode="stretch_width",
    styles={"max-width": "900px", "padding": "20px"},
)


# Final layout
template = pn.template.MaterialTemplate(
    title="CBAS",
    sidebar=[sidebar],
    main=[main_content],
    header_background="#4361ee",
)

template.servable()
