"""
CBAS Interactive GUI

A no-code interface for running Choice-Wide Behavioral Association Studies.
Launch with: panel serve app.py --show
Or install: pipx install pycbas[gui] && pycbas gui
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

import os

def _get_system_info():
    """Detect available RAM and CPU cores."""
    cores = os.cpu_count() or 1
    try:
        import psutil
        ram_gb = psutil.virtual_memory().total / (1024**3)
        available_gb = psutil.virtual_memory().available / (1024**3)
    except ImportError:
        ram_gb = None
        available_gb = None
    return {"cores": cores, "ram_gb": ram_gb, "available_gb": available_gb}

SYSTEM_INFO = _get_system_info() if not BROWSER_MODE else {"cores": 1, "ram_gb": None, "available_gb": None}

hv.extension("bokeh")
pn.extension("tabulator", sizing_mode="stretch_width", notifications=True)


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
    contingency = param.Integer(default=1, bounds=(0, 10), doc="Trial condition filter")
    block_aware = param.Boolean(default=False, doc="Sequences cannot span block/session boundaries")

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
                block_aware=self.block_aware,
            )
        else:
            self.result = run_cbas_correlative(
                self.subjects_data, self.covariate, params,
                contingency=self.contingency,
                encode_reward=self.encode_reward,
                block_aware=self.block_aware,
            )

    _observed_cache_key = param.Parameter(default=None)
    _observed_cache_val = param.Integer(default=0)

    def _count_observed_sequences(self):
        cache_key = (
            self.num_arms, self.seq_len_max, self.criterion,
            self.contingency, self.encode_reward, self.n_subjects,
        )
        if cache_key == self._observed_cache_key:
            return self._observed_cache_val

        from pycbas.io import extract_choice_stream, enumerate_sequences
        all_seqs = set()
        for subj_data in self.subjects_data:
            stream = extract_choice_stream(
                subj_data, self.contingency, self.num_arms,
                encode_reward=self.encode_reward)
            for seq_len in range(1, self.seq_len_max + 1):
                counts = enumerate_sequences(stream, seq_len, self.criterion)
                all_seqs.update(counts.keys())
        self._observed_cache_key = cache_key
        self._observed_cache_val = len(all_seqs)
        return self._observed_cache_val

    def get_resource_estimate(self):
        from pycbas import estimate_resources
        n_observed = None
        if self.data_loaded and self.subjects_data:
            n_observed = self._count_observed_sequences()
        return estimate_resources(
            num_arms=self.num_arms,
            seq_len_max=self.seq_len_max,
            n_observed=n_observed,
            resample_number=self.resample_number,
            encode_reward=self.encode_reward,
        )


app_state = CBASApp()


# =============================================================================
# UI Components
# =============================================================================

# --- Header ---
header = pn.pane.Markdown(
    """# pyCBAS
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
mode_detection_alert = pn.pane.Alert("", alert_type="info", visible=False)

def detect_mode(values):
    """Detect whether values look like group labels (comparative) or continuous scores (correlative)."""
    unique = set(values)
    if unique <= {0, 1, 0.0, 1.0}:
        return "Comparative"
    all_int = all(v == int(v) for v in values)
    if all_int and len(unique) == 2:
        return "Comparative"
    return "Correlative"

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


def set_detected_mode(values):
    """Auto-detect and set mode from info file values. Returns the detected mode string."""
    detected = detect_mode(values)
    mode_selector.value = detected
    app_state.mode = detected
    update_mode_explanation()
    mode_detection_alert.object = f"**Mode auto-detected: {detected}** based on your info file values. Change above if needed."
    mode_detection_alert.visible = True
    return detected


# --- Data upload ---
data_upload_type = pn.widgets.RadioButtonGroup(
    name="Data Format",
    options=["Local folder", "CSV per subject", "Single spreadsheet"],
    value="Local folder",
    button_type="default",
    button_style="outline",
)

# Local folder selector
folder_selector = pn.widgets.FileSelector(
    directory=str(Path.cwd()),
    only_files=False,
    file_pattern="*",
    name="Select a data folder",
)

folder_format_help = pn.pane.Markdown("""
**Navigate to your data folder** and select any file inside it (or the folder itself).
Subject files and group/score info will be detected automatically.

**Option A:** Folder with an info file
```
my_data/
├── subjectInfo.txt    ← *Info.txt: subject_id,label_or_score per line
├── subject0.txt       ← one file per subject (CSV: session,choice,reward,contingency)
└── ...
```

**Option B:** Folder with group names in filenames (no info file needed)
```
my_data/
├── control0.txt       ← prefix determines group membership
├── control1.txt
├── lesion0.txt
└── lesion1.txt
```
Recognized group keywords: control/ctrl/sham/wt (group 0), lesion/exp/ko/mutant (group 1).
""", styles={"font-size": "13px"})

folder_load_button = pn.widgets.Button(
    name="Load from this folder",
    button_type="primary",
    icon="folder-open",
    width=180,
)


def load_from_folder(event):
    selected = folder_selector.value
    if not selected:
        data_status.object = "Select something from the file browser above."
        data_status.alert_type = "warning"
        return

    # Reset previous state
    app_state.data_loaded = False
    app_state.result = None
    results_tabs.objects = []
    mode_detection_alert.visible = False

    folder_path = Path(selected[0])
    if not folder_path.is_dir():
        folder_path = folder_path.parent

    try:
        from pycbas import load_subject_data
        import re

        info_files = list(folder_path.glob("*Info.txt"))

        if info_files:
            # --- Mode 1: Info file present ---
            info_file = info_files[0]
            all_files = list(folder_path.glob("*.txt")) + list(folder_path.glob("*.csv"))
            subject_files_unsorted = [
                f for f in all_files if f != info_file and not f.name.endswith("Info.txt")
            ]

            def _numeric_key(p):
                digits = ""
                for ch in reversed(p.stem):
                    if ch.isdigit():
                        digits = ch + digits
                    else:
                        break
                return int(digits) if digits else 0

            subject_files = sorted(subject_files_unsorted, key=_numeric_key)

            # Try to parse: first detect if it has a CSV header
            info = {}
            info_df = None
            with open(info_file) as f:
                first_line = f.readline().strip()

            has_header = not first_line[0].isdigit() and "," in first_line

            if has_header:
                # Multi-column CSV with header (e.g. rat data: name,experiment,sex,genotype,lesion)
                info_df = pd.read_csv(info_file)
                # Row index = file index (row 0 -> an0.txt, etc.)
                # Look for a group column (lesion, group, label)
                group_col = None
                for col in ["lesion", "group", "label", "condition"]:
                    if col in info_df.columns:
                        group_col = col
                        break
                # Look for a filter column (genotype)
                filter_col = None
                for col in ["genotype", "geno", "genotype_filter"]:
                    if col in info_df.columns:
                        filter_col = col
                        break

                if group_col:
                    # Apply filter if present
                    if filter_col is not None:
                        mask = info_df[filter_col] == 0
                        valid_rows = info_df[mask]
                    else:
                        valid_rows = info_df

                    # Get rows with valid group labels (0 or 1)
                    valid_rows = valid_rows[valid_rows[group_col].isin([0, 1])]

                    for idx in valid_rows.index:
                        info[idx] = float(valid_rows.loc[idx, group_col])
                else:
                    # No group column found, try second column as score
                    cols = info_df.columns.tolist()
                    if len(cols) >= 2:
                        for idx, row in info_df.iterrows():
                            try:
                                info[idx] = float(row.iloc[1])
                            except (ValueError, TypeError):
                                continue
            else:
                # Simple format: subject_id, label_or_score per line
                with open(info_file) as f:
                    for line in f:
                        parts = line.strip().split(",")
                        if len(parts) >= 2:
                            try:
                                info[int(parts[0])] = float(parts[1])
                            except ValueError:
                                continue

            if not info:
                data_status.object = f"**Error:** info file `{info_file.name}` is empty or malformed."
                data_status.alert_type = "danger"
                return

            # Match subjects to files
            matched = {}
            if has_header and info_df is not None:
                # Row-indexed: row i -> file i (sorted by name)
                for subj_idx in sorted(info.keys()):
                    if subj_idx < len(subject_files):
                        matched[subj_idx] = subject_files[subj_idx]
            else:
                # ID-based matching from filename trailing digits
                for f in subject_files:
                    stem = f.stem
                    digits = ""
                    for ch in reversed(stem):
                        if ch.isdigit():
                            digits = ch + digits
                        else:
                            break
                    if digits:
                        subj_id = int(digits)
                        if subj_id in info:
                            matched[subj_id] = f

            if not matched:
                data_status.object = (
                    f"**Error:** found info file `{info_file.name}` with {len(info)} entries, "
                    "but could not match any subject data files. "
                    "Subject files should be named like `subject0.txt`, `fly12.txt`, etc."
                )
                data_status.alert_type = "danger"
                return

            subjects_data = []
            labels_or_scores = []
            for subj_id in sorted(matched.keys()):
                subjects_data.append(load_subject_data(matched[subj_id]))
                labels_or_scores.append(info[subj_id])

            source_desc = f"info file: `{info_file.name}`"

        else:
            # --- Mode 2: No info file, detect groups from filenames ---
            all_files = sorted(
                list(folder_path.glob("*.txt")) + list(folder_path.glob("*.csv"))
            )
            if not all_files:
                data_status.object = "**Error:** no data files found in that folder."
                data_status.alert_type = "danger"
                return

            # Extract prefix (letters before the trailing digits) for each file
            prefixes = {}
            for f in all_files:
                m = re.match(r'^([a-zA-Z]+)\d+$', f.stem)
                if m:
                    prefix = m.group(1)
                    prefixes.setdefault(prefix, []).append(f)

            if len(prefixes) < 2:
                data_status.object = (
                    "**Error:** no `*Info.txt` file found, and could not detect two groups "
                    "from filenames. Expected either an info file, or files named with "
                    "distinct prefixes like `control0.txt`/`lesion0.txt`."
                )
                data_status.alert_type = "danger"
                return

            # Assign group labels: sort prefixes alphabetically, group0=first, group1=second
            # If more than 2 prefixes, look for common comparative keywords
            group_keywords_0 = {"control", "ctrl", "sham", "wt", "wildtype"}
            group_keywords_1 = {"lesion", "experimental", "exp", "ko", "knockout", "mutant"}

            sorted_prefixes = sorted(prefixes.keys())
            if len(sorted_prefixes) == 2:
                grp0_prefix, grp1_prefix = sorted_prefixes
                # Swap if the second prefix looks more like a control
                if grp1_prefix.lower() in group_keywords_0 or grp0_prefix.lower() in group_keywords_1:
                    grp0_prefix, grp1_prefix = grp1_prefix, grp0_prefix
                elif grp0_prefix.lower() not in group_keywords_0 and grp1_prefix.lower() in group_keywords_1:
                    pass  # already correct
                group_map = {grp0_prefix: 0, grp1_prefix: 1}
            else:
                # Multiple prefixes: merge by keyword matching
                group_map = {}
                for p in sorted_prefixes:
                    p_lower = p.lower()
                    if any(kw in p_lower for kw in group_keywords_1):
                        group_map[p] = 1
                    else:
                        group_map[p] = 0

            subjects_data = []
            labels_or_scores = []
            for prefix in sorted_prefixes:
                label = group_map[prefix]
                for f in sorted(prefixes[prefix]):
                    subjects_data.append(load_subject_data(f))
                    labels_or_scores.append(float(label))

            group_counts = {}
            for prefix in sorted_prefixes:
                label = group_map[prefix]
                group_counts.setdefault(label, []).append(f"{prefix} ({len(prefixes[prefix])})")

            grp_desc = ", ".join(
                f"group {g}: {' + '.join(names)}"
                for g, names in sorted(group_counts.items())
            )
            source_desc = f"groups from filenames: {grp_desc}"

        # --- Common: set state and auto-detect params ---
        app_state.subjects_data = subjects_data

        set_detected_mode(labels_or_scores)
        if app_state.mode == "Comparative":
            app_state.group_labels = np.asarray(labels_or_scores, dtype=np.int32)
        else:
            app_state.covariate = np.asarray(labels_or_scores, dtype=np.float64)

        app_state.n_subjects = len(subjects_data)
        app_state.data_loaded = True

        # Auto-detect parameters from loaded data
        max_choice = max(int(arr[:, 1].max()) for arr in subjects_data)
        max_reward = max(int(arr[:, 2].max()) for arr in subjects_data)
        has_reward = max_reward > 0
        suggested_arms = max_choice + 1
        contingency_values = set()
        for arr in subjects_data:
            contingency_values.update(arr[:, 3].tolist())
        has_contingency = len(contingency_values) > 1 or (0 not in contingency_values)
        suggested_contingency = max(contingency_values) if has_contingency else 0

        trial_counts = []
        for arr in subjects_data:
            if has_contingency and suggested_contingency > 0:
                n = int((arr[:, 3] == suggested_contingency).sum())
            else:
                n = len(arr)
            trial_counts.append(n)
        suggested_criterion = min(trial_counts)

        num_arms_widget.value = suggested_arms
        app_state.num_arms = suggested_arms
        encode_reward_widget.value = has_reward
        app_state.encode_reward = has_reward
        contingency_widget.value = suggested_contingency
        app_state.contingency = suggested_contingency
        criterion_widget.value = suggested_criterion
        app_state.criterion = suggested_criterion

        # Auto-detect block_aware: enable if data has multiple blocks/sessions
        # (column 0) and contingency-filtered data spans multiple blocks
        suggested_block_aware = False
        for arr in subjects_data:
            unique_blocks = set(arr[:, 0].tolist())
            if len(unique_blocks) > 1:
                # Check if filtered data spans multiple blocks
                if suggested_contingency > 0:
                    filtered = arr[arr[:, 3] == suggested_contingency]
                else:
                    filtered = arr
                if len(filtered) > 0:
                    filtered_blocks = set(filtered[:, 0].tolist())
                    if len(filtered_blocks) > 1:
                        suggested_block_aware = True
                        break
        block_aware_widget.value = suggested_block_aware
        app_state.block_aware = suggested_block_aware

        data_status.object = (
            f"**Loaded {app_state.n_subjects} subjects** from `{folder_path.name}/` "
            f"({source_desc}). "
            f"Parameters auto-configured: {suggested_arms} arms, "
            f"{'reward encoded, ' if has_reward else ''}"
            f"criterion={suggested_criterion}"
            f"{f', contingency={suggested_contingency}' if has_contingency else ''}"
            f"{', block-aware' if suggested_block_aware else ''}."
        )
        data_status.alert_type = "success"
        update_resource_estimate()

    except Exception as e:
        data_status.object = f"**Error loading folder:** {str(e)}"
        data_status.alert_type = "danger"

folder_load_button.on_click(load_from_folder)

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

Optionally include `reward` (0 or 1) and `contingency` (trial condition integer) columns.
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

        has_group = "group" in df.columns
        has_score = "score" in df.columns
        if has_group and not has_score:
            labels = df.groupby("subject")["group"].first().loc[subjects].values
            set_detected_mode(labels.tolist())
            app_state.load_choice_streams(choice_streams, labels)
        elif has_score and not has_group:
            scores = df.groupby("subject")["score"].first().loc[subjects].values
            set_detected_mode(scores.tolist())
            app_state.load_choice_streams(choice_streams, scores)
        elif has_group and has_score:
            if app_state.mode == "Comparative":
                labels = df.groupby("subject")["group"].first().loc[subjects].values
                app_state.load_choice_streams(choice_streams, labels)
            else:
                scores = df.groupby("subject")["score"].first().loc[subjects].values
                app_state.load_choice_streams(choice_streams, scores)
        else:
            data_status.object = "**Error:** spreadsheet needs a 'group' or 'score' column."
            data_status.alert_type = "danger"
            return

        max_choice = max(s.max() for s in choice_streams)
        suggested_arms = int(max_choice) + 1
        num_arms_widget.value = suggested_arms
        app_state.num_arms = suggested_arms

        has_reward = "reward" in df.columns and df["reward"].max() > 0
        encode_reward_widget.value = has_reward
        app_state.encode_reward = has_reward

        if "contingency" in df.columns:
            cont_vals = df["contingency"].unique()
            if len(cont_vals) > 1 or (0 not in cont_vals):
                suggested_cont = int(df["contingency"].max())
                contingency_widget.value = suggested_cont
                app_state.contingency = suggested_cont

        min_len = min(len(s) for s in choice_streams)
        criterion_widget.value = min_len
        app_state.criterion = min_len

        data_status.object = (
            f"**Loaded {app_state.n_subjects} subjects** from spreadsheet. "
            f"Parameters auto-configured: {suggested_arms} arms, "
            f"{'reward encoded, ' if has_reward else ''}"
            f"criterion={min_len}."
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

        set_detected_mode(labels_values)
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
    description="Number of distinct choices in your task. With reward encoding enabled, effective alphabet becomes num_arms * 2.",
)
seq_len_max_widget = pn.widgets.IntInput(
    name="Max sequence length (L)",
    value=4, start=2, end=10, step=1,
    description="Maximum pattern length to test. All lengths 1 through L are evaluated. Larger L greatly increases the hypothesis space.",
)
criterion_widget = pn.widgets.IntInput(
    name="Trials per subject (criterion)",
    value=200, start=10, end=5000, step=10,
    description="Number of trials per subject used for sequence counting. Should not exceed the minimum trial count across subjects.",
)
resample_widget = pn.widgets.IntInput(
    name="Bootstrap resamples (M)",
    value=10000, start=100, end=50000, step=100,
    description="Number of bootstrap resamples for the null distribution. More gives tighter p-values but costs linearly in time and memory.",
)
encode_reward_widget = pn.widgets.Checkbox(
    name="Encode reward into symbols (doubles alphabet)",
    value=False,
)
encode_reward_tooltip = pn.widgets.TooltipIcon(
    value="When enabled, each trial's symbol becomes choice + reward * num_arms. Use for tasks where reward outcome is informative.",
)
contingency_widget = pn.widgets.IntInput(
    name="Contingency filter (trial condition)",
    value=1, start=0, end=10, step=1,
    description="Only trials matching this value in the contingency column are used. Filters by task condition, not session.",
)
block_aware_widget = pn.widgets.Checkbox(
    name="Block-aware (sequences cannot span sessions)",
    value=False,
)
block_aware_tooltip = pn.widgets.TooltipIcon(
    value="When enabled, sequences are counted within blocks only and cannot span block/session boundaries. Enable for multi-session experiments.",
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
        elif e.obj is block_aware_widget:
            app_state.block_aware = e.new
    update_resource_estimate()

for w in [num_arms_widget, seq_len_max_widget, criterion_widget,
          resample_widget, encode_reward_widget, contingency_widget,
          block_aware_widget]:
    w.param.watch(sync_params, "value")

resource_estimate_pane = pn.pane.Alert("", alert_type="light", visible=False)

def update_resource_estimate():
    if not app_state.data_loaded:
        resource_estimate_pane.visible = False
        return

    resource_estimate_pane.visible = True
    est = app_state.get_resource_estimate()
    alphabet = est["alphabet"]
    total_seq = est["total_sequences"]
    n_observed = est["observed_sequences"]
    mem_gb = est["memory_chunked_gb"]
    est_time = est["est_time_seconds"]

    available_gb = SYSTEM_INFO["available_gb"]
    cores = SYSTEM_INFO["cores"]

    if available_gb is not None:
        if mem_gb < available_gb * 0.3:
            verdict = "COMFORTABLE"
        elif mem_gb < available_gb * 0.7:
            verdict = "FITS (close other apps)"
        elif mem_gb < available_gb * 0.95:
            verdict = "TIGHT (may need reduced M)"
        else:
            verdict = "TOO LARGE (reduce M or sequence length)"
    else:
        if mem_gb < 1.0:
            verdict = "TRIVIAL"
        elif mem_gb < 8.0:
            verdict = "COMFORTABLE"
        elif mem_gb < 24.0:
            verdict = "FITS (close other apps)"
        else:
            verdict = "TOO LARGE (reduce M or sequence length)"

    color_map = {
        "TRIVIAL": "success",
        "COMFORTABLE": "success",
        "FITS (close other apps)": "warning",
        "TIGHT (may need reduced M)": "warning",
        "TOO LARGE (reduce M or sequence length)": "danger",
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

    sys_note = ""
    if available_gb is not None:
        sys_note = f"\n\n**System:** {available_gb:.0f} GB available, {cores} cores"

    if mem_gb < 0.1:
        mem_str = f"{mem_gb * 1024:.0f} MB"
    else:
        mem_str = f"{mem_gb:.1f} GB"

    if n_observed is not None:
        seq_str = f"**Sequences to test:** {n_observed:,} (of {total_seq:,} possible)"
    else:
        seq_str = f"**Hypothesis space:** {total_seq:,} sequences"

    resource_estimate_pane.object = (
        f"**Effective alphabet:** {alphabet} symbols | "
        f"{seq_str}\n\n"
        f"**Est. memory:** {mem_str} | "
        f"**Est. time:** ~{est_time:.0f}s | "
        f"**Verdict:** {verdict}"
        f"{browser_note}"
        f"{sys_note}"
    )
    resource_estimate_pane.alert_type = alert_type


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

    # --- Top sequences bar chart ---
    if n_sig > 0:
        top_seqs = make_top_sequences_plot(result)
        tabs.append(("Top Sequences", top_seqs))

    # --- k-convergence plot ---
    if result.k_history:
        k_plot = make_k_convergence_plot(result)
        tabs.append(("k-FWER Convergence", k_plot))

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
    unique_lens = sorted(set(seq_lengths))

    length_colors = ["#00e5ff", "#0099ff", "#0044dd", "#00aa44",
                     "#009922", "#006600", "#ccaa00", "#dd6600", "#cc0000"]
    color_map = {slen: length_colors[i % len(length_colors)]
                 for i, slen in enumerate(unique_lens)}

    # Rank sequences grouped by length (paper style)
    x_pos = np.zeros(n_seq)
    rank = 1
    for slen in unique_lens:
        for idx in np.where(seq_lengths == slen)[0]:
            x_pos[idx] = rank
            rank += 1

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
            seq_str = "-".join(str(x) for x in result.sequences[i])
            points_data.append({
                "x": x_pos[i],
                "neg_log_g": neg_log,
                "length": int(seq_lengths[i]),
                "direction": direction,
                "significant": bool(result.significant_mask[i]),
                "sequence": seq_str,
                "g_value": best_g,
                "color": color_map[int(seq_lengths[i])],
            })

    if not points_data:
        return pn.pane.Markdown("No valid g-values to plot.")

    df = pd.DataFrame(points_data)
    threshold = -np.log10(0.5)

    hover = HoverTool(tooltips=[
        ("Sequence", "@sequence"),
        ("g-value", "@g_value{0.0000}"),
        ("-log10(g)", "@neg_log_g{0.00}"),
        ("Direction", "@direction"),
        ("Length", "@length"),
    ])

    overlays = []
    for slen in unique_lens:
        slen_df = df[df["length"] == slen]
        if slen_df.empty:
            continue
        scatter = hv.Points(
            slen_df, kdims=["x", "neg_log_g"],
            vdims=["sequence", "g_value", "direction", "length", "color"],
        ).opts(color=color_map[slen], size=5, alpha=0.8, tools=[hover],
               line_color="black", line_width=0.3)
        overlays.append(scatter)

    hline = hv.HLine(threshold).opts(
        color="black", line_dash="dotted", line_width=1,
    )

    plot = hv.Overlay(overlays) * hline
    plot = plot.opts(
        width=700, height=400,
        xlabel="Sequence (ranked by length)", ylabel="-log10(g-value)",
        title="Manhattan Plot",
        logx=True,
    )

    legend_md = " | ".join(
        f'<span style="color:{color_map[slen]}">&#9679;</span> L={slen}'
        for slen in unique_lens
    )

    return pn.Column(
        pn.pane.HTML(f"<div style='text-align:center; font-size:12px;'>{legend_md}</div>"),
        pn.pane.HoloViews(plot, sizing_mode="stretch_width"),
    )


def make_top_sequences_plot(result, n_top=20):
    """Bar chart of the most significant sequences, colored by direction."""
    rows = []
    for i, seq in enumerate(result.sequences):
        if not result.significant_mask[i]:
            continue
        pos_g = result.g_values[i * 2]
        neg_g = result.g_values[i * 2 + 1]
        pos_t = result.test_stats[i * 2]
        neg_t = result.test_stats[i * 2 + 1]

        pos_sig = not np.isnan(pos_g) and pos_g < 0.5
        neg_sig = not np.isnan(neg_g) and neg_g < 0.5

        if pos_sig and neg_sig:
            if pos_g <= neg_g:
                g_val, t_val, direction = pos_g, pos_t, "positive"
            else:
                g_val, t_val, direction = neg_g, neg_t, "negative"
        elif pos_sig:
            g_val, t_val, direction = pos_g, pos_t, "positive"
        elif neg_sig:
            g_val, t_val, direction = neg_g, neg_t, "negative"
        else:
            continue

        seq_str = "-".join(str(x) for x in seq)
        rows.append({
            "sequence": seq_str,
            "t_stat": float(t_val) if not np.isnan(t_val) else 0.0,
            "g_value": g_val,
            "direction": direction,
            "length": len(seq),
        })

    if not rows:
        return pn.pane.Markdown("No significant sequences to display.")

    df = pd.DataFrame(rows)
    df = df.sort_values("g_value").head(n_top)
    # Use signed t-stat for bar direction
    df["signed_t"] = df.apply(
        lambda r: r["t_stat"] if r["direction"] == "positive" else -r["t_stat"], axis=1)
    df = df.sort_values("signed_t")

    if app_state.mode == "Correlative":
        pos_label = "Positive correlation"
        neg_label = "Negative correlation"
    else:
        pos_label = "Group 0 > Group 1"
        neg_label = "Group 1 > Group 0"

    pos_df = df[df["direction"] == "positive"]
    neg_df = df[df["direction"] == "negative"]

    bars_pos = hv.Bars(
        pos_df, kdims=["sequence"], vdims=["signed_t"],
    ).opts(color="#e76f51", alpha=0.85) if not pos_df.empty else hv.Bars([])

    bars_neg = hv.Bars(
        neg_df, kdims=["sequence"], vdims=["signed_t"],
    ).opts(color="#4361ee", alpha=0.85) if not neg_df.empty else hv.Bars([])

    plot = (bars_neg * bars_pos).opts(
        width=700, height=max(300, len(df) * 22),
        xlabel="Test statistic (t)", ylabel="",
        title=f"Top {len(df)} Significant Sequences",
        invert_axes=True,
        show_legend=False,
    )

    legend_html = (
        f'<div style="font-size:12px; text-align:center;">'
        f'<span style="color:#e76f51">&#9632;</span> {pos_label} &nbsp; '
        f'<span style="color:#4361ee">&#9632;</span> {neg_label}'
        f'</div>'
    )

    return pn.Column(
        pn.pane.HoloViews(plot, sizing_mode="stretch_width"),
        pn.pane.HTML(legend_html),
    )


def make_k_convergence_plot(result):
    history = result.k_history
    iterations = list(range(1, len(history) + 1))
    ks = [h["k"] for h in history]
    rejections = [h["rejections"] for h in history]

    k_curve = hv.Curve(
        list(zip(iterations, ks)), kdims=["Iteration"], vdims=["k"]
    ).opts(color="#4361ee", line_width=2, tools=["hover"])

    k_scatter = hv.Scatter(
        list(zip(iterations, ks)), kdims=["Iteration"], vdims=["k"]
    ).opts(color="#4361ee", size=8)

    rej_curve = hv.Curve(
        list(zip(iterations, rejections)), kdims=["Iteration"], vdims=["Rejections"]
    ).opts(color="#e76f51", line_width=2, tools=["hover"])

    rej_scatter = hv.Scatter(
        list(zip(iterations, rejections)), kdims=["Iteration"], vdims=["Rejections"]
    ).opts(color="#e76f51", size=8)

    k_plot = (k_curve * k_scatter).opts(
        width=350, height=250, ylabel="k", title="k convergence",
    )
    rej_plot = (rej_curve * rej_scatter).opts(
        width=350, height=250, ylabel="Rejections", title="Rejections per iteration",
    )

    explanation = pn.pane.Markdown(
        f"The k-FWER procedure iterates until the number of rejections stabilizes. "
        f"Converged at **k={result.k_final}** after **{len(history)}** iterations "
        f"with **{history[-1]['rejections']}** final rejections.",
        styles={"font-size": "13px"},
    )

    return pn.Column(
        explanation,
        pn.Row(
            pn.pane.HoloViews(k_plot),
            pn.pane.HoloViews(rej_plot),
        ),
    )


def make_sig_table(result):
    if result.n_significant == 0:
        return pn.pane.Markdown("No significant sequences found.")

    rows = []
    for i, seq in enumerate(result.sequences):
        if not result.significant_mask[i]:
            continue
        pos_g = result.g_values[i * 2]
        neg_g = result.g_values[i * 2 + 1]

        pos_sig = not np.isnan(pos_g) and pos_g < 0.5
        neg_sig = not np.isnan(neg_g) and neg_g < 0.5

        if pos_sig and neg_sig:
            g_val = min(pos_g, neg_g)
            if pos_g <= neg_g:
                direction = "positive" if app_state.mode == "Correlative" else "group0 > group1"
            else:
                direction = "negative" if app_state.mode == "Correlative" else "group1 > group0"
        elif pos_sig:
            direction = "positive" if app_state.mode == "Correlative" else "group0 > group1"
            g_val = pos_g
        elif neg_sig:
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

    if not rows:
        return pn.pane.Markdown("No significant sequences found.")

    df = pd.DataFrame(rows).sort_values("g-value").reset_index(drop=True)
    return pn.widgets.Tabulator(
        df, show_index=False, sizing_mode="stretch_width",
        page_size=25, pagination="remote",
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
    if data_upload_type.value == "Local folder":
        return pn.Column(
            folder_format_help,
            folder_selector,
            folder_load_button,
        )
    elif data_upload_type.value == "Single spreadsheet":
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
    make_step_indicator(1, "Load data"),
    make_step_indicator(2, "Confirm mode"),
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
    "For full performance, install locally: `pipx install pycbas[gui]` then run `pycbas gui`.",
    alert_type="warning",
) if BROWSER_MODE else None

# Main content
main_content = pn.Column(
    header,
    *([browser_banner] if browser_banner else []),
    pn.layout.Divider(),

    # Step 1: Data
    pn.pane.Markdown("## 1. Load your data"),
    pn.Row(data_upload_type, demo_button),
    data_upload_content,
    data_status,
    pn.layout.Divider(),

    # Step 2: Mode (auto-detected, with override)
    pn.pane.Markdown("## 2. Analysis mode"),
    mode_detection_alert,
    mode_selector,
    mode_explanation,
    pn.layout.Divider(),

    # Step 3: Parameters
    pn.pane.Markdown("## 3. Configure parameters"),
    pn.Row(
        pn.Column(num_arms_widget, seq_len_max_widget, pn.Row(encode_reward_widget, encode_reward_tooltip), width=300),
        pn.Column(criterion_widget, resample_widget, contingency_widget, pn.Row(block_aware_widget, block_aware_tooltip), width=300),
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
    title="pyCBAS",
    sidebar=[sidebar],
    main=[main_content],
    header_background="#4361ee",
)

template.servable()
