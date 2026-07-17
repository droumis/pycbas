"""
Shared I/O for CBAS analysis results.

Each analysis script writes:
  - results.json: summary numbers, params, timing (human-readable, used by reports/figures)
  - results.npz: full arrays (g_values, test_stats, etc.) for replotting

The JSON is the single source of truth for all numbers that appear in
figures, reports, and the validation summary.
"""

import json
import numpy as np
from pathlib import Path


def save_results_json(path, data):
    """Save results dict as JSON. Handles numpy types."""
    def convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=convert)


def load_results_json(path):
    """Load results JSON."""
    with open(path) as f:
        return json.load(f)


def compute_significance_summary(g_values, n_sequences, alpha=0.5):
    """Compute significance counts from g_values array."""
    n_sig = 0
    n_pos = 0
    n_neg = 0
    for i in range(n_sequences):
        pos_g = g_values[i * 2]
        neg_g = g_values[i * 2 + 1]
        pos_sig = not np.isnan(pos_g) and pos_g < alpha
        neg_sig = not np.isnan(neg_g) and neg_g < alpha
        if pos_sig or neg_sig:
            n_sig += 1
        if pos_sig:
            n_pos += 1
        if neg_sig:
            n_neg += 1
    return {
        "n_significant": n_sig,
        "n_positive": n_pos,
        "n_negative": n_neg,
        "fraction_significant": n_sig / n_sequences if n_sequences > 0 else 0,
    }
