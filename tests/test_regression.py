"""Golden-value regression test for the validated numerical path.

The Igor cross-validation lives in `test_cbas.py` and needs `igor_cbas/data`,
which is reference material that is not distributed. Those tests therefore skip
in CI, leaving the pipeline that reproduces Igor unguarded on every push.

This module closes that gap without shipping data. It builds a deterministic
synthetic cohort in code, runs the full comparative and correlative pipelines,
and compares every stage against values committed in
`tests/fixtures/regression_golden.json`.

That does not re-prove agreement with Igor. It proves the code path which was
validated against Igor is unchanged, which is the property CI can actually
enforce. Any edit to sequence enumeration, the studentized statistic, the
bootstrap, step-down, or k-FWER escalation moves at least one of these numbers
and fails the test.

Regenerate the fixture deliberately, never to make a red test green:

    pixi run python scripts/update_regression_golden.py
"""

import json
from pathlib import Path

import numpy as np
import pytest

from pycbas import (CBASParams, build_count_matrix, compute_test_stats,
                    bootstrap_test_stats, romano_wolf_stepdown,
                    run_cbas_comparative, run_cbas_correlative)

GOLDEN_PATH = Path(__file__).parent / "fixtures" / "regression_golden.json"

# Deterministic synthetic cohort. Two groups differing in how often they emit a
# planted motif, which guarantees a nonempty significant set and therefore
# exercises step-down and k-FWER escalation rather than just the early return.
N_PER_GROUP = 24
N_SESSIONS = 4
TRIALS_PER_SESSION = 60
NUM_ARMS = 3
MOTIF = (1, 2, 1)
SEED = 20260828


def make_cohort():
    """Build (subjects_data, group_labels, covariate) with no randomness escaping.

    Each subject is an array of (session, choice, reward, contingency) rows in the
    layout `load_subject_data` returns, so the pipeline sees exactly the shape it
    sees for real data. Group 1 inserts the motif more often, at a rate that also
    varies per subject so the studentized denominator is nonzero.
    """
    rng = np.random.default_rng(SEED)
    subjects, labels = [], []
    for group in (0, 1):
        for i in range(N_PER_GROUP):
            motif_p = 0.10 + 0.16 * group + 0.02 * (i % 4)
            rows = []
            for session in range(N_SESSIONS):
                t = 0
                while t < TRIALS_PER_SESSION:
                    if rng.random() < motif_p and t + len(MOTIF) <= TRIALS_PER_SESSION:
                        choices = MOTIF
                    else:
                        choices = (int(rng.integers(0, NUM_ARMS)),)
                    for c in choices:
                        reward = int(c == 2)
                        rows.append((session, c, reward, 2))
                        t += 1
            subjects.append(np.array(rows, dtype=np.int32))
            labels.append(group)
    labels = np.array(labels)
    # A covariate correlated with the motif rate, for the correlative path.
    covariate = np.array([0.10 + 0.16 * (l) + 0.02 * (i % 4)
                          for l, i in zip(labels, list(range(N_PER_GROUP)) * 2)])
    return subjects, labels, covariate


def encoded_motif():
    """MOTIF as it appears after reward encoding, so the lookup cannot silently miss.

    `build_count_matrix` with encode_reward=True stores choice + reward*num_arms,
    and reward here is deterministic from the choice, so the mapping is exact.
    """
    return tuple(c + int(c == 2) * NUM_ARMS for c in MOTIF)


def params():
    return CBASParams(num_arms=NUM_ARMS, seq_len_max=3, criterion=200,
                      resample_number=2000, alpha=0.5, gamma=0.05)


def compute_all():
    """Every quantity the golden file pins, as plain Python types."""
    subjects, labels, covariate = make_cohort()
    p = params()

    out = {}
    for block_aware in (False, True):
        seqs, counts = build_count_matrix(subjects, p, contingency=2,
                                          encode_reward=True, block_aware=block_aware)
        key = "block_aware" if block_aware else "pooled"
        out[key] = {
            "n_sequences": len(seqs),
            "count_total": int(counts.sum()),
            "count_checksum": float(np.round((counts * np.arange(1, counts.shape[1] + 1)).sum(), 6)),
            "row_sums_head": [int(v) for v in counts.sum(axis=1)[:6]],
            "motif_column_total": int(counts[:, seqs.index(encoded_motif())].sum()),
        }

    seqs, counts = build_count_matrix(subjects, p, contingency=2,
                                      encode_reward=True, block_aware=False)
    grp = [np.flatnonzero(labels == 0), np.flatnonzero(labels == 1)]

    stats = compute_test_stats(counts, grp)
    finite = np.isfinite(stats)
    out["test_stats"] = {
        "n_finite": int(finite.sum()),
        "max": float(np.round(np.nanmax(stats), 8)),
        "sum_finite": float(np.round(stats[finite].sum(), 6)),
    }

    null = bootstrap_test_stats(counts, grp, p, rng=np.random.default_rng(SEED))
    null_arr = null[0] if isinstance(null, tuple) else null
    out["bootstrap"] = {
        "shape": list(null_arr.shape),
        "n_finite": int(np.isfinite(null_arr).sum()),
        "mean_finite": float(np.round(np.nanmean(null_arr), 8)),
    }

    comp = run_cbas_comparative(subjects, labels, p, contingency=2,
                                encode_reward=True)
    out["comparative"] = {
        "n_sequences": len(comp.sequences),
        "n_significant": int(comp.n_significant),
        "k_final": int(comp.k_final),
        "g_value_sum": float(np.round(np.nansum(comp.g_values), 6)),
        "significant_sequences": sorted(
            list(s) for s, m in zip(comp.sequences, comp.significant_mask) if m
        )[:12],
    }

    corr = run_cbas_correlative(subjects, covariate, p, contingency=2,
                                encode_reward=True)
    out["correlative"] = {
        "n_sequences": len(corr.sequences),
        "n_significant": int(corr.n_significant),
        "k_final": int(corr.k_final),
        "g_value_sum": float(np.round(np.nansum(corr.g_values), 6)),
    }
    return out


@pytest.fixture(scope="module")
def actual():
    return compute_all()


@pytest.fixture(scope="module")
def golden():
    if not GOLDEN_PATH.exists():
        pytest.fail(f"missing golden fixture {GOLDEN_PATH}; "
                    "run scripts/update_regression_golden.py")
    return json.loads(GOLDEN_PATH.read_text())


def test_cohort_is_deterministic():
    """The synthetic cohort must not depend on global RNG state."""
    a = make_cohort()[0]
    b = make_cohort()[0]
    assert len(a) == len(b) == 2 * N_PER_GROUP
    for x, y in zip(a, b):
        assert np.array_equal(x, y)


@pytest.mark.parametrize("section", ["pooled", "block_aware"])
def test_count_matrix_unchanged(actual, golden, section):
    assert actual[section] == golden[section]


def test_test_stats_unchanged(actual, golden):
    assert actual["test_stats"] == golden["test_stats"]


def test_bootstrap_null_unchanged(actual, golden):
    assert actual["bootstrap"] == golden["bootstrap"]


def test_comparative_pipeline_unchanged(actual, golden):
    assert actual["comparative"] == golden["comparative"]


def test_correlative_pipeline_unchanged(actual, golden):
    assert actual["correlative"] == golden["correlative"]


def test_stepdown_monotonicity(actual):
    """g-values from step-down must be non-decreasing in rank, by construction."""
    subjects, labels, _ = make_cohort()
    p = params()
    _, counts = build_count_matrix(subjects, p, contingency=2, encode_reward=True)
    grp = [np.flatnonzero(labels == 0), np.flatnonzero(labels == 1)]
    stats = compute_test_stats(counts, grp)
    null = bootstrap_test_stats(counts, grp, p, rng=np.random.default_rng(SEED))
    null_arr = null[0] if isinstance(null, tuple) else null
    g = romano_wolf_stepdown(stats, null_arr, k=1)
    g = np.asarray(g[0] if isinstance(g, tuple) else g, dtype=float)
    order = np.argsort(np.where(np.isfinite(stats), -stats, np.inf))
    seq = g[order]
    seq = seq[np.isfinite(seq)]
    assert np.all(np.diff(seq) >= -1e-12), "step-down g-values must be monotone in rank"
