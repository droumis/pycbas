"""Tests for higher-order criteria.

Two layers. The synthetic tests run everywhere and pin the semantics that are
easy to get wrong. The reference test reproduces Kastner's Igor output exactly,
and skips when the unpublished lesion cohort is absent, following the same
pattern as the Igor cross-validation in test_cbas.py.
"""

from pathlib import Path

import numpy as np
import pytest

from pycbas.criterion import (perfect_run_starts, criterion_trial,
                              reached_criterion, criterion_trials_by_subject)
from pycbas.contingency import (assign_contingency_blocks,
                                load_subject_data_with_contingencies,
                                load_cohort_info)

DATA_DIR = Path(__file__).parent.parent / "data" / "rats_AllHipLesionData"
REFERENCE = Path(__file__).parent.parent / "igor_cbas" / "allTrialToPerfect.txt"


# ---------------------------------------------------------------------------
# Run counting
# ---------------------------------------------------------------------------

class TestPerfectRunStarts:
    def test_runs_overlap(self):
        """Seven rewarded trials contain four length-4 runs, not one."""
        starts = perfect_run_starts([np.ones(7, dtype=int)], order=4)
        assert list(starts) == [0, 1, 2, 3]

    def test_order_one_is_every_reward(self):
        rewards = np.array([1, 0, 1, 1, 0])
        starts = perfect_run_starts([rewards], order=1)
        assert list(starts) == [0, 2, 3]

    def test_runs_do_not_span_sessions(self):
        """Two sessions of two rewards each contain no length-3 run."""
        blocks = [np.ones(2, dtype=int), np.ones(2, dtype=int)]
        assert len(perfect_run_starts(blocks, order=3)) == 0

    def test_index_counts_trials_not_windows(self):
        """The position after a session advances by its trial count.

        Session one has 5 trials and can start two length-4 runs. The next
        session's first trial is therefore index 5, not index 2. Getting this
        wrong shifts every value by a few percent.
        """
        blocks = [np.ones(5, dtype=int), np.ones(4, dtype=int)]
        starts = perfect_run_starts(blocks, order=4)
        assert list(starts) == [0, 1, 5]

    def test_unrewarded_trial_breaks_the_run(self):
        rewards = np.array([1, 1, 0, 1, 1, 1, 1])
        starts = perfect_run_starts([rewards], order=3)
        assert list(starts) == [3, 4]

    def test_order_zero_rejected(self):
        with pytest.raises(ValueError):
            perfect_run_starts([np.ones(4, dtype=int)], order=0)


# ---------------------------------------------------------------------------
# Criterion
# ---------------------------------------------------------------------------

class TestCriterionTrial:
    def test_order_zero_returns_the_requested_index(self):
        blocks = [np.ones(500, dtype=int)]
        assert criterion_trial(blocks, order=0, count=200) == 200

    def test_order_zero_short_subject_has_no_cutoff(self):
        """Matches getZerothOrder, and the historical pycbas behaviour where
        min(criterion, len - seq_len) already admitted every window."""
        blocks = [np.ones(100, dtype=int)]
        assert criterion_trial(blocks, order=0, count=800) == np.inf

    def test_nth_run_index(self):
        starts = perfect_run_starts([np.ones(10, dtype=int)], order=4)
        assert criterion_trial([np.ones(10, dtype=int)], order=4, count=3) == starts[2]

    def test_not_reaching_count_gives_inf(self):
        blocks = [np.array([1, 1, 0, 1, 1])]
        assert criterion_trial(blocks, order=3, count=5) == np.inf

    def test_reached_criterion_helper(self):
        assert reached_criterion(412) is True
        assert reached_criterion(np.inf) is False

    def test_count_must_be_positive(self):
        with pytest.raises(ValueError):
            criterion_trial([np.ones(4, dtype=int)], order=1, count=0)

    def test_by_subject_preserves_inf(self):
        subjects = [[np.ones(10, dtype=int)], [np.zeros(10, dtype=int)]]
        out = criterion_trials_by_subject(subjects, order=2, count=1)
        assert out[0] == 0 and np.isinf(out[1])
        assert out.dtype == np.float64


# ---------------------------------------------------------------------------
# Contingency segmentation
# ---------------------------------------------------------------------------

class TestContingencyBlocks:
    def test_repeated_arm_pair_gets_a_new_block(self):
        """The central subtlety: a recurrence is a separate contingency.

        Keying on the arm pair rather than on block order would merge these two
        into one, which is what Igor's wConting avoids.
        """
        session = np.array([0, 0, 1, 1, 2, 2])
        centre = np.array([2, 2, 3, 3, 2, 2])
        left = np.array([1, 1, 2, 2, 1, 1])
        block, blocks = assign_contingency_blocks(session, centre, left)
        assert list(block) == [0, 0, 1, 1, 2, 2]
        assert len(blocks) == 3
        assert (blocks[0].centre, blocks[0].left_outer) == (2, 1)
        assert (blocks[2].centre, blocks[2].left_outer) == (2, 1)

    def test_exploration_is_flagged(self):
        block, blocks = assign_contingency_blocks(
            np.array([0, 1]), np.array([-1, 2]), np.array([-1, 1]))
        assert blocks[0].is_exploration
        assert not blocks[1].is_exploration
        assert blocks[0].right_outer is None

    def test_right_outer_is_implied(self):
        _, blocks = assign_contingency_blocks(
            np.array([0]), np.array([3]), np.array([1]))
        assert blocks[0].right_outer == 5


# ---------------------------------------------------------------------------
# Reference: exact reproduction of Kastner's Igor output
# ---------------------------------------------------------------------------

def _load_reference():
    """Parse allTrialToPerfect: one row per subject, one column per contingency.

    Column 0 is the exploration contingency and is always blank. A blank in any
    other column means the subject did not reach the criterion.
    """
    text = REFERENCE.read_bytes().decode("utf-8").replace("\r", "\n")
    rows = [l for l in text.split("\n") if l.strip(",").strip()]
    reference = {}
    for subject, line in enumerate(rows):
        for block, value in enumerate(line.split(",")[1:], start=1):
            value = value.strip()
            reference[(subject, block)] = (
                np.inf if value == "" or value.lower() in ("inf", "nan")
                else float(value))
    return reference


# Kastner's stated setting: 4th order, 100 runs of four rewarded choices.
REFERENCE_ORDER = 4
REFERENCE_COUNT = 100


@pytest.fixture(scope="module")
def computed():
    files = [f for f in DATA_DIR.glob("an*.txt") if f.stem != "anInfo"]
    files.sort(key=lambda p: int(p.stem[2:]))
    out = {}
    for path in files:
        subject = int(path.stem[2:])
        record = load_subject_data_with_contingencies(path)
        for block in record.alternation_blocks():
            out[(subject, block.block)] = criterion_trial(
                record.reward_blocks_for(block.block), REFERENCE_ORDER, REFERENCE_COUNT)
    return out


@pytest.mark.skipif(not DATA_DIR.is_dir() or not REFERENCE.exists(),
                    reason="unpublished lesion cohort or Igor reference not present")
class TestAgainstIgorReference:
    """Kastner's stated setting: 4th order, 100 runs of four rewarded choices."""

    def test_cohort_shape(self, computed):
        subjects = {s for s, _ in computed}
        assert len(subjects) == 222
        per_subject = {s: sum(1 for a, _ in computed if a == s) for s in subjects}
        assert set(per_subject.values()) == {6}, "every subject runs six contingencies"

    def test_every_value_matches_igor(self, computed):
        reference = _load_reference()
        shared = sorted(set(reference) & set(computed))
        assert len(shared) == 1332

        mismatches = [
            (k, reference[k], computed[k]) for k in shared
            if not ((np.isinf(reference[k]) and np.isinf(computed[k]))
                    or reference[k] == computed[k])
        ]
        assert not mismatches, (
            f"{len(mismatches)} of {len(shared)} disagree; first few: "
            + "; ".join(f"subject {a} block {b}: igor={r} ours={c}"
                        for (a, b), r, c in mismatches[:5]))

    def test_non_reachers_agree(self, computed):
        """Agreement on who failed matters as much as agreement on the values."""
        reference = _load_reference()
        shared = sorted(set(reference) & set(computed))
        igor_missed = {k for k in shared if np.isinf(reference[k])}
        ours_missed = {k for k in shared if np.isinf(computed[k])}
        assert igor_missed == ours_missed
        assert len(igor_missed) == 151

    def test_cohort_info_parses(self):
        info = load_cohort_info(DATA_DIR / "anInfo.txt")
        assert len(info) == 222
        lesion = [r["lesion"] for r in info]
        assert lesion.count("Control") == 111
        assert lesion.count("Hippocampal Lesion") == 98
        assert lesion.count(None) == 13, "blank means surgery with no CT evidence"
