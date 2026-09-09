"""Tests for higher-order criteria.

Two layers. The synthetic tests run everywhere and pin the semantics that are
easy to get wrong. The reference test reproduces the Igor reference output exactly,
and skips when the multi-contingency cohort is absent, following the same
pattern as the Igor cross-validation in test_cbas.py.
"""

from pathlib import Path

import numpy as np
import pytest

from pycbas import CBASParams, build_count_matrix, extract_choice_stream
from pycbas.io import enumerate_sequences
from pycbas.core import reward_blocks, subject_criteria
from pycbas.criterion import (perfect_run_starts, criterion_trial,
                              reached_criterion, criterion_trials_by_subject)
from pycbas.contingency import (assign_contingency_blocks,
                                load_subject_data_with_contingencies,
                                load_cohort_info)

from conftest import LESION_COHORT_DIR, CRITERION_REFERENCE  # noqa: E402


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

    def test_mid_session_change_raises(self):
        """The one case where this segmentation and Igor's wConting disagree.

        Igor gives the whole session to the new contingency; this splits it at the
        exact trial. Failing loudly is the point, because the divergence would
        otherwise be invisible in the output.
        """
        session = np.array([0, 0, 1, 1, 1])
        centre = np.array([2, 2, 2, 3, 3])
        left = np.array([1, 1, 1, 2, 2])
        with pytest.raises(ValueError, match="partway through session"):
            assign_contingency_blocks(session, centre, left)

    def test_mid_session_change_can_be_opted_into(self):
        session = np.array([0, 0, 1, 1, 1])
        centre = np.array([2, 2, 2, 3, 3])
        left = np.array([1, 1, 1, 2, 2])
        block, _ = assign_contingency_blocks(session, centre, left,
                                             allow_mid_session_change=True)
        assert list(block) == [0, 0, 0, 1, 1], "session 1 is split at the change"

    def test_boundary_change_does_not_raise(self):
        session = np.array([0, 0, 1, 1])
        centre = np.array([2, 2, 3, 3])
        left = np.array([1, 1, 2, 2])
        block, _ = assign_contingency_blocks(session, centre, left)
        assert list(block) == [0, 0, 1, 1]


# ---------------------------------------------------------------------------
# Pipeline wiring
# ---------------------------------------------------------------------------

def _subject(choices, rewards, sessions=None, contingency=2):
    """Build a subject array in the published 4-column layout."""
    n = len(choices)
    if sessions is None:
        sessions = np.zeros(n, dtype=int)
    return np.column_stack([sessions, choices, rewards,
                            np.full(n, contingency)]).astype(np.int32)


class TestPipelineWiring:
    def test_order_zero_matches_the_historical_path(self):
        """The default must reproduce a direct enumerate_sequences call."""
        rng = np.random.default_rng(0)
        choices = rng.integers(0, 3, 200)
        rewards = (choices == 2).astype(int)
        subj = _subject(choices, rewards)
        params = CBASParams(num_arms=3, seq_len_max=3, criterion=50)

        seqs, counts = build_count_matrix([subj], params, contingency=2)
        stream = extract_choice_stream(subj, 2, 3, encode_reward=True)
        expected = enumerate_sequences(stream, 2, 50)
        idx = seqs.index(next(iter(expected)))
        assert counts[0, idx] == expected[next(iter(expected))]

    def test_order_one_stops_at_the_nth_reward(self):
        rewards = np.array([0, 1, 0, 1, 0, 1, 0, 1])
        choices = np.arange(8) % 3
        subj = _subject(choices, rewards)
        params = CBASParams(num_arms=3, seq_len_max=1, criterion=3,
                            criterion_order=1)
        crit = subject_criteria([subj], params, contingency=2)
        assert crit[0] == 5, "third reward is at trial index 5"

    def test_order_one_yields_fewer_counts_than_order_zero(self):
        rng = np.random.default_rng(1)
        choices = rng.integers(0, 3, 300)
        rewards = (choices == 2).astype(int)
        subj = _subject(choices, rewards)
        low = CBASParams(num_arms=3, seq_len_max=2, criterion=10,
                         criterion_order=1)
        high = CBASParams(num_arms=3, seq_len_max=2, criterion=200,
                          criterion_order=0)
        _, few = build_count_matrix([subj], low, contingency=2)
        _, many = build_count_matrix([subj], high, contingency=2)
        assert few.sum() < many.sum()

    def test_subject_short_of_criterion_reports_inf(self):
        subj = _subject(np.zeros(20, dtype=int), np.zeros(20, dtype=int))
        params = CBASParams(num_arms=3, seq_len_max=1, criterion=5,
                            criterion_order=1)
        crit = subject_criteria([subj], params, contingency=2)
        assert np.isinf(crit[0])

    def test_inf_criterion_counts_every_window(self):
        """A subject that never reaches criterion is not truncated."""
        choices = np.zeros(30, dtype=int)
        subj = _subject(choices, np.zeros(30, dtype=int))
        unreachable = CBASParams(num_arms=3, seq_len_max=1, criterion=99,
                                 criterion_order=1)
        no_cutoff = CBASParams(num_arms=3, seq_len_max=1, criterion=10_000)
        _, a = build_count_matrix([subj], unreachable, contingency=2)
        _, b = build_count_matrix([subj], no_cutoff, contingency=2)
        assert a.sum() == b.sum() == 30

    def test_reward_blocks_follow_the_enumeration_structure(self):
        sessions = np.array([0, 0, 0, 1, 1])
        subj = _subject(np.arange(5) % 3, np.ones(5, dtype=int), sessions)
        pooled = reward_blocks(subj, 2, block_aware=False)
        split = reward_blocks(subj, 2, block_aware=True)
        assert len(pooled) == 1 and len(pooled[0]) == 5
        assert [len(b) for b in split] == [3, 2]

    def test_block_aware_criterion_respects_sessions(self):
        """A run cannot span sessions when the enumeration cannot either."""
        sessions = np.array([0, 0, 1, 1])
        subj = _subject(np.zeros(4, dtype=int), np.ones(4, dtype=int), sessions)
        params = CBASParams(num_arms=3, seq_len_max=1, criterion=1,
                            criterion_order=3)
        assert np.isinf(subject_criteria([subj], params, contingency=2,
                                         block_aware=True)[0])
        assert subject_criteria([subj], params, contingency=2,
                                block_aware=False)[0] == 0


# ---------------------------------------------------------------------------
# Reference: exact reproduction of the Igor reference output
# ---------------------------------------------------------------------------

def _load_reference(path):
    """Parse allTrialToPerfect: one row per subject, one column per contingency.

    Column 0 is the exploration contingency and is always blank. A blank in any
    other column means the subject did not reach the criterion.
    """
    text = path.read_bytes().decode("utf-8").replace("\r", "\n")
    rows = [l for l in text.split("\n") if l.strip(",").strip()]
    reference = {}
    for subject, line in enumerate(rows):
        for block, value in enumerate(line.split(",")[1:], start=1):
            value = value.strip()
            reference[(subject, block)] = (
                np.inf if value == "" or value.lower() in ("inf", "nan")
                else float(value))
    return reference


# The reference setting: 4th order, 100 runs of four rewarded choices.
REFERENCE_ORDER = 4
REFERENCE_COUNT = 100


@pytest.fixture(scope="module")
def computed(lesion_cohort_dir):
    files = [f for f in lesion_cohort_dir.glob("an*.txt") if f.stem != "anInfo"]
    files.sort(key=lambda p: int(p.stem[2:]))
    out = {}
    for path in files:
        subject = int(path.stem[2:])
        record = load_subject_data_with_contingencies(path)
        for block in record.alternation_blocks():
            out[(subject, block.block)] = criterion_trial(
                record.reward_blocks_for(block.block), REFERENCE_ORDER, REFERENCE_COUNT)
    return out


class TestAgainstIgorReference:
    """The reference setting: 4th order, 100 runs of four rewarded choices."""

    def test_cohort_shape(self, computed):
        subjects = {s for s, _ in computed}
        assert len(subjects) == 222
        per_subject = {s: sum(1 for a, _ in computed if a == s) for s in subjects}
        assert set(per_subject.values()) == {6}, "every subject runs six contingencies"

    def test_every_value_matches_igor(self, computed, criterion_reference_file):
        reference = _load_reference(criterion_reference_file)
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

    def test_non_reachers_agree(self, computed, criterion_reference_file):
        """Agreement on who failed matters as much as agreement on the values."""
        reference = _load_reference(criterion_reference_file)
        shared = sorted(set(reference) & set(computed))
        igor_missed = {k for k in shared if np.isinf(reference[k])}
        ours_missed = {k for k in shared if np.isinf(computed[k])}
        assert igor_missed == ours_missed
        assert len(igor_missed) == 151

    def test_cohort_info_parses(self, lesion_cohort_dir):
        info = load_cohort_info(lesion_cohort_dir / "anInfo.txt")
        assert len(info) == 222
        lesion = [r["lesion"] for r in info]
        assert lesion.count("Control") == 111
        assert lesion.count("Hippocampal Lesion") == 98
        assert lesion.count(None) == 13, "blank means surgery with no CT evidence"
