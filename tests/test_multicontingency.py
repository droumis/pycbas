"""Tests for counting several contingencies as separate hypothesis sets.

The published convention: each contingency is analysed separately, so the same arm
sequence under two contingencies is two hypotheses, and the multiplicity
correction runs over the concatenation. "If you have 100 sequences in the first
contingency and 200 in the second, you are evaluating 300 sequences total."

The synthetic tests run everywhere. The cohort tests use the unpublished lesion
data and skip without it, via the fixtures in conftest.py.
"""

import numpy as np
import pytest

from pycbas import CBASParams
from pycbas.contingency import (SubjectRecord, assign_contingency_blocks,
                                build_multicontingency_count_matrix,
                                shared_contingency_blocks,
                                load_cohort_with_contingencies,
                                load_subject_data_with_contingencies)
from pycbas.criterion import criterion_trial, as_enumeration_cutoff
from pycbas.io import enumerate_sequences_block_aware
from pycbas.pipeline import run_cbas_multicontingency


def make_record(sessions, choices, rewards, centres, lefts):
    """Build a SubjectRecord directly, bypassing the file format."""
    sessions = np.asarray(sessions, dtype=np.int64)
    block, blocks = assign_contingency_blocks(
        sessions, np.asarray(centres), np.asarray(lefts))
    return SubjectRecord(sessions, np.asarray(choices, dtype=np.int64),
                         np.asarray(rewards, dtype=np.int64), block, blocks)


def two_block_record(seed=0, n=40):
    """Explore block, then two alternation contingencies with distinct arms."""
    rng = np.random.default_rng(seed)
    sessions, choices, rewards, centres, lefts = [], [], [], [], []
    spec = [(-1, -1, 0), (2, 1, 1), (4, 3, 2)]   # exploration, arms (2,1), arms (4,3)
    for centre, left, session in spec:
        for _ in range(n):
            sessions.append(session)
            choices.append(int(rng.integers(0, 3)))
            rewards.append(int(rng.integers(0, 2)))
            centres.append(centre)
            lefts.append(left)
    return make_record(sessions, choices, rewards, centres, lefts)


# ---------------------------------------------------------------------------
# Block validation
# ---------------------------------------------------------------------------

class TestSharedBlocks:
    def test_aligned_cohort(self):
        records = [two_block_record(s) for s in range(3)]
        assert shared_contingency_blocks(records) == [1, 2]

    def test_exploration_is_excluded(self):
        assert 0 not in shared_contingency_blocks([two_block_record()])

    def test_misaligned_arms_raise(self):
        """Block i must denote the same arms for every subject."""
        a = two_block_record()
        b = make_record([0, 0, 1, 1], [0, 0, 0, 0], [1, 1, 1, 1],
                        [-1, -1, 3, 3], [-1, -1, 2, 2])   # block 1 is (3,2) here
        with pytest.raises(ValueError, match="not aligned across subjects"):
            shared_contingency_blocks([a, b])

    def test_missing_contingency_raises(self):
        """Zero would assert 'never produced it' rather than 'not measured'."""
        a = two_block_record()
        b = make_record([0, 0], [0, 0], [1, 1], [-1, -1], [-1, -1])
        with pytest.raises(ValueError, match="not every subject ran every contingency"):
            shared_contingency_blocks([a, b])

    def test_empty_input_raises(self):
        with pytest.raises(ValueError):
            shared_contingency_blocks([])


# ---------------------------------------------------------------------------
# Count matrix
# ---------------------------------------------------------------------------

class TestMultiContingencyCounts:
    @pytest.fixture
    def cohort(self):
        return [two_block_record(s) for s in range(5)]

    def test_columns_are_grouped_by_contingency(self, cohort):
        params = CBASParams(num_arms=3, seq_len_max=2, criterion=10_000)
        sequences, _ = build_multicontingency_count_matrix(cohort, params)
        blocks = [b for b, _ in sequences]
        assert blocks == sorted(blocks)

    def test_shape_and_label_format(self, cohort):
        params = CBASParams(num_arms=3, seq_len_max=2, criterion=10_000)
        sequences, counts = build_multicontingency_count_matrix(cohort, params)
        assert counts.shape == (len(cohort), len(sequences))
        block, seq = sequences[0]
        assert isinstance(block, int) and isinstance(seq, tuple)

    def test_per_contingency_columns_sum_to_the_total(self, cohort):
        """The convention's 100 + 200 = 300."""
        params = CBASParams(num_arms=3, seq_len_max=2, criterion=10_000)
        sequences, counts = build_multicontingency_count_matrix(cohort, params)
        from collections import Counter
        per_block = Counter(b for b, _ in sequences)
        assert sum(per_block.values()) == counts.shape[1]
        assert len(per_block) == 2

    def test_same_sequence_under_two_contingencies_is_two_columns(self, cohort):
        params = CBASParams(num_arms=3, seq_len_max=2, criterion=10_000)
        sequences, _ = build_multicontingency_count_matrix(cohort, params)
        from collections import Counter
        repeated = [s for s, n in Counter(s for _, s in sequences).items() if n > 1]
        assert repeated, "expected at least one sequence shared across contingencies"
        for seq in repeated:
            blocks = [b for b, s in sequences if s == seq]
            assert len(blocks) == len(set(blocks)), "one column per contingency"

    def test_matches_an_independent_per_block_build(self, cohort):
        """The whole point: each contingency counted as if it were alone."""
        params = CBASParams(num_arms=3, seq_len_max=2, criterion=10_000)
        sequences, counts = build_multicontingency_count_matrix(cohort, params)
        index = {key: i for i, key in enumerate(sequences)}

        for block in shared_contingency_blocks(cohort):
            for row, record in enumerate(cohort):
                streams = record.symbol_blocks_for(block, params.num_arms, True)
                cutoff = as_enumeration_cutoff(
                    criterion_trial(record.reward_blocks_for(block), 0,
                                    params.criterion),
                    sum(len(s) for s in streams))
                for seq_len in (1, 2):
                    expected = enumerate_sequences_block_aware(
                        streams, seq_len, cutoff)
                    for seq, n in expected.items():
                        assert counts[row, index[(block, seq)]] == n

    def test_blocks_argument_subsets(self, cohort):
        params = CBASParams(num_arms=3, seq_len_max=2, criterion=10_000)
        all_seqs, _ = build_multicontingency_count_matrix(cohort, params)
        one_seqs, one = build_multicontingency_count_matrix(cohort, params, blocks=[1])
        assert {b for b, _ in one_seqs} == {1}
        assert one.shape[1] == sum(1 for b, _ in all_seqs if b == 1)

    def test_unknown_block_raises(self, cohort):
        params = CBASParams(num_arms=3, seq_len_max=2, criterion=10_000)
        with pytest.raises(ValueError, match="not shared by all subjects"):
            build_multicontingency_count_matrix(cohort, params, blocks=[1, 99])

    def test_criterion_is_applied_per_contingency(self):
        """Each contingency is its own learning episode with its own cutoff.

        Block 1 is rewarded throughout and reaches the criterion early. Block 2 is
        unrewarded until its final trials, so a criterion evaluated over the whole
        stream rather than within the contingency would cut it in the wrong place.
        """
        n = 12
        sessions = [0] * n + [1] * n
        choices = [0] * (2 * n)
        rewards = [1] * n + [0] * (n - 3) + [1, 1, 1]
        centres = [2] * n + [4] * n
        lefts = [1] * n + [3] * n
        record = make_record(sessions, choices, rewards, centres, lefts)

        first = criterion_trial(record.reward_blocks_for(0), 1, 2)
        second = criterion_trial(record.reward_blocks_for(1), 1, 2)
        assert first == 1, "second reward of block 1 is at its own trial index 1"
        assert second == n - 2, "block 2's index restarts within the contingency"

    def test_count_matrix_uses_a_per_contingency_cutoff(self):
        """The cutoff must be local to each contingency, with exact counts.

        Block 1 is rewarded throughout, so an order-1 criterion of 2 lands at its
        trial index 1 and admits two windows. Block 2 is rewarded only at the end,
        so its own criterion lands at index 3 and admits four windows.

        A criterion evaluated over the concatenated stream instead would place the
        single cutoff at global index 1, leaving block 2 with no counts at all. So
        block 2 being populated is what distinguishes the two implementations, and
        asserting exact counts pins it.
        """
        record = make_record(
            sessions=[1, 1, 1, 1, 2, 2, 2, 2],
            choices=[0, 0, 0, 0, 1, 1, 1, 1],
            rewards=[1, 1, 1, 1, 0, 0, 1, 1],
            centres=[2, 2, 2, 2, 4, 4, 4, 4],
            lefts=[1, 1, 1, 1, 3, 3, 3, 3],
        )
        params = CBASParams(num_arms=3, seq_len_max=1, criterion=2,
                            criterion_order=1)
        sequences, counts = build_multicontingency_count_matrix([record], params)
        got = {key: int(counts[0, i]) for i, key in enumerate(sequences)}

        # symbol = choice + reward * num_arms
        assert got == {
            (0, (3,)): 2,        # block 1: rewarded choice 0, cutoff at index 1
            (1, (1,)): 2,        # block 2: unrewarded choice 1
            (1, (4,)): 2,        # block 2: rewarded choice 1, cutoff at index 3
        }, "cutoffs must be local to each contingency"


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

class TestPipeline:
    def test_runs_end_to_end(self):
        cohort = [two_block_record(s, n=60) for s in range(12)]
        labels = np.array([0, 1] * 6)
        params = CBASParams(num_arms=3, seq_len_max=2, criterion=10_000,
                            resample_number=200)
        result = run_cbas_multicontingency(cohort, labels, params)
        assert len(result.sequences) == len(result.significant_mask)
        assert result.k_final >= 1
        assert all(isinstance(b, int) and isinstance(s, tuple)
                   for b, s in result.sequences)

    def test_chunked_and_unchunked_agree(self):
        """The chunked path is the default here, so it needs pinning.

        Counting several contingencies multiplies the hypothesis space, which is
        exactly when the memory-efficient path matters, so a divergence between
        the two would show up first on multi-contingency data and nowhere else.
        Both default to the same seed, so agreement should be exact rather than
        approximate.
        """
        cohort = [two_block_record(s, n=60) for s in range(12)]
        labels = np.array([0, 1] * 6)
        params = CBASParams(num_arms=3, seq_len_max=2, criterion=10_000,
                            resample_number=300)
        chunked = run_cbas_multicontingency(cohort, labels, params, chunked=True)
        direct = run_cbas_multicontingency(cohort, labels, params, chunked=False)

        assert chunked.sequences == direct.sequences
        assert chunked.k_final == direct.k_final
        assert np.array_equal(chunked.significant_mask, direct.significant_mask)
        np.testing.assert_array_equal(np.isnan(chunked.g_values),
                                      np.isnan(direct.g_values))
        finite = ~np.isnan(chunked.g_values)
        np.testing.assert_allclose(chunked.g_values[finite],
                                   direct.g_values[finite])

    def test_single_block_subset_matches_hypothesis_count(self):
        cohort = [two_block_record(s, n=60) for s in range(12)]
        labels = np.array([0, 1] * 6)
        params = CBASParams(num_arms=3, seq_len_max=2, criterion=10_000,
                            resample_number=200)
        both = run_cbas_multicontingency(cohort, labels, params)
        one = run_cbas_multicontingency(cohort, labels, params, blocks=[1])
        assert len(one.sequences) < len(both.sequences)


# ---------------------------------------------------------------------------
# A multi-contingency cohort, supplied separately
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def cohort(lesion_cohort_dir):
    return load_cohort_with_contingencies(lesion_cohort_dir)


class TestLesionCohort:
    def test_loads_in_numeric_order(self, cohort):
        records, info = cohort
        # Consistency between the two loaders rather than a pinned cohort size: the
        # claim is that every info row produced a record, which is what makes the
        # positional correspondence between the two lists safe to rely on.
        assert len(records) == len(info) > 0

    def test_every_subject_runs_all_six_contingencies(self, cohort):
        records, _ = cohort
        assert shared_contingency_blocks(records) == [1, 2, 3, 4, 5, 6]

    def test_block_one_and_five_share_arms_but_stay_separate(self, cohort):
        """The recurrence the reference implementation does not merge."""
        records, _ = cohort
        blocks = {b.block: (b.centre, b.left_outer)
                  for b in records[0].alternation_blocks()}
        assert blocks[1] == blocks[5], "expected the arms-2-3-4 recurrence"
        params = CBASParams(num_arms=6, seq_len_max=2, criterion=100,
                            criterion_order=4)
        sequences, _ = build_multicontingency_count_matrix(
            records[:20], params, blocks=[1, 5])
        present = {b for b, _ in sequences}
        assert present == {1, 5}, "the recurrence must contribute its own columns"

    def test_counts_are_nonzero_and_shaped(self, cohort):
        records, _ = cohort
        params = CBASParams(num_arms=6, seq_len_max=2, criterion=100,
                            criterion_order=4)
        sequences, counts = build_multicontingency_count_matrix(
            records[:30], params)
        assert counts.shape == (30, len(sequences))
        assert counts.sum() > 0
        assert (counts.sum(axis=1) > 0).all(), "every subject contributes counts"


class TestHeaderDetection:
    """The export has shipped both with and without a header line.

    Unconditionally skipping the first line drops a real trial when no header is
    present, which shifts every later trial index and silently changes the
    criterion. Detection has to be based on the content of the line.
    """

    # Session 0 is exploration throughout and session 1 alternation throughout.
    # Mixing the two inside one session is a mid-session contingency change, which
    # assign_contingency_blocks rightly refuses.
    ROWS = ["0,2,1,,,0,72,1,1",
            "0,5,1,,,246,450,1,1",
            "1,4,0,3,2,609,735,1,1",
            "1,3,1,3,2,900,980,1,1"]

    def _load(self, tmp_path, text):
        fp = tmp_path / "an0.txt"
        fp.write_text(text)
        return load_subject_data_with_contingencies(fp)

    def test_headerless_file_keeps_every_trial(self, tmp_path):
        rec = self._load(tmp_path, "\n".join(self.ROWS) + "\n")
        assert len(rec) == 4, "first data row must not be eaten as a header"
        assert rec.choice[0] == 2
        assert rec.reward[0] == 1

    def test_header_is_still_skipped(self, tmp_path):
        header = "session,choice,reward,centre,left,t1,t2,f1,f2"
        rec = self._load(tmp_path, header + "\n" + "\n".join(self.ROWS) + "\n")
        assert len(rec) == 4
        assert rec.choice[0] == 2

    def test_both_forms_agree(self, tmp_path):
        bare = self._load(tmp_path, "\n".join(self.ROWS) + "\n")
        with_header = self._load(
            tmp_path, "anInfo\n" + "\n".join(self.ROWS) + "\n")
        assert np.array_equal(bare.choice, with_header.choice)
        assert np.array_equal(bare.reward, with_header.reward)
        assert np.array_equal(bare.session, with_header.session)

    def test_blank_choice_or_reward_still_dropped(self, tmp_path):
        rows = self.ROWS + ["1,,1,3,2,1100,1200,1,1",
                            "1,3,,3,2,1300,1400,1,1"]
        rec = self._load(tmp_path, "\n".join(rows) + "\n")
        assert len(rec) == 4, "blank choice/reward rows are dropped, matching Igor"
