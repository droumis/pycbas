"""Tests for `sequence_occurrences`, all on synthetic cohorts.

The central property is that the occurrence table is the count matrix unrolled. For
every subject and sequence the number of rows equals the count, and the positions
are exactly the windows the enumeration counted. A reference written here as plain
loops checks the positions, and the function's own check against the count matrix
checks the totals.
"""

from dataclasses import replace

import numpy as np
import pytest

from pycbas import (CBASParams, Cohort, Subject, build_count_matrix,
                    run_cbas_comparative, run_cbas_multicontingency,
                    sequence_occurrences)
from pycbas.contingency import build_multicontingency_count_matrix
from pycbas.core import reward_blocks
from pycbas.criterion import as_enumeration_cutoff, criterion_trial
from pycbas.io import extract_choice_stream, extract_choice_streams_by_block

from test_multicontingency import two_block_record


def subject(choices, rewards=None, sessions=None, conditions=None, id="s"):
    n = len(choices)
    return Subject(id=id, trials=np.column_stack([
        np.ones(n, dtype=int) if sessions is None else sessions,
        choices,
        np.zeros(n, dtype=int) if rewards is None else rewards,
        np.full(n, 2) if conditions is None else conditions,
    ]))


def irregular_cohort(seed=0, num_arms=3):
    """Subjects with the session layouts most likely to break a position map.

    Sessions of length one, session labels out of order and recurring after a gap,
    trials of another condition interleaved, and one subject too short for any
    sequence longer than two.
    """
    rng = np.random.default_rng(seed)
    layouts = [
        [5, 3, 1, 4] * 3 + [2] * 10,
        [1, 1, 1, 2, 2, 3, 3, 3, 3] * 6,
        [7] * 30 + [2] * 25 + [7] * 5,
        list(np.repeat([1, 2, 3, 4], 18)),
        [1, 1],
    ]
    subjects = []
    for i, sessions in enumerate(layouts):
        n = len(sessions)
        conditions = np.where(rng.random(n) < 0.15, 0, 2)
        subjects.append(subject(rng.integers(0, num_arms, n),
                                rng.integers(0, 2, n), np.asarray(sessions),
                                conditions, id=f"s{i}"))
    return Cohort(subjects)


def reference_starts(segments, length, cutoff):
    """Counted start positions per sequence, by the enumeration rules, as loops."""
    out, offset = {}, 0
    for seg in segments:
        for i in range(len(seg) - length + 1):
            if offset + i <= cutoff:
                out.setdefault(tuple(seg[i:i + length].tolist()), []).append(offset + i)
        offset += len(seg)
    return out


def as_table(occ):
    """{(subject_id, sequence): starts} from an occurrence table."""
    table = {}
    for sid, seq, start in zip(occ["subject_id"], occ["sequence"], occ["start"]):
        table.setdefault((sid, seq), []).append(int(start))
    return table


class TestKnownPositions:
    params = CBASParams(num_arms=3, seq_len_max=3, criterion=100)

    def cohort(self):
        return Cohort([subject([0, 0, 0, 0, 1, 0, 0, 1])])

    def test_overlapping_windows_each_count(self):
        occ, _ = sequence_occurrences(self.cohort(), self.params,
                                      [(0, 0, 0), (0, 0), (0, 1)],
                                      encode_reward=False)
        table = as_table(occ)
        assert table[("s", (0, 0, 0))] == [0, 1]
        assert table[("s", (0, 0))] == [0, 1, 2, 5]
        assert table[("s", (0, 1))] == [3, 6]

    def test_the_cutoff_bounds_the_start_inclusively(self):
        params = replace(self.params, criterion=5)
        occ, subjects = sequence_occurrences(self.cohort(), params, [(0, 0), (0, 1)],
                                             encode_reward=False)
        table = as_table(occ)
        assert table[("s", (0, 0))] == [0, 1, 2, 5]
        assert table[("s", (0, 1))] == [3]
        assert subjects["cutoff"].tolist() == [5]
        assert subjects["reached"].tolist() == [True]

    def test_an_unreached_criterion_counts_every_window(self):
        _, subjects = sequence_occurrences(self.cohort(), self.params, [(0,)],
                                           encode_reward=False)
        assert subjects["reached"].tolist() == [False]
        assert subjects["cutoff"].tolist() == [8]
        assert subjects["n_trials"].tolist() == [8]

    def test_a_cutoff_shorter_than_the_sequence(self):
        params = replace(self.params, criterion=1)
        occ, _ = sequence_occurrences(self.cohort(), params, [(0, 0, 0)],
                                      encode_reward=False)
        assert as_table(occ)[("s", (0, 0, 0))] == [0, 1]


@pytest.mark.parametrize("order,count", [(0, 40), (1, 8), (2, 3)])
@pytest.mark.parametrize("block_aware", [False, True])
@pytest.mark.parametrize("encode_reward", [False, True])
@pytest.mark.parametrize("contingency", [2, None])
def test_positions_are_exactly_the_counted_windows(order, count, block_aware,
                                                   encode_reward, contingency):
    cohort = irregular_cohort()
    params = CBASParams(num_arms=3, seq_len_max=3, criterion=count,
                        criterion_order=order)
    settings = dict(contingency=contingency, block_aware=block_aware,
                    encode_reward=encode_reward)
    matrix = build_count_matrix(cohort, params, **settings)
    occ, subjects = sequence_occurrences(cohort, params, matrix.sequences,
                                         counted=matrix, **settings)
    table = as_table(occ)

    for i, subj in enumerate(cohort):
        if block_aware:
            segments = extract_choice_streams_by_block(subj.trials, contingency, 3,
                                                       encode_reward)
        else:
            segments = [extract_choice_stream(subj.trials, contingency, 3,
                                              encode_reward)]
        n = sum(len(s) for s in segments)
        cutoff = as_enumeration_cutoff(criterion_trial(
            reward_blocks(subj.trials, contingency, block_aware), order, count), n)
        expected = {}
        for length in range(1, 4):
            expected.update(reference_starts(segments, length, cutoff))
        got = {seq: starts for (sid, seq), starts in table.items() if sid == subj.id}
        assert got == expected, subj.id
        assert subjects["cutoff"][i] == cutoff
        for j, seq in enumerate(matrix.sequences):
            assert len(got.get(seq, [])) == matrix.counts[i, j]

    # `row` points at the first trial of each window, in the raw trials.
    for sid, seq, row, session in zip(occ["subject_id"], occ["sequence"],
                                      occ["row"], occ["session"]):
        trials = cohort.by_id(sid).trials
        symbol = trials[row, 1] + (trials[row, 2] * 3 if encode_reward else 0)
        assert symbol == seq[0]
        assert session == trials[row, 0]
        if contingency is not None:
            assert trials[row, 3] == contingency


class TestMultiContingency:
    params = CBASParams(num_arms=3, seq_len_max=2, criterion=6, criterion_order=1,
                        resample_number=20)

    def cohort(self):
        return Cohort([two_block_record(seed, n=40) for seed in range(6)])

    def test_positions_are_local_to_each_contingency(self):
        cohort = self.cohort()
        matrix = build_multicontingency_count_matrix(cohort, self.params)
        occ, subjects = sequence_occurrences(cohort, self.params, matrix.sequences,
                                             counted=matrix)
        assert len(subjects["subject_id"]) == len(cohort) * 2
        assert set(subjects["block"]) == {1, 2}

        cutoff = {(s, b): c for s, b, c in zip(subjects["subject_id"],
                                               subjects["block"],
                                               subjects["cutoff"])}
        for sid, block, start, row in zip(occ["subject_id"], occ["block"],
                                          occ["start"], occ["row"]):
            assert start <= cutoff[(sid, block)]
            assert cohort.by_id(sid).trials[row, 3] == block
        assert occ["start"].min() == 0
        assert occ["start"].max() < 40

    def test_blocks_restrict_the_tables(self):
        cohort = self.cohort()
        matrix = build_multicontingency_count_matrix(cohort, self.params, blocks=[1])
        occ, subjects = sequence_occurrences(cohort, self.params, matrix.sequences,
                                             blocks=[1], counted=matrix)
        assert set(occ["block"]) == {1}
        assert set(subjects["block"]) == {1}

    def test_checked_against_a_pipeline_result(self):
        cohort = self.cohort()
        result = run_cbas_multicontingency(cohort, [0, 1] * 3, self.params)
        occ, _ = sequence_occurrences(cohort, self.params, result.sequences[:5],
                                      counted=result)
        assert len(occ["start"]) > 0


class TestTables:
    params = CBASParams(num_arms=3, seq_len_max=3, criterion=40)

    def test_non_users_are_in_subjects_but_not_occurrences(self):
        cohort = Cohort([subject([0, 1, 0, 1], id="a"), subject([2, 2, 2, 2], id="b")])
        occ, subjects = sequence_occurrences(cohort, self.params, [(0, 1)],
                                             encode_reward=False)
        assert set(occ["subject_id"]) == {"a"}
        assert subjects["subject_id"].tolist() == ["a", "b"]

    def test_columns_have_equal_lengths_and_documented_dtypes(self):
        cohort = irregular_cohort()
        matrix = build_count_matrix(cohort, self.params)
        occ, subjects = sequence_occurrences(cohort, self.params, matrix.sequences[:10])
        for table in (occ, subjects):
            assert len({len(column) for column in table.values()}) == 1
        for key in ("start", "row", "session", ):
            assert occ[key].dtype == np.int64
        for key in ("subject_id", "block", "sequence"):
            assert occ[key].dtype == object
        assert subjects["reached"].dtype == bool
        assert subjects["cutoff"].dtype == np.int64

    def test_reordering_the_cohort_changes_no_subject_rows(self):
        cohort = irregular_cohort()
        seqs = build_count_matrix(cohort, self.params).sequences[:20]
        before = as_table(sequence_occurrences(cohort, self.params, seqs)[0])
        after = as_table(sequence_occurrences(cohort.reorder([3, 0, 4, 2, 1]),
                                              self.params, seqs)[0])
        assert before == after

    def test_numpy_entries_are_accepted(self):
        cohort = irregular_cohort()
        seqs = build_count_matrix(cohort, self.params).sequences[:3]
        as_numpy = [np.asarray(s) for s in seqs]
        a = sequence_occurrences(cohort, self.params, seqs)[0]
        b = sequence_occurrences(cohort, self.params, as_numpy)[0]
        assert a["sequence"].tolist() == b["sequence"].tolist()


class TestVerification:
    """Mismatches are refused, on a hypothesis space small enough that the sequence
    list alone would often fail to notice them."""

    params = CBASParams(num_arms=2, seq_len_max=2, criterion=40, resample_number=20)

    def cohort(self):
        rng = np.random.default_rng(1)
        return Cohort([subject(rng.integers(0, 2, 60), rng.integers(0, 2, 60),
                               np.repeat([1, 2, 3], 20),
                               np.where(rng.random(60) < 0.2, 0, 2), id=f"s{i}")
                       for i in range(6)])

    def run(self, cohort):
        return run_cbas_comparative(cohort, [0, 0, 0, 1, 1, 1], self.params,
                                    encode_reward=False)

    def test_the_matching_result_passes_silently(self, recwarn):
        cohort = self.cohort()
        result = self.run(cohort)
        sequence_occurrences(cohort, self.params, result.sequences,
                             encode_reward=False, counted=result)
        assert not [w for w in recwarn if issubclass(w.category, UserWarning)]

    @pytest.mark.parametrize("change", [
        dict(params=dict(criterion=39)),
        dict(params=dict(criterion_order=1, criterion=10)),
        dict(settings=dict(block_aware=True)),
        dict(settings=dict(contingency=None)),
    ])
    def test_different_settings_are_refused(self, change):
        cohort = self.cohort()
        result = self.run(cohort)
        params = replace(self.params, **change.get("params", {}))
        settings = dict(encode_reward=False, **change.get("settings", {}))
        matrix = build_count_matrix(cohort, self.params, encode_reward=False)
        for counted in (result, matrix):
            with pytest.raises(ValueError, match="does not match"):
                sequence_occurrences(cohort, params, result.sequences[:2],
                                     counted=counted, **settings)

    def test_a_filtered_cohort_is_refused(self):
        cohort = self.cohort()
        result = self.run(cohort)
        matrix = build_count_matrix(cohort, self.params, encode_reward=False)
        fewer = Cohort(list(cohort)[1:])
        for counted in (result, matrix):
            with pytest.raises(ValueError, match="does not match"):
                sequence_occurrences(fewer, self.params, result.sequences[:2],
                                     encode_reward=False, counted=counted)

    def test_a_result_without_a_digest_warns_and_falls_back(self):
        cohort = self.cohort()
        result = replace(self.run(cohort), counts_digest=None)
        with pytest.warns(UserWarning, match="no counts_digest"):
            sequence_occurrences(cohort, self.params, result.sequences[:2],
                                 encode_reward=False, counted=result)

    def test_counted_must_be_a_result_or_a_matrix(self):
        with pytest.raises(TypeError, match="CBASResult or a CountMatrix"):
            sequence_occurrences(self.cohort(), self.params, [(0,)],
                                 encode_reward=False, counted=object())

    @pytest.mark.parametrize("sequences,match", [
        ([], "empty"),
        ([(0,), (0,)], "repeated"),
        ([(0,), (1, (0,))], "mixes"),
        ([(0, 1, 0)], "seq_len_max"),
        ([(3,)], "never counted"),
    ])
    def test_malformed_sequences_are_refused(self, sequences, match):
        with pytest.raises(ValueError, match=match):
            sequence_occurrences(self.cohort(), self.params, sequences,
                                 encode_reward=False)

    def test_blocks_need_paired_entries(self):
        with pytest.raises(ValueError, match="blocks applies"):
            sequence_occurrences(self.cohort(), self.params, [(0,)],
                                 encode_reward=False, blocks=[1])
