"""Tests for subject identity, cohort assembly, and the labelled count matrix.

These types exist so that nothing per-subject is aligned by position, so the tests
are mostly about what happens when something would have to be: a reordering, a
filter, a label supplied in a different order from the rows it describes.
"""

import numpy as np
import pytest

from pycbas import (CountMatrix, Cohort, Subject, default_group_coder, load_cohort,
                    load_subject, resolve_labels)


def trials(n=40, seed=0, condition=2, n_sessions=2):
    rng = np.random.default_rng(seed)
    sessions = np.repeat(np.arange(1, n_sessions + 1), n // n_sessions)
    return np.column_stack([
        sessions,
        rng.integers(0, 3, len(sessions)),
        rng.integers(0, 2, len(sessions)),
        np.full(len(sessions), condition),
    ]).astype(np.int64)


def subject(name, seed=0, **meta):
    return Subject(id=name, trials=trials(seed=seed), meta=meta)


def write_cohort(directory, names=("an1", "an2", "an10")):
    for i, name in enumerate(names):
        rows = ["1,%d,%d,2" % (i % 3, i % 2) for _ in range(30)]
        (directory / f"{name}.txt").write_text("\n".join(rows) + "\n")
    return directory


class TestSubject:
    def test_columns_are_views_of_one_array(self):
        s = subject("a")
        assert len(s) == len(s.trials)
        np.testing.assert_array_equal(s.session, s.trials[:, 0])
        np.testing.assert_array_equal(s.choice, s.trials[:, 1])
        np.testing.assert_array_equal(s.reward, s.trials[:, 2])
        np.testing.assert_array_equal(s.condition, s.trials[:, 3])

    def test_wrong_shape_is_refused_with_the_layout_named(self):
        with pytest.raises(ValueError, match="session, choice, reward, condition"):
            Subject(id="a", trials=np.zeros((10, 3)))

    def test_block_methods_refuse_a_single_contingency_subject(self):
        """The raw fourth column is not a validated block structure, so say so."""
        with pytest.raises(ValueError, match="no contingency block structure"):
            subject("a").alternation_blocks()

    def test_per_session_splits_do_not_span_sessions(self):
        s = subject("a", seed=3)
        blocks = s.reward_blocks_for(2)
        assert len(blocks) == len(np.unique(s.session))
        assert sum(len(b) for b in blocks) == len(s)

    def test_symbol_blocks_encode_reward_like_the_stream_does(self):
        s = subject("a", seed=4)
        plain = np.concatenate(s.symbol_blocks_for(2, num_arms=3, encode_reward=False))
        coded = np.concatenate(s.symbol_blocks_for(2, num_arms=3, encode_reward=True))
        np.testing.assert_array_equal(coded, plain + s.reward * 3)


class TestCohort:
    def test_duplicate_ids_are_refused(self):
        with pytest.raises(ValueError, match="unique"):
            Cohort([subject("a"), subject("a", seed=1)])

    def test_bare_arrays_are_refused_by_name(self):
        with pytest.raises(TypeError, match="load_subject"):
            Cohort([trials()])

    def test_addressable_by_position_and_by_id(self):
        cohort = Cohort([subject("a"), subject("b", seed=1)])
        assert cohort[0].id == "a"
        assert cohort.by_id("b").id == "b"
        with pytest.raises(KeyError):
            cohort.by_id("nobody")

    def test_an_id_is_not_a_position(self):
        """Ids are often numbers, so one syntax must not answer two questions."""
        cohort = Cohort([subject("2"), subject("0", seed=1), subject("1", seed=2)])
        assert cohort[0].id == "2", "indexing is positional"
        assert cohort.by_id("0").id == "0"
        with pytest.raises(TypeError, match="by_id"):
            cohort["0"]

    def test_ids_are_strings_however_they_were_given(self):
        """Cohort tables number their animals, so an int id must not lose the lookup."""
        cohort = Cohort([Subject(id=200, trials=trials()),
                         Subject(id=201, trials=trials(seed=1))])
        assert cohort.ids == ["200", "201"]
        assert cohort.by_id(200).id == "200"
        assert cohort.by_id("200").id == "200"

    def test_a_slice_is_a_cohort(self):
        cohort = Cohort([subject(n, seed=i) for i, n in enumerate("abcd")])
        assert isinstance(cohort[:2], Cohort)
        assert cohort[:2].ids == ["a", "b"]

    def test_filter_keeps_meta_with_its_subject(self):
        cohort = Cohort([subject("a", genotype="WT"), subject("b", seed=1, genotype="Mut"),
                         subject("c", seed=2, genotype="WT")])
        wt = cohort.filter(genotype="WT")
        assert wt.ids == ["a", "c"]
        assert [s.meta["genotype"] for s in wt] == ["WT", "WT"]

    def test_filter_accepts_several_values_and_a_predicate(self):
        cohort = Cohort([subject("a", g="WT"), subject("b", seed=1, g="Mut"),
                         subject("c", seed=2, g="KO")])
        assert cohort.filter(g={"WT", "KO"}).ids == ["a", "c"]
        assert cohort.filter(lambda s: s.id != "b").ids == ["a", "c"]

    def test_reorder_carries_identity(self):
        cohort = Cohort([subject(n, seed=i) for i, n in enumerate("abc")])
        moved = cohort.reorder([2, 0, 1])
        assert moved.ids == ["c", "a", "b"]
        np.testing.assert_array_equal(moved[0].trials, cohort.by_id("c").trials)

    def test_labels_from_a_meta_column(self):
        cohort = Cohort([subject("a", lesion=0), subject("b", seed=1, lesion=1)])
        np.testing.assert_array_equal(cohort.labels_from("lesion"), [0, 1])

    def test_labels_from_reports_unusable_values_instead_of_dropping(self):
        """Dropping here would change the cohort the caller thought it assembled."""
        cohort = Cohort([subject("a", lesion=0), subject("b", seed=1, lesion=None),
                         subject("c", seed=2, lesion=1)])
        with pytest.raises(ValueError, match=r"no usable 'lesion' value.*\['b'\]"):
            cohort.labels_from("lesion")

    def test_labels_from_takes_a_coder_for_string_values(self):
        cohort = Cohort([subject("a", lesion="control"), subject("b", seed=1, lesion="lesion")])
        coder = {"control": 0, "lesion": 1}.get
        np.testing.assert_array_equal(cohort.labels_from("lesion", coder), [0, 1])

    def test_covariate_from_rejects_a_non_number(self):
        cohort = Cohort([subject("a", score=1.5), subject("b", seed=1, score="high")])
        with pytest.raises(ValueError, match="not a number"):
            cohort.covariate_from("score")


class TestResolveLabels:
    """Labels follow ids, not positions, which is what makes a reordering harmless."""

    def test_a_mapping_follows_the_ids_given(self):
        got = resolve_labels({"a": 0, "b": 1, "c": 0}, ["c", "a", "b"])
        np.testing.assert_array_equal(got, [0, 0, 1])

    def test_a_sequence_is_read_in_cohort_order_then_reindexed(self):
        got = resolve_labels([0, 1, 0], ids=["c", "a", "b"], cohort_ids=["a", "b", "c"])
        np.testing.assert_array_equal(got, [0, 0, 1])

    def test_a_missing_id_is_named(self):
        with pytest.raises(ValueError, match=r"\['c'\]"):
            resolve_labels({"a": 0, "b": 1}, ["a", "b", "c"])

    def test_a_length_mismatch_is_named(self):
        with pytest.raises(ValueError, match="2 labels for 3 subjects"):
            resolve_labels([0, 1], ["a", "b", "c"])


class TestCountMatrix:
    def matrix(self):
        return CountMatrix(counts=np.arange(6, dtype=float).reshape(3, 2),
                           subject_ids=["a", "b", "c"],
                           sequences=[(0,), (1,)])

    def test_label_counts_must_match_the_array(self):
        with pytest.raises(ValueError, match="2 subject ids for 3 rows"):
            CountMatrix(np.zeros((3, 2)), ["a", "b"], [(0,), (1,)])
        with pytest.raises(ValueError, match="1 sequences for 2 columns"):
            CountMatrix(np.zeros((3, 2)), ["a", "b", "c"], [(0,)])

    def test_row_by_id(self):
        np.testing.assert_array_equal(self.matrix().row("b"), [2.0, 3.0])
        with pytest.raises(KeyError):
            self.matrix().row("nobody")

    def test_reorder_moves_rows_and_ids_together(self):
        moved = self.matrix().reorder([2, 0, 1])
        assert moved.subject_ids == ["c", "a", "b"]
        np.testing.assert_array_equal(moved.row("c"), [4.0, 5.0])

    def test_select_block_refuses_a_single_contingency_matrix(self):
        """Returning an empty selection here is the quiet wrong answer to avoid."""
        with pytest.raises(ValueError, match="no contingency blocks to select"):
            self.matrix().select_block(1)

    def test_select_block_keeps_one_block(self):
        m = CountMatrix(np.arange(9, dtype=float).reshape(3, 3), ["a", "b", "c"],
                        [(1, (0,)), (1, (2,)), (2, (0,))])
        one = m.select_block(1)
        assert one.sequences == [(1, (0,)), (1, (2,))]
        assert one.shape == (3, 2)
        with pytest.raises(ValueError, match="no columns for block 7"):
            m.select_block(7)

    def test_column_labels_decode(self):
        m = CountMatrix(np.zeros((1, 2)), ["a"], [(2,), (8,)])
        assert m.column_labels(num_arms=6) == ["3", "3*"]


class TestDefaultGroupCoder:
    """The vocabulary that decides group membership, shared by the library and the app.

    Untested it would be free to drift, and it is the one place where a cohort table's
    wording turns into a 0 or a 1.
    """

    @pytest.mark.parametrize("value", [0, "0", "control", "Control", "CTRL", "sham",
                                       "WT", "wildtype", "  control  ", False])
    def test_group_zero_words(self, value):
        assert default_group_coder(value) == 0

    @pytest.mark.parametrize("value", [1, "1", "lesion", "Hippocampal Lesion", "KO",
                                       "knockout", "mutant", "experimental", True])
    def test_group_one_words(self, value):
        assert default_group_coder(value) == 1

    @pytest.mark.parametrize("value", [None, "", "   ", "unknown", "genotype", 2, -1,
                                       0.5, "n/a"])
    def test_unusable_values_are_not_guessed_at(self, value):
        # None rather than a default group: a value that says nothing is not evidence.
        assert default_group_coder(value) is None

    def test_substring_fallback_places_a_descriptive_label(self):
        # What the fallback is for: real tables say "Hippocampal Lesion", not "lesion".
        assert default_group_coder("Hippocampal Lesion") == 1
        assert default_group_coder("sham control surgery") == 0

    def test_the_substring_fallback_reads_words_and_not_sentences(self):
        """A documented limitation, pinned so it is a known cost rather than a surprise.

        A table that spells out "no lesion evident" is read as lesion, because the
        fallback looks for the word. The cohorts this vocabulary was built for leave
        that field blank, which reads as unusable instead. A table that phrases its
        groups as sentences needs an explicit coder.
        """
        assert default_group_coder("no lesion evident") == 1
        assert default_group_coder("lesion-free") == 1


class TestLoaders:
    def test_load_subject_defaults_its_id_to_the_stem(self, tmp_path):
        write_cohort(tmp_path, ["an7"])
        s = load_subject(tmp_path / "an7.txt")
        assert s.id == "an7"
        assert len(s) == 30

    def test_load_cohort_orders_a_directory_numerically(self, tmp_path):
        write_cohort(tmp_path)
        assert load_cohort(tmp_path).ids == ["an1", "an2", "an10"]

    def test_load_cohort_keeps_an_explicit_file_order(self, tmp_path):
        write_cohort(tmp_path)
        paths = [tmp_path / "an10.txt", tmp_path / "an1.txt"]
        assert load_cohort(paths).ids == ["an10", "an1"]

    def test_load_cohort_takes_meta_by_id_or_in_file_order(self, tmp_path):
        write_cohort(tmp_path)
        by_id = load_cohort(tmp_path, meta={"an2": {"lesion": 1}})
        assert by_id.by_id("an2").meta == {"lesion": 1}
        assert by_id.by_id("an1").meta == {}
        in_order = load_cohort(tmp_path, meta=[{"lesion": 0}, {"lesion": 1}, {"lesion": 0}])
        np.testing.assert_array_equal(in_order.labels_from("lesion"), [0, 1, 0])

    def test_explicit_ids_disambiguate_repeated_filenames(self, tmp_path):
        write_cohort(tmp_path, ["an1"])
        other = tmp_path / "other"
        other.mkdir()
        write_cohort(other, ["an1"])
        paths = [tmp_path / "an1.txt", other / "an1.txt"]
        with pytest.raises(ValueError, match="unique"):
            load_cohort(paths)
        assert load_cohort(paths, ids=["first", "second"]).ids == ["first", "second"]

    def test_an_empty_source_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="no files matching"):
            load_cohort(tmp_path)
