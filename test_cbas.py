"""Tests for the core CBAS implementation."""

import numpy as np
import pytest
from pathlib import Path

from pycbas import (
    CBASParams,
    load_subject_data,
    extract_choice_stream,
    enumerate_sequences,
    build_count_matrix,
    compute_test_stats,
    bootstrap_test_stats,
    romano_wolf_stepdown,
    find_k_fwer,
    run_cbas_comparative,
)

DATA_DIR = Path(__file__).parent / "igor_cbas" / "data"


# --- Unit tests with synthetic data ---


class TestEnumerateSequences:
    def test_simple_stream(self):
        stream = np.array([0, 1, 2, 3, 4, 5])
        counts = enumerate_sequences(stream, seq_len=2, criterion=6)
        assert counts[(0, 1)] == 1
        assert counts[(1, 2)] == 1
        assert len(counts) == 5

    def test_repeated_sequence(self):
        stream = np.array([0, 1, 0, 1, 0, 1])
        counts = enumerate_sequences(stream, seq_len=2, criterion=6)
        assert counts[(0, 1)] == 3
        assert counts[(1, 0)] == 2

    def test_criterion_truncates(self):
        stream = np.array([0, 1, 2, 3, 4, 5, 6, 7])
        counts = enumerate_sequences(stream, seq_len=2, criterion=4)
        assert (4, 5) not in counts
        assert len(counts) == 3

    def test_length_3(self):
        stream = np.array([1, 2, 3, 1, 2, 3])
        counts = enumerate_sequences(stream, seq_len=3, criterion=6)
        assert counts[(1, 2, 3)] == 2
        assert counts[(2, 3, 1)] == 1


class TestExtractChoiceStream:
    def test_encoding(self):
        data = np.array([
            [0, 3, 1, 2],  # arm=3, reward=1, conting=2 -> 3 + 1*6 = 9
            [0, 2, 0, 2],  # arm=2, reward=0, conting=2 -> 2
            [0, 5, 1, 0],  # conting=0, excluded
        ], dtype=np.int32)
        stream = extract_choice_stream(data, contingency=2, num_arms=6)
        assert list(stream) == [9, 2]

    def test_filters_contingency(self):
        data = np.array([
            [0, 1, 0, 0],
            [0, 2, 0, 0],
            [1, 3, 1, 2],
            [1, 4, 0, 2],
        ], dtype=np.int32)
        stream = extract_choice_stream(data, contingency=2, num_arms=6)
        assert len(stream) == 2


class TestComputeTestStats:
    def test_clear_group_difference(self):
        rng = np.random.default_rng(0)
        n0, n1, n_seq = 20, 20, 5
        counts = np.zeros((40, n_seq))
        counts[:20] = rng.poisson(10, (n0, n_seq))
        counts[20:] = rng.poisson(2, (n1, n_seq))

        group_indices = [np.arange(20), np.arange(20, 40)]
        stats = compute_test_stats(counts, group_indices)

        assert stats.shape == (n_seq * 2,)
        for i in range(n_seq):
            assert not np.isnan(stats[i * 2])  # positive direction populated
            assert np.isnan(stats[i * 2 + 1])  # negative direction is NaN

    def test_no_difference_gives_nan(self):
        counts = np.ones((10, 3))
        group_indices = [np.arange(5), np.arange(5, 10)]
        stats = compute_test_stats(counts, group_indices)
        assert np.all(np.isnan(stats))

    def test_directionality(self):
        rng = np.random.default_rng(1)
        counts = np.zeros((10, 1))
        counts[:5, 0] = rng.normal(1.0, 0.3, 5)  # group 0 lower
        counts[5:, 0] = rng.normal(5.0, 0.3, 5)  # group 1 higher
        group_indices = [np.arange(5), np.arange(5, 10)]
        stats = compute_test_stats(counts, group_indices)
        assert np.isnan(stats[0])  # positive direction NaN (group0 < group1)
        assert not np.isnan(stats[1])  # negative direction populated


class TestRomanoWolfStepdown:
    def test_clear_signal_gets_low_pvalue(self):
        rng = np.random.default_rng(42)
        test_stats = np.array([5.0, 3.0, np.nan, np.nan])
        null_matrix = rng.normal(0, 1, (1000, 4))
        null_matrix[:, 2:] = np.nan

        p_values = romano_wolf_stepdown(test_stats, null_matrix, k=1)
        assert p_values[0] < 0.01
        assert p_values[1] < 0.05
        assert np.isnan(p_values[2])
        assert np.isnan(p_values[3])

    def test_monotonicity(self):
        rng = np.random.default_rng(7)
        test_stats = np.array([4.0, 2.0, 1.5, 0.5])
        null_matrix = rng.normal(0, 1, (500, 4))

        p_values = romano_wolf_stepdown(test_stats, null_matrix, k=1)
        valid = p_values[~np.isnan(p_values)]
        for i in range(len(valid) - 1):
            assert valid[i] <= valid[i + 1]

    def test_k_greater_than_1_is_more_liberal(self):
        rng = np.random.default_rng(99)
        test_stats = np.array([3.5, 2.5, 1.8, 1.2])
        null_matrix = rng.normal(0, 1, (500, 4))

        p_k1 = romano_wolf_stepdown(test_stats, null_matrix, k=1)
        p_k2 = romano_wolf_stepdown(test_stats, null_matrix, k=2)

        assert np.all(p_k2[~np.isnan(p_k2)] <= p_k1[~np.isnan(p_k1)] + 1e-10)


class TestFindKFWER:
    def test_converges(self):
        rng = np.random.default_rng(12)
        n_tests = 20
        test_stats = np.concatenate([
            rng.normal(4, 0.5, 5),
            rng.normal(1, 0.5, 15),
        ])
        null_matrix = rng.normal(0, 1, (500, n_tests))

        g_values, k_final = find_k_fwer(test_stats, null_matrix, alpha=0.5, gamma=0.05)
        assert k_final >= 1
        assert g_values.shape == (n_tests,)
        assert np.sum(g_values < 0.5) > 0


# --- Integration tests with real data (lightweight) ---


def load_rat_data():
    """Load all rat data files and return (subjects_data, group_labels)."""
    if not DATA_DIR.exists():
        pytest.skip("Igor data directory not found")

    subjects_data = []
    group_labels = []

    for f in sorted(DATA_DIR.glob("*.txt")):
        name = f.stem
        if "Control" in name:
            label = 0
        elif "Lesion" in name:
            label = 1
        else:
            continue
        data = load_subject_data(f)
        subjects_data.append(data)
        group_labels.append(label)

    return subjects_data, np.array(group_labels)


class TestLoadSubjectData:
    def test_loads_file(self):
        if not DATA_DIR.exists():
            pytest.skip("Igor data directory not found")
        f = sorted(DATA_DIR.glob("*.txt"))[0]
        data = load_subject_data(f)
        assert data.ndim == 2
        assert data.shape[1] == 4
        assert data[:, 1].max() <= 5  # 6 arms: 0-5
        assert set(np.unique(data[:, 2])).issubset({0, 1})

    def test_extract_stream_correct_length(self):
        if not DATA_DIR.exists():
            pytest.skip("Igor data directory not found")
        f = sorted(DATA_DIR.glob("femaleControl*"))[0]
        data = load_subject_data(f)
        stream = extract_choice_stream(data, contingency=2, num_arms=6)
        assert len(stream) >= 800
        assert stream.max() < 12  # 6 arms * 2 (reward/no reward)


class TestBuildCountMatrix:
    def test_small_subset(self):
        if not DATA_DIR.exists():
            pytest.skip("Igor data directory not found")
        subjects_data, group_labels = load_rat_data()
        params = CBASParams(seq_len_max=2, criterion=100)
        sequences, count_matrix = build_count_matrix(subjects_data[:6], params)

        assert count_matrix.shape[0] == 6
        assert count_matrix.shape[1] == len(sequences)
        assert count_matrix.shape[1] > 0
        assert np.all(count_matrix >= 0)
        for seq in sequences:
            assert 1 <= len(seq) <= 2


class TestIntegrationSmall:
    """Run CBAS on a small subset with reduced parameters to validate the pipeline."""

    def test_pipeline_runs(self):
        if not DATA_DIR.exists():
            pytest.skip("Igor data directory not found")
        subjects_data, group_labels = load_rat_data()

        params = CBASParams(
            seq_len_max=2,
            criterion=200,
            resample_number=200,
        )
        result = run_cbas_comparative(subjects_data, group_labels, params)

        assert len(result.sequences) > 0
        assert result.g_values.shape[0] == len(result.sequences) * 2
        assert result.k_final >= 1
        assert result.significant_mask.shape[0] == len(result.sequences)

    def test_finds_some_significant(self):
        """With real group differences, CBAS should find significant sequences."""
        if not DATA_DIR.exists():
            pytest.skip("Igor data directory not found")
        subjects_data, group_labels = load_rat_data()

        params = CBASParams(
            seq_len_max=3,
            criterion=400,
            resample_number=500,
        )
        result = run_cbas_comparative(subjects_data, group_labels, params)
        assert result.n_significant > 0, "Expected to find significant sequences in rat data"

    def test_no_signal_with_same_group(self):
        """Comparing a group to itself should yield few/no significant sequences."""
        if not DATA_DIR.exists():
            pytest.skip("Igor data directory not found")
        subjects_data, group_labels = load_rat_data()

        control_indices = np.where(group_labels == 0)[0]
        control_data = [subjects_data[i] for i in control_indices]
        n = len(control_data)
        fake_labels = np.array([0] * (n // 2) + [1] * (n - n // 2))

        params = CBASParams(
            seq_len_max=2,
            criterion=200,
            resample_number=300,
        )
        result = run_cbas_comparative(control_data, fake_labels, params)
        n_seq = len(result.sequences)
        frac_significant = result.n_significant / n_seq if n_seq > 0 else 0
        assert frac_significant < 0.15, (
            f"Expected few false positives but got {result.n_significant}/{n_seq}"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
