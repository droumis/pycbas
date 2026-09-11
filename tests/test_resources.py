"""`estimate_resources` has to agree with what the pipeline actually allocates.

The estimate is what a user consults before committing to a long run, so being
wrong in the pessimistic direction is not harmless: it talks people out of the
default configuration. Two errors of that kind are pinned here.

The unchunked figure was computed over the enumerable sequence space rather than
over the count matrix's columns. Sequence space grows as A**L and the observed
count does not, so on a sparse space the estimate was wrong by two orders of
magnitude and contradicted the halving that README.md, docs/guide.md and
docs/algorithm.md all describe.

The time estimate omitted M entirely, so a run at five times the resample count
was quoted the same number of seconds.

The reported gigabyte figures are rounded to two decimals as part of the public
contract, so the expectations below are rounded the same way rather than compared
against full precision.
"""

import numpy as np
import pytest

from pycbas import CBASParams, bootstrap_test_stats, compute_test_stats
from pycbas.resources import estimate_resources

PER_ENTRY = 9  # float64 magnitude + int8 direction
GB = 1024**3


def test_chunked_matches_the_real_allocation():
    """The chunked path allocates exactly M x n_valid of magnitude plus direction."""
    est = estimate_resources(6, 6, n_observed=16376, resample_number=10000)
    expected = round(10000 * 16376 * PER_ENTRY / GB, 2)
    assert est["memory_chunked_gb"] == expected


def test_unchunked_is_about_twice_chunked_not_a_multiple_of_sequence_space():
    """Chunking avoids one full-width copy, so peak roughly halves.

    The regression this pins: using the enumerable sequence space for the
    unchunked figure. At 6 arms with reward encoding and L=6 that space is over
    three million while the observed count is a few tens of thousands, so the bug
    inflated the estimate by more than a hundredfold.
    """
    est = estimate_resources(6, 6, n_observed=16376, resample_number=10000)
    ratio = est["memory_full_null_gb"] / est["memory_chunked_gb"]
    assert 1.9 < ratio < 2.1, ratio
    assert est["total_sequences"] > 100 * 16376, "precondition: the space is sparse"
    assert est["memory_full_null_gb"] < 10, est["memory_full_null_gb"]


def test_memory_and_time_are_both_linear_in_resample_number():
    small = estimate_resources(6, 6, n_observed=16376, resample_number=10000)
    large = estimate_resources(6, 6, n_observed=16376, resample_number=50000)
    assert large["memory_chunked_gb"] / small["memory_chunked_gb"] == pytest.approx(5, rel=1e-2)
    assert large["est_time_seconds"] / small["est_time_seconds"] == pytest.approx(5, rel=1e-2)


def test_time_is_linear_in_hypothesis_count():
    a = estimate_resources(6, 6, n_observed=8000, resample_number=10000)
    b = estimate_resources(6, 6, n_observed=16000, resample_number=10000)
    assert b["est_time_seconds"] / a["est_time_seconds"] == pytest.approx(2, rel=1e-2)


def test_time_keeps_its_calibration_point():
    """The one measured point the model is fitted to must come back unchanged."""
    est = estimate_resources(6, 6, n_observed=16376, resample_number=10000)
    assert est["est_time_seconds"] == pytest.approx(12.3, abs=0.05)


def test_sequence_space_is_the_fallback_when_the_count_is_unknown():
    """Without `n_observed` there is nothing to go on but the worst case."""
    est = estimate_resources(2, 4, resample_number=100000)
    S = sum(4**length for length in range(1, 5))  # A = 2 * 2 with reward encoding
    assert est["total_sequences"] == S
    assert est["memory_chunked_gb"] == round(100000 * S * PER_ENTRY / GB, 2)


def test_model_matches_the_arrays_the_pipeline_really_builds():
    """Check the formula's shape against real allocations, at small size.

    The claim is structural: the null and its directions are sized by the count
    matrix's columns, not by the enumerable space, and the step-down's submatrix by
    the statistics that are actually defined. Compared in bytes, because at this
    size the public gigabyte figures round to zero.
    """
    rng = np.random.default_rng(0)
    n0 = n1 = 9
    n_seq = 120
    counts = rng.integers(0, 4, size=(n0 + n1, n_seq)).astype(np.float64)
    grp = [np.arange(n0), np.arange(n0, n0 + n1)]
    M = 40

    null, dirs = bootstrap_test_stats(counts, grp, CBASParams(resample_number=M),
                                      rng=np.random.default_rng(0))
    assert null.shape == (M, n_seq), "null is sized by observed columns, not by A**L"
    assert dirs.shape == (M, n_seq)
    assert null.dtype == np.float64 and dirs.dtype == np.int8, "the 9 bytes an entry"

    stats = compute_test_stats(counts, grp)
    n_valid = int((~np.isnan(stats)).sum())
    assert 0 < n_valid <= n_seq

    real_chunked = M * n_valid * PER_ENTRY
    real_full = M * (n_seq + n_valid) * PER_ENTRY

    # Recompute the model in bytes, mirroring estimate_resources.
    model_chunked = M * n_seq * PER_ENTRY
    model_full = M * (n_seq + n_seq) * PER_ENTRY

    # The model assumes every sequence yields a defined statistic, so it is an
    # upper bound, and a tight one whenever most statistics are defined.
    assert model_chunked >= real_chunked
    assert model_full >= real_full
    assert model_chunked < 1.3 * real_chunked
    assert model_full < 1.3 * real_full


def test_hypothesis_sets_multiply_the_reported_ceiling():
    """A multi-contingency run counts the sequence space once per block.

    With the multiplier missing, the ceiling came out below the observed count and
    the GUI reported "3,542 of 1,884 possible", which cannot be read any way but
    wrong.
    """
    one = estimate_resources(6, 3, n_observed=3542, n_hypothesis_sets=1)
    three = estimate_resources(6, 3, n_observed=3542, n_hypothesis_sets=3)
    assert three["total_sequences"] == 3 * one["total_sequences"]
    assert one["observed_sequences"] > one["total_sequences"], "precondition"
    assert three["observed_sequences"] <= three["total_sequences"]


def test_hypothesis_sets_do_not_change_the_memory_when_the_count_is_known():
    """The allocation follows the observed count, not the ceiling."""
    a = estimate_resources(6, 3, n_observed=3542, n_hypothesis_sets=1)
    b = estimate_resources(6, 3, n_observed=3542, n_hypothesis_sets=3)
    assert a["memory_chunked_gb"] == b["memory_chunked_gb"]
    assert a["est_time_seconds"] == b["est_time_seconds"]
