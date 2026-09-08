"""The observed and bootstrap statistics must agree bitwise, not just closely.

The step-down compares an observed statistic against bootstrap statistics with
`>=`, and for rare sequences that comparison is an equality test on a large block
of resamples at once. A sequence occurring once in a handful of subjects of one
group has a statistic that is an exact small rational, and the bootstrap
reproduces the identical value whenever it draws the same pattern.

The bug these tests exist to prevent: `compute_test_stats` formed
`std(ddof=1) / sqrt(n)` and squared it while the bootstrap formed
`var / (n * (n - 1))`. Algebraically identical, two units in the last place
apart. In one real analysis that discarded a tenth of 10,000 resamples for a block
of tied sequences, moved their adjusted p-value by 0.105, and turned a block that
should have been rejected into false positives. The error is always
anti-conservative, so it cannot be dismissed as noise.

See pycbas/_moments.py for why raw sums rather than centred deviations.
"""

import sys

import numpy as np
import pytest

from pycbas import CBASParams, compute_test_stats
from pycbas._moments import (sem_from_sums, sigma_from_sums,
                             sigma_from_sums_scalar, tie_rtol_for)
from pycbas.bootstrap import (_bootstrap_chunk_into, _bootstrap_parallel,
                              bootstrap_test_stats)


def test_array_and_scalar_forms_agree_bitwise():
    """The two forms in _moments.py are the same arithmetic; keep them so."""
    rng = np.random.default_rng(0)
    n0, n1 = 55, 50
    for _ in range(200):
        x = rng.integers(0, 40, size=n0).astype(np.float64)
        y = rng.integers(0, 40, size=n1).astype(np.float64)
        args = (x.sum(), (x * x).sum(), float(n0),
                y.sum(), (y * y).sum(), float(n1))
        got = sigma_from_sums(*args)
        want = sigma_from_sums_scalar(*args)
        assert got == want, (args, got, want)


def test_observed_matches_bootstrap_bitwise_on_identity_resample():
    """A resample that draws exactly the real groups must reproduce the observed.

    This is the invariant the step-down's `>=` depends on. Equality here is
    required, not approximate agreement.
    """
    rng = np.random.default_rng(1)
    n0, n1 = 9, 7
    n_seq = 400
    # Small integer counts, many zeros: the regime that produces exact rational
    # statistics and therefore the ties that matter.
    counts = rng.integers(0, 3, size=(n0 + n1, n_seq)).astype(np.float64)
    counts[:, ::5] = 0.0
    counts[0, ::5] = 1.0
    grp = [np.arange(n0), np.arange(n0, n0 + n1)]

    observed = compute_test_stats(counts, grp)

    out = np.full((1, 2 * n_seq), -np.inf)
    out_dir = np.full((1, 2 * n_seq), np.int8(-1), dtype=np.int8)
    _bootstrap_chunk_into(
        np.ascontiguousarray(counts), grp[0].reshape(1, -1),
        grp[1].reshape(1, -1), n0, n1, np.zeros(n_seq),
        np.arange(2 * n_seq, dtype=np.int64), 0, 1, out, out_dir)

    defined = ~np.isnan(observed)
    assert defined.sum() > 0
    assert np.array_equal(observed[defined], out[0][defined])


def test_order_independence_under_permuted_resample():
    """The half of the fix that the identity resample cannot test.

    `test_observed_matches_bootstrap_bitwise_on_identity_resample` visits subjects
    in the same order `compute_test_stats` sums them, so it would still pass if the
    bootstrap went back to centred deviations. This one permutes the draw order,
    which is what a real resample does, and is the test that fails if raw sums are
    replaced by anything order-dependent.
    """
    rng = np.random.default_rng(3)
    n0, n1 = 55, 50
    n_seq = 300
    counts = np.zeros((n0 + n1, n_seq))
    # The regime that produces exact rational statistics: a handful of ones.
    for s in range(n_seq):
        counts[rng.choice(n0, size=rng.integers(1, 8), replace=False), s] = 1.0
    grp = [np.arange(n0), np.arange(n0, n0 + n1)]
    observed = compute_test_stats(counts, grp)
    defined = ~np.isnan(observed)

    for _ in range(15):
        b0 = rng.permutation(grp[0]).reshape(1, -1)
        b1 = rng.permutation(grp[1]).reshape(1, -1)
        out = np.full((1, 2 * n_seq), -np.inf)
        out_dir = np.full((1, 2 * n_seq), np.int8(-1), dtype=np.int8)
        _bootstrap_chunk_into(
            np.ascontiguousarray(counts), b0, b1, n0, n1, np.zeros(n_seq),
            np.arange(2 * n_seq, dtype=np.int64), 0, 1, out, out_dir)
        # Same multiset in a different order is the same statistic, so the
        # step-down's `>=` must hold with nothing to spare.
        assert np.array_equal(observed[defined], out[0][defined])


def test_non_chunked_kernel_matches_too():
    """`_bootstrap_parallel` got the same edit and needs the same guarantee."""
    rng = np.random.default_rng(4)
    n0, n1 = 12, 9
    n_seq = 150
    counts = rng.integers(0, 3, size=(n0 + n1, n_seq)).astype(np.float64)
    grp = [np.arange(n0), np.arange(n0, n0 + n1)]
    observed = compute_test_stats(counts, grp)

    params = CBASParams(resample_number=1)
    null, _dirs = bootstrap_test_stats(counts, grp, params,
                                       rng=np.random.default_rng(0))
    # magnitude only, so compare against |t| on the defined side
    mag = np.where(np.isnan(observed[0::2]), observed[1::2], observed[0::2])
    identity = _bootstrap_parallel(
        np.ascontiguousarray(counts), grp[0].reshape(1, -1),
        grp[1].reshape(1, -1), n0, n1, n_seq, 1, np.zeros(n_seq))[0][0]
    defined = ~np.isnan(mag)
    assert np.array_equal(mag[defined], identity[defined])


def test_kernels_are_not_compiled_with_fastmath():
    """`fastmath=True` would license reassociation and void the whole argument.

    Checks the compiled target options rather than the source text, so that
    mentioning the word in a comment does not pass or fail it.
    """
    for fn in (_bootstrap_chunk_into, _bootstrap_parallel,
               sigma_from_sums_scalar):
        opts = getattr(fn, "targetoptions", None)
        if opts is None:
            pytest.skip("numba not installed; kernels are plain Python")
        assert not opts.get("fastmath", False), fn


def test_compute_test_stats_reduces_over_contiguous_array():
    """numpy must accumulate in subject order, which needs C-contiguity.

    A Fortran-ordered input would switch numpy to a different summation and break
    agreement for non-integer matrices. `compute_test_stats` copies to C order for
    that reason; this pins the behaviour, not the implementation detail.
    """
    rng = np.random.default_rng(5)
    n0, n1, n_seq = 20, 18, 120
    counts = rng.integers(0, 4, size=(n0 + n1, n_seq)).astype(np.float64)
    grp = [np.arange(n0), np.arange(n0, n0 + n1)]
    c_order = compute_test_stats(np.ascontiguousarray(counts), grp)
    f_order = compute_test_stats(np.asfortranarray(counts), grp)
    both = ~np.isnan(c_order) & ~np.isnan(f_order)
    assert np.array_equal(c_order[both], f_order[both])


@pytest.mark.parametrize("scale", [1, 7, 113])
def test_variance_numerator_is_exact_for_integer_counts(scale):
    """`n * sq - sum * sum` must be an exact integer, which is what buys the
    order independence."""
    rng = np.random.default_rng(2)
    n = 41
    x = (rng.integers(0, 50, size=n) * scale).astype(np.float64)
    total, sq = x.sum(), (x * x).sum()
    numerator = n * sq - total * total
    assert numerator == float(round(numerator))
    # and independent of the order the subjects are visited in
    for _ in range(20):
        y = rng.permutation(x)
        assert y.sum() == total
        assert (y * y).sum() == sq


def test_tie_rtol_is_zero_for_integer_matrices_and_positive_otherwise():
    """Integer input is exact, so the tolerance must be exactly zero there."""
    rng = np.random.default_rng(6)
    counts = rng.integers(0, 500, size=(105, 200)).astype(np.float64)
    assert tie_rtol_for(counts) == 0.0
    assert tie_rtol_for(counts.astype(np.int64)) == 0.0
    # Detection is conservative by design: integers over a power of two are in
    # fact safe, but are not recognised, and get a harmless slack instead.
    assert 0.0 < tie_rtol_for(counts / 1024.0) < 1e-9

    rates = counts / counts.sum(axis=1, keepdims=True).clip(min=1)
    rtol = tie_rtol_for(rates)
    assert rtol > 0.0
    # derived, not tuned: far above one ULP, far below the spacing of distinct
    # statistics, which is order 1e-4 relative on real data
    assert 1e-15 < rtol < 1e-9


def test_rate_matrix_ties_are_counted_with_the_derived_tolerance():
    """The regime with no exactness guarantee must still not drop its own ties.

    A rate-normalised matrix loses order-independence, so an observed statistic
    and a null replication representing the same value can differ in the last
    places. With tie_rtol=0 the step-down discards them; with the derived value it
    does not. This pins that the derived tolerance is doing its job.
    """
    from pycbas.stepdown import _stepdown_core_directional

    rng = np.random.default_rng(7)
    n0, n1, n_seq = 30, 28, 240
    counts = np.zeros((n0 + n1, n_seq))
    for s in range(n_seq):
        counts[rng.choice(n0, size=rng.integers(1, 5), replace=False), s] = 1.0
    denom = counts.sum(axis=1, keepdims=True)
    denom[denom == 0] = 1.0
    rates = counts / denom
    grp = [np.arange(n0), np.arange(n0, n0 + n1)]

    observed = compute_test_stats(rates, grp)
    defined = ~np.isnan(observed)
    assert defined.sum() > 0

    # a permuted-order resample of the real groups: same values, same statistic
    out = np.full((1, 2 * n_seq), -np.inf)
    out_dir = np.full((1, 2 * n_seq), np.int8(-1), dtype=np.int8)
    _bootstrap_chunk_into(
        np.ascontiguousarray(rates), rng.permutation(grp[0]).reshape(1, -1),
        rng.permutation(grp[1]).reshape(1, -1), n0, n1, np.zeros(n_seq),
        np.arange(2 * n_seq, dtype=np.int64), 0, 1, out, out_dir)

    rtol = tie_rtol_for(rates)
    obs = observed[defined]
    nul = out[0][defined]
    strict_dropped = int((nul < obs).sum())
    tolerant_dropped = int((nul < obs - np.abs(obs) * rtol).sum())
    # the point of the tolerance: it recovers what strict comparison discards
    assert tolerant_dropped == 0
    assert strict_dropped > 0, "expected the rate matrix to expose the problem"


def test_paired_one_sample_form_is_order_independent():
    """`sem_from_sums` is the one-sample form, for a paired or sign-flip null.

    Flipping the sign of a subject whose paired difference is zero leaves the
    statistic unchanged, which is an exact tie. Integer differences must therefore
    give a bitwise identical standard error whatever the order.
    """
    rng = np.random.default_rng(8)
    n, n_seq = 50, 400
    D = np.zeros((n, n_seq))
    for s in range(n_seq):
        m = int(rng.integers(1, 6))
        D[rng.choice(n, size=m, replace=False), s] = rng.integers(-3, 4, size=m)
    total, ssq = D.sum(axis=0), (D * D).sum(axis=0)
    base = sem_from_sums(total, ssq, n)
    for _ in range(10):
        P = rng.permutation(n)
        Dp = D[P]
        assert np.array_equal(Dp.sum(axis=0), total)
        assert np.array_equal((Dp * Dp).sum(axis=0), ssq)
        assert np.array_equal(sem_from_sums(Dp.sum(axis=0),
                                           (Dp * Dp).sum(axis=0), n), base)
