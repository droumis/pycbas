"""One definition of the two-sample studentized statistic, shared by both paths.

The step-down compares an observed statistic against bootstrap statistics with
`>=`, and that comparison is decided by exact equality far more often than it
looks. Sequences occurring in only a few subjects of one group have statistics
that are exact small rationals, and the bootstrap reproduces the identical value
in a large block of resamples. A sequence seen once in a handful of subjects of one
group and never in the other is the common case: its statistic is then fixed by that
pattern alone, so every such sequence shares one value, and any resample drawing the
same pattern reproduces it. In one real analysis that block was a tenth of 10,000
resamples.

So the observed and bootstrap statistics must agree *bitwise* whenever they
represent the same value. Two properties give that, and both matter:

1. **One code path.** Both callers use `sigma_from_sums` below, so no difference of
   formula or operation order can creep in. Computing the observed value as
   `std(ddof=1) / sqrt(n)` and squaring it, while the bootstrap formed
   `var / (n * (n - 1))`, put the two, two units in the last place apart and
   silently discarded that whole block of resamples.

2. **Order-independent accumulation.** The statistic is built from the raw sum
   and the raw sum of squares rather than from centred deviations. For count
   data those are integers, integers below 2**53 are exact in float64, and exact
   integer addition is associative -- so both sums are bitwise identical no
   matter which order the subjects are visited in. Centred deviations are not:
   a resample sums the same six ones in a different position among the zeros and
   lands an ULP away.

   The variance numerator is then formed as `n * sum_sq - sum * sum`, an exact
   integer difference, *before* any division. Written as `sum_sq - sum * sum / n`
   the division would round first and the exactness would be lost.

The result is that the statistic is a deterministic function of `(sum, sum_sq, n)`
per group, so identical multiset statistics give identical floating-point output
and the step-down's comparison needs no tolerance.

The exactness argument holds for non-negative integer counts, which is what CBAS
counts are. **Integer input is the only regime with a guarantee.** For non-integer input, such
as a rate-normalised matrix, the two paths still share a formula and are far closer
than they were, but the sums are no longer order-independent and bitwise equality is
lost: on a rate matrix with permuted draw order, over half the entries differ in the
last place, and the strict comparison in the step-down is then anti-conservative by
roughly 0.3% of exact ties. Callers who need the guarantee should pass the integer
count matrix; `t` is scale-invariant, so where the normalising denominator is common
to all rows the rate and count matrices give the same statistic anyway.

The exact range is bounded by `n * sum_sq < 2**53`, so with counts at most `c` the
requirement is `c < sqrt(2**53) / n`: about 9.0e5 at n=105 and 1.7e6 at n=55.
Observed counts in this project peak near 800, so the margin is over three orders
of magnitude.

Two further things this relies on, both easy to break by accident. The numba
kernels must not be compiled with `fastmath=True`, which would license
reassociation of the accumulation loops. And `compute_test_stats` must reduce over
a C-contiguous array, so that numpy accumulates in subject order; both are asserted
in `tests/test_moments.py`.
"""

import numpy as np

from ._numba import njit

__all__ = ["sigma_from_sums", "sigma_from_sums_scalar", "sem_from_sums"]

# The two functions below are the same arithmetic in the same order, one for
# whole arrays and one scalar for use inside the numba bootstrap kernels. The only
# deliberate difference is the clamp: `np.maximum(a, 0.0)` normalises -0.0 to +0.0
# where the scalar branch does not. Unreachable for integer counts, where the
# numerator is a non-negative exact integer by Cauchy-Schwarz. Numba
# cannot compile the array form and numpy cannot vectorise the scalar form, so
# there are two. They must stay in step: `tests/test_moments.py` asserts they
# agree bitwise, which is the guard against editing one and not the other.


def sigma_from_sums(sum0, sq0, n0, sum1, sq1, n1):
    """Pooled standard error from raw sums, elementwise over arrays."""
    a0 = (n0 * sq0 - sum0 * sum0) / (n0 * n0 * (n0 - 1.0))
    a1 = (n1 * sq1 - sum1 * sum1) / (n1 * n1 * (n1 - 1.0))
    a0 = np.maximum(a0, 0.0)
    a1 = np.maximum(a1, 0.0)
    return np.sqrt(a0 + a1)


@njit(inline="always")
def sigma_from_sums_scalar(sum0, sq0, n0, sum1, sq1, n1):
    """Scalar form of `sigma_from_sums`, for the bootstrap kernels."""
    a0 = (n0 * sq0 - sum0 * sum0) / (n0 * n0 * (n0 - 1.0))
    a1 = (n1 * sq1 - sum1 * sum1) / (n1 * n1 * (n1 - 1.0))
    if a0 < 0.0:
        a0 = 0.0
    if a1 < 0.0:
        a1 = 0.0
    return np.sqrt(a0 + a1)


def sem_from_sums(total, ssq, n):
    """One-sample standard error of the mean, from raw sums.

    The paired/one-sample counterpart of `sigma_from_sums`, for a sign-flip or
    other one-sample null. Same reasoning and the same numerator formed before any
    division, so an observed statistic and a null replication that represent the
    same value agree bitwise. Paired differences of integer counts are integers, so
    the guarantee applies to them.

    The motivating case is a sign-flip null, where flipping a subject whose paired
    difference is exactly zero leaves the statistic unchanged. That makes an exact
    tie, and if the observed and null paths use different formulas the tie is
    decided by rounding. Where most differences are zero, as in sparse sequence
    counts, those tie blocks are the majority of the null.
    """
    a = (n * ssq - total * total) / (n * n * (n - 1.0))
    return np.sqrt(np.maximum(a, 0.0))
