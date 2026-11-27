"""CBAS parameter and result dataclasses."""

import hashlib
import json

import numpy as np
from dataclasses import dataclass


@dataclass
class CBASParams:
    """Configuration for a CBAS run.

    criterion and criterion_order together decide how much of each subject's
    stream is counted. criterion_order 0, the default, treats criterion as a
    fixed trial index and reproduces the behaviour pycbas has always had. Order 1
    stops each subject at its criterion-th reward, and order k at its
    criterion-th run of k consecutive rewarded choices, so subjects are matched on
    achievement rather than exposure and contribute different numbers of trials.
    See pycbas.criterion, and subject_criteria() for who fell short.
    """

    num_arms: int = 6
    seq_len_max: int = 6
    criterion: int = 800
    criterion_order: int = 0
    resample_number: int = 10_000
    alpha: float = 0.5
    gamma: float = 0.05
    centering: bool = False


@dataclass
class CBASResult:
    """Outcome of a CBAS run.

    `counts_digest` identifies the count matrix the result was computed from, by its
    counts, row ids and column labels. It lets a later step that recounts the same
    cohort, such as `sequence_occurrences`, confirm it counted exactly what the
    analysis counted. None for results built by hand.
    """

    sequences: list[tuple]
    test_stats: np.ndarray
    g_values: np.ndarray
    k_final: int
    significant_mask: np.ndarray
    k_history: list = None
    counts_digest: str | None = None

    @property
    def n_significant(self):
        return int(self.significant_mask.sum())


def _counts_digest(counts, subject_ids, sequences):
    """SHA-256 of a count matrix together with its row and column labels.

    Each part must mark its own end, or two different inputs could run together
    into the same bytes. The shape fixes the length of the counts, and the ids are
    a JSON array, which escapes any character an id holds and ends at its `]`.
    """
    counts = np.ascontiguousarray(counts, dtype="<f8")
    h = hashlib.sha256()
    h.update(repr(counts.shape).encode())
    h.update(counts.tobytes())
    h.update(json.dumps([str(i) for i in subject_ids]).encode())
    h.update(repr([_plain(entry) for entry in sequences]).encode())
    return h.hexdigest()


def _plain(value):
    """Builtin ints in place of numpy integers, whose repr differs across versions."""
    if isinstance(value, tuple):
        return tuple(_plain(v) for v in value)
    if isinstance(value, np.integer):
        return int(value)
    return value
