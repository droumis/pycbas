"""CBAS parameter and result dataclasses."""

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
    sequences: list[tuple]
    test_stats: np.ndarray
    g_values: np.ndarray
    k_final: int
    significant_mask: np.ndarray
    k_history: list = None

    @property
    def n_significant(self):
        return int(self.significant_mask.sum())
