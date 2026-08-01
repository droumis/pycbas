"""CBAS parameter and result dataclasses."""

import numpy as np
from dataclasses import dataclass


@dataclass
class CBASParams:
    num_arms: int = 6
    seq_len_max: int = 6
    criterion: int = 800
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
