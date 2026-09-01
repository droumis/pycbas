"""
pycbas -- Choice-Wide Behavioral Association Study.

A Python implementation of the CBAS algorithm for identifying behavioral
sequences that differ between experimental groups or correlate with a
continuous measure. Uses Romano-Wolf step-down for multiple comparison
correction and k-FWER iteration for false discovery proportion control.

Reference: Kastner et al., "Choice-Wide Behavioral Association Study"
(2026 preprint) https://www.biorxiv.org/content/10.1101/2024.02.26.582115v4
"""

from .params import CBASParams, CBASResult
from .io import (load_subject_data, extract_choice_stream, extract_choice_streams_by_block,
                 enumerate_sequences, enumerate_sequences_block_aware)
from .core import build_count_matrix, compute_test_stats, compute_test_stats_correlative
from .bootstrap import bootstrap_test_stats, bootstrap_test_stats_correlative
from .stepdown import (
    romano_wolf_stepdown,
    find_k_fwer,
    find_k_fwer_k1,
    find_k_fwer_chunked,
)
from .resources import estimate_resources, print_resource_estimate
from .pipeline import (run_cbas_comparative, run_cbas_correlative,
                       run_cbas_multicontingency)
from .contingency import (load_subject_data_with_contingencies,
                          load_cohort_with_contingencies,
                          shared_contingency_blocks,
                          build_multicontingency_count_matrix)
from .criterion import criterion_trial, reached_criterion
from .core import subject_criteria

__all__ = [
    "CBASParams",
    "CBASResult",
    "load_subject_data",
    "extract_choice_stream",
    "extract_choice_streams_by_block",
    "enumerate_sequences",
    "enumerate_sequences_block_aware",
    "build_count_matrix",
    "compute_test_stats",
    "compute_test_stats_correlative",
    "bootstrap_test_stats",
    "bootstrap_test_stats_correlative",
    "romano_wolf_stepdown",
    "find_k_fwer",
    "find_k_fwer_k1",
    "find_k_fwer_chunked",
    "run_cbas_comparative",
    "run_cbas_correlative",
    "run_cbas_multicontingency",
    "load_subject_data_with_contingencies",
    "load_cohort_with_contingencies",
    "shared_contingency_blocks",
    "build_multicontingency_count_matrix",
    "criterion_trial",
    "reached_criterion",
    "subject_criteria",
    "estimate_resources",
    "print_resource_estimate",
]
