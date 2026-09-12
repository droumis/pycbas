"""
pycbas -- Choice-Wide Behavioral Association Study.

A Python implementation of the CBAS algorithm for identifying behavioral
sequences that differ between experimental groups or correlate with a
continuous measure. Uses Romano-Wolf step-down for multiple comparison
correction and k-FWER iteration for false discovery proportion control.

Reference: Kastner et al., "Choice-Wide Behavioral Association Study"
Nature Communications (2026) https://www.nature.com/articles/s41467-026-76681-3
"""

#: Single source of truth for the version; `pyproject.toml` reads it from here via
#: hatch. Worth being able to check from Python: 0.2.0 changes numerical output at
#: ties, so "which version produced this result" is a question users will need to
#: answer.
__version__ = "0.2.0"


from .params import CBASParams, CBASResult
from .io import (load_subject_data, extract_choice_stream, extract_choice_streams_by_block,
                 enumerate_sequences, enumerate_sequences_block_aware,
                 split_sequence_entry, decode_symbol, decode_sequence)
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
                          build_multicontingency_count_matrix,
                          record_criteria)
from .criterion import criterion_trial, reached_criterion
from .core import subject_criteria, NonIntegerCountWarning

__all__ = [
    "__version__",
    "record_criteria",
    "CBASParams",
    "CBASResult",
    "load_subject_data",
    "extract_choice_stream",
    "extract_choice_streams_by_block",
    "enumerate_sequences",
    "enumerate_sequences_block_aware",
    "split_sequence_entry",
    "decode_symbol",
    "decode_sequence",
    "build_count_matrix",
    "compute_test_stats",
    "compute_test_stats_correlative",
    "NonIntegerCountWarning",
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
