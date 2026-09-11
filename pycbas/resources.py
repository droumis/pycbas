"""Resource estimation utilities."""


def estimate_resources(num_arms, seq_len_max, n_subjects=None, n_observed=None,
                       resample_number=10000, encode_reward=True,
                       n_hypothesis_sets=1):
    """Estimate memory and time requirements for a CBAS run.

    Args:
        num_arms: number of base symbols (before reward encoding)
        seq_len_max: maximum sequence length L
        n_subjects: number of subjects (for context only)
        n_observed: actual number of observed sequences (if known from data).
            Overrides worst-case S for memory/time estimates.
        resample_number: number of bootstrap resamples M (default 10,000)
        encode_reward: if True, effective alphabet is num_arms*2
        n_hypothesis_sets: how many times the sequence space is counted. One for an
            ordinary run; for a multi-contingency run it is the number of contingency
            blocks, because each block contributes its own column for the same
            sequence. Leaving it at one made the reported ceiling smaller than the
            observed count, so the estimate contradicted itself.

    Returns:
        dict with resource estimates and recommendation
    """
    A = num_arms * 2 if encode_reward else num_arms
    S = n_hypothesis_sets * sum(A**l for l in range(1, seq_len_max + 1))
    M = resample_number

    # Both the null matrix and the submatrix are sized by the *count matrix's*
    # columns, which is the number of sequences actually observed, not by the
    # enumerable sequence space S. S is only the fallback worst case for a caller
    # who has not counted yet. Using S here regardless made the unchunked estimate
    # wrong by the ratio S / n_observed, which is two orders of magnitude on a
    # sparse space, and contradicted the documented halving.
    n_cols = n_observed if n_observed is not None else S
    # One column per one-sided hypothesis, which is one per sequence rather than
    # two, because only the observed direction has a defined statistic.
    n_valid = n_cols

    # The step-down holds one entry per (resample, hypothesis): a float64 magnitude
    # and an int8 direction, so nine bytes, not sixteen.
    PER_ENTRY = 8 + 1
    # Chunked: only the sorted submatrix is allocated.
    chunked_bytes = M * n_valid * PER_ENTRY
    # Unchunked: the full-width null exists at the same time as the extracted
    # submatrix, so peak holds both. Hence roughly double the chunked figure, which
    # is what the guide and docs/algorithm.md describe.
    full_null_bytes = M * (n_cols + n_valid) * PER_ENTRY

    full_null_gb = full_null_bytes / (1024**3)
    chunked_gb = chunked_bytes / (1024**3)

    # Two terms, fitted to the two published validation runs at M=10,000: the human
    # analysis, 408 hypotheses in about 3 seconds, and the rat analysis, 16,376 in
    # about 12. A single point through the origin misses the fixed cost that
    # dominates a small run -- compiling the kernels, building the count matrix --
    # and reported "~0s" for a run the app then took several seconds to finish,
    # which is a contradiction a user can see in one screen.
    #
    # Only the variable term scales with M: the bootstrap fills M x n_valid entries
    # and the step-down scans them, while the fixed cost is paid once.
    REF_RESAMPLES = 10000
    FIXED_SECONDS = 2.8
    SECONDS_PER_HYPOTHESIS = (12.0 - FIXED_SECONDS) / 16376
    est_time = (FIXED_SECONDS
                + SECONDS_PER_HYPOTHESIS * n_valid * (M / REF_RESAMPLES))

    if chunked_gb < 1.0:
        verdict = "TRIVIAL"
    elif chunked_gb < 8.0:
        verdict = "COMFORTABLE"
    elif chunked_gb < 24.0:
        verdict = "FITS (close other apps)"
    elif chunked_gb < 48.0:
        verdict = "TIGHT (may need reduced M or chunking to disk)"
    else:
        verdict = "TOO LARGE (needs embedding or M reduction)"

    return {
        "alphabet": A,
        "seq_len_max": seq_len_max,
        "total_sequences": S,
        "observed_sequences": n_observed,
        "resample_number": M,
        "n_subjects": n_subjects,
        "memory_full_null_gb": round(full_null_gb, 2),
        "memory_chunked_gb": round(chunked_gb, 2),
        "est_time_seconds": round(est_time, 1),
        "recommendation": verdict,
    }


def print_resource_estimate(est):
    """Pretty-print output from estimate_resources()."""
    print("CBAS Resource Estimate")
    print(f"  Alphabet (A):          {est['alphabet']}")
    print(f"  Max length (L):        {est['seq_len_max']}")
    print(f"  Total sequences (S):   {est['total_sequences']:,}")
    if est["observed_sequences"] is not None:
        pct = est["observed_sequences"] / est["total_sequences"] * 100
        print(f"  Observed sequences:    {est['observed_sequences']:,} ({pct:.0f}%)")
    if est["n_subjects"] is not None:
        print(f"  Subjects:              {est['n_subjects']}")
    print(f"  Resamples (M):         {est['resample_number']:,}")
    print(f"  Memory (standard):     {est['memory_full_null_gb']:.1f} GB")
    print(f"  Memory (chunked):      {est['memory_chunked_gb']:.1f} GB")
    print(f"  Est. time (chunked):   {est['est_time_seconds']:.0f}s")
    print(f"  Verdict:               {est['recommendation']}")
