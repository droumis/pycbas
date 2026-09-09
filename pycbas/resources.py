"""Resource estimation utilities."""


def estimate_resources(num_arms, seq_len_max, n_subjects=None, n_observed=None,
                       resample_number=10000, encode_reward=True):
    """Estimate memory and time requirements for a CBAS run.

    Args:
        num_arms: number of base symbols (before reward encoding)
        seq_len_max: maximum sequence length L
        n_subjects: number of subjects (for context only)
        n_observed: actual number of observed sequences (if known from data).
            Overrides worst-case S for memory/time estimates.
        resample_number: number of bootstrap resamples M (default 10,000)
        encode_reward: if True, effective alphabet is num_arms*2

    Returns:
        dict with resource estimates and recommendation
    """
    A = num_arms * 2 if encode_reward else num_arms
    S = sum(A**l for l in range(1, seq_len_max + 1))
    M = resample_number

    n_valid = n_observed if n_observed is not None else S

    # The step-down holds one entry per (resample, hypothesis): a float64 magnitude
    # and an int8 direction, so nine bytes, not sixteen. And there is one column per
    # one-sided hypothesis, which is one per sequence rather than two, because only
    # the observed direction has a defined statistic. An earlier version doubled the
    # columns and omitted the directions, overestimating by about 1.8x.
    PER_ENTRY = 8 + 1
    # Chunked: only the sorted submatrix is allocated.
    chunked_bytes = M * n_valid * PER_ENTRY
    # Unchunked: the full-width null exists at the same time as the extracted
    # submatrix, so peak holds both.
    full_null_bytes = M * (S + n_valid) * PER_ENTRY

    full_null_gb = full_null_bytes / (1024**3)
    chunked_gb = chunked_bytes / (1024**3)

    # Calibrated on the published rat validation: 16,376 hypotheses with a defined
    # statistic, timed end to end. Linear in the hypothesis count, which is the
    # dominant term for both the bootstrap and the step-down.
    rat_cols = 16376
    rat_time = 12.3
    est_time = rat_time * n_valid / rat_cols

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
    print(f"CBAS Resource Estimate")
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
