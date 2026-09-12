"""Data loading and sequence enumeration."""

import numpy as np


def load_subject_data(filepath):
    """Load a single subject's data file. Returns one (n_trials, 4) int32 array of session, choice, reward, contingency."""
    rows = []
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split(",")
            session = int(parts[0])
            choice = int(parts[1])
            reward = int(parts[2])
            conting = int(parts[3]) if parts[3].strip() else 0
            rows.append((session, choice, reward, conting))
    arr = np.array(rows, dtype=np.int32)
    return arr


def extract_choice_stream(subject_data, contingency=2, num_arms=6, encode_reward=True):
    """Extract choice stream, optionally filtered by contingency.

    Args:
        subject_data: array with columns (session, symbol, reward, contingency)
        contingency: block type to filter on, or None to use all trials
        num_arms: number of base symbols (choices)
        encode_reward: if True, encode as symbol + reward*num_arms (doubles alphabet).
            Set False for tasks where outcome is deterministic from the symbol
            (e.g., 2AFC where symbol already encodes choice x stimulus_side).
    """
    if contingency is None:
        data = subject_data
    else:
        mask = subject_data[:, 3] == contingency
        data = subject_data[mask]
    if encode_reward:
        symbols = data[:, 1] + data[:, 2] * num_arms
    else:
        symbols = data[:, 1]
    return symbols


def split_sequence_entry(entry):
    """`(block, symbols)` for one entry of a `sequences` list.

    The multi-contingency pipeline keys columns by `(block, sequence)`, because the
    same arm sequence under two contingencies is two hypotheses; every other pipeline
    keys by the sequence alone. Callers that accept results from either should not
    have to know which produced them, and in particular must not take `len(entry)` as
    the sequence length: for a multi-contingency entry that is 2 for every hypothesis
    whatever its actual length.

    Returns:
        (block, symbols), where block is None for a bare sequence.
    """
    if len(entry) == 2 and isinstance(entry[1], tuple):
        return entry[0], entry[1]
    return None, tuple(entry)


def decode_symbol(sym, num_arms=6, encode_reward=True):
    """Invert the symbol encoding of `extract_choice_stream`.

    Args:
        sym: one encoded symbol.
        num_arms: number of base choices, the same value used to encode.
        encode_reward: whether reward was folded into the symbol.

    Returns:
        (choice, rewarded), where choice is 0-based and rewarded is a bool, or None
        when `encode_reward` is False, since the outcome is then not recoverable from
        the symbol rather than being known to be absent.

    Raises:
        ValueError: if `sym` cannot have been produced by this `num_arms`. Decoding
            with the wrong `num_arms` otherwise succeeds silently and reports the
            wrong arm, which is indistinguishable from a real result.
    """
    sym = int(sym)
    if sym < 0:
        raise ValueError(f"symbol {sym} is negative")
    if not encode_reward:
        return sym, None
    if sym >= 2 * num_arms:
        raise ValueError(
            f"symbol {sym} is out of range for num_arms={num_arms} with reward "
            f"encoding, which allows 0 to {2 * num_arms - 1}; check num_arms")
    return sym % num_arms, bool(sym // num_arms)


def decode_sequence(entry, num_arms=6, encode_reward=True, join=" "):
    """Readable label for one entry of a `sequences` list.

    Uses the published convention: arms are numbered from 1, and a trailing `*` marks
    a rewarded choice, so with `num_arms=6` the symbol 2 reads as `3` and the symbol 8
    reads as `3*`. A multi-contingency entry is prefixed with its contingency block.

    Args:
        entry: a sequence tuple, or a `(block, sequence)` pair.
        num_arms: number of base choices, the same value used to encode.
        encode_reward: whether reward was folded into the symbol.
        join: separator between symbols.
    """
    block, symbols = split_sequence_entry(entry)
    parts = []
    for sym in symbols:
        choice, rewarded = decode_symbol(sym, num_arms, encode_reward)
        parts.append(f"{choice + 1}{'*' if rewarded else ''}")
    body = join.join(parts)
    return f"c{block}: {body}" if block is not None else body


def extract_choice_streams_by_block(subject_data, contingency=2, num_arms=6, encode_reward=True):
    """Extract choice streams split by block/session boundaries.

    Returns a list of arrays, one per block, preserving temporal order.
    Sequences should not span across blocks.
    """
    if contingency is None:
        data = subject_data
    else:
        mask = subject_data[:, 3] == contingency
        data = subject_data[mask]

    if encode_reward:
        symbols = data[:, 1] + data[:, 2] * num_arms
    else:
        symbols = data[:, 1]

    blocks = data[:, 0]
    streams = []
    for block in np.unique(blocks):
        block_mask = blocks == block
        streams.append(symbols[block_mask])
    return streams


def enumerate_sequences(choice_stream, seq_len, criterion):
    """Find all subsequences of given length with start position <= criterion.

    Matches Igor's counting: sequences starting at positions 0..criterion
    (inclusive) are counted, using elements up to position criterion+seq_len-1.
    """
    max_start = min(criterion, len(choice_stream) - seq_len)
    counts = {}
    for i in range(max_start + 1):
        seq = tuple(choice_stream[i:i + seq_len].tolist())
        counts[seq] = counts.get(seq, 0) + 1
    return counts


def enumerate_sequences_block_aware(block_streams, seq_len, criterion):
    """Count sequences in a block-aware stream, matching Igor's approach.

    The criterion is applied to global position in the concatenated stream
    (start_position <= criterion, inclusive). A sequence is only counted if all
    its positions fall within the same block.

    Args:
        block_streams: list of symbol arrays, one per block
        seq_len: length of sequences to enumerate
        criterion: maximum global start position (inclusive)
    """
    counts = {}
    global_pos = 0
    for stream in block_streams:
        block_len = len(stream)
        for i in range(block_len - seq_len + 1):
            if global_pos + i > criterion:
                return counts
            seq = tuple(stream[i:i + seq_len].tolist())
            counts[seq] = counts.get(seq, 0) + 1
        global_pos += block_len
        if global_pos > criterion:
            break
    return counts
