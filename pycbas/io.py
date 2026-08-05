"""Data loading and sequence enumeration."""

import numpy as np


def load_subject_data(filepath):
    """Load a single subject's data file. Returns (session, choice, reward, contingency) arrays."""
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
