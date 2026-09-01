"""Loading and segmenting data that spans multiple contingencies.

A contingency in the spatial alternation task is the set of three arms that can
pay: a centre arm and two outer arms. Subjects run several contingencies in
sequence, and Kastner's convention is that each is analysed separately, so the
same arm sequence under two contingencies is two hypotheses rather than one.

Two facts about the format make this less obvious than it sounds.

The file identifies a contingency by centre arm and left outer arm, not by an
index. The right outer arm is implied. Critically, **the same arm pair can recur
later in the experiment**, and a recurrence is a separate contingency rather
than a continuation of the earlier one. In the hippocampal lesion cohort every
subject runs six contingency blocks drawn from only five distinct arm pairs,
because the arms 2-3-4 configuration appears twice. Keying on the arm pair
therefore silently merges two contingencies; identity has to come from a block
counter over sessions, which is what `assign_contingency_blocks` provides and
what Igor's `wConting` does.

Rows whose arm or reward field is blank are dropped rather than treated as gaps,
matching Igor, where a blank becomes NaN and `zapNaNs` removes it. Dropping
closes the stream, so a run of rewarded trials can span the removed trial. In
the lesion cohort there are 100 such rows and all of them fall inside
alternation contingencies rather than the exploration phase.

The exploration phase, where the contingency fields are blank and every arm can
pay, is block 0 and is excluded from analysis.

One known divergence from Igor, harmless on the lesion cohort but worth recording.
`assign_contingency_blocks` increments at the exact trial where the arm pair
changes, whereas Igor's `wConting` works session by session and increments at
most once per session, assigning a whole session to the new block even if the
change happened partway through. The two agree whenever a contingency change
coincides with a session boundary. In the hippocampal lesion cohort no session
contains a mid-session change, so the implementations cannot differ there, which
is consistent with the criterion values matching exactly. Data containing
mid-session changes would need this reconciled.
"""

from dataclasses import dataclass

import numpy as np

__all__ = [
    "ContingencyBlock",
    "SubjectRecord",
    "assign_contingency_blocks",
    "load_subject_data_with_contingencies",
    "load_cohort_info",
    "load_cohort_with_contingencies",
    "shared_contingency_blocks",
    "build_multicontingency_count_matrix",
]


@dataclass
class ContingencyBlock:
    """One contiguous run of sessions sharing a contingency.

    block: 0 for exploration, then 1, 2, ... in temporal order
    centre, left_outer: arm identities, or None during exploration
    sessions: session labels in this block, ascending
    """

    block: int
    centre: int | None
    left_outer: int | None
    sessions: np.ndarray

    @property
    def is_exploration(self):
        return self.centre is None

    @property
    def right_outer(self):
        """The implied third active arm, or None during exploration."""
        if self.centre is None:
            return None
        return 2 * self.centre - self.left_outer


@dataclass
class SubjectRecord:
    """One subject's trials, with contingency block membership per trial.

    session, choice, reward: per-trial arrays, blanks already dropped
    block: per-trial contingency block index
    blocks: the ContingencyBlock list, in temporal order
    """

    session: np.ndarray
    choice: np.ndarray
    reward: np.ndarray
    block: np.ndarray
    blocks: list

    def __len__(self):
        return len(self.choice)

    def alternation_blocks(self):
        """Blocks excluding exploration."""
        return [b for b in self.blocks if not b.is_exploration]

    def reward_blocks_for(self, block):
        """Per-session reward arrays for one contingency, in temporal order.

        This is the input `pycbas.criterion` expects. Sessions are kept separate
        so that runs cannot span a session boundary.
        """
        mask = self.block == block
        out = []
        for session in np.unique(self.session[mask]):
            out.append(self.reward[mask & (self.session == session)])
        return out

    def symbol_blocks_for(self, block, num_arms=6, encode_reward=True):
        """Per-session symbol arrays for one contingency, in temporal order."""
        mask = self.block == block
        out = []
        for session in np.unique(self.session[mask]):
            sel = mask & (self.session == session)
            symbols = self.choice[sel]
            if encode_reward:
                symbols = symbols + self.reward[sel] * num_arms
            out.append(symbols)
        return out


def assign_contingency_blocks(session, centre, left_outer,
                              allow_mid_session_change=False):
    """Number contingencies by contiguous block, the way Igor's wConting does.

    The counter increments whenever the (centre, left_outer) pair changes from
    one trial to the next, so a repeat of an earlier arm pair receives a new
    number rather than reusing the old one.

    Args:
        session: per-trial session labels
        centre, left_outer: per-trial arm identities, with -1 marking blank
        allow_mid_session_change: permit a contingency change partway through a
            session, accepting a known divergence from Igor. See below.

    Returns:
        (per-trial block index array, list of ContingencyBlock)

    Raises:
        ValueError: if a session contains more than one contingency and
            `allow_mid_session_change` is False.

    Why the mid-session case raises
    -------------------------------
    This function and Igor's `wConting` agree exactly as long as every
    contingency change coincides with a session boundary, which is true of every
    session in the hippocampal lesion cohort and is why the derived criterion
    values match the Igor reference on every subject-contingency pair.

    They disagree when a change happens partway through a session:

        this function   splits the session at the exact trial, so trials before
                        the change stay with the old contingency
        Igor            works session by session and increments at most once per
                        session, assigning the whole session to the new
                        contingency including trials that ran under the old one

    Neither is obviously wrong. Attributing trials to the contingency actually in
    force is the more literal reading, and Igor's is coarser. But agreement with
    the reference implementation is what makes results comparable, so the
    disagreement has to be resolved deliberately rather than absorbed silently.

    Failing loudly matters because the symptom is invisible: counts would differ
    from Igor's with no error, no warning, and nothing in the output to indicate
    which convention produced them. Pass `allow_mid_session_change=True` only if
    you have decided to accept this function's convention and do not need to
    match Igor on that data.
    """
    session = np.asarray(session)
    centre = np.asarray(centre)
    left_outer = np.asarray(left_outer)

    block_of_trial = np.empty(len(session), dtype=np.int64)
    blocks = []
    current = -1
    previous_key = None

    for i in range(len(session)):
        key = (int(centre[i]), int(left_outer[i]))
        if key != previous_key:
            current += 1
            c = None if key[0] < 0 else key[0]
            lo = None if key[1] < 0 else key[1]
            blocks.append(ContingencyBlock(current, c, lo, []))
            previous_key = key
        block_of_trial[i] = current
        blocks[current].sessions.append(session[i])

    if not allow_mid_session_change:
        # A session spanning more than one block means the contingency changed
        # partway through it, which is the one case where this function and Igor
        # disagree. See the note in the docstring above.
        offenders = [int(s) for s in np.unique(session)
                     if len(np.unique(block_of_trial[session == s])) > 1]
        if offenders:
            shown = ", ".join(str(s) for s in offenders[:5])
            more = f" and {len(offenders) - 5} more" if len(offenders) > 5 else ""
            raise ValueError(
                f"contingency changes partway through session(s) {shown}{more}. "
                "This function splits the session at the exact trial, while Igor's "
                "wConting assigns the whole session to the new contingency, so "
                "results would silently diverge from the reference implementation. "
                "Decide which convention you want, then pass "
                "allow_mid_session_change=True to accept this one."
            )

    for b in blocks:
        b.sessions = np.unique(np.asarray(b.sessions))
    return block_of_trial, blocks


def load_subject_data_with_contingencies(filepath, allow_mid_session_change=False):
    """Load one subject from the multi-contingency text format.

    The format has nine comma-separated columns: session, choice, reward, centre
    arm, left outer arm, two timestamps, and two flags. Only the first five are
    used; the algorithm never reads the rest.

    A one-line header may or may not be present: successive exports of the same
    cohort have differed. It is therefore detected rather than assumed, because
    unconditionally skipping the first line silently discards a real trial when no
    header is there, which shifts every subsequent trial index.

    Rows with a blank choice or reward are dropped, matching Igor.

    Args:
        filepath: path to one subject's file
        allow_mid_session_change: forwarded to `assign_contingency_blocks`, which
            raises by default when a contingency changes partway through a session.

    Returns:
        SubjectRecord
    """
    sessions, choices, rewards, centres, lefts = [], [], [], [], []
    with open(filepath) as fh:
        lines = fh.readlines()

    if lines:
        first_field = lines[0].split(",")[0].strip()
        try:
            int(first_field)
        except ValueError:
            lines = lines[1:]                   # a header, not a trial

    for line in lines:
        parts = line.rstrip("\n").split(",")
        if len(parts) < 5 or not parts[0].strip():
            continue
        if not parts[1].strip() or not parts[2].strip():
            continue                            # blank choice or reward: dropped
        sessions.append(int(parts[0]))
        choices.append(int(parts[1]))
        rewards.append(int(parts[2]))
        centres.append(int(parts[3]) if parts[3].strip() else -1)
        lefts.append(int(parts[4]) if parts[4].strip() else -1)

    session = np.array(sessions, dtype=np.int64)
    choice = np.array(choices, dtype=np.int64)
    reward = np.array(rewards, dtype=np.int64)
    block, blocks = assign_contingency_blocks(
        session, centres, lefts, allow_mid_session_change=allow_mid_session_change)
    return SubjectRecord(session, choice, reward, block, blocks)


def shared_contingency_blocks(records):
    """Validate that a block index means the same thing for every subject.

    Kastner's convention counts each contingency separately, so a column of the
    count matrix is identified by (block, sequence). That is only coherent if
    block 3 refers to the same three arms for every subject, which is what Igor
    checks when it prints "Contingencies not aligned".

    Two failure modes are rejected rather than papered over:

    Misalignment, where subjects disagree about which arms block i uses. Pooling
    those would compare unlike behaviours under one label.

    Incompleteness, where a subject never ran a contingency other subjects did.
    That subject has no data for those columns, and zero is the wrong fill because
    zero means "produced this sequence zero times" rather than "was not measured".
    Igor represents it as NaN and its statistic ignores NaN, giving a per-sequence
    n. pycbas has no NaN support in the statistic or the bootstrap yet, so this
    raises instead of silently writing zeros. Every subject in the hippocampal
    lesion cohort runs all six contingencies, so the case does not arise there.

    Returns:
        sorted list of block indices common to all subjects
    """
    if not records:
        raise ValueError("no subject records given")

    arms_by_block = {}
    block_sets = []
    for index, record in enumerate(records):
        blocks = {b.block: (b.centre, b.left_outer)
                  for b in record.alternation_blocks()}
        block_sets.append(set(blocks))
        for block, arms in blocks.items():
            if block not in arms_by_block:
                arms_by_block[block] = (arms, index)
            elif arms_by_block[block][0] != arms:
                first_arms, first_index = arms_by_block[block]
                raise ValueError(
                    f"contingency block {block} is not aligned across subjects: "
                    f"subject {first_index} has centre/left {first_arms} but "
                    f"subject {index} has {arms}. A block index must denote the "
                    "same arms for every subject before its sequences can be "
                    "pooled into one hypothesis."
                )

    common = set.intersection(*block_sets)
    union = set.union(*block_sets)
    if common != union:
        missing = sorted(union - common)
        offenders = [i for i, s in enumerate(block_sets) if not union <= s]
        raise ValueError(
            f"not every subject ran every contingency: block(s) {missing} are "
            f"absent for {len(offenders)} of {len(records)} subjects, first at "
            f"index {offenders[0]}. Those subjects have no data for the affected "
            "columns, and filling zero would assert they never produced those "
            "sequences. Representing it honestly needs NaN support in the "
            "statistic and the bootstrap, which pycbas does not have yet. Restrict "
            "to the blocks all subjects share by passing `blocks=`."
        )
    return sorted(common)


def build_multicontingency_count_matrix(records, params, blocks=None,
                                        encode_reward=True):
    """Count sequences separately per contingency and concatenate the columns.

    Each contingency is analysed as its own set of hypotheses, so the same arm
    sequence under two contingencies gives two columns. With 100 sequences in the
    first contingency and 200 in the second, the matrix has 300 columns and the
    multiplicity correction runs over all of them jointly.

    Counting is always session-respecting within a contingency, matching Igor,
    which enumerates per session and concatenates. The criterion is applied per
    subject *and* per contingency, since each contingency is its own learning
    episode, and its trial index is local to that contingency.

    Args:
        records: list of SubjectRecord
        params: CBASParams; `criterion_order` and `criterion` apply within each
            contingency
        blocks: contingency block indices to include, default all that every
            subject shares. Exploration is never included.
        encode_reward: encode reward into symbols

    Returns:
        (sequences, count_matrix) where `sequences` is a list of
        (block, sequence_tuple) pairs and `count_matrix` has shape
        (n_subjects, len(sequences)).
    """
    from .criterion import criterion_trial, as_enumeration_cutoff
    from .io import enumerate_sequences_block_aware

    available = shared_contingency_blocks(records)
    if blocks is None:
        blocks = available
    else:
        blocks = sorted(blocks)
        unknown = [b for b in blocks if b not in available]
        if unknown:
            raise ValueError(
                f"contingency block(s) {unknown} are not shared by all subjects; "
                f"available: {available}")

    order = getattr(params, "criterion_order", 0)

    # counts[subject][(block, sequence)] = n
    per_subject = []
    for record in records:
        counts = {}
        for block in blocks:
            streams = record.symbol_blocks_for(block, params.num_arms, encode_reward)
            n_trials = sum(len(s) for s in streams)
            cutoff = as_enumeration_cutoff(
                criterion_trial(record.reward_blocks_for(block), order,
                                params.criterion),
                n_trials)
            for seq_len in range(1, params.seq_len_max + 1):
                for seq, n in enumerate_sequences_block_aware(
                        streams, seq_len, cutoff).items():
                    counts[(block, seq)] = n
        per_subject.append(counts)

    totals = {}
    for counts in per_subject:
        for key, n in counts.items():
            totals[key] = totals.get(key, 0) + n

    # Ordered by contingency first, then as build_count_matrix orders within one:
    # descending total frequency, then length, then value. Keeping contingencies
    # contiguous makes the matrix readable and lets callers slice one out.
    sequences = sorted(totals, key=lambda k: (k[0], -totals[k], len(k[1]), k[1]))

    index = {key: i for i, key in enumerate(sequences)}
    count_matrix = np.zeros((len(records), len(sequences)), dtype=np.float64)
    for row, counts in enumerate(per_subject):
        for key, n in counts.items():
            count_matrix[row, index[key]] = n
    return sequences, count_matrix


def load_cohort_with_contingencies(directory, allow_mid_session_change=False):
    """Load every numbered subject file in a directory, plus the info table.

    Files are ordered numerically by the digits in their stem, so `an2` precedes
    `an10`, which is the order the info table rows correspond to.

    Returns:
        (records, info) where records is a list of SubjectRecord and info is the
        list of dicts from `load_cohort_info`, or None when no info file exists.
    """
    from pathlib import Path
    directory = Path(directory)
    files = [p for p in directory.glob("an*.txt") if p.stem.lower() != "aninfo"]
    if not files:
        raise ValueError(f"no subject files matching an*.txt in {directory}")
    files.sort(key=lambda p: int("".join(c for c in p.stem if c.isdigit())))

    records = [load_subject_data_with_contingencies(
        p, allow_mid_session_change=allow_mid_session_change) for p in files]

    info_path = directory / "anInfo.txt"
    info = load_cohort_info(info_path) if info_path.exists() else None
    if info is not None and len(info) != len(records):
        raise ValueError(
            f"{len(records)} subject files but {len(info)} rows in "
            f"{info_path.name}; they must correspond one to one")
    return records, info


def load_cohort_info(filepath):
    """Load the cohort info table accompanying the multi-contingency format.

    Two header lines, then one row per subject in the same order as the numbered
    subject files. Values are strings rather than the numeric codes used by the
    published format. A blank lesion field means the subject had lesion surgery
    but no lesion was evident on CT, which is neither control nor lesion.

    Returns:
        list of dicts keyed by the column names on the second header line.
    """
    with open(filepath) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    columns = [c.strip() for c in lines[1].split(",")]
    rows = []
    for line in lines[2:]:
        values = [v.strip() for v in line.split(",")]
        rows.append({c: (v if v else None) for c, v in zip(columns, values)})
    return rows
