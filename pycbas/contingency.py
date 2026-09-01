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

    The format carries a one-line header and nine comma-separated columns:
    session, choice, reward, centre arm, left outer arm, two timestamps, and two
    flags. Only the first five are used; the algorithm never reads the rest.

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

    for line in lines[1:]:                      # skip the header line
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
