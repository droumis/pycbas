"""Subject identity, cohort assembly, and the labelled count matrix.

Everything CBAS computes is per subject, and until these types existed the
correspondence between a subject's trials, its identity, and its group was
positional: separate lists that callers had to keep in step by hand. That is the one
mistake nothing downstream can detect, because a permuted list of the right length is
indistinguishable from a correct one, and per-subject count totals are identical
whenever a fixed criterion truncates every subject to the same number of windows.

So identity travels with the trials, group membership is resolved through that
identity rather than by position, and the count matrix carries the ids of its own
rows. A reordering anywhere then relabels nothing.

Group membership deliberately does **not** live on `Subject`. The same animal is
control in one comparison and, say, male in another, so a label belongs to an
analysis, not to the data. What lives on the subject is `meta`, the columns its
cohort file supplied, from which any analysis can derive its own grouping.
"""

from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

#: Column order of `Subject.trials`, matching the on-disk single-contingency format.
SESSION, CHOICE, REWARD, CONDITION = range(4)


@dataclass
class Subject:
    """One subject's trials, with the identity they were loaded under.

    `trials` is an (n_trials, 4) array of session, choice, reward, condition, the
    layout the text loaders produce. `condition` is whatever the file's fourth column
    held for single-contingency data, and the contingency block index for
    multi-contingency data, so `contingency=2` selects the same thing in both.

    `blocks` is present only for multi-contingency subjects, where it gives the arm
    identity and exploration flag of each condition value. Its absence is what
    `build_multicontingency_count_matrix` refuses on, rather than silently treating a
    raw condition column as a block structure it validated.
    """

    id: str
    trials: np.ndarray
    meta: dict = field(default_factory=dict)
    blocks: list | None = None

    def __post_init__(self):
        trials = np.asarray(self.trials)
        if trials.ndim != 2 or trials.shape[1] != 4:
            raise ValueError(
                f"subject {self.id!r}: trials must have shape (n_trials, 4) of "
                f"session, choice, reward, condition, got {trials.shape}")
        self.trials = trials

    def __len__(self):
        return len(self.trials)

    @property
    def session(self):
        return self.trials[:, SESSION]

    @property
    def choice(self):
        return self.trials[:, CHOICE]

    @property
    def reward(self):
        return self.trials[:, REWARD]

    @property
    def condition(self):
        return self.trials[:, CONDITION]

    @property
    def has_blocks(self):
        return self.blocks is not None

    def conditions(self):
        """Sorted condition values this subject has trials for."""
        return sorted(int(c) for c in np.unique(self.condition))

    def alternation_blocks(self):
        """Blocks excluding exploration. Multi-contingency subjects only."""
        self._require_blocks()
        return [b for b in self.blocks if not b.is_exploration]

    def reward_blocks_for(self, condition):
        """Per-session reward arrays for one condition, in temporal order.

        This is the input `pycbas.criterion` expects. Sessions are kept separate so
        that runs cannot span a session boundary.
        """
        return self._per_session(condition, self.reward)

    def symbol_blocks_for(self, condition, num_arms=6, encode_reward=True):
        """Per-session symbol arrays for one condition, in temporal order."""
        mask = self.condition == condition
        out = []
        for session in np.unique(self.session[mask]):
            sel = mask & (self.session == session)
            symbols = self.choice[sel]
            if encode_reward:
                symbols = symbols + self.reward[sel] * num_arms
            out.append(symbols)
        return out

    def _per_session(self, condition, values):
        mask = self.condition == condition
        return [values[mask & (self.session == s)]
                for s in np.unique(self.session[mask])]

    def _require_blocks(self):
        if self.blocks is None:
            raise ValueError(
                f"subject {self.id!r} has no contingency block structure; it was "
                f"loaded from the single-contingency format, where the fourth column "
                f"is a raw condition code rather than a validated block index")


class Cohort:
    """An ordered collection of `Subject`, addressable by id.

    Ordered because the count matrix's rows follow it, and addressable because
    everything aligned to those rows is resolved through ids rather than positions.
    Ids must be unique: two subjects answering to one name would make that resolution
    ambiguous, and the loaders derive ids from filenames, so a collision means two
    files of the same name from different directories. Pass explicit ids for that.
    """

    def __init__(self, subjects):
        subjects = list(subjects)
        if not subjects:
            raise ValueError("a cohort needs at least one subject")
        for s in subjects:
            if not isinstance(s, Subject):
                raise TypeError(
                    f"a Cohort holds Subject instances, got {type(s).__name__}; "
                    f"load with load_subject or load_cohort")
        ids = [s.id for s in subjects]
        duplicates = sorted({i for i in ids if ids.count(i) > 1})
        if duplicates:
            raise ValueError(
                f"subject ids must be unique, but {duplicates} appear more than "
                f"once; pass explicit ids to distinguish them")
        self.subjects = subjects

    def __len__(self):
        return len(self.subjects)

    def __iter__(self):
        return iter(self.subjects)

    def __getitem__(self, key):
        """Position, id, or slice. A slice gives another cohort."""
        if isinstance(key, slice):
            return Cohort(self.subjects[key])
        if isinstance(key, str):
            for s in self.subjects:
                if s.id == key:
                    return s
            raise KeyError(f"no subject with id {key!r}")
        return self.subjects[key]

    def __repr__(self):
        return f"Cohort({len(self)} subjects: {', '.join(self.ids[:4])}" + \
               (", ...)" if len(self) > 4 else ")")

    @property
    def ids(self):
        return [s.id for s in self.subjects]

    @property
    def has_blocks(self):
        return all(s.has_blocks for s in self.subjects)

    def reorder(self, order):
        """A cohort in the given order of positions. Identity follows each subject."""
        return Cohort([self.subjects[i] for i in order])

    def filter(self, predicate=None, **meta_equals):
        """Subjects passing a predicate and matching every `meta` value given.

        `cohort.filter(genotype="WT")` keeps the wild types; a value may also be a
        set or list to keep several. Filtering cannot misalign anything, which is the
        point of doing it here rather than over parallel lists.
        """
        kept = []
        for s in self.subjects:
            if predicate is not None and not predicate(s):
                continue
            if any(not _matches(s.meta.get(k), v) for k, v in meta_equals.items()):
                continue
            kept.append(s)
        return Cohort(kept)

    def meta_values(self, key):
        """The `meta` value of `key` for every subject, in cohort order."""
        return [s.meta.get(key) for s in self.subjects]

    def labels_from(self, key, coder=None):
        """Group labels from a `meta` column, as an array aligned to this cohort.

        `coder` maps a raw value to 0, 1, or None to exclude; without one the values
        must already be 0 and 1. Subjects the coder excludes are reported rather than
        dropped, since dropping them here would silently change the cohort a caller
        thought it had assembled.
        """
        coder = coder or _default_coder
        labels, unusable = [], []
        for s in self.subjects:
            value = coder(s.meta.get(key))
            if value is None:
                unusable.append(s.id)
            labels.append(value)
        if unusable:
            raise ValueError(
                f"{len(unusable)} of {len(self)} subjects have no usable {key!r} "
                f"value: {unusable[:6]}{'...' if len(unusable) > 6 else ''}. Drop "
                f"them with Cohort.filter, or pass a coder that handles them.")
        return np.asarray(labels, dtype=np.int64)

    def covariate_from(self, key):
        """A float covariate from a `meta` column, aligned to this cohort."""
        values = []
        for s in self.subjects:
            raw = s.meta.get(key)
            try:
                values.append(float(raw))
            except (TypeError, ValueError):
                raise ValueError(
                    f"subject {s.id!r} has {key!r} = {raw!r}, which is not a "
                    f"number") from None
        return np.asarray(values, dtype=np.float64)


def _matches(value, wanted):
    if isinstance(wanted, (set, frozenset, list, tuple)):
        return value in wanted
    return value == wanted


def _default_coder(value):
    if value is None:
        return None
    if isinstance(value, bool):
        return int(value)
    try:
        number = int(value)
    except (TypeError, ValueError):
        return None
    return number if number in (0, 1) else None


def resolve_labels(labels, ids, cohort_ids=None):
    """Group labels for `ids`, from a mapping, a cohort order, or a `meta` key.

    Accepts a `{id: label}` mapping, a sequence in cohort order, or the name of a
    `meta` column when `labels` is a string and a cohort is given. A sequence is
    converted to a mapping first, so the result follows `ids` even if those are not
    in cohort order. That is what keeps a reordering harmless: the labels are carried
    by identity, never by position.
    """
    if isinstance(labels, dict):
        missing = [i for i in ids if i not in labels]
        if missing:
            raise ValueError(
                f"no label for {len(missing)} of {len(ids)} subjects: "
                f"{missing[:6]}{'...' if len(missing) > 6 else ''}")
        return np.asarray([labels[i] for i in ids])

    values = np.asarray(labels)
    if values.ndim != 1:
        raise ValueError(f"labels must be one-dimensional, got shape {values.shape}")
    order = list(cohort_ids) if cohort_ids is not None else list(ids)
    if len(values) != len(order):
        raise ValueError(
            f"{len(values)} labels for {len(order)} subjects; they must correspond "
            f"one to one, in the same order")
    return resolve_labels(dict(zip(order, values.tolist())), ids)


@dataclass
class CountMatrix:
    """Sequence counts with the ids of its rows and the labels of its columns.

    `counts[i, j]` is how often subject `subject_ids[i]` used sequence
    `sequences[j]`. Both labels are carried rather than implied, so a caller can
    check the correspondence instead of trusting an argument order, and anything
    aligned to the rows can be resolved by id.
    """

    counts: np.ndarray
    subject_ids: list
    sequences: list

    def __post_init__(self):
        n_rows, n_cols = self.counts.shape
        if len(self.subject_ids) != n_rows:
            raise ValueError(
                f"{len(self.subject_ids)} subject ids for {n_rows} rows")
        if len(self.sequences) != n_cols:
            raise ValueError(
                f"{len(self.sequences)} sequences for {n_cols} columns")

    @property
    def shape(self):
        return self.counts.shape

    def row(self, subject_id):
        """One subject's counts, by id."""
        try:
            return self.counts[self.subject_ids.index(subject_id)]
        except ValueError:
            raise KeyError(f"no row for subject {subject_id!r}") from None

    def reorder(self, order):
        """Rows in the given order, ids following them so nothing is relabelled."""
        return replace(self, counts=self.counts[list(order)],
                       subject_ids=[self.subject_ids[i] for i in order])

    def select_block(self, block):
        """The columns of one contingency block.

        Raises on a single-contingency matrix rather than returning nothing. An empty
        selection there looks like a contingency with no sequences, which is the kind
        of quiet wrong answer this type exists to prevent.
        """
        from .io import split_sequence_entry
        if not self.sequences or split_sequence_entry(self.sequences[0])[0] is None:
            raise ValueError(
                "this count matrix has bare sequence entries, so it has no "
                "contingency blocks to select; it came from build_count_matrix "
                "rather than build_multicontingency_count_matrix")
        keep = [j for j, entry in enumerate(self.sequences)
                if split_sequence_entry(entry)[0] == block]
        if not keep:
            raise ValueError(f"no columns for block {block!r}")
        return replace(self, counts=self.counts[:, keep],
                       sequences=[self.sequences[j] for j in keep])

    def column_labels(self, num_arms=6, encode_reward=True, join=" "):
        """Readable labels for the columns, in the published convention."""
        from .io import decode_sequence
        return [decode_sequence(entry, num_arms, encode_reward, join)
                for entry in self.sequences]


def load_subject(filepath, id=None, meta=None):
    """One `Subject` from the single-contingency text format.

    The id defaults to the filename stem, which is the only identity the format
    carries. Pass one explicitly when the stem is not the name you want to see in
    results, or when two cohorts use the same filenames.
    """
    from .io import load_subject_data
    path = Path(filepath)
    return Subject(id=id if id is not None else path.stem,
                   trials=load_subject_data(path),
                   meta=dict(meta or {}))


def load_cohort(source, pattern="*.txt", meta=None, ids=None):
    """A `Cohort` from a directory or an explicit list of files.

    Files from a directory are ordered by the digits in their names, so `an2` precedes
    `an10`, and lexicographically when a name has no digits. An explicit list is used
    in the order given, since a caller who listed the files has already chosen one.

    Args:
        source: a directory, or an iterable of paths.
        pattern: glob used when `source` is a directory.
        meta: optional `{id: dict}` of per-subject metadata, or a list in file order.
        ids: optional explicit ids, in the same order as the files.
    """
    if isinstance(source, (str, Path)) and Path(source).is_dir():
        paths = sorted(Path(source).glob(pattern), key=_numeric_key)
        if not paths:
            raise ValueError(f"no files matching {pattern!r} in {source}")
    else:
        paths = [Path(p) for p in source]
        if not paths:
            raise ValueError("no subject files given")

    if ids is not None and len(ids) != len(paths):
        raise ValueError(f"{len(ids)} ids for {len(paths)} files")

    subjects = []
    for i, path in enumerate(paths):
        subject_id = ids[i] if ids is not None else path.stem
        if isinstance(meta, dict):
            subject_meta = meta.get(subject_id, {})
        elif meta is not None:
            if len(meta) != len(paths):
                raise ValueError(f"{len(meta)} meta rows for {len(paths)} files")
            subject_meta = meta[i]
        else:
            subject_meta = {}
        subjects.append(load_subject(path, id=subject_id, meta=subject_meta))
    return Cohort(subjects)


def _numeric_key(path):
    """Sort key that orders `an2` before `an10`, and names without digits by name."""
    digits = "".join(c for c in path.stem if c.isdigit())
    return (0, int(digits), path.stem) if digits else (1, 0, path.stem)
