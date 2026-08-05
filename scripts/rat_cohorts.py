"""Cohort selection for the rat spatial alternation dataset.

The rat data lives in `data/rats/` as `an{N}.txt` plus `anInfo.txt`. Row N of
`anInfo.txt` (0-indexed, after the header) describes `an{N}.txt`. This mapping
was established by matching stream contents against the 111 files released in
`igor_cbas/data/`, all of which match byte-exactly.

Column semantics, recovered from that match:

    name        animal identifier
    experiment  0 = male batch, 1 = female batch, 2 and 3 = later male batches
    sex         0 = male, 1 = female
    genotype    0 = the 111 animals released to the Igor repository
    lesion      0 = control, 1 = lesion, blank = lesion animal excluded

The published cohorts follow from these:

    initial      experiment in {0, 1}, genotype 0, lesion known -> 46 control, 39 lesion
    replication  experiment in {2, 3}, genotype 0, lesion known ->  9 control, 11 lesion

`genotype == 1` is a further 111 animals that were never released. Their status
is unresolved, so they are available here but not part of any published result.

Do not select rats by filename sort order. Doing so mixes the initial and
replication cohorts and silently includes the 6 excluded animals.
"""

from pathlib import Path

import numpy as np

from pycbas import load_subject_data

ROOT_DIR = Path(__file__).parent.parent
DEFAULT_DATA_DIR = ROOT_DIR / "data" / "rats"

# Published cohorts plus the unreleased genotype-1 set.
COHORTS = {
    "initial": {"experiment": (0, 1), "genotype": 0},
    "replication": {"experiment": (2, 3), "genotype": 0},
    "genotype1": {"experiment": (0, 1), "genotype": 1},
    "genotype1_late": {"experiment": (2, 3), "genotype": 1},
    "all_published": {"experiment": (0, 1, 2, 3), "genotype": 0},
}

# Expected control/lesion counts, asserted on load so that a silent change in
# the data directory cannot go unnoticed.
EXPECTED_COUNTS = {
    "initial": (46, 39),
    "replication": (9, 11),
    "all_published": (55, 50),
}


def read_info(data_dir=None):
    """Read anInfo.txt into a list of dicts, one per animal, in file order.

    Returns dicts with keys: index, name, experiment, sex, genotype, lesion,
    path. `lesion` is None where the field is blank.
    """
    data_dir = Path(data_dir) if data_dir is not None else DEFAULT_DATA_DIR
    info_path = data_dir / "anInfo.txt"
    if not info_path.exists():
        raise FileNotFoundError(f"missing {info_path}")

    lines = [ln for ln in info_path.read_text().splitlines() if ln.strip()]
    header = [h.strip() for h in lines[0].split(",")]
    expected = ["name", "experiment", "sex", "genotype", "lesion"]
    if header != expected:
        raise ValueError(f"unexpected anInfo.txt header {header}, expected {expected}")

    records = []
    for i, line in enumerate(lines[1:]):
        parts = [p.strip() for p in line.split(",")]
        path = data_dir / f"an{i}.txt"
        records.append({
            "index": i,
            "name": int(parts[0]),
            "experiment": int(parts[1]),
            "sex": int(parts[2]),
            "genotype": int(parts[3]),
            "lesion": _parse_lesion(parts[4]),
            "path": path,
        })

    missing = [r["path"].name for r in records if not r["path"].exists()]
    if missing:
        raise FileNotFoundError(
            f"anInfo.txt has {len(records)} rows but these files are missing: "
            f"{missing[:5]}{' ...' if len(missing) > 5 else ''}")

    stray = sorted(p.name for p in data_dir.glob("an*.txt")
                   if p.name != "anInfo.txt" and not _is_indexed(p, len(records)))
    if stray:
        raise ValueError(
            f"{len(stray)} an*.txt files have no anInfo.txt row: {stray[:5]}"
            f"{' ...' if len(stray) > 5 else ''}")

    return records


def _parse_lesion(token):
    """Parse the lesion field. Blank or 'NaN' marks an excluded lesion animal."""
    if token == "" or token.lower() == "nan":
        return None
    return int(token)


def _is_indexed(path, n_records):
    stem = path.stem
    if not stem.startswith("an") or not stem[2:].isdigit():
        return False
    return int(stem[2:]) < n_records


def select_cohort(cohort="initial", data_dir=None, require_lesion_known=True,
                  sex=None):
    """Return the anInfo records belonging to a cohort, without loading streams.

    Args:
        cohort: key of COHORTS, or a dict of column -> allowed value(s)
        data_dir: override the default data/rats directory
        require_lesion_known: drop animals whose lesion field is blank. These
            are lesion animals excluded from the published analysis, so keeping
            them changes the reported group sizes.
        sex: optionally restrict to 0 (male) or 1 (female)
    """
    spec = COHORTS[cohort] if isinstance(cohort, str) else dict(cohort)
    records = read_info(data_dir)

    def allowed(value, spec_value):
        if isinstance(spec_value, (tuple, list, set)):
            return value in spec_value
        return value == spec_value

    selected = [
        r for r in records
        if all(allowed(r[col], val) for col, val in spec.items())
    ]
    if require_lesion_known:
        selected = [r for r in selected if r["lesion"] is not None]
    if sex is not None:
        selected = [r for r in selected if r["sex"] == sex]

    if isinstance(cohort, str) and cohort in EXPECTED_COUNTS \
            and require_lesion_known and sex is None:
        n_ctrl = sum(1 for r in selected if r["lesion"] == 0)
        n_les = sum(1 for r in selected if r["lesion"] == 1)
        if (n_ctrl, n_les) != EXPECTED_COUNTS[cohort]:
            raise ValueError(
                f"cohort '{cohort}' resolved to {n_ctrl} control / {n_les} lesion, "
                f"expected {EXPECTED_COUNTS[cohort]}")

    return selected


def load_rat_cohort(cohort="initial", data_dir=None, require_lesion_known=True,
                    sex=None, groups=("control", "lesion")):
    """Load streams for a cohort.

    Returns:
        subjects_data: list of arrays from load_subject_data, controls first
        group_labels: ndarray, 0 = control, 1 = lesion
        records: the matching anInfo records, in the same order
    """
    selected = select_cohort(cohort, data_dir=data_dir,
                             require_lesion_known=require_lesion_known, sex=sex)

    wanted = set()
    if "control" in groups:
        wanted.add(0)
    if "lesion" in groups:
        wanted.add(1)
    selected = [r for r in selected if r["lesion"] in wanted]

    # Controls first, then lesion, each ordered by animal name for determinism.
    selected.sort(key=lambda r: (r["lesion"], r["name"]))

    subjects_data = [load_subject_data(r["path"]) for r in selected]
    group_labels = np.array([r["lesion"] for r in selected], dtype=int)
    return subjects_data, group_labels, selected


def describe(records):
    """One-line summary of a record list, for logging."""
    n_ctrl = sum(1 for r in records if r["lesion"] == 0)
    n_les = sum(1 for r in records if r["lesion"] == 1)
    n_male = sum(1 for r in records if r["sex"] == 0)
    n_female = sum(1 for r in records if r["sex"] == 1)
    return (f"{len(records)} animals: {n_ctrl} control, {n_les} lesion; "
            f"{n_male} male, {n_female} female")


if __name__ == "__main__":
    for name in ("initial", "replication", "all_published", "genotype1"):
        recs = select_cohort(name)
        print(f"{name:16s} {describe(recs)}")
        by_sex = {}
        for r in recs:
            by_sex.setdefault((r["sex"], r["lesion"]), 0)
            by_sex[(r["sex"], r["lesion"])] += 1
        detail = ", ".join(
            f"{'M' if s == 0 else 'F'}/{'ctrl' if l == 0 else 'les'}={n}"
            for (s, l), n in sorted(by_sex.items()))
        print(f"{'':16s}   {detail}")
    excluded = [r for r in read_info() if r["lesion"] is None]
    print(f"\nexcluded (blank lesion): {len(excluded)} animals, "
          f"experiments {sorted({r['experiment'] for r in excluded})}")
