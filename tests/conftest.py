"""Shared handling for tests that need data the repository cannot ship.

Two datasets are referenced by the suite but tracked by neither repo: the Igor
reference data under `igor_cbas/`, and the unpublished hippocampal lesion cohort
under `data/`. Tests depending on them skip when they are absent, which is what
makes the suite runnable by anyone who clones the project.

Skipping has one failure mode worth guarding: a typo in a path constant is
indistinguishable from "the data is not here", so a broken path would skip
silently and forever, including on a machine that does have the data.

The guard cannot key on whether some directory exists, because another developer
may well have their own `data/` directory holding unrelated data. It has to key on
intent instead. Set the environment variable to say "I have these datasets and I
expect these tests to run", and a missing path becomes a failure that names the
path it looked for:

    PYCBAS_REQUIRE_REFERENCE_DATA=1 pixi run python -m pytest tests/ -q

or equivalently `pixi run test-reference`. Nobody else is affected, because
without the variable the behaviour is an ordinary skip.
"""

import os

import pytest

REQUIRE_ENV = "PYCBAS_REQUIRE_REFERENCE_DATA"

REPO_ROOT = __import__("pathlib").Path(__file__).parent.parent

#: Igor's published reference data, used by the cross-validation in test_cbas.py.
IGOR_DATA_DIR = REPO_ROOT / "igor_cbas" / "data"

#: Kastner's unpublished multi-contingency lesion cohort.
LESION_COHORT_DIR = REPO_ROOT / "data" / "rats_AllHipLesionData"

#: Igor's allTrialToPerfect output for the 4th order / 100 criterion setting.
CRITERION_REFERENCE = REPO_ROOT / "igor_cbas" / "allTrialToPerfect.txt"


def strict_reference_data():
    """Whether missing reference data should fail rather than skip."""
    return os.environ.get(REQUIRE_ENV, "").strip() not in ("", "0", "false", "False")


def require_reference_path(path, description):
    """Skip, or fail under the strict flag, when an untracked path is absent.

    Naming the exact path in both messages is the point: a typo then shows up as
    "looked for <wrong path>" rather than as a bare skip.
    """
    if path.exists():
        return path

    detail = f"{description} not found at {path.relative_to(REPO_ROOT)}"
    if strict_reference_data():
        pytest.fail(
            f"{detail}. {REQUIRE_ENV} is set, so this is treated as a broken "
            "configuration rather than absent data. Fix the path, place the data, "
            f"or unset {REQUIRE_ENV}."
        )
    pytest.skip(f"{detail} (set {REQUIRE_ENV}=1 to make this a failure)")


@pytest.fixture(scope="session")
def igor_data_dir():
    return require_reference_path(IGOR_DATA_DIR, "Igor reference data")


@pytest.fixture(scope="session")
def lesion_cohort_dir():
    return require_reference_path(LESION_COHORT_DIR, "unpublished lesion cohort")


@pytest.fixture(scope="session")
def criterion_reference_file():
    return require_reference_path(CRITERION_REFERENCE,
                                  "Igor allTrialToPerfect reference")
