"""Regenerate tests/fixtures/regression_golden.json.

Run this only when a numerical change is intended, and record why in the commit
message. Regenerating to silence a failing test discards the guarantee the
fixture exists to provide.

    pixi run python scripts/update_regression_golden.py
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))
sys.path.insert(0, str(Path(__file__).parents[1] / "tests"))

from test_regression import GOLDEN_PATH, compute_all  # noqa: E402


def main():
    new = compute_all()
    old = json.loads(GOLDEN_PATH.read_text()) if GOLDEN_PATH.exists() else None

    if old == new:
        print(f"unchanged: {GOLDEN_PATH}")
        return 0

    if old is not None:
        print("changes:")
        for key in sorted(set(old) | set(new)):
            if old.get(key) != new.get(key):
                print(f"  {key}:\n    was {old.get(key)}\n    now {new.get(key)}")

    GOLDEN_PATH.parent.mkdir(parents=True, exist_ok=True)
    GOLDEN_PATH.write_text(json.dumps(new, indent=2, sort_keys=True) + "\n")
    print(f"\nwrote {GOLDEN_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
