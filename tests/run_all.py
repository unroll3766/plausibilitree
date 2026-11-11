"""Run every regression suite in tests/ and aggregate the exit codes.

    conda run -n ssf python tests/run_all.py [--keep]
"""
from __future__ import annotations

import argparse
import sys

import test_cross_artifact_invariants
import test_placement_reflection
import test_review_round_regression


SUITES = [
    ("placement & frame", test_placement_reflection.main),
    ("cross-artifact invariants", test_cross_artifact_invariants.main),
    ("review-round regressions (W21-W27)", test_review_round_regression.main),
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--keep", action="store_true",
                        help="keep fixture directories for debugging")
    args = parser.parse_args()
    failed = []
    for name, fn in SUITES:
        print(f"\n=== {name} ===")
        rc = fn(keep=args.keep)
        print(f"=== {name}: exit {rc} ===")
        if rc:
            failed.append(name)
    print()
    if failed:
        print(f"FAILED SUITES: {failed}")
        return 1
    print("all suites passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
