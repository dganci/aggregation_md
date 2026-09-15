"""Shared by the test_analysis_*.py files: the check counter, and the two helpers that write a
PLUMED file and run analyze_offline.py.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))

PASSED = 0
FAILURES = []


def check(name, condition):
    global PASSED
    print(f"[{'PASS' if condition else 'FAIL'}] {name}")
    if condition:
        PASSED += 1
    else:
        FAILURES.append(name)


def close(a, b, tol=1e-9):
    return bool(np.all(np.abs(np.asarray(a) - np.asarray(b)) <= tol))


def write_colvar(path, fields, rows):
    lines = ["#! FIELDS " + " ".join(fields)]
    lines += [" ".join(f"{v:.6f}" for v in row) for row in rows]
    path.write_text("\n".join(lines) + "\n")


def run_cli(*argv):
    return subprocess.run([sys.executable, str(SCRIPTS / "analyze_offline.py"), *argv],
                          capture_output=True, text=True)


def finish():
    """Prints the tally and exits non-zero on any failure."""
    print()
    print(f"{PASSED}/{PASSED + len(FAILURES)} tests passed")
    if FAILURES:
        for name in FAILURES:
            print(f"  failed: {name}")
        raise SystemExit(1)
