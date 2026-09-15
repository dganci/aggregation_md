#!/usr/bin/env python3
"""Regression tests for rebase_clocks."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from cvgen_estimators import rebase_clocks  # noqa: E402

PASSED = 0
FAILURES = []


def check(name, cond):
    global PASSED
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    if cond:
        PASSED += 1
    else:
        FAILURES.append(name)


def increasing(arrays):
    joined = np.concatenate(arrays)
    return bool(np.all(np.diff(joined) > 0))


one = np.arange(0.0, 100.0, 10.0)
out = rebase_clocks([one])
check("a single continuous file is returned unchanged", np.array_equal(out[0], one))

chunk = np.arange(0.0, 100.0, 2.0)
out = rebase_clocks([chunk.copy(), chunk.copy(), chunk.copy()])
check("restarted chunks are chained instead of discarded", increasing(out))
check("chunk 1 starts one frame after chunk 0 ends", out[1][0] == 100.0)
check("chunk 2 starts one frame after chunk 1 ends", out[2][0] == 200.0)
check("no frame is lost", sum(len(a) for a in out) == 3 * len(chunk))
check("the last frame carries the whole sampled time", out[-1][-1] == 298.0)

batched = np.concatenate([np.arange(0.0, 50.0, 10.0), np.arange(0.0, 50.0, 10.0)])
out = rebase_clocks([batched])
check("a mid-file restart is rebased too", increasing(out))
check("the second batch continues the first", out[0][5] == 50.0)
check("a mid-file restart loses no row", len(out[0]) == len(batched))

absolute = [np.array([0.0, 10.0, 20.0]), np.array([20.0, 30.0, 40.0])]
out = rebase_clocks([a.copy() for a in absolute])
check("absolute-time files are left alone", all(np.array_equal(a, b)
                                                for a, b in zip(out, absolute)))

out = rebase_clocks([np.array([]), np.array([5.0]), np.array([])])
check("empty and single-frame arrays are handled", len(out) == 3 and out[1][0] == 5.0)
check("nothing is returned for an empty list", rebase_clocks([]) == [])

src = np.array([0.0, 10.0, 0.0, 10.0])
out = rebase_clocks([src.copy()])
check("the input array is not mutated in place", np.array_equal(src, [0.0, 10.0, 0.0, 10.0]))

print()
print(f"{PASSED}/{PASSED + len(FAILURES)} tests passed")
if FAILURES:
    for name in FAILURES:
        print(f"  failed: {name}")
    raise SystemExit(1)
