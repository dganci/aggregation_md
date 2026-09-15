#!/usr/bin/env python3
"""cv_readiness.py end to end on a series with a known slowest timescale."""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from analysis_common import SCRIPTS, check, finish, write_colvar  # noqa: F401


with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    rng = np.random.default_rng(3)
    n = 40_000
    phi = np.exp(-1.0 / 500.0)
    slow = np.empty(n)
    slow[0] = 0.0
    for i in range(1, n):
        slow[i] = phi * slow[i - 1] + rng.normal(0.0, np.sqrt(1 - phi ** 2))
    t = np.arange(n) * 10.0
    noise = lambda: rng.normal(0.0, 0.3, n)   # noqa: E731 - three uses, one line
    write_colvar(root / "COLVAR",
                 ["time", "d_1_2", "cn_1_2", "cn_total", "rg1", "rg2", "rg_com"],
                 np.column_stack([t, 5.0 + slow, 1.0 + noise(), 1.0 + noise(),
                                  2.0 + noise(), 2.0 + noise(), 2.5 + 0.5 * slow]))
    out = root / "cv_readiness.txt"
    result = subprocess.run([sys.executable, str(SCRIPTS / "cv_readiness.py"),
                             "--colvar", str(root / "COLVAR"), "--n-prot", "2",
                             "--out", str(out)], capture_output=True, text=True)
    check("cv_readiness.py runs", result.returncode == 0)
    values = dict(line.split(None, 1) for line in out.read_text().splitlines()) \
        if out.is_file() else {}
    check("cv_readiness finds the plateau of an AR(1) series", values.get("plateau") == "1")
    its1_frames = float(values.get("its1_ps", 0)) / 10.0
    check("cv_readiness recovers the 500-frame timescale within a factor of two",
          250.0 < its1_frames < 1000.0)
    check("cv_readiness reports the run as many timescales long",
          float(values.get("time_over_its", 0)) > 40.0)

    write_colvar(root / "SHORT",
                 ["time", "d_1_2", "cn_1_2", "cn_total", "rg1", "rg2", "rg_com"],
                 np.column_stack([t[:1500], 5.0 + slow[:1500], 1.0 + noise()[:1500],
                                  1.0 + noise()[:1500], 2.0 + noise()[:1500],
                                  2.0 + noise()[:1500], 2.5 + 0.5 * slow[:1500]]))
    short = root / "short.txt"
    subprocess.run([sys.executable, str(SCRIPTS / "cv_readiness.py"),
                    "--colvar", str(root / "SHORT"), "--n-prot", "2", "--out", str(short)],
                   capture_output=True, text=True)
    short_values = dict(line.split(None, 1) for line in short.read_text().splitlines()) \
        if short.is_file() else {}
    check("a run three timescales long is reported as such, not as converged",
          float(short_values.get("time_over_its", 100)) < 10.0)

    write_colvar(root / "TINY",
                 ["time", "d_1_2", "cn_1_2", "cn_total", "rg1", "rg2", "rg_com"],
                 np.column_stack([t[:3], 5.0 + slow[:3], np.ones(3), np.ones(3),
                                  2.0 * np.ones(3), 2.0 * np.ones(3), 2.5 * np.ones(3)]))
    tiny = root / "tiny.txt"
    result = subprocess.run([sys.executable, str(SCRIPTS / "cv_readiness.py"),
                             "--colvar", str(root / "TINY"), "--n-prot", "2", "--out", str(tiny)],
                            capture_output=True, text=True)
    check("a run shorter than the equilibration window gives a note, not a crash",
          result.returncode == 2 and tiny.is_file() and tiny.read_text().startswith("note "))


finish()
