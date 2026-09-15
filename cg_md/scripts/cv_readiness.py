#!/usr/bin/env python3
"""Is this trajectory long enough to train a CV on?"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

from cv_quality import permutation_invariant_features
from report_io import column, load_colvars

MAX_SCAN_FRAMES = 500_000


def find_cvgen_estimators():
    """Imports cvgen_estimators from wherever this checkout keeps cg_cvgen."""
    here = Path(__file__).resolve().parent
    candidates = [os.environ.get("CG_CVGEN_SCRIPTS"),
                  here.parent.parent / "cg_cvgen" / "scripts",
                  Path(os.environ["CGMD_SRC"]) / "cg_cvgen" / "scripts" if "CGMD_SRC" in os.environ else None]
    for candidate in candidates:
        if candidate and (Path(candidate) / "cvgen_estimators.py").is_file():
            sys.path.insert(0, str(candidate))
            import cvgen_estimators  # noqa: E402
            return cvgen_estimators
    return None


def write_result(path, **values):
    lines = [f"{key} {value}" for key, value in values.items()]
    Path(path).write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--colvar", nargs="+", required=True, type=Path)
    ap.add_argument("--n-prot", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--equilibration-ps", type=float, default=50.0,
                    help="frames before this are dropped, as cg_cvgen does")
    ap.add_argument("--dt-ps", type=float, default=2.0,
                    help="frame spacing when the COLVAR has no time column")
    args = ap.parse_args()

    estimators = find_cvgen_estimators()
    if estimators is None:
        write_result(args.out, note="cvgen_estimators.py not found; set CG_CVGEN_SCRIPTS")
        print("cv_readiness: cvgen_estimators.py not found", file=sys.stderr)
        return 2

    fields, data = load_colvars([p for p in args.colvar if p.is_file()])
    if not data.size:
        write_result(args.out, note="no COLVAR rows")
        return 2

    time_ps = column(fields, data, "time")
    if time_ps is not None:
        data = data[time_ps > time_ps[0] + args.equilibration_ps]
        time_ps = column(fields, data, "time")
    dt_ps = float(np.median(np.diff(time_ps))) if time_ps is not None and len(time_ps) > 1 else args.dt_ps

    step = max(1, len(data) // MAX_SCAN_FRAMES)
    data = data[::step]
    dt_ps *= step

    names, x = permutation_invariant_features(fields, data, args.n_prot)
    if x.shape[0] < 100 or x.shape[1] == 0:
        write_result(args.out, note=f"{x.shape[0]} frames x {x.shape[1]} features is too little to scan")
        return 2

    x = (x - x.mean(axis=0)) / np.maximum(x.std(axis=0), 1e-12)
    rows = estimators.PlateauScan.sweep(x, dt_ps, len(x))
    if not rows:
        write_result(args.out, note="no lag on the ladder gave a finite timescale")
        return 2
    band = estimators.PlateauScan.plateau(rows)
    if len(band) < 3:
        band = []

    its1 = float(np.median([r["its1_ps"] for r in band])) if band else float(min(r["its1_ps"] for r in rows))
    total_time_ps = float(len(x) * dt_ps)
    write_result(args.out,
                 plateau=int(bool(band)),
                 its1_ps=f"{its1:.6g}",
                 lag_lo_ps=f"{band[0]['lag_ps']:.6g}" if band else 0,
                 lag_hi_ps=f"{band[-1]['lag_ps']:.6g}" if band else 0,
                 time_over_its=f"{total_time_ps / its1:.4f}" if its1 > 0 else 0,
                 frames=len(x),
                 dt_ps=f"{dt_ps:.6g}",
                 total_time_ps=f"{total_time_ps:.6g}",
                 features=len(names))
    print(f"cv_readiness: ITS-1 {its1 / 1000:.1f} ns, run {total_time_ps / 1000:.1f} ns "
          f"({total_time_ps / its1:.1f}x), plateau {'yes' if band else 'NO'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
