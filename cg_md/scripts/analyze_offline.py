#!/usr/bin/env python3
"""Offline analysis of aggregation, CV quality and metadynamics, from files."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import offline_aggregation
import offline_compare
import offline_cv
import offline_metad
from report_io import HAVE_MPL


def parse_args(argv=None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    def shared(parser):
        parser.add_argument("--out-dir", required=True, type=Path)
        parser.add_argument("--temperature-K", type=float, default=310.0)

    a = sub.add_parser("aggregation", help="oligomer sizes, lifetimes and assembly times")
    shared(a)
    a.add_argument("--colvar", nargs="+", required=True,
                   help="COLVAR paths or globs; several files are put on one clock")
    a.add_argument("--n-prot", type=int, required=True)
    a.add_argument("--contact-threshold", type=float, default=1.0,
                   help="cn_i_j at or above which two protomers count as touching")
    a.add_argument("--dt-ps", type=float, default=10.0,
                   help="frame spacing, used only when the COLVAR has no time column")
    a.add_argument("--max-lag", type=int, default=2000,
                   help="longest lag of the contact correlation function, in frames")
    a.add_argument("--max-cluster-frames", type=int, default=200_000)
    a.set_defaults(func=offline_aggregation.run)

    c = sub.add_parser("cv", help="quality control of a trained collective variable")
    shared(c)
    c.add_argument("--cv-dir", type=Path, default=None,
                   help="cg_cvgen output directory (its_scan.csv, lag_scores.csv, cv_params.json)")
    c.add_argument("--colvar", nargs="*", default=[],
                   help="biased COLVARs carrying the mycv.node-* columns")
    c.add_argument("--monitor", nargs="*", default=[],
                   help="COLVAR_monitor files, for what the CV corresponds to physically")
    c.add_argument("--observable", nargs="*", default=["cn_total", "rg_com"])
    c.set_defaults(func=offline_cv.run)

    m = sub.add_parser("metad", help="metadynamics convergence and reweighting")
    shared(m)
    m.add_argument("--hills", nargs="*", default=[])
    m.add_argument("--colvar", nargs="*", default=[])
    m.add_argument("--monitor", nargs="*", default=[])
    m.add_argument("--fes", type=Path, default=None, help="final sum_hills surface")
    m.add_argument("--fes-convergence", nargs="*", default=[],
                   help="the fes_*.dat series from sum_hills --stride")
    m.add_argument("--fes-window", type=float, default=30.0,
                   help="compare surfaces only where they are within this many kJ/mol of "
                        "the minimum; empty grid otherwise dominates the RMSD")
    m.add_argument("--observable", nargs="*", default=["cn_total", "rg_com"])
    m.add_argument("--bins", type=int, default=50)
    m.add_argument("--blocks", type=int, default=10,
                   help="blocks per walker for the error bar on the reweighted surface")
    m.set_defaults(func=offline_metad.run)

    r = sub.add_parser("compare", help="do independent replicas give the same size distribution")
    shared(r)
    r.add_argument("--summaries", nargs="+", required=True,
                   help="aggregation_summary.json of each replica (paths or globs)")
    r.add_argument("--max-jsd", type=float, default=0.05,
                   help="largest pairwise Jensen-Shannon divergence, in bits, still called agreement")
    r.set_defaults(func=offline_compare.run)

    return ap.parse_args(argv)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    status = args.func(args)
    if not HAVE_MPL:
        print("analyze_offline: matplotlib not available - the JSON is complete, figures skipped",
              file=sys.stderr)
    return status


if __name__ == "__main__":
    raise SystemExit(main())
