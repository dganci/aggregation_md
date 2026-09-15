#!/usr/bin/env python3
"""Offline analysis and figures for a cg_md run."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from report_adaptive import plot_adaptive_metrics
from report_health import plot_health
from report_metad import plot_metadynamics
from report_pmf import plot_pmf
from report_structure import plot_structure
from report_io import HAVE_MPL, load_colvars, read_jsonl, sampled_time_us


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Offline analysis and figures for a cg_md run.")
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--mode", default="production",
                    choices=["production", "adaptive", "metadynamics", "pmf"])
    ap.add_argument("--n-prot", type=int, required=True)
    ap.add_argument("--dt-ps", type=float, default=0.01)
    ap.add_argument("--plumed-stride", type=int, default=100)
    ap.add_argument("--temperature-K", type=float, default=300.0)
    ap.add_argument("--contact-threshold", type=float, default=1.0)
    ap.add_argument("--biasfactor", type=float, default=10.0)
    ap.add_argument("--pmf-wall-nm", type=float, default=6.0)
    ap.add_argument("--pmf-bound-cutoff-nm", type=float, default=0.0)
    ap.add_argument("--pmf-plateau-frac", type=float, default=0.75)
    ap.add_argument("--out-dir", type=Path, default=None)
    return ap.parse_args()


def colvar_paths(run_dir: Path, mode: str) -> list[Path]:
    """The unbiased COLVARs of a production or adaptive run, in time order."""
    paths = sorted(run_dir.glob("COLVAR_chunk_*.dat")) if mode == "adaptive" else [run_dir / "COLVAR"]
    return [p for p in paths if p.exists()]


def report_unbiased(run_dir: Path, out_dir: Path, produced: list[str], summary: dict,
                    args: argparse.Namespace) -> None:
    """Production and adaptive runs: one continuous COLVAR, plus the stop criteria."""
    paths = colvar_paths(run_dir, args.mode)
    summary["colvar_files"] = [p.name for p in paths]

    fields, data = load_colvars(paths)
    if data.size:
        summary["n_frames"] = int(len(data))
        span_us = sampled_time_us(fields, data)
        if span_us is not None:
            summary["sampled_time_us"] = span_us
        plot_structure(fields, data, args.n_prot, args.contact_threshold,
                       args.temperature_K, out_dir, produced, summary)
    else:
        summary["warning"] = "no COLVAR data found"

    if args.mode == "adaptive":
        plot_adaptive_metrics(read_jsonl(run_dir / "adaptive_sampling_metrics.jsonl"),
                              out_dir, produced, summary)


def report_biased(run_dir: Path, out_dir: Path, produced: list[str], summary: dict,
                  args: argparse.Namespace) -> None:
    """Metadynamics, and the PMF - which is metadynamics along one distance."""
    summary["biasfactor"] = args.biasfactor
    if args.mode == "pmf":
        summary["pmf_wall_nm"] = args.pmf_wall_nm
        plot_pmf(run_dir, out_dir, produced, summary, args.temperature_K, args)
    plot_metadynamics(run_dir, out_dir, produced, summary, args.temperature_K,
                      args.n_prot, args.contact_threshold)


def main() -> int:
    args = parse_args()

    run_dir: Path = args.run_dir
    if not run_dir.is_dir():
        print(f"analyze_run: no such run directory: {run_dir}", file=sys.stderr)
        return 1

    out_dir: Path = args.out_dir or (run_dir / "report")
    out_dir.mkdir(parents=True, exist_ok=True)

    produced: list[str] = []
    summary: dict = {
        "run_dir": str(run_dir),
        "mode": args.mode,
        "n_prot": args.n_prot,
        "temperature_K": args.temperature_K,
        "matplotlib": HAVE_MPL,
    }

    plot_health(read_jsonl(run_dir / "diagnostics.jsonl"), out_dir, produced)

    if args.mode in ("metadynamics", "pmf"):
        report_biased(run_dir, out_dir, produced, summary, args)
    else:
        report_unbiased(run_dir, out_dir, produced, summary, args)

    summary["figures"] = produced
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"analyze_run: wrote {len(produced)} figure(s) and summary.json to {out_dir}")
    if not HAVE_MPL:
        print("analyze_run: matplotlib not available - numeric summary written, figures skipped",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
