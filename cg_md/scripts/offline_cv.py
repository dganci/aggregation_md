"""The `analyze_offline.py cv` subcommand."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import cv_quality as cvq
import figures_cv as figs
from offline_io import CV_PREFIX, cv_columns, dump, expand, join_on_time
from report_io import HAVE_MPL, column, read_plumed


def read_walkers(patterns):
    """CV values per walker, kept separate so recrossings stay within a run."""
    records, names = [], []
    for path in expand(patterns):
        fields, data = read_plumed(path)
        nodes = cv_columns(fields)
        if not nodes or not data.size:
            continue
        names = names or nodes
        if nodes != names:
            print(f"analyze_offline: {path} has CV columns {nodes}, expected {names} - skipped",
                  file=sys.stderr)
            continue
        records.append((path, fields, data,
                        np.column_stack([column(fields, data, n) for n in nodes])))
    return records, names


def interpret(records, monitors, observables, names, summary):
    """Rank correlation of each CV component against the physical descriptors."""
    correlations = []
    for path, fields, data, values in records:
        beside = [m for m in monitors if m.parent == path.parent]
        if not beside:
            print(f"analyze_offline: no monitor file next to {path} - skipping the CV "
                  "interpretation for it", file=sys.stderr)
            continue
        summary["cv_interpretation_source"] = str(beside[0])
        mon_fields, mon = read_plumed(beside[0])
        if not mon.size or "time" not in mon_fields or "time" not in fields:
            continue

        index_cv, index_mon = join_on_time(column(fields, data, "time"),
                                           column(mon_fields, mon, "time"))
        if len(index_cv) < 10:
            print(f"analyze_offline: {path.name} and {beside[0].name} share only "
                  f"{len(index_cv)} time stamps - check --metad-print-stride against "
                  "--plumed-stride", file=sys.stderr)
            continue

        for observable in observables:
            series = column(mon_fields, mon, observable)
            if series is None:
                continue
            for k, name in enumerate(names):
                correlations.append({"component": name, "observable": observable,
                                     **cvq.correlation(values[index_cv, k], series[index_mon])})
        break
    return correlations


def run(args) -> int:
    out_dir = Path(args.out_dir)
    summary: dict = {}
    produced: list[str] = []

    cv_dir = Path(args.cv_dir) if args.cv_dir else None
    params = cvq.read_json(cv_dir / "cv_params.json") if cv_dir else {}
    manifest = cvq.read_json(cv_dir / "cv_manifest.json") if cv_dir else {}
    selected = manifest.get("selected_lag")

    plateau = cvq.plateau_summary(
        cvq.read_csv_rows(cv_dir / "its_scan.csv") if cv_dir else [], selected)
    lag_table = cvq.lag_score_summary(
        cvq.read_csv_rows(cv_dir / "lag_scores.csv") if cv_dir else [], selected)
    summary["cv_dir"] = str(cv_dir) if cv_dir else None
    summary["lag_plateau"] = plateau
    summary["lag_scores"] = lag_table
    summary["hill_widths"] = cvq.sigma_vs_spread(params)

    records, names = read_walkers(args.colvar)
    values = np.vstack([r[3] for r in records]) if records else None

    if records:
        summary["cv_frames"] = int(len(values))
        summary["walkers"] = len(records)
        summary["training_coverage"] = cvq.training_coverage(
            values, params.get("min", []), params.get("max", []), params.get("std", []), names)
        summary["grid_coverage"] = cvq.grid_coverage(
            values, params.get("grid_min", []), params.get("grid_max", []), names)
        summary["component_independence"] = cvq.component_independence(values)
        summary["recrossings"] = {
            names[k]: [cvq.basin_transitions(r[3][:, k]) for r in records]
            for k in range(values.shape[1])
        }
    else:
        summary["warning"] = f"no COLVAR with {CV_PREFIX}* columns matched --colvar"

    correlations = []
    monitors = expand(args.monitor)
    if records and monitors:
        correlations = interpret(records, monitors, args.observable, names, summary)
        summary["cv_interpretation"] = correlations
        summary["cv_interpretation_frames"] = correlations[0]["n"] if correlations else 0

    if HAVE_MPL:
        figs.plot_lag_selection(plateau, lag_table, out_dir, produced)
        if records:
            figs.plot_cv_coverage(None, values, summary["training_coverage"],
                                  summary["grid_coverage"], out_dir, produced)
        if correlations:
            figs.plot_cv_interpretation(correlations, out_dir, produced)

    summary["figures"] = produced
    dump(summary, out_dir, "cv_quality.json")
    return 0
