"""The `analyze_offline.py metad` subcommand."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import cv_quality as cvq
import metad_quality as mdq
import figures_metad as figs
from offline_io import cv_columns, dump, expand
from report_io import HAVE_MPL, column, read_plumed
from pmf_metrics import fes_layout


def read_hills(patterns, summary):
    """Hill height against time per walker, and how many each deposited."""
    curves, counts = [], {}
    for path in expand(patterns):
        fields, data = read_plumed(path)
        time_ps = column(fields, data, "time")
        height = column(fields, data, "height")
        if time_ps is None or height is None or not len(height):
            continue
        label = path.parent.name or path.name
        counts[label] = int(len(height))
        curves.append((label, time_ps / 1000.0, height))
        summary.setdefault("hill_decay", {})[label] = mdq.hill_decay(time_ps, height)
    return curves, counts


def read_surfaces(convergence_patterns, final_path, window_kJ, summary):
    """The sum_hills series and the final surface, reduced to two numbers."""
    surfaces = []
    for path in expand(convergence_patterns):
        fields, data = read_plumed(path)
        layout = fes_layout(fields) if data.size else None
        if layout:
            surfaces.append(data[:, layout[1]])
    rmsd = mdq.fes_rmsd(surfaces, window_kJ) if len(surfaces) > 1 else {}
    if rmsd:
        summary["fes_convergence"] = rmsd

    if final_path:
        if not Path(final_path).is_file():
            print(f"analyze_offline: no such file: {final_path}", file=sys.stderr)
        fields, data = read_plumed(Path(final_path))
        layout = fes_layout(fields) if data.size else None
        if layout:
            free = data[:, layout[1]]
            summary["fes"] = {"cv_names": layout[0], "bins": int(len(free)),
                              "depth_kJ_per_mol": float(np.nanmax(free) - np.nanmin(free))}
    return rmsd


def read_monitors(patterns, observables, temperature_K):
    """Pooled log-weights and observables from every walker's monitor file."""
    segments, values_by_name, weights = [], {}, []
    offset = 0
    for path in expand(patterns):
        fields, data = read_plumed(path)
        rbias = column(fields, data, "metad.rbias")
        if rbias is None or not data.size:
            continue
        segments.append((offset, offset + len(data)))
        offset += len(data)
        weights.append(mdq.log_weights_from_rbias(rbias, temperature_K))
        for name in observables:
            series = column(fields, data, name)
            if series is not None:
                values_by_name.setdefault(name, []).append(series)
    return segments, values_by_name, weights


def reweight(segments, values_by_name, weights, args, summary):
    """Free energy with a block error bar along each observable that is there."""
    if not weights:
        summary["reweighting"] = {
            "error": "no metad.rbias column in any monitor file - the run was made without "
                     "CALC_RCT, and a gridless METAD cannot compute c(t)"}
        return [], {}

    log_weights = np.concatenate(weights)
    ess = mdq.effective_sample_size(log_weights)
    summary["reweighting"] = dict(ess, walkers_pooled=len(segments))

    panels = []
    for name, chunks in values_by_name.items():
        if sum(len(c) for c in chunks) != len(log_weights):
            summary.setdefault("skipped_observables", []).append(name)
            continue
        centres, fes, sigma = mdq.block_free_energy(
            np.concatenate(chunks), log_weights, args.temperature_K,
            bins=args.bins, n_blocks=args.blocks, segments=segments)
        if centres.size:
            panels.append((name, centres, fes, sigma))
            summary.setdefault("reweighted_fes", {})[name] = {
                "centres": centres.tolist(),
                "fes_kJ_per_mol": fes.tolist(),
                "sigma_kJ_per_mol": sigma.tolist(),
                "max_sigma_kJ_per_mol": (float(np.nanmax(sigma))
                                         if np.any(np.isfinite(sigma)) else None),
            }
    return panels, ess


def read_bias_history(patterns, temperature_K, summary):
    """c(t), the recrossing count, and the running free-energy difference."""
    for path in expand(patterns):
        fields, data = read_plumed(path)
        bias = column(fields, data, "metad.bias")
        rbias = column(fields, data, "metad.rbias")
        if bias is None or rbias is None:
            continue
        ct = mdq.ct_curve(bias, rbias)
        summary["ct"] = {k: v for k, v in ct.items() if k != "curve"}

        nodes = cv_columns(fields)
        if not nodes:
            return None
        first = column(fields, data, nodes[0])
        summary["recrossings"] = cvq.basin_transitions(first)
        return mdq.running_free_energy_difference(
            first, mdq.log_weights_from_rbias(rbias, temperature_K),
            float(np.median(first)), temperature_K)
    return None


def run(args) -> int:
    out_dir = Path(args.out_dir)
    summary: dict = {"temperature_K": args.temperature_K}
    produced: list[str] = []

    decay_curves, hills_per_walker = read_hills(args.hills, summary)
    if hills_per_walker:
        summary["walker_balance"] = mdq.walker_balance(hills_per_walker)

    rmsd = read_surfaces(args.fes_convergence, args.fes, args.fes_window, summary)

    segments, values_by_name, weights = read_monitors(args.monitor, args.observable,
                                                      args.temperature_K)
    panels, ess = reweight(segments, values_by_name, weights, args, summary)
    running = read_bias_history(args.colvar, args.temperature_K, summary)

    if HAVE_MPL:
        if decay_curves:
            figs.plot_metad_convergence(decay_curves, rmsd, running, out_dir, produced)
        if panels:
            figs.plot_reweighted_with_error(panels, ess, out_dir, produced)
        if hills_per_walker:
            figs.plot_walker_balance(summary.get("walker_balance", {}), out_dir, produced)

    summary["figures"] = produced
    dump(summary, out_dir, "metad_quality.json")
    return 0
