"""The `analyze_offline.py aggregation` subcommand."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import aggregation_metrics as agg
import figures_aggregation as figs
from offline_io import dump, expand, frame_spacing
from report_io import HAVE_MPL, column, load_colvars


def run(args) -> int:
    paths = expand(args.colvar)
    if not paths:
        print("analyze_offline: no COLVAR matched --colvar", file=sys.stderr)
        return 1

    fields, data = load_colvars(paths)
    if not data.size:
        print("analyze_offline: the COLVARs held no data rows", file=sys.stderr)
        return 1

    pairs, cols = agg.pair_columns(fields, args.n_prot)
    if not cols:
        print(f"analyze_offline: no cn_i_j columns for {args.n_prot} protomers in "
              f"{', '.join(p.name for p in paths)}", file=sys.stderr)
        return 1

    time_ps = column(fields, data, "time")
    dt_ps = frame_spacing(time_ps, args.dt_ps)
    if time_ps is None:
        time_ps = np.arange(len(data), dtype=float) * dt_ps

    summary = {
        "files": [str(p) for p in paths],
        "n_prot": args.n_prot,
        "frames": int(len(data)),
        "frame_spacing_ps": dt_ps,
        "sampled_time_us": float((time_ps[-1] - time_ps[0] + dt_ps) / 1e6),
        "contact_threshold": args.contact_threshold,
        "pairs_found": len(cols),
        "pairs_expected": args.n_prot * (args.n_prot - 1) // 2,
        "temperature_K": args.temperature_K,
    }

    contact = data[:, cols] >= args.contact_threshold
    summary["contact_lifetimes"] = agg.contact_lifetimes(contact, dt_ps)
    curve = agg.contact_autocorrelation(contact, max_lag=min(args.max_lag, len(contact) // 4))
    tau_ps = agg.correlation_time(curve, dt_ps)
    summary["contact_lifetimes"]["intermittent_tau_ps"] = tau_ps

    step = max(1, len(data) // args.max_cluster_frames)
    frames = np.arange(0, len(data), step)
    sizes = agg.cluster_sizes(agg.cluster_labels(contact[frames], pairs, args.n_prot),
                              args.n_prot)
    largest = sizes.max(axis=1)

    counts = agg.size_histogram(sizes, args.n_prot)
    distribution = agg.size_distribution(counts)
    free_energy = agg.size_free_energy(counts, args.temperature_K)

    summary["clustering"] = {
        "frames_analysed": int(len(sizes)),
        "subsample_step": int(step),
        "largest_mean": float(largest.mean()),
        "largest_max": int(largest.max()),
        "clusters_mean": float(np.count_nonzero(sizes, axis=1).mean()),
    }
    summary["size_distribution"] = distribution
    summary["size_free_energy_kJ_per_mol"] = free_energy.tolist()
    summary["censoring"] = agg.censoring_report(largest, args.n_prot)
    summary["first_passage_ps"] = {
        str(target): agg.first_passage(largest, target, time_ps[frames])
        for target in range(2, args.n_prot + 1)
    }

    out_dir = Path(args.out_dir)
    produced: list[str] = []
    if HAVE_MPL:
        running_w = np.cumsum((sizes.astype(float) ** 2).sum(axis=1)) / (
            args.n_prot * np.arange(1, len(sizes) + 1))
        figs.plot_size_distribution(distribution, free_energy, out_dir, produced)
        figs.plot_aggregation_kinetics(time_ps[frames] / 1000.0, largest, running_w,
                                       out_dir, produced)
        complete = []
        for k in range(contact.shape[1]):
            run, censored = agg.contact_runs(contact[:, k])
            complete.append(run[~censored])
        lengths = np.concatenate(complete) if complete else np.empty(0)
        figs.plot_contact_lifetimes(lengths, dt_ps, curve, tau_ps, out_dir, produced)

    summary["figures"] = produced
    dump(summary, out_dir, "aggregation_summary.json")
    return 0
