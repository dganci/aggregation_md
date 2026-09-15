"""Figures 10-14: hill decay, the FES and its convergence, the CV trajectory, and the reweighted
free energies.
"""

from __future__ import annotations

import re

import numpy as np

from pmf_metrics import fes_layout
from report_io import HAVE_MPL, column, plt, read_plumed, save, thin
from report_metrics import reweighted_free_energy
from report_structure import plot_structure


def walker_number(path):
    """Sort key that orders walker10 after walker2, not after walker1."""
    digits = re.sub(r"\D", "", path.name)
    return (int(digits) if digits else 0, path.name)

def plot_metadynamics(run_dir, out_dir, produced, summary, temperature_K, n_prot, threshold):
    walkers = sorted((p for p in run_dir.glob("walker*") if p.is_dir()), key=walker_number)
    dirs = walkers if walkers else [run_dir]
    summary["walkers"] = len(dirs)

    hill_curves = []
    for d in dirs:
        fields, data = read_plumed(d / "HILLS")
        height = column(fields, data, "height")
        time = column(fields, data, "time")
        if height is not None and time is not None and len(height):
            hill_curves.append((d.name, time / 1000.0, height))
    if hill_curves:
        summary["hills_deposited"] = int(sum(len(h) for _, _, h in hill_curves))
        summary["final_hill_height_kJ_per_mol"] = float(np.mean([h[-1] for _, _, h in hill_curves]))
        if HAVE_MPL:
            fig, ax = plt.subplots(figsize=(8, 4))
            for name, t, h in hill_curves:
                t, h = thin(t, h)
                ax.plot(t, h, lw=0.8, label=name)
            ax.set_xlabel("time (ns)")
            ax.set_ylabel("hill height (kJ/mol)")
            ax.set_yscale("log")
            ax.set_title("Well-tempered hill height decay (must fall towards zero)")
            if len(hill_curves) > 1:
                ax.legend(fontsize=7)
            save(fig, out_dir, "10_hill_heights", produced)

    fields, fes = read_plumed(run_dir / "fes.dat")
    layout = fes_layout(fields) if fes.size else None
    if layout:
        cv_names, free_idx = layout
        free = fes[:, free_idx]
        summary["fes_depth_kJ_per_mol"] = float(np.nanmax(free) - np.nanmin(free))
        if HAVE_MPL and len(cv_names) == 1:
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.plot(fes[:, 0], free)
            ax.set_xlabel(cv_names[0])
            ax.set_ylabel("free energy (kJ/mol)")
            ax.set_title("Free-energy surface")
            save(fig, out_dir, "11_fes", produced)
        elif HAVE_MPL and len(cv_names) >= 2:
            xs, xi = np.unique(fes[:, 0], return_inverse=True)
            ys, yi = np.unique(fes[:, 1], return_inverse=True)
            if len(xs) * len(ys) == len(free):
                grid = np.full((len(ys), len(xs)), np.nan)
                grid[yi, xi] = free
                fig, ax = plt.subplots(figsize=(6.5, 5))
                mesh = ax.pcolormesh(xs, ys, grid, cmap="viridis_r", shading="auto")
                fig.colorbar(mesh, ax=ax, label="free energy (kJ/mol)")
                ax.set_xlabel(cv_names[0])
                ax.set_ylabel(cv_names[1])
                ax.set_title("Free-energy surface (first two CV components)")
                save(fig, out_dir, "11_fes", produced)

    conv_dir = run_dir / "fes_convergence"
    if conv_dir.is_dir():
        files = sorted(conv_dir.glob("fes_*.dat"),
                       key=lambda p: int(re.sub(r"\D", "", p.stem) or 0))
        surfaces = []
        for path in files:
            f, data = read_plumed(path)
            block_layout = fes_layout(f) if data.size else None
            if block_layout:
                surfaces.append(data[:, block_layout[1]])
        if len(surfaces) > 1:
            reference = surfaces[-1]
            deviations = [float(np.sqrt(np.nanmean((s - reference) ** 2)))
                          for s in surfaces if len(s) == len(reference)]
            summary["fes_convergence_rmsd_kJ_per_mol"] = deviations
            if HAVE_MPL:
                fig, ax = plt.subplots(figsize=(7, 4))
                ax.plot(range(1, len(deviations) + 1), deviations, marker="o")
                ax.set_xlabel("sum_hills block")
                ax.set_ylabel("RMSD to final FES (kJ/mol)")
                ax.set_title("FES convergence")
                save(fig, out_dir, "12_fes_convergence", produced)

    for d in dirs:
        fields, data = read_plumed(d / "COLVAR")
        if not data.size:
            continue
        nodes = [f for f in fields if f.startswith("mycv.node-")]
        time = column(fields, data, "time")
        if HAVE_MPL and nodes and time is not None:
            fig, axes = plt.subplots(len(nodes) + 1, 1, figsize=(10, 2.2 * (len(nodes) + 1)), sharex=True)
            axes = np.atleast_1d(axes)
            series = [column(fields, data, name) for name in nodes]
            bias = column(fields, data, "metad.bias")
            thinned = thin(time / 1000.0, *series, bias)
            for ax, name, values in zip(axes, nodes, thinned[1:]):
                ax.plot(thinned[0], values, lw=0.6)
                ax.set_ylabel(name)
            if bias is not None:
                axes[-1].plot(thinned[0], thinned[-1], lw=0.6, color="tab:red")
                axes[-1].set_ylabel("bias (kJ/mol)")
            axes[-1].set_xlabel("time (ns)")
            fig.suptitle(f"CV trajectory and deposited bias ({d.name})")
            save(fig, out_dir, f"13_cv_trajectory_{d.name}", produced)
        break

    monitors = [read_plumed(d / "COLVAR_monitor") for d in dirs]
    monitors = [(f, m) for f, m in monitors if m.size]
    monitor_fields, monitor = monitors[0] if monitors else (None, None)
    rbias = None
    if monitors and all(column(f, m, "metad.rbias") is not None for f, m in monitors):
        rbias = np.concatenate([column(f, m, "metad.rbias") for f, m in monitors])
        summary["reweighted_walkers"] = len(monitors)

    if monitor is not None and monitor.size:
        plot_structure(monitor_fields, monitor, n_prot, threshold, temperature_K,
                       out_dir, produced, summary, unbiased=False)

    if rbias is not None and monitor is not None and monitor.size and HAVE_MPL:
        observables = [("cn_total", "total contact number"),
                       ("rg_com", "assembly Rg from COMs (nm)")]
        available = [(k, t) for k, t in observables if k in monitor_fields]
        if available:
            fig, axes = plt.subplots(1, len(available), figsize=(5.5 * len(available), 4))
            axes = np.atleast_1d(axes)
            for ax, (key, title) in zip(axes, available):
                pooled = [column(f, m, key) for f, m in monitors]
                if any(series is None for series in pooled):
                    continue
                centers, fes_1d, _ = reweighted_free_energy(
                    np.concatenate(pooled), rbias, temperature_K)
                ax.plot(centers, fes_1d)
                ax.set_xlabel(title)
                ax.set_ylabel("free energy (kJ/mol)")
            fig.suptitle("Reweighted free energy along unbiased observables (c(t) reweighting)")
            save(fig, out_dir, "14_reweighted_fes", produced)
            summary["reweighting"] = "ok"
    elif monitor is not None and monitor.size:
        summary["reweighting"] = ("unavailable: no metad.rbias column - the run was made without "
                                  "CALC_RCT (a gridless METAD cannot compute c(t))")
