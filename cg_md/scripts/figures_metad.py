"""Figures 60-62: convergence, reweighted surfaces, walker balance."""

from __future__ import annotations

import numpy as np

from report_io import HAVE_MPL, plt, save, thin


def plot_metad_convergence(decay_curves, rmsd, running, out_dir, produced):
    """Hill decay, FES RMSD against the final surface, and a running dF."""
    if not HAVE_MPL:
        return
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))

    for name, time_ns, height in decay_curves:
        t, h = thin(time_ns, height)
        axes[0].plot(t, h, lw=0.8, label=name)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("time (ns)")
    axes[0].set_ylabel("hill height (kJ/mol)")
    axes[0].set_title("Deposition rate")
    if len(decay_curves) > 1:
        axes[0].legend(fontsize=7)

    if rmsd and "rmsd_kJ_per_mol" in rmsd:
        values = rmsd["rmsd_kJ_per_mol"]
        blocks = range(1, len(values) + 1)
        axes[1].plot(blocks, values, marker="o")
        axes[1].set_xticks(list(blocks))
        axes[1].set_title(f"RMSD to the final surface\n({rmsd['bins_compared']} bins within "
                          f"{rmsd['window_kJ_per_mol']:.0f} kJ/mol)")
    else:
        axes[1].set_title("RMSD to the final surface")
    axes[1].set_xlabel("sum_hills block")
    axes[1].set_ylabel("RMSD (kJ/mol)")

    if running is not None and len(running[0]):
        axes[2].plot(running[0], running[1], marker=".")
        axes[2].axhline(0.0, lw=0.6, color="k")
    axes[2].set_xlabel("frames used")
    axes[2].set_ylabel("dF between the two sides (kJ/mol)")
    axes[2].set_title("Running free-energy difference")

    fig.suptitle("Metadynamics convergence: three independent views")
    save(fig, out_dir, "60_metad_convergence", produced)

def plot_reweighted_with_error(panels, ess, out_dir, produced):
    """Reweighted free energies with block error bands, one panel per observable."""
    if not HAVE_MPL or not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(5.5 * len(panels), 4), squeeze=False)
    for ax, (label, centres, fes, sigma) in zip(axes[0], panels):
        ax.plot(centres, fes, color="tab:blue")
        ax.fill_between(centres, fes - sigma, fes + sigma, alpha=0.3, color="tab:blue")
        ax.set_xlabel(label)
        ax.set_ylabel("free energy (kJ/mol)")
    note = ""
    if ess and "ess_fraction" in ess:
        note = (f"  -  effective sample size {ess['ess']:.0f} of {ess['frames']} frames "
                f"({100 * ess['ess_fraction']:.2f}%)")
    fig.suptitle(f"Reweighted free energy, block error band{note}")
    save(fig, out_dir, "61_reweighted_with_error", produced)

def plot_walker_balance(balance, out_dir, produced):
    """Hills deposited per walker - a stalled walker shows up here and nowhere else."""
    if not HAVE_MPL or not balance or "per_walker" not in balance:
        return
    names = list(balance["per_walker"])
    counts = [balance["per_walker"][n] for n in names]
    fig, ax = plt.subplots(figsize=(max(4.5, 0.7 * len(names) + 2), 3.6))
    ax.bar(range(len(names)), counts)
    ax.axhline(float(np.mean(counts)), color="tab:red", lw=0.9, ls="--", label="mean")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("hills deposited")
    ax.set_title(f"Walker balance (min/max = {balance['min_over_max']:.2f})")
    ax.legend(fontsize=8)
    save(fig, out_dir, "62_walker_balance", produced)
