"""Figures 40-43: sizes, kinetics, lifetimes, replica agreement."""

from __future__ import annotations

import numpy as np

from report_io import HAVE_MPL, plt, save, thin


def plot_size_distribution(dist, free_energy, out_dir, produced):
    """Oligomer size distribution, both weightings, and on a free-energy axis."""
    if not HAVE_MPL or "error" in dist:
        return
    sizes = np.asarray(dist["sizes"])
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))

    axes[0].bar(sizes, dist["number_fraction"])
    axes[0].set_ylabel("fraction of oligomers")
    axes[0].set_title(f"Number-weighted, <s> = {dist['mean_size_number']:.2f}")

    axes[1].bar(sizes, dist["mass_fraction"], color="tab:orange")
    axes[1].set_ylabel("fraction of chains")
    axes[1].set_title(f"Mass-weighted, <s>_w = {dist['mean_size_weight']:.2f}")

    axes[2].plot(sizes, free_energy, marker="o")
    axes[2].axhline(0.0, lw=0.6, color="k")
    axes[2].set_ylabel("-kT ln (N_s / N_1)  (kJ/mol)")
    axes[2].set_title("Referenced to the free chain")

    for ax in axes:
        ax.set_xlabel("oligomer size")
        ax.set_xticks(sizes)
    fig.suptitle("Oligomer size distribution (ceiling at N is the box, not the physics)")
    save(fig, out_dir, "40_size_distribution", produced)

def plot_aggregation_kinetics(time_ns, largest, weight_avg, out_dir, produced):
    """Largest oligomer and the running weight-average size against time."""
    if not HAVE_MPL or largest is None or len(largest) == 0:
        return
    t, big, avg = thin(time_ns, largest, weight_avg)
    fig, ax = plt.subplots(figsize=(9, 3.8))
    ax.plot(t, big, lw=0.7, color="tab:red", label="largest oligomer")
    if avg is not None:
        ax.plot(t, avg, lw=1.2, color="tab:blue", label="running <s>_w")
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("chains")
    ax.legend(fontsize=8)
    ax.set_title("Assembly over time")
    save(fig, out_dir, "41_aggregation_kinetics", produced)

def plot_contact_lifetimes(run_lengths, dt_ps, autocorr, tau_ps, out_dir, produced):
    """Uninterrupted contact durations, and the intermittent correlation."""
    if not HAVE_MPL:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))

    if run_lengths is not None and len(run_lengths):
        axes[0].hist(np.asarray(run_lengths) * dt_ps / 1000.0, bins=40, log=True)
    axes[0].set_xlabel("uninterrupted contact duration (ns)")
    axes[0].set_ylabel("count")
    axes[0].set_title("Continuous lifetimes (any break ends a run)")

    if autocorr is not None and len(autocorr):
        lags = np.arange(len(autocorr)) * dt_ps / 1000.0
        axes[1].plot(lags, autocorr)
        axes[1].axhline(0.0, lw=0.6, color="k")
    axes[1].set_xlabel("lag (ns)")
    axes[1].set_ylabel("C(t) = <h(0)h(t)> / <h>")
    title = "Intermittent correlation"
    if tau_ps is not None and np.isfinite(tau_ps):
        title += f", tau = {tau_ps / 1000.0:.2f} ns"
    axes[1].set_title(title)
    save(fig, out_dir, "42_contact_lifetimes", produced)

def plot_replica_agreement(labels, number, pooled, pairs, out_dir, produced):
    """Every replica's size distribution, the pooled one with its error, the JSDs."""
    if not HAVE_MPL or len(labels) < 2:
        return
    sizes = np.asarray(pooled["sizes"])
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    width = 0.8 / len(labels)
    for k, (label, dist) in enumerate(zip(labels, number)):
        axes[0].bar(sizes + (k - 0.5 * (len(labels) - 1)) * width, dist, width, label=label)
    axes[0].errorbar(sizes, pooled["number_fraction_mean"], yerr=pooled["number_fraction_sem"],
                     fmt="k_", capsize=4, lw=1.4, label="pooled +/- SEM")
    axes[0].set_xlabel("oligomer size")
    axes[0].set_ylabel("fraction of oligomers")
    axes[0].set_xticks(sizes)
    axes[0].set_title("Size distribution per replica")
    axes[0].legend(fontsize=7)

    index = {label: k + 1 for k, label in enumerate(labels)}
    names = [f"{index[p['a']]}-{index[p['b']]}" for p in pairs]
    axes[1].bar(range(len(pairs)), [p["jsd_bits"] for p in pairs])
    axes[1].set_xticks(range(len(pairs)))
    axes[1].set_xticklabels(names, fontsize=8)
    axes[1].set_ylabel("Jensen-Shannon divergence (bits)")
    axes[1].set_title("Pairwise disagreement (replicas numbered as in the legend)")
    fig.suptitle("Do independent replicas give the same answer?")
    save(fig, out_dir, "43_replica_agreement", produced)
