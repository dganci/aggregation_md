"""Figures 01-06: the structural panels every run gets."""

from __future__ import annotations

import numpy as np

from report_io import HAVE_MPL, column, plt, sampled_time_us, save, thin
from report_metrics import (autocorrelation_time, block_average_error, cluster_timeseries,
                            contact_matrix, free_energy_2d)


def plot_structure(fields, data, n_prot, threshold, temperature_K, out_dir, produced, summary,
                   unbiased=True):
    time_ps = column(fields, data, "time")
    if time_ps is None:
        time_ps = np.arange(len(data), dtype=float)
    time_ns = time_ps / 1000.0

    cn = column(fields, data, "cn_total")
    rg = column(fields, data, "rg_com")
    largest, n_clusters, cluster_frames = cluster_timeseries(fields, data, n_prot, threshold)
    cluster_time_ns = time_ns[cluster_frames] if cluster_frames is not None else None

    summary["n_frames"] = int(len(data))
    span_us = sampled_time_us(fields, data)
    if span_us is not None:
        summary["sampled_time_us"] = span_us
    if not unbiased:
        summary["ensemble"] = "biased (the means below are not reweighted)"

    if cn is not None:
        summary["cn_total"] = {
            "mean": float(cn.mean()), "std": float(cn.std()),
            "min": float(cn.min()), "max": float(cn.max()),
            "autocorr_frames": autocorrelation_time(cn),
        }
    if rg is not None:
        summary["rg_com_nm"] = {
            "mean": float(rg.mean()), "std": float(rg.std()),
            "min": float(rg.min()), "max": float(rg.max()),
            "autocorr_frames": autocorrelation_time(rg),
        }
    if largest is not None:
        counts = np.bincount(largest, minlength=n_prot + 1)[1:]
        summary["largest_cluster"] = {
            "mean": float(largest.mean()),
            "max": int(largest.max()),
            "distribution": {str(i + 1): int(c) for i, c in enumerate(counts)},
            "fraction_fully_assembled": float(np.mean(largest == n_prot)),
            "frames_analysed": int(len(largest)),
            "subsampled": bool(cluster_frames is not None and len(largest) < len(data)),
        }

    if not HAVE_MPL:
        return

    plot_t, plot_cn, plot_rg = thin(time_ns, cn, rg)
    plot_ct, plot_lcc, plot_nc = thin(cluster_time_ns, largest, n_clusters)

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    if cn is not None:
        axes[0].plot(plot_t, plot_cn, lw=0.7)
        axes[0].set_ylabel("cn_total")
    if rg is not None:
        axes[1].plot(plot_t, plot_rg, lw=0.7, color="tab:green")
        axes[1].set_ylabel("Rg global (nm)")
    if largest is not None:
        axes[2].plot(plot_ct, plot_lcc, lw=0.7, color="tab:red", label="largest oligomer")
        axes[2].plot(plot_ct, plot_nc, lw=0.7, color="tab:blue", alpha=0.6,
                     label="n clusters (monomers included)")
        axes[2].set_ylabel("cluster size / count")
        axes[2].legend(fontsize=8)
    axes[2].set_xlabel("time (ns)")
    fig.suptitle("Aggregation order parameters")
    save(fig, out_dir, "01_order_parameters", produced)

    if largest is not None:
        fig, ax = plt.subplots(figsize=(6, 4))
        counts = np.bincount(largest, minlength=n_prot + 1)[1:]
        ax.bar(np.arange(1, n_prot + 1), counts / counts.sum())
        ax.set_xlabel("largest oligomer size")
        ax.set_ylabel("population")
        ax.set_title("Oligomer size distribution")
        save(fig, out_dir, "02_oligomer_distribution", produced)

    matrix = contact_matrix(fields, data, n_prot)
    if np.any(np.isfinite(matrix)):
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(matrix, cmap="viridis", origin="lower",
                       extent=(0.5, n_prot + 0.5, 0.5, n_prot + 0.5))
        fig.colorbar(im, ax=ax, label="mean contact number")
        ax.set_xlabel("protomer")
        ax.set_ylabel("protomer")
        ax.set_title("Mean inter-protomer contact map")
        save(fig, out_dir, "03_contact_map", produced)

    per_chain = [column(fields, data, f"rg{i}") for i in range(1, n_prot + 1)]
    per_chain = [c for c in per_chain if c is not None]
    if per_chain:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        thinned = thin(time_ns, *per_chain)
        for i, series in enumerate(thinned[1:]):
            axes[0].plot(thinned[0], series, lw=0.5, alpha=0.7, label=f"chain {i + 1}")
        axes[0].set_xlabel("time (ns)")
        axes[0].set_ylabel("Rg (nm)")
        axes[0].set_title("Per-chain radius of gyration")
        if len(per_chain) <= 12:
            axes[0].legend(fontsize=6, ncol=2)
        chain_labels = [str(i + 1) for i in range(len(per_chain))]
        try:
            axes[1].boxplot(per_chain, tick_labels=chain_labels)
        except TypeError:
            axes[1].boxplot(per_chain, labels=chain_labels)
        axes[1].set_xlabel("protomer")
        axes[1].set_ylabel("Rg (nm)")
        axes[1].set_title("Per-chain Rg distribution")
        save(fig, out_dir, "04_per_chain_rg", produced)
        summary["per_chain_rg_nm"] = [
            {"mean": float(s.mean()), "std": float(s.std())} for s in per_chain
        ]

    if unbiased and cn is not None and rg is not None and len(cn) > 50:
        fes, xe, ye = free_energy_2d(cn, rg, temperature_K)
        fig, ax = plt.subplots(figsize=(6.5, 5))
        mesh = ax.pcolormesh(xe, ye, fes.T, cmap="viridis_r", shading="auto")
        fig.colorbar(mesh, ax=ax, label="free energy (kJ/mol)")
        ax.set_xlabel("cn_total")
        ax.set_ylabel("Rg global (nm)")
        ax.set_title("-kT ln P (unbiased estimate)")
        save(fig, out_dir, "05_landscape_unbiased", produced)

    if cn is not None and len(cn) > 64:
        lengths, errors = block_average_error(cn)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        axes[0].semilogx(lengths, errors, marker="o")
        axes[0].set_xlabel("block length (frames)")
        axes[0].set_ylabel("standard error of mean")
        axes[0].set_title("Block averaging (plateau = true error bar)")
        running = np.cumsum(cn) / np.arange(1, len(cn) + 1)
        rt, rv = thin(time_ns, running)
        axes[1].plot(rt, rv)
        axes[1].set_xlabel("time (ns)")
        axes[1].set_ylabel("running mean cn_total")
        axes[1].set_title("Running average")
        save(fig, out_dir, "06_convergence", produced)
        if len(errors):
            worst = int(np.argmax(errors))
            summary.setdefault("cn_total", {})["block_error_max"] = float(errors[worst])
            summary["cn_total"]["block_error_at_frames"] = int(lengths[worst])
