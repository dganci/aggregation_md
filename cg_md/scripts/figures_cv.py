"""Figures 50-52: lag selection, CV coverage, CV interpretation."""

from __future__ import annotations

import numpy as np

from report_io import HAVE_MPL, plt, save


def plot_lag_selection(plateau, lag_table, out_dir, produced):
    """Implied timescales against lag, and the VAMP-2 scores beside them."""
    if not HAVE_MPL or "error" in plateau:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    lag_ns = np.asarray(plateau["lag_ns"], dtype=float)
    its_ns = np.asarray(plateau["its1_ns"], dtype=float)
    frames = np.asarray(plateau["lag_frames"], dtype=float)
    axes[0].plot(lag_ns, its_ns, marker="o", label="ITS-1")
    axes[0].plot(lag_ns, lag_ns, ls=":", color="grey", label="ITS = lag")

    band = set(plateau.get("plateau_lags") or [])
    if band:
        inside = lag_ns[np.isin(frames, list(band))]
        if inside.size:
            axes[0].axvspan(inside.min(), inside.max(), color="tab:green", alpha=0.15,
                            label="plateau")
    for lag in plateau.get("selected_lags") or []:
        picked = lag_ns[frames == lag]
        if picked.size:
            axes[0].axvline(picked[0], color="tab:green", lw=0.9)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("lag (ns)")
    axes[0].set_ylabel("implied timescale (ns)")
    axes[0].set_title("Plateau: where ITS stops depending on the lag")
    axes[0].legend(fontsize=8)

    if "lags" in lag_table and lag_table["lags"]:
        lags = [r["lag"] for r in lag_table["lags"]]
        axes[1].plot(lags, [r["vamp2_train"] for r in lag_table["lags"]],
                     marker="o", label="VAMP-2 train")
        axes[1].plot(lags, [r["vamp2_test"] for r in lag_table["lags"]],
                     marker="s", label="VAMP-2 test")
        for row in lag_table["lags"]:
            if not row["valid"]:
                axes[1].axvline(row["lag"], color="tab:red", alpha=0.15, lw=6)
        axes[1].set_xlabel("lag (frames)")
        axes[1].set_ylabel("VAMP-2")
        axes[1].set_title("Train vs test (red = failed validity)")
        axes[1].legend(fontsize=8)
    else:
        axes[1].set_axis_off()

    save(fig, out_dir, "50_cv_lag_selection", produced)

def plot_cv_coverage(train, produced_values, coverage, grid, out_dir, produced):
    """Training and production CV distributions, with the METAD grid drawn on."""
    if not HAVE_MPL or not coverage:
        return
    n = len(coverage)
    fig, axes = plt.subplots(1, n, figsize=(5.0 * n, 3.8), squeeze=False)
    for k, record in enumerate(coverage):
        ax = axes[0][k]
        if train is not None and len(train):
            ax.hist(train[:, k], bins=60, density=True, alpha=0.55, label="training")
        if produced_values is not None and len(produced_values):
            ax.hist(produced_values[:, k], bins=60, density=True, alpha=0.55, label="biased run")
        ax.axvline(record["train_min"], color="tab:blue", ls="--", lw=0.9,
                   label="training range")
        ax.axvline(record["train_max"], color="tab:blue", ls="--", lw=0.9)
        if k < len(grid):
            ax.axvline(grid[k]["grid_min"], color="k", ls=":", lw=0.9, label="METAD grid")
            ax.axvline(grid[k]["grid_max"], color="k", ls=":", lw=0.9)
        ax.set_xlabel(record["component"])
        ax.set_ylabel("density")
        ax.set_title(f"{100 * record['fraction_outside']:.1f}% outside the training range")
        ax.legend(fontsize=7)
    fig.suptitle("Where the CV was used against where it was fitted")
    save(fig, out_dir, "51_cv_coverage", produced)

def plot_cv_interpretation(correlations, out_dir, produced):
    """Rank correlation of every CV component with the physical descriptors."""
    if not HAVE_MPL or not correlations:
        return
    components = sorted({c["component"] for c in correlations})
    observables = sorted({c["observable"] for c in correlations})
    grid = np.full((len(components), len(observables)), np.nan)
    for record in correlations:
        grid[components.index(record["component"]),
             observables.index(record["observable"])] = record["spearman"]

    fig, ax = plt.subplots(figsize=(1.6 * len(observables) + 3, 0.8 * len(components) + 2.5))
    image = ax.imshow(grid, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    fig.colorbar(image, ax=ax, label="Spearman correlation")
    ax.set_xticks(range(len(observables)))
    ax.set_xticklabels(observables, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(components)))
    ax.set_yticklabels(components, fontsize=8)
    for i in range(len(components)):
        for j in range(len(observables)):
            if np.isfinite(grid[i, j]):
                ax.text(j, i, f"{grid[i, j]:+.2f}", ha="center", va="center", fontsize=8)
    ax.set_title("What the CV components track")
    save(fig, out_dir, "52_cv_interpretation", produced)
