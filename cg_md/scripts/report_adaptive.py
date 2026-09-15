"""Figures 07-08: the adaptive stop rule, chunk by chunk."""

from __future__ import annotations

import numpy as np

from report_io import HAVE_MPL, plt, save


def plot_adaptive_metrics(records, out_dir, produced, summary):
    if not records:
        return
    summary["adaptive_chunks"] = len(records)
    summary["adaptive_stopped"] = bool(records[-1].get("stop", False))
    summary["adaptive_last_metrics"] = records[-1]

    if not HAVE_MPL:
        return
    chunks = [r.get("last_chunk", i) for i, r in enumerate(records)]
    fig, axes = plt.subplots(2, 3, figsize=(15, 7))
    series = [
        ("cn_effective_samples", "independent samples of cn_total"),
        ("pattern_jsd_halves_bits", "JSD between halves (bits)"),
        ("pattern_growth_ratio", "pattern growth, whole / first half"),
        ("populated_contact_patterns", "patterns with >= 1% of frames"),
        ("total_time_us", "sampled time (us)"),
    ]
    for ax, (key, title) in zip(axes.ravel(), series):
        ax.plot(chunks, [r.get(key, np.nan) for r in records], marker="o")
        ax.set_xlabel("chunk")
        ax.set_title(title)
    ax = axes.ravel()[-1]
    ax.plot(chunks, [r.get("assembly_growth_events", np.nan) for r in records],
            marker="o", label="grew")
    ax.plot(chunks, [r.get("assembly_shrink_events", np.nan) for r in records],
            marker="s", label="shrank")
    ax.set_xlabel("chunk")
    ax.set_title("largest oligomer events (residence-filtered)")
    ax.legend(fontsize=8)
    fig.suptitle("Adaptive sampling: convergence gates per chunk")
    save(fig, out_dir, "07_adaptive_criteria", produced)

    checks = [k for k in records[-1].get("checks", {})]
    if checks:
        grid = np.array([[1.0 if r.get("checks", {}).get(k) else 0.0 for k in checks] for r in records])
        fig, ax = plt.subplots(figsize=(8, 0.4 * len(records) + 2))
        ax.imshow(grid, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
        ax.set_xticks(range(len(checks)))
        ax.set_xticklabels(checks, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(records)))
        ax.set_yticklabels([f"chunk {c}" for c in chunks], fontsize=7)
        ax.set_title("Stop-criteria satisfaction per chunk")
        save(fig, out_dir, "08_adaptive_checks", produced)
