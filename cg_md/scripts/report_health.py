"""Figure 00: temperature, pressure, density and energy per stage."""

from __future__ import annotations

from report_io import HAVE_MPL, plt, save


def plot_health(diag_records, out_dir, produced):
    """Temperature / pressure / density / total-energy mean and drift per stage."""
    if not HAVE_MPL or not diag_records:
        return
    labels = [r.get("segment", "?") for r in diag_records]
    keys = [
        ("temperature_K", "Temperature (K)"),
        ("pressure_bar", "Pressure (bar)"),
        ("density_kg_per_m3", "Density (kg/m3)"),
        ("total_energy_kJ_per_mol", "Total energy (kJ/mol)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for ax, (key, title) in zip(axes.ravel(), keys):
        means, stds, drifts, xs = [], [], [], []
        for i, rec in enumerate(diag_records):
            block = rec.get("energy", {}).get(key)
            if not block:
                continue
            xs.append(i)
            means.append(block["mean"])
            stds.append(block["stddev"])
            drifts.append(block["drift_per_ns"])
        if not xs:
            ax.set_axis_off()
            continue
        ax.errorbar(xs, means, yerr=stds, marker="o", capsize=3, label="mean +/- sd")
        ax.set_title(title)
        ax.set_xticks(xs)
        ax.set_xticklabels([labels[i] for i in xs], rotation=45, ha="right", fontsize=7)
        twin = ax.twinx()
        twin.bar(xs, drifts, alpha=0.25, color="tab:red")
        twin.set_ylabel("drift / ns", color="tab:red", fontsize=8)
        twin.axhline(0.0, color="tab:red", lw=0.6)
    fig.suptitle("Run health per stage (drift, not fluctuation size, is the red flag)")
    save(fig, out_dir, "00_run_health", produced)
