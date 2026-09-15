"""Figure 20: the dimerization PMF and its standard-state free energy."""

from __future__ import annotations

import sys

from pmf_metrics import dimerization_free_energy, fes_layout
from report_io import HAVE_MPL, plt, read_plumed, save


def plot_pmf(run_dir, out_dir, produced, summary, temperature_K, args):
    """Dimerization PMF: F(r) from sum_hills, w(r), and the standard-state dG."""
    fields, fes = read_plumed(run_dir / "fes.dat")
    layout = fes_layout(fields) if fes.size else None
    if not layout:
        summary["warning"] = "no fes.dat - run the fes stage first"
        return
    cv_names, free_idx = layout
    r, f = fes[:, 0], fes[:, free_idx]

    result = dimerization_free_energy(r, f, temperature_K, args.pmf_wall_nm,
                                      args.pmf_plateau_frac, args.pmf_bound_cutoff_nm)
    if "error" in result:
        summary["pmf_error"] = result["error"]
        print(f"analyze_run: {result['error']}", file=sys.stderr)
        return

    summary["pmf"] = {k: v for k, v in result.items()
                      if k not in ("r_nm", "w_kJ_per_mol")}
    dg = result["dG0_kJ_per_mol"]
    print(f"analyze_run: dG0 = {dg:.2f} kJ/mol ({dg/4.184:.2f} kcal/mol), "
          f"bound cut-off {result['bound_cutoff_nm']:.2f} nm "
          f"({result['bound_cutoff_from']})")

    if not HAVE_MPL:
        return

    rr, ww = result["r_nm"], result["w_kJ_per_mol"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    axes[0].plot(r, f - f.min(), label="F(r) from sum_hills", color="tab:blue")
    axes[0].plot(rr, ww, label="w(r) = F + 2kT ln r, referenced", color="tab:red")
    axes[0].axhline(0.0, lw=0.6, color="k")
    axes[0].axvline(result["bound_cutoff_nm"], ls="--", lw=0.9, color="tab:green",
                    label=f"bound cut-off {result['bound_cutoff_nm']:.2f} nm")
    axes[0].axvspan(*result["plateau_window_nm"], color="grey", alpha=0.15,
                    label="unbound reference")
    axes[0].axvline(args.pmf_wall_nm, ls=":", lw=0.9, color="k", label="wall")
    axes[0].set_xlabel("COM-COM distance (nm)")
    axes[0].set_ylabel("free energy (kJ/mol)")
    axes[0].set_title("PMF: raw vs Jacobian-corrected")
    axes[0].legend(fontsize=7)

    cut, vals = zip(*result["dG0_vs_cutoff"])
    axes[1].plot(cut, vals, marker="o")
    axes[1].axvline(result["bound_cutoff_nm"], ls="--", lw=0.9, color="tab:green")
    axes[1].set_xlabel("bound-state cut-off (nm)")
    axes[1].set_ylabel("dG0 (kJ/mol)")
    axes[1].set_title("Sensitivity to the bound-state definition")
    fig.suptitle(f"Dimerization free energy: dG0 = {dg:.2f} kJ/mol "
                 f"({dg/4.184:.2f} kcal/mol) at 1 M")
    save(fig, out_dir, "20_dimerization_pmf", produced)
