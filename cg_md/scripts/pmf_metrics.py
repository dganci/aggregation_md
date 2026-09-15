"""The dimerization PMF: from sum_hills' F(r) to a standard-state dG."""

from __future__ import annotations

import math

import numpy as np

from constants import KB_KJ_PER_MOL_K, trapezoid


def fes_layout(fields: list[str]) -> tuple[list[str], int] | None:
    """(CV column names, index of the free-energy column) for a sum_hills file."""
    if not fields:
        return None
    free_idx = next((i for i, f in enumerate(fields) if f.endswith(".free")), None)
    if free_idx is None:
        cvs = [f for f in fields if not f.startswith("der_")]
        free_idx = len(cvs) - 1
        if free_idx <= 0:
            return None
    return fields[:free_idx], free_idx


V_STANDARD_NM3 = 1.66053906717

def dimerization_free_energy(r_nm, fes_kJ, temperature_K, wall_nm,
                             plateau_frac=0.75, bound_cutoff_nm=0.0):
    """Standard-state association free energy from a COM-distance PMF."""
    kT = KB_KJ_PER_MOL_K * temperature_K
    r = np.asarray(r_nm, dtype=float)
    f = np.asarray(fes_kJ, dtype=float)

    keep = (r > 1e-9) & np.isfinite(f)
    r, f = r[keep], f[keep]
    if len(r) < 10:
        return {"error": "PMF has too few usable points"}

    w = f + 2.0 * kT * np.log(r)

    lo, hi = plateau_frac * wall_nm, 0.95 * wall_nm
    plateau = (r >= lo) & (r <= hi)
    if plateau.sum() < 3:
        return {"error": f"no unbound plateau sampled in [{lo:.2f}, {hi:.2f}] nm - "
                         f"increase --pmf-wall-nm or lower --pmf-plateau-frac"}
    w = w - w[plateau].mean()
    plateau_rms = float(np.std(w[plateau]))

    integrand = 4.0 * np.pi * r ** 2 * np.exp(-w / kT)

    def dG(rb):
        sel = r <= rb
        if sel.sum() < 2:
            return np.nan, np.nan
        ka = float(trapezoid(integrand[sel], r[sel]))
        if ka <= 0:
            return np.nan, np.nan
        return ka, float(-kT * math.log(ka / V_STANDARD_NM3))

    if bound_cutoff_nm > 0:
        rb, how = bound_cutoff_nm, "user"
    else:
        i_min = int(np.argmin(w))
        rb, how = lo, "plateau onset (no barrier found)"
        for i in range(i_min + 1, len(w) - 1):
            if w[i] >= w[i - 1] and w[i] >= w[i + 1] and w[i] > w[i_min] + 0.1 * kT:
                rb, how = float(r[i]), "first barrier past the minimum"
                break

    ka, dg = dG(rb)
    scan = [(float(c), dG(c)[1]) for c in np.linspace(0.5 * rb, min(hi, 1.5 * rb), 9)]

    return {
        "r_nm": r, "w_kJ_per_mol": w,
        "kT_kJ_per_mol": float(kT),
        "well_depth_kJ_per_mol": float(w.min()),
        "r_min_nm": float(r[int(np.argmin(w))]),
        "bound_cutoff_nm": float(rb),
        "bound_cutoff_from": how,
        "Ka_nm3": ka,
        "dG0_kJ_per_mol": dg,
        "dG0_kcal_per_mol": dg / 4.184 if np.isfinite(dg) else None,
        "dG0_vs_cutoff": scan,
        "plateau_rms_kJ_per_mol": plateau_rms,
        "plateau_window_nm": [float(lo), float(hi)],
    }
