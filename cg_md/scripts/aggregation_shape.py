"""Shape of an assembly from its gyration-tensor eigenvalues."""

from __future__ import annotations

import numpy as np


def shape_descriptors(eigenvalues):
    """Asphericity, acylindricity, anisotropy and NPMI from a gyration tensor."""
    lam = np.sort(np.asarray(eigenvalues, dtype=float))[::-1]
    if lam.size < 3:
        return {}
    rg2 = float(lam.sum())
    if rg2 <= 0:
        return {}
    b = float(lam[0] - 0.5 * (lam[1] + lam[2]))
    c = float(lam[1] - lam[2])
    i_min, i_max = rg2 - lam[0], rg2 - lam[2]
    return {
        "rg_nm": float(np.sqrt(rg2)),
        "asphericity": b,
        "acylindricity": c,
        "anisotropy": float((b ** 2 + 0.75 * c ** 2) / rg2 ** 2),
        "npmi": float(i_min / i_max) if i_max > 0 else float("nan"),
    }
