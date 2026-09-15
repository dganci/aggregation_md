"""NPMI and the other shape descriptors of the chain arrangement, from PLUMED's own COM-COM
distances rather than from the trajectory.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from aggregation_shape import shape_descriptors


def shape_descriptors_from_colvar(colvar_dir, n_prot):
    """NPMI of the chain arrangement, from PLUMED's own COM-COM distances."""
    files = sorted(Path(colvar_dir).glob("COLVAR_chunk_*.dat")) or \
            sorted(Path(colvar_dir).glob("COLVAR*.dat"))
    if not files:
        return None

    names, rows = None, []
    for f in files:
        for line in f.read_text().splitlines():
            if line.startswith("#! FIELDS"):
                names = line.split()[2:]
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if names and len(parts) == len(names):
                rows.append(dict(zip(names, map(float, parts))))
    if not rows:
        return None

    wanted = [f"d_{i}_{j}" for i in range(1, n_prot + 1) for j in range(i + 1, n_prot + 1)]
    if any(w not in rows[0] for w in wanted) or n_prot < 3:
        return None

    j_centre = np.eye(n_prot) - np.ones((n_prot, n_prot)) / n_prot
    out = np.empty(len(rows))
    for f, r in enumerate(rows):
        d2 = np.zeros((n_prot, n_prot))
        for i in range(n_prot):
            for j in range(i + 1, n_prot):
                d2[i, j] = d2[j, i] = r[f"d_{i+1}_{j+1}"] ** 2
        eig = np.linalg.eigvalsh(-0.5 * j_centre @ d2 @ j_centre)[::-1] / n_prot

        if eig[-1] < -1e-3 * max(abs(eig[0]), 1e-12):
            out[f] = np.nan
            continue
        shape = shape_descriptors(np.clip(eig[:3], 0.0, None))
        out[f] = shape.get("npmi", np.nan)
    return out
