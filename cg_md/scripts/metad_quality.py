"""Convergence and reweighting diagnostics for a well-tempered metadynamics run."""
from __future__ import annotations

import numpy as np

from metad_reweight import (block_free_energy, block_ranges, effective_sample_size,  # noqa: F401
                            log_weights_from_rbias, running_free_energy_difference)


def hill_decay(time_ps, height, tail_fraction=0.05):
    """How far the deposition rate fell, and when it halved."""
    time_ps = np.asarray(time_ps, dtype=float)
    height = np.asarray(height, dtype=float)
    keep = np.isfinite(time_ps) & np.isfinite(height)
    time_ps, height = time_ps[keep], height[keep]
    if height.size < 10:
        return {"error": "fewer than ten hills"}

    n_tail = max(int(tail_fraction * len(height)), 1)
    initial = float(np.median(height[:n_tail]))
    final = float(np.median(height[-n_tail:]))

    halved = None
    if initial > 0:
        below = np.flatnonzero(height <= 0.5 * initial)
        if below.size:
            halved = float(time_ps[below[0]])

    return {
        "hills": int(len(height)),
        "initial_height_kJ_per_mol": initial,
        "final_height_kJ_per_mol": final,
        "final_over_initial": float(final / initial) if initial > 0 else float("nan"),
        "half_height_time_ps": halved,
        "deposition_span_ps": float(time_ps[-1] - time_ps[0]) if len(time_ps) > 1 else 0.0,
    }

def fes_rmsd(surfaces, window_kJ=30.0):
    """RMSD of each surface to the last one, over the region that is populated."""
    usable = [np.asarray(s, dtype=float) for s in surfaces if np.size(s)]
    if len(usable) < 2:
        return {"error": "fewer than two comparable surfaces"}
    reference = usable[-1]
    same = [s for s in usable if s.shape == reference.shape]
    skipped = len(usable) - len(same)

    finite = np.isfinite(reference)
    if not np.any(finite):
        return {"error": "the reference surface has no finite bins"}
    mask = finite & (reference <= reference[finite].min() + window_kJ)
    if mask.sum() < 5:
        return {"error": f"fewer than five bins within {window_kJ} kJ/mol of the minimum"}

    rmsd = []
    for s in same:
        both = mask & np.isfinite(s)
        if both.sum() < 5:
            rmsd.append(float("nan"))
            continue
        diff = s[both] - reference[both]
        rmsd.append(float(np.sqrt(np.mean((diff - diff.mean()) ** 2))))

    return {
        "rmsd_kJ_per_mol": rmsd,
        "window_kJ_per_mol": float(window_kJ),
        "bins_compared": int(mask.sum()),
        "surfaces_skipped": skipped,
        "final_rmsd_kJ_per_mol": rmsd[-2] if len(rmsd) > 1 else float("nan"),
    }

def ct_curve(bias, rbias):
    """c(t) = bias - rbias, and whether it behaves the way it has to."""
    bias = np.asarray(bias, dtype=float)
    rbias = np.asarray(rbias, dtype=float)
    if bias.size == 0 or bias.shape != rbias.shape:
        return {}
    ct = bias - rbias
    steps = np.diff(ct)
    tail = max(len(ct) // 10, 2)
    return {
        "ct_final_kJ_per_mol": float(ct[-1]),
        "ct_decreasing_fraction": float(np.mean(steps < 0)) if steps.size else 0.0,
        "ct_tail_growth_kJ_per_mol": float(ct[-1] - ct[-tail]),
        "curve": ct,
    }

def walker_balance(hills_per_walker):
    """Whether every walker actually contributed hills."""
    counts = np.asarray(list(hills_per_walker.values()), dtype=float)
    if counts.size == 0:
        return {}
    return {
        "per_walker": {str(k): int(v) for k, v in hills_per_walker.items()},
        "walkers": int(counts.size),
        "total_hills": int(counts.sum()),
        "min_over_max": float(counts.min() / counts.max()) if counts.max() > 0 else 0.0,
        "empty_walkers": [str(k) for k, v in hills_per_walker.items() if v == 0],
    }
