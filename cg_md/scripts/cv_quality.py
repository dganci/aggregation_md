"""Quality control for a trained collective variable, after the fact."""
from __future__ import annotations

import numpy as np

from cv_features import (INVARIANT_SCALARS, PERMUTABLE_BLOCKS,  # noqa: F401
                         permutation_invariant_features)
from cv_stats import basin_transitions, component_independence, correlation, rank_data  # noqa: F401
from offline_io import read_csv_rows, read_json  # noqa: F401


def training_coverage(produced, train_min, train_max, train_std, names=None):
    """How far the production CV strayed outside the range it was trained on."""
    produced = np.atleast_2d(np.asarray(produced, dtype=float))
    lo_all = np.asarray(train_min, dtype=float).ravel()
    hi_all = np.asarray(train_max, dtype=float).ravel()
    sd_all = np.asarray(train_std, dtype=float).ravel()
    if produced.size == 0 or len(lo_all) != produced.shape[1] or len(hi_all) != produced.shape[1]:
        return []

    out = []
    for k in range(produced.shape[1]):
        lo, hi = float(lo_all[k]), float(hi_all[k])
        spread = float(sd_all[k]) if k < len(sd_all) else 0.0
        values = produced[:, k]
        excess = np.maximum(lo - values, values - hi)
        out.append({
            "component": names[k] if names and k < len(names) else f"cv{k}",
            "train_min": lo, "train_max": hi,
            "run_min": float(values.min()), "run_max": float(values.max()),
            "fraction_outside": float(np.mean(excess > 0.0)),
            "max_excess_in_train_sd": float(excess.max() / spread) if spread > 0 else float("nan"),
        })
    return out

def grid_coverage(produced, grid_min, grid_max, names=None):
    """Where the run sat inside the METAD grid, and whether it left it."""
    produced = np.atleast_2d(np.asarray(produced, dtype=float))
    grid_min = np.asarray(grid_min, dtype=float).ravel()
    grid_max = np.asarray(grid_max, dtype=float).ravel()
    if produced.size == 0 or len(grid_min) != produced.shape[1]:
        return []

    out = []
    for k in range(produced.shape[1]):
        values = produced[:, k]
        span = grid_max[k] - grid_min[k]
        used = values.max() - values.min()
        out.append({
            "component": names[k] if names and k < len(names) else f"cv{k}",
            "grid_min": float(grid_min[k]), "grid_max": float(grid_max[k]),
            "fraction_outside_grid": float(np.mean((values < grid_min[k]) | (values > grid_max[k]))),
            "fraction_of_grid_used": float(used / span) if span > 0 else float("nan"),
        })
    return out


def plateau_summary(its_rows, selected_lag=None):
    """What the implied-timescale scan says about the lag that was chosen."""
    if not its_rows:
        return {"error": "no its_scan.csv"}
    band = [r for r in its_rows if r.get("in_plateau")]
    chosen = [r for r in its_rows if r.get("selected")]
    out = {
        "lags_scanned": len(its_rows),
        "plateau_lags": [int(r["lag_frames"]) for r in band],
        "selected_lags": [int(r["lag_frames"]) for r in chosen],
        "lag_frames": [int(r["lag_frames"]) for r in its_rows],
        "its1_ns": [r["its1_ps"] / 1000.0 for r in its_rows],
        "lag_ns": [r["lag_ps"] / 1000.0 for r in its_rows],
    }
    if not band:
        out["warning"] = ("no implied-timescale plateau: the leading timescale never stopped "
                          "depending on the lag, so the trajectory is too short for its own "
                          "slowest mode")
    if selected_lag is not None:
        out["selected_lag"] = int(selected_lag)
        out["selected_lag_in_plateau"] = int(selected_lag) in out["plateau_lags"]
    return out

def lag_score_summary(rows, selected_lag=None):
    """The scoring table reduced to the numbers that decide whether to trust it."""
    if not rows:
        return {"error": "no lag_scores.csv"}
    table = []
    for row in rows:
        table.append({
            "lag": int(row.get("lag", 0)),
            "score": row.get("score"),
            "valid": bool(row.get("valid") in (1.0, "True", "true", True)),
            "vamp2_train": row.get("vamp2_train"),
            "vamp2_test": row.get("vamp2_test"),
            "overfit_gap": row.get("vamp_gap"),
            "ck_spectral_error": row.get("ck_spectral_error"),
            "spectral_gap": row.get("spectral_gap"),
        })
    out = {"lags": table, "any_valid": any(r["valid"] for r in table)}
    if not out["any_valid"]:
        out["warning"] = ("no lag passed cg_cvgen's validity checks, so the selection fell back "
                          "to raw VAMP-2 - see the warning in the training log")
    if selected_lag is not None:
        out["selected"] = next((r for r in table if r["lag"] == int(selected_lag)), None)
    return out

def sigma_vs_spread(params):
    """Whether any hill width had to be floored, from cv_params.json."""
    sigma = np.asarray(params.get("sigma", []), dtype=float)
    std = np.asarray(params.get("std", []), dtype=float)
    if sigma.size == 0 or sigma.size != std.size:
        return {}
    floored = [int(k) for k in range(len(sigma)) if sigma[k] > 0.5 * std[k] + 1e-9]
    return {
        "sigma": sigma.tolist(),
        "half_std": (0.5 * std).tolist(),
        "floored_components": floored,
    }
