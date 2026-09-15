"""Rank correlation, component independence and recrossings - the statistics the CV quality checks
are made of.
"""

from __future__ import annotations

import numpy as np


def rank_data(values):
    """Ranks with tied values averaged, i.e. scipy.stats.rankdata without scipy."""
    values = np.asarray(values, dtype=float)
    _, inverse, counts = np.unique(values, return_inverse=True, return_counts=True)
    starts = np.concatenate(([0], np.cumsum(counts)[:-1]))
    return (starts + 0.5 * (counts - 1))[inverse]

def correlation(x, y):
    """Pearson and Spearman coefficients, NaN when either series is constant."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 3:
        return {"pearson": float("nan"), "spearman": float("nan"), "n": int(keep.sum())}
    x, y = x[keep], y[keep]

    def pearson(a, b):
        a = a - a.mean()
        b = b - b.mean()
        denom = np.sqrt((a ** 2).sum() * (b ** 2).sum())
        return float(a @ b / denom) if denom > 0 else float("nan")

    return {"pearson": pearson(x, y),
            "spearman": pearson(rank_data(x), rank_data(y)),
            "n": int(len(x))}

def component_independence(values):
    """Largest off-diagonal correlation between CV components."""
    values = np.atleast_2d(np.asarray(values, dtype=float))
    if values.shape[1] < 2 or len(values) < 3:
        return {}
    corr = np.corrcoef(values.T)
    off = np.abs(corr - np.diag(np.diag(corr)))
    i, j = np.unravel_index(np.argmax(off), off.shape)
    return {"max_abs_offdiagonal": float(off[i, j]), "between": [int(i), int(j)]}

def basin_transitions(values, lo=None, hi=None):
    """How often the CV travelled from one end of its range to the other."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 3:
        return {"crossings": 0, "round_trips": 0}
    if lo is None:
        lo = float(np.percentile(values, 30))
    if hi is None:
        hi = float(np.percentile(values, 70))
    if hi <= lo:
        return {"crossings": 0, "round_trips": 0, "low": lo, "high": hi}

    state = np.zeros(len(values), dtype=np.int8)
    state[values <= lo] = -1
    state[values >= hi] = 1
    visited = state[state != 0]
    if visited.size < 2:
        return {"crossings": 0, "round_trips": 0, "low": lo, "high": hi}
    crossings = int(np.count_nonzero(visited[1:] != visited[:-1]))
    return {"crossings": crossings, "round_trips": crossings // 2, "low": lo, "high": hi}
