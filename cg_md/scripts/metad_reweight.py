"""c(t) reweighting of a well-tempered run: weights, effective sample size, and the block error on
a reweighted free energy.
"""

from __future__ import annotations

import numpy as np

from constants import KB_KJ_PER_MOL_K


def effective_sample_size(log_weights):
    """Kish effective sample size of a reweighted trajectory."""
    log_weights = np.asarray(log_weights, dtype=float)
    log_weights = log_weights[np.isfinite(log_weights)]
    if log_weights.size == 0:
        return {"error": "no finite weights"}

    w = np.exp(log_weights - log_weights.max())
    total = w.sum()
    ess = float(total ** 2 / (w ** 2).sum())
    return {
        "frames": int(len(w)),
        "ess": ess,
        "ess_fraction": ess / len(w),
        "max_weight_fraction": float(w.max() / total),
        "log_weight_span": float(log_weights.max() - log_weights.min()),
    }


def log_weights_from_rbias(rbias, temperature_K):
    """exp(+rbias/kT) as log-weights, the c(t) reweighting of a WT run."""
    return np.asarray(rbias, dtype=float) / (KB_KJ_PER_MOL_K * temperature_K)


def block_ranges(n_samples, n_blocks, segments=None):
    """Index ranges of contiguous blocks, never straddling a segment boundary."""
    if not segments:
        segments = [(0, int(n_samples))]
    out = []
    for start, stop in segments:
        if stop - start < n_blocks:
            continue
        edges = np.linspace(start, stop, n_blocks + 1).astype(int)
        out.extend(zip(edges[:-1], edges[1:]))
    return [(a, b) for a, b in out if b > a]


def block_free_energy(values, log_weights, temperature_K, bins=50, n_blocks=10,
                      value_range=None, segments=None):
    """Reweighted F(s) with an error bar, by block averaging over the weights."""
    values = np.asarray(values, dtype=float)
    log_weights = np.asarray(log_weights, dtype=float)
    if values.size < 2 * n_blocks or n_blocks < 2 or values.shape != log_weights.shape:
        return np.empty(0), np.empty(0), np.empty(0)

    bad = ~(np.isfinite(values) & np.isfinite(log_weights))
    finite_lw = log_weights[~bad]
    if finite_lw.size == 0:
        return np.empty(0), np.empty(0), np.empty(0)
    top = finite_lw.max()
    weights = np.where(bad, 0.0, np.exp(np.where(bad, top, log_weights) - top))
    values = np.where(bad, value_range[0] if value_range else 0.0, values)

    if value_range is None:
        good = values[~bad]
        value_range = (float(good.min()), float(good.max()))
    edges = np.linspace(value_range[0], value_range[1], bins + 1)
    width = edges[1] - edges[0]
    if width <= 0:
        return np.empty(0), np.empty(0), np.empty(0)

    hist, block_w = [], []
    for a, b in block_ranges(len(values), n_blocks, segments):
        counts, _ = np.histogram(values[a:b], bins=edges, weights=weights[a:b])
        total = weights[a:b].sum()
        if total <= 0:
            continue
        hist.append(counts / (total * width))
        block_w.append(total)
    if len(hist) < 2:
        return np.empty(0), np.empty(0), np.empty(0)

    hist = np.asarray(hist)
    block_w = np.asarray(block_w, dtype=float)
    w_sum = block_w.sum()
    w_sq = (block_w ** 2).sum()

    ave = (block_w[:, None] * hist).sum(axis=0) / w_sum
    ave2 = (block_w[:, None] * hist ** 2).sum(axis=0) / w_sum
    denom = w_sum ** 2 - w_sq
    var = (ave2 - ave ** 2) * (w_sum ** 2 / denom) if denom > 0 else np.zeros_like(ave)
    err = np.sqrt(np.maximum(var, 0.0) / len(hist))

    kT = KB_KJ_PER_MOL_K * temperature_K
    with np.errstate(divide="ignore", invalid="ignore"):
        fes = -kT * np.log(ave)
        sigma = kT * err / ave
    finite = np.isfinite(fes)
    if np.any(finite):
        fes = fes - fes[finite].min()
    fes = np.where(finite, fes, np.nan)
    sigma = np.where(finite, sigma, np.nan)
    return 0.5 * (edges[:-1] + edges[1:]), fes, sigma


def running_free_energy_difference(values, log_weights, split, temperature_K, n_points=20):
    """Free-energy difference between the two sides of `split`, versus time."""
    values = np.asarray(values, dtype=float)
    log_weights = np.asarray(log_weights, dtype=float)
    keep = np.isfinite(values) & np.isfinite(log_weights)
    values, log_weights = values[keep], log_weights[keep]
    if values.size < 2 * n_points:
        return np.empty(0), np.empty(0)

    weights = np.exp(log_weights - log_weights.max())
    left = np.cumsum(np.where(values < split, weights, 0.0))
    right = np.cumsum(np.where(values >= split, weights, 0.0))

    kT = KB_KJ_PER_MOL_K * temperature_K
    idx = np.linspace(len(values) // n_points, len(values), n_points).astype(int) - 1
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = -kT * np.log(left[idx] / right[idx])
    return idx + 1, np.where(np.isfinite(delta), delta, np.nan)
