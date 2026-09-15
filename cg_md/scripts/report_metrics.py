"""The numbers behind the report: clustering, free energies, error bars."""
from __future__ import annotations

import math

import numpy as np

from aggregation_clusters import cluster_labels, cluster_sizes, pair_columns
from constants import KB_KJ_PER_MOL_K
from report_io import column


def contact_matrix(fields, data, n_prot):
    """Mean pairwise contact number per protomer pair -> (n_prot, n_prot)."""
    m = np.full((n_prot, n_prot), np.nan)
    for i in range(1, n_prot + 1):
        for j in range(i + 1, n_prot + 1):
            values = column(fields, data, f"cn_{i}_{j}")
            if values is None:
                continue
            m[i - 1, j - 1] = m[j - 1, i - 1] = float(np.mean(values))
    return m

MAX_CLUSTER_FRAMES = 200_000

def cluster_timeseries(fields, data, n_prot, threshold):
    """Largest-oligomer size and number of clusters per frame."""
    pairs, cols = pair_columns(fields, n_prot)
    if not cols or data.size == 0:
        return None, None, None

    step = max(1, len(data) // MAX_CLUSTER_FRAMES)
    frame_index = np.arange(0, len(data), step)
    contact = data[np.ix_(frame_index, cols)] >= threshold

    sizes = cluster_sizes(cluster_labels(contact, pairs, n_prot), n_prot)
    return sizes.max(axis=1), np.count_nonzero(sizes, axis=1), frame_index

def free_energy_2d(x, y, temperature_K, bins=60):
    """-kT ln P(x, y) from an *unbiased* trajectory, in kJ/mol, min set to 0."""
    hist, xe, ye = np.histogram2d(x, y, bins=bins)
    kT = KB_KJ_PER_MOL_K * temperature_K
    with np.errstate(divide="ignore"):
        fes = -kT * np.log(hist / hist.sum())
    fes -= np.nanmin(fes[np.isfinite(fes)]) if np.any(np.isfinite(fes)) else 0.0
    return fes, xe, ye

def reweighted_free_energy(values, rbias, temperature_K, bins=60):
    """1D FES of `values` from a well-tempered metadynamics run."""
    kT = KB_KJ_PER_MOL_K * temperature_K
    logw = rbias / kT
    logw -= logw.max()
    weights = np.exp(logw)
    hist, edges = np.histogram(values, bins=bins, weights=weights)
    with np.errstate(divide="ignore"):
        fes = -kT * np.log(hist / hist.sum())
    finite = np.isfinite(fes)
    if np.any(finite):
        fes -= fes[finite].min()
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, fes, weights

def block_average_error(series, min_blocks=5):
    """Block-averaging curve: standard error of the mean vs block length."""
    n = len(series)
    lengths, errors = [], []
    size = 1
    while n // size >= min_blocks:
        blocks = series[: (n // size) * size].reshape(-1, size).mean(axis=1)
        lengths.append(size)
        errors.append(float(blocks.std(ddof=1) / math.sqrt(len(blocks))))
        size *= 2
    return np.asarray(lengths), np.asarray(errors)

def autocorrelation_time(series):
    """Integrated autocorrelation time in frames (initial-positive-sequence)."""
    x = np.asarray(series, dtype=float)
    n = len(x)
    if n < 8:
        return float("nan")
    x = x - x.mean()
    if not np.any(x):
        return float("nan")

    size = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(x, size)
    acf = np.fft.irfft(spectrum * np.conjugate(spectrum), size)[:n]
    if acf[0] <= 0:
        return float("nan")
    acf /= acf[0]

    cut = np.argmax(acf[1:] <= 0.0) + 1 if np.any(acf[1:] <= 0.0) else n
    return float(0.5 + acf[1:cut].sum())
