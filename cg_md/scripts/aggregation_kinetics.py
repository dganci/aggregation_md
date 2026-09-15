"""How long chains stay together, and when they first came together."""

from __future__ import annotations

import numpy as np

from constants import trapezoid


def contact_runs(series):
    """Lengths of the uninterrupted True runs in a boolean series."""
    h = np.asarray(series, dtype=bool)
    if h.size == 0 or not h.any():
        return np.empty(0, dtype=int), np.empty(0, dtype=bool)

    padded = np.concatenate(([False], h, [False]))
    edges = np.flatnonzero(np.diff(padded.astype(np.int8)))
    starts, ends = edges[0::2], edges[1::2]
    lengths = ends - starts

    censored = np.zeros(len(lengths), dtype=bool)
    if starts[0] == 0:
        censored[0] = True
    if ends[-1] == h.size:
        censored[-1] = True
    return lengths, censored


def contact_lifetimes(contact, dt_ps):
    """Mean uninterrupted contact duration, pooled over every protomer pair."""
    contact = np.asarray(contact, dtype=bool)
    if contact.ndim == 1:
        contact = contact[:, None]

    lengths, censored = [], 0
    for k in range(contact.shape[1]):
        run, cen = contact_runs(contact[:, k])
        lengths.append(run[~cen])
        censored += int(cen.sum())

    complete = np.concatenate(lengths) if lengths else np.empty(0, dtype=int)
    if complete.size == 0:
        return {"complete_runs": 0, "censored_runs": censored,
                "mean_lifetime_ps": None, "median_lifetime_ps": None,
                "max_lifetime_ps": None}
    return {
        "complete_runs": int(complete.size),
        "censored_runs": censored,
        "mean_lifetime_ps": float(complete.mean() * dt_ps),
        "median_lifetime_ps": float(np.median(complete) * dt_ps),
        "max_lifetime_ps": float(complete.max() * dt_ps),
    }


def contact_autocorrelation(contact, max_lag=None):
    """Intermittent contact correlation, pooled over protomer pairs."""
    contact = np.asarray(contact, dtype=float)
    if contact.ndim == 1:
        contact = contact[:, None]
    n, n_pairs = contact.shape
    if n < 4 or n_pairs == 0:
        return np.empty(0)

    lag = min(max_lag or n // 2, n - 1)
    size = 1 << (2 * n - 1).bit_length()
    pooled = np.zeros(lag + 1)
    for k in range(n_pairs):
        spectrum = np.fft.rfft(contact[:, k], size)
        raw = np.fft.irfft(spectrum * np.conjugate(spectrum), size)[: lag + 1]
        pooled += raw / (n - np.arange(lag + 1))
    pooled /= n_pairs

    occupancy = contact.mean(axis=0)
    mean_h = float(occupancy.mean())
    baseline = float((occupancy ** 2).mean())
    if mean_h <= 0.0:
        return np.zeros(lag + 1)
    if mean_h - baseline <= 1e-12:
        return np.ones(lag + 1)
    return (pooled - baseline) / (mean_h - baseline)


def correlation_time(curve, dt_ps):
    """Integral of a correlation function up to its first non-positive value."""
    c = np.asarray(curve, dtype=float)
    if c.size < 2 or c[0] <= 0:
        return float("nan")
    below = np.flatnonzero(c <= 0.0)
    cut = int(below[0]) if below.size else c.size
    if cut < 2:
        return float("nan")
    return float(trapezoid(c[:cut], dx=dt_ps))


def first_passage(largest, target, time_ps):
    """When the largest oligomer first reached `target` chains."""
    largest = np.asarray(largest)
    hit = np.flatnonzero(largest >= target)
    if hit.size == 0:
        return None
    return float(time_ps[hit[0]])
