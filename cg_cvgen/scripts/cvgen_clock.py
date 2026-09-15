"""Putting PLUMED's restarting clock onto one axis."""
import numpy as np


def rebase_clocks(times):
    """Puts a sequence of PLUMED time arrays onto one increasing clock."""
    out, last_end, spacing = [], None, 0.0

    for raw in times:
        t = np.asarray(raw, dtype=float).copy()
        if t.size == 0:
            out.append(t)
            continue

        gaps = np.diff(t)
        positive = gaps[gaps > 0.0]
        if positive.size:
            spacing = float(np.median(positive))

        breaks = np.flatnonzero(gaps < -0.5 * spacing) + 1 if t.size > 1 else np.empty(0, dtype=int)
        bounds = np.concatenate(([0], breaks, [t.size])).astype(int)

        for a, b in zip(bounds[:-1], bounds[1:]):
            if last_end is not None and t[a] < last_end - 0.5 * spacing:
                t[a:b] += (last_end + spacing) - t[a]
            last_end = float(t[b - 1])

        out.append(t)

    return out
