"""Reading what the pipeline writes, and putting figures on disk."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:  # pragma: no cover - environment dependent
    HAVE_MPL = False
    plt = None

from plumed_clock import Clock, load_colvars, read_plumed, sampled_time_us  # noqa: E402,F401


def read_jsonl(path: Path) -> list[dict]:
    if not path.exists():
        return []
    out = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return out

def column(fields: list[str], data: np.ndarray, name: str) -> np.ndarray | None:
    if name not in fields or data.size == 0:
        return None
    return data[:, fields.index(name)]

MAX_PLOT_POINTS = 50_000

def thin(*arrays):
    """Uniformly subsamples parallel arrays for plotting."""
    first = next((a for a in arrays if a is not None), None)
    if first is None or len(first) <= MAX_PLOT_POINTS:
        return arrays
    step = len(first) // MAX_PLOT_POINTS + 1
    return tuple(None if a is None else a[::step] for a in arrays)

def save(fig, out_dir: Path, name: str, produced: list[str]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{name}.png"
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    produced.append(path.name)
