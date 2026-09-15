"""Reading PLUMED files onto one increasing clock."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


class Clock:
    """Puts PLUMED's restarting clock back onto one increasing axis."""

    __slots__ = ("last", "spacing", "offset")

    def __init__(self):
        self.last = None
        self.spacing = 0.0
        self.offset = 0.0

    def rebase(self, t: float) -> float | None:
        if self.last is not None:
            if t + self.offset < self.last - 0.5 * self.spacing:
                self.offset = (self.last + self.spacing) - t
        t += self.offset
        if self.last is not None:
            if t <= self.last:
                return None
            if self.spacing <= 0.0:
                self.spacing = t - self.last
        self.last = t
        return t

def read_plumed(path: Path, clock: "Clock | None" = None) -> tuple[list[str], np.ndarray]:
    """Reads a PLUMED COLVAR/HILLS/FES file into (field names, data array)."""
    if not path.is_file():
        return [], np.empty((0, 0))

    fields: list[str] = []
    rows: list[list[float]] = []
    timed = False
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#!"):
                parts = line.split()
                if len(parts) >= 3 and parts[1] == "FIELDS":
                    if fields and parts[2:] != fields:
                        print(f"read_plumed: {path.name} changes columns part "
                              f"way through; keeping the {len(rows)} rows "
                              "before that", file=sys.stderr)
                        break
                    fields = parts[2:]
                    timed = bool(fields) and fields[0] == "time"
                    if timed and clock is None:
                        clock = Clock()
                continue
            if line.startswith("#"):
                continue
            parts = line.split()
            if fields and len(parts) != len(fields):
                continue
            try:
                row = [float(x) for x in parts]
            except ValueError:
                continue
            if timed:
                t = clock.rebase(row[0])
                if t is None:
                    continue
                row[0] = t
            rows.append(row)

    if not fields or not rows:
        return fields, np.empty((0, len(fields)))
    return fields, np.asarray(rows, dtype=float)

def load_colvars(paths: list[Path]) -> tuple[list[str], np.ndarray]:
    """Concatenates COLVAR files onto ONE continuous clock."""
    fields: list[str] = []
    chunks: list[np.ndarray] = []
    clock = Clock()

    for path in paths:
        f, data = read_plumed(path, clock)
        if not data.size:
            continue
        fields = f or fields
        chunks.append(data)

    if not chunks:
        return fields, np.empty((0, len(fields)))
    return fields, np.vstack(chunks)

def sampled_time_us(fields: list[str], data: np.ndarray) -> float | None:
    """Simulated time the table spans, in us."""
    if "time" not in fields or data.size == 0:
        return None
    t = data[:, fields.index("time")]
    if len(t) < 2:
        return 0.0
    return float((t[-1] - t[0] + np.median(np.diff(t))) / 1e6)
