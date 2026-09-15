"""Paths, time stamps and JSON for the offline toolkit."""

from __future__ import annotations

import csv
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np

CV_PREFIX = "mycv.node-"


def natural_key(text):
    """Sort key that puts walker10 after walker2 instead of after walker1."""
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", str(text))]


def expand(patterns):
    """Existing files named by a list of paths or shell globs, in walker order."""
    found = []
    for item in patterns or []:
        item = str(item)
        if any(c in item for c in "*?["):
            found.extend(glob.glob(item))
        elif Path(item).is_file():
            found.append(item)
        else:
            print(f"analyze_offline: no such file: {item}", file=sys.stderr)
    return sorted((Path(m) for m in found if Path(m).is_file()), key=natural_key)


def frame_spacing(time_ps, fallback):
    """Median spacing of a time column, or `fallback` when there is no column."""
    if time_ps is None or len(time_ps) < 2:
        return fallback
    steps = np.diff(time_ps)
    steps = steps[steps > 0]
    return float(np.median(steps)) if steps.size else fallback


def join_on_time(time_a, time_b):
    """Indices of the rows two PLUMED tables have in common."""
    key_a = np.round(np.asarray(time_a, dtype=float), 6)
    key_b = np.round(np.asarray(time_b, dtype=float), 6)
    _, index_a, index_b = np.intersect1d(key_a, key_b, return_indices=True)
    return index_a, index_b


def cv_columns(fields):
    """Names of the CV components in a biased COLVAR, in node order."""
    return sorted((f for f in fields if f.startswith(CV_PREFIX)), key=natural_key)


def plain(obj):
    """A structure json.dumps can write as valid JSON."""
    if isinstance(obj, dict):
        return {str(k): plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return plain(obj.tolist())
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        obj = float(obj)
    if isinstance(obj, float):
        return obj if np.isfinite(obj) else None
    if isinstance(obj, (str, int, bool)) or obj is None:
        return obj
    return str(obj)


def dump(summary, out_dir, name):
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / name
    path.write_text(json.dumps(plain(summary), indent=2))
    print(f"analyze_offline: wrote {path}")


def read_csv_rows(path):
    """Rows of a CSV as dicts with numbers parsed, or [] if it is not there."""
    path = Path(path)
    if not path.is_file():
        return []
    rows = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            parsed = {}
            for key, value in row.items():
                try:
                    parsed[key] = float(value)
                except (TypeError, ValueError):
                    parsed[key] = value
            rows.append(parsed)
    return rows

def read_json(path):
    path = Path(path)
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return {}
