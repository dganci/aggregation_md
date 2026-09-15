"""The feature vector cg_cvgen trains on, rebuilt from a COLVAR."""

from __future__ import annotations

import re

import numpy as np

from report_io import column

PERMUTABLE_BLOCKS = (
    ("sorted_d", re.compile(r"^d_\d+_\d+$")),
    ("sorted_cn", re.compile(r"^cn_\d+_\d+$")),
    ("sorted_rg", re.compile(r"^rg\d+$")),
)
INVARIANT_SCALARS = ("rg_com", "cn_total")


def permutation_invariant_features(fields, data, n_prot=None):
    """The feature matrix cg_cvgen trains on, from a COLVAR's columns."""
    if len(data) == 0:
        return [], np.empty((0, 0))
    names, columns = [], []
    for label, pattern in PERMUTABLE_BLOCKS:
        members = [f for f in fields if pattern.match(f)]
        if not members:
            continue
        if n_prot and label != "sorted_rg" and len(members) != n_prot * (n_prot - 1) // 2:
            raise ValueError(f"{len(members)} {label} columns for {n_prot} protomers; "
                             f"expected {n_prot * (n_prot - 1) // 2}")
        block = np.sort(np.column_stack([column(fields, data, m) for m in members]), axis=1)
        names.extend(f"{label}.{k + 1}" for k in range(block.shape[1]))
        columns.append(block)
    for name in INVARIANT_SCALARS:
        series = column(fields, data, name)
        if series is not None:
            names.append(name)
            columns.append(series[:, None])
    if not columns:
        return [], np.empty((len(data), 0))
    return names, np.hstack(columns).astype(np.float64)
