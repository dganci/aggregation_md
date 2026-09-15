"""Oligomers from the pairwise contact numbers: who is bound to whom, and how big the pieces are.
"""

from __future__ import annotations

import numpy as np

from constants import KB_KJ_PER_MOL_K


def pair_columns(fields, n_prot):
    """Indices of the cn_i_j columns present, and the protomer pair each is."""
    pairs, cols = [], []
    for i in range(1, n_prot + 1):
        for j in range(i + 1, n_prot + 1):
            name = f"cn_{i}_{j}"
            if name in fields:
                pairs.append((i - 1, j - 1))
                cols.append(fields.index(name))
    return pairs, cols


def cluster_labels(contact, pairs, n_prot):
    """Connected-component label of every protomer, frame by frame."""
    contact = np.asarray(contact, dtype=bool)
    labels = np.empty((len(contact), n_prot), dtype=np.int32)

    for f in range(len(contact)):
        parent = list(range(n_prot))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for k, active in enumerate(contact[f]):
            if not active:
                continue
            a, b = find(pairs[k][0]), find(pairs[k][1])
            if a != b:
                parent[a] = b

        for node in range(n_prot):
            labels[f, node] = find(node)

    return labels


def cluster_sizes(labels, n_prot):
    """Size of the oligomer each label stands for, frame by frame."""
    labels = np.asarray(labels)
    n_frames = len(labels)
    if n_frames == 0:
        return np.empty((0, n_prot), dtype=int)
    flat = labels + np.arange(n_frames)[:, None] * n_prot
    return np.bincount(flat.ravel(), minlength=n_frames * n_prot).reshape(n_frames, n_prot)


def size_histogram(sizes, n_prot):
    """Total number of oligomers of each size over the whole trajectory."""
    counts = np.bincount(np.asarray(sizes).ravel(), minlength=n_prot + 1).astype(float)
    counts[0] = 0.0
    return counts[: n_prot + 1]


def size_distribution(counts):
    """Number- and mass-weighted oligomer size distributions, and their means."""
    counts = np.asarray(counts, dtype=float)
    sizes = np.arange(len(counts), dtype=float)
    n_clusters = counts.sum()
    n_chains = (sizes * counts).sum()
    if n_clusters <= 0 or n_chains <= 0:
        return {"error": "no oligomers counted"}
    return {
        "sizes": sizes[1:].astype(int).tolist(),
        "counts": counts[1:].tolist(),
        "number_fraction": (counts[1:] / n_clusters).tolist(),
        "mass_fraction": (sizes[1:] * counts[1:] / n_chains).tolist(),
        "mean_size_number": float(n_chains / n_clusters),
        "mean_size_weight": float((sizes ** 2 * counts).sum() / n_chains),
    }


def size_free_energy(counts, temperature_K):
    """-kT ln( N_s / N_1 ), the oligomer size distribution on an energy axis."""
    counts = np.asarray(counts, dtype=float)
    kT = KB_KJ_PER_MOL_K * temperature_K
    if len(counts) < 2 or counts[1] <= 0:
        return np.full(max(len(counts) - 1, 0), np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        g = -kT * np.log(counts[1:] / counts[1])
    return np.where(counts[1:] > 0, g, np.nan)

def censoring_report(largest, n_prot):
    """How often the largest oligomer was every chain in the box."""
    largest = np.asarray(largest)
    if largest.size == 0:
        return {}
    return {
        "n_prot": int(n_prot),
        "fraction_at_n_prot": float(np.mean(largest >= n_prot)),
        "max_observed": int(largest.max()),
    }
