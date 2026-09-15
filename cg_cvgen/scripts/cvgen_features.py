"""The feature vector a CV is trained on, and the check that it is the whole one."""
import json
import re
import sys
from pathlib import Path

import numpy as np

PERMUTABLE_BLOCKS = (
    ("sorted_d", re.compile(r"^d_\d+_\d+$")),
    ("sorted_cn", re.compile(r"^cn_\d+_\d+$")),
    ("sorted_rg", re.compile(r"^rg\d+$")),
)
INVARIANT_SCALARS = ("rg_com", "cn_total")

def sort_permutable_blocks(df, feature_cols):
    """Replaces each permutable block by its row-wise ascending sort."""
    out = df.copy()
    cols, present = [], set(feature_cols)

    for label, pattern in PERMUTABLE_BLOCKS:
        members = [c for c in feature_cols if pattern.match(c)]
        if not members:
            continue
        values = np.sort(df[members].to_numpy(dtype=np.float64), axis=1)
        names = [f"{label}.{i + 1}" for i in range(values.shape[1])]
        for k, name in enumerate(names):
            out[name] = values[:, k]
        cols.extend(names)

    for name in INVARIANT_SCALARS:
        if name in present:
            cols.append(name)

    unclassified = [c for c in feature_cols
                    if c not in INVARIANT_SCALARS
                    and not any(p.match(c) for _, p in PERMUTABLE_BLOCKS)]
    if unclassified:
        raise ValueError(
            "Cannot make these features permutation-invariant because they belong to no known "
            f"block: {unclassified}. Either extend PERMUTABLE_BLOCKS or pass "
            "--no-permutation-invariant (accepting that the CV will distinguish states that "
            "differ only by how the chains happen to be numbered).")

    print(f"permutation-invariant features: {len(cols)} sorted columns")
    return out, cols


def check_against_cg_md_schema(paths, feature_cols):
    """Refuses to train on fewer features than cg_md said it emitted."""
    for colvar in paths:
        schema = Path(colvar).parent / "feature_schema.json"
        if not schema.exists():
            continue
        try:
            declared = json.loads(schema.read_text()).get("feature_cols") or []
        except (OSError, ValueError) as exc:
            print(f"warning: could not read {schema} ({exc}); "
                  f"feature selection is unchecked", file=sys.stderr)
            return
        if not declared:
            return
        missed = [c for c in declared if c not in set(feature_cols)]
        if missed:
            raise ValueError(
                f"feature_regex selected {len(feature_cols)} of the {len(declared)} "
                f"columns cg_md declared in {schema}. Missing: {missed}.\n"
                "Training would silently proceed on the smaller vector, and the resulting CV "
                "could not be fed the vector cg_md rebuilds at biasing time. Update "
                "--feature-regex to cover cg_md's current naming, or pass the columns you "
                "actually want explicitly.")
        return
