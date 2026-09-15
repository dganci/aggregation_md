#!/usr/bin/env python3
"""Inter-chain residue-residue contact map and shape descriptors."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from shape_from_colvar import shape_descriptors_from_colvar

try:
    import mdtraj as md
    HAVE_MDTRAJ = True
except Exception:
    HAVE_MDTRAJ = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False


def chain_residue_slices(traj, n_prot):
    """Residue index ranges of the n_prot protomers."""
    top = traj.topology
    protein = [r.index for r in top.residues if r.is_protein]
    if not protein:
        protein = [r.index for r in top.residues]

    chains = [[r.index for r in c.residues if r.index in set(protein)] for c in top.chains]
    chains = [c for c in chains if c]
    if len(chains) == n_prot:
        return chains

    if len(protein) % n_prot:
        raise ValueError(f"{len(protein)} protein residues do not divide into {n_prot} protomers; "
                         "check --n-prot against the trajectory")
    per = len(protein) // n_prot
    return [protein[i * per:(i + 1) * per] for i in range(n_prot)]


def residue_contact_map(traj, slices, cutoff_nm, stride):
    """Probability that residue i of one chain contacts residue j of another."""
    per = len(slices[0])
    if any(len(s) != per for s in slices):
        raise ValueError("protomers have different residue counts; a shared contact map is undefined")

    frames = range(0, traj.n_frames, max(1, stride))
    sub = traj[list(frames)]

    pairs, index = [], []
    for a in range(len(slices)):
        for b in range(a + 1, len(slices)):
            for i, ri in enumerate(slices[a]):
                for j, rj in enumerate(slices[b]):
                    pairs.append((ri, rj))
                    index.append((i, j))
    if not pairs:
        return None, 0

    dist, _ = md.compute_contacts(sub, contacts=np.asarray(pairs), scheme="closest", periodic=True)
    hit = dist < cutoff_nm

    m = np.zeros((per, per))
    for k, (i, j) in enumerate(index):
        c = hit[:, k].sum()
        m[i, j] += c
        m[j, i] += c

    n_chain_pairs = len(slices) * (len(slices) - 1) // 2
    m /= float(sub.n_frames * n_chain_pairs * 2)
    return m, sub.n_frames


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--traj", required=True, type=Path)
    ap.add_argument("--top", required=True, type=Path)
    ap.add_argument("--n-prot", required=True, type=int)
    ap.add_argument("--cutoff-nm", type=float, default=0.6,
                    help="contact cut-off; 0.6 nm suits Martini beads, 0.4 nm is the "
                         "all-atom value used in the reference protocol")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--label", default="")
    ap.add_argument("--colvar-dir", type=Path, default=None,
                    help="where to look for the COLVARs the NPMI is computed from "
                         "(default: next to the trajectory)")
    args = ap.parse_args()

    if not HAVE_MDTRAJ:
        print("contact_map: mdtraj not available - skipping the residue contact map", file=sys.stderr)
        return 0
    for p in (args.traj, args.top):
        if not p.is_file():
            print(f"contact_map: missing {p} - skipping", file=sys.stderr)
            return 0

    args.out_dir.mkdir(parents=True, exist_ok=True)
    traj = md.load(str(args.traj), top=str(args.top))
    slices = chain_residue_slices(traj, args.n_prot)

    summary = {"frames": int(traj.n_frames), "n_prot": args.n_prot,
               "cutoff_nm": args.cutoff_nm, "residues_per_protomer": len(slices[0])}

    m, used = residue_contact_map(traj, slices, args.cutoff_nm, args.stride)
    if m is not None:
        np.savetxt(args.out_dir / "residue_contact_map.dat", m, fmt="%.6f")
        summary["frames_analysed"] = int(used)
        flat = [(m[i, j], i + 1, j + 1) for i in range(m.shape[0]) for j in range(i, m.shape[1])]
        flat.sort(reverse=True)
        summary["top_contacts"] = [{"res_i": i, "res_j": j, "probability": round(float(p), 4)}
                                   for p, i, j in flat[:25] if p > 0]

        if HAVE_MPL:
            fig, ax = plt.subplots(figsize=(6.5, 5.5))
            im = ax.imshow(m, origin="lower", cmap="viridis",
                           extent=(0.5, m.shape[0] + 0.5, 0.5, m.shape[1] + 0.5))
            fig.colorbar(im, ax=ax, label="contact probability")
            ax.set_xlabel("residue (chain i)")
            ax.set_ylabel("residue (chain j)")
            ax.set_title(f"Inter-chain residue contact map{' - ' + args.label if args.label else ''}")
            fig.tight_layout()
            fig.savefig(args.out_dir / "30_residue_contact_map.png", dpi=150)
            plt.close(fig)

    try:
        npmi = shape_descriptors_from_colvar(args.colvar_dir or args.traj.parent, args.n_prot)
        if npmi is None:
            raise RuntimeError("no COLVAR with the d_i_j columns next to the trajectory; "
                               "NPMI is skipped rather than re-derived from coordinates, "
                               "which is not well defined under periodic boundaries")
        usable = int(np.count_nonzero(~np.isnan(npmi)))
        summary["npmi"] = {
            "mean": float(np.nanmean(npmi)) if usable else None,
            "std": float(np.nanstd(npmi)) if usable else None,
            "min": float(np.nanmin(npmi)) if usable else None,
            "max": float(np.nanmax(npmi)) if usable else None,
            "frames_used": usable,
            "frames_total": int(len(npmi)),
            "frames_not_embeddable": int(len(npmi) - usable),
        }
        if usable < len(npmi):
            print(f"note: {len(npmi) - usable} of {len(npmi)} frames "
                  f"({100 * (len(npmi) - usable) / len(npmi):.1f}%) have a non-embeddable "
                  f"COM distance matrix - the chains are more than half a box apart and the "
                  f"assembly has no shape to report. Those frames are NaN.")
        np.savetxt(args.out_dir / "npmi.dat", npmi, fmt="%.6f")
        if HAVE_MPL:
            fig, axes = plt.subplots(1, 2, figsize=(11, 4))
            axes[0].plot(npmi, lw=0.7)
            axes[0].set_xlabel("analysed frame")
            axes[0].set_ylabel("NPMI (I_min / I_max)")
            axes[0].set_title("Aggregate shape over time")
            axes[1].hist(npmi[~np.isnan(npmi)], bins=40)
            axes[1].set_xlabel("NPMI")
            axes[1].set_ylabel("count")
            axes[1].set_title("1 = spherical, small = elongated")
            fig.tight_layout()
            fig.savefig(args.out_dir / "31_shape_npmi.png", dpi=150)
            plt.close(fig)
    except Exception as e:
        summary["npmi_error"] = str(e)

    (args.out_dir / "contact_map_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"contact_map: wrote residue_contact_map.dat and summary to {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
