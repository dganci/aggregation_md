#!/usr/bin/env python3
"""Report non-bonded bead pairs that are too close in a coarse-grained protomer."""
import argparse
import itertools
import math
import sys


def read_beads(path) -> list:
    """[(resnum, atom name, xyz)] of the first model."""
    out = []
    for line in open(path):
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")) and line[16] in (" ", "A"):
            out.append((int(line[22:26]), line[12:16].strip(),
                        (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return out


def clashes(beads, cutoff, skip_neighbours) -> list:
    """Non-bonded pairs closer than cutoff, nearest first."""
    grid = {}
    for i, (_, _, x) in enumerate(beads):
        grid.setdefault(tuple(int(c // cutoff) for c in x), []).append(i)

    found = []
    for i, (ri, ni, xi) in enumerate(beads):
        key = tuple(int(c // cutoff) for c in xi)
        for d in itertools.product((-1, 0, 1), repeat=3):
            for j in grid.get((key[0] + d[0], key[1] + d[1], key[2] + d[2]), ()):
                if j <= i:
                    continue
                rj, nj, xj = beads[j]
                if abs(ri - rj) <= skip_neighbours:
                    continue
                dist = math.dist(xi, xj)
                if dist < cutoff:
                    found.append((dist, ri, ni, rj, nj))
    return sorted(found)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pdb", nargs="+")
    ap.add_argument("--cutoff", type=float, default=3.0, help="angstrom (default 3.0)")
    ap.add_argument("--skip-neighbours", type=int, default=1,
                    help="ignore pairs within N residues of each other (default 1)")
    ap.add_argument("--chain-size", type=int, default=0,
                    help="residues per chain, to label inter-chain clashes")
    ap.add_argument("--max-print", type=int, default=10)
    args = ap.parse_args()

    worst = 0
    for path in args.pdb:
        beads = read_beads(path)
        found = clashes(beads, args.cutoff, args.skip_neighbours)
        print(f"\n{path}\n  {len(beads)} beads, {len(found)} non-bonded pair(s) "
              f"< {args.cutoff:g} A")
        for dist, ri, ni, rj, nj in found[:args.max_print]:
            tag = ""
            if args.chain_size:
                ca, cb = (ri - 1) // args.chain_size, (rj - 1) // args.chain_size
                la, lb = ri - ca * args.chain_size, rj - cb * args.chain_size
                tag = (f"   [chain {chr(65+ca)} {la} <-> chain {chr(65+cb)} {lb}]"
                       + ("  INTER-CHAIN" if ca != cb else ""))
            print(f"    {dist:5.2f} A   res {ri} {ni}  <->  res {rj} {nj}{tag}")
        if len(found) > args.max_print:
            print(f"    ... and {len(found) - args.max_print} more")
        if found:
            print(f"  worst {found[0][0]:.2f} A"
                  + ("   <- severe, expect a huge initial Fmax" if found[0][0] < 2.0 else ""))
            worst = max(worst, len(found))
    return 1 if worst else 0


if __name__ == "__main__":
    sys.exit(main())
