#!/usr/bin/env python3
"""Writes a variant of martini_v3.0.0.itp with the protein-water interactions scaled by a factor
lambda_PW.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MARTINI = ROOT / "martini_v300"
BEAD_LIST = ROOT / "tools" / "martini3_protein_beads.txt"
WATER = "W"


def known_protein_beads() -> set:
    """Every bead type martinize2 -ff martini3IDP emits for these constructs."""
    if not BEAD_LIST.exists():
        return set()
    return {tok for line in BEAD_LIST.read_text().splitlines()
            if not line.startswith("#") for tok in line.split()}


def protein_bead_types(itp_paths) -> set:
    """Bead types that appear in the [ atoms ] section of the protomer topologies."""
    beads = set()
    for p in itp_paths:
        section = None
        for line in p.read_text().splitlines():
            s = line.strip()
            if s.startswith("["):
                section = s.strip("[] ").strip()
                continue
            if section == "atoms" and s and not s.startswith(";"):
                parts = s.split()
                if len(parts) >= 2:
                    beads.add(parts[1])
    return beads


def rescale(src: Path, dst: Path, lam: float, beads: set) -> int:
    out, in_nb, n = [], False, 0
    for line in src.read_text().split("\n"):
        s = line.strip()
        if s.startswith("["):
            in_nb = s.startswith("[ nonbond_params ]")
        if in_nb and s and not s.startswith((";", "[")):
            p = s.split()
            if len(p) >= 5 and WATER in (p[0], p[1]) and ({p[0], p[1]} - {WATER}) & beads:
                line = (f"{p[0]:>6}{p[1]:>6}  {p[2]} "
                        f"{float(p[3]):.6e} {float(p[4]) * lam:>15.6e}")
                n += 1
        out.append(line)
    header = (f"; interazioni proteina-acqua riscalate di lambda_PW = {lam}\n"
              f"; generato da tools/make_lambda_itp.py da {src.name}\n"
              f"; {n} coppie (bead proteico, {WATER}) modificate; acqua-acqua e ioni intatti\n"
              f"; bead proteici coperti: {' '.join(sorted(beads))}\n")
    dst.write_text(header + "\n".join(out))
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("lambda_pw", type=float, help="e.g. 1.10")
    ap.add_argument("--itp", nargs="*", help="topologies to read protein bead types from "
                                             "(on top of the shipped complete list)")
    ap.add_argument("--only-found", action="store_true",
                    help="use ONLY the beads seen in the topologies, ignoring "
                         f"{BEAD_LIST.name} (not recommended: see known_protein_beads)")
    args = ap.parse_args()

    if not 0.5 <= args.lambda_pw <= 2.0:
        sys.exit("[make_lambda_itp] lambda must lie in [0.5, 2.0]")

    src = MARTINI / "martini_v3.0.0.itp"
    if not src.exists():
        sys.exit(f"[make_lambda_itp] missing {src}")

    paths = [Path(x) for x in (args.itp or [])] + sorted((ROOT / "runs").rglob("*_0.itp"))
    if not paths and args.only_found:
        sys.exit("[make_lambda_itp] no protomer topology found: pass --itp, or make "
                 "a structure first with `--stage cg`")

    beads = protein_bead_types(paths)
    if not args.only_found:
        beads |= known_protein_beads()
    if not beads:
        sys.exit("[make_lambda_itp] no bead type read from the topologies")

    tag = int(round(args.lambda_pw * 100))
    dst = MARTINI / f"martini_v3.0.0_lpw{tag}.itp"
    n = rescale(src, dst, args.lambda_pw, beads)

    print(f"protein bead types : {len(beads)}  ({' '.join(sorted(beads))})")
    print(f"pairs rescaled     : {n}  (epsilon x {args.lambda_pw})")
    print(f"written            : {dst}")
    print(f"\nuse it with:  --martini-lambda-pw {args.lambda_pw}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
