#!/usr/bin/env python3
"""Derive --ss-string and --elastic-units for a protomer PDB, ready to paste into a cg_md
invocation.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

DSSP_TO_MARTINI = {
    "H": "H", "G": "H", "I": "H",
    "E": "E", "B": "E",
    "T": "C", "S": "C", "C": "C", " ": "C", "X": "C",
}
MARTINI_STRUCTURED = ("H", "E")


def pdb_residues(path) -> list:
    """[(resnum, chain)] of the first model, in file order, one entry per residue."""
    seen, out = set(), []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM") or line[16] not in (" ", "A"):
                continue
            key = (line[21], int(line[22:26]), line[26])
            if key not in seen:
                seen.add(key)
                out.append((int(line[22:26]), line[21]))
    return out


def run_mkdssp(pdb) -> str:
    exe = next((e for e in ("mkdssp", "dssp") if shutil.which(e)), None)
    if not exe:
        sys.exit("mkdssp/dssp not found on PATH. Run this inside the container, "
                 "or pass --dssp-file with a precomputed .dssp.")
    with tempfile.NamedTemporaryFile(suffix=".dssp", delete=False) as tmp:
        out = tmp.name
    try:
        proc = subprocess.run([exe, str(pdb), out], capture_output=True, text=True)
        if proc.returncode != 0:
            proc = subprocess.run([exe, "-i", str(pdb), "-o", out],
                                  capture_output=True, text=True)
        if proc.returncode != 0:
            sys.exit(f"{exe} failed:\n{proc.stderr.strip() or proc.stdout.strip()}")
        with open(out) as fh:
            return fh.read()
    finally:
        os.unlink(out)


def dssp_by_resnum(text) -> dict:
    """{resnum: raw DSSP code}."""
    out, started = {}, False
    for line in text.splitlines():
        if line.startswith("  #  RESIDUE"):
            started = True
            continue
        if not started or len(line) < 17 or line[13] == "!":
            continue
        num = line[5:10].strip()
        if num:
            out[int(num)] = line[16]
    return out


def martini_string(residues, raw) -> tuple:
    """(ss string aligned to `residues`, count of residues DSSP had nothing for)."""
    chars, missing = [], 0
    for num, _chain in residues:
        code = raw.get(num)
        if code is None:
            missing += 1
            chars.append("C")
        else:
            chars.append(DSSP_TO_MARTINI.get(code, "C"))
    return "".join(chars), missing


def structured_ranges(residues, ss, min_len) -> list:
    """Contiguous H/E stretches of at least min_len residues, as (lo, hi) resids."""
    out, start, prev = [], None, None
    for (num, _), code in zip(residues, ss):
        structured = code in MARTINI_STRUCTURED
        contiguous = prev is not None and num == prev + 1
        if structured and (start is None or not contiguous):
            if start is not None and prev - start + 1 >= min_len:
                out.append((start, prev))
            start = num
        elif not structured and start is not None:
            if prev - start + 1 >= min_len:
                out.append((start, prev))
            start = None
        prev = num
    if start is not None and prev - start + 1 >= min_len:
        out.append((start, prev))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdb", nargs="?", help="all-atom PDB (pdbs/AA/<protomer>.pdb)")
    ap.add_argument("--dssp-file", help="use this .dssp instead of running mkdssp")
    ap.add_argument("--seq-length", type=int,
                    help="expected residue count; checked against the PDB")
    ap.add_argument("--min-elastic", type=int, default=6,
                    help="shortest H/E stretch to restrain (default 6)")
    ap.add_argument("--no-elastic-suggestion", action="store_true")
    args = ap.parse_args()

    if not args.pdb and not args.dssp_file:
        ap.error("give a PDB, or --dssp-file together with --seq-length")

    if args.pdb:
        residues = pdb_residues(args.pdb)
    else:
        if not args.seq_length:
            ap.error("--dssp-file needs --seq-length")
        residues = [(i, "A") for i in range(1, args.seq_length + 1)]

    raw = dssp_by_resnum(args.dssp_file and open(args.dssp_file).read()
                         or run_mkdssp(args.pdb))
    ss, missing = martini_string(residues, raw)

    chains = sorted({c for _, c in residues})
    nums = [n for n, _ in residues]
    print(f"residues      {len(residues)}  (resid {min(nums)}-{max(nums)}, "
          f"chain{'s' if len(chains) > 1 else ''} {''.join(chains)})")
    if args.seq_length and args.seq_length != len(residues):
        print(f"WARNING --seq-length {args.seq_length} != {len(residues)} residues in "
              "the PDB. cg_md's check_inputs() will refuse this; fix --seq-length "
              "(and --ss-string, same length).")
    if missing:
        print(f"note          {missing} residue(s) absent from the DSSP output, "
              "written as C (coil)")
    counts = {c: ss.count(c) for c in ("H", "E", "C")}
    total = len(ss) or 1
    print("composition   " + "  ".join(
        f"{c} {counts[c]} ({100 * counts[c] / total:.0f}%)" for c in ("H", "E", "C")))
    print()
    print(f"--seq-length {len(residues)} --ss-string {ss}")

    if args.no_elastic_suggestion:
        return
    ranges = structured_ranges(residues, ss, args.min_elastic)
    print()
    if not ranges:
        print(f"no H/E stretch reaches {args.min_elastic} residues -> --no-elastic")
        print("  (an all-coil IDR has nothing to restrain; an elastic network here "
              "would freeze the disorder you are trying to sample)")
        return
    print("--elastic-units " + ",".join(f"{lo}:{hi}" for lo, hi in ranges))
    covered = sum(hi - lo + 1 for lo, hi in ranges)
    print(f"  {len(ranges)} structured stretch(es), {covered}/{len(residues)} residues "
          f"({100 * covered / total:.0f}%)")
    print("  resids are taken from THIS PDB. cg_md's default is 90:102, which is a "
          "leftover and almost certainly wrong for your construct -- always pass "
          "--elastic-units explicitly, or --no-elastic.")


if __name__ == "__main__":
    main()
