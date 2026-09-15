#!/usr/bin/env bash
# Regenerates tools/martini3_protein_beads.txt by coarse-graining every
# construct in pdbs/AA/ and taking the union of the bead types. Run it again
# when a construct is added: see the comment inside the file it writes.
#
#   tools/collect_bead_types.sh [docker-image]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${1:-dganci/aggregation-md:alma9-metatomic}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

cat > "$WORK/go.sh" <<'INNER'
#!/bin/bash
set -e
cd /work
for d in force_fields mappings scripts; do ln -sfn /proj/$d $d; done
for f in /proj/pdbs/AA/*_1.pdb; do
  c=$(basename "$f" .pdb)
  r=$(awk '/^ATOM/{print substr($0,22,6)}' "$f" | sort -u | wc -l)
  martinize2 -f "$f" -x "cg_$c.pdb" -o "top_$c.top" -name "$c" \
             -ff martini3IDP -ss "$(python3 -c "print('C'*$r)")" >"m_$c.log" 2>&1 \
    || { echo "martinize2 failed on $c; see m_$c.log" >&2; exit 1; }
done
INNER
chmod +x "$WORK/go.sh"

echo "coarse-graining every construct in pdbs/AA/ ..."
docker run --rm -v "$WORK:/work" -v "$ROOT:/proj" "$IMAGE" /work/go.sh

python3 - "$WORK" "$ROOT/tools/martini3_protein_beads.txt" <<'PY'
import sys, glob, os
def beads(p):
    sec, out = None, set()
    for line in open(p):
        s = line.strip()
        if s.startswith("["):
            sec = s.strip("[] ").strip(); continue
        if sec == "atoms" and s and not s.startswith(";"):
            f = s.split()
            if len(f) >= 2: out.add(f[1])
    return out
per = {os.path.basename(f).replace("_0.itp", ""): beads(f)
       for f in sorted(glob.glob(os.path.join(sys.argv[1], "*_0.itp")))}
if not per: sys.exit("no topology was produced")
union = sorted(set().union(*per.values()))
for k, v in per.items(): print(f"  {k:<12} {len(v):>2} types")
print(f"  {'union':<12} {len(union):>2} types")
dst = sys.argv[2]
head = [l for l in open(dst).read().split("\n") if l.startswith("#") or not l.strip()]
open(dst, "w").write("\n".join(head).rstrip("\n") + "\n\n" + " ".join(union) + "\n")
print(f"\nwritten {dst}")
PY
