#!/usr/bin/env bash
#
# Deep preflight for a batch manifest. Runs INSIDE the container.
#

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/manifest.sh
source "$HERE/lib/manifest.sh"
# shellcheck source=lib/rundir.sh
source "$HERE/lib/rundir.sh"

CG_MD="${CG_MD:-/data/build/cg_md}"
PROJECT_DIR="${PROJECT_DIR:-/data}"
NTOMP="${NTOMP:-1}"

MANIFEST="${1:-}"
[[ -n "$MANIFEST" && -f "$MANIFEST" ]] || {
  echo "usage: $0 <manifest.tsv>" >&2; exit 2; }
[[ -x "$CG_MD" ]] || { echo "cg_md not executable: $CG_MD" >&2; exit 2; }

manifest_load "$MANIFEST" || exit 2
names=("${MANIFEST_NAMES[@]}")
flagsets=("${MANIFEST_FLAGS[@]}")

echo "preflight: $MANIFEST  (${#names[@]} entries)  project: $PROJECT_DIR"

echo
echo "=== 1/2  configuration (--dry-run) ==="
rundir_resolve_all "$CG_MD" "$PROJECT_DIR" "$NTOMP" || exit 1

echo
echo "=== 2/2  structure and bead count (--stage cg) ==="
declare -A SEEN
cg_bad=0
for i in "${!names[@]}"; do
  name="${names[$i]}"
  flags_str="${flagsets[$i]}"
  protomer="$(manifest_flag "$flags_str" --protomer)"

  if [[ -z "$protomer" ]]; then
    printf '  FAIL  %-22s no --protomer in this entry\n' "$name"
    cg_bad=$((cg_bad + 1)); continue
  fi
  if [[ -n "${SEEN[$protomer]:-}" ]]; then
    printf '  ..    %-22s protomer %s already checked\n' "$name" "$protomer"
    continue
  fi
  SEEN["$protomer"]=1

  read -r -a flags <<< "$flags_str"
  if out=$("$CG_MD" --project-dir "$PROJECT_DIR" --ntomp "$NTOMP" "${flags[@]}" --stage cg 2>&1); then
    beads="$(grep -c '^ATOM' "$PROJECT_DIR/runs/${MANIFEST_RUNDIR[$name]}/prep/${protomer}_cg.pdb" 2>/dev/null || echo '?')"
    printf '  ok    %-22s protomer %-14s %s beads\n' "$name" "$protomer" "$beads"
  else
    printf '  FAIL  %-22s protomer %s\n' "$name" "$protomer"
    sed 's/^/          /' <<< "$(grep -E '^(ERROR|Coarse-grained|Input PDB|--)' <<< "$out" | head -8)"
    [[ -n "$(grep -E '^(ERROR|Coarse-grained|Input PDB)' <<< "$out")" ]] || \
      sed 's/^/          /' <<< "$(tail -6 <<< "$out")"
    cg_bad=$((cg_bad + 1))
  fi
done

echo
if [[ $cg_bad -ne 0 ]]; then
  echo "$cg_bad protomer(s) failed structure/bead checks; fix the manifest before running" >&2
  exit 1
fi
echo "preflight clean: ${#names[@]} entries, ${#SEEN[@]} distinct protomer(s)"
