#!/usr/bin/env bash
#
# Sequential batch runner for cg_md. Runs INSIDE the container.
#

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/manifest.sh
source "$HERE/lib/manifest.sh"
# shellcheck source=lib/rundir.sh
source "$HERE/lib/rundir.sh"
# shellcheck source=lib/logdir.sh
source "$HERE/lib/logdir.sh"
# shellcheck source=lib/resources.sh
source "$HERE/lib/resources.sh"
# shellcheck source=lib/entry.sh
source "$HERE/lib/gpu.sh"
source "$HERE/lib/entry.sh"

CG_MD="${CG_MD:-/data/build/cg_md}"
PROJECT_DIR="${PROJECT_DIR:-/data}"
LOG_ROOT="${LOG_ROOT:-/data/runs/_batch}"
PER_SIM_TIMEOUT="${PER_SIM_TIMEOUT:-0}"
# GPUs to spread the concurrent entries over. Empty (the default) means no
# pinning at all, which is right on a CPU-only machine.
#
# It matters because GROMACS, left alone, has every mdrun pick the first device
# it can see: four concurrent entries on a four-GPU node would all land on GPU 0
# and leave three idle while the accounting bills the whole node either way
# (B_H = T x N x R x C, R = 1.0, C = 32 - you pay for the node, not for the GPUs
# you used). Pinning one GPU per entry is the difference between 1/4 and 4/4 of
# what the allocation is charged for.
#
# Give a count ("4") or an explicit list ("0,1,2,3"). Entries take a free device
# under flock rather than by index: the pool starts a new entry whenever ANY slot
# frees, so an index-based assignment collides as soon as entries finish out of
# order - which they do, since they have different sizes and stop conditions.
GPUS="${GPUS:-}"
MIN_FREE_GB="${MIN_FREE_GB:-40}"
NTOMP="${NTOMP:-0}"
JOBS="${JOBS:-2}"
RETRIES="${RETRIES:-1}"

FATAL_PATTERNS=(
  "This run was started with"
  "Input PDB not found"
  "residue(s), but --seq-length"
  "but --phospho was not set"
  "Coarse-grained bead count mismatch"
  "Protein bead count mismatch"
  "must be <="
  "Invalid --"
  "unknown option"
)

MANIFEST="${1:-}"; shift || true
MODE="run"; ONLY=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --check) MODE="check"; shift ;;
    --only)  ONLY="${2:-}"; shift 2 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$MANIFEST" && -f "$MANIFEST" ]] || {
  echo "usage: $0 <manifest.tsv> [--check] [--only NAME]" >&2; exit 2; }
[[ -x "$CG_MD" ]] || { echo "cg_md not executable: $CG_MD" >&2; exit 2; }

logdir_ensure "$LOG_ROOT" || {
  echo "FATAL: cannot use $LOG_ROOT as the batch log directory." >&2
  echo "Check that the bind mount is present and writable from inside the container." >&2
  exit 3
}

SUMMARY="$LOG_ROOT/summary.tsv"
STATUS="$LOG_ROOT/status.txt"
[[ -f "$SUMMARY" ]] || printf 'name\trun_dir\tstatus\texit\tattempts\tstarted\tended\tminutes\tlog\n' > "$SUMMARY"

exec 9>"$LOG_ROOT/.batch.lock"
flock -n 9 || { echo "another run_batch.sh is already running (lock in $LOG_ROOT)" >&2; exit 2; }

say() {
  echo "$*"
  printf '%s  %s\n' "$(date -Is)" "$*" >> "$STATUS" 2>/dev/null || true
}
resources_detect
manifest_load "$MANIFEST" || exit 2
names=("${MANIFEST_NAMES[@]}")
flagsets=("${MANIFEST_FLAGS[@]}")

if [[ -n "$ONLY" ]]; then
  found=0
  for n in "${names[@]}"; do [[ "$n" == "$ONLY" ]] && found=1; done
  [[ $found -eq 1 ]] || { echo "--only $ONLY: no such entry in $MANIFEST" >&2; exit 2; }
fi

echo "manifest: $MANIFEST  (${#names[@]} entries)  project: $PROJECT_DIR"

echo
echo "=== preflight (--dry-run) ==="

if ! timeout --signal=TERM "$PER_SIM_TIMEOUT" true 2>/dev/null; then
  echo "ABORT: PER_SIM_TIMEOUT='$PER_SIM_TIMEOUT' is not a duration timeout(1) accepts." >&2
  echo "       Use 0 for no limit (the default), or a value like 48h, 90m, 3600." >&2
  exit 1
fi
if [[ "$PER_SIM_TIMEOUT" == "0" ]]; then
  echo "per-entry wall clock: none (PER_SIM_TIMEOUT=0); the adaptive stop condition"
  echo "                      and --adaptive-max-total-us decide when an entry ends"
else
  echo "per-entry wall clock: $PER_SIM_TIMEOUT  - an entry killed mid-chunk loses that"
  echo "                      whole chunk, because a SIGTERM checkpoint has no .gro"
fi
rundir_resolve_all "$CG_MD" "$PROJECT_DIR" "$NTOMP" || {
  echo "nothing was run" >&2; exit 1; }

if [[ "$MODE" == "check" ]]; then
  echo
  echo "configuration is clean; --check, stopping here"
  echo "for the structure and bead-count checks too, run:  ./tools/preflight.sh $MANIFEST"
  exit 0
fi
avail="$(df -Pk "$PROJECT_DIR" 2>/dev/null | awk 'NR==2 {printf "%d", $4/1048576}')"
if [[ -n "$avail" && "$avail" -lt "$MIN_FREE_GB" ]]; then
  say "ABORT: only ${avail} GB free under $PROJECT_DIR (need $MIN_FREE_GB)"
  exit 1
fi

if [[ "$JOBS" -gt 1 ]]; then
  dup="$(for i in "${!names[@]}"; do
           [[ -n "$ONLY" && "${names[$i]}" != "$ONLY" ]] && continue
           sed -n 's/.*--protomer \([^ ]*\).*/\1/p' <<< "${flagsets[$i]}"
         done | sort | uniq -d)"
  if [[ -n "$dup" ]]; then
    echo "ABORT: JOBS=$JOBS but these --protomer values appear in more than one entry:" >&2
    echo "$dup" | sed 's/^/         /' >&2
    echo "       Such entries share runs/relax_<protomer>/ and the .itp martinize2 writes" >&2
    echo "       into the project directory, so running them at the same time would have" >&2
    echo "       them overwrite each other. Re-run with JOBS=1." >&2
    exit 2
  fi
fi

# Il marcatore di avvio, e non e' cosmetico: e' l'unica prova, in status.txt,
# che il batch e' arrivato a partire. Viveva in lib/entry.sh, dove pero' girava
# al momento del `source` - prima che say() fosse definita, quattro righe piu'
# sotto - quindi falliva con "say: command not found" e la riga non veniva MAI
# scritta. Che e' esattamente la firma del 2026-08-03: nessuna traccia nei log
# di un batch che credi partito. Va emessa qui, dopo il preflight e il controllo
# dello spazio, cioe' quando la partenza e' un fatto e non un'intenzione.
echo
say "=== batch start $(date -Is) ==="

OUTCOME_DIR="$(mktemp -d)"
trap 'rm -rf "$OUTCOME_DIR"' EXIT

selected=()
for i in "${!names[@]}"; do
  [[ -n "$ONLY" && "${names[$i]}" != "$ONLY" ]] && continue
  selected+=("$i")
done

if [[ "$JOBS" -le 1 ]]; then
  for i in "${selected[@]}"; do run_entry "$i" || true; done
else
  if (sleep 0 & wait -n) >/dev/null 2>&1; then have_wait_n=1; else have_wait_n=0; fi

  running=0
  for i in "${selected[@]}"; do
    if [[ "$running" -ge "$JOBS" ]]; then
      if [[ "$have_wait_n" -eq 1 ]]; then
        wait -n 2>/dev/null || true
        running=$((running - 1))
      else
        wait
        running=0
      fi
    fi
    run_entry "$i" &
    running=$((running + 1))
  done
  wait
fi

ok=0; failed=0; skipped=0
for i in "${selected[@]}"; do
  case "$(cat "$OUTCOME_DIR/${names[$i]}" 2>/dev/null)" in
    done)            ok=$((ok + 1)) ;;
    skipped)         skipped=$((skipped + 1)) ;;
    timeout|failed)  failed=$((failed + 1)) ;;
    *)               failed=$((failed + 1)) ;;
  esac
done

echo
say "=== batch end $(date -Is) ==="
say "done $ok   failed $failed   skipped $skipped"
echo "summary: $SUMMARY"
[[ $failed -eq 0 ]]
