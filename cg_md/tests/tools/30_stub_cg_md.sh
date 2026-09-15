# shellcheck shell=bash
STUB_DIR="$TMP/stub"
mkdir -p "$STUB_DIR"
cat > "$TMP/cg_md" <<'STUB'
#!/usr/bin/env bash
# Stub cg_md. Resolves a run directory for --dry-run; otherwise acts out the
# scenario recorded for this entry.
protomer=""; dry=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --protomer) protomer="$2"; shift 2 ;;
    --dry-run)  dry=1; shift ;;
    *) shift ;;
  esac
done
if [[ $dry -eq 1 ]]; then
  echo "  run  = 5x${protomer}_stubhash"
  exit 0
fi
behaviour="$(cat "$STUB_DIR/$protomer" 2>/dev/null || echo ok)"
case "$behaviour" in
  ok)      echo "stub: finished cleanly"; exit 0 ;;
  fail)    echo "stub: something went wrong"; exit 1 ;;
  fatal)   echo "ERROR: Input PDB not found: /nope.pdb"; exit 1 ;;
  # The 2026-08-03 scenario: the job runs, and while it runs the batch log
  # directory is removed from under it. Its own writes keep working, because a
  # file descriptor outlives its directory entry.
  nuke_logdir)
      echo "stub: running"
      rm -rf "$LOG_ROOT"
      echo "stub: still writing to an unlinked log"
      exit 1 ;;
esac
exit 0
STUB
chmod +x "$TMP/cg_md"

printf 'one\t--protomer one\ntwo\t--protomer two\nthree\t--protomer three\n' > "$TMP/batch.tsv"

RUN_OUT="$TMP/run_out.txt"
RUN_STATUS=0
run_batch() {
  local logroot="$1"; shift
  CG_MD="$TMP/cg_md" PROJECT_DIR="$TMP/project" LOG_ROOT="$logroot" \
    STUB_DIR="$STUB_DIR" NTOMP=1 RETRIES="${RETRIES:-0}" MIN_FREE_GB=0 \
    JOBS="${JOBS:-1}" \
    bash "$ROOT/tools/run_batch.sh" "$TMP/batch.tsv" "$@" > "$RUN_OUT" 2>&1
  RUN_STATUS=$?
}

mkdir -p "$TMP/project/runs"
