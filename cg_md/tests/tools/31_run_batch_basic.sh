# shellcheck shell=bash
echo ok > "$STUB_DIR/one"; echo ok > "$STUB_DIR/two"; echo ok > "$STUB_DIR/three"
run_batch "$TMP/lr1"
check "a clean batch exits 0" "$RUN_STATUS" "0"
if grep -q "done 3   failed 0" "$RUN_OUT"; then ok "a clean batch reports 3 done"
else bad "a clean batch reports 3 done" "$(tail -3 "$RUN_OUT")"; fi
check "every entry gets its own log" "$(ls "$TMP/lr1"/*.log | wc -l | tr -d ' ')" "3"
check "the summary has one row per entry" "$(( $(wc -l < "$TMP/lr1/summary.tsv") - 1 ))" "3"

run_batch "$TMP/lr1"
check "re-running a finished batch exits 0" "$RUN_STATUS" "0"
if grep -q "skipped 3" "$RUN_OUT"; then ok "re-running skips entries already done"
else bad "re-running skips entries already done" "$(tail -3 "$RUN_OUT")"; fi

rm -rf "$TMP/project/runs"; mkdir -p "$TMP/project/runs"
echo fail > "$STUB_DIR/two"
run_batch "$TMP/lr2"
check "a batch with a failure exits 1" "$RUN_STATUS" "1"
if grep -q "done 2   failed 1" "$RUN_OUT"; then ok "a failing entry does not stop the batch"
else bad "a failing entry does not stop the batch" "$(tail -3 "$RUN_OUT")"; fi

rm -rf "$TMP/project/runs"; mkdir -p "$TMP/project/runs"
echo fatal > "$STUB_DIR/two"
RETRIES=2 run_batch "$TMP/lr3"
if grep -q "not retryable" "$RUN_OUT"; then ok "a configuration error is not retried"
else bad "a configuration error is not retried" "$(tail -5 "$RUN_OUT")"; fi
check "the attempt count records a single try" \
      "$(awk -F'\t' '$1=="two" {print $5}' "$TMP/lr3/summary.tsv")" "1"

cat > "$TMP/cg_md_seeded" <<'STUB'
#!/usr/bin/env bash
protomer=""; seed=""; dry=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --protomer) protomer="$2"; shift 2 ;;
    --seed)     seed="$2"; shift 2 ;;
    --dry-run)  dry=1; shift ;;
    *) shift ;;
  esac
done
[[ $dry -eq 1 ]] && { echo "  run  = 5x${protomer}_${seed}"; exit 0; }
echo "stub: finished cleanly"; exit 0
STUB
chmod +x "$TMP/cg_md_seeded"
printf 'p1\t--protomer shared --seed 1\np2\t--protomer shared --seed 2\np3\t--protomer other --seed 3\n' \
  > "$TMP/dup.tsv"

CG_MD="$TMP/cg_md_seeded" PROJECT_DIR="$TMP/project" LOG_ROOT="$TMP/dup1" \
  STUB_DIR="$STUB_DIR" NTOMP=1 RETRIES=0 MIN_FREE_GB=0 JOBS=2 \
  bash "$ROOT/tools/run_batch.sh" "$TMP/dup.tsv" > "$RUN_OUT" 2>&1
check "JOBS=2 refuses entries that share a --protomer" "$?" "2"
if grep -q "appear in more than one entry" "$RUN_OUT"; then
  ok "the refusal names the offending --protomer"
else
  bad "the refusal names the offending --protomer" "$(tail -4 "$RUN_OUT")"
fi

CG_MD="$TMP/cg_md_seeded" PROJECT_DIR="$TMP/project" LOG_ROOT="$TMP/dup2" \
  STUB_DIR="$STUB_DIR" NTOMP=1 RETRIES=0 MIN_FREE_GB=0 JOBS=1 \
  bash "$ROOT/tools/run_batch.sh" "$TMP/dup.tsv" > "$RUN_OUT" 2>&1
check "JOBS=1 runs the same manifest without complaint" "$?" "0"
