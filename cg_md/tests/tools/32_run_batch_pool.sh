# shellcheck shell=bash
probe_dir="$TMP/concurrent_probe"
mkdir -p "$probe_dir"
probe_pids=(); probe_fails=0
for _ in 1 2 3 4 5 6 7 8; do ( logdir_ensure "$probe_dir" 2>/dev/null ) & probe_pids+=("$!"); done
for p in "${probe_pids[@]}"; do wait "$p" || probe_fails=$((probe_fails + 1)); done
check "logdir_ensure survives eight concurrent callers" "$probe_fails" "0"
check "concurrent probes clean up after themselves" \
      "$(find "$probe_dir" -name '.probe.*' | wc -l | tr -d ' ')" "0"

rm -rf "$TMP/project/runs"; mkdir -p "$TMP/project/runs"
echo ok > "$STUB_DIR/one"; echo ok > "$STUB_DIR/two"; echo ok > "$STUB_DIR/three"
JOBS=2 run_batch "$TMP/par1"

check "the pool exits 0 when every entry succeeds" "$RUN_STATUS" "0"
missing=""
for e in one two three; do
  grep -q "stub: finished cleanly" "$TMP/par1/$e.log" 2>/dev/null || missing="$missing $e"
done
if [[ -z "$missing" ]]; then ok "every entry runs under JOBS=2"
else bad "every entry runs under JOBS=2" "no log for:$missing"; fi

check "the summary has exactly one row per entry" \
      "$(awk -F'\t' '$1=="one"||$1=="two"||$1=="three"' "$TMP/par1/summary.tsv" | wc -l | tr -d ' ')" "3"
check "the tally counts every success" \
      "$(grep -c 'done 3   failed 0   skipped 0' "$RUN_OUT")" "1"

rm -rf "$TMP/project/runs"; mkdir -p "$TMP/project/runs"
echo fail > "$STUB_DIR/two"
JOBS=2 run_batch "$TMP/par2"
check "a failure inside the pool is counted" \
      "$(grep -c 'done 2   failed 1   skipped 0' "$RUN_OUT")" "1"
check "the pool reports non-zero when an entry fails" "$RUN_STATUS" "1"
echo ok > "$STUB_DIR/two"

JOBS=2 run_batch "$TMP/par3"
check "completed entries are skipped on a concurrent re-run" \
      "$(grep -c 'done 1   failed 0   skipped 2' "$RUN_OUT")" "1"

rm -rf "$TMP/project/runs"; mkdir -p "$TMP/project/runs"
echo nuke_logdir > "$STUB_DIR/one"; echo ok > "$STUB_DIR/two"; echo ok > "$STUB_DIR/three"
run_batch "$TMP/lr4"

if grep -q "no longer readable" "$RUN_OUT"; then
  ok "a vanished log directory is diagnosed, not blamed on the entry"
else
  bad "a vanished log directory is diagnosed, not blamed on the entry" "$(tail -6 "$RUN_OUT")"
fi
if grep -q "stub: finished cleanly" "$TMP/lr4/two.log" 2>/dev/null; then
  ok "entries after the incident actually run, instead of failing instantly"
else
  bad "entries after the incident actually run, instead of failing instantly" \
      "two.log: $(cat "$TMP/lr4/two.log" 2>/dev/null || echo MISSING)"
fi
if [[ -f "$TMP/lr4/summary.tsv" ]] && grep -q "^two" "$TMP/lr4/summary.tsv"; then
  ok "the summary survives the log directory being removed"
else
  bad "the summary survives the log directory being removed" "no row for 'two'"
fi

rm -rf "$TMP/project/runs"; mkdir -p "$TMP/project/runs"
if [[ "$(id -u)" == "0" ]]; then
  echo "[SKIP] an unwritable log root aborts with status 3 (running as root)"
else
  mkdir -p "$TMP/lr5_parent"; chmod a-w "$TMP/lr5_parent"
  run_batch "$TMP/lr5_parent/logs"
  check "an unwritable log root aborts with status 3" "$RUN_STATUS" "3"
  chmod u+w "$TMP/lr5_parent"
fi
