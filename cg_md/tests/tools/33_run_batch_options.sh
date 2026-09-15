# shellcheck shell=bash
rm -rf "$TMP/project/runs"; mkdir -p "$TMP/project/runs"
echo ok > "$STUB_DIR/one"
run_batch "$TMP/lr6" --only one
check "--only exits 0" "$RUN_STATUS" "0"
check "--only runs a single entry" "$(ls "$TMP/lr6"/*.log | wc -l | tr -d ' ')" "1"

run_batch "$TMP/lr7" --only nosuchentry
check "--only rejects an unknown entry name" "$RUN_STATUS" "2"

rm -rf "$TMP/project/runs"; mkdir -p "$TMP/project/runs"
run_batch "$TMP/lr8" --check
check "--check exits 0 on a valid manifest" "$RUN_STATUS" "0"
check "--check runs nothing" "$(ls "$TMP/lr8"/*.log 2>/dev/null | wc -l | tr -d ' ')" "0"

printf 'a\t--protomer p1\n' > "$TMP/one.tsv"

CG_MD="$TMP/cg_md" PROJECT_DIR="$TMP/project" LOG_ROOT="$TMP/lrT1" \
  STUB_DIR="$STUB_DIR" NTOMP=1 MIN_FREE_GB=0 \
  bash "$ROOT/tools/run_batch.sh" "$TMP/one.tsv" > "$RUN_OUT" 2>&1
status=$?
check "the default PER_SIM_TIMEOUT runs a batch to completion" "$status" "0"
if grep -q "per-entry wall clock: none" "$RUN_OUT"; then ok "no wall clock by default, and it says so"
else bad "no wall clock by default, and it says so" "$(grep -i 'wall clock' "$RUN_OUT" | head -2)"; fi

CG_MD="$TMP/cg_md" PROJECT_DIR="$TMP/project" LOG_ROOT="$TMP/lrT2" \
  STUB_DIR="$STUB_DIR" NTOMP=1 MIN_FREE_GB=0 PER_SIM_TIMEOUT=30m \
  bash "$ROOT/tools/run_batch.sh" "$TMP/one.tsv" > "$RUN_OUT" 2>&1
status=$?
check "an explicit PER_SIM_TIMEOUT is honoured" "$status" "0"
if grep -q "per-entry wall clock: 30m" "$RUN_OUT"; then ok "an explicit wall clock is announced"
else bad "an explicit wall clock is announced" "$(grep -i 'wall clock' "$RUN_OUT" | head -2)"; fi

CG_MD="$TMP/cg_md" PROJECT_DIR="$TMP/project" LOG_ROOT="$TMP/lrT3" \
  STUB_DIR="$STUB_DIR" NTOMP=1 MIN_FREE_GB=0 PER_SIM_TIMEOUT=120hh \
  bash "$ROOT/tools/run_batch.sh" "$TMP/one.tsv" > "$RUN_OUT" 2>&1
status=$?
check "an invalid PER_SIM_TIMEOUT aborts the batch" "$status" "1"
check "an invalid PER_SIM_TIMEOUT runs no entry" \
  "$(ls "$TMP/lrT3"/*.log 2>/dev/null | wc -l | tr -d ' ')" "0"

printf 'a\t--protomer same\nb\t--protomer same\n' > "$TMP/collide.tsv"
CG_MD="$TMP/cg_md" PROJECT_DIR="$TMP/project" LOG_ROOT="$TMP/lr9" \
  STUB_DIR="$STUB_DIR" NTOMP=1 MIN_FREE_GB=0 \
  bash "$ROOT/tools/run_batch.sh" "$TMP/collide.tsv" > "$RUN_OUT" 2>&1
status=$?
check "colliding run directories stop the batch" "$status" "1"
if grep -q "same run directory" "$RUN_OUT"; then ok "the collision is explained"
else bad "the collision is explained" "$(tail -3 "$RUN_OUT")"; fi
