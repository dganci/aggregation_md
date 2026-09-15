# shellcheck shell=bash
logdir_ensure "$TMP/logs" >/dev/null 2>&1
check_true "logdir_ensure creates a usable directory" $?
check "logdir_ensure leaves no probe file behind" "$(ls -A "$TMP/logs" | wc -l | tr -d ' ')" "0"

logdir_ensure "$TMP/logs" >/dev/null 2>&1
check_true "logdir_ensure is idempotent" $?

mkdir -p "$TMP/readonly"
chmod a-w "$TMP/readonly"
if [[ "$(id -u)" == "0" ]]; then
  echo "[SKIP] logdir_ensure rejects an unwritable directory (running as root)"
else
  logdir_ensure "$TMP/readonly" >/dev/null 2>&1
  check_fail "logdir_ensure rejects an unwritable directory" $?
fi
chmod u+w "$TMP/readonly"

logdir_open_log 7 "$TMP/logs/a.log"
check_true "logdir_open_log opens a log for appending" $?
echo "first" >&7
logdir_close 7
logdir_open_log 7 "$TMP/logs/a.log"
echo "second" >&7
logdir_close 7
check "logdir_open_log appends rather than truncating" "$(wc -l < "$TMP/logs/a.log" | tr -d ' ')" "2"

logdir_open_log 7 "$TMP/no/such/dir/a.log" 2>/dev/null
check_fail "logdir_open_log fails on a missing directory" $?

logdir_append "$TMP/logs/rec.tsv" "a	b"
check_true "logdir_append writes a record" $?
logdir_append "$TMP/no/such/dir/rec.tsv" "x" 2>/dev/null
check_fail "logdir_append reports a lost record instead of dropping it" $?
