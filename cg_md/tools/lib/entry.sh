# shellcheck shell=bash
fatal_in() {
  local log="$1"
  [[ -r "$log" ]] || return 1
  tail -c 2000000 "$log" 2>/dev/null |
    grep -qF -f <(printf '%s\n' "${FATAL_PATTERNS[@]}")
}

run_entry() {
  local i="$1"
  local name="${names[$i]}"
  local rd="${MANIFEST_RUNDIR[$name]}"
  local done_marker="$PROJECT_DIR/runs/$rd/.batch_done"
  local log="$LOG_ROOT/$name.log"
  local progress="[$((i+1))/${#names[@]}]"
  local outcome="$OUTCOME_DIR/$name"

  if [[ -f "$done_marker" ]]; then
    say "$progress $name: already done, skipping"
    echo skipped > "$outcome"; return 0
  fi

  local started t0 attempt code
  started="$(date -Is)"; t0="$(date +%s)"
  attempt=0; code=1
  while :; do
    attempt=$((attempt + 1))

    logdir_ensure "$LOG_ROOT" || {
      echo "FATAL: $LOG_ROOT became unusable mid-batch; aborting rather than" >&2
      echo "       running the remaining entries with nowhere to record them." >&2
      echo logdir > "$outcome"; return 3
    }

    logdir_open_log 7 "$log" || {
      echo "FATAL: cannot open $log; aborting instead of blaming the entry." >&2
      echo logdir > "$outcome"; return 3
    }

    say "$progress $name  attempt $attempt  -> $log"
    local flags
    read -r -a flags <<< "${flagsets[$i]}"
    # One GPU per entry, held for as long as the entry runs. The subshell owns
    # the lock file descriptor, so the device is released when the entry ends -
    # including when it is killed - without any cleanup path to forget.
    if [[ -n "${GPUS:-}" ]]; then
      gpu_run 7 "$name" \
        timeout --signal=TERM --kill-after=5m "$PER_SIM_TIMEOUT" \
          "$CG_MD" --project-dir "$PROJECT_DIR" --ntomp "$NTOMP" "${flags[@]}"
    else
      timeout --signal=TERM --kill-after=5m "$PER_SIM_TIMEOUT" \
          "$CG_MD" --project-dir "$PROJECT_DIR" --ntomp "$NTOMP" "${flags[@]}" \
          >&7 2>&7
    fi
    code=$?
    logdir_close 7

    [[ $code -eq 0 ]] && break
    if [[ ! -r "$log" ]]; then
      echo "    $name: exit $code, and $log is no longer readable: the log directory"
      echo "    was removed while this entry ran. Retrying with a fresh log."
    elif fatal_in "$log"; then
      echo "    $name: not retryable (configuration or input error)"
      break
    fi
    [[ $attempt -gt $RETRIES ]] && break
    echo "    $name: exit $code, retrying once: cg_md resumes from the last complete chunk"
    sleep 30
  done

  local t1 ended status
  t1="$(date +%s)"; ended="$(date -Is)"
  if [[ $code -eq 0 ]]; then
    status="done"
    if ! { mkdir -p "$(dirname "$done_marker")" && : > "$done_marker"; } 2>/dev/null; then
      echo "    WARNING cannot write $done_marker; a restart will re-run this entry" >&2
    fi
  elif [[ $code -eq 124 || $code -eq 137 ]]; then
    status=timeout
  else
    status=failed
  fi
  logdir_append "$SUMMARY" "$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "$name" "$rd" "$status" "$code" "$attempt" "$started" "$ended" \
    "$(( (t1 - t0) / 60 ))" "$log")"
  say "    $name: $status (exit $code) after $(( (t1 - t0) / 60 )) min"
  if [[ "$status" != "done" && -r "$log" ]]; then
    grep -m2 -E "^ERROR:|Fatal error|^Coarse-grained bead count" "$log" | sed "s/^/    | $name | /"
  fi
  echo "$status" > "$outcome"
}
