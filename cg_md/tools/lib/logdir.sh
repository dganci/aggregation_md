# shellcheck shell=bash

logdir_ensure() {
  local dir="$1" probe

  if ! mkdir -p "$dir" 2>/dev/null; then
    echo "cannot create log directory: $dir" >&2
    return 1
  fi

  probe="$dir/.probe.${BASHPID:-$$}"
  if ! printf 'ok\n' > "$probe" 2>/dev/null; then
    echo "log directory is not writable: $dir" >&2
    return 1
  fi
  if [[ "$(cat "$probe" 2>/dev/null)" != "ok" ]]; then
    echo "log directory did not read back what was written: $dir" >&2
    rm -f "$probe" 2>/dev/null
    return 1
  fi
  rm -f "$probe" 2>/dev/null
  return 0
}

logdir_open_log() {
  local fd="$1" log="$2"
  eval "exec ${fd}>>\"\$log\"" 2>/dev/null || {
    echo "cannot open log file for appending: $log" >&2
    return 1
  }
  return 0
}

logdir_close() {
  eval "exec ${1}>&-" 2>/dev/null || true
}

logdir_append() {
  local file="$1" line="$2"
  printf '%s\n' "$line" >> "$file" 2>/dev/null && return 0
  echo "cannot append to $file" >&2
  return 1
}
