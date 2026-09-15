# shellcheck shell=bash

MANIFEST_NAMES=()
MANIFEST_FLAGS=()

manifest_load() {
  local file="$1" line name flags lineno=0 dupes

  [[ -f "$file" ]] || { echo "manifest not found: $file" >&2; return 2; }

  MANIFEST_NAMES=()
  MANIFEST_FLAGS=()

  while IFS= read -r line || [[ -n "$line" ]]; do
    lineno=$((lineno + 1))
    line="${line%%#*}"
    [[ -z "${line// }" ]] && continue

    if [[ "$line" != *$'\t'* ]]; then
      echo "manifest line $lineno: expected 'name<TAB>flags'" >&2
      return 2
    fi

    name="$(echo "${line%%$'\t'*}" | xargs)"
    flags="${line#*$'\t'}"

    if [[ -z "$name" ]]; then
      echo "manifest line $lineno: empty name" >&2
      return 2
    fi
    if [[ "$flags" == *"--project-dir"* ]]; then
      echo "manifest line $lineno ($name): remove --project-dir, the runner sets it" >&2
      return 2
    fi

    MANIFEST_NAMES+=("$name")
    MANIFEST_FLAGS+=("$flags")
  done < "$file"

  if [[ ${#MANIFEST_NAMES[@]} -eq 0 ]]; then
    echo "manifest is empty: $file" >&2
    return 2
  fi

  dupes="$(printf '%s\n' "${MANIFEST_NAMES[@]}" | sort | uniq -d)"
  if [[ -n "$dupes" ]]; then
    echo "duplicate entry names:"$'\n'"$dupes" >&2
    return 2
  fi

  return 0
}

manifest_flag() {
  local flags="$1" wanted="$2" tok take=0
  for tok in $flags; do
    if [[ $take -eq 1 ]]; then
      [[ "$tok" == --* ]] && return 0
      printf '%s' "$tok"
      return 0
    fi
    [[ "$tok" == "$wanted" ]] && take=1
  done
  return 0
}

manifest_has_flag() {
  local flags="$1" wanted="$2" tok
  for tok in $flags; do
    [[ "$tok" == "$wanted" ]] && return 0
  done
  return 1
}
