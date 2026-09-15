# shellcheck shell=bash

declare -A MANIFEST_RUNDIR=()

rundir_resolve_all() {
  local cg_md="$1" project_dir="$2" ntomp="$3"
  local i name rd out bad=0 collide d n

  MANIFEST_RUNDIR=()

  for i in "${!MANIFEST_NAMES[@]}"; do
    name="${MANIFEST_NAMES[$i]}"
    local -a flags
    read -r -a flags <<< "${MANIFEST_FLAGS[$i]}"

    if out=$("$cg_md" --project-dir "$project_dir" --ntomp "$ntomp" "${flags[@]}" --dry-run 2>&1); then
      rd="$(sed -n 's/^  run  *= \([^ ]*\).*/\1/p' <<< "$out" | head -1)"
      MANIFEST_RUNDIR["$name"]="${rd:-unknown}"
      printf '  ok    %-22s -> runs/%s\n' "$name" "${rd:-unknown}"
    else
      printf '  FAIL  %s\n' "$name"
      sed 's/^/          /' <<< "$(tail -4 <<< "$out")"
      bad=$((bad + 1))
    fi
  done

  collide="$(printf '%s\n' "${MANIFEST_RUNDIR[@]}" | sort | uniq -d)"
  if [[ -n "$collide" ]]; then
    echo >&2
    echo "ERROR: these entries resolve to the same run directory, so they are the same" >&2
    echo "run and would resume each other. Change a parameter, or --run-tag them:" >&2
    for d in $collide; do
      for n in "${!MANIFEST_RUNDIR[@]}"; do
        [[ "${MANIFEST_RUNDIR[$n]}" == "$d" ]] && echo "  $d  <- $n" >&2
      done
    done
    return 1
  fi

  if [[ $bad -ne 0 ]]; then
    echo >&2
    echo "$bad entr(y/ies) failed configuration checks" >&2
    return 1
  fi
  return 0
}
