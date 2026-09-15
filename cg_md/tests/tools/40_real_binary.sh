# shellcheck shell=bash
REAL_CG_MD="${REAL_CG_MD:-$ROOT/build/cg_md}"
if [[ ! -x "$REAL_CG_MD" ]]; then
  echo "[SKIP] --dry-run creates nothing (no build/cg_md)"
elif [[ ! -f "$ROOT/pdbs/AA/1-108_1.pdb" ]]; then
  echo "[SKIP] --dry-run creates nothing (no input PDB to name)"
else
  before="$(find "$ROOT/runs" -mindepth 1 2>/dev/null | wc -l | tr -d " ")"
  dry_rc=0
  for stage in all prepare adaptive relax em nvt npt production; do
    if ! "$REAL_CG_MD" --project-dir "$ROOT" --n-prot 3 --protomer 1-108_1 \
      --seq-length 108 --atoms-per-prot 236 --run-tag drynothing --no-elastic \
      --stage "$stage" --adaptive --dry-run >/dev/null 2>&1; then
      dry_rc=1
    fi
  done
  after="$(find "$ROOT/runs" -mindepth 1 2>/dev/null | wc -l | tr -d " ")"
  check "--dry-run succeeds for every stage" "$dry_rc" "0"
  if [[ "$before" == "$after" ]]; then
    ok "--dry-run creates nothing, for every stage"
  else
    bad "--dry-run creates nothing, for every stage" \
        "runs/ went from $before to $after entries: $(find "$ROOT/runs" -path "*drynothing*" | head -3 | tr "\n" " ")"
    rm -rf "$ROOT"/runs/*drynothing* 2>/dev/null
  fi
fi


if [[ ! -x "$REAL_CG_MD" ]]; then
  echo "[SKIP] run_all() skips prep when npt.gro/npt.cpt already exist (no build/cg_md)"
elif [[ ! -f "$ROOT/pdbs/AA/1-108_1.pdb" ]]; then
  echo "[SKIP] run_all() skips prep when npt.gro/npt.cpt already exist (no input PDB)"
else
  common=(--project-dir "$ROOT" --n-prot 3 --protomer 1-108_1 --seq-length 108
          --atoms-per-prot 236 --run-tag prepguard --adaptive --no-elastic --dry-run)

  out="$("$REAL_CG_MD" "${common[@]}" 2>&1)"
  rd="$(sed -n 's/^  run  *= \([^ ]*\).*/\1/p' <<< "$out" | head -1)"

  if [[ -z "$rd" ]]; then
    bad "run_all() skips prep when npt.gro/npt.cpt already exist" \
        "could not learn the run directory from --dry-run output"
  else
    rundir="$ROOT/runs/$rd"

    out_before="$("$REAL_CG_MD" "${common[@]}" 2>&1)"
    if grep -qE "'(martinize2|packmol)'" <<< "$out_before"; then
      ok "without npt.gro/npt.cpt, --dry-run still previews the prep pipeline"
    else
      bad "without npt.gro/npt.cpt, --dry-run still previews the prep pipeline" \
          "$(tail -5 <<< "$out_before")"
    fi

    mkdir -p "$rundir"
    : > "$rundir/npt.gro"
    : > "$rundir/npt.cpt"

    out_after="$("$REAL_CG_MD" "${common[@]}" 2>&1)"
    if grep -qE "'(martinize2|packmol)'" <<< "$out_after"; then
      bad "run_all() skips prep when npt.gro/npt.cpt already exist" \
          "prep commands still appear: $(grep -m1 -E "'(martinize2|packmol)'" <<< "$out_after")"
    else
      ok "run_all() skips prep when npt.gro/npt.cpt already exist"
    fi

    rm -rf "$rundir"
  fi
fi
