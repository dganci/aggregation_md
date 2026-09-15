# shellcheck shell=bash
cat > "$TMP/good.tsv" <<'EOF'
# a comment
#	an indented comment

alpha	--n-prot 5 --protomer p1 --relax --seed 1
beta	--n-prot 5 --protomer p2 --seed 2
EOF

manifest_load "$TMP/good.tsv"
check_true "manifest_load accepts a well-formed manifest" $?
check "manifest_load skips comments and blank lines" "${#MANIFEST_NAMES[@]}" "2"
check "manifest_load reads names"  "${MANIFEST_NAMES[0]}" "alpha"
check "manifest_load reads flags"  "${MANIFEST_FLAGS[1]}" "--n-prot 5 --protomer p2 --seed 2"

check "manifest_flag reads a flag value"        "$(manifest_flag "${MANIFEST_FLAGS[0]}" --protomer)" "p1"
check "manifest_flag on a missing flag is empty" "$(manifest_flag "${MANIFEST_FLAGS[0]}" --nope)" ""
check "manifest_flag does not treat a following flag as a value" \
      "$(manifest_flag "${MANIFEST_FLAGS[0]}" --relax)" ""
manifest_has_flag "${MANIFEST_FLAGS[0]}" --relax
check_true "manifest_has_flag finds a valueless flag" $?
manifest_has_flag "${MANIFEST_FLAGS[1]}" --relax
check_fail "manifest_has_flag rejects an absent flag" $?

printf 'alpha --no-tab-here\n' > "$TMP/notab.tsv"
manifest_load "$TMP/notab.tsv" 2>/dev/null
check_fail "manifest_load rejects a line without a tab" $?

printf 'alpha\t--a 1\nalpha\t--a 2\n' > "$TMP/dupe.tsv"
manifest_load "$TMP/dupe.tsv" 2>/dev/null
check_fail "manifest_load rejects duplicate entry names" $?

printf '\t--a 1\n' > "$TMP/noname.tsv"
manifest_load "$TMP/noname.tsv" 2>/dev/null
check_fail "manifest_load rejects an empty name" $?

printf 'alpha\t--project-dir /elsewhere\n' > "$TMP/projdir.tsv"
manifest_load "$TMP/projdir.tsv" 2>/dev/null
check_fail "manifest_load rejects an entry that overrides --project-dir" $?

printf '# only comments\n\n' > "$TMP/empty.tsv"
manifest_load "$TMP/empty.tsv" 2>/dev/null
check_fail "manifest_load rejects a manifest with no entries" $?

manifest_load "$TMP/does-not-exist.tsv" 2>/dev/null
check_fail "manifest_load rejects a missing file" $?
